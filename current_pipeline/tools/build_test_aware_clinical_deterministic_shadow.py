"""Apply clinical deterministic guardrails to test-aware forwarded candidates.

The scorer consumes the high-recall forwarding set, multi-assay event timeline,
hospital evidence, image mentions, host context, and provisional phenotype
packets. It never sums reads across tests and never treats phenotype or imaging
as organism confirmation.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools.build_index_event_clinical_timeline_shadow import sha256_file
from tools.organism_taxonomy_classifier import classify_organism, inferred_taxonomic_rank
from tools.pathogen_normalization import canonical_key


DEFAULT_POLICY = Path("rules/test_aware_clinical_deterministic_shadow_v1.json")
SCHEMA_VERSION = "test_aware_clinical_deterministic_shadow.v1"
DECISION_LEVELS = {
    "picked_shadow": "Level 2",
    "review_high_priority": "Level 3",
    "review_context_needed": "Level 4",
}


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot identify patient from {path}")
    return int(match.group(1))


def candidate_event_refs(candidate: dict[str, Any]) -> list[str]:
    refs = {
        str(item.get("case_review_id") or "").strip()
        for item in candidate.get("case_decision_refs") or []
        if isinstance(item, dict) and str(item.get("case_review_id") or "").strip()
    }
    direct = str(candidate.get("case_review_id") or "").strip()
    if direct:
        refs.add(direct)
    return sorted(refs)


def relevant_events(timeline: dict[str, Any], refs: list[str]) -> list[dict[str, Any]]:
    ref_set = set(refs)
    events = timeline.get("index_events") or []
    if not ref_set:
        return events
    return [event for event in events
            if str(event.get("case_review_id") or event.get("specimen_code") or "") in ref_set]


def image_mentions(
    timeline: dict[str, Any], candidate: dict[str, Any], refs: list[str], policy: dict[str, Any]
) -> list[dict[str, Any]]:
    key = canonical_key(candidate.get("organism_name"))
    statuses = set(policy["episode_evidence_statuses"])
    matches: list[dict[str, Any]] = []
    for event in relevant_events(timeline, refs):
        event_id = str(event.get("event_id") or "")
        for row in timeline.get("clinical_evidence") or []:
            if row.get("evidence_role") != "imaging_context":
                continue
            if (row.get("episode_evidence_status") or {}).get(event_id) not in statuses:
                continue
            for mention in row.get("diagnostic_mentions") or []:
                if canonical_key(mention.get("concept_name")) != key:
                    continue
                matches.append({
                    "evidence_id": row.get("evidence_id"),
                    "event_id": event_id,
                    "case_review_id": event.get("case_review_id"),
                    "concept_name": mention.get("concept_name"),
                    "assertion": mention.get("assertion"),
                    "source_text": row.get("source_text"),
                    "observed_at": row.get("observed_at"),
                    "reported_at": row.get("reported_at"),
                    "available_at_index": (row.get("available_at_index") or {}).get(event_id),
                    "episode_evidence_status": (row.get("episode_evidence_status") or {}).get(event_id),
                    "microbiologic_confirmation": False,
                    "source_ref": row.get("source_ref"),
                })
    unique: dict[tuple[Any, Any, Any], dict[str, Any]] = {}
    for item in matches:
        unique[(item.get("evidence_id"), item.get("event_id"), item.get("concept_name"))] = item
    return list(unique.values())


def host_evidence(timeline: dict[str, Any], refs: list[str]) -> list[dict[str, Any]]:
    output = []
    for event in relevant_events(timeline, refs):
        event_id = str(event.get("event_id") or "")
        for row in timeline.get("clinical_evidence") or []:
            if row.get("evidence_role") not in {"host_history", "host_medication"}:
                continue
            output.append({
                "evidence_id": row.get("evidence_id"),
                "event_id": event_id,
                "case_review_id": event.get("case_review_id"),
                "evidence_role": row.get("evidence_role"),
                "label": row.get("label"),
                "assertion": row.get("assertion"),
                "source_text": row.get("source_text"),
                "observed_at": row.get("observed_at"),
                "reported_at": row.get("reported_at"),
                "episode_evidence_status": (row.get("episode_evidence_status") or {}).get(event_id),
                "relative_to_treatment": row.get("relative_to_treatment"),
                "source_ref": row.get("source_ref"),
            })
    unique: dict[tuple[Any, Any], dict[str, Any]] = {}
    for item in output:
        unique[(item.get("evidence_id"), item.get("event_id"))] = item
    return list(unique.values())


def phenotype_context(
    packet: dict[str, Any], organism_name: Any, refs: list[str]
) -> dict[str, Any]:
    key = canonical_key(organism_name)
    ref_set = set(refs)
    contexts = []
    exact_micro = []
    event_ids = []
    for event in packet.get("events") or []:
        index = event.get("index_event") or {}
        case = str(index.get("case_review_id") or index.get("specimen_code") or "")
        if ref_set and case not in ref_set:
            continue
        for candidate in event.get("candidate_packets") or []:
            if str(candidate.get("canonical_key") or "") != key:
                continue
            event_ids.append(index.get("event_id"))
            contexts.extend(candidate.get("phenotype_context") or [])
            exact_micro.extend(candidate.get("exact_organism_workbook_microbiology") or [])
    return {
        "phenotype_status": packet.get("phenotype_status"),
        "production_link_verified": False,
        "matched_event_ids": event_ids,
        "routed_phenotype_context": contexts,
        "exact_organism_workbook_microbiology": exact_micro,
        "deterministic_effect": "none_provisional_link",
        "constraint": "phenotype is derived context, not independent organism proof",
    }


def direct_level(candidate: dict[str, Any]) -> int | None:
    if candidate.get("evidence_source") == "hospital_only":
        value = (candidate.get("hospital_evidence_detail") or {}).get("best_hospital_level")
    else:
        value = (candidate.get("hospital_profile") or {}).get("best_direct_level")
    match = re.search(r"(\d+)", str(value or ""))
    return int(match.group(1)) if match else None


def analytically_credible(candidate: dict[str, Any]) -> bool:
    analytical = candidate.get("analytical_profile") or {}
    rank = analytical.get("best_rank_in_retained_universe")
    try:
        rank = int(rank)
    except (TypeError, ValueError):
        return False
    return rank <= 3 or (bool(analytical.get("reproducibility_axis")) and rank <= 10)


def generic_identity(name: Any, taxonomy: dict[str, Any], policy: dict[str, Any]) -> bool:
    text = " ".join(str(name or "").casefold().split())
    if inferred_taxonomic_rank(name) == "genus_or_group_label":
        return True
    return any(term in text for term in policy["generic_identity_terms"])


def decide(
    candidate: dict[str, Any], taxonomy: dict[str, Any], images: list[dict[str, Any]],
    phenotype: dict[str, Any], policy: dict[str, Any],
) -> tuple[str, list[str], list[str]]:
    incoming = str(candidate.get("decision") or "review_context_needed")
    family = str(candidate.get("taxonomy_family") or taxonomy.get("primary_rule_family") or
                 "unmapped_or_uncertain")
    source = "hospital_only" if candidate.get("evidence_source") == "hospital_only" else "mngs_positive_reads"
    rules = ["CLIN-S1-PRESERVE-HIGH-RECALL-FORWARDING", "CLIN-S1-NO-CROSS-TEST-READ-SUM"]
    reasons: list[str] = []

    if incoming == "picked_shadow":
        if source == "hospital_only" and generic_identity(candidate.get("organism_name"), taxonomy, policy):
            rules.append("CLIN-G1-GENERIC-HOSPITAL-LABEL-NOT-SPECIES-PICK")
            reasons.append("Hospital evidence is retained, but a group/complex/spp./pending-identification label is not a formal species-level Picked result.")
            return "review_high_priority", rules, reasons
        rules.append("CLIN-P1-PRESERVE-EXISTING-PICK")
        reasons.append("The answer-blind analytical/hospital scorer already met its conservative Picked rule and no clinical guardrail invalidated it.")
        return "picked_shadow", rules, reasons

    if incoming == "review_high_priority":
        rules.append("CLIN-H1-PRESERVE-HIGH-REVIEW")
        reasons.append("Existing high-priority review status is preserved for clinical adjudication.")
        return "review_high_priority", rules, reasons

    credible = analytically_credible(candidate)
    level = direct_level(candidate)
    if source == "mngs_positive_reads" and credible and level is not None and level <= 3:
        rules.append("CLIN-H2-EXACT-HOSPITAL-SUPPORT")
        reasons.append("An analytically credible mNGS signal has exact/alias hospital-side Level 1-3 support; it is promoted for high-priority review, not automatically Picked.")
        return "review_high_priority", rules, reasons

    opportunistic = family in set(policy["opportunistic_families"])
    if source == "mngs_positive_reads" and credible and opportunistic and images:
        rules.extend(["CLIN-I1-EXACT-IMAGE-MENTION", "CLIN-I2-IMAGE-CONTEXT-NOT-CONFIRMATION"])
        reasons.append("The same infection episode includes an exact organism/disease imaging mention; imaging raises review priority but is not microbiological confirmation.")
        return "review_high_priority", rules, reasons
    if source == "mngs_positive_reads" and credible and opportunistic and bool(candidate.get("host_support")):
        rules.append("CLIN-HOST1-OPPORTUNISTIC-HOST-RISK")
        reasons.append("Host-risk context plus an analytically credible opportunistic signal warrants high-priority review but does not prove causation.")
        return "review_high_priority", rules, reasons

    if phenotype.get("routed_phenotype_context"):
        rules.append("CLIN-PHENO1-PROVISIONAL-CONTEXT-NO-TIER-EFFECT")
        reasons.append("Relevant phenotype evidence is retained for later review; provisional linkage prevents deterministic tier changes.")
    rules.append("CLIN-C1-RETAIN-CONTEXT")
    reasons.append("The candidate remains visible for clinical review, but current clinical evidence is insufficient for high-priority or Picked status.")
    return "review_context_needed", rules, reasons


def build_candidate(
    candidate: dict[str, Any], timeline: dict[str, Any], phenotype_packet: dict[str, Any],
    policy: dict[str, Any],
) -> dict[str, Any]:
    original_decision = str(candidate.get("decision") or "review_context_needed")
    routed_decision = str(candidate.get("history_routed_decision") or original_decision)
    routed_candidate = dict(candidate)
    routed_candidate["decision"] = routed_decision
    refs = candidate_event_refs(candidate)
    taxonomy = classify_organism(candidate.get("organism_name"),
                                 biological_class=candidate.get("category"))
    images = image_mentions(timeline, candidate, refs, policy)
    host = host_evidence(timeline, refs)
    phenotype = phenotype_context(phenotype_packet, candidate.get("organism_name"), refs)
    decision, rule_ids, reasons = decide(routed_candidate, taxonomy, images, phenotype, policy)
    if routed_decision != original_decision:
        rule_ids = ["CLIN-HIST1-HISTORY-ROUTE-TO-SCORER", *rule_ids]
        reasons = [
            "The versioned history route changed the scorer entry tier by one bounded level; "
            "consumed history remains provenance and cannot independently promote High to Picked.",
            *reasons,
        ]
    events = relevant_events(timeline, refs)
    return {
        "patient_id": str(candidate.get("patient_id") or timeline.get("patient_id")),
        "organism_name": candidate.get("organism_name"),
        "organism_key": canonical_key(candidate.get("organism_name")),
        "candidate_source": (
            "hospital_only" if candidate.get("evidence_source") == "hospital_only"
            else "mngs_positive_reads"
        ),
        "candidate_event_refs": refs,
        "linked_events": [{
            "event_id": event.get("event_id"),
            "case_review_id": event.get("case_review_id"),
            "specimen_site": event.get("specimen_site"),
            "collected_time": event.get("collected_time"),
            "seq_ids": event.get("seq_ids") or [],
        } for event in events],
        "taxonomy_profile": taxonomy,
        "incoming_decision": original_decision,
        "history_routed_decision": routed_decision,
        "incoming_level": candidate.get("integrated_level"),
        "analytic_route": candidate.get("analytic_route"),
        "history_adjusted_route": candidate.get("history_adjusted_route"),
        "history_route_changed": bool(candidate.get("history_route_changed")),
        "history_policy_family": candidate.get("history_policy_family"),
        "history_rule_ids": candidate.get("history_rule_ids") or [],
        "history_evidence_ids": candidate.get("history_evidence_ids") or [],
        "consumed_evidence_ids": candidate.get("consumed_evidence_ids") or [],
        "history_evidence_consumption": candidate.get("history_evidence_consumption") or [],
        "history_adjustment_blockers": candidate.get("history_adjustment_blockers") or [],
        "history_source_disagreement": candidate.get("history_source_disagreement"),
        "future_picked_gate_review_required": bool(
            candidate.get("future_picked_gate_review_required")
        ),
        "future_picked_gate_review_status": candidate.get(
            "future_picked_gate_review_status"
        ),
        "future_picked_gate_review_requirements": candidate.get(
            "future_picked_gate_review_requirements"
        ) or [],
        "clinical_decision": decision,
        "clinical_level": DECISION_LEVELS[decision],
        "formal_pick_allowed": decision == "picked_shadow",
        "rule_ids": rule_ids,
        "reasons": reasons,
        "analytical_profile": candidate.get("analytical_profile"),
        "hospital_profile": candidate.get("hospital_profile"),
        "hospital_evidence_detail": candidate.get("hospital_evidence_detail"),
        "host_support_flag_from_upstream": bool(candidate.get("host_support")),
        "host_evidence": host,
        "image_evidence": images,
        "image_constraint": "context_only_not_microbiologic_confirmation",
        "phenotype_evidence": phenotype,
        "illness_lab_and_treatment_response_effect": "none_for_organism_attribution",
        "source_candidate": candidate,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "patient_id", "organism_name", "candidate_source", "candidate_event_refs",
        "taxonomy_family", "incoming_decision", "history_routed_decision",
        "analytic_route", "history_adjusted_route", "history_route_changed",
        "history_policy_family",
        "clinical_decision", "clinical_level",
        "formal_pick_allowed", "best_rank", "reproducibility_axis", "direct_hospital_level",
        "host_support", "image_match_count", "phenotype_context_count",
        "consumed_evidence_ids", "future_picked_gate_review_required",
        "future_picked_gate_review_status", "rule_ids",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in rows:
            analytical = item.get("analytical_profile") or {}
            phenotype = item.get("phenotype_evidence") or {}
            source = item.get("source_candidate") or {}
            writer.writerow({
                "patient_id": item["patient_id"],
                "organism_name": item["organism_name"],
                "candidate_source": item["candidate_source"],
                "candidate_event_refs": "|".join(item["candidate_event_refs"]),
                "taxonomy_family": item["taxonomy_profile"].get("primary_rule_family"),
                "incoming_decision": item["incoming_decision"],
                "history_routed_decision": item["history_routed_decision"],
                "analytic_route": item.get("analytic_route"),
                "history_adjusted_route": item.get("history_adjusted_route"),
                "history_route_changed": item.get("history_route_changed"),
                "history_policy_family": item.get("history_policy_family"),
                "clinical_decision": item["clinical_decision"],
                "clinical_level": item["clinical_level"],
                "formal_pick_allowed": item["formal_pick_allowed"],
                "best_rank": analytical.get("best_rank_in_retained_universe"),
                "reproducibility_axis": analytical.get("reproducibility_axis"),
                "direct_hospital_level": direct_level(source),
                "host_support": item["host_support_flag_from_upstream"],
                "image_match_count": len(item["image_evidence"]),
                "phenotype_context_count": len(phenotype.get("routed_phenotype_context") or []),
                "consumed_evidence_ids": "|".join(item.get("consumed_evidence_ids") or []),
                "future_picked_gate_review_required": item.get(
                    "future_picked_gate_review_required"
                ),
                "future_picked_gate_review_status": item.get(
                    "future_picked_gate_review_status"
                ),
                "rule_ids": "|".join(item["rule_ids"]),
            })


def run(
    scorer_root: Path, timeline_root: Path, phenotype_packet_root: Path,
    output_dir: Path, *, policy_path: Path = DEFAULT_POLICY,
    history_route_root: Path | None = None,
    patients: set[int] | None = None,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    policy = read_json(policy_path)
    paths = sorted((scorer_root / "patient_outputs").glob(
        "NGS_patient_*_test_aware_deterministic_shadow.json"), key=patient_number)
    if patients is not None:
        paths = [path for path in paths if patient_number(path) in patients]
    if not paths:
        raise ValueError("No test-aware patient outputs found")
    output_dir.mkdir(parents=True, exist_ok=True)
    patient_dir = output_dir / "patient_outputs"
    patient_dir.mkdir()
    all_rows = []
    expected_forwarded_count = 0
    counts: Counter[str] = Counter()
    for source_path in paths:
        patient = patient_number(source_path)
        timeline_path = timeline_root / f"NGS_patient_{patient}_clinical_timeline_shadow.json"
        phenotype_path = phenotype_packet_root / f"patient_{patient}_phenotype_event_shadow.json"
        if not timeline_path.is_file() or not phenotype_path.is_file():
            raise FileNotFoundError(timeline_path if not timeline_path.is_file() else phenotype_path)
        source = read_json(source_path)
        history_path = None
        clinical_input = source
        if history_route_root is not None:
            history_path = (
                history_route_root / "patient_outputs"
                / f"NGS_patient_{patient}_history_adjusted_route_shadow.json"
            )
            if not history_path.is_file():
                raise FileNotFoundError(history_path)
            clinical_input = read_json(history_path)
        timeline = read_json(timeline_path)
        phenotype_packet = read_json(phenotype_path)
        forwarded = clinical_input.get("clinical_scorer_forward") or []
        expected_forwarded_count += len(forwarded)
        candidates = [build_candidate(item, timeline, phenotype_packet, policy)
                      for item in forwarded]
        if len(candidates) != len(forwarded):
            raise ValueError(f"Candidate loss in P{patient}")
        all_rows.extend(candidates)
        for item in candidates:
            counts[f"incoming:{item['incoming_decision']}"] += 1
            counts[f"history_routed:{item['history_routed_decision']}"] += 1
            counts[f"clinical:{item['clinical_decision']}"] += 1
            if item["history_route_changed"]:
                counts[
                    f"history_route:{item['analytic_route']}->{item['history_adjusted_route']}"
                ] += 1
            if item["incoming_decision"] != item["clinical_decision"]:
                counts[f"transition:{item['incoming_decision']}->{item['clinical_decision']}"] += 1
        payload = {
            "schema_version": SCHEMA_VERSION,
            "policy_id": policy.get("policy_id"),
            "answer_blind": True,
            "patient_id": str(patient),
            "candidate_count": len(candidates),
            "picked_shadow": [x for x in candidates if x["clinical_decision"] == "picked_shadow"],
            "review_high_priority": [x for x in candidates if x["clinical_decision"] == "review_high_priority"],
            "review_context_needed": [x for x in candidates if x["clinical_decision"] == "review_context_needed"],
            "all_forwarded_candidates": candidates,
            "no_forced_pick": not any(x["clinical_decision"] == "picked_shadow" for x in candidates),
            "source": {
                "test_aware_file": str(source_path.resolve()),
                "test_aware_sha256": sha256_file(source_path),
                "history_route_file": str(history_path.resolve()) if history_path else None,
                "history_route_sha256": sha256_file(history_path) if history_path else None,
                "timeline_file": str(timeline_path.resolve()),
                "timeline_sha256": sha256_file(timeline_path),
                "phenotype_packet_file": str(phenotype_path.resolve()),
                "phenotype_packet_sha256": sha256_file(phenotype_path),
            },
            "constraints": policy.get("constraints") or [],
        }
        write_json(patient_dir / f"NGS_patient_{patient}_test_aware_clinical_shadow.json", payload)
        counts["patients"] += 1

    if len(all_rows) != expected_forwarded_count:
        raise ValueError(
            f"Expected {expected_forwarded_count} forwarded candidates from the input files, "
            f"got {len(all_rows)}"
        )
    write_csv(output_dir / "clinical_decisions.csv", all_rows)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "policy_id": policy.get("policy_id"),
        "answer_blind": True,
        "candidate_count": len(all_rows),
        "expected_forwarded_count_from_inputs": expected_forwarded_count,
        "patient_count": len(paths),
        "counts": dict(sorted(counts.items())),
        "patients_without_picked": sorted({
            int(item["patient_id"]) for item in all_rows
            if not any(other["patient_id"] == item["patient_id"] and
                       other["clinical_decision"] == "picked_shadow" for other in all_rows)
        }),
        "source_scorer_root": str(scorer_root.resolve()),
        "history_route_root": str(history_route_root.resolve()) if history_route_root else None,
        "selected_patients": sorted(patients) if patients is not None else None,
        "timeline_root": str(timeline_root.resolve()),
        "phenotype_packet_root": str(phenotype_packet_root.resolve()),
        "policy_file": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "constraints": policy.get("constraints") or [],
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scorer_root", type=Path)
    parser.add_argument("timeline_root", type=Path)
    parser.add_argument("phenotype_packet_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--history-route-root", type=Path)
    parser.add_argument("--patients", type=int, nargs="*")
    args = parser.parse_args()
    selected = set(args.patients) if args.patients else None
    print(json.dumps(run(
        args.scorer_root, args.timeline_root, args.phenotype_packet_root,
        args.output_dir, policy_path=args.policy,
        history_route_root=args.history_route_root,
        patients=selected,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

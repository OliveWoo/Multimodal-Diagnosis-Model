"""Select answer-blind OBER cases that warrant a second Sol review."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

from tools.adjudicate_rag_v2_packets import validate_adjudication


NON_VISIBLE_ACTIONS = {"reject", "insufficient_evidence"}
FUNGAL_OR_VIRAL = {"fungal", "fungus", "viral", "virus"}
BACKGROUND_RULE_FAMILIES = {
    "commensal_virome_or_endogenous_element",
    "environmental_background_or_reagent",
}


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip()
    ]


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows),
        encoding="utf-8",
    )


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def normalize(value: Any) -> str:
    return str(value or "").strip().lower()


def escalation_reasons(packet: dict[str, Any], record: dict[str, Any]) -> list[str]:
    adjudication = as_dict(record.get("adjudication"))
    decision = as_dict(adjudication.get("decision"))
    knowledge = as_dict(adjudication.get("knowledge_card"))
    boundary = as_dict(record.get("boundary_guardrail"))
    snapshot = as_dict(boundary.get("evidence_snapshot"))
    exact = as_dict(snapshot.get("exact_species_evidence"))
    organism = as_dict(packet.get("organism"))
    taxonomy = as_dict(organism.get("taxonomy_profile"))
    local = as_dict(packet.get("local_evidence"))
    mngs = as_dict(local.get("mngs_evidence"))

    action = normalize(decision.get("recommended_action"))
    visibility = normalize(decision.get("clinician_visibility"))
    confidence = normalize(decision.get("confidence"))
    pathogenicity = normalize(knowledge.get("pulmonary_pathogenicity"))
    colonization = normalize(knowledge.get("respiratory_colonization_risk"))
    biological_class = normalize(
        taxonomy.get("biological_class") or organism.get("classification")
    )
    primary_family = normalize(taxonomy.get("primary_rule_family"))
    reads_tier = normalize(snapshot.get("reads_tier") or mngs.get("reads_tier"))
    signal_tier = normalize(mngs.get("mngs_signal_tier"))
    dominance = normalize(snapshot.get("dominance_tier") or mngs.get("dominance_tier"))
    strong_signal = (
        reads_tier in {"r2_moderate", "r3_high"}
        or signal_tier == "m1_strong"
        or dominance in {"d2_moderate", "d3_dominant"}
        or bool(exact.get("direct"))
    )

    reasons: list[str] = []
    if bool(boundary.get("clinician_visibility_allowed")) and visibility == "do_not_show":
        reasons.append("boundary_allows_visibility_but_luna_does_not_show")
    if bool(exact.get("direct")) and action in NON_VISIBLE_ACTIONS:
        reasons.append("direct_same_organism_evidence_but_non_visible_action")
    if pathogenicity == "established" and action == "reject":
        reasons.append("established_pathogen_rejected")
    if confidence == "low":
        reasons.append("low_confidence")
    if action in {"primary_pathogen", "secondary_pathogen"} and confidence == "low":
        reasons.append("positive_action_low_confidence_mismatch")

    taxonomy_uncertain = bool(organism.get("taxonomy_unknown")) or normalize(
        taxonomy.get("mapping_status")
    ) in {"provisional_unmapped", "unmapped", "unresolved", "ambiguous"}
    if taxonomy_uncertain and action != "reject":
        reasons.append("taxonomy_uncertainty_affects_non_reject_decision")

    if bool(exact.get("direct")) and colonization in {"high", "moderate"}:
        reasons.append("direct_evidence_with_colonization_conflict")

    if (
        biological_class in FUNGAL_OR_VIRAL
        and primary_family not in BACKGROUND_RULE_FAMILIES
        and strong_signal
        and action in NON_VISIBLE_ACTIONS
        and pathogenicity in {"uncertain", "recognized_but_uncommon", "rare_case_report_only"}
    ):
        reasons.append("strong_ambiguous_fungal_or_viral_candidate")
    return sorted(set(reasons))


def select(packets_path: Path, luna_dir: Path, output_dir: Path) -> dict[str, Any]:
    packets = read_jsonl(packets_path)
    packet_by_id = {str(packet.get("case_id") or ""): packet for packet in packets}
    if len(packet_by_id) != len(packets) or "" in packet_by_id:
        raise ValueError("Packets contain missing or duplicate case_id values")

    manifest = read_json(luna_dir / "run_manifest.json")
    records = read_jsonl(luna_dir / "adjudications.jsonl")
    record_by_id = {str(record.get("case_id") or ""): record for record in records}
    if set(record_by_id) != set(packet_by_id):
        raise ValueError(
            "Luna coverage must be complete before Sol escalation selection: "
            f"missing={sorted(set(packet_by_id) - set(record_by_id))}, "
            f"extra={sorted(set(record_by_id) - set(packet_by_id))}"
        )
    if int(manifest.get("completed_count") or 0) != len(packets):
        raise ValueError("Luna run manifest is incomplete")
    if int(manifest.get("failure_count") or 0) != 0:
        raise ValueError("Luna run manifest contains failures")
    if manifest.get("answer_source_read") is not False:
        raise ValueError("Luna run does not explicitly record answer_source_read=false")

    selected: list[dict[str, Any]] = []
    reason_counts: Counter[str] = Counter()
    for case_id in sorted(packet_by_id):
        packet = packet_by_id[case_id]
        record = record_by_id[case_id]
        if record.get("answer_source_read") is not False:
            raise ValueError(f"{case_id} does not record answer_source_read=false")
        validate_adjudication(as_dict(record.get("adjudication")), packet)
        reasons = escalation_reasons(packet, record)
        if not reasons:
            continue
        reason_counts.update(reasons)
        organism = as_dict(packet.get("organism"))
        decision = as_dict(as_dict(record.get("adjudication")).get("decision"))
        selected.append(
            {
                "case_id": case_id,
                "patient_id": packet.get("patient_id"),
                "organism_name": organism.get("display_name"),
                "luna_action": decision.get("recommended_action"),
                "luna_visibility": decision.get("clinician_visibility"),
                "luna_confidence": decision.get("confidence"),
                "escalation_reasons": reasons,
                "answer_source_read": False,
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    write_jsonl(output_dir / "sol_escalation_cases.jsonl", selected)
    (output_dir / "case_ids.txt").write_text(
        "".join(row["case_id"] + "\n" for row in selected), encoding="utf-8"
    )
    summary = {
        "schema_version": "ober_sol_escalation_selection_v1",
        "answer_blind": True,
        "ground_truth_read": False,
        "packet_count": len(packets),
        "luna_completed_count": len(records),
        "selected_count": len(selected),
        "selected_case_ids": [row["case_id"] for row in selected],
        "reason_counts": dict(sorted(reason_counts.items())),
        "policy": {
            "boundary_allowed_but_hidden": True,
            "direct_evidence_conflict": True,
            "established_pathogen_rejected": True,
            "low_confidence": True,
            "taxonomy_uncertainty": True,
            "strong_ambiguous_fungal_or_viral_candidate": True,
        },
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packets", type=Path, required=True)
    parser.add_argument("--luna-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            select(args.packets, args.luna_dir, args.output_dir),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

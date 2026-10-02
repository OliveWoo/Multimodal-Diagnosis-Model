"""Apply family-specific high-to-Picked promotion gates to clinical shadow v1.

The input clinical shadow remains immutable. This stage only asks whether an
existing review_high_priority candidate has enough independent axes for a
formal shadow Picked result. Candidates that fail a gate remain high priority.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from tools.build_index_event_clinical_timeline_shadow import sha256_file
from tools.pathogen_normalization import canonical_key


DEFAULT_POLICY = Path("rules/test_aware_clinical_promotion_v3.json")
SCHEMA_VERSION = "test_aware_clinical_deterministic_shadow.v4"
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


def as_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def direct_hospital_level(candidate: dict[str, Any]) -> int | None:
    if candidate.get("candidate_source") == "hospital_only":
        value = (candidate.get("hospital_evidence_detail") or {}).get("best_hospital_level")
    else:
        value = (candidate.get("hospital_profile") or {}).get("best_direct_level")
    match = re.search(r"(\d+)", str(value or ""))
    return int(match.group(1)) if match else None


def precise_identity(candidate: dict[str, Any], policy: dict[str, Any]) -> bool:
    taxonomy = candidate.get("taxonomy_profile") or {}
    text = " ".join(str(candidate.get("organism_name") or "").casefold().split())
    rank = str(taxonomy.get("taxonomic_rank") or "")
    mapping = str(taxonomy.get("mapping_status") or "")
    if rank not in {"species", "species_candidate"} or mapping == "unmapped":
        return False
    return not any(term in text for term in policy["generic_identity_terms"])


def specimen_context(candidate: dict[str, Any]) -> str:
    source = candidate.get("source_candidate") or {}
    explicit = str(source.get("specimen_context") or "").strip()
    if explicit:
        return explicit
    sites = " ".join(
        str(event.get("specimen_site") or "").casefold()
        for event in candidate.get("linked_events") or []
    )
    respiratory_terms = ("balf", "bronch", "sputum", "trache", "respiratory", "lung")
    sterile_terms = ("blood", "csf", "pleural", "tissue", "sterile")
    if any(term in sites for term in respiratory_terms):
        return "lower_respiratory"
    if any(term in sites for term in sterile_terms):
        return "sterile_or_systemic"
    return "other_or_unknown"


def exact_image_context(candidate: dict[str, Any]) -> bool:
    allowed = {"support", "supports", "supported", "compatible", "differential", "suspected"}
    return any(
        str(item.get("assertion") or "").casefold() in allowed
        for item in candidate.get("image_evidence") or []
    )


def verified_host_support_profile(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any]:
    """Separate a broad host-risk flag from evidence eligible for Picked promotion."""
    accepted_assertions = set(policy.get("pjp_verified_host_assertions") or [
        "verified_active_exposure",
        "verified_high_risk_host",
        "verified_immunosuppression",
    ])
    accepted_rows = []
    rejected_rows = []
    for row in candidate.get("host_evidence") or []:
        assertion = str(row.get("assertion") or "").strip()
        summary = {
            "evidence_id": row.get("evidence_id"),
            "evidence_role": row.get("evidence_role"),
            "label": row.get("label"),
            "assertion": assertion,
            "episode_evidence_status": row.get("episode_evidence_status"),
        }
        if assertion in accepted_assertions:
            accepted_rows.append(summary)
        else:
            rejected_rows.append(summary)

    phenotype = candidate.get("phenotype_evidence") or {}
    phenotype_verified = bool(phenotype.get("production_link_verified"))
    verified_phenotypes = []
    if phenotype_verified:
        for row in phenotype.get("routed_phenotype_context") or []:
            if (str(row.get("status") or "").upper() == "YES"
                    and str(row.get("confidence") or "").upper() == "HIGH"
                    and str(row.get("phenotype") or "") in {
                        "IMMUNOCOMPROMISED", "FUNGAL_INFECTION_RISK"
                    }):
                verified_phenotypes.append({
                    "source_id": row.get("source_id"),
                    "phenotype": row.get("phenotype"),
                    "reason": row.get("reason"),
                })
    return {
        "eligible": bool(accepted_rows or verified_phenotypes),
        "accepted_host_evidence": accepted_rows,
        "nonqualifying_host_evidence": rejected_rows,
        "phenotype_production_link_verified": phenotype_verified,
        "verified_phenotype_evidence": verified_phenotypes,
        "broad_host_flag_alone_is_insufficient": True,
    }


def _raw_key(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").casefold())


def _recursive_text(value: Any) -> str:
    if isinstance(value, dict):
        return " ".join(_recursive_text(item) for item in value.values())
    if isinstance(value, list):
        return " ".join(_recursive_text(item) for item in value)
    return str(value or "").casefold()


def _parse_datetime(value: Any) -> datetime | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return datetime.fromisoformat(text.replace("Z", "+00:00")).replace(tzinfo=None)
    except ValueError:
        return None


def candidate_event_times(candidate: dict[str, Any]) -> list[datetime]:
    values = [(candidate.get("source_candidate") or {}).get("collected_time")]
    values.extend(event.get("collected_time") for event in candidate.get("linked_events") or [])
    parsed = [_parse_datetime(value) for value in values]
    return list(dict.fromkeys(value for value in parsed if value is not None))


def candida_hospital_evidence_detail(candidate: dict[str, Any]) -> dict[str, Any]:
    if candidate.get("candidate_source") == "hospital_only":
        return candidate.get("hospital_evidence_detail") or {}
    return (candidate.get("hospital_profile") or {}).get("hospital_evidence_detail") or {}


def _module_evidence_rows(detail: dict[str, Any]) -> list[dict[str, Any]]:
    source_rows = [
        row for row in detail.get("source_observations") or [] if isinstance(row, dict)
    ]
    if source_rows:
        return source_rows
    output: list[dict[str, Any]] = []
    module_evidence = detail.get("module_evidence") or {}
    if isinstance(module_evidence, dict):
        for rows in module_evidence.values():
            output.extend(row for row in rows or [] if isinstance(row, dict))
    for key, rows in detail.items():
        if key in {
            "source_observations", "module_evidence", "module_level_summary",
            "source_observation_resolution",
        }:
            continue
        if isinstance(rows, list):
            output.extend(row for row in rows if isinstance(row, dict))
    return output


def candida_positive_assertion(row: dict[str, Any], policy: dict[str, Any]) -> str:
    config = policy["candida_invasive_evidence"]
    text = f" {_recursive_text(row)} "
    if any(term.casefold() in text for term in config["explicit_negative_terms"]):
        return "negative"
    status = str(row.get("quantitation_status") or "").casefold()
    if status.startswith("positive"):
        return "positive"
    if any(term.casefold() in text for term in config["explicit_positive_terms"]):
        return "positive"
    return "unknown"


def candida_invasive_specimen(row: dict[str, Any], policy: dict[str, Any]) -> bool:
    config = policy["candida_invasive_evidence"]
    category = _raw_key(row.get("specimen_category"))
    if category in {_raw_key(value) for value in config["sterile_specimen_categories"]}:
        return True
    text = " ".join(
        str(row.get(key) or "").casefold()
        for key in ("specimen_type", "specimen_category", "sample", "site")
    )
    return any(term.casefold() in text for term in config["sterile_specimen_terms"])


def candida_group_label(candidate: dict[str, Any], policy: dict[str, Any]) -> bool:
    name = str(candidate.get("organism_name") or "").casefold()
    return any(term.casefold() in name for term in policy["candida_group_identity_terms"])


def candida_evidence_profile(candidate: dict[str, Any], policy: dict[str, Any]) -> dict[str, Any]:
    detail = candida_hospital_evidence_detail(candidate)
    rows = _module_evidence_rows(detail)
    event_times = candidate_event_times(candidate)
    group_identity = candida_group_label(candidate, policy)
    window_hours = float(policy["candida_invasive_evidence"]["event_window_hours"])
    evaluated_rows: list[dict[str, Any]] = []
    for row in rows:
        assertion = candida_positive_assertion(row, policy)
        invasive_specimen = candida_invasive_specimen(row, policy)
        observed_at = _parse_datetime(row.get("collected_time"))
        if observed_at is None or not event_times:
            timing = "unknown"
        elif any(abs((observed_at - event_time).total_seconds()) <= window_hours * 3600
                 for event_time in event_times):
            timing = "within_event_window"
        else:
            timing = "outside_event_window"
        row_name = str(row.get("organism_name") or "").strip()
        row_name_text = row_name.casefold()
        identity_match = (
            not row_name
            or canonical_key(row_name) == canonical_key(candidate.get("organism_name"))
            or (group_identity and any(term in row_name_text for term in ("candida", "yeast")))
        )
        evaluated_rows.append({
            "assertion": assertion,
            "invasive_specimen": invasive_specimen,
            "timing": timing,
            "identity_match": identity_match,
            "observation_id": row.get("observation_id"),
            "organism_name": row_name or None,
            "specimen_type": row.get("specimen_type"),
            "specimen_category": row.get("specimen_category"),
            "collected_time": row.get("collected_time"),
            "raw_result": row.get("raw_result"),
            "quantitation_status": row.get("quantitation_status"),
        })
    qualifying = [
        row for row in evaluated_rows
        if row["assertion"] == "positive"
        and row["invasive_specimen"]
        and row["identity_match"]
        and row["timing"] == "within_event_window"
    ]
    unknown_timing = [
        row for row in evaluated_rows
        if row["assertion"] == "positive"
        and row["invasive_specimen"]
        and row["identity_match"]
        and row["timing"] == "unknown"
    ]
    return {
        "profile_version": "candida_invasive_evidence_v1",
        "assessment_target": "invasive_candida_infection_not_candida_pneumonia",
        "group_identity": group_identity,
        "detail_available": bool(detail),
        "evaluated_row_count": len(evaluated_rows),
        "qualifying_event_aligned_positive_count": len(qualifying),
        "positive_invasive_timing_unknown_count": len(unknown_timing),
        "exact_event_aligned_invasive_positive": bool(qualifying),
        "rows": evaluated_rows,
        "bdg_used_for_promotion": False,
        "missing_evidence_is_negative": False,
    }


def direct_evidence_timing_profile(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any]:
    """Verify that a direct hospital result is positive and near the index event."""
    if candidate.get("candidate_source") == "hospital_only":
        detail = candidate.get("hospital_evidence_detail") or {}
    else:
        detail = (candidate.get("hospital_profile") or {}).get(
            "hospital_evidence_detail"
        ) or {}
    rows = _module_evidence_rows(detail)
    event_times = candidate_event_times(candidate)
    window_hours = float(policy.get("direct_evidence_event_window_hours", 48))
    evaluated = []
    for row in rows:
        text = f" {_recursive_text(row)} "
        status = str(row.get("detection_status") or "").casefold()
        quantitation = str(row.get("quantitation_status") or "").casefold()
        negative = any(term in text for term in (
            " not detected", " no growth", " negative", " not isolated"
        ))
        positive = not negative and (
            status in {"detected", "positive"}
            or quantitation.startswith("positive")
            or any(term in text for term in (" isolated", " detected", " positive", " grew"))
        )
        observed_at = _parse_datetime(row.get("collected_time"))
        time_basis = "collected_time"
        if observed_at is None:
            observed_at = _parse_datetime(row.get("reported_time"))
            time_basis = "reported_time" if observed_at is not None else "unavailable"
        if observed_at is None or not event_times:
            timing = "unknown"
        elif any(abs((observed_at - event_time).total_seconds()) <= window_hours * 3600
                 for event_time in event_times):
            timing = "within_event_window"
        else:
            timing = "outside_event_window"
        evaluated.append({
            "observation_id": row.get("observation_id"),
            "positive": positive,
            "timing": timing,
            "time_basis": time_basis,
            "collected_time": row.get("collected_time"),
            "reported_time": row.get("reported_time"),
            "specimen_type": row.get("specimen_type"),
            "target": row.get("target") or row.get("organism_name"),
            "raw_result": row.get("raw_result"),
        })
    qualifying = [
        row for row in evaluated
        if row["positive"] and row["timing"] == "within_event_window"
    ]
    return {
        "event_window_hours": window_hours,
        "evaluated_row_count": len(evaluated),
        "event_aligned_positive_count": len(qualifying),
        "event_aligned_positive": bool(qualifying),
        "rows": evaluated,
        "missing_timing_is_negative": False,
    }


def build_patient_promotion_context(
    candidates: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    profiles = {
        canonical_key(candidate.get("organism_name")): candida_evidence_profile(candidate, policy)
        for candidate in candidates
        if ((candidate.get("taxonomy_profile") or {}).get("primary_rule_family") == "candida_or_yeast"
            or candida_group_label(candidate, policy))
    }
    exact_keys = sorted(
        key for key, profile in profiles.items()
        if profile["exact_event_aligned_invasive_positive"] and not profile["group_identity"]
    )
    return {
        "candida_profiles": profiles,
        "exact_event_aligned_invasive_candida_keys": exact_keys,
    }


def evidence_axes(candidate: dict[str, Any], policy: dict[str, Any]) -> dict[str, Any]:
    analytical = candidate.get("analytical_profile") or {}
    taxonomy = candidate.get("taxonomy_profile") or {}
    rank = as_int(analytical.get("best_rank_in_retained_universe"))
    selected_count = as_int(analytical.get("selected_positive_test_count")) or 0
    direct_level = direct_hospital_level(candidate)
    context = specimen_context(candidate)
    consumed_history = list(candidate.get("consumed_evidence_ids") or [])
    raw_host_support = bool(candidate.get("host_support_flag_from_upstream"))
    verified_host = verified_host_support_profile(candidate, policy)
    direct_timing = direct_evidence_timing_profile(candidate, policy)
    return {
        "taxonomy_family": str(taxonomy.get("primary_rule_family") or "unmapped_or_uncertain"),
        "precise_identity": precise_identity(candidate, policy),
        "specimen_context": context,
        "compatible_specimen": context in set(policy["compatible_specimen_contexts"]),
        "best_rank_in_retained_universe": rank,
        "selected_positive_test_count": selected_count,
        "reproducibility_axis": bool(analytical.get("reproducibility_axis")),
        "cross_molecule_selected": bool(analytical.get("cross_molecule_selected")),
        "direct_hospital_level": direct_level,
        "direct_hospital_level_1_2": direct_level is not None and direct_level <= 2,
        "direct_evidence_timing_profile": direct_timing,
        "host_support": raw_host_support,
        "verified_host_support_profile": verified_host,
        "consumed_history_evidence_ids": consumed_history,
        "host_support_eligible_for_promotion": (
            verified_host["eligible"] and not consumed_history
        ),
        "exact_image_context": exact_image_context(candidate),
        "phenotype_used_for_promotion": False,
    }


def promotion_decision(
    candidate: dict[str, Any], policy: dict[str, Any],
    patient_context: dict[str, Any] | None = None,
) -> tuple[str, dict[str, Any], list[str], list[str]]:
    current = str(candidate.get("clinical_decision") or "review_context_needed")
    axes = evidence_axes(candidate, policy)
    family = axes["taxonomy_family"]
    rank = axes["best_rank_in_retained_universe"]
    source = str(candidate.get("candidate_source") or "")
    threshold = policy["thresholds"]
    patient_context = patient_context or {}
    candida_related = family == "candida_or_yeast" or candida_group_label(candidate, policy)
    candida_profile = (patient_context.get("candida_profiles") or {}).get(
        canonical_key(candidate.get("organism_name"))
    ) or (candida_evidence_profile(candidate, policy) if candida_related else None)
    axes["candida_invasive_evidence"] = candida_profile
    rules = ["PROMO-S1-HIGH-ONLY", "PROMO-S2-NO-CROSS-TEST-READ-SUM"]
    if axes["consumed_history_evidence_ids"]:
        rules.append("PROMO-S3-CONSUMED-HISTORY-NOT-REUSED")
    blockers: list[str] = []

    if current != "review_high_priority":
        return current, {
            "evaluated": False,
            "outcome": "not_applicable",
            "route": "preserve_existing_tier",
            "axes": axes,
            "blockers": [],
        }, rules, ["Only existing high-priority candidates are evaluated by this promotion stage."]

    if family in set(policy["blocked_without_invasive_evidence_families"]):
        blockers.append("context-prone family has no implemented invasive or pathogen-specific promotion route")
    if family == "unmapped_or_uncertain" and not candida_related:
        blockers.append("taxonomy remains unmapped or uncertain")
    if not axes["precise_identity"] and not (candida_related and candida_profile["group_identity"]):
        blockers.append("identity is not an exact or species-candidate organism")
    if source == "mngs_positive_reads" and not axes["compatible_specimen"]:
        blockers.append("mNGS signal is not linked to a lower-respiratory or sterile/systemic event")

    route: str | None = None
    enough_repeats = (
        axes["reproducibility_axis"]
        and axes["selected_positive_test_count"] >= threshold["minimum_reproduced_tests"]
    )
    require_event_aligned_direct = bool(
        policy.get("require_event_aligned_direct_evidence", False)
    )
    direct_support_ready = bool(
        axes["direct_hospital_level_1_2"]
        and (
            not require_event_aligned_direct
            or axes["direct_evidence_timing_profile"]["event_aligned_positive"]
        )
    )

    if not blockers and candida_related:
        exact_keys = set(patient_context.get("exact_event_aligned_invasive_candida_keys") or [])
        if candida_profile["exact_event_aligned_invasive_positive"] and axes["precise_identity"]:
            route = "candida_exact_species_invasive_event"
            rules.extend(["PROMO-CAN1-EXACT-INVASIVE-POSITIVE", "PROMO-CAN2-EVENT-WINDOW-48H"])
        elif (candida_profile["group_identity"]
              and candida_profile["exact_event_aligned_invasive_positive"]
              and not exact_keys):
            route = "candida_group_invasive_event_species_unresolved"
            rules.extend(["PROMO-CAN3-GROUP-INVASIVE-POSITIVE", "PROMO-CAN4-SPECIES-UNRESOLVED"])
        elif candida_profile["group_identity"] and exact_keys:
            blockers.append("exact Candida species already represents the invasive event; group evidence is supporting only")
        elif candida_profile["positive_invasive_timing_unknown_count"]:
            blockers.append("positive invasive Candida evidence has unknown timing; missing timing is not treated as negative")
        else:
            blockers.append("no event-aligned positive blood, sterile-site, or tissue Candida evidence")

    if route is None and not blockers and source == "mngs_positive_reads" and family in set(policy["bacterial_families"]):
        allow_analytical_only = bool(
            policy.get("allow_analytical_only_bacterial_promotion", True)
        )
        if (allow_analytical_only and enough_repeats and rank is not None
                and rank <= threshold["bacterial_reproducible_max_rank"]):
            route = "bacterial_reproducible_top3"
            rules.append("PROMO-B1-REPRODUCIBLE-TOP3")
        elif (direct_support_ready and rank is not None
              and rank <= threshold["direct_supported_max_rank"]
              and (enough_repeats or rank <= threshold["bacterial_reproducible_max_rank"])):
            route = "bacterial_direct_support"
            rules.append("PROMO-B2-DIRECT-L1L2")

    if not blockers and family == "other_respiratory_virus":
        if source == "hospital_only" and direct_support_ready:
            route = "respiratory_assay_target"
            rules.append("PROMO-V1-DIRECT-RESPIRATORY-ASSAY")
        elif (source == "mngs_positive_reads" and direct_support_ready
              and rank is not None and rank <= threshold["direct_supported_max_rank"]):
            route = "respiratory_virus_direct_support"
            rules.append("PROMO-V2-MNGS-PLUS-DIRECT-ASSAY")

    if not blockers and canonical_key(candidate.get("organism_name")) == "pneumocystisjirovecii":
        base = enough_repeats and rank is not None and rank <= threshold["opportunistic_max_rank"]
        if (base and axes["host_support_eligible_for_promotion"]
                and rank <= threshold["bacterial_reproducible_max_rank"]):
            route = "pjp_host_risk_triangulation"
            rules.append("PROMO-PJP1-HOST-PLUS-REPRODUCIBLE-TOP3")
        elif base and axes["exact_image_context"] and axes["cross_molecule_selected"]:
            route = "pjp_image_cross_molecule_triangulation"
            rules.extend(["PROMO-PJP2-IMAGE-PLUS-CROSS-MOLECULE", "PROMO-I1-IMAGE-NOT-CONFIRMATION"])

    if not blockers and family == "mold_or_opportunistic_fungus":
        if (source == "mngs_positive_reads" and direct_support_ready
                and enough_repeats and rank is not None
                and rank <= threshold["opportunistic_max_rank"]):
            route = "opportunistic_fungus_direct_support"
            rules.append("PROMO-F1-DIRECT-L1L2-PLUS-REPRODUCIBILITY")

    if not blockers and family == "environmental_low_specificity":
        if (source == "mngs_positive_reads" and direct_support_ready
                and enough_repeats and rank is not None
                and rank <= threshold["environmental_max_rank"]):
            route = "low_specificity_direct_support"
            rules.append("PROMO-E1-STRICT-DIRECT-SUPPORT")

    if route:
        return "picked_shadow", {
            "evaluated": True,
            "outcome": "promoted_to_picked_shadow",
            "route": route,
            "axes": axes,
            "blockers": [],
        }, rules, [
            "The candidate passed identity and specimen gates plus its family-specific independent-evidence route."
        ]

    if (not blockers and not candida_related and require_event_aligned_direct
            and axes["direct_hospital_level_1_2"] and not direct_support_ready):
        blockers.append("direct hospital evidence is not positive and event-aligned within the configured window")
    if not blockers:
        blockers.append("available axes do not satisfy the family-specific Picked promotion route")
    rules.append("PROMO-H1-RETAIN-HIGH")
    return "review_high_priority", {
        "evaluated": True,
        "outcome": "retained_high_priority",
        "route": "no_promotion",
        "axes": axes,
        "blockers": blockers,
    }, rules, [
        "The candidate remains visible at high priority because its evidence does not meet the family-specific Picked gate."
    ]


def promote_candidate(
    candidate: dict[str, Any], policy: dict[str, Any],
    patient_context: dict[str, Any] | None = None,
) -> dict[str, Any]:
    output = copy.deepcopy(candidate)
    output["pre_promotion_decision"] = candidate.get("clinical_decision")
    decision, gate, rules, reasons = promotion_decision(candidate, policy, patient_context)
    output["clinical_decision"] = decision
    output["clinical_level"] = DECISION_LEVELS[decision]
    output["formal_pick_allowed"] = decision == "picked_shadow"
    output["promotion_gate"] = gate
    if decision == "picked_shadow":
        output["selection_role"] = (
            "possible_concurrent_pathogen"
            if gate.get("outcome") == "promoted_to_picked_shadow"
            else "existing_picked_pathogen"
        )
        output["selected_for_report"] = True
    elif decision == "review_high_priority":
        output["selection_role"] = "high_priority_review_candidate"
        output["selected_for_report"] = False
    else:
        output["selection_role"] = "context_review_candidate"
        output["selected_for_report"] = False
    output["rule_ids"] = list(candidate.get("rule_ids") or []) + rules
    output["reasons"] = list(candidate.get("reasons") or []) + reasons
    return output


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "patient_id", "organism_name", "candidate_source", "taxonomy_family",
        "incoming_decision", "pre_promotion_decision", "clinical_decision", "clinical_level",
        "promotion_outcome", "promotion_route", "promotion_blockers", "specimen_context",
        "precise_identity", "best_rank", "selected_positive_test_count", "reproducibility_axis",
        "cross_molecule_selected", "direct_hospital_level", "host_support",
        "exact_image_context", "formal_pick_allowed", "selection_role",
        "selected_for_report", "future_picked_gate_review_required",
        "future_picked_gate_review_status", "rule_ids",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in rows:
            gate = item.get("promotion_gate") or {}
            axes = gate.get("axes") or {}
            writer.writerow({
                "patient_id": item.get("patient_id"),
                "organism_name": item.get("organism_name"),
                "candidate_source": item.get("candidate_source"),
                "taxonomy_family": axes.get("taxonomy_family"),
                "incoming_decision": item.get("incoming_decision"),
                "pre_promotion_decision": item.get("pre_promotion_decision"),
                "clinical_decision": item.get("clinical_decision"),
                "clinical_level": item.get("clinical_level"),
                "promotion_outcome": gate.get("outcome"),
                "promotion_route": gate.get("route"),
                "promotion_blockers": "|".join(gate.get("blockers") or []),
                "specimen_context": axes.get("specimen_context"),
                "precise_identity": axes.get("precise_identity"),
                "best_rank": axes.get("best_rank_in_retained_universe"),
                "selected_positive_test_count": axes.get("selected_positive_test_count"),
                "reproducibility_axis": axes.get("reproducibility_axis"),
                "cross_molecule_selected": axes.get("cross_molecule_selected"),
                "direct_hospital_level": axes.get("direct_hospital_level"),
                "host_support": axes.get("host_support"),
                "exact_image_context": axes.get("exact_image_context"),
                "formal_pick_allowed": item.get("formal_pick_allowed"),
                "selection_role": item.get("selection_role"),
                "selected_for_report": item.get("selected_for_report"),
                "future_picked_gate_review_required": item.get(
                    "future_picked_gate_review_required"
                ),
                "future_picked_gate_review_status": item.get(
                    "future_picked_gate_review_status"
                ),
                "rule_ids": "|".join(item.get("rule_ids") or []),
            })


def run(
    input_root: Path, output_dir: Path, *, policy_path: Path = DEFAULT_POLICY,
    patients: set[int] | None = None,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    policy = read_json(policy_path)
    input_summary_path = input_root / "summary.json"
    input_summary = read_json(input_summary_path) if input_summary_path.is_file() else {}
    history_integrated_input = bool(input_summary.get("history_route_root"))
    paths = sorted((input_root / "patient_outputs").glob(
        "NGS_patient_*_test_aware_clinical_shadow.json"), key=patient_number)
    if patients is not None:
        paths = [path for path in paths if patient_number(path) in patients]
    if not paths:
        raise ValueError("No clinical shadow v1 patient outputs found")
    output_dir.mkdir(parents=True, exist_ok=True)
    patient_dir = output_dir / "patient_outputs"
    patient_dir.mkdir()
    all_rows: list[dict[str, Any]] = []
    expected_input_count = 0
    counts: Counter[str] = Counter()

    for source_path in paths:
        source = read_json(source_path)
        original = source.get("all_forwarded_candidates") or []
        expected_input_count += len(original)
        patient_context = build_patient_promotion_context(original, policy)
        candidates = [promote_candidate(item, policy, patient_context) for item in original]
        if len(candidates) != len(original):
            raise ValueError(f"Candidate loss in {source_path}")
        all_rows.extend(candidates)
        for item in candidates:
            gate = item["promotion_gate"]
            counts[f"clinical:{item['clinical_decision']}"] += 1
            counts[f"selection_role:{item['selection_role']}"] += 1
            counts[f"promotion:{gate['outcome']}"] += 1
            if gate["outcome"] == "promoted_to_picked_shadow":
                counts[f"promotion_route:{gate['route']}"] += 1
            for blocker in gate.get("blockers") or []:
                counts[f"blocker:{blocker}"] += 1

        patient = patient_number(source_path)
        payload = {
            "schema_version": SCHEMA_VERSION,
            "policy_id": policy["policy_id"],
            "answer_blind": True,
            "patient_id": str(patient),
            "candidate_count": len(candidates),
            "picked_shadow": [x for x in candidates if x["clinical_decision"] == "picked_shadow"],
            "existing_picked_pathogens": [
                x for x in candidates if x["selection_role"] == "existing_picked_pathogen"
            ],
            "possible_concurrent_pathogens": [
                x for x in candidates if x["selection_role"] == "possible_concurrent_pathogen"
            ],
            "review_high_priority": [x for x in candidates if x["clinical_decision"] == "review_high_priority"],
            "review_context_needed": [x for x in candidates if x["clinical_decision"] == "review_context_needed"],
            "all_forwarded_candidates": candidates,
            "no_forced_pick": not any(x["clinical_decision"] == "picked_shadow" for x in candidates),
            "source": {
                "clinical_shadow_v1_file": str(source_path.resolve()),
                "clinical_shadow_v1_sha256": sha256_file(source_path),
            },
            "constraints": policy.get("constraints") or [],
            "patient_promotion_context": patient_context,
        }
        write_json(patient_dir / f"NGS_patient_{patient}_test_aware_clinical_shadow.json", payload)

    if len(all_rows) != expected_input_count:
        raise ValueError(
            f"Expected {expected_input_count} candidates from the input files, "
            f"got {len(all_rows)}"
        )
    write_csv(output_dir / "clinical_decisions.csv", all_rows)
    write_csv(
        output_dir / "possible_concurrent_pathogens.csv",
        [x for x in all_rows if x["selection_role"] == "possible_concurrent_pathogen"],
    )
    patient_ids = sorted({int(item["patient_id"]) for item in all_rows})
    patients_with_picks = {
        int(item["patient_id"]) for item in all_rows if item["clinical_decision"] == "picked_shadow"
    }
    summary = {
        "schema_version": SCHEMA_VERSION,
        "policy_id": policy["policy_id"],
        "answer_blind": True,
        "candidate_count": len(all_rows),
        "expected_candidate_count_from_inputs": expected_input_count,
        "patient_count": len(paths),
        "counts": dict(sorted(counts.items())),
        "patients_without_picked": [x for x in patient_ids if x not in patients_with_picks],
        "source_clinical_shadow_root": str(input_root.resolve()),
        "history_integrated_input": history_integrated_input,
        "selected_patients": sorted(patients) if patients is not None else None,
        "policy_file": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "constraints": policy.get("constraints") or [],
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--patients", type=int, nargs="*")
    args = parser.parse_args()
    selected = set(args.patients) if args.patients else None
    print(json.dumps(run(
        args.input_root, args.output_dir, policy_path=args.policy, patients=selected,
    ),
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

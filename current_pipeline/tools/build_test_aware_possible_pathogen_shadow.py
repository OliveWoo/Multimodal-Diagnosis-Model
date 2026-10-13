"""Build a non-destructive possible-pathogen reporting shadow from clinical V4.

The source clinical tier is immutable. This stage adds a reporting role for
high-recall review without converting analytical evidence into a formal causal
claim. Benchmark answers are intentionally not accepted by this tool.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import re
from datetime import timedelta
from collections import Counter
from pathlib import Path
from typing import Any

from tools.build_index_event_clinical_timeline_shadow import sha256_file
from tools.pathogen_normalization import canonical_key
from tools import promote_test_aware_clinical_high_shadow as promotion_base
from tools.promote_test_aware_clinical_high_shadow import (
    as_int,
    patient_number,
    precise_identity,
    read_json,
    specimen_context,
    write_json,
)


DEFAULT_POLICY = Path("rules/test_aware_possible_pathogen_v1.json")
SCHEMA_VERSION = "test_aware_possible_pathogen_shadow.v1"
STRICT_ROLE = "strict_picked_report"
ANALYTICAL_DIRECT_ROLE = "analytical_or_direct_possible_pathogen"
LOW_SPECIFICITY_ROLE = "reproducible_low_specificity_possible_pathogen"
SYSTEMIC_ROLE = "possible_concurrent_systemic_pathogen"
NOT_SELECTED_ROLE = "not_selected_for_complete_report"
PJP_ROLE = "pneumocystis_possible_pathogen"
REACTIVATION_ROLE = "reactivation_or_shedding_possible_pathogen"
CANDIDA_ROLE = "candida_clinical_convergence_possible_pathogen"
STRONG_CROSS_MOLECULE_ROLE = "strong_cross_molecule_possible_pathogen"
ASPIRATION_ROLE = "aspiration_supported_possible_pathogen"
FALLBACK_ROLE = "fallback_best_available_possible_pathogen"
LOW_SPECIFICITY_REPEAT_ROLE = "analytically_compelling_low_specificity_repeat_possible"
MNGS_ONLY_ROLE = "mngs_only_possible_pathogen"
ASPERGILLUS_GROUP_ROLE = "aspergillus_group_mycology_convergence_possible"
RESPIRATORY_VIRUS_ROLE = "respiratory_virus_mngs_possible_pathogen"
OPPORTUNISTIC_FUNGUS_ROLE = (
    "opportunistic_fungus_direct_convergence_possible_pathogen"
)


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key == "extends":
            continue
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def load_policy(path: Path) -> tuple[dict[str, Any], Path | None]:
    policy = read_json(path)
    parent_value = policy.get("extends")
    if not parent_value:
        return policy, None
    parent_path = Path(str(parent_value))
    if not parent_path.is_absolute():
        parent_path = path.parent / parent_path
    parent = read_json(parent_path)
    return _deep_merge(parent, policy), parent_path


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().casefold() in {"1", "true", "yes"}


def _axes(candidate: dict[str, Any]) -> dict[str, Any]:
    return (candidate.get("promotion_gate") or {}).get("axes") or {}


def _positive_mngs_signals(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    """Return positive per-test mNGS signals without combining their reads."""
    profile = candidate.get("analytical_profile") or {}
    return [
        row for row in profile.get("per_test_signals") or []
        if isinstance(row, dict) and float(row.get("reads") or 0) > 0
    ]


def _has_positive_hospital_evidence(candidate: dict[str, Any]) -> bool:
    """Detect candidate-specific hospital corroboration, not file availability."""
    if str(candidate.get("candidate_source") or "").startswith("hospital_only"):
        return True
    if candidate.get("hospital_evidence_detail"):
        return True
    hospital = candidate.get("hospital_profile") or {}
    if hospital.get("hospital_evidence_detail"):
        return True
    if _bool(hospital.get("exact_or_alias_evidence")):
        return True
    return bool(_event_aligned_direct_rows(candidate))


def _has_clinical_context(candidate: dict[str, Any]) -> bool:
    """Detect available contextual evidence without treating it as organism proof."""
    phenotype = candidate.get("phenotype_evidence") or {}
    return bool(
        candidate.get("host_evidence")
        or candidate.get("image_evidence")
        or candidate.get("history_evidence_ids")
        or phenotype.get("routed_phenotype_context")
        or phenotype.get("verified_phenotype_evidence")
        or _bool(_axes(candidate).get("exact_image_context"))
    )


def candidate_decision_evidence_mode(candidate: dict[str, Any]) -> str:
    """Classify which positive/contextual evidence can support one candidate."""
    has_mngs = bool(_positive_mngs_signals(candidate))
    has_hospital = _has_positive_hospital_evidence(candidate)
    has_context = _has_clinical_context(candidate)
    if has_mngs and has_hospital:
        return "mngs_and_direct_hospital"
    if has_mngs and has_context:
        return "mngs_and_clinical_context"
    if has_mngs:
        return "mngs_only"
    if has_hospital:
        return "hospital_only"
    if has_context:
        return "clinical_context_only"
    return "no_positive_evidence"


def decision_evidence_profile(
    candidates: list[dict[str, Any]],
) -> dict[str, Any]:
    """Summarize evidence reaching the decision packet without guessing files."""
    has_mngs = any(_positive_mngs_signals(item) for item in candidates)
    has_hospital = any(_has_positive_hospital_evidence(item) for item in candidates)
    has_context = any(_has_clinical_context(item) for item in candidates)
    if has_mngs and has_hospital:
        mode = "multimodal_mngs_hospital"
    elif has_mngs and has_context:
        mode = "mngs_plus_clinical_context"
    elif has_mngs:
        mode = "mngs_only"
    elif has_hospital:
        mode = "hospital_only"
    else:
        mode = "insufficient_positive_evidence"
    positive_sources = []
    if has_mngs:
        positive_sources.append("mngs")
    if has_hospital:
        positive_sources.append("hospital_pathogen_specific")
    if has_context:
        positive_sources.append("clinical_context")
    return {
        "decision_evidence_mode": mode,
        "positive_or_contextual_sources": positive_sources,
        "mngs_positive_available": has_mngs,
        "hospital_pathogen_specific_positive_available": has_hospital,
        "clinical_context_available": has_context,
        "absence_semantics": (
            "A missing source here means it was not available as positive decision "
            "evidence; it must not be interpreted as a negative clinical result."
        ),
    }


def _event_aligned_direct_rows(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    profile = _axes(candidate).get("direct_evidence_timing_profile") or {}
    return [
        row for row in profile.get("rows") or []
        if _bool(row.get("positive")) and row.get("timing") == "within_event_window"
    ]


def _sterile_rows(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> list[dict[str, Any]]:
    terms = tuple(str(value).casefold() for value in policy["sterile_specimen_terms"])
    output = []
    for row in _event_aligned_direct_rows(candidate):
        text = " ".join(str(row.get(key) or "").casefold()
                        for key in ("specimen_type", "raw_result"))
        if any(term in text for term in terms):
            output.append(row)
    return output


def _host_context_available(
    candidate: dict[str, Any], axes: dict[str, Any]
) -> bool:
    if _bool(axes.get("host_support")):
        return True
    return str(candidate.get("history_policy_family") or "") in {
        "opportunistic_host_risk",
        "reactivation_interpretation",
    } and bool(candidate.get("history_evidence_ids")) and any(
        str(rule).startswith(("HIST-OPP4-", "HIST-REACT4-"))
        for rule in candidate.get("history_rule_ids") or []
    )


def _provisional_phenotype_host_context(
    candidate: dict[str, Any], config: dict[str, Any]
) -> dict[str, Any]:
    """Return tightly scoped phenotype host context for Possible-only routes."""
    if not config.get("enabled"):
        return {"available": False, "evidence": []}
    accepted_phenotypes = {
        str(value) for value in config.get("accepted_phenotypes") or []
    }
    accepted_statuses = {
        str(value) for value in config.get("accepted_statuses") or ["YES"]
    }
    accepted_confidences = {
        str(value) for value in config.get("accepted_confidences") or ["HIGH"]
    }
    accepted_temporal_relations = {
        str(value)
        for value in config.get("accepted_temporal_relations")
        or ["PRE_PNEUMONIA"]
    }
    accepted_qc_statuses = {
        str(value)
        for value in config.get("accepted_qc_statuses") or ["PASS"]
    }
    rows = []
    phenotype = candidate.get("phenotype_evidence") or {}
    for row in phenotype.get("routed_phenotype_context") or []:
        if not isinstance(row, dict):
            continue
        if (
            accepted_phenotypes
            and str(row.get("phenotype") or "") not in accepted_phenotypes
        ):
            continue
        if str(row.get("status") or "") not in accepted_statuses:
            continue
        if str(row.get("confidence") or "") not in accepted_confidences:
            continue
        if (
            str(row.get("temporal_relation") or "")
            not in accepted_temporal_relations
        ):
            continue
        if (
            str(row.get("phenotype_qc_status") or "PASS")
            not in accepted_qc_statuses
        ):
            continue
        rows.append({
            "source_id": row.get("source_id"),
            "phenotype": row.get("phenotype"),
            "status": row.get("status"),
            "confidence": row.get("confidence"),
            "temporal_relation": row.get("temporal_relation"),
            "production_link_verified": _bool(
                phenotype.get("production_link_verified")
            ),
        })
    return {
        "available": bool(rows),
        "evidence": rows,
        "possible_only": True,
        "production_link_verified": _bool(
            phenotype.get("production_link_verified")
        ),
    }


def _recorded_host_review_context(
    candidate: dict[str, Any], config: dict[str, Any]
) -> dict[str, Any]:
    """Resolve pre-sample host context for Possible-only review routes."""
    if not config.get("enabled"):
        return {"available": False, "evidence": [], "possible_only": True}
    accepted_roles = {
        str(value)
        for value in config.get("accepted_roles") or ["host_medication"]
    }
    verified_assertions = {
        str(value)
        for value in config.get("accepted_verified_assertions")
        or [
            "verified_active_exposure",
            "verified_high_risk_host",
            "verified_immunosuppression",
        ]
    }
    recorded_assertions = {
        str(value)
        for value in config.get("accepted_recorded_exposure_assertions")
        or ["recorded_pre_sample_exposure_not_administration_verified"]
    }
    verified_rows = []
    recorded_rows = []
    seen = set()
    for row in candidate.get("host_evidence") or []:
        if not isinstance(row, dict):
            continue
        role = str(row.get("evidence_role") or "")
        assertion = str(row.get("assertion") or "")
        if role not in accepted_roles:
            continue
        if assertion not in verified_assertions | recorded_assertions:
            continue
        identity = (
            str(row.get("label") or "").strip().casefold(),
            assertion,
            str((row.get("source_ref") or {}).get("source_claim_group") or ""),
        )
        if identity in seen:
            continue
        seen.add(identity)
        evidence = {
            "evidence_id": row.get("evidence_id"),
            "label": row.get("label"),
            "assertion": assertion,
            "observed_at": row.get("observed_at"),
            "episode_evidence_status": row.get("episode_evidence_status"),
            "source_ref": row.get("source_ref"),
        }
        if assertion in verified_assertions:
            verified_rows.append(evidence)
        else:
            recorded_rows.append(evidence)
    minimum_recorded = int(
        config.get("minimum_distinct_recorded_medications", 2)
    )
    recorded_labels = {
        str(row.get("label") or "").strip().casefold()
        for row in recorded_rows
        if str(row.get("label") or "").strip()
    }
    available = bool(verified_rows) or len(recorded_labels) >= minimum_recorded
    return {
        "available": available,
        "evidence": [*verified_rows, *recorded_rows],
        "verified_evidence_count": len(verified_rows),
        "recorded_distinct_medication_count": len(recorded_labels),
        "minimum_distinct_recorded_medications": minimum_recorded,
        "possible_only": True,
        "administration_verified": bool(verified_rows),
        "missing_dose_is_negative": False,
    }


def _candida_evidence_axes(
    candidate: dict[str, Any], axes: dict[str, Any]
) -> dict[str, Any]:
    profile = axes.get("candida_invasive_evidence") or {}
    lower_respiratory_rows = []
    invasive_rows = []
    for row in profile.get("rows") or []:
        if not isinstance(row, dict):
            continue
        if str(row.get("assertion") or "") != "positive":
            continue
        if str(row.get("timing") or "") != "within_event_window":
            continue
        if not _bool(row.get("identity_match")):
            continue
        if _bool(row.get("invasive_specimen")):
            invasive_rows.append(row)
        elif str(row.get("specimen_category") or "") == "Lower_Respiratory":
            lower_respiratory_rows.append(row)
    return {
        "exact_event_aligned_invasive_positive": bool(invasive_rows),
        "invasive_rows": invasive_rows,
        "exact_event_aligned_lower_respiratory_positive": bool(
            lower_respiratory_rows
        ),
        "lower_respiratory_rows": lower_respiratory_rows,
        "missing_evidence_is_negative": False,
    }


def _reactivation_status(candidate: dict[str, Any]) -> dict[str, Any]:
    profile = candidate.get("reactivation_interpretation") or {}
    explicit = _bool(profile.get("has_explicit_current_reactivation"))
    resolved = _bool(profile.get("has_resolved_or_asymptomatic_shedding"))
    quantitative = _bool(profile.get("has_quantitative_or_tissue_evidence"))
    syndrome = _bool(profile.get("has_compatible_pulmonary_syndrome"))
    if explicit and (quantitative or syndrome):
        status = "explicit_current_reactivation_with_clinical_convergence"
    elif explicit:
        status = "explicit_current_reactivation_without_full_convergence"
    elif resolved:
        status = "resolved_or_asymptomatic_shedding_documented"
    else:
        status = "reactivation_status_unknown"
    return {
        "status": status,
        "explicit_current_reactivation": explicit,
        "resolved_or_asymptomatic_shedding": resolved,
        "quantitative_or_tissue_evidence": quantitative,
        "compatible_pulmonary_syndrome": syndrome,
        "missing_reactivation_wording_is_negative": False,
    }


def _recursive_text(value: Any) -> str:
    if isinstance(value, dict):
        return " ".join(_recursive_text(item) for item in value.values())
    if isinstance(value, list):
        return " ".join(_recursive_text(item) for item in value)
    return str(value or "").casefold()


def _date_only_event_overlap(
    observed_value: Any,
    event_values: list[Any],
    *,
    window_hours: float,
) -> bool:
    """Treat a date-only timestamp as an interval, not as midnight certainty."""
    text = str(observed_value or "").strip()
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", text):
        return False
    start = promotion_base._parse_datetime(text)
    if start is None:
        return False
    end = start + timedelta(days=1)
    for raw_event in event_values:
        event = promotion_base._parse_datetime(raw_event)
        if event is None:
            continue
        if start <= event <= end:
            nearest_hours = 0.0
        else:
            nearest = start if event < start else end
            nearest_hours = abs((event - nearest).total_seconds()) / 3600
        if nearest_hours <= window_hours:
            return True
    return False


def _event_aligned_or_date_overlap_rows(
    candidate: dict[str, Any], config: dict[str, Any]
) -> list[dict[str, Any]]:
    rows = ((_axes(candidate).get("direct_evidence_timing_profile") or {}).get("rows") or [])
    event_values = [
        event.get("collected_time")
        for event in candidate.get("linked_events") or []
        if isinstance(event, dict) and event.get("collected_time")
    ]
    window_hours = float(config.get("event_window_hours", 48))
    output: list[dict[str, Any]] = []
    for row in rows:
        if not isinstance(row, dict) or not _bool(row.get("positive")):
            continue
        alignment = None
        if str(row.get("timing") or "") == "within_event_window":
            alignment = "exact_timestamp_within_window"
        elif config.get("allow_date_only_interval_overlap", False) and _date_only_event_overlap(
            row.get("collected_time"), event_values, window_hours=window_hours
        ):
            alignment = "date_only_interval_overlaps_window"
        if alignment:
            output.append({**row, "event_alignment_mode": alignment})
    return output


def _accepted_host_risk_rows(
    candidate: dict[str, Any], config: dict[str, Any]
) -> list[dict[str, Any]]:
    terms = tuple(str(value).casefold() for value in config.get("host_risk_terms") or [])
    assertions = set(config.get("accepted_host_assertions") or [
        "recorded_history",
        "verified_active_exposure",
        "verified_high_risk_host",
        "verified_immunosuppression",
    ])
    return [
        row for row in candidate.get("host_evidence") or []
        if isinstance(row, dict)
        and str(row.get("assertion") or "") in assertions
        and any(term in str(row.get("label") or "").casefold() for term in terms)
    ]


def resistance_interpretation(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any]:
    config = policy.get("resistance_interpretation") or {}
    key = canonical_key(candidate.get("organism_name"))
    target_keys = set(config.get("organism_keys") or [])
    if not config.get("enabled") or key not in target_keys:
        return {"status": "not_applicable"}
    detail = candidate.get("hospital_evidence_detail") or (
        (candidate.get("hospital_profile") or {}).get("hospital_evidence_detail")
    ) or {}
    text = f" {_recursive_text(detail)} "
    matched_terms = [
        term for term in config.get("carbapenem_resistance_terms") or []
        if str(term).casefold() in text
    ]
    status = "confirmed_carbapenem_resistant" if matched_terms else "unknown"
    return {
        "organism_identity": candidate.get("organism_name"),
        "resistance_target": "carbapenem_resistance",
        "status": status,
        "matched_terms": matched_terms,
        "identity_and_resistance_are_separate": True,
        "missing_resistance_evidence_is_negative": False,
    }


def _recall_first_v2_decision(
    candidate: dict[str, Any],
    policy: dict[str, Any],
    evaluated_axes: dict[str, Any],
) -> dict[str, Any] | None:
    if not policy.get("enable_recall_first_v2"):
        return None
    family = str(evaluated_axes["taxonomy_family"])
    context = str(evaluated_axes["specimen_context"])
    exact = bool(evaluated_axes["precise_identity"])
    rank = evaluated_axes["best_rank_in_retained_universe"]
    positive_tests = int(evaluated_axes["selected_positive_test_count"] or 0)
    reproducible = bool(evaluated_axes["reproducibility_axis"])
    cross_molecule = bool(evaluated_axes["cross_molecule_selected"])
    direct_level = as_int(evaluated_axes.get("direct_hospital_level"))
    host_context = _host_context_available(candidate, _axes(candidate))
    provisional_host_config = (
        policy.get("provisional_phenotype_host_context") or {}
    )
    provisional_host = _provisional_phenotype_host_context(
        candidate, provisional_host_config
    )
    recorded_host_config = policy.get("recorded_host_review_context") or {}
    recorded_host = _recorded_host_review_context(
        candidate, recorded_host_config
    )
    reactivation = _reactivation_status(candidate)
    evaluated_axes["recorded_host_context_available"] = host_context
    if provisional_host_config.get("enabled"):
        evaluated_axes["provisional_phenotype_host_context"] = provisional_host
    if recorded_host_config.get("enabled"):
        evaluated_axes["recorded_host_review_context"] = recorded_host
    evaluated_axes["reactivation_status"] = reactivation

    repeat_route = policy.get("low_specificity_repeat_possible") or {}
    analytical_profile = candidate.get("analytical_profile") or {}
    repeat_gate = analytical_profile.get("low_specificity_repeat_high_gate") or {}
    required_repeat_rule = str(
        repeat_route.get("required_scorer_rule_id") or ""
    )
    repeat_rule_present = required_repeat_rule in (
        set(candidate.get("rule_ids") or [])
        | set((candidate.get("source_candidate") or {}).get("rule_ids") or [])
    )
    low_specificity_repeat_possible = bool(
        repeat_route.get("enabled")
        and str(candidate.get("clinical_decision") or "")
        == repeat_route.get("required_clinical_decision", "review_high_priority")
        and family in set(repeat_route.get("families") or [])
        and exact
        and context == repeat_route.get("required_specimen_context", "lower_respiratory")
        and repeat_rule_present
        and _bool(repeat_gate.get("eligible"))
    )
    evaluated_axes["low_specificity_repeat_high_gate"] = repeat_gate
    if low_specificity_repeat_possible:
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": "analytically_compelling_low_specificity_repeat",
            "reporting_role": LOW_SPECIFICITY_REPEAT_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [
                "Technical branches support analytical reproducibility but are not independent biological specimens.",
                "This route creates Possible only; it cannot create Picked without an independent clinical axis.",
            ],
        }

    respiratory_virus = policy.get("respiratory_virus_mngs_possible") or {}
    taxonomy_profile = candidate.get("taxonomy_profile") or {}
    curated_group_label = bool(
        respiratory_virus.get("allow_curated_group_or_species_label", False)
        and str(taxonomy_profile.get("mapping_status") or "")
        == str(respiratory_virus.get("curated_mapping_status") or "exact_species")
        and str(taxonomy_profile.get("rule_source") or "")
        == str(respiratory_virus.get("curated_rule_source") or "curated_exact_species")
        and str(taxonomy_profile.get("taxonomic_rank") or "")
        == str(respiratory_virus.get("curated_taxonomic_rank") or "group_or_species_label")
    )
    curated_exact_species = bool(
        respiratory_virus.get("allow_curated_exact_species_identity", False)
        and str(taxonomy_profile.get("mapping_status") or "")
        == str(respiratory_virus.get("curated_mapping_status") or "exact_species")
        and str(taxonomy_profile.get("rule_source") or "")
        == str(respiratory_virus.get("curated_rule_source") or "curated_exact_species")
        and str(taxonomy_profile.get("taxonomic_rank") or "")
        in {"species", "species_candidate"}
    )
    accepted_rna_qc = set(
        respiratory_virus.get("accepted_positive_test_qc_statuses")
        or ["evaluable"]
    )
    positive_rna_signals = [
        row
        for row in (candidate.get("analytical_profile") or {}).get(
            "per_test_signals"
        ) or []
        if isinstance(row, dict)
        and float(row.get("reads") or 0) > 0
        and str(row.get("nucleic_type") or "").upper() == "RNA"
        and str(row.get("qc_status") or "") in accepted_rna_qc
    ]
    respiratory_virus_identity_ok = bool(
        curated_exact_species or curated_group_label
    )
    respiratory_virus_route = bool(
        respiratory_virus.get("enabled")
        and family
        == str(respiratory_virus.get("family") or "other_respiratory_virus")
        and str(candidate.get("clinical_decision") or "")
        == str(
            respiratory_virus.get(
                "required_clinical_decision", "review_high_priority"
            )
        )
        and context
        in set(
            respiratory_virus.get("compatible_specimen_contexts")
            or ["lower_respiratory"]
        )
        and rank is not None
        and rank <= int(respiratory_virus.get("maximum_rank", 3))
        and positive_tests
        >= int(respiratory_virus.get("minimum_positive_tests", 1))
        and respiratory_virus_identity_ok
        and (
            bool(positive_rna_signals)
            or not respiratory_virus.get("require_positive_rna_test", True)
        )
    )
    evaluated_axes["respiratory_virus_identity"] = {
        "precise_species_identity": exact,
        "curated_exact_species_identity": curated_exact_species,
        "curated_group_label": curated_group_label,
        "positive_rna_test_count": len(positive_rna_signals),
        "species_confirmation_required": curated_group_label,
    }
    if respiratory_virus_route:
        cautions = [
            "A lower-respiratory RNA mNGS signal supports Possible only; it cannot create Picked without independent confirmation.",
            "Missing targeted confirmation is unknown and is not treated as negative.",
        ]
        if curated_group_label:
            cautions.insert(
                0,
                "The source is a curated respiratory-virus group/species label; the report must preserve species-unresolved identity and request targeted confirmation.",
            )
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": (
                "curated_respiratory_virus_group_label_lower_respiratory_rna"
                if curated_group_label
                else "curated_exact_respiratory_virus_lower_respiratory_rna"
            ),
            "reporting_role": RESPIRATORY_VIRUS_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": cautions,
        }

    pjp = policy.get("pjp_possible") or {}
    is_pjp = (
        canonical_key(candidate.get("organism_name"))
        == str(pjp.get("canonical_organism_key") or "")
    )
    pjp_host_context = bool(
        host_context
        or (
            pjp.get("allow_provisional_phenotype_host_context", False)
            and provisional_host["available"]
        )
        or (
            pjp.get("allow_recorded_host_review_context", False)
            and recorded_host["available"]
        )
    )
    pjp_route = (
        is_pjp
        and exact
        and context in set(pjp.get("compatible_specimen_contexts") or [])
        and rank is not None
        and rank <= int(pjp.get("maximum_rank", 3))
        and positive_tests >= int(pjp.get("minimum_positive_tests", 2))
        and (reproducible or not pjp.get("require_reproducibility", True))
        and (
            pjp_host_context
            or not pjp.get("require_host_context", True)
        )
    )
    if pjp_route:
        provisional_host_used = bool(
            not host_context
            and pjp.get("allow_provisional_phenotype_host_context", False)
            and provisional_host["available"]
        )
        recorded_review_host_used = bool(
            not host_context
            and not provisional_host_used
            and pjp.get("allow_recorded_host_review_context", False)
            and recorded_host["available"]
        )
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": "pneumocystis_reproducible_top_rank_with_recorded_host_context",
            "reporting_role": PJP_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [
                (
                    "Policy-qualified provisional host context supports Possible only; it cannot create Picked."
                    if provisional_host_used
                    else (
                        "Pre-sample medication history supports Possible-only review; missing dose or unverified administration is unknown and cannot create Picked."
                        if recorded_review_host_used
                        else "Recorded host context supports Possible only; dose-unverified exposure cannot create Picked."
                    )
                ),
                "Missing PCR, beta-D-glucan, staining, or tissue evidence is unknown and not negative.",
            ],
        }

    virus = policy.get("reactivation_possible") or {}
    virus_context_ok = context in set(
        virus.get("compatible_specimen_contexts") or []
    )
    organism_key = canonical_key(candidate.get("organism_name"))
    candidate_specific_cross_keys = set(
        virus.get("candidate_specific_cross_molecule_organism_keys") or []
    )
    uses_candidate_specific_cross_route = (
        not candidate_specific_cross_keys
        or organism_key in candidate_specific_cross_keys
    )
    virus_minimum_tests = int(
        virus.get("explicit_current_reactivation_minimum_positive_tests", 2)
        if reactivation["explicit_current_reactivation"]
        else (
            virus.get("cross_molecule_minimum_positive_tests", 2)
            if uses_candidate_specific_cross_route
            else virus.get(
                "legacy_cross_molecule_minimum_positive_tests",
                virus.get("cross_molecule_minimum_positive_tests", 2),
            )
        )
    )
    event_aligned_direct_rows = _event_aligned_direct_rows(candidate)
    candidate_specific_current_support = bool(
        event_aligned_direct_rows
        or reactivation["explicit_current_reactivation"]
        or reactivation["quantitative_or_tissue_evidence"]
    )
    strong_virus = (
        cross_molecule
        and reproducible
        and rank is not None
        and rank <= int(virus.get("cross_molecule_maximum_rank", 3))
        and positive_tests >= virus_minimum_tests
    )
    virus_possible_host_context = bool(
        host_context
        or (
            virus.get("allow_provisional_phenotype_host_context", False)
            and provisional_host["available"]
        )
        or (
            virus.get("allow_recorded_host_review_context", False)
            and recorded_host["available"]
        )
    )
    virus_clinical_support = bool(
        virus_possible_host_context
        or reactivation["explicit_current_reactivation"]
        or reactivation["quantitative_or_tissue_evidence"]
        or reactivation["compatible_pulmonary_syndrome"]
        or bool(event_aligned_direct_rows)
        or _bool(_axes(candidate).get("exact_image_context"))
    )
    if (
        virus.get("require_clinical_support_for_cross_molecule", False)
        and uses_candidate_specific_cross_route
    ):
        strong_virus = strong_virus and virus_clinical_support
    if virus.get(
        "require_candidate_specific_current_support_for_cross_molecule", False
    ) and uses_candidate_specific_cross_route:
        strong_virus = strong_virus and candidate_specific_current_support
    host_supported_single_contexts = set(
        virus.get("host_supported_single_test_contexts") or []
    )
    host_supported_single_keys = set(
        virus.get("host_supported_single_test_organism_keys") or []
    )
    host_supported_single = (
        virus.get("enable_host_supported_single_test", True)
        and virus_possible_host_context
        and (
            not host_supported_single_contexts
            or context in host_supported_single_contexts
        )
        and (
            not host_supported_single_keys
            or organism_key in host_supported_single_keys
        )
        and rank is not None
        and rank <= int(
            virus.get("host_supported_single_test_maximum_rank", 1)
        )
        and positive_tests >= int(
            virus.get("host_supported_single_test_minimum_positive_tests", 1)
        )
        and positive_tests <= int(
            virus.get("host_supported_single_test_maximum_positive_tests", 1)
        )
    )
    direct_supported_single = (
        virus.get("enable_direct_supported_single_test", False)
        and bool(event_aligned_direct_rows)
        and rank is not None
        and rank <= int(virus.get("direct_supported_single_test_maximum_rank", 1))
        and positive_tests >= 1
    )
    host_supported_reproducible_contexts = set(
        virus.get("host_supported_reproducible_contexts") or []
    )
    host_supported_reproducible_keys = set(
        virus.get("host_supported_reproducible_organism_keys") or []
    )
    reproducible_host_context = (
        bool(recorded_host["available"])
        if virus.get(
            "host_supported_reproducible_require_recorded_host_review_context",
            False,
        )
        else virus_possible_host_context
    )
    host_supported_reproducible = bool(
        virus.get("enable_host_supported_reproducible", False)
        and reproducible_host_context
        and (
            not host_supported_reproducible_contexts
            or context in host_supported_reproducible_contexts
        )
        and (
            not host_supported_reproducible_keys
            or organism_key in host_supported_reproducible_keys
        )
        and rank is not None
        and rank
        <= int(virus.get("host_supported_reproducible_maximum_rank", 3))
        and positive_tests
        >= int(
            virus.get("host_supported_reproducible_minimum_positive_tests", 2)
        )
        and (
            reproducible
            or not virus.get(
                "host_supported_reproducible_require_reproducibility", True
            )
        )
    )
    virus_route = (
        family == virus.get("family")
        and exact
        and virus_context_ok
        and (
            strong_virus
            or host_supported_single
            or direct_supported_single
            or host_supported_reproducible
        )
    )
    if virus_route:
        evaluated_axes["reactivation_candidate_specific_current_support"] = (
            candidate_specific_current_support
        )
        if virus.get("enable_host_supported_reproducible", False):
            evaluated_axes["reactivation_reproducible_host_source"] = {
                "recorded_host_review_context_required": _bool(
                    virus.get(
                        "host_supported_reproducible_require_recorded_host_review_context"
                    )
                ),
                "recorded_host_review_context_available": bool(
                    recorded_host["available"]
                ),
                "provisional_phenotype_alone_is_sufficient": False,
            }
        if reactivation["explicit_current_reactivation"]:
            cautions = [
                "Candidate-matched current reactivation is documented; incidental shedding concern is reduced but not eliminated."
            ]
        else:
            cautions = [
                "Reactivation or shedding is not excluded; absent reactivation wording is unknown, not negative."
            ]
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": (
                "reactivation_family_cross_molecule_reproducible_signal"
                if strong_virus
                else (
                    "reactivation_family_event_aligned_direct_supported_signal"
                    if direct_supported_single
                    else (
                        "reactivation_family_reproducible_signal_with_recorded_host_context"
                        if host_supported_reproducible
                        else "reactivation_family_top_rank_with_recorded_host_context"
                    )
                )
            ),
            "reporting_role": REACTIVATION_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": cautions,
        }

    candida = policy.get("candida_possible") or {}
    candida_axes = (
        _candida_evidence_axes(candidate, _axes(candidate))
        if candida.get("enabled")
        else {
            "exact_event_aligned_invasive_positive": False,
            "exact_event_aligned_lower_respiratory_positive": False,
        }
    )
    if candida.get("enabled"):
        evaluated_axes["candida_clinical_evidence"] = candida_axes
    candida_invasive_route = bool(
        candida.get("enabled")
        and family == candida.get("family", "candida_or_yeast")
        and exact
        and candida_axes["exact_event_aligned_invasive_positive"]
    )
    candida_respiratory_route = bool(
        candida.get("enabled")
        and family == candida.get("family", "candida_or_yeast")
        and exact
        and context == "lower_respiratory"
        and rank is not None
        and rank <= int(candida.get("respiratory_maximum_rank", 2))
        and positive_tests
        >= int(candida.get("respiratory_minimum_positive_tests", 2))
        and (
            reproducible
            or not candida.get("respiratory_require_reproducibility", True)
        )
        and (
            cross_molecule
            or not candida.get("respiratory_require_cross_molecule", True)
        )
        and candida_axes[
            "exact_event_aligned_lower_respiratory_positive"
        ]
    )
    if candida_invasive_route or candida_respiratory_route:
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": (
                "candida_exact_event_aligned_invasive_evidence"
                if candida_invasive_route
                else "candida_strong_mngs_with_exact_event_aligned_lower_respiratory_culture"
            ),
            "reporting_role": CANDIDA_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [
                "This route flags possible clinically relevant Candida and does not prove Candida pneumonia.",
                "Species identity is not transferred between different Candida species.",
            ],
        }

    opportunistic_fungus = policy.get("opportunistic_fungus_possible") or {}
    fungus_direct_level_ok = bool(
        direct_level is not None
        and direct_level
        <= int(
            opportunistic_fungus.get("maximum_direct_hospital_level", 3)
        )
    )
    opportunistic_fungus_route = bool(
        opportunistic_fungus.get("enabled")
        and family
        == opportunistic_fungus.get(
            "family", "mold_or_opportunistic_fungus"
        )
        and str(candidate.get("clinical_decision") or "")
        in set(
            opportunistic_fungus.get("allowed_clinical_decisions")
            or ["review_high_priority", "review_context_needed"]
        )
        and exact
        and context
        in set(
            opportunistic_fungus.get("compatible_specimen_contexts")
            or ["lower_respiratory", "sterile_or_systemic"]
        )
        and rank is not None
        and rank <= int(opportunistic_fungus.get("maximum_rank", 3))
        and positive_tests
        >= int(opportunistic_fungus.get("minimum_positive_tests", 2))
        and (
            reproducible
            or not opportunistic_fungus.get("require_reproducibility", True)
        )
        and (
            cross_molecule
            or not opportunistic_fungus.get("require_cross_molecule", True)
        )
        and bool(event_aligned_direct_rows)
        and fungus_direct_level_ok
    )
    if opportunistic_fungus.get("enabled"):
        evaluated_axes["opportunistic_fungus_convergence"] = {
            "precise_identity": exact,
            "event_aligned_same_candidate_positive_count": len(
                event_aligned_direct_rows
            ),
            "direct_hospital_level_acceptable": fungus_direct_level_ok,
            "reproducible_mngs": reproducible,
            "cross_molecule_mngs": cross_molecule,
            "possible_only": True,
        }
    if opportunistic_fungus_route:
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": "opportunistic_fungus_reproducible_mngs_with_event_aligned_direct_detection",
            "reporting_role": OPPORTUNISTIC_FUNGUS_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [
                "Convergent exact-organism mNGS and event-aligned hospital evidence support Possible only; invasive disease is not proven.",
                "A respiratory fungal culture may represent colonization, so independent invasive-site or pathology evidence remains necessary for Picked.",
            ],
        }

    aspergillus = policy.get("aspergillus_group_possible") or {}
    taxonomy = candidate.get("taxonomy_profile") or {}
    organism_key = canonical_key(candidate.get("organism_name"))
    group_identity = bool(
        str(taxonomy.get("taxonomic_rank") or "") == "genus_or_group_label"
        or organism_key in set(aspergillus.get("group_organism_keys") or [])
    )
    aligned_rows = _event_aligned_or_date_overlap_rows(candidate, aspergillus)
    target_terms = tuple(
        str(value).casefold()
        for value in aspergillus.get("accepted_mycology_target_terms") or []
    )
    specimen_terms = tuple(
        str(value).casefold()
        for value in aspergillus.get("accepted_mycology_specimen_terms") or []
    )
    mycology_rows = [
        row for row in aligned_rows
        if any(term in _recursive_text(row) for term in target_terms)
        and any(
            term in str(row.get("specimen_type") or "").casefold()
            for term in specimen_terms
        )
    ]
    aspergillus_host_rows = _accepted_host_risk_rows(candidate, aspergillus)
    aspergillus_image_support = _bool(_axes(candidate).get("exact_image_context"))
    direct_level_ok = bool(
        direct_level is not None
        and direct_level <= int(aspergillus.get("maximum_direct_hospital_level", 3))
    )
    aspergillus_route = bool(
        aspergillus.get("enabled")
        and family == aspergillus.get("family", "mold_or_opportunistic_fungus")
        and group_identity
        and mycology_rows
        and direct_level_ok
        and (aspergillus_host_rows or aspergillus_image_support)
    )
    if aspergillus.get("enabled"):
        evaluated_axes["aspergillus_group_evidence"] = {
            "group_identity": group_identity,
            "mycology_rows": mycology_rows,
            "host_risk_rows": aspergillus_host_rows,
            "exact_image_context": aspergillus_image_support,
            "direct_hospital_level_acceptable": direct_level_ok,
            "species_resolution": "species_unresolved" if group_identity else "not_applicable",
            "missing_evidence_is_negative": False,
        }
    if aspergillus_route:
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": "aspergillus_group_mycology_plus_host_or_image_convergence",
            "reporting_role": ASPERGILLUS_GROUP_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [
                "The evidence supports possible invasive Aspergillus disease at group level; species remains unresolved.",
                "A positive galactomannan result is mycological support and is not equivalent to species culture or tissue invasion.",
                "Date-only timestamps are treated as intervals when they overlap the event window, not as exact collection times.",
            ],
        }

    cross = policy.get("strong_cross_molecule_possible") or {}
    cross_route = (
        family in set(cross.get("families") or [])
        and exact
        and context in set(cross.get("compatible_specimen_contexts") or [])
        and rank is not None
        and rank <= int(cross.get("maximum_rank", 4))
        and positive_tests >= int(cross.get("minimum_positive_tests", 3))
        and (reproducible or not cross.get("require_reproducibility", True))
        and (cross_molecule or not cross.get("require_cross_molecule", True))
    )
    if cross_route:
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": "strong_repeated_cross_molecule_signal_despite_source_caution",
            "reporting_role": STRONG_CROSS_MOLECULE_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [
                "Taxonomy or source caution blocks direct Picked but does not erase strong repeated DNA/RNA evidence.",
                "Missing culture or targeted PCR is unknown and not negative.",
            ],
        }

    aspiration = policy.get("aspiration_possible") or {}
    signals = [
        row
        for row in (candidate.get("analytical_profile") or {}).get(
            "per_test_signals"
        ) or []
        if isinstance(row, dict) and float(row.get("reads") or 0) > 0
    ]
    all_positive_tests_present = bool(signals) and len(signals) == positive_tests
    all_positive_tests_fully_evaluable = all_positive_tests_present and all(
        str(row.get("qc_status") or "") == "evaluable" for row in signals
    )
    aspiration_route = (
        family == aspiration.get("family")
        and exact
        and context in set(aspiration.get("compatible_specimen_contexts") or [])
        and str(candidate.get("history_policy_family") or "")
        == aspiration.get("required_history_policy_family")
        and aspiration.get("required_history_rule_id")
        in set(candidate.get("history_rule_ids") or [])
        and str(candidate.get("history_adjusted_route") or "")
        == aspiration.get("required_history_adjusted_route")
        and str(candidate.get("clinical_decision") or "")
        == aspiration.get("required_clinical_decision")
        and rank is not None
        and rank <= int(aspiration.get("maximum_rank", 3))
        and positive_tests >= int(aspiration.get("minimum_positive_tests", 3))
        and (reproducible or not aspiration.get("require_reproducibility", True))
        and (cross_molecule or not aspiration.get("require_cross_molecule", True))
        and (
            all_positive_tests_fully_evaluable
            or not aspiration.get(
                "require_all_positive_tests_fully_evaluable", True
            )
        )
    )
    if aspiration_route:
        evaluated_axes["aspiration_history_bridge"] = {
            "history_policy_family": candidate.get("history_policy_family"),
            "history_adjusted_route": candidate.get("history_adjusted_route"),
            "history_rule_ids": candidate.get("history_rule_ids") or [],
            "history_evidence_already_consumed": True,
            "all_positive_tests_fully_evaluable": (
                all_positive_tests_fully_evaluable
            ),
        }
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": "aspiration_history_priority_with_complete_cross_molecule_signal",
            "reporting_role": ASPIRATION_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [
                "Aspiration history was consumed upstream to establish Priority and is not counted again as organism proof.",
                "This route supports Possible only; independent same-organism or invasive evidence is still required for Picked.",
            ],
        }
    return None


def possible_pathogen_decision(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any]:
    axes = _axes(candidate)
    family = str(
        axes.get("taxonomy_family")
        or (candidate.get("taxonomy_profile") or {}).get("primary_rule_family")
        or "unmapped_or_uncertain"
    )
    decision = str(candidate.get("clinical_decision") or "")
    context = str(axes.get("specimen_context") or specimen_context(candidate))
    exact = _bool(axes.get("precise_identity"))
    if "precise_identity" not in axes:
        exact = precise_identity(candidate, policy)
    rank = as_int(axes.get("best_rank_in_retained_universe"))
    positive_tests = as_int(axes.get("selected_positive_test_count")) or 0
    reproducible = _bool(axes.get("reproducibility_axis"))
    cross_molecule = _bool(axes.get("cross_molecule_selected"))
    direct_level = as_int(axes.get("direct_hospital_level"))
    direct_rows = _event_aligned_direct_rows(candidate)
    sterile_rows = _sterile_rows(candidate, policy)
    unique_sterile_ids = sorted({
        str(row.get("observation_id")) for row in sterile_rows
        if row.get("observation_id")
    })
    thresholds = policy["thresholds"]
    compatible = context in set(policy["compatible_specimen_contexts"])

    evaluated_axes = {
        "taxonomy_family": family,
        "clinical_decision_unchanged": decision,
        "specimen_context": context,
        "precise_identity": exact,
        "best_rank_in_retained_universe": rank,
        "selected_positive_test_count": positive_tests,
        "reproducibility_axis": reproducible,
        "cross_molecule_selected": cross_molecule,
        "direct_hospital_level": direct_level,
        "event_aligned_direct_positive_count": len(direct_rows),
        "repeated_sterile_positive_count": len(unique_sterile_ids),
        "repeated_sterile_observation_ids": unique_sterile_ids,
        "repeated_sterile_independence_verified": False,
    }

    if decision == "picked_shadow":
        return {
            "outcome": "preserved_strict_picked",
            "route": "existing_clinical_picked",
            "reporting_role": STRICT_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [],
        }

    recall_first = _recall_first_v2_decision(
        candidate, policy, evaluated_axes
    )
    if recall_first is not None:
        return recall_first

    special = family in set(policy["special_family_guardrails"])
    if special:
        return {
            "outcome": "retained_in_original_review_tier",
            "route": "family_specific_guardrail_required",
            "reporting_role": NOT_SELECTED_ROLE,
            "selected_for_complete_report": False,
            "axes": evaluated_axes,
            "blockers": ["generic possible-pathogen routes are disabled for this family"],
            "cautions": [],
        }

    repeated_sterile = (
        exact
        and len(unique_sterile_ids) >= thresholds["repeated_sterile_minimum_observations"]
    )
    if repeated_sterile:
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": "repeated_event_aligned_sterile_site_detection",
            "reporting_role": SYSTEMIC_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [
                "Distinct observation IDs do not prove independently collected blood sets; verify source collection episodes.",
                "Contaminant-prone organisms still require clinical correlation.",
            ],
        }

    analytical_route = (
        decision == "review_high_priority"
        and family in set(policy["analytical_possible_families"])
        and exact and compatible and reproducible
        and rank is not None and rank <= thresholds["analytical_max_rank"]
        and positive_tests >= thresholds["analytical_minimum_positive_tests"]
    )
    direct_level_ok = (
        direct_level is not None
        and direct_level <= thresholds["direct_max_hospital_level"]
    )
    direct_without_level_ok = (
        direct_level is None
        and policy.get("allow_event_aligned_direct_positive_without_level", False)
    )
    direct_route = (
        decision == "review_high_priority"
        and family in set(policy["direct_possible_families"])
        and exact and compatible and reproducible
        and rank is not None and rank <= thresholds["direct_max_rank"]
        and positive_tests >= thresholds["direct_minimum_positive_tests"]
        and bool(direct_rows)
        and (direct_level_ok or direct_without_level_ok)
    )
    if analytical_route or direct_route:
        route = (
            "event_aligned_direct_detection_or_level3_culture"
            if direct_route else "reproducible_top_ranked_typical_or_hospital_bacterium"
        )
        cautions = []
        if direct_route and direct_level is None:
            cautions.append(
                "The candidate-scoped direct observation is positive and event aligned, but the generic hospital level is unavailable."
            )
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": route,
            "reporting_role": ANALYTICAL_DIRECT_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": cautions,
        }

    low_specificity_route = (
        decision in {"review_context_needed", "review_high_priority"}
        and family in set(policy["low_specificity_possible_families"])
        and exact and context == "lower_respiratory"
        and reproducible and cross_molecule
        and rank is not None and rank <= thresholds["low_specificity_max_rank"]
        and positive_tests >= thresholds["low_specificity_minimum_positive_tests"]
    )
    if low_specificity_route:
        return {
            "outcome": "selected_as_possible_pathogen",
            "route": "reproducible_cross_molecule_top_ranked_low_specificity_signal",
            "reporting_role": LOW_SPECIFICITY_ROLE,
            "selected_for_complete_report": True,
            "axes": evaluated_axes,
            "blockers": [],
            "cautions": [
                "This is a low-specificity respiratory signal and not a strict causative claim."
            ],
        }

    return {
        "outcome": "retained_in_original_review_tier",
        "route": "possible_pathogen_threshold_not_met",
        "reporting_role": NOT_SELECTED_ROLE,
        "selected_for_complete_report": False,
        "axes": evaluated_axes,
        "blockers": ["available generic axes do not meet a possible-pathogen route"],
        "cautions": [],
    }


def annotate_candidate(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any]:
    output = copy.deepcopy(candidate)
    immutable = {
        key: copy.deepcopy(output.get(key))
        for key in ("clinical_decision", "clinical_level", "formal_pick_allowed", "selection_role")
    }
    gate = possible_pathogen_decision(output, policy)
    output["possible_pathogen_gate"] = gate
    output["possible_reporting_role"] = gate["reporting_role"]
    output["selected_for_complete_report"] = gate["selected_for_complete_report"]
    output["resistance_interpretation"] = resistance_interpretation(output, policy)
    if policy.get("enable_recall_first_v2"):
        if gate["reporting_role"] == STRICT_ROLE:
            output["final_reporting_tier"] = "Picked"
        elif gate["reporting_role"] == NOT_SELECTED_ROLE:
            output["final_reporting_tier"] = "Context"
        else:
            output["final_reporting_tier"] = "Possible"
    for key, value in immutable.items():
        if output.get(key) != value:
            raise ValueError(f"Possible-pathogen stage changed immutable field {key}")
    return output


def _demote_from_primary_report(
    candidate: dict[str, Any], *, route: str, reason: str
) -> None:
    gate = candidate["possible_pathogen_gate"]
    gate["pre_compaction_decision"] = {
        "outcome": gate.get("outcome"),
        "route": gate.get("route"),
        "reporting_role": gate.get("reporting_role"),
        "selected_for_complete_report": gate.get("selected_for_complete_report"),
    }
    gate.update({
        "outcome": "retained_in_context_after_primary_report_compaction",
        "route": route,
        "reporting_role": NOT_SELECTED_ROLE,
        "selected_for_complete_report": False,
        "blockers": [reason],
        "cautions": [
            "The candidate remains visible in Context and is not deleted from the audit universe."
        ],
    })
    candidate["possible_reporting_role"] = NOT_SELECTED_ROLE
    candidate["selected_for_complete_report"] = False
    candidate["final_reporting_tier"] = "Context"


def apply_primary_report_compaction(
    candidates: list[dict[str, Any]], policy: dict[str, Any]
) -> list[dict[str, Any]]:
    config = policy.get("primary_report_compaction") or {}
    if not policy.get("enable_recall_first_v2") or not config.get("enabled"):
        return candidates
    analytical_route = str(config.get("analytical_only_route") or "")
    low_burden = config.get("low_normalized_burden") or {}
    inherited = config.get("inherited_taxonomy_without_support") or {}
    for candidate in candidates:
        gate = candidate.get("possible_pathogen_gate") or {}
        if (
            gate.get("reporting_role") != ANALYTICAL_DIRECT_ROLE
            or gate.get("route") != analytical_route
        ):
            continue
        axes = gate.get("axes") or {}
        signals = [
            row for row in (candidate.get("analytical_profile") or {}).get(
                "per_test_signals"
            ) or []
            if isinstance(row, dict) and float(row.get("reads") or 0) > 0
        ]
        rpms = [
            float(row["rpm_total"])
            for row in signals
            if row.get("rpm_total") is not None
        ]
        positive_tests = int(axes.get("selected_positive_test_count") or 0)
        all_signals_present = bool(signals) and len(signals) == positive_tests
        all_fully_evaluable = all_signals_present and all(
            str(row.get("qc_status") or "") == "evaluable"
            for row in signals
        )
        rpm_for_every_test = all_signals_present and len(rpms) == positive_tests
        low_burden_applies = (
            bool(rpms)
            and max(rpms) < float(
                low_burden.get(
                    "maximum_available_per_test_rpm_exclusive", 5.0
                )
            )
        )
        if low_burden.get("require_all_positive_tests_fully_evaluable", True):
            low_burden_applies = low_burden_applies and all_fully_evaluable
        if low_burden.get("require_rpm_for_every_positive_test", True):
            low_burden_applies = low_burden_applies and rpm_for_every_test
        if low_burden_applies:
            axes["primary_report_compaction"] = {
                "maximum_available_per_test_rpm": max(rpms),
                "positive_test_count": positive_tests,
                "all_positive_tests_fully_evaluable": all_fully_evaluable,
                "rpm_available_for_every_positive_test": rpm_for_every_test,
            }
            _demote_from_primary_report(
                candidate,
                route="analytical_only_low_normalized_burden_context",
                reason=(
                    "analytical-only signal has complete per-test QC and RPM coverage, "
                    "and every per-test RPM remains below the primary-report threshold"
                ),
            )
            continue

        taxonomy = candidate.get("taxonomy_profile") or {}
        mapping_status = str(taxonomy.get("mapping_status") or "")
        needs_host = bool(
            inherited.get(
                "require_recorded_host_context_when_not_exact", True
            )
        )
        host_context = _host_context_available(candidate, _axes(candidate))
        inherited_guardrail = (
            mapping_status != str(
                inherited.get("exact_mapping_status", "exact_species")
            )
            and needs_host
            and not host_context
        )
        if inherited_guardrail:
            axes["primary_report_compaction"] = {
                "taxonomy_mapping_status": mapping_status,
                "recorded_host_context_available": host_context,
            }
            _demote_from_primary_report(
                candidate,
                route="analytical_only_inherited_taxonomy_without_host_context",
                reason=(
                    "analytical-only primary-report route requires recorded host context "
                    "when the taxonomy family assignment is not exact-species mapped"
                ),
            )
    return candidates


def _fallback_features(
    candidate: dict[str, Any], policy: dict[str, Any],
    *, evidence_profile: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    config = policy.get("fallback_best_available") or {}
    gate = candidate.get("possible_pathogen_gate") or {}
    axes = gate.get("axes") or {}
    signals = _positive_mngs_signals(candidate)
    hospital_positive = _has_positive_hospital_evidence(candidate)
    if config.get("require_positive_evidence", False) and not (
        signals or hospital_positive
    ):
        return None
    if (
        not config.get("require_positive_evidence", False)
        and config.get("require_positive_test", True)
        and not signals
    ):
        return None
    family = str(axes.get("taxonomy_family") or "")
    taxonomy_profile = candidate.get("taxonomy_profile") or {}
    identity_exception_applied = bool(
        family in set(config.get("allow_curated_exact_group_labels_families") or [])
        and str(taxonomy_profile.get("mapping_status") or "") == "exact_species"
        and str(taxonomy_profile.get("rule_source") or "") == "curated_exact_species"
        and str(taxonomy_profile.get("taxonomic_rank") or "")
        == "group_or_species_label"
    )
    if (
        config.get("require_precise_identity", True)
        and not _bool(axes.get("precise_identity"))
        and not identity_exception_applied
    ):
        return None
    evaluable = [
        row for row in signals
        if str(row.get("qc_status") or "") == "evaluable"
    ]
    patient_mode = str(
        (evidence_profile or {}).get("decision_evidence_mode") or ""
    )
    partial_qc_rank_modes = set(
        config.get("partial_qc_rank_modes")
        or ["mngs_only", "mngs_plus_clinical_context"]
    )
    accepted_qc_statuses = {"evaluable"}
    if patient_mode in partial_qc_rank_modes:
        accepted_qc_statuses.update(
            str(value)
            for value in config.get("accepted_qc_statuses") or []
        )
    rankable = [
        row for row in signals
        if str(row.get("qc_status") or "") in accepted_qc_statuses
    ]
    best_evaluable_rank = min(
        (
            as_int(row.get("rank_in_retained_universe"))
            for row in evaluable
            if as_int(row.get("rank_in_retained_universe")) is not None
        ),
        default=999,
    )
    best_legacy_priority = min(
        (
            as_int(row.get("analytical_rank_priority"))
            for row in signals
            if as_int(row.get("analytical_rank_priority")) is not None
        ),
        default=999,
    )
    best_rankable_rank = min(
        (
            as_int(row.get("rank_in_retained_universe"))
            for row in rankable
            if as_int(row.get("rank_in_retained_universe")) is not None
        ),
        default=999,
    )
    max_evaluable_rpm = max(
        (
            float(row.get("rpm_total") or 0)
            for row in evaluable
            if row.get("rpm_total") is not None
        ),
        default=0.0,
    )
    decision_order = {
        "review_high_priority": 0,
        "review_context_needed": 1,
        "review_low_specificity": 2,
        "hold": 3,
    }
    context = str(axes.get("specimen_context") or "")
    tax_penalty = int(
        family in set(config.get("taxonomy_penalty_families") or [])
    )
    family_priority = int(
        (config.get("family_priority_order") or {}).get(
            family, config.get("default_family_priority", 0)
        )
    )
    decision_rank = decision_order.get(
        str(candidate.get("clinical_decision") or ""), 9
    )
    specimen_penalty = int(
        context not in set(policy.get("compatible_specimen_contexts") or [])
    )
    if config.get("family_aware_ordering_enabled", False):
        # Clinical tier remains first so family preference never promotes a
        # Context candidate over a High candidate.  QC and within-test rank
        # remain ahead of family: family is a tie-breaker for analytically
        # comparable candidates, not a substitute for the observed signal.
        core_rank = (
            decision_rank,
            specimen_penalty,
            0 if evaluable else 1 if rankable else 2,
            best_rankable_rank,
            family_priority,
            int(not _bool(axes.get("cross_molecule_selected"))),
            int(not _bool(axes.get("reproducibility_axis"))),
            as_int(axes.get("direct_hospital_level")) or 999,
            best_legacy_priority,
            tax_penalty,
        )
    else:
        core_rank = (
            decision_rank,
            0 if evaluable else 1 if rankable else 2,
            best_rankable_rank,
            int(not _bool(axes.get("cross_molecule_selected"))),
            int(not _bool(axes.get("reproducibility_axis"))),
            specimen_penalty,
            as_int(axes.get("direct_hospital_level")) or 999,
            best_legacy_priority,
            tax_penalty,
        )
    return {
        "core_rank": core_rank,
        "sort_key": (
            *core_rank,
            -int(axes.get("selected_positive_test_count") or 0),
            -max_evaluable_rpm,
            canonical_key(candidate.get("organism_name")),
        ),
        "best_evaluable_rank": (
            None if best_evaluable_rank == 999 else best_evaluable_rank
        ),
        "best_rankable_rank": (
            None if best_rankable_rank == 999 else best_rankable_rank
        ),
        "evaluable_positive_test_count": len(evaluable),
        "rankable_positive_test_count": len(rankable),
        "accepted_qc_statuses": sorted(accepted_qc_statuses),
        "hospital_positive_available": hospital_positive,
        "best_legacy_analytical_priority": (
            None if best_legacy_priority == 999 else best_legacy_priority
        ),
        "max_evaluable_rpm": max_evaluable_rpm,
        "taxonomy_penalty_applied": bool(tax_penalty),
        "family_priority": family_priority,
        "family_aware_ordering_enabled": bool(
            config.get("family_aware_ordering_enabled", False)
        ),
        "identity_exception_applied": identity_exception_applied,
    }


def _fallback_floor_failures(
    candidate: dict[str, Any],
    policy: dict[str, Any],
    features: dict[str, Any],
) -> list[str]:
    """Return answer-blind reasons why fallback should abstain."""
    config = policy.get("fallback_best_available") or {}
    floor = config.get("evidence_floor") or {}
    if not floor.get("abstain_enabled"):
        return []
    gate = candidate.get("possible_pathogen_gate") or {}
    axes = gate.get("axes") or {}
    family = str(axes.get("taxonomy_family") or "")
    decision = str(candidate.get("clinical_decision") or "")
    failures: list[str] = []
    if family in set(floor.get("excluded_families") or []):
        failures.append("fallback_family_excluded")
    accepted_decisions = set(floor.get("accepted_clinical_decisions") or [])
    if accepted_decisions and decision not in accepted_decisions:
        failures.append("fallback_clinical_tier_below_floor")
    if floor.get("require_compatible_specimen", False):
        accepted_contexts = set(
            floor.get("compatible_specimen_contexts")
            or policy.get("compatible_specimen_contexts")
            or []
        )
        if str(axes.get("specimen_context") or "") not in accepted_contexts:
            failures.append("fallback_specimen_context_not_compatible")

    hospital_positive = bool(features.get("hospital_positive_available"))
    rank = features.get("best_rankable_rank")
    if not hospital_positive:
        if floor.get("require_rankable_qc_for_mngs", False) and rank is None:
            failures.append("fallback_no_rankable_positive_mngs_test")
        maximum = (floor.get("maximum_rank_by_clinical_decision") or {}).get(decision)
        if maximum is not None and (rank is None or int(rank) > int(maximum)):
            failures.append("fallback_rank_below_evidence_floor")
    else:
        maximum_level = floor.get("hospital_only_maximum_direct_level")
        direct_level = as_int(axes.get("direct_hospital_level"))
        if direct_level is None:
            direct_level = 999
        if maximum_level is not None and direct_level > int(maximum_level):
            failures.append("fallback_hospital_evidence_below_direct_level_floor")

    if (
        family in set(floor.get("low_specificity_families") or [])
        and decision != "review_high_priority"
        and floor.get("low_specificity_context_requires_reproducibility", False)
        and not _bool(axes.get("reproducibility_axis"))
    ):
        failures.append("fallback_low_specificity_context_not_reproducible")
    return failures


def apply_best_available_fallback(
    candidates: list[dict[str, Any]], policy: dict[str, Any],
    *, evidence_profile: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    config = policy.get("fallback_best_available") or {}
    if not policy.get("enable_recall_first_v2") or not config.get("enabled"):
        return candidates
    if any(item.get("selected_for_complete_report") for item in candidates):
        return candidates
    profile = evidence_profile or decision_evidence_profile(candidates)
    ranked = []
    for item in candidates:
        features = _fallback_features(
            item, policy, evidence_profile=profile
        )
        gate = item.get("possible_pathogen_gate") or {}
        axes = gate.setdefault("axes", {})
        if features is None:
            axes["fallback_eligibility"] = {
                "eligible": False,
                "failures": ["base_fallback_requirements_not_met"],
            }
            continue
        failures = _fallback_floor_failures(item, policy, features)
        axes["fallback_eligibility"] = {
            "eligible": not failures,
            "failures": failures,
            "features": features,
        }
        if not failures:
            ranked.append((features["sort_key"], features, item))
    if not ranked:
        return candidates
    ranked.sort(key=lambda row: row[0])
    selected = [ranked[0]]
    maximum = int(config.get("maximum_candidates", 1))
    if (
        maximum > 1
        and len(ranked) > 1
        and config.get("include_runner_up_only_on_equal_core_rank", True)
        and ranked[1][1]["core_rank"] == ranked[0][1]["core_rank"]
    ):
        selected.append(ranked[1])
    for fallback_rank, (_, features, item) in enumerate(
        selected[:maximum], start=1
    ):
        gate = item["possible_pathogen_gate"]
        gate.update({
            "outcome": "selected_as_fallback_best_available",
            "route": (
                f"{profile['decision_evidence_mode']}_best_available_fallback"
            ),
            "reporting_role": FALLBACK_ROLE,
            "selected_for_complete_report": True,
            "blockers": [],
            "cautions": [
                "Low-confidence fallback because no candidate met Picked or ordinary Possible criteria.",
                "This role must not be interpreted as a formal causal claim.",
                "Missing corroborating evidence is unknown and is not treated as negative.",
            ],
        })
        gate["axes"]["fallback_rank"] = fallback_rank
        gate["axes"]["fallback_primary"] = fallback_rank == 1
        gate["axes"]["fallback_features"] = features
        item["possible_reporting_role"] = FALLBACK_ROLE
        item["selected_for_complete_report"] = True
        item["fallback_rank"] = fallback_rank
        item["fallback_primary"] = fallback_rank == 1
        item["final_reporting_confidence"] = "low"
        item["final_reporting_tier"] = "Fallback-Possible"
    return candidates


def apply_graceful_degradation(
    candidates: list[dict[str, Any]], policy: dict[str, Any],
    *, evidence_profile: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Keep sparse-input candidates reportable while applying confidence ceilings."""
    config = policy.get("graceful_degradation") or {}
    if not config.get("enabled"):
        return candidates
    patient_profile = evidence_profile or decision_evidence_profile(candidates)
    for item in candidates:
        candidate_mode = candidate_decision_evidence_mode(item)
        item["decision_evidence_profile"] = copy.deepcopy(patient_profile)
        item["candidate_decision_evidence_mode"] = candidate_mode
        gate = item["possible_pathogen_gate"]
        gate["axes"]["patient_decision_evidence_mode"] = patient_profile[
            "decision_evidence_mode"
        ]
        gate["axes"]["candidate_decision_evidence_mode"] = candidate_mode
        if candidate_mode != "mngs_only":
            continue
        if not item.get("selected_for_complete_report"):
            continue
        if item.get("possible_reporting_role") == FALLBACK_ROLE:
            continue
        gate["pre_graceful_degradation_decision"] = {
            "outcome": gate.get("outcome"),
            "route": gate.get("route"),
            "reporting_role": gate.get("reporting_role"),
            "selected_for_complete_report": gate.get(
                "selected_for_complete_report"
            ),
            "final_reporting_tier": item.get("final_reporting_tier"),
        }
        gate.update({
            "outcome": "selected_as_mngs_only_possible_pathogen",
            "route": "mngs_only_confidence_ceiling",
            "reporting_role": MNGS_ONLY_ROLE,
            "selected_for_complete_report": True,
            "blockers": [],
            "cautions": list(gate.get("cautions") or []) + [
                "Only positive mNGS evidence reached this candidate decision.",
                "Absence of Culture, PCR, history, or imaging support is unknown, not negative.",
                "The candidate is capped at Possible and requires clinical correlation before causal interpretation.",
            ],
        })
        item["possible_reporting_role"] = MNGS_ONLY_ROLE
        item["selected_for_complete_report"] = True
        item["final_reporting_tier"] = str(
            config.get("mngs_only_maximum_tier") or "Possible"
        )
        axes = gate["axes"]
        item["final_reporting_confidence"] = (
            "moderate"
            if _bool(axes.get("reproducibility_axis"))
            and _bool(axes.get("cross_molecule_selected"))
            else "low_to_moderate"
        )
    return candidates


def _csv_row(
    candidate: dict[str, Any], *, include_v2: bool = False
) -> dict[str, Any]:
    gate = candidate["possible_pathogen_gate"]
    axes = gate["axes"]
    row = {
        "patient_id": candidate.get("patient_id"),
        "organism_name": candidate.get("organism_name"),
        "taxonomy_family": axes.get("taxonomy_family"),
        "clinical_decision": candidate.get("clinical_decision"),
        "original_selection_role": candidate.get("selection_role"),
        "possible_reporting_role": candidate.get("possible_reporting_role"),
        "selected_for_complete_report": candidate.get("selected_for_complete_report"),
        "possible_outcome": gate.get("outcome"),
        "possible_route": gate.get("route"),
        "possible_blockers": "|".join(gate.get("blockers") or []),
        "possible_cautions": "|".join(gate.get("cautions") or []),
        "specimen_context": axes.get("specimen_context"),
        "precise_identity": axes.get("precise_identity"),
        "best_rank": axes.get("best_rank_in_retained_universe"),
        "selected_positive_test_count": axes.get("selected_positive_test_count"),
        "reproducibility_axis": axes.get("reproducibility_axis"),
        "cross_molecule_selected": axes.get("cross_molecule_selected"),
        "direct_hospital_level": axes.get("direct_hospital_level"),
        "event_aligned_direct_positive_count": axes.get("event_aligned_direct_positive_count"),
        "repeated_sterile_positive_count": axes.get("repeated_sterile_positive_count"),
        "identity_resolution": (
            candidate.get("taxonomy_profile") or {}
        ).get("identity_resolution"),
        "taxonomy_rule_source": (
            candidate.get("taxonomy_profile") or {}
        ).get("taxonomy_rule_source"),
        "resistance_target": (
            candidate.get("resistance_interpretation") or {}
        ).get("resistance_target"),
        "resistance_status": (
            candidate.get("resistance_interpretation") or {}
        ).get("status"),
    }
    if include_v2:
        profile = candidate.get("decision_evidence_profile") or {}
        row.update({
            "final_reporting_tier": candidate.get("final_reporting_tier"),
            "fallback_rank": axes.get("fallback_rank"),
            "fallback_primary": axes.get("fallback_primary"),
            "final_reporting_confidence": candidate.get(
                "final_reporting_confidence"
            ),
            "patient_decision_evidence_mode": profile.get(
                "decision_evidence_mode"
            ),
            "candidate_decision_evidence_mode": candidate.get(
                "candidate_decision_evidence_mode"
            ),
            "positive_or_contextual_sources": "|".join(
                profile.get("positive_or_contextual_sources") or []
            ),
        })
    return row


def _write_csv(
    path: Path,
    candidates: list[dict[str, Any]],
    *,
    include_v2: bool = False,
) -> None:
    rows = [
        _csv_row(candidate, include_v2=include_v2)
        for candidate in candidates
    ]
    fields = (
        list(_csv_row(candidates[0], include_v2=include_v2).keys())
        if candidates else [
        "patient_id", "organism_name", "taxonomy_family", "clinical_decision",
        "original_selection_role", "possible_reporting_role",
        "selected_for_complete_report", "possible_outcome", "possible_route",
        "possible_blockers", "possible_cautions", "specimen_context",
        "precise_identity", "best_rank", "selected_positive_test_count",
        "reproducibility_axis", "cross_molecule_selected", "direct_hospital_level",
        "event_aligned_direct_positive_count", "repeated_sterile_positive_count",
        "identity_resolution", "taxonomy_rule_source", "resistance_target",
        "resistance_status",
    ])
    if include_v2 and not candidates:
        fields.extend([
            "final_reporting_tier", "fallback_rank", "fallback_primary",
            "final_reporting_confidence", "patient_decision_evidence_mode",
            "candidate_decision_evidence_mode",
            "positive_or_contextual_sources",
        ])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(
    input_root: Path, output_dir: Path, *, policy_path: Path = DEFAULT_POLICY,
    patients: set[int] | None = None,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    policy, policy_parent_path = load_policy(policy_path)
    paths = sorted(
        (input_root / "patient_outputs").glob(
            "NGS_patient_*_test_aware_clinical_shadow.json"
        ),
        key=patient_number,
    )
    if patients is not None:
        paths = [path for path in paths if patient_number(path) in patients]
    if not paths:
        raise ValueError("No clinical V4 patient outputs found")

    output_dir.mkdir(parents=True, exist_ok=True)
    patient_dir = output_dir / "patient_outputs"
    patient_dir.mkdir()
    all_candidates: list[dict[str, Any]] = []
    abstention_rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()

    for source_path in paths:
        source = read_json(source_path)
        original = source.get("all_forwarded_candidates") or []
        candidates = [annotate_candidate(item, policy) for item in original]
        candidates = apply_primary_report_compaction(candidates, policy)
        evidence_profile = decision_evidence_profile(candidates)
        candidates = apply_graceful_degradation(
            candidates, policy, evidence_profile=evidence_profile
        )
        candidates = apply_best_available_fallback(
            candidates, policy, evidence_profile=evidence_profile
        )
        if policy.get("graceful_degradation", {}).get("enabled"):
            for item in candidates:
                item.setdefault(
                    "decision_evidence_profile", copy.deepcopy(evidence_profile)
                )
                item.setdefault(
                    "candidate_decision_evidence_mode",
                    candidate_decision_evidence_mode(item),
                )
        if len(candidates) != len(original):
            raise ValueError(f"Candidate loss in {source_path}")
        all_candidates.extend(candidates)
        for item in candidates:
            gate = item["possible_pathogen_gate"]
            counts[f"clinical:{item.get('clinical_decision')}"] += 1
            counts[f"reporting_role:{gate['reporting_role']}"] += 1
            counts[f"route:{gate['route']}"] += 1

        patient = patient_number(source_path)
        strict = [
            item for item in candidates
            if item["possible_reporting_role"] == STRICT_ROLE
        ]
        possible = [
            item for item in candidates
            if item["possible_reporting_role"] not in {STRICT_ROLE, NOT_SELECTED_ROLE}
        ]
        complete = [item for item in candidates if item["selected_for_complete_report"]]
        abstention = {
            "abstained": not bool(complete),
            "reason": (
                "no_candidate_met_picked_possible_or_fallback_evidence_floor"
                if not complete else "not_applicable"
            ),
            "candidate_count": len(candidates),
            "fallback_eligible_candidate_count": sum(
                bool(
                    (((item.get("possible_pathogen_gate") or {}).get("axes") or {})
                     .get("fallback_eligibility") or {}).get("eligible")
                )
                for item in candidates
            ),
        }
        if abstention["abstained"]:
            counts["patient:abstained"] += 1
            abstention_rows.append({"patient_id": patient, **abstention})
        payload = {
            "schema_version": (
                "test_aware_unified_reporting_shadow.v2"
                if policy.get("enable_recall_first_v2") else SCHEMA_VERSION
            ),
            "policy_id": policy["policy_id"],
            "answer_blind": True,
            "patient_id": str(patient),
            "candidate_count": len(candidates),
            "decision_evidence_profile": evidence_profile,
            "strict_picked": strict,
            "possible_pathogens": possible,
            "complete_report": complete,
            "abstention": abstention,
            "remaining_high": [
                item for item in candidates
                if item["clinical_decision"] == "review_high_priority"
                and not item["selected_for_complete_report"]
            ],
            "remaining_context": [
                item for item in candidates
                if item["clinical_decision"] == "review_context_needed"
                and not item["selected_for_complete_report"]
            ],
            "all_forwarded_candidates": candidates,
            "source": {
                "clinical_v4_file": str(source_path.resolve()),
                "clinical_v4_sha256": sha256_file(source_path),
            },
            "constraints": policy.get("constraints") or [],
        }
        write_json(
            patient_dir / f"NGS_patient_{patient}_test_aware_possible_pathogen_shadow.json",
            payload,
        )

    include_v2 = bool(policy.get("enable_recall_first_v2"))
    _write_csv(
        output_dir / "possible_pathogen_decisions.csv",
        all_candidates,
        include_v2=include_v2,
    )
    _write_csv(
        output_dir / "possible_pathogens.csv",
        [item for item in all_candidates
         if item["possible_reporting_role"] not in {STRICT_ROLE, NOT_SELECTED_ROLE}],
        include_v2=include_v2,
    )
    _write_csv(
        output_dir / "complete_report.csv",
        [item for item in all_candidates if item["selected_for_complete_report"]],
        include_v2=include_v2,
    )
    with (output_dir / "abstentions.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        fields = [
            "patient_id", "abstained", "reason", "candidate_count",
            "fallback_eligible_candidate_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(abstention_rows)
    summary = {
        "schema_version": (
            "test_aware_unified_reporting_shadow.v2"
            if policy.get("enable_recall_first_v2") else SCHEMA_VERSION
        ),
        "policy_id": policy["policy_id"],
        "answer_blind": True,
        "patient_count": len(paths),
        "candidate_count": len(all_candidates),
        "strict_picked_count": sum(
            item["possible_reporting_role"] == STRICT_ROLE
            for item in all_candidates
        ),
        "possible_pathogen_count": sum(
            item["possible_reporting_role"] not in {STRICT_ROLE, NOT_SELECTED_ROLE}
            for item in all_candidates
        ),
        "complete_report_count": sum(
            bool(item["selected_for_complete_report"]) for item in all_candidates
        ),
        "counts": dict(sorted(counts.items())),
        "source_clinical_v4_root": str(input_root.resolve()),
        "source_summary_sha256": (
            sha256_file(input_root / "summary.json")
            if (input_root / "summary.json").is_file() else None
        ),
        "policy_file": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "policy_extends": (
            str(policy_parent_path.resolve()) if policy_parent_path else None
        ),
        "policy_base_sha256": (
            sha256_file(policy_parent_path) if policy_parent_path else None
        ),
        "abstained_patient_count": len(abstention_rows),
        "selected_patients": sorted(patients) if patients is not None else None,
        "constraints": policy.get("constraints") or [],
    }
    if include_v2:
        summary["final_reporting_tier_counts"] = dict(Counter(
            str(item.get("final_reporting_tier") or "missing")
            for item in all_candidates
        ))
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
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

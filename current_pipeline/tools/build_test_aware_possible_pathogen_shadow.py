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
from collections import Counter
from pathlib import Path
from typing import Any

from tools.build_index_event_clinical_timeline_shadow import sha256_file
from tools.pathogen_normalization import canonical_key
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


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().casefold() in {"1", "true", "yes"}


def _axes(candidate: dict[str, Any]) -> dict[str, Any]:
    return (candidate.get("promotion_gate") or {}).get("axes") or {}


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
    host_context = _host_context_available(candidate, _axes(candidate))
    provisional_host_config = (
        policy.get("provisional_phenotype_host_context") or {}
    )
    provisional_host = _provisional_phenotype_host_context(
        candidate, provisional_host_config
    )
    reactivation = _reactivation_status(candidate)
    evaluated_axes["recorded_host_context_available"] = host_context
    if provisional_host_config.get("enabled"):
        evaluated_axes["provisional_phenotype_host_context"] = provisional_host
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

    pjp = policy.get("pjp_possible") or {}
    is_pjp = (
        canonical_key(candidate.get("organism_name"))
        == str(pjp.get("canonical_organism_key") or "")
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
            host_context
            or (
                pjp.get("allow_provisional_phenotype_host_context", False)
                and provisional_host["available"]
            )
            or not pjp.get("allow_recorded_host_context", True)
        )
    )
    if pjp_route:
        provisional_host_used = bool(
            not host_context
            and pjp.get("allow_provisional_phenotype_host_context", False)
            and provisional_host["available"]
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
                    else "Recorded host context supports Possible only; dose-unverified exposure cannot create Picked."
                ),
                "Missing PCR, beta-D-glucan, staining, or tissue evidence is unknown and not negative.",
            ],
        }

    virus = policy.get("reactivation_possible") or {}
    virus_context_ok = context in set(
        virus.get("compatible_specimen_contexts") or []
    )
    virus_minimum_tests = int(
        virus.get("explicit_current_reactivation_minimum_positive_tests", 2)
        if reactivation["explicit_current_reactivation"]
        else virus.get("cross_molecule_minimum_positive_tests", 2)
    )
    strong_virus = (
        cross_molecule
        and reproducible
        and rank is not None
        and rank <= int(virus.get("cross_molecule_maximum_rank", 3))
        and positive_tests >= virus_minimum_tests
    )
    virus_clinical_support = bool(
        host_context
        or (
            virus.get("allow_provisional_phenotype_host_context", False)
            and provisional_host["available"]
        )
        or reactivation["explicit_current_reactivation"]
        or reactivation["quantitative_or_tissue_evidence"]
        or reactivation["compatible_pulmonary_syndrome"]
        or bool(_event_aligned_direct_rows(candidate))
        or _bool(_axes(candidate).get("exact_image_context"))
    )
    if virus.get("require_clinical_support_for_cross_molecule", False):
        strong_virus = strong_virus and virus_clinical_support
    host_supported_single = (
        virus.get("enable_host_supported_single_test", True)
        and host_context
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
        and bool(_event_aligned_direct_rows(candidate))
        and rank is not None
        and rank <= int(virus.get("direct_supported_single_test_maximum_rank", 1))
        and positive_tests >= 1
    )
    virus_route = (
        family == virus.get("family")
        and exact
        and virus_context_ok
        and (strong_virus or host_supported_single or direct_supported_single)
    )
    if virus_route:
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
                    else "reactivation_family_top_rank_with_recorded_host_context"
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
    candidate: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any] | None:
    config = policy.get("fallback_best_available") or {}
    gate = candidate.get("possible_pathogen_gate") or {}
    axes = gate.get("axes") or {}
    profile = candidate.get("analytical_profile") or {}
    signals = [
        row for row in profile.get("per_test_signals") or []
        if isinstance(row, dict) and float(row.get("reads") or 0) > 0
    ]
    if config.get("require_positive_test", True) and not signals:
        return None
    if config.get("require_precise_identity", True) and not _bool(
        axes.get("precise_identity")
    ):
        return None
    evaluable = [
        row for row in signals
        if str(row.get("qc_status") or "") == "evaluable"
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
    family = str(axes.get("taxonomy_family") or "")
    tax_penalty = int(
        family in set(config.get("taxonomy_penalty_families") or [])
    )
    core_rank = (
        decision_order.get(str(candidate.get("clinical_decision") or ""), 9),
        int(not bool(evaluable)),
        best_evaluable_rank,
        int(not _bool(axes.get("cross_molecule_selected"))),
        int(not _bool(axes.get("reproducibility_axis"))),
        int(context not in set(policy.get("compatible_specimen_contexts") or [])),
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
        "evaluable_positive_test_count": len(evaluable),
        "best_legacy_analytical_priority": (
            None if best_legacy_priority == 999 else best_legacy_priority
        ),
        "max_evaluable_rpm": max_evaluable_rpm,
        "taxonomy_penalty_applied": bool(tax_penalty),
    }


def apply_best_available_fallback(
    candidates: list[dict[str, Any]], policy: dict[str, Any]
) -> list[dict[str, Any]]:
    config = policy.get("fallback_best_available") or {}
    if not policy.get("enable_recall_first_v2") or not config.get("enabled"):
        return candidates
    if any(item.get("selected_for_complete_report") for item in candidates):
        return candidates
    ranked = []
    for item in candidates:
        features = _fallback_features(item, policy)
        if features is not None:
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
            "route": "no_picked_or_possible_best_available_fallback",
            "reporting_role": FALLBACK_ROLE,
            "selected_for_complete_report": True,
            "blockers": [],
            "cautions": [
                "Low-confidence fallback because no candidate met Picked or ordinary Possible criteria.",
                "This role must not be interpreted as a formal causal claim.",
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
    }
    if include_v2:
        row.update({
            "final_reporting_tier": candidate.get("final_reporting_tier"),
            "fallback_rank": axes.get("fallback_rank"),
            "fallback_primary": axes.get("fallback_primary"),
            "final_reporting_confidence": candidate.get(
                "final_reporting_confidence"
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
    ])
    if include_v2 and not candidates:
        fields.extend([
            "final_reporting_tier", "fallback_rank", "fallback_primary",
            "final_reporting_confidence",
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
    policy = read_json(policy_path)
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
    counts: Counter[str] = Counter()

    for source_path in paths:
        source = read_json(source_path)
        original = source.get("all_forwarded_candidates") or []
        candidates = [annotate_candidate(item, policy) for item in original]
        candidates = apply_primary_report_compaction(candidates, policy)
        candidates = apply_best_available_fallback(candidates, policy)
        if len(candidates) != len(original):
            raise ValueError(f"Candidate loss in {source_path}")
        all_candidates.extend(candidates)
        for item in candidates:
            gate = item["possible_pathogen_gate"]
            counts[f"clinical:{item.get('clinical_decision')}"] += 1
            counts[f"reporting_role:{gate['reporting_role']}"] += 1
            counts[f"route:{gate['route']}"] += 1

        patient = patient_number(source_path)
        strict = [item for item in candidates if item["clinical_decision"] == "picked_shadow"]
        possible = [
            item for item in candidates
            if item["possible_reporting_role"] not in {STRICT_ROLE, NOT_SELECTED_ROLE}
        ]
        complete = [item for item in candidates if item["selected_for_complete_report"]]
        payload = {
            "schema_version": (
                "test_aware_unified_reporting_shadow.v2"
                if policy.get("enable_recall_first_v2") else SCHEMA_VERSION
            ),
            "policy_id": policy["policy_id"],
            "answer_blind": True,
            "patient_id": str(patient),
            "candidate_count": len(candidates),
            "strict_picked": strict,
            "possible_pathogens": possible,
            "complete_report": complete,
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
            item["clinical_decision"] == "picked_shadow" for item in all_candidates
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

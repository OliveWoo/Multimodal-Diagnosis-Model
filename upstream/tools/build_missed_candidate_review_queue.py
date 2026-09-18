from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any, Sequence

from tools import deterministic_mngs_max_scorer as scorer
from tools import mngs_common as mngs
from tools import organism_taxonomy_classifier as organism_taxonomy
from tools import pathogen_normalization as pathogen_names
from tools import pathogen_rule_lists


DEFAULT_DETERMINISTIC_SUFFIX = "mNGS_max_deterministic_resp_commensal_dominance_guardrail_opt_chosen_full"
DEFAULT_OUTPUT_SUFFIX = "mNGS_missed_candidate_review_queue"


RARE_OPPORTUNISTIC_GNB_PREFIXES = (
    "achromobacter",
    "chryseobacterium",
    "delftia",
    "elizabethkingia",
    "ralstonia",
    "sphingomonas",
    "phytobacter",
    "xanthomonas",
    "cupriavidus",
    "comamonas",
    "pandoraea",
    "brevundimonas",
    "ochrobactrum",
    "alcaligenes",
    "acidovorax",
)

RARE_OPPORTUNISTIC_GNB_EXACT_NAMES = {
    "brucellaintermedia",
    "enterobactersoli",
    "pseudomonaszhaodongensis",
}

ATYPICAL_PNEUMONIA_PREFIXES = (
    "legionella",
    "chlamydophila",
    "chlamydia",
    "mycoplasma",
    "coxiella",
)

RESPIRATORY_ANAEROBE_OR_ASPIRATION_PREFIXES = (
    "prevotella",
    "fusobacterium",
    "bacteroides",
    "phocaeicola",
    "segatella",
    "parabacteroides",
    "veillonella",
    "parvimonas",
    "porphyromonas",
    "peptostreptococcus",
    "finegoldia",
    "actinomyces",
    "capnocytophaga",
    "leptotrichia",
    "moraxellaosloensis",
    "mycoplasmasalivarium",
    "streptococcustaonis",
    "clostridium",
    "clostridioides",
)

STRICT_ASPIRATION_ANAEROBE_PREFIXES = (
    "prevotella",
    "fusobacterium",
    "bacteroides",
    "phocaeicola",
    "segatella",
    "parabacteroides",
    "veillonella",
    "parvimonas",
    "porphyromonas",
    "peptostreptococcus",
    "finegoldia",
    "actinomyces",
    "clostridium",
    "clostridioides",
)

BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_PREFIXES = (
    "gemella",
    "lacticaseibacillus",
    "ligilactobacillus",
    "limosilactobacillus",
    "rothia",
    "lactobacillus",
)

ATYPICAL_ORAL_ASSOCIATED_PREFIXES = (
    "capnocytophaga",
    "leptotrichia",
    "mycoplasmasalivarium",
    "mycoplasmaorale",
    "abiotrophia",
)

BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_EXACT_NAMES = {
    "moraxellaosloensis",
    "neisseriasubflava",
    "neisseriaflavescens",
    "neisseriamucosa",
    "neisseriasicca",
    "neisseriacinerea",
    "neisseriaperflava",
    "neisseriaelongata",
    "streptococcustaonis",
}

ATYPICAL_ORAL_ASSOCIATED_EXACT_NAMES = {
    "kingellaoralis",
    "campylobacterconcisus",
    "mycoplasmasalivarium",
}

HIGH_RISK_MOLD_OR_DIMORPHIC_PREFIXES = (
    "aspergillus",
    "mucorales",
    "rhizopus",
    "mucor",
    "rhizomucor",
    "lichtheimia",
    "cunninghamella",
    "scedosporium",
    "lomentospora",
    "fusarium",
    "talaromyces",
    "histoplasma",
    "cryptococcus",
)

YEAST_OPPORTUNISTIC_PREFIXES = (
    "candida",
    "nakaseomyces",
    "pichia",
    "trichosporon",
    "saccharomyces",
)

SPECIAL_BACKGROUND_PATHOGEN_PREFIXES = (
    "corynebacteriumstriatum",
    "corynebacteriumjeikeium",
    "cutibacteriumacnes",
)

ATYPICAL_ENTEROBACTERALES_PREFIXES = (
    "enterobacter",
    "citrobacter",
    "serratia",
    "morganella",
    "providencia",
    "raoultella",
    "kluyvera",
    "cronobacter",
    "hafnia",
    "pluralibacter",
)

OPPORTUNISTIC_HERPESVIRUS_EXACT_NAMES = {
    "humangammaherpesvirus4",  # EBV
    "humanbetaherpesvirus6",
    "humanbetaherpesvirus6a",
    "humanbetaherpesvirus6b",
    "humanbetaherpesvirus7",
}

REACTIVATION_HERPESVIRUS_EXACT_NAMES = {
    "humanbetaherpesvirus5",  # CMV
    "humanalphaherpesvirus1",  # HSV-1
    "humanalphaherpesvirus2",  # HSV-2
    "humanalphaherpesvirus3",  # VZV
}

HERPES_R4_NONDOMINANT_HOST_ONLY_REVIEW_RULE = "R-S4-HERPES-R4-NONDOMINANT-HOST-ONLY-REVIEW"
HERPES_HOST_R2R3_NONDOMINANT_REVIEW_RULE = "R-S4-HERPES-HOST-R2R3-NONDOMINANT-REVIEW"


def load_json(path: Path) -> Any:
    return scorer.load_json(path)


def write_json(path: Path, payload: Any) -> None:
    scorer.write_json(path, payload)


def output_path_for(patient_dir: Path, output_suffix: str, artifact_root: Path | None = None) -> Path:
    base = mngs.extract_patient_identifier(patient_dir)
    if artifact_root is not None:
        return artifact_root / patient_dir.name / mngs.SUMMARY_OUTPUT_DIR_NAME / f"{base}_{output_suffix}.json"
    return patient_dir / mngs.SUMMARY_OUTPUT_DIR_NAME / f"{base}_{output_suffix}.json"


def norm(value: Any) -> str:
    return mngs.normalize_organism_name(value)


def canonical_key(value: Any) -> str:
    return pathogen_names.canonical_key(value) or norm(value)


def is_rare_opportunistic_gnb_name(value: Any) -> bool:
    normalized = norm(value)
    return normalized in RARE_OPPORTUNISTIC_GNB_EXACT_NAMES or any(
        normalized.startswith(prefix) for prefix in RARE_OPPORTUNISTIC_GNB_PREFIXES
    )


def is_non_aureus_staphylococcus_name(value: Any) -> bool:
    normalized = norm(value)
    return normalized.startswith("staphylococcus") and not normalized.startswith("staphylococcusaureus")


def starts_with_any_name(value: Any, prefixes: Sequence[str]) -> bool:
    normalized = norm(value)
    return any(normalized.startswith(prefix) for prefix in prefixes)


def is_broad_or_atypical_oral_flora_name(value: Any) -> bool:
    normalized = norm(value)
    return (
        normalized in BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_EXACT_NAMES
        or normalized in ATYPICAL_ORAL_ASSOCIATED_EXACT_NAMES
        or starts_with_any_name(value, BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_PREFIXES)
        or starts_with_any_name(value, ATYPICAL_ORAL_ASSOCIATED_PREFIXES)
    )


def is_opportunistic_herpesvirus_name(value: Any) -> bool:
    return canonical_key(value) in OPPORTUNISTIC_HERPESVIRUS_EXACT_NAMES


def is_reactivation_herpesvirus_name(value: Any) -> bool:
    return canonical_key(value) in REACTIVATION_HERPESVIRUS_EXACT_NAMES


def is_host_vulnerable(ctx: dict[str, Any] | None) -> bool:
    if not isinstance(ctx, dict):
        return False
    return bool(
        ctx.get("host_level3_expansion")
        or str(ctx.get("host_vulnerability_tier") or "").upper() in {"V2", "V3"}
        or str(ctx.get("opportunistic_coverage_level") or "").upper() in {"O1", "O2", "O3"}
    )


def tier_rank(reads_tier: Any) -> int:
    order = {
        "R4_very_high": 4,
        "R3_high": 3,
        "R2_medium": 2,
        "R1_low": 1,
        "R0_trace": 0,
    }
    return order.get(str(reads_tier or ""), -1)


def is_common_or_background(candidate: dict[str, Any]) -> bool:
    name = candidate.get("organism_name")
    return bool(
        candidate.get("is_likely_colonizer_or_background")
        or mngs.is_background_pathogen_name(name)
        or mngs.is_oral_upper_airway_flora_name(name)
        or (
            mngs.is_yeast_like_background_name(name)
            and not scorer.is_rare_opportunistic_yeast_name(name)
        )
        or scorer.is_respiratory_commensal_flora_name(name)
        or scorer.is_coagulase_negative_staph_background_name(name)
    )


def clean_mngs_signal(candidate: dict[str, Any]) -> bool:
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads = mngs.to_float(candidate.get("reads"))
    rank_rule = str(candidate.get("rank_rule") or "").lower()
    if reads <= 0 or rank > 3:
        return False
    return rank <= 1 or "rk00" in rank_rule or "ntc=00" in rank_rule or "ntc" in rank_rule


def strong_lower_respiratory_aspiration_signal(candidate: dict[str, Any]) -> bool:
    name = candidate.get("organism_name")
    if not starts_with_any_name(name, STRICT_ASPIRATION_ANAEROBE_PREFIXES):
        return False
    if candidate.get("specimen_class") != "S2_lower_respiratory":
        return False
    if not clean_mngs_signal(candidate):
        return False
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads_rank = tier_rank(candidate.get("reads_tier"))
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    return rank <= 3 and (reads_rank >= 3 or percentile >= 0.85)


def strong_water_environment_watch_signal(candidate: dict[str, Any]) -> bool:
    name = candidate.get("organism_name")
    if not (scorer.is_water_environmental_gnb_name(name) or scorer.is_low_pulmonary_specificity_environmental_name(name)):
        return False
    if not clean_mngs_signal(candidate):
        return False
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads_rank = tier_rank(candidate.get("reads_tier"))
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    dominance = candidate.get("dominance_tier")
    specimen_cls = candidate.get("specimen_class")
    if dominance in {"D2_moderate", "D3_dominant"} or non_host_support_modules(candidate):
        return True
    if specimen_cls == "S2_lower_respiratory" and rank <= 3 and (reads_rank >= 3 or percentile >= 0.85):
        return True
    if specimen_cls == "S1_sterile_systemic" and rank <= 3 and (reads_rank >= 2 or percentile >= 0.75):
        return True
    return False


def support_modules(candidate: dict[str, Any]) -> list[str]:
    support = candidate.get("module_support_summary") or {}
    if not isinstance(support, dict):
        return []
    return sorted(key for key, value in support.items() if key != "mNGS" and value == "Support")


def non_host_support_modules(candidate: dict[str, Any]) -> list[str]:
    return [key for key in support_modules(candidate) if key != "host"]


def candida_related_respiratory_hospital_context(candidate: dict[str, Any]) -> bool:
    if not scorer.is_candida_or_generic_yeast_name(candidate.get("organism_name")):
        return False
    related = candidate.get("related_representative_hospital_support")
    if not isinstance(related, list) or not related:
        related = (candidate.get("key_evidence") or {}).get("related_representative_hospital_support")
    return isinstance(related, list) and bool(related)


def respiratory_candida_high_priority_context(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    if not scorer.is_candida_or_generic_yeast_name(candidate.get("organism_name")):
        return False
    if candidate.get("specimen_class") != "S2_lower_respiratory":
        return False
    if not is_host_vulnerable(ctx):
        return False
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads_rank = tier_rank(candidate.get("reads_tier"))
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    has_mngs_signal = rank <= 2 and (reads_rank >= 2 or percentile >= 0.90)
    has_hospital_context = (
        bool(non_host_support_modules(candidate))
        or bool((candidate.get("key_evidence") or {}).get("candida_invasive_hospital_support"))
        or bool((candidate.get("key_evidence") or {}).get("candida_related_invasive_hospital_context"))
        or candida_related_respiratory_hospital_context(candidate)
    )
    return has_mngs_signal and has_hospital_context


def herpes_r4_nondominant_host_only_review(candidate: dict[str, Any]) -> bool:
    rule_ids = {
        str(rule)
        for rule in [
            *((candidate.get("applied_rules") or []) if isinstance(candidate.get("applied_rules"), list) else []),
            (candidate.get("key_evidence") or {}).get("guardrail_rule", ""),
        ]
    }
    if HERPES_R4_NONDOMINANT_HOST_ONLY_REVIEW_RULE in rule_ids:
        return True
    return (
        is_reactivation_herpesvirus_name(candidate.get("organism_name"))
        and candidate.get("reads_tier") == "R4_very_high"
        and candidate.get("dominance_tier") in {"D0_not_top", "D1_low"}
        and not non_host_support_modules(candidate)
    )


def herpes_ranked_host_context_review(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    rule_ids = {
        str(rule)
        for rule in [
            *((candidate.get("applied_rules") or []) if isinstance(candidate.get("applied_rules"), list) else []),
            (candidate.get("key_evidence") or {}).get("guardrail_rule", ""),
        ]
    }
    if HERPES_HOST_R2R3_NONDOMINANT_REVIEW_RULE in rule_ids:
        return True
    return (
        is_reactivation_herpesvirus_name(candidate.get("organism_name"))
        and bool(ctx.get("host_level3_expansion"))
        and mngs.to_rank_int(candidate.get("rank_priority")) <= 2
        and candidate.get("reads_tier") in {"R2_medium", "R3_high"}
        and mngs.to_float(candidate.get("reads_percentile")) >= 0.85
        and candidate.get("dominance_tier") in {"D0_not_top", "D1_low"}
        and not non_host_support_modules(candidate)
    )


def review_reasons(candidate: dict[str, Any], ctx: dict[str, Any]) -> list[str]:
    if pathogen_rule_lists.is_impossible_infection_source(candidate.get("organism_name")):
        return []
    reasons: list[str] = []
    level = scorer.level_rank(candidate.get("integrated_causative_level"))
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads_tier = candidate.get("reads_tier")
    reads_rank = tier_rank(reads_tier)
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    dominance = candidate.get("dominance_tier")
    common_or_background = is_common_or_background(candidate)
    applied_rules = {str(rule) for rule in (candidate.get("applied_rules") or [])}

    if herpes_r4_nondominant_host_only_review(candidate):
        reasons.append("herpes_r4_nondominant_host_only_review_context")
    if herpes_ranked_host_context_review(candidate, ctx):
        reasons.append("herpes_ranked_host_context_review")
    if respiratory_candida_high_priority_context(candidate, ctx):
        reasons.append("respiratory_candida_high_priority_context")
    candida_profile = (candidate.get("key_evidence") or {}).get("candida_evidence_strength_profile")
    if isinstance(candida_profile, dict):
        best_candida_strength = str(candida_profile.get("best_strength") or "")
        if best_candida_strength == getattr(scorer, "CANDIDA_EVIDENCE_STRONG_INVASIVE", ""):
            reasons.append("candida_strong_invasive_evidence")
        elif best_candida_strength == getattr(scorer, "CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC", ""):
            reasons.append("candida_intermediate_systemic_or_invasive_clue")
        elif best_candida_strength == getattr(scorer, "CANDIDA_EVIDENCE_WEAK_NONINVASIVE", ""):
            reasons.append("candida_weak_noninvasive_or_colonization_evidence")
    if getattr(scorer, "CANDIDA_RESP_CONTEXT_NOT_PICKED_GUARDRAIL", "") in applied_rules:
        reasons.append("respiratory_candida_context_not_formal_picked")
    if "R-S4-RARE-YEAST-MNGS-ONLY-CONTEXT" in applied_rules:
        reasons.append("rare_opportunistic_yeast_mngs_only_context")
    if getattr(scorer, "CANDIDA_NONINVASIVE_HOSPITAL_ONLY_GUARDRAIL", "") in applied_rules:
        reasons.append("noninvasive_hospital_candida_not_formal_picked")
    if getattr(scorer, "NON_CORE_RESPIRATORY_VIRUS_FORMAL_PICK_GUARDRAIL", "") in applied_rules:
        reasons.append("non_core_respiratory_virus_context_not_formal_picked")
    if getattr(scorer, "MNGS_ONLY_NONDOMINANT_NON_TYPICAL_FORMAL_PICK_GUARDRAIL", "") in applied_rules:
        if reads_rank >= 3 or percentile >= 0.85:
            if (
                is_broad_or_atypical_oral_flora_name(candidate.get("organism_name"))
                and not non_host_support_modules(candidate)
                and dominance not in {"D2_moderate", "D3_dominant"}
            ):
                reasons.append("broad_or_atypical_oral_mngs_only_low_specificity")
            else:
                reasons.append("mngs_only_nondominant_non_typical_strong_context")
        else:
            reasons.append("mngs_only_nondominant_non_typical_low_specificity")
    if getattr(scorer, "COMMON_HOSPITAL_LOW_READ_NODOM_FORMAL_PICK_GUARDRAIL", "") in applied_rules:
        if scorer.related_representative_context_support(candidate):
            reasons.append("common_hospital_low_read_nondominant_related_context")
        else:
            reasons.append("common_hospital_low_read_nondominant_low_specificity")
    if level <= 2:
        reasons.append("full_ranked_level_1_2_not_in_formal_pick")
    if level == 3:
        reasons.append("borderline_level_3_not_in_formal_pick")
    if candidate.get("is_protected_pathogen"):
        reasons.append("protected_pathogen")
    if non_host_support_modules(candidate):
        reasons.append("same_organism_auxiliary_support")
    if rank <= 2 and reads_rank >= 2 and not common_or_background:
        reasons.append("high_rank_medium_or_high_reads_signal")
    if dominance in {"D2_moderate", "D3_dominant"} and reads_rank >= 1:
        reasons.append("dominant_signal")
    if candidate.get("specimen_class") == "S1_sterile_systemic" and rank <= 3:
        reasons.append("sterile_or_systemic_candidate")
    if ctx.get("host_level3_expansion") and rank <= 3 and reads_rank >= 1 and not common_or_background:
        reasons.append("host_vulnerable_ranked_signal")
    if rank <= 3 and percentile >= 0.75 and not common_or_background:
        reasons.append("high_percentile_uncommon_signal")
    if strong_lower_respiratory_aspiration_signal(candidate):
        reasons.append("lower_respiratory_aspiration_anaerobe_clean_high_signal_watchlist")
    if strong_water_environment_watch_signal(candidate):
        reasons.append("water_environmental_clean_high_signal_watchlist")

    return reasons


def candidate_summary(candidate: dict[str, Any], reasons: list[str] | None = None) -> dict[str, Any]:
    key_evidence = candidate.get("key_evidence") if isinstance(candidate.get("key_evidence"), dict) else {}
    classification = candidate.get("classification", "Unknown")
    taxonomy_profile = candidate.get("taxonomy_profile")
    if not isinstance(taxonomy_profile, dict):
        taxonomy_profile = organism_taxonomy.classify_organism(
            candidate.get("organism_name", ""),
            biological_class=classification,
        )
    output = {
        "organism_name": candidate.get("organism_name", ""),
        "classification": classification,
        "taxonomy_profile": taxonomy_profile,
        "integrated_causative_level": candidate.get("integrated_causative_level", "Unknown"),
        "mngs_signal_tier": candidate.get("mngs_signal_tier", "Unknown"),
        "rank_priority": candidate.get("rank_priority", "Unknown"),
        "rank_rule": candidate.get("rank_rule", ""),
        "reads": candidate.get("reads", 0),
        "reads_tier": candidate.get("reads_tier", "Unknown"),
        "reads_percentile": candidate.get("reads_percentile", 0),
        "dominance_tier": candidate.get("dominance_tier", "Unknown"),
        "specimen_class": candidate.get("specimen_class", "Unknown"),
        "specimen_alignment": candidate.get("specimen_alignment", "Unknown"),
        "source_category": candidate.get("source_category", "Unknown"),
        "is_protected_pathogen": bool(candidate.get("is_protected_pathogen")),
        "is_likely_colonizer_or_background": bool(candidate.get("is_likely_colonizer_or_background")),
        "support_modules": support_modules(candidate),
        "non_host_support_modules": non_host_support_modules(candidate),
        "guardrail_rule": key_evidence.get("guardrail_rule", ""),
        "formal_pick_exclusion_rule": candidate.get("formal_pick_exclusion_rule", {}),
        "related_representative_hospital_support": (
            candidate.get("related_representative_hospital_support")
            or key_evidence.get("related_representative_hospital_support")
            or []
        ),
        "candida_invasive_hospital_support": key_evidence.get("candida_invasive_hospital_support", []),
        "candida_related_invasive_hospital_context": key_evidence.get("candida_related_invasive_hospital_context", []),
        "candida_evidence_strength_profile": key_evidence.get("candida_evidence_strength_profile", {}),
        "candida_related_evidence_strength_profile": key_evidence.get("candida_related_evidence_strength_profile", {}),
        "applied_rules": list(candidate.get("applied_rules") or []),
    }
    if reasons is not None:
        output["review_reasons"] = reasons
    return output


def local_full_ranked_mngs_path(patient_dir: Path) -> Path | None:
    base = mngs.extract_patient_identifier(patient_dir)
    path = patient_dir / f"{base}_mNGS_ranked_candidates.json"
    return path if path.exists() else None


def supplemental_full_ranked_reasons(candidate: dict[str, Any], ctx: dict[str, Any] | None = None) -> list[str]:
    """Safety-only queue rules for candidates filtered out before max scoring.

    These rules intentionally do not make a candidate pickable. They only expose
    plausible but weak/rare signals to the LLM review layer. Keep thresholds
    narrower than deterministic scoring so common colonizers do not flood review.
    """
    reasons: list[str] = []
    name = candidate.get("organism_name")
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads = mngs.to_float(candidate.get("reads"))
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    specimen_cls = candidate.get("specimen_class")
    reads_tier = candidate.get("reads_tier")
    reads_rank = tier_rank(reads_tier)
    dominance = candidate.get("dominance_tier")
    host_vulnerable = is_host_vulnerable(ctx)
    has_aux_support = bool(non_host_support_modules(candidate))
    common_or_background = is_common_or_background(candidate)

    if common_or_background and not is_opportunistic_herpesvirus_name(name) and reads > 0:
        sterile_high_risk_signal = (
            specimen_cls == "S1_sterile_systemic"
            and rank <= 2
            and (reads_rank >= 1 or percentile >= 0.75 or has_aux_support)
        )
        respiratory_high_burden_signal = (
            specimen_cls == "S2_lower_respiratory"
            and rank <= 2
            and (reads_rank >= 3 or percentile >= 0.90)
        )
        respiratory_dominant_signal = (
            specimen_cls == "S2_lower_respiratory"
            and rank <= 3
            and reads_rank >= 1
            and dominance in {"D2_moderate", "D3_dominant"}
        )
        respiratory_auxiliary_signal = (
            specimen_cls == "S2_lower_respiratory"
            and rank <= 3
            and reads_rank >= 1
            and has_aux_support
        )
        if (
            sterile_high_risk_signal
            or respiratory_high_burden_signal
            or respiratory_dominant_signal
            or respiratory_auxiliary_signal
        ):
            reasons.append("full_ranked_common_background_high_burden_review")

    if is_rare_opportunistic_gnb_name(name) and rank <= 3 and reads > 0:
        if specimen_cls == "S2_lower_respiratory":
            reasons.append("full_ranked_rare_respiratory_opportunistic_gnb")
        elif specimen_cls == "S1_sterile_systemic" and percentile >= 0.75:
            reasons.append("full_ranked_rare_opportunistic_gnb_sterile_signal")

    if starts_with_any_name(name, ATYPICAL_ENTEROBACTERALES_PREFIXES) and rank <= 2:
        if reads_rank >= 2 or dominance in {"D2_moderate", "D3_dominant"} or (specimen_cls == "S1_sterile_systemic" and percentile >= 0.75):
            reasons.append("full_ranked_atypical_enterobacterales_strong_signal")

    if starts_with_any_name(name, ATYPICAL_PNEUMONIA_PREFIXES) and rank <= 4 and reads > 0:
        reasons.append("full_ranked_atypical_pneumonia_pathogen")

    if starts_with_any_name(name, HIGH_RISK_MOLD_OR_DIMORPHIC_PREFIXES) and rank <= 4 and reads > 0:
        if reads_rank >= 1 or host_vulnerable or has_aux_support:
            reasons.append("full_ranked_high_risk_mold_or_dimorphic_fungus")

    if starts_with_any_name(name, YEAST_OPPORTUNISTIC_PREFIXES) and rank <= 3 and reads > 0:
        if specimen_cls == "S1_sterile_systemic" and (reads_rank >= 1 or percentile >= 0.75):
            reasons.append("full_ranked_sterile_systemic_opportunistic_yeast")
        elif specimen_cls == "S2_lower_respiratory" and (has_aux_support or reads_rank >= 3 or (host_vulnerable and reads_rank >= 2)):
            reasons.append("full_ranked_respiratory_opportunistic_yeast_watchlist")

    if starts_with_any_name(name, STRICT_ASPIRATION_ANAEROBE_PREFIXES) and rank <= 3 and reads > 0:
        if specimen_cls == "S2_lower_respiratory" and clean_mngs_signal(candidate) and (reads_rank >= 3 or percentile >= 0.85):
            reasons.append("full_ranked_aspiration_anaerobe_or_oral_flora_clean_high_signal_watchlist")
        elif specimen_cls == "S2_lower_respiratory" and (
            reads_rank >= 3
            or (reads_rank >= 2 and (dominance in {"D2_moderate", "D3_dominant"} or percentile >= 0.75 or has_aux_support))
        ):
            reasons.append("full_ranked_aspiration_anaerobe_or_oral_flora_strong_signal")
        elif specimen_cls == "S1_sterile_systemic" and (reads_rank >= 2 or percentile >= 0.75):
            reasons.append("full_ranked_sterile_systemic_anaerobe_signal")

    if starts_with_any_name(name, SPECIAL_BACKGROUND_PATHOGEN_PREFIXES) and rank <= 2 and reads > 0:
        if specimen_cls == "S1_sterile_systemic" and (reads_rank >= 2 or percentile >= 0.75):
            reasons.append("full_ranked_special_background_pathogen_sterile_signal")
        elif specimen_cls == "S2_lower_respiratory" and (reads_rank >= 3 or has_aux_support):
            reasons.append("full_ranked_special_background_pathogen_respiratory_signal")

    if is_opportunistic_herpesvirus_name(name) and rank <= 2 and reads > 0:
        if (host_vulnerable and reads_rank >= 1) or reads_rank >= 3 or specimen_cls == "S1_sterile_systemic":
            reasons.append("full_ranked_opportunistic_herpesvirus_watchlist")

    if is_reactivation_herpesvirus_name(name) and rank <= 3 and reads > 0:
        if has_aux_support or dominance in {"D2_moderate", "D3_dominant"}:
            reasons.append("full_ranked_reactivation_herpesvirus_context_signal")
        elif reads_tier == "R4_very_high" and dominance in {"D0_not_top", "D1_low"}:
            reasons.append("full_ranked_reactivation_herpesvirus_r4_nondominant_context")
        elif host_vulnerable and rank <= 2 and reads_rank >= 2 and percentile >= 0.85 and dominance in {"D0_not_top", "D1_low"}:
            reasons.append("full_ranked_reactivation_herpesvirus_ranked_host_context")

    if (
        is_non_aureus_staphylococcus_name(name)
        and specimen_cls == "S1_sterile_systemic"
        and rank <= 2
        and (percentile >= 0.75 or reads >= 50)
    ):
        reasons.append("pulmonary_opt_removed_sterile_systemic_non_aureus_staphylococcus")

    if (
        scorer.is_coagulase_negative_staph_background_name(name)
        and specimen_cls == "S1_sterile_systemic"
        and rank <= 2
        and reads_tier in {"R4_very_high", "R3_high", "R2_medium"}
    ):
        reasons.append("pulmonary_opt_removed_sterile_systemic_skin_flora_watchlist")

    return reasons


def supplemental_full_ranked_candidates(
    patient_dir: Path,
    *,
    final_summary: dict[str, Any],
    existing_names: set[str],
    formal_picked_names: set[str],
    formal_picked_groups: set[str],
    host_context: dict[str, Any] | None = None,
) -> tuple[list[dict[str, Any]], Path | None, int]:
    full_ranked_path = local_full_ranked_mngs_path(patient_dir)
    if full_ranked_path is None:
        return [], None, 0
    ranked_payload = mngs.load_ranked_mngs_for_patient(full_ranked_path, patient_dir)
    if ranked_payload is None:
        return [], full_ranked_path, 0
    full_scored = scorer.score_payload(
        patient_dir=patient_dir,
        ranked_mngs=ranked_payload,
        final_summary=final_summary,
    )
    added: list[dict[str, Any]] = []
    for candidate in sorted(full_scored.get("pathogen_candidates") or [], key=scorer.candidate_sort_key):
        if not isinstance(candidate, dict):
            continue
        candidate_name = norm(candidate.get("organism_name"))
        if (
            not candidate_name
            or candidate_name in formal_picked_names
            or same_representative_group_already_picked(candidate.get("organism_name"), formal_picked_groups)
            or candidate_name in existing_names
        ):
            continue
        reasons = supplemental_full_ranked_reasons(candidate, host_context)
        if not reasons:
            continue
        item = candidate_summary(candidate, reasons)
        item["review_source"] = "local_full_ranked_safety_backfill"
        added.append(item)
        existing_names.add(candidate_name)
    return added, full_ranked_path, len(full_scored.get("pathogen_candidates") or [])


def first_nonempty(*values: Any) -> Any:
    for value in values:
        if value not in (None, "", [], {}):
            return value
    return None


def normalize_classification(value: Any) -> str:
    text = str(value or "Unknown").strip()
    mapping = {
        "bacteria": "Bacterial",
        "bacterial": "Bacterial",
        "fungi": "Fungal",
        "fungal": "Fungal",
        "virus": "Viral",
        "viral": "Viral",
        "parasite": "Parasitic",
        "parasitic": "Parasitic",
    }
    return mapping.get(text.lower(), text or "Unknown")


def representative_group_key(value: Any) -> str:
    return pathogen_names.representative_group_key(value)


def formal_picked_representative_groups(picked_items: Sequence[dict[str, Any]]) -> set[str]:
    groups: set[str] = set()
    for item in picked_items:
        if not isinstance(item, dict):
            continue
        group = representative_group_key(item.get("organism_name"))
        if group and not group.startswith("name:"):
            groups.add(group)
    return groups


def same_representative_group_already_picked(value: Any, formal_picked_groups: set[str]) -> bool:
    group = representative_group_key(value)
    return bool(
        group
        and not group.startswith("name:")
        and group in formal_picked_groups
        and pathogen_names.representative_group_suppression_allowed(value)
    )


def find_underlying_filmarray_summary(patient_dir: Path) -> Path | None:
    base = mngs.extract_patient_identifier(patient_dir)
    filename = f"{base}_final_summary_with_underlying_with_filmarray.json"
    for directory in (patient_dir / mngs.SUMMARY_OUTPUT_DIR_NAME, patient_dir):
        path = directory / filename
        if path.exists():
            return path
    return None


def add_hospital_candidate(
    candidates: dict[str, dict[str, Any]],
    *,
    organism_name: Any,
    classification: Any = None,
    level: Any = None,
    source_module: str,
    source_field: str,
    evidence: Any = None,
    pulmonary_causality_support: bool = True,
) -> None:
    name = str(organism_name or "").strip()
    if not name or norm(name) in {"no growth", "not detected", "none", "negative", "unknown"}:
        return
    key = norm(name)
    item = candidates.setdefault(
        key,
        {
            "organism_name": name,
            "classification": normalize_classification(classification),
            "taxonomy_profile": organism_taxonomy.classify_organism(
                name,
                biological_class=normalize_classification(classification),
            ),
            "integrated_causative_level": str(level or "Unknown"),
            "mngs_signal_tier": "Not_available",
            "rank_priority": "Not_available",
            "rank_rule": "hospital_detected_not_in_formal_pick",
            "reads": 0,
            "reads_tier": "Not_available",
            "reads_percentile": 0,
            "dominance_tier": "Not_available",
            "specimen_class": "Unknown",
            "specimen_alignment": "Unknown",
            "source_category": "hospital_detected_only",
            "is_protected_pathogen": False,
            "is_likely_colonizer_or_background": False,
            "support_modules": [],
            "non_host_support_modules": [],
            "guardrail_rule": "",
            "applied_rules": ["HOSPITAL-DETECTED-REVIEW"],
            "review_reasons": ["hospital_detected_not_in_formal_pick"],
            "hospital_review_source": [],
            "hospital_evidence": [],
        },
    )
    if classification and item.get("classification") in {"Unknown", ""}:
        item["classification"] = normalize_classification(classification)
    if level and item.get("integrated_causative_level") in {"Unknown", ""}:
        item["integrated_causative_level"] = str(level)
    if pulmonary_causality_support:
        if source_module not in item["support_modules"]:
            item["support_modules"].append(source_module)
        if source_module not in item["non_host_support_modules"]:
            item["non_host_support_modules"].append(source_module)
    else:
        context_modules = item.setdefault("context_only_modules", [])
        if source_module not in context_modules:
            context_modules.append(source_module)
        if scorer.NONPULMONARY_CULTURE_CONTEXT_RULE not in item["applied_rules"]:
            item["applied_rules"].append(scorer.NONPULMONARY_CULTURE_CONTEXT_RULE)
        if "nonpulmonary_hospital_evidence_context_only" not in item["review_reasons"]:
            item["review_reasons"].append("nonpulmonary_hospital_evidence_context_only")
    source_tag = f"{source_module}.{source_field}"
    if source_tag not in item["hospital_review_source"]:
        item["hospital_review_source"].append(source_tag)
    if evidence not in (None, "", [], {}) and evidence not in item["hospital_evidence"]:
        item["hospital_evidence"].append(evidence)


def hospital_detected_candidates(
    patient_dir: Path,
    *,
    final_summary: dict[str, Any],
) -> tuple[list[dict[str, Any]], Path | None]:
    candidates: dict[str, dict[str, Any]] = {}

    for evidence_item in final_summary.get("hospital_organism_evidence") or []:
        if not isinstance(evidence_item, dict):
            continue
        modules = evidence_item.get("evidence_modules") or []
        module = "+".join(str(x) for x in modules) if isinstance(modules, list) and modules else "hospital_organism_evidence"
        module_list = [str(x) for x in modules] if isinstance(modules, list) else []
        pulmonary_causality_support = any(
            source_module != "culture"
            or scorer.module_has_pulmonary_relevant_evidence(evidence_item, "culture")
            for source_module in module_list
        ) if module_list else True
        add_hospital_candidate(
            candidates,
            organism_name=evidence_item.get("organism_name") or evidence_item.get("name"),
            classification=evidence_item.get("classification"),
            level=evidence_item.get("best_hospital_level"),
            source_module=module,
            source_field="hospital_organism_evidence",
            evidence=evidence_item,
            pulmonary_causality_support=pulmonary_causality_support,
        )

    fallback_path = find_underlying_filmarray_summary(patient_dir)
    if fallback_path is not None:
        try:
            fallback = load_json(fallback_path)
        except Exception as exc:  # pragma: no cover - defensive CLI path
            print(f"warning {patient_dir.name}: failed to read fallback hospital summary {fallback_path}: {exc}")
            fallback = None
        if isinstance(fallback, dict):
            for item in fallback.get("pathogen_candidates") or []:
                if isinstance(item, dict):
                    add_hospital_candidate(
                        candidates,
                        organism_name=item.get("organism_name"),
                        classification=item.get("classification"),
                        level=item.get("integrated_causative_level"),
                        source_module="legacy_final_summary",
                        source_field="pathogen_candidates",
                        evidence={
                            "integrated_causative_level": item.get("integrated_causative_level"),
                            "module_level_summary": item.get("module_level_summary"),
                        },
                    )
            modules = fallback.get("module_summaries") or {}
            if isinstance(modules, dict):
                culture = modules.get("culture") or {}
                if isinstance(culture, dict):
                    for item in culture.get("isolates") or []:
                        if isinstance(item, dict):
                            add_hospital_candidate(
                                candidates,
                                organism_name=item.get("name") or item.get("organism_name"),
                                classification="Bacterial",
                                level=None,
                                source_module="legacy_culture",
                                source_field="module_summaries.culture.isolates",
                                evidence=item,
                            )
                filmarray = modules.get("filmarray") or modules.get("filmarray_gmtest") or {}
                if isinstance(filmarray, dict):
                    for field in ("detected_pathogens", "organism_labels"):
                        for item in filmarray.get(field) or []:
                            if isinstance(item, dict):
                                add_hospital_candidate(
                                    candidates,
                                    organism_name=first_nonempty(item.get("organism_name"), item.get("name")),
                                    classification=item.get("classification"),
                                    level=first_nonempty(item.get("causative_level"), item.get("level")),
                                    source_module="legacy_filmarray_gmtest",
                                    source_field=f"module_summaries.filmarray.{field}",
                                    evidence=item,
                                )
                            else:
                                add_hospital_candidate(
                                    candidates,
                                    organism_name=item,
                                    classification=None,
                                    level=None,
                                    source_module="legacy_filmarray_gmtest",
                                    source_field=f"module_summaries.filmarray.{field}",
                                    evidence=item,
                                )
    return list(candidates.values()), fallback_path


def build_queue_for_patient(
    patient_dir: Path,
    *,
    ranked_mngs_path: Path,
    deterministic_suffix: str,
    output_suffix: str,
    artifact_root: Path | None,
    overwrite: bool,
) -> Path | None:
    deterministic_path = output_path_for(patient_dir, deterministic_suffix)
    if not deterministic_path.exists():
        print(f"skip {patient_dir.name}: deterministic output missing {deterministic_path}")
        return None

    final_summary_path = mngs.find_final_summary_file(patient_dir, summary_mode="deterministic")
    if final_summary_path is None:
        print(f"skip {patient_dir.name}: deterministic final summary missing")
        return None

    ranked_payload = mngs.load_ranked_mngs_for_patient(ranked_mngs_path, patient_dir)
    if ranked_payload is None:
        print(f"skip {patient_dir.name}: ranked mNGS missing")
        return None

    destination = output_path_for(patient_dir, output_suffix, artifact_root=artifact_root)
    if destination.exists() and not overwrite:
        print(f"skip {patient_dir.name}: output exists {destination}")
        return None

    deterministic_payload = load_json(deterministic_path)
    final_summary = load_json(final_summary_path)
    full_scored = scorer.score_payload(
        patient_dir=patient_dir,
        ranked_mngs=ranked_payload,
        final_summary=final_summary,
    )

    formal_picked_items = [
        item
        for item in (deterministic_payload.get("best_available_summary") or {}).get("picked_pathogens") or []
        if isinstance(item, dict)
    ]
    formal_picked_names = {
        norm(item.get("organism_name"))
        for item in formal_picked_items
    }
    formal_picked_groups = formal_picked_representative_groups(formal_picked_items)

    review_queue: list[dict[str, Any]] = []
    skipped_picked: list[dict[str, Any]] = []
    skipped_same_representative_group: list[dict[str, Any]] = []
    excluded_low_risk: list[dict[str, Any]] = []
    excluded_impossible: list[dict[str, Any]] = []

    for candidate in sorted(full_scored.get("pathogen_candidates") or [], key=scorer.candidate_sort_key):
        if not isinstance(candidate, dict):
            continue
        candidate_name = norm(candidate.get("organism_name"))
        impossible_rule = pathogen_rule_lists.impossible_infection_source_rule(candidate.get("organism_name"))
        if impossible_rule:
            summary = candidate_summary(candidate, ["impossible_infection_source_excluded_from_llm_review"])
            summary["impossible_infection_source_rule"] = {
                "organism_name": impossible_rule.get("organism_name"),
                "action": impossible_rule.get("action"),
                "rationale_zh": impossible_rule.get("rationale_zh"),
                "notes_zh": impossible_rule.get("notes_zh"),
            }
            excluded_low_risk.append(summary)
            excluded_impossible.append(summary)
            continue
        if candidate_name in formal_picked_names:
            skipped_picked.append(candidate_summary(candidate))
            continue
        if same_representative_group_already_picked(candidate.get("organism_name"), formal_picked_groups):
            summary = candidate_summary(candidate, ["same_representative_group_already_picked"])
            summary["same_representative_group_key"] = representative_group_key(candidate.get("organism_name"))
            skipped_same_representative_group.append(summary)
            continue
        reasons = review_reasons(candidate, full_scored.get("host_context") or {})
        if reasons:
            review_queue.append(candidate_summary(candidate, reasons))
        else:
            excluded_low_risk.append(candidate_summary(candidate))

    hospital_candidates, legacy_hospital_summary_path = hospital_detected_candidates(
        patient_dir,
        final_summary=final_summary,
    )
    review_queue_by_name = {
        norm(item.get("organism_name")): item
        for item in review_queue
        if isinstance(item, dict)
    }
    supplemental_candidates, local_full_ranked_path, local_full_ranked_scored_count = supplemental_full_ranked_candidates(
        patient_dir,
        final_summary=final_summary,
        existing_names=set(review_queue_by_name),
        formal_picked_names=formal_picked_names,
        formal_picked_groups=formal_picked_groups,
        host_context=full_scored.get("host_context") or {},
    )
    for supplemental_candidate in supplemental_candidates:
        candidate_name = norm(supplemental_candidate.get("organism_name"))
        if not candidate_name:
            continue
        impossible_rule = pathogen_rule_lists.impossible_infection_source_rule(supplemental_candidate.get("organism_name"))
        if impossible_rule:
            supplemental_candidate = dict(supplemental_candidate)
            supplemental_candidate["review_reasons"] = ["impossible_infection_source_excluded_from_llm_review"]
            supplemental_candidate["impossible_infection_source_rule"] = {
                "organism_name": impossible_rule.get("organism_name"),
                "action": impossible_rule.get("action"),
                "rationale_zh": impossible_rule.get("rationale_zh"),
                "notes_zh": impossible_rule.get("notes_zh"),
            }
            excluded_impossible.append(supplemental_candidate)
            excluded_low_risk.append(supplemental_candidate)
            continue
        review_queue.append(supplemental_candidate)
        review_queue_by_name[candidate_name] = supplemental_candidate

    hospital_detected_added: list[dict[str, Any]] = []
    hospital_detected_merged: list[dict[str, Any]] = []
    hospital_detected_already_picked: list[dict[str, Any]] = []

    for hospital_candidate in sorted(hospital_candidates, key=lambda item: norm(item.get("organism_name"))):
        candidate_name = norm(hospital_candidate.get("organism_name"))
        if not candidate_name:
            continue
        impossible_rule = pathogen_rule_lists.impossible_infection_source_rule(hospital_candidate.get("organism_name"))
        if impossible_rule:
            hospital_candidate = dict(hospital_candidate)
            hospital_candidate["review_reasons"] = ["impossible_infection_source_excluded_from_llm_review"]
            hospital_candidate["impossible_infection_source_rule"] = {
                "organism_name": impossible_rule.get("organism_name"),
                "action": impossible_rule.get("action"),
                "rationale_zh": impossible_rule.get("rationale_zh"),
                "notes_zh": impossible_rule.get("notes_zh"),
            }
            excluded_impossible.append(hospital_candidate)
            excluded_low_risk.append(hospital_candidate)
            continue
        if candidate_name in formal_picked_names:
            hospital_detected_already_picked.append(hospital_candidate)
            continue
        if same_representative_group_already_picked(hospital_candidate.get("organism_name"), formal_picked_groups):
            hospital_candidate = dict(hospital_candidate)
            hospital_candidate["review_reasons"] = ["same_representative_group_already_picked"]
            hospital_candidate["same_representative_group_key"] = representative_group_key(
                hospital_candidate.get("organism_name")
            )
            skipped_same_representative_group.append(hospital_candidate)
            continue
        existing = review_queue_by_name.get(candidate_name)
        if existing is not None:
            reasons = existing.setdefault("review_reasons", [])
            for reason in hospital_candidate.get("review_reasons") or ["hospital_detected_not_in_formal_pick"]:
                if reason not in reasons:
                    reasons.append(reason)
            applied_rules = existing.setdefault("applied_rules", [])
            for rule in hospital_candidate.get("applied_rules") or []:
                if rule not in applied_rules:
                    applied_rules.append(rule)
            for key in ("hospital_review_source", "hospital_evidence"):
                values = existing.setdefault(key, [])
                for value in hospital_candidate.get(key) or []:
                    if value not in values:
                        values.append(value)
            for key in ("support_modules", "non_host_support_modules"):
                values = existing.setdefault(key, [])
                for value in hospital_candidate.get(key) or []:
                    if value not in values:
                        values.append(value)
            context_only_modules = existing.setdefault("context_only_modules", [])
            for value in hospital_candidate.get("context_only_modules") or []:
                if value not in context_only_modules:
                    context_only_modules.append(value)
            hospital_detected_merged.append(existing)
            continue
        review_queue.append(hospital_candidate)
        review_queue_by_name[candidate_name] = hospital_candidate
        hospital_detected_added.append(hospital_candidate)

    payload = {
        "review_queue_version": "missed_candidate_review_queue_v10_central_taxonomy_profile",
        "patient_id": mngs.extract_patient_identifier(patient_dir).replace("NGS_patient_", ""),
        "patient_dir": str(patient_dir),
        "source_files": {
            "formal_deterministic_output": str(deterministic_path),
            "deterministic_final_summary": str(final_summary_path),
            "legacy_underlying_filmarray_summary": str(legacy_hospital_summary_path) if legacy_hospital_summary_path else None,
            "ranked_mngs_for_deterministic_queue": str(ranked_mngs_path),
            "local_full_ranked_mngs_for_safety_backfill": str(local_full_ranked_path) if local_full_ranked_path else None,
        },
        "host_context": full_scored.get("host_context") or {},
        "formal_picked_count": len(formal_picked_names),
        "full_ranked_scored_count": len(full_scored.get("pathogen_candidates") or []),
        "local_full_ranked_scored_count": local_full_ranked_scored_count,
        "supplemental_full_ranked_added_count": len(supplemental_candidates),
        "hospital_detected_candidate_count": len(hospital_candidates),
        "hospital_detected_added_count": len(hospital_detected_added),
        "hospital_detected_merged_count": len(hospital_detected_merged),
        "hospital_detected_already_picked_count": len(hospital_detected_already_picked),
        "skipped_same_representative_group_count": len(skipped_same_representative_group),
        "review_queue_count": len(review_queue),
        "excluded_low_risk_count": len(excluded_low_risk),
        "excluded_impossible_infection_source_count": len(excluded_impossible),
        "excluded_impossible_infection_sources": excluded_impossible,
        "impossible_infection_source_policy": {
            "rule_file": str(pathogen_rule_lists.DEFAULT_IMPOSSIBLE_RULE_PATH),
            "organisms": pathogen_rule_lists.impossible_infection_source_names(),
            "policy_zh": "?? rules/impossible_infection_sources.json ????? deterministic max picked???? LLM missing queue?",
        },
        "review_queue_policy": [
            "This queue is a safety-review input only; it must not overwrite deterministic picked_pathogens.",
            "It reviews non-picked candidates from full ranked mNGS, especially level 1/2 misses, level 3 borderlines, protected pathogens, high-rank/high-read signals, dominance, sterile/systemic candidates, and host-vulnerable signals.",
            "Local full-ranked mNGS is also used as a safety backfill for rare opportunistic Gram-negative bacilli, atypical Enterobacterales, atypical pneumonia pathogens, high-risk fungi, selected opportunistic yeasts, aspiration anaerobes/oral flora with clean high signal, water/environmental organisms with clean high signal, special background pathogens, opportunistic herpesviruses, and high-burden common/background organisms that may have been removed by pulmonary opt.",
            "Hospital-detected organisms from deterministic hospital_organism_evidence are also queued when they are not already formal picked; legacy final_summary_with_underlying_with_filmarray is used as a fallback when deterministic evidence is unavailable.",
            "If a hospital or ranked candidate is already represented by a formal picked organism in the same representative group, it is kept as skipped/audit metadata rather than re-entering the LLM missing queue.",
            "The catch-all signal criteria are intentionally not limited to a fixed pathogen whitelist, reducing the chance of missing rare or newly recognized organisms.",
        ],
        "review_queue": review_queue,
        "supplemental_full_ranked_review_candidates": supplemental_candidates,
        "hospital_detected_review_candidates": hospital_detected_added,
        "hospital_detected_merged_candidates": hospital_detected_merged,
        "hospital_detected_already_picked": hospital_detected_already_picked,
        "skipped_already_picked": skipped_picked,
        "skipped_same_representative_group": skipped_same_representative_group,
        "excluded_low_risk": excluded_low_risk,
    }
    write_json(destination, payload)
    print(f"wrote {destination} review_queue_count={len(review_queue)}")
    return destination


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a deterministic-to-LLM missed-candidate review queue.")
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--patients", nargs="*")
    parser.add_argument("--ranked-mngs", type=Path, required=True)
    parser.add_argument("--deterministic-suffix", default=DEFAULT_DETERMINISTIC_SUFFIX)
    parser.add_argument("--output-suffix", default=DEFAULT_OUTPUT_SUFFIX)
    parser.add_argument(
        "--artifact-root",
        type=Path,
        help="Optional root for generated queue JSONs. Useful when the workspace drive is full.",
    )
    parser.add_argument("--summary-csv", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    requested = scorer.parse_patient_ids(args.patients)
    written: list[Path] = []

    patient_dirs = mngs.collect_patient_directories(
        [args.patient_root],
        summary_mode="deterministic",
        mngs_source="ranked_only",
    )
    if requested:
        patient_dirs = [
            path for path in patient_dirs if mngs.extract_patient_identifier(path) in requested
        ]

    for patient_dir in patient_dirs:
        out = build_queue_for_patient(
            patient_dir,
            ranked_mngs_path=args.ranked_mngs,
            deterministic_suffix=args.deterministic_suffix,
            output_suffix=args.output_suffix,
            artifact_root=args.artifact_root,
            overwrite=args.overwrite,
        )
        if out is not None:
            written.append(out)

    if args.summary_csv:
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.summary_csv.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "patient_id",
                    "formal_picked_count",
                    "full_ranked_scored_count",
                    "local_full_ranked_scored_count",
                    "supplemental_full_ranked_added_count",
                    "hospital_detected_candidate_count",
                    "hospital_detected_added_count",
                    "hospital_detected_merged_count",
                    "hospital_detected_already_picked_count",
                    "skipped_same_representative_group_count",
                    "review_queue_count",
                    "excluded_impossible_infection_source_count",
                    "queue_organisms",
                    "supplemental_full_ranked_added_organisms",
                    "hospital_detected_added_organisms",
                    "output_path",
                ],
            )
            writer.writeheader()
            for path in written:
                payload = load_json(path)
                writer.writerow(
                    {
                        "patient_id": payload.get("patient_id", ""),
                        "formal_picked_count": payload.get("formal_picked_count", 0),
                        "full_ranked_scored_count": payload.get("full_ranked_scored_count", 0),
                        "local_full_ranked_scored_count": payload.get("local_full_ranked_scored_count", 0),
                        "supplemental_full_ranked_added_count": payload.get("supplemental_full_ranked_added_count", 0),
                        "hospital_detected_candidate_count": payload.get("hospital_detected_candidate_count", 0),
                        "hospital_detected_added_count": payload.get("hospital_detected_added_count", 0),
                        "hospital_detected_merged_count": payload.get("hospital_detected_merged_count", 0),
                        "hospital_detected_already_picked_count": payload.get("hospital_detected_already_picked_count", 0),
                        "skipped_same_representative_group_count": payload.get("skipped_same_representative_group_count", 0),
                        "review_queue_count": payload.get("review_queue_count", 0),
                        "excluded_impossible_infection_source_count": payload.get("excluded_impossible_infection_source_count", 0),
                        "queue_organisms": "; ".join(
                            item.get("organism_name", "") for item in payload.get("review_queue") or []
                        ),
                        "supplemental_full_ranked_added_organisms": "; ".join(
                            item.get("organism_name", "")
                            for item in payload.get("supplemental_full_ranked_review_candidates") or []
                        ),
                        "hospital_detected_added_organisms": "; ".join(
                            item.get("organism_name", "") for item in payload.get("hospital_detected_review_candidates") or []
                        ),
                        "output_path": str(path),
                    }
                )
        print(f"wrote {args.summary_csv}")

    print(f"written_count={len(written)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

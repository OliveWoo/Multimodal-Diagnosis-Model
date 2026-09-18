from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import candidate_evidence_profile as evidence_profiles  # noqa: E402
from tools import mngs_common as mngs  # noqa: E402
from tools import organism_taxonomy_classifier as organism_taxonomy  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402
from tools import pathogen_rule_lists  # noqa: E402


DEFAULT_MODELLESS_OUTPUT_SUFFIX = "mNGS_max_deterministic_resp_commensal_dominance_guardrail_opt_chosen_full"
LEVEL_TEXT = {1: "Level 1", 2: "Level 2", 3: "Level 3", 4: "Level 4", 5: "Level 5"}
RANKED_RULE_VERSION = "mngs_deterministic_level_only_v20_central_taxonomy_profile"
DIRECT_HOSPITAL_MODULES = evidence_profiles.DIRECT_HOSPITAL_MODULES
LOW_SIGNAL_READ_TIERS = {"R0_trace", "R1_low"}
NON_DOMINANT_TIERS = {"D0_not_top", "D1_low"}
GENERIC_HOSPITAL_LABEL_GUARDRAIL = "R-S4-GENERIC-UNIDENTIFIED-HOSPITAL-LABEL"
L3_MNGS_ONLY_LOW_SIGNAL_GUARDRAIL = "R-S4-L3-MNGS-ONLY-HOST-LOW-NONDOMINANT"
SAME_REPRESENTATIVE_GROUP_GUARDRAIL = "R-S4-SAME-REPRESENTATIVE-GROUP-PICKED-OMIT"
RELATED_REPRESENTATIVE_SUPPORT_RULE = "R-S4-RELATED-REPRESENTATIVE-HOSPITAL-CONTEXT"
NONPULMONARY_CULTURE_CONTEXT_RULE = "R-S4-NONPULMONARY-CULTURE-CONTEXT-ONLY"
CANDIDA_RESP_CONTEXT_NOT_PICKED_GUARDRAIL = "R-S4-CANDIDA-RESP-CONTEXT-NOT-PICKED"
CANDIDA_NONINVASIVE_HOSPITAL_ONLY_GUARDRAIL = "R-S4-CANDIDA-HOSPITAL-NONINVASIVE-NOT-PICKED"
LOW_ACTIONABILITY_VIRUS_FORMAL_PICK_GUARDRAIL = "R-S4-LOW-ACTIONABILITY-VIRUS-PICKED-OMIT"
WATER_ENVIRONMENTAL_LOW_SPECIFICITY_FORMAL_PICK_GUARDRAIL = "R-S4-WATER-ENV-LOW-SPECIFICITY-PICKED-OMIT"
MNGS_R0_R1_HOST_ONLY_NONPROTECTED_FORMAL_PICK_GUARDRAIL = "R-S4-MNGS-R0R1-HOST-ONLY-NONPROTECTED-PICKED-OMIT"
WEAK_MOLD_HOST_ONLY_FORMAL_PICK_GUARDRAIL = "R-S4-WEAK-MOLD-HOST-ONLY-PICKED-OMIT"
NON_CORE_RESPIRATORY_VIRUS_FORMAL_PICK_GUARDRAIL = "R-S4-NON-CORE-RESP-VIRUS-PICKED-OMIT"
MNGS_ONLY_NONDOMINANT_NON_TYPICAL_FORMAL_PICK_GUARDRAIL = "R-S4-MNGS-ONLY-NONDOMINANT-NON-TYPICAL-PICKED-OMIT"
COMMON_HOSPITAL_LOW_READ_NODOM_FORMAL_PICK_GUARDRAIL = "R-S4-COMMON-HOSPITAL-R0R1-NONDOMINANT-PICKED-OMIT"
HERPES_HOST_R2R3_NONDOMINANT_REVIEW_RULE = "R-S4-HERPES-HOST-R2R3-NONDOMINANT-REVIEW"
ORAL_ASPIRATION_FLORA_L3_GUARDRAIL = "R-S4-ORAL-ASPIRATION-FLORA-L3"
ORAL_ASPIRATION_FLORA_L4_GUARDRAIL = "R-S4-ORAL-ASPIRATION-FLORA-MNGS-ONLY-NONDOMINANT-L4"
RESPIRATORY_TERMS = (
    "respiratory",
    "lung",
    "pneumonia",
    "bal",
    "balf",
    "bronchoalveolar",
    "eta",
    "endotracheal",
    "sputum",
    "tracheal",
    "呼吸道",
    "肺",
    "痰",
    "氣管",
)
STERILE_TERMS = (
    "blood",
    "csf",
    "pleural",
    "peritoneal",
    "ascites",
    "synovial",
    "tissue",
    "biopsy",
    "sterile",
    "血",
    "腦脊髓",
    "胸水",
    "腹水",
    "關節液",
    "組織",
)
LOWER_RESPIRATORY_TERMS = (
    "balf",
    "bal",
    "bronchoalveolar",
    "sputum",
    "eta",
    "endotracheal",
    "tracheal",
)
INVASIVE_CANDIDA_SPECIMEN_CATEGORIES = evidence_profiles.INVASIVE_CANDIDA_SPECIMEN_CATEGORIES
INVASIVE_CANDIDA_SPECIMEN_TERMS = evidence_profiles.INVASIVE_CANDIDA_SPECIMEN_TERMS
CANDIDA_EVIDENCE_STRONG_INVASIVE = evidence_profiles.CANDIDA_EVIDENCE_STRONG_INVASIVE
CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC = evidence_profiles.CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC
CANDIDA_EVIDENCE_WEAK_NONINVASIVE = evidence_profiles.CANDIDA_EVIDENCE_WEAK_NONINVASIVE
CANDIDA_INTERMEDIATE_EVIDENCE_TERMS = evidence_profiles.CANDIDA_INTERMEDIATE_EVIDENCE_TERMS
COAG_NEG_STAPH_BACKGROUND_PREFIXES = (
    "staphylococcuscapitis",
    "staphylococcusepidermidis",
    "staphylococcushaemolyticus",
    "staphylococcushominis",
    "staphylococcuswarneri",
    "staphylococcussimulans",
    "staphylococcuscohnii",
    "staphylococcuscaprae",
)
REACTIVATION_VIRUS_NAMES = {
    "humanbetaherpesvirus5",  # CMV
    "humanalphaherpesvirus1",  # HSV-1
    "humanalphaherpesvirus2",  # HSV-2
    "humanalphaherpesvirus3",  # VZV
}
LOW_SPECIFICITY_VIRUS_NAMES = {
    "torquetenovirus",
    "torquetenominivirus",
    "torquetenomidivirus",
}
LOW_ACTIONABILITY_VIRUS_PREFIXES = (
    "humanpapillomavirus",
    "papillomavirus",
    "ttv",
    "torqueteno",
    "anellovirus",
)
SKIN_ENVIRONMENTAL_YEAST_PREFIXES = (
    "malassezia",
)
RARE_OPPORTUNISTIC_YEAST_PREFIXES = (
    "trichosporon",
)
LOW_PULMONARY_SPECIFICITY_ATYPICAL_PREFIXES = (
    "metamycoplasmasalivarium",
    "metamycoplasmahominis",
    "mycoplasmasalivarium",
    "mycoplasmaorale",
)
SKIN_LOW_SPECIFICITY_BACTERIA_PREFIXES = (
    "cutibacterium",
    "propionibacterium",
)
SARS_COV_2_NAMES = {
    "severeacuterespiratorysyndromerelatedcoronavirus",
    "severeacuterespiratorysyndromecoronavirus2",
}
TYPICAL_PNEUMONIA_BACTERIA_PREFIXES = (
    "streptococcuspneumoniae",
    "haemophilusinfluenzae",
    "moraxellacatarrhalis",
    "staphylococcusaureus",
    "klebsiella",
    "escherichiacoli",
    "enterobacter",
    "citrobacter",
    "serratia",
    "proteus",
    "morganella",
    "pseudomonasaeruginosa",
    "acinetobacter",
    "stenotrophomonasmaltophilia",
    "burkholderia",
    "achromobacterxylosoxidans",
    "legionella",
)
ASPERGILLUS_MOLD_PREFIXES = (
    "aspergillus",
)
RESPIRATORY_COMMENSAL_FLORA_EXACT_NAMES = {
    "abiotrophiadefectiva",
    "neisseriasubflava",
    "neisseriaflavescens",
    "neisseriamucosa",
    "neisseriasicca",
    "neisseriacinerea",
    "neisseriaperflava",
    "neisseriaelongata",
}
RESPIRATORY_COMMENSAL_FLORA_PREFIXES = (
    "rothia",
    "gemella",
    "granulicatella",
    "prevotella",
    "veillonella",
    "lactiplantibacillus",
    "lacticaseibacillus",
    "limosilactobacillus",
    "ligilactobacillus",
    "lactobacillus",
    "alloscardovia",
    "parvimonas",
    "porphyromonas",
    "fusobacterium",
    "leptotrichia",
    "capnocytophaga",
    "actinomyces",
    "oribacterium",
    "slackia",
    "treponema",
    "bacteroides",
    "phocaeicola",
    "segatella",
    "blautia",
    "faecalibacterium",
    "parabacteroides",
    "alistipes",
    "odoribacter",
    "eubacterium",
)
ORAL_ASPIRATION_FLORA_EXACT_NAMES = {
    "campylobacterconcisus",
    "kingellaoralis",
    "moraxellaosloensis",
    "mycoplasmasalivarium",
    "neisseriacinerea",
    "neisseriaelongata",
    "neisseriaflavescens",
    "neisseriamucosa",
    "neisseriaperflava",
    "neisseriasicca",
    "neisseriasubflava",
    "rothiamucilaginosa",
    "streptococcusparasanguinis",
    "streptococcussalivarius",
    "streptococcustaonis",
    "viridansstreptococci",
}
ORAL_ASPIRATION_FLORA_PREFIXES = (
    "abiotrophia",
    "actinomyces",
    "bacteroides",
    "capnocytophaga",
    "finegoldia",
    "fusobacterium",
    "gemella",
    "lactobacillus",
    "lacticaseibacillus",
    "lactiplantibacillus",
    "leptotrichia",
    "ligilactobacillus",
    "limosilactobacillus",
    "parabacteroides",
    "parvimonas",
    "peptostreptococcus",
    "phocaeicola",
    "porphyromonas",
    "prevotella",
    "rothia",
    "segatella",
    "veillonella",
)
LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_PREFIXES = (
    "arthrobacter",
    "agrobacterium",
    "brachybacterium",
    "janibacter",
    "kocuria",
    "microbacterium",
    "micrococcus",
    "methylobacterium",
    "novosphingobium",
    "paracoccus",
    "polynucleobacter",
    "sphingobium",
    "sphingomonas",
    "undibacterium",
)
WATER_ENVIRONMENTAL_GNB_PREFIXES = (
    "achromobacter",
    "acidovorax",
    "brevundimonas",
    "cloacibacterium",
    "chryseobacterium",
    "comamonas",
    "delftia",
    "elizabethkingia",
    "flavobacterium",
    "paraburkholderia",
    "phytobacter",
    "ralstonia",
    "pseudoxanthomonas",
    "xanthomonas",
)
WATER_ENVIRONMENTAL_GNB_EXACT_NAMES = {
    "acinetobacterjohnsonii",
    "brucellaintermedia",
    "enterobactersoli",
    "pseudomonasalcaligenes",
    "pseudomonasceruminis",
    "pseudomonasjaponica",
    "pseudomonaskhazarica",
    "pseudomonaslactis",
    "pseudomonassoli",
    "pseudomonaszhaodongensis",
}
PICKED_SAFE_WATER_ENVIRONMENTAL_GNB_PREFIXES = (
    "achromobacter",
    "acidovorax",
    "brevundimonas",
    "cloacibacterium",
    "comamonas",
    "delftia",
    "elizabethkingia",
    "flavobacterium",
    "paraburkholderia",
    "phytobacter",
    "ralstonia",
    "pseudoxanthomonas",
    "xanthomonas",
)
PICKED_SAFE_LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_PREFIXES = (
    "arthrobacter",
    "agrobacterium",
    "brachybacterium",
    "janibacter",
    "kocuria",
    "microbacterium",
    "micrococcus",
    "methylobacterium",
    "novosphingobium",
    "paracoccus",
    "polynucleobacter",
    "sphingobium",
    "sphingomonas",
    "undibacterium",
)
PICKED_SAFE_WATER_ENVIRONMENTAL_GNB_EXACT_NAMES = {
    "acinetobacterjohnsonii",
    "brucellaintermedia",
    "enterobactersoli",
    "pseudomonasalcaligenes",
    "pseudomonasceruminis",
    "pseudomonasjaponica",
    "pseudomonaskhazarica",
    "pseudomonaslactis",
    "pseudomonassoli",
    "pseudomonaszhaodongensis",
}
PICKED_SAFE_MOLD_OPPORTUNISTIC_PREFIXES = (
    "aspergillus",
    "cunninghamella",
    "talaromyces",
    "schizophyllum",
    "rhizopus",
    "mucor",
    "lichtheimia",
    "fusarium",
)
PICKED_SAFE_HIGH_CONSEQUENCE_NAMES = {
    "pneumocystisjirovecii",
    "mycobacteriumtuberculosis",
    "legionellapneumophila",
    "nocardiathailandica",
    "burkholderiacenocepacia",
    "severeacuterespiratorysyndromerelatedcoronavirus",
}
PICKED_SAFE_COMMON_RESPIRATORY_NAMES = {
    "acinetobacterbaumannii",
    "klebsiellapneumoniae",
    "klebsiellavariicola",
    "pseudomonasaeruginosa",
    "serratiamarcescens",
    "stenotrophomonasmaltophilia",
    "escherichiacoli",
    "haemophilusinfluenzae",
    "staphylococcusaureus",
}
RESPIRATORY_VIRUS_PREFIXES = (
    "influenza",
    "parainfluenza",
    "humanrespirovirus",
    "humanrhinovirus",
    "rhinovirus",
    "respiratorysyncytialvirus",
    "severeacuterespiratorysyndrome",
    "sarscov2",
)


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def lower_text(value: Any) -> str:
    return str(value or "").strip().lower()


def rank_priority(candidate: dict[str, Any]) -> int:
    ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
    return mngs.to_rank_int(ranking.get("rank_priority"))


def reads_percentile(candidate: dict[str, Any]) -> float:
    ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
    return mngs.to_float(ranking.get("reads_percentile"))


def possibility(candidate: dict[str, Any]) -> str:
    ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
    return lower_text(ranking.get("possibility_level"))


def is_reactivation_virus_candidate(candidate: dict[str, Any]) -> bool:
    return mngs.normalize_organism_name(candidate.get("organism_name")) in REACTIVATION_VIRUS_NAMES



def is_low_specificity_virus_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return normalized in LOW_SPECIFICITY_VIRUS_NAMES or any(
        normalized.startswith(prefix) for prefix in LOW_ACTIONABILITY_VIRUS_PREFIXES
    )


def is_skin_environmental_yeast_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in SKIN_ENVIRONMENTAL_YEAST_PREFIXES)


def is_rare_opportunistic_yeast_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in RARE_OPPORTUNISTIC_YEAST_PREFIXES)


def is_low_pulmonary_specificity_atypical_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in LOW_PULMONARY_SPECIFICITY_ATYPICAL_PREFIXES)


def is_skin_low_specificity_bacteria_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in SKIN_LOW_SPECIFICITY_BACTERIA_PREFIXES)


def is_aspergillus_mold_candidate(candidate: dict[str, Any]) -> bool:
    normalized = mngs.normalize_organism_name(candidate.get("organism_name"))
    return any(normalized.startswith(prefix) for prefix in ASPERGILLUS_MOLD_PREFIXES)


def is_respiratory_commensal_flora_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return normalized in RESPIRATORY_COMMENSAL_FLORA_EXACT_NAMES or is_oral_aspiration_flora_name(value) or any(
        normalized.startswith(prefix) for prefix in RESPIRATORY_COMMENSAL_FLORA_PREFIXES
    )


def is_oral_aspiration_flora_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return normalized in ORAL_ASPIRATION_FLORA_EXACT_NAMES or mngs.is_oral_upper_airway_flora_name(value) or any(
        normalized.startswith(prefix) for prefix in ORAL_ASPIRATION_FLORA_PREFIXES
    )


def is_low_pulmonary_specificity_environmental_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_PREFIXES)


def is_water_environmental_gnb_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return normalized in WATER_ENVIRONMENTAL_GNB_EXACT_NAMES or any(
        normalized.startswith(prefix) for prefix in WATER_ENVIRONMENTAL_GNB_PREFIXES
    )


def is_picked_safe_water_environmental_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return (
        normalized in PICKED_SAFE_WATER_ENVIRONMENTAL_GNB_EXACT_NAMES
        or any(normalized.startswith(prefix) for prefix in PICKED_SAFE_WATER_ENVIRONMENTAL_GNB_PREFIXES)
        or any(
            normalized.startswith(prefix)
            for prefix in PICKED_SAFE_LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_PREFIXES
        )
    )


def is_picked_safe_mold_opportunistic_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in PICKED_SAFE_MOLD_OPPORTUNISTIC_PREFIXES)


def is_picked_safe_respiratory_virus_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in RESPIRATORY_VIRUS_PREFIXES) or normalized in REACTIVATION_VIRUS_NAMES


def is_core_answer_sensitive_respiratory_virus_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return (
        normalized in REACTIVATION_VIRUS_NAMES
        or normalized in SARS_COV_2_NAMES
        or normalized.startswith("influenza")
        or normalized.startswith("severeacuterespiratorysyndrome")
        or normalized.startswith("sarscov2")
    )


def is_non_core_respiratory_virus_name(value: Any) -> bool:
    return (
        is_picked_safe_respiratory_virus_name(value)
        and not is_core_answer_sensitive_respiratory_virus_name(value)
    )


def is_picked_safe_answer_sensitive_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return (
        normalized in PICKED_SAFE_HIGH_CONSEQUENCE_NAMES
        or normalized in PICKED_SAFE_COMMON_RESPIRATORY_NAMES
        or is_picked_safe_respiratory_virus_name(value)
    )


def pathogen_evidence_key(value: Any) -> str:
    """Clinical matching key for same-organism evidence across module labels."""
    key, _, _ = pathogen_names.alias_info(value)
    if not key:
        return ""
    normalized = mngs.normalize_organism_name(value)
    if key == "aspergillusspp" or normalized.startswith("aspergillus"):
        return "aspergillusspp"
    return key


def final_candidate_for_name(final_by_name: dict[str, dict[str, Any]], name: Any) -> dict[str, Any]:
    normalized = mngs.normalize_organism_name(name)
    candidate = final_by_name.get(normalized)
    if isinstance(candidate, dict):
        return candidate
    target = pathogen_evidence_key(name)
    if not target:
        return {}
    for key, value in final_by_name.items():
        if pathogen_evidence_key(key) == target and isinstance(value, dict):
            return value
        if isinstance(value, dict) and pathogen_evidence_key(value.get("organism_name") or value.get("name")) == target:
            return value
    return {}


def is_sars_cov2_candidate(candidate: dict[str, Any]) -> bool:
    return mngs.normalize_organism_name(candidate.get("organism_name")) in SARS_COV_2_NAMES


def is_typical_pneumonia_bacterium_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in TYPICAL_PNEUMONIA_BACTERIA_PREFIXES)


def clean_mngs_signal(candidate: dict[str, Any]) -> bool:
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads = mngs.to_float(candidate.get("reads"))
    rule = lower_text(candidate.get("rank_rule"))
    if reads <= 0 or rank > 3:
        return False
    return rank <= 1 or "ntc=00" in rule or "rk00" in rule or "ntc" in rule


def strong_sars_cov2_mngs_signal(candidate: dict[str, Any]) -> bool:
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    tier = candidate.get("reads_tier")
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    return clean_mngs_signal(candidate) and rank <= 3 and (tier in {"R4_very_high", "R3_high"} or percentile >= 0.85)


def typical_pneumonia_hospital_mngs_rule(candidate: dict[str, Any]) -> bool:
    if candidate.get("classification") != "Bacterial":
        return False
    if candidate.get("specimen_class") != "S2_lower_respiratory":
        return False
    if not is_typical_pneumonia_bacterium_name(candidate.get("organism_name")):
        return False
    if not clean_mngs_signal(candidate):
        return False
    support = candidate.get("module_support_summary") or {}
    return support.get("culture") == "Support" or support.get("filmarray_gmtest") == "Support"


def rare_opportunistic_yeast_supported_mngs_hospital_rule(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    if not is_rare_opportunistic_yeast_name(candidate.get("organism_name")):
        return False
    if candidate.get("specimen_class") != "S2_lower_respiratory":
        return False
    if not same_organism_auxiliary_support(candidate):
        return False
    if not host_support(ctx):
        return False
    reads_tier = candidate.get("reads_tier")
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    dominance = candidate.get("dominance_tier")
    return (
        reads_tier in {"R4_very_high", "R3_high"}
        or percentile >= 0.90
        or dominance in {"D2_moderate", "D3_dominant"}
    )


def same_organism_auxiliary_support(candidate: dict[str, Any]) -> bool:
    support = candidate.get("module_support_summary") or {}
    return any(
        support.get(key) == "Support"
        for key in ("culture", "filmarray_gmtest", "molecular_microbiology")
    )


def herpes_has_non_host_support(candidate: dict[str, Any]) -> bool:
    """Direct organism-specific testing support outside host vulnerability."""
    return same_organism_auxiliary_support(candidate)


def herpes_has_dominant_signal(candidate: dict[str, Any]) -> bool:
    return candidate.get("dominance_tier") in {"D2_moderate", "D3_dominant"}


def herpes_high_reads_nondominant_host_only(candidate: dict[str, Any]) -> bool:
    """High-burden herpes signal that should be reviewed, not auto-picked.

    R4 D0/D1 CMV/HSV/VZV can be clinically important, especially in vulnerable
    hosts, but without dominance or direct PCR/viral-load/hospital support it is
    safer as RAG-visible context than as formal picked pneumonia attribution.
    """
    return (
        is_reactivation_virus_candidate(candidate)
        and candidate.get("reads_tier") == "R4_very_high"
        and candidate.get("dominance_tier") in NON_DOMINANT_TIERS
        and not herpes_has_non_host_support(candidate)
    )


def herpes_ranked_host_context_review(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    """R2/R3 HSV/CMV/VZV host-only signal that should be RAG-visible, not picked."""
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    return (
        is_reactivation_virus_candidate(candidate)
        and ctx.get("host_level3_expansion")
        and rank <= 2
        and candidate.get("reads_tier") in {"R2_medium", "R3_high"}
        and percentile >= 0.85
        and candidate.get("dominance_tier") in NON_DOMINANT_TIERS
        and not herpes_has_non_host_support(candidate)
    )


def is_generic_unidentified_hospital_label(value: Any) -> bool:
    text = lower_text(value)
    if not text:
        return False
    compact = pathogen_names.raw_key(text)
    if "identificationtofollow" in compact:
        return True
    if compact in {
        "gbacilli",
        "gbacillus",
        "gnegativebacilli",
        "gnegativebacillus",
        "gpositivebacilli",
        "gpositivebacillus",
        "gramnegativebacilli",
        "gramnegativebacillus",
        "grampositivebacilli",
        "grampositivebacillus",
        "yeast",
        "yeasts",
    }:
        return True
    if ("gramnegative" in compact or "grampositive" in compact) and (
        "bacilli" in compact or "bacillus" in compact or "cocci" in compact or "coccus" in compact
    ):
        return True
    return False


def mark_formal_pick_exclusion(candidate: dict[str, Any], *, rule_id: str, reason: str) -> None:
    existing = candidate.get("formal_pick_exclusion_rule")
    if isinstance(existing, dict) and existing.get("rule_id") == rule_id:
        return
    candidate["formal_pick_exclusion_rule"] = {
        "rule_id": rule_id,
        "action": "exclude_from_formal_picked_pathogens",
        "reason": reason,
    }
    key_evidence = candidate.get("key_evidence")
    if not isinstance(key_evidence, dict):
        key_evidence = {}
        candidate["key_evidence"] = key_evidence
    if not key_evidence.get("guardrail_rule"):
        key_evidence["guardrail_rule"] = rule_id
    extra_rules = key_evidence.setdefault("formal_pick_exclusion_rules", [])
    if isinstance(extra_rules, list) and rule_id not in extra_rules:
        extra_rules.append(rule_id)
    applied_rules = candidate.setdefault("applied_rules", [])
    if isinstance(applied_rules, list) and rule_id not in applied_rules:
        applied_rules.append(rule_id)
    reasoning = candidate.setdefault("integrated_reasoning", [])
    if isinstance(reasoning, list):
        reasoning.append(reason)


def formal_pick_excluded(candidate: dict[str, Any]) -> bool:
    return isinstance(candidate.get("formal_pick_exclusion_rule"), dict)


def classify_source_category(value: Any) -> str:
    text = lower_text(value)
    if "virus" in text or "viral" in text:
        return "Viral"
    if "fung" in text:
        return "Fungal"
    if "bac" in text or "bacteria" in text:
        return "Bacterial"
    if "parasite" in text:
        return "Parasitic"
    return "Unknown"


def specimen_alignment(record: dict[str, Any]) -> str:
    haystack = " ".join(
        lower_text(record.get(key))
        for key in ("specimen_site", "specimen_name", "specimen_code", "seq_id")
    )
    if any(term in haystack for term in STERILE_TERMS):
        return "Sterile_or_Systemic"
    if any(term in haystack for term in RESPIRATORY_TERMS):
        return "Aligned"
    return "Unknown"


def specimen_class(record: dict[str, Any], alignment: str) -> str:
    haystack = " ".join(
        lower_text(record.get(key))
        for key in ("specimen_site", "specimen_name", "specimen_code", "seq_id")
    )
    if alignment == "Sterile_or_Systemic" or any(term in haystack for term in STERILE_TERMS):
        return "S1_sterile_systemic"
    if any(term in haystack for term in LOWER_RESPIRATORY_TERMS):
        return "S2_lower_respiratory"
    return "S3_unknown_or_low_value"


def is_sterile_or_systemic(record: dict[str, Any], alignment: str) -> bool:
    return alignment == "Sterile_or_Systemic"


def is_respiratory(record: dict[str, Any], alignment: str) -> bool:
    if alignment == "Sterile_or_Systemic":
        return False
    haystack = " ".join(
        lower_text(record.get(key))
        for key in ("specimen_site", "specimen_name", "specimen_code", "seq_id")
    )
    return alignment == "Aligned" or any(term in haystack for term in RESPIRATORY_TERMS)


def is_coagulase_negative_staph_background_name(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in COAG_NEG_STAPH_BACKGROUND_PREFIXES)


def host_context(final_summary: dict[str, Any]) -> dict[str, Any]:
    ctx = final_summary.get("host_context")
    if not isinstance(ctx, dict):
        ctx = {}
    module_host = ((final_summary.get("module_summaries") or {}).get("cbc_other_lab") or {})
    if not isinstance(module_host, dict):
        module_host = {}
    tier = ctx.get("host_vulnerability_tier") or module_host.get("host_vulnerability_tier") or "Unknown"
    opp = ctx.get("opportunistic_coverage_level") or module_host.get("opportunistic_coverage_level") or "Unknown"
    expanded = bool(
        ctx.get("expanded_candidate_policy")
        or module_host.get("expanded_candidate_policy")
        or tier in {"V2", "V3"}
        or opp in {"O1", "O2"}
    )
    return {
        "host_vulnerability_tier": tier,
        "opportunistic_coverage_level": opp,
        "host_level3_expansion": expanded,
    }


def final_candidate_map(final_summary: dict[str, Any]) -> dict[str, dict[str, Any]]:
    mapped: dict[str, dict[str, Any]] = {}
    for candidate in final_summary.get("pathogen_candidates") or []:
        if not isinstance(candidate, dict):
            continue
        key = mngs.normalize_organism_name(candidate.get("organism_name") or candidate.get("name"))
        if key:
            mapped[key] = candidate
    for evidence_item in final_summary.get("hospital_organism_evidence") or []:
        if not isinstance(evidence_item, dict):
            continue
        key = mngs.normalize_organism_name(evidence_item.get("organism_name") or evidence_item.get("name"))
        if not key:
            continue
        current = mapped.setdefault(
            key,
            {
                "organism_name": evidence_item.get("organism_name", ""),
                "classification": evidence_item.get("classification", "Unknown"),
                "evidence_modules": [],
                "module_level_summary": {
                    "culture": "Not_available",
                    "filmarray_gmtest": "Not_available",
                    "image": "Not_available",
                    "molecular_microbiology": "Not_available",
                    "cbc_other_lab": "Not_pathogen_specific",
                },
                "module_evidence": {},
            },
        )
        modules = evidence_item.get("evidence_modules")
        if isinstance(modules, list):
            existing = current.setdefault("evidence_modules", [])
            if isinstance(existing, list):
                for module in modules:
                    if module not in existing:
                        existing.append(module)
        levels = evidence_item.get("module_level_summary")
        if isinstance(levels, dict):
            current_levels = current.setdefault("module_level_summary", {})
            if isinstance(current_levels, dict):
                current_levels.update(levels)
        module_evidence = evidence_item.get("module_evidence")
        if isinstance(module_evidence, dict):
            current_evidence = current.setdefault("module_evidence", {})
            if isinstance(current_evidence, dict):
                current_evidence.update(module_evidence)
        if evidence_item.get("best_hospital_level"):
            current.setdefault("best_hospital_level", evidence_item.get("best_hospital_level"))
    return mapped


def final_candidate_level(final_by_name: dict[str, dict[str, Any]], name: Any) -> int:
    candidate = final_candidate_for_name(final_by_name, name)
    if not isinstance(candidate, dict) or not candidate:
        return 99
    return min(
        mngs._level_rank(candidate.get("integrated_causative_level")),  # pylint: disable=protected-access
        mngs._level_rank(candidate.get("best_hospital_level")),  # pylint: disable=protected-access
    )


def direct_hospital_support_modules(evidence_item: dict[str, Any], *, max_level_rank: int = 2) -> list[str]:
    """Direct organism tests that can stand alone when mNGS does not contain the organism."""
    levels = evidence_item.get("module_level_summary")
    if not isinstance(levels, dict):
        return []
    output: list[str] = []
    for module in DIRECT_HOSPITAL_MODULES:
        if level_rank(levels.get(module)) > max_level_rank:
            continue
        if module == "culture" and not module_has_pulmonary_relevant_evidence(evidence_item, module):
            continue
        if level_rank(levels.get(module)) <= max_level_rank:
            output.append(module)
    return output


def hospital_evidence_detail(evidence_item: dict[str, Any], modules: Sequence[str]) -> dict[str, Any]:
    module_evidence = evidence_item.get("module_evidence")
    if not isinstance(module_evidence, dict):
        module_evidence = {}
    details: dict[str, Any] = {
        "best_hospital_level": evidence_item.get("best_hospital_level", "Unknown"),
        "support_modules": list(modules),
        "module_level_summary": evidence_item.get("module_level_summary", {}),
    }
    for module in modules:
        values = module_evidence.get(module)
        if isinstance(values, list):
            details[module] = values
        elif values is not None:
            details[module] = [values]
    return details


def is_candida_or_generic_yeast_name(value: Any) -> bool:
    return pathogen_names.is_candida_or_generic_yeast(value)


def candida_evidence_row_text(row: dict[str, Any]) -> str:
    return evidence_profiles.candida_evidence_row_text(row)


def candida_evidence_row_strength(row: Any) -> tuple[str, str]:
    return evidence_profiles.candida_evidence_row_strength(row)


def candida_evidence_row_is_invasive(row: Any) -> bool:
    return evidence_profiles.candida_evidence_row_is_invasive(row)


def candida_support_details(
    evidence_item: dict[str, Any],
    *,
    strengths: set[str] | None = None,
) -> list[dict[str, Any]]:
    return evidence_profiles.candida_support_details(evidence_item, strengths=strengths)


def candida_invasive_support_details(evidence_item: dict[str, Any]) -> list[dict[str, Any]]:
    return candida_support_details(evidence_item, strengths={CANDIDA_EVIDENCE_STRONG_INVASIVE})


def candida_evidence_strength_profile(details: Sequence[dict[str, Any]]) -> dict[str, Any]:
    return evidence_profiles.candida_evidence_strength_profile(details)


def candida_hospital_evidence_strength_profile(evidence_item: dict[str, Any]) -> dict[str, Any]:
    return candida_evidence_strength_profile(candida_support_details(evidence_item))


def candida_invasive_hospital_support_details(
    final_by_name: dict[str, dict[str, Any]],
    name: Any,
    *,
    include_related_representative: bool = False,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    target_key = pathogen_names.canonical_key(name)
    target_group = pathogen_names.representative_group_key(name)
    for evidence_item in final_by_name.values():
        if not isinstance(evidence_item, dict):
            continue
        evidence_name = evidence_item.get("organism_name") or evidence_item.get("name")
        evidence_key = pathogen_names.canonical_key(evidence_name)
        exact = bool(target_key and evidence_key == target_key)
        related = bool(
            include_related_representative
            and target_group
            and target_group != "name:"
            and pathogen_names.representative_group_key(evidence_name) == target_group
        )
        if not exact and not related:
            continue
        for detail in candida_invasive_support_details(evidence_item):
            detail = dict(detail)
            detail["organism_name"] = evidence_name
            detail["support_type"] = "same_species" if exact else "related_representative_context"
            output.append(detail)
    return output


def candida_hospital_support_strength_details(
    final_by_name: dict[str, dict[str, Any]],
    name: Any,
    *,
    include_related_representative: bool = False,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    target_key = pathogen_names.canonical_key(name)
    target_group = pathogen_names.representative_group_key(name)
    for evidence_item in final_by_name.values():
        if not isinstance(evidence_item, dict):
            continue
        evidence_name = evidence_item.get("organism_name") or evidence_item.get("name")
        evidence_key = pathogen_names.canonical_key(evidence_name)
        exact = bool(target_key and evidence_key == target_key)
        related = bool(
            include_related_representative
            and target_group
            and target_group != "name:"
            and pathogen_names.representative_group_key(evidence_name) == target_group
        )
        if not exact and not related:
            continue
        for detail in candida_support_details(evidence_item):
            detail = dict(detail)
            detail["organism_name"] = evidence_name
            detail["support_type"] = "same_species" if exact else "related_representative_context"
            output.append(detail)
    return output


def hospital_only_candidate_from_evidence(
    evidence_item: dict[str, Any],
    ctx: dict[str, Any],
) -> dict[str, Any] | None:
    name = mngs.canonical_display_name(evidence_item.get("organism_name") or evidence_item.get("name"))
    if not name:
        return None
    best_level = evidence_item.get("best_hospital_level")
    best_level_rank = level_rank(best_level)
    impossible_rule = pathogen_rule_lists.impossible_infection_source_rule(name)
    protected_name = mngs.is_protected_pathogen_name(name)
    normalized_name = mngs.normalize_organism_name(name)
    allow_hospital_only_level3 = any(
        normalized_name.startswith(prefix) for prefix in ASPERGILLUS_MOLD_PREFIXES
    )
    direct_modules = direct_hospital_support_modules(evidence_item)
    if best_level_rank > 2:
        direct_modules = direct_hospital_support_modules(evidence_item, max_level_rank=3)
        if not (best_level_rank == 3 and protected_name and allow_hospital_only_level3 and direct_modules):
            return None
    if not direct_modules:
        return None
    classification = str(evidence_item.get("classification") or "Unknown")
    support = {
        "mNGS": "Not_detected_or_not_in_candidate_pool",
        "culture": "Support" if "culture" in direct_modules else "Not_available",
        "filmarray_gmtest": "Support" if "filmarray_gmtest" in direct_modules else "Not_available",
        "molecular_microbiology": "Support" if "molecular_microbiology" in direct_modules else "Not_available",
        "image": "Not_available",
        "host": "Support" if host_support(ctx) else "Not_available",
    }
    if impossible_rule:
        best_level = "Level 5"
        best_level_rank = 5
        direct_modules = []
    level_rule = f"R-HOSPITAL-ONLY-L{level_rank(best_level)}"
    candidate = {
        "organism_name": name,
        "classification": classification,
        "taxonomy_profile": organism_taxonomy.classify_organism(
            name,
            biological_class=classification,
        ),
        "integrated_causative_level": best_level,
        "mngs_signal_tier": "M0_not_detected_or_not_in_candidate_pool",
        "rank_priority": "Not_available",
        "rank_rule": "hospital_only",
        "reads": 0,
        "reads_tier": "Not_available",
        "reads_percentile": 0,
        "dominance_tier": "D0_not_top",
        "specimen_alignment": "Hospital_only",
        "specimen_class": "S3_unknown_or_low_value",
        "source_category": "Hospital_only",
        "evidence_source": "hospital_only",
        "is_protected_pathogen": protected_name,
        "protected_retention_condition": protected_name and best_level_rank <= 3 and (
            best_level_rank <= 2 or allow_hospital_only_level3
        ),
        "is_likely_colonizer_or_background": False,
        "module_support_summary": support,
        "applied_rules": ["D-S1-CANDIDATE_FROM_HOSPITAL_ONLY", level_rule],
        "key_evidence": {
            "best_record_rule": "hospital_only",
            "level_rule": level_rule,
            "guardrail_rule": "",
            "support_modules": direct_modules,
            "hospital_best_level": best_level,
            "hospital_evidence_modules": list(evidence_item.get("evidence_modules") or []),
            "hospital_module_level_summary": evidence_item.get("module_level_summary", {}),
            "hospital_module_evidence": evidence_item.get("module_evidence", {}),
            "hospital_evidence_trace": evidence_item.get("evidence_trace", []),
            "hospital_evidence_detail": hospital_evidence_detail(evidence_item, direct_modules),
        },
        "integrated_reasoning": [
            "deterministic scorer：此菌未出現在 mNGS ranked/chosen 候選池，但醫院端直接病原檢測"
            f"（{', '.join(direct_modules)}）支持 {best_level}，因此以 hospital_only 候選納入 final pick。"
        ],
        "possibility_level": "hospital_only",
    }
    if is_candida_or_generic_yeast_name(name):
        evidence_profile = candida_hospital_evidence_strength_profile(evidence_item)
        candidate["key_evidence"]["candida_evidence_strength_profile"] = evidence_profile
        invasive_support = candida_invasive_support_details(evidence_item)
        if invasive_support:
            candidate["key_evidence"]["candida_invasive_hospital_support"] = invasive_support
            candidate["applied_rules"].append("R-S4-CANDIDA-BLOOD-STERILE-CONTEXT")
            candidate["integrated_reasoning"].append(
                "Candida/yeast has blood or sterile-site hospital evidence, so it is not treated as respiratory Candida colonization. "
                "This preserves invasive/systemic Candida context but does not by itself prove Candida pneumonia."
            )
            candidate["interpretation_caution_zh"] = (
                "Candida/yeast 具有血液或無菌部位證據，不能用呼吸道定植邏輯直接排除；"
                "但仍需臨床/RAG 判斷是否為侵襲性念珠菌感染、血流感染或肺部感染主因。"
            )
        else:
            mark_formal_pick_exclusion(
                candidate,
                rule_id=CANDIDA_NONINVASIVE_HOSPITAL_ONLY_GUARDRAIL,
                reason=(
                    "Hospital-only Candida/yeast lacks same-species strong invasive evidence "
                    "(blood, sterile-site, tissue, pleural, abscess, or histopathology). "
                    "Respiratory/non-sterile Candida is retained as audit/review context rather than formal picked_pathogens."
                ),
            )
            candidate["interpretation_caution_zh"] = (
                "Candida/yeast 若沒有同菌種血液、無菌部位、組織、胸水、膿瘍或病理等強侵襲性證據，"
                "通常較符合定植、背景訊號或需進一步查證；"
                "不列入正式 picked pathogen。"
            )
    if impossible_rule:
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M5_impossible_infection_source"
        candidate["integrated_causative_level"] = "Level 5"
        candidate["key_evidence"]["guardrail_rule"] = "R-S4-IMPOSSIBLE-INFECTION-SOURCE"
        candidate["key_evidence"]["level_rule"] = "R-S5-L5"
        candidate["applied_rules"].extend(["R-S4-IMPOSSIBLE-INFECTION-SOURCE", "R-S3-M5-IMPOSSIBLE", "R-S5-L5"])
        candidate["impossible_infection_source_rule"] = {
            "organism_name": impossible_rule.get("organism_name"),
            "action": impossible_rule.get("action"),
            "rationale_zh": impossible_rule.get("rationale_zh"),
            "notes_zh": impossible_rule.get("notes_zh"),
        }
        candidate["integrated_reasoning"] = [
            "Hospital-side organism is listed in rules/impossible_infection_sources.json as not suitable for pulmonary infection-source output; force Level 5 and exclude from formal picked_pathogens."
        ]
    return candidate


def hospital_only_candidates(
    final_summary: dict[str, Any],
    existing_candidates: Sequence[dict[str, Any]],
    ctx: dict[str, Any],
) -> list[dict[str, Any]]:
    existing_names = {
        mngs.normalize_organism_name(candidate.get("organism_name"))
        for candidate in existing_candidates
        if isinstance(candidate, dict)
    }
    output: list[dict[str, Any]] = []
    for evidence_item in final_summary.get("hospital_organism_evidence") or []:
        if not isinstance(evidence_item, dict):
            continue
        normalized = mngs.normalize_organism_name(
            evidence_item.get("organism_name") or evidence_item.get("name")
        )
        if not normalized or normalized in existing_names:
            continue
        candidate = hospital_only_candidate_from_evidence(evidence_item, ctx)
        if candidate is not None:
            output.append(candidate)
            existing_names.add(normalized)
    return output


def has_same_organism_culture_support(final_by_name: dict[str, dict[str, Any]], name: Any) -> bool:
    candidate = final_candidate_for_name(final_by_name, name)
    if not isinstance(candidate, dict) or not candidate:
        return False
    module_evidence = candidate.get("module_evidence")
    if isinstance(module_evidence, dict):
        rows = module_evidence.get("culture")
        if isinstance(rows, list) and rows:
            return any(evidence_row_is_pulmonary_relevant(row) for row in rows)
    modules = candidate.get("evidence_modules")
    if isinstance(modules, list) and "culture" in {str(item).lower() for item in modules}:
        return True
    levels = candidate.get("module_level_summary")
    if isinstance(levels, dict):
        return mngs._level_rank(levels.get("culture")) <= 4  # pylint: disable=protected-access
    return False


def evidence_row_is_lower_respiratory(row: Any) -> bool:
    if not isinstance(row, dict):
        return False
    category = pathogen_names.raw_key(row.get("specimen_category"))
    if category == "lowerrespiratory":
        return True
    haystack = " ".join(
        lower_text(row.get(key))
        for key in (
            "specimen_type",
            "specimen_category",
            "sample",
            "sample_type",
            "source",
            "site",
            "body_site",
        )
    )
    return any(term in haystack for term in LOWER_RESPIRATORY_TERMS)


def evidence_row_is_pulmonary_relevant(row: Any) -> bool:
    """Whether a hospital result can support causality in a pulmonary model.

    Blood and other sterile-site evidence remains relevant because it can
    corroborate systemic or hematogenous infection. Explicit non-sterile
    nonpulmonary sources such as urine, Foley, or stool remain context only.
    """
    if not isinstance(row, dict):
        return False
    category = pathogen_names.raw_key(row.get("specimen_category"))
    if category in {"lowerrespiratory", "sterilesite"}:
        return True
    if category in {"nonsterileother", "urine", "urinary", "stool", "gastrointestinal"}:
        return False
    haystack = " ".join(
        lower_text(row.get(key))
        for key in (
            "specimen_type",
            "specimen_category",
            "sample",
            "sample_type",
            "source",
            "site",
            "body_site",
        )
    )
    return any(term in haystack for term in (*LOWER_RESPIRATORY_TERMS, *STERILE_TERMS))


def module_has_pulmonary_relevant_evidence(evidence_item: dict[str, Any], module: str) -> bool:
    module_evidence = evidence_item.get("module_evidence")
    if not isinstance(module_evidence, dict):
        return True
    rows = module_evidence.get(module)
    if not isinstance(rows, list) or not rows:
        return True
    return any(evidence_row_is_pulmonary_relevant(row) for row in rows)


def nonpulmonary_culture_context(evidence_item: dict[str, Any]) -> list[dict[str, Any]]:
    module_evidence = evidence_item.get("module_evidence")
    if not isinstance(module_evidence, dict):
        return []
    rows = module_evidence.get("culture")
    if not isinstance(rows, list) or not rows:
        return []
    return [dict(row) for row in rows if isinstance(row, dict) and not evidence_row_is_pulmonary_relevant(row)]


def has_same_organism_lower_respiratory_hospital_support(
    final_by_name: dict[str, dict[str, Any]],
    name: Any,
) -> bool:
    candidate = final_candidate_for_name(final_by_name, name)
    if not isinstance(candidate, dict) or not candidate:
        return False
    module_evidence = candidate.get("module_evidence")
    if not isinstance(module_evidence, dict):
        return False
    for module in DIRECT_HOSPITAL_MODULES:
        rows = module_evidence.get(module)
        if not isinstance(rows, list):
            continue
        for row in rows:
            if evidence_row_is_lower_respiratory(row):
                return True
    return False


def related_representative_hospital_support(
    final_by_name: dict[str, dict[str, Any]],
    name: Any,
) -> list[dict[str, Any]]:
    target_group = pathogen_names.representative_group_key(name)
    target_key = pathogen_names.canonical_key(name)
    if not target_group or not target_key:
        return []

    output: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for evidence_item in final_by_name.values():
        if not isinstance(evidence_item, dict):
            continue
        evidence_name = evidence_item.get("organism_name") or evidence_item.get("name")
        evidence_key = pathogen_names.canonical_key(evidence_name)
        if not evidence_key or evidence_key == target_key:
            continue
        if pathogen_names.representative_group_key(evidence_name) != target_group:
            continue
        levels = evidence_item.get("module_level_summary")
        if not isinstance(levels, dict):
            continue
        modules: list[dict[str, str]] = []
        for module in DIRECT_HOSPITAL_MODULES:
            level = levels.get(module)
            if level_rank(level) <= 4:
                modules.append({"module": module, "level": str(level)})
        if not modules:
            continue
        dedupe_key = (str(evidence_key), ";".join(item["module"] for item in modules))
        if dedupe_key in seen:
            continue
        seen.add(dedupe_key)
        output.append(
            {
                "organism_name": mngs.canonical_display_name(evidence_name),
                "canonical_key": evidence_key,
                "representative_group": target_group,
                "support_type": "related_representative_context",
                "support_modules": modules,
            }
        )
    return output


def related_typical_pneumonia_context_support(candidate: dict[str, Any]) -> bool:
    return (
        candidate.get("classification") == "Bacterial"
        and candidate.get("specimen_class") == "S2_lower_respiratory"
        and is_typical_pneumonia_bacterium_name(candidate.get("organism_name"))
        and clean_mngs_signal(candidate)
        and bool(candidate.get("related_representative_hospital_support"))
    )


def related_representative_context_support(candidate: dict[str, Any]) -> bool:
    if candidate.get("related_representative_hospital_support"):
        return True
    key_evidence = candidate.get("key_evidence")
    return isinstance(key_evidence, dict) and bool(key_evidence.get("related_representative_hospital_support"))


def mngs_ranked_host_only_without_nonhost_context(candidate: dict[str, Any]) -> bool:
    return (
        candidate.get("evidence_source", "mNGS_ranked") == "mNGS_ranked"
        and not same_organism_auxiliary_support(candidate)
        and not related_representative_context_support(candidate)
    )


def mngs_ranked_without_same_organism_auxiliary_support(candidate: dict[str, Any]) -> bool:
    return (
        candidate.get("evidence_source", "mNGS_ranked") == "mNGS_ranked"
        and not same_organism_auxiliary_support(candidate)
    )


def low_actionability_virus_formal_pick_omit(candidate: dict[str, Any]) -> bool:
    return (
        is_low_specificity_virus_name(candidate.get("organism_name"))
        and mngs_ranked_host_only_without_nonhost_context(candidate)
        and candidate.get("dominance_tier") in NON_DOMINANT_TIERS
    )


def water_environmental_low_specificity_formal_pick_omit(candidate: dict[str, Any]) -> bool:
    return (
        is_picked_safe_water_environmental_name(candidate.get("organism_name"))
        and mngs_ranked_host_only_without_nonhost_context(candidate)
    )


def mngs_r0_r1_host_only_nonprotected_formal_pick_omit(candidate: dict[str, Any]) -> bool:
    return (
        mngs_ranked_host_only_without_nonhost_context(candidate)
        and candidate.get("reads_tier") in LOW_SIGNAL_READ_TIERS
        and candidate.get("dominance_tier") in NON_DOMINANT_TIERS
        and not candidate.get("is_protected_pathogen")
        and not is_picked_safe_answer_sensitive_name(candidate.get("organism_name"))
    )


def common_hospital_low_read_nondominant_formal_pick_omit(candidate: dict[str, Any]) -> bool:
    return (
        mngs_ranked_without_same_organism_auxiliary_support(candidate)
        and mngs.normalize_organism_name(candidate.get("organism_name")) in PICKED_SAFE_COMMON_RESPIRATORY_NAMES
        and candidate.get("reads_tier") in LOW_SIGNAL_READ_TIERS
        and candidate.get("dominance_tier") in NON_DOMINANT_TIERS
    )


def weak_mold_host_only_formal_pick_omit(candidate: dict[str, Any]) -> bool:
    return (
        is_picked_safe_mold_opportunistic_name(candidate.get("organism_name"))
        and mngs_ranked_host_only_without_nonhost_context(candidate)
        and candidate.get("reads_tier") in {"R0_trace", "R1_low", "R2_medium"}
        and candidate.get("dominance_tier") in NON_DOMINANT_TIERS
        and mngs.normalize_organism_name(candidate.get("organism_name")) not in PICKED_SAFE_HIGH_CONSEQUENCE_NAMES
    )


def non_core_respiratory_virus_formal_pick_omit(candidate: dict[str, Any]) -> bool:
    return is_non_core_respiratory_virus_name(candidate.get("organism_name"))


def deferred_gi_urinary_colonizer_for_non_typical_guardrail(value: Any) -> bool:
    normalized = mngs.normalize_organism_name(value)
    return normalized.startswith("enterococcus") or normalized in {
        "clostridioidesdifficile",
        "clostridiumdifficile",
    }


def answer_sensitive_healthcare_opportunist_for_non_typical_guardrail(value: Any) -> bool:
    return mngs.normalize_organism_name(value).startswith("acinetobacter")


def mngs_only_nondominant_non_typical_formal_pick_omit(candidate: dict[str, Any]) -> bool:
    name = candidate.get("organism_name")
    return (
        mngs_ranked_host_only_without_nonhost_context(candidate)
        and candidate.get("dominance_tier") in NON_DOMINANT_TIERS
        and not candidate.get("is_protected_pathogen")
        and not is_picked_safe_answer_sensitive_name(name)
        and not deferred_gi_urinary_colonizer_for_non_typical_guardrail(name)
        and not answer_sensitive_healthcare_opportunist_for_non_typical_guardrail(name)
    )


def supported_high_burden_low_specificity_signal(candidate: dict[str, Any]) -> bool:
    """Allow a narrow picked-path exception for low-specificity taxa.

    This is intentionally independent of whether another stronger pathogen is
    already picked. CoNS/Corynebacterium-like organisms should not become
    primary by default, but a same-organism hospital signal plus high mNGS
    burden is enough to report them as secondary Level 3 candidates.
    """
    if not same_organism_auxiliary_support(candidate):
        return False
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads_tier = candidate.get("reads_tier")
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    dominance = candidate.get("dominance_tier")
    specimen_cls = candidate.get("specimen_class")
    if specimen_cls == "S2_lower_respiratory":
        return rank <= 2 and (
            reads_tier in {"R4_very_high", "R3_high"}
            or percentile >= 0.90
            or dominance in {"D2_moderate", "D3_dominant"}
        )
    if specimen_cls == "S1_sterile_systemic":
        return rank <= 2 and (
            reads_tier in {"R4_very_high", "R3_high", "R2_medium", "R1_low"}
            or percentile >= 0.75
        )
    return False


def molecular_support_categories(record: dict[str, Any], classification: str) -> set[str]:
    haystack = " ".join(
        lower_text(record.get(key))
        for key in ("specimen_site", "specimen_name", "specimen_code", "seq_id")
    )
    if any(term in haystack for term in LOWER_RESPIRATORY_TERMS):
        categories = {"Lower_Respiratory"}
        if classification == "Viral":
            categories.add("Upper_Respiratory")
        return categories
    if "csf" in haystack or "cerebrospinal" in haystack:
        return {"CNS", "Sterile_Site"}
    if any(term in haystack for term in STERILE_TERMS):
        return {"Sterile_Site"}
    if any(term in haystack for term in ("urine", "urinary")):
        return {"Urinary"}
    if any(term in haystack for term in ("stool", "feces", "faeces", "rectal")):
        return {"GI"}
    return {"Unknown"}


def has_compatible_molecular_support(
    final_by_name: dict[str, dict[str, Any]],
    name: Any,
    record: dict[str, Any],
    classification: str,
) -> bool:
    candidate = final_candidate_for_name(final_by_name, name)
    if not isinstance(candidate, dict) or not candidate:
        return False
    levels = candidate.get("module_level_summary")
    if not isinstance(levels, dict) or mngs._level_rank(  # pylint: disable=protected-access
        levels.get("molecular_microbiology")
    ) > 4:
        return False
    module_evidence = candidate.get("module_evidence")
    if not isinstance(module_evidence, dict):
        return True
    rows = module_evidence.get("molecular_microbiology")
    if not isinstance(rows, list) or not rows:
        return True
    accepted = molecular_support_categories(record, classification)
    for row in rows:
        if not isinstance(row, dict):
            continue
        category = str(row.get("specimen_category") or "Unknown")
        if category in accepted:
            return True
    return False


def host_support(ctx: dict[str, Any]) -> bool:
    return bool(
        ctx.get("host_level3_expansion")
        or ctx.get("host_vulnerability_tier") in {"V2", "V3"}
        or ctx.get("opportunistic_coverage_level") in {"O1", "O2"}
    )


def candida_resp_l3_condition(
    *,
    candidate: dict[str, Any],
    final_by_name: dict[str, dict[str, Any]],
    ctx: dict[str, Any],
) -> bool:
    name = candidate.get("organism_name")
    if final_candidate_level(final_by_name, name) == 3:
        return True
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads_tier = candidate.get("reads_tier")
    mngs_support = rank <= 2 or reads_tier in {"R4_very_high", "R3_high", "R2_medium"}
    host_or_culture_support = host_support(ctx) or has_same_organism_culture_support(final_by_name, name)
    if host_or_culture_support and mngs_support:
        return True
    return reads_tier == "R4_very_high"


def has_other_level12_strong_pathogen(candidates: list[dict[str, Any]]) -> bool:
    """Detect whether a case already has a non-commensal Level 1/2 explanation."""
    return any(
        level_rank(candidate.get("integrated_causative_level")) <= 2
        and candidate.get("is_likely_colonizer_or_background") is not True
        and candidate.get("_resp_commensal_guardrail_candidate") is not True
        for candidate in candidates
    )


def apply_resp_commensal_case_fallback(candidates: list[dict[str, Any]]) -> None:
    """Promote high-read respiratory commensals only when no stronger pathogen exists."""
    if not has_other_level12_strong_pathogen(candidates):
        for candidate in candidates:
            if candidate.get("_resp_commensal_high_reads_fallback") is not True:
                continue
            if is_oral_aspiration_flora_name(candidate.get("organism_name")):
                continue
            candidate["is_likely_colonizer_or_background"] = False
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["module_support_summary"]["mNGS"] = "M3_protected"
            candidate["integrated_causative_level"] = "Level 3"
            candidate["key_evidence"]["guardrail_rule"] = "R-S4-RESP-COMMENSAL-L3"
            candidate["key_evidence"]["level_rule"] = "R-S5-L3"
            candidate["applied_rules"] = [
                rule
                for rule in candidate.get("applied_rules", [])
                if rule not in {"R-S4-RESP-COMMENSAL-L4", "R-S5-L4"}
            ]
            for rule in ("R-S4-RESP-COMMENSAL-L3", "R-S3-M3", "R-S5-L3"):
                if rule not in candidate["applied_rules"]:
                    candidate["applied_rules"].append(rule)
            candidate["integrated_reasoning"].append(
                "Lower respiratory commensal/oral/gut flora has R3/R4 reads and no other "
                "non-commensal Level 1/2 pathogen in this case, so it is fallback-retained "
                "but capped at M3_protected / Level 3."
            )

    for candidate in candidates:
        candidate.pop("_resp_commensal_guardrail_candidate", None)
        candidate.pop("_resp_commensal_high_reads_fallback", None)


def dominance_by_record(record: dict[str, Any]) -> dict[str, str]:
    candidates = [item for item in record.get("candidates") or [] if isinstance(item, dict)]
    if not candidates:
        return {}
    sorted_candidates = sorted(
        candidates,
        key=lambda item: (
            -mngs.to_float(item.get("reads", item.get("sec_hit", 0))),
            mngs.normalize_organism_name(item.get("organism_name") or item.get("name")),
        ),
    )
    top_reads = mngs.to_float(sorted_candidates[0].get("reads", sorted_candidates[0].get("sec_hit", 0)))
    second_reads = (
        mngs.to_float(sorted_candidates[1].get("reads", sorted_candidates[1].get("sec_hit", 0)))
        if len(sorted_candidates) > 1
        else 0.0
    )
    ratio = top_reads / max(second_reads, 1.0)
    top_tier = "D3_dominant" if ratio >= 10 else "D2_moderate" if ratio >= 3 else "D1_low"
    output: dict[str, str] = {}
    for index, candidate in enumerate(sorted_candidates):
        key = mngs.normalize_organism_name(candidate.get("organism_name") or candidate.get("name"))
        output[key] = top_tier if index == 0 else "D0_not_top"
    return output


def m1_strong(candidate: dict[str, Any]) -> bool:
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    tier = candidate.get("reads_tier")
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    dominance = candidate.get("dominance_tier")
    return (
        (rank in {1, 2} and tier in {"R4_very_high", "R3_high"})
        or (rank in {1, 2} and percentile >= 0.85)
        or (
            dominance == "D3_dominant"
            and tier in {"R4_very_high", "R3_high", "R2_medium"}
            and rank <= 3
        )
    )


def m2_moderate(candidate: dict[str, Any]) -> bool:
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    tier = candidate.get("reads_tier")
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    poss = lower_text(candidate.get("possibility_level"))
    normalized = mngs.normalize_organism_name(candidate.get("organism_name"))

    # Level 2 should require a strong mNGS signal, same-organism hospital support,
    # or protected-pathogen retention. Broad percentile/possibility rules are
    # intentionally handled in M3 to avoid over-picking many non-top organisms.
    if rank in {1, 2} and tier == "R2_medium" and percentile >= 0.75:
        return True
    if (
        candidate.get("dominance_tier") in {"D2_moderate", "D3_dominant"}
        and rank <= 3
        and tier in {"R4_very_high", "R3_high", "R2_medium"}
    ):
        return True
    if (
        rank <= 3
        and tier in {"R4_very_high", "R3_high", "R2_medium", "R1_low"}
        and same_organism_auxiliary_support(candidate)
    ):
        return True
    if (
        candidate.get("is_protected_pathogen")
        and not is_reactivation_virus_candidate(candidate)
        and not is_aspergillus_mold_candidate(candidate)
        and rank <= 3
        and tier in {"R4_very_high", "R3_high", "R2_medium", "R1_low"}
        and poss != "low"
    ):
        return True
    if rank <= 4 and tier in {"R4_very_high", "R3_high"} and poss != "low":
        return True
    if (
        normalized == "pneumocystisjirovecii"
        and rank <= 3
        and poss != "low"
        and tier in {"R2_medium", "R1_low"}
    ):
        return True
    if (
        normalized.startswith(("mucorales", "cunninghamella", "rhizopus", "mucor", "rhizomucor"))
        and rank == 1
        and percentile >= 0.60
    ):
        return True
    return False


def reactivation_virus_m2(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    """CMV/HSV/VZV need stronger context than other protected pathogens."""
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    tier = candidate.get("reads_tier")
    poss = lower_text(candidate.get("possibility_level"))
    specimen_cls = candidate.get("specimen_class")
    if poss == "low":
        return False
    if same_organism_auxiliary_support(candidate) and rank <= 3 and tier in {
        "R4_very_high",
        "R3_high",
        "R2_medium",
        "R1_low",
    }:
        return True
    if herpes_has_dominant_signal(candidate) and rank <= 3 and tier in {
        "R4_very_high",
        "R3_high",
        "R2_medium",
    }:
        return True
    if rank <= 2 and tier == "R1_low" and ctx.get("host_level3_expansion") and specimen_cls == "S1_sterile_systemic":
        return True
    return False


def reactivation_virus_m3(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    tier = candidate.get("reads_tier")
    poss = lower_text(candidate.get("possibility_level"))
    specimen_cls = candidate.get("specimen_class")
    if poss == "low":
        return False
    if herpes_high_reads_nondominant_host_only(candidate):
        return False
    if herpes_has_dominant_signal(candidate) and rank <= 3 and tier in {
        "R4_very_high",
        "R3_high",
        "R2_medium",
        "R1_low",
    }:
        return True
    if same_organism_auxiliary_support(candidate) and rank <= 3 and tier in {
        "R4_very_high",
        "R3_high",
        "R2_medium",
        "R1_low",
    }:
        return True
    if rank <= 2 and tier == "R1_low":
        return True
    if rank == 3 and tier != "R0_trace" and ctx.get("host_level3_expansion"):
        return True
    if (
        rank <= 3
        and tier == "R0_trace"
        and specimen_cls in {"S1_sterile_systemic", "S2_lower_respiratory"}
        and ctx.get("host_level3_expansion")
    ):
        return True
    return False


def aspergillus_mold_m2(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    """Aspergillus needs stronger context than generic protected pathogens for Level 2."""
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    tier = candidate.get("reads_tier")
    poss = lower_text(candidate.get("possibility_level"))
    if poss == "low":
        return False
    if rank <= 2 and tier in {"R4_very_high", "R3_high", "R2_medium"}:
        return True
    if rank <= 2 and tier == "R1_low" and host_support(ctx):
        return True
    if rank <= 3 and tier in {"R4_very_high", "R3_high", "R2_medium"} and same_organism_auxiliary_support(candidate):
        return True
    if rank <= 3 and tier == "R1_low" and same_organism_auxiliary_support(candidate) and host_support(ctx):
        return True
    return False


def aspergillus_mold_m3(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    tier = candidate.get("reads_tier")
    poss = lower_text(candidate.get("possibility_level"))
    if poss == "low":
        return False
    if rank <= 2 and tier == "R1_low":
        return True
    if rank <= 3 and tier == "R2_medium":
        return True
    if rank <= 3 and tier == "R1_low" and host_support(ctx):
        return True
    if rank <= 3 and tier == "R0_trace" and (host_support(ctx) or same_organism_auxiliary_support(candidate)):
        return True
    if rank <= 4 and same_organism_auxiliary_support(candidate) and tier in {"R4_very_high", "R3_high", "R2_medium", "R1_low"}:
        return True
    return False


def m3_protected(candidate: dict[str, Any]) -> bool:
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    tier = candidate.get("reads_tier")
    percentile = mngs.to_float(candidate.get("reads_percentile"))
    poss = lower_text(candidate.get("possibility_level"))
    if candidate.get("protected_retention_condition") is True:
        return True
    if rank in {1, 2} and tier == "R2_medium":
        return True
    if rank in {1, 2} and tier in {"R1_low", "R0_trace"}:
        return True
    if rank <= 3 and percentile >= 0.60 and tier != "R0_trace":
        return True
    if rank <= 3 and poss in {"high", "medium"} and tier in {"R2_medium", "R1_low"}:
        return True
    if rank == 3 and tier in {"R2_medium", "R1_low"} and candidate.get("has_non_host_auxiliary_support") is True:
        return True
    if rank == 4 and tier == "R4_very_high":
        return True
    return False


def m4_weak(candidate: dict[str, Any]) -> bool:
    return (
        mngs.to_rank_int(candidate.get("rank_priority")) in {4, 5}
        or lower_text(candidate.get("possibility_level")) == "low"
        or candidate.get("reads_tier") in {"R1_low", "R0_trace"}
        or candidate.get("specimen_alignment") == "Non_aligned"
    )


def m4_host_level3_condition(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    """Allow host expansion for weak signals only when the weak candidate is still plausible."""
    if not ctx.get("host_level3_expansion"):
        return False
    if is_reactivation_virus_candidate(candidate):
        return False
    if is_aspergillus_mold_candidate(candidate):
        return aspergillus_mold_m3(candidate, ctx)
    if candidate.get("is_likely_colonizer_or_background"):
        return False
    if candidate.get("specimen_alignment") == "Non_aligned":
        return False
    if candidate.get("has_non_host_auxiliary_support") is True:
        return True
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    tier = candidate.get("reads_tier")
    possibility_level = lower_text(candidate.get("possibility_level"))
    if candidate.get("is_protected_pathogen") and rank <= 3:
        return True
    return rank <= 3 and tier == "R1_low" and possibility_level in {"high", "medium"}



def downgrade_candidate_to_level4(candidate: dict[str, Any], guardrail_rule: str, reasoning: str) -> None:
    candidate["is_likely_colonizer_or_background"] = True
    candidate["mngs_signal_tier"] = "M4_weak"
    candidate["module_support_summary"]["mNGS"] = "M4_weak"
    candidate["integrated_causative_level"] = "Level 4"
    candidate["key_evidence"]["guardrail_rule"] = guardrail_rule
    candidate["key_evidence"]["level_rule"] = "R-S5-L4"
    for rule in (guardrail_rule, "R-S3-M4", "R-S5-L4"):
        if rule not in candidate["applied_rules"]:
            candidate["applied_rules"].append(rule)
    candidate["integrated_reasoning"].append(reasoning)


def trace_protected_pathogen_should_not_be_picked(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    if not candidate.get("is_protected_pathogen"):
        return False
    if candidate.get("reads_tier") != "R0_trace":
        return False
    normalized = mngs.normalize_organism_name(candidate.get("organism_name"))
    if normalized == "pneumocystisjirovecii":
        return False
    if same_organism_auxiliary_support(candidate):
        return False
    if is_reactivation_virus_candidate(candidate):
        return False
    return not (host_support(ctx) and candidate.get("dominance_tier") in {"D2_moderate", "D3_dominant"})

def set_candidate_level_from_tier(candidate: dict[str, Any], ctx: dict[str, Any]) -> None:
    tier = candidate.get("mngs_signal_tier")
    alignment = candidate.get("specimen_alignment")
    reads = candidate.get("reads_tier")
    dominance = candidate.get("dominance_tier")
    if tier == "M1_strong":
        if alignment in {"Aligned", "Sterile_or_Systemic"} and (
            dominance == "D3_dominant" or reads in {"R4_very_high", "R3_high"}
        ):
            candidate["integrated_causative_level"] = "Level 1"
            candidate["key_evidence"]["level_rule"] = "R-S5-L1"
        else:
            candidate["integrated_causative_level"] = "Level 2"
            candidate["key_evidence"]["level_rule"] = "R-S5-L2"
    elif tier == "M2_moderate":
        candidate["integrated_causative_level"] = (
            "Level 2" if alignment in {"Aligned", "Sterile_or_Systemic"} else "Level 3"
        )
        candidate["key_evidence"]["level_rule"] = (
            "R-S5-L2" if candidate["integrated_causative_level"] == "Level 2" else "R-S5-L3"
        )
    elif tier == "M3_protected":
        candidate["integrated_causative_level"] = "Level 3"
        candidate["key_evidence"]["level_rule"] = "R-S5-L3"
    elif tier == "M4_weak":
        candidate["integrated_causative_level"] = (
            "Level 3" if m4_host_level3_condition(candidate, ctx) else "Level 4"
        )
        candidate["key_evidence"]["level_rule"] = (
            "R-S5-L3" if candidate["integrated_causative_level"] == "Level 3" else "R-S5-L4"
        )
    elif tier == "M5_background":
        candidate["integrated_causative_level"] = "Level 5"
        candidate["key_evidence"]["level_rule"] = "R-S5-L5"
    else:
        candidate["integrated_causative_level"] = "Unknown"
        candidate["key_evidence"]["level_rule"] = "R-S5-UNKNOWN"


def score_candidate(
    raw_candidate: dict[str, Any],
    record: dict[str, Any],
    dominance_map: dict[str, str],
    final_by_name: dict[str, dict[str, Any]],
    ctx: dict[str, Any],
) -> dict[str, Any]:
    raw_name = raw_candidate.get("organism_name") or raw_candidate.get("name")
    name = mngs.canonical_display_name(raw_name)
    normalized = mngs.normalize_organism_name(name)
    ranking = raw_candidate.get("ranking") if isinstance(raw_candidate.get("ranking"), dict) else {}
    reads_value = raw_candidate.get("reads", raw_candidate.get("sec_hit", 0))
    alignment = specimen_alignment(record)
    specimen_cls = specimen_class(record, alignment)
    rank = mngs.to_rank_int(ranking.get("rank_priority"))
    protected = mngs.is_protected_pathogen_name(name)
    final_candidate = final_candidate_for_name(final_by_name, name)
    module_levels = final_candidate.get("module_level_summary") if isinstance(final_candidate, dict) else {}
    module_support = {
        "mNGS": "Unknown",
        "culture": "Support" if has_same_organism_culture_support(final_by_name, name) else "Not_available",
        "filmarray_gmtest": "Support"
        if isinstance(module_levels, dict) and mngs._level_rank(module_levels.get("filmarray_gmtest")) <= 4  # pylint: disable=protected-access
        else "Not_available",
        "molecular_microbiology": "Support"
        if has_compatible_molecular_support(
            final_by_name,
            name,
            record,
            classify_source_category(raw_candidate.get("source_category")),
        )
        else "Not_available",
        "image": "Not_available",
        "host": "Support" if host_support(ctx) else "Not_available",
    }
    classification = classify_source_category(raw_candidate.get("source_category"))
    candidate = {
        "organism_name": name,
        "classification": classification,
        "taxonomy_profile": organism_taxonomy.classify_organism(
            name,
            biological_class=classification,
        ),
        "integrated_causative_level": "Unknown",
        "mngs_signal_tier": "Unknown",
        "rank_priority": str(rank) if rank != 999 else "Unknown",
        "rank_rule": str(ranking.get("rank_rule") or ""),
        "reads": int(mngs.to_float(reads_value)),
        "reads_tier": mngs.reads_tier(reads_value),
        "reads_percentile": mngs.normalize_number_value(ranking.get("reads_percentile", 0)),
        "dominance_tier": dominance_map.get(mngs.normalize_organism_name(raw_name), "D0_not_top"),
        "specimen_alignment": alignment,
        "specimen_class": specimen_cls,
        "source_category": raw_candidate.get("source_category", "Unknown"),
        "is_protected_pathogen": protected,
        "protected_retention_condition": False,
        "is_likely_colonizer_or_background": False,
        "module_support_summary": module_support,
        "applied_rules": ["D-S1-CANDIDATE_FROM_RANKED"],
        "key_evidence": {
            "best_record_rule": str(ranking.get("rank_rule") or ""),
            "level_rule": "",
            "guardrail_rule": "",
            "support_modules": [key for key, value in module_support.items() if key != "mNGS" and value == "Support"],
        },
        "integrated_reasoning": [],
        "possibility_level": str(ranking.get("possibility_level") or "Unknown"),
    }

    nonpulmonary_culture = (
        nonpulmonary_culture_context(final_candidate) if isinstance(final_candidate, dict) else []
    )
    if nonpulmonary_culture and module_support["culture"] != "Support":
        candidate["key_evidence"]["nonpulmonary_culture_context"] = nonpulmonary_culture
        candidate["key_evidence"]["context_support_rules"] = [NONPULMONARY_CULTURE_CONTEXT_RULE]
        candidate["applied_rules"].append(NONPULMONARY_CULTURE_CONTEXT_RULE)
        candidate["integrated_reasoning"].append(
            "Same-organism culture evidence is present only in a nonpulmonary source, so it is retained as patient context but does not support pulmonary causality."
        )

    candidate["has_auxiliary_support"] = any(
        value == "Support" for key, value in module_support.items() if key != "mNGS"
    )
    candidate["has_non_host_auxiliary_support"] = any(
        module_support.get(key) == "Support"
        for key in ("culture", "filmarray_gmtest", "molecular_microbiology", "image")
    )
    related_support = related_representative_hospital_support(final_by_name, name)
    if related_support:
        candidate["related_representative_hospital_support"] = related_support
        candidate["key_evidence"]["related_representative_hospital_support"] = related_support
        context_rules = candidate["key_evidence"].setdefault("context_support_rules", [])
        if isinstance(context_rules, list) and RELATED_REPRESENTATIVE_SUPPORT_RULE not in context_rules:
            context_rules.append(RELATED_REPRESENTATIVE_SUPPORT_RULE)
        if RELATED_REPRESENTATIVE_SUPPORT_RULE not in candidate["applied_rules"]:
            candidate["applied_rules"].append(RELATED_REPRESENTATIVE_SUPPORT_RULE)
        candidate["integrated_reasoning"].append(
            "Related representative-group hospital evidence is present but is not treated as same-organism direct support."
        )
    if protected:
        candidate["protected_retention_condition"] = (
            rank in {1, 2}
            or (rank == 3 and (candidate["has_non_host_auxiliary_support"] or host_support(ctx)))
            or (rank == 4 and candidate["reads_tier"] == "R4_very_high")
        )

    impossible_rule = pathogen_rule_lists.impossible_infection_source_rule(name)
    if impossible_rule:
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M5_impossible_infection_source"
        candidate["module_support_summary"]["mNGS"] = "M5_impossible_infection_source"
        candidate["integrated_causative_level"] = "Level 5"
        candidate["key_evidence"]["guardrail_rule"] = "R-S4-IMPOSSIBLE-INFECTION-SOURCE"
        candidate["key_evidence"]["level_rule"] = "R-S5-L5"
        candidate["applied_rules"].extend(["R-S4-IMPOSSIBLE-INFECTION-SOURCE", "R-S3-M5-IMPOSSIBLE", "R-S5-L5"])
        candidate["impossible_infection_source_rule"] = {
            "organism_name": impossible_rule.get("organism_name"),
            "action": impossible_rule.get("action"),
            "rationale_zh": impossible_rule.get("rationale_zh"),
            "notes_zh": impossible_rule.get("notes_zh"),
        }
        candidate["integrated_reasoning"].append(
            "This organism is listed in rules/impossible_infection_sources.json as not suitable for pulmonary infection-source output; "
            "force Level 5 and exclude from formal picked_pathogens."
        )
        candidate["interpretation_caution_zh"] = str(impossible_rule.get("rationale_zh") or "????????????????????????????")
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if is_sars_cov2_candidate(candidate):
        candidate["applied_rules"].append("R-S4-SARS_COV_2_PNEUMONIA_GUARDRAIL")
        candidate["key_evidence"]["guardrail_rule"] = "R-S4-SARS_COV_2_PNEUMONIA_GUARDRAIL"
        if same_organism_auxiliary_support(candidate):
            candidate["mngs_signal_tier"] = "M2_moderate"
            candidate["applied_rules"].append("R-S3-M2-SARS-HOSPITAL")
            set_candidate_level_from_tier(candidate, ctx)
            reasoning = (
                "SARS-CoV-2/COVID-19 has direct hospital support plus mNGS signal; "
                "retain as a pneumonia-relevant viral candidate while flagging possible current/recent infection context."
            )
        elif strong_sars_cov2_mngs_signal(candidate):
            candidate["mngs_signal_tier"] = "M2_moderate"
            candidate["applied_rules"].append("R-S3-M2-SARS-STRONG-MNGS")
            set_candidate_level_from_tier(candidate, ctx)
            reasoning = (
                "SARS-CoV-2/COVID-19 has high-clean mNGS signal; retain as picked or strong viral candidate."
            )
        else:
            candidate["mngs_signal_tier"] = "M4_weak"
            candidate["module_support_summary"]["mNGS"] = "M4_weak"
            candidate["integrated_causative_level"] = "Level 4"
            candidate["key_evidence"]["level_rule"] = "R-S5-L4"
            candidate["applied_rules"].extend(["R-S3-M4-SARS-LOW-SIGNAL", "R-S5-L4"])
            reasoning = (
                "SARS-CoV-2/COVID-19 has only low mNGS signal without direct hospital support; "
                "keep out of formal picked_pathogens and leave for safety review/watchlist if needed."
            )
        candidate["module_support_summary"]["mNGS"] = candidate["mngs_signal_tier"]
        if candidate["key_evidence"].get("level_rule") and candidate["key_evidence"].get("level_rule") not in candidate["applied_rules"]:
            candidate["applied_rules"].append(candidate["key_evidence"]["level_rule"])
        candidate["integrated_reasoning"].append(reasoning)
        candidate["interpretation_caution_zh"] = "COVID-19/SARS-CoV-2 可能代表目前或近期感染；若只有低 reads 且缺乏 PCR/抗原/臨床影像支持，不應單獨當作肺炎主因。"
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if is_low_specificity_virus_name(name):
        downgrade_candidate_to_level4(
            candidate,
            "R-S4-LOW-SPECIFICITY-VIRUS-L4",
            "Low-specificity virus such as Torque teno virus is not treated as a pulmonary infection source by mNGS alone; keep it out of formal picked_pathogens and leave only for contextual review if needed.",
        )
        candidate["interpretation_caution_zh"] = "低特異性病毒不適合單獨作為肺炎感染源；若需要解讀，應放在免疫狀態或背景病毒脈絡。"
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if trace_protected_pathogen_should_not_be_picked(candidate, ctx):
        downgrade_candidate_to_level4(
            candidate,
            "R-S4-TRACE-PROTECTED-WATCHLIST-PREFERRED",
            "Protected pathogen has only trace reads without direct same-organism hospital support or D2/D3 dominance with host risk; do not promote to formal picked_pathogens by protected status alone.",
        )
        candidate["interpretation_caution_zh"] = "受保護病原只有 trace reads 且缺乏同菌醫院端或優勢訊號支持；不直接作為正式 picked，較適合由 LLM missing/watchlist 或人工覆核。"
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if is_skin_environmental_yeast_name(name) and not same_organism_auxiliary_support(candidate):
        downgrade_candidate_to_level4(
            candidate,
            "R-S4-SKIN-ENV-YEAST-L4",
            "Skin/environment-associated yeast is mNGS-only without same-organism hospital support; avoid formal picked_pathogens despite rank/percentile signal.",
        )
        candidate["interpretation_caution_zh"] = "皮膚或環境相關 yeast 若只有 mNGS 訊號、沒有同菌醫院端支持，不宜作為肺炎正式 picked pathogen。"
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if is_low_pulmonary_specificity_atypical_name(name) and not same_organism_auxiliary_support(candidate):
        downgrade_candidate_to_level4(
            candidate,
            "R-S4-LOW-PULM-SPECIFICITY-ATYPICAL-L4",
            "Low-pulmonary-specificity atypical/mycoplasma-like organism lacks same-organism hospital support; rank alone is not enough for formal picked_pathogens.",
        )
        candidate["interpretation_caution_zh"] = "低肺部特異性的 atypical/mycoplasma-like 菌若缺乏同菌醫院端支持，不因 rank 或 percentile 單獨進正式 picked。"
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if (
        is_skin_low_specificity_bacteria_name(name)
        and candidate.get("specimen_class") != "S1_sterile_systemic"
        and not same_organism_auxiliary_support(candidate)
    ):
        downgrade_candidate_to_level4(
            candidate,
            "R-S4-SKIN-LOW-SPECIFICITY-BACTERIA-L4",
            "Skin-associated low-specificity bacterium is mNGS-only from a non-sterile specimen without same-organism hospital support; rank alone is not enough for formal picked_pathogens.",
        )
        candidate["interpretation_caution_zh"] = "皮膚相關低特異性菌若來自非無菌檢體且只有 mNGS 訊號，不因 rank 或 percentile 單獨進正式 picked。"
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if typical_pneumonia_hospital_mngs_rule(candidate):
        preserves_m1 = m1_strong(candidate)
        candidate["mngs_signal_tier"] = "M1_strong" if preserves_m1 else "M2_moderate"
        candidate["module_support_summary"]["mNGS"] = candidate["mngs_signal_tier"]
        candidate["key_evidence"]["guardrail_rule"] = "R-S4-TYPICAL_PNEUMONIA_HOSPITAL_MNGS_L2"
        candidate["applied_rules"].extend(
            [
                "R-S4-TYPICAL_PNEUMONIA_HOSPITAL_MNGS_L2",
                "R-S3-M1" if preserves_m1 else "R-S3-M2-HOSPITAL-MNGS",
            ]
        )
        set_candidate_level_from_tier(candidate, ctx)
        candidate["applied_rules"].append(candidate["key_evidence"]["level_rule"])
        candidate["integrated_reasoning"].append(
            "Typical pneumonia bacterium has lower-respiratory culture/FilmArray support plus clean mNGS signal; "
            + (
                "retain the existing M1_strong signal because an evidence-support rule must not downgrade it."
                if preserves_m1
                else "apply an M2_moderate floor so it is not lost because the hospital label and mNGS species label differ."
            )
        )
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    respiratory = is_respiratory(record, alignment)
    sterile = is_sterile_or_systemic(record, alignment)
    if is_rare_opportunistic_yeast_name(name) and respiratory and not sterile:
        if rare_opportunistic_yeast_supported_mngs_hospital_rule(candidate, ctx):
            candidate["is_likely_colonizer_or_background"] = False
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["module_support_summary"]["mNGS"] = "M3_protected"
            candidate["integrated_causative_level"] = "Level 3"
            candidate["key_evidence"]["guardrail_rule"] = "R-S4-RARE-YEAST-HOSPITAL-MNGS-L3"
            candidate["key_evidence"]["level_rule"] = "R-S5-L3"
            candidate["applied_rules"].extend(["R-S4-RARE-YEAST-HOSPITAL-MNGS-L3", "R-S3-M3", "R-S5-L3"])
            candidate["integrated_reasoning"].append(
                "Rare opportunistic yeast has same-organism lower-respiratory hospital support, host vulnerability, "
                "and high-burden or high-percentile clean mNGS signal. It is retained as a Level 3 formal candidate "
                "instead of being handled by the respiratory Candida colonization guardrail."
            )
            candidate["interpretation_caution_zh"] = (
                "Trichosporon 等 rare/opportunistic yeast 不能完全套用呼吸道 Candida 定植規則；"
                "若有同菌下呼吸道培養支持、免疫低下宿主與高 mNGS 訊號，可列為正式候選，但仍需臨床脈絡判讀。"
            )
        else:
            candidate["is_likely_colonizer_or_background"] = False
            candidate["mngs_signal_tier"] = "M4_weak"
            candidate["module_support_summary"]["mNGS"] = "M4_weak"
            candidate["integrated_causative_level"] = "Level 4"
            candidate["key_evidence"]["guardrail_rule"] = "R-S4-RARE-YEAST-MNGS-ONLY-CONTEXT"
            candidate["key_evidence"]["level_rule"] = "R-S5-L4"
            candidate["applied_rules"].extend(["R-S4-RARE-YEAST-MNGS-ONLY-CONTEXT", "R-S3-M4", "R-S5-L4"])
            candidate["integrated_reasoning"].append(
                "Rare opportunistic yeast lacks same-organism hospital support or strong enough contextual support; "
                "keep out of formal picked_pathogens and leave for review/RAG context if clinically needed."
            )
            candidate["interpretation_caution_zh"] = (
                "Trichosporon 等 rare/opportunistic yeast 若只有 mNGS 訊號，不能直接作為正式肺炎病原；"
                "應降為 review/RAG 脈絡或低優先候選。"
            )
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    candida_like = is_candida_or_generic_yeast_name(name)
    if candida_like and respiratory and not sterile:
        if candida_resp_l3_condition(candidate=candidate, final_by_name=final_by_name, ctx=ctx):
            exact_strength_details = candida_hospital_support_strength_details(final_by_name, name)
            related_strength_details = [
                detail
                for detail in candida_hospital_support_strength_details(
                    final_by_name,
                    name,
                    include_related_representative=True,
                )
                if detail.get("support_type") == "related_representative_context"
            ]
            exact_strength_profile = candida_evidence_strength_profile(exact_strength_details)
            related_strength_profile = candida_evidence_strength_profile(related_strength_details)
            exact_invasive_support = candida_invasive_hospital_support_details(final_by_name, name)
            related_invasive_context = candida_invasive_hospital_support_details(
                final_by_name,
                name,
                include_related_representative=True,
            )
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["module_support_summary"]["mNGS"] = "M3_protected"
            candidate["integrated_causative_level"] = "Level 3"
            candidate["key_evidence"]["guardrail_rule"] = "R-S4-CANDIDA-RESP-L3"
            candidate["key_evidence"]["level_rule"] = "R-S5-L3"
            candidate["key_evidence"]["candida_evidence_strength_profile"] = exact_strength_profile
            if related_strength_profile.get("best_strength") != "none":
                candidate["key_evidence"]["candida_related_evidence_strength_profile"] = related_strength_profile
            candidate["applied_rules"].extend(["R-S4-CANDIDA-RESP-L3", "R-S3-M3", "R-S5-L3"])
            if exact_invasive_support:
                candidate["key_evidence"]["candida_invasive_hospital_support"] = exact_invasive_support
                candidate["applied_rules"].append("R-S4-CANDIDA-BLOOD-STERILE-CONTEXT")
            elif related_invasive_context:
                candidate["key_evidence"]["candida_related_invasive_hospital_context"] = related_invasive_context
            candidate["integrated_reasoning"].append(
                "Respiratory Candida/yeast meets context-retention criteria and is capped at M3_protected / Level 3. "
                "Respiratory Candida alone is not formal picked unless same-species blood, sterile-site, tissue, pleural, abscess, or histopathology evidence is present."
            )
            candidate["interpretation_caution_zh"] = (
                "呼吸道 Candida/yeast 通常較像定植或背景訊號；若無同菌血液、無菌部位、組織、胸水、膿瘍或病理證據，"
                "不列入正式 picked pathogen，只保留為 review/RAG 脈絡。"
            )
            if not exact_invasive_support:
                mark_formal_pick_exclusion(
                    candidate,
                    rule_id=CANDIDA_RESP_CONTEXT_NOT_PICKED_GUARDRAIL,
                    reason=(
                        "Respiratory Candida/yeast is retained as Level 3 context but excluded from formal picked_pathogens "
                        "because there is no same-species blood, sterile-site, tissue, pleural, abscess, or histopathology evidence."
                    ),
                )
            candidate.pop("has_auxiliary_support", None)
            candidate.pop("has_non_host_auxiliary_support", None)
            return candidate
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M5_background"
        candidate["module_support_summary"]["mNGS"] = "M5_background"
        candidate["integrated_causative_level"] = "Level 5"
        candidate["key_evidence"]["guardrail_rule"] = "R-S4-BACKGROUND"
        candidate["key_evidence"]["level_rule"] = "R-S5-L5"
        candidate["applied_rules"].extend(["R-S4-BACKGROUND", "R-S3-M5", "R-S5-L5"])
        candidate["integrated_reasoning"].append(
            "呼吸道 Candida 未符合 Level 3 條件性保留規則，視為 background/colonizer。"
        )
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if is_coagulase_negative_staph_background_name(name) and not sterile:
        if supported_high_burden_low_specificity_signal(candidate):
            candidate["is_likely_colonizer_or_background"] = False
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["module_support_summary"]["mNGS"] = "M3_protected"
            candidate["integrated_causative_level"] = "Level 3"
            candidate["key_evidence"]["guardrail_rule"] = "R-S4-CONS-SKIN_FLORA_SUPPORTED_L3"
            candidate["key_evidence"]["level_rule"] = "R-S5-L3"
            candidate["applied_rules"].extend(["R-S4-CONS-SKIN_FLORA_SUPPORTED_L3", "R-S3-M3", "R-S5-L3"])
            candidate["integrated_reasoning"].append(
                "CoNS/skin flora has same-organism hospital support plus high-burden mNGS signal; "
                "report as a secondary Level 3 candidate rather than suppressing solely because another pathogen is present."
            )
            candidate.pop("has_auxiliary_support", None)
            candidate.pop("has_non_host_auxiliary_support", None)
            return candidate
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M4_weak"
        candidate["module_support_summary"]["mNGS"] = "M4_weak"
        candidate["integrated_causative_level"] = "Level 4"
        candidate["key_evidence"]["guardrail_rule"] = "R-S4-CONS-SKIN_FLORA_HARD"
        candidate["key_evidence"]["level_rule"] = "R-S5-L4"
        candidate["applied_rules"].extend(["R-S4-CONS-SKIN_FLORA_HARD", "R-S5-L4"])
        candidate["integrated_reasoning"].append(
            "Non-sterile respiratory CoNS/skin flora is capped at M4_weak / Level 4 unless same-organism hospital support and high-burden mNGS signal justify Level 3 review."
        )
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if specimen_cls == "S2_lower_respiratory" and is_oral_aspiration_flora_name(name) and not sterile:
        has_same_lower_respiratory_support = has_same_organism_lower_respiratory_hospital_support(
            final_by_name,
            name,
        )
        has_dominance = candidate.get("dominance_tier") in {"D2_moderate", "D3_dominant"}
        has_meaningful_reads = candidate.get("reads_tier") not in {"R0_trace", "Unknown", "Not_available"}
        candidate["_resp_commensal_guardrail_candidate"] = True
        if has_same_lower_respiratory_support or (has_dominance and has_meaningful_reads):
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["module_support_summary"]["mNGS"] = "M3_protected"
            candidate["integrated_causative_level"] = "Level 3"
            candidate["key_evidence"]["guardrail_rule"] = ORAL_ASPIRATION_FLORA_L3_GUARDRAIL
            candidate["key_evidence"]["level_rule"] = "R-S5-L3"
            candidate["applied_rules"].extend([ORAL_ASPIRATION_FLORA_L3_GUARDRAIL, "R-S3-M3", "R-S5-L3"])
            candidate["integrated_reasoning"].append(
                "Oral/aspiration flora has same-organism lower-respiratory hospital support or D2/D3 dominant "
                "mNGS signal, so it is retained but capped at M3_protected / Level 3."
            )
            candidate.pop("_resp_commensal_guardrail_candidate", None)
            candidate.pop("has_auxiliary_support", None)
            candidate.pop("has_non_host_auxiliary_support", None)
            return candidate
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M4_weak"
        candidate["module_support_summary"]["mNGS"] = "M4_weak"
        candidate["integrated_causative_level"] = "Level 4"
        candidate["key_evidence"]["guardrail_rule"] = ORAL_ASPIRATION_FLORA_L4_GUARDRAIL
        candidate["key_evidence"]["level_rule"] = "R-S5-L4"
        candidate["applied_rules"].extend([ORAL_ASPIRATION_FLORA_L4_GUARDRAIL, "R-S5-L4"])
        candidate["integrated_reasoning"].append(
            "Single oral/aspiration flora from lower-respiratory mNGS without same-organism lower-respiratory "
            "hospital support and without D2/D3 dominant signal is capped at M4_weak / Level 4. "
            "It should be evaluated as an aspiration pattern in review/RAG rather than auto-picked as a standalone pathogen."
        )
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if specimen_cls == "S2_lower_respiratory" and is_respiratory_commensal_flora_name(name) and not sterile:
        has_same_support = same_organism_auxiliary_support(candidate)
        has_high_reads = candidate.get("reads_tier") in {"R4_very_high", "R3_high"}
        has_dominance = candidate.get("dominance_tier") in {"D2_moderate", "D3_dominant"}
        candidate["_resp_commensal_guardrail_candidate"] = True
        if has_same_support or (has_high_reads and has_dominance):
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["module_support_summary"]["mNGS"] = "M3_protected"
            candidate["integrated_causative_level"] = "Level 3"
            candidate["key_evidence"]["guardrail_rule"] = "R-S4-RESP-COMMENSAL-L3"
            candidate["key_evidence"]["level_rule"] = "R-S5-L3"
            candidate["applied_rules"].extend(["R-S4-RESP-COMMENSAL-L3", "R-S3-M3", "R-S5-L3"])
            candidate["integrated_reasoning"].append(
                "Lower respiratory commensal/gut flora has same-organism culture/FilmArray support "
                "or R3/R4 reads with D2/D3 dominance, so it is retained but capped at "
                "M3_protected / Level 3."
            )
            candidate.pop("_resp_commensal_guardrail_candidate", None)
            candidate.pop("has_auxiliary_support", None)
            candidate.pop("has_non_host_auxiliary_support", None)
            return candidate
        if has_high_reads:
            candidate["_resp_commensal_high_reads_fallback"] = True
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M4_weak"
        candidate["module_support_summary"]["mNGS"] = "M4_weak"
        candidate["integrated_causative_level"] = "Level 4"
        candidate["key_evidence"]["guardrail_rule"] = "R-S4-RESP-COMMENSAL-L4"
        candidate["key_evidence"]["level_rule"] = "R-S5-L4"
        candidate["applied_rules"].extend(["R-S4-RESP-COMMENSAL-L4", "R-S5-L4"])
        candidate["integrated_reasoning"].append(
            "Lower respiratory commensal/gut flora without same-organism culture/FilmArray support "
            "and without D2/D3 dominant R3/R4 reads is capped at M4_weak / Level 4. "
            "If it has R3/R4 reads, it can be fallback-retained only when the case has no other "
            "non-commensal Level 1/2 pathogen."
        )
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    low_specificity_environmental = is_low_pulmonary_specificity_environmental_name(name)
    water_environmental_gnb = is_water_environmental_gnb_name(name)
    weak_sterile_low_specificity_environmental = (
        sterile
        and low_specificity_environmental
        and candidate.get("reads_tier") in {"R0_trace", "R1_low", "R2_medium"}
    )
    if (
        (low_specificity_environmental or water_environmental_gnb)
        and (not sterile or weak_sterile_low_specificity_environmental)
        and not same_organism_auxiliary_support(candidate)
        and candidate.get("dominance_tier") not in {"D2_moderate", "D3_dominant"}
    ):
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M4_weak"
        candidate["module_support_summary"]["mNGS"] = "M4_weak"
        candidate["integrated_causative_level"] = "Level 4"
        candidate["key_evidence"]["guardrail_rule"] = (
            "R-S4-WATER-ENV-GNB-L4"
            if water_environmental_gnb
            else "R-S4-LOW-PULM-SPECIFICITY-ENV-L4"
        )
        candidate["key_evidence"]["level_rule"] = "R-S5-L4"
        candidate["applied_rules"].extend([candidate["key_evidence"]["guardrail_rule"], "R-S5-L4"])
        candidate["integrated_reasoning"].append(
            "Low-pulmonary-specificity environmental or water-associated organism is mNGS-only "
            "without D2/D3 dominance or same-organism culture/FilmArray/molecular support, so it is "
            "downgraded to M4_weak / Level 4 and kept out of formal picked_pathogens. It remains "
            "available for safety review/watchlist if clinically needed."
        )
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if mngs.is_background_pathogen_name(name) and not protected:
        if supported_high_burden_low_specificity_signal(candidate):
            guardrail = "R-S4-BACKGROUND_SUPPORTED_L3"
            candidate["is_likely_colonizer_or_background"] = False
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["module_support_summary"]["mNGS"] = "M3_protected"
            candidate["integrated_causative_level"] = "Level 3"
            candidate["key_evidence"]["guardrail_rule"] = guardrail
            candidate["key_evidence"]["level_rule"] = "R-S5-L3"
            candidate["applied_rules"].extend([guardrail, "R-S3-M3", "R-S5-L3"])
            candidate["integrated_reasoning"].append(
                "Background-risk organism such as Corynebacterium has same-organism hospital support plus high-burden mNGS signal; "
                "retain as a secondary Level 3 candidate without requiring repeated culture."
            )
            candidate.pop("has_auxiliary_support", None)
            candidate.pop("has_non_host_auxiliary_support", None)
            return candidate
        guardrail = "R-S4-BACKGROUND"
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M5_background"
        candidate["module_support_summary"]["mNGS"] = candidate["mngs_signal_tier"]
        candidate["integrated_causative_level"] = "Level 5"
        candidate["key_evidence"]["guardrail_rule"] = guardrail
        candidate["key_evidence"]["level_rule"] = "R-S5-L5"
        candidate["applied_rules"].extend([guardrail, candidate["key_evidence"]["level_rule"]])
        candidate.pop("has_auxiliary_support", None)
        candidate.pop("has_non_host_auxiliary_support", None)
        return candidate

    if is_reactivation_virus_candidate(candidate):
        candidate["applied_rules"].append("R-S4-HERPES_REACTIVATION_GUARDRAIL")
        candidate["key_evidence"]["guardrail_rule"] = "R-S4-HERPES_REACTIVATION_GUARDRAIL"
        if herpes_high_reads_nondominant_host_only(candidate):
            downgrade_candidate_to_level4(
                candidate,
                "R-S4-HERPES-R4-NONDOMINANT-HOST-ONLY-REVIEW",
                "CMV/HSV/VZV has R4 very high reads but D0/D1 non-dominant pattern and no direct same-organism hospital/PCR support; keep out of formal picked_pathogens and send to review/RAG context instead of ignoring it.",
            )
            candidate["interpretation_caution_zh"] = (
                "CMV/HSV/VZV 高 reads 但非優勢且缺乏同菌 PCR/viral load/醫院端支持；不可忽略，"
                "但較適合進 review_context_needed 由 RAG/臨床脈絡判斷是否為 shedding/reactivation 或真正肺炎。"
            )
        elif herpes_ranked_host_context_review(candidate, ctx):
            candidate["is_likely_colonizer_or_background"] = False
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["module_support_summary"]["mNGS"] = "M3_protected"
            candidate["integrated_causative_level"] = "Level 3"
            candidate["key_evidence"]["level_rule"] = "R-S5-L3"
            candidate["applied_rules"].extend(["R-S3-M3-HERPES-HOST-CONTEXT", "R-S5-L3"])
            mark_formal_pick_exclusion(
                candidate,
                rule_id=HERPES_HOST_R2R3_NONDOMINANT_REVIEW_RULE,
                reason=(
                    "HSV/CMV/VZV has rank <=2, R2/R3 reads, high reads percentile, host vulnerability, "
                    "and D0/D1 non-dominant host-only support. Keep RAG-visible as Level 3 context, "
                    "but exclude from formal picked_pathogens without PCR/viral-load/direct hospital support."
                ),
            )
            candidate["interpretation_caution_zh"] = (
                "HSV/CMV/VZV 在免疫低下宿主有 rank/reads/percentile 支持時不可忽略；"
                "但若為 host-only、非優勢且沒有 PCR/viral load/院方同菌證據，"
                "應進 review_context_needed/RAG，而不是直接列為正式 picked pathogen。"
            )
        elif reactivation_virus_m2(candidate, ctx):
            candidate["mngs_signal_tier"] = "M2_moderate"
            candidate["applied_rules"].append("R-S3-M2-HERPES")
        elif reactivation_virus_m3(candidate, ctx):
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["applied_rules"].append("R-S3-M3-HERPES")
        elif m4_weak(candidate):
            candidate["mngs_signal_tier"] = "M4_weak"
            candidate["applied_rules"].append("R-S3-M4")
        else:
            candidate["mngs_signal_tier"] = "Unknown"
            candidate["applied_rules"].append("R-S3-UNKNOWN")
    elif is_aspergillus_mold_candidate(candidate):
        candidate["applied_rules"].append("R-S4-ASPERGILLUS_MOLD_GUARDRAIL")
        candidate["key_evidence"]["guardrail_rule"] = "R-S4-ASPERGILLUS_MOLD_GUARDRAIL"
        if aspergillus_mold_m2(candidate, ctx):
            candidate["mngs_signal_tier"] = "M2_moderate"
            candidate["applied_rules"].append("R-S3-M2-MOLD")
        elif aspergillus_mold_m3(candidate, ctx):
            candidate["mngs_signal_tier"] = "M3_protected"
            candidate["applied_rules"].append("R-S3-M3-MOLD")
        elif m4_weak(candidate):
            candidate["mngs_signal_tier"] = "M4_weak"
            candidate["applied_rules"].append("R-S3-M4")
        else:
            candidate["mngs_signal_tier"] = "Unknown"
            candidate["applied_rules"].append("R-S3-UNKNOWN")
    elif m1_strong(candidate):
        candidate["mngs_signal_tier"] = "M1_strong"
        candidate["applied_rules"].append("R-S3-M1")
    elif m2_moderate(candidate):
        candidate["mngs_signal_tier"] = "M2_moderate"
        candidate["applied_rules"].append("R-S3-M2")
    elif m3_protected(candidate):
        candidate["mngs_signal_tier"] = "M3_protected"
        candidate["applied_rules"].append("R-S3-M3")
    elif m4_weak(candidate):
        candidate["mngs_signal_tier"] = "M4_weak"
        candidate["applied_rules"].append("R-S3-M4")
    else:
        candidate["mngs_signal_tier"] = "Unknown"
        candidate["applied_rules"].append("R-S3-UNKNOWN")

    candidate["module_support_summary"]["mNGS"] = candidate["mngs_signal_tier"]
    set_candidate_level_from_tier(candidate, ctx)
    candidate["applied_rules"].append(candidate["key_evidence"]["level_rule"])
    candidate["integrated_reasoning"].append(
        f"deterministic scorer：rank_priority={candidate['rank_priority']}、reads_tier={candidate['reads_tier']}、"
        f"reads_percentile={candidate['reads_percentile']} -> {candidate['mngs_signal_tier']} / {candidate['integrated_causative_level']}。"
    )
    candidate.pop("has_auxiliary_support", None)
    candidate.pop("has_non_host_auxiliary_support", None)
    return candidate


def level_rank(value: Any) -> int:
    return mngs._level_rank(value)  # pylint: disable=protected-access


def candidate_sort_key(candidate: dict[str, Any]) -> tuple[Any, ...]:
    return (
        level_rank(candidate.get("integrated_causative_level")),
        mngs.to_rank_int(candidate.get("rank_priority")),
        -mngs.to_float(candidate.get("reads")),
        mngs.normalize_organism_name(candidate.get("organism_name")),
    )


def reads_tier_score(value: Any) -> int:
    return {
        "R4_very_high": 4,
        "R3_high": 3,
        "R2_medium": 2,
        "R1_low": 1,
        "R0_trace": 0,
    }.get(str(value or ""), -1)


def dominance_tier_score(value: Any) -> int:
    return {
        "D3_dominant": 3,
        "D2_moderate": 2,
        "D1_low": 1,
        "D0_not_top": 0,
    }.get(str(value or ""), -1)


def representative_pick_score(candidate: dict[str, Any]) -> tuple[Any, ...]:
    """Choose one formal representative while preserving related evidence.

    When species-level mNGS evidence exists, it must not be hidden by a broader
    hospital group/genus label. Prefer exact-species hospital + mNGS support,
    then species mNGS + related group/genus hospital support, and only then use
    the previous level/signal ordering. This keeps group evidence as support
    without replacing a more informative species-level output.
    """
    name = candidate.get("organism_name")
    is_specific_species = not pathogen_names.is_generic_representative_label(name)
    is_mngs_ranked = candidate.get("evidence_source", "mNGS_ranked") == "mNGS_ranked"
    if is_specific_species and is_mngs_ranked and same_organism_auxiliary_support(candidate):
        species_evidence_priority = 2
    elif is_specific_species and is_mngs_ranked and related_representative_context_support(candidate):
        species_evidence_priority = 1
    else:
        species_evidence_priority = 0
    return (
        species_evidence_priority,
        100 - level_rank(candidate.get("integrated_causative_level")) * 10,
        5 if is_mngs_ranked else 0,
        4 if same_organism_auxiliary_support(candidate) else 0,
        dominance_tier_score(candidate.get("dominance_tier")),
        reads_tier_score(candidate.get("reads_tier")),
        mngs.to_float(candidate.get("reads")),
        pathogen_names.representative_specificity_score(candidate.get("organism_name")),
        -mngs.to_rank_int(candidate.get("rank_priority")),
        str(candidate.get("organism_name") or ""),
    )


def candidate_support_modules(candidate: dict[str, Any]) -> list[str]:
    key_evidence = candidate.get("key_evidence")
    if isinstance(key_evidence, dict):
        values = key_evidence.get("support_modules")
        if isinstance(values, list):
            return [str(value) for value in values if str(value)]
    return []


def add_same_representative_group_evidence(
    retained: dict[str, Any],
    omitted: dict[str, Any],
    *,
    group_key: str,
    reason: str,
) -> None:
    key_evidence = retained.get("key_evidence")
    if not isinstance(key_evidence, dict):
        key_evidence = {}
        retained["key_evidence"] = key_evidence
    merged = key_evidence.setdefault("same_representative_group_merged_picked_candidates", [])
    if not isinstance(merged, list):
        merged = []
        key_evidence["same_representative_group_merged_picked_candidates"] = merged
    merged.append(
        {
            "organism_name": omitted.get("organism_name"),
            "classification": omitted.get("classification"),
            "observed_level": omitted.get("integrated_causative_level"),
            "evidence_source": omitted.get("evidence_source"),
            "reads": omitted.get("reads", 0),
            "rank_priority": omitted.get("rank_priority"),
            "support_modules": candidate_support_modules(omitted),
            "representative_group": group_key,
            "reason": reason,
        }
    )


def suppress_same_representative_group_candidate(
    omitted: dict[str, Any],
    *,
    group_key: str,
) -> bool:
    """Return True only when representative grouping should hide this item.

    Same-genus grouping is useful for broad labels such as "Aspergillus spp." or
    "Yeast", but different species should remain independently visible. Candida
    is stricter: generic Yeast/Candida labels may be folded into a species, but
    Candida species must not suppress each other.
    """
    name = omitted.get("organism_name")
    if group_key == "yeast:candida_or_generic_yeast":
        return pathogen_names.is_generic_representative_label(name)
    if group_key.startswith("alias:"):
        return True
    if group_key.startswith("genus:"):
        return True
    return False


def preferred_representative_for_group(group_key: str, group: list[dict[str, Any]]) -> dict[str, Any]:
    if group_key == "yeast:candida_or_generic_yeast":
        species = [
            candidate
            for candidate in group
            if not pathogen_names.is_generic_representative_label(candidate.get("organism_name"))
        ]
        if species:
            return max(species, key=representative_pick_score)
    return max(group, key=representative_pick_score)


def apply_same_representative_group_guardrail(
    candidates: list[dict[str, Any]],
    ctx: dict[str, Any],
) -> None:
    groups: dict[str, list[dict[str, Any]]] = {}
    for candidate in candidates:
        if not pickable(candidate, ctx):
            continue
        group_key = pathogen_names.representative_group_key(candidate.get("organism_name"))
        if not group_key or group_key == "name:":
            continue
        groups.setdefault(group_key, []).append(candidate)

    for group_key, group in groups.items():
        if len(group) <= 1:
            continue
        retained = preferred_representative_for_group(group_key, group)
        retained_name = str(retained.get("organism_name") or "")
        for omitted in group:
            if omitted is retained:
                continue
            if not suppress_same_representative_group_candidate(omitted, group_key=group_key):
                continue
            omitted_name = str(omitted.get("organism_name") or "")
            reason = (
                f"Same representative group ({group_key}) already has a stronger formal picked representative: "
                f"{retained_name}. Keep {omitted_name} as related/audit evidence rather than a separate picked pathogen."
            )
            omitted["same_representative_group_retained"] = retained_name
            omitted["same_representative_group_key"] = group_key
            mark_formal_pick_exclusion(
                omitted,
                rule_id=SAME_REPRESENTATIVE_GROUP_GUARDRAIL,
                reason=reason,
            )
            add_same_representative_group_evidence(retained, omitted, group_key=group_key, reason=reason)


def conservative_l3_mngs_only_low_signal(candidate: dict[str, Any]) -> bool:
    if level_rank(candidate.get("integrated_causative_level")) != 3:
        return False
    if candidate.get("evidence_source", "mNGS_ranked") != "mNGS_ranked":
        return False
    if same_organism_auxiliary_support(candidate):
        return False
    if related_typical_pneumonia_context_support(candidate):
        return False
    if candidate.get("reads_tier") not in LOW_SIGNAL_READ_TIERS:
        return False
    if candidate.get("dominance_tier") not in NON_DOMINANT_TIERS:
        return False
    return True


def apply_formal_pick_guardrails(candidates: list[dict[str, Any]], ctx: dict[str, Any]) -> None:
    for candidate in candidates:
        if candidate.get("evidence_source") == "hospital_only" and is_generic_unidentified_hospital_label(
            candidate.get("organism_name")
        ):
            mark_formal_pick_exclusion(
                candidate,
                rule_id=GENERIC_HOSPITAL_LABEL_GUARDRAIL,
                reason=(
                    "Generic unidentified hospital-side label is kept as culture/microbiology evidence "
                    "but excluded from formal picked_pathogens until species-level identification is available."
                ),
            )
            continue
        if conservative_l3_mngs_only_low_signal(candidate):
            mark_formal_pick_exclusion(
                candidate,
                rule_id=L3_MNGS_ONLY_LOW_SIGNAL_GUARDRAIL,
                reason=(
                    "Level 3 mNGS-only signal has no pathogen-specific hospital support, R0/R1 reads, "
                    "and D0/D1 non-dominant pattern; keep for audit/review context rather than formal picked_pathogens."
                ),
            )
        if low_actionability_virus_formal_pick_omit(candidate):
            mark_formal_pick_exclusion(
                candidate,
                rule_id=LOW_ACTIONABILITY_VIRUS_FORMAL_PICK_GUARDRAIL,
                reason=(
                    "Low-actionability viral signal such as HPV/TTV/Anellovirus is mNGS-only, "
                    "host-context, and non-dominant. Keep it as audit/context evidence rather than "
                    "a formal picked pulmonary pathogen."
                ),
            )
            candidate["interpretation_caution_zh"] = (
                "HPV、TTV/Anellovirus 類訊號通常不是肺炎致病源；若偵測到，較適合作為背景或免疫狀態脈絡，"
                "不列入正式 picked pathogen。"
            )
        if water_environmental_low_specificity_formal_pick_omit(candidate):
            mark_formal_pick_exclusion(
                candidate,
                rule_id=WATER_ENVIRONMENTAL_LOW_SPECIFICITY_FORMAL_PICK_GUARDRAIL,
                reason=(
                    "Water/environmental low-pulmonary-specificity organism lacks same-organism culture, "
                    "FilmArray, molecular, or related-representative hospital support. Even with a strong "
                    "mNGS burden, it should be reviewed as context rather than formal picked without local support."
                ),
            )
            candidate["interpretation_caution_zh"] = (
                "環境/水源相關低肺炎特異性菌若缺乏同菌醫院端支持，即使 mNGS reads 或 dominance 較高，"
                "也先不列入正式 picked；可保留給 review/RAG 判斷。"
            )
        if mngs_r0_r1_host_only_nonprotected_formal_pick_omit(candidate):
            mark_formal_pick_exclusion(
                candidate,
                rule_id=MNGS_R0_R1_HOST_ONLY_NONPROTECTED_FORMAL_PICK_GUARDRAIL,
                reason=(
                    "mNGS-only R0/R1 host-context signal is D0/D1 non-dominant and is not in protected, "
                    "answer-sensitive respiratory virus, or common hospital pneumonia pathogen groups. "
                    "Keep it out of formal picked_pathogens."
                ),
            )
        if common_hospital_low_read_nondominant_formal_pick_omit(candidate):
            mark_formal_pick_exclusion(
                candidate,
                rule_id=COMMON_HOSPITAL_LOW_READ_NODOM_FORMAL_PICK_GUARDRAIL,
                reason=(
                    "Common hospital pneumonia pathogen has only R0/R1 mNGS signal, D0/D1 non-dominant pattern, "
                    "and no direct same-organism culture, FilmArray, or molecular support. Keep it out of formal "
                    "picked_pathogens; retain as review/audit context depending on related hospital context."
                ),
            )
            candidate["interpretation_caution_zh"] = (
                "常見院內肺炎菌若只有低 reads mNGS、D0/D1 非優勢且沒有直接同菌院端支持，"
                "不直接列入正式 picked；若有 related group 院端脈絡可送 review，否則作為低特異性 audit。"
            )
        if weak_mold_host_only_formal_pick_omit(candidate):
            mark_formal_pick_exclusion(
                candidate,
                rule_id=WEAK_MOLD_HOST_ONLY_FORMAL_PICK_GUARDRAIL,
                reason=(
                    "Mold/opportunistic fungal signal is host-context only, R0/R1/R2, and D0/D1 non-dominant "
                    "without culture, GM, molecular, sterile-site, or related-representative support. "
                    "Keep it out of formal picked_pathogens and leave for review/RAG context if needed."
                ),
            )
            candidate["interpretation_caution_zh"] = (
                "弱 mold/opportunistic fungus 訊號若只有 mNGS、reads 不高且非優勢，沒有培養/GM/PCR/無菌部位支持時，"
                "不直接列為正式 picked pathogen。"
            )
        if non_core_respiratory_virus_formal_pick_omit(candidate):
            mark_formal_pick_exclusion(
                candidate,
                rule_id=NON_CORE_RESPIRATORY_VIRUS_FORMAL_PICK_GUARDRAIL,
                reason=(
                    "Non-core respiratory virus detection is clinically relevant but often reflects "
                    "recent infection, shedding, coinfection, or syndrome context rather than a formal "
                    "single-causative picked pathogen. Keep it RAG-visible/context instead of formal picked."
                ),
            )
            candidate["interpretation_caution_zh"] = (
                "非核心呼吸道病毒如 parainfluenza/rhinovirus/RSV/human respirovirus 可有臨床意義，"
                "但較需要結合症狀、時序、影像與共感染脈絡；先不列入正式 picked，改送 review/RAG 判斷。"
            )
        if mngs_only_nondominant_non_typical_formal_pick_omit(candidate):
            mark_formal_pick_exclusion(
                candidate,
                rule_id=MNGS_ONLY_NONDOMINANT_NON_TYPICAL_FORMAL_PICK_GUARDRAIL,
                reason=(
                    "mNGS-only D0/D1 non-dominant signal lacks same-organism hospital support and is outside "
                    "protected, common hospital pneumonia, core respiratory virus, deferred Enterococcus/C. difficile, "
                    "and Acinetobacter answer-sensitive groups. Keep out of formal picked_pathogens; use review tiering by signal strength."
                ),
            )
            candidate["interpretation_caution_zh"] = (
                "此菌為 mNGS-only、無同菌醫院端支持且 D0/D1 非優勢；又不是目前保護的強肺炎病原群。"
                "不直接列為正式 picked，可依 reads/percentile 強度保留給 review/RAG 或 low-specificity audit。"
            )
    apply_same_representative_group_guardrail(candidates, ctx)


def level3_pickable(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    """Level 3 is reported only when it has a protected/supported rationale.

    Host expansion permits more Level 3 candidates, but does not report every
    weak Level 3 signal. This keeps immunocompromised cases sensitive without
    letting non-specific environmental/background organisms dominate the report.
    """
    if candidate.get("is_protected_pathogen"):
        return True
    guardrail = candidate.get("key_evidence", {}).get("guardrail_rule")
    if guardrail == "R-S4-CANDIDA-RESP-L3":
        return bool((candidate.get("key_evidence") or {}).get("candida_invasive_hospital_support"))
    if guardrail == ORAL_ASPIRATION_FLORA_L3_GUARDRAIL:
        return True
    if guardrail == "R-S4-RESP-COMMENSAL-L3":
        return True
    if related_typical_pneumonia_context_support(candidate):
        return True
    if is_candida_or_generic_yeast_name(candidate.get("organism_name")):
        return False
    support = candidate.get("module_support_summary") or {}
    if any(
        support.get(key) == "Support"
        for key in ("culture", "filmarray_gmtest", "molecular_microbiology")
    ):
        return True
    if not ctx.get("host_level3_expansion"):
        return False
    rank = mngs.to_rank_int(candidate.get("rank_priority"))
    reads_tier = candidate.get("reads_tier")
    if rank <= 2 and reads_tier in {"R4_very_high", "R3_high", "R2_medium"}:
        return True
    if candidate.get("dominance_tier") in {"D2_moderate", "D3_dominant"}:
        return True
    return False


def pickable(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    if formal_pick_excluded(candidate):
        return False
    if candidate.get("impossible_infection_source_rule"):
        return False
    if candidate.get("is_likely_colonizer_or_background") is True:
        return False
    level = level_rank(candidate.get("integrated_causative_level"))
    if level <= 2:
        return True
    if level == 3:
        return level3_pickable(candidate, ctx)
    return False


def fallback_pickable(candidate: dict[str, Any], ctx: dict[str, Any]) -> bool:
    if formal_pick_excluded(candidate):
        return False
    if candidate.get("impossible_infection_source_rule"):
        return False
    if candidate.get("is_likely_colonizer_or_background") is True:
        return False
    if (
        mngs_ranked_host_only_without_nonhost_context(candidate)
        and candidate.get("reads_tier") in LOW_SIGNAL_READ_TIERS
        and candidate.get("dominance_tier") in NON_DOMINANT_TIERS
    ):
        return False
    return True


def picked_item(candidate: dict[str, Any], role: str) -> dict[str, Any]:
    key_evidence = candidate.get("key_evidence", {}) if isinstance(candidate.get("key_evidence"), dict) else {}
    caution_flags = ["Evidence_limited", "Not_recommended_as_sole_treatment_basis"]
    item = {
        "organism_name": candidate.get("organism_name", ""),
        "classification": candidate.get("classification", "Unknown"),
        "picked_role": role,
        "evidence_source": candidate.get("evidence_source", "mNGS_ranked"),
        "basis_level": candidate.get("integrated_causative_level", "Unknown"),
        "mngs_signal_tier": candidate.get("mngs_signal_tier", "Unknown"),
        "rank_priority": candidate.get("rank_priority", "Unknown"),
        "reads": candidate.get("reads", 0),
        "support_modules": list(key_evidence.get("support_modules") or []),
        "hospital_evidence_detail": key_evidence.get("hospital_evidence_detail", {}),
        "why_picked": list(candidate.get("integrated_reasoning") or []),
    }
    if key_evidence.get("related_representative_hospital_support"):
        item["related_representative_hospital_support"] = key_evidence.get(
            "related_representative_hospital_support"
        )
        caution_flags.append("Related_group_support_not_species_specific")
    if key_evidence.get("same_representative_group_merged_picked_candidates"):
        item["same_representative_group_merged_picked_candidates"] = key_evidence.get(
            "same_representative_group_merged_picked_candidates"
        )
        caution_flags.append("Same_representative_group_evidence_consolidated")
    if is_reactivation_virus_candidate(candidate):
        caution_flags.append("Viral_reactivation_possible")
        item["interpretation_caution_zh"] = "CMV/HSV/VZV 可能是病毒再活化；需結合宿主免疫狀態、檢體、病毒量趨勢、PCR/viral load 與影像，避免單獨當成肺炎主因。"
    if is_sars_cov2_candidate(candidate):
        caution_flags.append("COVID_current_or_recent_infection_context_needed")
        item["interpretation_caution_zh"] = "COVID-19/SARS-CoV-2 需結合 PCR/抗原、症狀與影像；低 reads mNGS-only 不應單獨當主病原。"
    if is_aspergillus_mold_candidate(candidate) and candidate.get("reads_tier") in {"R0_trace", "R1_low"}:
        caution_flags.append("Low_read_mold_watchlist_preferred")
        item["interpretation_caution_zh"] = "低 reads mold 較適合作為 watchlist；除非有 GM、影像、培養或高度免疫低下支持，否則不宜當主病原。"
    if candidate.get("interpretation_caution_zh") and not item.get("interpretation_caution_zh"):
        item["interpretation_caution_zh"] = candidate.get("interpretation_caution_zh")
    item["caution_flags"] = caution_flags
    return item

def build_best_available(candidates: list[dict[str, Any]], ctx: dict[str, Any]) -> dict[str, Any]:
    sorted_candidates = sorted(candidates, key=candidate_sort_key)
    picked_candidates = [candidate for candidate in sorted_candidates if pickable(candidate, ctx)]
    if not picked_candidates:
        fallback_level = 4 if any(level_rank(c.get("integrated_causative_level")) == 4 for c in sorted_candidates) else 5
        picked_candidates = [
            candidate
            for candidate in sorted_candidates
            if level_rank(candidate.get("integrated_causative_level")) == fallback_level
            and fallback_pickable(candidate, ctx)
        ]
        selection_mode = "BestAvailable"
    else:
        selection_mode = "Confirmed" if any(level_rank(c.get("integrated_causative_level")) <= 2 for c in picked_candidates) else "BestAvailable"

    items: list[dict[str, Any]] = []
    primary_assigned = False
    for candidate in picked_candidates:
        level = level_rank(candidate.get("integrated_causative_level"))
        if level <= 2 and not primary_assigned:
            role = "Primary"
            primary_assigned = True
        elif level == 3 and candidate.get("is_protected_pathogen"):
            role = "Protected_Level3"
        elif level >= 4:
            role = "BestAvailable"
        else:
            role = "Secondary"
        items.append(picked_item(candidate, role))

    return {
        "selection_mode": selection_mode if items else "No_high_priority_candidate",
        "picked_count": len(items),
        "picked_pathogens": items,
    }


def infection_likelihood(best: dict[str, Any]) -> str:
    ranks = [level_rank(item.get("basis_level")) for item in best.get("picked_pathogens") or []]
    if 1 in ranks:
        return "Severe"
    if 2 in ranks:
        return "Likely"
    if 3 in ranks:
        return "Possible"
    if ranks:
        return "Possible"
    return "Unknown"


def dominant_pathogen_type(best: dict[str, Any]) -> str:
    classes = {
        str(item.get("classification") or "Unknown")
        for item in best.get("picked_pathogens") or []
        if str(item.get("classification") or "Unknown") != "Unknown"
    }
    if not classes:
        return "Unknown"
    return classes.pop() if len(classes) == 1 else "Mixed"


def exclusion_reason_code(candidate: dict[str, Any]) -> str:
    if candidate.get("impossible_infection_source_rule"):
        return "impossible_infection_source"
    if formal_pick_excluded(candidate):
        rule = candidate.get("formal_pick_exclusion_rule") or {}
        return str(rule.get("rule_id") or "formal_pick_guardrail")
    if candidate.get("is_likely_colonizer_or_background"):
        return "colonizer_or_background"
    return "not_picked_low_level"


def excluded_candidates(candidates: list[dict[str, Any]], best: dict[str, Any]) -> list[dict[str, Any]]:
    picked_names = {
        mngs.normalize_organism_name(item.get("organism_name"))
        for item in best.get("picked_pathogens") or []
        if isinstance(item, dict)
    }
    excluded = []
    for candidate in sorted(candidates, key=candidate_sort_key):
        if mngs.normalize_organism_name(candidate.get("organism_name")) in picked_names:
            continue
        excluded.append(
            {
                "organism_name": candidate.get("organism_name", ""),
                "classification": candidate.get("classification", "Unknown"),
                "observed_level": candidate.get("integrated_causative_level", "Unknown"),
                "rank_priority": candidate.get("rank_priority", "Unknown"),
                "reads": candidate.get("reads", 0),
                "exclusion_reason_code": exclusion_reason_code(candidate),
                "impossible_infection_source_rule": candidate.get("impossible_infection_source_rule"),
                "formal_pick_exclusion_rule": candidate.get("formal_pick_exclusion_rule"),
                "same_representative_group_retained": candidate.get("same_representative_group_retained"),
                "same_representative_group_key": candidate.get("same_representative_group_key"),
            }
        )
    return excluded


def score_payload(
    *,
    patient_dir: Path,
    ranked_mngs: dict[str, Any],
    final_summary: dict[str, Any],
) -> dict[str, Any]:
    ctx = host_context(final_summary)
    final_by_name = final_candidate_map(final_summary)
    candidates: list[dict[str, Any]] = []
    for record in ranked_mngs.get("records") or []:
        if not isinstance(record, dict):
            continue
        dominance = dominance_by_record(record)
        for raw_candidate in record.get("candidates") or []:
            if not isinstance(raw_candidate, dict):
                continue
            candidates.append(score_candidate(raw_candidate, record, dominance, final_by_name, ctx))
    candidates = mngs._dedupe_candidates(candidates)  # pylint: disable=protected-access
    hospital_added = hospital_only_candidates(final_summary, candidates, ctx)
    candidates.extend(hospital_added)
    candidates = mngs._dedupe_candidates(candidates)  # pylint: disable=protected-access
    apply_resp_commensal_case_fallback(candidates)
    apply_formal_pick_guardrails(candidates, ctx)
    best = build_best_available(candidates, ctx)
    return {
        "rule_version": RANKED_RULE_VERSION,
        "best_available_summary": best,
        "final_infection_likelihood": infection_likelihood(best),
        "dominant_source": final_summary.get("dominant_source", "Unknown"),
        "dominant_pathogen_type": dominant_pathogen_type(best),
        "host_context": ctx,
        "pathogen_candidates": candidates,
        "excluded_candidates": excluded_candidates(candidates, best),
        "data_gaps": [],
        "cross_module_reasoning": [
            "deterministic scorer：候選池取自 ranked_mngs，並補入 Level 1/2 的直接醫院端"
            f" hospital_only 候選；本次共 {len(candidates)} 個去重候選，其中 hospital_only={len(hospital_added)}。",
            "Level / picked / excluded 由 Python 規則計算，未使用 LLM 判斷。",
        ],
        "source_files": {
            "patient_dir": str(patient_dir),
        },
    }


def output_path_for(patient_dir: Path, output_suffix: str) -> Path:
    base = mngs.extract_patient_identifier(patient_dir)
    return patient_dir / mngs.SUMMARY_OUTPUT_DIR_NAME / f"{base}_{output_suffix}.json"


def process_directory(
    patient_dir: Path,
    *,
    ranked_mngs_path: Path | None,
    summary_mode: str,
    summary_suffix: str | None,
    output_suffix: str,
    overwrite: bool,
) -> Path | None:
    if summary_suffix:
        patient_id = mngs.extract_patient_identifier(patient_dir)
        final_summary_path = (
            patient_dir
            / mngs.SUMMARY_OUTPUT_DIR_NAME
            / f"{patient_id}_{summary_suffix}.json"
        )
        if not final_summary_path.exists():
            final_summary_path = None
    else:
        final_summary_path = mngs.find_final_summary_file(patient_dir, summary_mode=summary_mode)
    if final_summary_path is None:
        print(f"skip {patient_dir.name}: final_summary missing")
        return None
    ranked_mngs = mngs.load_ranked_mngs_for_patient(ranked_mngs_path, patient_dir)
    if ranked_mngs is None:
        print(f"skip {patient_dir.name}: ranked_mngs missing")
        return None
    destination = output_path_for(patient_dir, output_suffix)
    if destination.exists() and not overwrite:
        print(f"skip {patient_dir.name}: output exists {destination}")
        return None
    final_summary = load_json(final_summary_path)
    payload = score_payload(patient_dir=patient_dir, ranked_mngs=ranked_mngs, final_summary=final_summary)
    payload["source_files"] = {
        "final_summary": str(final_summary_path),
        "ranked_mngs": str(ranked_mngs_path or "local/default ranked lookup"),
    }
    write_json(destination, payload)
    return destination


def parse_patient_ids(values: Sequence[str] | None) -> set[str]:
    output: set[str] = set()
    for value in values or []:
        text = str(value).strip()
        if not text:
            continue
        if text.startswith("NGS_patient_"):
            output.add(text.removesuffix("_json"))
        else:
            output.add(f"NGS_patient_{text}")
    return output


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run deterministic mNGS max scorer without LLM.")
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--patients", nargs="*")
    parser.add_argument("--ranked-mngs", type=Path, default=mngs.DEFAULT_RANKED_MNGS_PATH)
    parser.add_argument(
        "--summary-mode",
        choices=["deterministic", "full"],
        default="deterministic",
        help=(
            "Which final_summary variant to load. Default deterministic uses "
            "*_final_summary_with_filmarray_deterministic.json."
        ),
    )
    parser.add_argument(
        "--summary-suffix",
        help=(
            "Explicit final-summary suffix to load from each patient's summary_outputs directory. "
            "When provided, this takes precedence over --summary-mode."
        ),
    )
    parser.add_argument("--output-suffix", default=DEFAULT_MODELLESS_OUTPUT_SUFFIX)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    patient_dirs = mngs.collect_patient_directories(
        args.inputs,
        summary_mode=args.summary_mode,
        mngs_source="ranked_only",
    )
    requested = parse_patient_ids(args.patients)
    if requested:
        patient_dirs = [
            path for path in patient_dirs if mngs.extract_patient_identifier(path) in requested
        ]
    written = []
    for patient_dir in patient_dirs:
        output = process_directory(
            patient_dir,
            ranked_mngs_path=args.ranked_mngs,
            summary_mode=args.summary_mode,
            summary_suffix=args.summary_suffix,
            output_suffix=args.output_suffix,
            overwrite=args.overwrite,
        )
        if output is not None:
            written.append(str(output))
            print(f"wrote {output}")
    print(f"written_count={len(written)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

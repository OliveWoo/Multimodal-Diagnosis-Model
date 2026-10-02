"""Merge deterministic max output and LLM missed-candidate review without dropping fields."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import candidate_evidence_profile as evidence_profiles  # noqa: E402
from tools import organism_taxonomy_classifier as organism_taxonomy  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402

DEFAULT_MAX_SUFFIX = "mNGS_max_deterministic_resp_commensal_dominance_guardrail_opt_chosen_full"
DEFAULT_REVIEW_SUFFIX = "mNGS_missed_candidate_review"
DEFAULT_OUTPUT_SUFFIX = "mNGS_max_deterministic_with_missed_candidate_review"


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8-sig")


def ensure_taxonomy_profile(item: Any) -> None:
    """Attach routing metadata without changing tier or causality decisions."""
    if not isinstance(item, dict):
        return
    name = item.get("organism_name") or item.get("name")
    if not name:
        return
    profile = item.get("taxonomy_profile")
    if isinstance(profile, dict) and profile.get("profile_version") == organism_taxonomy.PROFILE_VERSION:
        return
    snapshot = item.get("evidence_snapshot") if isinstance(item.get("evidence_snapshot"), dict) else {}
    classification = item.get("classification") or snapshot.get("classification")
    item["taxonomy_profile"] = organism_taxonomy.classify_organism(
        name,
        biological_class=classification,
    )


def attach_taxonomy_profiles(deterministic_max: dict[str, Any], review: dict[str, Any]) -> None:
    """Backfill taxonomy profiles into both new and previously generated inputs."""
    for section in ("pathogen_candidates", "excluded_candidates"):
        for item in deterministic_max.get(section) or []:
            ensure_taxonomy_profile(item)
    best = deterministic_max.get("best_available_summary")
    if isinstance(best, dict):
        for section in (
            "picked_pathogens",
            "excluded_candidates",
            "best_available_pathogens",
            "representative_omitted_candidates",
        ):
            for item in best.get(section) or []:
                ensure_taxonomy_profile(item)
    for section in (
        "review_high_priority",
        "review_context_needed",
        "review_low_specificity",
        "review_omitted_with_reason",
    ):
        for item in review.get(section) or []:
            ensure_taxonomy_profile(item)


def patient_id(patient_dir: Path) -> str:
    name = patient_dir.name
    return name.removeprefix("NGS_patient_").removesuffix("_json")


def sort_key(patient_dir: Path) -> tuple[int, str]:
    value = patient_id(patient_dir)
    try:
        return int(value), value
    except ValueError:
        return 10**9, value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge deterministic max JSON and LLM missed-candidate review JSON per patient."
    )
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--patients", nargs="*")
    parser.add_argument("--max-suffix", default=DEFAULT_MAX_SUFFIX)
    parser.add_argument("--review-suffix", default=DEFAULT_REVIEW_SUFFIX)
    parser.add_argument("--output-suffix", default=DEFAULT_OUTPUT_SUFFIX)
    parser.add_argument("--summary-csv", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def normalize_requested(values: list[str] | None) -> set[str]:
    result: set[str] = set()
    for value in values or []:
        text = str(value).strip()
        if text.startswith("NGS_patient_"):
            text = text.removeprefix("NGS_patient_").removesuffix("_json")
        if text:
            result.add(text)
    return result


def normalize_name(value: Any) -> str:
    return pathogen_names.canonical_key(value)


def display_name(value: Any) -> str:
    return pathogen_names.clean_display_text(value)


def pathogen_alias_info(value: Any) -> tuple[str, str | None, str | None]:
    """Return canonical key, relationship type, and clinically useful note."""
    return pathogen_names.alias_info(value)


MERGEABLE_CLOSE_GENERA = {
    "elizabethkingia",
    "chryseobacterium",
    "sphingomonas",
    "ralstonia",
    "stenotrophomonas",
    "aspergillus",
}

LOW_SPECIFICITY_ENVIRONMENTAL_PREFIXES = (
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
    "sphingobium",
    "sphingomonas",
    "undibacterium",
    "acidovorax",
    "brevundimonas",
    "cloacibacterium",
    "comamonas",
    "flavobacterium",
    "ralstonia",
    "pseudoxanthomonas",
)

ACTIONABLE_RESPIRATORY_REVIEW_KEYS = {
    "haemophilusinfluenzae",
    "enterobactercloacaecomplex",
}

DIRECT_HOSPITAL_MODULES = set(evidence_profiles.DIRECT_HOSPITAL_MODULES)

NONPULMONARY_SOURCE_TERMS = (
    "foley",
    "urine",
    "urinary",
    "stool",
    "feces",
    "faeces",
    "rectal",
)

PULMONARY_SOURCE_TERMS = (
    "bal",
    "bronchoalveolar",
    "bronchial",
    "sputum",
    "tracheal",
    "endotracheal",
    "lower respiratory",
    "respiratory",
    "pleural",
    "lung",
)

CONTEXT_REVIEW_KEYS = {
    "acinetobacterbaumannii",
    "cronobactersakazakii",
    "enterobactersoli",
    "escherichiacoli",
    "mycobacteriumtuberculosis",
    "nocardiathailandica",
    "rhodococcusqingshengii",
    "schizophyllumcommune",
    "stenotrophomonasmaltophilia",
    "streptococcuspneumoniae",
    "talaromycespinophilus",
}

ASPIRATION_CONTEXT_KEYS = {
    "bacteroidesfragilis",
    "fusobacteriumnucleatum",
    "porphyromonasendodontalis",
}

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

ENVIRONMENTAL_WATER_SOIL_GNB_PREFIXES = LOW_SPECIFICITY_ENVIRONMENTAL_PREFIXES + (
    "achromobacter",
    "chryseobacterium",
    "delftia",
    "elizabethkingia",
    "paraburkholderia",
    "polynucleobacter",
    "phytobacter",
    "xanthomonas",
)

SKIN_AIRWAY_BACKGROUND_PREFIXES = (
    "cutibacterium",
    "lactiplantibacillus",
    "lacticaseibacillus",
    "ligilactobacillus",
    "limosilactobacillus",
    "lactobacillus",
)

COAG_NEG_STAPH_BACKGROUND_EXACT_NAMES = {
    "staphylococcuscapitis",
    "staphylococcusepidermidis",
    "staphylococcushaemolyticus",
    "staphylococcushominis",
    "staphylococcuswarneri",
    "staphylococcuscohnii",
}

RARE_OPPORTUNISTIC_YEAST_PREFIXES = ("trichosporon",)

CORE_ANSWER_SENSITIVE_RESPIRATORY_VIRUS_KEYS = {
    "severeacuterespiratorysyndromerelatedcoronavirus",
}

CORE_ANSWER_SENSITIVE_RESPIRATORY_VIRUS_PREFIXES = (
    "influenza",
    "sarscov2",
    "severeacuterespiratorysyndrome",
)

GENERIC_WEAK_REVIEW_DEMOTION_GROUPS = {
    "strict_aspiration_anaerobe",
    "broad_oral_upper_airway_commensal",
    "atypical_oral_associated",
    "other_oral_aspiration_flora",
    "environmental_water_soil_low_specificity_gnb",
    "skin_or_airway_background",
    "other_virus",
    "other_fungus",
    "other_bacterial_or_unclassified",
    "other_or_unclassified",
}

NON_CORE_RESPIRATORY_VIRUS_FORMAL_PICK_GUARDRAIL = "R-S4-NON-CORE-RESP-VIRUS-PICKED-OMIT"
MNGS_ONLY_NONDOMINANT_NON_TYPICAL_FORMAL_PICK_GUARDRAIL = "R-S4-MNGS-ONLY-NONDOMINANT-NON-TYPICAL-PICKED-OMIT"
COMMON_HOSPITAL_LOW_READ_NODOM_FORMAL_PICK_GUARDRAIL = "R-S4-COMMON-HOSPITAL-R0R1-NONDOMINANT-PICKED-OMIT"
CANDIDA_EVIDENCE_STRONG_INVASIVE = evidence_profiles.CANDIDA_EVIDENCE_STRONG_INVASIVE
CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC = evidence_profiles.CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC
CANDIDA_EVIDENCE_WEAK_NONINVASIVE = evidence_profiles.CANDIDA_EVIDENCE_WEAK_NONINVASIVE
REVIEW_TIER_CONVERGENCE_WEAK_CANDIDA_RULE = "RTCV2-WEAK-CANDIDA-NONINVASIVE"
REVIEW_TIER_CONVERGENCE_NONSTRIATUM_CORYNE_RULE = "RTCV2-NONSTRIATUM-CORYNE-WEAK"


def genus_of(value: Any) -> str:
    return pathogen_names.genus_name(value)


def starts_with_any(value: Any, prefixes: tuple[str, ...]) -> bool:
    name = normalize_name(value)
    return any(name.startswith(prefix) for prefix in prefixes)


def is_broad_or_atypical_oral_flora_name(value: Any) -> bool:
    name = normalize_name(value)
    return (
        name in BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_EXACT_NAMES
        or name in ATYPICAL_ORAL_ASSOCIATED_EXACT_NAMES
        or starts_with_any(value, BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_PREFIXES)
        or starts_with_any(value, ATYPICAL_ORAL_ASSOCIATED_PREFIXES)
    )


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def first_present(*values: Any) -> Any:
    for value in values:
        if value not in (None, "", [], {}):
            return value
    return ""


def is_strict_aspiration_anaerobe_name(value: Any) -> bool:
    return starts_with_any(value, STRICT_ASPIRATION_ANAEROBE_PREFIXES)


def oral_aspiration_flora_category(value: Any) -> str:
    name = normalize_name(value)
    if is_strict_aspiration_anaerobe_name(value):
        return "strict_aspiration_anaerobe"
    if name in BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_EXACT_NAMES or starts_with_any(
        value, BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_PREFIXES
    ):
        return "broad_oral_upper_airway_commensal"
    if name in ATYPICAL_ORAL_ASSOCIATED_EXACT_NAMES or starts_with_any(value, ATYPICAL_ORAL_ASSOCIATED_PREFIXES):
        return "atypical_oral_associated"
    if is_broad_or_atypical_oral_flora_name(value):
        return "other_oral_aspiration_flora"
    return "not_oral_aspiration_flora"


def is_rare_opportunistic_yeast_name(value: Any) -> bool:
    return starts_with_any(value, RARE_OPPORTUNISTIC_YEAST_PREFIXES)


def is_environmental_water_soil_gnb_name(value: Any) -> bool:
    return starts_with_any(value, ENVIRONMENTAL_WATER_SOIL_GNB_PREFIXES)


def is_skin_airway_background_name(value: Any) -> bool:
    name = normalize_name(value)
    return name in COAG_NEG_STAPH_BACKGROUND_EXACT_NAMES or starts_with_any(value, SKIN_AIRWAY_BACKGROUND_PREFIXES)


def is_core_answer_sensitive_respiratory_virus_name(value: Any) -> bool:
    name = normalize_name(value)
    key, _, _ = pathogen_alias_info(value)
    return (
        key in CORE_ANSWER_SENSITIVE_RESPIRATORY_VIRUS_KEYS
        or name in CORE_ANSWER_SENSITIVE_RESPIRATORY_VIRUS_KEYS
        or any(name.startswith(prefix) for prefix in CORE_ANSWER_SENSITIVE_RESPIRATORY_VIRUS_PREFIXES)
    )


def evidence_from_item(item: dict[str, Any]) -> dict[str, Any]:
    evidence = item.get("evidence_snapshot")
    return evidence if isinstance(evidence, dict) else item


def nested_text(value: Any) -> str:
    """Flatten evidence values for conservative specimen-source classification."""
    if isinstance(value, dict):
        return " ".join(nested_text(item) for item in value.values())
    if isinstance(value, list):
        return " ".join(nested_text(item) for item in value)
    return str(value or "").lower()


def hospital_source_guardrail_profile(evidence: dict[str, Any]) -> dict[str, Any]:
    """Identify hospital support that is exclusively urine/Foley/stool sourced.

    Unknown-source evidence is never classified as non-pulmonary-only. This
    keeps the guardrail conservative when source metadata are incomplete.
    """
    hospital_items = as_list(evidence.get("hospital_evidence"))
    direct_rows: list[dict[str, Any]] = []
    for hospital_item in hospital_items:
        if not isinstance(hospital_item, dict):
            continue
        module_evidence = as_dict(hospital_item.get("module_evidence"))
        for module in DIRECT_HOSPITAL_MODULES:
            for row in as_list(module_evidence.get(module)):
                text = nested_text(row)
                source_class = "unknown"
                if any(term in text for term in PULMONARY_SOURCE_TERMS):
                    source_class = "pulmonary_or_pleural"
                elif any(term in text for term in NONPULMONARY_SOURCE_TERMS):
                    source_class = "nonpulmonary_urine_foley_stool"
                direct_rows.append(
                    {
                        "module": module,
                        "source_class": source_class,
                        "source_text": text[:500],
                    }
                )
    classes = {row["source_class"] for row in direct_rows}
    return {
        "evaluated": bool(direct_rows),
        "nonpulmonary_only": bool(direct_rows)
        and classes == {"nonpulmonary_urine_foley_stool"},
        "source_classes": sorted(classes),
        "direct_evidence_rows": direct_rows,
    }


def strip_nonpulmonary_support_modules(evidence: dict[str, Any]) -> None:
    """Stop non-pulmonary hospital findings from acting as pneumonia support."""
    for field in ("non_host_support_modules", "support_modules"):
        values = evidence.get(field)
        if not isinstance(values, list):
            continue
        evidence.setdefault(f"original_{field}", list(values))
        evidence[field] = [
            value
            for value in values
            if str(value).lower() not in DIRECT_HOSPITAL_MODULES
        ]
    related = evidence.get("related_representative_hospital_support")
    if isinstance(related, list) and related:
        evidence.setdefault("original_related_representative_hospital_support", copy.deepcopy(related))
        evidence["related_representative_hospital_support"] = []


def evidence_rule_ids(evidence: dict[str, Any]) -> set[str]:
    rule_ids: set[str] = set()
    for field in ("applied_rules", "review_reasons"):
        values = evidence.get(field)
        if isinstance(values, list):
            rule_ids.update(str(value) for value in values if value)
    for field in ("guardrail_rule",):
        if evidence.get(field):
            rule_ids.add(str(evidence.get(field)))
    formal_rule = evidence.get("formal_pick_exclusion_rule")
    if isinstance(formal_rule, dict) and formal_rule.get("rule_id"):
        rule_ids.add(str(formal_rule.get("rule_id")))
    key_evidence = evidence.get("key_evidence")
    if isinstance(key_evidence, dict) and key_evidence.get("guardrail_rule"):
        rule_ids.add(str(key_evidence.get("guardrail_rule")))
    return rule_ids


def has_related_representative_context(evidence: dict[str, Any]) -> bool:
    related = evidence.get("related_representative_hospital_support")
    if isinstance(related, list) and related:
        return True
    key_evidence = evidence.get("key_evidence")
    if isinstance(key_evidence, dict):
        related = key_evidence.get("related_representative_hospital_support")
        return isinstance(related, list) and bool(related)
    return False


def has_non_host_support(evidence: dict[str, Any]) -> bool:
    support = evidence.get("non_host_support_modules")
    if isinstance(support, list) and support:
        return True
    support = evidence.get("support_modules")
    if isinstance(support, list):
        return any(str(item) != "host" for item in support)
    return False


def candida_best_evidence_strength(evidence: dict[str, Any]) -> str:
    return evidence_profiles.candida_best_evidence_strength(evidence)


def candida_has_invasive_context(evidence: dict[str, Any]) -> bool:
    return evidence_profiles.candida_has_invasive_context(evidence)


def candida_has_intermediate_context(evidence: dict[str, Any]) -> bool:
    return evidence_profiles.candida_has_intermediate_context(evidence)


def numeric(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def tier_score(value: Any) -> float:
    text = str(value or "")
    if "R4" in text or "M1" in text:
        return 40
    if "R3" in text or "M2" in text:
        return 30
    if "R2" in text or "M3" in text:
        return 20
    if "R1" in text or "M4" in text:
        return 10
    return 0


def level_score(value: Any) -> float:
    text = str(value or "")
    match = re.search(r"Level\\s*(\\d+)", text, flags=re.I)
    if not match:
        return 0
    return max(0, 60 - int(match.group(1)) * 10)


def entry_score(entry: dict[str, Any]) -> float:
    section_bonus = {
        "picked_pathogens": 10000,
        "possible_missed_pathogens": 5000,
        "review_high_priority": 5000,
        "review_context_needed": 3000,
        "watchlist_candidates": 1000,
        "review_low_specificity": 500,
        "review_low_specificity_audit": 500,
        "review_omitted_with_reason": 100,
    }.get(entry.get("section"), 0)
    item = entry.get("item") or {}
    evidence = evidence_from_item(item)
    rank = numeric(evidence.get("rank_priority"), default=99)
    reads = numeric(evidence.get("reads"))
    percentile = numeric(evidence.get("reads_percentile"))
    support = evidence.get("non_host_support_modules") or evidence.get("support_modules") or []
    if not isinstance(support, list):
        support = []
    return (
        section_bonus
        + pathogen_names.representative_specificity_score(item.get("organism_name") or item.get("name"))
        + (500 if reads > 0 else 0)
        + level_score(item.get("basis_level") or evidence.get("integrated_causative_level"))
        + tier_score(evidence.get("mngs_signal_tier"))
        + tier_score(evidence.get("reads_tier"))
        + min(reads, 100000) / 1000
        + percentile * 100
        + max(0, 80 - rank * 10)
        + len(support) * 10
    )


def is_weak_or_related_duplicate(entry: dict[str, Any]) -> bool:
    item = entry.get("item") or {}
    evidence = evidence_from_item(item)
    if entry.get("section") != "picked_pathogens":
        return True
    if item.get("evidence_source") == "hospital_only":
        return True
    if numeric(evidence.get("reads")) <= 0 and evidence.get("rank_priority") in {None, "", "Not_available"}:
        return True
    return False


def relationship_for_duplicate(loser_name: str, winner_name: str, loser: dict[str, Any]) -> tuple[str, str | None]:
    loser_key, rel, note = pathogen_alias_info(loser_name)
    winner_key, _, _ = pathogen_alias_info(winner_name)
    if loser_key and loser_key == winner_key and rel:
        return rel, note
    if pathogen_names.approved_group_member_match(loser_name, winner_name):
        return (
            "approved_group_member_support",
            "group/complex-level evidence supports the retained member species but does not confirm exact speciation",
        )
    loser_genus = genus_of(loser_name)
    winner_genus = genus_of(winner_name)
    if loser_genus and loser_genus == winner_genus:
        return "same_genus_related_evidence", f"same genus as retained representative: {winner_name}"
    return "related_evidence", None


def add_related_evidence(winner_item: dict[str, Any], loser_entry: dict[str, Any]) -> None:
    loser_item = loser_entry.get("item") or {}
    loser_name = display_name(
        loser_item.get("original_organism_name") or loser_item.get("organism_name") or loser_item.get("name")
    )
    winner_name = display_name(winner_item.get("organism_name") or winner_item.get("name"))
    rel, note = relationship_for_duplicate(loser_name, winner_name, loser_entry)
    related = winner_item.setdefault("related_evidence", [])
    if not isinstance(related, list):
        related = []
        winner_item["related_evidence"] = related
    payload = {
        "original_name": loser_name,
        "canonical_key": pathogen_names.canonical_key(loser_name),
        "source_section": loser_entry.get("section"),
        "relationship": rel,
    }
    if note:
        payload["note"] = note
    evidence = evidence_from_item(loser_item)
    snap = evidence_snapshot(evidence) if evidence else {}
    if snap:
        payload["evidence_snapshot"] = snap
    if loser_item.get("rationale_zh"):
        payload["rationale_zh"] = loser_item.get("rationale_zh")
    if loser_item.get("why_picked"):
        payload["why_picked"] = loser_item.get("why_picked")
    keys = {(entry.get("original_name"), entry.get("source_section"), entry.get("relationship")) for entry in related if isinstance(entry, dict)}
    key = (payload.get("original_name"), payload.get("source_section"), payload.get("relationship"))
    if key not in keys:
        related.append(payload)
    _, rel_type, alias_note = pathogen_alias_info(loser_name)
    if rel_type == "resistance_phenotype" and alias_note:
        winner_item["resistance_or_phenotype_note"] = alias_note
        winner_item["display_name"] = f"{winner_name} ({alias_note})"
    elif alias_note and not winner_item.get("display_name"):
        winner_item["display_name"] = alias_note


def is_candida_or_generic_yeast_name(value: Any) -> bool:
    return pathogen_names.is_candida_or_generic_yeast(value)


def is_reactivation_virus_name(value: Any) -> bool:
    key, _, _ = pathogen_alias_info(value)
    return key in {
        "humanbetaherpesvirus5",
        "humanalphaherpesvirus1",
        "humanalphaherpesvirus2",
        "humanalphaherpesvirus3",
        "humangammaherpesvirus4",
    }


ORGANISM_NAME_LABEL_PREFIX_TERMS = (
    "代表",
    "菌群",
    "相關菌",
    "oral",
    "aspiration",
    "flora",
    "representative",
    "pattern",
)


def clean_review_organism_name(value: Any) -> str:
    text = pathogen_names.clean_display_text(value)
    if not text:
        return ""
    text = re.sub(r"^\s*[-*•\d.)、\s]+", "", text).strip()
    for delimiter in ("：", ":"):
        if delimiter not in text:
            continue
        prefix, suffix = text.rsplit(delimiter, 1)
        if suffix.strip() and any(term in prefix.lower() for term in ORGANISM_NAME_LABEL_PREFIX_TERMS):
            text = suffix.strip()
            break
    return pathogen_names.clean_display_text(text)


def canonicalize_item_name(item: dict[str, Any]) -> None:
    name = clean_review_organism_name(item.get("organism_name") or item.get("name"))
    if name:
        item["organism_name"] = name
    key, rel_type, alias_note = pathogen_alias_info(name)
    if not key:
        return
    current = display_name(name)
    canonical = pathogen_names.display_name(name)
    item["canonical_key"] = key
    if current and canonical and current != canonical and not item.get("original_organism_name"):
        item["original_organism_name"] = current
    if canonical:
        item["organism_name"] = canonical
    if rel_type == "resistance_phenotype" and alias_note:
        item["resistance_or_phenotype_note"] = alias_note
        item["display_name"] = f"{canonical or current} ({alias_note})"
    elif alias_note:
        item.setdefault("alias_note", alias_note)


def add_interpretation_caution(item: dict[str, Any]) -> None:
    name = item.get("organism_name") or item.get("name")
    if is_reactivation_virus_name(name):
        item.setdefault(
            "interpretation_caution_zh",
            "CMV/HSV/EBV 可能是病毒再活化；需結合宿主免疫狀態、檢體、病毒量趨勢與影像，避免單獨當成肺炎主因。",
        )
    key, _, _ = pathogen_alias_info(name)
    if key == "severeacuterespiratorysyndromerelatedcoronavirus":
        item.setdefault(
            "interpretation_caution_zh",
            "COVID-19/SARS-CoV-2 需結合 PCR/抗原、症狀與影像；低 reads mNGS-only 不應單獨當主病原。",
        )


def grouping_key(entry: dict[str, Any]) -> str:
    item = entry.get("item") or {}
    name = item.get("organism_name") or item.get("name")
    return pathogen_names.representative_group_key(name)


def should_merge_group_loser(loser: dict[str, Any], winner: dict[str, Any]) -> bool:
    loser_item = loser.get("item") or {}
    winner_item = winner.get("item") or {}
    loser_name = loser_item.get("organism_name") or loser_item.get("name")
    winner_name = winner_item.get("organism_name") or winner_item.get("name")
    loser_key = pathogen_names.canonical_key(loser_name)
    winner_key = pathogen_names.canonical_key(winner_name)
    if loser_key and winner_key and loser_key == winner_key:
        return True
    return pathogen_names.representative_group_suppression_allowed(loser_name)


def reconcile_sections(deterministic_max: dict[str, Any], review: dict[str, Any]) -> dict[str, Any]:
    best = deterministic_max.get("best_available_summary")
    if not isinstance(best, dict):
        return {"enabled": True, "merged_count": 0, "merged_items": []}
    sections = {
        "picked_pathogens": best.get("picked_pathogens") or [],
        "possible_missed_pathogens": review.get("possible_missed_pathogens") or [],
        "watchlist_candidates": review.get("watchlist_candidates") or [],
        "review_high_priority": review.get("review_high_priority") or [],
        "review_context_needed": review.get("review_context_needed") or [],
        "review_low_specificity": review.get("review_low_specificity") or [],
        "review_low_specificity_audit": review.get("review_low_specificity_audit") or [],
        "review_omitted_with_reason": review.get("review_omitted_with_reason") or [],
    }
    entries: list[dict[str, Any]] = []
    for section, values in sections.items():
        if not isinstance(values, list):
            continue
        for index, item in enumerate(values):
            if isinstance(item, dict):
                canonicalize_item_name(item)
                add_interpretation_caution(item)
                entries.append({"section": section, "index": index, "item": item, "score": entry_score({"section": section, "item": item})})

    groups: dict[str, list[dict[str, Any]]] = {}
    for entry in entries:
        groups.setdefault(grouping_key(entry), []).append(entry)

    remove_ids: set[tuple[str, int]] = set()
    merged_items: list[dict[str, Any]] = []
    for key, group in groups.items():
        if len(group) <= 1:
            continue
        if key.startswith("genus:") and not any(is_weak_or_related_duplicate(entry) for entry in group):
            continue
        winner = max(group, key=lambda entry: entry["score"])
        winner_item = winner["item"]
        for loser in group:
            if loser is winner:
                continue
            if not should_merge_group_loser(loser, winner):
                continue
            remove_ids.add((loser["section"], loser["index"]))
            add_related_evidence(winner_item, loser)
            loser_item = loser["item"] or {}
            merged_items.append(
                {
                    "retained": winner_item.get("organism_name") or winner_item.get("name"),
                    "merged": (
                        loser_item.get("original_organism_name")
                        or loser_item.get("organism_name")
                        or loser_item.get("name")
                    ),
                    "from_section": loser["section"],
                    "into_section": winner["section"],
                    "group_key": key,
                }
            )

    for section, values in sections.items():
        if not isinstance(values, list):
            continue
        filtered = [
            item
            for index, item in enumerate(values)
            if (section, index) not in remove_ids
        ]
        if section == "picked_pathogens":
            best["picked_pathogens"] = filtered
            best["picked_count"] = len(filtered)
        else:
            review[section] = filtered

    return {
        "enabled": True,
        "policy": (
            "Merge-level output collapses same pathogen aliases, resistance phenotypes, disease/test labels, "
            "and selected close-genus duplicate signals into the strongest representative while preserving related evidence."
        ),
        "merged_count": len(merged_items),
        "merged_items": merged_items,
    }


GROUP_REPRESENTATIVE_REVIEW_SECTIONS = (
    "review_high_priority",
    "review_context_needed",
    "review_low_specificity",
    "review_low_specificity_audit",
    "review_omitted_with_reason",
)

GROUP_REPRESENTATIVE_TIER_PRIORITY = {
    "review_high_priority": 4,
    "review_context_needed": 3,
    "review_low_specificity": 2,
    "review_low_specificity_audit": 2,
    "review_omitted_with_reason": 0,
}


def group_member_representative_score(entry: dict[str, Any]) -> tuple[float, ...]:
    item = entry.get("item") or {}
    evidence = evidence_from_item(item)
    return (
        1.0 if numeric(evidence.get("reads")) > 0 else 0.0,
        tier_score(evidence.get("reads_tier")),
        tier_score(evidence.get("mngs_signal_tier")),
        numeric(evidence.get("reads_percentile")),
        numeric(evidence.get("reads")),
        -numeric(evidence.get("rank_priority"), default=99),
    )


def apply_species_preferred_group_representatives(review: dict[str, Any]) -> dict[str, Any]:
    """Replace a broad review group label with an explicit detected member.

    The broad hospital group evidence keeps its original tier and is attached
    to the selected mNGS member as group-level support. Other same-genus mNGS
    species remain visible only as related evidence. This changes display and
    evidence organization, not picked status or clinical tier.
    """
    entries: list[dict[str, Any]] = []
    for section in GROUP_REPRESENTATIVE_REVIEW_SECTIONS:
        values = review.get(section)
        if not isinstance(values, list):
            continue
        for index, item in enumerate(values):
            if isinstance(item, dict):
                canonicalize_item_name(item)
                entries.append({"section": section, "index": index, "item": item})

    replacements: dict[tuple[str, int], dict[str, Any]] = {}
    removals: set[tuple[str, int]] = set()
    transformations: list[dict[str, Any]] = []
    claimed_members: set[tuple[str, int]] = set()

    for group_entry in entries:
        group_item = group_entry["item"]
        group_name = group_item.get("organism_name") or group_item.get("name")
        approved_members = pathogen_names.approved_group_members(group_name)
        if not approved_members:
            continue

        group_priority = GROUP_REPRESENTATIVE_TIER_PRIORITY.get(group_entry["section"], 0)

        member_entries = [
            entry
            for entry in entries
            if (entry["section"], entry["index"]) not in claimed_members
            and GROUP_REPRESENTATIVE_TIER_PRIORITY.get(entry["section"], 0) <= group_priority
            and pathogen_names.canonical_key(
                entry["item"].get("organism_name") or entry["item"].get("name")
            )
            in approved_members
            and numeric(evidence_from_item(entry["item"]).get("reads")) > 0
        ]
        if not member_entries:
            continue

        winner = max(member_entries, key=group_member_representative_score)
        winner_item = copy.deepcopy(winner["item"])
        winner_name = winner_item.get("organism_name") or winner_item.get("name")
        winner_genus = genus_of(winner_name)
        group_section = group_entry["section"]

        winner_item["review_tier"] = group_section
        winner_item["review_tier_reason"] = (
            "explicit mNGS member species selected as the representative for group-level hospital evidence; "
            "the group result supports the candidate but does not confirm exact species"
        )
        winner_item["original_review_section"] = group_item.get("original_review_section") or group_section
        winner_item["representative_for_group"] = display_name(group_name)
        winner_item["group_level_support_only"] = True
        winner_item["species_confidence"] = "low_to_moderate"
        winner_item["group_evidence_snapshot"] = copy.deepcopy(evidence_from_item(group_item))
        winner_item["group_level_rationale_zh"] = group_item.get("rationale_zh")
        add_related_evidence(winner_item, group_entry)

        related_entries: list[dict[str, Any]] = []
        for entry in entries:
            entry_key = (entry["section"], entry["index"])
            if entry is group_entry or entry is winner or entry_key in claimed_members:
                continue
            entry_name = entry["item"].get("organism_name") or entry["item"].get("name")
            if (
                winner_genus
                and genus_of(entry_name) == winner_genus
                and GROUP_REPRESENTATIVE_TIER_PRIORITY.get(entry["section"], 0) <= group_priority
            ):
                related_entries.append(entry)
                add_related_evidence(winner_item, entry)

        group_key = (group_entry["section"], group_entry["index"])
        winner_key = (winner["section"], winner["index"])
        replacements[group_key] = winner_item
        removals.add(winner_key)
        claimed_members.add(winner_key)
        for entry in related_entries:
            related_key = (entry["section"], entry["index"])
            removals.add(related_key)
            claimed_members.add(related_key)

        transformations.append(
            {
                "group_label": display_name(group_name),
                "retained_representative": display_name(winner_name),
                "retained_tier": group_section,
                "member_reads": numeric(evidence_from_item(winner["item"]).get("reads")),
                "related_species": [
                    display_name(entry["item"].get("organism_name") or entry["item"].get("name"))
                    for entry in related_entries
                ],
                "species_confidence": "low_to_moderate",
            }
        )

    if replacements or removals:
        for section in GROUP_REPRESENTATIVE_REVIEW_SECTIONS:
            values = review.get(section)
            if not isinstance(values, list):
                continue
            rebuilt: list[dict[str, Any]] = []
            for index, item in enumerate(values):
                key = (section, index)
                if key in replacements:
                    rebuilt.append(replacements[key])
                elif key not in removals:
                    rebuilt.append(item)
            review[section] = rebuilt

    review["species_preferred_group_representative_policy"] = {
        "enabled": True,
        "placement": "after review-tier convergence and before final name reconciliation",
        "tier_policy": "preserve the broad group item's review tier; do not promote to picked",
        "species_policy": "retain an explicitly approved detected member as the visible representative; keep group and same-genus alternatives as related evidence",
        "transformation_count": len(transformations),
        "transformations": transformations,
    }
    return review


def weak_environmental_review_item(item: dict[str, Any], *, has_picked: bool) -> bool:
    name = item.get("organism_name") or item.get("name")
    if not starts_with_any(name, LOW_SPECIFICITY_ENVIRONMENTAL_PREFIXES):
        return False
    evidence = evidence_from_item(item)
    if evidence.get("source_category") != "hospital_detected_only":
        return False
    if numeric(evidence.get("reads")) > 0:
        return False
    haystack = json.dumps(evidence, ensure_ascii=False).lower()
    if any(term in haystack for term in ("blood", "csf", "tissue", "abscess", "pus", "pleural", "sterile")):
        return False
    if not any(term in haystack for term in ("sputum", "respiratory", "tracheal", "eta")):
        return False
    if any(term in haystack for term in ("heavy", "many", "moderate", "3+", "4+", "abundant", "predominant", "pure", "repeat")):
        return False
    return has_picked


def filter_weak_environmental_review_items(review: dict[str, Any], deterministic_max: dict[str, Any]) -> dict[str, Any]:
    best = deterministic_max.get("best_available_summary") if isinstance(deterministic_max, dict) else {}
    has_picked = bool(isinstance(best, dict) and (best.get("picked_pathogens") or []))
    omitted: list[dict[str, Any]] = []
    for section in ("possible_missed_pathogens", "watchlist_candidates"):
        values = review.get(section) or []
        if not isinstance(values, list):
            review[section] = []
            continue
        kept = []
        for item in values:
            if isinstance(item, dict) and weak_environmental_review_item(item, has_picked=has_picked):
                omitted.append(
                    {
                        "organism_name": item.get("organism_name"),
                        "from_section": section,
                        "reason": (
                            "low-specificity water/environmental organism detected only from weak non-sterile respiratory culture, "
                            "without mNGS signal or quantified/repeated support, and other more plausible pathogens are present"
                        ),
                    }
                )
                continue
            kept.append(item)
        review[section] = kept
    review["weak_environmental_review_items_omitted_count"] = len(omitted)
    if omitted:
        review["weak_environmental_review_items_omitted"] = omitted
    return review


def picked_name_sets(deterministic_max: dict[str, Any]) -> tuple[set[str], set[str]]:
    best = deterministic_max.get("best_available_summary") if isinstance(deterministic_max, dict) else {}
    picked = best.get("picked_pathogens") if isinstance(best, dict) else []
    keys: set[str] = set()
    representative_groups: set[str] = set()
    if not isinstance(picked, list):
        return keys, representative_groups
    for item in picked:
        if not isinstance(item, dict):
            continue
        name = item.get("organism_name") or item.get("name")
        key = pathogen_names.canonical_key(name)
        group = pathogen_names.representative_group_key(name)
        if key:
            keys.add(key)
        if group and not group.startswith("name:"):
            representative_groups.add(group)
    return keys, representative_groups


def review_item_tier(
    item: dict[str, Any],
    *,
    source_section: str,
    picked_keys: set[str],
    picked_genera: set[str],
) -> tuple[str, str]:
    name = item.get("organism_name") or item.get("name")
    key = pathogen_names.canonical_key(name)
    group = pathogen_names.representative_group_key(name)
    evidence = evidence_from_item(item)
    reads = numeric(evidence.get("reads"))
    dominance = str(evidence.get("dominance_tier") or "")
    confidence = str(item.get("confidence") or "")

    if key and key in picked_keys:
        return "review_omitted_with_reason", "same pathogen or alias is already retained in picked_pathogens"
    if (
        group
        and group in picked_genera
        and pathogen_names.representative_group_suppression_allowed(name)
    ):
        return "review_omitted_with_reason", "same representative alias/generic label already has a picked representative"

    if key == "corynebacteriumstriatum" and dominance == "D3_dominant" and reads >= 10000:
        return "review_high_priority", "D3 dominant, very high-burden lower-respiratory C. striatum safety signal"
    if source_section == "possible_missed_pathogens" and key in ACTIONABLE_RESPIRATORY_REVIEW_KEYS:
        return "review_high_priority", "actionable respiratory pathogen or respiratory panel target with non-picked evidence"
    if is_candida_or_generic_yeast_name(name) and "respiratory_candida_high_priority_context" in {
        str(reason) for reason in (evidence.get("review_reasons") or [])
    }:
        return (
            "review_high_priority",
            "BAL/lower-respiratory Candida mNGS signal plus highly vulnerable host and same-species or same-genus respiratory hospital context; prioritize RAG while keeping it out of formal picked without same-species strong invasive proof",
        )

    rule_ids = evidence_rule_ids(evidence)
    if is_candida_or_generic_yeast_name(name):
        if candida_has_invasive_context(evidence):
            return (
                "review_high_priority",
                "same-species Candida has strong invasive evidence such as blood, sterile-site, tissue, pleural, abscess, or histopathology; high-priority review is needed but this is not automatic Candida pneumonia attribution",
            )
        if candida_has_intermediate_context(evidence):
            return (
                "review_context_needed",
                "Candida has intermediate systemic or invasive-candidiasis clues; review clinical context, source, and pulmonary linkage before attribution",
            )
        if (
            "respiratory_candida_context_not_formal_picked" in rule_ids
            and (tier_score(evidence.get("reads_tier")) >= 20 or numeric(evidence.get("reads_percentile")) >= 0.85)
        ):
            return (
                "review_context_needed",
                "respiratory Candida/generic yeast has lower-respiratory mNGS or vulnerable-host context, but lacks same-species strong invasive evidence; RAG should assess colonization versus infection",
            )
        return (
            "review_low_specificity_audit",
            "respiratory or non-invasive Candida/generic yeast is low-specificity for pneumonia and should not enter RAG by default",
        )
    if (
        NON_CORE_RESPIRATORY_VIRUS_FORMAL_PICK_GUARDRAIL in rule_ids
        or "non_core_respiratory_virus_context_not_formal_picked" in rule_ids
    ):
        return (
            "review_context_needed",
            "non-core respiratory virus was moved out of formal picked but remains RAG-visible for symptom timing, panel/PCR, imaging, and coinfection context",
        )
    if (
        MNGS_ONLY_NONDOMINANT_NON_TYPICAL_FORMAL_PICK_GUARDRAIL in rule_ids
        or "mngs_only_nondominant_non_typical_strong_context" in rule_ids
        or "mngs_only_nondominant_non_typical_low_specificity" in rule_ids
        or "broad_or_atypical_oral_mngs_only_low_specificity" in rule_ids
    ):
        if (
            "broad_or_atypical_oral_mngs_only_low_specificity" in rule_ids
            or (
                is_broad_or_atypical_oral_flora_name(name)
                and not has_non_host_support(evidence)
                and str(evidence.get("dominance_tier") or "") not in {"D2_moderate", "D3_dominant"}
                and not has_related_representative_context(evidence)
            )
        ):
            return (
                "review_low_specificity_audit",
                "broad oral/upper-airway or atypical oral-associated organism has mNGS-only D0/D1 high-percentile signal but lacks same-organism support, D2/D3 dominance, sterile/systemic signal, or coherent aspiration pattern",
            )
        if (
            "mngs_only_nondominant_non_typical_strong_context" in rule_ids
            or tier_score(evidence.get("reads_tier")) >= 30
            or numeric(evidence.get("reads_percentile")) >= 0.85
        ):
            return (
                "review_context_needed",
                "mNGS-only D0/D1 non-dominant non-typical organism has high reads tier or percentile; review as context rather than formal picked",
            )
        return (
            "review_low_specificity_audit",
            "mNGS-only D0/D1 non-dominant non-typical organism lacks same-organism hospital support and does not have strong enough signal for RAG-visible context",
        )
    if (
        COMMON_HOSPITAL_LOW_READ_NODOM_FORMAL_PICK_GUARDRAIL in rule_ids
        or "common_hospital_low_read_nondominant_related_context" in rule_ids
        or "common_hospital_low_read_nondominant_low_specificity" in rule_ids
    ):
        if (
            "common_hospital_low_read_nondominant_related_context" in rule_ids
            or has_related_representative_context(evidence)
        ):
            return (
                "review_context_needed",
                "common hospital pneumonia pathogen has only R0/R1 non-dominant mNGS signal, but related representative hospital context merits RAG review",
            )
        return (
            "review_low_specificity_audit",
            "common hospital pneumonia pathogen has only R0/R1 non-dominant mNGS signal without direct same-organism or related hospital support; retain as audit rather than RAG-visible context",
        )

    if key == "corynebacteriumstriatum" and reads >= 1000:
        return "review_context_needed", "high-burden C. striatum needs airway-device, culture, imaging, and treatment-context review"
    if key in CONTEXT_REVIEW_KEYS:
        return "review_context_needed", "clinically plausible or high-consequence organism with incomplete pulmonary evidence"
    if key == "chryseobacteriumindologenes" and reads >= 10000:
        if has_non_host_support(evidence) or dominance in {"D2_moderate", "D3_dominant"}:
            return "review_context_needed", "water/environment-associated GNB has unusually strong local convergence"
        return (
            "review_low_specificity_audit",
            "high-read Chryseobacterium without same-organism support or D2/D3 dominance remains low pulmonary specificity",
        )
    if key == "elizabethkingiaanophelis" and confidence == "moderate":
        if has_non_host_support(evidence) or dominance in {"D2_moderate", "D3_dominant"}:
            return "review_context_needed", "potential nosocomial water-associated GNB has convergent support"
        return (
            "review_low_specificity_audit",
            "Elizabethkingia signal lacks same-organism support or D2/D3 dominance; retain as low-specificity audit",
        )
    if key in ASPIRATION_CONTEXT_KEYS and reads >= 1000:
        return "review_context_needed", "high-burden aspiration/anaerobe pattern needs clinical and imaging context"

    return (
        "review_low_specificity_audit",
        "low-specificity, colonizer/background/environmental, non-pulmonary culture-only, or weak mNGS signal; not sent to RAG",
    )


def annotate_review_item(item: dict[str, Any], *, tier: str, reason: str, source_section: str) -> dict[str, Any]:
    output = dict(item)
    canonicalize_item_name(output)
    output["review_tier"] = tier
    output["review_tier_reason"] = reason
    output["original_review_section"] = source_section
    return output


def triage_review_output(review: dict[str, Any], deterministic_max: dict[str, Any]) -> dict[str, Any]:
    picked_keys, picked_genera = picked_name_sets(deterministic_max)
    output = dict(review)
    has_direct_tiers = any(
        isinstance(review.get(key), list) and review.get(key)
        for key in (
            "review_high_priority",
            "review_context_needed",
            "review_low_specificity",
            "review_low_specificity_audit",
            "review_omitted_with_reason",
        )
    )
    if has_direct_tiers:
        low_items = review.get("review_low_specificity")
        if not isinstance(low_items, list):
            low_items = []
        legacy_low_items = review.get("review_low_specificity_audit")
        if isinstance(legacy_low_items, list):
            low_items = [*low_items, *legacy_low_items]
        output["review_high_priority"] = [
            annotate_review_item(item, tier="review_high_priority", reason=str(item.get("review_tier_reason") or "direct review tier"), source_section=str(item.get("original_review_section") or "review_high_priority"))
            for item in (review.get("review_high_priority") or [])
            if isinstance(item, dict)
        ]
        output["review_context_needed"] = [
            annotate_review_item(item, tier="review_context_needed", reason=str(item.get("review_tier_reason") or "direct review tier"), source_section=str(item.get("original_review_section") or "review_context_needed"))
            for item in (review.get("review_context_needed") or [])
            if isinstance(item, dict)
        ]
        output["review_low_specificity"] = [
            annotate_review_item(item, tier="review_low_specificity", reason=str(item.get("review_tier_reason") or "audit-only low-specificity tier"), source_section=str(item.get("original_review_section") or "review_low_specificity"))
            for item in low_items
            if isinstance(item, dict)
        ]
        direct_omitted = []
        for omitted_key in ("review_omitted_with_reason", "omitted_with_reason"):
            values = review.get(omitted_key)
            if isinstance(values, list):
                direct_omitted.extend(values)
        output["review_omitted_with_reason"] = [
            annotate_review_item(item, tier="review_omitted_with_reason", reason=str(item.get("review_tier_reason") or item.get("rationale_zh") or "omitted from RAG-visible review"), source_section=str(item.get("original_review_section") or "review_omitted_with_reason"))
            for item in direct_omitted
            if isinstance(item, dict)
        ]
        output.pop("possible_missed_pathogens", None)
        output.pop("watchlist_candidates", None)
        output.pop("review_low_specificity_audit", None)
        output.pop("omitted_with_reason", None)
        output["review_tiering_policy"] = {
            "enabled": True,
            "version": "direct_review_triage_high_context_only_v10_low_specificity_taxonomy",
            "rag_visible_tiers": ["review_high_priority", "review_context_needed"],
            "audit_only_tiers": ["review_low_specificity", "review_omitted_with_reason"],
            "counts_before_name_reconciliation": {
                "review_high_priority": len(output.get("review_high_priority") or []),
                "review_context_needed": len(output.get("review_context_needed") or []),
                "review_low_specificity": len(output.get("review_low_specificity") or []),
                "review_omitted_with_reason": len(output.get("review_omitted_with_reason") or []),
            },
            "policy_zh": (
                "review JSON 已直接輸出分層；merge 不再從 possible/watchlist 轉譯。"
                "最終 combined/RAG 輸出只保留 review_high_priority 與 review_context_needed。"
            ),
        }
        output["review_output_policy"] = (
            "Merged output uses direct review_high_priority and review_context_needed as RAG-visible missed-review tiers. "
            "review_low_specificity and review_omitted_with_reason are retained only as audit metadata."
        )
        return output

    original_counts = {
        "possible_missed_pathogens": len(review.get("possible_missed_pathogens") or [])
        if isinstance(review.get("possible_missed_pathogens"), list)
        else 0,
        "watchlist_candidates": len(review.get("watchlist_candidates") or [])
        if isinstance(review.get("watchlist_candidates"), list)
        else 0,
    }
    for key in (
        "possible_missed_pathogens",
        "watchlist_candidates",
        "review_high_priority",
        "review_context_needed",
        "review_low_specificity",
        "review_low_specificity_audit",
        "review_omitted_with_reason",
    ):
        output.pop(key, None)

    triaged: dict[str, list[dict[str, Any]]] = {
        "review_high_priority": [],
        "review_context_needed": [],
        "review_low_specificity": [],
        "review_omitted_with_reason": [],
    }
    for source_section in ("possible_missed_pathogens", "watchlist_candidates"):
        values = review.get(source_section) or []
        if not isinstance(values, list):
            continue
        for item in values:
            if not isinstance(item, dict):
                continue
            tier, reason = review_item_tier(
                item,
                source_section=source_section,
                picked_keys=picked_keys,
                picked_genera=picked_genera,
            )
            triaged[tier].append(
                annotate_review_item(item, tier=tier, reason=reason, source_section=source_section)
            )

    output.update(triaged)
    output["review_tiering_policy"] = {
        "enabled": True,
        "version": "review_triage_high_context_only_v1",
        "original_counts": original_counts,
        "rag_visible_tiers": ["review_high_priority", "review_context_needed"],
        "audit_only_tiers": ["review_low_specificity", "review_omitted_with_reason"],
        "counts_before_name_reconciliation": {key: len(value) for key, value in triaged.items()},
        "policy_zh": (
            "最終 combined/RAG 輸出只保留 review_high_priority 與 review_context_needed。"
            "review_low_specificity 與 review_omitted_with_reason 僅供追溯，不納入 combined precision/recall。"
        ),
    }
    output["review_output_policy"] = (
        "Merged output exposes review_high_priority and review_context_needed as the RAG-visible missed-review tiers. "
        "Low-specificity and duplicate/non-pulmonary items are retained only as audit metadata and are excluded from combined metrics."
    )
    return output


def apply_nonpulmonary_hospital_evidence_guardrail(review: dict[str, Any]) -> dict[str, Any]:
    """Prevent urine/Foley/stool-only evidence from supporting pneumonia tiers."""
    output = dict(review)
    sections = (
        "review_high_priority",
        "review_context_needed",
        "review_low_specificity",
        "review_omitted_with_reason",
    )
    processed: dict[str, list[dict[str, Any]]] = {section: [] for section in sections}
    affected: list[dict[str, Any]] = []

    for section in sections:
        for raw_item in as_list(review.get(section)):
            if not isinstance(raw_item, dict):
                continue
            item = copy.deepcopy(raw_item)
            evidence = evidence_from_item(item)
            profile = hospital_source_guardrail_profile(evidence)
            target = section
            action = "unchanged"
            if profile["nonpulmonary_only"]:
                strip_nonpulmonary_support_modules(evidence)
                reads = numeric(evidence.get("reads"))
                dominance = str(evidence.get("dominance_tier") or "")
                has_invasive = candida_has_invasive_context(evidence)
                if reads <= 0:
                    target = "review_omitted_with_reason"
                    action = "omit_nonpulmonary_hospital_only"
                elif section == "review_high_priority" and dominance in {
                    "",
                    "Not_available",
                    "D0_not_top",
                    "D1_low",
                } and not has_invasive:
                    target = "review_context_needed"
                    action = "downgrade_high_to_context_without_pulmonary_support"
                item["nonpulmonary_hospital_evidence_guardrail"] = {
                    **profile,
                    "action": action,
                    "original_tier": section,
                    "final_tier": target,
                    "reason": (
                        "Urine, Foley, or stool-only hospital evidence is retained for audit but cannot count as "
                        "same-organism pulmonary support."
                    ),
                }
                item["review_tier"] = target
                item["review_tier_reason"] = (
                    "hospital evidence is exclusively urine/Foley/stool sourced and is not pulmonary support"
                    if target == "review_omitted_with_reason"
                    else str(item.get("review_tier_reason") or "")
                )
                affected.append(
                    {
                        "organism_name": item.get("organism_name") or item.get("name"),
                        "original_tier": section,
                        "final_tier": target,
                        "action": action,
                    }
                )
            processed[target].append(item)

    for section, items in processed.items():
        output[section] = items
    policy = output.get("review_tiering_policy")
    if isinstance(policy, dict):
        policy["nonpulmonary_hospital_evidence_guardrail"] = {
            "enabled": True,
            "affected_count": len(affected),
            "items": affected,
        }
    output["nonpulmonary_hospital_evidence_guardrail"] = {
        "enabled": True,
        "policy": (
            "Urine/Foley/stool-only hospital evidence is audit context and cannot act as same-organism "
            "support for pulmonary causality."
        ),
        "affected_count": len(affected),
        "items": affected,
    }
    return output


def deterministic_candidate_index(deterministic_max: dict[str, Any]) -> dict[str, dict[str, Any]]:
    candidates = deterministic_max.get("pathogen_candidates") if isinstance(deterministic_max, dict) else []
    output: dict[str, dict[str, Any]] = {}
    if not isinstance(candidates, list):
        return output
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        name = candidate.get("organism_name") or candidate.get("name")
        for key in pathogen_names.expanded_name_keys(name):
            if key and key not in output:
                output[key] = candidate
    return output


def candidate_for_item(item: dict[str, Any], candidates: dict[str, dict[str, Any]]) -> dict[str, Any]:
    for key in pathogen_names.expanded_name_keys(item.get("organism_name") or item.get("name")):
        candidate = candidates.get(key)
        if candidate:
            return candidate
    return {}


def evidence_value(item: dict[str, Any], candidate: dict[str, Any], key: str) -> Any:
    evidence = evidence_from_item(item)
    candidate_key_evidence = as_dict(candidate.get("key_evidence"))
    return first_present(item.get(key), evidence.get(key), candidate.get(key), candidate_key_evidence.get(key))


def item_support_modules(item: dict[str, Any], candidate: dict[str, Any]) -> list[str]:
    evidence = evidence_from_item(item)
    modules: list[str] = []
    modules.extend(str(value) for value in as_list(item.get("support_modules")))
    modules.extend(str(value) for value in as_list(evidence.get("support_modules")))
    modules.extend(str(value) for value in as_list(candidate.get("support_modules")))
    module_summary = as_dict(candidate.get("module_support_summary"))
    modules.extend(str(key) for key, value in module_summary.items() if value == "Support")
    output: list[str] = []
    seen: set[str] = set()
    for module in modules:
        if not module or module in seen:
            continue
        seen.add(module)
        output.append(module)
    return output


def item_non_host_support_modules(item: dict[str, Any], candidate: dict[str, Any]) -> list[str]:
    evidence = evidence_from_item(item)
    modules: list[str] = []
    modules.extend(str(value) for value in as_list(evidence.get("non_host_support_modules")))
    modules.extend(module for module in item_support_modules(item, candidate) if module in DIRECT_HOSPITAL_MODULES)
    output: list[str] = []
    seen: set[str] = set()
    for module in modules:
        if not module or module == "host" or module in seen:
            continue
        seen.add(module)
        output.append(module)
    return output


def item_has_same_hospital_support(item: dict[str, Any], candidate: dict[str, Any]) -> bool:
    return bool(set(item_non_host_support_modules(item, candidate)) & DIRECT_HOSPITAL_MODULES)


def item_has_related_hospital_support(item: dict[str, Any], candidate: dict[str, Any]) -> bool:
    evidence = evidence_from_item(item)
    values = (
        item.get("same_genus_or_related_hospital_evidence"),
        item.get("related_representative_hospital_support"),
        evidence.get("related_representative_hospital_support"),
        candidate.get("related_representative_hospital_support"),
        as_dict(candidate.get("key_evidence")).get("related_representative_hospital_support"),
    )
    return any(bool(value) for value in values) or has_related_representative_context(evidence)


def lower_respiratory_evidence(item: dict[str, Any], candidate: dict[str, Any]) -> bool:
    specimen_class = str(evidence_value(item, candidate, "specimen_class") or "")
    specimen_alignment = str(evidence_value(item, candidate, "specimen_alignment") or "")
    text = f"{specimen_class} {specimen_alignment}".lower()
    return "s2_lower_respiratory" in text or "lower_respiratory" in text or "lower respiratory" in text


def clinical_ecology_group_for_convergence(item: dict[str, Any], candidate: dict[str, Any]) -> str:
    name = item.get("organism_name") or item.get("name")
    normalized = normalize_name(name)
    classification = str(evidence_value(item, candidate, "classification") or "").lower()

    if is_candida_or_generic_yeast_name(name):
        return "candida_or_generic_yeast"
    if is_rare_opportunistic_yeast_name(name):
        return "rare_opportunistic_yeast"
    if is_reactivation_virus_name(name):
        return "herpes_reactivation_virus"
    if is_core_answer_sensitive_respiratory_virus_name(name):
        return "core_respiratory_virus"
    oral_group = oral_aspiration_flora_category(name)
    if oral_group != "not_oral_aspiration_flora":
        return oral_group
    if is_environmental_water_soil_gnb_name(name):
        return "environmental_water_soil_low_specificity_gnb"
    if normalized.startswith("corynebacterium"):
        return "potentially_pathogenic_nondiphtherial_corynebacterium"
    if is_skin_airway_background_name(name):
        return "skin_or_airway_background"
    if any(token in normalized for token in ("klebsiella", "escherichia", "enterobacter", "serratia", "proteus", "cronobacter")):
        return "enterobacterales_or_hospital_gnb"
    if any(token in normalized for token in ("pseudomonas", "acinetobacter", "stenotrophomonas")):
        return "nonfermenter_hospital_gnb"
    if classification == "viral":
        return "other_virus"
    if classification == "fungal":
        return "other_fungus"
    if classification == "bacterial":
        return "other_bacterial_or_unclassified"
    return "other_or_unclassified"


def review_items_for_pattern_summary(
    review: dict[str, Any],
    deterministic_max: dict[str, Any],
    candidates: dict[str, dict[str, Any]],
) -> list[tuple[str, dict[str, Any], dict[str, Any]]]:
    best = deterministic_max.get("best_available_summary") if isinstance(deterministic_max, dict) else {}
    rows: list[tuple[str, dict[str, Any], dict[str, Any]]] = []
    for item in as_list(as_dict(best).get("picked_pathogens")):
        if isinstance(item, dict):
            rows.append(("picked", item, candidate_for_item(item, candidates)))
    for tier in ("review_high_priority", "review_context_needed", "review_low_specificity"):
        for item in as_list(review.get(tier)):
            if isinstance(item, dict):
                rows.append((tier, item, candidate_for_item(item, candidates)))
    return rows


def patient_syndrome_pattern_summary(
    review: dict[str, Any],
    deterministic_max: dict[str, Any],
    candidates: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    group_counts: dict[str, int] = {}
    for _, item, candidate in review_items_for_pattern_summary(review, deterministic_max, candidates):
        group = clinical_ecology_group_for_convergence(item, candidate)
        group_counts[group] = group_counts.get(group, 0) + 1
    oral_count = sum(
        group_counts.get(group, 0)
        for group in (
            "strict_aspiration_anaerobe",
            "broad_oral_upper_airway_commensal",
            "atypical_oral_associated",
            "other_oral_aspiration_flora",
        )
    )
    hospital_gnb_count = group_counts.get("enterobacterales_or_hospital_gnb", 0) + group_counts.get(
        "nonfermenter_hospital_gnb", 0
    )
    skin_count = group_counts.get("skin_or_airway_background", 0) + group_counts.get(
        "potentially_pathogenic_nondiphtherial_corynebacterium", 0
    )
    patterns: list[str] = []
    if oral_count >= 2 or group_counts.get("strict_aspiration_anaerobe", 0):
        patterns.append("aspiration_or_anaerobe_cluster")
    if hospital_gnb_count >= 2:
        patterns.append("hospital_gnb_cluster")
    if group_counts.get("herpes_reactivation_virus", 0):
        patterns.append("herpes_reactivation_context")
    if group_counts.get("candida_or_generic_yeast", 0) or group_counts.get("rare_opportunistic_yeast", 0):
        patterns.append("yeast_context")
    if skin_count:
        patterns.append("skin_airway_background_context")
    return {
        "patterns": patterns,
        "group_counts": group_counts,
        "oral_count": oral_count,
        "hospital_gnb_count": hospital_gnb_count,
        "skin_airway_count": skin_count,
    }


def convergence_pattern_strength(
    item: dict[str, Any],
    *,
    current_tier: str,
    candidate: dict[str, Any],
    patient_patterns: dict[str, Any],
) -> dict[str, str]:
    group = clinical_ecology_group_for_convergence(item, candidate)
    dominance = str(evidence_value(item, candidate, "dominance_tier") or "")
    reads_tier = str(evidence_value(item, candidate, "reads_tier") or "")
    reads = numeric(evidence_value(item, candidate, "reads"))
    rank = numeric(evidence_value(item, candidate, "rank_priority"), default=99)
    has_same = item_has_same_hospital_support(item, candidate)
    has_related = item_has_related_hospital_support(item, candidate)
    lower_resp = lower_respiratory_evidence(item, candidate)
    high_reads = tier_score(reads_tier) >= 30
    dominant = dominance in {"D2_moderate", "D3_dominant"}
    patterns = set(as_list(patient_patterns.get("patterns")))

    if current_tier == "picked":
        strength, role = "formal", "lead_or_selected_pathogen"
    elif group == "strict_aspiration_anaerobe":
        if has_same or dominant or ("aspiration_or_anaerobe_cluster" in patterns and lower_resp and high_reads):
            strength, role = "moderate", "syndrome_representative_context"
        else:
            strength, role = "weak", "single_anaerobe_context"
    elif group in {"broad_oral_upper_airway_commensal", "atypical_oral_associated", "other_oral_aspiration_flora"}:
        if "aspiration_or_anaerobe_cluster" in patterns and (has_same or dominant):
            strength, role = "moderate", "supporting_syndrome_member"
        else:
            strength, role = "weak", "background_oral_member"
    elif group == "potentially_pathogenic_nondiphtherial_corynebacterium":
        if has_same and (high_reads or dominant):
            strength, role = "moderate", "possible_device_or_airway_pathogen"
        elif high_reads or dominant:
            strength, role = "weak_to_moderate", "high-burden nondiphtherial Corynebacterium needs case-level review"
        else:
            strength, role = "weak", "background_corynebacterium_signal"
    elif group == "skin_or_airway_background":
        if has_same and (high_reads or dominant):
            strength, role = "moderate", "possible_device_or_airway_pathogen"
        elif has_related and high_reads:
            strength, role = "weak_to_moderate", "related_background_signal"
        else:
            strength, role = "weak", "background_skin_airway_signal"
    elif group == "environmental_water_soil_low_specificity_gnb":
        if has_same or dominant:
            strength, role = "moderate", "environmental_gnb_with_convergence"
        else:
            strength, role = "weak", "environmental_low_specificity_signal"
    elif group == "candida_or_generic_yeast":
        if candida_has_invasive_context(evidence_from_item(item)):
            strength, role = "strong", "candida_strong_invasive_context"
        elif candida_has_intermediate_context(evidence_from_item(item)) or (has_same and lower_resp and high_reads):
            strength, role = "moderate", "candida_respiratory_or_systemic_context"
        else:
            strength, role = "weak_to_moderate", "candida_colonization_vs_infection_context"
    elif group == "herpes_reactivation_virus":
        if dominant or has_same:
            strength, role = "moderate", "viral_reactivation_with_direct_support"
        else:
            strength, role = "weak_to_moderate", "viral_reactivation_context"
    elif group in {"enterobacterales_or_hospital_gnb", "nonfermenter_hospital_gnb", "core_respiratory_virus"}:
        if has_same or dominant:
            strength, role = "strong", "actionable_respiratory_pathogen_context"
        else:
            strength, role = "moderate", "actionable_but_incomplete_context"
    elif group in {"other_fungus", "other_bacterial_or_unclassified", "other_or_unclassified"}:
        if has_same or dominant:
            strength, role = "moderate", "non-core organism with direct convergence"
        elif high_reads and lower_resp and rank <= 3:
            strength, role = "weak_to_moderate", "high-rank lower-respiratory signal needs case-level review"
        else:
            strength, role = "weak", "weak non-core signal without convergence"
    elif group == "other_virus":
        if has_same and reads <= 0 and not dominant:
            strength, role = "weak", "panel-positive non-core virus without mNGS burden"
        elif has_same or dominant:
            strength, role = "moderate", "non-core virus with direct convergence"
        else:
            strength, role = "weak", "low-actionability or duplicate viral context"
    else:
        strength, role = "weak", "needs_rule_review"
    return {
        "clinical_ecology_group": group,
        "pattern_strength": strength,
        "pattern_role": role,
    }


def suggested_convergence_tier(current_tier: str, assessment: dict[str, str]) -> tuple[str, str]:
    if current_tier not in {"review_high_priority", "review_context_needed"}:
        return current_tier, "syndrome-pattern convergence only reviews RAG-visible advisory tiers"
    group = assessment["clinical_ecology_group"]
    strength = assessment["pattern_strength"]
    role = assessment["pattern_role"]
    if strength == "weak" and group in GENERIC_WEAK_REVIEW_DEMOTION_GROUPS:
        return (
            "review_low_specificity",
            f"syndrome-pattern convergence: {group} has weak local convergence ({role}); keep as audit rather than RAG-visible",
        )
    return current_tier, f"syndrome-pattern convergence kept tier: {group} / {strength} / {role}"


def apply_syndrome_pattern_convergence(review: dict[str, Any], deterministic_max: dict[str, Any]) -> dict[str, Any]:
    output = dict(review)
    for section in ("review_high_priority", "review_context_needed", "review_low_specificity", "review_omitted_with_reason"):
        if not isinstance(output.get(section), list):
            output[section] = []
    candidates = deterministic_candidate_index(deterministic_max)
    patient_patterns = patient_syndrome_pattern_summary(output, deterministic_max, candidates)
    counts_before = {
        section: len(output.get(section) or [])
        for section in ("review_high_priority", "review_context_needed", "review_low_specificity", "review_omitted_with_reason")
    }
    demoted_items: list[dict[str, Any]] = []
    for source_section in ("review_high_priority", "review_context_needed"):
        kept: list[dict[str, Any]] = []
        for item in as_list(output.get(source_section)):
            if not isinstance(item, dict):
                continue
            candidate = candidate_for_item(item, candidates)
            assessment = convergence_pattern_strength(
                item,
                current_tier=source_section,
                candidate=candidate,
                patient_patterns=patient_patterns,
            )
            target_tier, reason = suggested_convergence_tier(source_section, assessment)
            if target_tier == source_section:
                kept.append(item)
                continue
            moved = dict(item)
            moved["pre_convergence_review_tier"] = source_section
            moved["pre_convergence_review_tier_reason"] = item.get("review_tier_reason")
            moved["review_tier"] = target_tier
            moved["review_tier_reason"] = reason
            moved["syndrome_pattern_convergence"] = {
                "action": f"demote_{source_section}_to_{target_tier}",
                "previous_tier": source_section,
                "new_tier": target_tier,
                "clinical_ecology_group": assessment["clinical_ecology_group"],
                "pattern_strength": assessment["pattern_strength"],
                "pattern_role": assessment["pattern_role"],
                "patient_patterns": patient_patterns.get("patterns") or [],
            }
            output[target_tier].append(moved)
            demoted_items.append(
                {
                    "organism_name": moved.get("organism_name") or moved.get("name"),
                    "from_tier": source_section,
                    "to_tier": target_tier,
                    "clinical_ecology_group": assessment["clinical_ecology_group"],
                    "pattern_strength": assessment["pattern_strength"],
                    "pattern_role": assessment["pattern_role"],
                    "reason": reason,
                }
            )
        output[source_section] = kept
    counts_after = {
        section: len(output.get(section) or [])
        for section in ("review_high_priority", "review_context_needed", "review_low_specificity", "review_omitted_with_reason")
    }
    output["syndrome_pattern_convergence_policy"] = {
        "enabled": True,
        "version": "syndrome_pattern_convergence_v1_conservative_demotions",
        "placement": "after LLM review triage and before final name reconciliation",
        "promotions_enabled": False,
        "formal_picked_modified": False,
        "rag_visible_tiers": ["review_high_priority", "review_context_needed"],
        "audit_only_tiers": ["review_low_specificity", "review_omitted_with_reason"],
        "patient_patterns": patient_patterns,
        "counts_before": counts_before,
        "counts_after": counts_after,
        "demoted_count": len(demoted_items),
        "demoted_items": demoted_items,
        "policy_zh": (
            "在 LLM review 之後、merge 之前，將缺乏同菌院端支持、D2/D3 dominance 或 syndrome pattern 的"
            "低特異性 context/high 候選降到 review_low_specificity；不升級候選，也不修改 picked。"
        ),
    }
    policy = output.get("review_tiering_policy")
    if isinstance(policy, dict):
        policy["syndrome_pattern_convergence_applied"] = True
        policy["counts_before_syndrome_pattern_convergence"] = counts_before
        policy["counts_after_syndrome_pattern_convergence"] = counts_after
    output["review_output_policy"] = (
        str(output.get("review_output_policy") or "")
        + " Syndrome-pattern convergence demotes weak low-specificity advisory items to review_low_specificity before final name reconciliation."
    ).strip()
    return output


def review_tier_convergence_decision(
    item: dict[str, Any],
    candidate: dict[str, Any],
) -> tuple[str, str]:
    """Return a formal low-risk context demotion rule and reason, if any."""
    name = item.get("organism_name") or item.get("name")
    key = pathogen_names.canonical_key(name)
    group = clinical_ecology_group_for_convergence(item, candidate)
    dominance = str(evidence_value(item, candidate, "dominance_tier") or "")
    reads_tier = str(evidence_value(item, candidate, "reads_tier") or "")
    same_species_support = item_has_same_hospital_support(item, candidate)
    related_support = item_has_related_hospital_support(item, candidate)
    non_host_support = bool(item_non_host_support_modules(item, candidate))
    evidence = evidence_from_item(item)

    if group == "candida_or_generic_yeast":
        weak_signal = reads_tier.startswith(("R0", "R1", "R2"))
        no_invasive_or_systemic_rescue = not (
            evidence_profiles.candida_has_invasive_context(evidence)
            or evidence_profiles.candida_has_intermediate_context(evidence)
        )
        if (
            dominance.startswith(("D0", "D1"))
            and weak_signal
            and not same_species_support
            and not related_support
            and not non_host_support
            and no_invasive_or_systemic_rescue
        ):
            return (
                REVIEW_TIER_CONVERGENCE_WEAK_CANDIDA_RULE,
                "weak non-dominant respiratory Candida lacks exact-species hospital, invasive/systemic, or related convergent support",
            )

    if group == "potentially_pathogenic_nondiphtherial_corynebacterium" and key != "corynebacteriumstriatum":
        if dominance.startswith(("D0", "D1")) and not same_species_support and not non_host_support:
            return (
                REVIEW_TIER_CONVERGENCE_NONSTRIATUM_CORYNE_RULE,
                "non-striatum Corynebacterium lacks dominance and exact-species hospital support",
            )
    return "", ""


def apply_review_tier_convergence_v2(
    review: dict[str, Any],
    deterministic_max: dict[str, Any],
) -> dict[str, Any]:
    """Formally demote agreed weak context candidates without deleting audit evidence."""
    output = dict(review)
    for section in (
        "review_high_priority",
        "review_context_needed",
        "review_low_specificity",
        "review_omitted_with_reason",
    ):
        if not isinstance(output.get(section), list):
            output[section] = []

    candidates = deterministic_candidate_index(deterministic_max)
    counts_before = {
        section: len(output[section])
        for section in (
            "review_high_priority",
            "review_context_needed",
            "review_low_specificity",
            "review_omitted_with_reason",
        )
    }
    kept: list[dict[str, Any]] = []
    demoted_items: list[dict[str, Any]] = []
    for raw_item in output["review_context_needed"]:
        if not isinstance(raw_item, dict):
            continue
        item = copy.deepcopy(raw_item)
        candidate = candidate_for_item(item, candidates)
        rule_id, reason = review_tier_convergence_decision(item, candidate)
        if not rule_id:
            kept.append(item)
            continue
        item["pre_review_tier_convergence_tier"] = "review_context_needed"
        item["pre_review_tier_convergence_reason"] = item.get("review_tier_reason")
        item["review_tier"] = "review_low_specificity"
        item["review_tier_reason"] = reason
        item["review_tier_convergence"] = {
            "version": "review_tier_convergence_v2_agreed_low_risk",
            "rule_id": rule_id,
            "action": "demote_context_to_low",
            "previous_tier": "review_context_needed",
            "new_tier": "review_low_specificity",
            "reason": reason,
            "formal_picked_modified": False,
        }
        output["review_low_specificity"].append(item)
        demoted_items.append(
            {
                "organism_name": item.get("organism_name") or item.get("name"),
                "from_tier": "review_context_needed",
                "to_tier": "review_low_specificity",
                "rule_id": rule_id,
                "reason": reason,
            }
        )
    output["review_context_needed"] = kept

    counts_after = {
        section: len(output[section])
        for section in (
            "review_high_priority",
            "review_context_needed",
            "review_low_specificity",
            "review_omitted_with_reason",
        )
    }
    output["review_tier_convergence_policy"] = {
        "enabled": True,
        "version": "review_tier_convergence_v2_agreed_low_risk",
        "placement": "after syndrome convergence and before final name reconciliation",
        "formal_picked_modified": False,
        "promotions_enabled": False,
        "rag_visible_tiers": ["review_high_priority", "review_context_needed"],
        "audit_only_tiers": ["review_low_specificity", "review_omitted_with_reason"],
        "counts_before": counts_before,
        "counts_after": counts_after,
        "demoted_count": len(demoted_items),
        "demoted_items": demoted_items,
        "policy_zh": (
            "將缺乏收斂證據的弱 respiratory Candida，以及非 C. striatum、D0/D1、"
            "無同菌院端支持的 Corynebacterium，從 review_context_needed 降至 review_low_specificity；"
            "保留完整理由與證據，不修改 picked。"
        ),
    }
    policy = output.get("review_tiering_policy")
    if isinstance(policy, dict):
        policy["review_tier_convergence_v2_applied"] = True
        policy["counts_before_review_tier_convergence_v2"] = counts_before
        policy["counts_after_review_tier_convergence_v2"] = counts_after
    output["review_output_policy"] = (
        str(output.get("review_output_policy") or "")
        + " Review-tier convergence v2 formally moves agreed weak Candida and non-striatum Corynebacterium context items to audit-only low specificity."
    ).strip()
    return output


def apply_atypical_respiratory_serology_context_rescue(
    review: dict[str, Any],
) -> dict[str, Any]:
    """Restore answer-blind, context-only review for indirect atypical-pathogen serology."""
    output = dict(review)
    for section in (
        "review_high_priority",
        "review_context_needed",
        "review_low_specificity",
        "review_omitted_with_reason",
    ):
        if not isinstance(output.get(section), list):
            output[section] = []

    kept_low: list[dict[str, Any]] = []
    rescued: list[dict[str, Any]] = []
    for raw_item in output["review_low_specificity"]:
        if not isinstance(raw_item, dict):
            continue
        item = copy.deepcopy(raw_item)
        evidence = evidence_from_item(item)
        profile = evidence_profiles.atypical_respiratory_serology_profile(
            item.get("organism_name") or item.get("name"),
            as_list(evidence.get("hospital_evidence")),
        )
        original_section = str(
            first_present(
                item.get("original_review_section"),
                item.get("pre_review_tier_convergence_tier"),
            )
            or ""
        )
        if not (
            profile["context_rescue_match"]
            and original_section
            in {"review_high_priority", "review_context_needed", "possible_missed_pathogens"}
        ):
            kept_low.append(item)
            continue

        previous_reason = item.get("review_tier_reason")
        item["pre_atypical_serology_context_rescue_tier"] = "review_low_specificity"
        item["pre_atypical_serology_context_rescue_reason"] = previous_reason
        item["review_tier"] = "review_context_needed"
        item["review_tier_reason"] = f"{profile['rule_id']}: {profile['reason']}"
        item["atypical_serology_context_rescue"] = {
            "version": "atypical_respiratory_serology_context_rescue_v1",
            "rule_id": profile["rule_id"],
            "action": "restore_low_to_context",
            "previous_tier": "review_low_specificity",
            "new_tier": "review_context_needed",
            "maximum_tier": profile["maximum_tier"],
            "indirect_serology_found": profile["indirect_serology_found"],
            "direct_molecular_found": profile["direct_molecular_found"],
            "formal_picked_modified": False,
            "answer_blind": True,
        }
        output["review_context_needed"].append(item)
        rescued.append(
            {
                "organism_name": item.get("organism_name") or item.get("name"),
                "from_tier": "review_low_specificity",
                "to_tier": "review_context_needed",
                "rule_id": profile["rule_id"],
            }
        )

    output["review_low_specificity"] = kept_low
    output["atypical_serology_context_rescue_policy"] = {
        "enabled": True,
        "version": "atypical_respiratory_serology_context_rescue_v1",
        "formal_picked_modified": False,
        "answer_blind": True,
        "target_tier": "review_context_needed",
        "maximum_tier": "review_context_needed",
        "rescued_count": len(rescued),
        "rescued_items": rescued,
        "policy_zh": (
            "已知非典型呼吸道病原若有結構化間接 IgM/血清學證據，且 LLM 原始病例判斷為 context/high，"
            "可由 low 恢復至 review_context_needed；不得升至 high 或 picked。"
        ),
    }
    output["review_output_policy"] = (
        str(output.get("review_output_policy") or "")
        + " Atypical-respiratory indirect serology may restore an original case-fit review item from low to context only; it never promotes high or picked."
    ).strip()
    return output


def refresh_review_tiering_counts(review: dict[str, Any]) -> None:
    policy = review.get("review_tiering_policy")
    if not isinstance(policy, dict):
        return
    policy["counts_after_name_reconciliation"] = {
        "review_high_priority": len(review.get("review_high_priority") or [])
        if isinstance(review.get("review_high_priority"), list)
        else 0,
        "review_context_needed": len(review.get("review_context_needed") or [])
        if isinstance(review.get("review_context_needed"), list)
        else 0,
        "review_low_specificity": len(review.get("review_low_specificity") or [])
        if isinstance(review.get("review_low_specificity"), list)
        else 0,
        "review_low_specificity_audit": len(review.get("review_low_specificity_audit") or [])
        if isinstance(review.get("review_low_specificity_audit"), list)
        else 0,
        "review_omitted_with_reason": len(review.get("review_omitted_with_reason") or [])
        if isinstance(review.get("review_omitted_with_reason"), list)
        else 0,
    }


def read_optional_json(path: Path | None) -> Any:
    if path is None or not path.is_file():
        return None
    return read_json(path)


def queue_path_from_review(review: dict[str, Any]) -> Path | None:
    source_files = review.get("source_files")
    if not isinstance(source_files, dict):
        return None
    queue = source_files.get("queue")
    if not queue:
        return None
    return Path(str(queue))


def build_queue_evidence_index(queue_payload: Any) -> dict[str, dict[str, Any]]:
    if not isinstance(queue_payload, dict):
        return {}
    candidates = queue_payload.get("review_queue")
    if not isinstance(candidates, list):
        candidates = []
    index: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        for key in pathogen_names.expanded_name_keys(candidate.get("organism_name")):
            if key and key not in index:
                index[key] = candidate
    return index


def evidence_snapshot(candidate: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "classification",
        "integrated_causative_level",
        "mngs_signal_tier",
        "rank_priority",
        "rank_rule",
        "reads",
        "reads_tier",
        "reads_percentile",
        "dominance_tier",
        "specimen_class",
        "specimen_alignment",
        "source_category",
        "is_protected_pathogen",
        "is_likely_colonizer_or_background",
        "support_modules",
        "non_host_support_modules",
        "guardrail_rule",
        "applied_rules",
        "review_reasons",
        "review_source",
        "hospital_review_source",
        "hospital_evidence",
        "related_representative_hospital_support",
        "candida_invasive_hospital_support",
        "candida_related_invasive_hospital_context",
        "candida_evidence_strength_profile",
        "candida_related_evidence_strength_profile",
    ]
    return {key: candidate.get(key) for key in keys if candidate.get(key) not in (None, "", [], {})}


def enrich_review_items(items: Any, evidence_index: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    if not isinstance(items, list):
        return []
    enriched: list[dict[str, Any]] = []
    for item in items:
        if not isinstance(item, dict):
            continue
        output = dict(item)
        evidence = None
        for key in pathogen_names.expanded_name_keys(output.get("organism_name")):
            evidence = evidence_index.get(key)
            if evidence:
                break
        if evidence:
            queue_snapshot = evidence_snapshot(evidence)
            existing_snapshot = output.get("evidence_snapshot")
            if not isinstance(existing_snapshot, dict):
                existing_snapshot = {}
            output["evidence_snapshot"] = {**queue_snapshot, **existing_snapshot}
        enriched.append(output)
    return enriched


def sanitize_and_enrich_review(review: dict[str, Any]) -> dict[str, Any]:
    queue_payload = read_optional_json(queue_path_from_review(review))
    evidence_index = build_queue_evidence_index(queue_payload)
    output: dict[str, Any] = {}
    for key, value in review.items():
        if key == "not_recommended_for_upgrade":
            continue
        if key in {
            "possible_missed_pathogens",
            "review_high_priority",
            "review_context_needed",
            "review_low_specificity",
            "review_low_specificity_audit",
            "review_omitted_with_reason",
            "omitted_with_reason",
        }:
            output[key] = enrich_review_items(value, evidence_index)
        elif key == "watchlist_candidates":
            output[key] = enrich_review_items(value, evidence_index)
        else:
            output[key] = value
    not_recommended = review.get("not_recommended_for_upgrade")
    output["not_recommended_for_upgrade_omitted"] = True
    output["not_recommended_for_upgrade_omitted_count"] = len(not_recommended) if isinstance(not_recommended, list) else 0
    output["review_output_policy"] = (
        "Merged output keeps review_high_priority and review_context_needed as RAG-visible review tiers. "
        "review_low_specificity and omitted/not-recommended items are retained only as audit metadata. "
        "Each retained review item is enriched with review-queue evidence when available."
    )
    return output


def main() -> int:
    args = parse_args()
    requested = normalize_requested(args.patients)
    written: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    patient_dirs = sorted(args.patient_root.glob("NGS_patient_*_json"), key=sort_key)
    for patient_dir in patient_dirs:
        pid = patient_id(patient_dir)
        if requested and pid not in requested:
            continue
        summary_dir = patient_dir / "summary_outputs"
        max_path = summary_dir / f"NGS_patient_{pid}_{args.max_suffix}.json"
        review_path = summary_dir / f"NGS_patient_{pid}_{args.review_suffix}.json"
        output_path = summary_dir / f"NGS_patient_{pid}_{args.output_suffix}.json"
        missing = [str(path) for path in (max_path, review_path) if not path.is_file()]
        if missing:
            skipped.append({"patient_id": pid, "reason": "missing_input", "details": "; ".join(missing)})
            print(f"skip {patient_dir.name}: missing input")
            continue
        if output_path.exists() and not args.overwrite:
            skipped.append({"patient_id": pid, "reason": "output_exists", "details": str(output_path)})
            print(f"skip {patient_dir.name}: output exists")
            continue
        deterministic_max = copy.deepcopy(read_json(max_path))
        review = read_json(review_path)
        enriched_review = sanitize_and_enrich_review(review)
        enriched_review = filter_weak_environmental_review_items(enriched_review, deterministic_max)
        enriched_review = triage_review_output(enriched_review, deterministic_max)
        enriched_review = apply_nonpulmonary_hospital_evidence_guardrail(enriched_review)
        enriched_review = apply_syndrome_pattern_convergence(enriched_review, deterministic_max)
        enriched_review = apply_review_tier_convergence_v2(enriched_review, deterministic_max)
        enriched_review = apply_atypical_respiratory_serology_context_rescue(enriched_review)
        enriched_review = apply_species_preferred_group_representatives(enriched_review)
        attach_taxonomy_profiles(deterministic_max, enriched_review)
        reconciliation = reconcile_sections(deterministic_max, enriched_review)
        refresh_review_tiering_counts(enriched_review)
        payload = {
            "patient_id": pid,
            "report_type": "deterministic_max_with_llm_missed_candidate_review",
            "interpretation_note": "LLM review is advisory only and does not override deterministic picked_pathogens.",
            "deterministic_max": deterministic_max,
            "llm_missed_candidate_review": enriched_review,
            "final_name_reconciliation": reconciliation,
            "merged_source_files": {
                "deterministic_max": str(max_path),
                "llm_missed_candidate_review": str(review_path),
            },
        }
        write_json(output_path, payload)
        written.append(
            {
                "patient_id": pid,
                "deterministic_picked_count": len((deterministic_max.get("best_available_summary") or {}).get("picked_pathogens") or []),
                "review_high_priority_count": len(enriched_review.get("review_high_priority") or []),
                "review_context_needed_count": len(enriched_review.get("review_context_needed") or []),
                "review_low_specificity_count": len(enriched_review.get("review_low_specificity") or []),
                "review_omitted_with_reason_count": len(enriched_review.get("review_omitted_with_reason") or []),
                "output_file": str(output_path),
            }
        )
        print(f"wrote {output_path}")
    if args.summary_csv:
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.summary_csv.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "patient_id",
                    "deterministic_picked_count",
                    "review_high_priority_count",
                    "review_context_needed_count",
                    "review_low_specificity_count",
                    "review_omitted_with_reason_count",
                    "output_file",
                ],
            )
            writer.writeheader()
            writer.writerows(written)
    print(f"written_count={len(written)}")
    print(f"skipped_count={len(skipped)}")
    for item in skipped:
        print(f"skipped patient={item['patient_id']} reason={item['reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

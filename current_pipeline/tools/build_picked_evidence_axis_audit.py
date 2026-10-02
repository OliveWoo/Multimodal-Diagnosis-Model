"""Build a category-aware evidence-axis shadow audit for picked pathogens.

The audit never edits patient outputs. It separates shared evidence extraction
from category-specific recommendation rules, then uses the answer CSV only for
post-hoc evaluation. Answer labels are not inputs to recommendation logic.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import audit_picked_false_positives as picked_audit  # noqa: E402
from tools import build_syndrome_pattern_shadow_audit as syndrome  # noqa: E402
from tools import compare_aliasclean_baseline as compare  # noqa: E402
from tools import deterministic_mngs_max_scorer as scorer  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402


DIRECT_HOSPITAL_MODULES = {"culture", "filmarray_gmtest", "molecular_microbiology"}
MOLECULAR_MODULES = {"filmarray_gmtest", "molecular_microbiology"}
DOMINANT_TIERS = {"D2_moderate", "D3_dominant"}
HIGH_READ_TIERS = {"R3_high", "R4_very_high"}
LOW_READ_TIERS = {"R0_trace", "R1_low"}
VISIBLE_REVIEW_TIERS = {"review_high_priority", "review_context_needed"}

FIELDS = [
    "label",
    "patient_id",
    "organism_name",
    "classification",
    "clinical_category",
    "category_rule_id",
    "baseline_answer_hit",
    "baseline_match_type",
    "baseline_matched_answer",
    "strict_species_answer_hit",
    "strict_species_matched_answer",
    "basis_level",
    "picked_role",
    "evidence_source",
    "rank_priority",
    "reads",
    "reads_tier",
    "reads_percentile",
    "dominance_tier",
    "is_d2_d3_dominant",
    "is_r3_r4",
    "mngs_specimen_sites",
    "mngs_specimen_buckets",
    "exact_species_hospital_support",
    "exact_species_pulmonary_support",
    "exact_species_sterile_support",
    "exact_species_nonpulmonary_only",
    "exact_species_hospital_modules",
    "exact_species_hospital_sites",
    "exact_species_hospital_buckets",
    "related_genus_or_group_support",
    "host_vulnerability_tier",
    "oral_aspiration_cluster_count",
    "explicit_aspiration_syndrome_support",
    "low_specificity_flag",
    "stronger_competing_pathogen",
    "stronger_competing_reason",
    "evidence_axis_score",
    "suggested_shadow_tier",
    "suggested_shadow_reason",
    "guarded_shadow_tier",
    "guarded_rule_id",
    "guarded_reason",
]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a category-aware picked evidence-axis shadow audit."
    )
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--merged-suffix", required=True)
    parser.add_argument("--answer-csv", type=Path, required=True)
    parser.add_argument("--answer-column", default="auto")
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def resolve_source_path(value: Any) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    path = Path(text)
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path


def clean_name(value: Any) -> str:
    return pathogen_names.clean_display_text(value)


def canonical_key(value: Any) -> str:
    return pathogen_names.canonical_key(value)


def candidate_for_pick(
    picked: dict[str, Any], candidates: list[dict[str, Any]]
) -> dict[str, Any]:
    return picked_audit.best_candidate_for_pick(picked, candidates)


def classify_source_bucket(site: Any, category: Any = "") -> str:
    text = f"{site} {category}".lower()
    if "balf" in text or re.search(r"\bbal\b", text):
        return "BAL"
    if any(term in text for term in ("endotracheal", "tracheal aspirate", "ett", "endotracheal aspirate")):
        return "ET_aspirate"
    if "sputum" in text:
        return "sputum"
    if "blood" in text:
        return "blood"
    if any(term in text for term in ("pleural", "abscess", "tissue", "histolog", "cerebrospinal", "csf")):
        return "sterile_site_other"
    if "sterile_site" in text or "sterile site" in text:
        return "sterile_site_other"
    if any(term in text for term in ("foley", "urine", "urinary")):
        return "urine_or_foley"
    if any(term in text for term in ("stool", "fecal", "faecal", "gastrointestinal")):
        return "stool_or_gi"
    if any(term in text for term in ("nasopharyn", "oropharyn", "throat", "oral", "saliva")):
        return "upper_respiratory_or_oral"
    if any(term in text for term in ("bronch", "lower_respiratory", "lower respiratory")):
        return "lower_respiratory_other"
    return "unknown"


def evidence_site(record: dict[str, Any]) -> tuple[str, str]:
    site = next(
        (
            str(record.get(key) or "").strip()
            for key in (
                "specimen_type",
                "filmarray_sample_type",
                "sample_type",
                "specimen_site",
                "sample",
                "site",
            )
            if str(record.get(key) or "").strip()
        ),
        "Unknown",
    )
    category = str(record.get("specimen_category") or record.get("sample_category") or "")
    return site, classify_source_bucket(site, category)


def hospital_index(summary: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    output: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for item in as_list(summary.get("hospital_organism_evidence")):
        if not isinstance(item, dict):
            continue
        key = canonical_key(item.get("organism_name"))
        if key:
            output[key].append(item)
    return output


def hospital_profile(entries: list[dict[str, Any]]) -> dict[str, Any]:
    modules: set[str] = set()
    sites: set[str] = set()
    buckets: set[str] = set()
    for entry in entries:
        modules.update(str(value) for value in as_list(entry.get("evidence_modules")) if str(value))
        for module, records in as_dict(entry.get("module_evidence")).items():
            if records:
                modules.add(str(module))
            for record in as_list(records):
                if not isinstance(record, dict):
                    continue
                site, bucket = evidence_site(record)
                sites.add(site)
                buckets.add(bucket)
    direct = bool(modules & DIRECT_HOSPITAL_MODULES)
    pulmonary = bool(buckets & {"BAL", "ET_aspirate", "sputum", "lower_respiratory_other"})
    sterile = bool(buckets & {"blood", "sterile_site_other"})
    nonpulmonary_only = direct and not pulmonary and not sterile and bool(
        buckets & {"urine_or_foley", "stool_or_gi", "upper_respiratory_or_oral"}
    )
    return {
        "direct": direct,
        "pulmonary": pulmonary,
        "sterile": sterile,
        "nonpulmonary_only": nonpulmonary_only,
        "modules": sorted(modules),
        "sites": sorted(sites),
        "buckets": sorted(buckets),
    }


def load_final_summary(payload: dict[str, Any]) -> dict[str, Any]:
    deterministic = as_dict(payload.get("deterministic_max"))
    path = resolve_source_path(as_dict(deterministic.get("source_files")).get("final_summary"))
    if path and path.is_file():
        value = read_json(path)
        return as_dict(value)
    return {}


def load_ranked_patient(
    payload: dict[str, Any], pid: str, cache: dict[Path, dict[str, Any]]
) -> dict[str, Any]:
    deterministic = as_dict(payload.get("deterministic_max"))
    path = resolve_source_path(as_dict(deterministic.get("source_files")).get("ranked_mngs"))
    if not path or not path.is_file():
        return {}
    if path not in cache:
        value = read_json(path)
        cache[path] = as_dict(value)
    return as_dict(as_dict(cache[path].get("patients")).get(pid))


def mngs_sites_for_name(patient: dict[str, Any], name: str) -> tuple[list[str], list[str]]:
    key = canonical_key(name)
    sites: set[str] = set()
    buckets: set[str] = set()
    for record in as_list(patient.get("records")):
        if not isinstance(record, dict):
            continue
        if not any(
            canonical_key(candidate.get("organism_name")) == key
            for candidate in as_list(record.get("candidates"))
            if isinstance(candidate, dict)
        ):
            continue
        site = str(record.get("specimen_site") or "Unknown")
        sites.add(site)
        buckets.add(classify_source_bucket(site))
    return sorted(sites), sorted(buckets)


def classification_for(item: dict[str, Any], candidate: dict[str, Any]) -> str:
    return str(item.get("classification") or candidate.get("classification") or "Unknown")


def clinical_category(name: str, classification: str) -> str:
    key = canonical_key(name)
    genus = pathogen_names.genus_name(name)
    probe = {"organism_name": name, "canonical_key": key, "genus": genus}
    if picked_audit.is_low_actionability_virus(
        probe
    ):
        return "low_actionability_virus"
    if picked_audit.is_high_consequence_pathogen(probe):
        return "high-consequence/opportunistic respiratory pathogen"
    base = syndrome.clinical_ecology_group(name, classification)
    if key == "corynebacteriumstriatum":
        return "corynebacterium_striatum"
    if genus == "corynebacterium":
        return "other_nondiphtherial_corynebacterium"
    if genus == "enterococcus":
        return "enterococcus_gi_urinary_colonizer_prone"
    if base == "other_bacterial_or_unclassified" and picked_audit.is_common_hospital_respiratory_pathogen(
        {"organism_name": name, "canonical_key": key, "genus": genus}
    ):
        return "common_hospital_respiratory_pathogen"
    return base


def guarded_recommendation(row: dict[str, Any]) -> tuple[str, str, str]:
    """Apply only previously agreed, low-blast-radius category rules."""
    category = row["clinical_category"]
    dominant = row["is_d2_d3_dominant"] == "yes"
    high_reads = row["is_r3_r4"] == "yes"
    exact = row["exact_species_hospital_support"] == "yes"
    pulmonary = row["exact_species_pulmonary_support"] == "yes"
    sterile = row["exact_species_sterile_support"] == "yes"
    cluster = int(row["oral_aspiration_cluster_count"])
    syndrome_support = row["explicit_aspiration_syndrome_support"] == "yes"

    if category == "low_actionability_virus":
        return "review_context_needed", "PICK-GUARDED-LOW-ACTION-VIRUS", "Agreed low-actionability viral DNA guardrail."
    if category == "candida_or_generic_yeast" and not sterile:
        if pulmonary and (high_reads or dominant):
            return "review_high_priority", "PICK-GUARDED-CANDIDA-HIGH", "Respiratory Candida convergence is RAG-visible but lacks invasive exact-species evidence."
        return "review_context_needed", "PICK-GUARDED-CANDIDA-CONTEXT", "Candida without exact-species invasive evidence is not formal picked."
    if category == "environmental_water_soil_low_specificity_gnb" and not exact and not dominant:
        return "review_context_needed", "PICK-GUARDED-ENV-GNB", "Weak environmental GNB lacks exact-species hospital support and dominance."
    if category in {"broad_oral_upper_airway_commensal", "atypical_oral_associated"}:
        convergent_pattern = sterile or (pulmonary and (high_reads or dominant) and (cluster >= 2 or syndrome_support))
        if not convergent_pattern:
            return "review_context_needed", "PICK-GUARDED-ORAL", "Single low-specificity oral signal lacks pulmonary syndrome convergence."
    if category in {"mold/opportunistic fungus", "other_fungus"} and not exact and not dominant:
        if row["reads_tier"] in {"R0_trace", "R1_low", "R2_medium", "Not_available", "Unknown"}:
            return "review_context_needed", "PICK-GUARDED-WEAK-MOLD", "Weak mold/fungal signal lacks exact-species hospital support and dominance."
    return "keep_picked", "PICK-GUARDED-KEEP", "No agreed low-blast-radius demotion rule applies."


def related_support(item: dict[str, Any], candidate: dict[str, Any]) -> bool:
    key_evidence = as_dict(candidate.get("key_evidence"))
    return bool(
        item.get("same_genus_or_related_hospital_evidence")
        or item.get("related_representative_hospital_support")
        or candidate.get("related_representative_hospital_support")
        or key_evidence.get("related_representative_hospital_support")
    )


def host_tier(payload: dict[str, Any]) -> str:
    deterministic = as_dict(payload.get("deterministic_max"))
    return str(as_dict(deterministic.get("host_context")).get("host_vulnerability_tier") or "Unknown")


def merged_names(payload: dict[str, Any]) -> list[str]:
    return [
        clean_name(item.get("organism_name") or item.get("name"))
        for _, item in syndrome.merged_items(payload)
        if clean_name(item.get("organism_name") or item.get("name"))
    ]


def oral_cluster_count(payload: dict[str, Any]) -> int:
    keys = {
        canonical_key(name)
        for name in merged_names(payload)
        if picked_audit.oral_aspiration_category(name)
    }
    return len(keys)


def explicit_aspiration_support(payload: dict[str, Any]) -> bool:
    deterministic = as_dict(payload.get("deterministic_max"))
    text = json.dumps(
        {
            "cross_module_reasoning": deterministic.get("cross_module_reasoning"),
            "host_context": deterministic.get("host_context"),
            "best_available_summary": deterministic.get("best_available_summary"),
        },
        ensure_ascii=False,
    ).lower()
    positive_terms = (
        "aspiration pneumonia",
        "aspiration event",
        "lung abscess",
        "pulmonary abscess",
        "empyema",
        "necrotizing pneumonia",
        "吸入性肺炎",
        "肺膿瘍",
        "膿胸",
        "壞死性肺炎",
    )
    return any(term in text for term in positive_terms)


def evidence_axis_score(row: dict[str, Any]) -> int:
    score = 0
    if row["exact_species_sterile_support"] == "yes":
        score += 7
    elif row["exact_species_pulmonary_support"] == "yes":
        score += 5
    elif row["exact_species_hospital_support"] == "yes":
        score += 2
    if row["is_d2_d3_dominant"] == "yes":
        score += 4
    if row["reads_tier"] == "R4_very_high":
        score += 3
    elif row["reads_tier"] == "R3_high":
        score += 2
    elif row["reads_tier"] == "R2_medium":
        score += 1
    if row["clinical_category"] in {
        "core_respiratory_virus",
        "common_hospital_respiratory_pathogen",
        "enterobacterales_or_hospital_gnb",
        "nonfermenter_hospital_gnb",
        "high-consequence/opportunistic respiratory pathogen",
    }:
        score += 2
    if row["low_specificity_flag"] == "yes":
        score -= 2
    return score


def low_specificity(category: str) -> bool:
    return category in {
        "candida_or_generic_yeast",
        "herpes_reactivation_virus",
        "low_actionability_virus",
        "environmental_water_soil_low_specificity_gnb",
        "skin_or_airway_background",
        "other_nondiphtherial_corynebacterium",
        "broad_oral_upper_airway_commensal",
        "atypical_oral_associated",
        "enterococcus_gi_urinary_colonizer_prone",
    }


def add_competing_pathogens(rows: list[dict[str, Any]]) -> None:
    for row in rows:
        current = int(row["evidence_axis_score"])
        competitors: list[dict[str, Any]] = []
        for other in rows:
            if other is row:
                continue
            if pathogen_names.representative_group_key(other["organism_name"]) == pathogen_names.representative_group_key(
                row["organism_name"]
            ):
                continue
            other_score = int(other["evidence_axis_score"])
            convergent = (
                other["exact_species_pulmonary_support"] == "yes"
                or other["exact_species_sterile_support"] == "yes"
                or other["is_d2_d3_dominant"] == "yes"
            )
            if convergent and other_score >= current + 3:
                competitors.append(other)
            elif (
                row["low_specificity_flag"] == "yes"
                and other["low_specificity_flag"] == "no"
                and convergent
                and other_score >= current + 2
            ):
                competitors.append(other)
        competitors.sort(key=lambda value: int(value["evidence_axis_score"]), reverse=True)
        if competitors:
            best = competitors[0]
            row["stronger_competing_pathogen"] = best["organism_name"]
            row["stronger_competing_reason"] = (
                f"score {best['evidence_axis_score']} vs {current}; competitor has direct pulmonary/sterile support or D2/D3"
            )
        else:
            row["stronger_competing_pathogen"] = ""
            row["stronger_competing_reason"] = ""


def category_recommendation(row: dict[str, Any]) -> tuple[str, str, str]:
    category = row["clinical_category"]
    high_reads = row["is_r3_r4"] == "yes"
    dominant = row["is_d2_d3_dominant"] == "yes"
    exact = row["exact_species_hospital_support"] == "yes"
    pulmonary = row["exact_species_pulmonary_support"] == "yes"
    sterile = row["exact_species_sterile_support"] == "yes"
    nonpulmonary_only = row["exact_species_nonpulmonary_only"] == "yes"
    competitor = bool(row["stronger_competing_pathogen"])
    modules = set(str(row["exact_species_hospital_modules"]).split("; "))
    molecular = bool(modules & MOLECULAR_MODULES)
    cluster = int(row["oral_aspiration_cluster_count"])
    syndrome_support = row["explicit_aspiration_syndrome_support"] == "yes"
    vulnerable = row["host_vulnerability_tier"] in {"V1", "V2", "V3"}

    if category == "candida_or_generic_yeast":
        if sterile:
            return "keep_picked", "PICK-CANDIDA-01", "Exact-species blood/sterile-site evidence supports invasive Candida context."
        if pulmonary and (high_reads or dominant):
            return "review_high_priority", "PICK-CANDIDA-02", "Respiratory Candida converges locally but lacks exact-species invasive evidence."
        if high_reads and dominant and vulnerable:
            return "review_high_priority", "PICK-CANDIDA-03", "High, dominant mNGS signal in a vulnerable host needs review but not formal Candida attribution."
        return "review_context_needed", "PICK-CANDIDA-04", "Candida without exact-species invasive evidence is better treated as colonization-versus-infection context."

    if category == "rare_opportunistic_yeast":
        if sterile or (pulmonary and (high_reads or dominant)):
            return "keep_picked", "PICK-RARE-YEAST-01", "Rare yeast has exact-species invasive or convergent lower-respiratory support."
        if exact or (high_reads and dominant and vulnerable):
            return "review_high_priority", "PICK-RARE-YEAST-02", "Rare yeast is plausible but lacks enough convergent evidence for formal picked."
        return "review_context_needed", "PICK-RARE-YEAST-03", "Weak rare-yeast evidence requires case-level review."

    if category == "herpes_reactivation_virus":
        if molecular or (dominant and row["reads_tier"] == "R4_very_high"):
            return "keep_picked", "PICK-HERPES-01", "Direct molecular support or R4 plus D2/D3 supports formal retention."
        if high_reads or dominant or pulmonary:
            return "review_high_priority", "PICK-HERPES-02", "Important herpes signal lacks direct PCR/viral-load or dominant convergence."
        return "review_context_needed", "PICK-HERPES-03", "Non-dominant herpes detection is compatible with shedding/reactivation."

    if category == "core_respiratory_virus":
        if molecular or dominant or high_reads:
            return "keep_picked", "PICK-CORE-VIRUS-01", "Core respiratory virus has direct molecular, dominant, or high-burden support."
        return "review_high_priority", "PICK-CORE-VIRUS-02", "Core respiratory virus is clinically actionable but local evidence is incomplete."

    if category == "low_actionability_virus":
        return "review_context_needed", "PICK-LOW-ACTION-VIRUS-01", "Low-actionability viral DNA should not be a formal pulmonary pathogen without exceptional evidence."

    if category == "environmental_water_soil_low_specificity_gnb":
        if sterile or (pulmonary and (high_reads or dominant)):
            return "keep_picked", "PICK-ENV-GNB-01", "Environmental GNB has exact-species invasive or convergent lower-respiratory support."
        if exact or dominant or row["reads_tier"] == "R4_very_high":
            return "review_high_priority", "PICK-ENV-GNB-02", "Strong signal warrants review, but pulmonary causality remains uncertain."
        return "review_context_needed", "PICK-ENV-GNB-03", "Environmental GNB lacks direct or dominant pulmonary convergence."

    if category == "corynebacterium_striatum":
        if sterile or (pulmonary and (high_reads or dominant)):
            return "keep_picked", "PICK-CSTR-01", "C. striatum has exact-species invasive or convergent lower-respiratory support."
        if pulmonary or dominant or row["reads_tier"] == "R4_very_high":
            return "review_high_priority", "PICK-CSTR-02", "High-burden or directly supported C. striatum needs priority review."
        return "review_context_needed", "PICK-CSTR-03", "C. striatum lacks enough convergence for formal attribution."

    if category in {"skin_or_airway_background", "other_nondiphtherial_corynebacterium"}:
        if sterile:
            return "keep_picked", "PICK-SKIN-AIRWAY-01", "Exact-species sterile-site evidence supports formal retention."
        if pulmonary and (high_reads or dominant):
            return "review_high_priority", "PICK-SKIN-AIRWAY-02", "Lower-respiratory convergence warrants priority review but remains colonization-prone."
        return "review_context_needed", "PICK-SKIN-AIRWAY-03", "Skin/airway background organism lacks invasive or convergent pulmonary evidence."

    if category == "strict_aspiration_anaerobe":
        if sterile or (pulmonary and (high_reads or dominant)):
            return "keep_picked", "PICK-ANAEROBE-01", "Anaerobe has exact-species invasive or convergent lower-respiratory support."
        if (cluster >= 2 and (high_reads or dominant)) or syndrome_support:
            return "review_high_priority", "PICK-ANAEROBE-02", "Anaerobe participates in a strong aspiration/anaerobe syndrome pattern."
        return "review_context_needed", "PICK-ANAEROBE-03", "Single or weak anaerobe signal lacks syndrome-level convergence."

    if category in {"broad_oral_upper_airway_commensal", "atypical_oral_associated"}:
        if sterile:
            return "keep_picked", "PICK-ORAL-01", "Exact-species sterile-site evidence overcomes usual oral-flora low specificity."
        if pulmonary and (high_reads or dominant) and (cluster >= 2 or syndrome_support):
            return "review_high_priority", "PICK-ORAL-02", "Oral organism has lower-respiratory and syndrome-pattern convergence."
        return "review_context_needed", "PICK-ORAL-03", "Oral/upper-airway organism is not specific enough for formal picked alone."

    if category in {"mold/opportunistic fungus", "other_fungus"}:
        if sterile or (exact and (high_reads or dominant)):
            return "keep_picked", "PICK-MOLD-01", "Fungal signal has exact-species laboratory convergence."
        if exact or (high_reads and dominant):
            return "review_high_priority", "PICK-MOLD-02", "Fungal signal is plausible but needs invasive-disease adjudication."
        return "review_context_needed", "PICK-MOLD-03", "Weak mNGS-only fungal signal is insufficient for formal picked."

    if category == "high-consequence/opportunistic respiratory pathogen":
        if canonical_key(row["organism_name"]) == "pneumocystisjirovecii" and vulnerable:
            return "keep_picked", "PICK-HIGH-CONSEQUENCE-01", "PJP is retained in a vulnerable host despite possible low reads."
        if sterile or pulmonary or molecular or dominant or high_reads:
            return "keep_picked", "PICK-HIGH-CONSEQUENCE-02", "High-consequence pathogen has at least one strong local evidence axis."
        return "review_high_priority", "PICK-HIGH-CONSEQUENCE-03", "Trace high-consequence pathogen should not be discarded but lacks formal evidence."

    if category == "enterococcus_gi_urinary_colonizer_prone":
        if sterile or pulmonary:
            return "keep_picked", "PICK-ENTEROCOCCUS-01", "Exact-species sterile or lower-respiratory evidence supports retention."
        if nonpulmonary_only:
            return "review_context_needed", "PICK-ENTEROCOCCUS-02", "Only urinary/GI evidence is not direct support for pulmonary causality."
        if high_reads and dominant:
            return "review_high_priority", "PICK-ENTEROCOCCUS-03", "Strong mNGS convergence warrants review despite colonizer-prone ecology."
        return "review_context_needed", "PICK-ENTEROCOCCUS-04", "Enterococcus lacks pulmonary or invasive exact-species support."

    if category in {
        "enterobacterales_or_hospital_gnb",
        "nonfermenter_hospital_gnb",
        "common_hospital_respiratory_pathogen",
    }:
        if sterile or pulmonary or dominant:
            return "keep_picked", "PICK-HOSPITAL-GNB-01", "Typical hospital pathogen has direct or dominant convergence."
        if high_reads and not competitor:
            return "keep_picked", "PICK-HOSPITAL-GNB-02", "R3/R4 mNGS signal is retained when no clearly stronger convergent pathogen exists."
        if high_reads or exact:
            return "review_high_priority", "PICK-HOSPITAL-GNB-03", "Plausible hospital pathogen is strong but not the best-converged candidate."
        return "review_context_needed", "PICK-HOSPITAL-GNB-04", "Weak unsupported hospital-pathogen signal needs review rather than formal attribution."

    if sterile or pulmonary or dominant:
        return "keep_picked", "PICK-OTHER-01", "Exact invasive/pulmonary or dominant evidence supports formal retention."
    if high_reads and not competitor:
        return "review_high_priority", "PICK-OTHER-02", "High mNGS burden is relevant, but category-specific pulmonary evidence is incomplete."
    return "review_context_needed", "PICK-OTHER-03", "Unclassified organism lacks direct or dominant pulmonary convergence."


def exact_match_rows(rows: list[dict[str, Any]], answers: list[compare.NameEntry]) -> None:
    unmatched = set(range(len(answers)))
    for row in rows:
        row["strict_species_answer_hit"] = "no"
        row["strict_species_matched_answer"] = ""
        key = canonical_key(row["organism_name"])
        for index in list(unmatched):
            if key and key == answers[index].canonical:
                unmatched.remove(index)
                row["strict_species_answer_hit"] = "yes"
                row["strict_species_matched_answer"] = answers[index].raw
                break


def relaxed_match_rows(rows: list[dict[str, Any]], answers: list[compare.NameEntry]) -> None:
    entries: list[picked_audit.NameEntry] = []
    for index, row in enumerate(rows):
        entry = picked_audit.name_entry(f"row:{index}", row["organism_name"])
        if entry:
            entries.append(entry)
    matched, _ = picked_audit.pairwise_match(entries, answers)
    for index, row in enumerate(rows):
        info = matched.get(f"row:{index}")
        row["baseline_answer_hit"] = "yes" if info else "no"
        row["baseline_match_type"] = info["match_type"] if info else ""
        row["baseline_matched_answer"] = info["answer"] if info else ""


def ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def metrics(rows_by_patient: dict[str, list[dict[str, Any]]], answers: dict[str, list[compare.NameEntry]]) -> dict[str, Any]:
    answer_total = sum(len(value) for value in answers.values())
    all_rows = [row for rows in rows_by_patient.values() for row in rows]
    keep_rows = [row for row in all_rows if row["suggested_shadow_tier"] == "keep_picked"]
    keep_high_rows = [
        row
        for row in all_rows
        if row["suggested_shadow_tier"] in {"keep_picked", "review_high_priority"}
    ]
    guarded_keep_rows = [row for row in all_rows if row["guarded_shadow_tier"] == "keep_picked"]
    guarded_keep_high_rows = [
        row
        for row in all_rows
        if row["guarded_shadow_tier"] in {"keep_picked", "review_high_priority"}
    ]

    def count_hits(selected: list[dict[str, Any]], strict: bool = False) -> int:
        selected_ids = {id(row) for row in selected}
        total = 0
        for pid, patient_rows in rows_by_patient.items():
            retained = [row for row in patient_rows if id(row) in selected_ids]
            if strict:
                total += sum(1 for row in retained if row["strict_species_answer_hit"] == "yes")
            else:
                entries = []
                for index, row in enumerate(retained):
                    entry = picked_audit.name_entry(f"{pid}:{index}", row["organism_name"])
                    if entry:
                        entries.append(entry)
                matched, _ = picked_audit.pairwise_match(entries, answers.get(pid, []))
                total += len(matched)
        return total

    baseline_hits = count_hits(all_rows)
    keep_hits = count_hits(keep_rows)
    keep_high_hits = count_hits(keep_high_rows)
    guarded_keep_hits = count_hits(guarded_keep_rows)
    guarded_keep_high_hits = count_hits(guarded_keep_high_rows)
    strict_baseline_hits = sum(1 for row in all_rows if row["strict_species_answer_hit"] == "yes")
    strict_keep_hits = sum(1 for row in keep_rows if row["strict_species_answer_hit"] == "yes")
    return {
        "answer_organism_count": answer_total,
        "baseline_picked": {
            "output": len(all_rows),
            "matched": baseline_hits,
            "precision": ratio(baseline_hits, len(all_rows)),
            "recall": ratio(baseline_hits, answer_total),
        },
        "shadow_keep_picked_only": {
            "output": len(keep_rows),
            "matched": keep_hits,
            "precision": ratio(keep_hits, len(keep_rows)),
            "recall": ratio(keep_hits, answer_total),
        },
        "shadow_keep_picked_plus_high": {
            "output": len(keep_high_rows),
            "matched": keep_high_hits,
            "precision": ratio(keep_high_hits, len(keep_high_rows)),
            "recall": ratio(keep_high_hits, answer_total),
        },
        "guarded_shadow_keep_picked_only": {
            "output": len(guarded_keep_rows),
            "matched": guarded_keep_hits,
            "precision": ratio(guarded_keep_hits, len(guarded_keep_rows)),
            "recall": ratio(guarded_keep_hits, answer_total),
        },
        "guarded_shadow_keep_picked_plus_high": {
            "output": len(guarded_keep_high_rows),
            "matched": guarded_keep_high_hits,
            "precision": ratio(guarded_keep_high_hits, len(guarded_keep_high_rows)),
            "recall": ratio(guarded_keep_high_hits, answer_total),
        },
        "strict_species_baseline": {
            "output": len(all_rows),
            "matched": strict_baseline_hits,
            "precision": ratio(strict_baseline_hits, len(all_rows)),
            "recall": ratio(strict_baseline_hits, answer_total),
        },
        "strict_species_shadow_keep": {
            "output": len(keep_rows),
            "matched": strict_keep_hits,
            "precision": ratio(strict_keep_hits, len(keep_rows)),
            "recall": ratio(strict_keep_hits, answer_total),
        },
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    answers, raw_answers, answer_column, answer_row_count = compare.load_answers(
        args.answer_csv, args.answer_column
    )
    rows_by_patient: dict[str, list[dict[str, Any]]] = {}
    missing_patients: list[str] = []
    ranked_cache: dict[Path, dict[str, Any]] = {}

    for patient_dir in sorted(args.patient_root.glob("NGS_patient_*_json"), key=compare.patient_sort_key):
        pid = compare.patient_id(patient_dir)
        path = patient_dir / "summary_outputs" / f"NGS_patient_{pid}_{args.merged_suffix}.json"
        if not path.is_file():
            missing_patients.append(pid)
            continue
        payload = as_dict(read_json(path))
        deterministic = as_dict(payload.get("deterministic_max"))
        best = as_dict(deterministic.get("best_available_summary"))
        picks = [item for item in as_list(best.get("picked_pathogens")) if isinstance(item, dict)]
        candidates = [item for item in as_list(deterministic.get("pathogen_candidates")) if isinstance(item, dict)]
        final_summary = load_final_summary(payload)
        hospital = hospital_index(final_summary)
        ranked_patient = load_ranked_patient(payload, pid, ranked_cache)
        cluster_count = oral_cluster_count(payload)
        aspiration_support = explicit_aspiration_support(payload)
        patient_rows: list[dict[str, Any]] = []

        for picked in picks:
            name = clean_name(picked.get("organism_name") or picked.get("name"))
            if not name:
                continue
            candidate = candidate_for_pick(picked, candidates)
            classification = classification_for(picked, candidate)
            category = clinical_category(name, classification)
            profile = hospital_profile(hospital.get(canonical_key(name), []))
            mngs_sites, mngs_buckets = mngs_sites_for_name(ranked_patient, name)
            reads_tier = str(candidate.get("reads_tier") or picked.get("reads_tier") or "Unknown")
            dominance = str(candidate.get("dominance_tier") or picked.get("dominance_tier") or "Unknown")
            row = {
                "label": args.label,
                "patient_id": pid,
                "organism_name": name,
                "classification": classification,
                "clinical_category": category,
                "category_rule_id": "",
                "baseline_answer_hit": "no",
                "baseline_match_type": "",
                "baseline_matched_answer": "",
                "strict_species_answer_hit": "no",
                "strict_species_matched_answer": "",
                "basis_level": picked.get("basis_level") or candidate.get("integrated_causative_level") or "",
                "picked_role": picked.get("picked_role") or "",
                "evidence_source": picked.get("evidence_source") or candidate.get("evidence_source") or "mNGS_ranked",
                "rank_priority": picked.get("rank_priority") or candidate.get("rank_priority") or "",
                "reads": picked.get("reads") if picked.get("reads") is not None else candidate.get("reads", ""),
                "reads_tier": reads_tier,
                "reads_percentile": candidate.get("reads_percentile") if candidate.get("reads_percentile") is not None else "",
                "dominance_tier": dominance,
                "is_d2_d3_dominant": "yes" if dominance in DOMINANT_TIERS else "no",
                "is_r3_r4": "yes" if reads_tier in HIGH_READ_TIERS else "no",
                "mngs_specimen_sites": "; ".join(mngs_sites),
                "mngs_specimen_buckets": "; ".join(mngs_buckets),
                "exact_species_hospital_support": "yes" if profile["direct"] else "no",
                "exact_species_pulmonary_support": "yes" if profile["pulmonary"] else "no",
                "exact_species_sterile_support": "yes" if profile["sterile"] else "no",
                "exact_species_nonpulmonary_only": "yes" if profile["nonpulmonary_only"] else "no",
                "exact_species_hospital_modules": "; ".join(profile["modules"]),
                "exact_species_hospital_sites": "; ".join(profile["sites"]),
                "exact_species_hospital_buckets": "; ".join(profile["buckets"]),
                "related_genus_or_group_support": "yes" if related_support(picked, candidate) else "no",
                "host_vulnerability_tier": host_tier(payload),
                "oral_aspiration_cluster_count": cluster_count,
                "explicit_aspiration_syndrome_support": "yes" if aspiration_support else "no",
                "low_specificity_flag": "yes" if low_specificity(category) else "no",
                "stronger_competing_pathogen": "",
                "stronger_competing_reason": "",
                "evidence_axis_score": 0,
                "suggested_shadow_tier": "",
                "suggested_shadow_reason": "",
                "guarded_shadow_tier": "",
                "guarded_rule_id": "",
                "guarded_reason": "",
            }
            row["evidence_axis_score"] = evidence_axis_score(row)
            patient_rows.append(row)

        add_competing_pathogens(patient_rows)
        for row in patient_rows:
            tier, rule_id, reason = category_recommendation(row)
            row["suggested_shadow_tier"] = tier
            row["category_rule_id"] = rule_id
            row["suggested_shadow_reason"] = reason
            guarded_tier, guarded_rule_id, guarded_reason = guarded_recommendation(row)
            row["guarded_shadow_tier"] = guarded_tier
            row["guarded_rule_id"] = guarded_rule_id
            row["guarded_reason"] = guarded_reason
        relaxed_match_rows(patient_rows, answers.get(pid, []))
        exact_match_rows(patient_rows, answers.get(pid, []))
        rows_by_patient[pid] = patient_rows

    all_rows = [row for rows in rows_by_patient.values() for row in rows]
    summary = {
        "label": args.label,
        "patient_root": str(args.patient_root),
        "merged_suffix": args.merged_suffix,
        "answer_csv": str(args.answer_csv),
        "answer_column": answer_column,
        "answer_csv_row_count": answer_row_count,
        "answer_patient_count": len(answers),
        "missing_patients": missing_patients,
        "recommendation_is_answer_blind": True,
        "metrics": metrics(rows_by_patient, answers),
        "recommendation_counts": dict(Counter(row["suggested_shadow_tier"] for row in all_rows)),
        "guarded_recommendation_counts": dict(Counter(row["guarded_shadow_tier"] for row in all_rows)),
        "recommendation_answer_hits": {
            tier: {
                "matched": sum(
                    1
                    for row in all_rows
                    if row["suggested_shadow_tier"] == tier and row["baseline_answer_hit"] == "yes"
                ),
                "output": sum(1 for row in all_rows if row["suggested_shadow_tier"] == tier),
            }
            for tier in ("keep_picked", "review_high_priority", "review_context_needed")
        },
        "category_counts": dict(Counter(row["clinical_category"] for row in all_rows)),
        "category_answer_hits": {
            category: {
                "matched": sum(
                    1
                    for row in all_rows
                    if row["clinical_category"] == category and row["baseline_answer_hit"] == "yes"
                ),
                "output": sum(1 for row in all_rows if row["clinical_category"] == category),
            }
            for category in sorted({row["clinical_category"] for row in all_rows})
        },
        "suggested_demotions": [
            {
                key: row[key]
                for key in (
                    "patient_id",
                    "organism_name",
                    "clinical_category",
                    "baseline_answer_hit",
                    "baseline_matched_answer",
                    "suggested_shadow_tier",
                    "category_rule_id",
                    "suggested_shadow_reason",
                    "stronger_competing_pathogen",
                )
            }
            for row in all_rows
            if row["suggested_shadow_tier"] != "keep_picked"
        ],
        "guarded_demotions": [
            {
                key: row[key]
                for key in (
                    "patient_id",
                    "organism_name",
                    "clinical_category",
                    "baseline_answer_hit",
                    "baseline_matched_answer",
                    "guarded_shadow_tier",
                    "guarded_rule_id",
                    "guarded_reason",
                    "stronger_competing_pathogen",
                )
            }
            for row in all_rows
            if row["guarded_shadow_tier"] != "keep_picked"
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "picked_evidence_axis_audit.csv", all_rows, FIELDS)
    write_json(args.output_dir / "picked_evidence_axis_summary.json", summary)
    write_json(args.output_dir / "picked_evidence_axis_rows.json", all_rows)
    print(json.dumps(summary["metrics"], ensure_ascii=False, indent=2))
    print(f"recommendations={summary['recommendation_counts']}")
    print(f"guarded_recommendations={summary['guarded_recommendation_counts']}")
    print(f"missing_patients={len(missing_patients)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

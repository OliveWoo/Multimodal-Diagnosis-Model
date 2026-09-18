"""Shadow audit for false positives inside deterministic picked_pathogens.

This report is deliberately read-only. It flattens the current picked_pathogens,
matches them against the answer CSV, and simulates a small set of conservative
candidate guardrails before any deterministic max rule is changed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import pathogen_normalization as pathogen_names  # noqa: E402


DEFAULT_MAX_SUFFIX = "mNGS_max_deterministic_pulmonary_opt_guardrail_20260804_aliasclean_20260805_direct_triage"
NON_ANSWER = {"", "-", "na", "n/a", "none", "negative", "no pathogen", "nopatogen", "notinfectioncase"}
NON_HOST_SUPPORT = {"culture", "filmarray_gmtest", "molecular_microbiology"}
LOW_SIGNAL_READ_TIERS = {"R0_trace", "R1_low"}
NON_DOMINANT_TIERS = {"", "D0_not_top", "D1_low", "Not_available", "Unknown"}
RESP_COMMENSAL_RULES = {"R-S4-RESP-COMMENSAL-L3", "R-S4-CONS-SKIN_FLORA_SUPPORTED_L3"}
HERPES_KEYS = {
    "humanbetaherpesvirus5",
    "cytomegalovirus",
    "cmv",
    "humangammaherpesvirus4",
    "ebv",
    "humanalphaherpesvirus1",
    "humanalphaherpesvirus2",
    "hsv1",
    "hsv2",
}
LOW_ACTIONABILITY_VIRUS_PREFIXES = (
    "humanpapillomavirus",
    "papillomavirus",
    "ttv",
    "torqueteno",
    "anellovirus",
)
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
WATER_ENVIRONMENTAL_GNB_PREFIXES = (
    "achromobacter",
    "acidovorax",
    "brevundimonas",
    "cloacibacterium",
    "comamonas",
    "flavobacterium",
    "phytobacter",
    "ralstonia",
    "pseudoxanthomonas",
    "xanthomonas",
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
    "sphingobium",
    "sphingomonas",
    "undibacterium",
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
)
ORAL_UPPER_AIRWAY_COMMENSAL_PREFIXES = (
    "gemella",
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
ORAL_UPPER_AIRWAY_COMMENSAL_EXACT_NAMES = {
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
HIGH_CONSEQUENCE_KEYS = {
    "pneumocystisjirovecii",
    "mycobacteriumtuberculosis",
    "legionellapneumophila",
    "nocardiathailandica",
    "burkholderiacenocepacia",
    "severeacuterespiratorysyndromerelatedcoronavirus",
}
COMMON_HOSPITAL_RESPIRATORY_KEYS = {
    "acinetobacterbaumannii",
    "klebsiellapneumoniae",
    "klebsiellavariicola",
    "pseudomonasaeruginosa",
    "serratia marcescens".replace(" ", ""),
    "stenotrophomonasmaltophilia",
    "escherichiacoli",
    "haemophilusinfluenzae",
    "staphylococcusaureus",
}
SHADOW_PROFILES = {
    "conservative_low_signal_and_generic": {
        "description": "Only generic unidentified labels plus Level 3 mNGS-only low-read non-dominant signals.",
        "rules": [
            "generic_unidentified_hospital_label_shadow_omit",
            "mngs_only_level3_low_signal_shadow_omit",
        ],
    },
    "plus_same_representative_group": {
        "description": "Conservative profile plus one representative per alias/genus group.",
        "rules": [
            "generic_unidentified_hospital_label_shadow_omit",
            "mngs_only_level3_low_signal_shadow_omit",
            "same_representative_group_shadow_omit",
        ],
    },
    "plus_candida_yeast_nondominant": {
        "description": "Conservative profile plus non-dominant Candida/yeast Level 3 mNGS-only pruning.",
        "rules": [
            "generic_unidentified_hospital_label_shadow_omit",
            "mngs_only_level3_low_signal_shadow_omit",
            "candida_yeast_mngs_only_nondominant_shadow_omit",
        ],
    },
    "full_initial_shadow": {
        "description": "All currently flagged shadow rules.",
        "rules": [
            "same_representative_group_shadow_omit",
            "generic_unidentified_hospital_label_shadow_omit",
            "mngs_only_level3_low_signal_shadow_omit",
            "candida_yeast_mngs_only_nondominant_shadow_omit",
            "herpes_reactivation_low_signal_shadow_omit",
            "resp_commensal_level3_nondominant_shadow_omit",
            "low_read_mold_watchlist_shadow_omit",
        ],
    },
    "picked_safe_low_specificity_shadow": {
        "description": "Safest picked audit: low-actionability viruses and water/environmental low pulmonary-specificity organisms only.",
        "rules": [
            "low_actionability_virus_shadow_omit",
            "water_environmental_low_specificity_shadow_omit",
        ],
    },
    "picked_safe_plus_nonprotected_r1": {
        "description": "Safe profile plus mNGS-only R0/R1, D0/D1, host-only, non-protected organisms.",
        "rules": [
            "low_actionability_virus_shadow_omit",
            "water_environmental_low_specificity_shadow_omit",
            "mngs_r0_r1_host_only_nondominant_nonprotected_shadow_omit",
        ],
    },
    "picked_safe_plus_weak_mold": {
        "description": "Safe profile plus weak mold/fungal host-only non-dominant signals.",
        "rules": [
            "low_actionability_virus_shadow_omit",
            "water_environmental_low_specificity_shadow_omit",
            "mngs_r0_r1_host_only_nondominant_nonprotected_shadow_omit",
            "weak_mold_fungal_host_only_nondominant_shadow_omit",
        ],
    },
    "picked_mngs_only_nondominant_non_typical_shadow": {
        "description": "Exploratory picked FP profile: mNGS-only, no same-organism hospital support, D0/D1 non-dominant, and not a typical strong respiratory pathogen.",
        "rules": [
            "mngs_only_no_same_support_nondominant_non_typical_shadow_omit",
        ],
    },
    "picked_mngs_only_nondominant_non_typical_guarded_shadow": {
        "description": "Guarded version of the mNGS-only non-dominant profile; keeps deferred GI/urinary colonizer-prone organisms and Acinetobacter for separate review.",
        "rules": [
            "mngs_only_no_same_support_nondominant_non_typical_guarded_shadow_omit",
        ],
    },
    "picked_precision_070_non_typical_guarded_combo_shadow": {
        "description": "Current precision-target candidate: existing safe picked pruning plus non-core respiratory virus context demotion and guarded non-typical mNGS-only non-dominant pruning.",
        "rules": [
            "low_actionability_virus_shadow_omit",
            "water_environmental_low_specificity_shadow_omit",
            "mngs_r0_r1_host_only_nondominant_nonprotected_shadow_omit",
            "weak_mold_fungal_host_only_nondominant_shadow_omit",
            "non_core_respiratory_virus_context_shadow_omit",
            "mngs_only_no_same_support_nondominant_non_typical_guarded_shadow_omit",
        ],
    },
    "picked_target_precision_resp_virus_context_exploratory": {
        "description": "Exploratory target-precision profile: safe weak-signal profile plus non-SARS/non-influenza/non-herpes respiratory viruses moved from picked to context.",
        "rules": [
            "low_actionability_virus_shadow_omit",
            "water_environmental_low_specificity_shadow_omit",
            "mngs_r0_r1_host_only_nondominant_nonprotected_shadow_omit",
            "weak_mold_fungal_host_only_nondominant_shadow_omit",
            "non_core_respiratory_virus_context_shadow_omit",
        ],
    },
    "picked_candida_yeast_no_invasive_exploratory": {
        "description": "Exploratory only: removes Candida/generic yeast picked items without invasive Candida evidence.",
        "rules": [
            "candida_yeast_no_invasive_picked_shadow_omit",
        ],
    },
    "picked_enterococcus_host_only_exploratory": {
        "description": "Exploratory only: removes Enterococcus picked items when they are mNGS host-only, non-dominant, and lack non-host support.",
        "rules": [
            "enterococcus_mngs_host_only_nondominant_shadow_omit",
        ],
    },
    "picked_hospital_only_l2_exploratory": {
        "description": "Exploratory only: removes Level 2 hospital-only picked items without mNGS ranked support.",
        "rules": [
            "hospital_only_l2_no_mngs_shadow_omit",
        ],
    },
    "picked_aggressive_precision_shadow": {
        "description": "Exploratory aggressive profile combining all picked-FP pruning candidates; inspect hit loss before using.",
        "rules": [
            "low_actionability_virus_shadow_omit",
            "water_environmental_low_specificity_shadow_omit",
            "mngs_r0_r1_host_only_nondominant_nonprotected_shadow_omit",
            "mngs_r0_r1_host_only_nondominant_all_shadow_omit",
            "mngs_r2_host_only_nondominant_nonprotected_shadow_omit",
            "weak_mold_fungal_host_only_nondominant_shadow_omit",
            "non_core_respiratory_virus_context_shadow_omit",
            "candida_yeast_no_invasive_picked_shadow_omit",
            "enterococcus_mngs_host_only_nondominant_shadow_omit",
            "hospital_only_l2_no_mngs_shadow_omit",
        ],
    },
}

CSV_FIELDS = (
    "dataset",
    "patient_id",
    "organism_name",
    "canonical_key",
    "genus",
    "representative_group",
    "answer_patient",
    "baseline_answer_hit",
    "baseline_match_type",
    "baseline_matched_answer",
    "after_shadow_answer_hit",
    "after_shadow_match_type",
    "after_shadow_matched_answer",
    "shadow_action",
    "primary_shadow_rule",
    "all_shadow_rules",
    "shadow_reason",
    "basis_level",
    "mngs_signal_tier",
    "picked_role",
    "evidence_source",
    "rank_priority",
    "reads",
    "reads_tier",
    "reads_percentile",
    "dominance_tier",
    "specimen_class",
    "source_category",
    "clinical_ecology_group",
    "support_modules",
    "related_representative_hospital_support",
    "guardrail_rule",
    "level_rule",
    "applied_rules",
    "is_protected_pathogen",
    "is_likely_colonizer_or_background",
    "selection_mode",
)


@dataclass(frozen=True)
class NameEntry:
    row_id: str
    raw: str
    canonical: str
    genus: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a picked_pathogens false-positive shadow audit."
    )
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--max-suffix", default=DEFAULT_MAX_SUFFIX)
    parser.add_argument("--answer-csv", type=Path, required=True)
    parser.add_argument(
        "--answer-column",
        default="auto",
        help="Answer column in --answer-csv. Use auto to prefer answer, answers, then hospital_answer.",
    )
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--patient-id", action="append", default=[])
    parser.add_argument("--quiet", action="store_true", help="Write reports without printing full metrics JSON.")
    return parser.parse_args(argv)


def resolve_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})


def patient_id_from_dir(patient_dir: Path) -> str:
    text = patient_dir.name
    if text.startswith("NGS_patient_"):
        text = text.removeprefix("NGS_patient_")
    if text.endswith("_json"):
        text = text.removesuffix("_json")
    return text


def patient_sort_key(patient_dir: Path) -> tuple[int, str]:
    pid = patient_id_from_dir(patient_dir)
    return (int(pid), pid) if pid.isdigit() else (10**9, pid)


def patient_dirs(patient_root: Path, selected_ids: set[str]) -> list[Path]:
    dirs = [
        path
        for path in patient_root.iterdir()
        if path.is_dir() and path.name.startswith("NGS_patient_")
    ]
    if selected_ids:
        dirs = [path for path in dirs if patient_id_from_dir(path) in selected_ids]
    return sorted(dirs, key=patient_sort_key)


def output_path(patient_dir: Path, suffix: str) -> Path:
    pid = patient_id_from_dir(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{pid}_{suffix}.json"


def deterministic_payload(payload: Any) -> dict[str, Any]:
    if not isinstance(payload, dict):
        return {}
    nested = payload.get("deterministic_max")
    if isinstance(nested, dict):
        return nested
    return payload


def clean_text(value: Any) -> str:
    return pathogen_names.clean_display_text(value)


def canonical_key(value: Any) -> str:
    return pathogen_names.canonical_key(value)


def genus_name(value: Any) -> str:
    return pathogen_names.genus_name(value)


def name_entry(row_id: str, value: Any) -> NameEntry | None:
    raw = clean_text(value)
    if not raw:
        return None
    key = canonical_key(raw)
    if not key or key.isdigit():
        return None
    return NameEntry(row_id=row_id, raw=raw, canonical=key, genus=genus_name(raw))


def split_names(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text or text.strip().lower() in NON_ANSWER:
        return []
    names: list[str] = []
    for part in re.split(r"[;\n\r\uFF1B\u3001]+", text):
        name = clean_text(part)
        if not name or name.lower() in NON_ANSWER:
            continue
        names.append(name)
    return names


def resolve_answer_column(fieldnames: Sequence[str] | None, requested: str) -> str:
    fields = list(fieldnames or [])
    if requested != "auto":
        if requested not in fields:
            raise ValueError(f"Answer column {requested!r} not found in {fields}")
        return requested
    for candidate in ("answer", "answers", "hospital_answer"):
        if candidate in fields:
            return candidate
    raise ValueError(f"Cannot auto-detect answer column from CSV fields: {fields}")


def load_answers(path: Path, answer_column: str) -> tuple[dict[str, list[NameEntry]], dict[str, str], str, int]:
    answers: dict[str, list[NameEntry]] = {}
    raw_answers: dict[str, str] = {}
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        column = resolve_answer_column(reader.fieldnames, answer_column)
        row_count = 0
        for row in reader:
            row_count += 1
            pid = str(row.get("patient_id") or "").strip()
            if not pid:
                continue
            raw = str(row.get(column) or "").strip()
            entries: list[NameEntry] = []
            for idx, value in enumerate(split_names(raw)):
                entry = name_entry(f"{pid}:answer:{idx}", value)
                if entry:
                    entries.append(entry)
            raw_answers[pid] = raw
            if entries:
                answers[pid] = entries
        return answers, raw_answers, column, row_count


def names_match(output: NameEntry, answer: NameEntry) -> str:
    return pathogen_names.name_match_type(output.raw, answer.raw)


def pairwise_match(outputs: list[NameEntry], answers: list[NameEntry]) -> tuple[dict[str, dict[str, str]], list[str]]:
    unmatched = set(range(len(answers)))
    matched: dict[str, dict[str, str]] = {}
    for output in outputs:
        best_idx = None
        best_type = ""
        for idx in list(unmatched):
            match_type = names_match(output, answers[idx])
            if not match_type:
                continue
            if match_type == "exact_or_alias":
                best_idx = idx
                best_type = match_type
                break
            if best_idx is None or (
                match_type == "approved_group_member_match" and best_type == "genus_relaxed"
            ):
                best_idx = idx
                best_type = match_type
        if best_idx is None:
            continue
        unmatched.remove(best_idx)
        matched[output.row_id] = {
            "answer": answers[best_idx].raw,
            "match_type": best_type,
        }
    return matched, [answers[idx].raw for idx in sorted(unmatched)]


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def numeric(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def level_rank(value: Any) -> int:
    match = re.search(r"(\d+)", str(value or ""))
    return int(match.group(1)) if match else 99


def tier_score(value: Any) -> int:
    text = str(value or "")
    if text == "R4_very_high":
        return 40
    if text == "R3_high":
        return 30
    if text == "R2_medium":
        return 20
    if text == "R1_low":
        return 10
    if text == "R0_trace":
        return 0
    return -10


def dominance_score(value: Any) -> int:
    text = str(value or "")
    if text == "D3_dominant":
        return 80
    if text == "D2_moderate":
        return 60
    if text == "D1_low":
        return 20
    if text == "D0_not_top":
        return 0
    return -10


def support_modules_from(candidate: dict[str, Any], picked: dict[str, Any]) -> list[str]:
    values: list[str] = []
    values.extend(str(value) for value in as_list(picked.get("support_modules")))
    key_evidence = candidate.get("key_evidence") if isinstance(candidate.get("key_evidence"), dict) else {}
    values.extend(str(value) for value in as_list(key_evidence.get("support_modules")))
    module_support = candidate.get("module_support_summary")
    if isinstance(module_support, dict):
        for module, status in module_support.items():
            if status == "Support":
                values.append(str(module))
    output: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = value.strip()
        if text and text not in seen:
            output.append(text)
            seen.add(text)
    return output


def has_non_host_support(row: dict[str, Any]) -> bool:
    modules = {str(value).strip() for value in row.get("_support_modules_list", [])}
    return bool(modules & NON_HOST_SUPPORT)


def support_is_host_only(row: dict[str, Any]) -> bool:
    modules = {str(value).strip() for value in row.get("_support_modules_list", []) if str(value).strip()}
    return not modules or modules <= {"host"}


def has_related_representative_context(row: dict[str, Any]) -> bool:
    value = row.get("_related_representative_hospital_support")
    return isinstance(value, list) and bool(value)


def applied_rules_from(candidate: dict[str, Any]) -> list[str]:
    return [str(value) for value in as_list(candidate.get("applied_rules")) if str(value)]


def candidate_sort_score(row: dict[str, Any]) -> tuple[float, ...]:
    return (
        -level_rank(row.get("basis_level")),
        100.0 if has_non_host_support(row) else 0.0,
        float(dominance_score(row.get("dominance_tier"))),
        float(tier_score(row.get("reads_tier"))),
        math.log10(numeric(row.get("reads")) + 1.0),
        pathogen_names.representative_specificity_score(row.get("organism_name")),
    )


def best_candidate_for_pick(
    picked: dict[str, Any],
    candidates: list[dict[str, Any]],
) -> dict[str, Any]:
    name = picked.get("organism_name") or picked.get("name") or ""
    key = canonical_key(name)
    if not key:
        return {}
    scored: list[tuple[tuple[int, int, int, int, float], dict[str, Any]]] = []
    picked_reads = numeric(picked.get("reads"), default=-1)
    picked_level = level_rank(picked.get("basis_level"))
    picked_rank = str(picked.get("rank_priority") or "")
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        if canonical_key(candidate.get("organism_name")) != key:
            continue
        score = (
            1 if clean_text(candidate.get("organism_name")).lower() == clean_text(name).lower() else 0,
            1 if numeric(candidate.get("reads"), default=-2) == picked_reads else 0,
            1 if level_rank(candidate.get("integrated_causative_level")) == picked_level else 0,
            1 if str(candidate.get("rank_priority") or "") == picked_rank else 0,
            numeric(candidate.get("reads"), default=0.0),
        )
        scored.append((score, candidate))
    if not scored:
        return {}
    scored.sort(key=lambda item: item[0], reverse=True)
    return scored[0][1]


def build_picked_rows(
    *,
    label: str,
    patient_id: str,
    payload: dict[str, Any],
) -> list[dict[str, Any]]:
    deterministic = deterministic_payload(payload)
    best = deterministic.get("best_available_summary") if isinstance(deterministic, dict) else {}
    candidates = [item for item in as_list(deterministic.get("pathogen_candidates")) if isinstance(item, dict)]
    picked_items = [item for item in as_list((best or {}).get("picked_pathogens")) if isinstance(item, dict)]
    rows: list[dict[str, Any]] = []
    for index, picked in enumerate(picked_items):
        name = clean_text(picked.get("organism_name") or picked.get("name"))
        if not name:
            continue
        candidate = best_candidate_for_pick(picked, candidates)
        key_evidence = candidate.get("key_evidence") if isinstance(candidate.get("key_evidence"), dict) else {}
        modules = support_modules_from(candidate, picked)
        applied = applied_rules_from(candidate)
        related_support = key_evidence.get("related_representative_hospital_support")
        if not isinstance(related_support, list):
            related_support = picked.get("related_representative_hospital_support")
        if not isinstance(related_support, list):
            related_support = []
        row = {
            "dataset": label,
            "patient_id": patient_id,
            "row_id": f"{patient_id}:picked:{index}",
            "organism_name": name,
            "canonical_key": canonical_key(name),
            "genus": genus_name(name),
            "representative_group": pathogen_names.representative_group_key(name),
            "baseline_answer_hit": False,
            "baseline_match_type": "",
            "baseline_matched_answer": "",
            "after_shadow_answer_hit": False,
            "after_shadow_match_type": "",
            "after_shadow_matched_answer": "",
            "shadow_action": "keep",
            "primary_shadow_rule": "",
            "all_shadow_rules": "",
            "shadow_reason": "",
            "basis_level": picked.get("basis_level") or candidate.get("integrated_causative_level") or "",
            "mngs_signal_tier": picked.get("mngs_signal_tier") or candidate.get("mngs_signal_tier") or "",
            "picked_role": picked.get("picked_role") or "",
            "evidence_source": picked.get("evidence_source") or candidate.get("evidence_source") or "mNGS_ranked",
            "rank_priority": picked.get("rank_priority") or candidate.get("rank_priority") or "",
            "reads": picked.get("reads") if picked.get("reads") is not None else candidate.get("reads", ""),
            "reads_tier": candidate.get("reads_tier") or "",
            "reads_percentile": candidate.get("reads_percentile") if candidate.get("reads_percentile") is not None else "",
            "dominance_tier": candidate.get("dominance_tier") or "",
            "specimen_class": candidate.get("specimen_class") or "",
            "source_category": candidate.get("source_category") or "",
            "clinical_ecology_group": clinical_ecology_group(name),
            "support_modules": "; ".join(modules),
            "related_representative_hospital_support": json.dumps(related_support, ensure_ascii=False)
            if related_support
            else "",
            "guardrail_rule": key_evidence.get("guardrail_rule") or "",
            "level_rule": key_evidence.get("level_rule") or "",
            "applied_rules": "; ".join(applied),
            "is_protected_pathogen": candidate.get("is_protected_pathogen", ""),
            "is_likely_colonizer_or_background": candidate.get("is_likely_colonizer_or_background", ""),
            "selection_mode": (best or {}).get("selection_mode") or "",
            "_support_modules_list": modules,
            "_applied_rules_list": applied,
            "_related_representative_hospital_support": related_support,
            "_picked_index": index,
        }
        rows.append(row)
    return rows


def is_generic_hospital_label(name: str) -> bool:
    text = name.lower()
    return bool(
        re.search(r"identification\s+to\s+follow", text)
        or re.search(r"\bg\s*\([+-]\)", text)
        or re.search(r"\bgram[-\s]?(negative|positive)\b", text)
        or re.search(r"\b(yeast|bacteria|bacilli|cocci),?\s*identification", text)
    )


def is_herpes_reactivation(row: dict[str, Any]) -> bool:
    key = str(row.get("canonical_key") or "")
    raw = pathogen_names.raw_key(row.get("organism_name"))
    return (
        key in HERPES_KEYS
        or raw in HERPES_KEYS
        or "R-S4-HERPES_REACTIVATION_GUARDRAIL" in row.get("_applied_rules_list", [])
    )


def is_candida_or_yeast(row: dict[str, Any]) -> bool:
    return pathogen_names.is_candida_or_generic_yeast(row.get("organism_name"))


def normalized_name_keys(value: Any) -> set[str]:
    keys = {
        canonical_key(value),
        pathogen_names.raw_key(value),
        re.sub(r"[^a-z0-9]+", "", str(value or "").lower()),
    }
    keys.update(pathogen_names.expanded_name_keys(value))
    return {key for key in keys if key}


def row_name_keys(row: dict[str, Any]) -> set[str]:
    keys = {
        str(row.get("canonical_key") or ""),
        pathogen_names.raw_key(row.get("organism_name")),
    }
    keys.update(normalized_name_keys(row.get("organism_name")))
    return {key for key in keys if key}


def keys_start_with(keys: set[str], prefixes: tuple[str, ...]) -> bool:
    return any(key.startswith(prefix) for key in keys for prefix in prefixes)


def is_low_actionability_virus(row: dict[str, Any]) -> bool:
    return keys_start_with(row_name_keys(row), LOW_ACTIONABILITY_VIRUS_PREFIXES)


def is_respiratory_virus(row: dict[str, Any]) -> bool:
    keys = row_name_keys(row)
    return keys_start_with(keys, RESPIRATORY_VIRUS_PREFIXES) or is_herpes_reactivation(row)


def is_core_answer_sensitive_respiratory_virus(row: dict[str, Any]) -> bool:
    keys = row_name_keys(row)
    return (
        is_herpes_reactivation(row)
        or keys_start_with(keys, ("influenza", "severeacuterespiratorysyndrome", "sarscov2"))
        or "severeacuterespiratorysyndromerelatedcoronavirus" in keys
    )


def is_non_core_respiratory_virus(row: dict[str, Any]) -> bool:
    return (
        is_respiratory_virus(row)
        and not is_low_actionability_virus(row)
        and not is_core_answer_sensitive_respiratory_virus(row)
    )


def is_water_or_low_specificity_environmental(row: dict[str, Any]) -> bool:
    keys = row_name_keys(row)
    return (
        bool(keys & WATER_ENVIRONMENTAL_GNB_EXACT_NAMES)
        or keys_start_with(keys, WATER_ENVIRONMENTAL_GNB_PREFIXES)
        or keys_start_with(keys, LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_PREFIXES)
    )


def oral_aspiration_category(value: Any) -> str:
    keys = normalized_name_keys(value)
    if keys_start_with(keys, STRICT_ASPIRATION_ANAEROBE_PREFIXES):
        return "strict_aspiration_anaerobe"
    if keys & ORAL_UPPER_AIRWAY_COMMENSAL_EXACT_NAMES or keys_start_with(
        keys, ORAL_UPPER_AIRWAY_COMMENSAL_PREFIXES
    ):
        return "broad_oral_upper_airway_commensal"
    if keys & ATYPICAL_ORAL_ASSOCIATED_EXACT_NAMES or keys_start_with(
        keys, ATYPICAL_ORAL_ASSOCIATED_PREFIXES
    ):
        return "atypical_oral_associated"
    return ""


def is_high_consequence_pathogen(row: dict[str, Any]) -> bool:
    keys = row_name_keys(row)
    return bool(keys & HIGH_CONSEQUENCE_KEYS)


def is_common_hospital_respiratory_pathogen(row: dict[str, Any]) -> bool:
    keys = row_name_keys(row)
    return bool(keys & COMMON_HOSPITAL_RESPIRATORY_KEYS)


def is_typical_strong_respiratory_pathogen(row: dict[str, Any]) -> bool:
    return (
        is_high_consequence_pathogen(row)
        or is_common_hospital_respiratory_pathogen(row)
        or is_respiratory_virus(row)
        or str(row.get("is_protected_pathogen")).lower() == "true"
    )


def is_nonprotected_for_low_signal_prune(row: dict[str, Any]) -> bool:
    return not is_typical_strong_respiratory_pathogen(row)


def is_deferred_gi_urinary_colonizer_for_shadow(row: dict[str, Any]) -> bool:
    keys = row_name_keys(row)
    return str(row.get("genus") or "").lower() == "enterococcus" or bool(
        keys & {"clostridioidesdifficile", "clostridiumdifficile"}
    )


def is_answer_sensitive_healthcare_opportunist_for_shadow(row: dict[str, Any]) -> bool:
    return str(row.get("genus") or "").lower() == "acinetobacter"


def is_enterococcus(row: dict[str, Any]) -> bool:
    return str(row.get("genus") or "").lower() == "enterococcus"


def is_mold_or_opportunistic_fungus(row: dict[str, Any]) -> bool:
    genus = str(row.get("genus") or "").lower()
    return genus in {
        "aspergillus",
        "cunninghamella",
        "talaromyces",
        "schizophyllum",
        "rhizopus",
        "mucor",
        "lichtheimia",
        "fusarium",
    }


def has_explicit_invasive_candida_context(row: dict[str, Any]) -> bool:
    text = " ".join(
        str(row.get(field) or "")
        for field in (
            "guardrail_rule",
            "level_rule",
            "applied_rules",
            "related_representative_hospital_support",
            "shadow_reason",
        )
    ).lower()
    return any(term in text for term in ("invasive", "blood", "sterile", "pleural", "abscess", "tissue", "histolog"))


def clinical_ecology_group(name: Any) -> str:
    keys = normalized_name_keys(name)
    genus = genus_name(name)
    probe = {
        "organism_name": name,
        "canonical_key": canonical_key(name),
        "genus": genus,
    }
    if pathogen_names.is_candida_or_generic_yeast(name):
        return "Candida/generic yeast"
    if is_herpes_reactivation(probe):
        return "herpes/reactivation virus"
    if keys_start_with(keys, LOW_ACTIONABILITY_VIRUS_PREFIXES):
        return "low-actionability virus"
    if keys_start_with(keys, RESPIRATORY_VIRUS_PREFIXES):
        return "respiratory virus"
    if bool(keys & HIGH_CONSEQUENCE_KEYS):
        return "high-consequence/opportunistic respiratory pathogen"
    if (
        bool(keys & WATER_ENVIRONMENTAL_GNB_EXACT_NAMES)
        or keys_start_with(keys, WATER_ENVIRONMENTAL_GNB_PREFIXES)
        or keys_start_with(keys, LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_PREFIXES)
    ):
        return "water/environmental low-pulmonary-specificity organism"
    oral_category = oral_aspiration_category(name)
    if oral_category:
        return oral_category
    if genus == "enterococcus":
        return "Enterococcus / GI-urinary colonizer-prone"
    if genus in {"aspergillus", "cunninghamella", "talaromyces", "schizophyllum", "rhizopus", "mucor", "lichtheimia", "fusarium"}:
        return "mold/opportunistic fungus"
    if bool(keys & COMMON_HOSPITAL_RESPIRATORY_KEYS):
        return "common hospital respiratory pathogen"
    return "other bacterial/uncategorized"


def add_shadow_rules(rows: list[dict[str, Any]]) -> None:
    by_patient_group: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        group = str(row.get("representative_group") or "")
        if group:
            by_patient_group[(str(row.get("patient_id")), group)].append(row)

    duplicate_remove_ids: set[str] = set()
    for group_rows in by_patient_group.values():
        if len(group_rows) <= 1:
            continue
        keep = max(group_rows, key=candidate_sort_score)
        for row in group_rows:
            if row is not keep:
                duplicate_remove_ids.add(str(row.get("row_id")))

    for row in rows:
        rules: list[tuple[str, str]] = []
        level = level_rank(row.get("basis_level"))
        reads_tier = str(row.get("reads_tier") or "")
        dominance = str(row.get("dominance_tier") or "")
        mngs_ranked = str(row.get("evidence_source") or "") == "mNGS_ranked"
        hospital_only = str(row.get("evidence_source") or "") == "hospital_only"
        host_only = support_is_host_only(row)
        non_dominant = dominance in NON_DOMINANT_TIERS

        if str(row.get("row_id")) in duplicate_remove_ids:
            rules.append(
                (
                    "same_representative_group_shadow_omit",
                    "Same representative group already has a stronger picked item; keep as related evidence only.",
                )
            )
        if is_generic_hospital_label(str(row.get("organism_name") or "")):
            rules.append(
                (
                    "generic_unidentified_hospital_label_shadow_omit",
                    "Generic hospital label is not species-level; keep as source evidence rather than formal picked organism.",
                )
            )
        if (
            level == 3
            and mngs_ranked
            and host_only
            and not has_related_representative_context(row)
            and non_dominant
            and reads_tier in LOW_SIGNAL_READ_TIERS
        ):
            rules.append(
                (
                    "mngs_only_level3_low_signal_shadow_omit",
                    "Level 3 mNGS-only host-context signal is low-read and non-dominant.",
                )
            )
        if (
            is_candida_or_yeast(row)
            and level == 3
            and mngs_ranked
            and host_only
            and non_dominant
            and reads_tier in {"R0_trace", "R1_low", "R2_medium"}
        ):
            rules.append(
                (
                    "candida_yeast_mngs_only_nondominant_shadow_omit",
                    "Candida/yeast lower-respiratory mNGS-only non-dominant signal is low-specificity for pneumonia.",
                )
            )
        if (
            is_herpes_reactivation(row)
            and level == 3
            and mngs_ranked
            and host_only
            and non_dominant
            and reads_tier in LOW_SIGNAL_READ_TIERS
        ):
            rules.append(
                (
                    "herpes_reactivation_low_signal_shadow_omit",
                    "CMV/HSV/EBV-like reactivation virus is low-read, non-dominant, and host-context only.",
                )
            )
        if (
            level == 3
            and mngs_ranked
            and host_only
            and non_dominant
            and reads_tier in {"R0_trace", "R1_low", "R2_medium"}
            and bool(set(row.get("_applied_rules_list", [])) & RESP_COMMENSAL_RULES)
        ):
            rules.append(
                (
                    "resp_commensal_level3_nondominant_shadow_omit",
                    "Respiratory commensal/skin-flora Level 3 signal lacks non-host support and dominance.",
                )
            )
        if (
            "Low_read_mold_watchlist_preferred" in str(row.get("applied_rules") or "")
            or (
                str(row.get("genus") or "") == "aspergillus"
                and level == 3
                and mngs_ranked
                and host_only
                and non_dominant
                and reads_tier in LOW_SIGNAL_READ_TIERS
            )
        ):
            rules.append(
                (
                    "low_read_mold_watchlist_shadow_omit",
                    "Low-read mold without non-host support is better kept out of formal picked output.",
                )
            )
        if (
            is_low_actionability_virus(row)
            and mngs_ranked
            and host_only
            and non_dominant
            and not has_related_representative_context(row)
        ):
            rules.append(
                (
                    "low_actionability_virus_shadow_omit",
                    "Low-actionability viral DNA signal such as HPV/TTV is host-only and non-dominant; keep out of formal picked.",
                )
            )
        if is_non_core_respiratory_virus(row):
            rules.append(
                (
                    "non_core_respiratory_virus_context_shadow_omit",
                    "Exploratory respiratory virus audit: non-SARS, non-influenza, non-herpes virus should usually be RAG/context unless clinical evidence is strong.",
                )
            )
        if (
            is_water_or_low_specificity_environmental(row)
            and mngs_ranked
            and host_only
            and not has_related_representative_context(row)
        ):
            rules.append(
                (
                    "water_environmental_low_specificity_shadow_omit",
                    "Water/environmental low-pulmonary-specificity organism lacks non-host support; audit whether it should be RAG/context rather than picked.",
                )
            )
        if (
            mngs_ranked
            and host_only
            and not has_related_representative_context(row)
            and non_dominant
            and reads_tier in LOW_SIGNAL_READ_TIERS
            and is_nonprotected_for_low_signal_prune(row)
        ):
            rules.append(
                (
                    "mngs_r0_r1_host_only_nondominant_nonprotected_shadow_omit",
                    "mNGS-only R0/R1 non-dominant host-context signal is not in protected/common respiratory groups.",
                )
            )
        if (
            mngs_ranked
            and host_only
            and not has_related_representative_context(row)
            and non_dominant
            and reads_tier in LOW_SIGNAL_READ_TIERS
        ):
            rules.append(
                (
                    "mngs_r0_r1_host_only_nondominant_all_shadow_omit",
                    "Exploratory broad rule: any mNGS-only R0/R1 non-dominant host-context picked item.",
                )
            )
        if (
            mngs_ranked
            and host_only
            and not has_related_representative_context(row)
            and non_dominant
            and not is_typical_strong_respiratory_pathogen(row)
        ):
            rules.append(
                (
                    "mngs_only_no_same_support_nondominant_non_typical_shadow_omit",
                    "Exploratory picked FP rule: mNGS-only D0/D1 signal lacks same-organism hospital support and is outside typical strong respiratory pathogen groups.",
                )
            )
        if (
            mngs_ranked
            and host_only
            and not has_related_representative_context(row)
            and non_dominant
            and not is_typical_strong_respiratory_pathogen(row)
            and not is_deferred_gi_urinary_colonizer_for_shadow(row)
            and not is_answer_sensitive_healthcare_opportunist_for_shadow(row)
        ):
            rules.append(
                (
                    "mngs_only_no_same_support_nondominant_non_typical_guarded_shadow_omit",
                    "Guarded exploratory picked FP rule: mNGS-only D0/D1 signal outside typical strong respiratory pathogen groups, excluding deferred GI/urinary colonizer-prone organisms and Acinetobacter.",
                )
            )
        if (
            mngs_ranked
            and host_only
            and not has_related_representative_context(row)
            and non_dominant
            and reads_tier == "R2_medium"
            and is_nonprotected_for_low_signal_prune(row)
        ):
            rules.append(
                (
                    "mngs_r2_host_only_nondominant_nonprotected_shadow_omit",
                    "Exploratory rule: mNGS-only R2 non-dominant host-context signal outside protected/common respiratory groups.",
                )
            )
        if (
            is_mold_or_opportunistic_fungus(row)
            and mngs_ranked
            and host_only
            and not has_related_representative_context(row)
            and non_dominant
            and reads_tier in {"R0_trace", "R1_low", "R2_medium"}
            and not is_high_consequence_pathogen(row)
        ):
            rules.append(
                (
                    "weak_mold_fungal_host_only_nondominant_shadow_omit",
                    "Weak mold/opportunistic fungal signal is host-only, non-dominant, and lacks non-host support.",
                )
            )
        if is_candida_or_yeast(row) and not has_explicit_invasive_candida_context(row):
            rules.append(
                (
                    "candida_yeast_no_invasive_picked_shadow_omit",
                    "Exploratory Candida/yeast picked audit: no explicit blood, sterile-site, tissue, pleural, abscess, or invasive context.",
                )
            )
        if (
            is_enterococcus(row)
            and mngs_ranked
            and host_only
            and not has_related_representative_context(row)
            and non_dominant
        ):
            rules.append(
                (
                    "enterococcus_mngs_host_only_nondominant_shadow_omit",
                    "Exploratory Enterococcus audit: mNGS host-only non-dominant signal without same-organism non-host support.",
                )
            )
        if hospital_only and level == 2:
            rules.append(
                (
                    "hospital_only_l2_no_mngs_shadow_omit",
                    "Exploratory rule: Level 2 hospital-only picked item without mNGS ranked support.",
                )
            )

        if rules:
            row["shadow_action"] = "remove_from_formal_picked"
            row["primary_shadow_rule"] = rules[0][0]
            row["all_shadow_rules"] = "; ".join(rule for rule, _ in rules)
            row["shadow_reason"] = rules[0][1]


def output_entries(rows: list[dict[str, Any]]) -> list[NameEntry]:
    entries: list[NameEntry] = []
    for row in rows:
        entry = name_entry(str(row.get("row_id")), row.get("organism_name"))
        if entry:
            entries.append(entry)
    return entries


def apply_matches(
    *,
    rows: list[dict[str, Any]],
    answers: list[NameEntry],
    prefix: str,
) -> tuple[int, list[str]]:
    matched, missing = pairwise_match(output_entries(rows), answers)
    hit_count = 0
    for row in rows:
        info = matched.get(str(row.get("row_id")))
        if info:
            hit_count += 1
            row[f"{prefix}_answer_hit"] = True
            row[f"{prefix}_match_type"] = info["match_type"]
            row[f"{prefix}_matched_answer"] = info["answer"]
        else:
            row[f"{prefix}_answer_hit"] = False
            row[f"{prefix}_match_type"] = ""
            row[f"{prefix}_matched_answer"] = ""
    return hit_count, missing


def safe_ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def row_shadow_rules(row: dict[str, Any]) -> set[str]:
    return {
        value.strip()
        for value in str(row.get("all_shadow_rules") or "").split(";")
        if value.strip()
    }


def hit_count_for_rows(rows: list[dict[str, Any]], answer_entries: list[NameEntry]) -> tuple[int, list[str]]:
    entries = output_entries(rows)
    matched, missing = pairwise_match(entries, answer_entries)
    return len(matched), missing


def profile_metrics(
    rows_by_patient: dict[str, list[dict[str, Any]]],
    answers: dict[str, list[NameEntry]],
    profile_rules: set[str],
) -> dict[str, Any]:
    baseline_output = 0
    baseline_hits = 0
    profile_output = 0
    profile_hits = 0
    answer_total = 0
    removed_output = 0
    patient_rows: list[dict[str, Any]] = []
    removed_answer_patient_rows: list[dict[str, Any]] = []

    for pid, answer_entries in sorted(
        answers.items(), key=lambda item: (int(item[0]) if item[0].isdigit() else 10**9, item[0])
    ):
        rows = rows_by_patient.get(pid, [])
        removed = [row for row in rows if row_shadow_rules(row) & profile_rules]
        retained = [row for row in rows if row not in removed]
        base_hits, base_missing = hit_count_for_rows(rows, answer_entries)
        retained_hits, retained_missing = hit_count_for_rows(retained, answer_entries)
        baseline_output += len(rows)
        baseline_hits += base_hits
        profile_output += len(retained)
        profile_hits += retained_hits
        answer_total += len(answer_entries)
        removed_output += len(rows) - len(retained)
        for row in removed:
            rules = sorted(row_shadow_rules(row) & profile_rules)
            removed_answer_patient_rows.append(
                {
                    "patient_id": pid,
                    "organism_name": row.get("organism_name"),
                    "baseline_answer_hit": bool(row.get("baseline_answer_hit")),
                    "matched_answer": row.get("baseline_matched_answer"),
                    "rules": rules,
                    "clinical_ecology_group": row.get("clinical_ecology_group"),
                    "basis_level": row.get("basis_level"),
                    "evidence_source": row.get("evidence_source"),
                    "reads_tier": row.get("reads_tier"),
                    "dominance_tier": row.get("dominance_tier"),
                    "support_modules": row.get("support_modules"),
                }
            )
        patient_rows.append(
            {
                "patient_id": pid,
                "baseline_output": len(rows),
                "baseline_hits": base_hits,
                "profile_output": len(retained),
                "profile_hits": retained_hits,
                "removed_output": len(rows) - len(retained),
                "hit_delta": base_hits - retained_hits,
                "baseline_missing": base_missing,
                "profile_missing": retained_missing,
            }
        )

    return {
        "picked_output": profile_output,
        "picked_matched": profile_hits,
        "removed_output": removed_output,
        "removed_matched_delta": baseline_hits - profile_hits,
        "precision": safe_ratio(profile_hits, profile_output),
        "recall": safe_ratio(profile_hits, answer_total),
        "patient_rows_with_hit_loss": [
            row for row in patient_rows if row["hit_delta"] > 0
        ],
        "removed_answer_patient_rows": removed_answer_patient_rows,
        "removed_answer_hits": [
            row for row in removed_answer_patient_rows if row.get("baseline_answer_hit")
        ],
    }


def build_profile_metrics(
    rows_by_patient: dict[str, list[dict[str, Any]]],
    answers: dict[str, list[NameEntry]],
) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for profile_name, config in SHADOW_PROFILES.items():
        rules = {str(value) for value in config["rules"]}
        metrics = profile_metrics(rows_by_patient, answers, rules)
        metrics["description"] = config["description"]
        metrics["rules"] = sorted(rules)
        output[profile_name] = metrics
    return output


def summarize_metrics(rows_by_patient: dict[str, list[dict[str, Any]]], answers: dict[str, list[NameEntry]]) -> dict[str, Any]:
    baseline_output = 0
    baseline_hits = 0
    shadow_output = 0
    shadow_hits = 0
    answer_total = 0
    patient_rows: list[dict[str, Any]] = []

    for pid, answer_entries in sorted(answers.items(), key=lambda item: (int(item[0]) if item[0].isdigit() else 10**9, item[0])):
        rows = rows_by_patient.get(pid, [])
        retained = [row for row in rows if row.get("shadow_action") != "remove_from_formal_picked"]
        base_hit_count, base_missing = apply_matches(rows=rows, answers=answer_entries, prefix="baseline")
        shadow_hit_count, shadow_missing = apply_matches(rows=retained, answers=answer_entries, prefix="after_shadow")
        baseline_output += len(rows)
        baseline_hits += base_hit_count
        shadow_output += len(retained)
        shadow_hits += shadow_hit_count
        answer_total += len(answer_entries)
        patient_rows.append(
            {
                "patient_id": pid,
                "answers": "; ".join(entry.raw for entry in answer_entries),
                "answer_count": len(answer_entries),
                "baseline_output": len(rows),
                "baseline_hits": base_hit_count,
                "baseline_precision": safe_ratio(base_hit_count, len(rows)),
                "baseline_recall": safe_ratio(base_hit_count, len(answer_entries)),
                "shadow_output": len(retained),
                "shadow_hits": shadow_hit_count,
                "shadow_precision": safe_ratio(shadow_hit_count, len(retained)),
                "shadow_recall": safe_ratio(shadow_hit_count, len(answer_entries)),
                "removed_output": len(rows) - len(retained),
                "baseline_missing": base_missing,
                "shadow_missing": shadow_missing,
            }
        )
    return {
        "answer_patient_count": len(answers),
        "answer_organism_count": answer_total,
        "baseline": {
            "picked_output": baseline_output,
            "picked_matched": baseline_hits,
            "precision": safe_ratio(baseline_hits, baseline_output),
            "recall": safe_ratio(baseline_hits, answer_total),
        },
        "shadow": {
            "picked_output": shadow_output,
            "picked_matched": shadow_hits,
            "removed_output": baseline_output - shadow_output,
            "removed_matched_delta": baseline_hits - shadow_hits,
            "precision": safe_ratio(shadow_hits, shadow_output),
            "recall": safe_ratio(shadow_hits, answer_total),
        },
        "profiles": build_profile_metrics(rows_by_patient, answers),
        "patient_rows": patient_rows,
    }


def count_rows(rows: list[dict[str, Any]], key: str, *, only_fp: bool = False) -> dict[str, int]:
    counter: Counter[str] = Counter()
    for row in rows:
        if only_fp and row.get("baseline_answer_hit"):
            continue
        value = str(row.get(key) or "Unknown")
        counter[value] += 1
    return dict(counter.most_common())


def summarize(rows: list[dict[str, Any]], metrics: dict[str, Any], skipped: list[dict[str, str]]) -> dict[str, Any]:
    answer_rows = [row for row in rows if row.get("answer_patient")]
    answer_fp = [row for row in answer_rows if not row.get("baseline_answer_hit")]
    removed_rows = [row for row in rows if row.get("shadow_action") == "remove_from_formal_picked"]
    removed_answer_rows = [row for row in answer_rows if row.get("shadow_action") == "remove_from_formal_picked"]
    return {
        "processed_picked_rows": len(rows),
        "processed_answer_patient_picked_rows": len(answer_rows),
        "baseline_answer_patient_false_positive_rows": len(answer_fp),
        "shadow_removed_rows_all_patients": len(removed_rows),
        "shadow_removed_rows_answer_patients": len(removed_answer_rows),
        "skipped_outputs": skipped,
        "metrics": metrics,
        "fp_by_primary_shadow_rule": count_rows(answer_fp, "primary_shadow_rule"),
        "fp_by_genus": count_rows(answer_fp, "genus", only_fp=False),
        "fp_by_basis_level": count_rows(answer_fp, "basis_level", only_fp=False),
        "fp_by_evidence_source": count_rows(answer_fp, "evidence_source", only_fp=False),
        "fp_by_reads_tier": count_rows(answer_fp, "reads_tier", only_fp=False),
        "fp_by_dominance_tier": count_rows(answer_fp, "dominance_tier", only_fp=False),
        "fp_by_clinical_ecology_group": count_rows(answer_fp, "clinical_ecology_group", only_fp=False),
        "fp_by_guardrail_rule": count_rows(answer_fp, "guardrail_rule", only_fp=False),
        "removed_by_primary_shadow_rule_all": count_rows(removed_rows, "primary_shadow_rule"),
        "removed_by_primary_shadow_rule_answer_patients": count_rows(removed_answer_rows, "primary_shadow_rule"),
        "removed_by_clinical_ecology_group_answer_patients": count_rows(
            removed_answer_rows, "clinical_ecology_group"
        ),
        "removed_answer_hits": [
            {
                "patient_id": row.get("patient_id"),
                "organism_name": row.get("organism_name"),
                "matched_answer": row.get("baseline_matched_answer"),
                "primary_shadow_rule": row.get("primary_shadow_rule"),
                "clinical_ecology_group": row.get("clinical_ecology_group"),
            }
            for row in removed_answer_rows
            if row.get("baseline_answer_hit")
        ],
        "top_answer_patient_false_positive_rows": [
            {
                "patient_id": row.get("patient_id"),
                "organism_name": row.get("organism_name"),
                "basis_level": row.get("basis_level"),
                "reads": row.get("reads"),
                "reads_tier": row.get("reads_tier"),
                "dominance_tier": row.get("dominance_tier"),
                "evidence_source": row.get("evidence_source"),
                "guardrail_rule": row.get("guardrail_rule"),
                "primary_shadow_rule": row.get("primary_shadow_rule"),
            }
            for row in answer_fp[:50]
        ],
    }


def strip_internal_fields(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for row in rows:
        clean = {key: value for key, value in row.items() if not key.startswith("_") and key != "row_id"}
        output.append(clean)
    return output


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    patient_root = resolve_path(args.patient_root)
    answer_csv = resolve_path(args.answer_csv)
    output_csv = resolve_path(args.output_csv)
    output_json = resolve_path(args.output_json)
    selected_ids = {str(value).strip() for value in args.patient_id if str(value).strip()}

    answers, raw_answers, answer_column, answer_row_count = load_answers(answer_csv, args.answer_column)
    rows_by_patient: dict[str, list[dict[str, Any]]] = {}
    skipped: list[dict[str, str]] = []

    for patient_dir in patient_dirs(patient_root, selected_ids):
        pid = patient_id_from_dir(patient_dir)
        path = output_path(patient_dir, args.max_suffix)
        if not path.is_file():
            skipped.append({"patient_id": pid, "reason": "missing_output", "path": str(path)})
            continue
        payload = load_json(path)
        rows = build_picked_rows(label=args.label, patient_id=pid, payload=payload)
        for row in rows:
            row["answer_patient"] = pid in answers
        rows_by_patient[pid] = rows

    all_rows = [row for rows in rows_by_patient.values() for row in rows]
    add_shadow_rules(all_rows)

    metrics = summarize_metrics(rows_by_patient, answers)

    # Fill after-shadow fields for no-answer patients so the CSV remains explicit.
    for row in all_rows:
        if not row.get("answer_patient"):
            row["baseline_answer_hit"] = False
            row["after_shadow_answer_hit"] = False

    clean_rows = strip_internal_fields(all_rows)
    write_csv(output_csv, clean_rows)
    report = {
        "label": args.label,
        "patient_root": str(patient_root),
        "max_suffix": args.max_suffix,
        "answer_csv": str(answer_csv),
        "answer_column": answer_column,
        "answer_csv_row_count": answer_row_count,
        "answer_patient_count": len(answers),
        "raw_answers": raw_answers,
        "summary": summarize(all_rows, metrics, skipped),
        "rows": clean_rows,
    }
    write_json(output_json, report)
    if not args.quiet:
        print(json.dumps(report["summary"]["metrics"], ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

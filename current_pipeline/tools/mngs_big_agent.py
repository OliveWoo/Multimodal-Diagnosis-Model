from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import pathogen_normalization as pathogen_names
from utils import confirm_overwrite, sanitize_filename

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
FULL_PROMPT_PATH = Path("agents") / "prompts" / "mNGS_max_prompt.txt"
NO_FILMARRAY_PROMPT_PATH = Path("agents") / "prompts" / "mNGS_max_prompt_no_filmarray.txt"
PROMPT_PATH = FULL_PROMPT_PATH
DEFAULT_RANKED_MNGS_PATH = Path("mngs_candidate_microbes_ranked.json")
LOCAL_CHOSEN_RANKED_MNGS_NAME_TEMPLATE = "{base_name}_mNGS_ranked_candidates_chosen.json"
LOCAL_RANKED_MNGS_NAME_TEMPLATE = "{base_name}_mNGS_ranked_candidates.json"
MNGS_SOURCE_CHOICES = ("grouped", "all_rk_ntc", "ranked_only")
DEFAULT_MNGS_SOURCE = "grouped"
SUMMARY_DIR_NAME = "summary_outputs"
SUMMARY_OUTPUT_DIR_NAME = "summary_outputs"
DEFAULT_MODEL = "gpt-5"
OUTPUT_SUFFIX = ""
FAILED_DIR = Path("outputs") / "logs" / "failed"
NAME_SANITIZER = pathogen_names.NAME_SANITIZER
MATCH_ALIASES = pathogen_names.alias_mapping()
DISPLAY_NAME_BY_NORMALIZED = pathogen_names.display_mapping()
PROTECTED_EXACT_NAMES = {
    "pneumocystisjirovecii",
    "humanbetaherpesvirus5",
    "humanalphaherpesvirus1",
    "humanalphaherpesvirus2",
    "toxoplasmagondii",
}
PROTECTED_PREFIXES = (
    "nocardia",
    "legionella",
    "mycobacteriumtuberculosis",
    "mycobacteriumavium",
    "mycobacteriumkansasii",
    "mycobacteriumabscessus",
    "mycobacteroidesabscessus",
    "mycolicibacteriumabscessus",
    "mycolicibacterabscessus",
    "aspergillus",
    "mucorales",
    "rhizopus",
    "mucor",
    "rhizomucor",
    "lichtheimia",
    "cunninghamella",
    "saksenaea",
    "apophysomyces",
    "histoplasma",
    "cryptococcus",
)
ORAL_UPPER_AIRWAY_FLORA_EXACT_NAMES = {
    "streptococcussalivarius",
    "streptococcusparasanguinis",
    "viridansstreptococci",
    "rothiamucilaginosa",
}
ORAL_UPPER_AIRWAY_FLORA_PREFIXES = (
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
)
BACKGROUND_EXACT_NAMES = {
    "humangammaherpesvirus4",
    "humanbetaherpesvirus6",
    "humanbetaherpesvirus6a",
    "humanbetaherpesvirus6b",
    "humanbetaherpesvirus7",
    "staphylococcusepidermidis",
    "staphylococcushaemolyticus",
}
BACKGROUND_PREFIXES = ("corynebacterium",)
YEAST_LIKE_BACKGROUND_PREFIXES = ("candida", "trichosporon", "yeast")
SUMMARY_MODE_TO_SUFFIX = {
    "full": "with_filmarray",
    "deterministic": "with_filmarray_deterministic",
    "no_filmarray": "no_filmarray",
}
SUMMARY_MODE_COMPAT_SUFFIXES = {
    "full": ["with_filmarray", "with_underlying_with_filmarray"],
    "deterministic": ["with_filmarray_deterministic"],
    "no_filmarray": ["no_filmarray", "with_underlying_no_filmarray"],
}
SUMMARY_MODE_TO_PROMPT_PATH = {
    "full": FULL_PROMPT_PATH,
    "deterministic": FULL_PROMPT_PATH,
    "no_filmarray": NO_FILMARRAY_PROMPT_PATH,
}
SUMMARY_MODE_TO_OUTPUT_STEM = {
    "full": "mNGS_max_agent_full",
    "deterministic": "mNGS_max_agent_full_deterministic_summary",
    "no_filmarray": "mNGS_max_agent_no_filmarray",
}

try:
    import tiktoken
except ImportError:  # pragma: no cover - optional dependency
    tiktoken = None  # type: ignore[assignment]


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format=LOG_FORMAT)


def extract_patient_identifier(input_dir: Path) -> str:
    base_name = sanitize_filename(input_dir.name)
    if base_name.endswith("_json"):
        base_name = base_name[: -len("_json")] or base_name
    return base_name


def load_prompt(path: Path | None = None) -> str:
    prompt_path = path or PROMPT_PATH
    if not prompt_path.exists():
        raise FileNotFoundError(f"Prompt template not found: {prompt_path}")
    return prompt_path.read_text(encoding="utf-8")


def load_json(filepath: Path) -> dict:
    with filepath.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def find_final_summary_file(input_dir: Path, *, summary_mode: str = "auto") -> Path | None:
    if summary_mode not in {"auto", *SUMMARY_MODE_TO_SUFFIX.keys()}:
        raise ValueError(f"Unknown summary mode: {summary_mode}")

    base_name = extract_patient_identifier(input_dir)
    summary_dir = input_dir / SUMMARY_DIR_NAME
    search_dirs: list[Path] = [summary_dir, input_dir]

    if summary_mode == "auto":
        ordered_modes = ["full", "no_filmarray"]
    else:
        ordered_modes = [summary_mode]

    for mode in ordered_modes:
        for suffix in SUMMARY_MODE_COMPAT_SUFFIXES[mode]:
            filename = f"{base_name}_final_summary_{suffix}.json"
            for directory in search_dirs:
                path = directory / filename
                if path.exists():
                    return path

    legacy_filename = f"{base_name}_final_summary.json"
    for directory in search_dirs:
        path = directory / legacy_filename
        if path.exists():
            return path

    return None


def find_mngs_name_reads_file(input_dir: Path) -> Path | None:
    base_name = extract_patient_identifier(input_dir)
    candidates = [
        input_dir / f"{base_name}_mNGS_grouped.json",
        input_dir / f"{base_name}_mngs_name_reads.json",
        input_dir / f"{base_name}_mNGS_candidates_for_llm.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def find_all_rk_ntc_file(input_dir: Path) -> Path | None:
    base_name = extract_patient_identifier(input_dir)
    exact = input_dir / f"{base_name}_all_RK_NTC_microbes.json"
    if exact.exists():
        return exact
    matches = sorted(input_dir.glob(f"{base_name}_all_RK_NTC_microbes*.json"))
    return matches[0] if matches else None


def find_mngs_source_file(input_dir: Path, *, mngs_source: str = DEFAULT_MNGS_SOURCE) -> Path | None:
    if mngs_source == "all_rk_ntc":
        return find_all_rk_ntc_file(input_dir)
    if mngs_source == "ranked_only":
        return None
    if mngs_source == "grouped":
        return find_mngs_name_reads_file(input_dir)
    raise ValueError(f"Unsupported mNGS source: {mngs_source}")


def normalize_patient_id(input_dir: Path) -> str:
    patient_id = extract_patient_identifier(input_dir)
    if patient_id.startswith("NGS_patient_"):
        patient_id = patient_id.removeprefix("NGS_patient_")
    return patient_id


def find_ranked_mngs_file(
    input_dir: Path,
    ranked_path: Path | None = DEFAULT_RANKED_MNGS_PATH,
) -> Path | None:
    base_name = extract_patient_identifier(input_dir)
    local_candidates = [
        input_dir / LOCAL_CHOSEN_RANKED_MNGS_NAME_TEMPLATE.format(base_name=base_name),
        input_dir / LOCAL_RANKED_MNGS_NAME_TEMPLATE.format(base_name=base_name),
        input_dir / f"{base_name}_mngs_candidate_microbes_ranked.json",
    ]
    explicit_ranked_path = ranked_path is not None and ranked_path != DEFAULT_RANKED_MNGS_PATH
    candidates = []
    if explicit_ranked_path:
        candidates.append(ranked_path)
    candidates.extend(local_candidates)
    if ranked_path is not None and not explicit_ranked_path:
        candidates.append(ranked_path)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def load_ranked_mngs_for_patient(ranked_path: Path | None, input_dir: Path) -> dict | None:
    source_path = find_ranked_mngs_file(input_dir, ranked_path)
    if source_path is None:
        return None
    try:
        payload = load_json(source_path)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.warning("Unable to load ranked mNGS file %s: %s", source_path, exc)
        return None
    if isinstance(payload, dict) and isinstance(payload.get("records"), list):
        return deduplicate_ranked_mngs_payload(
            {
                "patient_id": str(payload.get("patient_id") or normalize_patient_id(input_dir)),
                "records": payload.get("records", []),
                "ranking_metadata": payload.get("ranking_metadata", {}),
            }
        )
    patients = payload.get("patients") if isinstance(payload, dict) else None
    if not isinstance(patients, dict):
        return None
    patient_id = normalize_patient_id(input_dir)
    ranked_patient = patients.get(patient_id) or patients.get(f"NGS_patient_{patient_id}")
    if not isinstance(ranked_patient, dict):
        return None
    return deduplicate_ranked_mngs_payload(
        {
            "patient_id": patient_id,
            "records": ranked_patient.get("records", []),
            "ranking_metadata": payload.get("ranking_metadata", {}),
        }
    )


def deduplicate_ranked_mngs_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Keep one deterministic best ranked mNGS record per organism before LLM use."""
    records = payload.get("records")
    if not isinstance(records, list):
        return payload

    flattened: list[tuple[int, dict[str, Any], dict[str, Any]]] = []
    exact_keys: list[str] = []
    for record_index, record in enumerate(records):
        if not isinstance(record, dict):
            continue
        specimen_code = str(record.get("specimen_code") or "")
        seq_id = str(record.get("seq_id") or "")
        for candidate in record.get("candidates") or []:
            if not isinstance(candidate, dict):
                continue
            normalized_name = normalize_organism_name(
                candidate.get("organism_name") or candidate.get("name")
            )
            if not normalized_name:
                continue
            item = dict(candidate)
            item.setdefault("source_specimen_code", specimen_code)
            item.setdefault("source_seq_id", seq_id)
            flattened.append((record_index, record, item))
            exact_keys.append(
                canonical_json(
                    {
                        "specimen_code": specimen_code,
                        "seq_id": seq_id,
                        "candidate": item,
                    }
                )
            )

    if not flattened:
        return payload

    selected: dict[str, tuple[int, dict[str, Any], dict[str, Any]]] = {}
    for entry in flattened:
        _, _, candidate = entry
        key = normalize_organism_name(candidate.get("organism_name") or candidate.get("name"))
        current = selected.get(key)
        if current is None or ranked_candidate_sort_key(entry) < ranked_candidate_sort_key(current):
            selected[key] = entry

    candidates_by_record: dict[int, list[dict[str, Any]]] = {}
    source_records: dict[int, dict[str, Any]] = {}
    for record_index, record, candidate in selected.values():
        candidates_by_record.setdefault(record_index, []).append(candidate)
        source_records[record_index] = record

    normalized_records: list[dict[str, Any]] = []
    for record_index in sorted(candidates_by_record):
        source_record = source_records[record_index]
        normalized_record = dict(source_record)
        candidates = sorted(
            candidates_by_record[record_index],
            key=lambda item: ranked_candidate_sort_key((record_index, source_record, item)),
        )
        normalized_record["candidates"] = candidates
        normalized_record["status"] = "parsed" if candidates else source_record.get("status", "parsed_no_match")
        normalized_records.append(normalized_record)

    exact_duplicate_count = len(exact_keys) - len(set(exact_keys))
    ranking_metadata = dict(payload.get("ranking_metadata") or {})
    ranking_metadata["deduplication"] = {
        "enabled": True,
        "method": "exact_duplicate_removed_then_best_record_per_normalized_organism",
        "input_candidate_rows": len(flattened),
        "output_candidate_rows": sum(len(record.get("candidates") or []) for record in normalized_records),
        "exact_duplicate_rows_removed": exact_duplicate_count,
        "same_organism_rows_removed": len(flattened) - len(selected),
        "best_record_sort": [
            "rank_priority ascending",
            "reads descending",
            "reads_percentile descending",
            "final_score descending",
            "specimen_code ascending",
            "organism_name ascending",
        ],
    }

    normalized_payload = dict(payload)
    normalized_payload["records"] = normalized_records
    normalized_payload["ranking_metadata"] = ranking_metadata
    return normalized_payload


def ranked_candidate_sort_key(entry: tuple[int, dict[str, Any], dict[str, Any]]) -> tuple[Any, ...]:
    record_index, record, candidate = entry
    ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
    priority = to_rank_int(ranking.get("rank_priority"))
    reads = to_float(candidate.get("reads", candidate.get("sec_hit", 0)))
    reads_percentile = to_float(ranking.get("reads_percentile"))
    final_score = to_float(ranking.get("final_score"))
    return (
        priority,
        -reads,
        -reads_percentile,
        -final_score,
        str(record.get("specimen_code") or ""),
        normalize_organism_name(candidate.get("organism_name") or candidate.get("name")),
        record_index,
    )


def to_rank_int(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 999


def to_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def normalize_organism_name(value: Any) -> str:
    return pathogen_names.canonical_key(value)


def canonical_display_name(value: Any) -> str:
    return pathogen_names.display_name(value)


def is_protected_pathogen_name(value: Any) -> bool:
    normalized = normalize_organism_name(value)
    return normalized in PROTECTED_EXACT_NAMES or any(
        normalized.startswith(prefix) for prefix in PROTECTED_PREFIXES
    )


def is_oral_upper_airway_flora_name(value: Any) -> bool:
    normalized = normalize_organism_name(value)
    return normalized in ORAL_UPPER_AIRWAY_FLORA_EXACT_NAMES or any(
        normalized.startswith(prefix) for prefix in ORAL_UPPER_AIRWAY_FLORA_PREFIXES
    )


def is_background_pathogen_name(value: Any) -> bool:
    normalized = normalize_organism_name(value)
    return normalized in BACKGROUND_EXACT_NAMES or any(
        normalized.startswith(prefix) for prefix in BACKGROUND_PREFIXES
    )


def is_yeast_like_background_name(value: Any) -> bool:
    normalized = normalize_organism_name(value)
    return any(normalized.startswith(prefix) for prefix in YEAST_LIKE_BACKGROUND_PREFIXES)


def reads_tier(value: Any) -> str:
    reads = to_float(value)
    if reads >= 10000:
        return "R4_very_high"
    if reads >= 1000:
        return "R3_high"
    if reads >= 100:
        return "R2_medium"
    if reads >= 10:
        return "R1_low"
    if reads >= 1:
        return "R0_trace"
    return "Unknown"


def _ranked_evidence_map(ranked_mngs: dict[str, Any] | None) -> dict[str, dict[str, Any]]:
    if not isinstance(ranked_mngs, dict):
        return {}
    selected: dict[str, tuple[int, dict[str, Any], dict[str, Any]]] = {}
    for record_index, record in enumerate(ranked_mngs.get("records") or []):
        if not isinstance(record, dict):
            continue
        for candidate in record.get("candidates") or []:
            if not isinstance(candidate, dict):
                continue
            normalized = normalize_organism_name(
                candidate.get("organism_name") or candidate.get("name")
            )
            if not normalized:
                continue
            entry = (record_index, record, candidate)
            current = selected.get(normalized)
            if current is None or ranked_candidate_sort_key(entry) < ranked_candidate_sort_key(current):
                selected[normalized] = entry
    return {key: entry[2] for key, entry in selected.items()}


def _ranked_dominance_map(ranked_mngs: dict[str, Any] | None) -> dict[str, str]:
    if not isinstance(ranked_mngs, dict):
        return {}
    selected: dict[str, tuple[tuple[Any, ...], str]] = {}
    for record_index, record in enumerate(ranked_mngs.get("records") or []):
        if not isinstance(record, dict):
            continue
        candidates = [item for item in record.get("candidates") or [] if isinstance(item, dict)]
        if not candidates:
            continue
        sorted_candidates = sorted(
            candidates,
            key=lambda item: (
                -to_float(item.get("reads", item.get("sec_hit", 0))),
                normalize_organism_name(item.get("organism_name") or item.get("name")),
            ),
        )
        top_reads = to_float(sorted_candidates[0].get("reads", sorted_candidates[0].get("sec_hit", 0)))
        second_reads = (
            to_float(sorted_candidates[1].get("reads", sorted_candidates[1].get("sec_hit", 0)))
            if len(sorted_candidates) > 1
            else 0.0
        )
        ratio = top_reads / max(second_reads, 1.0)
        top_tier = "D3_dominant" if ratio >= 10 else "D2_moderate" if ratio >= 3 else "D1_low"
        for index, candidate in enumerate(sorted_candidates):
            normalized = normalize_organism_name(
                candidate.get("organism_name") or candidate.get("name")
            )
            if not normalized:
                continue
            tier = top_tier if index == 0 else "D0_not_top"
            sort_key = ranked_candidate_sort_key((record_index, record, candidate))
            current = selected.get(normalized)
            if current is None or sort_key < current[0]:
                selected[normalized] = (sort_key, tier)
    return {key: value[1] for key, value in selected.items()}


def _final_summary_evidence_map(final_summary: dict[str, Any] | None) -> dict[str, dict[str, Any]]:
    if not isinstance(final_summary, dict):
        return {}
    mapped: dict[str, dict[str, Any]] = {}
    for candidate in final_summary.get("pathogen_candidates") or []:
        if not isinstance(candidate, dict):
            continue
        normalized = normalize_organism_name(candidate.get("organism_name") or candidate.get("name"))
        if normalized:
            mapped[normalized] = candidate
    for evidence_item in final_summary.get("hospital_organism_evidence") or []:
        if not isinstance(evidence_item, dict):
            continue
        normalized = normalize_organism_name(evidence_item.get("organism_name") or evidence_item.get("name"))
        if not normalized:
            continue
        current = mapped.setdefault(
            normalized,
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
            current_modules = current.setdefault("evidence_modules", [])
            if isinstance(current_modules, list):
                for module in modules:
                    if module not in current_modules:
                        current_modules.append(module)
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


def _final_summary_level_rank(
    final_summary_by_name: dict[str, dict[str, Any]],
    organism_name: Any,
) -> int:
    candidate = final_summary_by_name.get(normalize_organism_name(organism_name))
    if not isinstance(candidate, dict):
        return 99
    return min(
        _level_rank(candidate.get("integrated_causative_level")),
        _level_rank(candidate.get("best_hospital_level")),
    )


def _rank_priority_from_evidence(evidence: dict[str, Any] | None) -> int:
    if not isinstance(evidence, dict):
        return 999
    ranking = evidence.get("ranking") if isinstance(evidence.get("ranking"), dict) else {}
    return to_rank_int(ranking.get("rank_priority"))


def _possibility_from_evidence(evidence: dict[str, Any] | None) -> str:
    if not isinstance(evidence, dict):
        return ""
    ranking = evidence.get("ranking") if isinstance(evidence.get("ranking"), dict) else {}
    return str(ranking.get("possibility_level") or "").strip().lower()


def _reads_percentile_from_evidence(evidence: dict[str, Any] | None) -> float:
    if not isinstance(evidence, dict):
        return 0.0
    ranking = evidence.get("ranking") if isinstance(evidence.get("ranking"), dict) else {}
    return to_float(ranking.get("reads_percentile"))


def _level_rank(value: Any) -> int:
    text = str(value or "").strip().lower()
    for index in range(1, 6):
        if text in {f"level {index}", f"level{index}"}:
            return index
    return 99


def normalize_number_value(value: Any) -> int | float | str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if numeric.is_integer():
        return int(numeric)
    return numeric


def _candidate_guardrail_rule(candidate: dict[str, Any]) -> str:
    key_evidence = candidate.get("key_evidence")
    if isinstance(key_evidence, dict):
        return str(key_evidence.get("guardrail_rule") or "")
    return ""


def _candidate_has_guardrail(candidate: dict[str, Any]) -> bool:
    return _candidate_guardrail_rule(candidate) in {
        "R-S4-BACKGROUND",
        "R-S4-ORAL-UPPER_AIRWAY_HARD",
        "R-S4-HSV_LOW_SIGNAL_HARD",
    }


def _set_guardrail_rule(candidate: dict[str, Any], guardrail_rule: str, level_rule: str) -> None:
    key_evidence = candidate.setdefault("key_evidence", {})
    if isinstance(key_evidence, dict):
        key_evidence["guardrail_rule"] = guardrail_rule
        key_evidence["level_rule"] = level_rule
    applied_rules = candidate.get("applied_rules")
    if isinstance(applied_rules, list) and guardrail_rule not in applied_rules:
        applied_rules.append(guardrail_rule)


def _clear_candidate_guardrail(candidate: dict[str, Any]) -> None:
    key_evidence = candidate.get("key_evidence")
    if isinstance(key_evidence, dict):
        key_evidence["guardrail_rule"] = ""
    applied_rules = candidate.get("applied_rules")
    if isinstance(applied_rules, list):
        candidate["applied_rules"] = [
            rule
            for rule in applied_rules
            if rule
            not in {
                "R-S4-BACKGROUND",
                "R-S4-ORAL-UPPER_AIRWAY_HARD",
                "R-S4-HSV_LOW_SIGNAL_HARD",
            }
        ]


def _deterministic_guardrail_rule(
    candidate: dict[str, Any],
    final_summary_by_name: dict[str, dict[str, Any]],
) -> str:
    final_level = _final_summary_level_rank(final_summary_by_name, candidate.get("organism_name"))
    if final_level <= 2:
        return ""
    if is_oral_upper_airway_flora_name(candidate.get("organism_name")):
        return "R-S4-ORAL-UPPER_AIRWAY_HARD"
    if is_background_pathogen_name(candidate.get("organism_name")):
        return "R-S4-BACKGROUND"
    if is_yeast_like_background_name(candidate.get("organism_name")) and final_level > 2:
        return "R-S4-BACKGROUND"
    return ""


def _apply_guardrail(
    candidate: dict[str, Any],
    final_summary_by_name: dict[str, dict[str, Any]],
) -> None:
    guardrail = _deterministic_guardrail_rule(candidate, final_summary_by_name)
    if guardrail == "R-S4-ORAL-UPPER_AIRWAY_HARD":
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M4_weak"
        candidate["integrated_causative_level"] = "Level 4"
        module_support = candidate.get("module_support_summary")
        if isinstance(module_support, dict):
            module_support["mNGS"] = "M4_weak"
        _set_guardrail_rule(candidate, guardrail, "R-S5-L4")
        return
    if guardrail == "R-S4-BACKGROUND":
        candidate["is_likely_colonizer_or_background"] = True
        candidate["mngs_signal_tier"] = "M5_background"
        candidate["integrated_causative_level"] = "Level 5"
        module_support = candidate.get("module_support_summary")
        if isinstance(module_support, dict):
            module_support["mNGS"] = "M5_background"
        _set_guardrail_rule(candidate, guardrail, "R-S5-L5")
        return
    if _candidate_has_guardrail(candidate):
        _clear_candidate_guardrail(candidate)


def _matches_m1_strong(
    evidence: dict[str, Any] | None,
    candidate: dict[str, Any],
    dominance_tier: str,
) -> bool:
    rank_priority = _rank_priority_from_evidence(evidence)
    tier = str(candidate.get("reads_tier") or reads_tier(candidate.get("reads")))
    percentile = _reads_percentile_from_evidence(evidence)
    if rank_priority in {1, 2} and tier in {"R4_very_high", "R3_high"}:
        return True
    if rank_priority in {1, 2} and percentile >= 0.85:
        return True
    return dominance_tier == "D3_dominant" and tier in {
        "R4_very_high",
        "R3_high",
        "R2_medium",
    } and rank_priority <= 3


def _matches_m2_moderate(evidence: dict[str, Any] | None, candidate: dict[str, Any]) -> bool:
    rank_priority = _rank_priority_from_evidence(evidence)
    tier = str(candidate.get("reads_tier") or reads_tier(candidate.get("reads")))
    percentile = _reads_percentile_from_evidence(evidence)
    possibility = _possibility_from_evidence(evidence)
    normalized = normalize_organism_name(candidate.get("organism_name"))
    if rank_priority in {1, 2} and tier == "R2_medium":
        return True
    if rank_priority <= 3 and percentile >= 0.60:
        return True
    if rank_priority <= 3 and possibility in {"high", "medium"} and tier in {
        "R2_medium",
        "R1_low",
    }:
        return True
    if rank_priority <= 4 and tier in {"R4_very_high", "R3_high"}:
        return True
    if (
        normalized in {"pneumocystisjirovecii", "humanbetaherpesvirus5"}
        and rank_priority == 1
        and possibility == "high"
        and tier in {"R2_medium", "R1_low"}
    ):
        return True
    return False


def _candidate_has_auxiliary_support(candidate: dict[str, Any]) -> bool:
    module_support = candidate.get("module_support_summary")
    if not isinstance(module_support, dict):
        return False
    unavailable = {"", "not_available", "not available", "unknown", "negative", "no_support"}
    for key in ("culture", "filmarray_gmtest", "image", "host"):
        value = str(module_support.get(key) or "").strip().lower()
        if value and value not in unavailable:
            return True
    return False


def _protected_retention_condition(
    evidence: dict[str, Any] | None,
    candidate: dict[str, Any],
) -> bool:
    rank_priority = _rank_priority_from_evidence(evidence)
    if rank_priority in {1, 2}:
        return True
    if rank_priority == 3:
        return _candidate_has_auxiliary_support(candidate)
    return rank_priority == 4 and str(candidate.get("reads_tier") or "") == "R4_very_high"


def _sync_candidate_from_ranked_evidence(
    candidate: dict[str, Any],
    evidence_by_name: dict[str, dict[str, Any]],
    dominance_by_name: dict[str, str],
    final_summary_by_name: dict[str, dict[str, Any]],
) -> None:
    raw_name = candidate.get("organism_name") or candidate.get("name")
    normalized = normalize_organism_name(raw_name)
    if normalized:
        candidate["organism_name"] = canonical_display_name(raw_name)
    evidence = evidence_by_name.get(normalized)
    if isinstance(evidence, dict):
        ranking = evidence.get("ranking") if isinstance(evidence.get("ranking"), dict) else {}
        reads_value = evidence.get("reads", evidence.get("sec_hit", candidate.get("reads", 0)))
        candidate["reads"] = int(to_float(reads_value))
        candidate["reads_tier"] = reads_tier(reads_value)
        if ranking.get("rank_priority") is not None:
            candidate["rank_priority"] = str(to_rank_int(ranking.get("rank_priority")))
        if ranking.get("rank_rule") is not None:
            candidate["rank_rule"] = str(ranking.get("rank_rule") or "")
        if ranking.get("reads_percentile") is not None:
            candidate["reads_percentile"] = normalize_number_value(ranking.get("reads_percentile"))
        if evidence.get("source_category") is not None:
            candidate["source_category"] = evidence.get("source_category")
    dominance_tier = dominance_by_name.get(normalized)
    if dominance_tier:
        candidate["dominance_tier"] = dominance_tier

    protected = is_protected_pathogen_name(candidate.get("organism_name"))
    candidate["is_protected_pathogen"] = protected
    candidate["protected_retention_condition"] = (
        _protected_retention_condition(evidence, candidate) if protected else False
    )

    _apply_guardrail(candidate, final_summary_by_name)
    if _candidate_has_guardrail(candidate):
        return
    candidate["is_likely_colonizer_or_background"] = False
    if _matches_m1_strong(evidence, candidate, str(candidate.get("dominance_tier") or "")):
        candidate["mngs_signal_tier"] = "M1_strong"
        module_support = candidate.get("module_support_summary")
        if isinstance(module_support, dict):
            module_support["mNGS"] = "M1_strong"
        if candidate.get("specimen_alignment") in {"Aligned", "Sterile_or_Systemic"}:
            candidate["integrated_causative_level"] = "Level 2"
            key_evidence = candidate.get("key_evidence")
            if isinstance(key_evidence, dict):
                key_evidence["guardrail_rule"] = ""
                key_evidence["level_rule"] = "R-S5-L2"
        return
    if _matches_m2_moderate(evidence, candidate):
        candidate["mngs_signal_tier"] = "M2_moderate"
        module_support = candidate.get("module_support_summary")
        if isinstance(module_support, dict):
            module_support["mNGS"] = "M2_moderate"
        if candidate.get("specimen_alignment") in {"Aligned", "Sterile_or_Systemic"}:
            candidate["integrated_causative_level"] = "Level 2"
            key_evidence = candidate.get("key_evidence")
            if isinstance(key_evidence, dict):
                key_evidence["guardrail_rule"] = ""
                key_evidence["level_rule"] = "R-S5-L2"


def _candidate_output_sort_key(candidate: dict[str, Any]) -> tuple[Any, ...]:
    return (
        _level_rank(candidate.get("integrated_causative_level") or candidate.get("observed_level")),
        to_rank_int(candidate.get("rank_priority")),
        -to_float(candidate.get("reads")),
        normalize_organism_name(candidate.get("organism_name")),
    )


def _dedupe_candidates(candidates: list[Any]) -> list[dict[str, Any]]:
    selected: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        normalized = normalize_organism_name(candidate.get("organism_name"))
        if not normalized:
            continue
        current = selected.get(normalized)
        if current is None or _candidate_output_sort_key(candidate) < _candidate_output_sort_key(current):
            selected[normalized] = candidate
    return sorted(selected.values(), key=_candidate_output_sort_key)


def _candidate_is_pickable(candidate: dict[str, Any]) -> bool:
    if candidate.get("is_likely_colonizer_or_background") is True:
        return False
    level = _level_rank(candidate.get("integrated_causative_level"))
    if level <= 2:
        return True
    return level == 3 and is_protected_pathogen_name(candidate.get("organism_name"))


def _picked_item_from_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    return {
        "organism_name": candidate.get("organism_name", ""),
        "classification": candidate.get("classification", "Unknown"),
        "picked_role": "Secondary",
        "basis_level": candidate.get("integrated_causative_level", "Unknown"),
        "mngs_signal_tier": candidate.get("mngs_signal_tier", "Unknown"),
        "rank_priority": candidate.get("rank_priority", "Unknown"),
        "reads": candidate.get("reads", 0),
    }


def _picked_output_sort_key(item: dict[str, Any]) -> tuple[Any, ...]:
    return (
        _level_rank(item.get("basis_level")),
        to_rank_int(item.get("rank_priority")),
        -to_float(item.get("reads")),
        normalize_organism_name(item.get("organism_name")),
    )


def _assign_picked_roles(best: dict[str, Any]) -> None:
    picked = [item for item in best.get("picked_pathogens", []) if isinstance(item, dict)]
    primary_assigned = False
    for item in sorted(picked, key=_picked_output_sort_key):
        level = _level_rank(item.get("basis_level"))
        if level <= 2:
            if not primary_assigned:
                item["picked_role"] = "Primary"
                primary_assigned = True
            else:
                item["picked_role"] = "Secondary"
        elif level == 3:
            item["picked_role"] = (
                "Protected_Level3"
                if is_protected_pathogen_name(item.get("organism_name"))
                else "BestAvailable"
            )


def _sync_likelihood_from_picks(payload: dict[str, Any]) -> None:
    best = payload.get("best_available_summary")
    if not isinstance(best, dict):
        return
    picked = best.get("picked_pathogens", [])
    level_ranks = [_level_rank(item.get("basis_level")) for item in picked if isinstance(item, dict)]
    if 1 in level_ranks:
        payload["final_infection_likelihood"] = "Severe"
        best["selection_mode"] = "Confirmed"
    elif 2 in level_ranks:
        payload["final_infection_likelihood"] = "Likely"
        best["selection_mode"] = "Confirmed"
    elif 3 in level_ranks:
        payload["final_infection_likelihood"] = "Possible"
        best["selection_mode"] = "BestAvailable"
    elif not level_ranks:
        best["selection_mode"] = "No_high_priority_candidate"


def _sync_picked_pathogens(payload: dict[str, Any]) -> None:
    best = payload.get("best_available_summary")
    candidates = payload.get("pathogen_candidates")
    if not isinstance(best, dict) or not isinstance(candidates, list):
        return
    candidate_by_name = {
        normalize_organism_name(candidate.get("organism_name")): candidate
        for candidate in candidates
        if isinstance(candidate, dict)
    }
    picked = best.get("picked_pathogens")
    if not isinstance(picked, list):
        picked = []
    synced_picked: list[dict[str, Any]] = []
    seen_names: set[str] = set()
    for item in picked:
        if not isinstance(item, dict):
            continue
        normalized = normalize_organism_name(item.get("organism_name"))
        candidate = candidate_by_name.get(normalized)
        if candidate is not None:
            if not _candidate_is_pickable(candidate):
                continue
            item["organism_name"] = candidate.get("organism_name", item.get("organism_name", ""))
            item["classification"] = candidate.get("classification", item.get("classification", "Unknown"))
            item["basis_level"] = candidate.get("integrated_causative_level", item.get("basis_level", "Unknown"))
            item["mngs_signal_tier"] = candidate.get("mngs_signal_tier", item.get("mngs_signal_tier", "Unknown"))
            item["rank_priority"] = candidate.get("rank_priority", item.get("rank_priority", "Unknown"))
            item["reads"] = candidate.get("reads", item.get("reads", 0))
        else:
            item["organism_name"] = canonical_display_name(item.get("organism_name"))
        normalized = normalize_organism_name(item.get("organism_name"))
        if normalized and normalized not in seen_names:
            synced_picked.append(item)
            seen_names.add(normalized)

    for candidate in sorted(
        (candidate for candidate in candidates if isinstance(candidate, dict)),
        key=_candidate_output_sort_key,
    ):
        normalized = normalize_organism_name(candidate.get("organism_name"))
        if normalized in seen_names or not _candidate_is_pickable(candidate):
            continue
        synced_picked.append(_picked_item_from_candidate(candidate))
        seen_names.add(normalized)

    synced_picked.sort(key=_picked_output_sort_key)
    best["picked_pathogens"] = synced_picked
    best["picked_count"] = len(synced_picked)
    _assign_picked_roles(best)
    _sync_likelihood_from_picks(payload)


def _sync_excluded_candidates(payload: dict[str, Any]) -> None:
    excluded = payload.get("excluded_candidates")
    candidates = payload.get("pathogen_candidates")
    best = payload.get("best_available_summary")
    if not isinstance(excluded, list) or not isinstance(candidates, list):
        return
    candidate_by_name = {
        normalize_organism_name(candidate.get("organism_name")): candidate
        for candidate in candidates
        if isinstance(candidate, dict)
    }
    picked_names: set[str] = set()
    if isinstance(best, dict) and isinstance(best.get("picked_pathogens"), list):
        picked_names = {
            normalize_organism_name(item.get("organism_name"))
            for item in best["picked_pathogens"]
            if isinstance(item, dict)
        }

    synced_excluded: list[dict[str, Any]] = []
    seen_names: set[str] = set()
    for item in excluded:
        if not isinstance(item, dict):
            continue
        item["organism_name"] = canonical_display_name(item.get("organism_name"))
        normalized = normalize_organism_name(item.get("organism_name"))
        candidate = candidate_by_name.get(normalized)
        if normalized in picked_names or (candidate is not None and _candidate_is_pickable(candidate)):
            continue
        if candidate is not None:
            item["organism_name"] = candidate.get("organism_name", item.get("organism_name", ""))
            item["classification"] = candidate.get("classification", item.get("classification", "Unknown"))
            item["observed_level"] = candidate.get(
                "integrated_causative_level",
                item.get("observed_level", "Unknown"),
            )
            item["rank_priority"] = candidate.get("rank_priority", item.get("rank_priority", "Unknown"))
            item["reads"] = candidate.get("reads", item.get("reads", 0))
        if normalized and normalized not in seen_names:
            synced_excluded.append(item)
            seen_names.add(normalized)

    for candidate in sorted(
        (candidate for candidate in candidates if isinstance(candidate, dict)),
        key=_candidate_output_sort_key,
    ):
        normalized = normalize_organism_name(candidate.get("organism_name"))
        if (
            not normalized
            or normalized in picked_names
            or normalized in seen_names
            or _candidate_is_pickable(candidate)
        ):
            continue
        synced_excluded.append(
            {
                "organism_name": candidate.get("organism_name", ""),
                "classification": candidate.get("classification", "Unknown"),
                "observed_level": candidate.get("integrated_causative_level", "Unknown"),
                "rank_priority": candidate.get("rank_priority", "Unknown"),
                "reads": candidate.get("reads", 0),
                "exclusion_reason_code": "colonizer_or_background"
                if candidate.get("is_likely_colonizer_or_background") is True
                else "not_picked_low_level",
            }
        )
        seen_names.add(normalized)

    synced_excluded.sort(key=_candidate_output_sort_key)
    payload["excluded_candidates"] = synced_excluded


def postprocess_mngs_max_output(
    payload: dict[str, Any],
    *,
    ranked_mngs: dict[str, Any] | None = None,
    mngs_grouped: dict[str, Any] | None = None,  # reserved for future deterministic fixes
    final_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    del mngs_grouped
    if not isinstance(payload, dict):
        return payload
    if isinstance(final_summary, dict) and final_summary.get("dominant_source"):
        payload["dominant_source"] = final_summary.get("dominant_source")
    evidence_by_name = _ranked_evidence_map(ranked_mngs)
    dominance_by_name = _ranked_dominance_map(ranked_mngs)
    final_summary_by_name = _final_summary_evidence_map(final_summary)
    for candidate in payload.get("pathogen_candidates") or []:
        if isinstance(candidate, dict):
            _sync_candidate_from_ranked_evidence(
                candidate,
                evidence_by_name,
                dominance_by_name,
                final_summary_by_name,
            )
    candidates = payload.get("pathogen_candidates")
    if isinstance(candidates, list):
        payload["pathogen_candidates"] = _dedupe_candidates(candidates)
    for excluded in payload.get("excluded_candidates") or []:
        if isinstance(excluded, dict):
            excluded["organism_name"] = canonical_display_name(excluded.get("organism_name"))
    _sync_picked_pathogens(payload)
    _sync_excluded_candidates(payload)
    return payload


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def detect_summary_variant(base_name: str, summary_path: Path) -> str:
    name = summary_path.name
    for mode, suffixes in SUMMARY_MODE_COMPAT_SUFFIXES.items():
        for suffix in suffixes:
            expected = f"{base_name}_final_summary_{suffix}.json"
            if name == expected:
                return mode
    if name == f"{base_name}_final_summary.json":
        return "legacy"
    return "custom"


def looks_like_patient_dir(
    path: Path,
    *,
    summary_mode: str,
    mngs_source: str = DEFAULT_MNGS_SOURCE,
) -> bool:
    if mngs_source == "ranked_only":
        return (
            find_final_summary_file(path, summary_mode=summary_mode) is not None
            and find_ranked_mngs_file(path) is not None
        )
    return (
        find_final_summary_file(path, summary_mode=summary_mode) is not None
        and find_mngs_source_file(path, mngs_source=mngs_source) is not None
    )


def collect_patient_directories(
    entries: Sequence[Path],
    *,
    summary_mode: str,
    mngs_source: str = DEFAULT_MNGS_SOURCE,
) -> list[Path]:
    collected: list[Path] = []
    seen: set[Path] = set()

    for entry in entries:
        path = entry.expanduser()
        if not path.exists():
            logging.warning("Skipping non-existent path: %s", path)
            continue

        resolved = path.resolve()
        if resolved.is_file():
            logging.warning("Skipping file input (expecting directories): %s", resolved)
            continue
        if not resolved.is_dir():
            logging.warning("Skipping non-directory path: %s", resolved)
            continue

        candidates = [resolved]
        try:
            candidates.extend(child.resolve() for child in resolved.iterdir() if child.is_dir())
        except PermissionError:
            logging.warning("Unable to enumerate subdirectories for %s", resolved)

        for candidate in candidates:
            if candidate in seen:
                continue
            if looks_like_patient_dir(
                candidate,
                summary_mode=summary_mode,
                mngs_source=mngs_source,
            ):
                collected.append(candidate)
                seen.add(candidate)
            else:
                logging.debug("Ignoring directory without required inputs: %s", candidate)

    collected.sort()
    return collected


def build_prompt(
    template: str,
    *,
    final_summary: dict,
    mngs_grouped: dict | None,
    final_path: Path,
    mngs_path: Path | None,
    ranked_mngs: dict | None = None,
    ranked_path: Path | None = None,
) -> str:
    blocks: list[str] = [
        template.rstrip(),
        "",
        "### JSON payloads ###",
        f"#### final_summary ({final_path.name})",
        json.dumps(final_summary, ensure_ascii=False, indent=2),
        "",
    ]
    if mngs_grouped is not None and mngs_path is not None:
        blocks.extend(
            [
                f"#### mngs_name_reads ({mngs_path.name})",
                json.dumps(mngs_grouped, ensure_ascii=False, indent=2),
            ]
        )
    if ranked_mngs is not None:
        label = ranked_path.name if ranked_path is not None else "ranked_mngs"
        blocks.extend(
            [
                "",
                f"#### mngs_ranked_candidates ({label})",
                json.dumps(ranked_mngs, ensure_ascii=False, indent=2),
            ]
        )
    return "\n".join(blocks)


def _strip_json_markers(text: str) -> str:
    stripped = text.strip()
    if not stripped.startswith("```"):
        return stripped
    lines = stripped.splitlines()
    if len(lines) >= 3 and lines[-1].strip() == "```":
        return "\n".join(lines[1:-1]).strip()
    return stripped


def validate_json(output_text: str) -> dict | None:
    candidate = _strip_json_markers(output_text)
    try:
        return json.loads(candidate)
    except json.JSONDecodeError:
        logging.error("Failed to parse LLM output as JSON.")
        logging.debug("Raw output that failed to parse:\n%s", output_text)
        return None


def estimate_token_count(text: str, *, model: str) -> int | None:
    if tiktoken is None:
        return None
    try:
        encoding = tiktoken.encoding_for_model(model)
    except KeyError:
        encoding = tiktoken.get_encoding("cl100k_base")
    return len(encoding.encode(text))


def _extract_output_text(response: Any) -> str:
    if getattr(response, "output_text", None):
        text = response.output_text.strip()
        if text:
            return text

    segments: list[str] = []
    for item in getattr(response, "output", []) or []:
        for content in getattr(item, "content", []) or []:
            text = getattr(content, "text", None)
            if text:
                segments.append(str(text))
    return "".join(segments).strip()


def send_to_llm(prompt: str, *, model: str) -> tuple[str, dict[str, int | None] | None]:
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY environment variable is not set.")
    from openai import OpenAI  # type: ignore

    client = OpenAI(api_key=api_key)
    response = client.responses.create(
        model=model,
        temperature=1,
        input=prompt,
    )
    message = _extract_output_text(response)
    if not message:
        raise RuntimeError("OpenAI response did not contain any content.")
    usage_summary: dict[str, int | None] | None = None
    usage = getattr(response, "usage", None)
    if usage is not None:
        usage_summary = {
            "prompt_tokens": getattr(usage, "input_tokens", None),
            "completion_tokens": getattr(usage, "output_tokens", None),
            "total_tokens": getattr(usage, "total_tokens", None),
        }
    return message, usage_summary


def derive_output_path(
    input_dir: Path, *, summary_mode: str, output_suffix: str = OUTPUT_SUFFIX
) -> Path:
    base_name = extract_patient_identifier(input_dir)
    if summary_mode not in SUMMARY_MODE_TO_OUTPUT_STEM:
        raise ValueError(f"Unknown summary mode for output naming: {summary_mode}")
    if output_suffix:
        filename = f"{base_name}_{output_suffix}.json"
    else:
        filename = f"{base_name}_{SUMMARY_MODE_TO_OUTPUT_STEM[summary_mode]}.json"
    destination_dir = input_dir / SUMMARY_OUTPUT_DIR_NAME
    destination_dir.mkdir(parents=True, exist_ok=True)
    return destination_dir / filename


def save_result(payload: dict, destination: Path) -> Path:
    destination.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return destination


def save_failed_output(raw_output: str, input_dir: Path, *, output_suffix: str = OUTPUT_SUFFIX) -> Path:
    FAILED_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    identifier = sanitize_filename(f"{input_dir.name}_{output_suffix}")
    destination = FAILED_DIR / f"{identifier}_{timestamp}.txt"
    destination.write_text(raw_output, encoding="utf-8")
    return destination


def process_directory(
    input_dir: Path,
    *,
    model: str,
    prompt_path: Path | None = None,
    summary_mode: str = "auto",
    output_suffix: str = OUTPUT_SUFFIX,
    ranked_mngs_path: Path | None = DEFAULT_RANKED_MNGS_PATH,
    mngs_source: str = DEFAULT_MNGS_SOURCE,
    overwrite: bool = False,
) -> None:
    final_summary_path = find_final_summary_file(input_dir, summary_mode=summary_mode)
    if final_summary_path is None:
        logging.warning("Skipping %s; final_summary JSON not found.", input_dir)
        return

    mngs_path = None if mngs_source == "ranked_only" else find_mngs_source_file(
        input_dir,
        mngs_source=mngs_source,
    )
    if mngs_source != "ranked_only" and mngs_path is None:
        logging.warning("Skipping %s; mNGS source JSON not found for %s.", input_dir, mngs_source)
        return
    output_path = derive_output_path(
        input_dir,
        summary_mode=summary_mode,
        output_suffix=output_suffix,
    )
    if output_path.exists() and not overwrite and not confirm_overwrite(
        output_path, prompt=f"File {output_path} exists. Overwrite?"
    ):
        logging.info("Skipped %s; output already exists: %s", input_dir.name, output_path)
        return

    try:
        selected_prompt = prompt_path or SUMMARY_MODE_TO_PROMPT_PATH.get(summary_mode, PROMPT_PATH)
        template = load_prompt(selected_prompt)
        final_summary = load_json(final_summary_path)
        ranked_mngs = load_ranked_mngs_for_patient(ranked_mngs_path, input_dir)
        if mngs_source == "ranked_only":
            if ranked_mngs is None:
                logging.warning("Skipping %s; ranked mNGS JSON not found.", input_dir)
                return
            mngs_grouped = None
        else:
            mngs_grouped = load_json(mngs_path)
        prompt = build_prompt(
            template,
            final_summary=final_summary,
            mngs_grouped=mngs_grouped,
            ranked_mngs=ranked_mngs,
            final_path=final_summary_path,
            mngs_path=mngs_path,
            ranked_path=ranked_mngs_path,
        )
        estimated_tokens = estimate_token_count(prompt, model=model)
        if estimated_tokens is not None:
            logging.info("Estimated input tokens for %s: %d", input_dir.name, estimated_tokens)
        else:
            logging.debug("Token estimation unavailable (tiktoken missing) for %s", input_dir.name)

        raw_output, usage = send_to_llm(prompt, model=model)
        if usage is not None:
            logging.info(
                "LLM token usage - prompt: %s, completion: %s, total: %s",
                usage.get("prompt_tokens"),
                usage.get("completion_tokens"),
                usage.get("total_tokens"),
            )
        else:
            logging.debug("Token usage data not returned for %s", input_dir.name)

        parsed = validate_json(raw_output)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.exception("Failed while processing %s: %s", input_dir, exc)
        return

    if parsed is None:
        failed_suffix = output_suffix or SUMMARY_MODE_TO_OUTPUT_STEM[summary_mode]
        failed_path = save_failed_output(raw_output, input_dir, output_suffix=failed_suffix)
        logging.error("Stored unparsed output for %s at %s", input_dir.name, failed_path)
        return

    parsed = postprocess_mngs_max_output(
        parsed,
        ranked_mngs=ranked_mngs,
        mngs_grouped=mngs_grouped,
        final_summary=final_summary,
    )
    output_path = save_result(parsed, output_path)
    logging.info("Wrote %s", output_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine final_summary and mngs_name_reads JSON with the big prompt and send to GPT."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Directories containing patient data (patient folders or their parent directories).",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Model identifier to use when calling the GPT API (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--output-suffix",
        default=OUTPUT_SUFFIX,
        help=(
            "Output filename suffix override. "
            "Default uses clear mode names: "
            "'..._mNGS_max_agent_full.json' and "
            "'..._mNGS_max_agent_no_filmarray.json'. "
            "Set this option only when you want custom naming."
        ),
    )
    parser.add_argument(
        "--prompt",
        type=Path,
        help=(
            "Override prompt path. "
            f"Default depends on --summary-mode: full={FULL_PROMPT_PATH}, "
            f"no_filmarray={NO_FILMARRAY_PROMPT_PATH}."
        ),
    )
    parser.add_argument(
        "--ranked-mngs",
        type=Path,
        default=DEFAULT_RANKED_MNGS_PATH,
        help=(
            "Optional ranked mNGS JSON to include in the prompt "
            f"(default: {DEFAULT_RANKED_MNGS_PATH}; ignored if missing)."
        ),
    )
    parser.add_argument(
        "--mngs-source",
        choices=MNGS_SOURCE_CHOICES,
        default=DEFAULT_MNGS_SOURCE,
        help=(
            "mNGS supplemental source for max prompt. "
            "grouped uses *_mNGS_grouped.json; all_rk_ntc uses *_all_RK_NTC_microbes.json; "
            "ranked_only uses only *_mNGS_ranked_candidates.json plus final_summary."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    parser.add_argument(
        "--summary-mode",
        choices=["auto", "full", "deterministic", "no_filmarray"],
        default="auto",
        help=(
            "Choose which final_summary variant to load. "
            "auto runs both full and no_filmarray; deterministic uses "
            "*_final_summary_with_filmarray_deterministic.json."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing mNGS max output files without prompting.",
    )
    args = parser.parse_args()
    configure_logging(args.verbose)

    target_dirs = collect_patient_directories(
        args.inputs,
        summary_mode=args.summary_mode,
        mngs_source=args.mngs_source,
    )
    if not target_dirs:
        logging.error("No patient directories found in the provided inputs.")
        return

    logging.info("Found %d patient directories to process.", len(target_dirs))
    modes_to_run = ["full", "no_filmarray"] if args.summary_mode == "auto" else [args.summary_mode]
    for directory in target_dirs:
        for mode in modes_to_run:
            logging.info("Processing %s (summary-mode=%s)", directory, mode)
            process_directory(
                directory,
                model=args.model,
                prompt_path=args.prompt,
                summary_mode=mode,
                output_suffix=args.output_suffix,
                ranked_mngs_path=args.ranked_mngs,
                mngs_source=args.mngs_source,
                overwrite=args.overwrite,
            )


if __name__ == "__main__":
    main()

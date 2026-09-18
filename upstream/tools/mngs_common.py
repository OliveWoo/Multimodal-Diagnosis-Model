"""Shared deterministic mNGS I/O and candidate helpers.

This module contains the reusable, model-free utilities that were historically
embedded in ``mngs_big_agent``.  Keeping them here lets the deterministic v19
and v20 scorers, summary builder, chosen-ranked builders, and review tooling run
without retaining the retired LLM mNGS-max program.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import pathogen_normalization as pathogen_names
from utils import sanitize_filename


DEFAULT_RANKED_MNGS_PATH = Path("mngs_candidate_microbes_ranked.json")
LOCAL_CHOSEN_RANKED_MNGS_NAME_TEMPLATE = "{base_name}_mNGS_ranked_candidates_chosen.json"
LOCAL_RANKED_MNGS_NAME_TEMPLATE = "{base_name}_mNGS_ranked_candidates.json"
MNGS_SOURCE_CHOICES = ("grouped", "all_rk_ntc", "ranked_only")
DEFAULT_MNGS_SOURCE = "grouped"
SUMMARY_DIR_NAME = "summary_outputs"
SUMMARY_OUTPUT_DIR_NAME = "summary_outputs"

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


def extract_patient_identifier(input_dir: Path) -> str:
    base_name = sanitize_filename(input_dir.name)
    if base_name.endswith("_json"):
        base_name = base_name[: -len("_json")] or base_name
    return base_name


def load_json(filepath: Path) -> dict[str, Any]:
    with filepath.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def find_final_summary_file(input_dir: Path, *, summary_mode: str = "auto") -> Path | None:
    if summary_mode not in {"auto", *SUMMARY_MODE_TO_SUFFIX.keys()}:
        raise ValueError(f"Unknown summary mode: {summary_mode}")

    base_name = extract_patient_identifier(input_dir)
    search_dirs = [input_dir / SUMMARY_DIR_NAME, input_dir]
    ordered_modes = ["deterministic", "full"] if summary_mode == "auto" else [summary_mode]
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
    return next((candidate for candidate in candidates if candidate.exists()), None)


def find_all_rk_ntc_file(input_dir: Path) -> Path | None:
    base_name = extract_patient_identifier(input_dir)
    exact = input_dir / f"{base_name}_all_RK_NTC_microbes.json"
    if exact.exists():
        return exact
    matches = sorted(input_dir.glob(f"{base_name}_all_RK_NTC_microbes*.json"))
    return matches[0] if matches else None


def find_mngs_source_file(
    input_dir: Path,
    *,
    mngs_source: str = DEFAULT_MNGS_SOURCE,
) -> Path | None:
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
    candidates: list[Path] = []
    if explicit_ranked_path:
        candidates.append(ranked_path)
    candidates.extend(local_candidates)
    if ranked_path is not None and not explicit_ranked_path:
        candidates.append(ranked_path)
    return next((candidate for candidate in candidates if candidate.exists()), None)


def load_ranked_mngs_for_patient(
    ranked_path: Path | None,
    input_dir: Path,
) -> dict[str, Any] | None:
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


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


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


def normalize_number_value(value: Any) -> int | float | str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    return int(numeric) if numeric.is_integer() else numeric


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


def ranked_candidate_sort_key(
    entry: tuple[int, dict[str, Any], dict[str, Any]],
) -> tuple[Any, ...]:
    record_index, record, candidate = entry
    ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
    return (
        to_rank_int(ranking.get("rank_priority")),
        -to_float(candidate.get("reads", candidate.get("sec_hit", 0))),
        -to_float(ranking.get("reads_percentile")),
        -to_float(ranking.get("final_score")),
        str(record.get("specimen_code") or ""),
        normalize_organism_name(candidate.get("organism_name") or candidate.get("name")),
        record_index,
    )


def deduplicate_ranked_mngs_payload(payload: dict[str, Any]) -> dict[str, Any]:
    """Keep one deterministic best ranked mNGS record per organism."""
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
        candidate = entry[2]
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
        normalized_record["status"] = (
            "parsed" if candidates else source_record.get("status", "parsed_no_match")
        )
        normalized_records.append(normalized_record)

    ranking_metadata = dict(payload.get("ranking_metadata") or {})
    ranking_metadata["deduplication"] = {
        "enabled": True,
        "method": "exact_duplicate_removed_then_best_record_per_normalized_organism",
        "input_candidate_rows": len(flattened),
        "output_candidate_rows": sum(
            len(record.get("candidates") or []) for record in normalized_records
        ),
        "exact_duplicate_rows_removed": len(exact_keys) - len(set(exact_keys)),
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


def _level_rank(value: Any) -> int:
    text = str(value or "").strip().lower()
    for index in range(1, 6):
        if text in {f"level {index}", f"level{index}"}:
            return index
    return 99


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

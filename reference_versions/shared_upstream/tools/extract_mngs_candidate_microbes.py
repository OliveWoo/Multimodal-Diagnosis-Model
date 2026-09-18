from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

DEFAULT_MAPPING_PATH = Path("outputs/patient_mapping.json")
DEFAULT_OUTPUT_PATH = Path("outputs/mngs_candidate_microbes.json")
DEFAULT_SPECIMEN_EXCEL_DIR = Path(r"C:\Users\User\Desktop\patient_info")
TARGET_SUBFOLDERS = ("1.Bac", "2.Fungi", "3.Virus")
STATUS_PARSED = "parsed"
STATUS_PARSED_NO_MATCH = "parsed_no_match"
STATUS_TSV_NOT_FOUND = "tsv_not_found"
STATUS_CHIP_FOLDER_NOT_FOUND = "chip_folder_not_found"
STATUS_TSV_READ_ERROR = "tsv_read_error"
STATUS_SEQ_ID_NOT_FOUND_IN_EXCEL = "seq_id_not_found_in_excel"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    default_specimen_excel = _resolve_default_specimen_excel()
    parser = argparse.ArgumentParser(
        description=(
            "Extract candidate microbes from mNGS TSV files using patient_mapping.json and "
            "the authoritative seq_id from the specimen Excel. This step does not call any LLM."
        )
    )
    parser.add_argument(
        "--mapping-json",
        type=Path,
        default=DEFAULT_MAPPING_PATH,
        help=f"Path to patient_mapping.json (default: {DEFAULT_MAPPING_PATH})",
    )
    parser.add_argument(
        "--specimen-excel",
        type=Path,
        default=default_specimen_excel,
        help=f"Path to the specimen Excel file (default: {default_specimen_excel})",
    )
    parser.add_argument(
        "--root-dir",
        type=Path,
        required=True,
        help="Root directory containing chip folders.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help=f"Path to write the result JSON (default: {DEFAULT_OUTPUT_PATH})",
    )
    return parser.parse_args(argv)


def _load_mapping(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected mapping JSON to be a dict, got {type(payload)}")
    return payload


def _normalize_header(value: Any) -> str:
    return "".join(str(value).split()).lower()


def _resolve_default_specimen_excel() -> Path:
    xlsx_files = sorted(
        path
        for path in DEFAULT_SPECIMEN_EXCEL_DIR.glob("*.xlsx")
        if not path.name.startswith("~$")
    )
    if not xlsx_files:
        return DEFAULT_SPECIMEN_EXCEL_DIR / "specimen_lookup.xlsx"
    return xlsx_files[0]


def _detect_column(
    columns: Sequence[Any],
    accepted: set[str],
    label: str,
    fallback_index: int | None = None,
) -> Any:
    normalized_map = {_normalize_header(column): column for column in columns}
    for candidate in accepted:
        if candidate in normalized_map:
            return normalized_map[candidate]
    if fallback_index is not None and 0 <= fallback_index < len(columns):
        return columns[fallback_index]
    raise KeyError(f"Unable to find column for {label}. Checked: {sorted(accepted)}")


def _coerce_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    if text.endswith(".0") and text.replace(".", "", 1).isdigit():
        return text[:-2]
    return text


def _load_specimen_seq_map(path: Path) -> dict[str, list[str]]:
    frame = pd.read_excel(path, dtype=object)
    specimen_code_col = _detect_column(
        frame.columns,
        {"specimencode", "specimen_code", "samplecode", "sample_code"},
        "specimen code",
        fallback_index=9,
    )
    seq_id_col = _detect_column(
        frame.columns,
        {"seq_id", "seqid", "systemcode", "system_id"},
        "seq_id",
        fallback_index=10,
    )

    specimen_to_seq_ids: dict[str, list[str]] = defaultdict(list)
    for _, row in frame.iterrows():
        specimen_code = _coerce_text(row[specimen_code_col])
        seq_id = _coerce_text(row[seq_id_col])
        if not specimen_code or not seq_id:
            continue
        existing = specimen_to_seq_ids[specimen_code]
        if seq_id not in existing:
            existing.append(seq_id)
    return dict(specimen_to_seq_ids)


def _find_chip_folders(root_dir: Path, chip_id: str) -> list[Path]:
    matches: list[Path] = []
    needle = chip_id.strip().lower()
    for path in root_dir.rglob("*"):
        if not path.is_dir():
            continue
        if needle in path.name.lower():
            matches.append(path)
    return sorted(matches)


def _looks_like_system_code(value: Any) -> bool:
    normalized = _coerce_text(value).upper()
    return normalized.startswith(("23S", "24S", "25S", "26S"))


def _find_seq_tsv_files(chip_folder: Path, seq_id: str) -> list[Path]:
    matches: list[Path] = []
    needle = seq_id.strip().lower()
    for subfolder_name in TARGET_SUBFOLDERS:
        subfolder = chip_folder / subfolder_name
        if not subfolder.is_dir():
            continue
        for candidate in subfolder.rglob("*.tsv"):
            if needle in candidate.name.lower():
                matches.append(candidate)
    return sorted(matches)


def _find_seq_tsv_files_anywhere(root_dir: Path, seq_id: str) -> list[Path]:
    matches: list[Path] = []
    needle = seq_id.strip().lower()
    for candidate in root_dir.rglob("*.tsv"):
        if needle in candidate.name.lower():
            matches.append(candidate)
    return sorted(matches)


def _match_code_ntc(value: Any) -> bool:
    normalized = _coerce_text(value)
    if not normalized:
        return False
    try:
        return int(normalized) == 0
    except ValueError:
        return normalized == "0"


def _match_code_rk_ntc(value: Any) -> bool:
    normalized = _coerce_text(value).upper()
    if not normalized:
        return False
    if normalized.startswith("RK-"):
        normalized = normalized[3:]
    try:
        numeric = int(normalized)
    except ValueError:
        return normalized in {"RK-0", "RK-8", "0", "8"}
    return numeric in {0, 8}


def _coerce_sec_hit(value: Any) -> Any:
    text = _coerce_text(value)
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return text
    if number.is_integer():
        return int(number)
    return number


def _extract_candidates(tsv_path: Path) -> list[dict[str, Any]]:
    frame = pd.read_csv(tsv_path, sep="\t", dtype=object)
    required_columns = {"Code_NTC", "Code_RK_NTC", "Organism_name", "Sec.hit"}
    missing_columns = required_columns - set(frame.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise KeyError(f"Missing required columns in {tsv_path}: {missing}")

    matched_rows = frame[
        frame["Code_NTC"].apply(_match_code_ntc) | frame["Code_RK_NTC"].apply(_match_code_rk_ntc)
    ]

    candidates: list[dict[str, Any]] = []
    for _, row in matched_rows.iterrows():
        candidates.append(
            {
                "organism_name": _coerce_text(row["Organism_name"]),
                "sec_hit": _coerce_sec_hit(row["Sec.hit"]),
            }
        )
    return candidates


def _parse_record(chip_ids: Sequence[str], seq_id: str, root_dir: Path) -> dict[str, Any]:
    result = {"status": "", "candidates": []}

    chip_folders: list[Path] = []
    for chip_id in chip_ids:
        chip_folders.extend(_find_chip_folders(root_dir, chip_id))
    chip_folders = sorted({path for path in chip_folders})
    has_system_code = any(_looks_like_system_code(chip_id) for chip_id in chip_ids)
    if not chip_folders and not has_system_code:
        result["status"] = STATUS_CHIP_FOLDER_NOT_FOUND
        return result

    matched_tsv_files: list[Path] = []
    for chip_folder in chip_folders:
        matched_tsv_files.extend(_find_seq_tsv_files(chip_folder, seq_id))
    if not matched_tsv_files and has_system_code:
        matched_tsv_files.extend(_find_seq_tsv_files_anywhere(root_dir, seq_id))
    matched_tsv_files = sorted({path for path in matched_tsv_files})
    if not matched_tsv_files:
        result["status"] = STATUS_TSV_NOT_FOUND
        return result

    read_error = False
    candidates: list[dict[str, Any]] = []
    for tsv_path in matched_tsv_files:
        try:
            extracted = _extract_candidates(tsv_path)
        except Exception:
            read_error = True
            continue

        source_category = tsv_path.parent.name
        for item in extracted:
            candidates.append(
                {
                    "organism_name": item["organism_name"],
                    "sec_hit": item["sec_hit"],
                    "source_category": source_category,
                }
            )

    result["candidates"] = candidates
    if candidates:
        result["status"] = STATUS_PARSED
    elif read_error:
        result["status"] = STATUS_TSV_READ_ERROR
    else:
        result["status"] = STATUS_PARSED_NO_MATCH
    return result


def build_result(
    mapping: dict[str, Any],
    specimen_to_seq_ids: dict[str, list[str]],
    root_dir: Path,
) -> dict[str, Any]:
    patients: dict[str, Any] = {}
    for patient_id, patient_entry in mapping.items():
        original_records = patient_entry.get("records", [])
        chip_ids_by_specimen: dict[str, list[str]] = defaultdict(list)
        for record in original_records:
            specimen_code = _coerce_text(record.get("specimen_code"))
            chip_id = _coerce_text(record.get("chip_id"))
            if specimen_code and chip_id and chip_id not in chip_ids_by_specimen[specimen_code]:
                chip_ids_by_specimen[specimen_code].append(chip_id)

        records: list[dict[str, Any]] = []
        specimen_codes = patient_entry.get("specimen_codes", [])
        for specimen_code_raw in specimen_codes:
            specimen_code = _coerce_text(specimen_code_raw)
            seq_ids = specimen_to_seq_ids.get(specimen_code, [])
            chip_ids = chip_ids_by_specimen.get(specimen_code, [])

            if not seq_ids:
                records.append(
                    {
                        "specimen_code": specimen_code,
                        "chip_ids": chip_ids,
                        "seq_id": None,
                        "status": STATUS_SEQ_ID_NOT_FOUND_IN_EXCEL,
                        "candidates": [],
                    }
                )
                continue

            for seq_id in seq_ids:
                parsed = _parse_record(chip_ids=chip_ids, seq_id=seq_id, root_dir=root_dir)
                records.append(
                    {
                        "specimen_code": specimen_code,
                        "chip_ids": chip_ids,
                        "seq_id": seq_id,
                        **parsed,
                    }
                )

        patients[patient_id] = {"records": records}

    return {"patients": patients}


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    mapping = _load_mapping(args.mapping_json)
    specimen_to_seq_ids = _load_specimen_seq_map(args.specimen_excel)
    result = build_result(mapping, specimen_to_seq_ids, args.root_dir)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    print(f"Wrote result JSON to {args.output_json}")


if __name__ == "__main__":
    main()

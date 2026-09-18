from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd

DEFAULT_PATIENT_ROOT = Path(r"C:\Users\User\Desktop\patient_info\0310_data")
DEFAULT_MAPPING_JSON = Path("outputs/patient_mapping.json")
DEFAULT_OUTPUT_JSON = Path("outputs/mngs_max_level123_overlap_with_excel.json")
DEFAULT_SPECIMEN_EXCEL_DIR = Path(r"C:\Users\User\Desktop\patient_info")
TARGET_LEVELS = {"Level 1", "Level 2", "Level 3"}

VALUE_PATTERN = re.compile(r"^(?P<name>.+?)\s*\((?P<count>[^()]+)\)\s*$")
NAME_SANITIZER = re.compile(r"[^a-z0-9]+")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare Level 1-3 pathogens from mNGS_max_agent_full.json with suspected pathogens "
            "in the specimen Excel and report overlap per patient."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_PATIENT_ROOT,
        help="Path to a patient root directory or a single mNGS_max_agent_full.json file.",
    )
    parser.add_argument(
        "--mapping-json",
        type=Path,
        default=DEFAULT_MAPPING_JSON,
        help=f"Path to patient_mapping.json (default: {DEFAULT_MAPPING_JSON})",
    )
    parser.add_argument(
        "--specimen-excel",
        type=Path,
        default=None,
        help="Path to 日期及菌種對照.xlsx. Defaults to the first .xlsx in patient_info.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=DEFAULT_OUTPUT_JSON,
        help=f"Path to write comparison JSON (default: {DEFAULT_OUTPUT_JSON})",
    )
    return parser.parse_args(argv)


def _resolve_default_specimen_excel() -> Path:
    xlsx_files = sorted(
        path
        for path in DEFAULT_SPECIMEN_EXCEL_DIR.glob("*.xlsx")
        if not path.name.startswith("~$")
    )
    if not xlsx_files:
        return DEFAULT_SPECIMEN_EXCEL_DIR / "specimen_lookup.xlsx"
    return xlsx_files[0]


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def _normalize_header(value: Any) -> str:
    return "".join(str(value).split()).lower()


def _lookup_column(normalized_map: dict[str, str], candidates: Iterable[str]) -> str:
    for candidate in candidates:
        key = _normalize_header(candidate)
        if key in normalized_map:
            return normalized_map[key]
    raise KeyError(f"Unable to locate any of: {', '.join(candidates)}")


def _normalize_specimen_code(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    if text.endswith(".0") and text.replace(".", "", 1).isdigit():
        return text[:-2]
    return text


def _normalize_name(name: Any) -> str:
    text = str(name or "").strip()
    if not text:
        return ""
    # Remove parenthetical remarks so "Staphylococcus aureus (MRSA likely)" can
    # match Excel entries like "Staphylococcus aureus (123)".
    text = re.sub(r"\([^)]*\)", "", text).strip()
    return NAME_SANITIZER.sub("", text.lower())


def _parse_excel_pathogen_cell(value: Any) -> list[str]:
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text == "-":
        return []

    names: list[str] = []
    for raw in text.replace("\r", "\n").split("\n"):
        candidate = raw.strip()
        if not candidate or candidate == "-":
            continue
        match = VALUE_PATTERN.match(candidate)
        if match:
            candidate = match.group("name").strip()
        names.append(candidate)
    return names


def _collect_input_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        return [input_path]
    return sorted(input_path.rglob("*_mNGS_max_agent_full.json"))


def _extract_patient_id(path: Path) -> str:
    match = re.search(r"NGS_patient_(\d+)_json", str(path))
    if not match:
        raise ValueError(f"Unable to determine patient id from path: {path}")
    return match.group(1)


def _load_patient_specimens(path: Path) -> dict[str, list[str]]:
    payload = _load_json(path)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict patient mapping, got {type(payload)}")
    result: dict[str, list[str]] = {}
    for patient_id, patient_entry in payload.items():
        if not isinstance(patient_entry, dict):
            continue
        specimen_codes = patient_entry.get("specimen_codes", [])
        result[str(patient_id)] = [_normalize_specimen_code(code) for code in specimen_codes if _normalize_specimen_code(code)]
    return result


def _load_excel_suspected_lookup(path: Path) -> dict[str, list[dict[str, str]]]:
    frame = pd.read_excel(path, dtype=object)
    normalized_map = {_normalize_header(column): column for column in frame.columns}
    specimen_col = _lookup_column(normalized_map, ["檢體代碼", "specimen code", "specimen_code"])
    bac_col = _lookup_column(normalized_map, ["細菌-疑似致病菌"])
    fungi_col = _lookup_column(normalized_map, ["真菌-疑似致病菌"])
    virus_col = _lookup_column(normalized_map, ["病毒-疑似致病菌"])

    lookup: dict[str, list[dict[str, str]]] = {}
    for _, row in frame[[specimen_col, bac_col, fungi_col, virus_col]].iterrows():
        specimen_code = _normalize_specimen_code(row[specimen_col])
        if not specimen_code:
            continue

        entries = lookup.setdefault(specimen_code, [])
        for raw_name in _parse_excel_pathogen_cell(row[bac_col]):
            entries.append({"name": raw_name, "type": "bacteria"})
        for raw_name in _parse_excel_pathogen_cell(row[fungi_col]):
            entries.append({"name": raw_name, "type": "fungi"})
        for raw_name in _parse_excel_pathogen_cell(row[virus_col]):
            entries.append({"name": raw_name, "type": "virus"})
    return lookup


def _extract_level123_pathogens(payload: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = payload.get("pathogen_candidates", [])
    if not isinstance(candidates, list):
        return []

    extracted: list[dict[str, Any]] = []
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        level = str(candidate.get("integrated_causative_level") or "").strip()
        if level not in TARGET_LEVELS:
            continue
        classification = str(candidate.get("classification") or "").strip().lower()
        candidate_type = {
            "bacterial": "bacteria",
            "fungal": "fungi",
            "viral": "virus",
        }.get(classification, classification)
        extracted.append(
            {
                "name": candidate.get("organism_name"),
                "type": candidate_type,
                "level": level,
                "confidence": candidate.get("integrated_causative_confidence"),
                "normalized_name": _normalize_name(candidate.get("organism_name")),
            }
        )
    return extracted


def _dedupe_entries(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[tuple[str, str]] = set()
    deduped: list[dict[str, Any]] = []
    for entry in entries:
        key = (str(entry.get("normalized_name") or ""), str(entry.get("type") or ""))
        if not key[0] or key in seen:
            continue
        seen.add(key)
        deduped.append(entry)
    return deduped


def _build_overlap_payload(
    patient_id: str,
    mngs_entries: list[dict[str, Any]],
    excel_entries: list[dict[str, Any]],
    specimen_codes: list[str],
) -> dict[str, Any]:
    deduped_mngs = _dedupe_entries(mngs_entries)
    excel_with_norm = [
        {
            "name": entry["name"],
            "type": entry["type"],
            "normalized_name": _normalize_name(entry["name"]),
        }
        for entry in excel_entries
    ]
    deduped_excel = _dedupe_entries(excel_with_norm)

    excel_key_map = {
        (entry["normalized_name"], entry["type"]): entry
        for entry in deduped_excel
        if entry["normalized_name"]
    }
    overlaps: list[dict[str, Any]] = []
    for entry in deduped_mngs:
        key = (entry["normalized_name"], entry["type"])
        matched = excel_key_map.get(key)
        if matched is None:
            continue
        overlaps.append(
            {
                "name": entry["name"],
                "type": entry["type"],
                "level": entry["level"],
                "excel_name": matched["name"],
            }
        )

    llm_count = len(deduped_mngs)
    excel_count = len(deduped_excel)
    overlap_count = len(overlaps)
    return {
        "patient_id": patient_id,
        "specimen_codes": specimen_codes,
        "mngs_level_1_3": [
            {
                "name": entry["name"],
                "type": entry["type"],
                "level": entry["level"],
                "confidence": entry["confidence"],
            }
            for entry in deduped_mngs
        ],
        "excel_suspected_pathogens": [
            {
                "name": entry["name"],
                "type": entry["type"],
            }
            for entry in deduped_excel
        ],
        "overlaps": overlaps,
        "overlap_count": overlap_count,
        "mngs_level_1_3_count": llm_count,
        "excel_suspected_count": excel_count,
        "overlap_ratio_vs_mngs": overlap_count / llm_count if llm_count else 0.0,
        "overlap_ratio_vs_excel": overlap_count / excel_count if excel_count else 0.0,
        "jaccard": overlap_count / (llm_count + excel_count - overlap_count)
        if (llm_count + excel_count - overlap_count)
        else 0.0,
    }


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    specimen_excel = args.specimen_excel or _resolve_default_specimen_excel()

    patient_specimens = _load_patient_specimens(args.mapping_json)
    excel_lookup = _load_excel_suspected_lookup(specimen_excel)
    input_files = _collect_input_files(args.input)

    comparisons: dict[str, Any] = {}
    for file_path in input_files:
        patient_id = _extract_patient_id(file_path)
        payload = _load_json(file_path)
        if not isinstance(payload, dict):
            continue

        specimen_codes = patient_specimens.get(patient_id, [])
        excel_entries: list[dict[str, Any]] = []
        for specimen_code in specimen_codes:
            excel_entries.extend(excel_lookup.get(specimen_code, []))

        mngs_entries = _extract_level123_pathogens(payload)
        comparisons[patient_id] = _build_overlap_payload(
            patient_id=patient_id,
            mngs_entries=mngs_entries,
            excel_entries=excel_entries,
            specimen_codes=specimen_codes,
        )

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump({"patients": comparisons}, handle, ensure_ascii=False, indent=2)
    print(f"Wrote overlap report to {args.output_json}")


if __name__ == "__main__":
    main()

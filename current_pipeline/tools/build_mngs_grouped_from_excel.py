from __future__ import annotations

"""
Build specimen_lookup_summary_mNGS_grouped.json directly from a B Excel file.

The script:
1) reads the Excel (all sheets),
2) maps each 檢體代碼 to patient_mapping.json (with the BALF A/B exception),
3) groups pathogens into bacterial / viral / fungal / others, and
4) writes outputs/specimen_lookup_summary_mNGS_grouped.json (same schema as the existing file).
"""

import argparse
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from parsers.mNGS_parser import transform_mngs_payload
from tools import culture_bal_extractor as bal

BALF_LABELS = {"balf", "lower bal"}

SPECIMEN_CODE_HEADERS = {
    bal.SPECIMEN_CODE_COLUMN,
    "檢體代碼",
    "檢體代號",
    "specimen code",
    "specimen_code",
    "sample id",
    "sample_id",
}
SPECIMEN_TYPE_HEADERS = {
    "檢體類型",
    "檢體種類",
    "specimen type",
    "specimen_type",
    "sample type",
    "sample_type",
}
SPECIMEN_SITE_HEADERS = {
    "檢體部位",
    "檢體類型",
    "specimen site",
    "specimen_site",
    "sample site",
    "sample_site",
}
COLLECTED_TIME_HEADERS = {
    bal.COLLECTED_TIME_COLUMN,
    "採檢時間",
    "採檢日期",
    "採檢日期時間",
    "collection time",
    "collected time",
    "collected_time",
}
CATEGORY_HEADER_CANDIDATES = {
    "bacterial": {bal.BACTERIA_COLUMN, "細菌", "bacteria", "bacterial"},
    "viral": {bal.VIRUS_COLUMN, "病毒", "virus", "viral"},
    "fungal": {bal.FUNGUS_COLUMN, "真菌", "fungus", "fungal"},
}
CATEGORY_LABELS = {field[2]: field[3] for field in bal.SUSPECTED_FIELDS}

DEFAULT_PATIENT_MAPPING = Path("outputs") / "patient_mapping.json"
DEFAULT_OUTPUT = Path("outputs") / "specimen_lookup_summary_mNGS_grouped.json"


@dataclass
class SpecimenRecord:
    specimen_code: str
    base_code: str
    collected_time: str
    specimen_type: str
    specimen_site: str
    suspected_pathogens: list[dict[str, Any]]


def _normalize_header(label: Any) -> str:
    return "".join(str(label or "").split()).lower()


def _extract_base_code(code: str) -> str:
    return "".join(ch for ch in code if ch.isdigit())


def _is_balf(sample_type: str) -> bool:
    normalized = _normalize_header(sample_type)
    return any(normalized == _normalize_header(label) for label in BALF_LABELS)


def _pick_column(columns: Iterable[str], candidates: Iterable[str], label: str, *, required: bool = True) -> str | None:
    normalized_map = {_normalize_header(col): col for col in columns}
    for candidate in candidates:
        key = _normalize_header(candidate)
        if key in normalized_map:
            return normalized_map[key]
    if required:
        raise KeyError(f"Missing column {label}. Tried: {', '.join(map(str, candidates))}")
    logging.warning("Optional column %s not found; tried %s", label, ", ".join(map(str, candidates)))
    return None


def _resolve_columns(columns: Sequence[str]) -> dict[str, Any]:
    specimen_col = _pick_column(columns, SPECIMEN_CODE_HEADERS, "[檢體代碼]")
    collected_col = _pick_column(columns, COLLECTED_TIME_HEADERS, "[採檢時間]")
    specimen_type_col = _pick_column(columns, SPECIMEN_TYPE_HEADERS, "[檢體類型]", required=False)
    specimen_site_col = _pick_column(columns, SPECIMEN_SITE_HEADERS, "[檢體部位]", required=False)

    category_cols: dict[str, str | None] = {}
    for key, candidates in CATEGORY_HEADER_CANDIDATES.items():
        category_cols[key] = _pick_column(columns, candidates, f"{key} pathogens", required=False)

    return {
        "specimen_code": specimen_col,
        "collected_time": collected_col,
        "specimen_type": specimen_type_col,
        "specimen_site": specimen_site_col,
        "categories": category_cols,
    }


def _load_excel_rows(excel_path: Path) -> list[SpecimenRecord]:
    df = bal._read_excel_all_sheets(excel_path)
    if df is None or df.empty:
        return []

    column_map = _resolve_columns(list(df.columns))
    category_columns = column_map["categories"]

    records: list[SpecimenRecord] = []
    for _, row in df.iterrows():
        code = bal._normalize_specimen_code(row.get(column_map["specimen_code"]))
        if not code:
            continue

        specimen_type = bal._coerce_cell(row.get(column_map["specimen_type"])) if column_map["specimen_type"] else ""
        specimen_site = bal._coerce_cell(row.get(column_map["specimen_site"])) if column_map["specimen_site"] else ""
        suspected: list[dict[str, Any]] = []
        for category_key, column_name in category_columns.items():
            if not column_name:
                continue
            raw_value = bal._coerce_cell(row.get(column_name))
            parsed = bal._parse_suspected_entries(
                raw_value,
                category_key=category_key,
                category_label=CATEGORY_LABELS.get(category_key, category_key),
            )
            if parsed:
                suspected.extend(parsed)

        record = SpecimenRecord(
            specimen_code=code,
            base_code=_extract_base_code(code),
            collected_time=bal._format_datetime(row.get(column_map["collected_time"])),
            specimen_type=specimen_type,
            specimen_site=specimen_site,
            suspected_pathogens=suspected,
        )
        records.append(record)

    return records


def _build_patient_indexes(mapping: Mapping[str, Any]) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    exact: dict[str, set[str]] = {}
    base: dict[str, set[str]] = {}

    for patient_id, payload in mapping.items():
        codes = payload.get("specimen_codes") or []
        for code in codes:
            normalized = bal._normalize_specimen_code(code)
            if not normalized:
                continue
            exact.setdefault(normalized, set()).add(str(patient_id))
            base_code = _extract_base_code(normalized)
            if base_code:
                base.setdefault(base_code, set()).add(str(patient_id))

    return exact, base


def _resolve_patient_for_row(
    index: int,
    rows: Sequence[SpecimenRecord],
    exact_index: Mapping[str, set[str]],
    base_index: Mapping[str, set[str]],
) -> str | None:
    row = rows[index]
    patient_ids = exact_index.get(row.specimen_code)

    if not patient_ids:
        if not _is_balf(row.specimen_type) and index + 1 < len(rows):
            next_row = rows[index + 1]
            if next_row.base_code and next_row.base_code == row.base_code:
                patient_ids = exact_index.get(next_row.specimen_code) or base_index.get(next_row.base_code)
        if not patient_ids and row.base_code:
            patient_ids = base_index.get(row.base_code)

    if not patient_ids:
        return None
    if len(patient_ids) > 1:
        raise ValueError(f"Specimen code {row.specimen_code} matches multiple patients: {sorted(patient_ids)}")
    return next(iter(patient_ids))


def _build_patient_payload(rows: Sequence[SpecimenRecord], patient_mapping: Mapping[str, Any]) -> dict[str, list[dict[str, Any]]]:
    exact_index, base_index = _build_patient_indexes(patient_mapping)
    grouped: dict[str, list[dict[str, Any]]] = {}

    for idx, row in enumerate(rows):
        patient_id = _resolve_patient_for_row(idx, rows, exact_index, base_index)
        if not patient_id:
            continue

        entry: dict[str, Any] = {
            "specimen_code": row.specimen_code,
            "collected_time": row.collected_time,
        }
        if row.suspected_pathogens:
            entry["suspected_pathogens"] = row.suspected_pathogens
        if row.specimen_site:
            entry["specimen_site"] = row.specimen_site
        grouped.setdefault(str(patient_id), []).append(entry)

    return grouped


def _load_patient_mapping(path: Path) -> Mapping[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build specimen_lookup_summary_mNGS_grouped.json directly from the mNGS Excel.",
    )
    parser.add_argument("excel", type=Path, help="Path to the Excel file containing 檢體代碼 / 檢體類型 / pathogens.")
    parser.add_argument(
        "--patient-mapping",
        type=Path,
        default=DEFAULT_PATIENT_MAPPING,
        help=f"Path to patient_mapping.json (default: {DEFAULT_PATIENT_MAPPING}).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Where to write the grouped JSON (default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument("--indent", type=int, default=2, help="Pretty-print indentation for JSON output.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    excel_path = args.excel.expanduser().resolve()
    patient_mapping_path = args.patient_mapping.expanduser().resolve()
    output_path = args.output.expanduser().resolve()

    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    if not patient_mapping_path.exists():
        raise FileNotFoundError(f"patient_mapping.json not found: {patient_mapping_path}")

    rows = _load_excel_rows(excel_path)
    patient_mapping = _load_patient_mapping(patient_mapping_path)
    raw_payload = _build_patient_payload(rows, patient_mapping)
    grouped_payload = transform_mngs_payload(raw_payload)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(grouped_payload, ensure_ascii=False, indent=args.indent), encoding="utf-8")
    print(f"Wrote grouped mNGS summary to {output_path}")


if __name__ == "__main__":
    main()

"""
Tool for building a patient-to-specimen mapping by chaining three sources:
X = patient/specimen metadata (Excel), Y = 檢體代碼 <-> specimen_id/specimen_id_CaseReview,
Z = case review records keyed by specimen_id_CaseReview. The output schema is identical
to the previous implementation so downstream consumers remain compatible.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import json
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

import pandas as pd

from tools import pathogen_normalization as pathogen_names

# Column aliases are normalized (trim + lowercase) to make the parser resilient to
# header variations introduced by manual editing.
PATIENT_NO_HEADERS = {
    "patient no.",
    "patient no",
    "patient_no",
    "patient number",
    "patientnumber",
    "patient#",
    "patient id",
    "patientid",
    "病人編號",
    "病歷號",
}

SPECIMEN_CODE_HEADERS = {
    "檢體代碼",
    "檢體代號",
    "specimen code",
    "specimen_code",
    "specimen id",
    "specimen_id",
    "sample id",
    "sample_id",
}

SPECIMEN_ID_HEADERS = {
    "specimen_id",
    "specimen id",
    "specimenid",
}

CSV_SEQ_HEADERS = {
    "seq_id",
    "seq id",
    "sequence id",
    "sequence_id",
    "system code",
    "system_code",
    "系統代碼",
    "系統碼",
}
CSV_CHIP_HEADERS = {"chip_id", "chip id"}
BRIDGE_SPECIMEN_HEADERS = SPECIMEN_CODE_HEADERS | SPECIMEN_ID_HEADERS

CASE_REVIEW_ID_HEADERS = {
    "specimen_id_casereview",
    "specimen id casereview",
    "specimen_id_case_review",
    "specimen id case review",
    "specimenid_casereview",
    "specimenidcasereview",
    "specimen_id",
}
EXCEL_SUFFIXES = {".xlsx", ".xlsm", ".xls"}

COLLECTED_TIME_HEADERS = {
    "採檢時間",
    "採檢日期",
    "collection time",
    "collected time",
}

SUSPECTED_COLUMN_LABELS = {
    "bacterial": "細菌-\n疑似致病菌",
    "fungal": "真菌-\n疑似致病菌",
    "viral": "病毒-\n疑似致病菌",
}

CATEGORY_DISPLAY = {
    "bacterial": "細菌",
    "fungal": "真菌",
    "viral": "病毒",
}

NAME_SANITIZER = re.compile(r"[^a-z0-9]+")
VALUE_PATTERN = re.compile(r"^(?P<name>.+?)\s*\((?P<count>[^()]+)\)\s*$")
SPECIMEN_CODE_SPLIT_PATTERN = re.compile(r"[\r\n;,，、]+")


@dataclass(frozen=True)
class PatientSpecimenRow:
    patient_no: str
    specimen_code: str
    system_code: str | None = None


@dataclass(frozen=True)
class BridgeRow:
    specimen_code: str
    specimen_id: str
    specimen_id_case_review: str


@dataclass(frozen=True)
class CaseReviewRow:
    specimen_id_case_review: str
    specimen_id: str | None
    seq_id: str | None
    chip_id: str | None


def _normalize_header(header: str) -> str:
    return header.strip().lower()


def _coerce_value(value) -> str:
    if isinstance(value, pd.Series):
        for v in value:
            if not pd.isna(v):
                value = v
                break
        else:
            return ""
    if pd.isna(value):
        return ""
    return str(value).strip()


def _detect_column(columns: Sequence[str], accepted: Iterable[str], label: str) -> str:
    normalized_map = {_normalize_header(col): col for col in columns}
    for candidate in accepted:
        if candidate in normalized_map:
            return normalized_map[candidate]
    raise KeyError(f"Unable to find column for {label}. Checked: {', '.join(accepted)}")


def _load_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in EXCEL_SUFFIXES:
        return _read_excel_all_sheets(path)
    raise ValueError(
        f"Unsupported file type for {path}. Expected CSV or one of {sorted(EXCEL_SUFFIXES)}."
    )


def _read_excel_all_sheets(path: Path) -> pd.DataFrame:
    frames = pd.read_excel(path, sheet_name=None, engine="openpyxl")
    # Merge every non-empty sheet so multi-sheet workbooks are fully scanned.
    if isinstance(frames, dict):
        non_empty = [frame for frame in frames.values() if frame is not None and not frame.empty]
        if not non_empty:
            return pd.DataFrame()
        return pd.concat(non_empty, ignore_index=True)
    return frames


def _normalize_label(label: Any) -> str:
    return "".join(str(label).split()).lower()


def _normalize_specimen_code(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, (int, float)):
        if isinstance(value, float):
            if pd.isna(value):
                return ""
            if value.is_integer():
                value = int(value)
        return str(value)
    text = str(value).strip()
    if not text:
        return ""
    if text.endswith(".0") and text.replace(".", "").isdigit():
        return text[:-2]
    return text


def _coerce_cell(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, str):
        return value.strip()
    return str(value).strip()


def _split_specimen_codes(value: Any) -> List[str]:
    text = _coerce_cell(value)
    if not text:
        return []

    codes: List[str] = []
    for raw in SPECIMEN_CODE_SPLIT_PATTERN.split(text.replace("\r", "\n")):
        code = _normalize_specimen_code(raw)
        if code and code not in codes:
            codes.append(code)
    return codes


def _format_datetime(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, pd.Timestamp):
        if pd.isna(value):
            return ""
        return value.to_pydatetime().isoformat(sep=" ")
    if isinstance(value, datetime):
        return value.isoformat(sep=" ")
    text = str(value).strip()
    return text


def _normalize_pathogen_name(name: Any) -> str:
    return pathogen_names.canonical_key(name)


def _split_lines(value: str) -> List[str]:
    parts: List[str] = []
    for raw in value.replace("\r", "\n").split("\n"):
        candidate = raw.strip()
        if candidate and candidate != "-":
            parts.append(candidate)
    return parts


def _parse_suspected_field(value: Any) -> List[Dict[str, Any]]:
    text = _coerce_cell(value)
    if not text or text == "-":
        return []

    entries: List[Dict[str, Any]] = []
    for raw in _split_lines(text):
        match = VALUE_PATTERN.match(raw)
        if match:
            name = match.group("name").strip()
            count_text = match.group("count").strip()
        else:
            name = raw.strip()
            count_text = None

        entries.append(
            {
                "name": name,
                "count": _to_numeric_count(count_text),
                "count_text": count_text,
                "raw": raw,
                "normalized_name": _normalize_pathogen_name(name),
            }
        )
    return entries


def _to_numeric_count(text: str | None) -> int | float | str | None:
    if text is None:
        return None
    normalized = text.replace(",", "").strip()
    if not normalized:
        return None
    try:
        if "." in normalized:
            return float(normalized)
        return int(normalized)
    except ValueError:
        return text


def _lookup_column(normalized_map: Dict[str, str], candidates: Iterable[str]) -> str:
    for candidate in candidates:
        key = _normalize_label(candidate)
        if key in normalized_map:
            return normalized_map[key]
    raise KeyError(f"Unable to locate any of: {', '.join(candidates)}")


def _resolve_b_excel_columns(columns: Iterable[str]) -> Dict[str, str]:
    normalized_map = {_normalize_label(col): col for col in columns}
    resolved: Dict[str, str] = {}
    resolved["specimen_code"] = _lookup_column(normalized_map, SPECIMEN_CODE_HEADERS)
    resolved["collected_time"] = _lookup_column(normalized_map, COLLECTED_TIME_HEADERS)
    for key, label in SUSPECTED_COLUMN_LABELS.items():
        resolved[key] = _lookup_column(normalized_map, [label])
    return resolved


def _resolve_seq_excel_columns(columns: Iterable[str]) -> Dict[str, str]:
    normalized_map = {_normalize_label(col): col for col in columns}
    resolved: Dict[str, str] = {}
    resolved["specimen_code"] = _lookup_column(normalized_map, SPECIMEN_CODE_HEADERS)
    try:
        resolved["seq_id"] = _lookup_column(normalized_map, ["系統代碼", "系統碼", "system code", "system_code"])
    except KeyError:
        column_list = list(columns)
        if len(column_list) <= 10:
            raise
        resolved["seq_id"] = column_list[10]
    return resolved


def load_authoritative_seq_ids(seq_excel_path: Path) -> Dict[str, List[str]]:
    df = _read_excel_all_sheets(seq_excel_path)
    column_lookup = _resolve_seq_excel_columns(df.columns)
    specimen_col = column_lookup["specimen_code"]
    seq_col = column_lookup["seq_id"]

    mapping: Dict[str, List[str]] = {}
    for _, row in df[[specimen_col, seq_col]].iterrows():
        specimen_code = _normalize_specimen_code(row[specimen_col])
        seq_id = _coerce_value(row[seq_col])
        if not specimen_code or not seq_id:
            continue
        entries = mapping.setdefault(specimen_code, [])
        if seq_id not in entries:
            entries.append(seq_id)
    return mapping


def _build_specimen_patient_index(mapping: Dict[str, Dict[str, Any]]) -> Dict[str, List[str]]:
    index: Dict[str, List[str]] = {}
    for patient_id, payload in mapping.items():
        codes = payload.get("specimen_codes") or []
        for code in codes:
            normalized = _normalize_specimen_code(code)
            if not normalized:
                continue
            index.setdefault(normalized, []).append(str(patient_id))
    return index


def _load_culture_lookup(path: Path | None) -> Dict[str, set[str]]:
    if path is None:
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        logging.warning("Culture summary not found: %s", path)
        return {}
    except json.JSONDecodeError as exc:
        logging.warning("Culture summary is not valid JSON (%s): %s", path, exc)
        return {}

    lookup: Dict[str, set[str]] = {}
    for patient_id, entries in data.items():
        if not isinstance(entries, list):
            continue
        normalized_set: set[str] = set()
        for entry in entries:
            organism: Any = ""
            if isinstance(entry, dict):
                organism = entry.get("organism") or entry.get("organism_name")
            else:
                organism = entry
            normalized = _normalize_pathogen_name(organism)
            if normalized:
                normalized_set.add(normalized)
        if normalized_set:
            lookup[str(patient_id)] = normalized_set
    return lookup


def attach_suspected_pathogens(
    mapping: Dict[str, Dict[str, Any]],
    b_excel_path: Path,
    culture_lookup: Dict[str, set[str]] | None = None,
) -> None:
    if not b_excel_path.exists():
        raise FileNotFoundError(f"B Excel file not found: {b_excel_path}")

    df = _read_excel_all_sheets(b_excel_path)
    column_lookup = _resolve_b_excel_columns(df.columns)
    specimen_col = column_lookup["specimen_code"]
    collected_col = column_lookup["collected_time"]

    specimen_index = _build_specimen_patient_index(mapping)
    if not specimen_index:
        return

    for _, row in df.iterrows():
        normalized_code = _normalize_specimen_code(row.get(specimen_col))
        if not normalized_code:
            continue
        patient_ids = specimen_index.get(normalized_code)
        if not patient_ids:
            continue

        collected_time = _format_datetime(row.get(collected_col))
        code_display = _coerce_cell(row.get(specimen_col)) or normalized_code

        for category_key in ("bacterial", "fungal", "viral"):
            column_name = column_lookup[category_key]
            parsed_entries = _parse_suspected_field(row.get(column_name))
            if not parsed_entries:
                continue

            type_label = CATEGORY_DISPLAY[category_key]
            column_label = SUSPECTED_COLUMN_LABELS[category_key].replace("\n", "")

            for patient_id in patient_ids:
                patient_entry = mapping.setdefault(patient_id, {"specimen_codes": [], "records": []})
                target_list: List[Dict[str, Any]] = patient_entry.setdefault("suspected_pathogens", [])
                culture_set = culture_lookup.get(str(patient_id)) if culture_lookup else None

                for parsed in parsed_entries:
                    normalized_name = parsed["normalized_name"]
                    matched = bool(culture_set) and bool(normalized_name) and normalized_name in culture_set

                    payload = {
                        "specimen_code": code_display,
                        "collected_time": collected_time,
                        "category": category_key,
                        "pathogen_type": type_label,
                        "column_label": column_label,
                        "name": parsed["name"],
                        "count": parsed["count"],
                        "count_text": parsed["count_text"],
                        "raw_value": parsed["raw"],
                        "matched_culture": matched,
                    }
                    target_list.append(payload)
                    if matched:
                        patient_entry.setdefault("culture_overlaps", []).append(payload)


def load_patient_specimens(excel_path: Path) -> List[PatientSpecimenRow]:
    df = _read_excel_all_sheets(excel_path)
    patient_col = _detect_column(df.columns, PATIENT_NO_HEADERS, "[patient No.]")  # type: ignore
    specimen_col = _detect_column(df.columns, SPECIMEN_CODE_HEADERS, "[檢體代碼]")  # type: ignore
    system_code_col = _detect_column(df.columns, CSV_SEQ_HEADERS, "[seq_id/system code]")  # type: ignore

    records: List[PatientSpecimenRow] = []
    for _, row in df[[patient_col, specimen_col, system_code_col]].iterrows():
        patient_no = _coerce_value(row[patient_col])
        specimen_codes = _split_specimen_codes(row[specimen_col])
        system_code = _coerce_value(row[system_code_col]) or None
        if not patient_no or not specimen_codes:
            continue
        for specimen_code in specimen_codes:
            records.append(
                PatientSpecimenRow(
                    patient_no=patient_no,
                    specimen_code=specimen_code,
                    system_code=system_code,
                )
            )
    return records


def load_bridge_records(bridge_path: Path) -> Dict[str, List[BridgeRow]]:
    df = _load_table(bridge_path)
    specimen_code_col = _detect_column(df.columns, BRIDGE_SPECIMEN_HEADERS, "[檢體代碼/specimen_id]")  # type: ignore
    specimen_id_col = _detect_column(df.columns, SPECIMEN_ID_HEADERS, "[specimen_id]")  # type: ignore
    case_id_col = _detect_column(df.columns, CASE_REVIEW_ID_HEADERS, "[specimen_id_CaseReview]")  # type: ignore

    mapping: Dict[str, List[BridgeRow]] = {}
    for _, row in df[[specimen_code_col, specimen_id_col, case_id_col]].iterrows():
        specimen_code = _coerce_value(row[specimen_code_col])
        specimen_id = _coerce_value(row[specimen_id_col])
        case_id = _coerce_value(row[case_id_col])
        if not specimen_code or not specimen_id or not case_id:
            continue
        entries = mapping.setdefault(specimen_code, [])
        bridge_row = BridgeRow(
            specimen_code=specimen_code,
            specimen_id=specimen_id,
            specimen_id_case_review=case_id,
        )
        if bridge_row not in entries:
            entries.append(bridge_row)
    return mapping


def load_case_review_records(case_review_path: Path) -> Dict[str, List[CaseReviewRow]]:
    df = _load_table(case_review_path)
    case_id_col = _detect_column(df.columns, CASE_REVIEW_ID_HEADERS, "[specimen_id_CaseReview]")  # type: ignore
    seq_col = _detect_column(df.columns, CSV_SEQ_HEADERS, "[seq_id]")  # type: ignore
    chip_col = _detect_column(df.columns, CSV_CHIP_HEADERS, "[chip_id]")  # type: ignore
    specimen_id_col: str | None = None
    try:
        candidate = _detect_column(df.columns, SPECIMEN_ID_HEADERS, "[specimen_id]")  # type: ignore
        if candidate != case_id_col:
            specimen_id_col = candidate
    except KeyError:
        specimen_id_col = None

    subset = {case_id_col, seq_col, chip_col}
    if specimen_id_col:
        subset.add(specimen_id_col)

    mapping: Dict[str, List[CaseReviewRow]] = {}
    for _, row in df[list(subset)].iterrows():
        case_id = _coerce_value(row[case_id_col])
        if not case_id:
            continue
        seq_id = _coerce_value(row[seq_col]) or None
        chip_id = _coerce_value(row[chip_col]) or None
        specimen_id = _coerce_value(row[specimen_id_col]) if specimen_id_col else None
        mapping.setdefault(case_id, []).append(
            CaseReviewRow(
                specimen_id_case_review=case_id,
                specimen_id=specimen_id or None,
                seq_id=seq_id,
                chip_id=chip_id,
            )
        )
    return mapping


def build_patient_mapping(
    excel_path: Path,
    bridge_path: Path,
    case_review_path: Path,
    seq_excel_path: Path | None = None,
) -> Dict[str, Dict[str, Any]]:
    excel_records = load_patient_specimens(excel_path)
    bridge_records = load_bridge_records(bridge_path)
    case_review_records = load_case_review_records(case_review_path)
    authoritative_seq_ids = load_authoritative_seq_ids(seq_excel_path) if seq_excel_path else {}

    result: Dict[str, Dict[str, Any]] = {}
    for record in excel_records:
        patient_entry = result.setdefault(record.patient_no, {"specimen_codes": [], "records": []})

        if record.specimen_code not in patient_entry["specimen_codes"]:
            patient_entry["specimen_codes"].append(record.specimen_code)  # type: ignore

        normalized_specimen_code = _normalize_specimen_code(record.specimen_code)
        seq_ids_for_specimen = authoritative_seq_ids.get(normalized_specimen_code, [])
        preferred_seq_ids = [record.system_code] if record.system_code else seq_ids_for_specimen
        bridges = bridge_records.get(record.specimen_code)
        if not bridges:
            fallback_seq_ids = preferred_seq_ids or [None]
            for seq_id in fallback_seq_ids:
                patient_entry["records"].append(
                    {
                        "specimen_code": record.specimen_code,
                        "chip_id": record.system_code,
                        "seq_id": seq_id,
                    }
                )
            continue

        matched_any = False
        appended_pairs: set[tuple[str | None, str | None]] = set()
        for bridge_row in bridges:
            linked = case_review_records.get(bridge_row.specimen_id_case_review)
            if not linked:
                continue

            for case_row in linked:
                target_seq_ids = preferred_seq_ids or [case_row.seq_id]
                for seq_id in target_seq_ids:
                    effective_chip_id = record.system_code or case_row.chip_id
                    key = (effective_chip_id, seq_id)
                    if key in appended_pairs:
                        continue
                    patient_entry["records"].append(
                        {
                            "specimen_code": record.specimen_code,
                            "chip_id": effective_chip_id,
                            "seq_id": seq_id,
                        }
                    )
                    appended_pairs.add(key)
                matched_any = True

        if not matched_any:
            fallback_seq_ids = preferred_seq_ids or [None]
            for seq_id in fallback_seq_ids:
                patient_entry["records"].append(
                    {
                        "specimen_code": record.specimen_code,
                        "chip_id": record.system_code,
                        "seq_id": seq_id,
                    }
                )

    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a patient mapping that links "
            "X(2025-7-14-collect-coding-mNGS52.xlsx) -> Y(bridge 檢體代碼/specimen_id/specimen_id_CaseReview) "
            "-> Z(case review sequencer data)."
        )
    )
    parser.add_argument(
        "excel",
        type=Path,
        help="Path to the Excel file (X) containing patient_no and 檢體代碼 columns.",
    )
    parser.add_argument(
        "bridge",
        type=Path,
        help="Path to file Y that maps 檢體代碼 ↔ specimen_id ↔ specimen_id_CaseReview.",
    )
    parser.add_argument(
        "case_review",
        type=Path,
        help="Path to file Z that maps specimen_id_CaseReview to sequencing metadata.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional JSON output path. When omitted the mapping is only built in-memory.",
    )
    parser.add_argument(
        "--b-excel",
        type=Path,
        help="Optional path to the B Excel file that contains the 疑似致病菌 columns.",
    )
    parser.add_argument(
        "--seq-excel",
        type=Path,
        help=(
            "Optional path to the Excel file whose system-code column should be used as the "
            "authoritative seq_id source. When omitted, --b-excel is reused if available."
        ),
    )
    parser.add_argument(
        "--culture-summary",
        type=Path,
        help="Optional JSON file containing culture organisms (e.g., outputs/balf_culture_summary.json) "
        "to highlight overlaps.",
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="Indentation level for JSON output (default: 2).",
    )
    parser.add_argument(
        "--ensure-ascii",
        action="store_true",
        help="Force ASCII-only JSON output. Defaults to UTF-8 to keep Chinese headers readable.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    excel_path = args.excel.expanduser().resolve()
    bridge_path = args.bridge.expanduser().resolve()
    case_review_path = args.case_review.expanduser().resolve()
    b_excel_path = args.b_excel.expanduser().resolve() if args.b_excel else None
    seq_excel_path = args.seq_excel.expanduser().resolve() if args.seq_excel else b_excel_path
    culture_summary_path = args.culture_summary.expanduser().resolve() if args.culture_summary else None

    if not excel_path.exists():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    if not bridge_path.exists():
        raise FileNotFoundError(f"Bridge file not found: {bridge_path}")
    if not case_review_path.exists():
        raise FileNotFoundError(f"Case review file not found: {case_review_path}")
    if b_excel_path and not b_excel_path.exists():
        raise FileNotFoundError(f"B Excel file not found: {b_excel_path}")
    if seq_excel_path and not seq_excel_path.exists():
        raise FileNotFoundError(f"Seq Excel file not found: {seq_excel_path}")
    if culture_summary_path and not culture_summary_path.exists():
        logging.warning("Culture summary file does not exist yet: %s", culture_summary_path)

    mapping = build_patient_mapping(excel_path, bridge_path, case_review_path, seq_excel_path)
    if b_excel_path:
        culture_lookup = _load_culture_lookup(culture_summary_path)
        attach_suspected_pathogens(mapping, b_excel_path, culture_lookup)

    if args.output:
        output_path = args.output.expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", encoding="utf-8") as fh:
            json.dump(mapping, fh, indent=args.indent, ensure_ascii=args.ensure_ascii)
        print(f"Mapping written to {output_path}")
    else:
        print("Mapping built in memory. Use --output to persist as JSON.")


if __name__ == "__main__":
    main()

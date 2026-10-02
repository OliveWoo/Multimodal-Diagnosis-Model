from __future__ import annotations

import argparse
import copy
import json
from datetime import datetime, timedelta
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Union

import pandas as pd

from tools import pathogen_normalization as pathogen_names

BAL_LABELS = {"balf", "lower bal"}
NAME_SANITIZER = re.compile(r"[^a-z0-9]+")
VALUE_PATTERN = re.compile(r"^(?P<name>.+?)\s*\((?P<count>[^()]+)\)\s*$")
TIME_WINDOW = timedelta(days=1)
DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent / "culture_bal_config.json"
DEFAULT_PATIENT_MAPPING_PATH = Path("outputs") / "patient_mapping.json"
DEFAULT_BAL_OUTPUT_PATH: Path | None = None
DEFAULT_LOOKUP_OUTPUT_PATH = Path("outputs") / "specimen_lookup_summary.json"
DEFAULT_GT_OUTPUT_PATH = Path("outputs") / "possible_gt_matches.json"


def collect_bal_organisms(root_dir: Union[str, Path]) -> Dict[str, List[Dict[str, str]]]:
    """Scan culture JSON files and collect organisms for BALF-derived samples."""

    root = Path(root_dir).expanduser().resolve()
    results: Dict[str, List[Dict[str, str]]] = {}

    for patient_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        if not (patient_dir.name.startswith("NGS_patient_") and patient_dir.name.endswith("_json")):
            continue

        patient_id = patient_dir.name[len("NGS_patient_") : -len("_json")]
        culture_path = patient_dir / f"NGS_patient_{patient_id}_culture.json"
        if not culture_path.exists():
            continue

        try:
            payload = json.loads(culture_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue

        target_entries = _extract_bal_records(payload)
        if not target_entries:
            continue

        entries = results.setdefault(patient_id, [])
        entries.extend(
            [
                {
                    "sample": str(entry.get("sample", "")),
                    "organism": str(entry.get("organism", "")),
                    "collected_time": str(entry.get("collected_time", "")),
                    "reported_time": str(entry.get("reported_time", "")),
                    "source_file": str(culture_path),
                }
                for entry in target_entries
                if entry.get("organism")
            ]
        )

    return results


def _extract_bal_records(entries: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    matches: List[Dict[str, Any]] = []
    for entry in entries:
        sample = (entry.get("sample") or "").strip().lower()
        if sample in BAL_LABELS:
            matches.append(entry)
    return matches


def determine_possible_gt(
    specimen_lookup_source: Union[str, Path, Dict[str, Any]],
    balf_summary_source: Union[str, Path, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Determine possible ground truths by comparing lookup specimens with BALF cultures."""

    lookup_data = _load_json_payload(specimen_lookup_source)
    balf_data = _load_json_payload(balf_summary_source)

    normalized_culture: Dict[str, List[Dict[str, Any]]] = {}
    for patient_id, entries in balf_data.items():
        normalized_entries: List[Dict[str, Any]] = []
        for entry in entries:
            sample = str(entry.get("sample") or entry.get("specimen_type") or "").strip()
            if not sample:
                continue
            normalized_entries.append(
                {
                    "collected_time": _parse_datetime(entry.get("collected_time")),
                    "organism": str(entry.get("organism") or entry.get("organism_name") or ""),
                    "sample": sample,
                }
            )
        if normalized_entries:
            normalized_culture[str(patient_id)] = normalized_entries

    results: List[Dict[str, Any]] = []
    for patient_id, specimens in lookup_data.items():
        patient_key = str(patient_id)
        cultures = normalized_culture.get(patient_key)
        if not cultures:
            continue

        for specimen in specimens:
            specimen_time = _parse_datetime(specimen.get("collected_time"))
            suspected_lookup = _gather_suspected_pathogens(specimen)

            matched: List[str] = []
            if specimen_time and suspected_lookup:
                for culture in cultures:
                    if not _is_bal_sample(culture.get("sample")):
                        continue
                    culture_time = culture.get("collected_time")
                    if not culture_time:
                        continue
                    if abs(culture_time - specimen_time) > TIME_WINDOW:
                        continue
                    organism = culture.get("organism")
                    normalized_org = _normalize_name_for_compare(organism)
                    if not normalized_org:
                        continue
                    match_info = suspected_lookup.get(normalized_org)
                    if match_info:
                        label = f"{match_info['category']}: {match_info['name']}"
                        if label not in matched:
                            matched.append(label)

            results.append(
                {
                    "patient_id": patient_key,
                    "matched_pathogens": matched,
                }
            )

    return results


def _normalize_label(label: Any) -> str:
    return "".join(str(label).split()).lower()


def _resolve_columns(columns: Iterable[str]) -> Dict[str, str]:
    normalized_map = {_normalize_label(col): col for col in columns}
    required = [
        SPECIMEN_CODE_COLUMN,
        COLLECTED_TIME_COLUMN,
        BACTERIA_COLUMN,
        VIRUS_COLUMN,
        FUNGUS_COLUMN,
    ]
    resolved: Dict[str, str] = {}
    for target in required:
        key = _normalize_label(target)
        if key not in normalized_map:
            raise KeyError(f"Missing column '{target}' in B Excel file.")
        resolved[target] = normalized_map[key]
    return resolved


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


def _split_lines(text: str) -> List[str]:
    parts: List[str] = []
    normalized = text.replace("\r", "\n")
    for raw_line in normalized.split("\n"):
        for chunk in raw_line.split(";"):
            candidate = chunk.strip()
            if candidate and candidate != "-":
                parts.append(candidate)
    return parts


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


def _parse_suspected_entries(
    value: Any,
    *,
    category_key: str,
    category_label: str,
) -> List[Dict[str, Any]]:
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
                "category": category_label,
                "_normalized_name": _normalize_name_for_compare(name),
                "matched_culture": False,
            }
        )
    return entries


def _read_excel_all_sheets(path: Path) -> pd.DataFrame:
    frames = pd.read_excel(path, sheet_name=None, engine="openpyxl")
    # Combine every non-empty sheet so we don't miss duplicate specimens in other tabs.
    if isinstance(frames, dict):
        non_empty = [frame for frame in frames.values() if frame is not None and not frame.empty]
        if not non_empty:
            return pd.DataFrame()
        return pd.concat(non_empty, ignore_index=True)
    return frames


def _build_culture_lookup_from_bal(
    bal_results: Dict[str, List[Dict[str, Any]]] | None,
) -> Dict[str, set[str]]:
    if not bal_results:
        return {}
    lookup: Dict[str, set[str]] = {}
    for patient_id, entries in bal_results.items():
        normalized: set[str] = set()
        for entry in entries:
            organism = ""
            if isinstance(entry, dict):
                organism = entry.get("organism") or entry.get("organism_name") or ""
            else:
                organism = str(entry or "")
            norm = _normalize_name_for_compare(organism)
            if norm:
                normalized.add(norm)
        if normalized:
            lookup[str(patient_id)] = normalized
    return lookup


def _load_config(path: Path | None) -> Dict[str, Any]:
    if path is None:
        return {}
    path = path.expanduser().resolve()
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # pylint: disable=broad-exception-caught
        raise ValueError(f"Failed to read config file {path}: {exc}") from exc
    return data if isinstance(data, dict) else {}


def _pick_path(
    cli_value: Path | None,
    config_value: Any,
    default_path: Path | None = None,
) -> Path | None:
    if cli_value is not None:
        return cli_value.expanduser().resolve()
    if config_value:
        return Path(str(config_value)).expanduser().resolve()
    if default_path is not None:
        return default_path.expanduser().resolve()
    return None


def _write_json(path: Path, payload: Any, indent: int, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=indent), encoding="utf-8")
    print(f"Wrote {label} to {path}")


def _load_json_payload(source: Union[str, Path, Dict[str, Any]]) -> Dict[str, Any]:
    if isinstance(source, dict):
        return source
    path = Path(source).expanduser().resolve()
    return json.loads(path.read_text(encoding="utf-8"))


def _parse_datetime(value: Any) -> datetime | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    try:
        parsed = pd.to_datetime(value, errors="coerce")
    except Exception:
        return None
    if pd.isna(parsed):
        return None
    if isinstance(parsed, pd.Timestamp):
        return parsed.to_pydatetime()
    if isinstance(parsed, datetime):
        return parsed
    return None


def _is_bal_sample(sample: Any) -> bool:
    label = str(sample or "").strip().lower()
    return label in BAL_LABELS


def _split_pathogen_field(value: Any) -> List[str]:
    text = _coerce_cell(value)
    if not text or text == "-":
        return []
    parts = []
    for raw in text.replace("\r", "\n").split("\n"):
        candidate = _strip_count_suffix(raw)
        if candidate:
            parts.append(candidate)
    return parts


def _strip_count_suffix(value: str) -> str:
    text = value.strip()
    if not text:
        return ""
    if "(" in text:
        return text.split("(", 1)[0].strip()
    return text


def _gather_suspected_pathogens(specimen: Dict[str, Any]) -> Dict[str, Dict[str, str]]:
    """Return normalized pathogen names mapped to their original label and category."""

    mapping: Dict[str, Dict[str, str]] = {}

    detailed_entries = specimen.get("suspected_pathogens")
    if isinstance(detailed_entries, list):
        for entry in detailed_entries:
            name = entry.get("name")
            normalized = entry.get("normalized_name") or _normalize_name_for_compare(name)
            if not normalized:
                continue
            mapping[normalized] = {
                "name": str(name or ""),
                "category": str(
                    entry.get("pathogen_type") or entry.get("column_label") or entry.get("category") or ""
                ),
            }
        if mapping:
            return mapping

    sources = [
        ("bacteria_suspected_pathogens", "Bacteria"),
        ("virus_suspected_pathogens", "Virus"),
        ("fungus_suspected_pathogens", "Fungus"),
    ]

    for field, category in sources:
        for name in _split_pathogen_field(specimen.get(field)):
            normalized = _normalize_name_for_compare(name)
            if normalized and normalized not in mapping:
                mapping[normalized] = {"name": name.strip(), "category": category}

    return mapping


def _normalize_name_for_compare(name: Any) -> str:
    return pathogen_names.canonical_key(name)


__all__ = [
    "collect_bal_organisms",
    "collect_b_specimen_details",
    "determine_possible_gt",
]

SPECIMEN_CODE_COLUMN = "檢體代碼"
COLLECTED_TIME_COLUMN = "採檢時間"
BACTERIA_COLUMN = "細菌-\n疑似致病菌"
VIRUS_COLUMN = "病毒-\n疑似致病菌"
FUNGUS_COLUMN = "真菌-\n疑似致病菌"

SUSPECTED_FIELDS = [
    ("bacteria_suspected_pathogens", BACTERIA_COLUMN, "bacterial", "細菌"),
    ("fungus_suspected_pathogens", FUNGUS_COLUMN, "fungal", "真菌"),
    ("virus_suspected_pathogens", VIRUS_COLUMN, "viral", "病毒"),
]


def collect_b_specimen_details(
    patient_mapping_path: Union[str, Path],
    b_excel_path: Union[str, Path],
    culture_data: Dict[str, List[Dict[str, Any]]] | None = None,
) -> Dict[str, List[Dict[str, Any]]]:
    """Match specimen codes to the B Excel file and extract organism summaries."""

    patient_mapping = json.loads(Path(patient_mapping_path).read_text(encoding="utf-8"))
    df = _read_excel_all_sheets(b_excel_path)

    column_lookup = _resolve_columns(df.columns)
    specimen_col = column_lookup[SPECIMEN_CODE_COLUMN]
    collected_col = column_lookup[COLLECTED_TIME_COLUMN]
    culture_lookup = _build_culture_lookup_from_bal(culture_data)

    specimen_rows: Dict[str, List[Dict[str, Any]]] = {}
    for _, row in df.iterrows():
        code = _normalize_specimen_code(row.get(specimen_col))
        if not code:
            continue

        entry: Dict[str, Any] = {
            "specimen_code": code,
            "collected_time": _format_datetime(row.get(collected_col)),
        }

        suspected: List[Dict[str, Any]] = []
        for field_name, column_label, category_key, category_label in SUSPECTED_FIELDS:
            excel_column = column_lookup[column_label]
            raw_value = _coerce_cell(row.get(excel_column))
            parsed_entries = _parse_suspected_entries(
                raw_value,
                category_key=category_key,
                category_label=category_label,
            )
            if parsed_entries:
                suspected.extend(parsed_entries)

        if suspected:
            entry["suspected_pathogens"] = suspected

        specimen_rows.setdefault(code, []).append(entry)

    results: Dict[str, List[Dict[str, Any]]] = {}
    for patient_id, payload in patient_mapping.items():
        codes = payload.get("specimen_codes") or []
        matched_entries: List[Dict[str, Any]] = []
        culture_set = culture_lookup.get(str(patient_id))
        for code in codes:
            normalized = _normalize_specimen_code(code)
            if not normalized:
                continue
            matches = specimen_rows.get(normalized)
            if matches:
                for base_entry in matches:
                    entry = copy.deepcopy(base_entry)
                    suspected_entries = entry.get("suspected_pathogens") or []
                    overlaps: List[Dict[str, Any]] = []
                    if suspected_entries:
                        for suspect in suspected_entries:
                            normalized_name = suspect.pop("_normalized_name", None)
                            matched = (
                                bool(culture_set)
                                and bool(normalized_name)
                                and normalized_name in culture_set
                            )
                            suspect["matched_culture"] = matched
                            if matched:
                                overlaps.append(suspect)
                    if overlaps:
                        entry["culture_overlaps"] = overlaps
                    matched_entries.append(entry)
        if matched_entries:
            results[str(patient_id)] = matched_entries

    return results


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Collect organisms from BALF/Lower BAL culture reports."
    )
    parser.add_argument(
        "folder",
        type=Path,
        nargs="?",
        help="Root folder containing NGS_patient_[ID]_json subdirectories.",
    )
    parser.add_argument(
        "b_excel",
        type=Path,
        nargs="?",
        help="Path to the Excel file (B) that lists specimen codes and organism calls.",
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="Pretty-print indentation for STDOUT JSON (default: 2).",
    )
    parser.add_argument(
        "--patient-mapping",
        type=Path,
        help="Path to patient_mapping.json (default: outputs/patient_mapping.json).",
    )
    parser.add_argument(
        "--bal-output",
        type=Path,
        help="Output path for the BALF culture summary JSON.",
    )
    parser.add_argument(
        "--lookup-summary-output",
        type=Path,
        help="Output path for specimen/reference Excel matches.",
    )
    parser.add_argument(
        "--gt-output",
        type=Path,
        help="Output path for determine_possible_gt results.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help=(
            "Optional JSON config containing default paths. "
            "If omitted, the tool will use tools/culture_bal_config.json when it exists."
        ),
    )
    return parser


def _main() -> None:
    parser = _build_arg_parser()
    args = parser.parse_args()
    config = _load_config(args.config) if args.config else _load_config(DEFAULT_CONFIG_PATH)

    folder = _pick_path(args.folder, config.get("culture_root"))
    b_excel = _pick_path(args.b_excel, config.get("b_excel"))
    if folder is None or b_excel is None:
        parser.error("Culture folder and B Excel path must be provided via arguments or config.")

    patient_mapping = _pick_path(
        args.patient_mapping,
        config.get("patient_mapping"),
        DEFAULT_PATIENT_MAPPING_PATH,
    )
    bal_output = _pick_path(args.bal_output, config.get("bal_output"), DEFAULT_BAL_OUTPUT_PATH)
    lookup_output = _pick_path(
        args.lookup_summary_output,
        config.get("lookup_summary_output"),
        DEFAULT_LOOKUP_OUTPUT_PATH,
    )
    gt_output = _pick_path(args.gt_output, config.get("gt_output"), DEFAULT_GT_OUTPUT_PATH)
    if patient_mapping is None:
        parser.error("Patient mapping path is required.")
    if lookup_output is None:
        parser.error("Lookup summary output path is required.")

    bal_results = collect_bal_organisms(folder)
    if bal_output:
        _write_json(bal_output, bal_results, args.indent, label="BALF culture summary")
        bal_source: Union[str, Path, Dict[str, Any]] = bal_output
    else:
        bal_source = bal_results

    lookup_summary = collect_b_specimen_details(patient_mapping, b_excel, bal_results)
    _write_json(lookup_output, lookup_summary, args.indent, label="Reference Excel specimen summary")

    gt_results = determine_possible_gt(lookup_output, bal_source)
    _write_json(gt_output, gt_results, args.indent, label="possible GT matches")


if __name__ == "__main__":
    _main()

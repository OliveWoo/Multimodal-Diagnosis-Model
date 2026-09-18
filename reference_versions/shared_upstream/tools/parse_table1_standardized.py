from __future__ import annotations

import argparse
import copy
import json
import re
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.parse_table1_workbook import DEFAULT_SHEET_NAME  # noqa: E402
from tools.table1_sections import load_table1_grid, split_titled_sections  # noqa: E402
from tools.table1_standardizer_engine import build_standardized_payloads  # noqa: E402
from utils import sanitize_filename  # noqa: E402


STANDARD_SUFFIXES = (
    "underlying",
    "admission_diagnosis",
    "culture",
    "filmarray",
    "gm_test",
    "molecular_microbiology",
    "image",
    "CBC",
    "other_lab",
    "mNGS_grouped",
)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _is_blank(value: Any) -> bool:
    return value is None or value == "" or value == [] or value == {}


def _merge_dict(baseline: dict[str, Any], extracted: dict[str, Any]) -> dict[str, Any]:
    merged = copy.deepcopy(baseline)
    for key, value in extracted.items():
        previous = merged.get(key)
        if isinstance(previous, dict) and isinstance(value, dict):
            merged[key] = _merge_dict(previous, value)
        elif isinstance(previous, list) and isinstance(value, list):
            combined = copy.deepcopy(previous)
            fingerprints = {
                json.dumps(item, ensure_ascii=False, sort_keys=True) for item in combined
            }
            for item in value:
                fingerprint = json.dumps(item, ensure_ascii=False, sort_keys=True)
                if fingerprint not in fingerprints:
                    combined.append(copy.deepcopy(item))
                    fingerprints.add(fingerprint)
            merged[key] = combined
        elif not _is_blank(value):
            merged[key] = copy.deepcopy(value)
        elif key not in merged:
            merged[key] = copy.deepcopy(value)
    return merged


RECORD_IDENTITY_KEYS: dict[str, tuple[str, ...]] = {
    "admission_diagnosis": ("test", "diagnosis"),
    "CBC": ("reported_time", "item", "item_full", "value", "test"),
    "other_lab": ("reported_time", "item", "item_full", "value", "test"),
    "filmarray": ("test", "sample", "collected_time", "reported_time", "target", "result"),
    "gm_test": ("test", "sample", "collected_time", "reported_time", "target", "result", "value"),
    "molecular_microbiology": (
        "test",
        "sample",
        "collected_time",
        "reported_time",
        "target",
        "result",
        "value",
    ),
    "culture": (
        "test",
        "sample",
        "collected_time",
        "received_time",
        "reported_time",
        "organism",
        "status",
    ),
    "image": ("test", "exam_type", "collected_time", "reported_time", "findings"),
    "mNGS_grouped": ("specimen_code", "collected_time", "specimen_site"),
}


def _record_identity(suffix: str, record: Any) -> str:
    if not isinstance(record, dict):
        return json.dumps(record, ensure_ascii=False, sort_keys=True)
    keys = RECORD_IDENTITY_KEYS.get(suffix)
    if not keys:
        return json.dumps(record, ensure_ascii=False, sort_keys=True)
    values = [record.get(key) for key in keys]
    return json.dumps(values, ensure_ascii=False, sort_keys=True)


def _merge_record_lists(suffix: str, baseline: list[Any], extracted: list[Any]) -> list[Any]:
    if suffix == "underlying" and baseline and extracted:
        first_old = baseline[0] if isinstance(baseline[0], dict) else {}
        first_new = extracted[0] if isinstance(extracted[0], dict) else {}
        return [_merge_dict(first_old, first_new), *copy.deepcopy(baseline[1:]), *copy.deepcopy(extracted[1:])]

    if not baseline:
        return copy.deepcopy(extracted)
    if not extracted:
        return copy.deepcopy(baseline)

    merged = copy.deepcopy(baseline)
    baseline_positions = {
        _record_identity(suffix, record): index for index, record in enumerate(merged)
    }
    for record in extracted:
        identity = _record_identity(suffix, record)
        if identity in baseline_positions:
            index = baseline_positions[identity]
            if isinstance(merged[index], dict) and isinstance(record, dict):
                merged[index] = _merge_dict(merged[index], record)
            continue
        merged.append(copy.deepcopy(record))
    return merged


def _resolve_baseline_dir(path: Path, base_name: str) -> Path:
    resolved = path.expanduser().resolve()
    patient_dir = resolved / f"{base_name}_json"
    return patient_dir if patient_dir.is_dir() else resolved


def _merge_baseline_payloads(
    payloads: dict[str, Any], baseline_dir: Path, base_name: str
) -> dict[str, Any]:
    merged = copy.deepcopy(payloads)
    for suffix in STANDARD_SUFFIXES:
        baseline_path = baseline_dir / f"{base_name}_{suffix}.json"
        if not baseline_path.exists():
            continue
        baseline = _read_json(baseline_path)
        extracted = merged.get(suffix, [])
        if isinstance(baseline, list) and isinstance(extracted, list):
            merged[suffix] = _merge_record_lists(suffix, baseline, extracted)
        elif _is_blank(extracted):
            merged[suffix] = copy.deepcopy(baseline)
    return merged


def _build_source_mapping(base_name: str) -> dict[str, list[str]]:
    return {
        "CBC/otherLab": [f"{base_name}_CBC.json", f"{base_name}_other_lab.json"],
        "filmarray": [
            f"{base_name}_filmarray.json",
            f"{base_name}_gm_test.json",
            f"{base_name}_molecular_microbiology.json",
        ],
        "culture": [f"{base_name}_culture.json"],
        "image": [f"{base_name}_image.json"],
        "admission_diagnosis": [f"{base_name}_admission_diagnosis.json"],
        "underlying": [f"{base_name}_underlying.json"],
    }


def _patient_base_name(excel_path: Path, patient_id: str | None) -> str:
    if patient_id:
        normalized = patient_id.strip()
        return normalized if normalized.startswith("NGS_patient_") else f"NGS_patient_{normalized}"
    match = re.search(r"patient\s*[_ -]*(\d+)", excel_path.stem, re.IGNORECASE)
    if match:
        return f"NGS_patient_{int(match.group(1))}"
    return f"NGS_patient_{sanitize_filename(excel_path.stem)}"


def export_standardized(
    excel_path: Path,
    output_root: Path | None = None,
    *,
    patient_id: str | None = None,
    sheet_name: str = DEFAULT_SHEET_NAME,
    baseline_dir: Path | None = None,
) -> Path:
    excel_path = excel_path.expanduser().resolve()
    output_root = (output_root or (excel_path.parent / "table1_patient_info_standardized")).expanduser().resolve()
    base_name = _patient_base_name(excel_path, patient_id)
    patient_dir = output_root / f"{base_name}_json"
    patient_dir.mkdir(parents=True, exist_ok=True)
    (patient_dir / "agent_outputs").mkdir(exist_ok=True)
    (patient_dir / "summary_outputs").mkdir(exist_ok=True)

    grid = load_table1_grid(excel_path, sheet_name)
    image_start = next(
        (
            row.number
            for row in grid.iter_rows()
            if "Radiology / Special Examination Results" in row.text
            or "放射報告與特殊檢查記錄" in row.text
        ),
        None,
    )
    sections = split_titled_sections(grid, end_row=(image_start - 1 if image_start else None))
    payloads = build_standardized_payloads(grid, sections)
    if baseline_dir is not None:
        payloads = _merge_baseline_payloads(
            payloads,
            _resolve_baseline_dir(baseline_dir, base_name),
            base_name,
        )

    for suffix in STANDARD_SUFFIXES:
        _write_json(patient_dir / f"{base_name}_{suffix}.json", payloads.get(suffix, []))
    _write_json(patient_dir / "source_mapping.json", _build_source_mapping(base_name))
    return patient_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export a one-sheet Table 1 workbook to patient_info_standardized-style JSON files.",
    )
    parser.add_argument("target", type=Path, help="A Table 1 .xlsx/.xlsm file.")
    parser.add_argument("--output-root", type=Path, help="Output root for NGS_patient_x_json folders.")
    parser.add_argument("--patient-id", help="Patient id, e.g. 1 or NGS_patient_1.")
    parser.add_argument("--sheet", default=DEFAULT_SHEET_NAME, help="Source sheet name.")
    parser.add_argument(
        "--baseline-dir",
        type=Path,
        help=(
            "Optional existing NGS_patient_x_json folder (or its parent). "
            "Its mNGS/demographic data are preserved without modifying the baseline files."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    patient_dir = export_standardized(
        args.target,
        args.output_root,
        patient_id=args.patient_id,
        sheet_name=args.sheet,
        baseline_dir=args.baseline_dir,
    )
    print(f"[Table1 standardized] {args.target.name} -> {patient_dir}")
    for path in sorted(patient_dir.glob("*.json")):
        if path.name == "source_mapping.json":
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        print(f"  - {path.name}: {len(payload) if isinstance(payload, list) else 'object'}")


if __name__ == "__main__":
    main()

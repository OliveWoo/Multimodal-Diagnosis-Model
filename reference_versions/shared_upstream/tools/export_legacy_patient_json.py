from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable, Sequence


DEFAULT_PATIENT_EXPORT_ROOT = Path("outputs") / "excel_patient_exports"
DEFAULT_OUTPUT_ROOT = Path("outputs") / "patient_info_standardized"
LEGACY_INDEX_NAME = "legacy_patient_index.json"

CATEGORY_SUFFIXES = (
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


def _read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _find_category_path(patient_dir: Path, suffix: str) -> Path | None:
    matches = sorted(patient_dir.glob(f"*_{suffix}.json"))
    return matches[0] if matches else None


def _find_patient_manifest_path(patient_dir: Path) -> Path | None:
    matches = sorted(patient_dir.glob("*_patient_manifest.json"))
    return matches[0] if matches else None


def _default_payload_for_suffix(suffix: str) -> Any:
    if suffix in CATEGORY_SUFFIXES:
        return []
    return {}


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


def _load_patient_index(patient_export_root: Path) -> list[dict[str, Any]]:
    index_path = patient_export_root / "patient_index.json"
    if index_path.exists():
        payload = _read_json(index_path)
        if isinstance(payload, list):
            return [entry for entry in payload if isinstance(entry, dict)]

    entries: list[dict[str, Any]] = []
    for patient_dir in sorted(path for path in patient_export_root.iterdir() if path.is_dir()):
        entries.append(
            {
                "patient_key": patient_dir.name,
                "output_dir": str(patient_dir.resolve()),
            }
        )
    return entries


def _iter_patient_entries(patient_export_root: Path) -> Iterable[tuple[int, dict[str, Any], Path]]:
    for index, entry in enumerate(_load_patient_index(patient_export_root), start=1):
        patient_key = str(entry.get("patient_key") or "").strip()
        output_dir = str(entry.get("output_dir") or "").strip()
        patient_dir = Path(output_dir) if output_dir else patient_export_root / patient_key
        if not patient_dir.exists() or not patient_dir.is_dir():
            continue
        yield index, entry, patient_dir


def export_legacy_patient_json(patient_export_root: Path, output_root: Path) -> tuple[list[Path], Path]:
    output_root.mkdir(parents=True, exist_ok=True)
    written_dirs: list[Path] = []
    legacy_index: list[dict[str, Any]] = []

    for ordinal, entry, patient_dir in _iter_patient_entries(patient_export_root):
        legacy_name = f"NGS_patient_{ordinal}"
        legacy_dir = output_root / f"{legacy_name}_json"
        legacy_dir.mkdir(parents=True, exist_ok=True)
        (legacy_dir / "agent_outputs").mkdir(exist_ok=True)
        (legacy_dir / "summary_outputs").mkdir(exist_ok=True)

        for suffix in CATEGORY_SUFFIXES:
            source_path = _find_category_path(patient_dir, suffix)
            payload = _read_json(source_path) if source_path else _default_payload_for_suffix(suffix)
            target_path = legacy_dir / f"{legacy_name}_{suffix}.json"
            _write_json(target_path, payload)

        _write_json(legacy_dir / "source_mapping.json", _build_source_mapping(legacy_name))
        manifest_path = _find_patient_manifest_path(patient_dir)
        manifest = _read_json(manifest_path) if manifest_path else {}
        mngs_collected_dates = manifest.get("mngs_collected_dates", []) if isinstance(manifest, dict) else []

        legacy_index.append(
            {
                "legacy_patient_id": legacy_name,
                "legacy_dir": str(legacy_dir.resolve()),
                "patient_key": entry.get("patient_key", patient_dir.name),
                "case_code": entry.get("case_code", ""),
                "specimen_id": entry.get("specimen_id", ""),
                "hospital_id": entry.get("hospital_id", ""),
                "mngs_collected_dates": mngs_collected_dates if isinstance(mngs_collected_dates, list) else [],
                "source_dir": str(patient_dir.resolve()),
            }
        )
        written_dirs.append(legacy_dir)

    index_path = output_root / LEGACY_INDEX_NAME
    _write_json(index_path, legacy_index)
    return written_dirs, index_path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export per-patient JSON folders to the legacy NGS_patient_x_json layout.",
    )
    parser.add_argument(
        "--patient-export-root",
        type=Path,
        default=DEFAULT_PATIENT_EXPORT_ROOT,
        help=f"Path to the per-patient export root (default: {DEFAULT_PATIENT_EXPORT_ROOT}).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help=f"Where to write the legacy-format patient folders (default: {DEFAULT_OUTPUT_ROOT}).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    patient_export_root = args.patient_export_root.expanduser().resolve()
    output_root = args.output.expanduser().resolve()
    written_dirs, index_path = export_legacy_patient_json(patient_export_root, output_root)
    print(f"Exported {len(written_dirs)} legacy patient folders to {output_root}")
    print(f"Wrote legacy index: {index_path}")


if __name__ == "__main__":
    main()

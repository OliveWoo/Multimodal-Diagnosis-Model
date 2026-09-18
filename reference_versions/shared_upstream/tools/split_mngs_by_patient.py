from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping, Sequence

DEFAULT_PATIENT_ROOT = Path(r"C:\Users\User\Desktop\patient_info\data(all)")
OUTPUT_NAME_TEMPLATE = "NGS_patient_{patient_id}_mNGS_grouped.json"


def _load_payload(path: Path) -> Mapping[str, list[dict]]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise TypeError(f"Expected the grouped JSON to be a dict, got {type(data)}")
    return data


def _write_patient_file(
    patient_id: str,
    entries: list[dict],
    target_dir: Path,
) -> Path | None:
    folder = target_dir / f"NGS_patient_{patient_id}_json"
    if not folder.exists():
        return None
    output_path = folder / OUTPUT_NAME_TEMPLATE.format(patient_id=patient_id)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(entries, handle, ensure_ascii=False, indent=2)
    return output_path


def split_grouped_json(
    source_json: Path,
    patient_root: Path,
) -> tuple[list[Path], list[str]]:
    payload = _load_payload(source_json)
    written_files: list[Path] = []
    skipped_patients: list[str] = []
    for patient_id, entries in payload.items():
        if not isinstance(entries, list):
            raise TypeError(f"Expected list for patient {patient_id}, got {type(entries)}")
        output_path = _write_patient_file(patient_id, entries, patient_root)
        if output_path:
            written_files.append(output_path)
        else:
            skipped_patients.append(patient_id)
    return written_files, skipped_patients


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split the grouped mNGS JSON into per-patient files and store them in each patient's folder."
    )
    parser.add_argument(
        "grouped_json",
        type=Path,
        help="Path to the combined specimen_lookup_summary_mNGS_grouped.json file.",
    )
    parser.add_argument(
        "--patient-root",
        type=Path,
        default=DEFAULT_PATIENT_ROOT,
        help=f"Path containing NGS_patient_*_json folders (default: {DEFAULT_PATIENT_ROOT})",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    written, skipped = split_grouped_json(args.grouped_json, args.patient_root)
    print(f"Generated {len(written)} files:")
    for path in written:
        print(f" - {path}")
    if skipped:
        print(f"Skipped {len(skipped)} patients; missing folders: {', '.join(sorted(skipped, key=int))}")


if __name__ == "__main__":
    main()

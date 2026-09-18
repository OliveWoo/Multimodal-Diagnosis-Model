from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

DEFAULT_SOURCE_JSON = Path("outputs/mngs_candidate_microbes.json")
DEFAULT_COMBINED_OUTPUT = Path("outputs/mngs_candidate_microbes_for_llm.json")
DEFAULT_PATIENT_ROOT = Path(r"C:\Users\User\Desktop\patient_info\0310_data")
OUTPUT_NAME_TEMPLATE = "NGS_patient_{patient_id}_mNGS_candidates_for_llm.json"
CATEGORY_MAP = {
    "1.Bac": "bacteria",
    "2.Fungi": "fungi",
    "3.Virus": "virus",
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a simplified per-patient mNGS candidate JSON for LLM use and optionally "
            "write each patient file into existing NGS_patient_*_json folders."
        )
    )
    parser.add_argument(
        "--source-json",
        type=Path,
        default=DEFAULT_SOURCE_JSON,
        help=f"Path to mngs_candidate_microbes.json (default: {DEFAULT_SOURCE_JSON})",
    )
    parser.add_argument(
        "--combined-output",
        type=Path,
        default=DEFAULT_COMBINED_OUTPUT,
        help=f"Path to write the combined simplified JSON (default: {DEFAULT_COMBINED_OUTPUT})",
    )
    parser.add_argument(
        "--patient-root",
        type=Path,
        default=DEFAULT_PATIENT_ROOT,
        help=f"Root containing NGS_patient_*_json folders (default: {DEFAULT_PATIENT_ROOT})",
    )
    return parser.parse_args(argv)


def _load_source_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict payload, got {type(payload)}")
    return payload


def _simplify_patient_records(patient_payload: dict[str, Any]) -> list[dict[str, Any]]:
    records = patient_payload.get("records", [])
    simplified: list[dict[str, Any]] = []
    for record in records:
        for candidate in record.get("candidates", []):
            simplified.append(
                {
                    "name": candidate.get("organism_name"),
                    "reads": candidate.get("sec_hit"),
                    "type": CATEGORY_MAP.get(
                        str(candidate.get("source_category")),
                        candidate.get("source_category"),
                    ),
                }
            )
    return simplified


def build_combined_payload(source_payload: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    patients = source_payload.get("patients", {})
    if not isinstance(patients, dict):
        raise TypeError(f"Expected 'patients' to be a dict, got {type(patients)}")

    combined: dict[str, list[dict[str, Any]]] = {}
    for patient_id, patient_payload in patients.items():
        if not isinstance(patient_payload, dict):
            raise TypeError(f"Expected patient payload for {patient_id} to be a dict")
        combined[str(patient_id)] = _simplify_patient_records(patient_payload)
    return combined


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def write_patient_files(
    combined_payload: dict[str, list[dict[str, Any]]],
    patient_root: Path,
) -> tuple[list[Path], list[str]]:
    written_files: list[Path] = []
    skipped_patients: list[str] = []
    for patient_id, entries in combined_payload.items():
        folder = patient_root / f"NGS_patient_{patient_id}_json"
        if not folder.exists():
            skipped_patients.append(patient_id)
            continue
        output_path = folder / OUTPUT_NAME_TEMPLATE.format(patient_id=patient_id)
        _write_json(output_path, entries)
        written_files.append(output_path)
    return written_files, skipped_patients


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    source_payload = _load_source_json(args.source_json)
    combined_payload = build_combined_payload(source_payload)

    _write_json(args.combined_output, combined_payload)
    print(f"Wrote combined simplified JSON to {args.combined_output}")

    written_files, skipped_patients = write_patient_files(combined_payload, args.patient_root)
    print(f"Wrote {len(written_files)} per-patient files.")
    if skipped_patients:
        print(
            "Skipped patients with no folder: "
            + ", ".join(sorted(skipped_patients, key=lambda value: int(value)))
        )


if __name__ == "__main__":
    main()

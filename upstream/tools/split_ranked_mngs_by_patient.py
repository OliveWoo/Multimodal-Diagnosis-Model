from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

DEFAULT_SOURCE = Path("mngs_candidate_microbes_ranked.json")
DEFAULT_PATIENT_ROOT = Path("outputs") / "patient_info_filtered"
OUTPUT_NAME_TEMPLATE = "NGS_patient_{patient_id}_mNGS_ranked_candidates.json"


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected object JSON: {path}")
    return payload


def split_ranked_mngs(
    source: Path,
    patient_root: Path,
    *,
    overwrite: bool,
) -> dict[str, Any]:
    payload = load_json(source)
    patients = payload.get("patients")
    if not isinstance(patients, dict):
        raise ValueError(f"Missing patients object in {source}")
    ranking_metadata = payload.get("ranking_metadata", {})

    written: list[str] = []
    skipped_missing_folder: list[str] = []
    skipped_existing: list[str] = []

    for patient_id, patient_payload in sorted(patients.items(), key=sort_patient_key):
        if not isinstance(patient_payload, dict):
            continue
        folder = patient_root / f"NGS_patient_{patient_id}_json"
        if not folder.exists():
            skipped_missing_folder.append(str(patient_id))
            continue
        output_path = folder / OUTPUT_NAME_TEMPLATE.format(patient_id=patient_id)
        if output_path.exists() and not overwrite:
            skipped_existing.append(str(output_path))
            continue
        output_payload = {
            "patient_id": str(patient_id),
            "records": patient_payload.get("records", []),
            "ranking_metadata": ranking_metadata,
        }
        output_path.write_text(
            json.dumps(output_payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        written.append(str(output_path))

    return {
        "source": str(source),
        "patient_root": str(patient_root),
        "written_count": len(written),
        "skipped_missing_folder_count": len(skipped_missing_folder),
        "skipped_existing_count": len(skipped_existing),
        "written": written,
        "skipped_missing_folder": skipped_missing_folder,
        "skipped_existing": skipped_existing,
    }


def sort_patient_key(item: tuple[str, Any]) -> tuple[int, str]:
    value = str(item[0])
    try:
        return (0, f"{int(value):08d}")
    except ValueError:
        return (1, value)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split mngs_candidate_microbes_ranked.json into per-patient folders."
    )
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Skip files that already exist instead of overwriting them.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("outputs") / "split_ranked_mngs_by_patient_report.json",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report = split_ranked_mngs(
        args.source,
        args.patient_root,
        overwrite=not args.no_overwrite,
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"written={report['written_count']}")
    print(f"skipped_missing_folder={report['skipped_missing_folder_count']}")
    print(f"skipped_existing={report['skipped_existing_count']}")
    print(f"report={args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

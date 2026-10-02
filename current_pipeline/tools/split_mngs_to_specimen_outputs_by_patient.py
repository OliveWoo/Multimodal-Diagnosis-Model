from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence


DEFAULT_INPUT = Path("outputs") / "mngs_to_specimen_outputs_rich_v4.json"
DEFAULT_PATIENT_ROOT = Path("outputs") / "patient_info_rich_normalized"
DEFAULT_SUFFIX = "mngs_to_specimen_v4"
DEFAULT_REPORT = Path("outputs") / "mngs_to_specimen_outputs_rich_v4_split_report.json"


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def split_outputs(input_json: Path, patient_root: Path, suffix: str, report_json: Path) -> dict[str, Any]:
    payload = load_json(input_json)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict payload: {input_json}")

    written: list[dict[str, Any]] = []
    missing_dirs: list[str] = []
    for patient_id, specimen_outputs in sorted(payload.items(), key=lambda item: int(item[0]) if str(item[0]).isdigit() else 999999):
        patient_id_text = str(patient_id)
        patient_dir = patient_root / f"NGS_patient_{patient_id_text}_json"
        if not patient_dir.exists():
            missing_dirs.append(patient_id_text)
        summary_dir = patient_dir / "summary_outputs"
        output_path = summary_dir / f"NGS_patient_{patient_id_text}_{suffix}.json"
        write_json(output_path, specimen_outputs)
        written.append(
            {
                "patient_id": patient_id_text,
                "specimen_count": len(specimen_outputs) if isinstance(specimen_outputs, list) else 0,
                "output_path": str(output_path),
                "patient_dir_exists": patient_dir.exists(),
            }
        )

    report = {
        "input_json": str(input_json),
        "patient_root": str(patient_root),
        "suffix": suffix,
        "patient_count": len(written),
        "specimen_count": sum(item["specimen_count"] for item in written),
        "missing_patient_dirs": missing_dirs,
        "written": written,
    }
    write_json(report_json, report)
    return report


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split a combined mNGS-to-specimen output JSON into per-patient summary_outputs files."
    )
    parser.add_argument("--input-json", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--suffix", default=DEFAULT_SUFFIX)
    parser.add_argument("--report-json", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report = split_outputs(args.input_json, args.patient_root, args.suffix, args.report_json)
    print(f"patient_count={report['patient_count']}")
    print(f"specimen_count={report['specimen_count']}")
    print(f"report={args.report_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

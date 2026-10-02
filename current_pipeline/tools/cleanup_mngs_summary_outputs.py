from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Sequence


DEFAULT_PATIENT_ROOT = Path("outputs") / "patient_info_rich_normalized"
DEFAULT_REPORT = Path("outputs") / "summary_outputs_cleanup_report.json"

ARCHIVE_LLM_PATTERNS = (
    "*_mNGS_max_agent_full.json",
    "*_mNGS_max_agent_all_rk_ntc_opt_chosen_full.json",
)

DELETE_TEST_PATTERNS = (
    "*_tag_explained_once.json",
    "*_tag_explained_once2.json",
    "*_explanation_once.prompt.md",
    "*_explanation_once2.prompt.md",
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Clean mNGS summary_outputs by archiving old LLM max results and "
            "removing temporary once/once2 explanation files."
        )
    )
    parser.add_argument("patient_root", nargs="?", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--patients", nargs="*", help="Optional patient IDs to clean.")
    parser.add_argument("--archive-dir-name", default="archive_llm_baseline")
    parser.add_argument("--report-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite-archive", action="store_true")
    return parser.parse_args(argv)


def patient_id_from_dir(patient_dir: Path) -> str:
    parts = patient_dir.name.split("_")
    if len(parts) >= 3 and parts[0] == "NGS" and parts[1] == "patient":
        return parts[2]
    raise ValueError(f"Cannot parse patient ID from {patient_dir}")


def patient_dirs(patient_root: Path, patients: Sequence[str] | None) -> list[Path]:
    wanted = {str(item).removeprefix("NGS_patient_").removesuffix("_json") for item in patients or []}
    dirs = sorted(
        (path for path in patient_root.glob("NGS_patient_*_json") if path.is_dir()),
        key=lambda path: int(patient_id_from_dir(path)),
    )
    if wanted:
        dirs = [path for path in dirs if patient_id_from_dir(path) in wanted]
    return dirs


def unique_files(summary_dir: Path, patterns: Sequence[str]) -> list[Path]:
    files: dict[Path, None] = {}
    for pattern in patterns:
        for path in summary_dir.glob(pattern):
            if path.is_file():
                files[path] = None
    return sorted(files)


def archive_file(path: Path, archive_dir: Path, *, dry_run: bool, overwrite: bool) -> dict[str, str]:
    destination = archive_dir / path.name
    if destination.exists() and not overwrite:
        return {
            "action": "archive_skipped_exists",
            "source": str(path),
            "destination": str(destination),
        }
    if not dry_run:
        archive_dir.mkdir(parents=True, exist_ok=True)
        if destination.exists() and overwrite:
            destination.unlink()
        shutil.move(str(path), str(destination))
    return {
        "action": "archive",
        "source": str(path),
        "destination": str(destination),
    }


def delete_file(path: Path, *, dry_run: bool) -> dict[str, str]:
    if not dry_run:
        path.unlink()
    return {
        "action": "delete",
        "source": str(path),
    }


def clean_patient(
    patient_dir: Path,
    *,
    archive_dir_name: str,
    dry_run: bool,
    overwrite_archive: bool,
) -> dict[str, object]:
    summary_dir = patient_dir / "summary_outputs"
    result: dict[str, object] = {
        "patient_id": patient_id_from_dir(patient_dir),
        "summary_dir": str(summary_dir),
        "archived": [],
        "deleted": [],
        "skipped": [],
    }
    if not summary_dir.exists():
        result["skipped"] = ["summary_outputs_missing"]
        return result

    archive_dir = summary_dir / archive_dir_name
    for path in unique_files(summary_dir, ARCHIVE_LLM_PATTERNS):
        operation = archive_file(
            path,
            archive_dir,
            dry_run=dry_run,
            overwrite=overwrite_archive,
        )
        target = "skipped" if operation["action"] == "archive_skipped_exists" else "archived"
        result[target].append(operation)  # type: ignore[index, union-attr]

    for path in unique_files(summary_dir, DELETE_TEST_PATTERNS):
        result["deleted"].append(delete_file(path, dry_run=dry_run))  # type: ignore[index, union-attr]

    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    results = [
        clean_patient(
            patient_dir,
            archive_dir_name=args.archive_dir_name,
            dry_run=args.dry_run,
            overwrite_archive=args.overwrite_archive,
        )
        for patient_dir in patient_dirs(args.patient_root, args.patients)
    ]
    summary = {
        "patient_root": str(args.patient_root),
        "dry_run": args.dry_run,
        "archive_patterns": ARCHIVE_LLM_PATTERNS,
        "delete_patterns": DELETE_TEST_PATTERNS,
        "patients_processed": len(results),
        "archived_count": sum(len(item.get("archived", [])) for item in results),
        "deleted_count": sum(len(item.get("deleted", [])) for item in results),
        "skipped_count": sum(len(item.get("skipped", [])) for item in results),
    }
    report = {"summary": summary, "patients": results}
    args.report_output.parent.mkdir(parents=True, exist_ok=True)
    args.report_output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Wrote cleanup report: {args.report_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

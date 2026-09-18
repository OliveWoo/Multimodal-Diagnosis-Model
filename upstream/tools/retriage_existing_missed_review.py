from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import mngs_common as mngs  # noqa: E402
from tools.review_missed_mngs_candidates import (  # noqa: E402
    DEFAULT_MAX_SUFFIX,
    DEFAULT_MNGS_TO_SPECIMEN_SUFFIX,
    DEFAULT_QUEUE_SUFFIX,
    build_review_payload,
    load_json,
    normalize_direct_review_output,
    parse_patient_ids,
    path_for_suffix,
    write_json,
)


def optional_path(value: Any) -> Path | None:
    if not value:
        return None
    path = Path(str(value))
    return path if path.exists() else None


def source_or_default(source_files: dict[str, Any], key: str, default: Path) -> Path:
    path = optional_path(source_files.get(key))
    return path if path is not None else default


def process_patient(
    patient_dir: Path,
    *,
    input_review_suffix: str,
    output_review_suffix: str,
    queue_suffix: str,
    max_suffix: str,
    mngs_to_specimen_suffix: str,
    artifact_root: Path | None,
    overwrite: bool,
    skip_existing: bool,
) -> Path:
    input_path = path_for_suffix(patient_dir, input_review_suffix, artifact_root=artifact_root)
    output_path = path_for_suffix(patient_dir, output_review_suffix, artifact_root=artifact_root)
    if output_path.exists() and skip_existing:
        return output_path
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output_path}. Use --overwrite or --skip-existing.")
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    review = load_json(input_path)
    if not isinstance(review, dict):
        raise RuntimeError(f"Review JSON is not an object: {input_path}")

    source_files = review.get("source_files")
    if not isinstance(source_files, dict):
        source_files = {}

    queue_path = source_or_default(
        source_files,
        "queue",
        path_for_suffix(patient_dir, queue_suffix, artifact_root=artifact_root),
    )
    max_path = source_or_default(
        source_files,
        "deterministic_max",
        path_for_suffix(patient_dir, max_suffix),
    )
    summary_path = source_or_default(
        source_files,
        "deterministic_summary",
        path_for_suffix(patient_dir, "final_summary_with_filmarray_deterministic"),
    )
    mngs_to_specimen_path = optional_path(source_files.get("mngs_to_specimen"))
    if mngs_to_specimen_path is None:
        default_mngs_path = path_for_suffix(patient_dir, mngs_to_specimen_suffix)
        mngs_to_specimen_path = default_mngs_path if default_mngs_path.exists() else None

    payload = build_review_payload(
        patient_dir=patient_dir,
        queue_path=queue_path,
        max_path=max_path,
        summary_path=summary_path,
        mngs_to_specimen_path=mngs_to_specimen_path,
    )
    output = normalize_direct_review_output(review, payload)
    output["retriage_source_review"] = str(input_path)
    output["retriage_note"] = (
        "Existing LLM review was re-normalized with the current direct triage rules; "
        "no LLM call was made by retriage_existing_missed_review.py."
    )
    write_json(output_path, output)
    return output_path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Re-apply current direct missed-review triage rules to existing LLM review JSON files without calling an LLM."
    )
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--patients", nargs="*")
    parser.add_argument("--input-review-suffix", required=True)
    parser.add_argument("--output-review-suffix", required=True)
    parser.add_argument("--queue-suffix", default=DEFAULT_QUEUE_SUFFIX)
    parser.add_argument("--max-suffix", default=DEFAULT_MAX_SUFFIX)
    parser.add_argument("--mngs-to-specimen-suffix", default=DEFAULT_MNGS_TO_SPECIMEN_SUFFIX)
    parser.add_argument("--artifact-root", type=Path)
    parser.add_argument("--summary-mode", choices=["deterministic", "full"], default="deterministic")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    patient_dirs = mngs.collect_patient_directories(
        [args.patient_root],
        summary_mode=args.summary_mode,
        mngs_source="ranked_only",
    )
    requested = parse_patient_ids(args.patients)
    if requested:
        patient_dirs = [
            patient_dir for patient_dir in patient_dirs if mngs.extract_patient_identifier(patient_dir) in requested
        ]

    failures: list[str] = []
    written = 0
    for patient_dir in patient_dirs:
        try:
            output_path = process_patient(
                patient_dir,
                input_review_suffix=args.input_review_suffix,
                output_review_suffix=args.output_review_suffix,
                queue_suffix=args.queue_suffix,
                max_suffix=args.max_suffix,
                mngs_to_specimen_suffix=args.mngs_to_specimen_suffix,
                artifact_root=args.artifact_root,
                overwrite=args.overwrite,
                skip_existing=args.skip_existing,
            )
            written += 1
            print(f"wrote retriaged review {output_path}")
        except Exception as exc:  # noqa: BLE001
            failures.append(f"{patient_dir}: {exc}")
            print(f"failed {patient_dir}: {exc}", file=sys.stderr)

    print(f"processed {written} patient(s)")
    if failures:
        print("failures:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

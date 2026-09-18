#!/usr/bin/env python3
"""Run the candidate-evidence workflow with one command.

This wrapper keeps the three original processing scripts independently usable,
while providing a safer handoff entry point with input checks and a predictable
output layout.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build candidate evidence, then create either the full time-window "
            "dataset, per-test candidate evidence, or both."
        )
    )
    parser.add_argument(
        "--input-root",
        required=True,
        type=Path,
        help="Root containing NGS_patient_*_json directories.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        help="A new or empty directory for this run.",
    )
    parser.add_argument(
        "--mngs-root",
        type=Path,
        default=None,
        help="Optional root containing per-patient grouped mNGS JSON files.",
    )
    parser.add_argument(
        "--candidate-window-days",
        type=float,
        default=2,
        help="Days before and after the anchor used to select candidates.",
    )
    parser.add_argument(
        "--window-mode",
        choices=("elapsed", "date"),
        default="date",
        help="Use exact elapsed hours or inclusive calendar dates.",
    )
    parser.add_argument(
        "--anchor-source",
        choices=("auto", "input_dates", "mngs"),
        default="auto",
        help="Source used to choose the candidate-window anchor.",
    )
    parser.add_argument(
        "--include-mngs-pathogens",
        action="store_true",
        help="Also treat pathogens in the selected mNGS record as candidates.",
    )
    parser.add_argument(
        "--selected-root",
        type=Path,
        default=None,
        help="Directory containing per-patient selected-pathogens packages.",
    )
    parser.add_argument(
        "--candidate-source",
        choices=("window", "selected", "upstream_all", "union"),
        default="window",
        help="Choose candidates from the time window, frozen upstream output, or both.",
    )
    parser.add_argument(
        "--selected-missing-policy",
        choices=("error", "skip"),
        default="error",
        help="How to handle patients without a selected-pathogens package.",
    )
    parser.add_argument(
        "--rationale-views",
        action="store_true",
        help="Write a simple selected-pathogen table and traceable reasoning-chain JSON.",
    )
    parser.add_argument(
        "--mode",
        choices=("combined", "split", "both"),
        default="combined",
        help=(
            "combined keeps window records plus all-time candidate evidence; "
            "split writes only candidate-matching microbiology by test type."
        ),
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print generated JSON.",
    )
    parser.add_argument(
        "--write-empty-files",
        action="store_true",
        help="In combined mode, write [] for test files with no retained records.",
    )
    return parser.parse_args()


def require_directory(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_dir():
        raise SystemExit(f"{label} is not a directory: {resolved}")
    return resolved


def prepare_output_root(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved.exists() and not resolved.is_dir():
        raise SystemExit(f"output root is not a directory: {resolved}")
    if resolved.is_dir() and any(resolved.iterdir()):
        raise SystemExit(
            "output root must be new or empty so an earlier run is not mixed "
            f"with this run: {resolved}"
        )
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def run_checked(command: list[str]) -> None:
    print("\n$ " + " ".join(command), flush=True)
    subprocess.run(command, check=True)


def main() -> int:
    args = parse_args()
    input_root = require_directory(args.input_root, "input root")
    if args.candidate_window_days < 0:
        raise SystemExit("candidate window days must be zero or greater")
    if args.window_mode == "date" and not args.candidate_window_days.is_integer():
        raise SystemExit("date window mode requires a whole number of days")
    if not any(
        path.is_dir()
        and path.name.startswith("NGS_patient_")
        and path.name.endswith("_json")
        for path in input_root.iterdir()
    ):
        raise SystemExit(
            "input root contains no NGS_patient_<ID>_json directories: "
            f"{input_root}"
        )
    mngs_root = (
        require_directory(args.mngs_root, "mNGS root") if args.mngs_root else None
    )
    selected_root = (
        require_directory(args.selected_root, "selected root")
        if args.selected_root
        else None
    )
    if args.candidate_source != "window" and selected_root is None:
        raise SystemExit(f"candidate-source {args.candidate_source} requires --selected-root")
    if args.rationale_views and selected_root is None:
        raise SystemExit("--rationale-views requires --selected-root")
    output_root = prepare_output_root(args.output_root)

    candidate_root = output_root / "candidate_evidence"
    candidate_command = [
        sys.executable,
        str(SCRIPT_DIR / "build_candidate_window_evidence.py"),
        "--input-root",
        str(input_root),
        "--output-root",
        str(candidate_root),
        "--candidate-window-days",
        str(args.candidate_window_days),
        "--window-mode",
        args.window_mode,
        "--anchor-source",
        args.anchor_source,
        "--candidate-source",
        args.candidate_source,
        "--selected-missing-policy",
        args.selected_missing_policy,
    ]
    if mngs_root:
        candidate_command.extend(["--mngs-root", str(mngs_root)])
    if args.include_mngs_pathogens:
        candidate_command.append("--include-mngs-pathogens")
    if selected_root:
        candidate_command.extend(["--selected-root", str(selected_root)])
    if args.pretty:
        candidate_command.append("--pretty")
    run_checked(candidate_command)

    if args.mode in {"combined", "both"}:
        combined_root = output_root / "window_plus_candidate_evidence"
        combined_command = [
            sys.executable,
            str(SCRIPT_DIR / "build_mngs_window_plus_candidate_evidence.py"),
            "--input-root",
            str(input_root),
            "--candidate-root",
            str(candidate_root),
            "--output-root",
            str(combined_root),
        ]
        if args.pretty:
            combined_command.append("--pretty")
        if args.write_empty_files:
            combined_command.append("--write-empty-files")
        run_checked(combined_command)

    if args.mode in {"split", "both"}:
        split_root = output_root / "candidate_evidence_by_test"
        split_command = [
            sys.executable,
            str(SCRIPT_DIR / "split_candidate_evidence_by_test.py"),
            "--input-root",
            str(input_root),
            "--candidate-root",
            str(candidate_root),
            "--output-root",
            str(split_root),
        ]
        if args.pretty:
            split_command.append("--pretty")
        run_checked(split_command)

    if args.rationale_views:
        rationale_command = [
            sys.executable,
            str(SCRIPT_DIR / "build_rationale_views.py"),
            "--candidate-root",
            str(candidate_root),
            "--output-root",
            str(output_root / "rationale_views"),
        ]
        if args.pretty:
            rationale_command.append("--pretty")
        run_checked(rationale_command)

    print(f"\nPipeline completed. Output: {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

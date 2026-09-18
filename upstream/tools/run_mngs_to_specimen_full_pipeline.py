from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Sequence

DEFAULT_RAW_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs.json"
DEFAULT_PREOPT_REPORT = Path("outputs") / "mngs_to_specimen_vs_excel_report_preopt.json"
DEFAULT_OPTIMIZED_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs.json"
DEFAULT_FINAL_REPORT = Path("outputs") / "mngs_to_specimen_vs_excel_report.json"
DEFAULT_OPTIMIZATION_CONFIG = Path("outputs") / "mngs_low_colonizer_optimization.json"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "One-command pipeline: generate mNGS-to-specimen output, compare with Excel, "
            "optimize low_colonizer under recall threshold, and write final comparison JSON."
        )
    )
    parser.add_argument("--model", default="gpt-5", help="OpenAI model name (default: gpt-5)")
    parser.add_argument("--patient-id", help="Optional patient id for quick test run")
    parser.add_argument("--specimen-excel", type=Path, default=None, help="Path to specimen Excel")
    parser.add_argument(
        "--min-recall",
        type=float,
        default=0.95,
        help="Recall floor for low_colonizer optimization (default: 0.95)",
    )
    parser.add_argument(
        "--scope",
        choices=("by_specimen", "by_patient"),
        default="by_specimen",
        help="Scope used to enforce min recall during optimization",
    )
    parser.add_argument(
        "--min-output-only-count",
        type=int,
        default=1,
        help=(
            "Only consider species with at least this many output_only occurrences "
            "during low_colonizer optimization (default: 1)."
        ),
    )
    parser.add_argument(
        "--allow-remove-matched",
        action="store_true",
        help=(
            "Allow optimizer to remove species that also have matched occurrences, "
            "while still enforcing min recall."
        ),
    )
    parser.add_argument(
        "--opt-granularity",
        choices=("specimen_name", "name"),
        default="specimen_name",
        help=(
            "Low_colonizer optimization unit. specimen_name is recommended because "
            "it can remove false positives only in specific specimens."
        ),
    )
    parser.add_argument(
        "--opt-include-medium",
        action="store_true",
        help=(
            "Also optimize medium tier (in addition to low_colonizer). "
            "Use this for stronger output_only reduction while keeping min_recall."
        ),
    )
    parser.add_argument("--raw-output-json", type=Path, default=DEFAULT_RAW_OUTPUT)
    parser.add_argument("--preopt-report-json", type=Path, default=DEFAULT_PREOPT_REPORT)
    parser.add_argument("--optimized-output-json", type=Path, default=DEFAULT_OPTIMIZED_OUTPUT)
    parser.add_argument("--final-report-json", type=Path, default=DEFAULT_FINAL_REPORT)
    parser.add_argument("--optimization-json", type=Path, default=DEFAULT_OPTIMIZATION_CONFIG)
    parser.add_argument(
        "--skip-optimize",
        action="store_true",
        help="Skip optimization step and keep final report from raw output",
    )
    return parser.parse_args(argv)


def _run(cmd: list[str]) -> None:
    print("$ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def _same_path(lhs: Path, rhs: Path) -> bool:
    return lhs.resolve(strict=False) == rhs.resolve(strict=False)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    py = sys.executable

    run_agent_cmd = [
        py,
        "-m",
        "tools.run_mngs_to_specimen_agent",
        "--model",
        args.model,
        "--output",
        str(args.raw_output_json),
    ]
    if args.patient_id:
        run_agent_cmd.extend(["--patient-id", str(args.patient_id)])
    _run(run_agent_cmd)

    compare_raw_cmd = [
        py,
        "-m",
        "tools.compare_mngs_to_specimen_with_excel",
        "--input-json",
        str(args.raw_output_json),
        "--output-json",
        str(args.preopt_report_json),
    ]
    if args.specimen_excel:
        compare_raw_cmd.extend(["--specimen-excel", str(args.specimen_excel)])
    _run(compare_raw_cmd)

    if args.skip_optimize:
        compare_final_cmd = [
            py,
            "-m",
            "tools.compare_mngs_to_specimen_with_excel",
            "--input-json",
            str(args.raw_output_json),
            "--output-json",
            str(args.final_report_json),
        ]
        if args.specimen_excel:
            compare_final_cmd.extend(["--specimen-excel", str(args.specimen_excel)])
        _run(compare_final_cmd)
        print(f"Final comparison JSON: {args.final_report_json}")
        return

    optimize_cmd = [
        py,
        "-m",
        "tools.optimize_low_colonizer_precision",
        "--input-json",
        str(args.raw_output_json),
        "--min-recall",
        str(args.min_recall),
        "--scope",
        args.scope,
        "--view",
        "all_levels",
        "--granularity",
        args.opt_granularity,
        "--min-output-only-count",
        str(args.min_output_only_count),
        "--output-json",
        str(args.optimized_output_json),
        "--report-json",
        str(args.final_report_json),
        "--config-json",
        str(args.optimization_json),
    ]
    if args.allow_remove_matched:
        optimize_cmd.append("--allow-remove-matched")
    if args.opt_include_medium:
        optimize_cmd.append("--include-medium")
    if args.specimen_excel:
        optimize_cmd.extend(["--specimen-excel", str(args.specimen_excel)])
    _run(optimize_cmd)

    if _same_path(args.raw_output_json, args.optimized_output_json):
        print(f"Final optimized output JSON: {args.optimized_output_json}")
    else:
        print(f"Raw output JSON: {args.raw_output_json}")
        print(f"Optimized output JSON: {args.optimized_output_json}")
    print(f"Raw report JSON: {args.preopt_report_json}")
    print(f"Final report JSON: {args.final_report_json}")
    print(f"Optimization summary JSON: {args.optimization_json}")


if __name__ == "__main__":
    main()

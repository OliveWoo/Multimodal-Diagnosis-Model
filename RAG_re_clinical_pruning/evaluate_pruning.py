#!/usr/bin/env python
"""Evaluate locked clinical-pruning artifacts without mutating generation output."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from rag_re_clinical_pruning.evaluation import (
    EvaluationError,
    evaluate_directory,
    write_json,
    write_markdown,
)


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_PROTOCOL = PROJECT_ROOT / "config" / "evaluation_protocol_v1.json"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Score rag_re_clinical_pruning.output.v1 artifacts against frozen gold. "
            "Recall-safety checks are reporting-only and never change predictions."
        )
    )
    parser.add_argument("--outputs", required=True, type=Path, help="Pruning artifact directory.")
    parser.add_argument("--gold", required=True, type=Path, help="Frozen clinical gold JSON.")
    parser.add_argument(
        "--protocol",
        type=Path,
        default=DEFAULT_PROTOCOL,
        help=f"Frozen evaluation protocol (default: {DEFAULT_PROTOCOL}).",
    )
    parser.add_argument("--json-out", required=True, type=Path, help="Evaluation JSON output.")
    parser.add_argument(
        "--markdown-out", required=True, type=Path, help="Evaluation Markdown output."
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress the success summary.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        report = evaluate_directory(
            args.outputs,
            args.gold,
            protocol_path=args.protocol,
        )
        json_path = write_json(report, args.json_out)
        markdown_path = write_markdown(report, args.markdown_out)
    except (EvaluationError, OSError) as exc:
        parser.exit(2, f"error: {exc}\n")
    if not args.quiet:
        missing = report["artifacts"]["missing_labeled_artifact_count"]
        primary = report["primary_arm"]
        safety = report["recall_safety_summary"].get(primary, {}).get("status", "NA")
        print(
            "clinical pruning evaluation complete: "
            f"labeled={report['gold']['labeled_case_count']}, "
            f"artifacts={report['artifacts']['artifact_count']}, "
            f"missing={missing}, primary={primary}, safety={safety}"
        )
        print(f"json: {json_path.resolve()}")
        print(f"markdown: {markdown_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


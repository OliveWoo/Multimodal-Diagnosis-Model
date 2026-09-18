#!/usr/bin/env python
"""Evaluate a completed post-L5 A/B/C replay without feeding gold to generation."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from rag_re_clinical_pruning_abc.evaluation import (
    EvaluationError,
    evaluate_directory,
    write_json,
    write_markdown,
)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Gold-isolated post-L5 A/B/C evaluation")
    parser.add_argument("--outputs", required=True, help="Completed ABC replay directory")
    parser.add_argument("--gold", required=True, help="Clinical gold JSON")
    parser.add_argument(
        "--protocol",
        default=str(root / "config" / "evaluation_protocol_v1.json"),
        help="Frozen evaluation protocol JSON",
    )
    parser.add_argument("--json-out", required=True, help="Output JSON report")
    parser.add_argument("--markdown-out", required=True, help="Output Markdown report")
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=None,
        help="Override only for smoke/tests; omit for the frozen protocol value",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = evaluate_directory(
            args.outputs,
            args.gold,
            args.protocol,
            bootstrap_replicates=args.bootstrap_replicates,
        )
        write_json(report, args.json_out)
        write_markdown(report, args.markdown_out)
    except EvaluationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    primary = report["scopes"][report["primary_scope"]]["arms"][report["primary_anchor"]]
    print(
        json.dumps(
            {
                "status": "complete",
                "patient_count": report["gold"]["labeled_case_count"],
                "artifact_count": report["artifacts"]["artifact_count"],
                "arm_count": report["arm_count"],
                "primary_anchor": report["primary_anchor"],
                "primary_auto_precision": primary["auto"]["precision"],
                "primary_auto_recall": primary["auto"]["recall"],
                "json_report": str(Path(args.json_out).resolve()),
                "markdown_report": str(Path(args.markdown_out).resolve()),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


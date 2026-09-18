"""Command-line interface for the independent clinical replay."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from .engine import load_config, run_patient, write_json
from .io_adapter import discover_jsons, load_json


def _pid(value: dict[str, Any], path: Path) -> str:
    patient = value.get("patient_id")
    if patient is None:
        raise ValueError(f"Missing patient_id: {path}")
    patient = str(patient).strip()
    if not patient:
        raise ValueError(f"Empty patient_id: {path}")
    return patient


def _index_inputs(directory: str | Path, *, expected_schema: str | None = None) -> dict[str, tuple[Path, dict[str, Any]]]:
    index: dict[str, tuple[Path, dict[str, Any]]] = {}
    for path in discover_jsons(directory):
        payload = load_json(path)
        if expected_schema and payload.get("schema_version") != expected_schema:
            continue
        patient = _pid(payload, path)
        if patient in index:
            raise ValueError(
                f"Duplicate patient_id {patient}: {index[patient][0]} and {path}"
            )
        index[patient] = (path, payload)
    return index


def replay_directory(
    frozen_outputs: str | Path,
    merge_inputs: str | Path,
    output: str | Path,
    config_path: str | Path,
    *,
    force: bool = False,
) -> dict[str, Any]:
    frozen = _index_inputs(frozen_outputs, expected_schema="rag_re.output.v1")
    merges = _index_inputs(merge_inputs)
    if not frozen:
        raise ValueError(f"No rag_re.output.v1 artifacts found in {frozen_outputs}")
    missing = sorted(set(frozen) - set(merges), key=lambda item: (len(item), item))
    if missing:
        raise ValueError(f"Missing merge inputs for patients: {', '.join(missing)}")

    config = load_config(config_path)
    destination = Path(output)
    destination.mkdir(parents=True, exist_ok=True)
    existing = list(destination.glob("*.json"))
    if existing and not force:
        raise FileExistsError(
            f"Output directory already contains JSON artifacts: {destination}; "
            "use a new directory or pass --force"
        )

    written: list[str] = []
    for patient in sorted(frozen, key=lambda item: (len(item), item)):
        frozen_path, frozen_payload = frozen[patient]
        merge_path, merge_payload = merges[patient]
        result = run_patient(
            frozen_payload,
            merge_payload,
            config,
            frozen_path=frozen_path,
            merge_path=merge_path,
        )
        out_path = destination / f"NGS_patient_{patient}_clinical_v1.json"
        write_json(out_path, result)
        written.append(str(out_path))

    manifest = {
        "schema_version": "rag_re_clinical.batch_manifest.v1",
        "run_status": "complete",
        "patient_count": len(written),
        "patients": sorted(frozen, key=lambda item: (len(item), item)),
        "artifacts": written,
        "config": str(Path(config_path)),
        "frozen_outputs": str(Path(frozen_outputs)),
        "merge_inputs": str(Path(merge_inputs)),
    }
    write_json(destination / "manifest.json", manifest)
    return manifest


def _evaluate(outputs: str | Path, gold: str | Path, report_json: str | Path, report_markdown: str | Path) -> dict[str, Any]:
    # Delayed import is intentional: replay generation never imports the gold reader.
    from .evaluation import evaluate_directory, write_markdown

    report = evaluate_directory(outputs, gold)
    write_json(report_json, report)
    write_markdown(report, report_markdown)
    return report


def _parser() -> argparse.ArgumentParser:
    default_config = Path(__file__).resolve().parents[1] / "config" / "clinical_v1.json"
    parser = argparse.ArgumentParser(
        description="Run clinical, baseline-preserving experiments beside RAG_re."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    replay = sub.add_parser("replay", help="Generate clinical artifacts without reading gold")
    replay.add_argument("--frozen-outputs", required=True)
    replay.add_argument("--merge-inputs", required=True)
    replay.add_argument("--output", required=True)
    replay.add_argument("--config", default=str(default_config))
    replay.add_argument("--force", action="store_true")

    evaluate = sub.add_parser("evaluate", help="Evaluate already generated artifacts")
    evaluate.add_argument("--outputs", required=True)
    evaluate.add_argument("--gold", required=True)
    evaluate.add_argument("--report-json", required=True)
    evaluate.add_argument("--report-markdown", required=True)

    all_command = sub.add_parser("run-all", help="Replay, then evaluate in a separate phase")
    all_command.add_argument("--frozen-outputs", required=True)
    all_command.add_argument("--merge-inputs", required=True)
    all_command.add_argument("--gold", required=True)
    all_command.add_argument("--output", required=True)
    all_command.add_argument("--report-json", required=True)
    all_command.add_argument("--report-markdown", required=True)
    all_command.add_argument("--config", default=str(default_config))
    all_command.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "replay":
            manifest = replay_directory(
                args.frozen_outputs,
                args.merge_inputs,
                args.output,
                args.config,
                force=args.force,
            )
            print(json.dumps(manifest, ensure_ascii=False, indent=2))
            return 0
        if args.command == "evaluate":
            report = _evaluate(
                args.outputs,
                args.gold,
                args.report_json,
                args.report_markdown,
            )
            print(json.dumps(report.get("summary", report), ensure_ascii=False, indent=2))
            return 0
        manifest = replay_directory(
            args.frozen_outputs,
            args.merge_inputs,
            args.output,
            args.config,
            force=args.force,
        )
        report = _evaluate(
            args.output,
            args.gold,
            args.report_json,
            args.report_markdown,
        )
        print(
            json.dumps(
                {"replay": manifest, "evaluation_summary": report.get("summary", report)},
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0
    except Exception as exc:  # CLI boundary: keep traceback opt-in, return non-zero.
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

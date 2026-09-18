from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

from .config import config_hash, load_config, load_env_file, validate_config
from .engine import RagReEngine
from .evaluation import evaluate_outputs, load_gold, load_outputs, render_markdown
from .io_utils import load_json, sha256_file, stable_hash, write_json_atomic
from .versioning import (
    EVALUATOR_FINGERPRINT,
    EVALUATOR_VERSION,
    PIPELINE_FINGERPRINT,
)


def _path(value: str | None) -> Path | None:
    return Path(value).resolve() if value else None


def _load_runtime_config(
    args: argparse.Namespace, *, base: Dict[str, Any] | None = None
) -> Dict[str, Any]:
    env_file = _path(getattr(args, "env_file", None))
    if env_file is not None:
        load_env_file(env_file)
    return load_config(_path(getattr(args, "config", None)), base=base)


def _filter_gold(gold: Dict[str, Any], split: str | None) -> Dict[str, Any]:
    if not split:
        return gold
    filtered = dict(gold)
    filtered["cases"] = [row for row in gold["cases"] if row.get("split") == split]
    if not filtered["cases"]:
        raise ValueError(f"Gold file has no cases with split={split!r}")
    return filtered


def _default_cache(output_path: Path, explicit: str | None) -> Path:
    return _path(explicit) or output_path.parent / "_rag_re_cache"


def command_run(args: argparse.Namespace) -> int:
    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    config = _load_runtime_config(args)
    engine = RagReEngine(
        config,
        cache_dir=_default_cache(output_path, args.cache_dir),
        skip_literature=args.skip_literature,
    )
    result = engine.run_file(input_path)
    write_json_atomic(output_path, result)
    print(
        json.dumps(
            {
                "output": str(output_path),
                "patient_id": result["patient_id"],
                "run_status": result["run_status"],
                "candidate_count": len(result["candidates"]),
            },
            ensure_ascii=False,
        )
    )
    return 0 if result["run_status"] == "complete" else 2


def _valid_existing_output(
    path: Path,
    *,
    input_path: Path,
    config: Dict[str, Any],
    skip_literature: bool,
) -> bool:
    if not path.exists():
        return False
    try:
        payload = load_json(path)
    except Exception:
        return False
    return bool(
        isinstance(payload, dict)
        and payload.get("schema_version") == "rag_re.output.v1"
        and payload.get("run_status") == "complete"
        and payload.get("pipeline_fingerprint") == PIPELINE_FINGERPRINT
        and (payload.get("input_meta") or {}).get("input_sha256") == sha256_file(input_path)
        and payload.get("config_hash") == config_hash(config)
        and (payload.get("execution") or {}).get("skip_literature") is skip_literature
    )


def command_batch(args: argparse.Namespace) -> int:
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Batch input directory does not exist: {input_dir}")
    cache_dir = _path(args.cache_dir) or output_dir / "_rag_re_cache"
    if output_dir == input_dir or input_dir in output_dir.parents:
        raise ValueError("Batch output directory must not be inside the input directory")
    if cache_dir == input_dir or input_dir in cache_dir.parents:
        raise ValueError("Batch cache directory must not be inside the input directory")
    paths = sorted(input_dir.rglob(args.pattern) if args.recursive else input_dir.glob(args.pattern))
    if not paths:
        raise FileNotFoundError(
            f"No input files matched pattern {args.pattern!r} under {input_dir}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    config = _load_runtime_config(args)
    engine = RagReEngine(
        config,
        cache_dir=cache_dir,
        skip_literature=args.skip_literature,
    )
    completed = skipped = 0
    errors: List[Dict[str, str]] = []
    used_targets: Dict[Path, Path] = {}
    for input_path in paths:
        target = output_dir / f"{input_path.stem}_rag_re_out.json"
        if target in used_targets and used_targets[target] != input_path:
            errors.append(
                {
                    "input": str(input_path),
                    "error": f"Output name collision with {used_targets[target]}",
                }
            )
            continue
        used_targets[target] = input_path
        if args.resume and _valid_existing_output(
            target,
            input_path=input_path,
            config=config,
            skip_literature=args.skip_literature,
        ):
            skipped += 1
            continue
        try:
            result = engine.run_file(input_path)
            write_json_atomic(target, result)
            completed += 1
            if result["run_status"] != "complete":
                errors.append(
                    {
                        "input": str(input_path),
                        "output": str(target),
                        "error": "Run completed with candidate-level errors; inspect output.errors.",
                    }
                )
        except Exception as exc:
            errors.append(
                {"input": str(input_path), "error": f"{type(exc).__name__}: {exc}"}
            )
            if args.fail_fast:
                break
    manifest = {
        "schema_version": "rag_re.batch.v1",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "matched_files": len(paths),
        "completed": completed,
        "skipped": skipped,
        "error_count": len(errors),
        "errors": errors,
    }
    write_json_atomic(output_dir / "batch_manifest.json", manifest)
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0 if not errors else 2


def command_replay(args: argparse.Namespace) -> int:
    source_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    source = load_json(source_path)
    if not isinstance(source, dict):
        raise ValueError("Replay input must be a JSON object")
    config = _load_runtime_config(args, base=source.get("config_snapshot"))
    # No provider is called during replay. skip_literature avoids requiring an API key.
    engine = RagReEngine(config, skip_literature=True)
    result = engine.replay(source)
    write_json_atomic(output_path, result)
    print(
        json.dumps(
            {
                "output": str(output_path),
                "patient_id": result["patient_id"],
                "run_status": result["run_status"],
            }
        )
    )
    return 0 if result["run_status"] == "complete" else 2


def command_evaluate(args: argparse.Namespace) -> int:
    predictions = Path(args.predictions).resolve()
    gold_path = Path(args.gold).resolve()
    output_path = Path(args.output).resolve()
    outputs = load_outputs(predictions)
    if not outputs:
        raise ValueError(f"No rag_re.output.v1 files found under {predictions}")
    gold = _filter_gold(load_gold(gold_path), args.split)
    report = evaluate_outputs(
        outputs,
        gold,
        bootstrap_iterations=args.bootstrap_iterations,
        bootstrap_seed=args.bootstrap_seed,
        allow_missing_outputs=args.exploratory_allow_missing,
        allow_partial_outputs=args.exploratory_allow_partial,
        allow_mixed_provenance=args.exploratory_allow_mixed_provenance,
        analysis_mode="bc_only" if args.bc_only else "full_abc",
    )
    write_json_atomic(output_path, report)
    markdown_path = _path(args.markdown) or output_path.with_suffix(".md")
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_path.write_text(render_markdown(report), encoding="utf-8", newline="\n")
    print(
        json.dumps(
            {
                "output": str(output_path),
                "markdown": str(markdown_path),
                "evaluated_outputs": report["evaluated_output_count"],
            },
            ensure_ascii=False,
        )
    )
    return 0


def command_sweep(args: argparse.Namespace) -> int:
    source_outputs = load_outputs(Path(args.predictions).resolve())
    if not source_outputs:
        raise ValueError("No rag_re.output.v1 artifacts found for threshold sweep")
    raw_gold = load_gold(Path(args.gold).resolve())
    if not 0 <= args.recall_floor <= 1:
        raise ValueError("--recall-floor must be within 0..1")
    if not (
        raw_gold.get("development_only")
        or args.split == "development"
        or args.allow_nondevelopment_sweep
    ):
        raise ValueError(
            "Threshold sweep is restricted to a development split. Mark the gold file "
            "development_only=true, pass --split development, or explicitly use the unsafe "
            "--allow-nondevelopment-sweep override."
        )
    gold = _filter_gold(raw_gold, args.split)
    noncomplete = [
        str(payload.get("patient_id"))
        for payload in source_outputs
        if payload.get("run_status") != "complete"
    ]
    if noncomplete:
        raise ValueError(
            "Threshold sweep requires complete provider artifacts; invalid patient IDs: "
            + ", ".join(noncomplete)
        )
    if any(
        (payload.get("execution") or {}).get("literature_attempted") is not True
        for payload in source_outputs
    ):
        raise ValueError("Threshold sweep requires Module A to have been attempted")
    if any(
        not (payload.get("execution") or {}).get("literature_evidence_available")
        for payload in source_outputs
    ):
        raise ValueError("Threshold sweep requires cached Module-A evidence, not B/C-only runs")
    judgeable_votes = sum(
        judgment.get("status") == "ok"
        for payload in source_outputs
        for candidate in payload.get("candidates", [])
        for judgment in (candidate.get("literature_evidence") or {}).get("judgments", [])
        if isinstance(judgment, dict)
    )
    if judgeable_votes == 0:
        raise ValueError("Threshold sweep has zero judgeable Module-A article votes")
    base_config = _load_runtime_config(
        args, base=source_outputs[0].get("config_snapshot")
    )
    thresholds = [int(value.strip()) for value in args.thresholds.split(",") if value.strip()]
    if not thresholds:
        raise ValueError("At least one threshold is required")
    rows = []
    for threshold in thresholds:
        config = copy.deepcopy(base_config)
        config["module_a_literature"]["min_support"] = threshold
        validate_config(config)
        engine = RagReEngine(config, skip_literature=True)
        replayed = [engine.replay(payload) for payload in source_outputs]
        report = evaluate_outputs(replayed, gold, bootstrap_iterations=0)
        experiment = args.experiment
        if experiment not in report["experiments"]:
            raise ValueError(f"Unknown experiment for sweep: {experiment}")
        metrics = report["experiments"][experiment]["metrics"]
        rows.append({"min_support": threshold, **metrics})

    scorable = [row for row in rows if row.get("f2") is not None]
    if not scorable:
        raise ValueError(
            "Threshold calibration is undefined because the selected gold/experiment "
            "has no scorable positive-label metrics."
        )
    eligible = [
        row
        for row in scorable
        if row.get("end_to_end_recall") is not None
        and row["end_to_end_recall"] >= args.recall_floor
    ]
    if eligible:
        selected = max(
            eligible,
            key=lambda row: (
                row.get("agreement_precision") if row.get("agreement_precision") is not None else -1,
                row.get("f2") if row.get("f2") is not None else -1,
                row["min_support"],
            ),
        )
        selection_rule = "highest agreement precision subject to recall floor; ties by F2 then larger k"
        selection_status = "selected"
        failure_reason = None
        exploratory_best_f2 = None
    else:
        selected = None
        exploratory_best_f2 = max(
            scorable,
            key=lambda row: (
                row.get("f2") if row.get("f2") is not None else -1,
                row["min_support"],
            ),
        )
        max_ceiling = max(
            (
                row.get("candidate_pool_recall_ceiling")
                for row in scorable
                if row.get("candidate_pool_recall_ceiling") is not None
            ),
            default=None,
        )
        failure_reason = (
            "requested recall floor exceeds the observed candidate-pool recall ceiling"
            if max_ceiling is not None and args.recall_floor > max_ceiling
            else "no tested threshold met the pre-specified recall floor"
        )
        selection_rule = "no confirmatory threshold selected"
        selection_status = "failed_no_threshold_meets_recall_floor"
    result = {
        "schema_version": "rag_re.threshold_sweep.v1",
        "development_only": True,
        "experiment": args.experiment,
        "recall_floor": args.recall_floor,
        "selection_status": selection_status,
        "selection_failure_reason": failure_reason,
        "selection_rule": selection_rule,
        "selected_threshold": selected["min_support"] if selected else None,
        "selected_metrics": selected,
        "exploratory_best_f2_when_selection_failed": exploratory_best_f2,
        "thresholds": rows,
        "provenance": {
            "gold_sha256": gold.get("source_sha256"),
            "gold_case_ids": sorted(str(row["patient_id"]) for row in gold["cases"]),
            "gold_split": args.split,
            "evaluator_version": EVALUATOR_VERSION,
            "evaluator_fingerprint": EVALUATOR_FINGERPRINT,
            "base_replay_config_hash": config_hash(base_config),
            "source_pipeline_fingerprints": sorted(
                {str(payload.get("pipeline_fingerprint")) for payload in source_outputs}
            ),
            "source_prompt_sha256": sorted(
                {str(payload.get("prompt_sha256")) for payload in source_outputs}
            ),
            "source_evidence_config_hashes": sorted(
                {
                    str((payload.get("evidence_provenance") or {}).get("source_config_hash"))
                    for payload in source_outputs
                }
            ),
            "source_artifacts": sorted(
                [
                    {
                        "patient_id": str(payload.get("patient_id")),
                        "artifact_hash": payload.get("artifact_sha256_basis")
                        or stable_hash(
                            {
                                key: value
                                for key, value in payload.items()
                                if not str(key).startswith("_")
                            }
                        ),
                        "input_sha256": (payload.get("input_meta") or {}).get(
                            "input_sha256"
                        ),
                    }
                    for payload in source_outputs
                ],
                key=lambda row: (row["patient_id"], row["artifact_hash"]),
            ),
        },
        "warning": (
            "Freeze the selected threshold, prompt, model, aliases, candidate universe, and "
            "publication cutoff before evaluating an untouched holdout."
        ),
    }
    output = Path(args.output).resolve()
    write_json_atomic(output, result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if selected else 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rag-re",
        description="Simple, replayable A/B/C PubMed ablation pipeline for pathogen candidates.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_runtime_flags(command: argparse.ArgumentParser, *, env: bool = True) -> None:
        command.add_argument("--config", help="JSON config override; defaults are built in.")
        if env:
            command.add_argument("--env-file", help="Optional NAME=VALUE environment file.")

    run = subparsers.add_parser("run", help="Run one infection bundle.")
    run.add_argument("--input", required=True)
    run.add_argument("--output", required=True)
    run.add_argument("--cache-dir")
    run.add_argument("--skip-literature", action="store_true", help="Run E0/B/C without network or LLM.")
    add_runtime_flags(run)
    run.set_defaults(func=command_run)

    batch = subparsers.add_parser("batch", help="Run a folder of infection bundles.")
    batch.add_argument("--input-dir", required=True)
    batch.add_argument("--output-dir", required=True)
    batch.add_argument("--pattern", default="*.json")
    batch.add_argument("--recursive", action="store_true")
    batch.add_argument("--resume", action="store_true")
    batch.add_argument("--fail-fast", action="store_true")
    batch.add_argument("--cache-dir")
    batch.add_argument("--skip-literature", action="store_true")
    add_runtime_flags(batch)
    batch.set_defaults(func=command_batch)

    replay = subparsers.add_parser("replay", help="Reapply thresholds/rules to cached evidence.")
    replay.add_argument("--input", required=True)
    replay.add_argument("--output", required=True)
    add_runtime_flags(replay)
    replay.set_defaults(func=command_replay)

    evaluate = subparsers.add_parser("evaluate", help="Evaluate all ablation arms.")
    evaluate.add_argument("--predictions", required=True)
    evaluate.add_argument("--gold", required=True)
    evaluate.add_argument("--output", required=True)
    evaluate.add_argument("--markdown")
    evaluate.add_argument("--split", help="Evaluate only cases with this frozen split label.")
    evaluate.add_argument("--bootstrap-iterations", type=int, default=1000)
    evaluate.add_argument("--bootstrap-seed", type=int, default=20260810)
    evaluate.add_argument("--exploratory-allow-missing", action="store_true")
    evaluate.add_argument("--exploratory-allow-partial", action="store_true")
    evaluate.add_argument("--exploratory-allow-mixed-provenance", action="store_true")
    evaluate.add_argument(
        "--bc-only",
        action="store_true",
        help="Explicitly evaluate only E0/B/C/BC/E0+C when Module A was skipped.",
    )
    evaluate.set_defaults(func=command_evaluate)

    sweep = subparsers.add_parser("sweep", help="Sweep Module-A support thresholds offline.")
    sweep.add_argument("--predictions", required=True)
    sweep.add_argument("--gold", required=True)
    sweep.add_argument("--output", required=True)
    sweep.add_argument("--thresholds", default="1,2,3,4,5,6,7,8,9,10")
    sweep.add_argument("--experiment", default="E0_OR_C_OR_A_AND_B")
    sweep.add_argument("--recall-floor", type=float, default=0.90)
    sweep.add_argument("--split", help="Use only this development split from the gold file.")
    sweep.add_argument("--allow-nondevelopment-sweep", action="store_true")
    add_runtime_flags(sweep)
    sweep.set_defaults(func=command_sweep)
    return parser


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

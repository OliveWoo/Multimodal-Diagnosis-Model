"""Command line interface for label-blind post-L5 factorial replay."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from .engine import OUTPUT_SCHEMA, PIPELINE_VERSION, load_config, run_patient
from .source_adapter import (
    CLINICAL_SCHEMA,
    PRUNING_SCHEMA,
    RAW_SCHEMA,
    SourceValidationError,
    canonical_json_sha256,
    discover_patient_artifacts,
    load_and_join_patient,
    load_pruning_manifest_anchors,
    sha256_file,
)
from .states import ARM_ORDER, BASE_ARM, FAMILIES, PREDICATES

BATCH_SCHEMA = "rag_re_clinical_pruning_abc.batch_manifest.v1"


def _patient_sort_key(patient_id: str) -> tuple[int, int | str]:
    return (0, int(patient_id)) if patient_id.isdigit() else (1, patient_id.casefold())


def _ensure_output_is_separate(output: Path, inputs: list[Path]) -> None:
    resolved_output = output.resolve()
    for source in inputs:
        resolved_source = source.resolve()
        if resolved_output == resolved_source or resolved_source in resolved_output.parents:
            raise SourceValidationError(
                f"Output must not equal or be nested under an input directory: {output}"
            )
    if output.exists():
        if not output.is_dir():
            raise SourceValidationError(f"Output path exists and is not a directory: {output}")
        if any(output.iterdir()):
            raise SourceValidationError(f"Output directory must be absent or empty: {output}")


def _write_json(path: Path, payload: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _aggregate(outputs: list[dict[str, Any]]) -> dict[str, Any]:
    counts = Counter()
    l5_states = Counter()
    signal_counts = {family: {key: Counter() for key in ("A", "B", "C")} for family in FAMILIES}
    predicate_counts = {
        family: {predicate: Counter() for predicate in PREDICATES} for family in FAMILIES
    }
    arm_counts = {
        arm: Counter(
            {
                "auto_positive_count": 0,
                "review_inclusive_count": 0,
                "manual_review_count": 0,
                "pruned_count": 0,
                "rescued_count": 0,
                "no_rescue_count": 0,
                "changed_from_l5_count": 0,
                "added_non_e0_auto_count": 0,
                "restored_e0_prune_auto_count": 0,
                "restored_e0_prune_review_count": 0,
                "promoted_e0_manual_auto_count": 0,
            }
        )
        for arm in ARM_ORDER
    }
    for output in outputs:
        counts.update(output["counts"])
        for candidate in output["candidates"]:
            l5_states[candidate["l5_state"]] += 1
            for family in FAMILIES:
                for key in ("A", "B", "C"):
                    signal_counts[family][key][_tri_label(candidate["signals"][family][key])] += 1
                for predicate in PREDICATES:
                    predicate_counts[family][predicate][
                        _tri_label(candidate["predicates"][family][predicate])
                    ] += 1
        for arm in ARM_ORDER:
            summary = output["arms"][arm]
            target = arm_counts[arm]
            for key in (
                "auto_positive_count",
                "review_inclusive_count",
                "manual_review_count",
                "pruned_count",
                "rescued_count",
                "no_rescue_count",
            ):
                target[key] += summary[key]
            for key, value in summary["transition_summary"].items():
                if key.endswith("_count"):
                    target[key] += value
    return {
        "counts": dict(counts),
        "l5_state_counts": dict(sorted(l5_states.items())),
        "signal_counts": {
            family: {key: dict(counter) for key, counter in values.items()}
            for family, values in signal_counts.items()
        },
        "predicate_counts": {
            family: {key: dict(counter) for key, counter in values.items()}
            for family, values in predicate_counts.items()
        },
        "arm_counts": {arm: dict(values) for arm, values in arm_counts.items()},
    }


def _tri_label(value: bool | None) -> str:
    return "true" if value is True else "false" if value is False else "unknown"


def replay_batch(args: argparse.Namespace) -> int:
    pruning_dir = Path(args.pruning_input)
    clinical_dir = Path(args.clinical_input)
    raw_dir = Path(args.raw_input)
    output_dir = Path(args.output)
    _ensure_output_is_separate(output_dir, [pruning_dir, clinical_dir, raw_dir])
    config, config_meta = load_config(args.config)
    pruning_index = discover_patient_artifacts(pruning_dir, PRUNING_SCHEMA)
    clinical_index = discover_patient_artifacts(clinical_dir, CLINICAL_SCHEMA)
    raw_index = discover_patient_artifacts(raw_dir, RAW_SCHEMA)
    sets = (set(pruning_index), set(clinical_index), set(raw_index))
    if not sets[0] == sets[1] == sets[2]:
        raise SourceValidationError(
            "Patient sets differ across pruning, clinical, and raw inputs: "
            f"pruning-only={sorted(sets[0] - sets[1] - sets[2])}, "
            f"clinical-only={sorted(sets[1] - sets[0] - sets[2])}, "
            f"raw-only={sorted(sets[2] - sets[0] - sets[1])}"
        )
    pruning_anchors, pruning_manifest_meta = load_pruning_manifest_anchors(
        pruning_dir, pruning_index, clinical_index
    )
    # Validate and decide the complete batch before writing the first artifact. A
    # late patient-level hash/schema error therefore cannot leave a plausible
    # looking half-run in the requested output directory.
    decided: list[tuple[str, dict[str, Any]]] = []
    for patient_id in sorted(sets[0], key=_patient_sort_key):
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", patient_id):
            raise SourceValidationError(f"Unsafe patient ID for output filename: {patient_id!r}")
        joined = load_and_join_patient(
            pruning_index[patient_id],
            clinical_index[patient_id],
            raw_index[patient_id],
            anchor=pruning_anchors[patient_id],
        )
        decided.append((patient_id, run_patient(joined, config, config_meta=config_meta)))

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[dict[str, Any]] = []
    artifacts: list[dict[str, Any]] = []
    for patient_id, output in decided:
        filename = f"NGS_patient_{patient_id}_post_l5_abc_v1.json"
        path = output_dir / filename
        _write_json(path, output)
        outputs.append(output)
        artifacts.append(
            {
                "patient_id": patient_id,
                "path": str(path.resolve()),
                "sha256": sha256_file(path),
                "decision_sha256": output["decision_sha256"],
                "candidate_count": output["counts"]["candidate_count"],
                "e0_count": output["counts"]["e0_count"],
                "source_hashes": {
                    key: value["sha256"] for key, value in output["sources"].items()
                },
            }
        )
    aggregate = _aggregate(outputs)
    manifest: dict[str, Any] = {
        "schema_version": BATCH_SCHEMA,
        "pipeline_version": PIPELINE_VERSION,
        "output_schema_version": OUTPUT_SCHEMA,
        "run_status": "complete",
        "study_status": config["study_status"],
        "primary_arm": BASE_ARM,
        "new_arm_primary_eligibility": False,
        "arm_order": list(ARM_ORDER),
        "sources": {
            "pruning_directory": str(pruning_dir.resolve()),
            "clinical_directory": str(clinical_dir.resolve()),
            "raw_directory": str(raw_dir.resolve()),
            "locked_read_only": True,
            "pruning_manifest": pruning_manifest_meta,
        },
        "config": config_meta,
        "generation_contract": {
            "reads_gold": False,
            "new_api_or_llm_calls": False,
            "exact_patient_and_candidate_join": True,
        },
        "patient_count": len(outputs),
        **aggregate,
        "artifacts": artifacts,
    }
    manifest["decision_set_sha256"] = canonical_json_sha256(
        [artifact["decision_sha256"] for artifact in artifacts]
    )
    _write_json(output_dir / "manifest.json", manifest)
    print(
        json.dumps(
            {
                "status": "complete",
                "output": str(output_dir.resolve()),
                "patient_count": manifest["patient_count"],
                "candidate_count": manifest["counts"]["candidate_count"],
                "arm_count": len(ARM_ORDER),
                "primary_arm": BASE_ARM,
                "decision_set_sha256": manifest["decision_set_sha256"],
            },
            ensure_ascii=False,
        )
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Replay frozen L5 with a label-blind corrected/raw A/B/C factorial matrix."
    )
    parser.add_argument("--pruning-input", required=True, help="Frozen pruning output directory")
    parser.add_argument("--clinical-input", required=True, help="Frozen clinical output directory")
    parser.add_argument("--raw-input", required=True, help="Frozen raw RAG_re output directory")
    parser.add_argument("--output", required=True, help="New, empty output directory")
    parser.add_argument(
        "--config",
        default=str(project_root / "config" / "post_l5_abc_policy_v1.json"),
        help="Frozen post-L5 policy JSON",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return replay_batch(args)
    except SourceValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

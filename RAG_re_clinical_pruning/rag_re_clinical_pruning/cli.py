from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Sequence

from .engine import (
    BATCH_MANIFEST_SCHEMA_VERSION,
    OUTPUT_SCHEMA_VERSION,
    PIPELINE_VERSION,
    PRIMARY_ARM,
    load_policy_config,
    run_patient,
)
from .source_adapter import index_source_directory, sha256_file


JSONDict = dict[str, Any]
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "pruning_policy_v1.json"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Replay the locked rag_re_clinical.output.v1 artifacts through the "
            "independent clinical pruning policy."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Locked clinical artifact JSON or directory.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="New/empty output directory; the source directory is never modified.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"Frozen pruning policy (default: {DEFAULT_CONFIG}).",
    )
    parser.add_argument(
        "--quiet", action="store_true", help="Suppress the success summary."
    )
    return parser


def _ensure_separate_output(source: Path, output: Path) -> None:
    source_resolved = source.resolve()
    output_resolved = output.resolve()
    protected_dir = source_resolved if source_resolved.is_dir() else source_resolved.parent
    if output_resolved == protected_dir or protected_dir in output_resolved.parents:
        raise ValueError(
            "Output must not be the locked source directory or one of its descendants"
        )
    if output_resolved.exists() and not output_resolved.is_dir():
        raise ValueError(f"Output exists and is not a directory: {output_resolved}")
    if output_resolved.exists() and any(output_resolved.glob("*.json")):
        raise ValueError(
            f"Output already contains JSON artifacts; choose a new/empty directory: "
            f"{output_resolved}"
        )


def _safe_patient_token(patient_id: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", patient_id).strip("._-")
    return token or "unknown"


def _artifact_filename(patient_id: str) -> str:
    return f"NGS_patient_{_safe_patient_token(patient_id)}_clinical_pruning_v1.json"


def _write_json_atomic(path: Path, value: Any) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    temporary.replace(path)


def run_batch(source: Path, output: Path, config: Path) -> JSONDict:
    _ensure_separate_output(source, output)
    policy = load_policy_config(config)
    indexed = index_source_directory(source)

    generated: list[tuple[str, Path, JSONDict]] = []
    used_filenames: set[str] = set()
    for patient_id, (source_path, projected) in indexed.items():
        filename = _artifact_filename(patient_id)
        folded = filename.casefold()
        if folded in used_filenames:
            raise ValueError(f"Output filename collision for patient_id={patient_id!r}")
        used_filenames.add(folded)
        generated.append((patient_id, output.resolve() / filename, run_patient(projected, policy)))
        if sha256_file(source_path) != projected["source"]["sha256"]:
            raise ValueError(f"Locked source changed during replay: {source_path}")

    if sha256_file(config.resolve()) != policy.sha256:
        raise ValueError(f"Config changed during replay: {config.resolve()}")

    output.mkdir(parents=True, exist_ok=True)
    artifacts: list[JSONDict] = []
    for patient_id, artifact_path, artifact in generated:
        _write_json_atomic(artifact_path, artifact)
        artifacts.append(
            {
                "patient_id": patient_id,
                "path": str(artifact_path),
                "sha256": sha256_file(artifact_path),
                "decision_sha256": artifact["decision_sha256"],
                "source_path": artifact["source"]["path"],
                "source_sha256": artifact["source"]["sha256"],
                "candidate_count": artifact["counts"]["candidate_count"],
                "e0_count": artifact["counts"]["e0_count"],
                "non_e0_count": artifact["counts"]["non_e0_count"],
            }
        )

    manifest: JSONDict = {
        "schema_version": BATCH_MANIFEST_SCHEMA_VERSION,
        "pipeline_version": PIPELINE_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "run_status": "complete",
        "source": {
            "path": str(source.resolve()),
            "schema_version": "rag_re_clinical.output.v1",
            "locked_read_only": True,
        },
        "config": {
            "path": policy.path,
            "sha256": policy.sha256,
            "schema_version": policy.raw["schema_version"],
            "policy_version": policy.raw["policy_version"],
        },
        "generation_contract": {
            "reads_gold": False,
            "configured_recall_safety_used_for_generation": False,
        },
        "primary_arm": PRIMARY_ARM,
        "patient_count": len(artifacts),
        "candidate_count": sum(row["candidate_count"] for row in artifacts),
        "e0_count": sum(row["e0_count"] for row in artifacts),
        "non_e0_count": sum(row["non_e0_count"] for row in artifacts),
        "artifacts": artifacts,
    }
    _write_json_atomic(output.resolve() / "manifest.json", manifest)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        manifest = run_batch(args.input, args.output, args.config)
    except (FileNotFoundError, OSError, ValueError, AssertionError) as exc:
        parser.exit(2, f"error: {exc}\n")
    if not args.quiet:
        print(
            "clinical pruning replay complete: "
            f"patients={manifest['patient_count']}, "
            f"candidates={manifest['candidate_count']}, "
            f"E0={manifest['e0_count']}, nonE0={manifest['non_e0_count']}, "
            f"primary={manifest['primary_arm']}"
        )
        print(f"manifest: {(args.output.resolve() / 'manifest.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

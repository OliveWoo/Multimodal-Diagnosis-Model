"""Run standardized patient JSON through the offline deterministic v20 pipeline.

The pipeline stages normalized patient files into a new output directory, builds
ranked mNGS candidates, creates a with-FilmArray deterministic clinical summary,
and runs the v20 deterministic scorer.  It never calls an LLM or the network.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from tools import build_deterministic_summary as summary_builder
from tools import build_ranked_mngs_from_grouped_folder as ranked_builder
from tools import deterministic_mngs_max_scorer as scorer
from tools import mngs_common as mngs


PIPELINE_VERSION = "standardized_json_to_deterministic_v20_v1"
FINAL_RESULTS_SCHEMA = "research.deterministic_v20_results.v1"
RUN_MANIFEST_SCHEMA = "research.deterministic_v20_run_manifest.v1"
PATIENT_DIR_RE = re.compile(r"^NGS_patient_(?P<id>\d+)(?:_json)?$")
PATIENT_FILE_RE = re.compile(r"^NGS_patient_(?P<id>\d+)_.*\.json$", re.IGNORECASE)
SOURCE_KINDS = ("auto", "mngs_grouped", "all_rk_ntc")


class PipelineInputError(ValueError):
    """Raised when standardized input cannot be mapped to one patient safely."""


@dataclass(frozen=True)
class PatientSource:
    patient_id: str
    files: tuple[Path, ...]
    agent_output_files: tuple[Path, ...]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def normalize_requested_patient(value: str) -> str:
    text = str(value).strip()
    match = re.search(r"(?:NGS_patient_)?(?P<id>\d+)(?:_json)?$", text)
    if not match:
        raise PipelineInputError(f"Invalid patient selector: {value}")
    return match.group("id")


def _direct_patient_files(directory: Path, patient_id: str) -> tuple[Path, ...]:
    prefix = f"NGS_patient_{patient_id}_"
    return tuple(
        sorted(
            path
            for path in directory.iterdir()
            if path.is_file()
            and path.suffix.lower() == ".json"
            and path.name.startswith(prefix)
        )
    )


def _agent_output_files(directory: Path) -> tuple[Path, ...]:
    agent_dir = directory / summary_builder.AGENT_OUTPUT_DIR_NAME
    if not agent_dir.is_dir():
        return ()
    return tuple(sorted(path for path in agent_dir.iterdir() if path.is_file() and path.suffix.lower() == ".json"))


def discover_patient_sources(
    input_root: Path,
    requested_patients: Sequence[str] | None = None,
) -> list[PatientSource]:
    input_root = input_root.resolve()
    if not input_root.is_dir():
        raise PipelineInputError(f"Input root is not a directory: {input_root}")

    requested = {
        normalize_requested_patient(value)
        for value in (requested_patients or [])
    }
    discovered: dict[str, PatientSource] = {}

    root_match = PATIENT_DIR_RE.fullmatch(input_root.name)
    patient_directories = [input_root] if root_match else [
        child
        for child in sorted(input_root.iterdir())
        if child.is_dir() and PATIENT_DIR_RE.fullmatch(child.name)
    ]
    for directory in patient_directories:
        match = PATIENT_DIR_RE.fullmatch(directory.name)
        assert match is not None
        patient_id = match.group("id")
        if requested and patient_id not in requested:
            continue
        files = _direct_patient_files(directory, patient_id)
        if not files:
            raise PipelineInputError(f"No standardized JSON files found in {directory}")
        discovered[patient_id] = PatientSource(
            patient_id=patient_id,
            files=files,
            agent_output_files=_agent_output_files(directory),
        )

    if not root_match:
        flat_groups: dict[str, list[Path]] = {}
        for path in sorted(input_root.iterdir()):
            if not path.is_file():
                continue
            match = PATIENT_FILE_RE.fullmatch(path.name)
            if match:
                flat_groups.setdefault(match.group("id"), []).append(path)
        for patient_id, files in flat_groups.items():
            if requested and patient_id not in requested:
                continue
            if patient_id in discovered:
                raise PipelineInputError(
                    f"Patient {patient_id} is present in both a patient directory and flat files"
                )
            discovered[patient_id] = PatientSource(
                patient_id=patient_id,
                files=tuple(files),
                agent_output_files=(),
            )

    if requested:
        missing = sorted(requested - set(discovered), key=int)
        if missing:
            raise PipelineInputError(f"Requested patients not found: {', '.join(missing)}")
    if not discovered:
        raise PipelineInputError(f"No NGS_patient_<id> standardized JSON inputs found in {input_root}")
    return [discovered[key] for key in sorted(discovered, key=int)]


def stage_patient_source(source: PatientSource, patients_root: Path) -> tuple[Path, list[dict[str, Any]]]:
    patient_dir = patients_root / f"NGS_patient_{source.patient_id}_json"
    patient_dir.mkdir(parents=True, exist_ok=False)
    agent_dir = patient_dir / summary_builder.AGENT_OUTPUT_DIR_NAME
    agent_dir.mkdir()

    input_rows: list[dict[str, Any]] = []
    for original in source.files:
        destination = patient_dir / original.name
        if destination.exists():
            raise PipelineInputError(f"Duplicate input filename for patient {source.patient_id}: {original.name}")
        shutil.copy2(original, destination)
        input_rows.append(
            {
                "source": str(original.resolve()),
                "staged": str(destination.relative_to(patients_root.parent).as_posix()),
                "sha256": sha256_file(original),
                "bytes": original.stat().st_size,
            }
        )
    for original in source.agent_output_files:
        destination = agent_dir / original.name
        if destination.exists():
            raise PipelineInputError(f"Duplicate agent output for patient {source.patient_id}: {original.name}")
        shutil.copy2(original, destination)
        input_rows.append(
            {
                "source": str(original.resolve()),
                "staged": str(destination.relative_to(patients_root.parent).as_posix()),
                "sha256": sha256_file(original),
                "bytes": original.stat().st_size,
            }
        )
    return patient_dir, input_rows


def _raw_mngs_source(patient_dir: Path, patient_id: str, source_kind: str) -> tuple[str, Path] | None:
    found = {
        kind: ranked_builder.find_source_file(patient_dir, patient_id, kind)
        for kind in ("mngs_grouped", "all_rk_ntc")
    }
    if source_kind != "auto":
        path = found[source_kind]
        return (source_kind, path) if path is not None else None
    available = [(kind, path) for kind, path in found.items() if path is not None]
    if len(available) > 1:
        raise PipelineInputError(
            f"Patient {patient_id} has both mNGS_grouped and all_RK_NTC inputs; "
            "select --source-kind explicitly"
        )
    return available[0] if available else None


def build_ranked_inputs(
    patient_dirs: Sequence[Path],
    *,
    source_kind: str,
    output_root: Path,
) -> tuple[Path, dict[str, str]]:
    patients: dict[str, dict[str, Any]] = {}
    source_modes: dict[str, str] = {}
    for patient_dir in patient_dirs:
        patient_id = mngs.normalize_patient_id(patient_dir)
        raw_source = _raw_mngs_source(patient_dir, patient_id, source_kind)
        local_ranked = patient_dir / f"NGS_patient_{patient_id}_mNGS_ranked_candidates.json"
        if raw_source is not None:
            selected_kind, source_path = raw_source
            payload = ranked_builder.build_ranked_payload(source_path, patient_id)
            if local_ranked.exists():
                local_ranked.unlink()
            ranked_builder.write_json(local_ranked, payload)
            source_modes[patient_id] = selected_kind
        else:
            payload = mngs.load_ranked_mngs_for_patient(None, patient_dir)
            if payload is None:
                expected = source_kind if source_kind != "auto" else "mNGS_grouped/all_RK_NTC/ranked"
                raise PipelineInputError(f"Patient {patient_id} has no usable {expected} mNGS input")
            source_modes[patient_id] = "existing_ranked"
        if not payload.get("records"):
            raise PipelineInputError(f"Patient {patient_id} produced no ranked mNGS records")
        patients[patient_id] = payload

    combined_path = output_root / "ranked_mngs.json"
    write_json(
        combined_path,
        {
            "patients": patients,
            "ranking_metadata": {
                "pipeline_version": PIPELINE_VERSION,
                "source_modes": source_modes,
            },
        },
    )
    return combined_path, source_modes


def run_pipeline(
    *,
    input_root: Path,
    output_root: Path,
    source_kind: str = "auto",
    requested_patients: Sequence[str] | None = None,
) -> dict[str, Any]:
    if source_kind not in SOURCE_KINDS:
        raise PipelineInputError(f"Unknown source kind: {source_kind}")
    input_root = input_root.resolve()
    output_root = output_root.resolve()
    if output_root.exists():
        raise FileExistsError(f"Output already exists; choose a new path: {output_root}")
    if output_root == input_root:
        raise PipelineInputError("Output root must differ from input root")

    patient_sources = discover_patient_sources(input_root, requested_patients)
    output_root.mkdir(parents=True, exist_ok=False)
    patients_root = output_root / "patients"
    patients_root.mkdir()

    patient_dirs: list[Path] = []
    input_manifest: dict[str, list[dict[str, Any]]] = {}
    try:
        for source in patient_sources:
            patient_dir, input_rows = stage_patient_source(source, patients_root)
            patient_dirs.append(patient_dir)
            input_manifest[source.patient_id] = input_rows

        ranked_path, source_modes = build_ranked_inputs(
            patient_dirs,
            source_kind=source_kind,
            output_root=output_root,
        )

        result_rows: dict[str, dict[str, Any]] = {}
        output_manifest: dict[str, dict[str, str]] = {}
        for patient_dir in patient_dirs:
            patient_id = mngs.normalize_patient_id(patient_dir)
            summary_path = summary_builder.process_patient(
                patient_dir,
                output_suffix=summary_builder.DEFAULT_OUTPUT_SUFFIX,
                overwrite=False,
                skip_existing=False,
            )
            assert summary_path is not None
            scorer_path = scorer.process_directory(
                patient_dir,
                ranked_mngs_path=ranked_path,
                summary_mode="deterministic",
                summary_suffix=None,
                output_suffix=scorer.DEFAULT_MODELLESS_OUTPUT_SUFFIX,
                overwrite=False,
            )
            if scorer_path is None:
                raise RuntimeError(f"v20 scorer did not produce output for patient {patient_id}")
            payload = scorer.load_json(scorer_path)
            best = payload.get("best_available_summary") or {}
            result_rows[patient_id] = {
                "picked_count": best.get("picked_count", 0),
                "picked_pathogens": best.get("picked_pathogens", []),
                "pathogen_candidates": payload.get("pathogen_candidates", []),
                "excluded_candidates": payload.get("excluded_candidates", []),
                "final_infection_likelihood": payload.get("final_infection_likelihood"),
                "dominant_source": payload.get("dominant_source"),
                "dominant_pathogen_type": payload.get("dominant_pathogen_type"),
                "rule_version": payload.get("rule_version"),
                "scorer_output": str(scorer_path.relative_to(output_root).as_posix()),
            }
            output_manifest[patient_id] = {
                "summary": str(summary_path.relative_to(output_root).as_posix()),
                "summary_sha256": sha256_file(summary_path),
                "scorer": str(scorer_path.relative_to(output_root).as_posix()),
                "scorer_sha256": sha256_file(scorer_path),
            }

        final_results_path = output_root / "final_results.json"
        write_json(
            final_results_path,
            {
                "schema": FINAL_RESULTS_SCHEMA,
                "pipeline_version": PIPELINE_VERSION,
                "patient_count": len(result_rows),
                "patients": result_rows,
            },
        )
        manifest = {
            "schema": RUN_MANIFEST_SCHEMA,
            "pipeline_version": PIPELINE_VERSION,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "mode": "offline_deterministic_with_filmarray_v20",
            "network_or_llm_used": False,
            "input_root": str(input_root),
            "output_root": str(output_root),
            "requested_source_kind": source_kind,
            "resolved_source_modes": source_modes,
            "patient_count": len(result_rows),
            "inputs": input_manifest,
            "outputs": output_manifest,
            "combined_ranked": {
                "path": str(ranked_path.relative_to(output_root).as_posix()),
                "sha256": sha256_file(ranked_path),
            },
            "final_results": {
                "path": final_results_path.name,
                "sha256": sha256_file(final_results_path),
            },
        }
        write_json(output_root / "run_manifest.json", manifest)
        return manifest
    except Exception as exc:
        failure_path = output_root / "run_failed.json"
        if not failure_path.exists():
            write_json(
                failure_path,
                {
                    "schema": "research.deterministic_v20_run_failure.v1",
                    "pipeline_version": PIPELINE_VERSION,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                },
            )
        raise


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run standardized patient JSON through ranked mNGS, with-FilmArray "
            "deterministic summary, and the deterministic v20 scorer."
        )
    )
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument(
        "--source-kind",
        choices=SOURCE_KINDS,
        default="auto",
        help="mNGS input kind. auto refuses ambiguous patients containing both source types.",
    )
    parser.add_argument(
        "--patients",
        nargs="*",
        help="Optional patient IDs, for example 3 10 or NGS_patient_3_json.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        manifest = run_pipeline(
            input_root=args.input_root,
            output_root=args.output_root,
            source_kind=args.source_kind,
            requested_patients=args.patients,
        )
    except (OSError, PipelineInputError, RuntimeError, TypeError, KeyError) as exc:
        print(f"Pipeline failed: {exc}", file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": "complete",
                "patient_count": manifest["patient_count"],
                "output_root": manifest["output_root"],
                "final_results": manifest["final_results"]["path"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

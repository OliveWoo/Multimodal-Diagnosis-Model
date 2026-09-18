"""Run the project flowchart from standardized patient JSON to audited delivery.

This module is deliberately an orchestrator.  It does not select pathogens or
reimplement any research rule; it validates handoffs, calls the existing stage
programs in their required order, and records the exact commands and outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parent
UPSTREAM_ROOT = REPO_ROOT / "upstream"
if str(UPSTREAM_ROOT) not in sys.path:
    sys.path.insert(0, str(UPSTREAM_ROOT))

from tools.run_deterministic_v20_pipeline import (  # noqa: E402
    PipelineInputError,
    discover_patient_sources,
)


CONFIG_SCHEMA = "research.full_pipeline_config.v1"
PLAN_SCHEMA = "research.full_pipeline_plan.v1"
MANIFEST_SCHEMA = "research.full_pipeline_run_manifest.v1"
MERGED_SUFFIX = "mNGS_max_deterministic_with_missed_candidate_review"
SAFE_COHORT_RE = re.compile(r"^[A-Za-z0-9_-]+$")


class FullPipelineError(ValueError):
    """Raised when the workflow cannot be wired without an ambiguous handoff."""


@dataclass(frozen=True)
class CohortSpec:
    name: str
    input_root: Path
    hospital: str
    dataset: str
    source_kind: str
    patient_ids: tuple[str, ...]


@dataclass(frozen=True)
class Stage:
    name: str
    program: str
    command: tuple[str, ...] | None
    cwd: Path
    inputs: tuple[Path, ...]
    outputs: tuple[Path, ...]
    cohort: str | None = None
    external: bool = False
    action: str | None = None
    action_data: Mapping[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "cohort": self.cohort,
            "program": self.program,
            "execution": "external_provider" if self.external else "local",
            "command": list(self.command) if self.command is not None else None,
            "cwd": str(self.cwd),
            "inputs": [str(path) for path in self.inputs],
            "outputs": [str(path) for path in self.outputs],
            "internal_handoff_action": self.action,
        }


@dataclass(frozen=True)
class WorkflowPlan:
    config_path: Path
    output_root: Path
    cohorts: tuple[CohortSpec, ...]
    stages: tuple[Stage, ...]
    env_file: Path | None

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": PLAN_SCHEMA,
            "config": str(self.config_path),
            "output_root": str(self.output_root),
            "cohorts": [
                {
                    "name": item.name,
                    "input_root": str(item.input_root),
                    "hospital": item.hospital,
                    "dataset": item.dataset,
                    "source_kind": item.source_kind,
                    "patient_ids": list(item.patient_ids),
                }
                for item in self.cohorts
            ],
            "stage_count": len(self.stages),
            "external_stage_count": sum(stage.external for stage in self.stages),
            "stages": [stage.as_dict() for stage in self.stages],
        }


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FullPipelineError(f"Could not read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise FullPipelineError(f"Expected a JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def config_path(base: Path, value: Any, *, default: Path | None = None) -> Path:
    if value in (None, ""):
        if default is None:
            raise FullPipelineError("A required path is absent from the workflow config")
        return default.resolve()
    path = Path(str(value)).expanduser()
    return (base / path).resolve() if not path.is_absolute() else path.resolve()


def nonempty_text(value: Any, label: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise FullPipelineError(f"{label} must not be empty")
    return text


def require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise FullPipelineError(f"{label} must be an object")
    return value


def discover_cohorts(config: Mapping[str, Any], base: Path) -> tuple[CohortSpec, ...]:
    raw_cohorts = require_mapping(config.get("cohorts"), "cohorts")
    if not raw_cohorts:
        raise FullPipelineError("cohorts must contain at least one cohort")
    cohorts: list[CohortSpec] = []
    for name, raw_value in raw_cohorts.items():
        if not SAFE_COHORT_RE.fullmatch(str(name)):
            raise FullPipelineError(f"Unsafe cohort name: {name!r}")
        raw = require_mapping(raw_value, f"cohorts.{name}")
        input_root = config_path(base, raw.get("input_root"))
        requested = [str(value) for value in (raw.get("patients") or [])]
        try:
            sources = discover_patient_sources(input_root, requested)
        except (OSError, PipelineInputError) as exc:
            raise FullPipelineError(f"Invalid input for cohort {name}: {exc}") from exc
        source_kind = str(raw.get("source_kind") or "auto")
        if source_kind not in {"auto", "mngs_grouped", "all_rk_ntc"}:
            raise FullPipelineError(f"Invalid source_kind for cohort {name}: {source_kind}")
        cohorts.append(
            CohortSpec(
                name=str(name),
                input_root=input_root,
                hospital=nonempty_text(raw.get("hospital"), f"cohorts.{name}.hospital"),
                dataset=nonempty_text(raw.get("dataset"), f"cohorts.{name}.dataset"),
                source_kind=source_kind,
                patient_ids=tuple(source.patient_id for source in sources),
            )
        )
    return tuple(cohorts)


def add_optional(command: list[str], flag: str, value: Any) -> None:
    if value not in (None, ""):
        command.extend([flag, str(value)])


def build_plan(config_file: Path, output_root: Path) -> WorkflowPlan:
    config_file = config_file.resolve()
    config = load_json(config_file)
    if config.get("schema_version") != CONFIG_SCHEMA:
        raise FullPipelineError(
            f"schema_version must be {CONFIG_SCHEMA!r}, got {config.get('schema_version')!r}"
        )
    base = config_file.parent
    output_root = output_root.resolve()
    cohorts = discover_cohorts(config, base)
    python = sys.executable
    env_file = (
        config_path(base, config.get("env_file")) if config.get("env_file") else None
    )

    review = require_mapping(config.get("missed_review") or {}, "missed_review")
    rag = require_mapping(config.get("rag_re") or {}, "rag_re")
    casefit = require_mapping(config.get("casefit") or {}, "casefit")
    evidence = require_mapping(config.get("evidence") or {}, "evidence")
    rationale = require_mapping(config.get("rationale") or {}, "rationale")

    rag_config = config_path(
        base,
        rag.get("config"),
        default=REPO_ROOT / "RAG_re" / "config" / "formal_merge_full_a.json",
    )
    casefit_config = config_path(
        base,
        casefit.get("config"),
        default=REPO_ROOT / "RAG_re_casefit" / "config" / "default.json",
    )
    clinical_config = config_path(
        base,
        config.get("clinical_config"),
        default=REPO_ROOT / "RAG_re_clinical" / "config" / "clinical_v1.json",
    )
    prompt = config_path(
        base,
        rationale.get("prompt"),
        default=(
            REPO_ROOT
            / "mngs_candidate_evidence_pipeline_20260902"
            / "mngs_candidate_evidence_pipeline"
            / "prompts"
            / "llm_reason_generator_zh_v1_citation_scope.txt"
        ),
    )

    stages: list[Stage] = []
    rationale_runs_root = output_root / "07_initial_rationale_runs"
    initial_audit_root = output_root / "08_initial_rationale_audit"
    runtime_root = output_root / "_runtime"
    initial_delivery_config = runtime_root / "delivery_initial.json"
    final_decisions_root = output_root / "09_r5_final_decisions"
    final_rationale_file = output_root / "11_r5_accepted_rationales.json"
    final_delivery_config = runtime_root / "delivery_final.json"
    delivery_root = output_root / "12_delivery"

    cohort_paths: dict[str, dict[str, Path]] = {}
    for cohort in cohorts:
        root = output_root / "cohorts" / cohort.name
        deterministic = root / "01_deterministic_v20"
        patient_root = deterministic / "patients"
        ranked = deterministic / "ranked_mngs.json"
        merged = root / "02_merged"
        selected = root / "03_selected_pathogens"
        evidence_run = root / "run"
        reasoning_chains = evidence_run / "rationale_views" / "traceable_reasoning_chains"
        rag_output = root / "05_rag_re"
        casefit_output = root / "06_casefit"
        rationale_run = rationale_runs_root / cohort.name
        r5_inputs_batch = root / "10_r5_rationale_inputs"
        r5_inputs = r5_inputs_batch / cohort.name
        r5_rationale_run = output_root / "10_r5_rationale_runs" / cohort.name
        cohort_paths[cohort.name] = {
            "root": root,
            "deterministic": deterministic,
            "patient_root": patient_root,
            "ranked": ranked,
            "merged": merged,
            "selected": selected,
            "evidence_run": evidence_run,
            "reasoning_chains": reasoning_chains,
            "rag": rag_output,
            "casefit": casefit_output,
            "rationale_run": rationale_run,
            "r5_inputs_batch": r5_inputs_batch,
            "r5_inputs": r5_inputs,
            "r5_rationale_run": r5_rationale_run,
        }

        deterministic_command = [
            python,
            "-B",
            str(REPO_ROOT / "run_deterministic_v20.py"),
            "--input-root",
            str(cohort.input_root),
            "--output-root",
            str(deterministic),
            "--source-kind",
            cohort.source_kind,
        ]
        if (require_mapping(config["cohorts"][cohort.name], f"cohorts.{cohort.name}")).get("patients"):
            deterministic_command.extend(["--patients", *cohort.patient_ids])
        stages.append(
            Stage(
                name=f"{cohort.name}:deterministic_v20",
                cohort=cohort.name,
                program="run_deterministic_v20.py",
                command=tuple(deterministic_command),
                cwd=REPO_ROOT,
                inputs=(cohort.input_root,),
                outputs=(deterministic / "run_manifest.json", ranked, patient_root),
            )
        )

        queue_csv = root / "review_queue_summary.csv"
        stages.append(
            Stage(
                name=f"{cohort.name}:missed_candidate_queue",
                cohort=cohort.name,
                program="upstream.tools.build_missed_candidate_review_queue",
                command=(
                    python,
                    "-B",
                    "-m",
                    "tools.build_missed_candidate_review_queue",
                    str(patient_root),
                    "--ranked-mngs",
                    str(ranked),
                    "--summary-csv",
                    str(queue_csv),
                ),
                cwd=UPSTREAM_ROOT,
                inputs=(patient_root, ranked),
                outputs=(queue_csv,),
            )
        )

        review_command = [
            python,
            "-B",
            "-m",
            "tools.review_missed_mngs_candidates",
            str(patient_root),
            "--summary-mode",
            "deterministic",
            "--model",
            str(review.get("model") or "gpt-5"),
            "--max-workers",
            str(int(review.get("max_workers", 1))),
        ]
        add_optional(review_command, "--reasoning-effort", review.get("reasoning_effort", "medium"))
        add_optional(review_command, "--temperature", review.get("temperature"))
        if review.get("keep_prompt", True):
            review_command.append("--keep-prompt")
        if review.get("keep_raw", True):
            review_command.append("--keep-raw")
        review_outputs = tuple(
            patient_root
            / f"NGS_patient_{patient_id}_json"
            / "summary_outputs"
            / f"NGS_patient_{patient_id}_mNGS_missed_candidate_review.json"
            for patient_id in cohort.patient_ids
        )
        stages.append(
            Stage(
                name=f"{cohort.name}:missed_candidate_review",
                cohort=cohort.name,
                program="upstream.tools.review_missed_mngs_candidates",
                command=tuple(review_command),
                cwd=UPSTREAM_ROOT,
                inputs=(patient_root,),
                outputs=review_outputs,
                external=True,
            )
        )

        merge_csv = root / "merge_summary.csv"
        stages.append(
            Stage(
                name=f"{cohort.name}:merge_deterministic_and_review",
                cohort=cohort.name,
                program="upstream.tools.merge_deterministic_max_with_missed_review",
                command=(
                    python,
                    "-B",
                    "-m",
                    "tools.merge_deterministic_max_with_missed_review",
                    str(patient_root),
                    "--summary-csv",
                    str(merge_csv),
                ),
                cwd=UPSTREAM_ROOT,
                inputs=(patient_root, *review_outputs),
                outputs=(merge_csv,),
            )
        )
        stages.append(
            Stage(
                name=f"{cohort.name}:flatten_merged_handoff",
                cohort=cohort.name,
                program="internal handoff adapter (copy only)",
                command=None,
                cwd=REPO_ROOT,
                inputs=(patient_root,),
                outputs=(merged,),
                action="flatten_merged",
                action_data={
                    "patient_root": str(patient_root),
                    "output_root": str(merged),
                    "patient_ids": list(cohort.patient_ids),
                },
            )
        )

        stages.append(
            Stage(
                name=f"{cohort.name}:extract_selected_pathogens",
                cohort=cohort.name,
                program="extract_selected_pathogens.py",
                command=(
                    python,
                    "-B",
                    str(
                        REPO_ROOT
                        / "mngs_candidate_evidence_pipeline_20260902"
                        / "mngs_candidate_evidence_pipeline"
                        / "extract_selected_pathogens.py"
                    ),
                    "--input-root",
                    str(merged),
                    "--output-root",
                    str(selected),
                    "--hospital",
                    cohort.hospital,
                    "--dataset",
                    cohort.dataset,
                    "--pipeline-version",
                    "deterministic_v20_plus_missed_review",
                    "--pretty",
                ),
                cwd=REPO_ROOT,
                inputs=(merged,),
                outputs=(selected / "selected_pathogens_manifest.json",),
            )
        )

        evidence_command = [
            python,
            "-B",
            str(
                REPO_ROOT
                / "mngs_candidate_evidence_pipeline_20260902"
                / "mngs_candidate_evidence_pipeline"
                / "run_pipeline.py"
            ),
            "--input-root",
            str(patient_root),
            "--output-root",
            str(evidence_run),
            "--selected-root",
            str(selected),
            "--candidate-source",
            "upstream_all",
            "--candidate-window-days",
            str(evidence.get("candidate_window_days", 2)),
            "--window-mode",
            str(evidence.get("window_mode") or "date"),
            "--anchor-source",
            str(evidence.get("anchor_source") or "auto"),
            "--mode",
            "combined",
            "--rationale-views",
            "--pretty",
        ]
        stages.append(
            Stage(
                name=f"{cohort.name}:candidate_evidence",
                cohort=cohort.name,
                program="candidate evidence run_pipeline.py",
                command=tuple(evidence_command),
                cwd=REPO_ROOT,
                inputs=(patient_root, selected),
                outputs=(reasoning_chains,),
            )
        )

        rag_command = [
            python,
            "-B",
            str(REPO_ROOT / "RAG_re" / "run_rag_re.py"),
            "batch",
            "--input-dir",
            str(merged),
            "--output-dir",
            str(rag_output),
            "--cache-dir",
            str(root / "_cache" / "rag_re"),
            "--config",
            str(rag_config),
            "--fail-fast",
        ]
        if rag.get("skip_literature", False):
            rag_command.append("--skip-literature")
        stages.append(
            Stage(
                name=f"{cohort.name}:rag_re",
                cohort=cohort.name,
                program="RAG_re/run_rag_re.py batch",
                command=tuple(rag_command),
                cwd=REPO_ROOT / "RAG_re",
                inputs=(merged, rag_config),
                outputs=(rag_output / "batch_manifest.json",),
                external=not bool(rag.get("skip_literature", False)),
            )
        )

        stages.append(
            Stage(
                name=f"{cohort.name}:casefit",
                cohort=cohort.name,
                program="RAG_re_casefit/run_casefit_batch.py",
                command=(
                    python,
                    "-B",
                    str(REPO_ROOT / "RAG_re_casefit" / "run_casefit_batch.py"),
                    "--merge-dir",
                    str(merged),
                    "--frozen-rag-dir",
                    str(rag_output),
                    "--output-dir",
                    str(casefit_output),
                    "--cache-dir",
                    str(root / "_cache" / "casefit"),
                    "--config",
                    str(casefit_config),
                    "--rag-re-config",
                    str(rag_config),
                    "--rag-re-cache-dir",
                    str(root / "_cache" / "rag_re"),
                ),
                cwd=REPO_ROOT / "RAG_re_casefit",
                inputs=(merged, rag_output, casefit_config, rag_config),
                outputs=(casefit_output / "manifest.json",),
                external=True,
            )
        )

        rationale_command = [
            python,
            "-B",
            str(
                REPO_ROOT
                / "mngs_candidate_evidence_pipeline_20260902"
                / "mngs_candidate_evidence_pipeline"
                / "generate_llm_rationales.py"
            ),
            "--input-root",
            str(reasoning_chains),
            "--output-root",
            str(rationale_run),
            "--patients",
            ",".join(cohort.patient_ids),
            "--prompt",
            str(prompt),
            "--model",
            str(rationale.get("model") or "gpt-5.6-luna"),
            "--reasoning-effort",
            str(rationale.get("reasoning_effort") or "low"),
            "--max-output-tokens",
            str(int(rationale.get("max_output_tokens", 4000))),
            "--max-excluded-candidates",
            str(int(rationale.get("max_excluded_candidates", 5))),
            "--retries",
            str(int(rationale.get("retries", 1))),
            "--pretty",
        ]
        stages.append(
            Stage(
                name=f"{cohort.name}:initial_rationales",
                cohort=cohort.name,
                program="generate_llm_rationales.py",
                command=tuple(rationale_command),
                cwd=REPO_ROOT,
                inputs=(reasoning_chains, prompt),
                outputs=(rationale_run / "llm_rationale_manifest.json",),
                external=True,
            )
        )

    audit_command = [
        python,
        "-B",
        str(
            REPO_ROOT
            / "mngs_candidate_evidence_pipeline_20260902"
            / "mngs_candidate_evidence_pipeline"
            / "audit_llm_rationale_batch.py"
        ),
        "--run-root",
        str(rationale_runs_root),
        "--report-root",
        str(initial_audit_root),
        "--cohorts",
        *[cohort.name for cohort in cohorts],
    ]
    for cohort in cohorts:
        audit_command.extend(
            ["--input-dir", f"{cohort.name}={cohort_paths[cohort.name]['reasoning_chains']}"]
        )
    stages.append(
        Stage(
            name="audit_initial_rationales",
            program="audit_llm_rationale_batch.py",
            command=tuple(audit_command),
            cwd=REPO_ROOT,
            inputs=tuple(
                path
                for cohort in cohorts
                for path in (
                    cohort_paths[cohort.name]["rationale_run"],
                    cohort_paths[cohort.name]["reasoning_chains"],
                )
            ),
            outputs=(
                initial_audit_root / "batch_audit.json",
                initial_audit_root / "all_patient_rationales.json",
            ),
        )
    )
    stages.append(
        Stage(
            name="write_initial_delivery_config",
            program="internal handoff adapter (configuration only)",
            command=None,
            cwd=REPO_ROOT,
            inputs=(initial_audit_root / "all_patient_rationales.json",),
            outputs=(initial_delivery_config,),
            action="write_delivery_config",
            action_data={
                "path": str(initial_delivery_config),
                "accepted_rationales": str(initial_audit_root / "all_patient_rationales.json"),
                "final_decisions_root": str(final_decisions_root),
                "cohorts": {
                    cohort.name: {
                        "selected_root": str(cohort_paths[cohort.name]["selected"]),
                        "casefit_root": str(cohort_paths[cohort.name]["casefit"]),
                    }
                    for cohort in cohorts
                },
            },
        )
    )
    stages.append(
        Stage(
            name="r5_final_decisions",
            program="OBER_patient_delivery/generate_r5_final_decisions.py",
            command=(
                python,
                "-B",
                str(REPO_ROOT / "OBER_patient_delivery" / "generate_r5_final_decisions.py"),
                "--workspace",
                str(REPO_ROOT),
                "--delivery-config",
                str(initial_delivery_config),
                "--clinical-config",
                str(clinical_config),
                "--output",
                str(final_decisions_root),
            ),
            cwd=REPO_ROOT,
            inputs=(initial_delivery_config, clinical_config),
            outputs=(final_decisions_root / "manifest.json",),
        )
    )

    regenerated_args: list[str] = []
    for cohort in cohorts:
        paths = cohort_paths[cohort.name]
        stages.append(
            Stage(
                name=f"{cohort.name}:prepare_r5_rationale_inputs",
                cohort=cohort.name,
                program="OBER_patient_delivery/prepare_r5_rationale_inputs.py",
                command=(
                    python,
                    "-B",
                    str(REPO_ROOT / "OBER_patient_delivery" / "prepare_r5_rationale_inputs.py"),
                    "--config",
                    str(initial_delivery_config),
                    "--final-decisions",
                    str(final_decisions_root),
                    "--output",
                    str(paths["r5_inputs_batch"]),
                    "--cohorts",
                    cohort.name,
                    "--patients",
                    ",".join(cohort.patient_ids),
                ),
                cwd=REPO_ROOT,
                inputs=(initial_delivery_config, final_decisions_root, paths["reasoning_chains"]),
                outputs=(paths["r5_inputs_batch"] / "manifest.json", paths["r5_inputs"]),
            )
        )
        stages.append(
            Stage(
                name=f"{cohort.name}:r5_rationales",
                cohort=cohort.name,
                program="generate_llm_rationales.py",
                command=(
                    python,
                    "-B",
                    str(
                        REPO_ROOT
                        / "mngs_candidate_evidence_pipeline_20260902"
                        / "mngs_candidate_evidence_pipeline"
                        / "generate_llm_rationales.py"
                    ),
                    "--input-root",
                    str(paths["r5_inputs"]),
                    "--output-root",
                    str(paths["r5_rationale_run"]),
                    "--patients",
                    ",".join(cohort.patient_ids),
                    "--prompt",
                    str(prompt),
                    "--model",
                    str(rationale.get("model") or "gpt-5.6-luna"),
                    "--reasoning-effort",
                    str(rationale.get("reasoning_effort") or "low"),
                    "--max-output-tokens",
                    str(int(rationale.get("max_output_tokens", 4000))),
                    "--max-excluded-candidates",
                    str(int(rationale.get("max_excluded_candidates", 5))),
                    "--retries",
                    str(int(rationale.get("retries", 1))),
                    "--pretty",
                ),
                cwd=REPO_ROOT,
                inputs=(paths["r5_inputs"], prompt),
                outputs=(paths["r5_rationale_run"] / "llm_rationale_manifest.json",),
                external=True,
            )
        )
        regenerated_args.extend(
            ["--regenerated-cohort", f"{cohort.name}={paths['r5_rationale_run']}"]
        )

    stages.append(
        Stage(
            name="assemble_r5_rationales",
            program="OBER_patient_delivery/assemble_r5_rationales.py",
            command=(
                python,
                "-B",
                str(REPO_ROOT / "OBER_patient_delivery" / "assemble_r5_rationales.py"),
                "--config",
                str(initial_delivery_config),
                "--final-decisions",
                str(final_decisions_root),
                *regenerated_args,
                "--output",
                str(final_rationale_file),
            ),
            cwd=REPO_ROOT,
            inputs=(
                initial_delivery_config,
                final_decisions_root,
                *tuple(cohort_paths[cohort.name]["r5_rationale_run"] for cohort in cohorts),
            ),
            outputs=(final_rationale_file,),
        )
    )
    stages.append(
        Stage(
            name="write_final_delivery_config",
            program="internal handoff adapter (configuration only)",
            command=None,
            cwd=REPO_ROOT,
            inputs=(final_rationale_file, final_decisions_root),
            outputs=(final_delivery_config,),
            action="write_delivery_config",
            action_data={
                "path": str(final_delivery_config),
                "accepted_rationales": str(final_rationale_file),
                "final_decisions_root": str(final_decisions_root),
                "cohorts": {
                    cohort.name: {
                        "selected_root": str(cohort_paths[cohort.name]["selected"]),
                        "casefit_root": str(cohort_paths[cohort.name]["casefit"]),
                    }
                    for cohort in cohorts
                },
            },
        )
    )
    stages.append(
        Stage(
            name="export_delivery",
            program="OBER_patient_delivery/export_delivery.mjs",
            command=(
                "node",
                str(REPO_ROOT / "OBER_patient_delivery" / "export_delivery.mjs"),
                "--config",
                str(final_delivery_config),
                "--output",
                str(delivery_root),
            ),
            cwd=REPO_ROOT,
            inputs=(final_delivery_config, final_rationale_file, final_decisions_root),
            outputs=(delivery_root / "manifest.json", delivery_root / "final_selection_list.json"),
        )
    )
    stages.append(
        Stage(
            name="audit_delivery",
            program="OBER_patient_delivery/audit_delivery.mjs",
            command=(
                "node",
                str(REPO_ROOT / "OBER_patient_delivery" / "audit_delivery.mjs"),
                "--output",
                str(delivery_root),
            ),
            cwd=REPO_ROOT,
            inputs=(delivery_root / "manifest.json",),
            outputs=(delivery_root / "manifest.json",),
        )
    )

    referenced_files = [
        rag_config,
        casefit_config,
        clinical_config,
        prompt,
        *(
            [env_file]
            if env_file is not None
            else []
        ),
    ]
    missing = [path for path in referenced_files if not path.is_file()]
    if missing:
        raise FullPipelineError("Missing configured files: " + ", ".join(map(str, missing)))
    return WorkflowPlan(
        config_path=config_file,
        output_root=output_root,
        cohorts=cohorts,
        stages=tuple(stages),
        env_file=env_file,
    )


def load_env(path: Path | None) -> dict[str, str]:
    env = dict(os.environ)
    if path is None:
        return env
    for raw_line in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        if key:
            env.setdefault(key, value)
    return env


def ensure_new_directory(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Output already exists; choose a new path: {path}")
    path.mkdir(parents=True, exist_ok=False)


def flatten_merged(data: Mapping[str, Any]) -> None:
    patient_root = Path(str(data["patient_root"])).resolve()
    output_root = Path(str(data["output_root"])).resolve()
    patient_ids = [str(value) for value in data["patient_ids"]]
    if output_root.exists():
        raise FileExistsError(f"Merged handoff already exists: {output_root}")
    output_root.mkdir(parents=True, exist_ok=False)
    for patient_id in patient_ids:
        source = (
            patient_root
            / f"NGS_patient_{patient_id}_json"
            / "summary_outputs"
            / f"NGS_patient_{patient_id}_{MERGED_SUFFIX}.json"
        )
        if not source.is_file():
            raise FullPipelineError(f"Merged handoff source is missing: {source}")
        destination = output_root / source.name
        shutil.copy2(source, destination)
        if sha256_file(source) != sha256_file(destination):
            raise FullPipelineError(f"Merged handoff copy hash mismatch: {source}")


def write_delivery_config(data: Mapping[str, Any]) -> None:
    destination = Path(str(data["path"])).resolve()
    if destination.exists():
        raise FileExistsError(f"Delivery config already exists: {destination}")
    final_root = Path(str(data["final_decisions_root"])).resolve()
    cohorts = {
        name: {
            "selected_root": str(Path(str(values["selected_root"])).resolve()),
            "casefit_root": str(Path(str(values["casefit_root"])).resolve()),
            "final_decision_root": str((final_root / name).resolve()),
        }
        for name, values in require_mapping(data["cohorts"], "delivery cohorts").items()
    }
    write_json(
        destination,
        {
            "schema_version": "ober.delivery_config.v1",
            "batch_name": destination.stem,
            "decision_stage": "workflow_final",
            "accepted_rationales": str(Path(str(data["accepted_rationales"])).resolve()),
            "selection_stage_note": (
                "Generated by run_full_pipeline.py; computational final pending clinical review."
            ),
            "cohorts": cohorts,
        },
    )


def execute_internal(stage: Stage) -> None:
    if stage.action == "flatten_merged":
        flatten_merged(stage.action_data)
    elif stage.action == "write_delivery_config":
        write_delivery_config(stage.action_data)
    else:
        raise FullPipelineError(f"Unknown internal handoff action: {stage.action}")


def verify_stage_outputs(stage: Stage) -> None:
    missing = [path for path in stage.outputs if not path.exists()]
    if missing:
        raise FullPipelineError(
            f"Stage {stage.name} returned successfully but outputs are missing: "
            + ", ".join(map(str, missing))
        )


def execute_plan(plan: WorkflowPlan) -> dict[str, Any]:
    child_env = load_env(plan.env_file)
    child_env["PYTHONDONTWRITEBYTECODE"] = "1"
    if any(stage.external for stage in plan.stages) and not child_env.get("OPENAI_API_KEY"):
        raise FullPipelineError(
            "OPENAI_API_KEY is required by review/RAG/casefit/rationale stages; "
            "set it in the environment or config env_file"
        )
    if shutil.which("node") is None:
        raise FullPipelineError("Node.js is required for the delivery export and audit stages")
    ensure_new_directory(plan.output_root)

    manifest: dict[str, Any] = {
        "schema_version": MANIFEST_SCHEMA,
        "status": "running",
        "started_at_utc": utc_now(),
        "config": str(plan.config_path),
        "config_sha256": sha256_file(plan.config_path),
        "output_root": str(plan.output_root),
        "cohorts": [cohort.name for cohort in plan.cohorts],
        "stages": [],
    }
    running_path = plan.output_root / "RUNNING.json"
    write_json(running_path, manifest)
    try:
        for index, stage in enumerate(plan.stages, start=1):
            row: dict[str, Any] = {
                "index": index,
                **stage.as_dict(),
                "status": "running",
                "started_at_utc": utc_now(),
            }
            manifest["stages"].append(row)
            write_json(running_path, manifest)
            print(f"[{index}/{len(plan.stages)}] {stage.name}", flush=True)
            if stage.command is None:
                execute_internal(stage)
            else:
                subprocess.run(
                    list(stage.command),
                    cwd=stage.cwd,
                    env=child_env,
                    check=True,
                )
            verify_stage_outputs(stage)
            row["status"] = "complete"
            row["finished_at_utc"] = utc_now()
            row["output_files"] = [
                {
                    "path": str(path),
                    "sha256": sha256_file(path) if path.is_file() else None,
                }
                for path in stage.outputs
            ]
            write_json(running_path, manifest)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["failed_at_utc"] = utc_now()
        manifest["error_type"] = type(exc).__name__
        manifest["error"] = str(exc)
        if manifest["stages"] and manifest["stages"][-1]["status"] == "running":
            manifest["stages"][-1]["status"] = "failed"
        write_json(plan.output_root / "run_failed.json", manifest)
        raise

    manifest["status"] = "complete"
    manifest["finished_at_utc"] = utc_now()
    manifest["delivery_root"] = str(plan.output_root / "12_delivery")
    final_manifest = plan.output_root / "run_manifest.json"
    write_json(final_manifest, manifest)
    running_path.unlink()
    return manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Orchestrate the existing project programs from standardized patient JSON "
            "through deterministic v20, review, RAG/casefit, R5 rationale, and delivery audit."
        )
    )
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument(
        "--plan",
        action="store_true",
        help="Validate inputs and print the exact stage plan without creating outputs or calling providers.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        plan = build_plan(args.config, args.output_root)
        if args.plan:
            print(json.dumps(plan.as_dict(), ensure_ascii=False, indent=2))
            return 0
        manifest = execute_plan(plan)
    except (OSError, FullPipelineError, subprocess.CalledProcessError) as exc:
        print(f"Full pipeline failed: {exc}", file=sys.stderr)
        return 1
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "cohorts": manifest["cohorts"],
                "stage_count": len(manifest["stages"]),
                "delivery_root": manifest["delivery_root"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

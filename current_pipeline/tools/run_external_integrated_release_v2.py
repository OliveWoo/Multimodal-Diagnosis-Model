"""Run the v2.6-compatible deterministic pipeline on an answer-free cohort.

This entry point is intentionally separate from the fixed 33-patient KH
release runner.  It reuses already generated hospital summaries, applies the
current taxonomy overlays and deterministic policies, and never loads
benchmark answers.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from tools import apply_test_aware_history_route_shadow as history
from tools import build_compact_multi_assay_entry_shadow as compact
from tools import build_index_event_clinical_timeline_shadow as timeline
from tools import build_mngs_only_external_input as mngs_only_input
from tools import build_phenotype_event_packets as phenotype_packets
from tools import build_test_aware_candidate_phenotype_shadow as candidate_phenotype
from tools import build_test_aware_clinical_deterministic_shadow as clinical
from tools import build_test_aware_deterministic_shadow_scorer as scorer
from tools import build_test_aware_possible_pathogen_shadow as reporting
from tools import promote_test_aware_clinical_high_shadow_v9 as promotion
from tools import rehydrate_clinical_hospital_evidence_shadow as rehydrate
from tools import hydrate_scorer_evidence_predecision as predecision_hydrate
from tools.run_kh_integrated_decision_pipeline import build_integrated_results
from tools import run_kh_integrated_release_v2 as kh_release


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "outputs/runs/2026-09-30_51patients_integrated_release_v2_6_heldout_freeze_v1"
)
DEFAULT_INVENTORY = (
    ROOT / "outputs/runs/2026-09-30_51patients_multi_assay_inventory_no_zero_v1"
)
DEFAULT_PATIENT_ROOT = (
    ROOT / "outputs/extracted_workbooks_pipeline_input_v3_51patients/patients"
)
DEFAULT_PHENOTYPE_ROOT = ROOT / "outputs/phenotype_decision_evidence_v1"
DEFAULT_SUMMARY_SUFFIX = "final_summary_51patients_no_zero_20260930"
EXTERNAL_SUMMARY_SUFFIX = "release_input_external_51_v1"
DEFAULT_DEVELOPMENT_FREEZE = (
    ROOT
    / "outputs/runs/2026-09-30_KH_integrated_release_v2_6_taxonomy_manifest_complete"
    / "answer_blind_freeze.json"
)
SCHEMA_VERSION = "external_integrated_release_shadow.v1"


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def load_taxonomy_completion_manifest(path: Path | None) -> dict[str, Any]:
    """Require a complete run-scoped taxonomy review before scoring."""

    if path is None:
        raise ValueError(
            "A taxonomy completion manifest is required. Run the organism-only "
            "NCBI/LLM/validator workflow before the deterministic scorer."
        )
    if not path.is_file():
        raise FileNotFoundError(f"Missing taxonomy completion manifest: {path}")
    payload = read_json(path)
    if payload.get("schema_version") != "organism_taxonomy_completion_manifest_v1":
        raise ValueError("Unsupported taxonomy completion manifest schema")
    if payload.get("passed") is not True:
        raise ValueError("Taxonomy completion manifest did not pass")
    if payload.get("answer_blind") is not True or payload.get("ground_truth_loaded") is not False:
        raise ValueError("Taxonomy completion manifest is not answer-blind")
    if payload.get("unreviewed_keys"):
        raise ValueError("Taxonomy completion manifest still contains unreviewed organisms")

    staged_overlay = Path(str(payload.get("staged_overlay") or ""))
    if not staged_overlay.is_file():
        raise FileNotFoundError(f"Missing staged taxonomy overlay: {staged_overlay}")
    expected_hash = str(payload.get("staged_overlay_sha256") or "")
    if kh_release.sha256_file(staged_overlay) != expected_hash:
        raise ValueError("Staged taxonomy overlay hash does not match completion manifest")

    active_overlay = str(os.environ.get("MNGS_TAXONOMY_AUTO_OVERLAY") or "").strip()
    if not active_overlay:
        raise ValueError(
            "MNGS_TAXONOMY_AUTO_OVERLAY must point to the validated run-scoped overlay"
        )
    if Path(active_overlay).resolve() != staged_overlay.resolve():
        raise ValueError("Active taxonomy overlay differs from the validated staged overlay")
    return payload


def patient_id_from_name(value: str) -> int:
    match = re.search(r"NGS_patient_(\d+)", value)
    if not match:
        raise ValueError(f"Cannot identify patient number from {value!r}")
    return int(match.group(1))


def patient_ids_from_root(patient_root: Path) -> set[int]:
    return {
        patient_id_from_name(path.name)
        for path in patient_root.glob("NGS_patient_*_json")
        if path.is_dir()
    }


def cohort_partition(
    expected_ids: set[int], development_freeze: Path
) -> dict[str, Any]:
    payload = read_json(development_freeze)
    development_ids = {
        int(value)
        for value in (
            payload.get("expected_patient_ids") or payload.get("patient_ids") or []
        )
    }
    if not development_ids:
        raise ValueError(
            "Development freeze has neither expected_patient_ids nor patient_ids: "
            f"{development_freeze}"
        )
    overlap_ids = expected_ids & development_ids
    heldout_ids = expected_ids - development_ids
    missing_development_ids = development_ids - expected_ids
    if not heldout_ids:
        raise ValueError("No held-out patients remain after development-cohort exclusion")
    return {
        "development_freeze": str(development_freeze.resolve()),
        "development_freeze_sha256": kh_release.sha256_file(development_freeze),
        "development_patient_ids": sorted(development_ids),
        "development_overlap_patient_ids": sorted(overlap_ids),
        "heldout_patient_ids": sorted(heldout_ids),
        "missing_development_patient_ids": sorted(missing_development_ids),
        "development_overlap_count": len(overlap_ids),
        "heldout_patient_count": len(heldout_ids),
        "patient_sets_disjoint": not bool(overlap_ids & heldout_ids),
    }


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8-sig")
        return
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def freeze_heldout_outputs(
    output_root: Path,
    heldout_ids: set[int],
    decision_patient_ids: set[int] | None = None,
) -> dict[str, Any]:
    subset_root = output_root / "heldout_subset"
    integrated_source = output_root / "integrated_decisions.csv"
    report_source = (
        output_root / "10_reporting_v3_route_a" / "complete_report.csv"
    )
    integrated_rows = [
        row
        for row in _read_csv(integrated_source)
        if int(row["patient_id"]) in heldout_ids
    ]
    report_rows = [
        row
        for row in _read_csv(report_source)
        if int(row["patient_id"]) in heldout_ids
    ]
    integrated_path = subset_root / "heldout_integrated_decisions.csv"
    report_path = subset_root / "heldout_complete_report.csv"
    _write_csv(integrated_path, integrated_rows)
    _write_csv(report_path, report_rows)
    observed_ids = {int(row["patient_id"]) for row in integrated_rows}
    decision_ids = (
        set(decision_patient_ids) & heldout_ids
        if decision_patient_ids is not None
        else set(observed_ids)
    )
    tier_counts = Counter(
        row.get("final_reporting_tier") or "" for row in integrated_rows
    )
    summary = {
        "answer_blind": True,
        "heldout_patient_ids": sorted(heldout_ids),
        "heldout_patient_count": len(heldout_ids),
        "patients_with_decision_packets": sorted(decision_ids),
        "patients_with_integrated_candidates": sorted(observed_ids),
        "patients_with_no_integrated_candidates": sorted(heldout_ids - observed_ids),
        "decision_patient_coverage_complete": decision_ids == heldout_ids,
        "integrated_candidate_coverage_complete": observed_ids == heldout_ids,
        # A patient may validly have no reportable candidate.  Cohort
        # completeness therefore means that every patient reached the
        # deterministic decision stage, not that the pipeline was forced to
        # emit at least one candidate for every patient.
        "patient_coverage_complete": decision_ids == heldout_ids,
        "integrated_candidate_count": len(integrated_rows),
        "complete_report_candidate_count": len(report_rows),
        "final_tier_counts": dict(sorted(tier_counts.items())),
        "integrated_decisions": str(integrated_path.resolve()),
        "integrated_decisions_sha256": kh_release.sha256_file(integrated_path),
        "complete_report": str(report_path.resolve()),
        "complete_report_sha256": kh_release.sha256_file(report_path),
        "answers_loaded": False,
    }
    write_json(subset_root / "summary.json", summary)
    return summary


def refresh_heldout_metadata(output_root: Path) -> dict[str, Any]:
    """Repair derived held-out summaries without rerunning or changing decisions."""

    freeze_path = output_root / "answer_blind_freeze.json"
    manifest_path = output_root / "run_manifest.json"
    summary_path = output_root / "summary.json"
    freeze = read_json(freeze_path)
    partition = freeze["cohort_partition"]
    heldout = freeze_heldout_outputs(
        output_root,
        {int(value) for value in partition["heldout_patient_ids"]},
        {int(value) for value in freeze.get("decision_patient_ids") or []},
    )
    freeze["heldout_subset"] = heldout
    freeze["passed"] = bool(
        freeze.get("passed") and heldout["patient_coverage_complete"]
    )
    write_json(freeze_path, freeze)
    manifest = read_json(manifest_path)
    manifest["heldout_subset"] = heldout
    manifest["answer_blind_freeze_sha256"] = kh_release.sha256_file(freeze_path)
    write_json(manifest_path, manifest)
    summary = read_json(summary_path)
    summary["heldout_subset"] = heldout
    write_json(summary_path, summary)
    return heldout


def inventory_patient_ids(inventory_root: Path) -> set[int]:
    return {
        patient_id_from_name(path.name)
        for path in (inventory_root / "patients").glob(
            "NGS_patient_*_multi_assay_mngs_evidence.json"
        )
    }


def freeze_existing_summaries(
    patient_root: Path,
    inventory_root: Path,
    output_dir: Path,
    *,
    source_suffix: str,
) -> dict[str, Any]:
    patient_ids = patient_ids_from_root(patient_root)
    inventory_ids = inventory_patient_ids(inventory_root)
    if not patient_ids:
        raise ValueError(f"No patient directories found under {patient_root}")
    if patient_ids != inventory_ids:
        raise ValueError(
            "Patient/inventory mismatch: "
            f"patient_only={sorted(patient_ids - inventory_ids)}, "
            f"inventory_only={sorted(inventory_ids - patient_ids)}"
        )
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Summary freeze output already exists: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    files = []
    for patient_id in sorted(patient_ids):
        source = (
            patient_root
            / f"NGS_patient_{patient_id}_json"
            / "summary_outputs"
            / f"NGS_patient_{patient_id}_{source_suffix}.json"
        )
        if not source.is_file():
            raise FileNotFoundError(source)
        payload = read_json(source)
        if not isinstance(payload, dict):
            raise ValueError(f"Expected JSON object: {source}")
        destination = output_dir / (
            f"NGS_patient_{patient_id}_{EXTERNAL_SUMMARY_SUFFIX}.json"
        )
        shutil.copy2(source, destination)
        files.append({
            "patient_id": patient_id,
            "source_path": str(source.resolve()),
            "source_sha256": kh_release.sha256_file(source),
            "frozen_path": str(destination.resolve()),
            "frozen_sha256": kh_release.sha256_file(destination),
        })
    summary = {
        "schema_version": "external_source_summary_freeze.v1",
        "answer_blind": True,
        "patient_count": len(files),
        "source_patient_root": str(patient_root.resolve()),
        "source_summary_suffix": source_suffix,
        "frozen_summary_suffix": EXTERNAL_SUMMARY_SUFFIX,
        "patient_files": files,
        "constraints": [
            "Existing cohort summaries are copied without reinterpretation.",
            "No benchmark answers are loaded.",
            "This differs from the fixed KH host-only overlay freeze.",
        ],
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def build_provisional_linkage_manifest(
    patient_ids: set[int], phenotype_root: Path, output_path: Path
) -> dict[str, Any]:
    episode_policy = read_json(kh_release.EPISODE_IDENTITY_FILE)
    clusters = [
        sorted(int(value) for value in item.get("members") or [])
        for item in episode_policy.get("clusters") or []
    ]
    counts: Counter[str] = Counter()
    patients = []
    for patient_id in sorted(patient_ids):
        ledger = (
            phenotype_root
            / f"patient_{patient_id}_phenotype_decision_evidence_v1.json"
        )
        available = ledger.is_file()
        status = (
            "dataset_contract_same_number_provisional"
            if available
            else "phenotype_file_missing"
        )
        counts[status] += 1
        cluster = next(
            (members for members in clusters if patient_id in members), []
        )
        patients.append({
            "pipeline_patient_number": patient_id,
            "proposed_phenotype_patient_number": patient_id,
            "phenotype_file_available": available,
            "local_linkage_status": status,
            "patient_identity_assessment": status,
            "infection_episode_assessment": (
                "index_infection_episode_unresolved"
                if available else "phenotype_file_missing"
            ),
            "phenotype_duplicate_export_cluster": cluster,
            "production_link_verified": False,
        })
    payload = {
        "schema_version": "external_phenotype_linkage_manifest.v1",
        "scope": {
            "patient_count": len(patient_ids),
            "contains_direct_identifiers": False,
            "production_link_verified": False,
            "linkage_basis": "same-number dataset contract; source-holder verification absent",
        },
        "status_counts": dict(sorted(counts.items())),
        "patients": patients,
    }
    write_json(output_path, payload)
    return payload


def component(path: Path, role: str) -> dict[str, Any]:
    return {
        "role": role,
        "path": str(path.resolve()),
        "sha256": kh_release.sha256_file(path),
    }


def prepare_optional_input_roots(
    *,
    output_root: Path,
    inventory_root: Path,
    patient_root: Path,
    phenotype_root: Path,
    source_summary_suffix: str,
    mngs_only: bool,
) -> tuple[Path, Path, str, dict[str, Any] | None]:
    """Resolve clinical inputs, generating explicit unknowns for mNGS-only runs."""
    if not mngs_only:
        return patient_root, phenotype_root, source_summary_suffix, None

    generated_root = output_root / "00_mngs_only_input"
    generated_manifest = mngs_only_input.run(inventory_root, generated_root)
    return (
        generated_root,
        generated_root / "phenotype_unavailable",
        mngs_only_input.SUMMARY_SUFFIX,
        generated_manifest,
    )


def run(
    output_root: Path = DEFAULT_OUTPUT,
    *,
    inventory_root: Path = DEFAULT_INVENTORY,
    patient_root: Path = DEFAULT_PATIENT_ROOT,
    phenotype_root: Path = DEFAULT_PHENOTYPE_ROOT,
    source_summary_suffix: str = DEFAULT_SUMMARY_SUFFIX,
    upstream_extraction_model: str = "not_declared",
    development_freeze: Path = DEFAULT_DEVELOPMENT_FREEZE,
    compatibility_target: str = "KH integrated release v2.6",
    analytical_policy: Path = kh_release.ANALYTICAL_POLICY,
    promotion_policy: Path = kh_release.PROMOTION_POLICY,
    reporting_policy: Path = kh_release.REPORTING_POLICY,
    taxonomy_completion_manifest: Path | None = None,
    mngs_only: bool = False,
    predecision_evidence_graph: bool = False,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty release: {output_root}")
    taxonomy_completion = load_taxonomy_completion_manifest(
        taxonomy_completion_manifest
    )
    effective_taxonomy_overlay = Path(taxonomy_completion["staged_overlay"])
    (
        patient_root,
        phenotype_root,
        source_summary_suffix,
        generated_mngs_only_manifest,
    ) = prepare_optional_input_roots(
        output_root=output_root,
        inventory_root=inventory_root,
        patient_root=patient_root,
        phenotype_root=phenotype_root,
        source_summary_suffix=source_summary_suffix,
        mngs_only=mngs_only,
    )
    input_evidence_mode = (
        "mngs_only_optional_clinical_absent"
        if mngs_only
        else "mngs_plus_provided_optional_clinical"
    )
    required = [
        inventory_root / "inventory_summary.json",
        patient_root,
        phenotype_root,
        kh_release.COMPACT_POLICY,
        analytical_policy,
        kh_release.HISTORY_POLICY,
        kh_release.CLINICAL_POLICY,
        promotion_policy,
        reporting_policy,
        kh_release.IDENTITY_ALIAS_FILE,
        kh_release.TAXONOMY_BASE_FILE,
        kh_release.TAXONOMY_GAP_FILE,
        effective_taxonomy_overlay,
        kh_release.TAXONOMY_REVIEWED_FILE,
        kh_release.EPISODE_IDENTITY_FILE,
        development_freeze,
        taxonomy_completion_manifest,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required external-shadow inputs: {missing}")
    reporting_policy_payload = read_json(reporting_policy)
    output_root.mkdir(parents=True, exist_ok=True)
    stage = {
        name: output_root / name
        for name in (
            "00_source_summary_freeze",
            "01_compact_taxonomy",
            "02_analytical_scorer",
            "02b_predecision_evidence_graph",
            "03_clinical_timeline",
            "04_candidate_phenotype",
            "05_phenotype_event_packets",
            "06_history_v4c",
            "07_clinical_scorer",
            "08_evidence_rehydrated",
            "09_promotion_v8",
            "10_reporting_v3_route_a",
        )
    }
    expected_ids = patient_ids_from_root(patient_root)
    partition = cohort_partition(expected_ids, development_freeze)
    stage_summaries: dict[str, Any] = {}
    stage_summaries["source_summary_freeze"] = freeze_existing_summaries(
        patient_root,
        inventory_root,
        stage["00_source_summary_freeze"],
        source_suffix=source_summary_suffix,
    )
    linkage_path = output_root / "phenotype_linkage_manifest.json"
    linkage = build_provisional_linkage_manifest(
        expected_ids, phenotype_root, linkage_path
    )
    stage_summaries["compact_taxonomy"] = compact.build_shadow(
        inventory_root, stage["01_compact_taxonomy"], kh_release.COMPACT_POLICY
    )
    stage_summaries["analytical_scorer"] = scorer.build_shadow(
        stage["01_compact_taxonomy"],
        patient_root,
        stage["02_analytical_scorer"],
        summary_suffix=EXTERNAL_SUMMARY_SUFFIX,
        summary_root=stage["00_source_summary_freeze"],
        policy_path=analytical_policy,
    )
    decision_scorer_root = stage["02_analytical_scorer"]
    if predecision_evidence_graph:
        stage_summaries["predecision_evidence_graph"] = predecision_hydrate.run(
            stage["02_analytical_scorer"],
            stage["02b_predecision_evidence_graph"],
        )
        decision_scorer_root = stage["02b_predecision_evidence_graph"]
    stage_summaries["clinical_timeline"] = timeline.build_shadow(
        patient_root,
        stage["03_clinical_timeline"],
        multi_assay_root=inventory_root / "patients",
    )
    stage_summaries["candidate_phenotype"] = candidate_phenotype.run(
        stage["02_analytical_scorer"],
        phenotype_root,
        stage["04_candidate_phenotype"],
        linkage_mode="same_number_shadow",
    )
    stage_summaries["phenotype_event_packets"] = phenotype_packets.run(
        stage["03_clinical_timeline"],
        phenotype_root,
        linkage_path,
        stage["04_candidate_phenotype"],
        stage["05_phenotype_event_packets"],
    )
    stage_summaries["history_v4c"] = history.run(
        decision_scorer_root,
        patient_root,
        stage["06_history_v4c"],
        policy_path=kh_release.HISTORY_POLICY,
        phenotype_root=stage["04_candidate_phenotype"],
    )
    if predecision_evidence_graph:
        stage_summaries["evidence_consumption_ledger"] = (
            predecision_hydrate.sync_history_consumption_ledger(
                stage["02b_predecision_evidence_graph"],
                stage["06_history_v4c"],
            )
        )
    stage_summaries["clinical_scorer"] = clinical.run(
        decision_scorer_root,
        stage["03_clinical_timeline"],
        stage["05_phenotype_event_packets"],
        stage["07_clinical_scorer"],
        policy_path=kh_release.CLINICAL_POLICY,
        history_route_root=stage["06_history_v4c"],
    )
    stage_summaries["evidence_rehydrated"] = rehydrate.run(
        stage["07_clinical_scorer"],
        decision_scorer_root,
        stage["08_evidence_rehydrated"],
    )
    stage_summaries["promotion_v8"] = promotion.run(
        stage["08_evidence_rehydrated"],
        stage["09_promotion_v8"],
        policy_path=promotion_policy,
    )
    stage_summaries["reporting_v3_route_a"] = reporting.run(
        stage["09_promotion_v8"],
        stage["10_reporting_v3_route_a"],
        policy_path=reporting_policy,
    )
    stage_summaries["integrated_output"] = build_integrated_results(
        stage["10_reporting_v3_route_a"], output_root
    )
    decision_ids = {
        patient_id_from_name(path.name)
        for path in (stage["02_analytical_scorer"] / "patient_outputs").glob(
            "NGS_patient_*_test_aware_deterministic_shadow.json"
        )
    }
    heldout_summary = freeze_heldout_outputs(
        output_root,
        set(partition["heldout_patient_ids"]),
        decision_ids,
    )
    patient_count_assertions = {
        name: summary.get("patient_count") == len(expected_ids)
        for name, summary in stage_summaries.items()
        if isinstance(summary, dict) and "patient_count" in summary
    }
    episode_validation = kh_release.validate_episode_separation(decision_ids)
    freeze = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "answers_loaded": False,
        "expected_patient_ids": sorted(expected_ids),
        "decision_patient_ids": sorted(decision_ids),
        "patient_count": len(decision_ids),
        "cohort_partition": partition,
        "heldout_subset": heldout_summary,
        "patient_count_assertions": patient_count_assertions,
        "episode_separation": episode_validation,
        "decision_artifacts": {
            "integrated_decisions_sha256": kh_release.sha256_file(
                output_root / "integrated_decisions.csv"
            ),
            "complete_report_sha256": kh_release.sha256_file(
                stage["10_reporting_v3_route_a"] / "complete_report.csv"
            ),
            "reporting_summary_sha256": kh_release.sha256_file(
                stage["10_reporting_v3_route_a"] / "summary.json"
            ),
        },
        "passed": (
            decision_ids == expected_ids
            and all(patient_count_assertions.values())
            and episode_validation["passed"]
            and partition["patient_sets_disjoint"]
            and heldout_summary["patient_coverage_complete"]
        ),
    }
    write_json(output_root / "answer_blind_freeze.json", freeze)
    if not freeze["passed"]:
        raise RuntimeError("External answer-blind completeness assertions failed")

    components = [
        component(kh_release.COMPACT_POLICY, "compact_entry_policy"),
        component(analytical_policy, "analytical_scorer_policy"),
        component(kh_release.HISTORY_POLICY, "history_v4c_policy"),
        component(kh_release.CLINICAL_POLICY, "clinical_scorer_policy"),
        component(promotion_policy, "promotion_policy"),
        component(reporting_policy, "reporting_policy"),
        component(kh_release.IDENTITY_ALIAS_FILE, "production_identity_aliases"),
        component(kh_release.TAXONOMY_BASE_FILE, "central_taxonomy_base"),
        component(kh_release.TAXONOMY_GAP_FILE, "taxonomy_gap_resolution_overlay"),
        component(effective_taxonomy_overlay, "run_scoped_auto_taxonomy_overlay"),
        component(kh_release.TAXONOMY_REVIEWED_FILE, "reviewed_taxonomy_overlay"),
        component(taxonomy_completion_manifest, "taxonomy_completion_manifest"),
        component(kh_release.EPISODE_IDENTITY_FILE, "episode_separation_policy"),
    ]
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "release_id": output_root.name,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "answer_free_external_shadow",
        "compatibility_target": compatibility_target,
        "answer_blind_decision_stages": True,
        "answer_blind_freeze_sha256": kh_release.sha256_file(
            output_root / "answer_blind_freeze.json"
        ),
        "upstream_extraction_model": upstream_extraction_model,
        "decision_model": "deterministic_only",
        "input": {
            "evidence_mode": input_evidence_mode,
            "inventory_root": str(inventory_root.resolve()),
            "inventory_summary_sha256": kh_release.sha256_file(
                inventory_root / "inventory_summary.json"
            ),
            "patient_root": str(patient_root.resolve()),
            "source_summary_suffix": source_summary_suffix,
            "phenotype_root": str(phenotype_root.resolve()),
            "phenotype_linkage_manifest": str(linkage_path.resolve()),
            "phenotype_linkage_manifest_sha256": kh_release.sha256_file(
                linkage_path
            ),
            "analytical_policy": str(analytical_policy.resolve()),
            "analytical_policy_sha256": kh_release.sha256_file(analytical_policy),
            "promotion_policy": str(promotion_policy.resolve()),
            "promotion_policy_sha256": kh_release.sha256_file(promotion_policy),
            "reporting_policy": str(reporting_policy.resolve()),
            "reporting_policy_sha256": kh_release.sha256_file(reporting_policy),
            "predecision_evidence_graph": predecision_evidence_graph,
            "taxonomy_completion_manifest": str(
                taxonomy_completion_manifest.resolve()
            ),
            "taxonomy_completion_manifest_sha256": kh_release.sha256_file(
                taxonomy_completion_manifest
            ),
            "run_scoped_taxonomy_overlay": str(
                effective_taxonomy_overlay.resolve()
            ),
            "run_scoped_taxonomy_overlay_sha256": kh_release.sha256_file(
                effective_taxonomy_overlay
            ),
            "generated_mngs_only_input_manifest": (
                generated_mngs_only_manifest
                if generated_mngs_only_manifest is not None
                else None
            ),
        },
        "components": components,
        "stage_roots": {
            name: str(path.resolve())
            for name, path in stage.items()
            if path.exists()
        },
        "patient_count": len(decision_ids),
        "candidate_count": stage_summaries["integrated_output"][
            "candidate_count"
        ],
        "cohort_partition": partition,
        "heldout_subset": heldout_summary,
        "phenotype_linkage_status_counts": linkage["status_counts"],
        "episode_separation": episode_validation,
        "identity_aliases_active": True,
        "taxonomy_auto_overlay_active": True,
        "taxonomy_review_complete": True,
        "taxonomy_review_status_counts": taxonomy_completion.get(
            "status_counts", {}
        ),
        "route_a_possible_active": True,
        "route_b_fallback_active": True,
        "graceful_degradation_reporting_active": bool(
            (reporting_policy_payload.get("graceful_degradation") or {}).get(
                "enabled"
            )
        ),
        "predecision_evidence_graph_active": predecision_evidence_graph,
        "ober_included": False,
        "evaluation_attached_posthoc": False,
        "constraints": [
            "No benchmark answers were loaded or evaluated.",
            "Phenotype links are provisional and cannot independently create Picked.",
            "Missing phenotype files remain explicit missingness.",
            (
                "Clinical, imaging, phenotype, and other hospital evidence was not provided; absence is unknown, not negative."
                if mngs_only
                else "Existing external-cohort hospital summaries are frozen as-is."
            ),
            "OBER and Luna are not decision stages; the declared LLM model was upstream extraction only.",
        ],
    }
    write_json(output_root / "run_manifest.json", manifest)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "release_id": manifest["release_id"],
        "answer_blind_freeze_passed": freeze["passed"],
        "patient_count": len(decision_ids),
        "candidate_count": manifest["candidate_count"],
        "input_evidence_mode": input_evidence_mode,
        "predecision_evidence_graph_active": predecision_evidence_graph,
        "cohort_partition": partition,
        "heldout_subset": heldout_summary,
        "stage_summaries": stage_summaries,
        "evaluation": None,
        "ober_included": False,
    }
    write_json(output_root / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--inventory-root", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--phenotype-root", type=Path, default=DEFAULT_PHENOTYPE_ROOT)
    parser.add_argument("--source-summary-suffix", default=DEFAULT_SUMMARY_SUFFIX)
    parser.add_argument(
        "--mngs-only",
        action="store_true",
        help=(
            "Use the multi-assay inventory as the only required cohort data. "
            "Missing clinical and phenotype inputs are generated as explicit unknowns."
        ),
    )
    parser.add_argument("--upstream-extraction-model", default="not_declared")
    parser.add_argument(
        "--development-freeze", type=Path, default=DEFAULT_DEVELOPMENT_FREEZE
    )
    parser.add_argument(
        "--refresh-heldout-only",
        action="store_true",
        help="Refresh derived held-out CSV/summary metadata without rerunning decisions.",
    )
    parser.add_argument(
        "--compatibility-target", default="KH integrated release v2.6"
    )
    parser.add_argument(
        "--analytical-policy",
        type=Path,
        default=kh_release.ANALYTICAL_POLICY,
        help="Analytical scorer policy; defaults to the frozen production policy.",
    )
    parser.add_argument(
        "--promotion-policy",
        type=Path,
        default=kh_release.PROMOTION_POLICY,
        help="Clinical promotion policy; defaults to the frozen production policy.",
    )
    parser.add_argument(
        "--reporting-policy",
        type=Path,
        default=kh_release.REPORTING_POLICY,
        help="Reporting policy; defaults to the frozen production policy.",
    )
    parser.add_argument(
        "--predecision-evidence-graph",
        action="store_true",
        help=(
            "Hydrate canonical candidate evidence and initialize the evidence "
            "ledger before history and clinical decisions. This is an opt-in "
            "vNext shadow and does not alter the frozen default pipeline."
        ),
    )
    parser.add_argument(
        "--taxonomy-completion-manifest",
        type=Path,
        help=(
            "Passed organism_taxonomy_completion_manifest_v1. The staged overlay "
            "must also be active through MNGS_TAXONOMY_AUTO_OVERLAY."
        ),
    )
    args = parser.parse_args()
    if args.refresh_heldout_only:
        print(json.dumps(
            refresh_heldout_metadata(args.output_root),
            ensure_ascii=False,
            indent=2,
        ))
        return
    print(json.dumps(run(
        args.output_root,
        inventory_root=args.inventory_root,
        patient_root=args.patient_root,
        phenotype_root=args.phenotype_root,
        source_summary_suffix=args.source_summary_suffix,
        upstream_extraction_model=args.upstream_extraction_model,
        development_freeze=args.development_freeze,
        compatibility_target=args.compatibility_target,
        analytical_policy=args.analytical_policy,
        promotion_policy=args.promotion_policy,
        reporting_policy=args.reporting_policy,
        taxonomy_completion_manifest=args.taxonomy_completion_manifest,
        mngs_only=args.mngs_only,
        predecision_evidence_graph=args.predecision_evidence_graph,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

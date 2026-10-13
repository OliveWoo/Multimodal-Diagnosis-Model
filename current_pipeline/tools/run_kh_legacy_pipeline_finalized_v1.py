"""Run the frozen KH legacy-compatible pipeline release.

This release keeps the old selected-DNA input and v20-compatible scorer, then
applies the locked taxonomy v2.6, history v4C, HAP/airway-device Priority
overlay, and aspiration-convergence Picked gate.  Final decisions are frozen
before benchmark answers are loaded.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

from tools.apply_healthcare_device_opportunist_history_route import (
    run as run_healthcare_device_route,
)
from tools.apply_legacy_history_picked_promotion import run as run_picked_gate
from tools.apply_test_aware_history_route_shadow import run as run_history_route
from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.legacy_primary_selection import (
    build_primary_selections,
    complete_report_rows,
    doctor_report_rows,
)
from tools.legacy_taxonomy_snapshot import materialize_locked_overlays
from tools.run_kh_ablation_d_history_only import (
    answer_status,
    build_bridge,
    read_inventory,
    read_json,
    sha256_file,
    write_csv,
    write_json,
)
from tools.run_kh_ablation_e_taxonomy_only import score_mode
from tools.run_kh_ablation_h_legacy_taxonomy_history import load_taxonomy_profiles


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = REPO_ROOT / "rules/kh_legacy_pipeline_finalized_v1.json"
DEFAULT_OUTPUT = (
    REPO_ROOT / "outputs/runs/2026-10-07_KH_legacy_pipeline_finalized_v1"
)
SCHEMA_VERSION = "kh_legacy_pipeline_finalized.v1"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def resolve_repo_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def require_hash(path: Path, expected: str, label: str) -> str:
    actual = sha256_file(path)
    if actual.casefold() != expected.casefold():
        raise ValueError(
            f"Locked {label} hash mismatch: expected {expected}, got {actual}: {path}"
        )
    return actual


def group_labels(
    rows: list[dict[str, str]], allowed_labels: set[str]
) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = {}
    for row in rows:
        if row.get("final_label") not in allowed_labels:
            continue
        patient = int(row["patient_id"])
        organism = str(row.get("organism_name") or "")
        if organism and organism not in grouped.setdefault(patient, []):
            grouped[patient].append(organism)
    return grouped


def run(contract_path: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")

    contract = read_json(contract_path)
    if contract.get("release_id") != "kh_legacy_pipeline_finalized_v1":
        raise ValueError("Unexpected release contract")

    paths = contract["paths"]
    locked = contract["locked_sha256"]
    a_root = resolve_repo_path(paths["frozen_a_root"])
    patient_root = resolve_repo_path(paths["patient_root"])
    history_policy = resolve_repo_path(paths["history_policy"])
    healthcare_policy = resolve_repo_path(paths["healthcare_device_policy"])
    picked_policy = resolve_repo_path(paths["picked_policy"])
    answers_path = resolve_repo_path(paths["benchmark_answers"])

    baseline_contract_path = a_root / "baseline_contract.json"
    require_hash(
        baseline_contract_path, locked["frozen_a_contract"], "frozen A contract"
    )
    require_hash(history_policy, locked["history_policy"], "history policy")
    require_hash(
        healthcare_policy,
        locked["healthcare_device_policy"],
        "healthcare-device policy",
    )
    require_hash(picked_policy, locked["picked_policy"], "Picked policy")

    baseline_contract = read_json(baseline_contract_path)
    patient_ids = [int(value) for value in baseline_contract["patient_selection"]]
    inventory = read_inventory(a_root / baseline_contract["artifact_inventory_file"])
    output_root.mkdir(parents=True)
    locked_taxonomy_overlays = materialize_locked_overlays(
        repo_root=REPO_ROOT,
        output_root=output_root,
        contract=contract,
    )

    taxonomy_root = output_root / "01_legacy_scorer_taxonomy_v2_6"
    _, taxonomy_summary = score_mode(
        mode="current-overlays",
        patient_ids=patient_ids,
        inventory=inventory,
        patient_root=patient_root,
        output_dir=taxonomy_root,
        taxonomy_overlay_paths=locked_taxonomy_overlays,
    )
    if taxonomy_summary["source_hash_mismatches"]:
        raise ValueError("Frozen A source hash validation failed")

    taxonomy_snapshot = taxonomy_summary["taxonomy_snapshot"]
    if taxonomy_snapshot["base_rules_sha256"] != locked["taxonomy_base_rules"]:
        raise ValueError("Locked taxonomy base rules changed")
    overlay_hashes = {
        Path(path).name: digest
        for path, digest in taxonomy_snapshot["loaded_overlay_sha256"].items()
    }
    expected_overlays = {
        "organism_taxonomy_gap_resolution_v1.json": locked[
            "taxonomy_gap_resolution_overlay"
        ],
        "organism_taxonomy_auto_v1.json": locked["taxonomy_auto_overlay"],
        "organism_taxonomy_reviewed_v1.json": locked[
            "taxonomy_reviewed_overlay"
        ],
    }
    if overlay_hashes != expected_overlays:
        raise ValueError(
            f"Locked taxonomy overlays changed: {overlay_hashes} != {expected_overlays}"
        )

    taxonomy_profiles = load_taxonomy_profiles(taxonomy_root / "patient_outputs")
    bridge_root = output_root / "02_legacy_route_bridge"
    bridge_summary = build_bridge(
        a_root=a_root,
        patient_ids=patient_ids,
        inventory=inventory,
        output_dir=bridge_root,
        taxonomy_profiles=taxonomy_profiles,
        taxonomy_mode="locked_taxonomy_v2_6",
    )
    if bridge_summary["source_hash_mismatches"]:
        raise ValueError("Frozen A bridge source hash validation failed")

    history_root = output_root / "03_history_route_v4c"
    history_summary = run_history_route(
        bridge_root,
        patient_root,
        history_root,
        policy_path=history_policy,
        phenotype_root=None,
        patients=set(patient_ids),
    )
    if history_summary["resolved_policy_sha256"] != locked["resolved_history_policy"]:
        raise ValueError("Locked resolved history policy changed")

    healthcare_root = output_root / "04_healthcare_device_priority_route"
    healthcare_summary = run_healthcare_device_route(
        history_root, healthcare_root, healthcare_policy
    )
    picked_root = output_root / "05_aspiration_picked_gate"
    picked_summary = run_picked_gate(healthcare_root, picked_root, picked_policy)

    label_map = contract["output_labels"]
    final_rows = read_csv(picked_root / "promotion_decisions.csv")
    unknown_decisions = sorted(
        {row["final_decision"] for row in final_rows} - set(label_map)
    )
    if unknown_decisions:
        raise ValueError(f"Unmapped final decisions: {unknown_decisions}")
    for row in final_rows:
        row["final_label"] = label_map[row["final_decision"]]
    if any(row["final_label"] == "Possible" for row in final_rows):
        raise ValueError("The finalized legacy release must not emit Possible")

    final_root = output_root / "06_frozen_decisions"
    final_path = final_root / "final_decisions.csv"
    write_csv(final_path, final_rows)
    primary_rows = build_primary_selections(final_rows, patient_ids)
    primary_path = final_root / "patient_primary_output.csv"
    write_csv(primary_path, primary_rows)
    complete_fallback_rows = complete_report_rows(final_rows, primary_rows)
    doctor_report = doctor_report_rows(final_rows, primary_rows)
    doctor_report_path = final_root / "doctor_visible_report.csv"
    write_csv(doctor_report_path, doctor_report)
    label_counts = dict(sorted(Counter(row["final_label"] for row in final_rows).items()))
    observed = {
        "patient_count": len({int(row["patient_id"]) for row in final_rows}),
        "candidate_count": len(final_rows),
        "label_counts": label_counts,
    }
    if observed != contract["expected_development_output"]:
        raise ValueError(
            "Finalized release regression mismatch: "
            f"{observed} != {contract['expected_development_output']}"
        )
    primary_selection_counts = dict(
        sorted(
            Counter(
                row["primary_selection_role"] for row in primary_rows
            ).items()
        )
    )

    # This file is deliberately written before the benchmark is opened.
    decision_freeze = {
        "schema_version": SCHEMA_VERSION,
        "release_id": contract["release_id"],
        "status": "decisions_frozen_before_benchmark_load",
        "answers_loaded_during_decision_generation": False,
        "cohort_role": contract["cohort_role"],
        "independent_validation_claim": contract["independent_validation_claim"],
        "components": contract["components"],
        "constraints": contract["constraints"],
        "observed_output": observed,
        "artifacts": {
            "release_contract": str(contract_path.resolve()),
            "release_contract_sha256": sha256_file(contract_path),
            "taxonomy_summary_sha256": sha256_file(taxonomy_root / "summary.json"),
            "bridge_summary_sha256": sha256_file(bridge_root / "summary.json"),
            "history_summary_sha256": sha256_file(history_root / "summary.json"),
            "healthcare_summary_sha256": sha256_file(
                healthcare_root / "summary.json"
            ),
            "picked_summary_sha256": sha256_file(picked_root / "summary.json"),
            "final_decisions_sha256": sha256_file(final_path),
            "patient_primary_output_sha256": sha256_file(primary_path),
            "doctor_visible_report_sha256": sha256_file(doctor_report_path),
            "runner_sha256": sha256_file(Path(__file__)),
        },
    }
    freeze_path = final_root / "decision_freeze.json"
    write_json(freeze_path, decision_freeze)

    # Benchmark access is intentionally below the persisted decision freeze.
    answers = load_answer_rows(answers_path)
    metric_sets = {
        "strict_picked": {"Picked"},
        "picked_or_priority": {"Picked", "Priority"},
        "complete_visible": {"Picked", "Priority", "Context"},
    }
    metrics = {
        name: evaluate_set(group_labels(final_rows, labels), answers)
        for name, labels in metric_sets.items()
    }
    metrics["picked_plus_mandatory_fallback"] = evaluate_set(
        group_labels(
            complete_fallback_rows,
            {"Picked", "Priority", "Context", "Audit"},
        ),
        answers,
    )
    visible_rows: list[dict[str, Any]] = []
    for row in final_rows:
        if row["final_label"] == "Audit":
            continue
        visible_rows.append(
            {
                "patient_id": row["patient_id"],
                "organism_name": row["organism_name"],
                "final_label": row["final_label"],
                "posthoc_answer_status": answer_status(
                    int(row["patient_id"]), row["organism_name"], answers
                ),
            }
        )
    evaluation_root = output_root / "07_postfreeze_evaluation"
    write_csv(evaluation_root / "visible_candidate_status.csv", visible_rows)
    evaluation = {
        "schema_version": SCHEMA_VERSION,
        "scope": "Post-freeze benchmark evaluation; not independent validation.",
        "decision_freeze_sha256": sha256_file(freeze_path),
        "benchmark_path": str(answers_path.resolve()),
        "benchmark_sha256": sha256_file(answers_path),
        "metrics": metrics,
        "primary_selection_counts": primary_selection_counts,
    }
    write_json(evaluation_root / "metrics.json", evaluation)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "release_id": contract["release_id"],
        "status": "complete_frozen_development_release",
        "decision_generation_answer_blind": True,
        "benchmark_evaluation_is_postfreeze": True,
        "independent_validation_claim": False,
        "observed_output": observed,
        "stage_counts": {
            "taxonomy": taxonomy_summary["counts"],
            "history": history_summary["counts"],
            "healthcare_device_route": healthcare_summary["counts"],
            "picked_gate": picked_summary["counts"],
        },
        "metrics": metrics,
        "primary_selection_counts": primary_selection_counts,
        "decision_freeze": str(freeze_path.resolve()),
        "decision_freeze_sha256": sha256_file(freeze_path),
        "final_decisions": str(final_path.resolve()),
        "final_decisions_sha256": sha256_file(final_path),
        "patient_primary_output": str(primary_path.resolve()),
        "patient_primary_output_sha256": sha256_file(primary_path),
        "doctor_visible_report": str(doctor_report_path.resolve()),
        "doctor_visible_report_sha256": sha256_file(doctor_report_path),
    }
    write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run(args.contract, args.output), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

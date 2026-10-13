"""Run KH ablation H: old selected-DNA/scorer with taxonomy v2.6 and history v4C.

The candidate generation and scorer remain legacy-v20-compatible. Current
taxonomy profiles are created answer-blind first, then supplied to the legacy
compatibility bridge before the v4C history route is applied. Benchmark labels
are loaded only after all candidate and history decisions have been written.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools.apply_test_aware_history_route_shadow import run as run_history_route
from tools.apply_legacy_history_picked_promotion import run as run_picked_promotion
from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.pathogen_normalization import canonical_key
from tools.run_kh_ablation_d_history_only import (
    after_rows,
    build_bridge,
    group_decisions,
    per_patient_changes,
    read_inventory,
    read_json,
    sha256_file,
    write_csv,
    write_json,
)
from tools.run_kh_ablation_e_taxonomy_only import score_mode


DEFAULT_A_ROOT = Path("outputs/runs/2026-09-24_KH_ablation_A_frozen_old_baseline_v1")
DEFAULT_PATIENT_ROOT = Path("outputs/patient_info_KH_0728_2Days")
DEFAULT_POLICY = Path("rules/test_aware_history_route_v4c_combined.json")
DEFAULT_PROMOTION_POLICY = Path("rules/legacy_history_picked_promotion_v1.json")
DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics/"
    "kh_answers_clinical_revision_20260918.csv"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-10-07_KH_ablation_H3_legacy_taxonomy_history_picked_v1"
)
SCHEMA_VERSION = "kh_ablation_h_legacy_taxonomy_history.v1"


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot identify patient from {path}")
    return int(match.group(1))


def load_taxonomy_profiles(patient_output: Path) -> dict[tuple[int, str], dict[str, Any]]:
    profiles: dict[tuple[int, str], dict[str, Any]] = {}
    for path in sorted(patient_output.glob("*.json")):
        patient = patient_number(path)
        payload = read_json(path)
        for candidate in payload.get("pathogen_candidates") or []:
            if not isinstance(candidate, dict):
                continue
            key = canonical_key(candidate.get("organism_name"))
            profile = candidate.get("taxonomy_profile")
            if key and isinstance(profile, dict) and profile:
                profiles[(patient, key)] = profile
    return profiles


def grouped_picked(rows: list[dict[str, Any]]) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = {}
    for row in rows:
        if not bool(row.get("picked")):
            continue
        grouped.setdefault(int(row["patient_id"]), []).append(str(row["organism_name"]))
    return grouped


def grouped_promoted_picked(path: Path) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if str(row.get("final_decision")) != "picked_shadow":
                continue
            grouped.setdefault(int(row["patient_id"]), []).append(str(row["organism_name"]))
    return grouped


def run(
    a_root: Path,
    patient_root: Path,
    policy_path: Path,
    promotion_policy_path: Path,
    answers_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")

    contract_path = a_root / "baseline_contract.json"
    contract = read_json(contract_path)
    patient_ids = [int(value) for value in contract["patient_selection"]]
    inventory = read_inventory(a_root / "artifact_inventory.csv")
    output_root.mkdir(parents=True)

    taxonomy_root = output_root / "01_legacy_scorer_taxonomy_v2_6"
    taxonomy_rows, taxonomy_summary = score_mode(
        mode="current-overlays",
        patient_ids=patient_ids,
        inventory=inventory,
        patient_root=patient_root,
        output_dir=taxonomy_root,
    )
    if taxonomy_summary["source_hash_mismatches"]:
        raise ValueError("Frozen A source hash validation failed")
    taxonomy_profiles = load_taxonomy_profiles(taxonomy_root / "patient_outputs")

    bridge_root = output_root / "02_legacy_route_bridge_taxonomy_v2_6"
    bridge_summary = build_bridge(
        a_root=a_root,
        patient_ids=patient_ids,
        inventory=inventory,
        output_dir=bridge_root,
        taxonomy_profiles=taxonomy_profiles,
        taxonomy_mode="current_taxonomy_v2_6",
    )
    if bridge_summary["source_hash_mismatches"]:
        raise ValueError("Frozen A bridge source hash validation failed")

    history_root = output_root / "03_history_route_v4c"
    history_summary = run_history_route(
        bridge_root,
        patient_root,
        history_root,
        policy_path=policy_path,
        phenotype_root=None,
        patients=set(patient_ids),
    )

    promotion_root = output_root / "04_history_supported_picked_gate"
    promotion_summary = run_picked_promotion(
        history_root, promotion_root, promotion_policy_path
    )

    # Only evaluation below this line can read benchmark labels.
    answers = load_answer_rows(answers_path)
    with (bridge_root / "baseline_candidate_decisions.csv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        before_history = list(csv.DictReader(handle))
    after_history = after_rows(history_root)
    evaluation_root = output_root / "05_evaluation"
    write_csv(evaluation_root / "final_candidate_decisions.csv", after_history)
    changes, summaries = per_patient_changes(before_history, after_history, answers)
    write_csv(evaluation_root / "per_patient_candidate_changes.csv", changes)
    write_csv(
        evaluation_root / "per_patient_summary.csv",
        summaries,
        [
            "patient_id", "before_visible_count", "after_visible_count", "tp_gained",
            "tp_lost", "fp_added", "fp_removed", "unlabeled_added", "unlabeled_removed",
        ],
    )

    strict = {"picked_shadow"}
    visible = {"picked_shadow", "review_high_priority", "review_context_needed"}
    metrics = {
        "E_taxonomy_strict_picked": evaluate_set(grouped_picked(taxonomy_rows), answers),
        "H_before_history_strict_picked": evaluate_set(
            group_decisions(before_history, strict), answers
        ),
        "H_after_history_strict_picked": evaluate_set(
            group_decisions(after_history, strict), answers
        ),
        "H3_after_history_promotion_strict_picked": evaluate_set(
            grouped_promoted_picked(promotion_root / "promotion_decisions.csv"), answers
        ),
        "H_before_history_complete_visible": evaluate_set(
            group_decisions(before_history, visible), answers
        ),
        "H_after_history_complete_visible": evaluate_set(
            group_decisions(after_history, visible), answers
        ),
    }
    expected = contract["frozen_metrics"]
    validation = {
        "taxonomy_scorer_picked_matches_A": (
            metrics["E_taxonomy_strict_picked"]["matched"] == expected["picked"]["matched"]
            and metrics["E_taxonomy_strict_picked"]["predicted_in_labeled_patients"]
            == expected["picked"]["outputs_in_labeled_patients"]
        ),
        "bridge_strict_matches_E": (
            metrics["H_before_history_strict_picked"] == metrics["E_taxonomy_strict_picked"]
        ),
        "bridge_visible_matches_A": (
            metrics["H_before_history_complete_visible"]["matched"]
            == expected["combined_picked_high_context"]["matched"]
            and metrics["H_before_history_complete_visible"]["predicted_in_labeled_patients"]
            == expected["combined_picked_high_context"]["outputs_in_labeled_patients"]
        ),
    }
    if not all(validation.values()):
        raise ValueError(f"H compatibility validation failed: {validation}")

    evaluation = {
        "schema_version": "kh_ablation_h_evaluation.v1",
        "scope": "Post-hoc benchmark evaluation after answer-blind taxonomy and history outputs.",
        "metrics": metrics,
        "validation": validation,
        "route_change_count": len(changes),
        "route_change_counts": dict(Counter(row["visibility_change"] for row in changes)),
        "history_rule_counts": dict(Counter(row["history_policy_family"] for row in changes)),
    }
    write_json(evaluation_root / "metrics.json", evaluation)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "experiment_id": "KH_ABLATION_H3_LEGACY_TAXONOMY_HISTORY_PICKED_V1",
        "status": "complete_answer_blind_shadow",
        "answer_blind_generation": True,
        "benchmark_evaluation_is_posthoc": True,
        "factors": {
            "mngs": contract["factors"]["mngs"],
            "scorer": (
                f"{contract['factors']['scorer']}+legacy_history_picked_promotion_v1"
            ),
            "taxonomy": "current_versioned_gap_auto_reviewed_overlays_v2_6",
            "history": "test_aware_history_route_v4c_combined",
        },
        "common_contract": {
            "patient_count": len(patient_ids),
            "patient_selection": [str(value) for value in patient_ids],
            "A_contract_path": str(contract_path.resolve()),
            "A_contract_sha256": sha256_file(contract_path),
            "answer_path": str(answers_path.resolve()),
            "answer_sha256": sha256_file(answers_path),
        },
        "taxonomy": {
            "snapshot": taxonomy_summary["taxonomy_snapshot"],
            "profile_count": len(taxonomy_profiles),
        },
        "history_policy": {
            "path": str(policy_path.resolve()),
            "sha256": sha256_file(policy_path),
            "policy_id": history_summary["policy_id"],
            "resolved_policy_sha256": history_summary["resolved_policy_sha256"],
        },
        "picked_promotion_policy": {
            "path": str(promotion_policy_path.resolve()),
            "sha256": sha256_file(promotion_policy_path),
            "policy_id": promotion_summary["policy_id"],
            "counts": promotion_summary["counts"],
        },
        "validation": validation,
        "counts": history_summary["counts"],
        "metrics": metrics,
        "route_change_counts": evaluation["route_change_counts"],
        "history_rule_counts": evaluation["history_rule_counts"],
        "outputs": {
            "taxonomy_summary_sha256": sha256_file(taxonomy_root / "summary.json"),
            "bridge_summary_sha256": sha256_file(bridge_root / "summary.json"),
            "history_summary_sha256": sha256_file(history_root / "summary.json"),
            "promotion_summary_sha256": sha256_file(promotion_root / "summary.json"),
            "promotion_decisions_sha256": sha256_file(
                promotion_root / "promotion_decisions.csv"
            ),
            "candidate_changes_sha256": sha256_file(
                evaluation_root / "per_patient_candidate_changes.csv"
            ),
            "metrics_sha256": sha256_file(evaluation_root / "metrics.json"),
        },
        "limitations": [
            "The mNGS source remains one selected-DNA analytical test per organism.",
            "Taxonomy and history are current, but the scorer interface remains legacy-compatible.",
            "History routing moves at most one level; a separate frozen gate may promote a qualifying Context-to-Priority candidate to Picked.",
            "Frozen A Luna review tiers are retained as starting routes, but no new Luna call is made.",
            "OBER/R5 and current Possible reporting are excluded from this mechanistic ablation.",
        ],
    }
    write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-root", type=Path, default=DEFAULT_A_ROOT)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--promotion-policy", type=Path, default=DEFAULT_PROMOTION_POLICY)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(
        json.dumps(
            run(
                args.a_root,
                args.patient_root,
                args.policy,
                args.promotion_policy,
                args.answers,
                args.output,
            ),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

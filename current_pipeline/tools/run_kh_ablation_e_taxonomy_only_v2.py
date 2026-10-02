"""Run refreshed KH ablation E against the immutable frozen-A control.

The historical A taxonomy file is no longer the working taxonomy file.  This
runner therefore does not attempt to recreate A with today's code.  It reads
the hash-checked preserved A decisions as the control and runs only the
treatment (the legacy scorer with the current taxonomy overlays).  This keeps
the ablation honest while allowing the current taxonomy snapshot to be tested.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

from tools import run_kh_ablation_e_taxonomy_only as e1
from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.run_legacy_scorer_on_multi_assay_shadow import score_rows


DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-10-01_KH_ablation_E_legacy_v20_taxonomy_v2_6_v2"
)
SCHEMA_VERSION = "kh_ablation_e_legacy_taxonomy_only.v2"


def preserved_control_rows(
    patient_ids: list[int],
    inventory: dict[tuple[int, str], dict[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    hash_mismatches: list[dict[str, Any]] = []
    for patient in patient_ids:
        record = inventory[(patient, "deterministic")]
        path = Path(record["path"])
        actual = e1.sha256_file(path)
        if actual.lower() != str(record["sha256"]).lower():
            hash_mismatches.append(
                {
                    "patient_id": patient,
                    "role": "deterministic",
                    "expected": record["sha256"],
                    "actual": actual,
                }
            )
        score = e1.read_json(path)
        patient_rows = score_rows(patient, score)
        candidates = {
            e1.canonical_key(item.get("organism_name")): item
            for item in score.get("pathogen_candidates") or []
            if isinstance(item, dict)
        }
        for row in patient_rows:
            candidate = candidates.get(e1.canonical_key(row.get("organism_name")), {})
            profile = candidate.get("taxonomy_profile") or {}
            row.update(
                {
                    "taxonomy_family": profile.get("primary_rule_family"),
                    "taxonomy_mapping_status": profile.get("mapping_status"),
                    "taxonomy_taxid": profile.get("taxid"),
                    "taxonomy_display_name": profile.get("display_name"),
                }
            )
        rows.extend(patient_rows)
    return rows, hash_mismatches


def run(
    a_root: Path,
    patient_root: Path,
    answers_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")

    contract_path = a_root / "baseline_contract.json"
    inventory_path = a_root / "artifact_inventory.csv"
    contract = e1.read_json(contract_path)
    patient_ids = [int(value) for value in contract["patient_selection"]]
    inventory = e1.read_inventory(inventory_path)
    answers = load_answer_rows(answers_path)
    output_root.mkdir(parents=True)

    control_rows, control_hash_mismatches = preserved_control_rows(patient_ids, inventory)
    if control_hash_mismatches:
        raise ValueError(f"Frozen A hash validation failed: {control_hash_mismatches}")
    control_dir = output_root / "control_frozen_A"
    control_dir.mkdir()
    e1.write_csv(control_dir / "candidate_decisions.csv", control_rows)

    treatment_rows, treatment_summary = e1.score_mode(
        mode="current-overlays",
        patient_ids=patient_ids,
        inventory=inventory,
        patient_root=patient_root,
        output_dir=output_root / "treatment_taxonomy_v2_6",
    )
    if treatment_summary["source_hash_mismatches"]:
        raise ValueError("Frozen A input hash validation failed in treatment")

    control_picked = e1.grouped_picked(control_rows)
    treatment_picked = e1.grouped_picked(treatment_rows)
    changes = e1.compare_rows(control_rows, treatment_rows, answers)
    evaluation_dir = output_root / "evaluation"
    e1.write_csv(evaluation_dir / "per_patient_candidate_changes.csv", changes)
    metrics = {
        "control_frozen_A_picked": evaluate_set(control_picked, answers),
        "treatment_current_taxonomy_picked": evaluate_set(treatment_picked, answers),
    }
    expected = contract["frozen_metrics"]["picked"]
    control_matches_contract = all(
        metrics["control_frozen_A_picked"][key] == expected[expected_key]
        for key, expected_key in (
            ("matched", "matched"),
            ("predicted_in_labeled_patients", "outputs_in_labeled_patients"),
            ("answer_count", "answer_count"),
        )
    )
    if not control_matches_contract:
        raise ValueError("Preserved control rows do not reproduce frozen A metrics")

    evaluation = {
        "schema_version": "kh_ablation_e_evaluation.v2",
        "scope": "Post-hoc evaluation after the treatment scorer output was frozen.",
        "metrics": metrics,
        "change_counts": dict(Counter(row["change"] for row in changes)),
        "changed_candidate_count": len(changes),
    }
    e1.write_json(evaluation_dir / "metrics.json", evaluation)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "experiment_id": "KH_ABLATION_E_LEGACY_V20_TAXONOMY_V2_6_V2",
        "status": "complete_frozen",
        "answer_blind_generation": True,
        "benchmark_evaluation_is_posthoc": True,
        "factors": {
            "mngs": contract["factors"]["mngs"],
            "scorer": contract["factors"]["scorer"],
            "taxonomy": "current_v2_6_versioned_gap_auto_reviewed_alias_overlays",
            "history": contract["factors"]["history"],
        },
        "common_contract": {
            "patient_count": len(patient_ids),
            "patient_selection": [str(value) for value in patient_ids],
            "A_contract_path": str(contract_path.resolve()),
            "A_contract_sha256": e1.sha256_file(contract_path),
            "answer_path": str(answers_path.resolve()),
            "answer_sha256": e1.sha256_file(answers_path),
        },
        "control_validation": {
            "mode": "immutable_preserved_A_outputs",
            "source_hash_mismatch_count": len(control_hash_mismatches),
            "metrics_match_frozen_A_contract": control_matches_contract,
            "reason_not_rerun": (
                "The historical taxonomy file is no longer the working file; "
                "A is explicitly frozen as an output comparator."
            ),
        },
        "implementation": {
            "runner_path": str(Path(__file__).resolve()),
            "runner_sha256": e1.sha256_file(Path(__file__)),
            "legacy_scorer_path": str(Path(e1.legacy.__file__).resolve()),
            "legacy_scorer_sha256": e1.sha256_file(Path(e1.legacy.__file__)),
        },
        "metrics": metrics,
        "change_counts": evaluation["change_counts"],
        "taxonomy_snapshots": {
            "control": contract["input_and_reference_artifacts"]["taxonomy_rules_snapshot"],
            "treatment": treatment_summary["taxonomy_snapshot"],
        },
        "outputs": {
            "control_candidate_sha256": e1.sha256_file(control_dir / "candidate_decisions.csv"),
            "treatment_candidate_sha256": e1.sha256_file(
                output_root / "treatment_taxonomy_v2_6/candidate_decisions.csv"
            ),
            "change_table_sha256": e1.sha256_file(
                evaluation_dir / "per_patient_candidate_changes.csv"
            ),
            "metrics_sha256": e1.sha256_file(evaluation_dir / "metrics.json"),
        },
        "limitations": [
            "This isolates current taxonomy inside the frozen legacy-v20-compatible scorer.",
            "The immutable A output, rather than a rerun with mutated source files, is the control.",
            "New mNGS, new history, current scorer, Possible reporting, and OBER are excluded.",
        ],
    }
    e1.write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-root", type=Path, default=e1.DEFAULT_A_ROOT)
    parser.add_argument("--patient-root", type=Path, default=e1.DEFAULT_PATIENT_ROOT)
    parser.add_argument("--answers", type=Path, default=e1.DEFAULT_ANSWERS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.a_root, args.patient_root, args.answers, args.output),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

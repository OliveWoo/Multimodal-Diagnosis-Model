"""Run fixed-policy RPM sensitivity and patient-cluster robustness analyses.

The sensitivity variants differ only in the predeclared per-test RPM threshold.
The leave-one-patient-out and bootstrap analyses do not retune the policy and
must not be described as external validation.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import random
from datetime import datetime
from pathlib import Path
from typing import Any

from tools.build_test_aware_possible_pathogen_shadow import run as run_shadow
from tools.evaluate_multi_assay_candidate_entry import load_answer_rows
from tools.evaluate_test_aware_possible_pathogen_shadow import (
    group_names,
    run as run_evaluation,
)
from tools.recalculate_kh_benchmark_metrics import match_type


DEFAULT_INPUT = Path(
    "outputs/runs/2026-09-28_KH_ablation_F_full_non_ober_v1/promotion_v4"
)
DEFAULT_POLICY = Path("rules/test_aware_unified_reporting_v2.json")
DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics/"
    "kh_answers_clinical_revision_20260918.csv"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-29_KH_unified_reporting_v2_2_robustness"
)


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0]) if rows else ["status"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def metric(tp: int, predicted: int, answers: int) -> dict[str, Any]:
    precision = tp / predicted if predicted else 0.0
    recall = tp / answers if answers else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall else 0.0
    )
    return {
        "matched": tp,
        "predicted": predicted,
        "answers": answers,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def patient_contributions(
    shadow_root: Path, answer_path: Path
) -> list[dict[str, Any]]:
    answers = load_answer_rows(answer_path)
    rows = read_csv(shadow_root / "possible_pathogen_decisions.csv")
    truthy = {"1", "true", "yes"}
    predictions = group_names(
        rows,
        lambda row: row["selected_for_complete_report"].casefold() in truthy,
    )
    output = []
    for patient_id, gold in sorted(answers.items()):
        names = predictions.get(patient_id, [])
        matched = sum(
            any(match_type(name, answer, genus_relaxed=False) for name in names)
            for answer in gold
        )
        output.append({
            "patient_id": patient_id,
            "matched": matched,
            "predicted": len(names),
            "answers": len(gold),
        })
    return output


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def earliest_collection_time(path: Path) -> datetime | None:
    payload = read_json(path)
    values: list[datetime] = []
    for candidate in payload.get("all_forwarded_candidates") or []:
        for event in candidate.get("linked_events") or []:
            raw = str(event.get("collected_time") or "").strip()
            if not raw:
                continue
            try:
                values.append(datetime.fromisoformat(raw))
            except ValueError:
                continue
    return min(values) if values else None


def run(
    input_root: Path,
    base_policy_path: Path,
    answer_path: Path,
    output_root: Path,
    *,
    thresholds: list[float],
    bootstrap_replicates: int,
    seed: int,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    policy_dir = output_root / "policies"
    policy_dir.mkdir()
    base_policy = read_json(base_policy_path)
    sensitivity_rows = []
    reference_shadow: Path | None = None

    for threshold in thresholds:
        label = str(threshold).replace(".", "p")
        policy = copy.deepcopy(base_policy)
        policy["policy_id"] = f"{base_policy['policy_id']}_rpm_{label}"
        policy["primary_report_compaction"]["low_normalized_burden"][
            "maximum_available_per_test_rpm_exclusive"
        ] = threshold
        policy_path = policy_dir / f"policy_rpm_{label}.json"
        write_json(policy_path, policy)
        shadow_root = output_root / f"rpm_{label}_shadow"
        evaluation_root = output_root / f"rpm_{label}_evaluation"
        run_shadow(input_root, shadow_root, policy_path=policy_path)
        result = run_evaluation(shadow_root, answer_path, evaluation_root)
        metrics = result["metrics"]["complete_report"]
        sensitivity_rows.append({
            "rpm_threshold_exclusive": threshold,
            "matched": metrics["matched"],
            "predicted_in_labeled_patients": metrics[
                "predicted_in_labeled_patients"
            ],
            "precision": metrics["precision"],
            "recall": metrics["recall"],
            "f1": metrics["f1"],
            "total_candidates_all_patients": metrics[
                "total_candidates_all_patients"
            ],
        })
        if threshold == 5:
            reference_shadow = shadow_root

    if reference_shadow is None:
        raise ValueError("Threshold list must include the fixed reference value 5")
    write_csv(output_root / "rpm_sensitivity.csv", sensitivity_rows)

    contributions = patient_contributions(reference_shadow, answer_path)
    totals = {
        key: sum(int(row[key]) for row in contributions)
        for key in ("matched", "predicted", "answers")
    }
    reference_metric = metric(
        totals["matched"], totals["predicted"], totals["answers"]
    )
    lopo_rows = []
    for row in contributions:
        result = metric(
            totals["matched"] - int(row["matched"]),
            totals["predicted"] - int(row["predicted"]),
            totals["answers"] - int(row["answers"]),
        )
        lopo_rows.append({
            "excluded_patient_id": row["patient_id"],
            "excluded_matched": row["matched"],
            "excluded_predictions": row["predicted"],
            "excluded_answers": row["answers"],
            **result,
            "precision_delta_from_full": result["precision"]
            - reference_metric["precision"],
            "recall_delta_from_full": result["recall"]
            - reference_metric["recall"],
            "f1_delta_from_full": result["f1"] - reference_metric["f1"],
        })
    write_csv(output_root / "leave_one_patient_out.csv", lopo_rows)

    rng = random.Random(seed)
    bootstrap_rows = []
    for replicate in range(1, bootstrap_replicates + 1):
        sampled = [rng.choice(contributions) for _ in contributions]
        result = metric(
            sum(int(row["matched"]) for row in sampled),
            sum(int(row["predicted"]) for row in sampled),
            sum(int(row["answers"]) for row in sampled),
        )
        bootstrap_rows.append({"replicate": replicate, **result})
    write_csv(output_root / "patient_cluster_bootstrap.csv", bootstrap_rows)

    date_by_patient: dict[int, datetime] = {}
    for path in (reference_shadow / "patient_outputs").glob("*.json"):
        found = path.name.split("NGS_patient_", 1)[-1].split("_", 1)[0]
        if not found.isdigit():
            continue
        timestamp = earliest_collection_time(path)
        if timestamp is not None:
            date_by_patient[int(found)] = timestamp
    dated = [
        row for row in contributions if int(row["patient_id"]) in date_by_patient
    ]
    dated.sort(key=lambda row: date_by_patient[int(row["patient_id"])])
    midpoint = len(dated) // 2
    temporal_rows = []
    for name, group in (("earlier", dated[:midpoint]), ("later", dated[midpoint:])):
        result = metric(
            sum(int(row["matched"]) for row in group),
            sum(int(row["predicted"]) for row in group),
            sum(int(row["answers"]) for row in group),
        )
        temporal_rows.append({
            "time_group": name,
            "patient_count": len(group),
            "first_collection_time": (
                date_by_patient[int(group[0]["patient_id"])].isoformat()
                if group else ""
            ),
            "last_collection_time": (
                date_by_patient[int(group[-1]["patient_id"])].isoformat()
                if group else ""
            ),
            **result,
        })
    write_csv(output_root / "temporal_subgroup_audit.csv", temporal_rows)

    bootstrap_summary = {}
    for key in ("precision", "recall", "f1"):
        values = [float(row[key]) for row in bootstrap_rows]
        bootstrap_summary[key] = {
            "lower_2_5_percentile": quantile(values, 0.025),
            "median": quantile(values, 0.5),
            "upper_97_5_percentile": quantile(values, 0.975),
        }
    lopo_summary = {
        key: {
            "minimum": min(float(row[key]) for row in lopo_rows),
            "maximum": max(float(row[key]) for row in lopo_rows),
        }
        for key in ("precision", "recall", "f1")
    }
    summary = {
        "scope": (
            "Internal robustness analysis of one cohort; no rule is selected "
            "from the best sensitivity result and this is not external validation."
        ),
        "fixed_reference_threshold": 5,
        "reference_complete_report": reference_metric,
        "rpm_sensitivity": sensitivity_rows,
        "leave_one_patient_out": lopo_summary,
        "patient_cluster_bootstrap": {
            "replicates": bootstrap_replicates,
            "seed": seed,
            **bootstrap_summary,
        },
        "temporal_subgroups": temporal_rows,
        "external_validation_status": (
            "not_performed_no_independent_compatible_per_test_DNA_RNA_cohort_"
            "was_supplied_to_this_run"
        ),
        "limitations": [
            "All rules were developed with access to the KH cohort and its benchmark labels.",
            "Leave-one-patient-out removes one patient after the policy is fixed; it does not retrain or provide independent validation.",
            "Bootstrap intervals quantify patient-sampling uncertainty, not dataset-shift robustness.",
            "The temporal split is a subgroup audit because both periods informed policy development.",
        ],
    }
    write_json(output_root / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--thresholds", type=float, nargs="+", default=[3, 5, 10])
    parser.add_argument("--bootstrap-replicates", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260929)
    args = parser.parse_args()
    print(json.dumps(run(
        args.input,
        args.policy,
        args.answers,
        args.output,
        thresholds=args.thresholds,
        bootstrap_replicates=args.bootstrap_replicates,
        seed=args.seed,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

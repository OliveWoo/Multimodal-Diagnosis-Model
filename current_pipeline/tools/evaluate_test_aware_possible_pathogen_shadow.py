"""Post-hoc evaluation of strict and possible-pathogen report sets.

The answer contract may contain both positive and explicit-no-pathogen cases.
Explicit negatives stay in the precision denominator; older positive-only
contracts remain supported for backwards compatibility.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable

from tools.build_test_aware_possible_pathogen_shadow import (
    ANALYTICAL_DIRECT_ROLE,
    STRICT_ROLE,
)
from tools.recalculate_kh_benchmark_metrics import (
    match_type,
    pairwise_match,
    split_names,
    unique_names,
)


DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics/"
    "kh_answers_clinical_revision_20260918.csv"
)
TRUTHY = {"1", "true", "yes"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path, fieldnames: list[str], rows: list[dict[str, Any]]
) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_answer_contract(path: Path) -> dict[int, dict[str, Any]]:
    """Load a status-aware contract, while preserving legacy CSV behavior."""
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"Answer contract is empty: {path}")
    has_status = "label_status" in rows[0]
    labels: dict[int, dict[str, Any]] = {}
    for row in rows:
        patient_id = int(row["patient_id"])
        if patient_id in labels:
            raise ValueError(f"Duplicate answer-contract patient: P{patient_id}")
        answers = unique_names(split_names(row.get("answers", "")))
        if not has_status:
            # Preserve historical behavior only when no explicit status exists.
            if not answers:
                continue
            status = "positive"
        else:
            status = row.get("label_status", "").strip()
            if status not in {"positive", "explicit_no_pathogen"}:
                raise ValueError(
                    f"Unsupported label_status for P{patient_id}: {status!r}"
                )
            if status == "positive" and not answers:
                raise ValueError(
                    f"Positive patient P{patient_id} has no answer organisms"
                )
            if status == "explicit_no_pathogen" and answers:
                raise ValueError(
                    f"No-pathogen patient P{patient_id} has answer organisms"
                )
        labels[patient_id] = {
            "status": status,
            "answers": answers,
            "cohort": row.get("cohort", "").strip(),
            "row": row,
        }
    if not labels:
        raise ValueError(f"Answer contract has no evaluable patients: {path}")
    return labels


def group_names(
    rows: list[dict[str, str]], predicate: Callable[[dict[str, str]], bool]
) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        if predicate(row):
            grouped[int(row["patient_id"])].append(row["organism_name"])
    return {patient: unique_names(names) for patient, names in grouped.items()}


def has_match(names: list[str], answer: str) -> bool:
    return any(match_type(name, answer, genus_relaxed=False) for name in names)


def calculate_metrics(
    predictions: dict[int, list[str]],
    labels: dict[int, dict[str, Any]],
    patient_ids: set[int],
) -> tuple[dict[str, Any], dict[int, dict[str, Any]]]:
    matched = predicted = answer_count = 0
    audits: dict[int, dict[str, Any]] = {}
    for patient_id in sorted(patient_ids):
        output = predictions.get(patient_id, [])
        answers = labels[patient_id]["answers"]
        audit = pairwise_match(output, answers, genus_relaxed=False)
        audits[patient_id] = audit
        matched += len(audit["matched"])
        predicted += len(output)
        answer_count += len(answers)
    false_positive = predicted - matched
    false_negative = answer_count - matched
    precision = matched / predicted if predicted else 0.0
    recall = matched / answer_count if answer_count else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    negative_ids = {
        patient_id
        for patient_id in patient_ids
        if labels[patient_id]["status"] == "explicit_no_pathogen"
    }
    return {
        "true_positive": matched,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "predicted": predicted,
        "answer_count": answer_count,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "patient_count": len(patient_ids),
        "positive_patient_count": sum(
            labels[patient_id]["status"] == "positive"
            for patient_id in patient_ids
        ),
        "explicit_no_pathogen_patient_count": len(negative_ids),
        "explicit_no_pathogen_patients_with_predictions": sorted(
            patient_id
            for patient_id in negative_ids
            if predictions.get(patient_id)
        ),
    }, audits


def build_strata(labels: dict[int, dict[str, Any]]) -> dict[str, set[int]]:
    by_cohort: dict[str, set[int]] = defaultdict(set)
    for patient_id, label in labels.items():
        if label["cohort"]:
            by_cohort[label["cohort"]].add(patient_id)
    strata: dict[str, set[int]] = {"all_labeled": set(labels)}
    for cohort in ("KH33", "two_hospital_A11", "two_hospital_B11"):
        if by_cohort.get(cohort):
            strata[cohort] = by_cohort[cohort]
    if by_cohort.get("two_hospital_A11") or by_cohort.get("two_hospital_B11"):
        strata["two_hospital_complete"] = (
            by_cohort.get("two_hospital_A11", set())
            | by_cohort.get("two_hospital_B11", set())
        )
    if by_cohort.get("KH33") or by_cohort.get("two_hospital_A11"):
        strata["without_B11"] = (
            by_cohort.get("KH33", set())
            | by_cohort.get("two_hospital_A11", set())
        )
    return strata


def run(shadow_root: Path, answer_path: Path, output_dir: Path) -> dict[str, Any]:
    labels = load_answer_contract(answer_path)
    patient_ids = set(labels)
    positive_ids = {
        patient_id
        for patient_id, label in labels.items()
        if label["status"] == "positive"
    }
    rows = read_csv(shadow_root / "possible_pathogen_decisions.csv")
    groups = {
        "strict_picked": group_names(
            rows, lambda row: row["possible_reporting_role"] == STRICT_ROLE
        ),
        "strict_plus_analytical_or_direct_possible": group_names(
            rows,
            lambda row: row["possible_reporting_role"]
            in {STRICT_ROLE, ANALYTICAL_DIRECT_ROLE},
        ),
        "complete_report": group_names(
            rows,
            lambda row: row["selected_for_complete_report"].casefold() in TRUTHY,
        ),
        "all_forwarded_visible": group_names(rows, lambda row: True),
    }

    metrics: dict[str, Any] = {}
    legacy_metrics: dict[str, Any] = {}
    endpoint_audits: dict[str, dict[int, dict[str, Any]]] = {}
    for endpoint, predictions in groups.items():
        metrics[endpoint], endpoint_audits[endpoint] = calculate_metrics(
            predictions, labels, patient_ids
        )
        legacy_metrics[endpoint], _ = calculate_metrics(
            predictions, labels, positive_ids
        )

    strata = build_strata(labels)
    stratified_metrics = {
        stratum: {
            endpoint: calculate_metrics(predictions, labels, ids)[0]
            for endpoint, predictions in groups.items()
        }
        for stratum, ids in strata.items()
    }

    routes: list[dict[str, Any]] = []
    for patient_id in sorted(positive_ids):
        patient_rows = [
            row for row in rows if int(row["patient_id"]) == patient_id
        ]
        for answer in labels[patient_id]["answers"]:
            matched = next(
                (
                    row
                    for row in patient_rows
                    if match_type(row["organism_name"], answer, genus_relaxed=False)
                ),
                None,
            )
            routes.append({
                "patient_id": patient_id,
                "cohort": labels[patient_id]["cohort"],
                "answer": answer,
                "clinical_decision": (
                    matched["clinical_decision"]
                    if matched
                    else "absent_from_forwarded_candidates"
                ),
                "possible_reporting_role": (
                    matched["possible_reporting_role"]
                    if matched
                    else "absent_from_forwarded_candidates"
                ),
                "possible_route": matched["possible_route"] if matched else "",
                "selected_for_complete_report": (
                    matched["selected_for_complete_report"] if matched else "False"
                ),
            })

    prediction_audit: list[dict[str, Any]] = []
    explicit_negative_predictions: list[dict[str, Any]] = []
    for row in rows:
        if row["selected_for_complete_report"].casefold() not in TRUTHY:
            continue
        patient_id = int(row["patient_id"])
        label = labels.get(patient_id)
        matched_answer = next(
            (
                answer
                for answer in (label["answers"] if label else [])
                if match_type(row["organism_name"], answer, genus_relaxed=False)
            ),
            None,
        )
        if label is None:
            outcome = "unlabeled_patient"
        elif label["status"] == "explicit_no_pathogen":
            outcome = "unmatched_in_explicit_no_pathogen_patient"
            explicit_negative_predictions.append({
                "patient_id": patient_id,
                "cohort": label["cohort"],
                "organism_name": row["organism_name"],
                "possible_reporting_role": row["possible_reporting_role"],
            })
        elif matched_answer:
            outcome = "matched_benchmark_answer"
        else:
            outcome = "unmatched_in_labeled_patient"
        prediction_audit.append({
            "patient_id": patient_id,
            "cohort": label["cohort"] if label else "",
            "label_status": label["status"] if label else "unlabeled",
            "organism_name": row["organism_name"],
            "taxonomy_family": row["taxonomy_family"],
            "clinical_decision": row["clinical_decision"],
            "possible_reporting_role": row["possible_reporting_role"],
            "possible_route": row["possible_route"],
            "posthoc_outcome": outcome,
            "matched_answer": matched_answer or "",
        })

    missed_answers: list[dict[str, Any]] = []
    for endpoint, audits in endpoint_audits.items():
        for patient_id, audit in sorted(audits.items()):
            for answer in audit["answer_only"]:
                missed_answers.append({
                    "endpoint": endpoint,
                    "patient_id": patient_id,
                    "cohort": labels[patient_id]["cohort"],
                    "answer": answer,
                })

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        output_dir / "answer_routes.csv",
        [
            "patient_id", "cohort", "answer", "clinical_decision",
            "possible_reporting_role", "possible_route",
            "selected_for_complete_report",
        ],
        routes,
    )
    write_csv(
        output_dir / "prediction_audit.csv",
        [
            "patient_id", "cohort", "label_status", "organism_name",
            "taxonomy_family", "clinical_decision", "possible_reporting_role",
            "possible_route", "posthoc_outcome", "matched_answer",
        ],
        prediction_audit,
    )
    write_csv(
        output_dir / "explicit_no_pathogen_predictions.csv",
        ["patient_id", "cohort", "organism_name", "possible_reporting_role"],
        explicit_negative_predictions,
    )
    write_csv(
        output_dir / "missed_answers.csv",
        ["endpoint", "patient_id", "cohort", "answer"],
        missed_answers,
    )
    metric_rows: list[dict[str, Any]] = []
    for stratum, endpoint_values in stratified_metrics.items():
        for endpoint, value in endpoint_values.items():
            metric_rows.append({"stratum": stratum, "endpoint": endpoint, **value})
    write_csv(
        output_dir / "stratified_metrics.csv",
        [
            "stratum", "endpoint", "true_positive", "false_positive",
            "false_negative", "predicted", "answer_count", "precision",
            "recall", "f1", "patient_count", "positive_patient_count",
            "explicit_no_pathogen_patient_count",
            "explicit_no_pathogen_patients_with_predictions",
        ],
        metric_rows,
    )

    result = {
        "schema_version": "test_aware_posthoc_evaluation.v2",
        "scope": (
            "Post-hoc evaluation only; the reporting policy does not read "
            "benchmark answers."
        ),
        "matching": (
            "Exact name, approved alias, or approved group member; "
            "no genus-relaxed matching."
        ),
        "shadow_root": str(shadow_root.resolve()),
        "answer_path": str(answer_path.resolve()),
        "metrics": metrics,
        "legacy_positive_patients_only": legacy_metrics,
        "stratified_metrics": stratified_metrics,
        "answer_contract": {
            "patient_count": len(patient_ids),
            "positive_patient_count": len(positive_ids),
            "explicit_no_pathogen_patient_count": len(patient_ids - positive_ids),
            "answer_count": sum(len(label["answers"]) for label in labels.values()),
        },
        "answer_role_counts": dict(
            Counter(row["possible_reporting_role"] for row in routes)
        ),
        "complete_report_posthoc_outcome_counts": dict(
            Counter(row["posthoc_outcome"] for row in prediction_audit)
        ),
        "possible_only_posthoc_outcome_counts": dict(
            Counter(
                row["posthoc_outcome"]
                for row in prediction_audit
                if row["possible_reporting_role"] != STRICT_ROLE
            )
        ),
        "candidate_count": len(rows),
        "notes": [
            "metrics includes explicit-no-pathogen patients in the precision denominator.",
            "legacy_positive_patients_only is retained only for historical comparison.",
        ],
    }
    (output_dir / "metrics.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shadow_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    args = parser.parse_args()
    print(json.dumps(
        run(args.shadow_root, args.answers, args.output_dir),
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()

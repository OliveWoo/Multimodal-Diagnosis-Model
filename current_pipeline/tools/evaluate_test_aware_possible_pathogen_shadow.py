"""Post-hoc evaluation of strict and possible-pathogen report sets."""

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
from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.recalculate_kh_benchmark_metrics import match_type, unique_names


DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics/kh_answers_clinical_revision_20260918.csv"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


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


def run(shadow_root: Path, answer_path: Path, output_dir: Path) -> dict[str, Any]:
    answers = load_answer_rows(answer_path)
    rows = read_csv(shadow_root / "possible_pathogen_decisions.csv")
    truthy = {"1", "true", "yes"}
    groups = {
        "strict_picked": group_names(
            rows, lambda row: row["possible_reporting_role"] == STRICT_ROLE
        ),
        "strict_plus_analytical_or_direct_possible": group_names(
            rows, lambda row: row["possible_reporting_role"] in {
                STRICT_ROLE, ANALYTICAL_DIRECT_ROLE,
            },
        ),
        "complete_report": group_names(
            rows,
            lambda row: row["selected_for_complete_report"].casefold() in truthy,
        ),
        "all_forwarded_visible": group_names(rows, lambda row: True),
    }
    metrics = {name: evaluate_set(values, answers) for name, values in groups.items()}

    routes = []
    for patient, gold in sorted(answers.items()):
        patient_rows = [row for row in rows if int(row["patient_id"]) == patient]
        for answer in gold:
            matched = next(
                (row for row in patient_rows
                 if match_type(row["organism_name"], answer, genus_relaxed=False)),
                None,
            )
            routes.append({
                "patient_id": patient,
                "answer": answer,
                "clinical_decision": (
                    matched["clinical_decision"] if matched else "absent_from_forwarded_candidates"
                ),
                "possible_reporting_role": (
                    matched["possible_reporting_role"] if matched else "absent_from_forwarded_candidates"
                ),
                "possible_route": matched["possible_route"] if matched else "",
                "selected_for_complete_report": (
                    matched["selected_for_complete_report"] if matched else "False"
                ),
            })

    prediction_audit = []
    for row in rows:
        if row["selected_for_complete_report"].casefold() not in truthy:
            continue
        patient = int(row["patient_id"])
        gold = answers.get(patient)
        matched_answer = next(
            (
                answer for answer in (gold or [])
                if match_type(row["organism_name"], answer, genus_relaxed=False)
            ),
            None,
        )
        if gold is None:
            outcome = "unlabeled_patient"
        elif matched_answer:
            outcome = "matched_benchmark_answer"
        else:
            outcome = "unmatched_in_labeled_patient"
        prediction_audit.append({
            "patient_id": patient,
            "organism_name": row["organism_name"],
            "taxonomy_family": row["taxonomy_family"],
            "clinical_decision": row["clinical_decision"],
            "possible_reporting_role": row["possible_reporting_role"],
            "possible_route": row["possible_route"],
            "posthoc_outcome": outcome,
            "matched_answer": matched_answer or "",
        })

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "answer_routes.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        fields = [
            "patient_id", "answer", "clinical_decision", "possible_reporting_role",
            "possible_route", "selected_for_complete_report",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(routes)
    with (output_dir / "prediction_audit.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        fields = [
            "patient_id", "organism_name", "taxonomy_family", "clinical_decision",
            "possible_reporting_role", "possible_route", "posthoc_outcome",
            "matched_answer",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(prediction_audit)
    result = {
        "scope": "Post-hoc evaluation only; the possible-pathogen policy does not read benchmark answers.",
        "matching": "Exact name, approved alias, or approved group member; no genus-relaxed matching.",
        "shadow_root": str(shadow_root.resolve()),
        "answer_path": str(answer_path.resolve()),
        "metrics": metrics,
        "answer_role_counts": dict(Counter(
            row["possible_reporting_role"] for row in routes
        )),
        "complete_report_posthoc_outcome_counts": dict(Counter(
            row["posthoc_outcome"] for row in prediction_audit
        )),
        "possible_only_posthoc_outcome_counts": dict(Counter(
            row["posthoc_outcome"] for row in prediction_audit
            if row["possible_reporting_role"] != STRICT_ROLE
        )),
        "candidate_count": len(rows),
    }
    (output_dir / "metrics.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
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

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools.recalculate_kh_benchmark_metrics import match_type, pairwise_match, split_names, unique_names


DEFAULT_ANSWERS = Path("outputs/runs/2026-09-18_KH_answer_revision_metrics/kh_answers_clinical_revision_20260918.csv")
DEFAULT_COMPARISON = Path("outputs/runs/2026-09-23_KH_multi_assay_positive_reads_shadow/selected_dna_comparison.csv")
DEFAULT_OLD_POOL = Path("outputs/runs/2026-09-22_KH_multi_assay_screening_shadow_taxonomy_v2/shadow_candidate_pool.csv")
DEFAULT_NEW_DECISIONS = Path("outputs/runs/2026-09-23_KH_independent_multi_assay_entry_shadow/all_decisions.csv")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def group_names(rows: list[dict[str, str]], *, predicate=lambda row: True) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        if predicate(row):
            grouped[int(row["patient_id"])].append(row["organism_name"])
    return {patient_id: unique_names(names) for patient_id, names in grouped.items()}


def load_answer_rows(path: Path) -> dict[int, list[str]]:
    return {
        int(row["patient_id"]): split_names(row["answers"])
        for row in read_csv(path)
        if split_names(row["answers"])
    }


def evaluate_set(predictions: dict[int, list[str]], answers: dict[int, list[str]]) -> dict[str, Any]:
    matched = predicted = answer_count = 0
    for patient_id, gold in answers.items():
        output = predictions.get(patient_id, [])
        result = pairwise_match(output, gold, genus_relaxed=False)
        matched += len(result["matched"])
        predicted += len(output)
        answer_count += len(gold)
    precision = matched / predicted if predicted else 0.0
    recall = matched / answer_count if answer_count else 0.0
    return {
        "matched": matched,
        "predicted_in_labeled_patients": predicted,
        "answer_count": answer_count,
        "precision": precision,
        "recall": recall,
        "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0,
        "patients_with_answers": len(answers),
        "total_candidates_all_patients": sum(len(names) for names in predictions.values()),
    }


def has_answer(names: list[str], answer: str) -> bool:
    return any(match_type(name, answer, genus_relaxed=False) for name in names)


def build_audit(
    answer_path: Path, comparison_path: Path, old_pool_path: Path, new_decisions_path: Path, output_dir: Path
) -> dict[str, Any]:
    answers = load_answer_rows(answer_path)
    comparison = read_csv(comparison_path)
    old_pool = read_csv(old_pool_path)
    decisions = read_csv(new_decisions_path)
    old_baseline = group_names(
        comparison, predicate=lambda row: row["comparison_status"] == "retained_from_selected_dna"
    )
    old_expanded = group_names(old_pool)
    scorer = group_names(decisions, predicate=lambda row: row["disposition"] == "scorer_entry")
    review = group_names(decisions, predicate=lambda row: row["disposition"] == "clinical_review")
    qc = group_names(decisions, predicate=lambda row: row["disposition"] == "qc_only")
    forwarded = group_names(decisions, predicate=lambda row: row["disposition"] != "qc_only")
    all_positive = group_names(decisions)
    selected_positive = group_names(
        decisions, predicate=lambda row: int(row["selected_positive_test_count"]) > 0
    )
    groups = {
        "old_selected_dna_baseline": old_baseline,
        "old_expanded_shadow_pool": old_expanded,
        "new_scorer_entry": scorer,
        "new_clinical_review_only": review,
        "new_scorer_or_review": forwarded,
        "new_all_positive_union": all_positive,
    }
    metrics = {name: evaluate_set(predictions, answers) for name, predictions in groups.items()}

    answer_routes: list[dict[str, Any]] = []
    route_order = {"scorer_entry": 0, "clinical_review": 1, "qc_only": 2}
    for patient_id, gold in sorted(answers.items()):
        for organism in gold:
            matching_rows = [
                row for row in decisions
                if int(row["patient_id"]) == patient_id
                and match_type(row["organism_name"], organism, genus_relaxed=False)
            ]
            matching_rows.sort(key=lambda row: route_order[row["disposition"]])
            chosen = matching_rows[0] if matching_rows else None
            new_route = chosen["disposition"] if chosen else "absent_from_positive_union"
            answer_routes.append(
                {
                    "patient_id": patient_id,
                    "answer": organism,
                    "new_route": new_route,
                    "matching_organism": chosen["organism_name"] if chosen else "",
                    "case_review_id": chosen["case_review_id"] if chosen else "",
                    "taxonomy_family": chosen["taxonomy_family"] if chosen else "",
                    "rule_ids": chosen["rule_ids"] if chosen else "",
                    "source_status": (
                        "selected_positive"
                        if has_answer(selected_positive.get(patient_id, []), organism)
                        else "filtered_positive_only"
                        if has_answer(all_positive.get(patient_id, []), organism)
                        else "absent_from_positive_union"
                    ),
                    "old_selected_dna_baseline_hit": has_answer(old_baseline.get(patient_id, []), organism),
                    "old_expanded_shadow_pool_hit": has_answer(old_expanded.get(patient_id, []), organism),
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "answer_routes.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(answer_routes[0]))
        writer.writeheader()
        writer.writerows(answer_routes)
    scorer_hits: list[dict[str, Any]] = []
    for patient_id, names in sorted(scorer.items()):
        gold = answers.get(patient_id)
        matched_outputs = (
            {item["output"]: item["answer"] for item in pairwise_match(names, gold, genus_relaxed=False)["matched"]}
            if gold else {}
        )
        for name in names:
            scorer_hits.append(
                {
                    "patient_id": patient_id,
                    "organism_name": name,
                    "answer_status": "matched" if name in matched_outputs else "unmatched" if gold else "unlabeled_patient",
                    "matched_answer": matched_outputs.get(name, ""),
                }
            )
    with (output_dir / "scorer_candidate_hits.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["patient_id", "organism_name", "answer_status", "matched_answer"])
        writer.writeheader()
        writer.writerows(scorer_hits)
    summary = {
        "scope": "Candidate-entry comparison, not final model Picked accuracy",
        "matching": "Exact name, approved alias, or approved group member; no genus-relaxed matches",
        "unlabeled_patients_excluded_from_precision": True,
        "answer_path": str(answer_path.resolve()),
        "new_decisions_path": str(new_decisions_path.resolve()),
        "metrics": metrics,
        "answer_route_counts": dict(Counter(row["new_route"] for row in answer_routes)),
        "answer_source_counts": dict(Counter(row["source_status"] for row in answer_routes)),
    }
    (output_dir / "metrics.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Answer audit for old and independent DNA/RNA candidate-entry pools")
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--comparison", type=Path, default=DEFAULT_COMPARISON)
    parser.add_argument("--old-pool", type=Path, default=DEFAULT_OLD_POOL)
    parser.add_argument("--new-decisions", type=Path, default=DEFAULT_NEW_DECISIONS)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_audit(args.answers, args.comparison, args.old_pool, args.new_decisions, args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

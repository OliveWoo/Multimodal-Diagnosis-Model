"""Post-hoc benchmark evaluation for the answer-blind clinical shadow."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable

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


def run(clinical_root: Path, answer_path: Path, output_dir: Path) -> dict[str, Any]:
    answers = load_answer_rows(answer_path)
    rows = read_csv(clinical_root / "clinical_decisions.csv")
    groups = {
        "incoming_picked_shadow": group_names(rows, lambda row: row["incoming_decision"] == "picked_shadow"),
        "clinical_picked_shadow": group_names(rows, lambda row: row["clinical_decision"] == "picked_shadow"),
        "clinical_picked_plus_high_review": group_names(
            rows, lambda row: row["clinical_decision"] in {"picked_shadow", "review_high_priority"}),
        "clinical_all_forwarded_visible": group_names(rows, lambda row: True),
    }
    metrics = {name: evaluate_set(values, answers) for name, values in groups.items()}
    routes = []
    tier_order = (
        ("picked_shadow", groups["clinical_picked_shadow"]),
        ("review_high_priority", group_names(
            rows, lambda row: row["clinical_decision"] == "review_high_priority")),
        ("review_context_needed", group_names(
            rows, lambda row: row["clinical_decision"] == "review_context_needed")),
    )
    for patient, gold in sorted(answers.items()):
        for answer in gold:
            route = "absent_from_forwarded_candidates"
            for name, predictions in tier_order:
                if has_match(predictions.get(patient, []), answer):
                    route = name
                    break
            routes.append({"patient_id": patient, "answer": answer, "clinical_route": route})

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "answer_routes.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["patient_id", "answer", "clinical_route"])
        writer.writeheader()
        writer.writerows(routes)
    result = {
        "scope": "Post-hoc evaluation only; the clinical shadow and rules were built without benchmark answers.",
        "matching": "Exact name, approved alias, or approved group member; no genus-relaxed matching.",
        "clinical_root": str(clinical_root.resolve()),
        "answer_path": str(answer_path.resolve()),
        "metrics": metrics,
        "answer_route_counts": dict(Counter(row["clinical_route"] for row in routes)),
        "candidate_count": len(rows),
    }
    (output_dir / "metrics.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("clinical_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    args = parser.parse_args()
    print(json.dumps(run(args.clinical_root, args.answers, args.output_dir),
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

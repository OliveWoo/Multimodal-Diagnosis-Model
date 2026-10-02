from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable

from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.recalculate_kh_benchmark_metrics import match_type, pairwise_match, unique_names


DEFAULT_ANSWERS = Path("outputs/runs/2026-09-18_KH_answer_revision_metrics/kh_answers_clinical_revision_20260918.csv")
DEFAULT_OLD_SUFFIX = "mNGS_max_deterministic_selected_dna_20260918"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def group_names(rows: list[dict[str, str]], predicate: Callable[[dict[str, str]], bool]) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        if predicate(row):
            grouped[int(row["patient_id"])].append(row["organism_name"])
    return {patient_id: unique_names(names) for patient_id, names in grouped.items()}


def merge_groups(*groups: dict[int, list[str]]) -> dict[int, list[str]]:
    patient_ids = {patient_id for group in groups for patient_id in group}
    return {
        patient_id: unique_names([name for group in groups for name in group.get(patient_id, [])])
        for patient_id in patient_ids
    }


def old_deterministic_picked(patient_root: Path, suffix: str) -> dict[int, list[str]]:
    output: dict[int, list[str]] = {}
    for patient_dir in patient_root.glob("NGS_patient_*_json"):
        try:
            patient_id = int(patient_dir.name.split("_")[2])
        except (IndexError, ValueError):
            continue
        path = patient_dir / "summary_outputs" / f"NGS_patient_{patient_id}_{suffix}.json"
        if not path.exists():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        names = [
            row.get("organism_name")
            for row in ((payload.get("best_available_summary") or {}).get("picked_pathogens") or [])
            if isinstance(row, dict) and row.get("organism_name")
        ]
        output[patient_id] = unique_names(names)
    return output


def has_match(names: list[str], answer: str) -> bool:
    return any(match_type(name, answer, genus_relaxed=False) for name in names)


def build_evaluation(
    shadow_dir: Path,
    patient_root: Path,
    answer_path: Path,
    output_dir: Path,
    *,
    old_suffix: str = DEFAULT_OLD_SUFFIX,
) -> dict[str, Any]:
    answers = load_answer_rows(answer_path)
    candidate_rows = read_csv(shadow_dir / "patient_organism_decisions.csv")
    hospital_rows = read_csv(shadow_dir / "hospital_only_decisions.csv")
    mngs_picked = group_names(candidate_rows, lambda row: row["decision"] == "picked_shadow")
    hospital_picked = group_names(hospital_rows, lambda row: row["decision"] == "picked_shadow")
    combined_picked = merge_groups(mngs_picked, hospital_picked)
    all_high = merge_groups(
        group_names(candidate_rows, lambda row: row["decision"] == "review_high_priority"),
        group_names(hospital_rows, lambda row: row["decision"] == "review_high_priority"),
    )
    context = merge_groups(
        group_names(candidate_rows, lambda row: row["decision"] == "review_context_needed"),
        group_names(hospital_rows, lambda row: row["decision"] == "review_context_needed"),
    )
    retained_low = group_names(candidate_rows, lambda row: row["decision"] == "review_low_specificity")
    all_routed = merge_groups(combined_picked, all_high, context, retained_low)
    picked_high = merge_groups(combined_picked, all_high)
    picked_high_context = merge_groups(picked_high, context)
    old_picked = old_deterministic_picked(patient_root, old_suffix)
    groups = {
        "old_selected_dna_deterministic_picked": old_picked,
        "test_aware_mngs_picked_shadow": mngs_picked,
        "test_aware_combined_picked_shadow": combined_picked,
        "test_aware_picked_plus_high_review": picked_high,
        "test_aware_clinical_scorer_forward": picked_high_context,
        "test_aware_picked_high_context_visible": picked_high_context,
        "test_aware_all_positive_routed_candidates": all_routed,
    }
    metrics = {name: evaluate_set(predictions, answers) for name, predictions in groups.items()}

    answer_routes = []
    route_sets = [
        ("picked_shadow", combined_picked),
        ("review_high_priority", all_high),
        ("review_context_needed", context),
        ("retained_low_not_forwarded", retained_low),
    ]
    for patient_id, gold in sorted(answers.items()):
        for answer in gold:
            route = "absent_from_positive_routed_candidates"
            for route_name, names_by_patient in route_sets:
                if has_match(names_by_patient.get(patient_id, []), answer):
                    route = route_name
                    break
            answer_routes.append({"patient_id": patient_id, "answer": answer, "test_aware_route": route})

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "answer_routes.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(answer_routes[0]))
        writer.writeheader()
        writer.writerows(answer_routes)
    picked_rows = []
    for patient_id, names in sorted(combined_picked.items()):
        gold = answers.get(patient_id, [])
        matched = {row["output"]: row["answer"] for row in pairwise_match(names, gold, genus_relaxed=False)["matched"]}
        for name in names:
            picked_rows.append({
                "patient_id": patient_id,
                "organism_name": name,
                "answer_status": "matched" if name in matched else "unmatched" if gold else "unlabeled_patient",
                "matched_answer": matched.get(name, ""),
            })
    with (output_dir / "picked_shadow_audit.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["patient_id", "organism_name", "answer_status", "matched_answer"])
        writer.writeheader()
        writer.writerows(picked_rows)

    summary = {
        "scope": "Post-hoc evaluation of an answer-blind test-aware deterministic shadow; not production diagnostic accuracy.",
        "matching": "Exact name, approved alias, or approved group member; no genus-relaxed matches.",
        "unlabeled_patients_excluded_from_precision": True,
        "shadow_dir": str(shadow_dir.resolve()),
        "answer_path": str(answer_path.resolve()),
        "old_deterministic_suffix": old_suffix,
        "metrics": metrics,
        "answer_route_counts": dict(Counter(row["test_aware_route"] for row in answer_routes)),
        "patients_without_combined_picked": sorted(set(answers) - set(combined_picked)),
    }
    (output_dir / "metrics.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate the test-aware deterministic shadow after scoring")
    parser.add_argument("shadow_dir", type=Path)
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--old-suffix", default=DEFAULT_OLD_SUFFIX)
    args = parser.parse_args()
    print(json.dumps(build_evaluation(
        args.shadow_dir,
        args.patient_root,
        args.answers,
        args.output_dir,
        old_suffix=args.old_suffix,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

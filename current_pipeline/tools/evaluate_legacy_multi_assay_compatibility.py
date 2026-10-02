"""Compare legacy-scorer/new-mNGS compatibility results with both baselines."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.evaluate_test_aware_deterministic_shadow import old_deterministic_picked
from tools.recalculate_kh_benchmark_metrics import match_type, unique_names


DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics/kh_answers_clinical_revision_20260918.csv"
)
DEFAULT_OLD_SUFFIX = "mNGS_max_deterministic_selected_dna_20260918"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def group(rows: list[dict[str, str]], predicate) -> dict[int, list[str]]:
    values: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        if predicate(row):
            values[int(row["patient_id"])].append(row["organism_name"])
    return {patient: unique_names(names) for patient, names in values.items()}


def level_rank(value: Any) -> int:
    text = str(value or "")
    for number in range(1, 6):
        if str(number) in text:
            return number
    return 99


def has_match(names: list[str], answer: str) -> bool:
    return any(match_type(name, answer, genus_relaxed=False) for name in names)


def set_changes(
    old: dict[int, list[str]], new: dict[int, list[str]]
) -> list[dict[str, Any]]:
    rows = []
    for patient in sorted(set(old) | set(new)):
        old_names = old.get(patient, [])
        new_names = new.get(patient, [])
        for name in new_names:
            if not has_match(old_names, name):
                rows.append({"patient_id": patient, "change": "added_by_multi_assay", "organism_name": name})
        for name in old_names:
            if not has_match(new_names, name):
                rows.append({"patient_id": patient, "change": "removed_or_replaced", "organism_name": name})
    return rows


def run(
    compatibility_root: Path, patient_root: Path, clinical_root: Path,
    answer_path: Path, output_dir: Path,
) -> dict[str, Any]:
    answers = load_answer_rows(answer_path)
    legacy_rows = read_csv(compatibility_root / "legacy_multi_assay_candidate_decisions.csv")
    clinical_rows = read_csv(clinical_root / "clinical_decisions.csv")
    groups = {
        "old_selected_dna_legacy_picked": old_deterministic_picked(patient_root, DEFAULT_OLD_SUFFIX),
        "new_multi_assay_through_legacy_picked": group(legacy_rows, lambda row: row["picked"].lower() == "true"),
        "new_multi_assay_through_legacy_level1_3_visible": group(
            legacy_rows, lambda row: level_rank(row["integrated_level"]) <= 3),
        "new_multi_assay_through_legacy_all_scored": group(legacy_rows, lambda row: True),
        "new_test_aware_clinical_picked": group(
            clinical_rows, lambda row: row["clinical_decision"] == "picked_shadow"),
        "new_test_aware_clinical_picked_plus_high": group(
            clinical_rows, lambda row: row["clinical_decision"] in {"picked_shadow", "review_high_priority"}),
        "new_test_aware_all_forwarded": group(clinical_rows, lambda row: True),
    }
    metrics = {name: evaluate_set(values, answers) for name, values in groups.items()}
    picked_changes = set_changes(
        groups["old_selected_dna_legacy_picked"],
        groups["new_multi_assay_through_legacy_picked"],
    )
    routes = []
    legacy_picked = groups["new_multi_assay_through_legacy_picked"]
    legacy_visible = groups["new_multi_assay_through_legacy_level1_3_visible"]
    legacy_all = groups["new_multi_assay_through_legacy_all_scored"]
    for patient, gold in sorted(answers.items()):
        for answer in gold:
            if has_match(legacy_picked.get(patient, []), answer):
                route = "picked"
            elif has_match(legacy_visible.get(patient, []), answer):
                route = "level1_3_not_picked"
            elif has_match(legacy_all.get(patient, []), answer):
                route = "level4_5_or_excluded"
            else:
                route = "absent"
            routes.append({"patient_id": patient, "answer": answer, "legacy_multi_assay_route": route})

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "answer_routes.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["patient_id", "answer", "legacy_multi_assay_route"])
        writer.writeheader()
        writer.writerows(routes)
    with (output_dir / "picked_changes_vs_selected_dna.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=["patient_id", "change", "organism_name"])
        writer.writeheader()
        writer.writerows(picked_changes)
    result = {
        "scope": "Post-hoc compatibility evaluation; neither scorer read benchmark answers.",
        "compatibility_interpretation": (
            "This tests the current legacy-shaped deterministic scorer on per-test adapted multi-assay data; "
            "it is not a bit-for-bit historical replay."
        ),
        "metrics": metrics,
        "answer_route_counts": dict(Counter(row["legacy_multi_assay_route"] for row in routes)),
        "picked_change_counts": dict(Counter(row["change"] for row in picked_changes)),
        "compatibility_root": str(compatibility_root.resolve()),
        "clinical_root": str(clinical_root.resolve()),
        "answer_path": str(answer_path.resolve()),
    }
    (output_dir / "metrics.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("compatibility_root", type=Path)
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("clinical_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    args = parser.parse_args()
    print(json.dumps(run(
        args.compatibility_root, args.patient_root, args.clinical_root,
        args.answers, args.output_dir,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

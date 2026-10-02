"""Build the patient-organism delta between frozen A visible and current F report.

This is a post-hoc evaluation utility.  It reads the already frozen outputs and
benchmark answers; it does not change either pipeline or any decision.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from tools.compare_aliasclean_baseline import load_sections
from tools.recalculate_kh_benchmark_metrics import match_type, pairwise_match, split_names, unique_names


ROOT = Path(__file__).resolve().parents[1]
OLD_ROOT = ROOT / "outputs" / "patient_info_KH_0728_2Days"
ANSWERS = (
    ROOT
    / "outputs"
    / "runs"
    / "2026-09-18_KH_answer_revision_metrics"
    / "kh_answers_clinical_revision_20260918.csv"
)
NEW_COMPLETE = (
    ROOT
    / "outputs"
    / "runs"
    / "2026-10-01_KH_ablation_F_current_full_v1"
    / "10_reporting_v3_route_a"
    / "complete_report.csv"
)
NEW_PATIENT_OUTPUT_ROOT = (
    ROOT
    / "outputs"
    / "runs"
    / "2026-10-01_KH_ablation_F_current_full_v1"
    / "10_reporting_v3_route_a"
    / "patient_outputs"
)
OUTPUT_DIR = (
    ROOT
    / "outputs"
    / "runs"
    / "2026-10-01_KH_old_A_visible_vs_new_F_complete_delta_v1"
)
OLD_SUFFIX = "mNGS_max_merged_selected_dna_20260918"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def old_visible(patient_id: int) -> list[str]:
    patient_dir = OLD_ROOT / f"NGS_patient_{patient_id}_json"
    sections = load_sections(patient_dir, OLD_SUFFIX)
    if sections is None:
        raise FileNotFoundError(f"Missing frozen A output for patient {patient_id}")
    names = []
    for section in ("picked", "possible", "watch"):
        names.extend(
            str(item.get("organism_name") or item.get("name") or "").strip()
            for item in sections[section]
        )
    return unique_names(name for name in names if name)


def old_output_path(patient_id: int) -> Path:
    return (
        OLD_ROOT
        / f"NGS_patient_{patient_id}_json"
        / "summary_outputs"
        / f"NGS_patient_{patient_id}_{OLD_SUFFIX}.json"
    )


def old_tiered_display(patient_id: int) -> str:
    sections = load_sections(OLD_ROOT / f"NGS_patient_{patient_id}_json", OLD_SUFFIX)
    if sections is None:
        raise FileNotFoundError(f"Missing frozen A output for patient {patient_id}")

    def names(section: str) -> str:
        values = unique_names(
            str(item.get("organism_name") or item.get("name") or "").strip()
            for item in sections[section]
            if str(item.get("organism_name") or item.get("name") or "").strip()
        )
        return "；".join(values) if values else "—"

    return "\n".join(
        (
            f"Picked：{names('picked')}",
            f"Luna High：{names('possible')}",
            f"Luna Context：{names('watch')}",
        )
    )


def new_complete_by_patient() -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in read_csv(NEW_COMPLETE):
        grouped[int(row["patient_id"])].append(row["organism_name"])
    return {patient: unique_names(names) for patient, names in grouped.items()}


def new_rows_by_patient() -> dict[int, list[dict[str, str]]]:
    grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in read_csv(NEW_COMPLETE):
        grouped[int(row["patient_id"])].append(row)
    return dict(grouped)


def new_tiered_display(rows: list[dict[str, str]]) -> str:
    def names(tier: str) -> str:
        values = unique_names(
            row["organism_name"]
            for row in rows
            if row.get("final_reporting_tier") == tier
        )
        return "；".join(values) if values else "—"

    return "\n".join(
        (
            f"Picked：{names('Picked')}",
            f"Possible：{names('Possible')}",
            f"Fallback-Possible：{names('Fallback-Possible')}",
        )
    )


def has_match(outputs: list[str], answer: str) -> bool:
    return any(match_type(output, answer, genus_relaxed=False) for output in outputs)


def subtract_names(left: list[str], right: list[str]) -> list[str]:
    """Return left-only names, treating aliases/group members as equivalent."""
    unmatched_left = set(range(len(left)))
    unmatched_right = set(range(len(right)))
    for desired in ("exact_or_alias", "approved_group_member_match"):
        for left_index in sorted(unmatched_left):
            right_index = next(
                (
                    idx
                    for idx in sorted(unmatched_right)
                    if match_type(left[left_index], right[idx], genus_relaxed=False)
                    == desired
                ),
                None,
            )
            if right_index is not None:
                unmatched_left.remove(left_index)
                unmatched_right.remove(right_index)
    return [left[index] for index in sorted(unmatched_left)]


def main() -> None:
    answer_rows = read_csv(ANSWERS)
    answers = {
        int(row["patient_id"]): split_names(row.get("answers") or "")
        for row in answer_rows
        if split_names(row.get("answers") or "")
    }
    current = new_complete_by_patient()
    current_rows = new_rows_by_patient()

    answer_deltas: list[dict[str, Any]] = []
    prediction_deltas: list[dict[str, Any]] = []
    old_tp = old_fp = new_tp = new_fp = 0

    for patient_id, gold in sorted(answers.items()):
        old_outputs = old_visible(patient_id)
        new_outputs = current.get(patient_id, [])
        old_eval = pairwise_match(old_outputs, gold, genus_relaxed=False)
        new_eval = pairwise_match(new_outputs, gold, genus_relaxed=False)
        old_tp += len(old_eval["matched"])
        old_fp += len(old_eval["output_only"])
        new_tp += len(new_eval["matched"])
        new_fp += len(new_eval["output_only"])

        for answer in gold:
            in_old = has_match(old_outputs, answer)
            in_new = has_match(new_outputs, answer)
            if in_old == in_new:
                status = "retained_tp" if in_old else "missed_by_both"
            elif in_new:
                status = "tp_gained_in_new"
            else:
                status = "tp_lost_in_new"
            answer_deltas.append(
                {
                    "patient_id": patient_id,
                    "answer": answer,
                    "old_visible_hit": in_old,
                    "new_complete_hit": in_new,
                    "status": status,
                }
            )

        old_only = subtract_names(old_outputs, new_outputs)
        new_only = subtract_names(new_outputs, old_outputs)
        for name in old_only:
            matching_answer = next(
                (answer for answer in gold if has_match([name], answer)), None
            )
            if matching_answer and has_match(new_outputs, matching_answer):
                role = "redundant_answer_equivalent_removed"
            elif matching_answer:
                role = "TP_lost"
            else:
                role = "FP_removed"
            prediction_deltas.append(
                {
                    "patient_id": patient_id,
                    "organism_name": name,
                    "change": "removed_from_new_report",
                    "benchmark_role": role,
                }
            )
        for name in new_only:
            matching_answer = next(
                (answer for answer in gold if has_match([name], answer)), None
            )
            if matching_answer and has_match(old_outputs, matching_answer):
                role = "redundant_answer_equivalent_added"
            elif matching_answer:
                role = "TP_gained"
            else:
                role = "FP_added"
            prediction_deltas.append(
                {
                    "patient_id": patient_id,
                    "organism_name": name,
                    "change": "added_to_new_report",
                    "benchmark_role": role,
                }
            )

    summary = {
        "scope": "30 labeled patients / 60 answer organisms; post-hoc comparison only",
        "matching": "exact or approved alias/group-member; no genus-relaxed matching",
        "old_endpoint": "frozen A Picked + Luna High + Luna Context",
        "new_endpoint": "current F Picked + Possible + mutually exclusive Fallback-Possible",
        "old": {"tp": old_tp, "fp": old_fp, "predictions": old_tp + old_fp},
        "new": {"tp": new_tp, "fp": new_fp, "predictions": new_tp + new_fp},
        "answer_delta_counts": {
            status: sum(1 for row in answer_deltas if row["status"] == status)
            for status in ("retained_tp", "tp_gained_in_new", "tp_lost_in_new", "missed_by_both")
        },
        "prediction_delta_counts": {
            label: sum(1 for row in prediction_deltas if row["change"] == label)
            for label in ("removed_from_new_report", "added_to_new_report")
        },
    }

    patient_comparison = []
    for row in answer_rows:
        patient_id = int(row["patient_id"])
        answer_values = split_names(row.get("answers") or "")
        patient_comparison.append(
            {
                "patient_id": f"P{patient_id}",
                "old_output": old_tiered_display(patient_id),
                "new_output": new_tiered_display(current_rows.get(patient_id, [])),
                "answers": "\n".join(answer_values) if answer_values else "—（無 benchmark 答案）",
                "old_final_source": str(old_output_path(patient_id).resolve()),
                "new_final_source": str(
                    (
                        NEW_PATIENT_OUTPUT_ROOT
                        / f"NGS_patient_{patient_id}_test_aware_possible_pathogen_shadow.json"
                    ).resolve()
                ),
            }
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(
        OUTPUT_DIR / "answer_level_delta.csv",
        answer_deltas,
        ["patient_id", "answer", "old_visible_hit", "new_complete_hit", "status"],
    )
    write_csv(
        OUTPUT_DIR / "prediction_level_delta.csv",
        prediction_deltas,
        ["patient_id", "organism_name", "change", "benchmark_role"],
    )
    (OUTPUT_DIR / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (OUTPUT_DIR / "patient_output_comparison.json").write_text(
        json.dumps(patient_comparison, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

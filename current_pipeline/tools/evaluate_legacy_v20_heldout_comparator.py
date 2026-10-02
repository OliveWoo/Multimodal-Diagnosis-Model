"""Evaluate the locked legacy-v20-compatible scorer on the 18 new KH patients.

This is deliberately separate from the historical A result.  The original A
pipeline depended on a manual selected-DNA input and an unavailable historical
execution environment, so it cannot be replayed bit-for-bit.  This evaluator
uses the retained legacy-v20-compatible scorer output produced from the same
per-test source used by the 51-patient release.

The primary evaluation keeps explicit no-pathogen patients in the precision
denominator.  A second, explicitly non-final pulmonary sensitivity removes
only the answer that the source workbook itself marks as urine-only.  Plasma
answers remain unresolved until blinded physician adjudication.
"""

from __future__ import annotations

import argparse
import copy
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

from tools.evaluate_external_kh_frozen_release import (
    TRUTHY,
    calculate_metrics,
    load_answer_contract,
    read_csv,
    read_json,
    sha256_file,
    verify_postfreeze_contract,
    write_csv,
)
from tools.recalculate_kh_benchmark_metrics import pairwise_match, unique_names


DEFAULT_LEGACY_ROOT = Path(
    "outputs/runs/2026-10-01_51patients_legacy_v20_compatible_v1"
)
DEFAULT_EVALUATION_ROOT = Path(
    "outputs/runs/2026-09-30_51patients_integrated_release_v2_6_heldout_freeze_v3_evaluation"
)
DEFAULT_ANSWER_PATH = DEFAULT_EVALUATION_ROOT / "answer_contract" / (
    "kh_18patient_heldout_answers_postfreeze_20260930.csv"
)
DEFAULT_ANSWER_MANIFEST_PATH = (
    DEFAULT_EVALUATION_ROOT / "answer_contract" / "answer_contract_manifest.json"
)
DEFAULT_NEW_METRICS_PATH = (
    DEFAULT_EVALUATION_ROOT / "all_labeled_v1" / "heldout" / "metrics.json"
)
DEFAULT_OUTPUT_DIR = Path(
    "outputs/runs/2026-10-01_18patients_legacy_v20_compatible_evaluation_v2"
)


def level_rank(value: Any) -> int:
    text = str(value or "")
    for rank in range(1, 6):
        if str(rank) in text:
            return rank
    return 99


def group_legacy_names(
    rows: list[dict[str, str]],
    patient_ids: set[int],
    predicate: Callable[[dict[str, str]], bool],
) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        patient_id = int(row["patient_id"])
        if patient_id in patient_ids and predicate(row):
            grouped[patient_id].append(row["organism_name"])
    return {patient_id: unique_names(names) for patient_id, names in grouped.items()}


def build_pulmonary_minimum_sensitivity(
    labels: dict[int, dict[str, Any]],
) -> tuple[dict[int, dict[str, Any]], list[dict[str, Any]]]:
    """Remove only explicit urine-only truth; do not adjudicate plasma answers."""

    sensitivity = copy.deepcopy(labels)
    changes: list[dict[str, Any]] = []
    patient = sensitivity.get(31)
    if patient is None:
        raise ValueError("Expected P31 in the 18-patient answer contract")
    raw_value = str(patient["row"].get("raw_workbook_values", ""))
    if "urine" not in raw_value.casefold():
        raise ValueError("P31 is no longer explicitly marked urine-only in source truth")
    removed = list(patient["answers"])
    patient["answers"] = []
    patient["status"] = "explicit_no_pathogen"
    changes.append({
        "patient_id": 31,
        "action": "exclude_explicit_nonpulmonary_answer",
        "removed_answers": "; ".join(removed),
        "source_text": raw_value,
        "reason": "Source workbook explicitly marks the diagnosis as urine-only.",
    })
    return sensitivity, changes


def collect_audit_rows(
    endpoint_predictions: dict[str, dict[int, list[str]]],
    labels: dict[int, dict[str, Any]],
    patient_ids: set[int],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    prediction_audit: list[dict[str, Any]] = []
    missed_answers: list[dict[str, Any]] = []
    negative_predictions: list[dict[str, Any]] = []
    for endpoint, grouped in endpoint_predictions.items():
        for patient_id in sorted(patient_ids):
            label = labels[patient_id]
            audit = pairwise_match(
                grouped.get(patient_id, []), label["answers"], genus_relaxed=False
            )
            matched_by_output = {item["output"]: item for item in audit["matched"]}
            for organism in grouped.get(patient_id, []):
                matched_item = matched_by_output.get(organism)
                prediction_audit.append({
                    "endpoint": endpoint,
                    "patient_id": patient_id,
                    "label_status": label["status"],
                    "organism_name": organism,
                    "outcome": "matched_answer" if matched_item else "false_positive",
                    "matched_answer": matched_item["answer"] if matched_item else "",
                    "match_type": matched_item["match_type"] if matched_item else "",
                })
                if label["status"] == "explicit_no_pathogen":
                    negative_predictions.append({
                        "endpoint": endpoint,
                        "patient_id": patient_id,
                        "organism_name": organism,
                    })
            for answer in audit["answer_only"]:
                missed_answers.append({
                    "endpoint": endpoint,
                    "patient_id": patient_id,
                    "answer": answer,
                })
    return prediction_audit, missed_answers, negative_predictions


def patient_summary_rows(
    predictions: dict[int, list[str]], labels: dict[int, dict[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for patient_id in sorted(labels):
        output = predictions.get(patient_id, [])
        audit = pairwise_match(output, labels[patient_id]["answers"], genus_relaxed=False)
        rows.append({
            "patient_id": patient_id,
            "label_status": labels[patient_id]["status"],
            "answers": "; ".join(labels[patient_id]["answers"]),
            "legacy_picked": "; ".join(output),
            "matched": "; ".join(item["output"] for item in audit["matched"]),
            "false_positive": "; ".join(audit["output_only"]),
            "missed_answer": "; ".join(audit["answer_only"]),
        })
    return rows


def run(
    legacy_root: Path,
    answer_path: Path,
    answer_manifest_path: Path,
    new_metrics_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty directory: {output_dir}")

    answer_manifest = verify_postfreeze_contract(answer_path, answer_manifest_path)
    labels = load_answer_contract(answer_path)
    patient_ids = set(labels)
    decision_path = legacy_root / "legacy_multi_assay_candidate_decisions.csv"
    legacy_summary_path = legacy_root / "summary.json"
    rows = read_csv(decision_path)
    available_ids = {int(row["patient_id"]) for row in rows}
    missing = patient_ids - available_ids
    if missing:
        raise ValueError(f"Legacy-compatible decisions miss patients: {sorted(missing)}")

    endpoints = {
        "legacy_picked": lambda row: row.get("picked", "").casefold() in TRUTHY,
        "legacy_level1_3_visible": lambda row: level_rank(row.get("integrated_level")) <= 3,
        "legacy_all_scored": lambda row: True,
    }
    predictions = {
        endpoint: group_legacy_names(rows, patient_ids, predicate)
        for endpoint, predicate in endpoints.items()
    }
    primary_metrics: dict[str, Any] = {}
    for endpoint, grouped in predictions.items():
        primary_metrics[endpoint], _ = calculate_metrics(grouped, labels, patient_ids)

    pulmonary_labels, pulmonary_changes = build_pulmonary_minimum_sensitivity(labels)
    pulmonary_metrics: dict[str, Any] = {}
    for endpoint, grouped in predictions.items():
        pulmonary_metrics[endpoint], _ = calculate_metrics(
            grouped, pulmonary_labels, patient_ids
        )

    prediction_audit, missed_answers, negative_predictions = collect_audit_rows(
        predictions, labels, patient_ids
    )
    patient_rows = patient_summary_rows(predictions["legacy_picked"], labels)

    new_metrics_document = read_json(new_metrics_path)
    new_metrics = new_metrics_document[
        "all_labeled_patients_including_explicit_no_pathogen"
    ]
    new_decision_path = Path(new_metrics_document["decision_path"])
    new_rows = read_csv(new_decision_path)
    new_predictions = {
        "strict_picked": group_legacy_names(
            new_rows,
            patient_ids,
            lambda row: row.get("possible_reporting_role") == "strict_picked_report",
        ),
        "complete_report": group_legacy_names(
            new_rows,
            patient_ids,
            lambda row: row.get("selected_for_complete_report", "").casefold() in TRUTHY,
        ),
    }
    new_pulmonary_metrics: dict[str, Any] = {}
    for endpoint, grouped in new_predictions.items():
        calculated, _ = calculate_metrics(grouped, pulmonary_labels, patient_ids)
        new_pulmonary_metrics[endpoint] = calculated
        primary_recalculated, _ = calculate_metrics(grouped, labels, patient_ids)
        if primary_recalculated != new_metrics[endpoint]:
            raise ValueError(
                f"Recalculated v2.6 {endpoint} metrics do not match the frozen evaluation"
            )
    side_by_side = [
        {
            "comparison_role": "formal_strict_output",
            "pipeline": "legacy_v20_compatible",
            "endpoint": "legacy_picked",
            **primary_metrics["legacy_picked"],
        },
        {
            "comparison_role": "formal_strict_output",
            "pipeline": "integrated_v2_6",
            "endpoint": "strict_picked",
            **new_metrics["strict_picked"],
        },
        {
            "comparison_role": "exploratory_non_equivalent_broad_output",
            "pipeline": "legacy_v20_compatible",
            "endpoint": "legacy_level1_3_visible",
            **primary_metrics["legacy_level1_3_visible"],
        },
        {
            "comparison_role": "exploratory_non_equivalent_broad_output",
            "pipeline": "integrated_v2_6",
            "endpoint": "complete_report",
            **new_metrics["complete_report"],
        },
    ]

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        output_dir / "prediction_audit.csv",
        [
            "endpoint", "patient_id", "label_status", "organism_name",
            "outcome", "matched_answer", "match_type",
        ],
        prediction_audit,
    )
    write_csv(
        output_dir / "missed_answers.csv",
        ["endpoint", "patient_id", "answer"],
        missed_answers,
    )
    write_csv(
        output_dir / "explicit_no_pathogen_predictions.csv",
        ["endpoint", "patient_id", "organism_name"],
        negative_predictions,
    )
    write_csv(
        output_dir / "patient_legacy_picked_summary.csv",
        [
            "patient_id", "label_status", "answers", "legacy_picked", "matched",
            "false_positive", "missed_answer",
        ],
        patient_rows,
    )
    write_csv(
        output_dir / "side_by_side_metrics.csv",
        [
            "comparison_role", "pipeline", "endpoint", "true_positive",
            "false_positive", "false_negative", "predicted", "answer_count",
            "precision", "recall", "f1", "patient_count", "positive_patient_count",
            "explicit_no_pathogen_patient_count",
            "explicit_no_pathogen_patients_with_predictions",
        ],
        side_by_side,
    )
    write_csv(
        output_dir / "pulmonary_endpoint_minimum_changes.csv",
        ["patient_id", "action", "removed_answers", "source_text", "reason"],
        pulmonary_changes,
    )
    pulmonary_side_by_side = [
        {
            "pipeline": "legacy_v20_compatible",
            "endpoint": "legacy_picked",
            **pulmonary_metrics["legacy_picked"],
        },
        {
            "pipeline": "integrated_v2_6",
            "endpoint": "strict_picked",
            **new_pulmonary_metrics["strict_picked"],
        },
        {
            "pipeline": "integrated_v2_6",
            "endpoint": "complete_report",
            **new_pulmonary_metrics["complete_report"],
        },
    ]
    write_csv(
        output_dir / "pulmonary_minimum_sensitivity_metrics.csv",
        [
            "pipeline", "endpoint", "true_positive", "false_positive",
            "false_negative", "predicted", "answer_count", "precision", "recall",
            "f1", "patient_count", "positive_patient_count",
            "explicit_no_pathogen_patient_count",
            "explicit_no_pathogen_patients_with_predictions",
        ],
        pulmonary_side_by_side,
    )
    write_csv(
        output_dir / "pulmonary_endpoint_adjudication_queue.csv",
        ["patient_id", "source_text", "status", "required_action"],
        [
            {
                "patient_id": 42,
                "source_text": labels[42]["row"].get("raw_workbook_values", ""),
                "status": "unresolved_plasma_linkage",
                "required_action": "Blinded physician decides whether each plasma organism belongs to the index pneumonia episode.",
            },
            {
                "patient_id": 44,
                "source_text": labels[44]["row"].get("raw_workbook_values", ""),
                "status": "unresolved_plasma_linkage",
                "required_action": "Blinded physician decides whether the plasma organism belongs to the index pneumonia episode.",
            },
        ],
    )

    result = {
        "schema_version": "kh_legacy_v20_heldout_comparator.v1",
        "scope": (
            "Post-unblinding locked-policy algorithmic replay on the 18 held-out patients. "
            "The scorer did not read answers, but execution was not frozen before unblinding."
        ),
        "historical_A_warning": (
            "This is the closest rerunnable legacy-v20-compatible comparator, not a "
            "bit-for-bit replay of historical A. Historical A used manually selected DNA "
            "inputs and an unavailable execution environment."
        ),
        "primary_truth_contract": (
            "All clinical pathogens from the physician workbook, including explicit "
            "no-pathogen patients in the precision denominator."
        ),
        "matching": (
            "Exact name, approved alias, or approved group-member match; no genus-relaxed matching."
        ),
        "legacy_root": str(legacy_root.resolve()),
        "legacy_decision_path": str(decision_path.resolve()),
        "legacy_decision_sha256": sha256_file(decision_path),
        "legacy_summary_path": str(legacy_summary_path.resolve()),
        "legacy_summary_sha256": sha256_file(legacy_summary_path),
        "answer_path": str(answer_path.resolve()),
        "answer_sha256": sha256_file(answer_path),
        "answer_manifest_path": str(answer_manifest_path.resolve()),
        "answer_manifest_decision_freeze_sha256": (
            answer_manifest.get("source_contracts") or {}
        ).get("decision_freeze_sha256"),
        "legacy_all_clinical_pathogen_metrics": primary_metrics,
        "pulmonary_minimum_sensitivity": {
            "status": "not_final_primary_endpoint",
            "definition": (
                "Removes only P31 Enterococcus faecium because the source explicitly says urine; "
                "P42/P44 plasma labels remain included pending blinded physician adjudication."
            ),
            "legacy_metrics": pulmonary_metrics,
            "new_v2_6_metrics": new_pulmonary_metrics,
        },
        "new_v2_6_all_clinical_pathogen_metrics": {
            "strict_picked": new_metrics["strict_picked"],
            "complete_report": new_metrics["complete_report"],
        },
        "comparison_cautions": [
            "Legacy Picked versus v2.6 strict Picked is the formal strict-output comparison.",
            "Legacy Level 1-3 visible versus v2.6 complete report is exploratory because the report roles are not semantically identical.",
            "Do not tune either pipeline on these 18 labels and then reuse the same metrics as external validation.",
        ],
    }
    (output_dir / "metrics.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-root", type=Path, default=DEFAULT_LEGACY_ROOT)
    parser.add_argument("--answer-path", type=Path, default=DEFAULT_ANSWER_PATH)
    parser.add_argument(
        "--answer-manifest-path", type=Path, default=DEFAULT_ANSWER_MANIFEST_PATH
    )
    parser.add_argument("--new-metrics-path", type=Path, default=DEFAULT_NEW_METRICS_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    print(json.dumps(run(
        args.legacy_root,
        args.answer_path,
        args.answer_manifest_path,
        args.new_metrics_path,
        args.output_dir,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

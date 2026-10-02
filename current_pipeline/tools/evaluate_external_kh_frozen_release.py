"""Evaluate a frozen KH release after loading the post-freeze answer contract.

This evaluator keeps explicit no-pathogen patients in the precision denominator.
It also reports the legacy positive-patient-only view for comparison with earlier
KH evaluations that silently omitted rows whose answer was ``-``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

from tools.build_test_aware_possible_pathogen_shadow import (
    ANALYTICAL_DIRECT_ROLE,
    FALLBACK_ROLE,
    STRICT_ROLE,
)
from tools.recalculate_kh_benchmark_metrics import pairwise_match, split_names, unique_names


TRUTHY = {"1", "true", "yes"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_answer_contract(path: Path) -> dict[int, dict[str, Any]]:
    result: dict[int, dict[str, Any]] = {}
    for row in read_csv(path):
        patient_id = int(row["patient_id"])
        if patient_id in result:
            raise ValueError(f"Duplicate answer-contract patient: P{patient_id}")
        status = row.get("label_status", "").strip()
        answers = unique_names(split_names(row.get("answers", "")))
        if status not in {"positive", "explicit_no_pathogen"}:
            raise ValueError(f"Unsupported label_status for P{patient_id}: {status!r}")
        if status == "positive" and not answers:
            raise ValueError(f"Positive patient P{patient_id} has no answer organisms")
        if status == "explicit_no_pathogen" and answers:
            raise ValueError(f"No-pathogen patient P{patient_id} has answer organisms")
        result[patient_id] = {"status": status, "answers": answers, "row": row}
    if not result:
        raise ValueError("Answer contract is empty")
    return result


def verify_postfreeze_contract(answer_path: Path, manifest_path: Path) -> dict[str, Any]:
    manifest = read_json(manifest_path)
    if manifest.get("scope") != "post-freeze evaluation only":
        raise ValueError("Answer manifest is not marked post-freeze evaluation only")
    if manifest.get("freeze_passed_before_unblinding") is not True:
        raise ValueError("Answer manifest does not confirm a pre-unblinding freeze")
    if manifest.get("answers_loaded_into_decision_pipeline") is not False:
        raise ValueError("Answer manifest permits answers in the decision pipeline")
    output_contracts = manifest.get("outputs") or {}
    expected_hashes = {
        str(value).casefold()
        for key, value in output_contracts.items()
        if str(key).endswith("_sha256") and value
    }
    actual_hash = sha256_file(answer_path).casefold()
    if actual_hash not in expected_hashes:
        raise ValueError("Answer CSV hash is not registered in the answer manifest")
    return manifest


def group_names(
    rows: list[dict[str, str]],
    allowed_ids: set[int],
    predicate: Callable[[dict[str, str]], bool],
) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        patient_id = int(row["patient_id"])
        if patient_id in allowed_ids and predicate(row):
            grouped[patient_id].append(row["organism_name"])
    return {patient_id: unique_names(names) for patient_id, names in grouped.items()}


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
    negative_with_predictions = sorted(
        patient_id for patient_id in negative_ids if predictions.get(patient_id)
    )
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
            labels[patient_id]["status"] == "positive" for patient_id in patient_ids
        ),
        "explicit_no_pathogen_patient_count": len(negative_ids),
        "explicit_no_pathogen_patients_with_predictions": negative_with_predictions,
    }, audits


def run(
    shadow_root: Path,
    answer_path: Path,
    answer_manifest_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty directory: {output_dir}")
    manifest = verify_postfreeze_contract(answer_path, answer_manifest_path)
    labels = load_answer_contract(answer_path)
    patient_ids = set(labels)
    positive_ids = {
        patient_id
        for patient_id, label in labels.items()
        if label["status"] == "positive"
    }
    decision_path = shadow_root / "possible_pathogen_decisions.csv"
    rows = read_csv(decision_path)
    available_ids = {int(row["patient_id"]) for row in rows}
    missing = patient_ids - available_ids
    if missing:
        raise ValueError(f"Frozen decisions are missing labeled patients: {sorted(missing)}")

    endpoints = {
        "strict_picked": lambda row: row["possible_reporting_role"] == STRICT_ROLE,
        "strict_plus_analytical_or_direct_possible": lambda row: (
            row["possible_reporting_role"] in {STRICT_ROLE, ANALYTICAL_DIRECT_ROLE}
        ),
        "picked_plus_all_possible": lambda row: (
            row["selected_for_complete_report"].casefold() in TRUTHY
            and row["possible_reporting_role"] != FALLBACK_ROLE
        ),
        "complete_report": lambda row: (
            row["selected_for_complete_report"].casefold() in TRUTHY
        ),
        "all_forwarded_visible": lambda row: True,
    }
    predictions = {
        endpoint: group_names(rows, patient_ids, predicate)
        for endpoint, predicate in endpoints.items()
    }

    all_metrics: dict[str, Any] = {}
    legacy_metrics: dict[str, Any] = {}
    all_audits: dict[str, dict[int, dict[str, Any]]] = {}
    for endpoint, grouped in predictions.items():
        all_metrics[endpoint], all_audits[endpoint] = calculate_metrics(
            grouped, labels, patient_ids
        )
        legacy_metrics[endpoint], _ = calculate_metrics(grouped, labels, positive_ids)

    prediction_audit: list[dict[str, Any]] = []
    missed_answers: list[dict[str, Any]] = []
    explicit_negative_predictions: list[dict[str, Any]] = []
    for endpoint, grouped in predictions.items():
        for patient_id in sorted(patient_ids):
            label = labels[patient_id]
            audit = all_audits[endpoint][patient_id]
            matched_by_output = {
                item["output"]: item for item in audit["matched"]
            }
            for organism in grouped.get(patient_id, []):
                matched_item = matched_by_output.get(organism)
                outcome = "matched_answer" if matched_item else "false_positive"
                prediction_audit.append({
                    "endpoint": endpoint,
                    "patient_id": patient_id,
                    "label_status": label["status"],
                    "organism_name": organism,
                    "outcome": outcome,
                    "matched_answer": matched_item["answer"] if matched_item else "",
                    "match_type": matched_item["match_type"] if matched_item else "",
                })
                if label["status"] == "explicit_no_pathogen":
                    explicit_negative_predictions.append({
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
        explicit_negative_predictions,
    )
    result = {
        "schema_version": "kh_external_frozen_release_evaluation.v1",
        "scope": "Post-freeze evaluation only; answers were not available to the decision pipeline.",
        "matching": "Exact name, approved alias, or approved group member; no genus-relaxed matching.",
        "shadow_root": str(shadow_root.resolve()),
        "decision_path": str(decision_path.resolve()),
        "decision_sha256": sha256_file(decision_path),
        "answer_path": str(answer_path.resolve()),
        "answer_sha256": sha256_file(answer_path),
        "answer_manifest_path": str(answer_manifest_path.resolve()),
        "decision_freeze_sha256": (
            (manifest.get("source_contracts") or {}).get("decision_freeze_sha256")
        ),
        "all_labeled_patients_including_explicit_no_pathogen": all_metrics,
        "legacy_positive_patients_only": legacy_metrics,
        "notes": [
            "Use all_labeled_patients_including_explicit_no_pathogen for the primary precision estimate.",
            "legacy_positive_patients_only is retained only for comparison with earlier KH reports.",
        ],
    }
    (output_dir / "metrics.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shadow_root", type=Path)
    parser.add_argument("answer_path", type=Path)
    parser.add_argument("answer_manifest_path", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(
        run(
            args.shadow_root,
            args.answer_path,
            args.answer_manifest_path,
            args.output_dir,
        ),
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()

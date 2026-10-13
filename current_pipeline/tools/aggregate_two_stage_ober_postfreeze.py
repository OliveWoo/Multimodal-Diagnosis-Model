"""Aggregate disjoint cohort-level two-stage OBER post-freeze evaluations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from tools.evaluate_test_aware_possible_pathogen_shadow import calculate_metrics
from tools.recalculate_kh_benchmark_metrics import match_type, split_names, unique_names


ENDPOINT_COLUMNS = {
    "deterministic_strict_picked": ("deterministic_picked",),
    "deterministic_picked_plus_mandatory_fallback": (
        "deterministic_picked",
        "mandatory_fallback_addition",
    ),
    "strict_picked_plus_ober_accepted": (
        "deterministic_picked",
        "ober_accepted_addition",
    ),
    "complete_legacy_with_fallback_and_ober": ("complete_legacy_output",),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def merge_names(*values: str) -> list[str]:
    output: list[str] = []
    for value in values:
        for name in split_names(value):
            if not any(match_type(existing, name, genus_relaxed=False) for existing in output):
                output.append(name)
    return unique_names(output)


def parse_source(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Expected COHORT=PATH")
    cohort, raw_path = value.split("=", 1)
    if not cohort.strip() or not raw_path.strip():
        raise argparse.ArgumentTypeError("Expected non-empty COHORT=PATH")
    return cohort.strip(), Path(raw_path.strip())


def run(sources: list[tuple[str, Path]], output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_dir}")
    source_summaries: dict[str, dict[str, Any]] = {}
    rows: list[dict[str, str]] = []
    source_artifacts: dict[str, Any] = {}
    seen_patients: dict[int, str] = {}
    for cohort, root in sources:
        summary_path = root / "summary.json"
        patient_path = root / "patient_complete_legacy_output.csv"
        summary = json.loads(summary_path.read_text(encoding="utf-8-sig"))
        if summary.get("status") != "complete_postfreeze_benchmark":
            raise ValueError(f"Incomplete cohort evaluation: {cohort}")
        if summary.get("cohort") != cohort:
            raise ValueError(f"Cohort name mismatch for {cohort}: {summary.get('cohort')}")
        cohort_rows = read_csv(patient_path)
        if len(cohort_rows) != int(summary["patient_count"]):
            raise ValueError(f"Patient count mismatch for {cohort}")
        for row in cohort_rows:
            patient_id = int(row["patient_id"])
            if patient_id in seen_patients:
                raise ValueError(
                    f"Patient P{patient_id} overlaps {seen_patients[patient_id]} and {cohort}"
                )
            seen_patients[patient_id] = cohort
            if row.get("cohort") != cohort:
                raise ValueError(f"Patient row cohort mismatch for P{patient_id}")
            rows.append(row)
        source_summaries[cohort] = summary
        source_artifacts[cohort] = {
            "root": str(root.resolve()),
            "summary_sha256": sha256_file(summary_path),
            "patient_output_sha256": sha256_file(patient_path),
        }

    labels: dict[int, dict[str, Any]] = {}
    predictions: dict[str, dict[int, list[str]]] = {
        endpoint: {} for endpoint in ENDPOINT_COLUMNS
    }
    cohort_patients: dict[str, set[int]] = {}
    for row in rows:
        patient_id = int(row["patient_id"])
        cohort = row["cohort"]
        status = row["label_status"]
        answers = merge_names(row.get("reference_answers", ""))
        if status not in {"positive", "explicit_no_pathogen"}:
            raise ValueError(f"Unsupported label status for P{patient_id}: {status}")
        if status == "positive" and not answers:
            raise ValueError(f"Positive patient P{patient_id} has no answers")
        if status == "explicit_no_pathogen" and answers:
            raise ValueError(f"No-pathogen patient P{patient_id} has answers")
        labels[patient_id] = {"status": status, "answers": answers, "cohort": cohort}
        cohort_patients.setdefault(cohort, set()).add(patient_id)
        for endpoint, columns in ENDPOINT_COLUMNS.items():
            predictions[endpoint][patient_id] = merge_names(
                *(row.get(column, "") for column in columns)
            )

    required = {"KH33", "KH18", "two_hospital_A11", "two_hospital_B11"}
    if set(cohort_patients) != required:
        raise ValueError(
            f"Expected cohorts {sorted(required)}, found {sorted(cohort_patients)}"
        )
    strata = {
        "KH33": cohort_patients["KH33"],
        "KH18": cohort_patients["KH18"],
        "two_hospital_A11": cohort_patients["two_hospital_A11"],
        "two_hospital_B11": cohort_patients["two_hospital_B11"],
        "KH_complete_51": cohort_patients["KH33"] | cohort_patients["KH18"],
        "two_hospital_complete_22": (
            cohort_patients["two_hospital_A11"]
            | cohort_patients["two_hospital_B11"]
        ),
        "original_55": (
            cohort_patients["KH33"]
            | cohort_patients["two_hospital_A11"]
            | cohort_patients["two_hospital_B11"]
        ),
        "original_44_without_B11": (
            cohort_patients["KH33"] | cohort_patients["two_hospital_A11"]
        ),
        "all_62_without_B11": (
            cohort_patients["KH33"]
            | cohort_patients["KH18"]
            | cohort_patients["two_hospital_A11"]
        ),
        "all_73": set(labels),
    }

    metrics: dict[str, dict[str, Any]] = {}
    metric_rows: list[dict[str, Any]] = []
    audits: dict[str, dict[str, dict[int, dict[str, Any]]]] = {}
    for stratum, patient_ids in strata.items():
        metrics[stratum] = {}
        audits[stratum] = {}
        for endpoint, endpoint_predictions in predictions.items():
            values, audit = calculate_metrics(endpoint_predictions, labels, patient_ids)
            metrics[stratum][endpoint] = values
            audits[stratum][endpoint] = audit
            metric_rows.append({"stratum": stratum, "endpoint": endpoint, **values})

    for cohort, summary in source_summaries.items():
        for endpoint in ENDPOINT_COLUMNS:
            direct = summary["metrics"][endpoint]
            combined = metrics[cohort][endpoint]
            for key in ("true_positive", "false_positive", "false_negative", "predicted", "answer_count"):
                if direct[key] != combined[key]:
                    raise ValueError(
                        f"Aggregate mismatch for {cohort}/{endpoint}/{key}: "
                        f"{direct[key]} != {combined[key]}"
                    )

    complete_endpoint = "complete_legacy_with_fallback_and_ober"
    complete_audit = audits["all_73"][complete_endpoint]
    prediction_rows: list[dict[str, Any]] = []
    missed_rows: list[dict[str, Any]] = []
    for patient_id in sorted(labels):
        audit = complete_audit[patient_id]
        matched_by_output = {item["output"]: item for item in audit["matched"]}
        for organism in predictions[complete_endpoint].get(patient_id, []):
            match = matched_by_output.get(organism)
            prediction_rows.append(
                {
                    "patient_id": patient_id,
                    "cohort": labels[patient_id]["cohort"],
                    "label_status": labels[patient_id]["status"],
                    "organism_name": organism,
                    "outcome": "matched_answer" if match else "false_positive",
                    "matched_answer": match["answer"] if match else "",
                    "match_type": match["match_type"] if match else "",
                }
            )
        for answer in audit["answer_only"]:
            missed_rows.append(
                {
                    "patient_id": patient_id,
                    "cohort": labels[patient_id]["cohort"],
                    "answer": answer,
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "metrics_table.csv", metric_rows, list(metric_rows[0]))
    write_csv(
        output_dir / "patient_complete_legacy_output.csv",
        sorted(rows, key=lambda row: int(row["patient_id"])),
        list(rows[0]),
    )
    write_csv(
        output_dir / "complete_endpoint_prediction_audit.csv",
        prediction_rows,
        list(prediction_rows[0]),
    )
    write_csv(
        output_dir / "complete_endpoint_missed_answers.csv",
        missed_rows,
        ["patient_id", "cohort", "answer"],
    )
    result = {
        "schema_version": "legacy_two_stage_ober_aggregate_postfreeze.v1",
        "status": "complete",
        "recommended_endpoint": complete_endpoint,
        "matching": (
            "Exact name, approved alias, or approved group member; "
            "no genus-relaxed matching."
        ),
        "patient_count": len(labels),
        "positive_patient_count": sum(
            label["status"] == "positive" for label in labels.values()
        ),
        "explicit_no_pathogen_patient_count": sum(
            label["status"] == "explicit_no_pathogen" for label in labels.values()
        ),
        "answer_count": sum(len(label["answers"]) for label in labels.values()),
        "source_artifacts": source_artifacts,
        "strata_patient_counts": {
            name: len(patient_ids) for name, patient_ids in strata.items()
        },
        "metrics": metrics,
        "complete_endpoint_false_positive_count": len(prediction_rows)
        - sum(row["outcome"] == "matched_answer" for row in prediction_rows),
        "complete_endpoint_missed_answer_count": len(missed_rows),
        "notes": [
            "All source OBER decisions were frozen before answer contracts were opened.",
            "Explicit no-pathogen patients remain in the precision denominator.",
            "OBER accepts only primary_pathogen and secondary_pathogen actions.",
            "Context-only and insufficient-evidence decisions are not formal positives.",
            "A11 and B11 are selected-DNA legacy compatibility replays, not untouched external validation.",
        ],
    }
    write_json(output_dir / "summary.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        action="append",
        type=parse_source,
        required=True,
        help="Cohort evaluation in COHORT=PATH form; repeat four times.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.source, args.output_dir)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

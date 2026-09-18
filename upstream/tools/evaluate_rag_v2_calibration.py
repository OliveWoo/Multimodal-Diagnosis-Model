"""Evaluate RAG v2 adjudication results against private calibration labels."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Sequence


SHOW_VISIBILITY = {"show_primary", "show_secondary", "show_context"}
DO_NOT_SHOW = {"do_not_show"}
KEEP_ACTIONS = {"primary_pathogen", "secondary_pathogen", "context_only"}
REJECT_ACTIONS = {"reject", "insufficient_evidence"}


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def label_rows(labels_path: Path) -> list[dict[str, Any]]:
    payload = read_json(labels_path)
    rows = payload.get("cases") if isinstance(payload, dict) else []
    return [row for row in rows if isinstance(row, dict)]


def adjudication_rows(adjudications_path: Path) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for record in read_jsonl(adjudications_path):
        case_id = str(record.get("case_id") or "")
        adjudication = record.get("adjudication")
        if case_id and isinstance(adjudication, dict):
            rows[case_id] = record
    return rows


def decision_fields(record: dict[str, Any] | None) -> tuple[str, str, str]:
    if not record:
        return "", "", ""
    adjudication = record.get("adjudication") if isinstance(record.get("adjudication"), dict) else {}
    decision = adjudication.get("decision") if isinstance(adjudication.get("decision"), dict) else {}
    return (
        str(decision.get("recommended_action") or ""),
        str(decision.get("clinician_visibility") or ""),
        str(decision.get("confidence") or ""),
    )


def is_clinician_visible(visibility: str) -> bool:
    return visibility in SHOW_VISIBILITY


def is_success(label_type: str, action: str, visibility: str) -> bool:
    if label_type == "positive_control":
        return is_clinician_visible(visibility) or action in KEEP_ACTIONS
    if label_type == "negative_control":
        return (visibility in DO_NOT_SHOW) or (action in REJECT_ACTIONS)
    return False


def evaluate(labels: Sequence[dict[str, Any]], adjudications: dict[str, dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    decision_counts: Counter[str] = Counter()
    category_counts: Counter[str] = Counter()
    category_success_counts: Counter[str] = Counter()

    for label in labels:
        case_id = str(label.get("case_id") or "")
        label_type = str(label.get("label_type") or "")
        category = str(label.get("pathogen_category") or "")
        record = adjudications.get(case_id)
        action, visibility, confidence = decision_fields(record)
        visible = is_clinician_visible(visibility)
        success = bool(record) and is_success(label_type, action, visibility)
        false_retention = label_type == "negative_control" and visible
        false_rejection = label_type == "positive_control" and not visible

        counts["total"] += 1
        counts[label_type] += 1
        counts[f"{label_type}:completed"] += 1 if record else 0
        counts[f"{label_type}:success"] += 1 if success else 0
        counts[f"{label_type}:visible"] += 1 if visible else 0
        counts["success"] += 1 if success else 0
        counts["completed"] += 1 if record else 0
        if false_retention:
            counts["negative_control:false_retention"] += 1
        if false_rejection:
            counts["positive_control:false_rejection"] += 1
        decision_counts[f"{label_type}|{action}|{visibility}"] += 1
        category_counts[f"{label_type}|{category}"] += 1
        if success:
            category_success_counts[f"{label_type}|{category}"] += 1

        rows.append(
            {
                "case_id": case_id,
                "label_type": label_type,
                "expected_rag_behavior": label.get("expected_rag_behavior"),
                "dataset": label.get("dataset"),
                "patient_id": label.get("patient_id"),
                "organism_name": label.get("organism_name"),
                "pathogen_category": category,
                "source_tier": label.get("source_tier"),
                "recommended_action": action,
                "clinician_visibility": visibility,
                "confidence": confidence,
                "clinician_visible": visible,
                "calibration_success": success,
                "false_retention_negative": false_retention,
                "false_rejection_positive": false_rejection,
                "patient_answers_private": label.get("patient_answers_private"),
            }
        )

    positive_total = counts["positive_control"]
    negative_total = counts["negative_control"]
    summary = {
        "schema_version": "rag_calibration_evaluation_v1.0",
        "case_count": counts["total"],
        "completed_count": counts["completed"],
        "success_count": counts["success"],
        "calibration_accuracy": round(counts["success"] / counts["total"], 4) if counts["total"] else None,
        "positive_controls": {
            "total": positive_total,
            "completed": counts["positive_control:completed"],
            "retained_or_visible": counts["positive_control:success"],
            "false_rejection": counts["positive_control:false_rejection"],
            "retention_rate": round(counts["positive_control:success"] / positive_total, 4) if positive_total else None,
        },
        "negative_controls": {
            "total": negative_total,
            "completed": counts["negative_control:completed"],
            "rejected_or_not_visible": counts["negative_control:success"],
            "false_retention": counts["negative_control:false_retention"],
            "rejection_rate": round(counts["negative_control:success"] / negative_total, 4) if negative_total else None,
            "false_retention_rate": round(counts["negative_control:false_retention"] / negative_total, 4) if negative_total else None,
        },
        "decision_counts": dict(sorted(decision_counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
        "category_success_counts": dict(sorted(category_success_counts.items())),
    }
    return rows, summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("labels", type=Path)
    parser.add_argument("adjudications", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    labels = label_rows(args.labels)
    adjudications = adjudication_rows(args.adjudications)
    rows, summary = evaluate(labels, adjudications)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json(args.output_dir / "calibration_evaluation_summary.json", summary)
    write_csv(
        args.output_dir / "calibration_evaluation_cases.csv",
        rows,
        (
            "case_id",
            "label_type",
            "expected_rag_behavior",
            "dataset",
            "patient_id",
            "organism_name",
            "pathogen_category",
            "source_tier",
            "recommended_action",
            "clinician_visibility",
            "confidence",
            "clinician_visible",
            "calibration_success",
            "false_retention_negative",
            "false_rejection_positive",
            "patient_answers_private",
        ),
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

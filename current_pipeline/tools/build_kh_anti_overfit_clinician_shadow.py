from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any


TIER_ORDER = {
    "High-confidence candidate": 0,
    "Provisional-Possible": 1,
    "Possible": 2,
    "Doctor-visible Fallback": 3,
    "Context": 4,
    "Audit": 5,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: Any) -> bool:
    return str(value).strip().casefold() in {"true", "1", "yes"}


def profile_row(row: dict[str, str], policy: dict[str, Any]) -> dict[str, Any]:
    current = row["final_reporting_tier"]
    promotion = row.get("promotion_route", "")
    family = row.get("taxonomy_family", "")
    reason = ""
    if family == "unmapped_or_uncertain":
        target = "Audit"
        doctor_visible = False
        reason = "taxonomy unresolved; audit only"
    elif current == "Fallback-Possible":
        target = "Doctor-visible Fallback"
        doctor_visible = True
        reason = "shown in a separate low-confidence supplement; excluded from the primary pulmonary-pathogen positive endpoint"
    elif current == "Picked" and promotion in set(policy["picked_routes_requiring_provisional_demotion"]):
        target = "Provisional-Possible"
        doctor_visible = True
        reason = "promotion depends on analytical-only, imaging-only, or group-to-member evidence"
    elif current == "Picked":
        target = "High-confidence candidate"
        doctor_visible = True
        reason = "preserved high-confidence candidate; still pending external validation"
    elif current == "Possible":
        target = "Possible"
        doctor_visible = True
        reason = "retained recall layer with route-specific caution"
    else:
        target = "Context"
        doctor_visible = False
        reason = "internal review only"
    return {
        "patient_id": row["patient_id"],
        "organism_name": row["organism_name"],
        "current_tier": current,
        "anti_overfit_tier": target,
        "doctor_visible": doctor_visible,
        "taxonomy_family": family,
        "promotion_route": promotion,
        "reporting_route": row.get("reporting_route", ""),
        "reason": reason,
        "benchmark_status_posthoc_only": row.get("benchmark_status", ""),
    }


def metric(rows: list[dict[str, Any]], accepted: set[str]) -> dict[str, Any]:
    selected = [row for row in rows if row["anti_overfit_tier"] in accepted]
    labeled = [row for row in selected if row["benchmark_status_posthoc_only"] in {"reported_TP", "reported_FP"}]
    tp = sum(row["benchmark_status_posthoc_only"] == "reported_TP" for row in labeled)
    fp = sum(row["benchmark_status_posthoc_only"] == "reported_FP" for row in labeled)
    answers = 60
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / answers
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "tp": tp,
        "fp": fp,
        "fn": answers - tp,
        "predicted_in_labeled_patients": tp + fp,
        "unlabeled_predictions": sum(row["benchmark_status_posthoc_only"] == "unlabeled_patient" for row in selected),
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def run(matrix: Path, policy_path: Path, endpoint_path: Path, output: Path) -> dict[str, Any]:
    source = read_csv(matrix)
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    endpoint = json.loads(endpoint_path.read_text(encoding="utf-8"))
    rows = [profile_row(row, policy) for row in source]
    rows.sort(key=lambda row: (int(row["patient_id"]), TIER_ORDER[row["anti_overfit_tier"]], row["organism_name"].casefold()))
    output.mkdir(parents=True, exist_ok=False)
    write_csv(output / "anti_overfit_patient_organism_profile.csv", rows)
    changes = [row for row in rows if row["current_tier"] not in {"Context"} and row["anti_overfit_tier"] not in {"High-confidence candidate", "Possible"}]
    write_csv(output / "tier_changes_for_review.csv", changes)
    counts = Counter(row["anti_overfit_tier"] for row in rows)
    summary = {
        "schema_version": "kh_anti_overfit_clinician_shadow.v1",
        "status": "research_shadow_not_production_validated",
        "primary_endpoint": endpoint["primary_question"],
        "source_rows": len(source),
        "tier_counts": dict(counts),
        "changed_noncontext_rows": len(changes),
        "development_posthoc_metrics": {
            "high_confidence_only": metric(rows, {"High-confidence candidate"}),
            "primary_pulmonary_report_excluding_fallback": metric(rows, {"High-confidence candidate", "Provisional-Possible", "Possible"}),
            "doctor_visible_including_fallback_descriptive_only": metric(rows, {"High-confidence candidate", "Provisional-Possible", "Possible", "Doctor-visible Fallback"}),
        },
        "warnings": [
            "Metrics are post-hoc descriptions of the reused 33-patient development cohort.",
            "This profile reduces unsupported confidence claims; it does not prove that overfitting has been removed.",
            "Formal replacement still requires a new untouched cohort under the pulmonary endpoint."
        ]
    }
    (output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build conservative anti-overfit clinician-output shadow.")
    parser.add_argument("--matrix", type=Path, default=Path("outputs/runs/2026-10-01_KH_organism_rule_review_v3/patient_organism_rule_matrix.csv"))
    parser.add_argument("--policy", type=Path, default=Path("rules/anti_overfit_clinician_output_profile_v1.json"))
    parser.add_argument("--endpoint", type=Path, default=Path("rules/pulmonary_pathogen_evaluation_target_v1.json"))
    parser.add_argument("--output", type=Path, default=Path("outputs/runs/2026-10-01_KH_anti_overfit_clinician_shadow_v1"))
    args = parser.parse_args()
    print(json.dumps(run(args.matrix, args.policy, args.endpoint, args.output), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

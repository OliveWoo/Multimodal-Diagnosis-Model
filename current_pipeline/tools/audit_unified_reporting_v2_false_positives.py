"""Build a post-hoc FP/TP contrast audit for unified reporting v2.

This tool is deliberately answer-aware and read-only with respect to scorer
decisions. It must never be imported into the answer-blind reporting builder.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0]) if rows else ["patient_id", "organism_name"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def patient_number(path: Path) -> int:
    found = re.search(r"NGS_patient_(\d+)_", path.name)
    if not found:
        raise ValueError(f"Cannot identify patient from {path}")
    return int(found.group(1))


def as_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def as_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def audit_focus(role: str, route: str) -> str:
    if role == "strict_picked_report":
        return "upstream_strict_scorer_or_promotion_audit"
    if role == "fallback_best_available_possible_pathogen":
        return "forced_fallback_keep_separate_from_primary_confidence"
    if route == "reproducible_top_ranked_typical_or_hospital_bacterium":
        return "analytical_only_possible_precision_guardrail"
    if route == "event_aligned_direct_detection_or_level3_culture":
        return "direct_evidence_identity_timing_and_quantity_audit"
    if route == "reproducible_cross_molecule_top_ranked_low_specificity_signal":
        return "source_or_colonization_prone_family_guardrail"
    return "candidate_specific_evidence_audit"


def candidate_row(
    candidate: dict[str, Any],
    *,
    outcome: str,
    matched_answer: str,
) -> dict[str, Any]:
    taxonomy = candidate.get("taxonomy_profile") or {}
    gate = candidate.get("possible_pathogen_gate") or {}
    axes = gate.get("axes") or {}
    promotion = candidate.get("promotion_gate") or {}
    promotion_axes = promotion.get("axes") or {}
    analytical = candidate.get("analytical_profile") or {}
    signals = [
        row
        for row in analytical.get("per_test_signals") or []
        if isinstance(row, dict) and float(row.get("reads") or 0) > 0
    ]
    qc_counts = Counter(str(row.get("qc_status") or "missing") for row in signals)
    nucleic_counts = Counter(
        str(row.get("nucleic_type") or "missing").upper() for row in signals
    )
    rpms = [
        value
        for value in (as_float(row.get("rpm_total")) for row in signals)
        if value is not None
    ]
    reads = [float(row.get("reads") or 0) for row in signals]
    colonization = candidate.get("colonization_interpretation") or {}
    reactivation = candidate.get("reactivation_interpretation") or {}
    role = str(candidate.get("possible_reporting_role") or "")
    route = str(gate.get("route") or "")
    return {
        "patient_id": candidate.get("patient_id"),
        "organism_name": candidate.get("organism_name"),
        "posthoc_outcome": outcome,
        "matched_answer": matched_answer,
        "final_reporting_tier": candidate.get("final_reporting_tier"),
        "possible_reporting_role": role,
        "possible_route": route,
        "audit_focus": audit_focus(role, route),
        "clinical_decision": candidate.get("clinical_decision"),
        "clinical_level": candidate.get("clinical_level"),
        "formal_pick_allowed": candidate.get("formal_pick_allowed"),
        "taxonomy_family": axes.get("taxonomy_family"),
        "taxonomy_mapping_status": taxonomy.get("mapping_status"),
        "taxonomy_confidence": taxonomy.get("classification_confidence"),
        "precise_identity": axes.get("precise_identity"),
        "specimen_context": axes.get("specimen_context"),
        "best_rank": axes.get("best_rank_in_retained_universe"),
        "selected_positive_test_count": axes.get("selected_positive_test_count"),
        "reproducibility_axis": axes.get("reproducibility_axis"),
        "cross_molecule_selected": axes.get("cross_molecule_selected"),
        "dna_positive_test_count": nucleic_counts.get("DNA", 0),
        "rna_positive_test_count": nucleic_counts.get("RNA", 0),
        "fully_evaluable_positive_test_count": qc_counts.get("evaluable", 0),
        "partial_evaluable_positive_test_count": qc_counts.get(
            "partial_evaluable", 0
        ),
        "maximum_per_test_reads": max(reads, default=0),
        "maximum_available_per_test_rpm": max(rpms, default=""),
        "direct_hospital_level": axes.get("direct_hospital_level"),
        "event_aligned_direct_positive_count": axes.get(
            "event_aligned_direct_positive_count"
        ),
        "repeated_sterile_positive_count": axes.get(
            "repeated_sterile_positive_count"
        ),
        "host_support": promotion_axes.get("host_support"),
        "host_support_eligible_for_promotion": promotion_axes.get(
            "host_support_eligible_for_promotion"
        ),
        "history_route_changed": candidate.get("history_route_changed"),
        "history_policy_family": candidate.get("history_policy_family"),
        "explicit_colonization": colonization.get("has_explicit_colonization"),
        "explicit_current_not_infection": colonization.get(
            "has_explicit_current_not_infection"
        ),
        "explicit_current_reactivation": reactivation.get(
            "has_explicit_current_reactivation"
        ),
        "fallback_primary": candidate.get("fallback_primary"),
        "final_reporting_confidence": candidate.get("final_reporting_confidence"),
        "possible_blockers": "|".join(gate.get("blockers") or []),
        "possible_cautions": "|".join(gate.get("cautions") or []),
        "promotion_blockers": "|".join(promotion.get("blockers") or []),
        "direct_evidence_rows_json": json.dumps(
            ((promotion_axes.get("direct_evidence_timing_profile") or {}).get(
                "rows"
            ) or []),
            ensure_ascii=False,
            sort_keys=True,
        ),
        "hospital_evidence_detail_json": json.dumps(
            candidate.get("hospital_evidence_detail") or {},
            ensure_ascii=False,
            sort_keys=True,
        ),
        "per_test_signals_json": json.dumps(
            signals, ensure_ascii=False, sort_keys=True
        ),
    }


def run(shadow_root: Path, evaluation_root: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    audit_path = evaluation_root / "prediction_audit.csv"
    audit_rows = read_csv(audit_path)
    labels = {
        (str(row["patient_id"]), row["organism_name"]): (
            row["posthoc_outcome"], row.get("matched_answer") or ""
        )
        for row in audit_rows
    }
    rows: list[dict[str, Any]] = []
    paths = sorted(
        (shadow_root / "patient_outputs").glob(
            "NGS_patient_*_test_aware_possible_pathogen_shadow.json"
        ),
        key=patient_number,
    )
    for path in paths:
        payload = read_json(path)
        for candidate in payload.get("complete_report") or []:
            key = (str(candidate.get("patient_id")), candidate.get("organism_name"))
            outcome, answer = labels.get(key, ("missing_evaluation_row", ""))
            rows.append(candidate_row(candidate, outcome=outcome, matched_answer=answer))

    false_positives = [
        row for row in rows
        if row["posthoc_outcome"] == "unmatched_in_labeled_patient"
    ]
    true_positives = [
        row for row in rows if row["posthoc_outcome"] == "matched_benchmark_answer"
    ]
    unlabeled = [row for row in rows if row["posthoc_outcome"] == "unlabeled_patient"]

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "false_positive_audit.csv", false_positives)
    write_csv(output_dir / "true_positive_contrast.csv", true_positives)
    write_csv(output_dir / "unlabeled_report_candidates.csv", unlabeled)

    strict_fp = sum(
        row["possible_reporting_role"] == "strict_picked_report"
        for row in false_positives
    )
    report_layer_fp = len(false_positives) - strict_fp
    summary = {
        "scope": (
            "Post-hoc answer-aware audit only; no scorer or reporting decision "
            "is changed by this tool."
        ),
        "shadow_root": str(shadow_root.resolve()),
        "evaluation_root": str(evaluation_root.resolve()),
        "complete_report_rows": len(rows),
        "matched_answer_rows": len(true_positives),
        "unmatched_labeled_rows": len(false_positives),
        "unlabeled_rows": len(unlabeled),
        "strict_picked_false_positives": strict_fp,
        "report_layer_false_positives": report_layer_fp,
        "theoretical_labeled_prediction_floor_if_only_report_layer_fp_removed": (
            len(true_positives) + strict_fp
        ),
        "false_positive_role_counts": dict(sorted(Counter(
            row["possible_reporting_role"] for row in false_positives
        ).items())),
        "false_positive_route_counts": dict(sorted(Counter(
            row["possible_route"] for row in false_positives
        ).items())),
        "false_positive_family_counts": dict(sorted(Counter(
            row["taxonomy_family"] for row in false_positives
        ).items())),
        "true_positive_role_counts": dict(sorted(Counter(
            row["possible_reporting_role"] for row in true_positives
        ).items())),
        "constraints": [
            "Strict Picked false positives originate upstream and cannot be removed by a Possible-only guardrail.",
            "Any proposed guardrail must be tested against true_positive_contrast.csv before implementation.",
            "The benchmark is development data; answer-aware patterns are hypotheses, not production rules.",
        ],
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shadow_root", type=Path)
    parser.add_argument("evaluation_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(
        args.shadow_root, args.evaluation_root, args.output_dir
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

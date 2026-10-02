"""Post-hoc value audit for same-specimen technical-repeat evidence.

The routing variants are defined without benchmark answers.  Physician answers
are joined only after the variants are fixed, to quantify their retrospective
precision/recall trade-offs.  No frozen pipeline output is overwritten.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable

from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.evaluate_test_aware_deterministic_shadow import group_names, merge_groups
from tools.recalculate_kh_benchmark_metrics import match_type, unique_names


DEFAULT_REPEAT_AUDIT = Path(
    "outputs/runs/2026-09-30_KH_mngs_selection_protocol_semantics_audit_v1"
    "/technical_repeat_scorer_counterfactual.csv"
)
DEFAULT_SCORER = Path(
    "outputs/runs/2026-09-25_KH_test_aware_deterministic_taxonomy_auto_v1"
)
DEFAULT_INTEGRATED = Path(
    "outputs/runs/2026-09-29_KH_integrated_decision_pipeline_v1"
    "/integrated_decisions.csv"
)
DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics"
    "/kh_answers_clinical_revision_20260918.csv"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-30_KH_technical_repeat_value_audit_v1"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def as_bool(value: Any) -> bool:
    return str(value or "").strip().lower() == "true"


def candidate_key(row: dict[str, Any]) -> tuple[str, str]:
    return str(row.get("patient_id") or ""), str(row.get("organism_name") or "")


def benchmark_status(
    patient_id: str,
    organism_name: str,
    answers: dict[int, list[str]],
) -> tuple[str, str]:
    gold = answers.get(int(patient_id), [])
    for answer in gold:
        if match_type(organism_name, answer, genus_relaxed=False):
            return "matched", answer
    return ("unmatched", "") if gold else ("unlabeled_patient", "")


def variant_candidate_rows(
    current_rows: list[dict[str, str]],
    changed_repeat_rows: dict[tuple[str, str], dict[str, str]],
    mode: str,
) -> list[dict[str, str]]:
    output = copy.deepcopy(current_rows)
    for row in output:
        change = changed_repeat_rows.get(candidate_key(row))
        if change is None:
            continue
        if mode == "context_only":
            # Technical consistency may retain visibility but cannot create
            # generic High or Picked.
            if row.get("decision") in {"picked_shadow", "review_high_priority"}:
                row["decision"] = "review_context_needed"
                row["integrated_level"] = "Level 4"
        elif mode == "audit_only":
            row["decision"] = change["no_technical_repeat_decision"]
            row["integrated_level"] = change["no_technical_repeat_level"]
        elif mode != "current_full":
            raise ValueError(f"Unsupported technical-repeat mode: {mode}")
    return output


def scorer_metrics(
    candidate_rows: list[dict[str, str]],
    hospital_rows: list[dict[str, str]],
    answers: dict[int, list[str]],
) -> dict[str, Any]:
    picked = merge_groups(
        group_names(candidate_rows, lambda row: row["decision"] == "picked_shadow"),
        group_names(hospital_rows, lambda row: row["decision"] == "picked_shadow"),
    )
    high = merge_groups(
        group_names(
            candidate_rows, lambda row: row["decision"] == "review_high_priority"
        ),
        group_names(
            hospital_rows, lambda row: row["decision"] == "review_high_priority"
        ),
    )
    context = merge_groups(
        group_names(
            candidate_rows, lambda row: row["decision"] == "review_context_needed"
        ),
        group_names(
            hospital_rows, lambda row: row["decision"] == "review_context_needed"
        ),
    )
    low = group_names(
        candidate_rows, lambda row: row["decision"] == "review_low_specificity"
    )
    return {
        "picked": evaluate_set(picked, answers),
        "picked_plus_high": evaluate_set(merge_groups(picked, high), answers),
        "forwarded_picked_high_context": evaluate_set(
            merge_groups(picked, high, context), answers
        ),
        "all_positive_routed": evaluate_set(
            merge_groups(picked, high, context, low), answers
        ),
    }


def report_predictions(rows: list[dict[str, str]]) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        if as_bool(row.get("selected_for_complete_report")):
            grouped[int(row["patient_id"])].append(row["organism_name"])
    return {
        patient_id: unique_names(names) for patient_id, names in grouped.items()
    }


def audit(
    repeat_audit_path: Path,
    scorer_dir: Path,
    integrated_path: Path,
    answers_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    repeat_rows = read_csv(repeat_audit_path)
    case_candidate_rows = read_csv(scorer_dir / "case_candidate_decisions.csv")
    patient_candidate_rows = read_csv(scorer_dir / "patient_organism_decisions.csv")
    hospital_rows = read_csv(scorer_dir / "hospital_only_decisions.csv")
    integrated_rows = read_csv(integrated_path)
    answers = load_answer_rows(answers_path)

    changed_rows = [row for row in repeat_rows if as_bool(row["decision_or_level_changed"])]
    changed_by_patient_organism = {candidate_key(row): row for row in changed_rows}
    integrated_index = {candidate_key(row): row for row in integrated_rows}

    value_rows = []
    for row in repeat_rows:
        final = integrated_index.get(candidate_key(row))
        status, matched_answer = benchmark_status(
            row["patient_id"], row["organism_name"], answers
        )
        value_rows.append(
            {
                "patient_id": row["patient_id"],
                "case_review_id": row["case_review_id"],
                "organism_name": row["organism_name"],
                "repeat_group_type": row["repeat_group_type"],
                "selected_test_count": row["selected_test_count"],
                "best_rank_in_retained_universe": row[
                    "best_rank_in_retained_universe"
                ],
                "current_scorer_decision": row["current_decision"],
                "audit_only_scorer_decision": row[
                    "no_technical_repeat_decision"
                ],
                "scorer_route_depends_on_technical_repeat": row[
                    "decision_or_level_changed"
                ],
                "current_final_reporting_tier": (
                    final.get("final_reporting_tier") if final else "absent"
                ),
                "selected_for_complete_report": (
                    final.get("selected_for_complete_report") if final else "False"
                ),
                "reporting_role": final.get("reporting_role") if final else "",
                "benchmark_status": status,
                "matched_answer": matched_answer,
                "interpretation": (
                    "post_hoc benchmark annotation; not used to define the routing variants"
                ),
            }
        )

    scorer_variants = {}
    for mode in ("current_full", "context_only", "audit_only"):
        rows = variant_candidate_rows(
            patient_candidate_rows, changed_by_patient_organism, mode
        )
        scorer_variants[mode] = scorer_metrics(rows, hospital_rows, answers)

    current_report_predictions = report_predictions(integrated_rows)
    current_report_metrics = evaluate_set(current_report_predictions, answers)
    route_dependent_complete_keys = {
        candidate_key(row)
        for row in value_rows
        if as_bool(row["scorer_route_depends_on_technical_repeat"])
        and as_bool(row["selected_for_complete_report"])
    }
    pjp_key = ("9", "Pneumocystis jirovecii")

    audit_only_report_rows = [
        row
        for row in integrated_rows
        if candidate_key(row) not in route_dependent_complete_keys
    ]
    context_pjp_report_rows = [
        row
        for row in integrated_rows
        if candidate_key(row) not in (route_dependent_complete_keys - {pjp_key})
    ]
    report_variants = {
        "current_full": current_report_metrics,
        "context_only_with_pjp_family_bridge": evaluate_set(
            report_predictions(context_pjp_report_rows), answers
        ),
        "audit_only_no_repeat_reporting_support": evaluate_set(
            report_predictions(audit_only_report_rows), answers
        ),
    }

    changed_value_rows = [
        row for row in value_rows if as_bool(row["scorer_route_depends_on_technical_repeat"])
    ]
    complete_value_rows = [
        row for row in value_rows if as_bool(row["selected_for_complete_report"])
    ]
    changed_complete_value_rows = [
        row
        for row in changed_value_rows
        if as_bool(row["selected_for_complete_report"])
    ]

    assertions = {
        "technical_repeat_only_entries_are_48": len(value_rows) == 48,
        "scorer_route_dependent_entries_are_28": len(changed_value_rows) == 28,
        "route_dependent_benchmark_matches_are_2": (
            sum(row["benchmark_status"] == "matched" for row in changed_value_rows)
            == 2
        ),
        "route_dependent_complete_report_entry_is_only_p9_pjp": (
            {(row["patient_id"], row["organism_name"]) for row in changed_complete_value_rows}
            == {pjp_key}
        ),
        "complete_report_technical_repeat_entries_are_8": len(complete_value_rows)
        == 8,
        "complete_report_technical_repeat_matches_are_6": (
            sum(row["benchmark_status"] == "matched" for row in complete_value_rows)
            == 6
        ),
        "context_only_preserves_forwarding_recall": (
            scorer_variants["context_only"]["forwarded_picked_high_context"][
                "matched"
            ]
            == 59
        ),
        "audit_only_loses_one_forwarded_answer": (
            scorer_variants["audit_only"]["forwarded_picked_high_context"][
                "matched"
            ]
            == 58
        ),
        "context_pjp_report_preserves_55_complete_matches": (
            report_variants["context_only_with_pjp_family_bridge"]["matched"]
            == 55
        ),
        "audit_only_report_loses_p9_pjp": (
            report_variants["audit_only_no_repeat_reporting_support"]["matched"]
            == 54
        ),
    }

    result = {
        "schema_version": "kh_technical_repeat_value_audit.v1",
        "scope": (
            "Post-hoc value audit. Routing variants were fixed before answer joining; "
            "metrics are retrospective and not external diagnostic validation."
        ),
        "inputs": {
            "repeat_counterfactual": str(repeat_audit_path.resolve()),
            "scorer_dir": str(scorer_dir.resolve()),
            "integrated_decisions": str(integrated_path.resolve()),
            "answers": str(answers_path.resolve()),
        },
        "counts": {
            "technical_repeat_only_entries": len(value_rows),
            "scorer_route_dependent_entries": len(changed_value_rows),
            "scorer_route_dependent_benchmark_status": dict(
                sorted(Counter(row["benchmark_status"] for row in changed_value_rows).items())
            ),
            "technical_repeat_entries_selected_for_complete_report": len(
                complete_value_rows
            ),
            "complete_report_benchmark_status": dict(
                sorted(Counter(row["benchmark_status"] for row in complete_value_rows).items())
            ),
            "route_dependent_complete_report_entries": [
                {
                    "patient_id": row["patient_id"],
                    "organism_name": row["organism_name"],
                    "final_reporting_tier": row["current_final_reporting_tier"],
                    "benchmark_status": row["benchmark_status"],
                }
                for row in changed_complete_value_rows
            ],
        },
        "scorer_variant_metrics": scorer_variants,
        "complete_report_variant_metrics": report_variants,
        "recommended_policy": {
            "generic_high_or_picked": (
                "same-molecule technical repeat cannot create generic High or Picked"
            ),
            "context_visibility": (
                "retain technical repeat as one weak same-specimen consistency axis so positive candidates remain visible"
            ),
            "pjp_possible": (
                "permit a family-specific Possible route only with exact PJP identity, compatible specimen, top-three rank, at least two technical branches, and recorded host context"
            ),
            "lab_mngs_priority": (
                "use only as an intra-tier tie-breaker; it does not change tier or metrics"
            ),
            "next_validation": (
                "implement and rerun the full answer-blind context-only plus PJP-family bridge shadow before adopting"
            ),
        },
        "assertions": assertions,
        "assertion_pass_count": sum(assertions.values()),
        "assertion_total_count": len(assertions),
        "all_assertions_pass": all(assertions.values()),
    }
    write_csv(
        output_dir / "technical_repeat_candidate_value.csv",
        value_rows,
        list(value_rows[0]),
    )
    write_json(output_dir / "summary.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze the retrospective value of technical-repeat routing"
    )
    parser.add_argument("--repeat-audit", type=Path, default=DEFAULT_REPEAT_AUDIT)
    parser.add_argument("--scorer-dir", type=Path, default=DEFAULT_SCORER)
    parser.add_argument("--integrated", type=Path, default=DEFAULT_INTEGRATED)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit(
        args.repeat_audit,
        args.scorer_dir,
        args.integrated,
        args.answers,
        args.output_dir,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

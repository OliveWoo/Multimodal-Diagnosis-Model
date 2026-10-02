"""Shadow two generalized low-specificity technical-repeat routes.

The answer-blind decision table is written and hashed before benchmark answers
are loaded.  Route A can move an analytically compelling low-specificity
candidate from Context to scorer High.  Route B preserves a best-available
low-specificity candidate as Fallback-Possible plus high review priority; it
does not claim High analytical confidence.  Neither route can create Picked.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any

from tools.analyze_kh_family_aware_technical_repeat_shadow import (
    DEFAULT_ANSWERS,
    DEFAULT_INTEGRATED,
    DEFAULT_SCORER,
    as_bool,
    benchmark_status,
    key,
    read_csv,
    report_predictions,
    scorer_metrics,
    sha256_file,
    write_csv,
    write_json,
)
from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows


DEFAULT_BASE_POLICY = Path(
    "outputs/runs/2026-09-30_KH_family_aware_technical_repeat_high_shadow_v1"
    "/answer_blind_policy_decisions.csv"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-30_KH_low_specificity_repeat_routes_shadow_v2"
)
LOW_SPECIFICITY_FAMILIES = {
    "environmental_low_specificity",
    "skin_airway_colonizer_prone",
}
DEFAULT_RPM_THRESHOLD = 10.0
SENSITIVITY_THRESHOLDS = (5.0, 10.0, 50.0)


def _split_ints(value: Any) -> list[int]:
    output: list[int] = []
    for item in str(value or "").split("|"):
        try:
            output.append(int(float(item)))
        except (TypeError, ValueError):
            continue
    return output


def _split_floats(value: Any) -> list[float]:
    output: list[float] = []
    for item in str(value or "").split("|"):
        if not item.strip():
            continue
        try:
            output.append(float(item))
        except (TypeError, ValueError):
            continue
    return output


def classify_routes(
    base_row: dict[str, Any],
    integrated_row: dict[str, Any],
    *,
    rpm_threshold: float = DEFAULT_RPM_THRESHOLD,
) -> dict[str, Any]:
    family = str(
        integrated_row.get("taxonomy_family")
        or base_row.get("taxonomy_family")
        or ""
    )
    mapping_status = str(integrated_row.get("taxonomy_mapping_status") or "")
    ranks = _split_ints(base_row.get("branch_ranks"))
    rpms = _split_floats(base_row.get("branch_rpms"))
    base_strength = as_bool(base_row.get("rank_only_strength_eligible"))
    low_specificity = family in LOW_SPECIFICITY_FAMILIES
    all_rank_one = len(ranks) >= 2 and all(rank == 1 for rank in ranks)
    all_branches_normalized = len(rpms) >= 2 and len(rpms) == len(ranks)
    minimum_rpm = min(rpms) if rpms else None
    burden_gate = bool(
        all_branches_normalized
        and minimum_rpm is not None
        and minimum_rpm >= rpm_threshold
    )
    route_a = bool(
        base_strength
        and low_specificity
        and mapping_status == "exact_species"
        and all_rank_one
        and burden_gate
        and str(base_row.get("repeat_group_type") or "")
        != "explicit_original_library_rerun"
    )
    route_b = bool(
        base_strength
        and low_specificity
        and not route_a
        and str(integrated_row.get("final_reporting_tier") or "")
        == "Fallback-Possible"
        and str(integrated_row.get("reporting_route") or "")
        == "no_picked_or_possible_best_available_fallback"
    )
    if route_a:
        reason = "analytically_compelling_low_specificity_repeat"
    elif route_b:
        reason = "best_available_low_specificity_fallback"
    elif not base_strength:
        reason = "base_technical_strength_gate_not_met"
    elif not low_specificity:
        reason = "not_target_low_specificity_family"
    elif mapping_status != "exact_species":
        reason = "identity_not_exact_species"
    elif not all_rank_one:
        reason = "not_all_selected_branches_rank_one"
    elif not all_branches_normalized:
        reason = "normalization_incomplete"
    elif not burden_gate:
        reason = "normalized_burden_below_shadow_threshold"
    else:
        reason = "no_generalized_route"
    return {
        "taxonomy_family": family,
        "taxonomy_mapping_status": mapping_status,
        "branch_ranks_parsed": "|".join(map(str, ranks)),
        "branch_rpms_parsed": "|".join(str(value) for value in rpms),
        "all_selected_branches_rank_one": all_rank_one,
        "all_selected_branches_have_rpm": all_branches_normalized,
        "minimum_selected_branch_rpm": "" if minimum_rpm is None else minimum_rpm,
        "rpm_threshold": rpm_threshold,
        "route_a_scorer_high_eligible": route_a,
        "route_b_fallback_review_eligible": route_b,
        "route_reason": reason,
        "picked_created": False,
    }


def _apply_route_a_high(
    candidate_rows: list[dict[str, str]], eligible: set[tuple[str, str]]
) -> list[dict[str, str]]:
    output = copy.deepcopy(candidate_rows)
    for row in output:
        if key(row) in eligible and row.get("decision") != "picked_shadow":
            row["decision"] = "review_high_priority"
            row["integrated_level"] = "Level 3"
    return output


def _apply_route_a_possible(
    integrated_rows: list[dict[str, str]], eligible: set[tuple[str, str]]
) -> list[dict[str, str]]:
    output = copy.deepcopy(integrated_rows)
    for row in output:
        if key(row) in eligible and not as_bool(row.get("selected_for_complete_report")):
            row["selected_for_complete_report"] = "True"
            row["final_reporting_tier"] = "Possible"
            row["reporting_role"] = "analytically_compelling_possible_shadow"
            row["reporting_route"] = "low_specificity_repeat_route_a_shadow"
    patients_with_formal_report_candidate = {
        str(row.get("patient_id") or "")
        for row in output
        if as_bool(row.get("selected_for_complete_report"))
        and str(row.get("final_reporting_tier") or "") in {"Picked", "Possible"}
    }
    for row in output:
        if (
            str(row.get("patient_id") or "") in patients_with_formal_report_candidate
            and as_bool(row.get("selected_for_complete_report"))
            and str(row.get("final_reporting_tier") or "") == "Fallback-Possible"
        ):
            row["selected_for_complete_report"] = "False"
            row["final_reporting_tier"] = "Context"
            row["reporting_role"] = "not_selected_for_complete_report"
            row["reporting_route"] = "fallback_suppressed_by_picked_or_possible_shadow"
    return output


def _metric_delta(after: dict[str, Any], before: dict[str, Any]) -> dict[str, Any]:
    before_tp = int(before.get("matched") or 0)
    after_tp = int(after.get("matched") or 0)
    before_predictions = int(before.get("predicted_in_labeled_patients") or 0)
    after_predictions = int(after.get("predicted_in_labeled_patients") or 0)
    return {
        "tp_delta": after_tp - before_tp,
        "fp_delta": (after_predictions - after_tp) - (before_predictions - before_tp),
        "prediction_delta": after_predictions - before_predictions,
        "precision_delta": float(after.get("precision") or 0)
        - float(before.get("precision") or 0),
        "recall_delta": float(after.get("recall") or 0)
        - float(before.get("recall") or 0),
        "f1_delta": float(after.get("f1") or 0) - float(before.get("f1") or 0),
    }


def audit(
    base_policy_path: Path,
    scorer_dir: Path,
    integrated_path: Path,
    answers_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    base_rows = read_csv(base_policy_path)
    candidate_rows = read_csv(scorer_dir / "patient_organism_decisions.csv")
    hospital_rows = read_csv(scorer_dir / "hospital_only_decisions.csv")
    integrated_rows = read_csv(integrated_path)
    integrated_index = {key(row): row for row in integrated_rows}

    policy = {
        "schema_version": "kh_low_specificity_repeat_routes_shadow.v2",
        "mode": "answer_blind_non_mutating_shadow",
        "route_a": {
            "result": "scorer High; never Picked",
            "families": sorted(LOW_SPECIFICITY_FAMILIES),
            "requirements": [
                "lower-respiratory species-level technical-repeat strength gate",
                "exact-species taxonomy mapping",
                "at least two selected branches and every branch rank equals 1",
                "every selected branch has normalized RPM",
                f"minimum selected-branch RPM >= {DEFAULT_RPM_THRESHOLD:g}",
                "not an explicit same-library rerun",
            ],
        },
        "route_b": {
            "result": "Fallback-Possible plus high review priority; not scorer High",
            "requirements": [
                "same answer-blind technical strength gate",
                "low-specificity family",
                "existing no-Picked/no-Possible best-available fallback route",
                "route A not met",
            ],
        },
        "sensitivity_thresholds": list(SENSITIVITY_THRESHOLDS),
        "benchmark_answers_used_by_policy": False,
        "laboratory_caution": (
            "same-run 12000g/unannotated branches are analytically separate branches, "
            "but biological independence is not assumed."
        ),
    }
    write_json(output_dir / "policy.json", policy)

    decision_rows: list[dict[str, Any]] = []
    for base in base_rows:
        integrated = integrated_index.get(key(base), {})
        route = classify_routes(base, integrated)
        decision_rows.append(
            {
                "patient_id": base.get("patient_id", ""),
                "case_review_id": base.get("case_review_id", ""),
                "organism_name": base.get("organism_name", ""),
                "organism_key": base.get("organism_key", ""),
                "repeat_group_type": base.get("repeat_group_type", ""),
                "base_rank_only_strength_eligible": base.get(
                    "rank_only_strength_eligible", ""
                ),
                "current_scorer_decision": base.get("current_scorer_decision", ""),
                "current_reporting_tier": integrated.get("final_reporting_tier", ""),
                "current_selected_for_complete_report": integrated.get(
                    "selected_for_complete_report", ""
                ),
                **route,
                "answer_blind": True,
            }
        )
    fields = list(decision_rows[0])
    decision_path = output_dir / "answer_blind_policy_decisions.csv"
    write_csv(decision_path, decision_rows, fields)
    decision_hash = sha256_file(decision_path)

    route_a_keys = {
        key(row)
        for row in decision_rows
        if as_bool(row["route_a_scorer_high_eligible"])
    }
    route_b_keys = {
        key(row)
        for row in decision_rows
        if as_bool(row["route_b_fallback_review_eligible"])
    }
    route_a_scorer_rows = _apply_route_a_high(candidate_rows, route_a_keys)
    route_a_report_rows = _apply_route_a_possible(integrated_rows, route_a_keys)

    # Answers are joined only after the answer-blind table has been frozen.
    answers = load_answer_rows(answers_path)
    benchmark_rows: list[dict[str, Any]] = []
    for row in decision_rows:
        status, answer = benchmark_status(
            str(row["patient_id"]), str(row["organism_name"]), answers
        )
        benchmark_rows.append(
            {
                **row,
                "benchmark_status": status,
                "matched_answer": answer,
                "benchmark_join_stage": "post_policy_freeze",
            }
        )
    write_csv(
        output_dir / "benchmark_comparison.csv",
        benchmark_rows,
        list(benchmark_rows[0]),
    )

    current_scorer = scorer_metrics(candidate_rows, hospital_rows, answers)
    route_a_scorer = scorer_metrics(route_a_scorer_rows, hospital_rows, answers)
    current_report = evaluate_set(report_predictions(integrated_rows), answers)
    scorer_only_report = evaluate_set(report_predictions(integrated_rows), answers)
    route_a_possible_report = evaluate_set(
        report_predictions(route_a_report_rows), answers
    )

    sensitivity: dict[str, Any] = {}
    for threshold in SENSITIVITY_THRESHOLDS:
        eligible = {
            key(base)
            for base in base_rows
            if classify_routes(
                base, integrated_index.get(key(base), {}), rpm_threshold=threshold
            )["route_a_scorer_high_eligible"]
        }
        sensitivity[f"rpm_{threshold:g}"] = {
            "eligible_count": len(eligible),
            "eligible_candidates": [
                {"patient_id": patient_id, "organism_name": organism_name}
                for patient_id, organism_name in sorted(eligible)
            ],
        }

    selected_route_rows = [
        row
        for row in benchmark_rows
        if as_bool(row["route_a_scorer_high_eligible"])
        or as_bool(row["route_b_fallback_review_eligible"])
    ]
    assertions = {
        "answer_blind_policy_has_no_benchmark_columns": not any(
            "benchmark" in field or "answer" in field
            for field in fields
            if field != "answer_blind"
        ),
        "all_base_repeat_candidates_audited": len(decision_rows) == len(base_rows) == 48,
        "no_route_can_create_picked": not any(
            as_bool(row["picked_created"]) for row in decision_rows
        ),
        "route_a_is_stable_across_rpm_sensitivity": len(
            {
                tuple(
                    (item["patient_id"], item["organism_name"])
                    for item in value["eligible_candidates"]
                )
                for value in sensitivity.values()
            }
        )
        == 1,
        "route_b_does_not_change_complete_report_membership": (
            all(
                as_bool(integrated_index[item].get("selected_for_complete_report"))
                for item in route_b_keys
            )
        ),
    }

    result = {
        "schema_version": "kh_low_specificity_repeat_routes_shadow.v2",
        "scope": (
            "Answer-blind non-mutating shadow with post-freeze retrospective benchmark "
            "evaluation; not external validation."
        ),
        "inputs": {
            "base_policy": str(base_policy_path.resolve()),
            "scorer": str(scorer_dir.resolve()),
            "integrated": str(integrated_path.resolve()),
            "answers_posthoc_only": str(answers_path.resolve()),
        },
        "policy": policy,
        "answer_blind_policy_decisions_sha256": decision_hash,
        "counts": {
            "technical_repeat_candidates": len(decision_rows),
            "route_a_scorer_high_eligible": len(route_a_keys),
            "route_b_fallback_review_eligible": len(route_b_keys),
            "selected_routes_posthoc": [
                {
                    "patient_id": row["patient_id"],
                    "organism_name": row["organism_name"],
                    "route_reason": row["route_reason"],
                    "benchmark_status": row["benchmark_status"],
                }
                for row in selected_route_rows
            ],
        },
        "scorer_metrics": {
            "current": current_scorer,
            "route_a_high": route_a_scorer,
            "picked_plus_high_delta": _metric_delta(
                route_a_scorer["picked_plus_high"],
                current_scorer["picked_plus_high"],
            ),
        },
        "complete_report_metrics": {
            "current": current_report,
            "route_a_scorer_high_only": scorer_only_report,
            "route_a_high_plus_possible_reporting_sensitivity": route_a_possible_report,
            "high_plus_possible_delta": _metric_delta(
                route_a_possible_report, current_report
            ),
            "route_b_membership_delta": {
                "tp_delta": 0,
                "fp_delta": 0,
                "reason": "The generalized best-available fallback is already active for its eligible case.",
            },
        },
        "rpm_sensitivity": sensitivity,
        "interpretation": {
            "route_a": (
                "A robust analytically compelling route can send qualifying exact-species "
                "low-specificity candidates to scorer High without creating Picked."
            ),
            "route_b": (
                "Best-available fallback is a reporting visibility and review-priority route, "
                "not a claim of High analytical confidence."
            ),
            "reporting": (
                "Scorer High alone does not automatically enter the complete report. The "
                "Possible-report scenario is shown separately and recomputes the mutually "
                "exclusive best-available fallback for that patient."
            ),
        },
        "assertions": assertions,
        "all_assertions_pass": all(assertions.values()),
    }
    write_json(output_dir / "summary.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-policy", type=Path, default=DEFAULT_BASE_POLICY)
    parser.add_argument("--scorer-dir", type=Path, default=DEFAULT_SCORER)
    parser.add_argument("--integrated", type=Path, default=DEFAULT_INTEGRATED)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit(
        args.base_policy,
        args.scorer_dir,
        args.integrated,
        args.answers,
        args.output_dir,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

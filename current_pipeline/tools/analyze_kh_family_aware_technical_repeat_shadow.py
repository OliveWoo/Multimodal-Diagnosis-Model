"""Answer-blind family-aware technical-repeat High shadow.

The policy is frozen and written before benchmark answers are loaded.  It does
not overwrite the scorer, history, promotion, or integrated reporting runs.
Benchmark labels are joined only after the answer-blind decision table has been
written and hashed.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.evaluate_test_aware_deterministic_shadow import group_names, merge_groups
from tools.recalculate_kh_benchmark_metrics import match_type, unique_names


DEFAULT_REPEAT_AUDIT = Path(
    "outputs/runs/2026-09-30_KH_mngs_selection_protocol_semantics_audit_v1"
    "/technical_repeat_scorer_counterfactual.csv"
)
DEFAULT_COMPACT = Path(
    "outputs/runs/2026-09-25_KH_compact_multi_assay_entry_taxonomy_auto_v1"
    "/patient_packets"
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
    "outputs/runs/2026-09-30_KH_family_aware_technical_repeat_high_shadow_v1"
)

HIGH_FAMILIES = {
    "typical_respiratory_pathogen",
    "hospital_or_nonfermenter_gnb",
}
OPPORTUNISTIC_FAMILIES = {
    "high_consequence_opportunistic",
    "mold_or_opportunistic_fungus",
}
RNA_RESPIRATORY_FAMILIES = {"other_respiratory_virus"}
BLOCKED_CONTEXT_FAMILIES = {
    "candida_or_yeast",
    "herpes_or_reactivation_virus",
    "oral_aspiration_or_anaerobe",
    "skin_airway_colonizer_prone",
    "environmental_low_specificity",
    "gi_urinary_or_nonpulmonary_prone",
    "commensal_virome_or_endogenous_element",
    "unmapped_or_uncertain",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
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


def as_bool(value: Any) -> bool:
    return str(value or "").strip().lower() == "true"


def as_int(value: Any, default: int = 999) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def key(row: dict[str, Any]) -> tuple[str, str]:
    return str(row.get("patient_id") or ""), str(row.get("organism_name") or "")


def exact_species_identity(name: str) -> bool:
    normalized = " ".join(str(name or "").strip().lower().split())
    ambiguous = (" sp.", " species", " group", " complex", "unclassified")
    return bool(normalized) and not any(token in normalized for token in ambiguous)


def compact_entry_index(packet_dir: Path) -> dict[tuple[str, str, str], dict[str, Any]]:
    output: dict[tuple[str, str, str], dict[str, Any]] = {}
    for path in sorted(packet_dir.glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        patient_id = str(payload["patient_id"])
        for case in payload.get("cases") or []:
            case_id = str(case.get("case_review_id") or "")
            for entry in case.get("scorer_entries") or []:
                output[(patient_id, case_id, str(entry.get("organism_key") or ""))] = entry
    return output


def policy_definition() -> dict[str, Any]:
    return {
        "schema_version": "kh_family_aware_technical_repeat_high_shadow.v1",
        "mode": "answer_blind_non_mutating_shadow",
        "target": "same-specimen, same-molecule technical-repeat candidates",
        "technical_strength_gate": {
            "minimum_selected_branches": 2,
            "all_selected_branch_ranks_at_most": 3,
            "exclude_repeat_group_types": ["explicit_original_library_rerun"],
            "require_species_level_identity": True,
            "require_compatible_specimen": "lower_respiratory",
            "rpm_threshold": None,
            "qc_missingness": "unknown_not_negative",
        },
        "family_routes": {
            "ordinary_high": sorted(HIGH_FAMILIES),
            "opportunistic_high_requires_history_priority": sorted(
                OPPORTUNISTIC_FAMILIES
            ),
            "rna_respiratory_high_requires_all_rna_branches": sorted(
                RNA_RESPIRATORY_FAMILIES
            ),
            "blocked_from_repeat_only_high": sorted(BLOCKED_CONTEXT_FAMILIES),
        },
        "picked_gate": (
            "Technical repeat alone never creates Picked; a separately approved "
            "independent evidence axis remains required."
        ),
        "same_library_rerun": "audit_only",
        "lab_mngs_priority": "intra_tier_tie_break_only",
        "benchmark_answers_used_by_policy": False,
    }


def family_route(
    *, family: str, history_route: str, all_rna: bool
) -> tuple[bool, str]:
    if family in HIGH_FAMILIES:
        return True, "ordinary_respiratory_or_hospital_family"
    if family in OPPORTUNISTIC_FAMILIES:
        if history_route == "Priority":
            return True, "opportunistic_family_with_history_priority"
        return False, "opportunistic_family_without_history_priority"
    if family in RNA_RESPIRATORY_FAMILIES:
        if all_rna:
            return True, "rna_respiratory_family_with_rna_repeat"
        return False, "rna_respiratory_family_without_all_rna_branches"
    if family in BLOCKED_CONTEXT_FAMILIES:
        return False, "context_or_low_specificity_family_block"
    return False, "family_not_approved_for_repeat_only_high"


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
    return {
        "picked": evaluate_set(picked, answers),
        "picked_plus_high": evaluate_set(merge_groups(picked, high), answers),
        "forwarded_picked_high_context": evaluate_set(
            merge_groups(picked, high, context), answers
        ),
    }


def report_predictions(rows: list[dict[str, str]]) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        if as_bool(row.get("selected_for_complete_report")):
            grouped[int(row["patient_id"])].append(row["organism_name"])
    return {patient_id: unique_names(names) for patient_id, names in grouped.items()}


def benchmark_status(
    patient_id: str, organism_name: str, answers: dict[int, list[str]]
) -> tuple[str, str]:
    gold = answers.get(int(patient_id), [])
    for answer in gold:
        if match_type(organism_name, answer, genus_relaxed=False):
            return "matched", answer
    return ("unmatched", "") if gold else ("unlabeled_patient", "")


def apply_context_only_baseline(
    rows: list[dict[str, str]], route_dependent: set[tuple[str, str]]
) -> list[dict[str, str]]:
    output = copy.deepcopy(rows)
    for row in output:
        if key(row) not in route_dependent:
            continue
        if row["decision"] in {"picked_shadow", "review_high_priority"}:
            row["decision"] = "review_context_needed"
            row["integrated_level"] = "Level 4"
    return output


def apply_shadow_high(
    context_rows: list[dict[str, str]],
    policy_rows: list[dict[str, Any]],
    *,
    open_rank_only: bool = False,
) -> list[dict[str, str]]:
    output = copy.deepcopy(context_rows)
    row_index = {key(row): row for row in output}
    for policy_row in policy_rows:
        eligible = (
            as_bool(policy_row["rank_only_strength_eligible"])
            if open_rank_only
            else as_bool(policy_row["family_aware_high_eligible"])
        )
        if not eligible:
            continue
        row = row_index.get(key(policy_row))
        if row is None or row["decision"] == "picked_shadow":
            continue
        row["decision"] = "review_high_priority"
        row["integrated_level"] = "Level 3"
    return output


def audit(
    repeat_audit_path: Path,
    compact_dir: Path,
    scorer_dir: Path,
    integrated_path: Path,
    answers_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    repeat_rows = read_csv(repeat_audit_path)
    entries = compact_entry_index(compact_dir)
    candidate_rows = read_csv(scorer_dir / "patient_organism_decisions.csv")
    hospital_rows = read_csv(scorer_dir / "hospital_only_decisions.csv")
    integrated_rows = read_csv(integrated_path)
    integrated_index = {key(row): row for row in integrated_rows}

    policy = policy_definition()
    write_json(output_dir / "policy.json", policy)

    decision_rows: list[dict[str, Any]] = []
    for row in repeat_rows:
        entry = entries[
            (
                str(row["patient_id"]),
                str(row["case_review_id"]),
                str(row["organism_key"]),
            )
        ]
        integrated = integrated_index.get(key(row), {})
        signals = entry.get("selected_test_signals") or []
        ranks = [as_int(signal.get("rank_in_retained_universe")) for signal in signals]
        nucleic_types = [str(signal.get("nucleic_type") or "") for signal in signals]
        qcs = [str(signal.get("qc_status") or "") for signal in signals]
        rpms = [signal.get("rpm_total") for signal in signals]
        species_level = exact_species_identity(row["organism_name"])
        not_same_library = row["repeat_group_type"] != "explicit_original_library_rerun"
        all_top_three = len(signals) >= 2 and all(rank <= 3 for rank in ranks)
        compatible_specimen = integrated.get("specimen_context") == "lower_respiratory"
        strength_eligible = bool(
            len(signals) >= 2
            and not_same_library
            and all_top_three
            and species_level
            and compatible_specimen
        )
        family = str(integrated.get("taxonomy_family") or row.get("taxonomy_family") or "")
        history_route = str(integrated.get("history_adjusted_route") or "")
        all_rna = bool(nucleic_types) and all(value == "RNA" for value in nucleic_types)
        route_allowed, route_reason = family_route(
            family=family, history_route=history_route, all_rna=all_rna
        )
        eligible = strength_eligible and route_allowed
        if not strength_eligible:
            if not not_same_library:
                final_reason = "same_library_rerun_audit_only"
            elif not all_top_three:
                final_reason = "not_all_technical_branches_top_three"
            elif not species_level:
                final_reason = "identity_not_species_level"
            elif not compatible_specimen:
                final_reason = "specimen_not_lower_respiratory"
            else:
                final_reason = "technical_strength_gate_not_met"
        else:
            final_reason = route_reason

        context_decision = (
            "review_context_needed"
            if as_bool(row["decision_or_level_changed"])
            and row["current_decision"] in {"picked_shadow", "review_high_priority"}
            else row["current_decision"]
        )
        shadow_decision = context_decision
        if eligible and shadow_decision != "picked_shadow":
            shadow_decision = "review_high_priority"
        decision_rows.append(
            {
                "patient_id": row["patient_id"],
                "case_review_id": row["case_review_id"],
                "organism_name": row["organism_name"],
                "organism_key": row["organism_key"],
                "taxonomy_family": family,
                "specimen_context": integrated.get("specimen_context", ""),
                "history_adjusted_route": history_route,
                "repeat_group_type": row["repeat_group_type"],
                "selected_branch_count": len(signals),
                "branch_ranks": "|".join(str(value) for value in ranks),
                "branch_nucleic_types": "|".join(nucleic_types),
                "branch_qc_statuses": "|".join(qcs),
                "branch_rpms": "|".join("" if value is None else str(value) for value in rpms),
                "species_level_identity": species_level,
                "all_branches_top_three": all_top_three,
                "compatible_lower_respiratory_specimen": compatible_specimen,
                "rank_only_strength_eligible": strength_eligible,
                "family_route_allowed": route_allowed,
                "family_aware_high_eligible": eligible,
                "policy_reason": final_reason,
                "current_scorer_decision": row["current_decision"],
                "context_only_scorer_decision": context_decision,
                "family_aware_shadow_decision": shadow_decision,
                "family_aware_changes_current": shadow_decision != row["current_decision"],
                "family_aware_changes_context_only": shadow_decision != context_decision,
                "picked_created_by_repeat": False,
                "answer_blind": True,
            }
        )

    policy_fields = list(decision_rows[0])
    decision_path = output_dir / "answer_blind_policy_decisions.csv"
    write_csv(decision_path, decision_rows, policy_fields)
    decision_hash = sha256_file(decision_path)

    route_dependent = {
        key(row) for row in repeat_rows if as_bool(row["decision_or_level_changed"])
    }
    context_rows = apply_context_only_baseline(candidate_rows, route_dependent)
    family_rows = apply_shadow_high(context_rows, decision_rows)
    rank_only_rows = apply_shadow_high(
        context_rows, decision_rows, open_rank_only=True
    )

    # Benchmark data is intentionally loaded only after policy decisions are frozen.
    answers = load_answer_rows(answers_path)
    annotated_rows = []
    for row in decision_rows:
        status, answer = benchmark_status(
            str(row["patient_id"]), str(row["organism_name"]), answers
        )
        annotated_rows.append(
            {
                **row,
                "benchmark_status": status,
                "matched_answer": answer,
                "benchmark_join_stage": "post_policy_freeze",
            }
        )
    write_csv(
        output_dir / "benchmark_comparison.csv",
        annotated_rows,
        list(annotated_rows[0]),
    )

    current_report = report_predictions(integrated_rows)
    route_dependent_report_keys = {
        key(row)
        for row in integrated_rows
        if key(row) in route_dependent and as_bool(row["selected_for_complete_report"])
    }
    family_eligible_keys = {
        key(row) for row in decision_rows if as_bool(row["family_aware_high_eligible"])
    }
    context_report_rows = [
        row for row in integrated_rows if key(row) not in route_dependent_report_keys
    ]
    # Existing reporting routes are retained only when the generalized family-aware
    # policy re-establishes their scorer High.  New High never auto-implies Possible.
    family_report_rows = [
        row
        for row in integrated_rows
        if key(row) not in (route_dependent_report_keys - family_eligible_keys)
    ]

    metrics = {
        "current_full": scorer_metrics(candidate_rows, hospital_rows, answers),
        "context_only": scorer_metrics(context_rows, hospital_rows, answers),
        "family_aware_repeat_high": scorer_metrics(
            family_rows, hospital_rows, answers
        ),
        "rank_only_open_repeat_high_negative_control": scorer_metrics(
            rank_only_rows, hospital_rows, answers
        ),
    }
    report_metrics = {
        "current_full": evaluate_set(current_report, answers),
        "context_only": evaluate_set(report_predictions(context_report_rows), answers),
        "family_aware_repeat_high": evaluate_set(
            report_predictions(family_report_rows), answers
        ),
    }

    changed_current = [
        row for row in annotated_rows if as_bool(row["family_aware_changes_current"])
    ]
    restored_from_context = [
        row
        for row in annotated_rows
        if as_bool(row["family_aware_changes_context_only"])
    ]
    rank_only_new = [
        row
        for row in annotated_rows
        if as_bool(row["rank_only_strength_eligible"])
        and row["context_only_scorer_decision"]
        not in {"picked_shadow", "review_high_priority"}
    ]
    assertions = {
        "answer_blind_policy_has_no_benchmark_columns": not any(
            "benchmark" in field or "answer" in field
            for field in policy_fields
            if field != "answer_blind"
        ),
        "all_48_repeat_candidates_audited": len(decision_rows) == 48,
        "technical_repeat_never_creates_picked": not any(
            as_bool(row["picked_created_by_repeat"]) for row in decision_rows
        ),
        "same_library_rerun_not_high_eligible": not any(
            row["repeat_group_type"] == "explicit_original_library_rerun"
            and as_bool(row["family_aware_high_eligible"])
            for row in decision_rows
        ),
        "family_aware_shadow_preserves_current_scorer_metrics": (
            metrics["family_aware_repeat_high"] == metrics["current_full"]
        ),
        "family_aware_shadow_preserves_current_complete_report_metrics": (
            report_metrics["family_aware_repeat_high"]
            == report_metrics["current_full"]
        ),
        "no_family_aware_current_tier_changes": not changed_current,
        "no_repeat_only_picked_change": (
            metrics["family_aware_repeat_high"]["picked"]
            == metrics["current_full"]["picked"]
        ),
    }

    result = {
        "schema_version": "kh_family_aware_technical_repeat_high_shadow.v1",
        "scope": (
            "Non-mutating answer-blind policy shadow followed by retrospective "
            "benchmark evaluation; not external clinical validation."
        ),
        "inputs": {
            "repeat_audit": str(repeat_audit_path.resolve()),
            "compact_packets": str(compact_dir.resolve()),
            "scorer": str(scorer_dir.resolve()),
            "integrated": str(integrated_path.resolve()),
            "answers_posthoc_only": str(answers_path.resolve()),
        },
        "policy": policy,
        "answer_blind_policy_decisions_sha256": decision_hash,
        "counts": {
            "technical_repeat_candidates": len(decision_rows),
            "rank_only_strength_eligible": sum(
                as_bool(row["rank_only_strength_eligible"])
                for row in decision_rows
            ),
            "family_aware_high_eligible": sum(
                as_bool(row["family_aware_high_eligible"])
                for row in decision_rows
            ),
            "family_aware_eligible_by_route": dict(
                sorted(
                    Counter(
                        row["policy_reason"]
                        for row in decision_rows
                        if as_bool(row["family_aware_high_eligible"])
                    ).items()
                )
            ),
            "blocked_by_reason": dict(
                sorted(
                    Counter(
                        row["policy_reason"]
                        for row in decision_rows
                        if not as_bool(row["family_aware_high_eligible"])
                    ).items()
                )
            ),
            "changes_vs_current": len(changed_current),
            "restored_high_vs_context_only": len(restored_from_context),
            "rank_only_new_high_vs_context_only": len(rank_only_new),
            "restored_high_candidates": [
                {
                    "patient_id": row["patient_id"],
                    "organism_name": row["organism_name"],
                    "benchmark_status": row["benchmark_status"],
                    "policy_reason": row["policy_reason"],
                }
                for row in restored_from_context
            ],
        },
        "scorer_metrics": metrics,
        "complete_report_metrics": report_metrics,
        "interpretation": {
            "current_comparison": (
                "The family-aware shadow is parity-preserving against the current "
                "scorer and complete report; it makes the permissive repeat route explicit."
            ),
            "context_only_comparison": (
                "The generalized opportunistic-family plus host-support route restores "
                "P9 PJP without an organism-specific answer rule."
            ),
            "picked": (
                "No Picked result is created by technical repeat; independent evidence "
                "remains required."
            ),
            "adoption": (
                "Safe to retain as a shadow policy; production replacement still requires "
                "laboratory confirmation of technical-branch independence and external validation."
            ),
        },
        "assertions": assertions,
        "assertion_pass_count": sum(assertions.values()),
        "assertion_total_count": len(assertions),
        "all_assertions_pass": all(assertions.values()),
    }
    write_json(output_dir / "summary.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the KH family-aware technical-repeat High shadow"
    )
    parser.add_argument("--repeat-audit", type=Path, default=DEFAULT_REPEAT_AUDIT)
    parser.add_argument("--compact-dir", type=Path, default=DEFAULT_COMPACT)
    parser.add_argument("--scorer-dir", type=Path, default=DEFAULT_SCORER)
    parser.add_argument("--integrated", type=Path, default=DEFAULT_INTEGRATED)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit(
        args.repeat_audit,
        args.compact_dir,
        args.scorer_dir,
        args.integrated,
        args.answers,
        args.output_dir,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

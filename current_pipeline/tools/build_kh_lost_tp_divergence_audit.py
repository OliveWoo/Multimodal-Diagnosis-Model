"""Explain true positives picked by legacy B but missed by strict new scorer C.

This is a post-hoc, answer-aware audit. It never changes scorer decisions.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from tools import deterministic_mngs_max_scorer as legacy
from tools.evaluate_multi_assay_candidate_entry import load_answer_rows
from tools.pathogen_normalization import canonical_key
from tools.recalculate_kh_benchmark_metrics import match_type
from tools.run_legacy_scorer_on_multi_assay_shadow import (
    build_legacy_ranked_payload,
    configure_taxonomy_mode,
    picked_names,
)


DEFAULT_B_ROOT = Path(
    "outputs/runs/2026-09-28_KH_ablation_B_legacy_v20_multi_assay_stable_v1"
)
DEFAULT_C_ROOT = Path("outputs/runs/2026-09-28_KH_ablation_C_new_scorer_v1")
DEFAULT_COMPACT_ROOT = Path(
    "outputs/runs/2026-09-23_KH_compact_multi_assay_entry_source_contract_v5"
)
DEFAULT_PATIENT_ROOT = Path("outputs/patient_info_KH_0728_2Days")
DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics/"
    "kh_answers_clinical_revision_20260918.csv"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-29_KH_B_to_C_lost_tp_divergence_v1"
)
DEFAULT_SUMMARY_SUFFIX = "final_summary_selected_dna_20260918"
EXPECTED_LOST_TP = 22


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


def same_organism(left: Any, right: Any) -> bool:
    if not left or not right:
        return False
    if canonical_key(left) == canonical_key(right):
        return True
    return bool(
        match_type(str(left), str(right), genus_relaxed=False)
        or match_type(str(right), str(left), genus_relaxed=False)
    )


def find_candidate(
    candidates: Iterable[dict[str, Any]], organism_name: str
) -> dict[str, Any] | None:
    items = [item for item in candidates if isinstance(item, dict)]
    wanted = canonical_key(organism_name)
    for item in items:
        if canonical_key(item.get("organism_name")) == wanted:
            return item
    for item in items:
        if same_organism(item.get("organism_name"), organism_name):
            return item
    return None


def best_answer_match(organism_name: str, answers: list[str]) -> str | None:
    return next(
        (answer for answer in answers if same_organism(organism_name, answer)),
        None,
    )


def old_pick_basis(candidate: dict[str, Any]) -> str:
    source = str(candidate.get("evidence_source") or "")
    signal = str(candidate.get("mngs_signal_tier") or "")
    modules = candidate.get("module_support_summary") or {}
    evidence = candidate.get("key_evidence") or {}
    applied_rules = candidate.get("applied_rules") or []
    supported = bool(
        modules.get("supporting_modules")
        or modules.get("direct_support_modules")
        or evidence.get("support_modules")
        or candidate.get("hospital_support")
    )
    if source == "hospital_only" or any(
        "HOSPITAL_ONLY" in str(rule) for rule in applied_rules
    ):
        return "hospital_only"
    if signal.startswith("M1"):
        return "M1_plus_hospital" if supported else "M1_signal_only"
    if signal.startswith("M2"):
        return "M2_plus_hospital" if supported else "M2_signal_only"
    if signal.startswith("M3"):
        return "M3_protected_plus_hospital" if supported else "M3_protected"
    return "other_legacy_basis"


def first_divergence_stage(
    analytical: dict[str, Any] | None,
    clinical: dict[str, Any] | None,
) -> str:
    if analytical is None:
        return "candidate_absent_from_new_scorer"
    if not analytical.get("forward_to_clinical_scorer"):
        return "analytical_not_forwarded"
    if str(analytical.get("decision")) != "picked_shadow":
        return "analytical_not_picked"
    if clinical is None:
        return "clinical_candidate_missing"
    if str(clinical.get("clinical_decision")) != "picked_shadow":
        return "clinical_guardrail_demotion"
    return "after_clinical_or_output_difference"


def explanation_class(
    old_candidate: dict[str, Any],
    analytical: dict[str, Any] | None,
    clinical: dict[str, Any] | None,
    promotion: dict[str, Any] | None,
) -> str:
    stage = first_divergence_stage(analytical, clinical)
    if stage == "candidate_absent_from_new_scorer":
        return "source_or_candidate_gap"
    if stage == "analytical_not_forwarded":
        return "new_analytical_threshold_not_met"
    if stage == "clinical_guardrail_demotion":
        return "clinical_guardrail_overrode_analytical_pick"

    profile = (analytical or {}).get("analytical_profile") or {}
    hospital = (analytical or {}).get("hospital_profile") or {}
    family = str((analytical or {}).get("taxonomy_family") or "")
    direct = hospital.get("best_direct_level")
    cross_molecule = bool(profile.get("cross_molecule_selected"))
    repeated = bool(profile.get("reproducibility_axis"))
    basis = old_pick_basis(old_candidate)
    guarded = {
        "environmental_low_specificity",
        "respiratory_colonizer",
        "commensal_opportunist",
        "fungal_colonizer",
    }
    if family in guarded and not direct and not cross_molecule:
        return "taxonomy_family_requires_independent_evidence"
    if basis in {"M1_signal_only", "M2_signal_only"} and not direct:
        return "legacy_best_test_signal_not_independent_clinical_axis"
    if promotion and (promotion.get("promotion_gate") or {}).get("blockers"):
        return "promotion_independent_axis_not_satisfied"
    if repeated and not direct and not cross_molecule:
        return "same_specimen_repeat_counts_as_one_analytical_axis"
    return "new_multi_axis_pick_threshold_not_met"


def short_explanation(
    old_candidate: dict[str, Any],
    analytical: dict[str, Any] | None,
    clinical: dict[str, Any] | None,
    promotion: dict[str, Any] | None,
) -> str:
    basis = old_pick_basis(old_candidate)
    family = str((analytical or {}).get("taxonomy_family") or "unknown")
    stage = first_divergence_stage(analytical, clinical)
    reason = explanation_class(old_candidate, analytical, clinical, promotion)
    return (
        f"Legacy reached Picked from {basis} using the best per-test signal; "
        f"the new scorer stopped at {stage}, taxonomy={family}, "
        f"primary explanation={reason}."
    )


def reconstruct_legacy_scores(
    compact_root: Path,
    patient_root: Path,
    summary_suffix: str,
) -> dict[int, dict[str, Any]]:
    configure_taxonomy_mode("base-only")
    scores: dict[int, dict[str, Any]] = {}
    packet_paths = sorted(
        (compact_root / "patient_packets").glob(
            "NGS_patient_*_compact_entry_shadow.json"
        ),
        key=patient_number,
    )
    if not packet_paths:
        raise ValueError(f"No compact packets under {compact_root}")
    for packet_path in packet_paths:
        patient_id = patient_number(packet_path)
        patient_dir = patient_root / f"NGS_patient_{patient_id}_json"
        summary_path = (
            patient_dir
            / "summary_outputs"
            / f"NGS_patient_{patient_id}_{summary_suffix}.json"
        )
        ranked, provenance = build_legacy_ranked_payload(read_json(packet_path))
        score = legacy.score_payload(
            patient_dir=patient_dir,
            ranked_mngs=ranked,
            final_summary=read_json(summary_path),
        )
        for candidate in score.get("pathogen_candidates") or []:
            if isinstance(candidate, dict):
                candidate["compatibility_source_tests"] = provenance.get(
                    canonical_key(candidate.get("organism_name")), []
                )
        scores[patient_id] = score
    return scores


def frozen_b_picks(path: Path) -> dict[int, set[str]]:
    grouped: dict[int, set[str]] = {}
    for row in read_csv(path):
        patient_id = int(row["patient_id"])
        grouped.setdefault(patient_id, set())
        if str(row.get("picked") or "").casefold() in {"true", "1", "yes"}:
            grouped[patient_id].add(canonical_key(row.get("organism_name")))
    return grouped


def verify_legacy_reconstruction(
    scores: dict[int, dict[str, Any]],
    frozen: dict[int, set[str]],
) -> None:
    mismatches: list[str] = []
    for patient_id in sorted(set(scores) | set(frozen)):
        rebuilt = picked_names(scores.get(patient_id, {}))
        expected = frozen.get(patient_id, set())
        if rebuilt != expected:
            mismatches.append(
                f"P{patient_id}: rebuilt={sorted(rebuilt)} frozen={sorted(expected)}"
            )
    if mismatches:
        raise AssertionError(
            "Legacy reconstruction differs from frozen B:\n" + "\n".join(mismatches)
        )


def patient_stage_payload(
    root: Path,
    subdir: str,
    patient_id: int,
    suffix: str,
) -> dict[str, Any]:
    return read_json(
        root
        / subdir
        / "patient_outputs"
        / f"NGS_patient_{patient_id}_{suffix}.json"
    )


def stage_candidates(
    payload: dict[str, Any],
    *,
    analytical: bool = False,
) -> list[dict[str, Any]]:
    if analytical:
        return [
            item
            for item in (
                (payload.get("patient_organism_decisions") or [])
                + (payload.get("hospital_only_decisions") or [])
            )
            if isinstance(item, dict)
        ]
    return [
        item
        for item in payload.get("all_forwarded_candidates") or []
        if isinstance(item, dict)
    ]


def compact_list(value: Any) -> str:
    if not value:
        return ""
    if isinstance(value, list):
        return "|".join(str(item) for item in value)
    return str(value)


def proposed_disposition(
    analytical: dict[str, Any] | None,
    possible: dict[str, Any] | None,
) -> tuple[str, str, str]:
    role = str((possible or {}).get("possible_reporting_role") or "")
    selected = bool((possible or {}).get("selected_for_complete_report"))
    family = str((analytical or {}).get("taxonomy_family") or "")
    if selected and role == "strict_picked_report":
        return (
            "Picked",
            "already_resolved_by_existing_answer_blind_promotion",
            "no_incremental_risk_from_this_audit",
        )
    if selected:
        return (
            "Possible",
            "already_resolved_by_existing_answer_blind_possible_route",
            "existing_route_requires_separate_18_fp_family_audit",
        )
    if family == "high_consequence_opportunistic":
        return (
            "Context_pending_verified_host_or_direct_axis",
            "source_contract_or_verified_evidence_review_required",
            "unknown_until_rule_is_tested_on_all_same_family_candidates",
        )
    return (
        "Context",
        "no_safe_generic_upgrade_identified",
        "high_if_family_guardrail_is_blanket_relaxed",
    )


def candidate_record(
    patient_id: int,
    answer: str,
    old_candidate: dict[str, Any],
    analytical: dict[str, Any] | None,
    clinical: dict[str, Any] | None,
    promotion: dict[str, Any] | None,
    possible: dict[str, Any] | None,
) -> dict[str, Any]:
    analytical_profile = (analytical or {}).get("analytical_profile") or {}
    hospital_profile = (analytical or {}).get("hospital_profile") or {}
    promotion_gate = (promotion or {}).get("promotion_gate") or {}
    possible_gate = (possible or {}).get("possible_pathogen_gate") or {}
    disposition, generic_resolution, fp_risk = proposed_disposition(
        analytical, possible
    )
    return {
        "patient_id": patient_id,
        "answer_name": answer,
        "legacy_organism_name": old_candidate.get("organism_name"),
        "new_organism_name": (analytical or {}).get("organism_name"),
        "legacy_pick_basis": old_pick_basis(old_candidate),
        "legacy_integrated_level": old_candidate.get("integrated_causative_level"),
        "legacy_mngs_signal_tier": old_candidate.get("mngs_signal_tier"),
        "legacy_rank_priority": old_candidate.get("rank_priority"),
        "legacy_reads": old_candidate.get("reads"),
        "legacy_reads_tier": old_candidate.get("reads_tier"),
        "legacy_reads_percentile": old_candidate.get("reads_percentile"),
        "legacy_dominance_tier": old_candidate.get("dominance_tier"),
        "legacy_specimen_class": old_candidate.get("specimen_class"),
        "legacy_applied_rules": compact_list(old_candidate.get("applied_rules")),
        "legacy_support_modules": compact_list(
            (old_candidate.get("key_evidence") or {}).get("support_modules")
        ),
        "legacy_key_evidence": json.dumps(
            old_candidate.get("key_evidence") or {},
            ensure_ascii=False,
            sort_keys=True,
        ),
        "legacy_source_test_count": len(
            old_candidate.get("compatibility_source_tests") or []
        ),
        "new_taxonomy_family": (analytical or {}).get("taxonomy_family"),
        "new_specimen_context": (analytical or {}).get("specimen_context"),
        "new_analytical_decision": (analytical or {}).get("decision"),
        "new_analytical_level": (analytical or {}).get("integrated_level"),
        "new_forwarded_to_clinical": bool(
            (analytical or {}).get("forward_to_clinical_scorer")
        ),
        "new_best_rank": analytical_profile.get(
            "best_rank_in_retained_universe"
        ),
        "new_selected_positive_test_count": analytical_profile.get(
            "selected_positive_test_count"
        ),
        "new_cross_molecule": bool(
            analytical_profile.get("cross_molecule_selected")
        ),
        "new_reproducibility_axis": bool(
            analytical_profile.get("reproducibility_axis")
        ),
        "new_per_test_signals_json": json.dumps(
            analytical_profile.get("per_test_signals") or [],
            ensure_ascii=False,
            sort_keys=True,
        ),
        "new_direct_hospital_level": hospital_profile.get("best_direct_level"),
        "new_direct_modules_level_1_2": compact_list(
            hospital_profile.get("direct_support_modules_level_1_2")
        ),
        "new_direct_modules_level_1_3": compact_list(
            hospital_profile.get("direct_support_modules_level_1_3")
        ),
        "new_hospital_evidence_detail_json": json.dumps(
            hospital_profile.get("hospital_evidence_detail") or {},
            ensure_ascii=False,
            sort_keys=True,
        ),
        "new_analytical_rules": compact_list(
            (analytical or {}).get("rule_ids")
        ),
        "new_analytical_reasons": compact_list(
            (analytical or {}).get("reasons")
        ),
        "new_clinical_decision": (clinical or {}).get("clinical_decision"),
        "new_clinical_level": (clinical or {}).get("clinical_level"),
        "new_history_routed_decision": (clinical or {}).get(
            "history_routed_decision"
        ),
        "new_history_rule_ids": compact_list(
            (clinical or {}).get("history_rule_ids")
        ),
        "new_history_evidence_ids": compact_list(
            (clinical or {}).get("history_evidence_ids")
        ),
        "new_consumed_evidence_ids": compact_list(
            (clinical or {}).get("consumed_evidence_ids")
        ),
        "new_host_evidence_count": len(
            (clinical or {}).get("host_evidence") or []
        ),
        "new_image_evidence_count": len(
            (clinical or {}).get("image_evidence") or []
        ),
        "new_clinical_rules": compact_list((clinical or {}).get("rule_ids")),
        "new_clinical_reasons": compact_list((clinical or {}).get("reasons")),
        "new_promotion_decision": (promotion or {}).get("clinical_decision"),
        "new_promotion_outcome": promotion_gate.get("outcome"),
        "new_promotion_route": promotion_gate.get("route"),
        "new_promotion_blockers": compact_list(promotion_gate.get("blockers")),
        "new_promotion_axes_json": json.dumps(
            promotion_gate.get("axes") or {},
            ensure_ascii=False,
            sort_keys=True,
        ),
        "new_possible_role": (possible or {}).get("possible_reporting_role"),
        "new_possible_selected": bool(
            (possible or {}).get("selected_for_complete_report")
        ),
        "new_possible_route": possible_gate.get("route"),
        "new_possible_blockers": compact_list(possible_gate.get("blockers")),
        "first_divergence_stage": first_divergence_stage(
            analytical, clinical
        ),
        "explanation_class": explanation_class(
            old_candidate, analytical, clinical, promotion
        ),
        "proposed_disposition": disposition,
        "generic_answer_blind_resolution": generic_resolution,
        "fp_risk_if_generalized": fp_risk,
        "one_line_explanation": short_explanation(
            old_candidate, analytical, clinical, promotion
        ),
    }


def aggregate_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    dimensions = [
        "legacy_pick_basis",
        "new_taxonomy_family",
        "first_divergence_stage",
        "new_clinical_decision",
        "new_promotion_decision",
        "new_possible_role",
        "explanation_class",
        "proposed_disposition",
    ]
    rows: list[dict[str, Any]] = []
    for dimension in dimensions:
        counts = Counter(str(record.get(dimension) or "missing") for record in records)
        rows.extend(
            {
                "dimension": dimension,
                "value": value,
                "count": count,
            }
            for value, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        )
    return rows


def report_markdown(
    records: list[dict[str, Any]],
    aggregates: list[dict[str, Any]],
) -> str:
    lines = [
        "# KH legacy B to new C lost-TP divergence audit",
        "",
        "This is an answer-aware post-hoc audit. It explains scorer divergence",
        "and is not an input to either scorer.",
        "",
        f"Lost true positives audited: {len(records)}",
        "",
        "## Aggregate counts",
        "",
        "| Dimension | Value | Count |",
        "|---|---|---:|",
    ]
    for row in aggregates:
        lines.append(
            f"| {row['dimension']} | {row['value']} | {row['count']} |"
        )
    lines.extend(
        [
            "",
            "## Per-case divergence",
            "",
            "| Patient | Answer | Legacy basis | New family | First divergence | "
            "Clinical tier | Possible role | Proposed disposition | "
            "Primary explanation |",
            "|---:|---|---|---|---|---|---|---|---|",
        ]
    )
    for row in records:
        lines.append(
            "| P{patient_id} | {answer_name} | {legacy_pick_basis} | "
            "{new_taxonomy_family} | {first_divergence_stage} | "
            "{new_clinical_decision} | {new_possible_role} | "
            "{proposed_disposition} | "
            "{explanation_class} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "- Legacy B is reconstructed with base-only taxonomy and checked against "
            "the frozen B candidate CSV before the audit is accepted.",
            "- A benchmark answer is used only to select the lost-TP review set.",
            "- No patient-specific exception or scorer threshold is created here.",
            "- The next step is to test general rules against both the 22 recovered "
            "TPs and the 18 FPs removed by C.",
            "",
        ]
    )
    return "\n".join(lines)


def build_audit(
    b_root: Path,
    c_root: Path,
    compact_root: Path,
    patient_root: Path,
    answers_path: Path,
    output_dir: Path,
    *,
    summary_suffix: str = DEFAULT_SUMMARY_SUFFIX,
    expected_lost_tp: int = EXPECTED_LOST_TP,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")

    scores = reconstruct_legacy_scores(
        compact_root, patient_root, summary_suffix
    )
    frozen = frozen_b_picks(
        b_root / "legacy_multi_assay_candidate_decisions.csv"
    )
    verify_legacy_reconstruction(scores, frozen)
    answers = load_answer_rows(answers_path)
    records: list[dict[str, Any]] = []

    for patient_id in sorted(answers):
        analytical_payload = patient_stage_payload(
            c_root,
            "scorer",
            patient_id,
            "test_aware_deterministic_shadow",
        )
        clinical_payload = patient_stage_payload(
            c_root,
            "clinical",
            patient_id,
            "test_aware_clinical_shadow",
        )
        promotion_payload = patient_stage_payload(
            c_root,
            "promotion_v4",
            patient_id,
            "test_aware_clinical_shadow",
        )
        possible_payload = patient_stage_payload(
            c_root,
            "possible_pathogen_v1",
            patient_id,
            "test_aware_possible_pathogen_shadow",
        )
        analytical_candidates = stage_candidates(
            analytical_payload, analytical=True
        )
        clinical_candidates = stage_candidates(clinical_payload)
        promotion_candidates = stage_candidates(promotion_payload)
        possible_candidates = stage_candidates(possible_payload)
        strict_new_picks = [
            item
            for item in clinical_candidates
            if item.get("clinical_decision") == "picked_shadow"
        ]
        old_score = scores[patient_id]
        old_picked_keys = picked_names(old_score)
        for old_candidate in old_score.get("pathogen_candidates") or []:
            if not isinstance(old_candidate, dict):
                continue
            old_name = str(old_candidate.get("organism_name") or "")
            if canonical_key(old_name) not in old_picked_keys:
                continue
            answer = best_answer_match(old_name, answers[patient_id])
            if not answer:
                continue
            if find_candidate(strict_new_picks, old_name):
                continue
            analytical = find_candidate(analytical_candidates, old_name)
            clinical = find_candidate(clinical_candidates, old_name)
            promotion = find_candidate(promotion_candidates, old_name)
            possible = find_candidate(possible_candidates, old_name)
            records.append(
                candidate_record(
                    patient_id,
                    answer,
                    old_candidate,
                    analytical,
                    clinical,
                    promotion,
                    possible,
                )
            )

    records.sort(
        key=lambda item: (
            int(item["patient_id"]),
            canonical_key(item["legacy_organism_name"]),
        )
    )
    if len(records) != expected_lost_tp:
        raise AssertionError(
            f"Expected {expected_lost_tp} lost TPs, found {len(records)}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    aggregates = aggregate_rows(records)
    summary = {
        "schema_version": "kh_lost_tp_divergence_audit.v1",
        "answer_aware_posthoc_audit": True,
        "scorer_rules_changed": False,
        "legacy_reconstruction_verified_against_frozen_b": True,
        "lost_tp_count": len(records),
        "patient_count": len({row["patient_id"] for row in records}),
        "legacy_pick_basis": dict(
            Counter(row["legacy_pick_basis"] for row in records)
        ),
        "first_divergence_stage": dict(
            Counter(row["first_divergence_stage"] for row in records)
        ),
        "explanation_class": dict(
            Counter(row["explanation_class"] for row in records)
        ),
        "new_taxonomy_family": dict(
            Counter(str(row["new_taxonomy_family"]) for row in records)
        ),
        "new_possible_selected_count": sum(
            bool(row["new_possible_selected"]) for row in records
        ),
        "proposed_disposition": dict(
            Counter(row["proposed_disposition"] for row in records)
        ),
        "sources": {
            "b_root": str(b_root.resolve()),
            "c_root": str(c_root.resolve()),
            "compact_root": str(compact_root.resolve()),
            "patient_root": str(patient_root.resolve()),
            "answers": str(answers_path.resolve()),
        },
    }
    write_csv(output_dir / "lost_tp_divergence.csv", records)
    write_json(output_dir / "lost_tp_divergence.json", records)
    write_csv(output_dir / "aggregate_counts.csv", aggregates)
    write_json(output_dir / "summary.json", summary)
    (output_dir / "lost_tp_divergence_report.md").write_text(
        report_markdown(records, aggregates),
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--b-root", type=Path, default=DEFAULT_B_ROOT)
    parser.add_argument("--c-root", type=Path, default=DEFAULT_C_ROOT)
    parser.add_argument("--compact-root", type=Path, default=DEFAULT_COMPACT_ROOT)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-suffix", default=DEFAULT_SUMMARY_SUFFIX)
    parser.add_argument("--expected-lost-tp", type=int, default=EXPECTED_LOST_TP)
    args = parser.parse_args()
    summary = build_audit(
        args.b_root,
        args.c_root,
        args.compact_root,
        args.patient_root,
        args.answers,
        args.output_dir,
        summary_suffix=args.summary_suffix,
        expected_lost_tp=args.expected_lost_tp,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

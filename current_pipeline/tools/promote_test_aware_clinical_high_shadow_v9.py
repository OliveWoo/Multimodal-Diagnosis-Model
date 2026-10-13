"""Add pathogen-specific strong-evidence Picked gates on top of promotion v8.

V9 preserves every v4-v8 route.  New routes are intentionally narrow: CMV
requires pathogen-specific clinical convergence, and an Aspergillus group label
requires direct invasive proof before it can be Picked.  Missing evidence is
unknown and never satisfies either route.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from tools import promote_test_aware_clinical_high_shadow as base
from tools import promote_test_aware_clinical_high_shadow_v7 as v7
from tools import promote_test_aware_clinical_high_shadow_v8 as v8
from tools.pathogen_normalization import canonical_key


DEFAULT_POLICY = Path(
    "rules/test_aware_clinical_promotion_v9_pathogen_specific_convergence.json"
)


def _reactivation_profile(candidate: dict[str, Any]) -> dict[str, bool]:
    profile = candidate.get("reactivation_interpretation") or {}
    return {
        "explicit_current": bool(profile.get("has_explicit_current_reactivation")),
        "quantitative_or_tissue": bool(profile.get("has_quantitative_or_tissue_evidence")),
        "compatible_syndrome": bool(profile.get("has_compatible_pulmonary_syndrome")),
    }


def _cmv_picked_gate(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> tuple[bool, dict[str, Any]]:
    config = policy.get("cmv_picked_convergence") or {}
    axes = base.evidence_axes(candidate, policy)
    taxonomy = candidate.get("taxonomy_profile") or {}
    reactivation = _reactivation_profile(candidate)
    key = canonical_key(candidate.get("organism_name"))
    rank = axes.get("best_rank_in_retained_universe")
    selected = int(axes.get("selected_positive_test_count") or 0)
    direct_level = axes.get("direct_hospital_level")
    direct_aligned = bool(
        (axes.get("direct_evidence_timing_profile") or {}).get(
            "event_aligned_positive"
        )
    )
    common = bool(
        config.get("enabled")
        and key in set(config.get("organism_keys") or [])
        and str(taxonomy.get("mapping_status") or "") == "exact_species"
        and candidate.get("candidate_source") == "mngs_positive_reads"
        and axes.get("specimen_context")
        in set(config.get("compatible_specimen_contexts") or [])
        and rank is not None
        and direct_level is not None
        and direct_level <= int(config.get("maximum_direct_hospital_level", 2))
        and direct_aligned
        and reactivation["explicit_current"]
    )
    analytical_convergence = bool(
        common
        and rank <= int(config.get("maximum_rank", 3))
        and selected >= int(config.get("minimum_positive_tests", 2))
        and axes.get("reproducibility_axis")
        and (
            axes.get("cross_molecule_selected")
            or not config.get("require_cross_molecule", True)
        )
        and (
            reactivation["quantitative_or_tissue"]
            or reactivation["compatible_syndrome"]
        )
    )
    single_with_invasive_evidence = bool(
        common
        and config.get("allow_single_test_with_quantitative_or_tissue_evidence")
        and rank <= int(config.get("single_test_maximum_rank", 1))
        and selected >= 1
        and reactivation["quantitative_or_tissue"]
    )
    return analytical_convergence or single_with_invasive_evidence, {
        **axes,
        "reactivation_interpretation": reactivation,
        "analytical_convergence": analytical_convergence,
        "single_test_with_invasive_evidence": single_with_invasive_evidence,
        "missing_confirmation_is_negative": False,
    }


def _aspergillus_group_picked_gate(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> tuple[bool, dict[str, Any]]:
    config = policy.get("aspergillus_group_picked_convergence") or {}
    axes = base.evidence_axes(candidate, policy)
    taxonomy = candidate.get("taxonomy_profile") or {}
    detail = candidate.get("hospital_evidence_detail") or (
        (candidate.get("hospital_profile") or {}).get("hospital_evidence_detail")
    ) or {}
    text = f" {base._recursive_text(detail)} "
    proof_terms = [
        str(term).casefold()
        for term in config.get("invasive_proof_terms") or []
        if str(term).casefold() in text
    ]
    excluded_only_terms = [
        str(term).casefold()
        for term in config.get("insufficient_alone_terms") or []
        if str(term).casefold() in text
    ]
    key = canonical_key(candidate.get("organism_name"))
    direct_level = axes.get("direct_hospital_level")
    direct_aligned = bool(
        (axes.get("direct_evidence_timing_profile") or {}).get(
            "event_aligned_positive"
        )
    )
    qualifies = bool(
        config.get("enabled")
        and key in set(config.get("group_organism_keys") or [])
        and str(taxonomy.get("taxonomic_rank") or "") == "genus_or_group_label"
        and direct_level is not None
        and direct_level <= int(config.get("maximum_direct_hospital_level", 2))
        and direct_aligned
        and proof_terms
        and not (
            excluded_only_terms
            and not any("histopath" in term or "tissue" in term for term in proof_terms)
        )
    )
    return qualifies, {
        **axes,
        "species_resolution": "species_unresolved",
        "matched_invasive_proof_terms": proof_terms,
        "matched_insufficient_alone_terms": excluded_only_terms,
        "missing_invasive_proof_is_negative": False,
    }


def promotion_decision(
    candidate: dict[str, Any],
    policy: dict[str, Any],
    patient_context: dict[str, Any] | None = None,
) -> tuple[str, dict[str, Any], list[str], list[str]]:
    decision, gate, rules, reasons = v8.promotion_decision(
        candidate, policy, patient_context
    )
    if decision == "picked_shadow":
        return decision, gate, rules, reasons

    cmv_qualifies, cmv_axes = _cmv_picked_gate(candidate, policy)
    if cmv_qualifies:
        return (
            "picked_shadow",
            {
                "evaluated": True,
                "outcome": "promoted_to_picked_shadow",
                "route": "cmv_pathogen_specific_clinical_convergence",
                "axes": cmv_axes,
                "blockers": [],
                "cautions": [
                    "CMV identity, reactivation, and organ-disease attribution remain separate claims; this route requires evidence for all configured axes."
                ],
            },
            rules + [
                "PROMO-CMV1-EXACT-IDENTITY",
                "PROMO-CMV2-EVENT-ALIGNED-DIRECT-SUPPORT",
                "PROMO-CMV3-CURRENT-REACTIVATION-CONVERGENCE",
            ],
            reasons + [
                "Exact CMV mNGS evidence converges with event-aligned pathogen-specific clinical evidence and current reactivation or invasive-disease support."
            ],
        )

    aspergillus_qualifies, aspergillus_axes = _aspergillus_group_picked_gate(
        candidate, policy
    )
    if aspergillus_qualifies:
        return (
            "picked_shadow",
            {
                "evaluated": True,
                "outcome": "promoted_to_picked_shadow",
                "route": "aspergillus_group_direct_invasive_proof",
                "axes": aspergillus_axes,
                "blockers": [],
                "cautions": [
                    "The pathogen is reported as Aspergillus spp. (species unresolved); the rule does not infer a species name."
                ],
            },
            rules + [
                "PROMO-ASP1-GROUP-IDENTITY-EXPLICIT",
                "PROMO-ASP2-EVENT-ALIGNED-INVASIVE-PROOF",
                "PROMO-ASP3-SPECIES-UNRESOLVED",
            ],
            reasons + [
                "Event-aligned direct invasive evidence supports Aspergillus disease while species identity remains unresolved."
            ],
        )

    return decision, gate, rules, reasons


def promote_candidate(
    candidate: dict[str, Any],
    policy: dict[str, Any],
    patient_context: dict[str, Any] | None = None,
) -> dict[str, Any]:
    original = base.promotion_decision
    try:
        base.promotion_decision = promotion_decision
        return base.promote_candidate(candidate, policy, patient_context)
    finally:
        base.promotion_decision = original


def run(
    input_root: Path,
    output_dir: Path,
    *,
    policy_path: Path = DEFAULT_POLICY,
    patients: set[int] | None = None,
    rpm_threshold: float | None = None,
    qc_mode: str | None = None,
) -> dict[str, Any]:
    policy = v7._effective_policy(
        policy_path, rpm_threshold=rpm_threshold, qc_mode=qc_mode
    )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    effective_path = output_dir.with_name(f"{output_dir.name}_effective_policy.json")
    if effective_path.exists():
        raise FileExistsError(f"Refusing to overwrite {effective_path}")
    effective_text = json.dumps(policy, ensure_ascii=False, indent=2) + "\n"
    effective_path.write_text(effective_text, encoding="utf-8")

    original_decision = base.promotion_decision
    original_context = base.build_patient_promotion_context
    try:
        base.promotion_decision = promotion_decision
        base.build_patient_promotion_context = v8.build_patient_promotion_context
        summary = base.run(
            input_root,
            output_dir,
            policy_path=effective_path,
            patients=patients,
        )
    finally:
        base.promotion_decision = original_decision
        base.build_patient_promotion_context = original_context
    (output_dir / "effective_policy.json").write_text(
        effective_text, encoding="utf-8"
    )
    effective_path.unlink()
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--patients", nargs="*", type=int)
    parser.add_argument("--rpm-threshold", type=float)
    parser.add_argument("--qc-mode", choices=["mixed", "strict"])
    args = parser.parse_args()
    print(json.dumps(run(
        args.input_root,
        args.output_dir,
        policy_path=args.policy,
        patients=set(args.patients) if args.patients else None,
        rpm_threshold=args.rpm_threshold,
        qc_mode=args.qc_mode,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

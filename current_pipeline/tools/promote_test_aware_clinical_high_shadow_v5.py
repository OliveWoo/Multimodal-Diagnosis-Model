"""Apply the v5 targeted-respiratory-assay bridge on top of frozen promotion v4.

The frozen v4 implementation remains byte-stable.  This opt-in shadow adds a
single narrow route for an exact-species bacterial candidate when a positive,
event-aligned targeted respiratory observation was preserved but the generic
hospital-level interface did not expose it as Level 1/2.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from tools import promote_test_aware_clinical_high_shadow as base


DEFAULT_POLICY = Path("rules/test_aware_clinical_promotion_v5_targeted_assay_bridge.json")
_V4_PROMOTION_DECISION = base.promotion_decision


def promotion_decision(
    candidate: dict[str, Any], policy: dict[str, Any],
    patient_context: dict[str, Any] | None = None,
) -> tuple[str, dict[str, Any], list[str], list[str]]:
    decision, gate, rules, reasons = _V4_PROMOTION_DECISION(
        candidate, policy, patient_context
    )
    bridge = policy.get("targeted_respiratory_assay_bridge") or {}
    axes = gate.get("axes") or {}
    taxonomy = candidate.get("taxonomy_profile") or {}
    hospital = candidate.get("hospital_profile") or {}
    mapping_status = str(taxonomy.get("mapping_status") or "")
    exact_or_alias = bool(hospital.get("exact_or_alias_evidence"))
    axes["taxonomy_mapping_status"] = mapping_status
    axes["exact_or_alias_hospital_evidence"] = exact_or_alias

    if (
        not bridge.get("enabled", False)
        or str(candidate.get("clinical_decision") or "") != "review_high_priority"
        or decision != "review_high_priority"
        or str(candidate.get("candidate_source") or "") != "mngs_positive_reads"
    ):
        return decision, gate, rules, reasons

    rank = axes.get("best_rank_in_retained_universe")
    bridge_ready = bool(
        axes.get("taxonomy_family") in set(bridge.get("allowed_families") or [])
        and mapping_status in set(bridge.get("allowed_mapping_statuses") or [])
        and exact_or_alias
        and axes.get("compatible_specimen")
        and (axes.get("direct_evidence_timing_profile") or {}).get(
            "event_aligned_positive"
        )
        and axes.get("reproducibility_axis")
        and axes.get("cross_molecule_selected")
        and int(axes.get("selected_positive_test_count") or 0)
            >= int(bridge.get("minimum_positive_tests", 3))
        and rank is not None
        and int(rank) <= int(bridge.get("maximum_rank", 3))
    )
    if not bridge_ready:
        return decision, gate, rules, reasons

    return (
        "picked_shadow",
        {
            "evaluated": True,
            "outcome": "promoted_to_picked_shadow",
            "route": "bacterial_targeted_respiratory_assay_bridge",
            "axes": axes,
            "blockers": [],
        },
        rules + [
            "PROMO-B3-TARGETED-RESPIRATORY-ASSAY",
            "PROMO-B4-EXACT-ALIAS-EVENT-ALIGNED",
            "PROMO-B5-CROSS-MOLECULE-TOP3",
        ],
        [
            "An exact-species candidate has a positive event-aligned targeted "
            "respiratory observation plus complete top-three cross-molecule mNGS support."
        ],
    )


def promote_candidate(
    candidate: dict[str, Any], policy: dict[str, Any],
    patient_context: dict[str, Any] | None = None,
) -> dict[str, Any]:
    original = base.promotion_decision
    try:
        base.promotion_decision = promotion_decision
        return base.promote_candidate(candidate, policy, patient_context)
    finally:
        base.promotion_decision = original


def run(
    input_root: Path, output_dir: Path, *, policy_path: Path = DEFAULT_POLICY,
    patients: set[int] | None = None,
) -> dict[str, Any]:
    original = base.promotion_decision
    try:
        base.promotion_decision = promotion_decision
        return base.run(
            input_root, output_dir, policy_path=policy_path, patients=patients
        )
    finally:
        base.promotion_decision = original


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--patients", type=int, nargs="*")
    args = parser.parse_args()
    selected = set(args.patients) if args.patients else None
    print(json.dumps(
        run(args.input_root, args.output_dir, policy_path=args.policy, patients=selected),
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()

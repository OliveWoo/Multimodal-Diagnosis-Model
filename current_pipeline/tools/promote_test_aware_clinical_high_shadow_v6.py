"""Apply targeted-assay v5 plus a narrow analytically-compelling Picked route.

The frozen v4 implementation and the opt-in v5 implementation remain
unchanged.  V6 adds a separately labelled analytical route that requires an
exact-species, top-three, cross-molecule signal in at least three tests, full
per-test normalization coverage, a minimum count of tests above a fixed RPM
threshold, and no explicit colonization/current-not-infection counterevidence.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

from tools import promote_test_aware_clinical_high_shadow as base
from tools import promote_test_aware_clinical_high_shadow_v5 as v5


DEFAULT_POLICY = Path(
    "rules/test_aware_clinical_promotion_v6_analytically_compelling_rpm5.json"
)


def _colonization_profile(candidate: dict[str, Any]) -> dict[str, Any]:
    direct = candidate.get("colonization_interpretation")
    if isinstance(direct, dict):
        return direct
    source = candidate.get("source_candidate") or {}
    nested = source.get("colonization_interpretation")
    return nested if isinstance(nested, dict) else {}


def _burden_profile(
    candidate: dict[str, Any], route: dict[str, Any]
) -> dict[str, Any]:
    analytical = candidate.get("analytical_profile") or {}
    signals = analytical.get("per_test_signals") or []
    selected_count = int(analytical.get("selected_positive_test_count") or 0)
    threshold = float(route.get("per_test_rpm_threshold", 5.0))
    accepted_qc = set(route.get("accepted_qc_statuses") or [])
    usable: list[dict[str, Any]] = []
    all_normalized = len(signals) == selected_count and selected_count > 0
    all_qc_accepted = all_normalized

    for signal in signals:
        rpm = signal.get("rpm_total")
        normalized = bool(signal.get("normalization_available")) and rpm is not None
        qc_status = str(signal.get("qc_status") or "")
        if not normalized:
            all_normalized = False
        if accepted_qc and qc_status not in accepted_qc:
            all_qc_accepted = False
        if normalized:
            usable.append({
                "seq_id": signal.get("seq_id"),
                "nucleic_type": signal.get("nucleic_type"),
                "qc_status": qc_status,
                "rpm_total": float(rpm),
                "meets_rpm_threshold": float(rpm) >= threshold,
            })

    return {
        "per_test_rpm_threshold": threshold,
        "selected_positive_test_count": selected_count,
        "per_test_signal_count": len(signals),
        "all_selected_tests_normalized": all_normalized,
        "all_selected_tests_qc_accepted": all_qc_accepted,
        "accepted_qc_statuses": sorted(accepted_qc),
        "tests_meeting_rpm_threshold": sum(
            1 for signal in usable if signal["meets_rpm_threshold"]
        ),
        "per_test_rpm": usable,
    }


def promotion_decision(
    candidate: dict[str, Any], policy: dict[str, Any],
    patient_context: dict[str, Any] | None = None,
) -> tuple[str, dict[str, Any], list[str], list[str]]:
    decision, gate, rules, reasons = v5.promotion_decision(
        candidate, policy, patient_context
    )
    route = policy.get("analytically_compelling_picked") or {}
    axes = gate.get("axes") or {}
    taxonomy = candidate.get("taxonomy_profile") or {}
    analytical = candidate.get("analytical_profile") or {}
    colonization = _colonization_profile(candidate)
    burden = _burden_profile(candidate, route)
    explicit_counterevidence = bool(
        colonization.get("has_explicit_colonization")
        or colonization.get("has_explicit_current_not_infection")
    )
    axes["analytically_compelling_burden_profile"] = burden
    axes["explicit_colonization_or_current_not_infection"] = (
        explicit_counterevidence
    )

    if (
        not route.get("enabled", False)
        or str(candidate.get("clinical_decision") or "")
            != "review_high_priority"
        or decision != "review_high_priority"
        or str(candidate.get("candidate_source") or "") != "mngs_positive_reads"
    ):
        return decision, gate, rules, reasons

    rank = analytical.get("best_rank_in_retained_universe")
    ready = bool(
        str(taxonomy.get("mapping_status") or "")
            in set(route.get("allowed_mapping_statuses") or [])
        and str(taxonomy.get("primary_rule_family") or "")
            in set(route.get("allowed_families") or [])
        and axes.get("compatible_specimen")
        and int(analytical.get("selected_positive_test_count") or 0)
            >= int(route.get("minimum_positive_tests", 3))
        and bool(analytical.get("reproducibility_axis"))
        and bool(analytical.get("cross_molecule_selected"))
        and rank is not None
        and int(rank) <= int(route.get("maximum_rank", 3))
        and burden["all_selected_tests_normalized"]
        and burden["all_selected_tests_qc_accepted"]
        and burden["tests_meeting_rpm_threshold"]
            >= int(route.get("minimum_tests_at_or_above_rpm_threshold", 2))
        and not explicit_counterevidence
    )
    if not ready:
        return decision, gate, rules, reasons

    threshold = burden["per_test_rpm_threshold"]
    return (
        "picked_shadow",
        {
            "evaluated": True,
            "outcome": "promoted_to_picked_shadow",
            "route": "bacterial_analytically_compelling_cross_molecule",
            "axes": axes,
            "blockers": [],
        },
        rules + [
            "PROMO-B6-EXACT-SPECIES-ANALYTICALLY-COMPELLING",
            "PROMO-B7-CROSS-MOLECULE-THREE-TESTS",
            "PROMO-B8-PER-TEST-RPM-THRESHOLD",
            "PROMO-B9-NO-EXPLICIT-COLONIZATION-COUNTEREVIDENCE",
        ],
        [
            "Exact-species top-three mNGS evidence is reproduced across at "
            "least three tests and DNA/RNA, with at least two per-test RPM "
            f"values at or above {threshold:g} and no explicit colonization "
            "or current-not-infection counterevidence."
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


def _effective_policy(
    policy_path: Path, *, rpm_threshold: float | None,
    qc_mode: str | None,
) -> dict[str, Any]:
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    policy = copy.deepcopy(policy)
    route = policy.setdefault("analytically_compelling_picked", {})
    if rpm_threshold is not None:
        route["per_test_rpm_threshold"] = float(rpm_threshold)
    if qc_mode == "full_only":
        route["accepted_qc_statuses"] = ["evaluable"]
    elif qc_mode == "partial_or_full":
        route["accepted_qc_statuses"] = ["evaluable", "partial_evaluable"]
    policy["effective_sensitivity_overrides"] = {
        "rpm_threshold": route.get("per_test_rpm_threshold"),
        "qc_mode": qc_mode or "policy_default",
    }
    return policy


def run(
    input_root: Path, output_dir: Path, *, policy_path: Path = DEFAULT_POLICY,
    patients: set[int] | None = None, rpm_threshold: float | None = None,
    qc_mode: str | None = None,
) -> dict[str, Any]:
    policy = _effective_policy(
        policy_path, rpm_threshold=rpm_threshold, qc_mode=qc_mode
    )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    effective_path = output_dir.with_name(
        f"{output_dir.name}_effective_policy.json"
    )
    if effective_path.exists():
        raise FileExistsError(f"Refusing to overwrite {effective_path}")
    effective_text = json.dumps(policy, ensure_ascii=False, indent=2) + "\n"
    effective_path.write_text(effective_text, encoding="utf-8")

    original = base.promotion_decision
    try:
        base.promotion_decision = promotion_decision
        summary = base.run(
            input_root, output_dir, policy_path=effective_path, patients=patients
        )
    finally:
        base.promotion_decision = original
    stored_effective_path = output_dir / "effective_policy.json"
    stored_effective_path.write_text(effective_text, encoding="utf-8")
    effective_path.unlink()
    summary["policy_file"] = str(stored_effective_path.resolve())
    summary["effective_policy_sha256"] = hashlib.sha256(
        stored_effective_path.read_bytes()
    ).hexdigest()
    (output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--patients", type=int, nargs="*")
    parser.add_argument("--rpm-threshold", type=float)
    parser.add_argument(
        "--qc-mode", choices=("partial_or_full", "full_only")
    )
    args = parser.parse_args()
    selected = set(args.patients) if args.patients else None
    print(json.dumps(
        run(
            args.input_root,
            args.output_dir,
            policy_path=args.policy,
            patients=selected,
            rpm_threshold=args.rpm_threshold,
            qc_mode=args.qc_mode,
        ),
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()

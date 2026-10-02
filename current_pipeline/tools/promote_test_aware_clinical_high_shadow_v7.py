"""Add a narrow cross-site convergence route on top of promotion v6.

V7 keeps every v4-v6 route unchanged.  Its only additional route is for an
exact, colonizer-prone organism that would otherwise remain Context/High but
has convergent evidence from two independent clinical axes:

* strong, repeated lower-respiratory mNGS evidence across DNA and RNA; and
* an exact, positive, event-aligned sterile-site culture observation.

The route is deliberately answer blind.  A provisional model contaminant
label is retained as a caution, while explicit clinical colonization or
current-not-infection evidence remains a hard blocker.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

from tools import promote_test_aware_clinical_high_shadow as base
from tools import promote_test_aware_clinical_high_shadow_v6 as v6
from tools.pathogen_normalization import canonical_key


DEFAULT_POLICY = Path(
    "rules/test_aware_clinical_promotion_v7_cross_site_convergence.json"
)


def _positive_assertion(row: dict[str, Any]) -> bool:
    text = f" {base._recursive_text(row)} "
    if any(term in text for term in (
        " not detected", " no growth", " negative", " not isolated"
    )):
        return False
    status = str(row.get("detection_status") or "").casefold()
    quantitation = str(row.get("quantitation_status") or "").casefold()
    return bool(
        status in {"detected", "positive"}
        or quantitation.startswith("positive")
        or any(term in text for term in (
            " isolated", " detected", " positive", " grew", " growth"
        ))
    )


def _sterile_site(row: dict[str, Any], route: dict[str, Any]) -> bool:
    allowed_categories = {
        base._raw_key(value)
        for value in route.get("sterile_specimen_categories") or []
    }
    category = base._raw_key(row.get("specimen_category"))
    if category and category in allowed_categories:
        return True
    text = " ".join(
        str(row.get(key) or "").casefold()
        for key in ("specimen_type", "specimen_category", "sample", "site")
    )
    return any(
        str(term).casefold() in text
        for term in route.get("sterile_specimen_terms") or []
    )


def _pure_observation_ids(
    detail: dict[str, Any], route: dict[str, Any]
) -> set[str]:
    accepted = {
        str(value).casefold()
        for value in route.get("accepted_growth_purity") or []
    }
    output: set[str] = set()
    module_evidence = detail.get("module_evidence") or {}
    if not isinstance(module_evidence, dict):
        return output
    for rows in module_evidence.values():
        for row in rows or []:
            if not isinstance(row, dict):
                continue
            if (
                str(row.get("growth_purity") or "").casefold() in accepted
                and _sterile_site(row, route)
            ):
                output.update(
                    str(value)
                    for value in row.get("evidence_observation_ids") or []
                    if value not in (None, "")
                )
    return output


def cross_site_convergence_profile(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any]:
    route = policy.get("colonizer_prone_cross_site_convergence") or {}
    hospital = candidate.get("hospital_profile") or {}
    detail = hospital.get("hospital_evidence_detail") or {}
    rows = base._module_evidence_rows(detail)
    event_times = base.candidate_event_times(candidate)
    window_hours = float(
        route.get(
            "event_window_hours",
            policy.get("direct_evidence_event_window_hours", 48),
        )
    )
    accepted_purity = {
        str(value).casefold()
        for value in route.get("accepted_growth_purity") or []
    }
    pure_ids = _pure_observation_ids(detail, route)
    target_key = canonical_key(candidate.get("organism_name"))
    evaluated: list[dict[str, Any]] = []

    for row in rows:
        row_name = row.get("organism_name") or row.get("target")
        identity_match = bool(
            row_name and canonical_key(row_name) == target_key
        )
        observed_at = base._parse_datetime(row.get("collected_time"))
        if observed_at is None or not event_times:
            timing = "unknown"
            nearest_delta_hours = None
        else:
            deltas = [
                abs((observed_at - event_time).total_seconds()) / 3600
                for event_time in event_times
            ]
            nearest_delta_hours = min(deltas)
            timing = (
                "within_event_window"
                if nearest_delta_hours <= window_hours
                else "outside_event_window"
            )
        observation_id = str(row.get("observation_id") or "")
        purity = str(row.get("growth_purity") or "").casefold()
        isolated_pure = bool(
            purity in accepted_purity
            or (observation_id and observation_id in pure_ids)
        )
        evaluated.append({
            "observation_id": observation_id or None,
            "organism_name": row_name,
            "identity_match": identity_match,
            "positive": _positive_assertion(row),
            "sterile_site": _sterile_site(row, route),
            "isolated_pure": isolated_pure,
            "timing": timing,
            "nearest_event_delta_hours": nearest_delta_hours,
            "specimen_type": row.get("specimen_type"),
            "specimen_category": row.get("specimen_category"),
            "collected_time": row.get("collected_time"),
            "reported_time": row.get("reported_time"),
            "raw_result": row.get("raw_result"),
        })

    require_pure = bool(route.get("require_isolated_pure", True))
    qualifying = [
        row for row in evaluated
        if row["identity_match"]
        and row["positive"]
        and row["sterile_site"]
        and row["timing"] == "within_event_window"
        and (row["isolated_pure"] or not require_pure)
    ]
    return {
        "profile_version": "colonizer_prone_cross_site_convergence_v1",
        "event_window_hours": window_hours,
        "evaluated_row_count": len(evaluated),
        "qualifying_observation_count": len(qualifying),
        "exact_event_aligned_sterile_positive": bool(qualifying),
        "require_isolated_pure": require_pure,
        "provisional_contaminant_label_preserved_as_caution": bool(
            (candidate.get("source_candidate") or {})
            .get("colonization_interpretation", {})
            .get("has_provisional_contaminant_label")
        ),
        "rows": evaluated,
        "missing_timing_is_negative": False,
    }


def _cross_site_burden_profile(
    candidate: dict[str, Any], route: dict[str, Any]
) -> dict[str, Any]:
    analytical = candidate.get("analytical_profile") or {}
    signals = analytical.get("per_test_signals") or []
    threshold = float(route.get("per_test_rpm_threshold", 5.0))
    accepted_qc = set(route.get("accepted_qc_statuses") or [])
    normalized = []
    qualifying = []
    for signal in signals:
        rpm = signal.get("rpm_total")
        if not signal.get("normalization_available") or rpm is None:
            continue
        row = {
            "seq_id": signal.get("seq_id"),
            "nucleic_type": signal.get("nucleic_type"),
            "qc_status": signal.get("qc_status"),
            "rpm_total": float(rpm),
        }
        normalized.append(row)
        if (
            (not accepted_qc or str(signal.get("qc_status") or "") in accepted_qc)
            and float(rpm) >= threshold
        ):
            qualifying.append(row)
    return {
        "per_test_rpm_threshold": threshold,
        "accepted_qc_statuses": sorted(accepted_qc),
        "normalization_available_test_count": len(normalized),
        "tests_meeting_rpm_and_qc_threshold": len(qualifying),
        "normalized_tests": normalized,
        "qualifying_tests": qualifying,
        "reads_summed_across_tests": False,
    }


def promotion_decision(
    candidate: dict[str, Any], policy: dict[str, Any],
    patient_context: dict[str, Any] | None = None,
) -> tuple[str, dict[str, Any], list[str], list[str]]:
    decision, gate, rules, reasons = v6.promotion_decision(
        candidate, policy, patient_context
    )
    route = policy.get("colonizer_prone_cross_site_convergence") or {}
    current = str(candidate.get("clinical_decision") or "")
    allowed_incoming = set(
        route.get("allowed_incoming_decisions")
        or ["review_context_needed", "review_high_priority"]
    )
    if (
        not route.get("enabled", False)
        or current not in allowed_incoming
        or decision != current
        or str(candidate.get("candidate_source") or "") != "mngs_positive_reads"
    ):
        return decision, gate, rules, reasons

    axes = base.evidence_axes(candidate, policy)
    taxonomy = candidate.get("taxonomy_profile") or {}
    analytical = candidate.get("analytical_profile") or {}
    colonization = v6._colonization_profile(candidate)
    burden = _cross_site_burden_profile(candidate, route)
    convergence = cross_site_convergence_profile(candidate, policy)
    rank = analytical.get("best_rank_in_retained_universe")
    explicit_counterevidence = bool(
        colonization.get("has_explicit_colonization")
        or colonization.get("has_explicit_current_not_infection")
    )
    axes["cross_site_convergence_profile"] = convergence
    axes["cross_site_burden_profile"] = burden
    axes["explicit_colonization_or_current_not_infection"] = (
        explicit_counterevidence
    )

    ready = bool(
        str(taxonomy.get("mapping_status") or "")
            in set(route.get("allowed_mapping_statuses") or [])
        and str(taxonomy.get("primary_rule_family") or "")
            in set(route.get("allowed_families") or [])
        and axes.get("specimen_context")
            in set(route.get("allowed_mngs_specimen_contexts") or [])
        and int(analytical.get("selected_positive_test_count") or 0)
            >= int(route.get("minimum_positive_tests", 3))
        and bool(analytical.get("reproducibility_axis"))
        and bool(analytical.get("cross_molecule_selected"))
        and rank is not None
        and int(rank) <= int(route.get("maximum_rank", 3))
        and burden["normalization_available_test_count"]
            >= int(route.get("minimum_normalized_tests", 1))
        and burden["tests_meeting_rpm_and_qc_threshold"]
            >= int(route.get("minimum_tests_at_or_above_rpm_threshold", 1))
        and convergence["exact_event_aligned_sterile_positive"]
        and not explicit_counterevidence
    )
    if not ready:
        return decision, gate, rules, reasons

    caution = convergence[
        "provisional_contaminant_label_preserved_as_caution"
    ]
    route_rules = [
        "PROMO-CS1-CONTEXT-OR-HIGH-CROSS-SITE-REVIEW",
        "PROMO-CS2-EXACT-LOWER-RESPIRATORY-THREE-TEST-CROSS-MOLECULE",
        "PROMO-CS3-EVENT-ALIGNED-EXACT-STERILE-CULTURE",
        "PROMO-CS4-NO-EXPLICIT-COLONIZATION-COUNTEREVIDENCE",
    ]
    if caution:
        route_rules.append(
            "PROMO-CS5-PROVISIONAL-CONTAMINANT-LABEL-RETAINED-AS-CAUTION"
        )
    return (
        "picked_shadow",
        {
            "evaluated": True,
            "outcome": "promoted_to_picked_shadow",
            "route": "colonizer_prone_cross_site_convergence",
            "axes": axes,
            "blockers": [],
            "cautions": ([
                "A provisional contaminant label remains visible, but it is "
                "not treated as explicit clinical colonization because exact "
                "lower-respiratory mNGS and event-aligned sterile-site culture "
                "converge."
            ] if caution else []),
        },
        rules + route_rules,
        [
            "Exact top-three lower-respiratory mNGS evidence is reproduced "
            "across at least three tests and DNA/RNA, includes a qualifying "
            "per-test normalized burden, and converges with an exact positive "
            "sterile-site culture collected within the event window."
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
    policy = v6._effective_policy(
        policy_path, rpm_threshold=rpm_threshold, qc_mode=qc_mode
    )
    route = policy.setdefault("colonizer_prone_cross_site_convergence", {})
    if rpm_threshold is not None:
        route["per_test_rpm_threshold"] = float(rpm_threshold)
    if qc_mode == "full_only":
        route["accepted_qc_statuses"] = ["evaluable"]
    elif qc_mode == "partial_or_full":
        route["accepted_qc_statuses"] = ["evaluable", "partial_evaluable"]
    policy["effective_sensitivity_overrides"]["cross_site_convergence"] = {
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

"""Consolidate an event-aligned respiratory-virus group assay into an exact mNGS identity.

V8 preserves every v4-v7 route.  Its only additional behavior is a narrow,
answer-blind representation rule: an exact respiratory-virus mNGS candidate
may represent a broad FilmArray group target when the names have an explicitly
approved group/member relationship and the positive lower-respiratory assay is
within the configured event window.  The assay is retained as group-level
support and is never described as exact subtype confirmation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from tools import promote_test_aware_clinical_high_shadow as base
from tools import promote_test_aware_clinical_high_shadow_v7 as v7
from tools.pathogen_normalization import (
    approved_group_member_match,
    canonical_key,
)


DEFAULT_POLICY = Path(
    "rules/test_aware_clinical_promotion_v8_respiratory_group_convergence.json"
)
_BASE_BUILD_PATIENT_CONTEXT = base.build_patient_promotion_context


def _positive_assertion(row: dict[str, Any]) -> bool:
    text = f" {base._recursive_text(row)} "
    if any(term in text for term in (
        " not detected", " negative", " invalid", " equivocal"
    )):
        return False
    status = str(row.get("detection_status") or "").casefold()
    return bool(
        status in {"detected", "positive"}
        or any(term in text for term in (" detected", " positive", " reactive"))
    )


def _lower_respiratory(row: dict[str, Any]) -> bool:
    text = " ".join(
        str(row.get(key) or "").casefold()
        for key in ("specimen_type", "sample", "site", "filmarray_sample_type")
    )
    return any(term in text for term in (
        "bal", "bronchoalveolar", "endotracheal", "tracheal", "sputum",
        "lower respiratory",
    ))


def _hospital_detail(candidate: dict[str, Any]) -> dict[str, Any]:
    if candidate.get("candidate_source") == "hospital_only":
        return candidate.get("hospital_evidence_detail") or {}
    return (candidate.get("hospital_profile") or {}).get(
        "hospital_evidence_detail"
    ) or {}


def respiratory_group_convergence_profile(
    exact_candidate: dict[str, Any],
    group_candidate: dict[str, Any],
    policy: dict[str, Any],
) -> dict[str, Any]:
    route = policy.get("respiratory_virus_group_convergence") or {}
    detail = _hospital_detail(group_candidate)
    rows = base._module_evidence_rows(detail)
    event_times = base.candidate_event_times(exact_candidate)
    window_hours = float(
        route.get(
            "event_window_hours",
            policy.get("direct_evidence_event_window_hours", 48),
        )
    )
    evaluated: list[dict[str, Any]] = []
    for row in rows:
        observed_at = base._parse_datetime(row.get("collected_time"))
        time_basis = "collected_time"
        if observed_at is None:
            observed_at = base._parse_datetime(row.get("reported_time"))
            time_basis = "reported_time" if observed_at is not None else "unavailable"
        if observed_at is None or not event_times:
            timing = "unknown"
            nearest_delta_hours = None
        else:
            deltas = [
                abs((observed_at - event).total_seconds()) / 3600
                for event in event_times
            ]
            nearest_delta_hours = min(deltas)
            timing = (
                "within_event_window"
                if nearest_delta_hours <= window_hours
                else "outside_event_window"
            )
        evaluated.append({
            "observation_id": row.get("observation_id"),
            "target": row.get("target") or row.get("organism_name"),
            "positive": _positive_assertion(row),
            "lower_respiratory": _lower_respiratory(row),
            "timing": timing,
            "time_basis": time_basis,
            "nearest_event_delta_hours": nearest_delta_hours,
            "collected_time": row.get("collected_time"),
            "reported_time": row.get("reported_time"),
            "specimen_type": row.get("specimen_type"),
            "raw_result": row.get("raw_result"),
        })

    exact_name = exact_candidate.get("organism_name")
    group_name = group_candidate.get("organism_name")
    analytical = exact_candidate.get("analytical_profile") or {}
    taxonomy = exact_candidate.get("taxonomy_profile") or {}
    rank = base.as_int(analytical.get("best_rank_in_retained_universe"))
    selected_count = base.as_int(analytical.get("selected_positive_test_count")) or 0
    group_level = base.direct_hospital_level(group_candidate)
    verified_exact_keys = {
        canonical_key(value)
        for value in route.get("verified_exact_member_names") or []
    }
    identity_verified = bool(
        str(taxonomy.get("mapping_status") or "") == "exact_species"
        or canonical_key(exact_name) in verified_exact_keys
    )
    qualifying_rows = [
        row for row in evaluated
        if row["positive"]
        and row["lower_respiratory"]
        and row["timing"] == "within_event_window"
    ]
    qualifies = bool(
        approved_group_member_match(exact_name, group_name)
        and identity_verified
        and str((taxonomy.get("primary_rule_family") or ""))
            in set(route.get("allowed_families") or ["other_respiratory_virus"])
        and exact_candidate.get("candidate_source") == "mngs_positive_reads"
        and group_candidate.get("candidate_source") == "hospital_only"
        and base.specimen_context(exact_candidate)
            in set(route.get("allowed_mngs_specimen_contexts") or ["lower_respiratory"])
        and rank is not None
        and rank <= int(route.get("maximum_rank", 3))
        and selected_count >= int(route.get("minimum_positive_tests", 1))
        and group_level is not None
        and group_level <= int(route.get("maximum_group_assay_level", 2))
        and qualifying_rows
    )
    return {
        "profile_version": "respiratory_virus_group_convergence_v1",
        "qualifies": qualifies,
        "exact_organism_name": exact_name,
        "exact_organism_key": canonical_key(exact_name),
        "group_assay_name": group_name,
        "group_assay_key": canonical_key(group_name),
        "approved_group_member_match": approved_group_member_match(
            exact_name, group_name
        ),
        "identity_verified": identity_verified,
        "identity_scope": "group_level_support_not_exact_species_confirmation",
        "best_rank": rank,
        "selected_positive_test_count": selected_count,
        "group_assay_level": group_level,
        "event_window_hours": window_hours,
        "qualifying_observation_count": len(qualifying_rows),
        "rows": evaluated,
        "missing_timing_is_negative": False,
    }


def build_patient_promotion_context(
    candidates: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    context = _BASE_BUILD_PATIENT_CONTEXT(candidates, policy)
    route = policy.get("respiratory_virus_group_convergence") or {}
    links: list[dict[str, Any]] = []
    if route.get("enabled", False):
        exact_candidates = [
            item for item in candidates
            if item.get("candidate_source") == "mngs_positive_reads"
        ]
        group_candidates = [
            item for item in candidates
            if item.get("candidate_source") == "hospital_only"
        ]
        for exact in exact_candidates:
            for group in group_candidates:
                if not approved_group_member_match(
                    exact.get("organism_name"), group.get("organism_name")
                ):
                    continue
                profile = respiratory_group_convergence_profile(
                    exact, group, policy
                )
                if profile["qualifies"]:
                    links.append(profile)
    context["respiratory_virus_group_convergence_links"] = links
    return context


def _candidate_link(
    candidate: dict[str, Any],
    patient_context: dict[str, Any] | None,
    *,
    role: str,
) -> dict[str, Any] | None:
    key = canonical_key(candidate.get("organism_name"))
    field = "exact_organism_key" if role == "exact" else "group_assay_key"
    for link in (patient_context or {}).get(
        "respiratory_virus_group_convergence_links"
    ) or []:
        if link.get(field) == key:
            return link
    return None


def promotion_decision(
    candidate: dict[str, Any],
    policy: dict[str, Any],
    patient_context: dict[str, Any] | None = None,
) -> tuple[str, dict[str, Any], list[str], list[str]]:
    decision, gate, rules, reasons = v7.promotion_decision(
        candidate, policy, patient_context
    )
    route = policy.get("respiratory_virus_group_convergence") or {}
    if not route.get("enabled", False):
        return decision, gate, rules, reasons

    exact_link = _candidate_link(candidate, patient_context, role="exact")
    if exact_link and decision != "picked_shadow":
        axes = base.evidence_axes(candidate, policy)
        axes["respiratory_virus_group_convergence"] = exact_link
        return (
            "picked_shadow",
            {
                "evaluated": True,
                "outcome": "promoted_to_picked_shadow",
                "route": "respiratory_virus_exact_mngs_plus_group_assay",
                "axes": axes,
                "blockers": [],
                "cautions": [
                    "The FilmArray target supplies group-level confirmation only; exact species identity comes from mNGS and is not claimed as exact panel subtyping."
                ],
            },
            rules + [
                "PROMO-RV1-APPROVED-GROUP-MEMBER-RELATIONSHIP",
                "PROMO-RV2-EXACT-TOP3-MNGS-MEMBER",
                "PROMO-RV3-EVENT-ALIGNED-LOWER-RESPIRATORY-GROUP-ASSAY",
                "PROMO-RV4-GROUP-ASSAY-NOT-EXACT-SUBTYPE-CONFIRMATION",
            ],
            reasons + [
                "An exact top-three respiratory-virus mNGS member converges with an approved, positive lower-respiratory group assay within the event window."
            ],
        )

    group_link = _candidate_link(candidate, patient_context, role="group")
    if group_link:
        axes = base.evidence_axes(candidate, policy)
        axes["respiratory_virus_group_convergence"] = group_link
        return (
            "review_context_needed",
            {
                "evaluated": True,
                "outcome": "represented_by_exact_mngs_member",
                "route": "broad_group_assay_consolidated_under_exact_mngs_member",
                "axes": axes,
                "blockers": [],
                "cautions": [
                    "The broad assay observation remains in audit evidence but is not reported as a second pathogen identity."
                ],
            },
            rules + ["PROMO-RV5-SUPPRESS-DUPLICATE-BROAD-REPORT-IDENTITY"],
            reasons + [
                "The broad panel target is represented by its event-aligned exact mNGS member while retaining the original assay wording in audit evidence."
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
    effective_path = output_dir.with_name(
        f"{output_dir.name}_effective_policy.json"
    )
    if effective_path.exists():
        raise FileExistsError(f"Refusing to overwrite {effective_path}")
    effective_text = json.dumps(policy, ensure_ascii=False, indent=2) + "\n"
    effective_path.write_text(effective_text, encoding="utf-8")

    original_decision = base.promotion_decision
    original_context = base.build_patient_promotion_context
    try:
        base.promotion_decision = promotion_decision
        base.build_patient_promotion_context = build_patient_promotion_context
        summary = base.run(
            input_root,
            output_dir,
            policy_path=effective_path,
            patients=patients,
        )
    finally:
        base.promotion_decision = original_decision
        base.build_patient_promotion_context = original_context
    stored_effective_path = output_dir / "effective_policy.json"
    stored_effective_path.write_text(effective_text, encoding="utf-8")
    effective_path.unlink()
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--patients", type=int, nargs="*")
    args = parser.parse_args()
    selected = set(args.patients) if args.patients else None
    print(json.dumps(run(
        args.input_root,
        args.output_dir,
        policy_path=args.policy,
        patients=selected,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

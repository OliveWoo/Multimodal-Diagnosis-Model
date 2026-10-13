"""Apply a bounded HAP/airway-device Context-to-Priority shadow overlay."""

from __future__ import annotations

import argparse
import copy
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

from tools.pathogen_normalization import canonical_key
from tools.run_kh_ablation_d_history_only import read_json, sha256_file, write_json


DEFAULT_POLICY = Path("rules/healthcare_device_opportunist_history_route_v2.json")
SCHEMA_VERSION = "healthcare_device_opportunist_history_route.v1"


def _taxonomy_family(candidate: dict[str, Any]) -> str:
    profile = candidate.get("taxonomy_profile") or {}
    return str(profile.get("primary_rule_family") or candidate.get("taxonomy_family") or "")


def _taxonomy_traits(candidate: dict[str, Any]) -> set[str]:
    profile = candidate.get("taxonomy_profile") or {}
    return {str(value) for value in profile.get("clinical_traits") or []}


def _specimen_context(candidate: dict[str, Any]) -> str:
    value = str(candidate.get("specimen_context") or "")
    if value:
        return value
    specimen_class = str(candidate.get("specimen_class") or "").casefold()
    return "lower_respiratory" if "lower_respiratory" in specimen_class else "unknown"


def _analytical_axes(candidate: dict[str, Any]) -> tuple[int | None, int]:
    profile = candidate.get("analytical_profile") or {}
    rank = profile.get("best_rank_in_retained_universe")
    try:
        rank = int(rank) if rank is not None else None
    except (TypeError, ValueError):
        rank = None
    try:
        positive_tests = int(profile.get("selected_positive_test_count") or 0)
    except (TypeError, ValueError):
        positive_tests = 0
    return rank, positive_tests


def _direct_colonization_block(candidate: dict[str, Any]) -> bool:
    review = candidate.get("colonization_interpretation") or {}
    return bool(
        review.get("has_direct_nonsterile_match")
        and (
            review.get("has_explicit_colonization")
            or review.get("has_explicit_current_not_infection")
        )
    )


def _has_exact_hospital_support(candidate: dict[str, Any]) -> bool:
    if candidate.get("evidence_source") == "hospital_only":
        detail = candidate.get("hospital_evidence_detail") or {}
        return bool(detail.get("best_hospital_level") or detail.get("source_observations"))
    profile = candidate.get("hospital_profile") or {}
    return bool(
        profile.get("exact_or_alias_evidence")
        or profile.get("best_direct_level")
        or (profile.get("hospital_evidence_detail") or {}).get("source_observations")
    )


def candidate_eligibility(
    candidate: dict[str, Any], patient_payload: dict[str, Any], policy: dict[str, Any]
) -> tuple[bool, list[str], list[str], list[dict[str, Any]]]:
    transition = policy["eligible_transition"]
    rank, positive_tests = _analytical_axes(candidate)
    traits = _taxonomy_traits(candidate)
    hap_profile = patient_payload.get("hap_mdr_history_profile") or {}
    colonization_blocked = _direct_colonization_block(candidate)
    checks: list[tuple[bool, str]] = [
        (
            str(candidate.get("history_adjusted_route")) == str(transition["from_route"]),
            "eligible_from_route",
        ),
        (
            _taxonomy_family(candidate) in set(policy["eligible_taxonomy_families"]),
            "eligible_taxonomy_family",
        ),
        (
            set(policy["required_all_taxonomy_traits"]).issubset(traits),
            "required_taxonomy_traits",
        ),
        (
            str(hap_profile.get("support_level"))
            in set(policy["required_history_support_levels"]),
            "required_hap_or_airway_device_history",
        ),
        (
            _specimen_context(candidate) in set(policy["required_specimen_contexts"]),
            "compatible_specimen_context",
        ),
        (
            rank is not None
            and rank <= int(policy["max_best_rank_in_retained_universe"]),
            "analytical_rank_threshold",
        ),
        (
            positive_tests >= int(policy["min_selected_positive_test_count"]),
            "selected_positive_test_threshold",
        ),
    ]
    if policy.get("require_positive_mngs_signal"):
        checks.append((positive_tests > 0, "positive_mngs_signal"))
    if policy.get("require_upstream_route_unchanged"):
        checks.append((not bool(candidate.get("history_route_changed")), "upstream_route_unchanged"))
    if policy.get("block_when_direct_colonization_or_not_infection"):
        checks.append((not colonization_blocked, "no_direct_colonization_or_not_infection"))

    passed = [reason for ok, reason in checks if ok]
    blockers = [reason for ok, reason in checks if not ok]
    categories = set(policy.get("consumable_history_categories") or [])
    evidence = [
        row
        for row in hap_profile.get("matched_evidence") or []
        if set(row.get("categories") or []) & categories
    ]
    if not evidence:
        blockers.append("no_consumable_hap_or_airway_device_evidence")
    return not blockers, passed, blockers, evidence


def _priority_sort_key(candidate: dict[str, Any]) -> tuple[Any, ...]:
    rank, positive_tests = _analytical_axes(candidate)
    analytical = candidate.get("analytical_profile") or {}
    reads = candidate.get("reads")
    try:
        reads_value = int(reads or 0)
    except (TypeError, ValueError):
        reads_value = 0
    return (
        -int(_has_exact_hospital_support(candidate)),
        rank if rank is not None else 10**9,
        -positive_tests,
        -int(bool(analytical.get("reproducibility_axis"))),
        -int(bool(analytical.get("cross_molecule_selected"))),
        -reads_value,
        canonical_key(candidate.get("organism_name")),
    )


def apply_patient(
    payload: dict[str, Any], policy: dict[str, Any]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    output = copy.deepcopy(payload)
    candidates = output.get("candidates") or []
    evaluated: list[tuple[int, bool, list[str], list[str], list[dict[str, Any]]]] = []
    for index, candidate in enumerate(candidates):
        eligible, passed, blockers, evidence = candidate_eligibility(candidate, output, policy)
        evaluated.append((index, eligible, passed, blockers, evidence))

    eligible_indices = [index for index, eligible, *_ in evaluated if eligible]
    eligible_indices.sort(key=lambda index: _priority_sort_key(candidates[index]))
    cap_exempt = {
        index
        for index in eligible_indices
        if policy.get("cap_exempt_when_exact_hospital_support")
        and _has_exact_hospital_support(candidates[index])
    }
    capped_candidates = [index for index in eligible_indices if index not in cap_exempt]
    retained = cap_exempt | set(
        capped_candidates[: int(policy["max_promoted_representatives_per_patient"])]
    )
    audit_rows: list[dict[str, Any]] = []
    transition = policy["eligible_transition"]

    for index, eligible, passed, blockers, evidence in evaluated:
        candidate = candidates[index]
        promoted = eligible and index in retained
        if eligible and not promoted:
            blockers = [*blockers, "patient_representative_cap"]
        gate = {
            "policy_id": policy["policy_id"],
            "eligible_before_patient_cap": eligible,
            "promoted_to_priority": promoted,
            "patient_cap_exempt_exact_hospital_support": index in cap_exempt,
            "passed_checks": passed,
            "blockers": blockers,
            "consumed_history_evidence_ids": (
                [str(row.get("evidence_id")) for row in evidence if row.get("evidence_id")]
                if promoted
                else []
            ),
        }
        candidate["healthcare_device_opportunist_route_gate"] = gate
        before_route = str(candidate.get("history_adjusted_route") or "")
        before_decision = str(candidate.get("history_routed_decision") or "")
        if promoted:
            evidence_ids = gate["consumed_history_evidence_ids"]
            candidate["history_adjusted_route"] = transition["to_route"]
            candidate["history_route_changed"] = True
            candidate["history_policy_family"] = "healthcare_device_opportunist"
            candidate["history_rule_ids"] = list(dict.fromkeys([
                *(candidate.get("history_rule_ids") or []), policy["rule_id"]
            ]))
            candidate["history_evidence_ids"] = list(dict.fromkeys([
                *(candidate.get("history_evidence_ids") or []), *evidence_ids
            ]))
            candidate["consumed_evidence_ids"] = list(dict.fromkeys([
                *(candidate.get("consumed_evidence_ids") or []), *evidence_ids
            ]))
            candidate["history_evidence_consumption"] = [
                *(candidate.get("history_evidence_consumption") or []),
                {
                    "policy_id": policy["policy_id"],
                    "rule_id": policy["rule_id"],
                    "evidence_ids": evidence_ids,
                    "reuse_as_independent_support_allowed": False,
                },
            ]
            candidate["upstream_history_adjustment_blockers"] = list(
                candidate.get("history_adjustment_blockers") or []
            )
            candidate["history_adjustment_blockers"] = []
            candidate["history_routed_decision"] = transition["final_decision"]
            candidate["downstream_route"] = {
                "scorer_review_tier": "Priority",
                "projected_final_tier": "Priority",
                "picked_allowed_from_history_alone": False,
                "picked_promotion_blockers": [
                    "HAP or airway-device setting does not independently identify a pathogen"
                ],
                "note": "Priority review only; Picked requires independent exact-organism evidence.",
            }
        audit_rows.append({
            "patient_id": output.get("patient_id"),
            "organism_name": candidate.get("organism_name"),
            "taxonomy_family": _taxonomy_family(candidate),
            "history_support_level": (output.get("hap_mdr_history_profile") or {}).get("support_level"),
            "before_route": before_route,
            "after_route": candidate.get("history_adjusted_route"),
            "before_decision": before_decision,
            "after_decision": candidate.get("history_routed_decision"),
            "eligible_before_patient_cap": eligible,
            "promoted_to_priority": promoted,
            "patient_cap_exempt_exact_hospital_support": index in cap_exempt,
            "passed_checks": "|".join(passed),
            "blockers": "|".join(blockers),
        })

    output["schema_version"] = SCHEMA_VERSION
    output["healthcare_device_opportunist_policy_id"] = policy["policy_id"]
    output["healthcare_device_opportunist_route_changes"] = [
        row["organism_name"] for row in audit_rows if row["promoted_to_priority"]
    ]
    return output, audit_rows


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["patient_id", "organism_name"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(
    history_root: Path, output_root: Path, policy_path: Path = DEFAULT_POLICY
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    policy = read_json(policy_path)
    rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    patient_output = output_root / "patient_outputs"
    for source in sorted((history_root / "patient_outputs").glob("*.json")):
        payload = read_json(source)
        adjusted, audit_rows = apply_patient(payload, policy)
        adjusted["source_health_history_output"] = {
            "path": str(source.resolve()),
            "sha256": sha256_file(source),
        }
        target = patient_output / source.name.replace(
            "_history_adjusted_route_shadow.json",
            "_healthcare_device_route_shadow.json",
        )
        write_json(target, adjusted)
        rows.extend(audit_rows)
        counts["patients"] += 1
        counts["candidates"] += len(audit_rows)
        counts["promoted_to_priority"] += sum(
            bool(row["promoted_to_priority"]) for row in audit_rows
        )
    _write_csv(output_root / "route_decisions.csv", rows)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "policy_id": policy["policy_id"],
        "answer_blind": True,
        "source_history_root": str(history_root.resolve()),
        "source_history_summary_sha256": sha256_file(history_root / "summary.json"),
        "policy_path": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "counts": dict(sorted(counts.items())),
        "constraints": policy.get("constraints") or [],
    }
    write_json(output_root / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("history_root", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    print(json.dumps(run(args.history_root, args.output_root, args.policy), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

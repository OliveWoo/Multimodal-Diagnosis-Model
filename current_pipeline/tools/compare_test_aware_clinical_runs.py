from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

from . import pathogen_normalization


TIER_LABELS = {
    "picked_shadow": "Picked",
    "review_high_priority": "High",
    "review_context_needed": "Context",
}


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def clinical_summary(promotion_summary: dict[str, Any]) -> dict[str, Any]:
    root = promotion_summary.get("source_clinical_shadow_root")
    if not root:
        return {}
    path = Path(str(root)) / "summary.json"
    return load_json(path) if path.is_file() else {}


def candidate_key(patient_id: str, candidate: dict[str, Any]) -> tuple[str, str]:
    organism_key = str(candidate.get("organism_key") or "").strip()
    if not organism_key:
        name = str(candidate.get("organism_name") or "").strip()
        organism_key = pathogen_normalization.canonical_key(name) or name.casefold()
    if not organism_key:
        raise ValueError(f"Candidate without organism identity for patient {patient_id}")
    return patient_id, organism_key


def load_candidates(root: Path) -> dict[tuple[str, str], dict[str, Any]]:
    indexed: dict[tuple[str, str], dict[str, Any]] = {}
    patient_root = root / "patient_outputs"
    for path in sorted(patient_root.glob("*.json")):
        payload = load_json(path)
        patient_id = str(payload.get("patient_id") or "").strip()
        if not patient_id:
            raise ValueError(f"Missing patient_id: {path}")
        candidates = payload.get("all_forwarded_candidates") or []
        if not isinstance(candidates, list):
            raise ValueError(f"Expected all_forwarded_candidates list: {path}")
        for candidate in candidates:
            if not isinstance(candidate, dict):
                continue
            key = candidate_key(patient_id, candidate)
            if key in indexed:
                raise ValueError(f"Duplicate candidate key: patient={key[0]} organism={key[1]}")
            indexed[key] = candidate
    return indexed


def tier(candidate: dict[str, Any] | None) -> str:
    if candidate is None:
        return "Not forwarded"
    decision = str(candidate.get("clinical_decision") or "")
    return TIER_LABELS.get(decision, decision or "Unknown")


def policy_checks(
    control_summary: dict[str, Any],
    candidate_summary: dict[str, Any],
) -> dict[str, Any]:
    control_clinical = clinical_summary(control_summary)
    candidate_clinical = clinical_summary(candidate_summary)
    checks = {
        "promotion_policy_sha256_equal": control_summary.get("policy_sha256")
        == candidate_summary.get("policy_sha256"),
        "clinical_policy_sha256_equal": control_clinical.get("policy_sha256")
        == candidate_clinical.get("policy_sha256"),
        "source_scorer_root_equal": control_clinical.get("source_scorer_root")
        == candidate_clinical.get("source_scorer_root"),
        "timeline_root_equal": control_clinical.get("timeline_root")
        == candidate_clinical.get("timeline_root"),
        "phenotype_packet_root_equal": control_clinical.get("phenotype_packet_root")
        == candidate_clinical.get("phenotype_packet_root"),
    }
    checks["all_equal"] = all(checks.values())
    return checks


def double_count_guard_pass(candidate: dict[str, Any]) -> bool:
    consumed = candidate.get("consumed_evidence_ids") or []
    if not consumed:
        return True
    gate = candidate.get("promotion_gate") or {}
    if not gate.get("evaluated"):
        return True
    axes = gate.get("axes") or {}
    return axes.get("host_support_eligible_for_promotion") is False


def route_effect_state(candidate: dict[str, Any] | None) -> tuple[Any, ...]:
    if candidate is None or not candidate.get("history_route_changed"):
        return (False,)
    return (
        True,
        candidate.get("analytic_route"),
        candidate.get("history_adjusted_route"),
        candidate.get("history_policy_family"),
        tuple(candidate.get("history_rule_ids") or []),
        tuple(candidate.get("consumed_evidence_ids") or []),
    )


def compare_runs(control_root: Path, candidate_root: Path) -> dict[str, Any]:
    control_summary = load_json(control_root / "summary.json")
    candidate_summary = load_json(candidate_root / "summary.json")
    checks = policy_checks(control_summary, candidate_summary)
    control = load_candidates(control_root)
    candidate = load_candidates(candidate_root)

    rows: list[dict[str, Any]] = []
    all_keys = sorted(
        set(control) | set(candidate),
        key=lambda value: (int(value[0]) if value[0].isdigit() else 10**9, value),
    )
    for key in all_keys:
        before = control.get(key)
        after = candidate.get(key)
        after_route_changed = bool(after and after.get("history_route_changed"))
        incremental_route_change = route_effect_state(before) != route_effect_state(after)
        if (
            before is not None
            and after is not None
            and tier(before) == tier(after)
            and not incremental_route_change
        ):
            continue
        chosen = after or before or {}
        gate = (after or {}).get("promotion_gate") or {}
        rows.append(
            {
                "patient_id": key[0],
                "organism_key": key[1],
                "organism_name": chosen.get("organism_name") or key[1],
                "before_final_tier": tier(before),
                "after_final_tier": tier(after),
                "analytic_route": (after or {}).get("analytic_route"),
                "history_adjusted_route": (after or {}).get("history_adjusted_route"),
                "history_route_changed": after_route_changed,
                "incremental_route_effect_changed": incremental_route_change,
                "history_rule_ids": (after or {}).get("history_rule_ids") or [],
                "history_evidence_ids": (after or {}).get("history_evidence_ids") or [],
                "consumed_evidence_ids": (after or {}).get("consumed_evidence_ids") or [],
                "promotion_outcome": gate.get("outcome"),
                "promotion_route": gate.get("route"),
                "promotion_blockers": gate.get("blockers") or [],
                "double_count_guard_pass": double_count_guard_pass(after or {}),
            }
        )

    before_tiers = Counter(tier(item) for item in control.values())
    after_tiers = Counter(tier(item) for item in candidate.values())
    added = set(candidate) - set(control)
    removed = set(control) - set(candidate)
    common_tier_changes = [
        key for key in set(control) & set(candidate) if tier(control[key]) != tier(candidate[key])
    ]
    history_changed = [item for item in candidate.values() if item.get("history_route_changed")]
    incremental_route_changes = [
        key for key in set(control) | set(candidate)
        if route_effect_state(control.get(key)) != route_effect_state(candidate.get(key))
    ]
    double_count_violations = [
        row for row in rows if row["consumed_evidence_ids"] and not row["double_count_guard_pass"]
    ]
    summary = {
        "control_candidate_count": len(control),
        "candidate_candidate_count": len(candidate),
        "candidate_added_count": len(added),
        "candidate_removed_count": len(removed),
        "existing_candidate_tier_change_count": len(common_tier_changes),
        "changed_patient_count": len({row["patient_id"] for row in rows}),
        "history_route_changed_count": len(history_changed),
        "incremental_route_effect_change_count": len(incremental_route_changes),
        "history_changed_ending_picked_count": sum(tier(item) == "Picked" for item in history_changed),
        "picked_gained_count": sum(
            tier(candidate[key]) == "Picked" and tier(control.get(key)) != "Picked"
            for key in candidate
        ),
        "picked_lost_count": sum(
            tier(control[key]) == "Picked" and tier(candidate.get(key)) != "Picked"
            for key in control
        ),
        "double_count_guard_violation_count": len(double_count_violations),
        "tier_counts_before": dict(sorted(before_tiers.items())),
        "tier_counts_after": dict(sorted(after_tiers.items())),
        "tier_count_delta": {
            label: after_tiers[label] - before_tiers[label]
            for label in sorted(set(before_tiers) | set(after_tiers))
        },
    }
    return {
        "schema_version": "test_aware_clinical_run_comparison.v1",
        "answer_blind": True,
        "control_root": str(control_root.resolve()),
        "candidate_root": str(candidate_root.resolve()),
        "policy_and_input_checks": checks,
        "summary": summary,
        "changes": rows,
    }


def list_text(values: list[Any]) -> str:
    return " | ".join(str(value) for value in values)


def write_csv(path: Path, changes: list[dict[str, Any]]) -> None:
    fieldnames = [
        "patient_id",
        "organism_name",
        "before_final_tier",
        "after_final_tier",
        "analytic_route",
        "history_adjusted_route",
        "history_route_changed",
        "incremental_route_effect_changed",
        "history_rule_ids",
        "history_evidence_ids",
        "consumed_evidence_ids",
        "promotion_outcome",
        "promotion_route",
        "promotion_blockers",
        "double_count_guard_pass",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for change in changes:
            row = {field: change.get(field) for field in fieldnames}
            for field in (
                "history_rule_ids",
                "history_evidence_ids",
                "consumed_evidence_ids",
                "promotion_blockers",
            ):
                row[field] = list_text(row[field] or [])
            writer.writerow(row)


def write_markdown(path: Path, report: dict[str, Any]) -> None:
    summary = report["summary"]
    checks = report["policy_and_input_checks"]
    lines = [
        "# Test-aware clinical run comparison",
        "",
        "This comparison is answer-blind and does not read benchmark labels.",
        "",
        "## Comparability checks",
        "",
        f"- Same frozen inputs and policy snapshots: **{checks['all_equal']}**",
        f"- Same promotion policy SHA-256: {checks['promotion_policy_sha256_equal']}",
        f"- Same clinical policy SHA-256: {checks['clinical_policy_sha256_equal']}",
        f"- Same source scorer, timeline, and phenotype packets: "
        f"{checks['source_scorer_root_equal'] and checks['timeline_root_equal'] and checks['phenotype_packet_root_equal']}",
        "",
        "## Summary",
        "",
        f"- Candidates: {summary['control_candidate_count']} -> {summary['candidate_candidate_count']}",
        f"- Added / removed: {summary['candidate_added_count']} / {summary['candidate_removed_count']}",
        f"- Existing candidate tier changes: {summary['existing_candidate_tier_change_count']}",
        f"- Patients changed: {summary['changed_patient_count']}",
        f"- History route changes: {summary['history_route_changed_count']}",
        f"- Incremental route-effect changes versus control: "
        f"{summary['incremental_route_effect_change_count']}",
        f"- History-adjusted candidates ending as Picked: {summary['history_changed_ending_picked_count']}",
        f"- Picked gained / lost: {summary['picked_gained_count']} / {summary['picked_lost_count']}",
        f"- Consumed-history double-count violations: {summary['double_count_guard_violation_count']}",
        "",
        "## Tier counts",
        "",
        "| Tier | Before | After | Delta |",
        "|---|---:|---:|---:|",
    ]
    labels = sorted(
        set(summary["tier_counts_before"]) | set(summary["tier_counts_after"]),
        key=lambda label: {"Picked": 0, "High": 1, "Context": 2}.get(label, 9),
    )
    for label in labels:
        lines.append(
            f"| {label} | {summary['tier_counts_before'].get(label, 0)} | "
            f"{summary['tier_counts_after'].get(label, 0)} | "
            f"{summary['tier_count_delta'].get(label, 0):+d} |"
        )
    lines.extend(
        [
            "",
            "## Candidate changes",
            "",
            "| Patient | Organism | Before | After | Route | Promotion |",
            "|---:|---|---|---|---|---|",
        ]
    )
    for row in report["changes"]:
        route = f"{row['analytic_route']} -> {row['history_adjusted_route']}"
        promotion = row["promotion_outcome"] or "not evaluated"
        lines.append(
            f"| {row['patient_id']} | {row['organism_name']} | "
            f"{row['before_final_tier']} | {row['after_final_tier']} | "
            f"{route} | {promotion} |"
        )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare two test-aware clinical promotion outputs while verifying that "
            "the frozen inputs and policy snapshots are identical."
        )
    )
    parser.add_argument("control_root", type=Path)
    parser.add_argument("candidate_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--allow-policy-mismatch",
        action="store_true",
        help="Write a report even when frozen input or policy checks differ.",
    )
    args = parser.parse_args()

    report = compare_runs(args.control_root, args.candidate_root)
    if not report["policy_and_input_checks"]["all_equal"] and not args.allow_policy_mismatch:
        raise SystemExit(
            "Refusing non-isolated comparison: policy or frozen input snapshots differ. "
            "Use --allow-policy-mismatch only for an explicitly non-causal audit."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "comparison.json"
    csv_path = args.output_dir / "candidate_changes.csv"
    markdown_path = args.output_dir / "comparison.md"
    json_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_csv(csv_path, report["changes"])
    write_markdown(markdown_path, report)
    print(json.dumps(report["summary"], ensure_ascii=False, indent=2))
    print(f"json={json_path}")
    print(f"csv={csv_path}")
    print(f"markdown={markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

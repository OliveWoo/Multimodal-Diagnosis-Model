"""Apply an answer-blind, bounded Priority-to-Picked gate to legacy history output."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any


DEFAULT_POLICY = Path("rules/legacy_history_picked_promotion_v1.json")
SCHEMA_VERSION = "legacy_history_picked_promotion.v1"


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot identify patient from {path}")
    return int(match.group(1))


def as_int(value: Any) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def promotion_reasons(
    candidate: dict[str, Any],
    patient_payload: dict[str, Any],
    policy: dict[str, Any],
) -> tuple[bool, list[str], list[str]]:
    passed: list[str] = []
    blockers: list[str] = []
    transition = policy["eligible_transition"]

    checks = [
        (
            str(candidate.get("analytic_route")) == str(transition["from_route"]),
            "eligible_from_route",
        ),
        (
            str(candidate.get("history_adjusted_route")) == str(transition["to_route"]),
            "eligible_to_route",
        ),
        (
            str(candidate.get("history_policy_family"))
            == str(transition["history_policy_family"]),
            "eligible_history_policy_family",
        ),
        (
            str(candidate.get("taxonomy_family"))
            in set(policy["eligible_taxonomy_families"]),
            "eligible_taxonomy_family",
        ),
        (
            str(candidate.get("specimen_context"))
            in set(policy["required_specimen_contexts"]),
            "compatible_specimen_context",
        ),
    ]
    if policy.get("require_history_route_changed"):
        checks.append((bool(candidate.get("history_route_changed")), "history_route_changed"))

    aspiration = patient_payload.get("aspiration_history_profile") or {}
    checks.append(
        (
            str(aspiration.get("support_level"))
            == str(policy["required_history_support_level"]),
            "required_history_support_level",
        )
    )
    if policy.get("require_explicit_aspiration_event"):
        checks.append(
            (bool(aspiration.get("has_explicit_aspiration_event")), "explicit_aspiration_event")
        )

    analytical = candidate.get("analytical_profile") or {}
    rank = as_int(analytical.get("best_rank_in_retained_universe"))
    checks.extend(
        [
            (
                rank is not None
                and rank <= int(policy["max_best_rank_in_retained_universe"]),
                "analytical_rank_threshold",
            ),
            (
                int(analytical.get("selected_positive_test_count") or 0)
                >= int(policy["min_selected_positive_test_count"]),
                "selected_positive_test_threshold",
            ),
        ]
    )

    if policy.get("block_when_history_adjustment_blockers_present"):
        checks.append(
            (not bool(candidate.get("history_adjustment_blockers")), "no_history_blockers")
        )
    if policy.get("block_when_direct_colonization_or_not_infection"):
        colonization = candidate.get("colonization_interpretation") or {}
        checks.append(
            (
                not bool(colonization.get("has_explicit_colonization"))
                and not bool(colonization.get("has_explicit_current_not_infection")),
                "no_direct_colonization_or_not_infection",
            )
        )

    for ok, reason in checks:
        (passed if ok else blockers).append(reason)
    return not blockers, passed, blockers


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["patient_id", "organism_name", "promoted_to_picked"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def routed_decision(candidate: dict[str, Any]) -> str:
    """Return the latest upstream decision while preserving legacy inputs."""
    return str(
        candidate.get("history_routed_decision")
        or candidate.get("decision")
        or "hold"
    )


def run(history_root: Path, output_root: Path, policy_path: Path = DEFAULT_POLICY) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    policy = read_json(policy_path)
    patient_output = output_root / "patient_outputs"
    patient_output.mkdir(parents=True)
    counts: Counter[str] = Counter()
    audit_rows: list[dict[str, Any]] = []

    for source_path in sorted((history_root / "patient_outputs").glob("*.json")):
        patient = patient_number(source_path)
        payload = read_json(source_path)
        output_candidates: list[dict[str, Any]] = []
        for candidate in payload.get("candidates") or []:
            candidate = dict(candidate)
            upstream_decision = routed_decision(candidate)
            already_picked = upstream_decision == "picked_shadow"
            eligible, passed, blockers = promotion_reasons(candidate, payload, policy)
            promoted = bool(eligible and not already_picked)
            final_decision = (
                "picked_shadow" if already_picked or promoted else upstream_decision
            )
            candidate["legacy_history_picked_gate"] = {
                "policy_id": policy["policy_id"],
                "answer_blind": True,
                "eligible": eligible,
                "promoted_to_picked": promoted,
                "passed_checks": passed,
                "blockers": blockers,
                "consumed_history_evidence_ids": (
                    list(candidate.get("consumed_evidence_ids") or []) if promoted else []
                ),
            }
            candidate["final_decision"] = final_decision
            output_candidates.append(candidate)
            counts["candidates"] += 1
            counts[f"final:{final_decision}"] += 1
            if promoted:
                counts["promoted_to_picked"] += 1
            audit_rows.append(
                {
                    "patient_id": patient,
                    "organism_name": candidate.get("organism_name"),
                    "baseline_decision": candidate.get("decision"),
                    "history_routed_decision": candidate.get("history_routed_decision"),
                    "analytic_route": candidate.get("analytic_route"),
                    "history_adjusted_route": candidate.get("history_adjusted_route"),
                    "taxonomy_family": candidate.get("taxonomy_family"),
                    "specimen_context": candidate.get("specimen_context"),
                    "best_rank": (candidate.get("analytical_profile") or {}).get(
                        "best_rank_in_retained_universe"
                    ),
                    "promoted_to_picked": promoted,
                    "final_decision": final_decision,
                    "passed_checks": "|".join(passed),
                    "blockers": "|".join(blockers),
                }
            )

        output = {
            "schema_version": SCHEMA_VERSION,
            "answer_blind": True,
            "patient_id": str(patient),
            "policy_id": policy["policy_id"],
            "source_history_output": str(source_path.resolve()),
            "source_history_sha256": sha256_file(source_path),
            "candidates": output_candidates,
        }
        write_json(
            patient_output / f"NGS_patient_{patient}_legacy_history_picked_gate.json", output
        )
        counts["patients"] += 1

    write_csv(output_root / "promotion_decisions.csv", audit_rows)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "answer_blind": True,
        "policy_id": policy["policy_id"],
        "policy_path": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "source_history_root": str(history_root.resolve()),
        "counts": dict(sorted(counts.items())),
    }
    write_json(output_root / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history_root", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.history_root, args.output_root, args.policy),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

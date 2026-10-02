"""Shadow-audit a staged, human-reviewed taxonomy overlay before applying it."""

from __future__ import annotations

import argparse
import copy
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any
from unittest.mock import patch

from tools import organism_taxonomy_classifier as taxonomy
from tools import pathogen_normalization as pathogen_names
from tools.build_independent_multi_assay_screening_shadow import build_shadow
from tools.build_manual_style_multi_assay_screening_shadow import (
    DEFAULT_POLICY,
    route_candidate,
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _combined_rules(staged_overlay: Path) -> tuple[dict[str, Any], list[str]]:
    overlay = json.loads(staged_overlay.read_text(encoding="utf-8-sig"))
    if not isinstance(overlay, dict) or not isinstance(overlay.get("exact_profiles"), dict):
        raise ValueError("Staged overlay must contain an exact_profiles object")
    combined = copy.deepcopy(taxonomy.rules_payload())
    changed_keys: list[str] = []
    for key, profile in overlay["exact_profiles"].items():
        existing = combined["exact_profiles"].get(key)
        if existing is not None and existing != profile:
            raise ValueError(f"Staged overlay conflicts with an existing exact profile: {key}")
        if existing is None:
            changed_keys.append(key)
        combined["exact_profiles"][key] = profile
    taxonomy._validate_rules(combined)
    return combined, sorted(changed_keys)


def _is_reviewed_candidate(row: dict[str, str], reviewed_keys: set[str]) -> bool:
    """Match reviewed profiles through the same alias normalization as taxonomy.

    Candidate-entry rows retain the source organism key for provenance, while
    exact taxonomy profiles are keyed by the canonical alias.  Comparing only
    the source key incorrectly labels a reviewed alias (for example,
    ``Metamycoplasma hominis`` -> ``mycoplasmahominis``) as an unrelated route
    change.
    """

    source_key = str(row.get("organism_key") or "")
    canonical_key = pathogen_names.canonical_key(row.get("organism_name"))
    return source_key in reviewed_keys or canonical_key in reviewed_keys


def audit(
    staged_overlay: Path,
    inventory_dir: Path,
    baseline_dir: Path,
    output_dir: Path,
    policy_path: Path = DEFAULT_POLICY,
) -> dict[str, Any]:
    combined, reviewed_keys = _combined_rules(staged_overlay)
    baseline_path = baseline_dir / "all_decisions.csv"
    if not baseline_path.exists():
        raise FileNotFoundError(f"Missing baseline decisions: {baseline_path}")

    shadow_dir = output_dir / "shadow"
    with patch.object(taxonomy, "rules_payload", return_value=combined):
        shadow_summary = build_shadow(
            inventory_dir,
            shadow_dir,
            policy_path,
            router=route_candidate,
        )

    key_fields = ("patient_id", "case_review_id", "organism_key")
    baseline_rows = _read_csv(baseline_path)
    shadow_rows = _read_csv(shadow_dir / "all_decisions.csv")
    baseline_by_key = {tuple(row[field] for field in key_fields): row for row in baseline_rows}
    shadow_by_key = {tuple(row[field] for field in key_fields): row for row in shadow_rows}
    if baseline_by_key.keys() != shadow_by_key.keys():
        added = sorted(shadow_by_key.keys() - baseline_by_key.keys())
        removed = sorted(baseline_by_key.keys() - shadow_by_key.keys())
        raise ValueError(f"Candidate identity changed during taxonomy audit: added={added}, removed={removed}")

    changes: list[dict[str, Any]] = []
    watched = (
        "taxonomy_family",
        "taxonomy_mapping_status",
        "disposition",
        "screening_tier",
        "rule_ids",
    )
    reviewed_key_set = set(reviewed_keys)
    for key in sorted(baseline_by_key):
        old = baseline_by_key[key]
        new = shadow_by_key[key]
        if not any(old[field] != new[field] for field in watched):
            continue
        changes.append(
            {
                "patient_id": key[0],
                "case_review_id": key[1],
                "organism_name": new["organism_name"],
                "organism_key": key[2],
                "reviewed_key": _is_reviewed_candidate(new, reviewed_key_set),
                "old_taxonomy_family": old["taxonomy_family"],
                "new_taxonomy_family": new["taxonomy_family"],
                "old_mapping_status": old["taxonomy_mapping_status"],
                "new_mapping_status": new["taxonomy_mapping_status"],
                "old_disposition": old["disposition"],
                "new_disposition": new["disposition"],
                "old_screening_tier": old["screening_tier"],
                "new_screening_tier": new["screening_tier"],
                "old_rule_ids": old["rule_ids"],
                "new_rule_ids": new["rule_ids"],
            }
        )

    fields = list(changes[0]) if changes else [
        "patient_id",
        "case_review_id",
        "organism_name",
        "organism_key",
        "reviewed_key",
        "old_taxonomy_family",
        "new_taxonomy_family",
        "old_mapping_status",
        "new_mapping_status",
        "old_disposition",
        "new_disposition",
        "old_screening_tier",
        "new_screening_tier",
        "old_rule_ids",
        "new_rule_ids",
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "route_changes.csv", changes, fields)
    unexpected = [row for row in changes if not row["reviewed_key"]]
    summary = {
        "mode": "staged_reviewed_taxonomy_shadow",
        "staged_overlay": str(staged_overlay.resolve()),
        "reviewed_exact_profiles": reviewed_keys,
        "candidate_count": len(shadow_rows),
        "route_changes": len(changes),
        "unexpected_nonreviewed_changes": len(unexpected),
        "changes_by_disposition": dict(
            sorted(Counter(f"{row['old_disposition']} -> {row['new_disposition']}" for row in changes).items())
        ),
        "shadow_summary": shadow_summary,
        "safe_to_consider_apply": not unexpected,
    }
    (output_dir / "audit_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit a staged reviewed taxonomy overlay against a frozen candidate-entry run.")
    parser.add_argument("staged_overlay", type=Path)
    parser.add_argument("inventory_dir", type=Path)
    parser.add_argument("baseline_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    print(
        json.dumps(
            audit(args.staged_overlay, args.inventory_dir, args.baseline_dir, args.output_dir, args.policy),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

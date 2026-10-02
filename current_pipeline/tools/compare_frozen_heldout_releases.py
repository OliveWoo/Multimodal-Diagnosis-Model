"""Compare two frozen held-out releases without reading outcome labels."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def key(row: dict[str, str]) -> tuple[str, str]:
    return row["patient_id"], row["organism_name"].strip().casefold()


def compare(reference: Path, current: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    ref_freeze = read_json(reference / "answer_blind_freeze.json")
    cur_freeze = read_json(current / "answer_blind_freeze.json")
    ref_ids = ref_freeze["heldout_subset"]["heldout_patient_ids"]
    cur_ids = cur_freeze["heldout_subset"]["heldout_patient_ids"]
    if ref_ids != cur_ids:
        raise ValueError("Held-out patient sets differ")
    ref_path = reference / "heldout_subset/heldout_integrated_decisions.csv"
    cur_path = current / "heldout_subset/heldout_integrated_decisions.csv"
    ref_rows = {key(row): row for row in read_csv(ref_path)}
    cur_rows = {key(row): row for row in read_csv(cur_path)}
    fields = [
        "taxonomy_family",
        "taxonomy_mapping_status",
        "analytical_route",
        "history_adjusted_route",
        "post_promotion_clinical_decision",
        "promotion_outcome",
        "final_reporting_tier",
        "reporting_role",
        "selected_for_complete_report",
    ]
    changes: list[dict[str, Any]] = []
    for item_key in sorted(ref_rows.keys() | cur_rows.keys()):
        before = ref_rows.get(item_key)
        after = cur_rows.get(item_key)
        patient_id, organism_key = item_key
        if before is None:
            changes.append({
                "patient_id": patient_id,
                "organism_name": after["organism_name"],
                "change_type": "added",
                "changed_fields": "",
                "before_tier": "",
                "after_tier": after["final_reporting_tier"],
            })
            continue
        if after is None:
            changes.append({
                "patient_id": patient_id,
                "organism_name": before["organism_name"],
                "change_type": "removed",
                "changed_fields": "",
                "before_tier": before["final_reporting_tier"],
                "after_tier": "",
            })
            continue
        changed = [field for field in fields if before.get(field) != after.get(field)]
        if changed:
            changes.append({
                "patient_id": patient_id,
                "organism_name": after["organism_name"],
                "change_type": "field_change",
                "changed_fields": "|".join(changed),
                "before_tier": before["final_reporting_tier"],
                "after_tier": after["final_reporting_tier"],
            })
    ref_report = reference / "heldout_subset/heldout_complete_report.csv"
    cur_report = current / "heldout_subset/heldout_complete_report.csv"
    summary = {
        "schema_version": "frozen_heldout_release_comparison.v1",
        "answer_blind": True,
        "heldout_patient_ids": ref_ids,
        "heldout_patient_count": len(ref_ids),
        "reference_root": str(reference.resolve()),
        "current_root": str(current.resolve()),
        "reference_integrated_candidate_count": len(ref_rows),
        "current_integrated_candidate_count": len(cur_rows),
        "change_count": len(changes),
        "change_type_counts": dict(Counter(row["change_type"] for row in changes)),
        "complete_report_identical": sha256_file(ref_report) == sha256_file(cur_report),
        "reference_complete_report_sha256": sha256_file(ref_report),
        "current_complete_report_sha256": sha256_file(cur_report),
        "reference_answers_loaded": ref_freeze.get("answers_loaded"),
        "current_answers_loaded": cur_freeze.get("answers_loaded"),
    }
    return summary, changes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("current", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    summary, changes = compare(args.reference, args.current)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    with (args.output_dir / "candidate_changes.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        fields = [
            "patient_id",
            "organism_name",
            "change_type",
            "changed_fields",
            "before_tier",
            "after_tier",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(changes)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

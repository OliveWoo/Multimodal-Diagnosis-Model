from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import pathogen_normalization as pathogen_names  # noqa: E402


def canonicalize_name(name: str) -> str:
    return pathogen_names.canonical_key(name)


def load_best(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    best = (payload or {}).get("best_available_summary", {}) or {}
    picked = best.get("picked_pathogens", []) or []
    return {
        "selection_mode": best.get("selection_mode"),
        "picked_count": best.get("picked_count"),
        "picked_items": [
            {
                "organism_name": (item.get("organism_name") or "").strip(),
                "basis_level": item.get("basis_level"),
                "picked_role": item.get("picked_role"),
                "basis_confidence": item.get("basis_confidence"),
            }
            for item in picked
            if (item.get("organism_name") or "").strip()
        ],
    }


def build_norm_index(items: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    index: dict[str, list[dict[str, Any]]] = {}
    for item in items:
        key = canonicalize_name(item["organism_name"])
        index.setdefault(key, []).append(item)
    return index


def compare_pair(no_filmarray: dict[str, Any], full: dict[str, Any]) -> dict[str, Any]:
    no_items = no_filmarray["picked_items"]
    full_items = full["picked_items"]

    no_exact = {x["organism_name"]: x for x in no_items}
    full_exact = {x["organism_name"]: x for x in full_items}

    no_exact_names = set(no_exact.keys())
    full_exact_names = set(full_exact.keys())

    exact_added_in_full = sorted(full_exact_names - no_exact_names)
    exact_removed_in_full = sorted(no_exact_names - full_exact_names)

    no_norm = build_norm_index(no_items)
    full_norm = build_norm_index(full_items)
    no_norm_keys = set(no_norm.keys())
    full_norm_keys = set(full_norm.keys())

    normalized_added_in_full = sorted(full_norm_keys - no_norm_keys)
    normalized_removed_in_full = sorted(no_norm_keys - full_norm_keys)

    naming_variant_matches: list[dict[str, Any]] = []
    for key in sorted(no_norm_keys & full_norm_keys):
        no_names = sorted({x["organism_name"] for x in no_norm[key]})
        full_names = sorted({x["organism_name"] for x in full_norm[key]})
        if no_names != full_names:
            naming_variant_matches.append(
                {
                    "canonical_name": key,
                    "no_filmarray_names": no_names,
                    "full_names": full_names,
                }
            )

    normalized_level_role_changes: list[dict[str, Any]] = []
    for key in sorted(no_norm_keys & full_norm_keys):
        no_group = no_norm[key]
        full_group = full_norm[key]
        if len(no_group) != 1 or len(full_group) != 1:
            continue
        n = no_group[0]
        f = full_group[0]
        if n["basis_level"] != f["basis_level"] or n["picked_role"] != f["picked_role"]:
            normalized_level_role_changes.append(
                {
                    "canonical_name": key,
                    "no_filmarray_name": n["organism_name"],
                    "full_name": f["organism_name"],
                    "no_filmarray_level": n["basis_level"],
                    "full_level": f["basis_level"],
                    "no_filmarray_role": n["picked_role"],
                    "full_role": f["picked_role"],
                }
            )

    return {
        "exact_diff": {
            "added_in_full": exact_added_in_full,
            "removed_in_full": exact_removed_in_full,
        },
        "name_aware_diff": {
            "normalized_added_in_full": normalized_added_in_full,
            "normalized_removed_in_full": normalized_removed_in_full,
            "naming_variant_matches": naming_variant_matches,
            "normalized_level_role_changes": normalized_level_role_changes,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare best_available_summary between no_filmarray and full with name-aware matching."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("outputs") / "patient_info_filtered",
        help="Root directory containing NGS_patient_*_json folders.",
    )
    parser.add_argument(
        "--out-json",
        type=Path,
        default=Path("outputs") / "best_available_no_filmarray_vs_full_ALL_nameaware_comparison.json",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("outputs") / "best_available_no_filmarray_vs_full_ALL_nameaware_comparison.csv",
    )
    args = parser.parse_args()

    rows: list[dict[str, Any]] = []
    for patient_dir in sorted(args.root.glob("NGS_patient_*_json")):
        if not patient_dir.is_dir():
            continue
        patient_id = patient_dir.name.replace("NGS_patient_", "").replace("_json", "")
        summary_dir = patient_dir / "summary_outputs"
        no_file = summary_dir / f"NGS_patient_{patient_id}_mNGS_max_agent_no_filmarray.json"
        full_file = summary_dir / f"NGS_patient_{patient_id}_mNGS_max_agent_full.json"
        if not (no_file.exists() and full_file.exists()):
            continue

        no_best = load_best(no_file)
        full_best = load_best(full_file)
        diff = compare_pair(no_best, full_best)

        rows.append(
            {
                "patient_id": patient_id,
                "no_filmarray_file": str(no_file),
                "full_file": str(full_file),
                "no_filmarray": no_best,
                "full": full_best,
                "diff": diff,
            }
        )

    report = {
        "status": "ok",
        "patient_count": len(rows),
        "patients": rows,
        "notes": {
            "name_aware_matching": "Uses canonicalized pathogen names to reduce false diffs from naming variants.",
            "example": "Acinetobacter baumannii complex ~= Acinetobacter baumannii",
        },
    }
    args.out_json.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "patient_id",
                "no_filmarray_picked_count",
                "full_picked_count",
                "exact_added_in_full",
                "exact_removed_in_full",
                "normalized_added_in_full",
                "normalized_removed_in_full",
                "naming_variant_matches",
                "normalized_level_role_changes",
            ]
        )
        for row in rows:
            diff = row["diff"]
            writer.writerow(
                [
                    row["patient_id"],
                    row["no_filmarray"]["picked_count"],
                    row["full"]["picked_count"],
                    "; ".join(diff["exact_diff"]["added_in_full"]),
                    "; ".join(diff["exact_diff"]["removed_in_full"]),
                    "; ".join(diff["name_aware_diff"]["normalized_added_in_full"]),
                    "; ".join(diff["name_aware_diff"]["normalized_removed_in_full"]),
                    json.dumps(diff["name_aware_diff"]["naming_variant_matches"], ensure_ascii=False),
                    json.dumps(diff["name_aware_diff"]["normalized_level_role_changes"], ensure_ascii=False),
                ]
            )

    print(f"JSON: {args.out_json}")
    print(f"CSV:  {args.out_csv}")
    print(f"Compared patients: {len(rows)}")


if __name__ == "__main__":
    main()

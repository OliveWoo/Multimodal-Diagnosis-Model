#!/usr/bin/env python3
"""Merge unchanged validated rationales with regenerated R5 rationales."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def selected_from_result(result: dict[str, Any]) -> list[tuple[str, int]]:
    rows = (result.get("generated_rationale") or {}).get("selected_pathogen_explanations") or []
    return sorted(((row["organism"], int(row["rank"])) for row in rows), key=lambda row: row[1])


def selected_from_decision(decision: dict[str, Any]) -> list[tuple[str, int]]:
    return sorted(((row["organism"], int(row["rank"])) for row in decision.get("selected") or []), key=lambda row: row[1])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--final-decisions", required=True, type=Path)
    regenerated_group = parser.add_mutually_exclusive_group(required=True)
    regenerated_group.add_argument(
        "--regenerated",
        type=Path,
        help="One regenerated batch; patient IDs must be unique across cohorts.",
    )
    regenerated_group.add_argument(
        "--regenerated-cohort",
        action="append",
        default=[],
        metavar="COHORT=PATH",
        help="Regenerated batch for one cohort; repeat for multi-cohort runs.",
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    config_file = args.config.resolve()
    config_root = config_file.parent
    config = read_json(config_file)
    base_file = (config_root / config["accepted_rationales"]).resolve()
    base = read_json(base_file)
    entries: dict[tuple[str, str], dict[str, Any]] = {}
    for row in base.get("patients") or []:
        patient_id = str(row["result"]["patient_id"])
        entries[(row["cohort"], patient_id)] = row

    # Apply any already-approved per-patient rationale migrations first.
    for identity, relative_path in (config.get("rationale_overrides") or {}).items():
        cohort, patient_id = identity.split("/", 1)
        result_file = (config_root / relative_path).resolve()
        result = read_json(result_file)
        if str(result["patient_id"]) != patient_id:
            raise ValueError(f"rationale override patient mismatch: {identity}")
        entries[(cohort, patient_id)] = {
            "cohort": cohort,
            "output_file": str(result_file),
            "result": result,
        }

    regenerated_batches: list[tuple[str | None, Path]] = []
    if args.regenerated is not None:
        regenerated_batches.append((None, args.regenerated.resolve()))
    else:
        seen_cohorts: set[str] = set()
        for value in args.regenerated_cohort:
            if "=" not in value:
                raise ValueError("--regenerated-cohort must be COHORT=PATH")
            cohort, path_value = value.split("=", 1)
            if not cohort or cohort in seen_cohorts:
                raise ValueError(f"duplicate or empty regenerated cohort: {cohort!r}")
            if cohort not in config.get("cohorts", {}):
                raise ValueError(f"regenerated cohort is absent from config: {cohort}")
            seen_cohorts.add(cohort)
            regenerated_batches.append((cohort, Path(path_value).resolve()))

    regenerated_count = 0
    regenerated_provenance: list[dict[str, Any]] = []
    replaced: set[tuple[str, str]] = set()
    for cohort_hint, regenerated_root in regenerated_batches:
        regenerated_manifest = read_json(regenerated_root / "llm_rationale_manifest.json")
        if regenerated_manifest.get("failures"):
            raise ValueError(f"regenerated rationale batch contains failures: {regenerated_root}")
        batch_count = 0
        for row in regenerated_manifest.get("patients") or []:
            result_file = Path(row["output_file"]).resolve()
            result = read_json(result_file)
            patient_id = str(result["patient_id"])
            if cohort_hint is None:
                matches = [key for key in entries if key[1] == patient_id]
                if len(matches) != 1:
                    raise ValueError(f"cannot uniquely map regenerated P{patient_id} to cohort")
                cohort = matches[0][0]
            else:
                cohort = cohort_hint
                if (cohort, patient_id) not in entries:
                    raise ValueError(f"regenerated patient is absent from base collection: {cohort}/P{patient_id}")
            key = (cohort, patient_id)
            if key in replaced:
                raise ValueError(f"duplicate regenerated rationale: {cohort}/P{patient_id}")
            replaced.add(key)
            entries[key] = {
                "cohort": cohort,
                "output_file": str(result_file),
                "result": result,
            }
            regenerated_count += 1
            batch_count += 1
        regenerated_provenance.append(
            {
                "cohort": cohort_hint,
                "path": str(regenerated_root),
                "patients": batch_count,
            }
        )

    decisions_root = args.final_decisions.resolve()
    expected_files = list(decisions_root.glob("*/*_final_decision.json"))
    if len(entries) != len(expected_files):
        raise ValueError(f"patient count mismatch: rationales={len(entries)}, decisions={len(expected_files)}")
    for (cohort, patient_id), row in entries.items():
        decision_file = decisions_root / cohort / f"P{patient_id}_final_decision.json"
        decision = read_json(decision_file)
        result = row["result"]
        if str(result["patient_id"]) != patient_id:
            raise ValueError(f"result patient mismatch: {cohort}/P{patient_id}")
        if selected_from_result(result) != selected_from_decision(decision):
            raise ValueError(
                f"selection mismatch {cohort}/P{patient_id}: "
                f"rationale={selected_from_result(result)}, decision={selected_from_decision(decision)}"
            )
        if json.dumps(result, ensure_ascii=False, sort_keys=True) != json.dumps(
            read_json(Path(row["output_file"])), ensure_ascii=False, sort_keys=True
        ):
            raise ValueError(f"embedded rationale differs from file: {cohort}/P{patient_id}")

    ordered = [
        entries[key]
        for key in sorted(entries, key=lambda item: (item[0], int(item[1])))
    ]
    collection = {
        "status": "structurally_validated",
        "clinical_review_status": "not_expert_validated",
        "schema_version": "ober.accepted_rationales.r5.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "provenance": {
            "base_collection": str(base_file),
            "regenerated_batches": regenerated_provenance,
            "final_decisions_root": str(decisions_root),
            "regenerated_patients": regenerated_count,
            "unchanged_or_migrated_patients": len(ordered) - regenerated_count,
        },
        "patients": ordered,
    }
    write_json(args.output.resolve(), collection)
    print(
        json.dumps(
            {
                "patients": len(ordered),
                "regenerated": regenerated_count,
                "selected_reconciled": sum(len(selected_from_result(row["result"])) for row in ordered),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

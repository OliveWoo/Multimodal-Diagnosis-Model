"""Fail closed unless every unmapped organism has a reviewed taxonomy outcome.

This validator is deliberately organism-only.  It verifies the artifacts made
by the NCBI resolver, two LLM reviewers, deterministic validator, and route
audit without reading patient answers or using patient-level evidence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from tools import pathogen_normalization


ALLOWED_STATUSES = {"auto_approved", "provisional_unmapped", "rejected"}


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            payload = json.loads(line)
            if not isinstance(payload, dict):
                raise ValueError(f"Queue row {line_number} is not an object")
            rows.append(payload)
    return rows


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical(value: Any) -> str:
    return pathogen_normalization.canonical_key(value)


def validate(
    baseline_csv: Path,
    queue_path: Path,
    decision_dir: Path,
    staged_overlay: Path,
    audit_summary_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    baseline_rows = _read_csv(baseline_csv)
    baseline_unmapped = {
        _canonical(row.get("organism_name"))
        for row in baseline_rows
        if row.get("taxonomy_mapping_status") == "unmapped"
    }
    baseline_unmapped.discard("")

    queue = _read_jsonl(queue_path)
    queue_keys = {str(row.get("canonical_key") or "") for row in queue}
    queue_keys.discard("")

    run_summary_path = decision_dir / "run_summary.json"
    run_summary = _read_json(run_summary_path)
    decision_paths = sorted((decision_dir / "decisions").glob("*.json"))
    decisions = [_read_json(path) for path in decision_paths]
    decision_by_key = {
        str(row.get("canonical_key") or ""): row
        for row in decisions
        if row.get("canonical_key")
    }
    decision_keys = set(decision_by_key)
    invalid_statuses = {
        key: row.get("status")
        for key, row in decision_by_key.items()
        if row.get("status") not in ALLOWED_STATUSES
    }

    overlay = _read_json(staged_overlay)
    overlay_keys = set((overlay.get("exact_profiles") or {}).keys())
    approved_keys = {
        key
        for key, row in decision_by_key.items()
        if row.get("status") == "auto_approved"
    }
    nonapproved_keys = decision_keys - approved_keys
    audit = _read_json(audit_summary_path)

    errors: list[str] = []
    if baseline_unmapped != queue_keys:
        errors.append(
            "Queue does not exactly cover baseline unmapped organisms: "
            f"missing={sorted(baseline_unmapped - queue_keys)}, "
            f"unexpected={sorted(queue_keys - baseline_unmapped)}"
        )
    if run_summary.get("mode") != "executed":
        errors.append("Taxonomy auto-decision was not executed")
    if run_summary.get("completion_status") != "complete":
        errors.append(f"Taxonomy auto-decision is not complete: {run_summary.get('completion_status')!r}")
    if int(run_summary.get("failures") or 0) != 0:
        errors.append(f"Taxonomy auto-decision reports failures={run_summary.get('failures')}")
    if run_summary.get("pending_keys"):
        errors.append(f"Taxonomy auto-decision has pending keys: {run_summary.get('pending_keys')}")
    if int(run_summary.get("remaining_eligible_not_selected") or 0) != 0:
        errors.append("Not every eligible organism was selected for review")
    if decision_keys != queue_keys:
        errors.append(
            "Decision set does not exactly cover the queue: "
            f"missing={sorted(queue_keys - decision_keys)}, "
            f"unexpected={sorted(decision_keys - queue_keys)}"
        )
    if invalid_statuses:
        errors.append(f"Invalid decision statuses: {invalid_statuses}")
    if not approved_keys.issubset(overlay_keys):
        errors.append(f"Approved profiles missing from staged overlay: {sorted(approved_keys - overlay_keys)}")
    if nonapproved_keys & overlay_keys:
        errors.append(
            "Held/rejected organisms unexpectedly appear in the staged overlay: "
            f"{sorted(nonapproved_keys & overlay_keys)}"
        )
    if audit.get("safe_to_consider_apply") is not True:
        errors.append("Route-impact audit did not pass")
    if Path(str(audit.get("staged_overlay") or "")).resolve() != staged_overlay.resolve():
        errors.append("Route-impact audit refers to a different staged overlay")

    status_counts = Counter(str(row.get("status") or "unknown") for row in decisions)
    manifest = {
        "schema_version": "organism_taxonomy_completion_manifest_v1",
        "passed": not errors,
        "answer_blind": True,
        "patient_level_fields_used": False,
        "ground_truth_loaded": False,
        "baseline_csv": str(baseline_csv.resolve()),
        "baseline_sha256": _sha256(baseline_csv),
        "queue_path": str(queue_path.resolve()),
        "queue_sha256": _sha256(queue_path),
        "run_summary": str(run_summary_path.resolve()),
        "run_summary_sha256": _sha256(run_summary_path),
        "staged_overlay": str(staged_overlay.resolve()),
        "staged_overlay_sha256": _sha256(staged_overlay),
        "audit_summary": str(audit_summary_path.resolve()),
        "audit_summary_sha256": _sha256(audit_summary_path),
        "baseline_unmapped_count": len(baseline_unmapped),
        "queued_count": len(queue_keys),
        "decided_count": len(decision_keys),
        "status_counts": dict(sorted(status_counts.items())),
        "auto_approved_count": len(approved_keys),
        "held_count": len(nonapproved_keys),
        "unreviewed_keys": sorted(baseline_unmapped - decision_keys),
        "errors": errors,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if errors:
        raise RuntimeError("Taxonomy completion validation failed: " + "; ".join(errors))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-csv", type=Path, required=True)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--decision-dir", type=Path, required=True)
    parser.add_argument("--staged-overlay", type=Path, required=True)
    parser.add_argument("--audit-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            validate(
                args.baseline_csv,
                args.queue,
                args.decision_dir,
                args.staged_overlay,
                args.audit_summary,
                args.output,
            ),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

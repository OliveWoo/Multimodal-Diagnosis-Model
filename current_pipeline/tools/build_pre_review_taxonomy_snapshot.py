"""Reconstruct the pre-double-review auto-taxonomy overlay.

The current overlay contains the profiles that existed before the 97-name
review plus the auto-approved review additions.  The decision summary is the
authoritative list of newly approved keys, so removing exactly those keys
reconstructs the frozen comparator without consulting patient outcomes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CURRENT = ROOT / "rules/organism_taxonomy_auto_v1.json"
DEFAULT_DECISIONS = (
    ROOT
    / "outputs/runs/2026-09-30_KH_taxonomy_auto_decision_v2_4/decision_summary.csv"
)
DEFAULT_OUTPUT = (
    ROOT
    / "outputs/runs/2026-09-30_taxonomy_pre_double_review_snapshot_v2_4"
    / "organism_taxonomy_auto_v2_4_frozen.json"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build(current: Path, decisions: Path, output: Path) -> dict[str, Any]:
    payload = json.loads(current.read_text(encoding="utf-8-sig"))
    with decisions.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    approved = {
        row["canonical_key"]
        for row in rows
        if row.get("status") == "auto_approved"
    }
    profiles = dict(payload.get("exact_profiles") or {})
    missing = sorted(approved - profiles.keys())
    if missing:
        raise ValueError(f"Approved keys missing from current overlay: {missing}")
    retained = {key: value for key, value in profiles.items() if key not in approved}
    if len(approved) != 21 or len(profiles) != 32 or len(retained) != 11:
        raise ValueError(
            "Expected 32 current profiles - 21 review additions = 11 comparator profiles; "
            f"got current={len(profiles)}, approved={len(approved)}, retained={len(retained)}"
        )
    payload["schema_version"] = "organism_taxonomy_auto_v2_4_frozen"
    payload["purpose"] = (
        "Pre-double-review comparator reconstructed answer-blind by removing "
        "the 21 keys recorded as auto_approved in the 97-name decision summary."
    )
    payload["exact_profiles"] = retained
    payload["comparator_provenance"] = {
        "source_overlay": str(current.resolve()),
        "source_overlay_sha256": sha256_file(current),
        "decision_summary": str(decisions.resolve()),
        "decision_summary_sha256": sha256_file(decisions),
        "removed_auto_approved_keys": sorted(approved),
        "removed_count": len(approved),
        "retained_count": len(retained),
        "patient_answers_used": False,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return {
        "output": str(output.resolve()),
        "output_sha256": sha256_file(output),
        "current_profile_count": len(profiles),
        "removed_auto_approved_count": len(approved),
        "retained_profile_count": len(retained),
        "patient_answers_used": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current", type=Path, default=DEFAULT_CURRENT)
    parser.add_argument("--decisions", type=Path, default=DEFAULT_DECISIONS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(
        build(args.current, args.decisions, args.output),
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()

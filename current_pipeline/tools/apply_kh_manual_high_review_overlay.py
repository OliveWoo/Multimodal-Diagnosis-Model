"""Apply explicit human High-review flags without changing deterministic tiers."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "outputs/runs/2026-09-29_KH_integrated_decision_pipeline_v1/integrated_decisions.csv"
DEFAULT_POLICY = ROOT / "rules/kh_manual_high_review_v1.json"
DEFAULT_OUTPUT = ROOT / "outputs/runs/2026-09-30_KH_manual_high_review_overlay_v1"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _deterministic_tier(row: dict[str, str]) -> str:
    return {
        "picked_shadow": "Picked",
        "high_priority_review": "High",
        "review_context_needed": "Context",
    }.get(row.get("post_promotion_clinical_decision", ""), "Unknown")


def apply(rows: list[dict[str, str]], policy: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    review_items = {
        (str(item["patient_id"]), str(item["organism_name"]).casefold()): item
        for item in policy.get("review_items", [])
    }
    matched: set[tuple[str, str]] = set()
    output: list[dict[str, Any]] = []
    queue: list[dict[str, Any]] = []
    for source in rows:
        row: dict[str, Any] = dict(source)
        key = (str(row.get("patient_id") or ""), str(row.get("organism_name") or "").casefold())
        item = review_items.get(key)
        deterministic_tier = _deterministic_tier(source)
        row["deterministic_review_tier"] = deterministic_tier
        row["manual_high_review_requested"] = bool(item)
        row["manual_review_tier"] = item.get("manual_review_tier", "") if item else ""
        row["effective_review_tier"] = item.get("manual_review_tier", deterministic_tier) if item else deterministic_tier
        row["manual_review_decision_source"] = item.get("decision_source", "") if item else ""
        row["manual_review_reason"] = item.get("reason", "") if item else ""
        row["manual_review_can_create_picked"] = False
        if item:
            matched.add(key)
            queue.append(
                {
                    "patient_id": row["patient_id"],
                    "organism_name": row["organism_name"],
                    "deterministic_review_tier": deterministic_tier,
                    "deterministic_reporting_tier": row.get("final_reporting_tier", ""),
                    "manual_review_tier": item["manual_review_tier"],
                    "taxonomy_family": row.get("taxonomy_family", ""),
                    "taxonomy_mapping_status": row.get("taxonomy_mapping_status", ""),
                    "best_rank": row.get("best_rank", ""),
                    "selected_positive_test_count": row.get("selected_positive_test_count", ""),
                    "cross_molecule_selected": row.get("cross_molecule_selected", ""),
                    "direct_hospital_level": row.get("direct_hospital_level", ""),
                    "reason": item["reason"],
                    "decision_source": item["decision_source"],
                    "automatic_picked_blocked": True,
                }
            )
        output.append(row)

    requested = set(review_items)
    missing = sorted(requested - matched)
    assertions = {
        "all_requested_candidates_matched": not missing,
        "row_count_preserved": len(output) == len(rows),
        "deterministic_fields_preserved": all(
            all(out.get(field) == src.get(field) for field in src)
            for out, src in zip(output, rows)
        ),
        "no_manual_item_can_create_picked": all(not row["manual_review_can_create_picked"] for row in output),
    }
    summary = {
        "schema_version": "kh_manual_high_review_overlay.v1",
        "role": "manual_review_queue_only",
        "input_candidate_count": len(rows),
        "requested_high_review_count": len(requested),
        "matched_high_review_count": len(matched),
        "missing_requested_candidates": [list(value) for value in missing],
        "queue": queue,
        "metric_policy": "Do not attribute manual High flags to deterministic model performance.",
        "assertions": assertions,
        "all_assertions_pass": all(assertions.values()),
    }
    return output, queue, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    rows = _read_csv(args.input)
    policy = json.loads(args.policy.read_text(encoding="utf-8"))
    output, queue, summary = apply(rows, policy)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "integrated_decisions_with_manual_review.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    with (args.output_dir / "manual_high_review_queue.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(queue[0]))
        writer.writeheader()
        writer.writerows(queue)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

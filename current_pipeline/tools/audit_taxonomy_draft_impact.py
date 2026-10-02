"""Audit unapproved organism drafts and their hypothetical screening impact."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from unittest.mock import patch

from tools import organism_taxonomy_classifier as taxonomy
from tools.build_test_aware_mngs_screening_shadow import DEFAULT_POLICY, build_shadow
from tools.draft_organism_taxonomy_reviews import _validate_proposal, load_queue, prompt_for


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def audit(
    queue_path: Path,
    draft_dir: Path,
    literature_path: Path,
    inventory_dir: Path,
    baseline_dir: Path,
    output_dir: Path,
    policy_path: Path = DEFAULT_POLICY,
) -> dict[str, object]:
    queue = load_queue(queue_path)
    literature = json.loads(literature_path.read_text(encoding="utf-8"))
    baseline_path = baseline_dir / "screening_shadow.csv"
    baseline = _read_csv(baseline_path)
    baseline_by_key = {(row["patient_id"], row["organism_key"]): row for row in baseline}
    candidate_counts = Counter(row["organism_key"] for row in baseline if row["comparison_status"] == "new_from_multi_assay_union")

    review_rows: list[dict[str, object]] = []
    hypothetical_rules = copy.deepcopy(taxonomy.rules_payload())
    for item in queue:
        key = item["canonical_key"]
        draft_path = draft_dir / f"{key}.json"
        draft = json.loads(draft_path.read_text(encoding="utf-8"))
        proposal = draft["proposal"]
        messages = prompt_for(item, literature.get(key))
        input_hash = hashlib.sha256(json.dumps(messages, ensure_ascii=False, sort_keys=True).encode("utf-8")).hexdigest()
        if draft.get("canonical_key") != key or draft.get("input_sha256") != input_hash:
            raise ValueError(f"Draft input mismatch: {key}")
        if draft.get("review_status") != "llm_draft_unverified":
            raise ValueError(f"Unexpected review status: {key}")
        available_pmids = _validate_proposal(proposal, item, literature.get(key))
        if sorted(available_pmids) != draft.get("retrieved_pmids"):
            raise ValueError(f"Retrieved literature mismatch: {key}")
        family = proposal["family_hypothesis"]
        review_rows.append({
            "organism_name": item["organism_name"],
            "canonical_key": key,
            "current_family": item["current_family"],
            "draft_family": family,
            "draft_confidence": proposal["confidence"],
            "identity_assessment": proposal["identity_assessment"],
            "candidate_patients": candidate_counts[key],
            "retrieved_pmid_count": len(available_pmids),
            "supporting_pmids": "|".join(proposal["supporting_pmids"]),
            "human_infection_evidence_as_drafted": proposal["human_infection_evidence"],
            "review_status": "human_review_required",
        })
        if family == "unmapped_or_uncertain":
            continue
        biological_class = next((value for value in item.get("input_classes", []) if value != "unknown"), "unknown")
        ncbi = item.get("ncbi") or {}
        hypothetical_rules["exact_profiles"][key] = {
            "rule_id": f"HYPOTHETICAL-UNAPPROVED-{key.upper()}",
            "biological_class": biological_class,
            "primary_rule_family": family,
            "clinical_traits": [],
            "literature_review_recommended": False,
            "taxid": ncbi.get("taxid"),
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    review_fields = list(review_rows[0]) if review_rows else []
    _write_csv(output_dir / "draft_review_queue.csv", review_rows, review_fields)
    with patch.object(taxonomy, "rules_payload", return_value=hypothetical_rules):
        shadow = build_shadow(inventory_dir, output_dir / "hypothetical_shadow", policy_path)

    hypothetical = _read_csv(output_dir / "hypothetical_shadow" / "screening_shadow.csv")
    changes: list[dict[str, object]] = []
    for row in hypothetical:
        key = (row["patient_id"], row["organism_key"])
        old = baseline_by_key.pop(key, None)
        if old is None:
            raise ValueError(f"New candidate appeared in hypothetical shadow: {key}")
        if row["screening_action"] == old["screening_action"]:
            continue
        changes.append({
            "patient_id": row["patient_id"],
            "organism_name": row["organism_name"],
            "organism_key": row["organism_key"],
            "molecule_support": row["molecule_support"],
            "detected_test_count": row["detected_test_count"],
            "baseline_family": old["primary_rule_family"],
            "hypothetical_family": row["primary_rule_family"],
            "baseline_action": old["screening_action"],
            "hypothetical_action": row["screening_action"],
        })
    if baseline_by_key:
        raise ValueError(f"Baseline candidates disappeared: {len(baseline_by_key)}")
    change_fields = list(changes[0]) if changes else [
        "patient_id", "organism_name", "organism_key", "molecule_support", "detected_test_count",
        "baseline_family", "hypothetical_family", "baseline_action", "hypothetical_action",
    ]
    _write_csv(output_dir / "patient_route_changes.csv", changes, change_fields)
    summary = {
        "mode": "unapproved_draft_impact_only",
        "drafts_checked": len(review_rows),
        "proposed_non_unmapped": sum(row["draft_family"] != "unmapped_or_uncertain" for row in review_rows),
        "without_retrieved_literature": sum(row["retrieved_pmid_count"] == 0 for row in review_rows),
        "patient_route_changes": len(changes),
        "changes_by_route": dict(sorted(Counter(f"{row['baseline_action']} -> {row['hypothetical_action']}" for row in changes).items())),
        "shadow_candidate_pool": shadow["counts"]["shadow_candidate_pool"],
        "baseline_candidate_pool": sum(row["screening_action"] in {"retain_baseline", "forward_priority_review", "forward_context_review"} for row in baseline),
        "output_dir": str(output_dir.resolve()),
    }
    (output_dir / "audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Review unapproved taxonomy drafts with a hypothetical patient routing shadow.")
    parser.add_argument("queue", type=Path)
    parser.add_argument("draft_dir", type=Path)
    parser.add_argument("literature_json", type=Path)
    parser.add_argument("inventory_dir", type=Path)
    parser.add_argument("baseline_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    print(json.dumps(audit(args.queue, args.draft_dir, args.literature_json, args.inventory_dir, args.baseline_dir, args.output_dir, args.policy), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

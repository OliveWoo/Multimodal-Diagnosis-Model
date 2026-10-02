from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from tools.build_independent_multi_assay_screening_shadow import (
    build_shadow,
    positive_selected_observations,
)


DEFAULT_POLICY = Path("rules/manual_style_multi_assay_entry_v2.json")


def route_candidate(candidate: dict[str, Any], taxonomy: dict[str, Any], policy: dict[str, Any]) -> dict[str, Any]:
    """Keep plausible positive candidates for scoring; tiers are not diagnoses."""
    selected = positive_selected_observations(candidate)
    selected_ids = {str(item.get("seq_id") or "") for item in selected}
    if "" in selected_ids:
        raise ValueError("Selected positive observation lacks seq_id")
    dna = any(item.get("nucleic_type") == "DNA" for item in selected)
    rna = any(item.get("nucleic_type") == "RNA" for item in selected)
    ranks = [
        int(item["rank_by_reads_in_retained_universe"])
        for item in selected
        if item.get("rank_by_reads_in_retained_universe") is not None
    ]
    best_rank = min(ranks) if ranks else None
    top_rank = best_rank is not None and best_rank <= int(policy["rank_in_retained_universe_at_most"])
    repeated_test = len(selected_ids) >= 2
    cross_molecule = dna and rna
    family = str(taxonomy.get("primary_rule_family") or "unmapped_or_uncertain")
    mapping_status = str(taxonomy.get("mapping_status") or "unmapped")
    biological_class = str(taxonomy.get("biological_class") or "unknown")
    rules: list[str] = []
    tier = ""

    if not selected:
        disposition = "qc_only"
        rules.append("NO_SELECTED_POSITIVE_TEST")
    elif mapping_status == "unmapped" or family == "unmapped_or_uncertain":
        disposition = "clinical_review"
        rules.append("TAXONOMY_REVIEW_REQUIRED")
    elif family in policy["excluded_background_families"]:
        disposition = "qc_only"
        rules.append("COMMENSAL_VIROME_QC")
    elif family in policy["priority_families"]:
        disposition = "scorer_entry"
        tier = "high" if top_rank and repeated_test else "medium" if top_rank or repeated_test else "low_colonizer"
        if rna and not dna and biological_class not in {"virus", "viral"} and tier == "high":
            tier = "medium"
            rules.append("RNA_ONLY_NONVIRAL_TIER_CAP")
        rules.append("PRIORITY_FAMILY_ENTRY")
    elif family in policy["conditional_families"] or family in policy["context_families"]:
        disposition = "scorer_entry"
        tier = "medium" if top_rank and repeated_test else "low_colonizer"
        rules.append("CONTEXT_FAMILY_ENTRY")
    elif family in policy["background_families"]:
        if top_rank or repeated_test:
            disposition = "scorer_entry"
            tier = "low_colonizer"
            rules.append("BACKGROUND_WITH_ANALYTICAL_SUPPORT")
        else:
            disposition = "qc_only"
            rules.append("WEAK_BACKGROUND_QC")
    else:
        disposition = "clinical_review"
        rules.append("UNRECOGNIZED_FAMILY_REVIEW")

    if tier:
        rules.append("ENTRY_TIER_" + tier.upper())
    if top_rank:
        rules.append("WITHIN_TEST_TOP_RANK")
    if repeated_test:
        rules.append("MULTIPLE_SELECTED_TESTS_SAME_CASE")
    if cross_molecule:
        rules.append("DNA_RNA_SELECTED_CONCORDANCE")
    if selected and any(item.get("test_qc_status") == "partial_evaluable" for item in selected):
        rules.append("PARTIAL_TEST_QC")
    return {
        "disposition": disposition,
        "screening_tier": tier,
        "rule_ids": rules,
        "selected_positive_test_count": len(selected_ids),
        "best_rank_in_retained_universe_selected": best_rank,
        "dna_selected": dna,
        "rna_selected": rna,
        "technical_repeat": repeated_test,
        "cross_molecule_selected": cross_molecule,
        "analytical_rank_scope": "within_test_retained_universe_only",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build manual-style multi-assay candidate-entry shadow")
    parser.add_argument("inventory_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    print(json.dumps(build_shadow(args.inventory_dir, args.output_dir, args.policy, router=route_candidate), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

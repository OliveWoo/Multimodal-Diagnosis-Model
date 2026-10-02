from __future__ import annotations

from tools.run_kh_ablation_d_history_only import bridge_candidate, review_tiers
from tools.run_kh_ablation_e_taxonomy_only import compare_rows


def test_history_bridge_preserves_old_taxonomy_and_single_test_axis() -> None:
    candidate = {
        "organism_name": "Bacteroides fragilis",
        "taxonomy_profile": {"primary_rule_family": "oral_aspiration_or_anaerobe"},
        "rank_priority": "3",
        "reads": 3551,
        "evidence_source": "mNGS_ranked",
        "specimen_class": "S2_lower_respiratory",
    }
    result = bridge_candidate(candidate, decision="review_context_needed", review_item=None)
    assert result["taxonomy_family"] == "oral_aspiration_or_anaerobe"
    assert result["analytical_profile"]["selected_positive_test_count"] == 1
    assert result["analytical_profile"]["reproducibility_axis"] is False
    assert result["analytical_profile"]["cross_molecule_selected"] is False


def test_review_tier_priority_prefers_high_over_lower_sections() -> None:
    review = {
        "review_low_specificity": [{"organism_name": "Example species"}],
        "review_context_needed": [{"organism_name": "Example species"}],
        "review_high_priority": [{"organism_name": "Example species"}],
    }
    tiers = review_tiers(review)
    assert tiers["examplespecies"][0] == "review_high_priority"


def test_taxonomy_comparison_does_not_invent_picked_change() -> None:
    control = [
        {
            "patient_id": 1,
            "organism_name": "Example species",
            "picked": False,
            "integrated_level": "Level 4",
            "taxonomy_family": "unmapped_or_uncertain",
            "taxonomy_mapping_status": "unmapped",
            "taxonomy_taxid": "",
            "taxonomy_display_name": "Example species",
            "formal_pick_exclusion_rule": "",
        }
    ]
    treatment = [
        {
            **control[0],
            "taxonomy_family": "oral_aspiration_or_anaerobe",
            "taxonomy_mapping_status": "exact_species",
        }
    ]
    changes = compare_rows(control, treatment, {1: ["Other species"]})
    assert len(changes) == 1
    assert changes[0]["change"] == "taxonomy_metadata_changed"
    assert changes[0]["before_picked"] is False
    assert changes[0]["after_picked"] is False

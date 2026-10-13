from tools.run_kh_ablation_d_history_only import bridge_candidate


def _candidate() -> dict:
    return {
        "organism_name": "Example organism",
        "reads": 120,
        "rank_priority": 3,
        "specimen_class": "S2_BALF",
        "taxonomy_profile": {
            "primary_rule_family": "unmapped_or_uncertain",
            "mapping_status": "unmapped",
        },
    }


def test_bridge_candidate_preserves_embedded_taxonomy_by_default() -> None:
    result = bridge_candidate(
        _candidate(), decision="review_context_needed", review_item=None
    )

    assert result["taxonomy_family"] == "unmapped_or_uncertain"
    assert result["ablation_bridge"]["taxonomy_mode"] == "embedded_old"
    assert result["ablation_bridge"]["old_taxonomy_profile_preserved"] is True


def test_bridge_candidate_can_receive_answer_blind_taxonomy_override() -> None:
    result = bridge_candidate(
        _candidate(),
        decision="review_context_needed",
        review_item=None,
        taxonomy_profile_override={
            "primary_rule_family": "oral_aspiration_or_anaerobe",
            "mapping_status": "exact",
            "taxid": 123,
        },
        taxonomy_mode="current_taxonomy_v2_6",
    )

    assert result["taxonomy_family"] == "oral_aspiration_or_anaerobe"
    assert result["taxonomy_profile"]["taxid"] == 123
    assert result["ablation_bridge"]["taxonomy_mode"] == "current_taxonomy_v2_6"
    assert result["ablation_bridge"]["old_taxonomy_profile_preserved"] is False

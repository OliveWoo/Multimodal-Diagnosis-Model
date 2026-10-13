from __future__ import annotations

from tools.build_two_hospital_legacy_ober_bridge import BRIDGE_SUFFIX, TIER_MAP


def test_bridge_has_fixed_answer_blind_tier_mapping() -> None:
    assert BRIDGE_SUFFIX == "legacy_external_ober_bridge_v1"
    assert TIER_MAP == {
        "Priority": "review_high_priority",
        "Context": "review_context_needed",
    }
    assert "Picked" not in TIER_MAP
    assert "Audit" not in TIER_MAP

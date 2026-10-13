from tools.apply_legacy_history_picked_promotion import promotion_reasons, routed_decision


POLICY = {
    "eligible_transition": {
        "from_route": "Context",
        "to_route": "Priority",
        "history_policy_family": "aspiration",
    },
    "eligible_taxonomy_families": ["oral_aspiration_or_anaerobe"],
    "required_history_support_level": "durable_family_specific_risk",
    "require_explicit_aspiration_event": True,
    "required_specimen_contexts": ["lower_respiratory"],
    "max_best_rank_in_retained_universe": 3,
    "min_selected_positive_test_count": 1,
    "require_history_route_changed": True,
    "block_when_history_adjustment_blockers_present": True,
    "block_when_direct_colonization_or_not_infection": True,
}


def candidate() -> dict:
    return {
        "analytic_route": "Context",
        "history_adjusted_route": "Priority",
        "history_policy_family": "aspiration",
        "taxonomy_family": "oral_aspiration_or_anaerobe",
        "specimen_context": "lower_respiratory",
        "history_route_changed": True,
        "analytical_profile": {
            "best_rank_in_retained_universe": 3,
            "selected_positive_test_count": 1,
        },
        "history_adjustment_blockers": [],
        "colonization_interpretation": {
            "has_explicit_colonization": False,
            "has_explicit_current_not_infection": False,
        },
    }


def patient_payload() -> dict:
    return {
        "aspiration_history_profile": {
            "support_level": "durable_family_specific_risk",
            "has_explicit_aspiration_event": True,
        }
    }


def test_full_aspiration_convergence_is_eligible_for_picked() -> None:
    eligible, _, blockers = promotion_reasons(candidate(), patient_payload(), POLICY)
    assert eligible is True
    assert blockers == []


def test_audit_to_context_is_not_promoted() -> None:
    value = candidate()
    value["analytic_route"] = "Audit"
    value["history_adjusted_route"] = "Context"
    eligible, _, blockers = promotion_reasons(value, patient_payload(), POLICY)
    assert eligible is False
    assert "eligible_from_route" in blockers


def test_explicit_colonization_blocks_promotion() -> None:
    value = candidate()
    value["colonization_interpretation"]["has_explicit_colonization"] = True
    eligible, _, blockers = promotion_reasons(value, patient_payload(), POLICY)
    assert eligible is False
    assert "no_direct_colonization_or_not_infection" in blockers


def test_downstream_gate_preserves_latest_history_routed_decision() -> None:
    value = {
        "decision": "review_context_needed",
        "history_routed_decision": "review_high_priority",
    }
    assert routed_decision(value) == "review_high_priority"


def test_downstream_gate_accepts_legacy_decision_without_history_field() -> None:
    assert routed_decision({"decision": "review_context_needed"}) == "review_context_needed"

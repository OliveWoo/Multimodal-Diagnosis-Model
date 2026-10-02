from __future__ import annotations

from tools.build_kh_lost_tp_divergence_audit import (
    explanation_class,
    first_divergence_stage,
    old_pick_basis,
    same_organism,
)


def test_old_pick_basis_reads_legacy_signal_and_support_modules() -> None:
    assert old_pick_basis(
        {
            "mngs_signal_tier": "M1_strong",
            "key_evidence": {"support_modules": []},
        }
    ) == "M1_signal_only"
    assert old_pick_basis(
        {
            "mngs_signal_tier": "M3_protected",
            "key_evidence": {"support_modules": ["culture"]},
        }
    ) == "M3_protected_plus_hospital"


def test_first_divergence_is_analytical_when_candidate_is_context() -> None:
    analytical = {
        "decision": "review_context_needed",
        "forward_to_clinical_scorer": True,
    }
    assert first_divergence_stage(analytical, {}) == "analytical_not_picked"


def test_low_specificity_without_independent_axis_is_explained() -> None:
    old = {
        "mngs_signal_tier": "M1_strong",
        "key_evidence": {"support_modules": []},
    }
    analytical = {
        "decision": "review_context_needed",
        "forward_to_clinical_scorer": True,
        "taxonomy_family": "environmental_low_specificity",
        "analytical_profile": {
            "cross_molecule_selected": False,
            "reproducibility_axis": True,
        },
        "hospital_profile": {"best_direct_level": None},
    }
    assert explanation_class(old, analytical, {}, None) == (
        "taxonomy_family_requires_independent_evidence"
    )


def test_known_aliases_match_for_audit_selection() -> None:
    assert same_organism("HSV-1", "Human alphaherpesvirus 1")

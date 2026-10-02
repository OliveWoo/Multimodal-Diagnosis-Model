from __future__ import annotations

import json
import unittest
from pathlib import Path

from tools.build_test_aware_possible_pathogen_shadow import (
    ASPIRATION_ROLE,
    CANDIDA_ROLE,
    FALLBACK_ROLE,
    LOW_SPECIFICITY_REPEAT_ROLE,
    PJP_ROLE,
    REACTIVATION_ROLE,
    STRONG_CROSS_MOLECULE_ROLE,
    annotate_candidate,
    apply_best_available_fallback,
    apply_primary_report_compaction,
)


POLICY = json.loads(
    Path("rules/test_aware_unified_reporting_v2.json").read_text(
        encoding="utf-8"
    )
)
ROUTE_A_POLICY = json.loads(
    Path("rules/test_aware_unified_reporting_v3_route_a.json").read_text(
        encoding="utf-8"
    )
)
SPECIAL_FAMILY_V4_POLICY = json.loads(
    Path(
        "rules/test_aware_unified_reporting_v4_special_family_possible_shadow.json"
    ).read_text(encoding="utf-8")
)


def candidate(
    organism: str,
    family: str,
    *,
    decision: str = "review_context_needed",
    rank: int = 1,
    count: int = 2,
    cross: bool = False,
    reproducible: bool = True,
    context: str = "lower_respiratory",
    host: bool = False,
    signals: list[dict] | None = None,
) -> dict:
    return {
        "patient_id": "1",
        "organism_name": organism,
        "clinical_decision": decision,
        "clinical_level": "Level 3" if decision == "review_high_priority" else "Level 4",
        "formal_pick_allowed": False,
        "selection_role": (
            "high_priority_review_candidate"
            if decision == "review_high_priority"
            else "context_review_candidate"
        ),
        "taxonomy_profile": {
            "primary_rule_family": family,
            "taxonomic_rank": "species",
            "mapping_status": "exact_species",
        },
        "source_candidate": {"specimen_context": context},
        "analytical_profile": {
            "per_test_signals": signals or [{
                "seq_id": "S1",
                "nucleic_type": "DNA",
                "reads": 100,
                "rank_in_retained_universe": rank,
                "analytical_rank_priority": 1,
                "rpm_total": 10,
                "qc_status": "evaluable",
            }],
        },
        "host_evidence": (
            [{"evidence_id": "H1", "assertion": "recorded_history"}]
            if host else []
        ),
        "promotion_gate": {
            "axes": {
                "taxonomy_family": family,
                "precise_identity": True,
                "specimen_context": context,
                "best_rank_in_retained_universe": rank,
                "selected_positive_test_count": count,
                "reproducibility_axis": reproducible,
                "cross_molecule_selected": cross,
                "direct_hospital_level": None,
                "host_support": host,
                "direct_evidence_timing_profile": {"rows": []},
            }
        },
    }


def test_pjp_recorded_host_context_creates_possible_not_picked() -> None:
    item = candidate(
        "Pneumocystis jirovecii",
        "high_consequence_opportunistic",
        decision="review_high_priority",
        rank=2,
        count=2,
        host=True,
    )
    result = annotate_candidate(item, POLICY)
    assert result["possible_reporting_role"] == PJP_ROLE
    assert result["final_reporting_tier"] == "Possible"
    assert result["clinical_decision"] == "review_high_priority"


def test_scorer_verified_route_a_creates_possible_not_picked() -> None:
    item = candidate(
        "Corynebacterium striatum",
        "skin_airway_colonizer_prone",
        decision="review_high_priority",
        rank=1,
        count=2,
        cross=False,
    )
    item["rule_ids"] = ["TA-S4-LOW-SPECIFICITY-REPEAT-HIGH"]
    item["source_candidate"] = {
        "specimen_context": "lower_respiratory",
        "rule_ids": ["TA-S4-LOW-SPECIFICITY-REPEAT-HIGH"],
    }
    item["analytical_profile"]["low_specificity_repeat_high_gate"] = {
        "eligible": True,
        "technical_repeat_group_type": (
            "same_run_parallel_12000g_and_unannotated"
        ),
        "minimum_selected_branch_rpm": 205.6,
    }
    result = annotate_candidate(item, ROUTE_A_POLICY)
    assert result["possible_reporting_role"] == LOW_SPECIFICITY_REPEAT_ROLE
    assert result["final_reporting_tier"] == "Possible"
    assert result["clinical_decision"] == "review_high_priority"
    assert not result["formal_pick_allowed"]


def test_reactivation_cross_molecule_signal_creates_possible_with_caution() -> None:
    item = candidate(
        "Human alphaherpesvirus 1",
        "herpes_or_reactivation_virus",
        rank=1,
        count=6,
        cross=True,
    )
    result = annotate_candidate(item, POLICY)
    assert result["possible_reporting_role"] == REACTIVATION_ROLE
    assert "not excluded" in " ".join(
        result["possible_pathogen_gate"]["cautions"]
    )


def test_explicit_current_reactivation_reduces_incidental_caution() -> None:
    item = candidate(
        "Human alphaherpesvirus 1",
        "herpes_or_reactivation_virus",
        rank=1,
        count=2,
        cross=True,
    )
    item["reactivation_interpretation"] = {
        "has_explicit_current_reactivation": True,
        "has_compatible_pulmonary_syndrome": True,
    }
    result = annotate_candidate(item, POLICY)
    assert "concern is reduced" in " ".join(
        result["possible_pathogen_gate"]["cautions"]
    )


def test_recorded_host_context_supports_exactly_one_top_rank_virus_test() -> None:
    item = candidate(
        "Human alphaherpesvirus 3",
        "herpes_or_reactivation_virus",
        rank=1,
        count=1,
        host=True,
    )
    result = annotate_candidate(item, POLICY)
    assert result["possible_reporting_role"] == REACTIVATION_ROLE
    assert result["possible_pathogen_gate"]["route"] == (
        "reactivation_family_top_rank_with_recorded_host_context"
    )


def test_recorded_host_context_does_not_mislabel_multiple_tests_as_single() -> None:
    item = candidate(
        "Human alphaherpesvirus 1",
        "herpes_or_reactivation_virus",
        rank=1,
        count=3,
        cross=False,
        host=True,
    )
    result = annotate_candidate(item, POLICY)
    assert result["possible_reporting_role"] != REACTIVATION_ROLE
    assert result["final_reporting_tier"] == "Context"


def test_v4_pjp_accepts_high_confidence_pre_pneumonia_phenotype_for_possible_only() -> None:
    item = candidate(
        "Pneumocystis jirovecii",
        "high_consequence_opportunistic",
        rank=2,
        count=3,
        cross=True,
    )
    item["phenotype_evidence"] = {
        "production_link_verified": False,
        "routed_phenotype_context": [{
            "source_id": "PHENO-1",
            "phenotype": "IMMUNOCOMPROMISED",
            "status": "YES",
            "confidence": "HIGH",
            "temporal_relation": "PRE_PNEUMONIA",
            "phenotype_qc_status": "REVIEW_REQUIRED",
        }],
    }
    result = annotate_candidate(item, SPECIAL_FAMILY_V4_POLICY)
    assert result["possible_reporting_role"] == PJP_ROLE
    assert result["final_reporting_tier"] == "Possible"
    assert not result["formal_pick_allowed"]


def test_v4_pjp_rejects_uncertain_phenotype_context() -> None:
    item = candidate(
        "Pneumocystis jirovecii",
        "high_consequence_opportunistic",
        rank=2,
        count=3,
        cross=True,
    )
    item["phenotype_evidence"] = {
        "production_link_verified": False,
        "routed_phenotype_context": [{
            "source_id": "PHENO-1",
            "phenotype": "IMMUNOCOMPROMISED",
            "status": "NO",
            "confidence": "LOW",
            "temporal_relation": "UNCLEAR",
            "phenotype_qc_status": "PASS",
        }],
    }
    result = annotate_candidate(item, SPECIAL_FAMILY_V4_POLICY)
    assert result["final_reporting_tier"] == "Context"


def test_v4_reactivation_requires_support_with_three_cross_molecule_tests() -> None:
    unsupported = candidate(
        "Human betaherpesvirus 5",
        "herpes_or_reactivation_virus",
        rank=2,
        count=3,
        cross=True,
    )
    supported = candidate(
        "Human betaherpesvirus 5",
        "herpes_or_reactivation_virus",
        rank=2,
        count=3,
        cross=True,
        host=True,
    )
    assert annotate_candidate(
        unsupported, SPECIAL_FAMILY_V4_POLICY
    )["final_reporting_tier"] == "Context"
    assert annotate_candidate(
        supported, SPECIAL_FAMILY_V4_POLICY
    )["possible_reporting_role"] == REACTIVATION_ROLE


def test_v4_disables_host_only_single_test_reactivation_route() -> None:
    item = candidate(
        "Human alphaherpesvirus 1",
        "herpes_or_reactivation_virus",
        rank=1,
        count=1,
        host=True,
    )
    result = annotate_candidate(item, SPECIAL_FAMILY_V4_POLICY)
    assert result["final_reporting_tier"] == "Context"


def test_v4_single_virus_test_requires_event_aligned_direct_support() -> None:
    item = candidate(
        "Human alphaherpesvirus 3",
        "herpes_or_reactivation_virus",
        rank=1,
        count=1,
        host=True,
        context="blood_or_systemic",
    )
    item["promotion_gate"]["axes"]["direct_evidence_timing_profile"] = {
        "rows": [{
            "observation_id": "PCR-1",
            "positive": True,
            "timing": "within_event_window",
            "specimen_type": "C.S.F.",
        }]
    }
    result = annotate_candidate(item, SPECIAL_FAMILY_V4_POLICY)
    assert result["possible_reporting_role"] == REACTIVATION_ROLE
    assert result["possible_pathogen_gate"]["route"] == (
        "reactivation_family_event_aligned_direct_supported_signal"
    )


def test_v4_candida_requires_clinical_convergence_not_mngs_alone() -> None:
    item = candidate(
        "Candida albicans",
        "candida_or_yeast",
        rank=1,
        count=3,
        cross=True,
    )
    assert annotate_candidate(
        item, SPECIAL_FAMILY_V4_POLICY
    )["final_reporting_tier"] == "Context"
    item["promotion_gate"]["axes"]["candida_invasive_evidence"] = {
        "rows": [{
            "assertion": "positive",
            "timing": "within_event_window",
            "identity_match": True,
            "invasive_specimen": False,
            "specimen_category": "Lower_Respiratory",
            "specimen_type": "Lower BAL",
            "observation_id": "CUL-1",
        }]
    }
    result = annotate_candidate(item, SPECIAL_FAMILY_V4_POLICY)
    assert result["possible_reporting_role"] == CANDIDA_ROLE
    assert result["final_reporting_tier"] == "Possible"


def test_v4_event_aligned_invasive_candida_is_possible() -> None:
    item = candidate(
        "Candida tropicalis",
        "candida_or_yeast",
        rank=8,
        count=1,
        cross=False,
    )
    item["promotion_gate"]["axes"]["candida_invasive_evidence"] = {
        "rows": [{
            "assertion": "positive",
            "timing": "within_event_window",
            "identity_match": True,
            "invasive_specimen": True,
            "specimen_category": "Sterile_Site",
            "specimen_type": "Blood",
            "observation_id": "CUL-2",
        }]
    }
    result = annotate_candidate(item, SPECIAL_FAMILY_V4_POLICY)
    assert result["possible_reporting_role"] == CANDIDA_ROLE


def test_rank_four_gi_source_cross_molecule_signal_creates_possible() -> None:
    item = candidate(
        "Enterococcus faecalis",
        "gi_urinary_or_nonpulmonary_prone",
        rank=4,
        count=6,
        cross=True,
    )
    result = annotate_candidate(item, POLICY)
    assert result["possible_reporting_role"] == STRONG_CROSS_MOLECULE_ROLE
    assert result["final_reporting_tier"] == "Possible"


def test_history_promoted_aspiration_candidate_with_complete_signal_is_possible() -> None:
    item = candidate(
        "Bacteroides fragilis",
        "oral_aspiration_or_anaerobe",
        decision="review_high_priority",
        rank=3,
        count=3,
        cross=True,
        signals=[
            {
                "seq_id": f"S{index}",
                "nucleic_type": "DNA" if index < 3 else "RNA",
                "reads": 3000,
                "rank_in_retained_universe": 3,
                "analytical_rank_priority": 3,
                "rpm_total": 80.0,
                "qc_status": "evaluable",
            }
            for index in range(1, 4)
        ],
    )
    item["history_policy_family"] = "aspiration"
    item["history_adjusted_route"] = "Priority"
    item["history_rule_ids"] = ["HIST-ASP3-CONTEXT-TO-PRIORITY"]
    item["history_evidence_ids"] = ["H1"]
    result = annotate_candidate(item, POLICY)
    assert result["possible_reporting_role"] == ASPIRATION_ROLE
    assert result["final_reporting_tier"] == "Possible"
    assert result["clinical_decision"] == "review_high_priority"


def test_aspiration_signal_without_history_priority_remains_context() -> None:
    item = candidate(
        "Fusobacterium nucleatum",
        "oral_aspiration_or_anaerobe",
        decision="review_context_needed",
        rank=2,
        count=3,
        cross=True,
        signals=[
            {
                "seq_id": f"S{index}",
                "nucleic_type": "DNA" if index < 3 else "RNA",
                "reads": 3000,
                "rank_in_retained_universe": 2,
                "analytical_rank_priority": 3,
                "rpm_total": 25.0,
                "qc_status": "evaluable",
            }
            for index in range(1, 4)
        ],
    )
    item["history_policy_family"] = "aspiration"
    item["history_adjusted_route"] = "Context"
    item["history_rule_ids"] = []
    result = annotate_candidate(item, POLICY)
    assert result["possible_reporting_role"] != ASPIRATION_ROLE
    assert result["final_reporting_tier"] == "Context"


def test_low_rpm_repeated_analytical_only_candidate_moves_to_context() -> None:
    item = candidate(
        "Pseudomonas aeruginosa",
        "hospital_or_nonfermenter_gnb",
        decision="review_high_priority",
        rank=3,
        count=3,
        cross=True,
        signals=[
            {
                "seq_id": f"S{index}",
                "nucleic_type": "DNA" if index < 3 else "RNA",
                "reads": 100,
                "rank_in_retained_universe": 3,
                "analytical_rank_priority": 1,
                "rpm_total": 2.7,
                "qc_status": "evaluable",
            }
            for index in range(1, 4)
        ],
    )
    annotated = annotate_candidate(item, POLICY)
    assert annotated["final_reporting_tier"] == "Possible"
    result = apply_primary_report_compaction([annotated], POLICY)[0]
    assert result["final_reporting_tier"] == "Context"
    assert result["possible_pathogen_gate"]["pre_compaction_decision"][
        "reporting_role"
    ] == "analytical_or_direct_possible_pathogen"


def test_incomplete_qc_low_rpm_analytical_candidate_is_not_compacted() -> None:
    item = candidate(
        "Stenotrophomonas maltophilia",
        "hospital_or_nonfermenter_gnb",
        decision="review_high_priority",
        rank=3,
        count=2,
        cross=False,
        signals=[
            {
                "seq_id": "S1",
                "nucleic_type": "RNA",
                "reads": 100,
                "rank_in_retained_universe": 3,
                "analytical_rank_priority": 1,
                "rpm_total": 2.0,
                "qc_status": "partial_evaluable",
            },
            {
                "seq_id": "S2",
                "nucleic_type": "RNA",
                "reads": 80,
                "rank_in_retained_universe": 3,
                "analytical_rank_priority": 1,
                "rpm_total": None,
                "qc_status": "partial_evaluable",
            },
        ],
    )
    result = apply_primary_report_compaction(
        [annotate_candidate(item, POLICY)], POLICY
    )[0]
    assert result["final_reporting_tier"] == "Possible"


def test_two_fully_evaluable_low_rpm_tests_are_compacted_without_count_penalty() -> None:
    item = candidate(
        "Pseudomonas aeruginosa",
        "hospital_or_nonfermenter_gnb",
        decision="review_high_priority",
        rank=3,
        count=2,
        cross=True,
        signals=[
            {
                "seq_id": f"S{index}",
                "nucleic_type": "DNA" if index == 1 else "RNA",
                "reads": 100,
                "rank_in_retained_universe": 3,
                "analytical_rank_priority": 1,
                "rpm_total": 2.0,
                "qc_status": "evaluable",
            }
            for index in range(1, 3)
        ],
    )
    result = apply_primary_report_compaction(
        [annotate_candidate(item, POLICY)], POLICY
    )[0]
    assert result["final_reporting_tier"] == "Context"


def test_inherited_taxonomy_without_host_moves_analytical_only_to_context() -> None:
    item = candidate(
        "Serratia nevei",
        "hospital_or_nonfermenter_gnb",
        decision="review_high_priority",
        rank=2,
        count=3,
        cross=True,
        signals=[{
            "seq_id": "S1",
            "nucleic_type": "DNA",
            "reads": 9000,
            "rank_in_retained_universe": 2,
            "analytical_rank_priority": 1,
            "rpm_total": 200.0,
            "qc_status": "evaluable",
        }],
    )
    item["taxonomy_profile"]["mapping_status"] = "genus_inherited"
    result = apply_primary_report_compaction(
        [annotate_candidate(item, POLICY)], POLICY
    )[0]
    assert result["final_reporting_tier"] == "Context"


def test_inherited_taxonomy_with_host_stays_possible() -> None:
    item = candidate(
        "Klebsiella pneumoniae",
        "typical_respiratory_pathogen",
        decision="review_high_priority",
        rank=3,
        count=3,
        cross=True,
        host=True,
        signals=[{
            "seq_id": "S1",
            "nucleic_type": "DNA",
            "reads": 15000,
            "rank_in_retained_universe": 3,
            "analytical_rank_priority": 1,
            "rpm_total": 500.0,
            "qc_status": "partial_evaluable",
        }],
    )
    item["taxonomy_profile"]["mapping_status"] = "genus_inherited"
    result = apply_primary_report_compaction(
        [annotate_candidate(item, POLICY)], POLICY
    )[0]
    assert result["final_reporting_tier"] == "Possible"


def test_fallback_prefers_best_fully_evaluable_test_and_keeps_low_confidence() -> None:
    acinetobacter = candidate(
        "Acinetobacter ursingii",
        "environmental_low_specificity",
        signals=[
            {
                "seq_id": "S1", "nucleic_type": "DNA", "reads": 614,
                "rank_in_retained_universe": 1, "analytical_rank_priority": 2,
                "rpm_total": 13.4, "qc_status": "evaluable",
            },
            {
                "seq_id": "S2", "nucleic_type": "DNA", "reads": 70,
                "rank_in_retained_universe": 3, "analytical_rank_priority": 2,
                "rpm_total": None, "qc_status": "partial_evaluable",
            },
        ],
    )
    hsv = candidate(
        "Human alphaherpesvirus 1",
        "herpes_or_reactivation_virus",
        host=False,
        signals=[
            {
                "seq_id": "S1", "nucleic_type": "DNA", "reads": 73,
                "rank_in_retained_universe": 2, "analytical_rank_priority": 1,
                "rpm_total": 1.6, "qc_status": "evaluable",
            },
            {
                "seq_id": "S2", "nucleic_type": "DNA", "reads": 517,
                "rank_in_retained_universe": 1, "analytical_rank_priority": 1,
                "rpm_total": None, "qc_status": "partial_evaluable",
            },
        ],
    )
    candidates = [
        annotate_candidate(acinetobacter, POLICY),
        annotate_candidate(hsv, POLICY),
    ]
    result = apply_best_available_fallback(candidates, POLICY)
    selected = [item for item in result if item["selected_for_complete_report"]]
    assert [item["organism_name"] for item in selected] == [
        "Acinetobacter ursingii"
    ]
    assert selected[0]["possible_reporting_role"] == FALLBACK_ROLE
    assert selected[0]["final_reporting_tier"] == "Fallback-Possible"
    assert selected[0]["final_reporting_confidence"] == "low"


def load_tests(
    loader: unittest.TestLoader,
    tests: unittest.TestSuite,
    pattern: str | None,
) -> unittest.TestSuite:
    del loader, tests, pattern
    suite = unittest.TestSuite()
    for name, value in sorted(globals().items()):
        if name.startswith("test_") and callable(value):
            suite.addTest(unittest.FunctionTestCase(value, description=name))
    return suite

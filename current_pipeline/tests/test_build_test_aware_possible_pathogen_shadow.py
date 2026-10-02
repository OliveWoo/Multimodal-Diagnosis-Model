import copy
import json
import unittest
from pathlib import Path

from tools.build_test_aware_possible_pathogen_shadow import (
    ANALYTICAL_DIRECT_ROLE,
    DEFAULT_POLICY,
    LOW_SPECIFICITY_ROLE,
    NOT_SELECTED_ROLE,
    STRICT_ROLE,
    SYSTEMIC_ROLE,
    annotate_candidate,
)


def candidate(
    *, family="hospital_or_nonfermenter_gnb", decision="review_high_priority",
    rank=2, count=2, reproducible=True, cross=False, direct_level=None,
    direct_rows=None, context="lower_respiratory",
):
    role = (
        "existing_picked_pathogen" if decision == "picked_shadow"
        else "high_priority_review_candidate" if decision == "review_high_priority"
        else "context_review_candidate"
    )
    return {
        "patient_id": "1",
        "organism_name": "Pseudomonas aeruginosa",
        "clinical_decision": decision,
        "clinical_level": "Level 2" if decision == "picked_shadow" else "Level 3",
        "formal_pick_allowed": decision == "picked_shadow",
        "selection_role": role,
        "taxonomy_profile": {
            "primary_rule_family": family,
            "taxonomic_rank": "species",
            "mapping_status": "exact_species",
        },
        "source_candidate": {"specimen_context": context},
        "promotion_gate": {
            "axes": {
                "taxonomy_family": family,
                "precise_identity": True,
                "specimen_context": context,
                "best_rank_in_retained_universe": rank,
                "selected_positive_test_count": count,
                "reproducibility_axis": reproducible,
                "cross_molecule_selected": cross,
                "direct_hospital_level": direct_level,
                "direct_evidence_timing_profile": {
                    "rows": direct_rows or [],
                },
            }
        },
    }


class TestPossiblePathogenShadow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.policy = json.loads(Path(DEFAULT_POLICY).read_text(encoding="utf-8"))

    def test_strict_pick_is_preserved_without_mutating_clinical_fields(self):
        source = candidate(decision="picked_shadow")
        result = annotate_candidate(source, self.policy)
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual("existing_picked_pathogen", result["selection_role"])
        self.assertEqual(STRICT_ROLE, result["possible_reporting_role"])
        self.assertTrue(result["selected_for_complete_report"])

    def test_reproducible_top_three_typical_bacterium_is_possible(self):
        result = annotate_candidate(candidate(), self.policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])
        self.assertEqual(ANALYTICAL_DIRECT_ROLE, result["possible_reporting_role"])
        self.assertFalse(result["formal_pick_allowed"])

    def test_event_aligned_direct_positive_without_level_is_possible(self):
        rows = [{
            "observation_id": "ASSAY-1",
            "positive": True,
            "timing": "within_event_window",
            "specimen_type": "BAL",
        }]
        result = annotate_candidate(
            candidate(rank=8, count=3, direct_rows=rows), self.policy
        )
        self.assertEqual(ANALYTICAL_DIRECT_ROLE, result["possible_reporting_role"])
        self.assertIn("generic hospital level is unavailable", " ".join(
            result["possible_pathogen_gate"]["cautions"]
        ))

    def test_level_three_direct_culture_is_possible(self):
        rows = [{
            "observation_id": "CUL-1",
            "positive": True,
            "timing": "within_event_window",
            "specimen_type": "Lower BAL",
        }]
        result = annotate_candidate(
            candidate(
                family="environmental_low_specificity", rank=1, count=3,
                cross=True, direct_level=3, direct_rows=rows,
            ),
            self.policy,
        )
        self.assertEqual(ANALYTICAL_DIRECT_ROLE, result["possible_reporting_role"])

    def test_context_top_rank_cross_molecule_low_specificity_is_possible(self):
        result = annotate_candidate(
            candidate(
                family="skin_airway_colonizer_prone",
                decision="review_context_needed", rank=1, count=3, cross=True,
            ),
            self.policy,
        )
        self.assertEqual(LOW_SPECIFICITY_ROLE, result["possible_reporting_role"])
        self.assertEqual("review_context_needed", result["clinical_decision"])

    def test_high_top_rank_cross_molecule_low_specificity_stays_possible(self):
        result = annotate_candidate(
            candidate(
                family="skin_airway_colonizer_prone",
                decision="review_high_priority", rank=1, count=3, cross=True,
            ),
            self.policy,
        )
        self.assertEqual(LOW_SPECIFICITY_ROLE, result["possible_reporting_role"])
        self.assertEqual("review_high_priority", result["clinical_decision"])

    def test_insufficient_context_signal_stays_unselected(self):
        result = annotate_candidate(
            candidate(
                family="skin_airway_colonizer_prone",
                decision="review_context_needed", rank=2, count=2, cross=False,
            ),
            self.policy,
        )
        self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_special_families_cannot_use_generic_route(self):
        for family in self.policy["special_family_guardrails"]:
            with self.subTest(family=family):
                item = candidate(family=family, rank=1, count=6, cross=True)
                result = annotate_candidate(item, self.policy)
                self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_repeated_blood_observations_trigger_systemic_review_with_caution(self):
        rows = [
            {
                "observation_id": "BC-1", "positive": True,
                "timing": "within_event_window", "specimen_type": "Blood",
            },
            {
                "observation_id": "BC-2", "positive": True,
                "timing": "within_event_window", "specimen_type": "Blood",
            },
        ]
        source = candidate(
            family="skin_airway_colonizer_prone",
            decision="review_context_needed", rank=1, count=3, cross=True,
            direct_rows=rows,
        )
        snapshot = copy.deepcopy(source)
        result = annotate_candidate(source, self.policy)
        self.assertEqual(SYSTEMIC_ROLE, result["possible_reporting_role"])
        self.assertFalse(
            result["possible_pathogen_gate"]["axes"]["repeated_sterile_independence_verified"]
        )
        self.assertEqual(snapshot, source)


if __name__ == "__main__":
    unittest.main()

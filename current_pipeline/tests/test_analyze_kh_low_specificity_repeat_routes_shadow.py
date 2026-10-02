from __future__ import annotations

import unittest

from tools.analyze_kh_low_specificity_repeat_routes_shadow import (
    _apply_route_a_possible,
    classify_routes,
)


def base(*, ranks: str, rpms: str) -> dict[str, str]:
    return {
        "rank_only_strength_eligible": "True",
        "branch_ranks": ranks,
        "branch_rpms": rpms,
        "repeat_group_type": "same_run_parallel_12000g_and_unannotated",
    }


class LowSpecificityRepeatRoutesShadowTest(unittest.TestCase):
    def test_compelling_exact_species_repeat_enters_route_a(self) -> None:
        result = classify_routes(
            base(ranks="1|1", rpms="1072.5|205.6"),
            {
                "taxonomy_family": "skin_airway_colonizer_prone",
                "taxonomy_mapping_status": "exact_species",
                "final_reporting_tier": "Context",
            },
        )
        self.assertTrue(result["route_a_scorer_high_eligible"])
        self.assertFalse(result["route_b_fallback_review_eligible"])

    def test_best_available_incomplete_normalization_is_route_b_not_high(self) -> None:
        result = classify_routes(
            base(ranks="1|3", rpms="13.48|"),
            {
                "taxonomy_family": "environmental_low_specificity",
                "taxonomy_mapping_status": "genus_inherited",
                "final_reporting_tier": "Fallback-Possible",
                "reporting_route": "no_picked_or_possible_best_available_fallback",
            },
        )
        self.assertFalse(result["route_a_scorer_high_eligible"])
        self.assertTrue(result["route_b_fallback_review_eligible"])

    def test_low_rpm_environmental_repeat_is_not_high(self) -> None:
        result = classify_routes(
            base(ranks="1|1", rpms="1.98|1.45"),
            {
                "taxonomy_family": "environmental_low_specificity",
                "taxonomy_mapping_status": "exact_species",
                "final_reporting_tier": "Context",
            },
        )
        self.assertFalse(result["route_a_scorer_high_eligible"])
        self.assertFalse(result["route_b_fallback_review_eligible"])
        self.assertEqual(
            "normalized_burden_below_shadow_threshold", result["route_reason"]
        )

    def test_possible_suppresses_same_patient_fallback(self) -> None:
        rows = [
            {
                "patient_id": "37",
                "organism_name": "Corynebacterium striatum",
                "selected_for_complete_report": "False",
                "final_reporting_tier": "Context",
            },
            {
                "patient_id": "37",
                "organism_name": "Candida albicans",
                "selected_for_complete_report": "True",
                "final_reporting_tier": "Fallback-Possible",
            },
        ]
        updated = _apply_route_a_possible(
            rows, {("37", "Corynebacterium striatum")}
        )
        self.assertEqual("Possible", updated[0]["final_reporting_tier"])
        self.assertEqual("False", updated[1]["selected_for_complete_report"])
        self.assertEqual("Context", updated[1]["final_reporting_tier"])


if __name__ == "__main__":
    unittest.main()

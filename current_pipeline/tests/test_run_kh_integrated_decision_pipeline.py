import unittest

from tools.run_kh_integrated_decision_pipeline import (
    _compact_candidate,
    compare_candidate_views,
)


class IntegratedDecisionPipelineTest(unittest.TestCase):
    def test_compare_candidate_views_passes_for_identical_semantics(self):
        current = [{
            "patient_id": "1",
            "organism_name": "Example bacterium",
            "candidate_source": "mngs_positive_reads",
            "clinical_decision": "picked_shadow",
        }]
        result = compare_candidate_views(
            current,
            [dict(current[0])],
            fields=("clinical_decision",),
            include_source=True,
            stage="promotion",
        )
        self.assertTrue(result["passed"])
        self.assertEqual(result["differences"], [])

    def test_compare_candidate_views_reports_changed_field(self):
        current = [{
            "patient_id": "1",
            "organism_name": "Example bacterium",
            "candidate_source": "mngs_positive_reads",
            "clinical_decision": "review_high_priority",
        }]
        reference = [{**current[0], "clinical_decision": "picked_shadow"}]
        result = compare_candidate_views(
            current,
            reference,
            fields=("clinical_decision",),
            include_source=True,
            stage="promotion",
        )
        self.assertFalse(result["passed"])
        self.assertEqual(result["field_difference_count"], 1)
        self.assertEqual(result["differences"][0]["field"], "clinical_decision")

    def test_compact_candidate_keeps_route_and_audit_without_changing_tier(self):
        candidate = {
            "patient_id": "9",
            "organism_name": "Pneumocystis jirovecii",
            "candidate_source": "mngs_positive_reads",
            "analytic_route": "Context",
            "history_adjusted_route": "Priority",
            "pre_promotion_decision": "review_high_priority",
            "clinical_decision": "review_high_priority",
            "final_reporting_tier": "Possible",
            "selected_for_complete_report": True,
            "possible_reporting_role": "pneumocystis_possible_pathogen",
            "history_rule_ids": ["HIST-OPP4-EXAMPLE"],
            "rule_ids": ["CLIN-S1-EXAMPLE"],
            "promotion_gate": {
                "outcome": "retained_high_priority",
                "route": "no_promotion",
                "blockers": ["requires direct evidence for Picked"],
                "axes": {"specimen_context": "lower_respiratory"},
            },
            "possible_pathogen_gate": {
                "outcome": "included_as_possible_pathogen",
                "route": "pneumocystis_host_supported",
                "blockers": [],
                "cautions": ["direct confirmation unavailable"],
                "axes": {},
            },
        }
        result = _compact_candidate(candidate)
        self.assertEqual(result["decision_path"], {
            "analytical_route": "Context",
            "history_adjusted_route": "Priority",
            "pre_promotion_clinical_decision": "review_high_priority",
            "post_promotion_clinical_decision": "review_high_priority",
            "final_reporting_tier": "Possible",
        })
        self.assertTrue(result["reporting"]["selected_for_complete_report"])
        self.assertEqual(
            result["audit"]["history_rule_ids"], ["HIST-OPP4-EXAMPLE"]
        )


if __name__ == "__main__":
    unittest.main()

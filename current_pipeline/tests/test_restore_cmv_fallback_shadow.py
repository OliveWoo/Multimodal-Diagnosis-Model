# All fixtures below are synthetic, not patient exports.
import copy
import json
import unittest

from tools.restore_cmv_fallback_shadow import DEFAULT_POLICY, finalize_patient


class PreserveCmvFallbackTests(unittest.TestCase):
    def setUp(self):
        self.policy = json.loads(DEFAULT_POLICY.read_text(encoding="utf-8"))
        self.cmv = {"organism_name": "Human betaherpesvirus 5", "reporting": {
            "tier": "Fallback-Possible", "selected_for_complete_report": True}}
        self.patient = {"patient_id": "999", "all_candidates": [self.cmv]}
        self.row = {"patient_id": "999", "organism_name": "Human betaherpesvirus 5", "decision": {
            "shadow_reporting_tier": "Context", "selected_as_pulmonary_candidate": False,
            "axes": {"detection": "D1", "activity": "A0", "organ_support": "O0"}}}

    def test_existing_fallback_restored_and_counts_positive_without_changing_axes(self):
        original = copy.deepcopy((self.patient, self.row))
        rows, report = finalize_patient(self.patient, [self.row], self.policy)
        d = rows[0]["decision"]
        self.assertEqual(d["shadow_reporting_tier"], "Fallback-Possible")
        self.assertTrue(d["included_in_positive_prediction_endpoint"])
        self.assertEqual(d["axes"], self.row["decision"]["axes"])
        self.assertFalse(d["fallback_means_clinical_disease_confirmed"])
        self.assertTrue(report["has_at_least_one_reported_candidate"])
        self.assertEqual((self.patient, self.row), original)

    def test_nonviral_picked_blocks_restoring_viral_fallback(self):
        self.cmv["reporting"].update(tier="Context", selected_for_complete_report=False)
        self.patient["all_candidates"].append({"organism_name": "Escherichia coli", "reporting": {
            "tier": "Picked", "selected_for_complete_report": True}})
        rows, report = finalize_patient(self.patient, [self.row], self.policy)
        self.assertFalse(rows[0]["decision"]["fallback_preserved"])
        self.assertEqual(len(report["final_report_after_shadow"]), 1)

    def test_conflicting_existing_fallback_and_ordinary_selection_fail_closed(self):
        self.patient["all_candidates"].append({"organism_name": "Escherichia coli", "reporting": {
            "tier": "Picked", "selected_for_complete_report": True}})
        # The CMV is removed by the base shadow, so the final report remains valid;
        # it must not be restored just because the viral subset has no ordinary pick.
        rows, report = finalize_patient(self.patient, [self.row], self.policy)
        self.assertFalse(rows[0]["decision"]["fallback_preserved"])
        self.assertTrue(report["has_picked_or_ordinary_possible"])

    def test_new_ordinary_cmv_possible_is_not_downgraded_to_fallback(self):
        self.row["decision"].update(shadow_reporting_tier="Possible", selected_as_pulmonary_candidate=True)
        rows, report = finalize_patient(self.patient, [self.row], self.policy)
        self.assertEqual(rows[0]["decision"]["shadow_reporting_tier"], "Possible")
        self.assertFalse(rows[0]["decision"]["fallback_preserved"])

    def test_no_prior_fallback_is_not_invented(self):
        self.cmv["reporting"].update(tier="Context", selected_for_complete_report=False)
        rows, report = finalize_patient(self.patient, [self.row], self.policy)
        self.assertFalse(rows[0]["decision"]["included_in_positive_prediction_endpoint"])
        self.assertFalse(report["has_at_least_one_reported_candidate"])

    def test_failed_detection_or_explicit_shedding_not_blindly_restored(self):
        for change in [{"detection": "D0"}, {"activity": "A3"}]:
            self.row["decision"]["axes"] = {"detection": "D1", "activity": "A0", "organ_support": "O0", **change}
            rows, _ = finalize_patient(self.patient, [self.row], self.policy)
            self.assertFalse(rows[0]["decision"]["fallback_preserved"])

    def test_other_pathogen_fallback_is_preserved(self):
        self.cmv["reporting"].update(tier="Context", selected_for_complete_report=False)
        self.patient["all_candidates"].append({"organism_name": "Pneumocystis jirovecii", "reporting": {
            "tier": "Fallback-Possible", "selected_for_complete_report": True}})
        rows, report = finalize_patient(self.patient, [self.row], self.policy)
        self.assertFalse(rows[0]["decision"]["fallback_preserved"])
        self.assertEqual(report["final_report_after_shadow"][0]["organism_name"], "Pneumocystis jirovecii")

    def test_wrong_patient_is_rejected(self):
        self.row["patient_id"] = "different"
        with self.assertRaises(ValueError):
            finalize_patient(self.patient, [self.row], self.policy)


if __name__ == "__main__":
    unittest.main()

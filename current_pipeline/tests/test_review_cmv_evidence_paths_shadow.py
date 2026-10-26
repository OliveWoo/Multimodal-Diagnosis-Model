# All fixtures below are synthetic, not patient exports.
import copy
import json
import unittest
from pathlib import Path

from tools.review_cmv_evidence_paths_shadow import (
    DEFAULT_POLICY, apply_clarifications, counterfactual_independence, decide,
)
from tools.viral_evidence_axes_shadow import DEFAULT_POLICY as AXES_POLICY


class CmvReportingPathTests(unittest.TestCase):
    def setUp(self):
        self.policy = json.loads(DEFAULT_POLICY.read_text(encoding="utf-8"))
        self.axes = json.loads(AXES_POLICY.read_text(encoding="utf-8"))
        self.candidate = {"organism_name": "Human betaherpesvirus 5",
            "taxonomy_profile": {"mapping_status": "exact_species"},
            "analytical_profile": {"per_test_signals": [
                {"reads": 10, "rpm_total": 1, "rank_in_retained_universe": 4,
                 "nucleic_type": "DNA", "qc_status": "evaluable"}]}}
        self.bound = {"cohort": "SYNTHETIC", "patient_id": "999", "organism_key": "humanbetaherpesvirus5",
                      "episode_id": "X:E1", "claims": []}

    def claim(self, **changes):
        return {"evidence_id": "resp-assay-1", "patient_id": "999", "organism_key": "humanbetaherpesvirus5",
            "episode_id": "X:E1", "specimen": "BALF", "collected_time": "2024-01-10",
            "source_file": "source.json", "source_row": "[0]", "source_sha256": "hash",
            "review_complete": True, "event_aligned": True, "claim_type": "detection",
            "identity_match": "exact_or_alias", "assay_method": "PCR",
            "independence_class": "independent_clinical_assay", "biological_evidence_group": "resp-group",
            **changes}

    def evaluate(self, claims=(), candidate=None):
        bound = copy.deepcopy(self.bound)
        bound["claims"] = list(claims)
        return decide(bound, candidate or self.candidate, self.policy, self.axes)

    def test_independent_respiratory_detection_supports_low_confidence_possible(self):
        result = self.evaluate([self.claim()])
        self.assertEqual(result["shadow_reporting_tier"], "Possible")
        self.assertEqual((result["axes"]["activity"], result["axes"]["organ_support"]), ("A0", "O0"))
        self.assertFalse(result["reporting_selection_is_clinical_diagnosis"])

    def test_single_dna_branch_is_not_blocked(self):
        self.assertEqual(self.evaluate([self.claim()])["shadow_reporting_tier"], "Possible")

    def test_rank_reads_rpm_and_technical_branch_count_do_not_change_disease_route(self):
        for rank, reads, rpm, count in [(1, 1, .001, 1), (4, 100, 10, 2), (99, 1e6, 1e5, 12)]:
            candidate = copy.deepcopy(self.candidate)
            signal = candidate["analytical_profile"]["per_test_signals"][0]
            signal.update(reads=reads, rpm_total=rpm, rank_in_retained_universe=rank)
            candidate["analytical_profile"]["per_test_signals"] = [copy.deepcopy(signal) for _ in range(count)]
            self.assertEqual(self.evaluate([self.claim()], candidate)["shadow_reporting_tier"], "Possible")

    def test_mngs_only_large_signal_stays_visible_unknown(self):
        candidate = copy.deepcopy(self.candidate)
        candidate["analytical_profile"]["per_test_signals"] *= 20
        result = self.evaluate(candidate=candidate)
        self.assertEqual(result["shadow_reporting_tier"], "Context")
        self.assertTrue(result["doctor_visible"])

    def test_no_credible_mngs_does_not_use_corroboration_route(self):
        for change in [{"reads": 0}, {"qc_status": "failed"}]:
            candidate = copy.deepcopy(self.candidate)
            candidate["analytical_profile"]["per_test_signals"][0].update(change)
            self.assertEqual(self.evaluate([self.claim()], candidate)["shadow_reporting_tier"], "Context")

    def test_hospital_method_unknown_is_not_independent(self):
        result = self.evaluate([self.claim(assay_method=None, independence_class="not_established")])
        self.assertEqual(result["shadow_reporting_tier"], "Context")
        self.assertIn("respiratory_hospital_assay_independence_unknown", result["gaps"])

    def test_blood_only_does_not_create_pulmonary_possible(self):
        self.assertEqual(self.evaluate([self.claim(specimen="Blood")])["shadow_reporting_tier"], "Context")

    def test_different_biological_respiratory_specimen_can_corroborate(self):
        claim = self.claim(specimen="Sputum (Endo)", assay_method=None,
                           independence_class="independent_biological_specimen")
        self.assertEqual(self.evaluate([claim])["shadow_reporting_tier"], "Possible")

    def test_group_identity_is_not_exact_cmv_support(self):
        self.assertEqual(self.evaluate([self.claim(identity_match="group_only")])["shadow_reporting_tier"], "Context")

    def test_wrong_patient_episode_organism_or_unverified_is_not_used(self):
        for changes in [{"patient_id": "other"}, {"episode_id": "X:E2"},
                        {"organism_key": "humanalphaherpesvirus1"}, {"review_complete": False},
                        {"event_aligned": False}, {"collected_time": None}]:
            self.assertEqual(self.evaluate([self.claim(**changes)])["shadow_reporting_tier"], "Context")

    def test_generic_pneumonia_and_host_risk_do_not_create_attribution(self):
        claim = self.claim(claim_type="compatible_pulmonary_syndrome", specimen="clinical_context")
        self.candidate["host_evidence"] = [{"label": "immunocompromised"}]
        self.assertEqual(self.evaluate([claim])["shadow_reporting_tier"], "Context")

    def test_other_pathogen_does_not_veto_supported_cmv(self):
        self.candidate["other_selected_pathogens"] = ["Strong bacterial pathogen"]
        result = self.evaluate([self.claim()])
        self.assertEqual(result["shadow_reporting_tier"], "Possible")
        self.assertFalse(result["co_pathogen_veto_applied"])

    def test_numeric_value_does_not_create_activity_or_picked(self):
        result = self.evaluate([self.claim(numeric_value=9999999, unit="IU/mL")])
        self.assertEqual((result["axes"]["activity"], result["shadow_reporting_tier"]), ("A0", "Possible"))

    def test_reviewed_lung_attribution_route_remains_available(self):
        for kind, basis, tier in [("reviewed_probable_pulmonary_disease", "lower_respiratory_assay_interpreted", "Possible"),
                                  ("reviewed_proven_pulmonary_disease", "lung_specific_pathology", "Picked")]:
            claim = self.claim(claim_type=kind, organ="lung", attribution_basis=basis,
                               clinical_interpretation="Source-reviewed lung-specific interpretation")
            self.assertEqual(self.evaluate([claim])["shadow_reporting_tier"], tier)

    def test_explicit_shedding_requires_review_not_automatic_selection(self):
        claim = self.claim(claim_type="resolved_or_asymptomatic_shedding", evidence_id="review-1")
        self.assertEqual(self.evaluate([self.claim(), claim])["shadow_reporting_tier"], "Context")

    def test_fallback_and_reference_answer_are_not_decision_inputs(self):
        self.bound.update(current_tier="Fallback-Possible", answer_contains_cmv=True)
        result = decide(self.bound, self.candidate, self.policy, self.axes)
        self.assertEqual(result["shadow_reporting_tier"], "Context")
        self.assertFalse(result["benchmark_answers_used"])

    def note(self, **changes):
        return {"cohort": "SYNTHETIC", "patient_id": "999", "organism_key": "humanbetaherpesvirus5",
            "episode_id": "X:E1", "evidence_id": "resp-assay-1", "original_source_sha256": "hash",
            "original_source_row": "[0]", "assertion": "assay_independent_of_mngs", "clarification_id": "note-1",
            **changes}

    def test_user_independence_clarification_preserves_unknown_method_and_original(self):
        self.bound["claims"] = [self.claim(assay_method=None, independence_class="not_established", numeric_value=777, unit=None)]
        original = copy.deepcopy(self.bound)
        enriched, notes = apply_clarifications(self.bound, [self.note()])
        self.assertEqual(notes, ["note-1"])
        self.assertEqual(self.bound, original)
        self.assertIsNone(enriched["claims"][0]["assay_method"])
        self.assertIsNone(enriched["claims"][0]["unit"])
        self.assertEqual(decide(enriched, self.candidate, self.policy, self.axes)["shadow_reporting_tier"], "Possible")

    def test_clarification_does_not_spill_into_another_patient(self):
        self.bound["claims"] = [self.claim()]
        enriched, notes = apply_clarifications(self.bound, [self.note(patient_id="different")])
        self.assertEqual(enriched, self.bound)
        self.assertEqual(notes, [])

    def test_clarification_source_mismatch_fails_closed(self):
        self.bound["claims"] = [self.claim()]
        with self.assertRaises(ValueError):
            apply_clarifications(self.bound, [self.note(original_source_sha256="different")])

    def test_counterfactual_is_explicit_and_not_applied(self):
        self.bound["claims"] = [self.claim(assay_method=None, independence_class="not_established")]
        original = copy.deepcopy(self.bound)
        result = counterfactual_independence(self.bound, self.candidate, self.policy, self.axes)
        self.assertTrue(result[0]["counterfactual_only_not_applied"])
        self.assertEqual(result[0]["shadow_tier_if_confirmed"], "Possible")
        self.assertEqual(self.bound, original)

    def test_policy_has_no_patient_specific_gate(self):
        text = DEFAULT_POLICY.read_text(encoding="utf-8")
        self.assertNotIn('"patient_id"', text)
        self.assertNotIn('"maximum_rank"', text)
        self.assertNotIn('"minimum_positive_tests"', text)


if __name__ == "__main__":
    unittest.main()

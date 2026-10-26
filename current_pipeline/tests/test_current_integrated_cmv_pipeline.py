# All fixtures below are synthetic, not patient exports.
import copy
import json
import unittest

from tools import current_pipeline_defaults as defaults
from tools import run_kh_integrated_decision_pipeline as decisions
from tools import run_kh_integrated_release_v2 as release
from tools import run_external_integrated_release_v2 as external
from tools import build_test_aware_possible_pathogen_shadow as reporting
from tools.run_current_integrated_pipeline import apply_patient_policy, remap_clarifications


class CurrentIntegratedCmvTests(unittest.TestCase):
    def setUp(self):
        self.policy = json.loads(defaults.CMV_POLICY.read_text(encoding="utf-8"))
        from tools.viral_evidence_axes_shadow import DEFAULT_POLICY
        self.axes = json.loads(DEFAULT_POLICY.read_text(encoding="utf-8"))
        self.base, _ = reporting.load_policy(defaults.REPORTING_POLICY)
        self.candidate = {"patient_id": "999", "organism_name": "Human betaherpesvirus 5",
            "clinical_decision": "review_context_needed", "clinical_level": "Level 4", "formal_pick_allowed": False,
            "candidate_source": "mngs_positive_reads",
            "taxonomy_profile": {"mapping_status": "exact_species", "primary_rule_family": "herpes_or_reactivation_virus"},
            "analytical_profile": {"per_test_signals": [{"reads": 10, "qc_status": "evaluable",
                "rank_in_retained_universe": 1, "analytical_rank_priority": 1, "nucleic_type": "DNA", "rpm_total": 1}]},
            "possible_reporting_role": reporting.NOT_SELECTED_ROLE, "selected_for_complete_report": False,
            "final_reporting_tier": "Context", "possible_pathogen_gate": {
                "outcome": "not_selected", "route": "family_guardrail", "reporting_role": reporting.NOT_SELECTED_ROLE,
                "selected_for_complete_report": False, "axes": {"precise_identity": True,
                    "taxonomy_family": "herpes_or_reactivation_virus", "specimen_context": "blood_or_systemic",
                    "selected_positive_test_count": 1, "best_rank_in_retained_universe": 1}, "blockers": [], "cautions": []}}
        self.bound = {"cohort": "TEST", "patient_id": "999", "organism_key": "humanbetaherpesvirus5",
                      "episode_id": "X:E1", "claims": []}
        self.packet = {"patient_id": "999", "policy_id": "v7", "all_forwarded_candidates": [self.candidate],
                       "decision_evidence_profile": {"decision_evidence_mode": "mngs_only"}}

    def claim(self, **extra):
        return {"evidence_id": "new-path-id", "patient_id": "999", "organism_key": "humanbetaherpesvirus5",
            "episode_id": "X:E1", "claim_type": "detection", "specimen": "BALF", "collected_time": "2024-01-01",
            "source_file": "source.json", "source_row": "[0]", "source_sha256": "hash", "review_complete": True,
            "event_aligned": True, "identity_match": "exact_or_alias", "assay_method": None,
            "independence_class": "independent_clinical_assay", **extra}

    def apply(self):
        return apply_patient_policy(self.packet, {"humanbetaherpesvirus5": self.bound}, [], self.policy, self.axes, self.base)

    def test_public_lower_level_defaults_are_unified(self):
        self.assertEqual(release.ANALYTICAL_POLICY, defaults.ANALYTICAL_POLICY)
        self.assertEqual(release.PROMOTION_POLICY, decisions.DEFAULT_PROMOTION_POLICY)
        self.assertEqual(release.REPORTING_POLICY, decisions.DEFAULT_REPORTING_POLICY)
        self.assertEqual(external.run.__kwdefaults__["reporting_policy"], defaults.REPORTING_POLICY)
        self.assertEqual(decisions.promotion_current.__name__, "tools.promote_test_aware_clinical_high_shadow_v9")

    def test_source_id_remapped_by_content_row_and_episode_not_patient_exception(self):
        self.bound["claims"] = [self.claim()]
        note = {"clarification_id": "N1", "cohort": "TEST", "patient_id": "999", "organism_key": "humanbetaherpesvirus5",
            "episode_id": "X:E1", "evidence_id": "old-path-id", "original_source_sha256": "hash", "original_source_row": "[0]"}
        result = remap_clarifications([self.bound], [note])
        self.assertEqual(result[0]["evidence_id"], "new-path-id")
        self.assertEqual(result[0]["prior_evidence_id"], "old-path-id")
        self.assertEqual(note["evidence_id"], "old-path-id")
        for change in [{"original_source_sha256": "changed"}, {"episode_id": "different"}]:
            with self.assertRaises(ValueError):
                remap_clarifications([self.bound], [{**note, **change}])

    def test_duplicate_source_confirmation_resolution_rejected(self):
        self.bound["claims"] = [self.claim(), self.claim()]
        note = {"clarification_id": "N1", "cohort": "TEST", "patient_id": "999", "organism_key": "humanbetaherpesvirus5",
            "episode_id": "X:E1", "evidence_id": "old", "original_source_sha256": "hash", "original_source_row": "[0]"}
        with self.assertRaises(ValueError):
            remap_clarifications([self.bound], [note])

    def test_independent_respiratory_support_replaces_old_rank_gate(self):
        self.bound["claims"] = [self.claim()]
        self.candidate["analytical_profile"]["per_test_signals"][0]["rank_in_retained_universe"] = 20
        result, _ = self.apply()
        cmv = result["all_forwarded_candidates"][0]
        self.assertEqual(cmv["final_reporting_tier"], "Possible")
        self.assertEqual(cmv["clinical_decision"], "review_context_needed")
        self.assertEqual(cmv["viral_evidence"]["cmv_decision"]["axes"]["activity"], "A0")
        self.assertEqual(len(result["complete_report"]), 1)

    def test_fallback_recomputed_from_whole_patient_and_counted_positive(self):
        result, _ = self.apply()
        self.assertEqual(result["complete_report"][0]["final_reporting_tier"], "Fallback-Possible")
        self.assertTrue(result["complete_report"][0]["viral_evidence"]["counts_as_positive_prediction"])

    def test_non_cmv_picked_prevents_extra_cmv_fallback(self):
        bacterium = copy.deepcopy(self.candidate)
        bacterium.update(organism_name="Escherichia coli", final_reporting_tier="Picked",
                         selected_for_complete_report=True, possible_reporting_role=reporting.STRICT_ROLE)
        bacterium["possible_pathogen_gate"].update(reporting_role=reporting.STRICT_ROLE, selected_for_complete_report=True)
        self.packet["all_forwarded_candidates"].append(bacterium)
        result, _ = self.apply()
        self.assertEqual([c["organism_name"] for c in result["complete_report"]], ["Escherichia coli"])
        self.assertEqual(result["all_forwarded_candidates"][1], bacterium)

    def test_missing_method_stays_unknown_no_confirmation_invented(self):
        self.bound["claims"] = [self.claim(independence_class="not_established")]
        bacterium = copy.deepcopy(self.candidate)
        bacterium.update(organism_name="Escherichia coli", final_reporting_tier="Picked", selected_for_complete_report=True)
        self.packet["all_forwarded_candidates"].append(bacterium)
        result, _ = self.apply()
        self.assertEqual(result["all_forwarded_candidates"][0]["final_reporting_tier"], "Context")

    def test_original_packet_and_upstream_data_not_mutated(self):
        original = copy.deepcopy(self.packet)
        result, _ = self.apply()
        self.assertEqual(self.packet, original)
        self.assertEqual(result["all_forwarded_candidates"][0]["analytical_profile"], self.candidate["analytical_profile"])

    def test_cmv_without_binding_fails_closed(self):
        with self.assertRaises(ValueError):
            apply_patient_policy(self.packet, {}, [], self.policy, self.axes, self.base)


if __name__ == "__main__":
    unittest.main()

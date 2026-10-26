# All fixtures below are synthetic, not patient exports.
import copy
import unittest
from pathlib import Path

from tools.project_viral_evidence_patient_outputs import attach_annotations, without_annotations, pointer_value


class ProjectViralEvidenceTest(unittest.TestCase):
    def setUp(self):
        self.viral = {"patient_id": "901", "organism_name": "Human betaherpesvirus 5",
                      "reporting": {"tier": "Possible", "selected_for_complete_report": True},
                      "decision_path": {"final_reporting_tier": "Possible"}}
        self.bound = {"patient_id": "901", "current_tier": "Possible", "current_selected": True,
                      "episode_id": "SYNTHETIC:P901:mngs:1", "claims": [],
                      "frozen_hospital_observation_reconciliation": {"requested": 2, "resolved": 2},
                      "axes": {"source_mapping_status": "completed_for_available_sources",
                               "detection": "D2", "activity": "A0", "organ_support": "O1",
                               "gaps": ["activity_unknown"], "unique_supporting_evidence_count": 3,
                               "evidence_consumption_ledger": {}}}
        # Only the hash of this stable fixture file is needed by the projector.
        self.binding_file = Path(__file__)

    def annotate(self, payload):
        return attach_annotations(payload, {("901", "humanbetaherpesvirus5"): (0, self.bound)}, binding_file=self.binding_file)

    def test_all_views_annotated_without_changing_original_fields(self):
        other = {"patient_id": "901", "organism_name": "Stenotrophomonas maltophilia",
                 "reporting": {"tier": "Picked", "selected_for_complete_report": True}}
        payload = {"patient_id": "901", "final_report": [copy.deepcopy(self.viral), copy.deepcopy(other)],
                   "all_candidates": [copy.deepcopy(self.viral), copy.deepcopy(other)]}
        original = copy.deepcopy(payload)
        output, count = self.annotate(payload)
        self.assertEqual(count, 2)
        self.assertEqual(output["final_report"][0]["viral_evidence_shadow"], output["all_candidates"][0]["viral_evidence_shadow"])
        self.assertEqual(output["all_candidates"][1], other)
        self.assertEqual(without_annotations(output), original)
        self.assertEqual(payload, original)

    def test_stale_binding_cannot_attach_to_new_decision(self):
        payload = {"patient_id": "901", "all_candidates": [copy.deepcopy(self.viral)]}
        payload["all_candidates"][0]["reporting"]["tier"] = "Picked"
        with self.assertRaises(ValueError):
            self.annotate(payload)

    def test_patient_mismatch_is_rejected(self):
        payload = {"patient_id": "901", "all_candidates": [copy.deepcopy(self.viral)]}
        payload["all_candidates"][0]["patient_id"] = "902"
        with self.assertRaises(ValueError):
            self.annotate(payload)

    def test_empty_patient_retained_with_empty_lists(self):
        payload = {"patient_id": "903", "candidate_count": 0, "all_candidates": [], "final_report": []}
        output, count = self.annotate(payload)
        self.assertEqual(count, 0)
        self.assertEqual(without_annotations(output), payload)

    def test_full_reporting_packet_views_are_annotated(self):
        full = {"patient_id": "901", "organism_name": self.viral["organism_name"],
                "final_reporting_tier": "Possible", "selected_for_complete_report": True,
                "clinical_decision": "review_high_priority"}
        payload = {"patient_id": "901", "all_forwarded_candidates": [copy.deepcopy(full)],
                   "complete_report": [copy.deepcopy(full)]}
        output, count = self.annotate(payload)
        self.assertEqual(count, 2)
        self.assertEqual(without_annotations(output), payload)

    def test_source_pointer_can_resolve_nested_candidate(self):
        payload = [{"claims": [{"evidence_id": "same-source"}]}]
        self.assertEqual(pointer_value(payload, "[0].claims[0].evidence_id"), "same-source")


if __name__ == "__main__":
    unittest.main()

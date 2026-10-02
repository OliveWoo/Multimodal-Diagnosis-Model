import json
import tempfile
import unittest
from pathlib import Path

from tools.phenotype_event_context import candidate_event_view, read_patient_packet
from tools.build_rag_v2_inputs import build_case
from tools.review_missed_mngs_candidates import build_review_payload


class TestPhenotypeEventContext(unittest.TestCase):
    def test_old_and_new_phenotype_inputs_cannot_be_combined(self):
        root = Path("unused")
        with self.assertRaisesRegex(ValueError, "not both"):
            build_review_payload(patient_dir=root, queue_path=root, max_path=root,
                                 summary_path=root, mngs_to_specimen_path=None,
                                 phenotype_root=root, phenotype_event_root=root)
        with self.assertRaisesRegex(ValueError, "not both"):
            build_case(label="shadow", merged_suffix="merge", patient_id="9",
                       tier="review_context_needed", review_item={}, deterministic={},
                       review={}, summary={}, phenotype_root=root, phenotype_event_root=root)

    def test_candidate_view_excludes_parser_lab_values_and_other_organisms(self):
        packet = {
            "patient_id": 9,
            "decision_effect": "none_shadow_only",
            "phenotype_status": "provisional_link",
            "source": {"source_workbook": "patient 9.xlsx"},
            "events": [{
                "index_event": {"event_id": "P9:mngs:1", "collected_time": "2024-04-30 10:40"},
                "linkage": {"infection_episode_assessment": "index_infection_episode_unresolved",
                            "production_link_verified": False},
                "candidate_packets": [
                    {"organism_name": "Pneumocystis jirovecii", "canonical_key": "pneumocystisjirovecii",
                     "phenotype_context": [{"phenotype": "IMMUNOCOMPROMISED"}],
                     "exact_organism_workbook_microbiology": [], "model_use": "provisional_shadow_review_only_event_unverified"},
                    {"organism_name": "CMV", "canonical_key": "humanbetaherpesvirus5",
                     "phenotype_context": [], "exact_organism_workbook_microbiology": [],
                     "model_use": "provisional_shadow_review_only_event_unverified"},
                ],
                "source_records": {
                    "Clinical_Facts": [{"event_window_relation": "within_window", "description": "pneumonia"},
                                       {"event_window_relation": "outside_window", "description": "old event"}],
                    "Clinical_Timeline": [],
                    "Lab_Results": [{"test_name": "WBC", "value": "5.0"}],
                },
                "lab_packet_rows_temporally_retrieved": 1,
            }],
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "patient_9_phenotype_event_shadow.json"
            path.write_text(json.dumps(packet), encoding="utf-8")
            loaded = read_patient_packet(Path(directory), "P9")
            view = candidate_event_view(loaded, "Pneumocystis jirovecii")
            self.assertEqual(len(view["event_views"]), 1)
            event = view["event_views"][0]
            self.assertEqual(event["phenotype_context"][0]["phenotype"], "IMMUNOCOMPROMISED")
            self.assertEqual(len(event["dated_clinical_facts"]), 1)
            self.assertFalse(event["lab_values_supplied_to_model"])
            self.assertNotIn("Lab_Results", json.dumps(view))
            self.assertNotIn("CMV", json.dumps(view))
            with self.assertRaises(FileNotFoundError):
                read_patient_packet(Path(directory), "P8")


if __name__ == "__main__":
    unittest.main()

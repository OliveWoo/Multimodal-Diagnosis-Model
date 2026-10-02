import json
import tempfile
import unittest
from pathlib import Path

from tools.build_index_event_clinical_timeline_shadow import (
    availability_at_collection, build_patient, build_shadow, episode_window_relation, temporal_relation,
)


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value), encoding="utf-8")


class TestClinicalTimelineShadow(unittest.TestCase):
    def test_time_precision_does_not_invent_order(self):
        self.assertEqual(temporal_relation("2024-04", "2024-04-30 10:40:00"), "overlaps_index_precision_uncertain")
        self.assertEqual(temporal_relation("2024-03", "2024-04-30 10:40:00"), "before_index")
        self.assertEqual(temporal_relation("2024-05-01 23:59", "2024-04-30 10:40:00"), "after_index")
        self.assertEqual(temporal_relation(None, "2024-04-30 10:40:00"), "unknown")

    def test_collection_cutoff_uses_report_availability_not_only_exam_time(self):
        index = "2024-04-30 10:40:00"
        self.assertEqual(availability_at_collection("imaging_context", "2024-04-29 09:00", "2024-05-01 12:00", index), "post_index_not_available")
        self.assertEqual(availability_at_collection("imaging_context", "2024-04-29 09:00", "2024-04-30", index), "timing_uncertain")
        self.assertEqual(availability_at_collection("imaging_context", "2024-04-29 09:00", "2024-04-29 15:00", index), "confirmed_available_before_index")
        self.assertEqual(availability_at_collection("illness_lab", None, "2024-04-29", index), "confirmed_available_before_index")
        self.assertEqual(availability_at_collection("illness_lab", None, "2024-05-01", index), "post_index_not_available")
        self.assertEqual(availability_at_collection("host_medication", "2024-03", None, index), "preexisting_context_documentation_unverified")
        self.assertEqual(availability_at_collection("host_medication", None, None, index), "timing_unknown")

    def test_episode_window_keeps_later_exams_and_respects_date_precision(self):
        index = "2024-04-30 10:40:00"
        self.assertEqual(episode_window_relation("2024-05-01 23:59", index), "within_window")
        self.assertEqual(episode_window_relation("2024-05-02 10:40:00", index), "within_window")
        self.assertEqual(episode_window_relation("2024-05-02 15:20:42", index), "outside_window")
        self.assertEqual(episode_window_relation("2024-05-02", index), "boundary_uncertain")
        self.assertEqual(episode_window_relation("2024-05-04", index), "outside_window")
        self.assertEqual(episode_window_relation(None, index), "unknown")

    def test_source_roles_event_separation_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            patient = root / "NGS_patient_1_json"
            patient.mkdir()
            write_json(patient / "NGS_patient_1_mNGS_ranked_candidates.json", {
                "patient_id": "1",
                "records": [
                    {"specimen_code": "A", "specimen_site": "Blood", "collected_time": "2024-04-30 10:40:00", "seq_id": "DNA-A"},
                    {"specimen_code": "B", "specimen_site": "BALF", "collected_time": "2024-05-01 10:40:00", "seq_id": "DNA-B"},
                ],
            })
            write_json(patient / "NGS_patient_1_underlying.json", [
                {"underlying_diseases": ["Dermatomyositis"]}
            ])
            write_json(patient / "NGS_patient_1_admission_diagnosis.json", [
                {"notes": ["add Prednisolone (2024/03-), Methotrexate (2024/03-)"]}
            ])
            write_json(patient / "NGS_patient_1_CBC.json", [
                {"test": "CBC", "item": "WBC", "item_full": "WBC (x1000/uL)", "value": "6.65", "reported_time": "2024-05-02"}
            ])
            write_json(patient / "NGS_patient_1_image.json", [
                {"exam_type": "CT", "collected_time": "2024-05-01 23:59", "reported_time": "2024-05-02 15:20:42", "findings": [
                    "Differential diagnosis: PJP (Pneumocystis jirovecii pneumonia)."
                ]}
            ])
            packet = build_patient(patient)
            self.assertEqual(packet["schema_version"], "index_event_clinical_timeline_shadow.v3")
            self.assertEqual(len(packet["index_events"]), 2)
            self.assertEqual([event["specimen_site"] for event in packet["index_events"]], ["Blood", "BALF"])
            self.assertEqual(packet["decision_effect"], "none_shadow_only")
            self.assertEqual(packet["primary_analysis"], "retrospective_infection_episode_window")
            self.assertEqual(packet["episode_window_hours_each_side"], 48)
            self.assertFalse(packet["primary_view_is_at_collection_prediction"])
            self.assertFalse(packet["source_availability"]["other_lab"])
            self.assertEqual(packet["clinical_data_status"]["cbc"]["status"], "cbc_available_differential_missing")
            medications = [row for row in packet["clinical_evidence"] if row["evidence_role"] == "host_medication"]
            self.assertEqual({row["label"] for row in medications}, {"prednisolone", "methotrexate"})
            self.assertEqual(len({row["source_ref"]["source_claim_group"] for row in medications}), 1)
            self.assertTrue(all(row["details"]["dose"] == "Unknown" for row in medications))
            self.assertTrue(all(set(row["available_at_index"].values()) == {"preexisting_context_documentation_unverified"} for row in medications))
            lab = next(row for row in packet["clinical_evidence"] if row["evidence_role"] == "illness_lab")
            self.assertEqual(lab["time_basis"], "reported_time_only_observation_time_unknown")
            self.assertEqual(set(lab["observation_relation_to_index"].values()), {"unknown"})
            self.assertEqual(set(lab["available_at_index"].values()), {"post_index_not_available"})
            self.assertEqual(lab["episode_evidence_status"]["P1:mngs:2"], "report_only_observation_time_unverified")
            image = next(row for row in packet["clinical_evidence"] if row["evidence_role"] == "imaging_context")
            self.assertEqual(image["observation_relation_to_index"]["P1:mngs:1"], "after_index")
            self.assertEqual(image["observation_relation_to_index"]["P1:mngs:2"], "after_index")
            self.assertEqual(set(image["available_at_index"].values()), {"post_index_not_available"})
            self.assertEqual(set(image["episode_evidence_status"].values()), {"episode_observation"})
            self.assertEqual(image["observation_window_relation_to_index"]["P1:mngs:1"], "within_window")
            self.assertEqual(image["report_window_relation_to_index"]["P1:mngs:1"], "outside_window")
            self.assertFalse(image["microbiologic_confirmation"])
            self.assertIn("Pneumocystis jirovecii", [item["concept_name"] for item in image["diagnostic_mentions"]])
            self.assertEqual(image["diagnostic_mentions"][0]["assertion"], "differential")
            self.assertEqual(len(image["source_ref"]["sha256"]), 64)
            output = root / "shadow"
            summary = build_shadow(root, output)
            self.assertEqual(summary["counts"]["patients"], 1)
            self.assertEqual(summary["counts"]["episode:episode_observation"], 2)
            with self.assertRaises(FileExistsError):
                build_shadow(root, output)

    def test_multi_assay_cases_replace_selected_dna_event_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            patient = root / "NGS_patient_9_json"
            patient.mkdir()
            write_json(patient / "NGS_patient_9_mNGS_ranked_candidates.json", {
                "records": [{
                    "specimen_code": "A", "specimen_site": "Blood",
                    "collected_time": "2024-02-22 10:55:00", "seq_id": "A",
                }],
            })
            multi = root / "NGS_patient_9_multi_assay_mngs_evidence.json"
            write_json(multi, {
                "patient_id": "9",
                "cases": [
                    {
                        "case_review_id": "A", "specimen_code": "A",
                        "specimen_site": "Blood", "collected_time": "2024-02-22 10:55:00",
                        "tests": [{"seq_id": "DNA-A"}],
                    },
                    {
                        "case_review_id": "B", "specimen_code": "B",
                        "specimen_site": "BALF", "collected_time": "2024-02-22 10:55:00",
                        "tests": [{"seq_id": "DNA-B"}, {"seq_id": "RNA-B"}],
                    },
                ],
            })

            packet = build_patient(patient, multi)

        self.assertEqual("index_event_clinical_timeline_shadow.v4", packet["schema_version"])
        self.assertEqual("multi_assay_case_inventory", packet["index_event_source"])
        self.assertEqual(["A", "B"], [item["case_review_id"] for item in packet["index_events"]])
        self.assertEqual(["DNA-B", "RNA-B"], packet["index_events"][1]["seq_ids"])


if __name__ == "__main__":
    unittest.main()

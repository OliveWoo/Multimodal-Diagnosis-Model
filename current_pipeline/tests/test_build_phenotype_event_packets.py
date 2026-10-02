import json
import tempfile
import unittest
from pathlib import Path

from tools.build_phenotype_event_packets import build_packet, patient_number, run


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


class TestPhenotypeEventPackets(unittest.TestCase):
    def test_patient_number_accepts_trailing_data_quality_notes(self):
        self.assertEqual(patient_number("patient 31 - filmarray not the same day"), 31)
        self.assertEqual(patient_number("patient 47- multi filmarray"), 47)
        self.assertEqual(patient_number("patient_9"), 9)
        self.assertEqual(patient_number("P16"), 16)
        with self.assertRaises(ValueError):
            patient_number("filmarray note for patient 31")

    def fixtures(self, root):
        timeline = {
            "patient_id": "9",
            "index_events": [{"event_id": "P9:mngs:1", "collected_time": "2024-04-30 10:40:00", "specimen_site": "BALF"}],
            "clinical_evidence": [{"evidence_id": "P9:E1", "evidence_role": "host_medication",
                                   "source_text": "Prednisolone and Methotrexate for dermatomyositis"}],
        }
        ledger = {
            "patient_id": "patient 9",
            "source": {"workbook_file": "patient 9.xlsx", "patient_qc_status": "SYSTEM_REVIEW_REQUIRED"},
            "records": {
                "Phenotypes_Long": [{"id": "Phenotypes_Long:2", "row": {
                    "phenotype": "IMMUNOCOMPROMISED", "status": "YES", "confidence": "HIGH",
                    "phenotype_qc_status": "PASS", "source_date": "2024-04-29", "_source_row": 2,
                    "evidence": "Prednisolone and Methotrexate for dermatomyositis"}}],
                "Evidence": [{"id": "Evidence:2", "row": {
                    "phenotype": "IMMUNOCOMPROMISED", "source_date": "2024-04-29",
                    "evidence": "Prednisolone and Methotrexate for dermatomyositis"}}],
                "Clinical_Facts": [], "Clinical_Timeline": [],
                "Microbiology": [{"id": "Microbiology:2", "row": {
                    "organism": "Pneumocystis jirovecii", "date": "2024-05-01", "specimen": "BALF"}},
                    {"id": "Microbiology:3", "row": {
                        "organism": "Pneumocystis jirovecii", "date": "2024-06-01", "specimen": "BALF"}}],
                "Lab_Results": [{"id": "Lab_Results:2", "row": {
                    "date": "2024-04-30", "test_name": "WBC", "value": "5.0"}},
                    {"id": "Lab_Results:3", "row": {
                        "date": "2024-06-01", "test_name": "WBC", "value": "4.0"}}],
            },
        }
        candidate = {
            "patient_id": "9", "linkage": {"production_link_verified": False},
            "all_phenotype_overview": [],
            "candidate_contexts": [{"organism_name": "Pneumocystis jirovecii",
                                    "canonical_key": "pneumocystisjirovecii",
                                    "expanded_phenotypes": [{
                                        "source_id": "Phenotypes_Long:2", "phenotype": "IMMUNOCOMPROMISED",
                                        "corroborating_evidence": [{"source_id": "Evidence:2"}]}]}],
        }
        link = {"pipeline_patient_number": 9,
                "patient_identity_assessment": "local_strong_same_patient_support",
                "infection_episode_assessment": "index_infection_episode_unresolved"}
        timeline_path = root / "timeline" / "NGS_patient_9_clinical_timeline_shadow.json"
        ledger_path = root / "phenotype" / "patient_9_phenotype_decision_evidence_v1.json"
        candidate_path = root / "candidates" / "patient_9_candidate_phenotype_shadow_v1.json"
        write(timeline_path, timeline)
        write(ledger_path, ledger)
        write(candidate_path, candidate)
        return timeline, ledger, link, candidate, timeline_path, ledger_path, candidate_path

    def test_provisional_packet_aligns_dates_and_marks_derived_overlap(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.fixtures(Path(directory))
            packet = build_packet(*args[:4], timeline_path=args[4], ledger_path=args[5], candidate_path=args[6])
            event = packet["events"][0]
            self.assertFalse(event["linkage"]["production_link_verified"])
            self.assertEqual(event["candidate_count"], 1)
            self.assertEqual(event["lab_source_total_in_ledger"], 2)
            self.assertEqual(event["lab_packet_rows_temporally_retrieved"], 1)
            phenotype = event["candidate_packets"][0]["phenotype_context"][0]
            self.assertEqual(phenotype["event_alignment"]["event_window_relation"], "within_window")
            self.assertEqual(phenotype["event_alignment"]["possible_same_kh_evidence_ids"], ["P9:E1"])
            self.assertFalse(phenotype["independent_corroboration"])
            self.assertEqual(phenotype["source_evidence_alignment"][0]["source_id"], "Evidence:2")
            self.assertEqual(event["candidate_packets"][0]["exact_organism_workbook_microbiology"][0]["specimen"], "BALF")
            self.assertEqual(len(event["candidate_packets"][0]["exact_organism_workbook_microbiology"]), 1)
            self.assertEqual(len(event["source_records"]["Microbiology"]), 2)
            self.assertEqual(event["linkage"]["infection_episode_assessment"], "index_infection_episode_unresolved")

    def test_missing_phenotype_stays_missing_and_output_does_not_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            timeline, _, link, candidate, timeline_path, _, candidate_path = self.fixtures(root)
            packet = build_packet(timeline, None, link, candidate, timeline_path=timeline_path,
                                  ledger_path=None, candidate_path=candidate_path)
            self.assertEqual(packet["phenotype_status"], "phenotype_file_missing")
            self.assertEqual(packet["events"][0]["candidate_packets"][0]["model_use"], "unavailable")
            write(root / "audit.json", {"patients": [link]})
            output = root / "result"
            summary = run(root / "timeline", root / "phenotype", root / "audit.json", root / "candidates", output)
            self.assertEqual(summary["patient_count"], 1)
            with self.assertRaises(FileExistsError):
                run(root / "timeline", root / "phenotype", root / "audit.json", root / "candidates", output)

    def test_candidate_event_refs_prevent_cross_event_attachment(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            timeline, ledger, link, candidate, timeline_path, ledger_path, candidate_path = self.fixtures(root)
            timeline["index_events"] = [
                {"event_id": "P9:mngs:1", "case_review_id": "CASE-A",
                 "collected_time": "2024-04-30 10:40:00", "specimen_site": "Blood"},
                {"event_id": "P9:mngs:2", "case_review_id": "CASE-B",
                 "collected_time": "2024-04-30 10:40:00", "specimen_site": "BALF"},
            ]
            candidate["candidate_contexts"][0]["candidate_event_refs"] = ["CASE-B"]
            packet = build_packet(timeline, ledger, link, candidate, timeline_path=timeline_path,
                                  ledger_path=ledger_path, candidate_path=candidate_path)
            self.assertEqual(packet["events"][0]["candidate_count"], 0)
            self.assertEqual(packet["events"][1]["candidate_count"], 1)
            self.assertEqual(
                packet["events"][1]["candidate_packets"][0]["candidate_event_refs"],
                ["CASE-B"],
            )


if __name__ == "__main__":
    unittest.main()

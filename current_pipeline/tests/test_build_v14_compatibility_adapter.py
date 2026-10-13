import json
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

from tools.build_v14_compatibility_adapter import (
    adapt_patient,
    alignment,
    is_explicit_positive,
)


class V14CompatibilityAdapterTests(unittest.TestCase):
    def test_alignment_uses_two_calendar_date_window(self):
        events = [("event", datetime(2024, 1, 3, 12, 0))]
        self.assertTrue(alignment({"collected_time": "2024-01-01 12:00"}, events)["passed"])
        boundary = alignment({"collected_time": "2024-01-01 00:01"}, events)
        self.assertTrue(boundary["passed"])
        self.assertEqual(boundary["calendar_day_delta"], -2)
        self.assertEqual(boundary["reason"], "within_two_calendar_dates")
        self.assertFalse(alignment({"collected_time": "2023-12-31 23:59"}, events)["passed"])
        self.assertEqual(alignment({}, events)["reason"], "time_unavailable")

    def test_positive_result_must_be_explicit(self):
        self.assertTrue(is_explicit_positive({"quantitation_status": "positive_detected_quantity_unknown"}))
        self.assertTrue(is_explicit_positive({"result": "Positive"}))
        self.assertFalse(is_explicit_positive({"result": "Negative"}))
        self.assertFalse(is_explicit_positive({"result": ""}))

    def test_canonical_event_aligned_addition_and_non_regression(self):
        with tempfile.TemporaryDirectory() as directory:
            patient_dir = Path(directory)
            culture_path = patient_dir / "NGS_patient_3_culture.json"
            culture_path.write_text(json.dumps({
                "record_id": "culture-1",
                "source": {"source_file": "patient_3.xlsx", "source_sheet": "Lab", "row": 10},
            }), encoding="utf-8")
            molecular_path = patient_dir / "NGS_patient_3_molecular_microbiology.json"
            molecular_path.write_text(json.dumps({
                "record_id": "molecular-1",
                "source": {"source_file": "patient_3.xlsx", "source_sheet": "Lab", "row": 11},
            }), encoding="utf-8")

            baseline = {"hospital_organism_evidence": [
                {
                    "organism_name": "Klebsiella pneumoniae",
                    "classification": "Bacterial",
                    "evidence_modules": ["culture"],
                    "module_level_summary": {"culture": "Level 2"},
                    "module_evidence": {"culture": [{"legacy_verified": True}]},
                    "best_hospital_level": "Level 2",
                },
                {
                    "organism_name": "COVID-19 Ag rapid test",
                    "classification": "Unknown",
                    "evidence_modules": ["molecular_microbiology"],
                    "module_level_summary": {"molecular_microbiology": "Level 1"},
                    "module_evidence": {},
                    "best_hospital_level": "Level 1",
                },
            ]}
            bridge = {"hospital_organism_evidence": [
                {
                    "organism_name": "Escherichia coli",
                    "classification": "Bacterial",
                    "module_level_summary": {"culture": "Level 2"},
                    "module_evidence": {"culture": [{
                        "specimen_type": "BAL",
                        "specimen_category": "Lower_Respiratory",
                        "quantitation_status": "positive_detected_quantity_unknown",
                        "collected_time": "2024-01-03 12:00",
                        "reported_time": "2024-01-04 12:00",
                        "evidence_observation_ids": ["CUL-001"],
                    }]},
                },
                {
                    "organism_name": "COVID-19 Ag rapid test",
                    "classification": "Unknown",
                    "module_level_summary": {"molecular_microbiology": "Level 1"},
                    "module_evidence": {"molecular_microbiology": [{
                        "specimen_type": "BAL",
                        "specimen_category": "Lower_Respiratory",
                        "result": "Positive",
                        "collected_time": "2024-01-03 12:00",
                        "evidence_observation_ids": ["MOL-001"],
                    }]},
                },
            ]}

            output, audit, summary = adapt_patient(
                patient_id=3,
                patient_dir=patient_dir,
                bridge_packet=bridge,
                events=[("event", datetime(2024, 1, 3, 12, 0))],
                baseline_summary=baseline,
                v14_summary={},
                allow_v14=True,
            )

            names = {item["organism_name"] for item in output["hospital_organism_evidence"]}
            self.assertIn("Klebsiella pneumoniae", names)
            self.assertIn("Escherichia coli", names)
            self.assertNotIn("COVID-19 Ag rapid test", names)
            klebsiella = next(item for item in output["hospital_organism_evidence"] if item["organism_name"] == "Klebsiella pneumoniae")
            self.assertEqual(klebsiella["best_hospital_level"], "Level 2")
            ecoli = next(item for item in output["hospital_organism_evidence"] if item["organism_name"] == "Escherichia coli")
            source_refs = ecoli["module_evidence"]["culture"][0]["source_references"]
            self.assertEqual(source_refs[0]["source_row"], 10)
            self.assertTrue(ecoli["module_evidence"]["culture"][0]["event_alignment"]["passed"])
            self.assertTrue(summary["non_regression_level_passed"])
            self.assertTrue(any(row["reason"] == "taxonomy_unmapped_or_uncertain" for row in audit))


if __name__ == "__main__":
    unittest.main()

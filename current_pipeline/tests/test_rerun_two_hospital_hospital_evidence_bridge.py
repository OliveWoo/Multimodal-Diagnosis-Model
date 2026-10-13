import json
import tempfile
import unittest
from pathlib import Path

from tools import rerun_two_hospital_hospital_evidence_bridge as bridge


class TwoHospitalHospitalEvidenceBridgeTest(unittest.TestCase):
    def test_runs_without_mngs_and_preserves_timing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            patient = root / "patients" / "NGS_patient_54_json"
            patient.mkdir(parents=True)
            (patient / "NGS_patient_54_culture.json").write_text(
                json.dumps([{
                    "sample": "Sputum",
                    "organism": "Klebsiella pneumoniae",
                    "result": "Isolated",
                    "collected_time": "2026-10-01",
                    "reported_time": "2026-10-02",
                }]),
                encoding="utf-8",
            )
            output = root / "output"

            summary = bridge.run(root / "patients", output)

            self.assertEqual(1, summary["patient_count"])
            self.assertEqual(1, summary["hospital_organism_evidence_count"])
            self.assertFalse(summary["mngs_scoring_performed"])
            packet = json.loads(
                (output / "patient_packets" / "NGS_patient_54_hospital_evidence_bridge.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual("pending_mngs_anchor", packet["event_alignment_status"])
            self.assertEqual(
                "Klebsiella pneumoniae",
                packet["hospital_organism_evidence"][0]["organism_name"],
            )


if __name__ == "__main__":
    unittest.main()

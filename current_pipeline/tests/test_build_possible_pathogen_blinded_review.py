import json
import tempfile
import unittest
from pathlib import Path

from tools.build_possible_pathogen_blinded_review import (
    forbidden_paths,
    make_case,
    run,
)


def candidate():
    return {
        "patient_id": "10",
        "organism_name": "Staphylococcus aureus",
        "candidate_source": "mngs_positive_reads",
        "clinical_decision": "review_high_priority",
        "possible_reporting_role": "analytical_or_direct_possible_pathogen",
        "taxonomy_profile": {
            "biological_class": "bacterium",
            "mapping_status": "exact_species",
            "primary_rule_family": "typical_respiratory_pathogen",
        },
        "analytical_profile": {
            "best_rank_in_retained_universe": 4,
            "selected_positive_test_count": 3,
            "dna_selected": True,
            "rna_selected": True,
            "cross_molecule_selected": True,
            "reproducibility_axis": True,
            "per_test_signals": [
                {"seq_id": "DNA-1", "nucleic_type": "DNA", "reads": 89, "rank_in_retained_universe": 5},
                {"seq_id": "RNA-1", "nucleic_type": "RNA", "reads": 280, "rank_in_retained_universe": 6},
            ],
        },
        "hospital_profile": {
            "hospital_evidence_detail": {
                "source_observations": [{
                    "observation_id": "ASSAY-1",
                    "specimen_type": "BAL",
                    "target": "Staphylococcus aureus",
                    "raw_result": "Detected; 10^5 copy/mL",
                }]
            }
        },
        "source_candidate": {
            "specimen_context": "lower_respiratory",
            "collected_time": "2024-02-26",
            "case_review_id": "PRIVATE-PATIENT-ID",
        },
        "possible_pathogen_gate": {
            "route": "event_aligned_direct_detection_or_level3_culture",
            "cautions": [],
            "axes": {"direct_hospital_level": 3},
        },
        "promotion_gate": {
            "axes": {
                "direct_evidence_timing_profile": {
                    "event_window_hours": 48,
                    "event_aligned_positive_count": 1,
                    "rows": [{
                        "observation_id": "ASSAY-1", "positive": True,
                        "timing": "within_event_window", "specimen_type": "BAL",
                    }],
                }
            }
        },
        "host_evidence": [],
        "image_evidence": [],
        "phenotype_evidence": {},
        "reasons": [],
        "rule_ids": [],
    }


class TestPossiblePathogenBlindedReview(unittest.TestCase):
    def test_case_excludes_patient_number_and_forbidden_fields(self):
        case = make_case("BLIND-001", candidate(), [], "abc")
        serialized = json.dumps(case, ensure_ascii=False)
        self.assertEqual("BLIND-001", case["patient_id"])
        self.assertNotIn('"patient_id": "10"', serialized)
        self.assertNotIn("PRIVATE-PATIENT-ID", serialized)
        self.assertEqual([], forbidden_paths(case))
        self.assertEqual(2, len(case["local_evidence"]["mngs_per_test_evidence"]["per_test_signals"]))

    def test_forbidden_benchmark_key_is_detected(self):
        self.assertEqual(["$.benchmark_answer"], forbidden_paths({"benchmark_answer": "x"}))

    def test_run_separates_private_case_map(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "input"
            patient_dir = root / "patient_outputs"
            patient_dir.mkdir(parents=True)
            payload = {
                "possible_pathogens": [candidate()],
                "strict_picked": [],
            }
            source = patient_dir / "NGS_patient_10_test_aware_possible_pathogen_shadow.json"
            source.write_text(json.dumps(payload), encoding="utf-8")
            (root / "summary.json").write_text("{}", encoding="utf-8")
            output = Path(temp) / "output"
            manifest = run(root, output)
            self.assertEqual(1, manifest["case_count"])
            production = (output / "production_queue.jsonl").read_text(encoding="utf-8")
            private = (output / "private_do_not_send_to_model" / "case_identity_map.json").read_text(encoding="utf-8")
            self.assertNotIn('"patient_id": "10"', production)
            self.assertIn('"patient_id": "10"', private)


if __name__ == "__main__":
    unittest.main()

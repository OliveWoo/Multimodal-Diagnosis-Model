from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.build_two_hospital_v2_6_compatibility_input import (
    build_multi_assay_document,
    infer_molecule,
    parse_args,
    supplied_multi_assay_document,
)


class TwoHospitalCompatibilityInputTest(unittest.TestCase):
    def test_patient_subset_cli_is_explicit(self) -> None:
        args = parse_args(["--patient-ids", "54, 67,90"])
        self.assertEqual("54, 67,90", args.patient_ids)

    def test_supplied_full_multi_assay_is_loaded_without_reconstruction(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "NGS_patient_55_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
            payload = {
                "patient_id": 55,
                "cases": [{"tests": [{"seq_id": "DNA1", "nucleic_type": "DNA"}]}],
                "organisms": [{
                    "name": "Example bacterium",
                    "measurements": [{"seq_id": "DNA1", "reads": 5}],
                }],
            }
            path.write_text(json.dumps(payload), encoding="utf-8")
            observed, observed_path = supplied_multi_assay_document(root, 55)
        self.assertEqual(payload, observed)
        self.assertEqual(path, observed_path)

    def test_molecule_inference_is_explicit_and_conservative(self) -> None:
        self.assertEqual(("DNA", "inferred_from_source_path"), infer_molecule("sample_DNA.zip"))
        self.assertEqual(("RNA", "inferred_from_source_path"), infer_molecule("sample_RNA.zip"))
        self.assertEqual(("UNKNOWN", "source_unknown"), infer_molecule("sample.zip"))
        self.assertEqual(("UNKNOWN", "source_unknown"), infer_molecule("DNA_RNA_mixed.zip"))

    def test_adapter_preserves_reads_and_unknowns_without_inventing_rpm(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            patient_dir = Path(temporary)
            (patient_dir / "NGS_patient_55_mNGS_ranked_candidates.json").write_text(
                json.dumps({
                    "records": [{
                        "seq_id": "SEQ1",
                        "specimen_code": "SEQ1",
                        "collected_time": "2024-01-01 00:00:00",
                        "specimen_site": "BALF",
                    }]
                }),
                encoding="utf-8",
            )
            rows = [{
                "specimen_id": "CASE1",
                "seq_id": "SEQ1",
                "source_category": "1.Bac",
                "source_file": "legacy_without_molecule.zip::1.Bac/result.tsv",
                "organism_name": "Streptococcus pneumoniae",
                "selection_reasons": ["Code_RK_NTC=RK0"],
                "row_index": 3,
                "row": {"Sec.hit": "19", "taxid": "1313", "Code_RK_NTC": "RK-0"},
            }]

            document = build_multi_assay_document(patient_dir, 55, rows)

        test = document["cases"][0]["tests"][0]
        measurement = document["organisms"][0]["measurements"][0]
        self.assertEqual("UNKNOWN", test["nucleic_type"])
        self.assertEqual("source_unknown", test["molecule_status"])
        self.assertIsNone(test["parameters"]["input_reads"])
        self.assertIsNone(test["parameters"]["host_ratio"])
        self.assertEqual(19.0, measurement["reads"])
        self.assertTrue(measurement["selected"])
        self.assertEqual("1313", measurement["source_taxid"])

    def test_duplicate_export_paths_do_not_become_independent_observations(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            patient_dir = Path(temporary)
            rows = []
            for source_file in ("old_folder/result.tsv", "new_folder/result.tsv"):
                rows.append({
                    "specimen_id": "CASE1",
                    "seq_id": "SEQ1",
                    "source_category": "1.Bac",
                    "source_file": source_file,
                    "organism_name": "Streptococcus pneumoniae",
                    "selection_reasons": ["Code_RK_NTC=RK0"],
                    "row_index": 3,
                    "row": {"Sec.hit": "19", "taxid": "1313", "Code_RK_NTC": "RK-0"},
                })

            document = build_multi_assay_document(patient_dir, 55, rows)

        measurements = document["organisms"][0]["measurements"]
        self.assertEqual(1, len(measurements))
        self.assertEqual(2, measurements[0]["duplicate_export_row_count"])
        self.assertEqual(2, len(measurements[0]["duplicate_export_source_files"]))
        self.assertFalse(measurements[0]["duplicate_export_reads_conflict"])
        self.assertFalse(measurements[0]["duplicate_export_taxid_conflict"])


if __name__ == "__main__":
    unittest.main()

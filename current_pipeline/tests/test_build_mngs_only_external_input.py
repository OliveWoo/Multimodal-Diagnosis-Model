from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.build_mngs_only_external_input import SUMMARY_SUFFIX, run


class BuildMngsOnlyExternalInputTest(unittest.TestCase):
    def test_builds_empty_clinical_root_without_copying_sources(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory = root / "inventory"
            patients = inventory / "patients"
            patients.mkdir(parents=True)
            (patients / "NGS_patient_58_multi_assay_mngs_evidence.json").write_text(
                json.dumps({"patient_id": "58", "cases": [{"tests": []}]}),
                encoding="utf-8",
            )
            output = root / "mngs_only"

            manifest = run(inventory, output)

            self.assertEqual([58], manifest["patient_ids"])
            self.assertEqual(0, manifest["clinical_files_copied"])
            summary_path = (
                output / "NGS_patient_58_json" / "summary_outputs"
                / f"NGS_patient_58_{SUMMARY_SUFFIX}.json"
            )
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual("mngs_only", summary["decision_input_mode"])
            self.assertEqual([], summary["hospital_organism_evidence"])
            self.assertEqual([], summary["pathogen_candidates"])
            self.assertFalse(summary["host_context"]["opportunistic_host_support"])
            self.assertEqual({}, summary["source_files"])
            self.assertTrue((output / "phenotype_unavailable").is_dir())
            self.assertEqual(
                str((output / "phenotype_unavailable").resolve()),
                manifest["phenotype_root"],
            )

    def test_refuses_to_overwrite_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory = root / "inventory" / "patients"
            inventory.mkdir(parents=True)
            (inventory / "NGS_patient_1_multi_assay_mngs_evidence.json").write_text(
                json.dumps({"patient_id": 1}), encoding="utf-8"
            )
            output = root / "output"
            output.mkdir()
            (output / "keep.txt").write_text("keep", encoding="utf-8")

            with self.assertRaises(FileExistsError):
                run(root / "inventory", output)


if __name__ == "__main__":
    unittest.main()

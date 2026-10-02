import json
import csv
import tempfile
import unittest
from pathlib import Path

from tools import run_external_integrated_release_v2 as release


class ExternalIntegratedReleaseV2Test(unittest.TestCase):
    def test_cohort_partition_excludes_development_patients(self):
        with tempfile.TemporaryDirectory() as tmp:
            freeze = Path(tmp) / "freeze.json"
            freeze.write_text(
                json.dumps({"expected_patient_ids": [1, 2, 3]}),
                encoding="utf-8",
            )
            result = release.cohort_partition({1, 2, 3, 4, 5}, freeze)
            self.assertEqual(result["development_overlap_patient_ids"], [1, 2, 3])
            self.assertEqual(result["heldout_patient_ids"], [4, 5])
            self.assertTrue(result["patient_sets_disjoint"])

    def test_cohort_partition_accepts_integrated_release_patient_ids_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            freeze = Path(tmp) / "freeze.json"
            freeze.write_text(
                json.dumps({"patient_ids": [3, 4]}), encoding="utf-8"
            )
            result = release.cohort_partition({3, 4, 40}, freeze)
            self.assertEqual(result["heldout_patient_ids"], [40])

    def test_freeze_heldout_outputs_filters_before_any_answer_use(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reporting = root / "10_reporting_v3_route_a"
            reporting.mkdir(parents=True)
            fields = ["patient_id", "organism_name", "final_reporting_tier"]
            for path in (
                root / "integrated_decisions.csv",
                reporting / "complete_report.csv",
            ):
                with path.open("w", encoding="utf-8-sig", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows([
                        {"patient_id": "1", "organism_name": "A", "final_reporting_tier": "Picked"},
                        {"patient_id": "4", "organism_name": "B", "final_reporting_tier": "Context"},
                        {"patient_id": "5", "organism_name": "C", "final_reporting_tier": "Possible"},
                    ])
            result = release.freeze_heldout_outputs(root, {4, 5})
            self.assertEqual(result["heldout_patient_count"], 2)
            self.assertEqual(result["integrated_candidate_count"], 2)
            self.assertEqual(
                result["final_tier_counts"], {"Context": 1, "Possible": 1}
            )
            self.assertTrue(result["patient_coverage_complete"])
            self.assertFalse(result["answers_loaded"])

    def test_freeze_existing_summaries_requires_matching_patient_sets(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            patients = root / "patients"
            inventory = root / "inventory"
            (patients / "NGS_patient_1_json").mkdir(parents=True)
            (inventory / "patients").mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, "Patient/inventory mismatch"):
                release.freeze_existing_summaries(
                    patients,
                    inventory,
                    root / "output",
                    source_suffix="source",
                )

    def test_freeze_existing_summaries_copies_validated_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            patients = root / "patients"
            inventory = root / "inventory"
            for patient_id in (1, 2):
                summary_dir = (
                    patients / f"NGS_patient_{patient_id}_json" / "summary_outputs"
                )
                summary_dir.mkdir(parents=True)
                (summary_dir / f"NGS_patient_{patient_id}_source.json").write_text(
                    json.dumps({"patient_id": patient_id}), encoding="utf-8"
                )
                inventory_path = (
                    inventory
                    / "patients"
                    / f"NGS_patient_{patient_id}_multi_assay_mngs_evidence.json"
                )
                inventory_path.parent.mkdir(parents=True, exist_ok=True)
                inventory_path.write_text("{}", encoding="utf-8")
            output = root / "output"
            summary = release.freeze_existing_summaries(
                patients, inventory, output, source_suffix="source"
            )
            self.assertEqual(summary["patient_count"], 2)
            self.assertTrue(
                (output / "NGS_patient_1_release_input_external_51_v1.json").is_file()
            )

    def test_provisional_linkage_manifest_keeps_missingness_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            phenotype = root / "phenotype"
            phenotype.mkdir()
            (phenotype / "patient_9_phenotype_decision_evidence_v1.json").write_text(
                "{}", encoding="utf-8"
            )
            payload = release.build_provisional_linkage_manifest(
                {9, 14}, phenotype, root / "linkage.json"
            )
            by_id = {
                item["pipeline_patient_number"]: item
                for item in payload["patients"]
            }
            self.assertFalse(by_id[9]["production_link_verified"])
            self.assertEqual(by_id[9]["phenotype_duplicate_export_cluster"], [9, 11])
            self.assertEqual(
                by_id[14]["local_linkage_status"], "phenotype_file_missing"
            )


if __name__ == "__main__":
    unittest.main()

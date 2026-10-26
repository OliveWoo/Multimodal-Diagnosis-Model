"""Synthetic configuration and publication boundaries, no clinical exports."""
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from tools.public_integrated_config import cohort_for, cohort_runs, load_config
from tools.evaluate_current_integrated_cmv import verify_freeze


class PublicIntegratedConfigTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        for name in ("inventory/patients", "patients", "phenotype", "reference/SYNTHETIC/patient_results"):
            (self.root / name).mkdir(parents=True)
        self.write("inventory/patients/NGS_patient_901_multi_assay_mngs_evidence.json", {"cases": []})
        self.write("reference/SYNTHETIC/patient_results/patient_901_integrated_decision.json", {"patient_id": "901"})
        self.write("source.json", {"answer_blind_decision_stages": True, "input": {
            "inventory_root": str(self.root / "inventory"), "patient_root": str(self.root / "patients"),
            "phenotype_root": str(self.root / "phenotype"), "source_summary_suffix": "synthetic"}})
        self.write("overlay.json", {"exact_organisms": {}})
        self.write("completion.json", {"staged_overlay": str(self.root / "overlay.json"),
            "staged_overlay_sha256": hashlib.sha256((self.root / "overlay.json").read_bytes()).hexdigest()})
        self.write("episodes.json", {"clusters": []})
        self.write("development.json", {"expected_patient_ids": [901]})
        self.write("clarifications.json", {"clarifications": []})
        self.config = {"schema_version": "public_integrated_run.v1", "source_manifests": {"SYNTHETIC": "source.json"},
            "taxonomy_completion_manifest": "completion.json", "episode_identity_file": "episodes.json",
            "source_clarifications": "clarifications.json", "reference_root": "reference",
            "development_freeze": "development.json"}

    def write(self, name, value):
        path = self.root / name
        path.write_text(json.dumps(value), encoding="utf-8")
        return path

    def load(self, config=None):
        return load_config(self.write("run.private.json", config or self.config))

    def test_explicit_private_paths_resolve_without_embedded_cohort_roster(self):
        config = self.load()
        self.assertEqual(cohort_runs(config), {"SYNTHETIC": "SYNTHETIC"})
        self.assertEqual(cohort_for(config, "SYNTHETIC", "901"), "SYNTHETIC")
        self.assertTrue(config["source_manifests"]["SYNTHETIC"].is_absolute())
        self.assertIn(str(self.root / "clarifications.json"), config["_input_hashes"])
        self.assertIn(str(self.root / "overlay.json"), config["_input_hashes"])

    def test_answer_configuration_rejected_at_inference_boundary(self):
        for key in ("answers", "answer_contract", "labels", "unexpected"):
            with self.assertRaises(ValueError):
                self.load({**self.config, key: "must_not_be_read.csv"})

    def test_missing_private_source_fails_before_any_run(self):
        with self.assertRaises(FileNotFoundError):
            self.load({**self.config, "source_clarifications": "missing.private.json"})

    def test_partition_covers_inventory_exactly(self):
        config = self.load({**self.config, "partitions": {"SYNTHETIC": {"SUBSET": [901]}}})
        self.assertEqual(cohort_for(config, "SYNTHETIC", "901"), "SUBSET")
        self.assertEqual(cohort_runs(config), {"SUBSET": "SYNTHETIC"})
        with self.assertRaises(ValueError):
            cohort_for(config, "SYNTHETIC", "902")

    def test_partition_overlap_missing_and_extra_members_fail_closed(self):
        for groups in ({"FIRST": [901], "SECOND": [901]}, {"FIRST": [902]}, {"FIRST": [901, 902]}, {"FIRST": []}):
            with self.assertRaises(ValueError):
                self.load({**self.config, "partitions": {"SYNTHETIC": groups}})

    def test_unsafe_directory_names_rejected(self):
        for name in ("../escape", "run/name", ".", ""):
            with self.assertRaises(ValueError):
                self.load({**self.config, "source_manifests": {name: "source.json"}})

    def test_missing_baseline_is_not_silently_skipped(self):
        (self.root / "reference/SYNTHETIC/patient_results/patient_901_integrated_decision.json").rename(self.root / "baseline.json")
        with self.assertRaises(FileNotFoundError):
            self.load()

    def test_taxonomy_overlay_hash_must_match(self):
        self.write("overlay.json", {"changed": True})
        with self.assertRaises(ValueError):
            self.load()

    def test_source_manifest_must_declare_answer_blind_stages(self):
        source = json.loads((self.root / "source.json").read_text())
        source["answer_blind_decision_stages"] = False
        self.write("source.json", source)
        with self.assertRaises(ValueError):
            self.load()

    def test_decision_freeze_required_before_evaluation(self):
        self.write("summary.json", {"completed": True})
        self.write("answer_blind_freeze.json", {"passed": False, "answers_loaded": False})
        with self.assertRaises(ValueError):
            verify_freeze(self.root)
        self.write("answer_blind_freeze.json", {"passed": True, "answers_loaded": False,
            "summary_sha256": "not-the-current-hash", "artifacts_sha256": {}})
        with self.assertRaises(ValueError):
            verify_freeze(self.root)


if __name__ == "__main__":
    unittest.main()

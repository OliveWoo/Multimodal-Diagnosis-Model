"""Offline tests for batch completeness and patient-pathogen scoring."""
import copy
from pathlib import Path
import unittest
from unittest.mock import patch

import astra_batch as batch
import direct_raw_runner as r
from evaluate_astra_benchmark import metrics, score_row
from input_identity import sha256_file


class AstraBatchTests(unittest.TestCase):
    def setUp(self):
        self.summary = {"failures": [], "success_patient_ids": ["P3", "P4"]}
        self.config = {
            "model": batch.MODEL, "reasoning_effort": "medium", "max_output_tokens": 4000,
            "store": False, "tools": [],
            "runner_sha256": sha256_file(r.SCRIPT_DIR / "direct_raw_runner.py"),
            "output_schema_sha256": r.sha256_text(r.stable_json(r.OUTPUT_SCHEMA)),
        }

    def verify(self, summary=None, config=None):
        with patch.object(r, "load_json", side_effect=[summary or self.summary, config or self.config]):
            batch.verify_execution(Path("unused_test_path"), ["P3", "P4"])

    def test_complete_matching_execution_passes(self):
        self.verify()

    def test_missing_patient_rejected(self):
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            self.verify(summary={"failures": [], "success_patient_ids": ["P3"]})

    def test_failed_execution_rejected(self):
        with self.assertRaisesRegex(ValueError, "failed"):
            self.verify(summary={"failures": [{"patient_id": "P4"}], "success_patient_ids": ["P3", "P4"]})

    def test_changed_model_cap_tools_or_effort_rejected(self):
        for field, value in [("model", "gpt-5.6-sol"), ("max_output_tokens", 8000),
                             ("tools", [{"type": "web_search"}]), ("reasoning_effort", "high")]:
            with self.subTest(field=field):
                config = copy.deepcopy(self.config)
                config[field] = value
                with self.assertRaisesRegex(ValueError, "Different request"):
                    self.verify(config=config)

    def test_patient_pairs_not_globally_deduplicated(self):
        rows = [score_row("3", "kmuh_answer_positive_30", {"organism_a"}, {"organism_a"}),
                score_row("4", "kmuh_answer_positive_30", {"organism_b"}, {"organism_a"})]
        total = metrics(rows)["all_41"]
        self.assertEqual((total["tp"], total["fp"], total["fn"]), (1, 1, 1))
        self.assertEqual(total["precision"], .5)
        self.assertEqual(total["recall"], .5)

    def test_different_candida_species_remain_different(self):
        row = score_row("3", "kmuh_answer_positive_30", {"candida albicans"}, {"candida tropicalis"})
        self.assertEqual(row["tp"], [])
        self.assertEqual(row["fp"], ["candida tropicalis"])
        self.assertEqual(row["fn"], ["candida albicans"])

    def test_no_prediction_retains_false_negatives(self):
        row = score_row("3", "kmuh_answer_positive_30", {"organism_a"}, set())
        self.assertEqual(row["fn"], ["organism_a"])
        self.assertEqual(row["precision"], 0)
        self.assertEqual(row["recall"], 0)


if __name__ == "__main__":
    unittest.main()

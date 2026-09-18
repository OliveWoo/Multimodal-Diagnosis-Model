import unittest

from rag_re.evaluation import evaluate_outputs


def output(patient_id, candidates, predicted, baseline=None, abstained=None):
    return {
        "schema_version": "rag_re.output.v1",
        "run_status": "complete",
        "pipeline_version": "test",
        "pipeline_fingerprint": "test-pipeline",
        "patient_id": str(patient_id),
        "config_hash": "test-config",
        "prompt_sha256": "test-prompt",
        "evidence_provenance": {"source_config_hash": "test-evidence"},
        "execution": {"literature_attempted": True},
        "input_meta": {"candidate_pool_provenance": "original_candidate_pool"},
        "candidates": [{"organism_name": value} for value in candidates],
        "experiments": {
            "E0_BASELINE": {
                "semantics": "clinical_baseline",
                "expression": "baseline",
                "predicted_pathogens": baseline or [],
                "abstained_pathogens": [],
            },
            "A_LITERATURE": {
                "semantics": "screening_trigger",
                "expression": "A",
                "predicted_pathogens": predicted,
                "abstained_pathogens": abstained or [],
            },
        },
    }


class EvaluationTests(unittest.TestCase):
    def test_full_and_conditional_recall_are_separate(self):
        outputs = [
            output("1", ["A", "B", "C"], ["A", "C"], baseline=["A"]),
            output("2", ["X"], [], baseline=[]),
            output("3", ["X"], ["X"], baseline=[]),
        ]
        gold = {
            "schema_version": "rag_re.gold.v1",
            "label_semantics": "infection_source_positive_list",
            "complete_candidate_negative_labels": False,
            "cases": [
                {"patient_id": "1", "pathogens": ["A", "B"]},
                {"patient_id": "2", "pathogens": ["D"]},
                {"patient_id": "3", "pathogens": []},
            ],
            "aliases": {},
            "source_sha256": "test",
        }
        report = evaluate_outputs(outputs, gold, bootstrap_iterations=0)
        metrics = report["experiments"]["A_LITERATURE"]["metrics"]
        self.assertAlmostEqual(metrics["agreement_precision"], 1 / 3)
        self.assertAlmostEqual(metrics["end_to_end_recall"], 1 / 3)
        self.assertAlmostEqual(metrics["candidate_pool_recall_ceiling"], 2 / 3)
        self.assertAlmostEqual(metrics["conditional_classifier_recall"], 1 / 2)
        self.assertEqual(metrics["negative_case_specificity"], 0.0)
        self.assertAlmostEqual(metrics["false_positives_per_case"], 2 / 3)

    def test_missing_output_counts_as_false_negative(self):
        outputs = [output("1", ["A"], ["A"])]
        gold = {
            "schema_version": "rag_re.gold.v1",
            "label_semantics": "infection_source_positive_list",
            "complete_candidate_negative_labels": True,
            "cases": [
                {"patient_id": "1", "pathogens": ["A"]},
                {"patient_id": "2", "pathogens": ["B"]},
            ],
            "aliases": {},
            "source_sha256": "test",
        }
        report = evaluate_outputs(
            outputs,
            gold,
            bootstrap_iterations=0,
            allow_missing_outputs=True,
        )
        metrics = report["experiments"]["A_LITERATURE"]["metrics"]
        self.assertEqual(metrics["end_to_end_recall"], 0.5)
        self.assertEqual(metrics["candidate_pool_recall_ceiling"], 1.0)
        self.assertEqual(metrics["conditional_classifier_recall"], 1.0)
        self.assertEqual(metrics["output_case_coverage"], 0.5)
        self.assertEqual(report["missing_output_patient_ids"], ["2"])

    def test_confirmatory_rejects_missing_and_out_of_pool_predictions(self):
        outputs = [output("1", ["A"], ["B"])]
        gold = {
            "schema_version": "rag_re.gold.v1",
            "label_semantics": "infection_source_positive_list",
            "complete_candidate_negative_labels": True,
            "cases": [{"patient_id": "1", "pathogens": ["A"]}],
            "aliases": {},
            "source_sha256": "test",
        }
        with self.assertRaisesRegex(ValueError, "outside its frozen candidate pool"):
            evaluate_outputs(outputs, gold, bootstrap_iterations=0)

    def test_full_abc_rejects_skipped_a_and_bc_only_is_explicit(self):
        artifact = output("1", ["A"], [], baseline=[])
        artifact["execution"]["literature_attempted"] = False
        artifact["experiments"].update(
            {
                "B_MNGS": {
                    "semantics": "screening_trigger",
                    "expression": "B",
                    "predicted_pathogens": ["A"],
                    "abstained_pathogens": [],
                },
                "C_DIRECT_SUPPORT": {
                    "semantics": "screening_trigger",
                    "expression": "C",
                    "predicted_pathogens": [],
                    "abstained_pathogens": [],
                },
                "BC_OR": {
                    "semantics": "screening_trigger",
                    "expression": "B OR C",
                    "predicted_pathogens": ["A"],
                    "abstained_pathogens": [],
                },
                "E0_OR_C": {
                    "semantics": "baseline_preserving_cumulative",
                    "expression": "E0 OR C",
                    "predicted_pathogens": [],
                    "abstained_pathogens": [],
                },
            }
        )
        gold = {
            "schema_version": "rag_re.gold.v1",
            "label_semantics": "infection_source_positive_list",
            "complete_candidate_negative_labels": True,
            "cases": [{"patient_id": "1", "pathogens": ["A"]}],
            "aliases": {},
            "source_sha256": "test",
        }
        with self.assertRaisesRegex(ValueError, "requires Module A"):
            evaluate_outputs([artifact], gold, bootstrap_iterations=0)
        report = evaluate_outputs(
            [artifact], gold, bootstrap_iterations=0, analysis_mode="bc_only"
        )
        self.assertEqual(report["analysis_mode"], "bc_only")
        self.assertNotIn("A_LITERATURE", report["experiments"])
        self.assertEqual(report["experiments"]["B_MNGS"]["metrics"]["tp"], 1)

    def test_cumulative_e0_sequence_has_paired_delta_rows(self):
        artifact = output("1", ["A", "B"], ["A"], baseline=[])
        for name, predicted in (
            ("E0_OR_A", ["A"]),
            ("E0_OR_A_OR_B", ["A", "B"]),
            ("E0_OR_A_OR_B_OR_C", ["A", "B"]),
        ):
            artifact["experiments"][name] = {
                "semantics": "baseline_preserving_cumulative",
                "expression": name,
                "predicted_pathogens": predicted,
                "abstained_pathogens": [],
            }
        gold = {
            "schema_version": "rag_re.gold.v1",
            "label_semantics": "infection_source_positive_list",
            "complete_candidate_negative_labels": True,
            "cases": [{"patient_id": "1", "pathogens": ["A", "B"]}],
            "aliases": {},
            "source_sha256": "test",
        }
        report = evaluate_outputs([artifact], gold, bootstrap_iterations=10)
        self.assertIn("ADD_A_TO_E0", report["incremental_deltas"])
        self.assertIn("ADD_B_AFTER_E0_A", report["incremental_deltas"])
        self.assertIn("ADD_C_AFTER_E0_AB", report["incremental_deltas"])
        self.assertTrue(report["evaluator_fingerprint"])

    def test_precision_pruning_arms_have_paired_delta_rows(self):
        artifact = output("1", ["A", "B"], ["A"], baseline=["A", "B"])
        for name, predicted in (
            ("E0_AND_B", ["A"]),
            ("E0_AND_C", []),
            ("E0_AND_B_OR_C", ["A"]),
            ("E0_AND_B_AND_C", []),
        ):
            artifact["experiments"][name] = {
                "semantics": "precision_pruning_policy",
                "expression": name,
                "predicted_pathogens": predicted,
                "abstained_pathogens": [],
            }
        gold = {
            "schema_version": "rag_re.gold.v1",
            "label_semantics": "infection_source_positive_list",
            "complete_candidate_negative_labels": True,
            "cases": [{"patient_id": "1", "pathogens": ["A"]}],
            "aliases": {},
            "source_sha256": "test",
        }
        report = evaluate_outputs([artifact], gold, bootstrap_iterations=10)
        self.assertIn("PRUNE_WITH_B_VS_E0", report["incremental_deltas"])
        self.assertIn("PRUNE_WITH_C_VS_E0", report["incremental_deltas"])
        self.assertIn("PRUNE_WITH_B_OR_C_VS_E0", report["incremental_deltas"])
        self.assertIn("PRUNE_WITH_B_AND_C_VS_E0", report["incremental_deltas"])


if __name__ == "__main__":
    unittest.main()

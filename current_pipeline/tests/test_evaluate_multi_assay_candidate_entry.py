import unittest

from tools.evaluate_multi_assay_candidate_entry import evaluate_set


class TestCandidateMetrics(unittest.TestCase):
    def test_precision_uses_only_labeled_patients_and_recall_uses_all_answers(self):
        answers = {1: ["Klebsiella pneumoniae", "Pneumocystis jirovecii"]}
        predictions = {
            1: ["Klebsiella pneumoniae", "Candida albicans"],
            2: ["Pneumocystis jirovecii"],
        }
        result = evaluate_set(predictions, answers)
        self.assertEqual(result["matched"], 1)
        self.assertEqual(result["predicted_in_labeled_patients"], 2)
        self.assertEqual(result["answer_count"], 2)
        self.assertEqual(result["precision"], 0.5)
        self.assertEqual(result["recall"], 0.5)
        self.assertEqual(result["total_candidates_all_patients"], 3)


if __name__ == "__main__":
    unittest.main()

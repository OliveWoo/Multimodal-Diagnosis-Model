import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from generate_r5_final_decisions import compose_decision, evaluate_r5


def signal(**overrides):
    value = {
        "candidate_organism": "Organism B",
        "canonical_organism": "Organism B",
        "review_tier": "review_context_needed",
        "A1_strong": 0,
        "A1_partial": 0,
        "A1_match": 0,
        "B_strict": False,
        "B_absolute_axis": False,
        "B_relative_axis": False,
        "B_site_aligned": False,
        "C_grade": "C0",
        "D_state": "D_PASS",
        "clinical_features_source": "test",
        "casefit_file": "test.json",
        "casefit_sha256": "x",
        "latest_merge_file": "merge.json",
        "latest_merge_sha256": "y",
    }
    value.update(overrides)
    return value


class R5FinalizerTests(unittest.TestCase):
    def test_c3_direct_path(self):
        result = evaluate_r5(signal(C_grade="C3"))
        self.assertTrue(result["selected"])
        self.assertEqual(result["rule_path"], "C3_direct")

    def test_strong_site_absolute_path(self):
        result = evaluate_r5(
            signal(A1_strong=2, A1_match=2, B_site_aligned=True, B_absolute_axis=True)
        )
        self.assertTrue(result["selected"])
        self.assertEqual(result["rule_path"], "A1_site_B")

    def test_match_site_relative_path(self):
        result = evaluate_r5(
            signal(A1_partial=7, A1_match=7, B_site_aligned=True, B_relative_axis=True)
        )
        self.assertTrue(result["selected"])

    def test_site_is_required_for_non_c3(self):
        result = evaluate_r5(
            signal(A1_strong=2, A1_match=2, B_site_aligned=None, B_absolute_axis=True)
        )
        self.assertFalse(result["selected"])

    def test_picked_is_immutable_and_rescue_is_appended(self):
        rescued = signal(
            A1_strong=2,
            A1_match=2,
            B_site_aligned=True,
            B_relative_axis=True,
        )
        decision = compose_decision(
            cohort="test",
            patient_id="1",
            candidate_names=["Organism A", "Organism B", "Organism C"],
            picked=[{"canonical_name": "Organism A", "rank": 1}],
            signals=[rescued],
            source_snapshot={},
        )
        self.assertEqual(
            [(item["organism"], item["rank"]) for item in decision["selected"]],
            [("Organism A", 1), ("Organism B", 2)],
        )
        self.assertEqual(decision["counts"]["R5_rescued"], 1)
        self.assertFalse(decision["validation"]["gold_answer_used_in_this_run"])

    def test_candidate_missing_casefit_stays_unselected(self):
        decision = compose_decision(
            cohort="test",
            patient_id="1",
            candidate_names=["Organism A"],
            picked=[],
            signals=[],
            source_snapshot={},
        )
        self.assertEqual(decision["selected"], [])
        self.assertEqual(decision["candidates"][0]["decision"], "not_selected")


if __name__ == "__main__":
    unittest.main()

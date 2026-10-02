from __future__ import annotations

import unittest

from tools.freeze_kh_ablation_a import answer_names, f1, patient_sort_key


class FreezeKhAblationATests(unittest.TestCase):
    def test_answer_names_handles_empty_and_multiple_answers(self) -> None:
        self.assertEqual(answer_names("-"), [])
        self.assertEqual(answer_names(""), [])
        self.assertEqual(answer_names("A; B ;C"), ["A", "B", "C"])

    def test_f1(self) -> None:
        self.assertAlmostEqual(f1(48 / 62, 48 / 60), 0.7868852459016393)
        self.assertEqual(f1(0.0, 0.0), 0.0)

    def test_patient_sort_key_is_numeric(self) -> None:
        self.assertEqual(sorted(["10", "3", "9"], key=patient_sort_key), ["3", "9", "10"])


if __name__ == "__main__":
    unittest.main()

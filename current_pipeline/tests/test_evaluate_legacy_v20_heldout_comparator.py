from __future__ import annotations

import unittest

from tools.evaluate_legacy_v20_heldout_comparator import (
    build_pulmonary_minimum_sensitivity,
    level_rank,
)


class LegacyHeldoutComparatorTest(unittest.TestCase):
    def test_level_rank(self) -> None:
        self.assertEqual(level_rank("Level 1"), 1)
        self.assertEqual(level_rank("Level 3"), 3)
        self.assertEqual(level_rank(""), 99)

    def test_pulmonary_sensitivity_removes_only_explicit_urine_answer(self) -> None:
        labels = {
            31: {
                "status": "positive",
                "answers": ["Enterococcus faecium"],
                "row": {"raw_workbook_values": "Enterococcus faecium (VRE)-urine"},
            },
            42: {
                "status": "positive",
                "answers": ["Enterococcus faecium"],
                "row": {"raw_workbook_values": "Enterococcus faecium (VRE)_plasma"},
            },
        }
        sensitivity, changes = build_pulmonary_minimum_sensitivity(labels)
        self.assertEqual(sensitivity[31]["status"], "explicit_no_pathogen")
        self.assertEqual(sensitivity[31]["answers"], [])
        self.assertEqual(sensitivity[42]["answers"], ["Enterococcus faecium"])
        self.assertEqual(labels[31]["answers"], ["Enterococcus faecium"])
        self.assertEqual(len(changes), 1)


if __name__ == "__main__":
    unittest.main()

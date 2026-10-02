from __future__ import annotations

import unittest

from tools.evaluate_test_aware_deterministic_shadow import merge_groups


class TestEvaluateTestAwareDeterministicShadow(unittest.TestCase):
    def test_merge_groups_deduplicates_patient_organisms(self):
        merged = merge_groups(
            {8: ["Candida albicans", "Staphylococcus aureus"]},
            {8: ["Candida albicans"], 9: ["Pneumocystis jirovecii"]},
        )
        self.assertEqual(2, len(merged[8]))
        self.assertEqual(["Pneumocystis jirovecii"], merged[9])


if __name__ == "__main__":
    unittest.main()

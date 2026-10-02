import unittest

from tools.evaluate_legacy_multi_assay_compatibility import set_changes


class TestLegacyMultiAssayCompatibilityEvaluation(unittest.TestCase):
    def test_set_changes_uses_approved_alias_matching(self):
        old = {1: ["HSV-1", "Candida albicans"]}
        new = {1: ["Human alphaherpesvirus 1", "Pseudomonas aeruginosa"]}
        self.assertEqual(set_changes(old, new), [
            {"patient_id": 1, "change": "added_by_multi_assay", "organism_name": "Pseudomonas aeruginosa"},
            {"patient_id": 1, "change": "removed_or_replaced", "organism_name": "Candida albicans"},
        ])


if __name__ == "__main__":
    unittest.main()

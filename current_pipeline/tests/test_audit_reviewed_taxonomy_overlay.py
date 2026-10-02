import unittest

from tools.audit_reviewed_taxonomy_overlay import _is_reviewed_candidate


class TestAuditReviewedTaxonomyOverlay(unittest.TestCase):
    def test_reviewed_candidate_matches_canonical_alias(self):
        row = {
            "organism_name": "Metamycoplasma hominis",
            "organism_key": "metamycoplasmahominis",
        }

        self.assertTrue(_is_reviewed_candidate(row, {"mycoplasmahominis"}))

    def test_unrelated_candidate_is_not_marked_reviewed(self):
        row = {
            "organism_name": "Staphylococcus aureus",
            "organism_key": "staphylococcusaureus",
        }

        self.assertFalse(_is_reviewed_candidate(row, {"mycoplasmahominis"}))


if __name__ == "__main__":
    unittest.main()

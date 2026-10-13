import tempfile
import unittest
from pathlib import Path

from tools.audit_reviewed_taxonomy_overlay import (
    _candidate_identity,
    _is_reviewed_candidate,
    _reset_owned_shadow_dir,
)


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

    def test_missing_case_review_id_serializations_have_same_identity(self):
        base = {"patient_id": "107", "organism_key": "pneumocystisjirovecii"}

        self.assertEqual(
            _candidate_identity({**base, "case_review_id": None}),
            _candidate_identity({**base, "case_review_id": "None"}),
        )
        self.assertEqual(
            _candidate_identity({**base, "case_review_id": "None"}),
            _candidate_identity({**base, "case_review_id": ""}),
        )

    def test_distinct_case_review_ids_remain_distinct(self):
        base = {"patient_id": "107", "organism_key": "pneumocystisjirovecii"}

        self.assertNotEqual(
            _candidate_identity({**base, "case_review_id": "CASE-A"}),
            _candidate_identity({**base, "case_review_id": "CASE-B"}),
        )

    def test_reset_owned_shadow_dir_removes_only_shadow_child(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "audit"
            shadow = output / "shadow"
            sibling = output / "keep.txt"
            shadow.mkdir(parents=True)
            (shadow / "stale.txt").write_text("stale", encoding="utf-8")
            sibling.write_text("keep", encoding="utf-8")

            reset_path = _reset_owned_shadow_dir(output)

            self.assertEqual(reset_path, shadow.resolve())
            self.assertFalse(shadow.exists())
            self.assertTrue(sibling.exists())


if __name__ == "__main__":
    unittest.main()

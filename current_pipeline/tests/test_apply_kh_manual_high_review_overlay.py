import unittest

from tools.apply_kh_manual_high_review_overlay import apply


class ManualHighReviewOverlayTest(unittest.TestCase):
    def test_manual_high_preserves_deterministic_fields_and_cannot_pick(self):
        rows = [{
            "patient_id": "6",
            "organism_name": "Acinetobacter ursingii",
            "post_promotion_clinical_decision": "review_context_needed",
            "final_reporting_tier": "Fallback-Possible",
        }]
        policy = {"review_items": [{
            "patient_id": "6",
            "organism_name": "Acinetobacter ursingii",
            "manual_review_tier": "High",
            "decision_source": "test",
            "reason": "review",
        }]}
        output, queue, summary = apply(rows, policy)
        self.assertEqual(output[0]["effective_review_tier"], "High")
        self.assertEqual(output[0]["post_promotion_clinical_decision"], "review_context_needed")
        self.assertFalse(output[0]["manual_review_can_create_picked"])
        self.assertEqual(len(queue), 1)
        self.assertTrue(summary["all_assertions_pass"])


if __name__ == "__main__":
    unittest.main()

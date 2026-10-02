import unittest

from tools.build_kh_taxonomy_completion_audit import _clinical_tier, _front_route


class TaxonomyCompletionAuditHelpersTest(unittest.TestCase):
    def test_route_labels(self):
        self.assertEqual(_front_route("forward_priority_review"), "Priority")
        self.assertEqual(_front_route("hold_not_forwarded"), "Hold")
        self.assertEqual(_clinical_tier("picked_shadow"), "Picked")
        self.assertEqual(_clinical_tier("high_priority_review"), "High")


if __name__ == "__main__":
    unittest.main()

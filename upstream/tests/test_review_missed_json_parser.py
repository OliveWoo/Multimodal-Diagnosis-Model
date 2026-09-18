from __future__ import annotations

import unittest

from tools import review_missed_mngs_candidates as review


class ReviewMissedJsonParserTests(unittest.TestCase):
    def test_repairs_trailing_commas_outside_strings(self) -> None:
        payload = review.parse_json_response(
            '{"items": [1, 2,], "nested": {"value": 3,}, '
            '"literal": "keep comma, }",}'
        )

        self.assertEqual([1, 2], payload["items"])
        self.assertEqual({"value": 3}, payload["nested"])
        self.assertEqual("keep comma, }", payload["literal"])

    def test_repairs_unterminated_container_with_trailing_comma(self) -> None:
        payload = review.parse_json_response('{"items": [1, 2],')

        self.assertEqual([1, 2], payload["items"])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.audit_kh_taxonomy_end_to_end import (
    DEFAULT_AUTO_DECISIONS,
    DEFAULT_COMPACT,
    DEFAULT_INTEGRATED,
    DEFAULT_SCORER,
    DEFAULT_UNION,
    audit,
)


class TaxonomyEndToEndAuditTests(unittest.TestCase):
    def test_frozen_chain(self) -> None:
        if not DEFAULT_UNION.exists():
            self.skipTest("Frozen KH taxonomy inputs are unavailable")
        with tempfile.TemporaryDirectory() as temporary:
            result = audit(
                DEFAULT_UNION,
                DEFAULT_COMPACT,
                DEFAULT_SCORER,
                DEFAULT_INTEGRATED,
                DEFAULT_AUTO_DECISIONS,
                Path(temporary) / "audit",
            )
        self.assertTrue(result["all_required_assertions_pass"])
        self.assertEqual(result["counts"]["case_organism_rows"], 782)
        self.assertEqual(result["counts"]["scorer_rows"], 444)
        self.assertEqual(result["counts"]["family_replay_drift_rows"], 0)
        self.assertEqual(result["counts"]["canonical_collision_count"], 0)


if __name__ == "__main__":
    unittest.main()

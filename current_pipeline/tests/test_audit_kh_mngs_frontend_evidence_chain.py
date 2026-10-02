from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.audit_kh_mngs_frontend_evidence_chain import (
    DEFAULT_COMPACT,
    DEFAULT_INVENTORY,
    DEFAULT_SCORER,
    as_bool,
    as_float,
    audit,
)


class MngsFrontendEvidenceChainAuditTests(unittest.TestCase):
    def test_value_parsers_preserve_missingness(self) -> None:
        self.assertIsNone(as_float(""))
        self.assertEqual(as_float("12.5"), 12.5)
        self.assertTrue(as_bool("True"))
        self.assertFalse(as_bool("0"))

    def test_frozen_frontend_chain_passes_all_assertions(self) -> None:
        if not DEFAULT_INVENTORY.is_dir():
            self.skipTest("Frozen KH inventory is not available")
        with tempfile.TemporaryDirectory() as temp_dir:
            result = audit(
                DEFAULT_INVENTORY,
                DEFAULT_COMPACT,
                DEFAULT_SCORER,
                Path(temp_dir) / "audit",
            )
        self.assertTrue(result["all_assertions_pass"])
        self.assertEqual(result["assertion_pass_count"], 27)
        self.assertEqual(result["counts"]["tests"], 117)
        self.assertEqual(result["counts"]["scorer_entries"], 428)
        self.assertEqual(result["counts"]["scorer_selected_test_signals"], 856)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.analyze_kh_technical_repeat_value import (
    DEFAULT_ANSWERS,
    DEFAULT_INTEGRATED,
    DEFAULT_REPEAT_AUDIT,
    DEFAULT_SCORER,
    audit,
)


class TechnicalRepeatValueAuditTests(unittest.TestCase):
    def test_frozen_value_audit(self) -> None:
        if not DEFAULT_REPEAT_AUDIT.exists():
            self.skipTest("Frozen technical-repeat audit is not available")
        with tempfile.TemporaryDirectory() as temporary:
            result = audit(
                DEFAULT_REPEAT_AUDIT,
                DEFAULT_SCORER,
                DEFAULT_INTEGRATED,
                DEFAULT_ANSWERS,
                Path(temporary) / "audit",
            )
        self.assertTrue(result["all_assertions_pass"])
        self.assertEqual(result["counts"]["technical_repeat_only_entries"], 48)
        self.assertEqual(result["counts"]["scorer_route_dependent_entries"], 28)
        self.assertEqual(
            result["complete_report_variant_metrics"]
            ["audit_only_no_repeat_reporting_support"]["matched"],
            54,
        )


if __name__ == "__main__":
    unittest.main()

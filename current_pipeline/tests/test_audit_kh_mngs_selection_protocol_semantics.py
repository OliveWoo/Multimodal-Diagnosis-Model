from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.audit_kh_mngs_selection_protocol_semantics import (
    DEFAULT_COMPACT,
    DEFAULT_INVENTORY,
    DEFAULT_PATIENT_ROOT,
    DEFAULT_POLICY,
    DEFAULT_SCORER,
    audit,
    expected_priority,
    repeat_group_type,
)


class MngsSelectionProtocolSemanticsAuditTests(unittest.TestCase):
    def test_selection_priority_mapping(self) -> None:
        self.assertEqual(expected_priority([]), (None, "unranked"))
        self.assertEqual(expected_priority(["Code_NTC=00"]), (5, "NTC00"))
        self.assertEqual(expected_priority(["Code_RK_NTC=RK8"]), (4, "RK8"))
        self.assertEqual(expected_priority(["Code_RK_NTC=RK00"]), (3, "RK00"))
        self.assertEqual(
            expected_priority(["Code_NTC=00", "Code_RK_NTC=RK00"]),
            (1, "NTC00+RK00"),
        )

    def test_repeat_group_classification(self) -> None:
        paired = [
            {"condition": "-", "start_time": "T1"},
            {"condition": "12000g/10min，adapter 2倍稀釋", "start_time": "T1"},
        ]
        rerun = [
            {"condition": "原library重新上機", "start_time": "T1"},
            {"condition": "原library重新上機", "start_time": "T2"},
        ]
        self.assertEqual(
            repeat_group_type(paired), "same_run_parallel_12000g_and_unannotated"
        )
        self.assertEqual(repeat_group_type(rerun), "explicit_original_library_rerun")

    def test_frozen_semantics_audit_passes(self) -> None:
        if not DEFAULT_INVENTORY.is_dir():
            self.skipTest("Frozen KH inventory is not available")
        with tempfile.TemporaryDirectory() as temporary:
            result = audit(
                DEFAULT_INVENTORY,
                DEFAULT_COMPACT,
                DEFAULT_SCORER,
                DEFAULT_PATIENT_ROOT,
                DEFAULT_POLICY,
                Path(temporary) / "audit",
            )
        self.assertTrue(result["all_assertions_pass"])
        self.assertEqual(result["counts"]["tests"], 117)
        self.assertEqual(result["counts"]["technical_repeat_only_scorer_entries"], 48)
        self.assertEqual(result["counts"]["technical_repeat_scorer_route_changes"], 28)


if __name__ == "__main__":
    unittest.main()

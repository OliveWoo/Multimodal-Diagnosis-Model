from __future__ import annotations

import unittest

from tools import audit_evidence_preservation as audit


class AuditEvidencePreservationTests(unittest.TestCase):
    def test_shadow_suffix_is_applied_to_agent_patterns_only(self) -> None:
        self.assertEqual(
            ("*culture_agent_evidence_v2_full33_20260917.json",),
            audit.tagged_agent_patterns(
                ("*culture_agent.json",), "evidence_v2_full33_20260917"
            ),
        )
        self.assertEqual(
            ("*culture_agent.json",),
            audit.tagged_agent_patterns(("*culture_agent.json",), None),
        )

    def test_invalid_suffix_path_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            audit.tagged_agent_patterns(("*culture_agent.json",), "../baseline")


if __name__ == "__main__":
    unittest.main()

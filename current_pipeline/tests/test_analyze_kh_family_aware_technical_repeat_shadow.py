from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.analyze_kh_family_aware_technical_repeat_shadow import (
    DEFAULT_ANSWERS,
    DEFAULT_COMPACT,
    DEFAULT_INTEGRATED,
    DEFAULT_REPEAT_AUDIT,
    DEFAULT_SCORER,
    audit,
    exact_species_identity,
    family_route,
)


class FamilyAwareTechnicalRepeatShadowTests(unittest.TestCase):
    def test_family_routes_are_general_not_pjp_specific(self) -> None:
        self.assertEqual(
            family_route(
                family="high_consequence_opportunistic",
                history_route="Priority",
                all_rna=False,
            ),
            (True, "opportunistic_family_with_history_priority"),
        )
        self.assertFalse(
            family_route(
                family="high_consequence_opportunistic",
                history_route="Context",
                all_rna=False,
            )[0]
        )
        self.assertFalse(
            family_route(
                family="skin_airway_colonizer_prone",
                history_route="Priority",
                all_rna=False,
            )[0]
        )

    def test_species_level_identity_rejects_ambiguous_labels(self) -> None:
        self.assertTrue(exact_species_identity("Pneumocystis jirovecii"))
        self.assertFalse(exact_species_identity("Candida sp."))
        self.assertFalse(exact_species_identity("Acinetobacter group"))

    def test_frozen_shadow(self) -> None:
        if not DEFAULT_REPEAT_AUDIT.exists():
            self.skipTest("Frozen KH technical-repeat inputs are unavailable")
        with tempfile.TemporaryDirectory() as temporary:
            result = audit(
                DEFAULT_REPEAT_AUDIT,
                DEFAULT_COMPACT,
                DEFAULT_SCORER,
                DEFAULT_INTEGRATED,
                DEFAULT_ANSWERS,
                Path(temporary) / "shadow",
            )
        self.assertTrue(result["all_assertions_pass"])
        self.assertEqual(result["counts"]["technical_repeat_candidates"], 48)
        self.assertEqual(result["counts"]["changes_vs_current"], 0)
        self.assertEqual(
            result["scorer_metrics"]["family_aware_repeat_high"]["picked"],
            result["scorer_metrics"]["current_full"]["picked"],
        )


if __name__ == "__main__":
    unittest.main()

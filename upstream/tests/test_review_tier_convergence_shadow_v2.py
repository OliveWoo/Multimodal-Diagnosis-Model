from __future__ import annotations

import unittest

from tools import build_review_tier_convergence_shadow_v2 as shadow
from tools import compare_aliasclean_baseline as compare


class ReviewTierConvergenceShadowV2Test(unittest.TestCase):
    def test_distinct_candida_species_is_related_not_strict(self) -> None:
        answers = [compare.name_entry("Candida albicans")]
        strict, strict_answer, related, related_answer, relationship = shadow.strict_and_related_match(
            "Candida tropicalis", [answer for answer in answers if answer]
        )

        self.assertEqual("no", strict)
        self.assertEqual("", strict_answer)
        self.assertEqual("yes", related)
        self.assertEqual("Candida albicans", related_answer)
        self.assertEqual("distinct_species_same_genus", relationship)

    def test_alias_is_a_strict_match(self) -> None:
        answers = [compare.name_entry("COVID-19")]
        strict, strict_answer, related, _, relationship = shadow.strict_and_related_match(
            "SARS-CoV-2", [answer for answer in answers if answer]
        )

        self.assertEqual("yes", strict)
        self.assertEqual("COVID-19", strict_answer)
        self.assertEqual("no", related)
        self.assertEqual("exact_or_alias", relationship)

    def test_stronger_distinct_species_picked_requires_direct_support(self) -> None:
        picked = [
            {
                "organism_name": "Serratia marcescens",
                "support_modules": ["culture", "filmarray_gmtest"],
            }
        ]

        representative, supported = shadow.stronger_same_genus_picked("Serratia nevei", picked)

        self.assertEqual("Serratia marcescens", representative)
        self.assertEqual("yes", supported)


if __name__ == "__main__":
    unittest.main()

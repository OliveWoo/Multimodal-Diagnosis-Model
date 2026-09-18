from __future__ import annotations

import unittest

from tools import compare_aliasclean_baseline as compare
from tools import deterministic_mngs_max_scorer as scorer
from tools import merge_deterministic_max_with_missed_review as merger
from tools import pathogen_normalization as names


def candidate(
    organism_name: str,
    *,
    level: str,
    evidence_source: str,
    reads: int = 0,
    rank: str = "1",
    support: dict[str, str] | None = None,
    related: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    item: dict[str, object] = {
        "organism_name": organism_name,
        "integrated_causative_level": level,
        "evidence_source": evidence_source,
        "reads": reads,
        "reads_tier": "R2_medium" if reads else "Not_available",
        "reads_percentile": 0.8 if reads else 0,
        "dominance_tier": "D0_not_top",
        "rank_priority": rank,
        "module_support_summary": support or {},
        "key_evidence": {"support_modules": []},
    }
    if related:
        item["related_representative_hospital_support"] = related
    return item


class PathogenIdentityAndRepresentativeTests(unittest.TestCase):
    def test_approved_group_member_matches_are_formal_but_not_aliases(self) -> None:
        self.assertNotEqual(
            names.canonical_key("Klebsiella pneumoniae group"),
            names.canonical_key("Klebsiella pneumoniae"),
        )
        self.assertEqual(
            "approved_group_member_match",
            names.name_match_type("Klebsiella pneumoniae group", "Klebsiella pneumoniae"),
        )
        self.assertEqual(
            "approved_group_member_match",
            names.name_match_type(
                "Acinetobacter calcoaceticus-baumannii complex",
                "Acinetobacter baumannii",
            ),
        )
        self.assertNotEqual(
            names.canonical_key("Acinetobacter calcoaceticus-baumannii complex"),
            names.canonical_key("Acinetobacter baumannii"),
        )
        self.assertEqual(
            "approved_group_member_match",
            names.name_match_type("Enterobacter cloacae complex", "Enterobacter hormaechei"),
        )
        self.assertEqual(
            "",
            names.name_match_type(
                "Enterobacter cloacae complex",
                "Enterobacter quasihormaechei",
                allow_genus_relaxed=False,
            ),
        )

    def test_distinct_species_and_generic_aspergillus_are_not_formal_matches(self) -> None:
        self.assertEqual(
            "",
            names.name_match_type(
                "Klebsiella variicola",
                "Klebsiella pneumoniae",
                allow_genus_relaxed=False,
            ),
        )
        self.assertEqual(
            "",
            names.name_match_type(
                "Candida albicans",
                "Candida tropicalis",
                allow_genus_relaxed=False,
            ),
        )
        self.assertEqual(
            "",
            names.name_match_type(
                "Aspergillus spp.",
                "Aspergillus flavus",
                allow_genus_relaxed=False,
            ),
        )

    def test_pairwise_comparison_prefers_exact_over_group_and_group_over_genus(self) -> None:
        outputs = [
            compare.name_entry("Klebsiella pneumoniae group"),
            compare.name_entry("Klebsiella variicola"),
        ]
        answers = [
            compare.name_entry("Klebsiella pneumoniae"),
            compare.name_entry("Klebsiella variicola"),
        ]
        result = compare.pairwise_compare(
            [item for item in outputs if item],
            [item for item in answers if item],
        )
        self.assertEqual(
            ["approved_group_member_match", "exact_or_alias"],
            [item["match_type"] for item in result["matched"]],
        )

    def test_species_mngs_plus_genus_hospital_evidence_beats_generic_label(self) -> None:
        generic = candidate(
            "Aspergillus spp.",
            level="Level 1",
            evidence_source="hospital_only",
            support={"filmarray_gmtest": "Support"},
        )
        flavus = candidate(
            "Aspergillus flavus",
            level="Level 2",
            evidence_source="mNGS_ranked",
            reads=224,
            related=[{"organism_name": "Aspergillus spp."}],
        )
        self.assertIs(
            flavus,
            scorer.preferred_representative_for_group(
                "genus:aspergillus", [generic, flavus]
            ),
        )

    def test_exact_species_hospital_and_mngs_beats_generic_label(self) -> None:
        generic = candidate(
            "Aspergillus spp.",
            level="Level 1",
            evidence_source="hospital_only",
            support={"filmarray_gmtest": "Support"},
        )
        terreus = candidate(
            "Aspergillus terreus",
            level="Level 2",
            evidence_source="mNGS_ranked",
            reads=450,
            support={"culture": "Support"},
        )
        self.assertIs(
            terreus,
            scorer.preferred_representative_for_group(
                "genus:aspergillus", [generic, terreus]
            ),
        )

    def test_review_group_uses_approved_detected_member_as_representative(self) -> None:
        review = {
            "review_high_priority": [
                {
                    "organism_name": "Enterobacter cloacae complex",
                    "review_tier": "review_high_priority",
                    "evidence_snapshot": {
                        "reads": 0,
                        "support_modules": ["filmarray_gmtest"],
                    },
                }
            ],
            "review_context_needed": [],
            "review_low_specificity": [
                {
                    "organism_name": "Enterobacter hormaechei",
                    "evidence_snapshot": {
                        "reads": 10,
                        "reads_tier": "R0_trace",
                        "rank_priority": "3",
                    },
                },
                {
                    "organism_name": "Enterobacter quasihormaechei",
                    "evidence_snapshot": {
                        "reads": 82,
                        "reads_tier": "R1_low",
                        "rank_priority": "1",
                    },
                },
            ],
            "review_omitted_with_reason": [],
        }

        output = merger.apply_species_preferred_group_representatives(review)

        self.assertEqual(1, len(output["review_high_priority"]))
        retained = output["review_high_priority"][0]
        self.assertEqual("Enterobacter hormaechei", retained["organism_name"])
        self.assertEqual("Enterobacter cloacae complex", retained["representative_for_group"])
        self.assertEqual("review_high_priority", retained["review_tier"])
        self.assertTrue(retained["group_level_support_only"])
        self.assertEqual([], output["review_low_specificity"])
        self.assertEqual(
            {"Enterobacter cloacae complex", "Enterobacter quasihormaechei"},
            {item["original_name"] for item in retained["related_evidence"]},
        )

    def test_omitted_group_does_not_demote_visible_member(self) -> None:
        review = {
            "review_high_priority": [],
            "review_context_needed": [
                {
                    "organism_name": "Klebsiella pneumoniae",
                    "review_tier": "review_context_needed",
                    "evidence_snapshot": {"reads": 57, "reads_tier": "R1_low"},
                }
            ],
            "review_low_specificity": [],
            "review_omitted_with_reason": [
                {
                    "organism_name": "Klebsiella pneumoniae group",
                    "review_tier": "review_omitted_with_reason",
                    "evidence_snapshot": {"reads": 0},
                }
            ],
        }

        output = merger.apply_species_preferred_group_representatives(review)

        self.assertEqual(
            ["Klebsiella pneumoniae"],
            [item["organism_name"] for item in output["review_context_needed"]],
        )
        self.assertEqual(0, output["species_preferred_group_representative_policy"]["transformation_count"])


if __name__ == "__main__":
    unittest.main()

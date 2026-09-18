from __future__ import annotations

import unittest

from tools import merge_deterministic_max_with_missed_review as merge
from tools import organism_taxonomy_classifier as taxonomy


class OrganismTaxonomyClassifierTests(unittest.TestCase):
    def test_pjp_exact_profile_is_identity_not_causality(self) -> None:
        profile = taxonomy.classify_organism("PJP")
        self.assertEqual("pneumocystisjirovecii", profile["canonical_key"])
        self.assertEqual(42068, profile["taxid"])
        self.assertEqual("fungus", profile["biological_class"])
        self.assertEqual("high_consequence_opportunistic", profile["primary_rule_family"])
        self.assertEqual("exact_species", profile["mapping_status"])
        self.assertTrue(profile["routing_only"])
        self.assertFalse(profile["biological_class_conflict"])

    def test_candida_inherits_genus_profile(self) -> None:
        profile = taxonomy.classify_organism("Candida dubliniensis", biological_class="Fungal")
        self.assertEqual("candida_or_yeast", profile["primary_rule_family"])
        self.assertEqual("genus_inherited", profile["mapping_status"])
        self.assertIn("species_identity_matters", profile["clinical_traits"])

    def test_pseudomonas_has_multiple_clinical_roles(self) -> None:
        profile = taxonomy.classify_organism("Pseudomonas aeruginosa")
        self.assertEqual("hospital_or_nonfermenter_gnb", profile["primary_rule_family"])
        self.assertNotIn("environmental_low_specificity", profile["secondary_rule_families"])
        self.assertIn("environmental_association", profile["clinical_traits"])
        self.assertIn("colonization_possible", profile["clinical_traits"])

    def test_c_striatum_is_not_treated_as_impossible_pathogen(self) -> None:
        profile = taxonomy.classify_organism("Corynebacterium striatum")
        self.assertEqual("skin_airway_colonizer_prone", profile["primary_rule_family"])
        self.assertIn("pulmonary_pathogen_possible", profile["clinical_traits"])

    def test_unknown_species_uses_safe_fallback_and_keeps_supplied_class(self) -> None:
        profile = taxonomy.classify_organism("Examplebacter rareii", biological_class="Bacterial")
        self.assertEqual("unmapped_or_uncertain", profile["primary_rule_family"])
        self.assertEqual("bacterium", profile["biological_class"])
        self.assertEqual("low", profile["classification_confidence"])
        self.assertTrue(profile["needs_literature_review"])
        self.assertIn("do_not_auto_pick_from_taxonomy", profile["clinical_traits"])

    def test_rules_are_answer_blind(self) -> None:
        self.assertEqual([], taxonomy._find_forbidden_rule_keys(taxonomy.rules_payload()))

    def test_curated_identity_reports_conflicting_supplied_class(self) -> None:
        profile = taxonomy.classify_organism("PJP", biological_class="Bacterial")
        self.assertEqual("fungus", profile["biological_class"])
        self.assertTrue(profile["biological_class_conflict"])

    def test_merge_backfills_profiles_without_changing_tiers(self) -> None:
        deterministic = {
            "pathogen_candidates": [{"organism_name": "Pneumocystis jirovecii", "classification": "Fungal"}],
            "best_available_summary": {
                "picked_pathogens": [{"organism_name": "Staphylococcus aureus", "classification": "Bacterial"}]
            },
        }
        review = {
            "review_high_priority": [{"organism_name": "Candida albicans", "classification": "Fungal"}],
            "review_context_needed": [],
            "review_low_specificity": [],
            "review_omitted_with_reason": [],
        }
        merge.attach_taxonomy_profiles(deterministic, review)
        self.assertEqual(
            "high_consequence_opportunistic",
            deterministic["pathogen_candidates"][0]["taxonomy_profile"]["primary_rule_family"],
        )
        self.assertEqual(
            "candida_or_yeast",
            review["review_high_priority"][0]["taxonomy_profile"]["primary_rule_family"],
        )
        self.assertEqual(1, len(review["review_high_priority"]))


if __name__ == "__main__":
    unittest.main()

import unittest

from tools.organism_taxonomy_classifier import classify_organism, rules_payload


class TaxonomyGapResolutionTest(unittest.TestCase):
    def test_overlay_loaded(self):
        self.assertIn("rules/organism_taxonomy_gap_resolution_v1.json", rules_payload()["overlay_sources"])

    def test_torque_teno_is_commensal_virome(self):
        profile = classify_organism("Torque teno virus 29", biological_class="viral")
        self.assertEqual(profile["primary_rule_family"], "commensal_virome_or_endogenous_element")
        self.assertEqual(profile["mapping_status"], "key_prefix_inferred")

    def test_oral_streptococcus_is_not_generic_typical_pathogen(self):
        profile = classify_organism("Streptococcus salivarius", biological_class="bacterial")
        self.assertEqual(profile["primary_rule_family"], "oral_aspiration_or_anaerobe")
        self.assertEqual(profile["mapping_status"], "exact_species")

    def test_environmental_water_bacterium(self):
        profile = classify_organism("Polynucleobacter corsicus", biological_class="bacterial")
        self.assertEqual(profile["primary_rule_family"], "environmental_low_specificity")

    def test_skin_yeast(self):
        profile = classify_organism("Malassezia globosa", biological_class="fungal")
        self.assertEqual(profile["primary_rule_family"], "skin_airway_colonizer_prone")

    def test_opportunistic_black_yeast(self):
        profile = classify_organism("Exophiala oligosperma", biological_class="fungal")
        self.assertEqual(profile["primary_rule_family"], "mold_or_opportunistic_fungus")

    def test_environmental_trichoderma_species_overrides_genus(self):
        profile = classify_organism("Trichoderma atroviride", biological_class="fungal")
        self.assertEqual(profile["primary_rule_family"], "environmental_low_specificity")
        self.assertEqual(profile["mapping_status"], "exact_species")


if __name__ == "__main__":
    unittest.main()

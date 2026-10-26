import unittest
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

from tools.organism_taxonomy_classifier import classify_organism, rules_payload
from tools import organism_taxonomy_classifier as taxonomy


class TaxonomyGapResolutionTest(unittest.TestCase):
    def test_overlay_loaded(self):
        self.assertIn("rules/organism_taxonomy_gap_resolution_v1_20261020.json", rules_payload()["overlay_sources"])

    def test_polyomavirus_is_not_routed_as_herpesvirus(self):
        profile = classify_organism("Betapolyomavirus hominis", biological_class="viral")
        self.assertEqual(profile["primary_rule_family"], "polyomavirus_persistence_or_shedding")

    def test_reviewed_arthrobacter_identity_does_not_establish_causality(self):
        profile = classify_organism("Arthrobacter woluwensis", biological_class="bacterial")
        self.assertEqual(profile["primary_rule_family"], "environmental_low_specificity")
        self.assertEqual(profile["mapping_status"], "exact_species")

    def test_private_external_overlay_keeps_traceable_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            overlay = Path(temporary) / "synthetic_overlay.json"
            overlay.write_text(json.dumps({"exact_profiles": {}}), encoding="utf-8")
            try:
                with patch.object(taxonomy, "RULE_OVERLAY_PATHS", (overlay,)):
                    rules_payload.cache_clear()
                    self.assertIn(overlay.resolve().as_posix(), rules_payload()["overlay_sources"])
            finally:
                rules_payload.cache_clear()

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

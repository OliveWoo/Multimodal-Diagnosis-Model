from __future__ import annotations

import unittest

from tools import import_phenotype_context as importer


class PhenotypeContextImportTests(unittest.TestCase):
    def test_contract_defines_exactly_22_unique_phenotypes(self) -> None:
        self.assertEqual(22, len(importer.EXPECTED_PHENOTYPES))
        self.assertEqual(22, len(set(importer.EXPECTED_PHENOTYPES)))

    def test_passed_pre_pneumonia_yes_is_etiologic_prior(self) -> None:
        row = {
            "status": "YES",
            "confidence": "HIGH",
            "temporal_relation": "PRE_PNEUMONIA",
            "phenotype_qc_status": "PASS",
        }
        self.assertEqual("etiologic_prior", importer.phenotype_model_role(row))

    def test_after_onset_is_not_etiologic_prior(self) -> None:
        row = {
            "status": "YES",
            "confidence": "HIGH",
            "temporal_relation": "AFTER_ONSET",
            "phenotype_qc_status": "PASS",
        }
        self.assertEqual("outcome_or_context", importer.phenotype_model_role(row))

    def test_no_is_not_treated_as_strong_negative(self) -> None:
        row = {
            "status": "NO",
            "confidence": "LOW",
            "temporal_relation": "UNCLEAR",
            "phenotype_qc_status": "PASS",
        }
        self.assertEqual("not_established", importer.phenotype_model_role(row))

    def test_cbc_change_is_not_antimicrobial_response(self) -> None:
        context = importer.treatment_context(
            clinical_facts=[],
            clinical_timeline=[],
            labs=[{"test_name": "WBC", "value": 4.2, "date": "2026-01-01"}],
        )
        self.assertFalse(context["explicit_antimicrobial_response_available"])
        self.assertEqual(1, len(context["inflammatory_and_immune_lab_series"]))

    def test_explicit_drug_linked_improvement_is_retained(self) -> None:
        context = importer.treatment_context(
            clinical_facts=[
                {
                    "fact_type": "TREATMENT",
                    "description": "Clinical improvement after meropenem treatment",
                    "evidence": "Clinical improvement after meropenem treatment",
                }
            ],
            clinical_timeline=[],
            labs=[],
        )
        self.assertTrue(context["explicit_antimicrobial_response_available"])
        self.assertEqual(1, len(context["explicit_antimicrobial_response_mentions"]))

    def test_death_note_no_response_to_stimulation_is_not_treatment_response(self) -> None:
        text = "There was no light reflex nor other response to stimulation."
        self.assertFalse(importer.is_explicit_response(text))


if __name__ == "__main__":
    unittest.main()

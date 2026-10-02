from __future__ import annotations

import unittest

from tools import adjudicate_rag_v2_packets as adjudicate


def packet() -> dict:
    return {
        "case_id": "P30_Bacteroides_fragilis",
        "literature_evidence": {
            "knowledge_card_key": "adult_pulmonary:bacteroidesfragilis",
            "articles": [{"pmid": "123", "doi": "10.1/example"}],
        },
    }


def result() -> dict:
    return {
        "schema_version": "rag_adjudication_v2.0",
        "case_id": "P30_Bacteroides_fragilis",
        "knowledge_card": {
            "knowledge_card_key": "adult_pulmonary:bacteroidesfragilis",
            "pulmonary_pathogenicity": "recognized_but_uncommon",
            "respiratory_colonization_risk": "moderate",
            "typical_supporting_evidence": ["lung abscess"],
        },
        "retrieval_summary": {
            "search_completed": True,
            "sources": [{
                "citation_id": "PMID:123",
                "title": "Example",
                "year": 2020,
                "pmid_or_doi": "123",
                "study_type": "case_series",
                "quality": "moderate",
                "case_similarity": "moderate",
                "direction": "supports_pathogenicity",
                "key_point": "Example",
            }],
            "weighted_evidence": {
                "supports_pathogenicity": 1,
                "supports_colonization_or_background": 0,
                "uncertainty": 1,
            },
        },
        "patient_fit": {
            "supporting_local_evidence": ["lower respiratory"],
            "contradicting_local_evidence": ["D0"],
            "missing_key_evidence": ["culture"],
            "competing_explanation_strength": "moderate",
        },
        "decision": {
            "recommended_action": "context_only",
            "clinician_visibility": "show_context",
            "confidence": "low",
            "one_sentence_reason": "Evidence is incomplete.",
        },
        "audit": {"human_review_required": True, "limitations": ["Sparse literature"]},
    }


def test_aware_packet(
    *, name: str = "Acinetobacter baumannii", category: str = "bacterium"
) -> dict:
    value = packet()
    value["organism"] = {
        "display_name": name,
        "pathogen_category": category,
    }
    value["current_model_state"] = {
        "strict_competing_pathogens": [{
            "organism_name": "Staphylococcus aureus",
            "selection_role": "existing_picked_pathogen",
        }],
    }
    value["local_evidence"] = {
        "index_event": {
            "specimen_context": "lower_respiratory",
            "linked_events": [{"specimen_site": "BALF"}],
        },
        "mngs_per_test_evidence": {
            "best_rank_in_retained_universe": 2,
            "rank_band": "top_3",
            "selected_positive_test_count": 3,
            "dna_selected": True,
            "rna_selected": True,
            "cross_molecule_selected": True,
            "technical_repeat": True,
            "reproducibility_axis": True,
            "normalization_available_test_count": 2,
            "fully_evaluable_test_count": 1,
            "per_test_signals": [
                {
                    "seq_id": "DNA-1",
                    "nucleic_type": "DNA",
                    "reads": 206,
                    "rank_in_retained_universe": 4,
                    "rpm_total": 2.7,
                    "qc_status": "evaluable",
                    "normalization_available": True,
                },
                {
                    "seq_id": "DNA-2",
                    "nucleic_type": "DNA",
                    "reads": 1191,
                    "rank_in_retained_universe": 2,
                    "rpm_total": 31.5,
                    "qc_status": "evaluable",
                    "normalization_available": True,
                },
                {
                    "seq_id": "RNA-1",
                    "nucleic_type": "RNA",
                    "reads": 0,
                    "rank_in_retained_universe": 1,
                    "rpm_total": None,
                    "qc_status": "filtered",
                    "normalization_available": False,
                },
            ],
            "reads_sum_across_tests": None,
        },
        "same_organism_hospital_observations": [],
        "direct_evidence_timing": {
            "event_window_hours": 48,
            "event_aligned_positive_count": 0,
            "rows": [],
        },
        "host_context": {},
        "image_context": {},
        "phenotype_context": {},
    }
    return value


class AdjudicateRagV2Tests(unittest.TestCase):
    def test_valid_result_passes(self) -> None:
        adjudicate.validate_adjudication(result(), packet())

    def test_invented_citation_fails(self) -> None:
        value = result()
        value["retrieval_summary"]["sources"][0]["pmid_or_doi"] = "999"
        with self.assertRaises(ValueError):
            adjudicate.validate_adjudication(value, packet())

    def test_prefixed_composite_citation_passes(self) -> None:
        value = result()
        value["retrieval_summary"]["sources"][0]["pmid_or_doi"] = "PMID:123; DOI:10.1/example"
        adjudicate.validate_adjudication(value, packet())

    def test_case_id_mismatch_fails(self) -> None:
        value = result()
        value["case_id"] = "P99_X"
        with self.assertRaises(ValueError):
            adjudicate.validate_adjudication(value, packet())

    def test_test_aware_profile_uses_positive_tests_without_summing_reads(self) -> None:
        profile = adjudicate.test_aware_mngs_profile(test_aware_packet())
        self.assertEqual(profile["schema"], "per_test_v1")
        self.assertTrue(profile["lower_respiratory"])
        self.assertEqual(profile["best_rank_in_retained_universe"], 2)
        self.assertEqual(len(profile["per_test_signals"]), 2)
        self.assertIsNone(profile["reads_sum_across_tests"])

    def test_typical_reproducible_top_rank_signal_allows_context_review(self) -> None:
        boundary = adjudicate.patient_fit_boundary(test_aware_packet(), "test-aware-v2")
        self.assertTrue(boundary["minimum_patient_support_met"])
        self.assertTrue(boundary["clinician_visibility_allowed"])
        self.assertIn(
            "TOP3_REPRODUCIBLE_PER_TEST_LOWER_RESPIRATORY_SIGNAL",
            boundary["generalized_evidence_safeguards"],
        )
        self.assertFalse(boundary["strong_competing_pathogen_present"])
        self.assertEqual(
            boundary["competing_pathogens_with_strength_unknown"],
            ["Staphylococcus aureus"],
        )

    def test_event_aligned_filmarray_is_direct_lower_respiratory_molecular_evidence(self) -> None:
        value = test_aware_packet(name="Haemophilus influenzae")
        value["local_evidence"]["same_organism_hospital_observations"] = [{
            "observation_id": "ASSAY-004",
            "test_type": "filmarray",
            "panel_or_assay": "BioFire FilmArray Pneumonia Panel",
            "specimen_type": "BAL",
            "target": "Haemophilus influenzae",
            "raw_result": "Detected; 10^4 copy/mL",
            "detection_status": "detected",
        }]
        value["local_evidence"]["direct_evidence_timing"] = {
            "event_aligned_positive_count": 1,
            "rows": [{
                "observation_id": "ASSAY-004",
                "positive": True,
                "timing": "within_event_window",
                "specimen_type": "BAL",
                "raw_result": "Detected; 10^4 copy/mL",
            }],
        }
        exact = adjudicate.exact_species_evidence_profile(value)
        self.assertTrue(exact["direct_molecular"])
        self.assertTrue(exact["event_aligned_positive"])
        self.assertTrue(exact["event_aligned_lower_respiratory"])

    def test_event_aligned_blood_culture_is_invasive_evidence(self) -> None:
        value = test_aware_packet(name="Corynebacterium striatum", category="skin/airway flora")
        value["local_evidence"]["same_organism_hospital_observations"] = [{
            "observation_id": "CUL-002",
            "specimen_type": "Blood",
            "specimen_category": "Sterile_Site",
            "organism_name": "Corynebacterium striatum",
            "raw_result": "Isolated",
        }]
        value["local_evidence"]["direct_evidence_timing"] = {
            "event_aligned_positive_count": 1,
            "rows": [{
                "observation_id": "CUL-002",
                "positive": True,
                "timing": "within_event_window",
                "specimen_type": "Blood",
                "raw_result": "Isolated",
            }],
        }
        boundary = adjudicate.patient_fit_boundary(value, "test-aware-v2")
        self.assertIn(
            "EVENT_ALIGNED_EXACT_SPECIES_INVASIVE_EVIDENCE",
            boundary["generalized_evidence_safeguards"],
        )
        self.assertTrue(boundary["clinician_visibility_allowed"])

    def test_prompt_contains_derived_test_aware_summary_without_mutating_packet(self) -> None:
        value = test_aware_packet()
        prompt_text = adjudicate.prompt_for_packet(value, "test-aware-v2")
        self.assertIn('"derived_test_aware_evidence"', prompt_text)
        self.assertIn('"reads_sum_across_tests": null', prompt_text)
        self.assertNotIn("derived_test_aware_evidence", value)

    def test_haemophilus_influenzae_is_not_misclassified_as_influenza_virus(self) -> None:
        value = test_aware_packet(name="Haemophilus influenzae", category="bacterium")
        self.assertEqual(
            adjudicate.organism_category(value),
            "hospital_or_typical_bacterial_pathogen",
        )

    def test_staphylococcus_aureus_is_typical_bacterial_pathogen(self) -> None:
        value = test_aware_packet(name="Staphylococcus aureus", category="bacterium")
        self.assertEqual(
            adjudicate.organism_category(value),
            "hospital_or_typical_bacterial_pathogen",
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import unittest

from tools import normalized_agent_fallback as fallback
from tools import deterministic_mngs_max_scorer as scorer
from tools import pathogen_normalization as pathogen_names


class NormalizedAgentFallbackTests(unittest.TestCase):
    def test_culture_fallback_preserves_every_record_and_all_positive_links(self) -> None:
        payload = fallback.culture_agent_from_source(
            [
                {
                    "sample": "Lower BAL",
                    "collected_time": "2024-07-04 11:30",
                    "reported_time": "2024-07-06 08:00",
                    "organism": "Trichosporon asahii",
                    "status": "Isolated",
                    "colony_count": "50x10^3 CFU/mL",
                },
                {
                    "sample": "Lower BAL",
                    "collected_time": "2024-07-04 11:30",
                    "reported_time": "2024-07-06 08:00",
                    "organism": "Trichosporon asahii",
                    "status": "Isolated",
                    "colony_count": ">10^5 CFU/mL",
                },
                {
                    "sample": "Blood",
                    "collected_time": "2024-07-04 12:00",
                    "reported_time": "2024-07-09 08:00",
                    "status": "No growth",
                },
            ]
        )
        self.assertEqual(3, len(payload["culture_observations"]))
        self.assertTrue(
            payload["source_record_contract"]["one_observation_per_source_record"]
        )
        label = payload["organism_labels"][0]
        self.assertEqual(["CUL-001", "CUL-002"], label["evidence_observation_ids"])
        self.assertEqual("Q4", label["key_evidence"]["quantity_tier"])
        self.assertEqual(
            "no_growth", payload["culture_observations"][2]["growth_purity"]
        )

    def test_culture_fallback_preserves_strong_respiratory_bacteria_and_candida_caution(self) -> None:
        payload = fallback.culture_agent_from_source(
            [
                {
                    "sample": "Lower BAL",
                    "organism": "Candida albicans",
                    "status": "Isolated",
                    "colony_count": "10x10^3 CFU/mL",
                },
                {
                    "sample": "Lower BAL",
                    "organism": "Stenotrophomonas maltophilia",
                    "status": "Isolated",
                    "colony_count": ">10^5 CFU/mL",
                },
                {
                    "sample": "Sputum",
                    "organism": "Klebsiella pneumoniae",
                    "status": "Isolated",
                },
            ]
        )
        levels = {item["organism_name"]: item["causative_level"] for item in payload["organism_labels"]}
        self.assertEqual("Level 3", levels["Candida albicans"])
        self.assertEqual("Level 2", levels["Stenotrophomonas maltophilia"])
        self.assertEqual("Level 2", levels["Klebsiella pneumoniae"])
        self.assertIn("missing_quantity", payload["data_gaps"])

    def test_culture_fallback_recognizes_baumannii_complex_without_claiming_exact_species(self) -> None:
        payload = fallback.culture_agent_from_source(
            [
                {
                    "sample": "Endotracheal aspirate",
                    "organism": "Acinetobacter baumannii complex",
                    "status": "Isolated",
                }
            ]
        )
        label = payload["organism_labels"][0]
        self.assertEqual("Acinetobacter baumannii complex", label["organism_name"])
        self.assertEqual("Level 2", label["causative_level"])

    def test_exact_candida_from_blood_is_level1_without_cfu(self) -> None:
        payload = fallback.culture_agent_from_source(
            [
                {
                    "sample": "Blood",
                    "organism": "Candida tropicalis",
                    "status": "Isolated",
                },
                {
                    "sample": "Blood",
                    "organism": "Yeast",
                    "status": "Isolated",
                },
            ]
        )
        levels = {
            item["organism_name"]: item["causative_level"]
            for item in payload["organism_labels"]
        }
        self.assertEqual("Level 1", levels["Candida tropicalis"])
        self.assertEqual("Level 2", levels["Yeast"])

    def test_exact_colonizer_prone_species_from_blood_is_provisional_level3(self) -> None:
        payload = fallback.culture_agent_from_source(
            [
                {
                    "sample": "Blood",
                    "organism": "Corynebacterium striatum",
                    "status": "Isolated",
                }
            ]
        )
        label = payload["organism_labels"][0]
        self.assertEqual("Level 3", label["causative_level"])
        self.assertEqual("Unknown", label["key_evidence"]["contaminant_flag"])
        self.assertIn("R-S1-02A", label["applied_rules"])

    def test_filmarray_fallback_excludes_negative_and_amr_from_organism_labels(self) -> None:
        payload = fallback.filmarray_agent_from_source(
            [
                {
                    "sample": "Lower BAL",
                    "target": "Klebsiella pneumoniae group",
                    "result": "Detected",
                    "value": "10^5 copy/mL",
                },
                {
                    "sample": "Lower BAL",
                    "target": "Human Rhinovirus/Enterovirus",
                    "result": "Detected",
                },
                {"sample": "Lower BAL", "target": "CTX-M", "result": "Detected"},
                {
                    "sample": "Lower BAL",
                    "target": "Pseudomonas aeruginosa",
                    "result": "Not Detected",
                },
            ]
        )
        levels = {item["organism_name"]: item["causative_level"] for item in payload["organism_labels"]}
        self.assertEqual("Level 3", levels["Klebsiella pneumoniae group"])
        self.assertEqual("Level 2", levels["Human Rhinovirus/Enterovirus"])
        self.assertNotIn("CTX-M", levels)
        self.assertNotIn("Pseudomonas aeruginosa", levels)
        self.assertEqual(["CTX-M"], [item["gene"] for item in payload["resistance_findings"]])
        self.assertEqual(4, len(payload["assay_observations"]))
        self.assertTrue(
            payload["source_record_contract"]["one_observation_per_source_record"]
        )
        klebsiella = next(
            item for item in payload["organism_labels"]
            if item["organism_name"] == "Klebsiella pneumoniae group"
        )
        self.assertEqual(["ASSAY-001"], klebsiella["evidence_observation_ids"])
        self.assertEqual(
            ["ASSAY-003"],
            payload["resistance_findings"][0]["evidence_observation_ids"],
        )
        self.assertEqual(
            "not_detected", payload["assay_observations"][3]["detection_status"]
        )

    def test_combined_filmarray_gm_fallback_preserves_all_rows_and_gm_levels(self) -> None:
        payload = fallback.filmarray_gm_agent_from_sources(
            [{
                "sample": "Lower BAL",
                "target": "Parainfluenza Virus",
                "result": "Detected",
                "reported_time": "2024-02-26 13:42",
            }],
            [
                {
                    "sample": "Lower BAL",
                    "target": "GM test (Index)",
                    "value": "1.695",
                    "reported_time": "2024-02-26",
                },
                {
                    "sample": "Serum",
                    "target": "GM test (Index)",
                    "value": "0.111",
                    "reported_time": "2024-02-26",
                },
                {
                    "sample": "Urine",
                    "target": "L.pne Ag",
                    "result": "Negative",
                },
            ],
        )
        self.assertEqual(4, len(payload["assay_observations"]))
        self.assertTrue(
            payload["source_record_contract"]["one_observation_per_source_record"]
        )
        levels = {
            item["organism_name"]: item["causative_level"]
            for item in payload["organism_labels"]
        }
        # With only one qualitative target in this unit fixture, panel type is
        # intentionally unknown; the same target is Level 2 in a full
        # lower-respiratory Pneumonia Panel record containing semiquant bins.
        self.assertEqual("Level 3", levels["Parainfluenza Virus"])
        self.assertEqual("Level 1", levels["Aspergillus spp."])
        aspergillus = next(
            item for item in payload["organism_labels"]
            if item["organism_name"] == "Aspergillus spp."
        )
        self.assertEqual(["ASSAY-002"], aspergillus["evidence_observation_ids"])
        self.assertEqual("Both", payload["assay_summary"]["gm_test_sample_type"])
        self.assertNotIn("L.pne", levels)

    def test_combined_fallback_keeps_positive_non_gm_antigen_conservatively(self) -> None:
        payload = fallback.filmarray_gm_agent_from_sources(
            [],
            [{
                "sample": "Urine",
                "target": "Legionella pneumophila Ag",
                "result": "Positive",
            }],
        )
        label = payload["organism_labels"][0]
        self.assertEqual("Legionella pneumophila", label["organism_name"])
        self.assertEqual("Level 3", label["causative_level"])
        self.assertEqual("antigen", label["key_evidence"]["non_gm_test_type"])
        self.assertEqual(["ASSAY-001"], label["evidence_observation_ids"])

    def test_hospital_support_floor_does_not_downgrade_existing_m1_signal(self) -> None:
        raw_candidate = {
            "organism_name": "Klebsiella pneumoniae",
            "source_category": "1.Bac",
            "reads": 2116,
            "ranking": {
                "rank_priority": 1,
                "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                "reads_percentile": 1.0,
                "possibility_level": "high",
            },
        }
        record = {
            "specimen_site": "BALF",
            "candidates": [raw_candidate],
        }
        final_by_name = {
            pathogen_names.canonical_key("Klebsiella pneumoniae"): {
                "organism_name": "Klebsiella pneumoniae",
                "best_hospital_level": "Level 2",
                "module_level_summary": {
                    "culture": "Level 2",
                    "filmarray_gmtest": "Not_available",
                    "molecular_microbiology": "Not_available",
                    "image": "Not_available",
                },
                "module_evidence": {
                    "culture": [{"specimen_type": "Lower BAL"}],
                },
            }
        }
        candidate = scorer.score_candidate(
            raw_candidate,
            record,
            {pathogen_names.canonical_key("Klebsiella pneumoniae"): "D1_low"},
            final_by_name,
            {
                "host_vulnerability_tier": "Unknown",
                "opportunistic_coverage_level": "Unknown",
                "host_level3_expansion": False,
            },
        )
        self.assertEqual("M1_strong", candidate["mngs_signal_tier"])
        self.assertEqual("Level 1", candidate["integrated_causative_level"])

    def test_bcid_resistance_target_does_not_imply_pneumonia_panel(self) -> None:
        payload = fallback.filmarray_agent_from_source(
            [
                {
                    "sample": "Blood culture",
                    "target": "Enterococcus faecium",
                    "result": "Detected",
                },
                {"sample": "Blood culture", "target": "vanA", "result": "Detected"},
            ]
        )
        self.assertEqual("BioFire BCID2", payload["assay_summary"]["filmarray_panel"])
        self.assertEqual(
            "Level 2",
            next(
                item["causative_level"]
                for item in payload["organism_labels"]
                if item["organism_name"] == "Enterococcus faecium"
            ),
        )

    def test_molecular_fallback_preserves_every_row_and_only_labels_positives(self) -> None:
        payload = fallback.molecular_microbiology_agent_from_source(
            [
                {
                    "test": "Non-panel Molecular Test",
                    "sample": "BALF",
                    "target": "Pneumocystis jirovecii",
                    "result": "Detected",
                },
                {
                    "test": "Non-panel Molecular Test",
                    "sample": "BALF",
                    "target": "Pseudomonas aeruginosa",
                    "result": "Not Detected",
                },
                {
                    "test": "Non-panel Molecular Test",
                    "sample": "Blood",
                    "target": "CMV",
                    "result": "Invalid",
                },
                {
                    "test": "Resistance PCR",
                    "sample": "BALF",
                    "target": "blaKPC",
                    "result": "Detected",
                },
            ]
        )
        self.assertEqual(4, len(payload["molecular_observations"]))
        self.assertTrue(
            payload["source_record_contract"]["one_observation_per_source_record"]
        )
        self.assertEqual(
            ["Pneumocystis jirovecii"],
            [item["organism_name"] for item in payload["organism_labels"]],
        )
        self.assertEqual(
            ["detected_positive", "not_detected", "invalid", "detected_positive"],
            [
                item["result_interpretation"]
                for item in payload["molecular_observations"]
            ],
        )
        self.assertEqual(
            ["blaKPC"], [item["gene"] for item in payload["resistance_findings"]]
        )
        self.assertEqual(
            ["MOL-001"], payload["organism_labels"][0]["evidence_observation_ids"]
        )

    def test_molecular_fallback_uses_specimen_specific_levels(self) -> None:
        payload = fallback.molecular_microbiology_agent_from_source(
            [
                {"sample": "BALF", "target": "CMV", "result": "Detected"},
                {"sample": "Blood", "target": "EBV", "result": "Detected"},
                {
                    "sample": "CSF",
                    "target": "Herpes simplex virus 1",
                    "result": "Detected",
                },
                {
                    "sample": "NPS",
                    "target": "Influenza A virus",
                    "result": "Detected",
                },
            ]
        )
        levels = {
            item["organism_name"]: item["causative_level"]
            for item in payload["organism_labels"]
        }
        self.assertEqual("Level 3", levels["CMV"])
        self.assertEqual("Level 2", levels["EBV"])
        self.assertEqual("Level 2", levels["Herpes simplex virus 1"])
        self.assertEqual("Level 3", levels["Influenza A virus"])
        self.assertEqual(
            "CNS", payload["molecular_observations"][2]["specimen_category"]
        )


if __name__ == "__main__":
    unittest.main()

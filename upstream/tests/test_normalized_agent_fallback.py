from __future__ import annotations

import unittest

from tools import normalized_agent_fallback as fallback
from tools import deterministic_mngs_max_scorer as scorer
from tools import pathogen_normalization as pathogen_names


class NormalizedAgentFallbackTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()

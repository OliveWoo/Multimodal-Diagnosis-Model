import unittest
from unittest.mock import patch

from tools import organism_taxonomy_classifier as taxonomy
from tools.run_legacy_scorer_on_multi_assay_shadow import (
    DEFAULT_TAXONOMY_OVERLAY_PATHS,
    build_legacy_ranked_payload,
    configure_taxonomy_mode,
)


class TestLegacyMultiAssayAdapter(unittest.TestCase):
    def test_builds_one_record_per_test_without_summing_reads(self):
        packet = {
            "patient_id": "9",
            "cases": [{
                "case_review_id": "CASE-B", "specimen_code": "CASE-B",
                "specimen_site": "BALF", "collected_time": "2024-01-01",
                "scorer_entries": [
                    {"organism_name": "Pneumocystis jirovecii", "category": "fungal",
                     "rule_ids": ["PRIORITY_FAMILY_ENTRY"], "selected_test_signals": [
                         {"seq_id": "DNA-1", "nucleic_type": "DNA", "reads": 79,
                          "analytical_rank_priority": 1, "rank_in_retained_universe": 1},
                         {"seq_id": "RNA-1", "nucleic_type": "RNA", "reads": 63,
                          "analytical_rank_priority": 1, "rank_in_retained_universe": 2},
                     ]},
                    {"organism_name": "Candida albicans", "category": "fungal",
                     "rule_ids": [], "selected_test_signals": [
                         {"seq_id": "DNA-1", "nucleic_type": "DNA", "reads": 20,
                          "analytical_rank_priority": 3, "rank_in_retained_universe": 2},
                     ]},
                ],
            }],
        }
        ranked, provenance = build_legacy_ranked_payload(packet)
        self.assertEqual(len(ranked["records"]), 2)
        self.assertEqual(sum(len(row["candidates"]) for row in ranked["records"]), 3)
        dna = next(row for row in ranked["records"] if row["seq_id"] == "DNA-1")
        pjp = next(row for row in dna["candidates"] if row["organism_name"] == "Pneumocystis jirovecii")
        self.assertEqual(pjp["reads"], 79)
        self.assertEqual(pjp["ranking"]["reads_percentile"], 1.0)
        self.assertEqual(len(provenance["pneumocystisjirovecii"]), 2)
        self.assertTrue(ranked["ranking_metadata"]["no_cross_test_read_sum"])

    def test_rejects_zero_read_signal(self):
        packet = {"patient_id": "1", "cases": [{
            "case_review_id": "A", "scorer_entries": [{
                "organism_name": "X y", "category": "bacterial",
                "selected_test_signals": [{"seq_id": "S", "reads": 0}],
            }],
        }]}
        with self.assertRaises(ValueError):
            build_legacy_ranked_payload(packet)

    def test_base_only_taxonomy_ignores_later_overlays(self):
        try:
            base = configure_taxonomy_mode("base-only")
            self.assertEqual(base["loaded_overlay_paths"], [])
            self.assertEqual(taxonomy.RULE_OVERLAY_PATHS, ())
            with patch.object(taxonomy, "RULE_OVERLAY_PATHS", DEFAULT_TAXONOMY_OVERLAY_PATHS):
                taxonomy.rules_payload.cache_clear()
                current = taxonomy.rules_payload()
                self.assertEqual(
                    current.get("overlay_sources"),
                    [
                        path.relative_to(taxonomy.REPO_ROOT).as_posix()
                        for path in DEFAULT_TAXONOMY_OVERLAY_PATHS
                        if path.exists()
                    ],
                )
        finally:
            taxonomy.RULE_OVERLAY_PATHS = DEFAULT_TAXONOMY_OVERLAY_PATHS
            taxonomy.rules_payload.cache_clear()

    def test_unknown_taxonomy_mode_is_rejected(self):
        with self.assertRaises(ValueError):
            configure_taxonomy_mode("latest-ish")


if __name__ == "__main__":
    unittest.main()

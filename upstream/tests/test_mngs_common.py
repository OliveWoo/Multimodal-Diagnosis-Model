from __future__ import annotations

import unittest

from tools import mngs_common


class RankedDeduplicationTests(unittest.TestCase):
    def test_keeps_best_record_per_normalized_organism(self) -> None:
        payload = {
            "patient_id": "1",
            "records": [
                {
                    "specimen_code": "late",
                    "candidates": [
                        {
                            "organism_name": "Pseudomonas aeruginosa",
                            "reads": 10,
                            "ranking": {"rank_priority": 3},
                        }
                    ],
                },
                {
                    "specimen_code": "best",
                    "candidates": [
                        {
                            "organism_name": "Pseudomonas aeruginosa",
                            "reads": 100,
                            "ranking": {"rank_priority": 1},
                        }
                    ],
                },
            ],
        }

        result = mngs_common.deduplicate_ranked_mngs_payload(payload)

        self.assertEqual(len(result["records"]), 1)
        self.assertEqual(result["records"][0]["specimen_code"], "best")
        self.assertEqual(result["records"][0]["candidates"][0]["reads"], 100)
        self.assertEqual(
            result["ranking_metadata"]["deduplication"]["same_organism_rows_removed"],
            1,
        )

    def test_candidate_deduplication_prefers_stronger_level(self) -> None:
        candidates = [
            {
                "organism_name": "Candida albicans",
                "integrated_causative_level": "Level 4",
                "rank_priority": "1",
                "reads": 100,
            },
            {
                "organism_name": "Candida albicans",
                "integrated_causative_level": "Level 2",
                "rank_priority": "3",
                "reads": 10,
            },
        ]

        result = mngs_common._dedupe_candidates(candidates)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["integrated_causative_level"], "Level 2")


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import unittest

from tools import reconstruct_pipeline_run_manifest as provenance


class ReconstructedRunManifestTests(unittest.TestCase):
    def test_ranking_formula_prefers_explicit_field(self) -> None:
        payload = {
            "ranking_metadata": {
                "score_formula": "final_score = code_score * 0.7 + reads_percentile * 0.3",
                "notes": ["ignored"],
            }
        }
        self.assertEqual(
            "final_score = code_score * 0.7 + reads_percentile * 0.3",
            provenance.ranking_formula(payload),
        )

    def test_ranking_formula_reads_historical_note(self) -> None:
        payload = {
            "ranking_metadata": {
                "notes": ["final_score = 0.60 * code_score + 0.40 * reads_percentile"]
            }
        }
        self.assertEqual(
            "final_score = 0.60 * code_score + 0.40 * reads_percentile",
            provenance.ranking_formula(payload),
        )

    def test_ranking_formula_is_inferred_from_all_stored_scores(self) -> None:
        payload = {
            "patients": {
                "1": {
                    "records": [
                        {
                            "candidates": [
                                {"ranking": {"code_score": 0.65, "reads_percentile": 1.0, "final_score": 0.755}},
                                {"ranking": {"code_score": 0.65, "reads_percentile": 0.7222, "final_score": 0.6717}},
                                {"ranking": {"code_score": 0.35, "reads_percentile": 0.8333, "final_score": 0.495}},
                            ]
                        }
                    ]
                }
            }
        }
        observation = provenance.ranking_formula_observation(payload)
        self.assertEqual(
            "final_score = 0.70 * code_score + 0.30 * reads_percentile",
            observation["formula"],
        )
        self.assertEqual("inferred_from_stored_candidate_scores", observation["evidence_source"])
        self.assertEqual(3, observation["informative_candidate_count"])

    def test_review_settings_are_read_from_saved_request(self) -> None:
        payload = {
            "llm_request": {
                "model": "gpt-5.6-luna",
                "reasoning_effort": "medium",
                "temperature_requested": None,
                "temperature_sent": None,
            }
        }
        self.assertEqual("gpt-5.6-luna", provenance.review_settings(payload)["model"])

    def test_merge_semantics_do_not_claim_review_overrides_picked(self) -> None:
        merged = {"final_name_reconciliation": {"enabled": True}}
        review = {
            "syndrome_pattern_convergence_policy": {
                "version": "v1",
                "promotions_enabled": False,
            }
        }
        result = provenance.merge_semantics(merged, review)
        self.assertFalse(result["deterministic_picked_advisory_override"])
        self.assertTrue(result["name_reconciliation_enabled"])


if __name__ == "__main__":
    unittest.main()

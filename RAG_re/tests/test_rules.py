import copy
import unittest

from rag_re.config import DEFAULT_CONFIG
from rag_re.rules import (
    experiment_values,
    module_a_from_judgments,
    module_b_mngs,
    module_c_direct_support,
    tri_and,
    tri_or,
)


class RuleTests(unittest.TestCase):
    def setUp(self):
        self.config = copy.deepcopy(DEFAULT_CONFIG)

    @staticmethod
    def judgment(pmid: str, verdict: str = "support", **overrides):
        row = {
            "pmid": pmid,
            "status": "ok",
            "verdict": verdict,
            "organism_scope": "exact_species",
            "human_clinical_evidence": True,
            "target_site_match": True,
            "evidence_span_valid": True,
        }
        row.update(overrides)
        return row

    def test_module_a_threshold_boundary_and_validation(self):
        retrieval = {"status": "ok"}
        judgments = [self.judgment(str(index)) for index in range(3)] + [
            self.judgment("4", "unclear"),
            self.judgment("5", "irrelevant"),
        ]
        result = module_a_from_judgments(
            retrieval, judgments, self.config["module_a_literature"]
        )
        self.assertIs(result["positive"], True)
        self.assertEqual(result["support_count"], 3)
        self.assertEqual(result["judgeable_count"], 5)

        invalid = list(judgments)
        invalid[0] = self.judgment("0", organism_scope="genus_only")
        result = module_a_from_judgments(
            retrieval, invalid, self.config["module_a_literature"]
        )
        self.assertIs(result["positive"], False)
        self.assertEqual(result["support_count"], 2)

    def test_module_a_insufficient_is_abstention(self):
        result = module_a_from_judgments(
            {"status": "ok"},
            [self.judgment("1"), self.judgment("2")],
            self.config["module_a_literature"],
        )
        self.assertIsNone(result["positive"])
        self.assertEqual(result["status"], "insufficient_articles")

    def test_api_failure_is_not_zero_support(self):
        result = module_a_from_judgments(
            {"status": "error"}, [], self.config["module_a_literature"]
        )
        self.assertIsNone(result["positive"])
        self.assertEqual(result["status"], "error")

    def test_module_b_uses_only_structured_mngs_features(self):
        candidate = {
            "mngs_signal_tier": "M2_moderate",
            "rank_priority": "1",
            "reads_percentile": 0.9,
            "specimen_alignment": "Aligned",
        }
        result = module_b_mngs(candidate, self.config["module_b_mngs"])
        self.assertIs(result["positive"], True)

        candidate.pop("specimen_alignment")
        candidate.pop("specimen_class", None)
        result = module_b_mngs(candidate, self.config["module_b_mngs"])
        self.assertIsNone(result["positive"])
        self.assertEqual(result["status"], "insufficient_input")
        candidate["specimen_alignment"] = "Not_aligned"
        result = module_b_mngs(candidate, self.config["module_b_mngs"])
        self.assertIs(result["positive"], False)

        candidate["specimen_alignment"] = "Sterile_or_Systemic"
        candidate["specimen_class"] = "S1_sterile_systemic"
        result = module_b_mngs(candidate, self.config["module_b_mngs"])
        self.assertIs(result["positive"], True)

    def test_module_c_excludes_host_and_mngs(self):
        candidate = {
            "specimen_alignment": "Aligned",
            "module_support_summary": {
                "mNGS": "M1_strong",
                "host": "Support",
                "culture": "Not_available",
                "filmarray_gmtest": "Support",
            }
        }
        result = module_c_direct_support(
            candidate, self.config["module_c_direct_support"]
        )
        self.assertIs(result["positive"], True)
        self.assertEqual([row["module"] for row in result["supports"]], ["filmarray_gmtest"])

    def test_module_c_hospital_only_uses_direct_evidence_site_and_rejects_pending(self):
        candidate = {
            "specimen_alignment": "Hospital_only",
            "specimen_class": "S3_unknown_or_low_value",
            "module_support_summary": {"culture": "Support"},
            "key_evidence": {
                "hospital_module_evidence": {
                    "culture": [
                        {
                            "specimen_type": "Sputum",
                            "specimen_category": "Lower_Respiratory",
                        }
                    ]
                }
            },
        }
        result = module_c_direct_support(
            candidate,
            self.config["module_c_direct_support"],
            context={"site_code": "lower_respiratory"},
        )
        self.assertIs(result["positive"], True)

        contradictory = dict(candidate)
        contradictory["specimen_class"] = "S2_lower_respiratory"
        contradictory["key_evidence"] = {
            "hospital_module_evidence": {
                "culture": [{"specimen_type": "Urine", "result": "Detected"}]
            }
        }
        result = module_c_direct_support(
            contradictory,
            self.config["module_c_direct_support"],
            context={"site_code": "lower_respiratory"},
        )
        self.assertIs(result["positive"], False)
        self.assertIn("did not match", result["site_alignment_basis"])

        candidate = {
            "specimen_alignment": "Aligned",
            "pcr": "Pending",
            "culture": {"status": "negative", "date": "2026-01-01"},
        }
        result = module_c_direct_support(
            candidate, self.config["module_c_direct_support"]
        )
        self.assertIs(result["positive"], False)

        candidate = {
            "module_support_summary": {"culture": "Support"},
        }
        result = module_c_direct_support(
            candidate, self.config["module_c_direct_support"]
        )
        self.assertIsNone(result["positive"])
        self.assertEqual(result["status"], "insufficient_input")

    def test_three_valued_experiment_truth_table(self):
        self.assertIs(tri_or(None, False), None)
        self.assertIs(tri_or(None, True), True)
        self.assertIs(tri_and(None, True), None)
        self.assertIs(tri_and(None, False), False)
        values = experiment_values(baseline=False, a=True, b=True, c=False)
        self.assertIs(values["ABC_MAJORITY"], True)
        self.assertIs(values["GATED_C_OR_A_AND_B"], True)
        self.assertIs(values["E0_OR_C_OR_A_AND_B"], True)
        preserved = experiment_values(baseline=True, a=False, b=False, c=False)
        self.assertIs(preserved["E0_OR_C_OR_A_AND_B"], True)

        pruning = experiment_values(baseline=True, a=None, b=True, c=False)
        self.assertIs(pruning["E0_AND_B"], True)
        self.assertIs(pruning["E0_AND_C"], False)
        self.assertIs(pruning["E0_AND_B_OR_C"], True)
        self.assertIs(pruning["E0_AND_B_AND_C"], False)

        outside_e0 = experiment_values(baseline=False, a=True, b=True, c=True)
        self.assertIs(outside_e0["E0_AND_B"], False)
        self.assertIs(outside_e0["E0_AND_C"], False)
        self.assertIs(outside_e0["E0_AND_B_OR_C"], False)
        self.assertIs(outside_e0["E0_AND_B_AND_C"], False)

        unknown_confirmation = experiment_values(
            baseline=True, a=False, b=None, c=False
        )
        self.assertIsNone(unknown_confirmation["E0_AND_B"])
        self.assertIsNone(unknown_confirmation["E0_AND_B_OR_C"])


if __name__ == "__main__":
    unittest.main()

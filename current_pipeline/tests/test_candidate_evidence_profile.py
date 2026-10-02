from __future__ import annotations

import unittest

from tools import candidate_evidence_profile as profiles
from tools import deterministic_mngs_max_scorer as scorer
from tools import merge_deterministic_max_with_missed_review as merge
from tools import review_missed_mngs_candidates as review


class CandidateEvidenceProfileTests(unittest.TestCase):
    def test_candida_blood_is_strong_invasive(self) -> None:
        strength, _ = profiles.candida_evidence_row_strength(
            {"specimen_type": "Blood", "specimen_category": "Blood"}
        )
        self.assertEqual(profiles.CANDIDA_EVIDENCE_STRONG_INVASIVE, strength)

    def test_candida_bal_is_weak_noninvasive(self) -> None:
        strength, _ = profiles.candida_evidence_row_strength(
            {"specimen_type": "BAL", "specimen_category": "Lower_Respiratory"}
        )
        self.assertEqual(profiles.CANDIDA_EVIDENCE_WEAK_NONINVASIVE, strength)

    def test_candida_bdg_is_intermediate(self) -> None:
        strength, _ = profiles.candida_evidence_row_strength(
            {"test": "Serum beta-D-glucan", "specimen_type": "Serum"}
        )
        self.assertEqual(profiles.CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC, strength)

    def test_profile_is_derived_from_nested_module_evidence(self) -> None:
        evidence = {
            "key_evidence": {
                "hospital_module_evidence": {
                    "culture": [
                        {"specimen_type": "Blood", "specimen_category": "Blood"}
                    ]
                }
            }
        }
        profile = profiles.candida_evidence_profile(evidence)
        self.assertEqual("strong_invasive", profile["best_strength"])

    def test_negative_flags_use_exact_support_and_candida_profile(self) -> None:
        snapshot = {
            "dominance_tier": "D1_low",
            "reads_tier": "R1_low",
            "specimen_class": "S2_lower_respiratory",
            "support_modules": ["host"],
            "non_host_support_modules": [],
        }
        hospital = {
            "same_organism_or_alias": [],
            "same_genus_or_related": [{"organism_name": "Candida albicans"}],
        }
        flags = profiles.build_negative_evidence(
            snapshot,
            hospital,
            [],
            organism_name="Candida tropicalis",
        )
        codes = [item["code"] for item in flags]
        self.assertIn("NO_SAME_ORGANISM_HOSPITAL_SUPPORT", codes)
        self.assertIn("NON_DOMINANT_MNGS", codes)
        self.assertIn("LOW_MNGS_SIGNAL", codes)
        self.assertIn("NO_NON_HOST_CONVERGENT_SUPPORT", codes)
        self.assertIn("CANDIDA_NO_STRONG_INVASIVE_EVIDENCE", codes)

    def test_all_pipeline_wrappers_share_central_strength(self) -> None:
        evidence = {
            "candida_evidence_strength_profile": {
                "best_strength": "intermediate_systemic_or_invasive_clue"
            }
        }
        self.assertEqual(
            profiles.CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC,
            review.candida_best_evidence_strength(evidence),
        )
        self.assertTrue(review.candida_has_intermediate_context(evidence))
        self.assertTrue(merge.candida_has_intermediate_context(evidence))
        self.assertEqual(
            profiles.candida_evidence_row_strength({"specimen_type": "BAL"}),
            scorer.candida_evidence_row_strength({"specimen_type": "BAL"}),
        )

    def test_review_tier_convergence_demotes_weak_candida(self) -> None:
        item = {
            "organism_name": "Candida albicans",
            "evidence_snapshot": {
                "classification": "Fungal",
                "reads_tier": "R1_low",
                "dominance_tier": "D0_not_top",
                "support_modules": ["host"],
                "non_host_support_modules": [],
            },
        }
        rule_id, _ = merge.review_tier_convergence_decision(item, {})
        self.assertEqual("RTCV2-WEAK-CANDIDA-NONINVASIVE", rule_id)

    def test_review_tier_convergence_keeps_invasive_candida(self) -> None:
        item = {
            "organism_name": "Candida albicans",
            "evidence_snapshot": {
                "classification": "Fungal",
                "reads_tier": "R1_low",
                "dominance_tier": "D0_not_top",
                "support_modules": ["host"],
                "non_host_support_modules": [],
                "candida_evidence_strength_profile": {"best_strength": "strong_invasive"},
            },
        }
        rule_id, _ = merge.review_tier_convergence_decision(item, {})
        self.assertEqual("", rule_id)

    def test_review_tier_convergence_keeps_c_striatum(self) -> None:
        item = {
            "organism_name": "Corynebacterium striatum",
            "evidence_snapshot": {
                "classification": "Bacterial",
                "reads_tier": "R3_high",
                "dominance_tier": "D0_not_top",
                "support_modules": ["host"],
                "non_host_support_modules": [],
            },
        }
        rule_id, _ = merge.review_tier_convergence_decision(item, {})
        self.assertEqual("", rule_id)

    def test_review_tier_convergence_demotes_non_striatum_coryne(self) -> None:
        item = {
            "organism_name": "Corynebacterium simulans",
            "evidence_snapshot": {
                "classification": "Bacterial",
                "reads_tier": "R3_high",
                "dominance_tier": "D0_not_top",
                "support_modules": ["host"],
                "non_host_support_modules": [],
            },
        }
        rule_id, _ = merge.review_tier_convergence_decision(item, {})
        self.assertEqual("RTCV2-NONSTRIATUM-CORYNE-WEAK", rule_id)

    def test_atypical_respiratory_igm_is_context_only_profile(self) -> None:
        profile = profiles.atypical_respiratory_serology_profile(
            "Chlamydophila pneumoniae",
            [
                {
                    "module_evidence": {
                        "filmarray_gmtest": [
                            {
                                "non_gm_test_type": "serology_IgM",
                                "evidence_type": "indirect_serology",
                            }
                        ]
                    }
                }
            ],
        )
        self.assertTrue(profile["context_rescue_match"])
        self.assertEqual("review_context_needed", profile["maximum_tier"])

    def test_atypical_respiratory_direct_pcr_is_not_indirect_only_rescue(self) -> None:
        profile = profiles.atypical_respiratory_serology_profile(
            "Chlamydia pneumoniae",
            [
                {
                    "module_evidence": {
                        "filmarray_gmtest": [
                            {
                                "non_gm_test_type": "serology_IgM",
                                "evidence_type": "indirect_serology",
                            }
                        ],
                        "molecular_microbiology": [
                            {"test": "respiratory PCR detected"}
                        ],
                    }
                }
            ],
        )
        self.assertTrue(profile["direct_molecular_found"])
        self.assertFalse(profile["context_rescue_match"])

    def test_review_normalizer_keeps_case_fit_atypical_igm_in_context(self) -> None:
        hospital_item = {
            "organism_name": "Chlamydophila pneumoniae",
            "module_evidence": {
                "filmarray_gmtest": [
                    {
                        "non_gm_test_type": "serology_IgM",
                        "evidence_type": "indirect_serology",
                    }
                ]
            },
        }
        tier, reason = review.direct_review_tier_for_item(
            {"organism_name": "Chlamydophila pneumoniae", "confidence": "moderate"},
            source_section="review_context_needed",
            evidence={"organism_name": "Chlamydophila pneumoniae", "reads": 0},
            payload={
                "hospital_side_summary": {
                    "cross_module_reasoning": "Severe pneumonia with bilateral infiltrates.",
                    "hospital_organism_evidence": [hospital_item],
                }
            },
            picked_keys=set(),
            picked_genera=set(),
        )
        self.assertEqual("review_context_needed", tier)
        self.assertIn("R-CTX-ATYPICAL-RESP-INDIRECT-SEROLOGY", reason)

    def test_merge_restores_existing_case_fit_atypical_igm_to_context(self) -> None:
        review_payload = {
            "review_high_priority": [],
            "review_context_needed": [],
            "review_low_specificity": [
                {
                    "organism_name": "Chlamydophila pneumoniae",
                    "review_tier": "review_low_specificity",
                    "original_review_section": "review_context_needed",
                    "evidence_snapshot": {
                        "hospital_evidence": [
                            {
                                "module_evidence": {
                                    "filmarray_gmtest": [
                                        {
                                            "non_gm_test_type": "serology_IgM",
                                            "evidence_type": "indirect_serology",
                                        }
                                    ]
                                }
                            }
                        ]
                    },
                }
            ],
            "review_omitted_with_reason": [],
        }
        output = merge.apply_atypical_respiratory_serology_context_rescue(review_payload)
        self.assertEqual(1, len(output["review_context_needed"]))
        self.assertEqual(0, len(output["review_low_specificity"]))
        rescued = output["review_context_needed"][0]
        self.assertEqual("review_context_needed", rescued["review_tier"])
        self.assertFalse(rescued["atypical_serology_context_rescue"]["formal_picked_modified"])


if __name__ == "__main__":
    unittest.main()

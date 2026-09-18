from __future__ import annotations

import copy
import unittest

from rag_re_clinical.clinical_rules import (
    evaluate_b_strict,
    evaluate_convergence,
    evaluate_modules,
    evaluate_temporal_coherence,
    grade_direct_evidence,
    normalize_site,
    specimen_colonization_gate,
)


def candidate(
    name: str = "Klebsiella pneumoniae",
    *,
    site: str = "lower_respiratory",
    alignment: str = "Aligned",
    key_evidence=None,
    support_modules=None,
):
    return {
        "organism_name": name,
        "canonical_organism_name": name,
        "source_context": {
            "site_code": site,
            "specimen_site_code": site,
            "specimen_alignment": alignment,
        },
        "rule_input": {
            "organism_name": name,
            "classification": "Fungal" if "Candida" in name else "Bacterial",
            "mngs_signal_tier": "M2_moderate",
            "reads_tier": "R3_high",
            "dominance_tier": "D2",
            "rank_priority": "1",
            "reads_percentile": 0.90,
            "specimen_alignment": alignment,
            "specimen_class": "S2_lower_respiratory" if site == "lower_respiratory" else site,
            "key_evidence": key_evidence or {},
            "support_modules": support_modules,
        },
        "modules": {"A": {"positive": False, "status": "ok"}},
    }


def direct_evidence(module, record):
    return {"hospital_module_evidence": {module: [record]}}


class SiteNormalizationTests(unittest.TestCase):
    def test_common_sites(self):
        self.assertEqual(normalize_site("Lower BAL"), "lower_respiratory")
        self.assertEqual(normalize_site("Endotracheal aspirate"), "lower_respiratory")
        self.assertEqual(normalize_site("Blood"), "bloodstream")
        self.assertEqual(normalize_site("Foley catheter"), "urinary")
        self.assertIsNone(normalize_site("Not_available"))


class SpecimenGateTests(unittest.TestCase):
    def test_structured_mismatch_beats_generic_aligned(self):
        row = candidate(
            key_evidence=direct_evidence(
                "culture",
                {"status": "positive", "specimen_type": "Urine"},
            )
        )
        result = specimen_colonization_gate({}, row)
        self.assertEqual(result["state"], "D_BLOCK")
        self.assertFalse(result["site_aligned"])
        self.assertTrue(result["direct_site_priority_applied"])

    def test_resp_candida_is_guarded(self):
        row = candidate("Candida albicans")
        result = specimen_colonization_gate({}, row)
        self.assertEqual(result["state"], "D_GUARDED")

    def test_explicit_false_background_flag_does_not_create_guard(self):
        row = candidate()
        row["rule_input"]["is_likely_colonizer_or_background"] = False
        row["rule_input"]["applied_rules"] = [
            "D-S1-CANDIDATE_FROM_HOSPITAL_ONLY",
            "R-HOSPITAL-ONLY-L1",
        ]
        result = specimen_colonization_gate({}, row)
        self.assertEqual(result["state"], "D_PASS")

    def test_negative_or_unknown_contaminant_flags_do_not_create_guard(self):
        for value in ("No", "Unknown", False):
            with self.subTest(value=value):
                row = candidate()
                merged = {
                    "effective_candidate": {
                        "organism_name": "Klebsiella pneumoniae",
                        "key_evidence": {
                            "hospital_module_evidence": {
                                "culture": [
                                    {
                                        "status": "positive",
                                        "specimen_type": "BAL",
                                        "contaminant_flag": value,
                                    }
                                ]
                            }
                        },
                    }
                }
                result = specimen_colonization_gate(merged, row)
                self.assertEqual(result["state"], "D_PASS")
                self.assertEqual(result["guard_reasons"], [])

    def test_benign_positive_guardrail_rule_is_not_treated_as_colonization(self):
        row = candidate("Haemophilus influenzae")
        row["rule_input"]["key_evidence"] = {
            "guardrail_rule": "R-S4-TYPICAL_PNEUMONIA_HOSPITAL_MNGS_L2"
        }
        result = specimen_colonization_gate({}, row)
        self.assertEqual(result["state"], "D_PASS")
        self.assertEqual(result["guard_reasons"], [])

    def test_exact_matched_raw_background_flag_is_not_lost_by_frozen_allowlist(self):
        row = candidate("Corynebacterium striatum")
        merged = {
            "match_scope": "exact_species",
            "effective_candidate": {
                "organism_name": "Corynebacterium striatum",
                "is_likely_colonizer_or_background": True,
                "formal_pick_exclusion_rule": "LOW_DOMINANCE_BACKGROUND_NOT_PICKED",
            },
            "deterministic_candidate": {
                "organism_name": "Corynebacterium striatum",
                "is_likely_colonizer_or_background": True,
            },
        }
        result = specimen_colonization_gate(merged, row)
        self.assertEqual(result["state"], "D_GUARDED")
        self.assertIn("upstream colonizer/background flag", result["guard_reasons"])


class StrictMngsTests(unittest.TestCase):
    def test_strict_two_axes_and_site_but_ntc_rpm_are_partial(self):
        row = candidate()
        result = evaluate_b_strict({}, row)
        self.assertTrue(result["positive"])
        self.assertEqual(result["status"], "partial")
        self.assertEqual(
            set(result["missing_quality_features"]), {"ntc_control_ratio", "rpm"}
        )

    def test_missing_site_is_unknown_not_false(self):
        row = candidate()
        row["source_context"] = {"site_code": "unknown", "specimen_site_code": "unknown"}
        row["rule_input"]["specimen_alignment"] = "Not_available"
        row["rule_input"]["specimen_class"] = "Not_available"
        result = evaluate_b_strict({}, row)
        self.assertIsNone(result["positive"])
        self.assertEqual(result["status"], "missing")

    def test_unrelated_modules_do_not_enter_b_strict(self):
        row = candidate()
        baseline = evaluate_b_strict({}, row)
        row["modules"]["A"]["positive"] = True
        row["modules"]["C"] = {"positive": True, "grade": "C3"}
        row["host_response"] = "Support"
        self.assertEqual(evaluate_b_strict({}, row)["positive"], baseline["positive"])


class DirectEvidenceGradeTests(unittest.TestCase):
    def test_aligned_structured_bal_record_is_c2(self):
        row = candidate(
            key_evidence=direct_evidence(
                "filmarray_gmtest",
                {
                    "organism_name": "Klebsiella pneumoniae",
                    "detection_status": "detected",
                    "filmarray_sample_type": "BAL",
                    "evidence_type": "direct_detection",
                },
            )
        )
        result = grade_direct_evidence({}, row)
        self.assertEqual(result["grade"], "C2")
        self.assertTrue(result["positive"])

    def test_combined_filmarray_gm_record_without_exact_target_is_only_c1(self):
        row = candidate(
            "Aspergillus spp.",
            key_evidence=direct_evidence(
                "filmarray_gmtest",
                {
                    "detection_status": "detected",
                    "filmarray_sample_type": "BAL",
                    "filmarray_panel": "Unknown",
                    "non_gm_test_type": "none",
                },
            ),
        )
        result = grade_direct_evidence({}, row)
        self.assertEqual(result["grade"], "C1")
        self.assertTrue(result["manual_review"])

    def test_aligned_sterile_record_is_c3(self):
        row = candidate(site="bloodstream")
        row["rule_input"]["key_evidence"] = direct_evidence(
            "culture",
            {
                "specimen_type": "Blood",
                "specimen_category": "Sterile_Site",
                "growth_purity": "isolated_pure",
            },
        )
        result = grade_direct_evidence({}, row)
        self.assertEqual(result["grade"], "C3")

    def test_sterile_culture_flagged_contaminant_is_manual_c1(self):
        row = candidate(site="bloodstream")
        row["rule_input"]["key_evidence"] = direct_evidence(
            "culture",
            {
                "specimen_type": "Blood",
                "specimen_category": "Sterile_Site",
                "status": "positive",
                "contaminant_flag": "Yes",
            },
        )
        result = grade_direct_evidence({}, row)
        self.assertEqual(result["grade"], "C1")
        self.assertTrue(result["manual_review"])

    def test_support_modules_alone_is_only_c1(self):
        row = candidate()
        row["rule_input"]["key_evidence"] = {"support_modules": ["culture"]}
        result = grade_direct_evidence({}, row)
        self.assertEqual(result["grade"], "C1")
        self.assertIsNone(result["positive"])

    def test_resp_candida_culture_only_is_c1_manual(self):
        row = candidate(
            "Candida albicans",
            key_evidence=direct_evidence(
                "culture",
                {
                    "specimen_type": "Lower BAL",
                    "specimen_category": "Lower_Respiratory",
                    "quantity_tier": "Q4",
                },
            ),
        )
        result = grade_direct_evidence({}, row)
        self.assertEqual(result["grade"], "C1")
        self.assertTrue(result["manual_review"])
        self.assertTrue(result["respiratory_candida_culture_only"])

    def test_related_representative_is_not_exact(self):
        row = candidate("Candida tropicalis")
        row["rule_input"]["key_evidence"] = {
            "related_representative_hospital_support": [
                {
                    "organism_name": "Candida albicans",
                    "support_type": "related_representative_context",
                    "hospital_module_evidence": {
                        "culture": [
                            {"status": "positive", "specimen_type": "Lower BAL"}
                        ]
                    },
                }
            ]
        }
        result = grade_direct_evidence({}, row)
        self.assertNotIn(result["grade"], {"C2", "C3"})
        self.assertFalse(result["related_representative_is_exact"])

    def test_pending_not_done_and_missing_are_unknown_not_cneg(self):
        for status in ("pending", "not_done", "missing"):
            with self.subTest(status=status):
                row = candidate(
                    key_evidence=direct_evidence(
                        "pcr", {"status": status, "specimen_type": "BAL"}
                    )
                )
                result = grade_direct_evidence({}, row)
                self.assertEqual(result["grade"], "C0")
                self.assertEqual(result["status"], "unknown")
                self.assertIsNone(result["positive"])

    def test_explicit_aligned_negative_is_cneg(self):
        row = candidate(
            key_evidence=direct_evidence(
                "pcr", {"result": "not_detected", "specimen_type": "BAL"}
            )
        )
        result = grade_direct_evidence({}, row)
        self.assertEqual(result["grade"], "CNEG")
        self.assertFalse(result["positive"])

    def test_mismatched_structured_site_cannot_be_promoted_by_aligned(self):
        row = candidate(
            key_evidence=direct_evidence(
                "culture", {"status": "positive", "specimen_type": "Urine"}
            )
        )
        result = grade_direct_evidence({}, row)
        self.assertEqual(result["grade"], "C1")
        self.assertTrue(result["manual_review"])


class TemporalAndConvergenceTests(unittest.TestCase):
    def test_temporal_coherence_and_post_treatment_flag(self):
        row = candidate()
        row["rule_input"]["key_evidence"] = direct_evidence(
            "pcr",
            {
                "status": "positive",
                "specimen_type": "BAL",
                "collection_date": "2026-08-11",
            },
        )
        context = {
            "syndrome_onset": "2026-08-10",
            "antimicrobial_start": "2026-08-09",
        }
        result = evaluate_temporal_coherence(context, row)
        self.assertEqual(result["state"], "F_COHERENT")
        self.assertTrue(result["post_treatment_specimen"])

    def test_missing_temporal_data_is_unknown(self):
        result = evaluate_temporal_coherence({}, candidate())
        self.assertEqual(result["state"], "F_UNKNOWN")
        self.assertIsNone(result["positive"])

    def test_e2_requires_d_pass_b_and_c2(self):
        result = evaluate_convergence(
            {"state": "D_PASS"},
            {"positive": True},
            {"grade": "C2"},
        )
        self.assertEqual(result["tier"], "E2_CONVERGENT_HIGH")
        self.assertTrue(result["auto_rescue"])
        self.assertTrue(result["e0_immutable"])

    def test_guarded_c2_is_manual(self):
        result = evaluate_convergence(
            {"state": "D_GUARDED"},
            {"positive": True},
            {"grade": "C2"},
        )
        self.assertEqual(result["disposition"], "manual_review")
        self.assertFalse(result["auto_rescue"])

    def test_full_entrypoint_is_pure(self):
        row = candidate(
            key_evidence=direct_evidence(
                "filmarray_gmtest",
                {"detection_status": "detected", "filmarray_sample_type": "BAL"},
            )
        )
        before = copy.deepcopy(row)
        result = evaluate_modules({}, row, {})
        self.assertEqual(row, before)
        self.assertEqual(set(result), {"D", "B_STRICT", "C", "E", "F", "contract"})
        self.assertTrue(result["contract"]["e0_immutable"])


if __name__ == "__main__":
    unittest.main()

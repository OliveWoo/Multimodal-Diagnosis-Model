from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.build_test_aware_deterministic_shadow_scorer import (
    analytical_profile,
    decide_entry,
    expected_scorer_entry_count,
    forward_to_clinical_scorer,
    hospital_only_shadow,
    hospital_profile,
    hospital_observation_index,
)


POLICY = json.loads(Path("rules/test_aware_deterministic_shadow_v2.json").read_text(encoding="utf-8"))
ROUTE_A_POLICY = json.loads(
    Path("rules/test_aware_deterministic_shadow_v3_route_a.json").read_text(
        encoding="utf-8"
    )
)


def entry(*, family="typical_respiratory_pathogen", rank=1, count=2, cross=True, reads=20):
    signals = [
        {
            "seq_id": f"S{index}",
            "nucleic_type": "DNA" if index == 1 else "RNA",
            "reads": reads,
            "rank_in_retained_universe": rank,
            "analytical_rank_priority": 1,
            "rpm_total": None,
            "qc_status": "evaluable",
            "normalization_available": False,
        }
        for index in range(1, count + 1)
    ]
    return {
        "organism_name": "Haemophilus influenzae",
        "organism_key": "haemophilusinfluenzae",
        "taxonomy_family": family,
        "screening_tier": "medium",
        "specimen_context": "lower_respiratory",
        "selected_test_signals": signals,
        "cross_molecule_selected": cross,
        "dna_selected": True,
        "rna_selected": cross,
    }


class TestAwareDeterministicShadowTests(unittest.TestCase):
    def test_raw_molecular_observations_are_reconstructed_for_evidence_links(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "molecular.json"
            source.write_text(
                json.dumps(
                    [
                        {
                            "sample": "BALF",
                            "target": "Pneumocystis jirovecii",
                            "result": "Detected",
                        },
                        {
                            "sample": "BALF",
                            "target": "CMV",
                            "result": "Not Detected",
                        },
                    ]
                ),
                encoding="utf-8",
            )
            index = hospital_observation_index(
                {"source_files": {"raw_molecular_microbiology": str(source)}}
            )
            self.assertEqual({"MOL-001", "MOL-002"}, set(index))
            self.assertEqual(
                "detected_positive", index["MOL-001"]["result_interpretation"]
            )
            self.assertEqual(
                "not_detected", index["MOL-002"]["result_interpretation"]
            )

    def test_expected_entry_count_comes_from_compact_summary(self):
        self.assertEqual(420, expected_scorer_entry_count({"counts": {"scorer_entry": 420}}))
        with self.assertRaises(ValueError):
            expected_scorer_entry_count({"counts": {}})

    def test_positive_signals_remain_separate_and_are_not_summed(self):
        profile = analytical_profile(entry(reads=40), POLICY)
        self.assertEqual([40, 40], [row["reads"] for row in profile["per_test_signals"]])
        self.assertIsNone(profile["reads_sum_across_tests"])
        self.assertEqual(1, sum([profile["reproducibility_axis"]]))

    def test_zero_reads_are_rejected(self):
        with self.assertRaises(ValueError):
            analytical_profile(entry(reads=0), POLICY)

    def test_no_positive_mngs_signal_does_not_block_hospital_only_candidate(self):
        final_summary = {
            "hospital_organism_evidence": [{
                "organism_name": "Streptococcus pneumoniae",
                "classification": "Bacterial",
                "best_hospital_level": "Level 2",
                "evidence_modules": ["culture"],
                "module_level_summary": {"culture": "Level 2"},
                "module_evidence": {
                    "culture": [{
                        "specimen_type": "BALF",
                        "specimen_category": "Lower_Respiratory",
                        "quantitation_status": "positive_detected_quantity_unknown",
                        "evidence_observation_ids": ["CUL-HOSP-ONLY"],
                    }]
                },
            }]
        }

        result = hospital_only_shadow(
            final_summary,
            patient_candidates=[],
            host={},
            observation_index={
                "CUL-HOSP-ONLY": {
                    "observation_id": "CUL-HOSP-ONLY",
                    "organism_name": "Streptococcus pneumoniae",
                    "specimen_type": "BALF",
                    "specimen_category": "Lower_Respiratory",
                    "collected_time": "2024-11-14 09:30",
                    "raw_result": "Isolated",
                }
            },
        )

        self.assertEqual(1, len(result))
        self.assertEqual("Streptococcus pneumoniae", result[0]["organism_name"])
        self.assertEqual("hospital_only", result[0]["evidence_source"])
        self.assertTrue(result[0]["forward_to_clinical_scorer"])

    def test_priority_top1_reproducible_candidate_can_be_shadow_picked(self):
        result = decide_entry(entry(), final_by_name={}, host={}, policy=POLICY)
        self.assertEqual("picked_shadow", result["decision"])
        self.assertIn("TA-S4-ANALYTICAL-TOP1-REPRODUCIBLE-PICK", result["rule_ids"])

    def test_candida_is_not_analytical_only_picked(self):
        item = entry(family="candida_or_yeast")
        item["organism_name"] = "Candida albicans"
        result = decide_entry(item, final_by_name={}, host={}, policy=POLICY)
        self.assertEqual("review_context_needed", result["decision"])
        self.assertFalse(result["formal_pick_allowed"])

    def test_opportunistic_candidate_requires_more_than_analytical_signal(self):
        item = entry(family="high_consequence_opportunistic")
        item["organism_name"] = "Pneumocystis jirovecii"
        result = decide_entry(item, final_by_name={}, host={}, policy=POLICY)
        self.assertNotEqual("picked_shadow", result["decision"])

    def test_single_test_priority_is_high_review_not_forced_pick(self):
        result = decide_entry(entry(count=1, cross=False), final_by_name={}, host={}, policy=POLICY)
        self.assertEqual("review_high_priority", result["decision"])
        self.assertFalse(result["formal_pick_allowed"])

    def test_reproducible_top_ten_opportunist_is_forwarded_as_context(self):
        item = entry(family="high_consequence_opportunistic", rank=4, count=3)
        item["organism_name"] = "Pneumocystis jirovecii"
        result = decide_entry(item, final_by_name={}, host={}, policy=POLICY)
        self.assertEqual("review_context_needed", result["decision"])
        self.assertFalse(result["formal_pick_allowed"])

    def test_top_three_low_specificity_signal_is_not_discarded(self):
        item = entry(family="gi_urinary_or_nonpulmonary_prone", rank=2, count=1, cross=False)
        item["organism_name"] = "Enterococcus faecium"
        result = decide_entry(item, final_by_name={}, host={}, policy=POLICY)
        self.assertEqual("review_context_needed", result["decision"])
        self.assertFalse(result["formal_pick_allowed"])

    def test_repeated_below_top_ten_signal_remains_low(self):
        item = entry(family="environmental_low_specificity", rank=11, count=2)
        result = decide_entry(item, final_by_name={}, host={}, policy=POLICY)
        self.assertEqual("review_low_specificity", result["decision"])

    def test_route_a_moves_compelling_low_specificity_repeat_to_high_only(self):
        item = entry(
            family="skin_airway_colonizer_prone", rank=1, count=2, cross=False
        )
        item.update({
            "organism_name": "Corynebacterium striatum",
            "taxonomy_mapping_status": "exact_species",
        })
        for index, signal in enumerate(item["selected_test_signals"]):
            signal.update({
                "rpm_total": 205.0 + index,
                "normalization_available": True,
                "protocol_condition": "-" if index == 0 else "12000g/10min",
                "test_start_time": "2024/10/17 23:59",
            })
        result = decide_entry(
            item, final_by_name={}, host={}, policy=ROUTE_A_POLICY
        )
        self.assertEqual("review_high_priority", result["decision"])
        self.assertFalse(result["formal_pick_allowed"])
        self.assertIn(
            "TA-S4-LOW-SPECIFICITY-REPEAT-HIGH", result["rule_ids"]
        )
        self.assertTrue(
            result["analytical_profile"]["low_specificity_repeat_high_gate"][
                "eligible"
            ]
        )

    def test_route_a_rejects_explicit_same_library_rerun(self):
        item = entry(
            family="environmental_low_specificity", rank=1, count=2, cross=False
        )
        item["taxonomy_mapping_status"] = "exact_species"
        for signal in item["selected_test_signals"]:
            signal.update({
                "rpm_total": 100.0,
                "normalization_available": True,
                "protocol_condition": "原library重新上機",
                "test_start_time": "2024/07/17 12:00",
            })
        result = decide_entry(
            item, final_by_name={}, host={}, policy=ROUTE_A_POLICY
        )
        self.assertNotIn(
            "TA-S4-LOW-SPECIFICITY-REPEAT-HIGH", result["rule_ids"]
        )

    def test_route_a_does_not_replace_existing_cross_molecule_route(self):
        item = entry(
            family="skin_airway_colonizer_prone", rank=1, count=4, cross=True
        )
        item["taxonomy_mapping_status"] = "exact_species"
        for signal in item["selected_test_signals"]:
            signal.update({
                "rpm_total": 100.0,
                "normalization_available": True,
                "protocol_condition": "-",
                "test_start_time": "2024/08/27 12:00",
            })
        result = decide_entry(
            item, final_by_name={}, host={}, policy=ROUTE_A_POLICY
        )
        self.assertNotIn(
            "TA-S4-LOW-SPECIFICITY-REPEAT-HIGH", result["rule_ids"]
        )
        self.assertEqual("review_context_needed", result["decision"])

    def test_forwarding_is_separate_from_picked(self):
        self.assertTrue(forward_to_clinical_scorer("review_context_needed", POLICY))
        self.assertFalse(forward_to_clinical_scorer("review_low_specificity", POLICY))

    def test_hospital_profile_preserves_source_observation_detail(self):
        evidence = {
            "organism_name": "Candida dubliniensis",
            "classification": "Fungal",
            "best_hospital_level": "Level 2",
            "evidence_modules": ["culture"],
            "module_level_summary": {"culture": "Level 2"},
            "module_evidence": {
                "culture": [{
                    "specimen_type": "Blood",
                    "specimen_category": "Sterile_Site",
                    "quantitation_status": "positive_detected_quantity_unknown",
                    "evidence_observation_ids": ["CUL-001"],
                }]
            },
        }
        profile = hospital_profile(
            "Candida dubliniensis",
            {"candidadubliniensis": evidence},
            {"CUL-001": {
                "observation_id": "CUL-001",
                "organism_name": "Candida dubliniensis",
                "specimen_type": "Blood",
                "specimen_category": "Sterile_Site",
                "collected_time": "2024-11-13 21:37",
                "raw_result": "Isolated",
            }},
        )
        detail = profile["hospital_evidence_detail"]
        self.assertEqual("CUL-001", detail["source_observations"][0]["observation_id"])
        self.assertEqual(1, detail["source_observation_resolution"]["resolved"])

    def test_parainfluenza_group_is_related_context_not_exact_hpiv3_evidence(self):
        final_by_name = {
            "parainfluenzavirus": {
                "organism_name": "Parainfluenza Virus",
                "best_hospital_level": "Level 2",
                "module_level_summary": {
                    "filmarray_gmtest": "Level 2",
                },
                "module_evidence": {
                    "filmarray_gmtest": [{
                        "filmarray_sample_type": "BAL",
                        "detection_status": "detected",
                        "evidence_observation_ids": ["ASSAY-024"],
                    }],
                },
            }
        }
        profile = hospital_profile("Human respirovirus 3", final_by_name)
        self.assertFalse(profile["exact_or_alias_evidence"])
        self.assertIsNone(profile["best_direct_level"])
        self.assertEqual(1, len(profile["related_representative_context"]))
        related = profile["related_representative_context"][0]
        self.assertEqual("approved_group_member_context", related["support_type"])
        self.assertEqual(
            "group_level_support_not_exact_species_confirmation",
            related["identity_scope"],
        )


if __name__ == "__main__":
    unittest.main()

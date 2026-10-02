from __future__ import annotations

import unittest

from tools import deterministic_mngs_max_scorer as scorer
from tools import build_missed_candidate_review_queue as queue_builder


def culture_candidate(specimen_type: str, specimen_category: str) -> dict:
    return {
        "organism_name": "Enterococcus faecium (VRE)",
        "evidence_modules": ["culture"],
        "module_level_summary": {"culture": "Level 3"},
        "module_evidence": {
            "culture": [
                {
                    "specimen_type": specimen_type,
                    "specimen_category": specimen_category,
                }
            ]
        },
    }


class PulmonaryCultureSupportTests(unittest.TestCase):
    def test_explicit_opportunistic_host_false_overrides_acute_v3(self) -> None:
        context = scorer.host_context(
            {
                "host_context": {
                    "host_vulnerability_tier": "V3",
                    "opportunistic_coverage_level": "O2",
                    "expanded_candidate_policy": True,
                    "opportunistic_host_support": False,
                    "opportunistic_host_review_context": False,
                    "immunocompromise_tier": "H1",
                    "acute_instability_tier": "A3",
                }
            }
        )

        self.assertFalse(scorer.host_support(context))
        self.assertFalse(queue_builder.is_host_vulnerable(context))

    def test_incomplete_immunosuppressant_evidence_remains_review_context(self) -> None:
        context = scorer.host_context(
            {
                "host_context": {
                    "host_vulnerability_tier": "V3",
                    "opportunistic_host_support": False,
                    "opportunistic_host_review_context": True,
                }
            }
        )

        self.assertFalse(scorer.host_support(context))
        self.assertTrue(queue_builder.is_host_vulnerable(context))

    def test_foley_culture_is_context_only(self) -> None:
        candidate = culture_candidate("Foley catheter", "NonSterile_Other")
        final_by_name = {"enterococcusfaecium": candidate}

        self.assertFalse(
            scorer.has_same_organism_culture_support(final_by_name, "Enterococcus faecium")
        )
        self.assertNotIn(
            "culture",
            scorer.direct_hospital_support_modules(candidate, max_level_rank=3),
        )

    def test_lower_bal_culture_supports_pulmonary_candidate(self) -> None:
        candidate = culture_candidate("Lower BAL", "Lower_Respiratory")
        final_by_name = {"enterococcusfaecium": candidate}

        self.assertTrue(
            scorer.has_same_organism_culture_support(final_by_name, "Enterococcus faecium")
        )
        self.assertIn(
            "culture",
            scorer.direct_hospital_support_modules(candidate, max_level_rank=3),
        )

    def test_blood_culture_remains_relevant(self) -> None:
        candidate = culture_candidate("Blood", "Sterile_Site")
        final_by_name = {"enterococcusfaecium": candidate}

        self.assertTrue(
            scorer.has_same_organism_culture_support(final_by_name, "Enterococcus faecium")
        )

    def test_queue_keeps_foley_evidence_as_context_not_support(self) -> None:
        candidates: dict[str, dict] = {}
        evidence = culture_candidate("Foley catheter", "NonSterile_Other")

        queue_builder.add_hospital_candidate(
            candidates,
            organism_name="Enterococcus faecium (VRE)",
            classification="Bacterial",
            level="Level 3",
            source_module="culture",
            source_field="hospital_organism_evidence",
            evidence=evidence,
            pulmonary_causality_support=False,
        )

        item = next(iter(candidates.values()))
        self.assertEqual(item["support_modules"], [])
        self.assertEqual(item["non_host_support_modules"], [])
        self.assertEqual(item["context_only_modules"], ["culture"])
        self.assertEqual(
            "gi_urinary_or_nonpulmonary_prone",
            item["taxonomy_profile"]["primary_rule_family"],
        )

    def _hospital_evidence(
        self,
        organism_name: str,
        *,
        classification: str,
        observation_count: int,
    ) -> dict:
        return {
            "organism_name": organism_name,
            "classification": classification,
            "evidence_modules": ["culture"],
            "module_level_summary": {"culture": "Level 2"},
            "module_evidence": {
                "culture": [
                    {
                        "specimen_type": "Lower BAL",
                        "specimen_category": "Lower_Respiratory",
                        "quantity_tier": "Q3",
                        "quantitation_status": "reported",
                        "observation_count": observation_count,
                    }
                ]
            },
        }

    def _score(
        self,
        organism_name: str,
        *,
        source_category: str,
        reads: int,
        rank: int,
        percentile: float,
        specimen_name: str,
        dominance: str,
        final_by_name: dict | None = None,
        host_expansion: bool = False,
    ) -> dict:
        raw_candidate = {
            "organism_name": organism_name,
            "source_category": source_category,
            "reads": reads,
            "ranking": {
                "rank_priority": rank,
                "reads_percentile": percentile,
                "rank_rule": "Code_NTC=00 + Code_RK_NTC=RK00",
                "possibility_level": "high",
            },
        }
        return scorer.score_candidate(
            raw_candidate,
            {"specimen_name": specimen_name},
            {scorer.mngs.normalize_organism_name(organism_name): dominance},
            final_by_name or {},
            {
                "host_level3_expansion": host_expansion,
                "opportunistic_host_support": host_expansion,
            },
        )

    def test_blood_only_cmv_host_context_is_review_not_picked(self) -> None:
        candidate = self._score(
            "CMV",
            source_category="Virus",
            reads=15,
            rank=1,
            percentile=0.67,
            specimen_name="Blood",
            dominance="D0_not_top",
            host_expansion=True,
        )

        self.assertEqual("Level 3", candidate["integrated_causative_level"])
        self.assertTrue(scorer.formal_pick_excluded(candidate))
        self.assertFalse(scorer.pickable(candidate, {"host_level3_expansion": True}))

    def test_nondominant_hsv_without_direct_support_has_level4_fallback(self) -> None:
        candidate = self._score(
            "HSV-1",
            source_category="Virus",
            reads=1640,
            rank=1,
            percentile=0.92,
            specimen_name="Lower BAL",
            dominance="D0_not_top",
        )

        self.assertEqual("Level 4", candidate["integrated_causative_level"])
        self.assertEqual("M4_weak", candidate["mngs_signal_tier"])
        self.assertIn(scorer.HERPES_NONDIRECT_FALLBACK_RULE, candidate["applied_rules"])

    def test_repeated_quantified_trichosporon_culture_can_replace_host_requirement(self) -> None:
        evidence = self._hospital_evidence(
            "Trichosporon asahii",
            classification="Fungal",
            observation_count=3,
        )
        candidate = self._score(
            "Trichosporon asahii",
            source_category="Fungi",
            reads=4223,
            rank=1,
            percentile=0.93,
            specimen_name="Lower BAL",
            dominance="D1_low",
            final_by_name={"trichosporonasahii": evidence},
        )

        self.assertEqual("Level 3", candidate["integrated_causative_level"])
        self.assertTrue(scorer.pickable(candidate, {"host_level3_expansion": False}))
        self.assertIn(scorer.RARE_YEAST_REPEATED_LR_CULTURE_RULE, candidate["applied_rules"])

    def test_exact_species_repeated_culture_beats_broad_group_context(self) -> None:
        evidence = self._hospital_evidence(
            "Acinetobacter baumannii",
            classification="Bacterial",
            observation_count=4,
        )
        candidate = self._score(
            "Acinetobacter baumannii",
            source_category="Bacteria",
            reads=696,
            rank=5,
            percentile=0.90,
            specimen_name="Lower BAL",
            dominance="D0_not_top",
            final_by_name={"acinetobacterbaumannii": evidence},
        )

        self.assertEqual("Level 3", candidate["integrated_causative_level"])
        self.assertTrue(scorer.pickable(candidate, {"host_level3_expansion": False}))
        self.assertIn(scorer.TYPICAL_EXACT_REPEATED_LR_CULTURE_RULE, candidate["applied_rules"])

    def test_no_level4_or_level5_candidate_is_forced_into_picked(self) -> None:
        candidate = {
            "organism_name": "Context organism",
            "classification": "Bacterial",
            "integrated_causative_level": "Level 4",
            "mngs_signal_tier": "M4_weak",
            "key_evidence": {},
        }

        best = scorer.build_best_available([candidate], {"host_level3_expansion": False})

        self.assertEqual("No_high_priority_candidate", best["selection_mode"])
        self.assertEqual(0, best["picked_count"])
        self.assertEqual([], best["picked_pathogens"])

    def test_image_respiratory_virus_context_builds_nonpicked_candidate(self) -> None:
        evidence = {
            "organism_name": "SARS-CoV-2",
            "classification": "Viral",
            "best_hospital_level": "Level 4",
            "evidence_modules": ["image"],
            "module_level_summary": {"image": "Level 4"},
            "module_evidence": {
                "image": [
                    {
                        "context_candidate_eligible": True,
                        "microbiologic_confirmation": False,
                        "source_sentence": "Covid-19 pneumonia.",
                    }
                ]
            },
        }

        candidate = scorer.hospital_only_candidate_from_evidence(evidence, {})

        self.assertIsNotNone(candidate)
        self.assertEqual("Level 4", candidate["integrated_causative_level"])
        self.assertEqual("Context_only", candidate["module_support_summary"]["image"])
        self.assertTrue(scorer.formal_pick_excluded(candidate))
        self.assertFalse(scorer.pickable(candidate, {}))


if __name__ == "__main__":
    unittest.main()

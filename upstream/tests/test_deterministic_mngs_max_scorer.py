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


if __name__ == "__main__":
    unittest.main()

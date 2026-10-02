import unittest

from tools.rehydrate_clinical_hospital_evidence_shadow import hydrate_candidate


class ClinicalEvidenceRehydrationTests(unittest.TestCase):
    def setUp(self):
        self.final_by_name = {
            "candidaalbicans": {
                "organism_name": "Candida albicans",
                "classification": "Fungal",
                "best_hospital_level": "Level 2",
                "module_level_summary": {"culture": "Level 2"},
                "module_evidence": {
                    "culture": [{
                        "evidence_observation_ids": ["CUL-001"],
                        "specimen_category": "Sterile_Site",
                    }]
                },
            }
        }
        self.observations = {
            "CUL-001": {
                "observation_id": "CUL-001",
                "organism_name": "Candida albicans",
                "raw_result": "Isolated",
                "collected_time": "2024-01-01 10:00",
                "specimen_type": "Blood",
                "specimen_category": "Sterile_Site",
            }
        }

    def candidate(self):
        return {
            "organism_name": "Candida albicans",
            "organism_key": "candidaalbicans",
            "candidate_source": "mngs",
            "clinical_decision": "review_high_priority",
            "formal_pick_allowed": False,
            "analytic_route": "Priority",
            "history_adjusted_route": "Priority",
            "history_route_changed": False,
            "history_rule_ids": [],
            "history_evidence_ids": [],
            "consumed_evidence_ids": [],
            "hospital_profile": {
                "exact_or_alias_evidence": True,
                "matched_hospital_name": "Candida albicans",
                "best_direct_level": "Level 2",
                "direct_support_modules_level_1_2": ["culture"],
                "direct_support_modules_level_1_3": ["culture"],
                "module_level_summary": {"culture": "Level 2"},
                "related_representative_context": [],
            },
        }

    def test_rehydrates_source_rows_without_changing_decision(self):
        before = self.candidate()
        result, audit = hydrate_candidate(
            before, self.final_by_name, self.observations
        )
        detail = result["hospital_profile"]["hospital_evidence_detail"]
        self.assertEqual(detail["source_observations"][0]["raw_result"], "Isolated")
        self.assertEqual(audit["requested"], 1)
        self.assertEqual(audit["resolved"], 1)
        self.assertEqual(result["clinical_decision"], before["clinical_decision"])
        self.assertFalse(
            result["evidence_rehydration"]["semantic_decision_fields_changed"]
        )

    def test_semantic_hospital_profile_change_is_rejected(self):
        item = self.candidate()
        item["hospital_profile"]["best_direct_level"] = "Level 1"
        with self.assertRaisesRegex(ValueError, "semantics changed"):
            hydrate_candidate(item, self.final_by_name, self.observations)

    def test_hospital_only_detail_is_rehydrated(self):
        item = {
            "organism_name": "Candida albicans",
            "candidate_source": "hospital_only",
            "clinical_decision": "review_high_priority",
            "hospital_evidence_detail": {
                "module_evidence": {
                    "culture": [{"evidence_observation_ids": ["CUL-001"]}]
                }
            },
        }
        result, audit = hydrate_candidate(item, self.final_by_name, self.observations)
        self.assertEqual(audit["resolved"], 1)
        self.assertEqual(
            result["hospital_evidence_detail"]["source_observations"][0]["observation_id"],
            "CUL-001",
        )


if __name__ == "__main__":
    unittest.main()

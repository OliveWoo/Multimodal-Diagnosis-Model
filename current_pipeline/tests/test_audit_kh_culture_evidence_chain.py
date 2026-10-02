import unittest

from tools.audit_kh_culture_evidence_chain import repair_agent_provenance


class CultureEvidenceChainAuditTests(unittest.TestCase):
    def test_repair_adds_missing_source_observation_without_changing_level(self):
        source = [
            {
                "sample": "Lower BAL",
                "collected_time": "2024-07-04 11:30",
                "reported_time": "2024-07-06 08:00",
                "organism": "Trichosporon asahii",
                "status": "Isolated",
                "colony_count": "50x10^3 CFU/mL",
            },
            {
                "sample": "Lower BAL",
                "collected_time": "2024-07-04 11:30",
                "reported_time": "2024-07-06 08:00",
                "organism": "Trichosporon asahii",
                "status": "Isolated",
                "colony_count": ">10^5 CFU/mL",
            },
        ]
        current = {
            "culture_observations": [{
                "observation_id": "CUL-001",
                "organism_name": "Trichosporon asahii",
            }],
            "organism_labels": [{
                "organism_name": "Trichosporon asahii",
                "classification": "Fungal",
                "causative_level": "Level 2",
                "applied_rules": ["CURRENT-RULE"],
                "evidence_observation_ids": ["CUL-001"],
                "key_evidence": {"quantity_tier": "Q3"},
            }],
        }
        repaired = repair_agent_provenance(current, source)
        self.assertEqual(2, len(repaired["culture_observations"]))
        label = repaired["organism_labels"][0]
        self.assertEqual("Level 2", label["causative_level"])
        self.assertEqual(["CURRENT-RULE"], label["applied_rules"])
        self.assertEqual(["CUL-001", "CUL-002"], label["evidence_observation_ids"])
        self.assertEqual("Q4", label["key_evidence"]["quantity_tier"])
        self.assertFalse(label["provenance_repair"]["semantic_level_changed"])


if __name__ == "__main__":
    unittest.main()

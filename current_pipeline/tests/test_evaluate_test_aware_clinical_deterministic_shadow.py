import unittest

from tools.evaluate_test_aware_clinical_deterministic_shadow import group_names


class TestEvaluateClinicalShadow(unittest.TestCase):
    def test_group_names_deduplicates_patient_organisms(self):
        rows = [
            {"patient_id": "9", "organism_name": "Pneumocystis jirovecii",
             "clinical_decision": "review_high_priority"},
            {"patient_id": "9", "organism_name": "Pneumocystis jirovecii",
             "clinical_decision": "review_high_priority"},
            {"patient_id": "9", "organism_name": "Candida albicans",
             "clinical_decision": "review_context_needed"},
        ]
        result = group_names(rows, lambda row: row["clinical_decision"] == "review_high_priority")
        self.assertEqual(result, {9: ["Pneumocystis jirovecii"]})


if __name__ == "__main__":
    unittest.main()

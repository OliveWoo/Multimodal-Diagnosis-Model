from __future__ import annotations

import unittest

from tools import evidence_preservation as preservation


class EvidencePreservationTests(unittest.TestCase):
    def test_p14_pjp_differential_is_preserved_as_context_not_confirmation(self) -> None:
        raw = [
            {
                "test": "Radiology / Special Examination",
                "exam_type": "CT of the chest",
                "collected_time": "2024-05-01 23:59",
                "findings": [
                    "ARDS and pneumonia at bilateral lung.",
                    "Differential diagnosis: viral pneumonitis, PJP (Pneumocystis jirovecii pneumonia).",
                ],
            }
        ]
        output = preservation.enrich_image_agent_output({"organism_hints": []}, raw)
        pjp = next(item for item in output["diagnostic_mentions"] if item["concept_name"] == "Pneumocystis jirovecii")
        self.assertEqual("differential", pjp["assertion"])
        self.assertEqual("imaging_context_only", pjp["evidence_role"])
        self.assertFalse(pjp["microbiologic_confirmation"])
        self.assertEqual([], output["organism_hints"])

    def test_negated_imaging_mention_is_not_promoted(self) -> None:
        mentions = preservation.extract_imaging_diagnostic_mentions(
            [{"exam_type": "Chest CT", "findings": ["No evidence of pulmonary tuberculosis."]}]
        )
        self.assertEqual(1, len(mentions))
        self.assertEqual("excluded", mentions[0]["assertion"])
        self.assertFalse(mentions[0]["microbiologic_confirmation"])

    def test_p9_medications_are_pre_sample_intermediate_evidence(self) -> None:
        admission = [
            {
                "test": "past medical history",
                "diagnosis": "Dermatomyositis",
                "notes": [
                    "progressive left finger swelling, add Prednisolone (2023/11-), Methotrexate (2023/12-)"
                ],
            }
        ]
        output = preservation.enrich_host_agent_output(
            None,
            raw_cbc=[{"test": "CBC", "item": "WBC", "value": "6.65", "reported_time": "2024-02-24"}],
            raw_underlying=[{"underlying_diseases": ["Dermatomyositis"]}],
            raw_admission=admission,
            sample_time="2024-02-22 10:55:00",
        )
        profile = output["structured_host_evidence"]["medication_profile"]
        self.assertEqual(2, profile["risk_medication_count"])
        self.assertEqual("intermediate", profile["host_risk_evidence_strength"])
        self.assertEqual(
            {"pre_sample_risk_medication"},
            {item["role"] for item in profile["medications"]},
        )
        self.assertEqual(
            {"prednisolone": "2023-11", "methotrexate": "2023-12"},
            {item["medication_name"]: item["start_date"] for item in profile["medications"]},
        )
        self.assertIn("systemic_corticosteroid_dose_missing", profile["data_gaps"])
        self.assertEqual("Unknown", output["host_state"]["host_vulnerability_tier"])
        self.assertFalse(output["host_state"]["expanded_candidate_policy"])
        self.assertIn("immunosuppressant_exposure", output["host_state"]["key_host_flags"])

    def test_cbc_available_is_distinct_from_differential_missing(self) -> None:
        status = preservation.cbc_data_status(
            [
                {"test": "CBC", "item": "WBC", "value": "6.65"},
                {"test": "CBC", "item": "PLT", "value": "57"},
            ]
        )
        self.assertEqual("cbc_available_differential_missing", status["status"])
        self.assertTrue(status["cbc_available"])
        self.assertFalse(status["differential_available"])

    def test_susceptibility_text_is_not_treated_as_exposure(self) -> None:
        profile = preservation.extract_host_risk_medications(
            {"culture_susceptibility": [{"result": "Methotrexate resistant MIC 8"}]},
            sample_time="2024-01-01",
        )
        self.assertEqual("susceptibility_only", profile["medications"][0]["role"])
        self.assertEqual(0, profile["risk_medication_count"])
        self.assertEqual("none", profile["host_risk_evidence_strength"])

    def test_source_index_preserves_critical_values_and_repeat_count(self) -> None:
        source = [
            {
                "organism_name": "Candida albicans",
                "specimen_type": "Lower BAL",
                "collected_time": "2024-01-01 08:00",
                "quantity": ">10^5 CFU/mL",
                "unit": "CFU/mL",
            },
            {
                "organism_name": "Candida albicans",
                "specimen_type": "Lower BAL",
                "collected_time": "2024-01-02 08:00",
                "result": "isolated",
            },
        ]
        index = preservation.build_source_evidence_index({"culture": source})
        section = index["sections"][0]
        self.assertEqual(2, section["row_count"])
        self.assertGreaterEqual(section["field_counts"]["specimen_or_source"], 2)
        self.assertGreaterEqual(section["field_counts"]["date_or_time"], 2)
        self.assertGreaterEqual(section["field_counts"]["quantity_or_result"], 2)
        self.assertEqual(
            [{"normalized_name": "candidaalbicans", "occurrence_count": 2}],
            section["repeated_organism_or_target"],
        )


if __name__ == "__main__":
    unittest.main()

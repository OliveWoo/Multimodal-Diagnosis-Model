from __future__ import annotations

import unittest
import json
import tempfile
from pathlib import Path

from tools import evidence_preservation as preservation


class EvidencePreservationTests(unittest.TestCase):
    def test_patient_mngs_sample_time_falls_back_to_per_test_inventory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            patient_dir = Path(tmp)
            path = patient_dir / "NGS_patient_9_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
            path.write_text(
                json.dumps(
                    {
                        "cases": [
                            {"collected_time": "2024-05-03 09:00"},
                            {"collected_time": "2024-05-01 10:30"},
                        ]
                    }
                ),
                encoding="utf-8",
            )

            self.assertEqual(
                preservation.patient_mngs_sample_time(patient_dir),
                "2024-05-01 10:30",
            )

    def test_category_wrapper_forwards_other_lab_to_host_enrichment(self) -> None:
        output = preservation.enrich_agent_payload(
            "cbc_other_lab",
            {},
            {
                "CBC": [{"item": "WBC", "value": "8.0", "unit": "x1000/uL"}],
                "OtherLab": [{"item": "CRP", "value": "12", "unit": "mg/L"}],
                "Underlying": [],
                "AdmissionDiagnosis": [],
            },
        )

        metadata = output["evidence_preservation"]
        self.assertTrue(metadata["raw_other_lab_available"])
        self.assertEqual(metadata["lab_observation_count"], 2)
        self.assertEqual(
            output["structured_host_evidence"]["lab_source_record_contract"]
            ["other_lab_source_record_count"],
            1,
        )

    def test_lab_observations_are_source_complete_and_flag_impossible_values(self) -> None:
        observations = preservation.lab_observations_from_sources(
            [
                {
                    "item": "WBC",
                    "item_full": "WBC (x1000/uL)",
                    "value": "15.2",
                    "reported_time": "2026-09-01",
                }
            ],
            [
                {
                    "item": "pH",
                    "item_full": "pH",
                    "value": "57",
                    "reported_time": "2026-09-01",
                }
            ],
        )
        self.assertEqual(2, len(observations))
        self.assertEqual("x1000/uL", observations[0]["unit"])
        self.assertEqual("CBC", observations[0]["source_section"])
        self.assertEqual("OtherLab", observations[1]["source_section"])
        self.assertIn("physiologically_implausible_ph", observations[1]["qc_flags"])
        self.assertFalse(observations[1]["eligible_for_acute_tiering"])

    def test_acute_profile_applies_fixed_thresholds_and_excludes_bad_ph(self) -> None:
        profile = preservation.acute_instability_profile(
            [{"item": "PLT", "item_full": "PLT (x1000/uL)", "value": "57"}],
            [
                {"item": "CRP", "item_full": "CRP (mg/dL)", "value": "16"},
                {"item": "pH", "value": "57"},
            ],
        )
        self.assertEqual("A2", profile["acute_instability_tier"])
        self.assertIn("crp_ge_150_mg_l", {item["rule"] for item in profile["triggers"]})
        self.assertEqual(1, len(profile["invalid_for_tiering"]))
        self.assertEqual("none", profile["opportunistic_host_support_effect"])

    def test_enrichment_replaces_agent_acute_tier_and_recomputes_vulnerability(self) -> None:
        output = preservation.enrich_host_agent_output(
            {
                "host_state": {
                    "immunocompromise_tier": "H0",
                    "acute_instability_tier": "A1",
                    "host_vulnerability_tier": "V1",
                    "expanded_candidate_policy": False,
                }
            },
            raw_cbc=[
                {"item": "WBC", "item_full": "WBC (x1000/uL)", "value": "8.9"},
                {"item": "PLT", "item_full": "PLT (x1000/uL)", "value": "131"},
            ],
            raw_other_lab=[
                {"item": "Creatinine", "item_full": "Creatinine (mg/dL)", "value": "161.87"}
            ],
            raw_underlying=[{"age": "40", "underlying_diseases": []}],
            raw_admission=[],
        )
        state = output["host_state"]
        self.assertEqual("A0", state["acute_instability_tier"])
        self.assertEqual("V0", state["host_vulnerability_tier"])
        self.assertFalse(state["opportunistic_host_support"])
        self.assertTrue(output["evidence_preservation"]["tier_changed_by_preservation"])

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

    def test_image_enrichment_drops_only_undated_duplicate_findings(self) -> None:
        payload = {
            "observed_findings": [
                {
                    "finding_text": "Pneumonia in both lungs",
                    "exam_type": "Chest AP",
                    "collected_time": "2024-05-01 09:00",
                    "reported_time": "2024-05-01 10:00",
                },
                {
                    "finding_text": "Pneumonia in both lungs",
                    "exam_type": "Chest AP",
                    "collected_time": "Unknown",
                    "reported_time": "Unknown",
                },
                {
                    "finding_text": "Pneumonia in both lungs",
                    "exam_type": "Chest AP",
                    "collected_time": "2024-05-02 09:00",
                    "reported_time": "2024-05-02 10:00",
                },
            ]
        }

        output = preservation.enrich_image_agent_output(payload, [])

        self.assertEqual(2, len(output["observed_findings"]))
        metadata = output["evidence_preservation"]
        self.assertEqual(1, metadata["observed_findings_removed_as_duplicates"])

    def test_empty_filmarray_and_gm_are_unknown_not_low_risk(self) -> None:
        output = preservation.enrich_agent_payload(
            "filmarray_gmtest",
            {"infection_likelihood": "None", "resistance_risk": "Low"},
            {"FilmArray": [], "GMTest": []},
        )

        self.assertEqual("Unknown", output["infection_likelihood"])
        self.assertEqual("Unknown", output["resistance_risk"])
        self.assertIn("resistance_evidence_not_available", output["data_gaps"])

    def test_empty_molecular_input_is_unknown_with_explicit_gaps(self) -> None:
        output = preservation.enrich_agent_payload(
            "molecular_microbiology",
            {"infection_likelihood": "None", "data_gaps": []},
            {"MolecularMicrobiology": []},
        )

        self.assertEqual("Unknown", output["infection_likelihood"])
        self.assertEqual([], output["molecular_observations"])
        self.assertIn("molecular_microbiology_data_not_available", output["data_gaps"])

    def test_preliminary_culture_identity_links_to_same_event_final_name(self) -> None:
        raw = [
            {
                "sample": "Blood",
                "collected_time": "2023-10-09 16:11",
                "received_time": "2023-10-10 08:55",
                "organism": "Escherichia coli",
                "status": "Isolated",
            },
            {
                "sample": "Blood",
                "collected_time": "2023-10-09 16:11",
                "received_time": "2023-10-10 08:55",
                "organism": "G(-) Bacilli, identification to follow.",
                "status": "Preliminary report",
            },
        ]
        payload = {
            "probable_pathogen_type": "Mixed",
            "culture_observations": [
                {
                    "observation_id": "CUL-001",
                    "organism_name": "Escherichia coli",
                    "specimen_type": "Blood",
                    "collected_time": "2023-10-09 16:11",
                },
                {
                    "observation_id": "CUL-002",
                    "organism_name": "G(-) Bacilli, identification to follow.",
                    "specimen_type": "Blood",
                    "collected_time": "2023-10-09 16:11",
                },
            ],
            "organism_labels": [
                {
                    "organism_name": "Escherichia coli",
                    "classification": "Bacterial",
                    "causative_level": "Level 2",
                    "evidence_observation_ids": ["CUL-001"],
                },
                {
                    "organism_name": "G(-) Bacilli, identification to follow.",
                    "classification": "Unknown",
                    "causative_level": "Level 2",
                    "evidence_observation_ids": ["CUL-002"],
                },
            ],
        }

        output = preservation.enrich_culture_agent_output(payload, raw)

        self.assertEqual(["Escherichia coli"], [x["organism_name"] for x in output["organism_labels"]])
        self.assertEqual("Bacterial", output["probable_pathogen_type"])
        self.assertEqual(
            ["CUL-001", "CUL-002"],
            output["organism_labels"][0]["evidence_observation_ids"],
        )
        self.assertEqual(1, output["evidence_preservation"]["preliminary_final_link_count"])

    def test_p9_medications_are_pre_sample_combination_evidence(self) -> None:
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
        self.assertEqual("strong", profile["host_risk_evidence_strength"])
        self.assertEqual(
            {"pre_sample_risk_medication"},
            {item["role"] for item in profile["medications"]},
        )
        self.assertEqual(
            {"prednisolone": "2023-11", "methotrexate": "2023-12"},
            {item["medication_name"]: item["start_date"] for item in profile["medications"]},
        )
        self.assertIn("systemic_corticosteroid_dose_missing", profile["data_gaps"])
        self.assertEqual("H3", output["host_state"]["immunocompromise_tier"])
        self.assertEqual("V3", output["host_state"]["host_vulnerability_tier"])
        self.assertTrue(output["host_state"]["expanded_candidate_policy"])
        self.assertTrue(output["host_state"]["opportunistic_host_support"])
        self.assertIn("immunosuppressant_exposure", output["host_state"]["key_host_flags"])

    def test_long_term_steroid_with_unknown_dose_is_context_not_host_support(self) -> None:
        output = preservation.enrich_host_agent_output(
            {"host_state": {"immunocompromise_tier": "H1"}},
            raw_cbc=[{"item": "WBC", "value": "8.0"}],
            raw_underlying=[
                {
                    "age": "99",
                    "underlying_diseases": [
                        "Suspect asthma, under prednisolone since 2023/10"
                    ],
                }
            ],
            raw_admission=[],
            sample_time="2024-05-09 13:30:00",
        )
        state = output["host_state"]
        self.assertEqual("H2", state["immunocompromise_tier"])
        self.assertFalse(state["opportunistic_host_support"])
        self.assertEqual("insufficient", state["opportunistic_host_support_status"])

    def test_parathyroid_autotransplant_does_not_trigger_transplant_h3(self) -> None:
        profile = preservation.deterministic_host_risk_profile(
            [
                {
                    "age": "70",
                    "underlying_diseases": [
                        "post parathyroidectomy and parathyroid autotransplantation"
                    ],
                }
            ],
            [],
        )
        self.assertEqual("H1", profile["immunocompromise_tier"])
        self.assertNotIn(
            "transplant_history", {item["rule"] for item in profile["all_triggers"]}
        )

    def test_medication_dose_and_end_date_stay_with_correct_drug(self) -> None:
        profile = preservation.extract_host_risk_medications(
            {
                "underlying": [
                    {
                        "note": (
                            "Prednisolone 10-15mg/day (2023/11/29-) and "
                            "Methotrexate (2023/12/27-)"
                        )
                    },
                    {
                        "note": (
                            "Decitabine 20 mg/m2 on "
                            "2024/01/17-2024/01/21"
                        )
                    },
                ]
            },
            sample_time="2024-02-25 11:30:00",
        )
        by_name = {item["medication_name"]: item for item in profile["medications"]}
        self.assertEqual("10-15mg/day", by_name["prednisolone"]["dose"])
        self.assertEqual("Unknown", by_name["methotrexate"]["dose"])
        self.assertEqual("2024-01-21", by_name["decitabine"]["end_date"])
        self.assertEqual("ended_before_sample", by_name["decitabine"]["temporal_status"])

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

    def test_seg_and_band_are_recognized_as_cbc_differential(self) -> None:
        status = preservation.cbc_data_status(
            [
                {"test": "CBC", "item": "WBC", "value": "8.1"},
                {"test": "CBC", "item": "Seg", "value": "72", "unit": "%"},
                {"test": "CBC", "item": "Band", "value": "3", "unit": "%"},
            ]
        )
        self.assertEqual("cbc_and_differential_available_alc_missing", status["status"])
        self.assertTrue(status["differential_available"])
        self.assertFalse(status["absolute_lymphocyte_count_available"])

    def test_wbc_and_lymphocyte_percent_mark_alc_as_derivable_not_observed(self) -> None:
        status = preservation.cbc_data_status(
            [
                {"test": "CBC", "item": "WBC", "value": "4.0", "unit": "10^3/uL"},
                {"test": "CBC", "item": "Lymphocyte %", "value": "10", "unit": "%"},
            ]
        )
        self.assertTrue(status["alc_potentially_derivable"])
        self.assertFalse(status["absolute_lymphocyte_count_available"])

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

    def test_source_record_index_preserves_complete_rows_and_relationships(self) -> None:
        source = [
            {
                "organism_name": "Candida albicans",
                "specimen_type": "BALF",
                "collected_time": "2024-01-01 08:00",
                "result": "isolated",
                "quantity": ">10^5",
                "unit": "CFU/mL",
                "free_text": "mixed yeast morphology",
            },
            {
                "organism_name": "Candida albicans",
                "specimen_type": "Blood",
                "collected_time": "2024-01-03 08:00",
                "result": "not detected",
            },
        ]
        index = preservation.build_source_record_index({"culture": source})
        section = index["sections"][0]
        self.assertTrue(index["lossless_for_normalized_records"])
        self.assertEqual(2, section["record_count"])
        self.assertEqual("culture[0]", section["records"][0]["record_path"])
        self.assertEqual(source[0], section["records"][0]["record"])
        self.assertEqual(source[1], section["records"][1]["record"])

    def test_duplicate_medication_mentions_do_not_inflate_distinct_count(self) -> None:
        profile = preservation.extract_host_risk_medications(
            {
                "underlying": [{"note": "Prednisolone since 2023/01-"}],
                "admission": [{"note": "Prednisolone since 2023/01-"}],
            },
            sample_time="2024-01-01",
        )
        self.assertEqual(2, profile["risk_medication_observation_count"])
        self.assertEqual(1, profile["risk_medication_count"])
        self.assertEqual(["prednisolone"], profile["distinct_risk_medications"])


if __name__ == "__main__":
    unittest.main()

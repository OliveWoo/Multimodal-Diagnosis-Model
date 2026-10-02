from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.build_multi_assay_mngs_evidence import build_patient


class MultiAssayMngsEvidenceTest(unittest.TestCase):
    def test_build_patient_keeps_tests_separate_and_builds_two_unions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "NGS_patient_1_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
            source.write_text(
                json.dumps(
                    {
                        "patient_id": 1,
                        "cases": [
                            {
                                "case_review_id": "case-a",
                                "specimen_code": "case-a",
                                "collected_time": "2024-01-01 12:00:00",
                                "specimen_site": "BALF",
                                "tests": [
                                    {
                                        "seq_id": "dna-1",
                                        "nucleic_type": "DNA",
                                        "parameters": {"input_reads": 1_000_000, "host_ratio": 90, "condition": "-"},
                                        "raw_result_status": {
                                            key: "available" for key in ("bacterial", "viral", "fungal", "others")
                                        },
                                    },
                                    {
                                        "seq_id": "rna-1",
                                        "nucleic_type": "RNA",
                                        "parameters": {"input_reads": 2_000_000, "host_ratio": 80, "condition": "adapter 10x"},
                                        "raw_result_status": {
                                            key: "available" for key in ("bacterial", "viral", "fungal", "others")
                                        },
                                    },
                                ],
                            }
                        ],
                        "organisms": [
                            {
                                "category": "fungal",
                                "name": "Aspergillus testii",
                                "measurements": [
                                    {
                                        "case_review_id": "case-a",
                                        "seq_id": "dna-1",
                                        "reads": 10,
                                        "detected": True,
                                        "selected": True,
                                        "selection_tags": ["Code_NTC=00", "Code_RK_NTC=RK0"],
                                        "raw_available": True,
                                    },
                                    {
                                        "case_review_id": "case-a",
                                        "seq_id": "rna-1",
                                        "reads": 20,
                                        "detected": True,
                                        "selected": True,
                                        "selection_tags": ["Code_NTC=00"],
                                        "raw_available": True,
                                    },
                                ],
                            },
                            {
                                "category": "bacterial",
                                "name": "Filtered bacterium",
                                "measurements": [
                                    {
                                        "case_review_id": "case-a",
                                        "seq_id": "dna-1",
                                        "reads": 5,
                                        "detected": True,
                                        "selected": False,
                                        "selection_tags": [],
                                        "raw_available": True,
                                    }
                                ],
                            },
                        ],
                    },
                    ensure_ascii=False,
                ),
                encoding="utf-8",
            )

            patient, tests, candidates = build_patient(source)

        self.assertEqual(len(tests), 2)
        self.assertTrue(all(row["qc_status"] == "evaluable" for row in tests))
        selected = patient["cases"][0]["candidate_union"]["selected_union"]
        filtered = patient["cases"][0]["candidate_union"]["detected_but_filtered_union"]
        self.assertEqual([item["organism_name"] for item in selected], ["Aspergillus testii"])
        self.assertEqual([item["organism_name"] for item in filtered], ["Filtered bacterium"])
        self.assertTrue(selected[0]["cross_molecule_support"])
        self.assertIsNone(selected[0]["reads_sum_across_tests"])
        self.assertEqual([item["reads"] for item in selected[0]["observations"]], [10.0, 20.0])
        self.assertEqual(selected[0]["observations"][0]["rpm_total"], 10.0)
        self.assertEqual(len(candidates), 2)

    def test_missing_raw_category_is_not_treated_as_negative(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "NGS_patient_2_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
            source.write_text(
                json.dumps(
                    {
                        "patient_id": 2,
                        "cases": [
                            {
                                "case_review_id": "case-b",
                                "specimen_code": "case-b",
                                "specimen_site": "BALF",
                                "tests": [
                                    {
                                        "seq_id": "dna-2",
                                        "nucleic_type": "DNA",
                                        "parameters": {"input_reads": 1000, "host_ratio": 50, "condition": "-"},
                                        "raw_result_status": {
                                            "bacterial": "available",
                                            "viral": "missing",
                                            "fungal": "available",
                                            "others": "available",
                                        },
                                    }
                                ],
                            }
                        ],
                        "organisms": [
                            {
                                "category": "viral",
                                "name": "Example virus",
                                "measurements": [
                                    {
                                        "case_review_id": "case-b",
                                        "seq_id": "dna-2",
                                        "reads": None,
                                        "detected": False,
                                        "selected": False,
                                        "selection_tags": [],
                                        "raw_available": False,
                                    }
                                ],
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            patient, tests, candidates = build_patient(source)

        self.assertEqual(tests[0]["qc_status"], "partial_evaluable")
        self.assertEqual(candidates, [])
        self.assertEqual(patient["cases"][0]["candidate_union"]["selected_union"], [])

    def test_zero_read_detected_flag_is_audit_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "NGS_patient_3_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
            source.write_text(
                json.dumps({
                    "patient_id": 3,
                    "cases": [{
                        "case_review_id": "case-c",
                        "specimen_code": "case-c",
                        "specimen_site": "BALF",
                        "tests": [
                            {"seq_id": seq_id, "nucleic_type": "DNA", "parameters": {"input_reads": 1000, "host_ratio": 50},
                             "raw_result_status": {key: "available" for key in ("bacterial", "viral", "fungal", "others")}}
                            for seq_id in ("dna-1", "dna-2")
                        ],
                    }],
                    "organisms": [
                        {"category": "bacterial", "name": "Example bacterium", "measurements": [
                            {"case_review_id": "case-c", "seq_id": "dna-1", "reads": 0, "detected": True,
                             "selected": False, "selection_tags": [], "raw_available": True},
                            {"case_review_id": "case-c", "seq_id": "dna-2", "reads": 5, "detected": True,
                             "selected": True, "selection_tags": ["Code_NTC=00"], "raw_available": True},
                        ]},
                        {"category": "bacterial", "name": "Zero only bacterium", "measurements": [
                            {"case_review_id": "case-c", "seq_id": "dna-1", "reads": 0, "detected": True,
                             "selected": False, "selection_tags": [], "raw_available": True},
                        ]},
                    ],
                }),
                encoding="utf-8",
            )

            patient, _, candidates = build_patient(source)

        self.assertEqual([row["organism_name"] for row in candidates], ["Example bacterium"])
        selected = patient["cases"][0]["candidate_union"]["selected_union"]
        self.assertEqual([row["seq_id"] for row in selected[0]["observations"]], ["dna-2"])
        self.assertEqual(selected[0]["detected_test_count"], 1)
        self.assertEqual(selected[0]["max_reads_single_test"], 5)
        zero_audit = [row for row in patient["observations"] if row["reads"] == 0]
        self.assertEqual(len(zero_audit), 2)
        self.assertTrue(all(row["source_case_detected"] for row in zero_audit))
        self.assertTrue(all(not row["observed_in_current_test"] for row in zero_audit))
        self.assertTrue(all(row["selection_status"] == "case_manifest_only" for row in zero_audit))
        self.assertTrue(
            all("case_manifest_presence_without_current_test_reads" in row["observation_warnings"] for row in zero_audit)
        )

    def test_source_missing_fields_have_explicit_nonnegative_semantics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "NGS_patient_4_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
            source.write_text(
                json.dumps({
                    "patient_id": 4,
                    "cases": [{
                        "case_review_id": "case-d",
                        "specimen_code": "case-d",
                        "collected_time": "2024-01-01 12:00:00",
                        "specimen_site": "BALF",
                        "tests": [{
                            "seq_id": "rna-missing",
                            "nucleic_type": "RNA",
                            "report_time": "-",
                            "parameters": {"input_reads": None, "host_ratio": None},
                            "raw_result_status": {
                                key: "missing" for key in ("bacterial", "viral", "fungal", "others")
                            },
                        }],
                    }],
                    "organisms": [],
                }),
                encoding="utf-8",
            )

            patient, tests, candidates = build_patient(source)

        self.assertEqual(candidates, [])
        self.assertEqual(patient["schema_version"], "multi_assay_mngs_evidence.v2")
        self.assertEqual(tests[0]["normalization_status"], "source_denominators_missing")
        self.assertEqual(tests[0]["normalization_missingness"], "confirmed_absent_in_source")
        self.assertIn("source_missing_input_reads", tests[0]["qc_issues"])
        self.assertIn("source_missing_host_ratio", tests[0]["qc_issues"])
        self.assertEqual(tests[0]["report_time_status"], "source_missing_confirmed")
        self.assertEqual(tests[0]["decision_time_anchor"], "specimen_collected_time")
        self.assertEqual(
            tests[0]["raw_data_interpretation"],
            "no_corresponding_raw_file_unevaluable_not_negative",
        )


if __name__ == "__main__":
    unittest.main()

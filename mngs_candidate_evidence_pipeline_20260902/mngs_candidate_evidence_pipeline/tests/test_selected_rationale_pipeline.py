from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from build_candidate_window_evidence import TimedRecord, find_mngs_anchor
from generate_llm_rationales import validate_output


class SelectedRationalePipelineTest(unittest.TestCase):
    def run_command(self, *args: str) -> None:
        subprocess.run([sys.executable, *args], cwd=ROOT, check=True)

    def test_extract_and_build_both_views(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            temp = Path(temporary)
            selected = temp / "selected"
            output = temp / "output"
            self.run_command(
                "extract_selected_pathogens.py",
                "--input-root",
                "examples/upstream",
                "--output-root",
                str(selected),
                "--hospital",
                "Synthetic Hospital",
                "--dataset",
                "synthetic_demo",
                "--pipeline-version",
                "synthetic_v1",
                "--pretty",
            )
            self.run_command(
                "run_pipeline.py",
                "--input-root",
                "examples/input",
                "--mngs-root",
                "examples/mngs",
                "--selected-root",
                str(selected),
                "--candidate-source",
                "upstream_all",
                "--output-root",
                str(output),
                "--anchor-source",
                "mngs",
                "--window-mode",
                "date",
                "--mode",
                "both",
                "--rationale-views",
                "--pretty",
            )

            with (output / "rationale_views" / "human_readable_selected_pathogens.csv").open(
                "r", encoding="utf-8-sig", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["選擇菌種"], "Klebsiella pneumoniae")
            self.assertEqual(rows[0]["檢驗編號"], "SYNTHETIC-001")

            readable_chain = json.loads(
                (
                    output
                    / "rationale_views"
                    / "traceable_reasoning_chains"
                    / "NGS_patient_1_explainable_reasoning_zh.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(
                readable_chain["一眼看懂"]["最終入選病原"],
                ["Klebsiella pneumoniae"],
            )
            self.assertIn("步驟2至5_逐候選病原判讀", readable_chain)

            chain = json.loads(
                (
                    output
                    / "rationale_views"
                    / "technical_traces"
                    / "NGS_patient_1_technical_trace.json"
                ).read_text(encoding="utf-8")
            )
            statuses = {
                item["organism"]: item["selection_status"]
                for item in chain["candidate_diagnoses"]
            }
            self.assertEqual(statuses["Klebsiella pneumoniae"], "picked")
            self.assertEqual(statuses["Candida albicans"], "not_picked")
            self.assertEqual(
                chain["traceability_status"]["external_literature_evidence"],
                "pending_ober",
            )

    def test_mngs_anchor_never_crosses_patient_ids(self) -> None:
        patient_record = TimedRecord(
            record_type="culture",
            source_file="patient_1.json",
            index=0,
            record={},
            record_time=datetime(2024, 1, 10, 8, 0),
            time_field="collected_time",
            microbe_name=None,
            normalized_name=None,
            result=None,
            is_positive=False,
            is_candidate_eligible=False,
        )
        other_patient_mngs = {
            "source_patient_id": 2,
            "record_time": datetime(2024, 1, 10, 8, 30),
        }
        anchor, record, source = find_mngs_anchor(
            [patient_record], [other_patient_mngs], patient_id=1
        )
        self.assertIsNone(anchor)
        self.assertIsNone(record)
        self.assertEqual(source, "none")

    def test_patient_result_explicitly_keeps_no_pick_patient(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            temp = Path(temporary)
            candidate_root = temp / "candidate"
            output_root = temp / "output"
            candidate_root.mkdir()
            bundle = {
                "patient_id": 2,
                "selected_input": {"hospital": "Synthetic Hospital"},
                "candidates": [
                    {
                        "name": "Candida albicans",
                        "selection_status": "not_picked",
                        "selection_provenance": {
                            "decision_reasons": ["Insufficient evidence"]
                        },
                    }
                ],
            }
            (candidate_root / "NGS_patient_2_candidate_window_2d_all_time_evidence.json").write_text(
                json.dumps(bundle), encoding="utf-8"
            )
            self.run_command(
                "build_rationale_views.py",
                "--candidate-root",
                str(candidate_root),
                "--output-root",
                str(output_root),
                "--pretty",
            )
            with (output_root / "human_readable_patient_results.csv").open(
                "r", encoding="utf-8-sig", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["選擇菌種"], "無入選病原")
            self.assertEqual(rows[0]["選擇狀態"], "none_picked")
            readable_chain = json.loads(
                (
                    output_root
                    / "traceable_reasoning_chains"
                    / "NGS_patient_2_explainable_reasoning_zh.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(readable_chain["一眼看懂"]["最終入選病原"], ["無入選病原"])

    def test_llm_reason_validation_rejects_changed_selection(self) -> None:
        payload = {
            "patient_id": "1",
            "frozen_decision": {"no_selected_pathogen": False},
            "candidates_available_to_reason_generator": [
                {
                    "organism": "Klebsiella pneumoniae",
                    "frozen_decision": "選入",
                    "frozen_rank": 1,
                    "evidence": [{"evidence_id": "C01-E01"}],
                }
            ],
        }
        changed = {
            "patient_id": "1",
            "no_selected_pathogen": False,
            "patient_result_summary": "Changed",
            "patient_level_evidence_ids": [],
            "selected_pathogen_explanations": [
                {
                    "organism": "Escherichia coli",
                    "rank": 1,
                    "evidence_ids": ["C01-E01"],
                }
            ],
        }
        with self.assertRaisesRegex(ValueError, "selected pathogens/ranks changed"):
            validate_output(changed, payload, {"C01-E01": "source.json#0"})

        correct_item = {"organism": "Klebsiella pneumoniae", "rank": 1, "evidence_ids": ["C01-E01"]}
        changed["selected_pathogen_explanations"] = [correct_item]
        validate_output(changed, payload, {"C01-E01": "source.json#0"})
        changed["selected_pathogen_explanations"] = [correct_item, dict(correct_item)]
        with self.assertRaisesRegex(ValueError, "duplicate organisms"):
            validate_output(changed, payload, {"C01-E01": "source.json#0"})

        changed["selected_pathogen_explanations"] = [
            {**correct_item, "limitations": [{"text": "another candidate", "evidence_ids": ["C02-E01"]}]}
        ]
        with self.assertRaisesRegex(ValueError, "another candidate's evidence"):
            validate_output(changed, payload, {"C01-E01": "source.json#0", "C02-E01": "source.json#1"})

        changed["selected_pathogen_explanations"] = [correct_item]
        changed["patient_level_evidence_ids"] = ["C02-E01"]
        validate_output(changed, payload, {"C01-E01": "source.json#0", "C02-E01": "source.json#1"})


if __name__ == "__main__":
    unittest.main()

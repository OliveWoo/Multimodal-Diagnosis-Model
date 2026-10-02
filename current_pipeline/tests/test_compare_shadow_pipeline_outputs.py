from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.compare_shadow_pipeline_outputs import compare_patient


def _write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _path(patient_dir: Path, suffix: str) -> Path:
    return patient_dir / "summary_outputs" / f"NGS_patient_1_{suffix}.json"


class CompareShadowPipelineOutputsTests(unittest.TestCase):
    def test_compare_patient_reports_membership_and_level_changes(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            patient_dir = Path(temp_dir) / "NGS_patient_1_json"
            baseline_scorer = {
                "best_available_summary": {
                    "picked_pathogens": [
                        {"organism_name": "Alpha bacterium", "basis_level": "M1"}
                    ]
                },
                "pathogen_candidates": [
                    {
                        "organism_name": "Alpha bacterium",
                        "integrated_causative_level": "M1",
                    },
                    {
                        "organism_name": "Beta bacterium",
                        "integrated_causative_level": "M3",
                    },
                ],
            }
            candidate_scorer = {
                "best_available_summary": {
                    "picked_pathogens": [
                        {"organism_name": "Beta bacterium", "basis_level": "M2"}
                    ]
                },
                "pathogen_candidates": [
                    {
                        "organism_name": "Alpha bacterium",
                        "integrated_causative_level": "M2",
                    },
                    {
                        "organism_name": "Beta bacterium",
                        "integrated_causative_level": "M2",
                    },
                ],
            }
            baseline_queue = {
                "review_queue": [{"organism_name": "Gamma bacterium"}]
            }
            candidate_queue = {
                "review_queue": [{"organism_name": "Delta bacterium"}]
            }
            _write(_path(patient_dir, "baseline_scorer"), baseline_scorer)
            _write(_path(patient_dir, "candidate_scorer"), candidate_scorer)
            _write(_path(patient_dir, "baseline_queue"), baseline_queue)
            _write(_path(patient_dir, "candidate_queue"), candidate_queue)

            result = compare_patient(
                patient_dir,
                "baseline_scorer",
                "candidate_scorer",
                "baseline_queue",
                "candidate_queue",
            )

            self.assertEqual(result["missing_outputs"], [])
            self.assertEqual(result["picked_removed"], ["Alpha bacterium"])
            self.assertEqual(result["picked_added"], ["Beta bacterium"])
            self.assertEqual(result["candidate_removed"], [])
            self.assertEqual(result["candidate_added"], [])
            self.assertEqual(
                result["level_changes"],
                [
                    {
                        "organism": "Alpha bacterium",
                        "before": "M1",
                        "after": "M2",
                    },
                    {
                        "organism": "Beta bacterium",
                        "before": "M3",
                        "after": "M2",
                    },
                ],
            )
            self.assertEqual(result["queue_removed"], ["Gamma bacterium"])
            self.assertEqual(result["queue_added"], ["Delta bacterium"])

    def test_compare_patient_reports_missing_files(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            patient_dir = Path(temp_dir) / "NGS_patient_1_json"
            patient_dir.mkdir()

            result = compare_patient(
                patient_dir,
                "baseline_scorer",
                "candidate_scorer",
                "baseline_queue",
                "candidate_queue",
            )

            self.assertEqual(
                result["missing_outputs"],
                [
                    "baseline_scorer",
                    "candidate_scorer",
                    "baseline_queue",
                    "candidate_queue",
                ],
            )

    def test_compare_patient_reports_candidate_membership_changes(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            patient_dir = Path(temp_dir) / "NGS_patient_1_json"
            baseline_scorer = {
                "best_available_summary": {"picked_pathogens": []},
                "pathogen_candidates": [{"organism_name": "Alpha bacterium"}],
            }
            candidate_scorer = {
                "best_available_summary": {"picked_pathogens": []},
                "pathogen_candidates": [{"organism_name": "Beta bacterium"}],
            }
            queue = {"review_queue": []}
            _write(_path(patient_dir, "baseline_scorer"), baseline_scorer)
            _write(_path(patient_dir, "candidate_scorer"), candidate_scorer)
            _write(_path(patient_dir, "baseline_queue"), queue)
            _write(_path(patient_dir, "candidate_queue"), queue)

            result = compare_patient(
                patient_dir,
                "baseline_scorer",
                "candidate_scorer",
                "baseline_queue",
                "candidate_queue",
            )

            self.assertEqual(result["candidate_removed"], ["Alpha bacterium"])
            self.assertEqual(result["candidate_added"], ["Beta bacterium"])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools import run_deterministic_v20_pipeline as pipeline


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")


def make_flat_patient(root: Path, patient_id: str = "54", *, both_sources: bool = False) -> None:
    prefix = f"NGS_patient_{patient_id}"
    grouped = [
        {
            "specimen_code": f"SYN-{patient_id}",
            "collected_time": "2025-01-02",
            "specimen_site": "BAL",
            "pathogens": {
                "bacterial": [
                    {
                        "name": "Synthetic bacterium",
                        "reads": 120,
                        "rk_ntc": ["Code_NTC=00", "Code_RK_NTC=RK0"],
                    }
                ],
                "viral": [],
                "fungal": [],
                "others": [],
            },
        }
    ]
    write_json(root / f"{prefix}_mNGS_grouped.json", grouped)
    if both_sources:
        write_json(root / f"{prefix}_all_RK_NTC_microbes.json", grouped)
    for suffix, payload in {
        "admission_diagnosis": {},
        "CBC": {},
        "culture": [],
        "filmarray": [],
        "gm_test": [],
        "image": [],
        "other_lab": {},
        "underlying": {},
    }.items():
        write_json(root / f"{prefix}_{suffix}.json", payload)


class DeterministicV20PipelineTests(unittest.TestCase):
    def test_flat_standardized_json_runs_to_v20_final_results(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "source"
            source.mkdir()
            make_flat_patient(source)
            output = root / "run"

            manifest = pipeline.run_pipeline(
                input_root=source,
                output_root=output,
                source_kind="auto",
            )

            self.assertEqual(1, manifest["patient_count"])
            self.assertFalse(manifest["network_or_llm_used"])
            self.assertEqual("mngs_grouped", manifest["resolved_source_modes"]["54"])
            final = json.loads((output / "final_results.json").read_text(encoding="utf-8"))
            patient = final["patients"]["54"]
            self.assertEqual(1, patient["picked_count"])
            self.assertEqual("Synthetic bacterium", patient["picked_pathogens"][0]["organism_name"])
            self.assertEqual(
                "mngs_deterministic_level_only_v20_central_taxonomy_profile",
                patient["rule_version"],
            )
            patient_root = output / "patients" / "NGS_patient_54_json"
            self.assertTrue((patient_root / "agent_outputs").is_dir())
            self.assertTrue((patient_root / "summary_outputs").is_dir())

    def test_patient_directory_without_agent_outputs_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "source" / "NGS_patient_8_json"
            source.mkdir(parents=True)
            make_flat_patient(source, "8")
            output = root / "run"

            pipeline.run_pipeline(
                input_root=source.parent,
                output_root=output,
            )

            self.assertTrue((output / "run_manifest.json").is_file())
            self.assertTrue(
                (
                    output
                    / "patients"
                    / "NGS_patient_8_json"
                    / "summary_outputs"
                    / "NGS_patient_8_final_summary_with_filmarray_deterministic.json"
                ).is_file()
            )

    def test_auto_refuses_ambiguous_mngs_source(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "source"
            source.mkdir()
            make_flat_patient(source, both_sources=True)

            with self.assertRaisesRegex(pipeline.PipelineInputError, "both mNGS_grouped and all_RK_NTC"):
                pipeline.run_pipeline(
                    input_root=source,
                    output_root=root / "run",
                    source_kind="auto",
                )
            failure = json.loads((root / "run" / "run_failed.json").read_text(encoding="utf-8"))
            self.assertEqual("PipelineInputError", failure["error_type"])

    def test_existing_output_is_never_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "source"
            source.mkdir()
            make_flat_patient(source)
            output = root / "run"
            output.mkdir()
            marker = output / "keep.txt"
            marker.write_text("keep", encoding="utf-8")

            with self.assertRaises(FileExistsError):
                pipeline.run_pipeline(input_root=source, output_root=output)
            self.assertEqual("keep", marker.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()

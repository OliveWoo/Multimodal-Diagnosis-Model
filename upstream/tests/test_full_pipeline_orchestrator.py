from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import run_full_pipeline as pipeline


class FullPipelineOrchestratorTests(unittest.TestCase):
    def test_template_builds_complete_plan_without_creating_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "planned_run"
            plan = pipeline.build_plan(
                pipeline.REPO_ROOT / "examples" / "full_pipeline_config.template.json",
                output,
            )

            self.assertFalse(output.exists())
            self.assertEqual(("KH", "two_hospital"), tuple(c.name for c in plan.cohorts))
            names = [stage.name for stage in plan.stages]
            self.assertIn("KH:deterministic_v20", names)
            self.assertIn("audit_initial_rationales", names)
            self.assertIn("r5_final_decisions", names)
            self.assertIn("assemble_r5_rationales", names)
            self.assertEqual("audit_delivery", names[-1])
            self.assertTrue(any(stage.external for stage in plan.stages))

    def test_overlapping_patient_ids_are_cohort_qualified_for_r5(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = pipeline.REPO_ROOT / "examples" / "frontend" / "cohorts" / "two_hospital"
            for cohort in ("alpha", "beta"):
                shutil.copytree(source, root / cohort)
            config = {
                "schema_version": pipeline.CONFIG_SCHEMA,
                "cohorts": {
                    name: {
                        "input_root": name,
                        "hospital": name,
                        "dataset": f"{name}_dataset",
                        "source_kind": "mngs_grouped",
                    }
                    for name in ("alpha", "beta")
                },
            }
            config_file = root / "config.json"
            config_file.write_text(json.dumps(config), encoding="utf-8")

            plan = pipeline.build_plan(config_file, root / "output")
            assemble = next(stage for stage in plan.stages if stage.name == "assemble_r5_rationales")
            command = list(assemble.command or ())

            self.assertEqual(("54",), plan.cohorts[0].patient_ids)
            self.assertEqual(("54",), plan.cohorts[1].patient_ids)
            self.assertIn("alpha=" + str(root / "output" / "10_r5_rationale_runs" / "alpha"), command)
            self.assertIn("beta=" + str(root / "output" / "10_r5_rationale_runs" / "beta"), command)

    def test_flatten_merged_is_copy_only_and_hash_preserving(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            patient_root = root / "patients"
            source = (
                patient_root
                / "NGS_patient_7_json"
                / "summary_outputs"
                / f"NGS_patient_7_{pipeline.MERGED_SUFFIX}.json"
            )
            source.parent.mkdir(parents=True)
            source.write_text('{"patient_id":"7"}\n', encoding="utf-8")
            output = root / "merged"

            pipeline.flatten_merged(
                {
                    "patient_root": str(patient_root),
                    "output_root": str(output),
                    "patient_ids": ["7"],
                }
            )

            copied = output / source.name
            self.assertEqual(source.read_bytes(), copied.read_bytes())
            self.assertEqual(pipeline.sha256_file(source), pipeline.sha256_file(copied))


if __name__ == "__main__":
    unittest.main()

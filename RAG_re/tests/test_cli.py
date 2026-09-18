import copy
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

from rag_re.cli import _valid_existing_output, command_batch, command_sweep
from rag_re.config import DEFAULT_CONFIG, config_hash
from rag_re.io_utils import sha256_file
from rag_re.versioning import PIPELINE_FINGERPRINT


class CliTests(unittest.TestCase):
    def test_resume_retries_partial_artifact(self):
        config = copy.deepcopy(DEFAULT_CONFIG)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            input_path = root / "input.json"
            output_path = root / "output.json"
            input_path.write_text("{}", encoding="utf-8")
            payload = {
                "schema_version": "rag_re.output.v1",
                "run_status": "partial_with_errors",
                "pipeline_fingerprint": PIPELINE_FINGERPRINT,
                "input_meta": {"input_sha256": sha256_file(input_path)},
                "config_hash": config_hash(config),
                "execution": {"skip_literature": False},
            }
            output_path.write_text(json.dumps(payload), encoding="utf-8")
            self.assertFalse(
                _valid_existing_output(
                    output_path,
                    input_path=input_path,
                    config=config,
                    skip_literature=False,
                )
            )
            payload["run_status"] = "complete"
            output_path.write_text(json.dumps(payload), encoding="utf-8")
            self.assertTrue(
                _valid_existing_output(
                    output_path,
                    input_path=input_path,
                    config=config,
                    skip_literature=False,
                )
            )

    def test_sweep_emits_no_confirmatory_threshold_when_floor_fails(self):
        source = {
            "schema_version": "rag_re.output.v1",
            "run_status": "complete",
            "patient_id": "1",
            "execution": {
                "literature_evidence_available": True,
                "literature_attempted": True,
            },
            "config_snapshot": copy.deepcopy(DEFAULT_CONFIG),
            "candidates": [
                {
                    "literature_evidence": {
                        "judgments": [{"status": "ok"}]
                    }
                }
            ],
        }
        gold = {
            "development_only": True,
            "cases": [{"patient_id": "1", "pathogens": ["A"], "split": "development"}],
        }

        class FakeEngine:
            def __init__(self, config, skip_literature):
                self.config = config

            def replay(self, payload):
                return payload

        def fake_evaluate(outputs, loaded_gold, bootstrap_iterations):
            return {
                "experiments": {
                    "E0_OR_C_OR_A_AND_B": {
                        "metrics": {
                            "agreement_precision": 0.8,
                            "end_to_end_recall": 0.5,
                            "candidate_pool_recall_ceiling": 0.5,
                            "f2": 0.55,
                        }
                    }
                }
            }

        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "sweep.json"
            args = Namespace(
                predictions="unused",
                gold="unused",
                output=str(output),
                thresholds="1,2",
                experiment="E0_OR_C_OR_A_AND_B",
                recall_floor=0.9,
                split=None,
                allow_nondevelopment_sweep=False,
                config=None,
                env_file=None,
            )
            with (
                patch("rag_re.cli.load_outputs", return_value=[source]),
                patch("rag_re.cli.load_gold", return_value=gold),
                patch("rag_re.cli._load_runtime_config", return_value=copy.deepcopy(DEFAULT_CONFIG)),
                patch("rag_re.cli.RagReEngine", FakeEngine),
                patch("rag_re.cli.evaluate_outputs", side_effect=fake_evaluate),
                patch("builtins.print"),
            ):
                code = command_sweep(args)
            payload = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(code, 2)
        self.assertEqual(payload["selection_status"], "failed_no_threshold_meets_recall_floor")
        self.assertIsNone(payload["selected_threshold"])

    def test_batch_rejects_output_or_cache_inside_input_tree(self):
        with tempfile.TemporaryDirectory() as temp:
            input_dir = Path(temp) / "inputs"
            input_dir.mkdir()
            args = Namespace(
                input_dir=str(input_dir),
                output_dir=str(input_dir / "outputs"),
                cache_dir=None,
                pattern="*.json",
                recursive=True,
                config=None,
                env_file=None,
                skip_literature=True,
                resume=False,
                fail_fast=False,
            )
            with self.assertRaisesRegex(ValueError, "output directory"):
                command_batch(args)


if __name__ == "__main__":
    unittest.main()

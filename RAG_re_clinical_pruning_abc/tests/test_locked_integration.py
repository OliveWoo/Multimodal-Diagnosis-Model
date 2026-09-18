from __future__ import annotations

import argparse
import json
import tempfile
import unittest
from pathlib import Path

from rag_re_clinical_pruning_abc.cli import replay_batch
from rag_re_clinical_pruning_abc.evaluation import EvaluationError, evaluate_directory


WORKSPACE = Path(__file__).resolve().parents[2]
PROJECT = WORKSPACE / "RAG_re_clinical_pruning_abc"


class LockedIntegrationTests(unittest.TestCase):
    def _replay(self, output: Path) -> None:
        args = argparse.Namespace(
            pruning_input=str(
                WORKSPACE
                / "RAG_re_clinical_pruning"
                / "runs"
                / "pruning_v1_33_20260811_final_locked"
            ),
            clinical_input=str(
                WORKSPACE / "RAG_re_clinical" / "runs" / "clinical_v1_33_20260811_locked"
            ),
            raw_input=str(
                WORKSPACE
                / "RAG_re"
                / "runs"
                / "KH_0728_2Days_latest_20260811_luna_precision_pruning_replay"
            ),
            output=str(output),
            config=str(PROJECT / "config" / "post_l5_abc_policy_v1.json"),
        )
        self.assertEqual(replay_batch(args), 0)

    def test_locked_33_replay_and_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "abc"
            self._replay(output)
            report = evaluate_directory(
                output,
                WORKSPACE / "RAG_re" / "gold" / "KMUH_clinical_pathogen_gold_20260811.json",
                PROJECT / "config" / "evaluation_protocol_v1.json",
                bootstrap_replicates=0,
            )
            self.assertEqual(report["artifacts"]["artifact_count"], 33)
            self.assertTrue(report["artifacts"]["batch_manifest_validation"]["validated"])
            primary = report["scopes"][report["primary_scope"]]["arms"]
            sensitivity = report["scopes"][report["sensitivity_scope"]]["arms"]
            self.assertEqual(
                (primary["L5_FROZEN"]["auto"]["tp"], primary["L5_FROZEN"]["auto"]["fp"]),
                (42, 31),
            )
            self.assertAlmostEqual(primary["L5_FROZEN"]["auto"]["precision"], 42 / 73)
            self.assertAlmostEqual(primary["L5_FROZEN"]["auto"]["recall"], 42 / 58)
            self.assertEqual(
                (
                    primary["S_CORR_RSC_A"]["attribution"]["added_non_e0_auto"]["tp"],
                    primary["S_CORR_RSC_A"]["attribution"]["added_non_e0_auto"]["fp"],
                ),
                (0, 6),
            )
            self.assertEqual(
                (
                    primary["X_CORR_RSC_UNGATED_AB"]["attribution"]["added_non_e0_auto"]["tp"],
                    primary["X_CORR_RSC_UNGATED_AB"]["attribution"]["added_non_e0_auto"]["fp"],
                ),
                (3, 7),
            )
            self.assertEqual(
                (
                    primary["X_RAW_RSC_UNGATED_AB"]["attribution"]["added_non_e0_auto"]["tp"],
                    primary["X_RAW_RSC_UNGATED_AB"]["attribution"]["added_non_e0_auto"]["fp"],
                ),
                (5, 13),
            )
            self.assertEqual(
                (
                    primary["S_RAW_M_RST_C"]["attribution"]["promoted_e0_manual_auto"]["tp"],
                    primary["S_RAW_M_RST_C"]["attribution"]["promoted_e0_manual_auto"]["fp"],
                ),
                (2, 5),
            )
            self.assertEqual(
                (sensitivity["L5_FROZEN"]["auto"]["tp"], sensitivity["L5_FROZEN"]["auto"]["fp"]),
                (44, 25),
            )

    def test_evaluator_rejects_manifest_hash_tamper(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "abc"
            self._replay(output)
            manifest_path = output / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["artifacts"][0]["sha256"] = "0" * 64
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaises(EvaluationError):
                evaluate_directory(
                    output,
                    WORKSPACE / "RAG_re" / "gold" / "KMUH_clinical_pathogen_gold_20260811.json",
                    PROJECT / "config" / "evaluation_protocol_v1.json",
                    bootstrap_replicates=0,
                )


if __name__ == "__main__":
    unittest.main()


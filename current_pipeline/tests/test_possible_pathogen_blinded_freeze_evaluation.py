import csv
import json
import tempfile
import unittest
from pathlib import Path

from tests.test_adjudicate_rag_v2_packets import packet as base_packet
from tests.test_adjudicate_rag_v2_packets import result as base_result
from tools.evaluate_possible_pathogen_blinded_adjudications import run as evaluate
from tools.freeze_possible_pathogen_blinded_adjudications import run as freeze


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def write_jsonl(path: Path, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")


def blind_fixture(root: Path):
    packet = base_packet()
    packet["case_id"] = "BLIND-001"
    result = base_result()
    result["case_id"] = "BLIND-001"
    packets_path = root / "packets.jsonl"
    adjudication_dir = root / "adjudication"
    write_jsonl(packets_path, [packet])
    write_jsonl(adjudication_dir / "adjudications.jsonl", [{
        "case_id": "BLIND-001",
        "adjudication": result,
        "answer_source_read": False,
    }])
    write_jsonl(adjudication_dir / "failures.jsonl", [])
    write_json(adjudication_dir / "run_manifest.json", {
        "completed_count": 1,
        "failure_count": 0,
        "answer_source_read": False,
        "backend": "test",
        "model": "test",
        "reasoning_effort": "test",
    })
    return packets_path, adjudication_dir


class TestPossiblePathogenBlindFreezeEvaluation(unittest.TestCase):
    def test_complete_blind_run_freezes_and_evaluates(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            packets, adjudication = blind_fixture(root)
            frozen = root / "frozen"
            manifest = freeze(packets, adjudication, frozen)
            self.assertTrue(manifest["ready_for_posthoc_unblinding"])

            private_map = root / "private_map.json"
            write_json(private_map, {"rows": [{
                "blind_case_id": "BLIND-001",
                "patient_id": "1",
                "organism_name": "Bacteroides fragilis",
                "possible_reporting_role": "analytical_or_direct_possible_pathogen",
            }]})
            shadow = root / "shadow"
            shadow.mkdir()
            with (shadow / "possible_pathogen_decisions.csv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "patient_id", "organism_name", "clinical_decision",
                    "possible_reporting_role",
                ])
                writer.writeheader()
                writer.writerow({
                    "patient_id": 1,
                    "organism_name": "Bacteroides fragilis",
                    "clinical_decision": "review_high_priority",
                    "possible_reporting_role": "analytical_or_direct_possible_pathogen",
                })
            answers = root / "answers.csv"
            answers.write_text(
                "patient_id,answers\n1,Bacteroides fragilis\n", encoding="utf-8"
            )
            evaluation = evaluate(
                frozen, private_map, shadow, answers, root / "evaluation"
            )
            self.assertEqual(1, evaluation["adjudication_count"])
            self.assertEqual(1, evaluation["ober_visible_possible_count"])
            self.assertEqual(
                1, evaluation["metrics"]["strict_plus_ober_visible_possible"]["matched"]
            )

    def test_freeze_rejects_failure_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            packets, adjudication = blind_fixture(root)
            write_jsonl(adjudication / "failures.jsonl", [{"case_id": "BLIND-001"}])
            with self.assertRaises(ValueError):
                freeze(packets, adjudication, root / "frozen")


if __name__ == "__main__":
    unittest.main()

import csv
import json
import tempfile
import unittest
from pathlib import Path

from tools.validate_taxonomy_review_completion import validate


class ValidateTaxonomyReviewCompletionTest(unittest.TestCase):
    def _fixture(self, root: Path, *, write_decision: bool) -> dict[str, Path]:
        baseline = root / "baseline.csv"
        with baseline.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["organism_name", "taxonomy_mapping_status"],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "organism_name": "Example organism",
                    "taxonomy_mapping_status": "unmapped",
                }
            )

        queue = root / "queue.jsonl"
        queue.write_text(
            json.dumps(
                {
                    "organism_name": "Example organism",
                    "canonical_key": "exampleorganism",
                }
            )
            + "\n",
            encoding="utf-8",
        )

        decisions = root / "decisions"
        (decisions / "decisions").mkdir(parents=True)
        if write_decision:
            (decisions / "decisions" / "exampleorganism.json").write_text(
                json.dumps(
                    {
                        "canonical_key": "exampleorganism",
                        "status": "provisional_unmapped",
                    }
                ),
                encoding="utf-8",
            )
        (decisions / "run_summary.json").write_text(
            json.dumps(
                {
                    "mode": "executed",
                    "completion_status": "complete" if write_decision else "partial_failure",
                    "failures": 0 if write_decision else 1,
                    "pending_keys": [] if write_decision else ["exampleorganism"],
                    "remaining_eligible_not_selected": 0,
                }
            ),
            encoding="utf-8",
        )

        overlay = root / "overlay.json"
        overlay.write_text(
            json.dumps(
                {
                    "schema_version": "organism_taxonomy_auto_v1",
                    "exact_profiles": {},
                }
            ),
            encoding="utf-8",
        )
        audit = root / "audit.json"
        audit.write_text(
            json.dumps(
                {
                    "safe_to_consider_apply": True,
                    "staged_overlay": str(overlay.resolve()),
                }
            ),
            encoding="utf-8",
        )
        return {
            "baseline": baseline,
            "queue": queue,
            "decisions": decisions,
            "overlay": overlay,
            "audit": audit,
            "output": root / "completion.json",
        }

    def test_passes_when_every_unmapped_name_has_a_terminal_decision(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._fixture(Path(tmp), write_decision=True)
            result = validate(
                paths["baseline"],
                paths["queue"],
                paths["decisions"],
                paths["overlay"],
                paths["audit"],
                paths["output"],
            )
            self.assertTrue(result["passed"])
            self.assertEqual(result["held_count"], 1)
            self.assertEqual(result["unreviewed_keys"], [])

    def test_fails_closed_when_a_queue_item_has_no_decision(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._fixture(Path(tmp), write_decision=False)
            with self.assertRaisesRegex(RuntimeError, "validation failed"):
                validate(
                    paths["baseline"],
                    paths["queue"],
                    paths["decisions"],
                    paths["overlay"],
                    paths["audit"],
                    paths["output"],
                )
            manifest = json.loads(paths["output"].read_text(encoding="utf-8"))
            self.assertFalse(manifest["passed"])
            self.assertEqual(manifest["unreviewed_keys"], ["exampleorganism"])


if __name__ == "__main__":
    unittest.main()

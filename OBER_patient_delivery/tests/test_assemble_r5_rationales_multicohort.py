from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "assemble_r5_rationales.py"


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")


def rationale(patient_id: str, organism: str) -> dict:
    return {
        "patient_id": patient_id,
        "generated_rationale": {
            "selected_pathogen_explanations": [{"organism": organism, "rank": 1}]
        },
    }


class AssembleR5RationalesMultiCohortTests(unittest.TestCase):
    def test_explicit_cohort_batches_allow_overlapping_patient_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            base_rows = []
            for cohort, organism in (("alpha", "Organism A"), ("beta", "Organism B")):
                old_file = root / "old" / cohort / "p1.json"
                old = rationale("1", f"Old {cohort}")
                write_json(old_file, old)
                base_rows.append(
                    {"cohort": cohort, "output_file": str(old_file), "result": old}
                )

                run_root = root / "runs" / cohort
                result_file = run_root / "patient_rationales" / "NGS_patient_1_llm_rationale.json"
                result = rationale("1", organism)
                write_json(result_file, result)
                write_json(
                    run_root / "llm_rationale_manifest.json",
                    {
                        "failures": [],
                        "patients": [{"patient_id": "1", "output_file": str(result_file)}],
                    },
                )
                write_json(
                    root / "decisions" / cohort / "P1_final_decision.json",
                    {"patient_id": "1", "selected": [{"organism": organism, "rank": 1}]},
                )

            write_json(
                root / "accepted.json",
                {"status": "structurally_validated", "patients": base_rows},
            )
            write_json(
                root / "config.json",
                {
                    "accepted_rationales": "accepted.json",
                    "cohorts": {"alpha": {}, "beta": {}},
                },
            )
            output = root / "assembled.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(SCRIPT),
                    "--config",
                    str(root / "config.json"),
                    "--final-decisions",
                    str(root / "decisions"),
                    "--regenerated-cohort",
                    f"alpha={root / 'runs' / 'alpha'}",
                    "--regenerated-cohort",
                    f"beta={root / 'runs' / 'beta'}",
                    "--output",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )

            self.assertEqual(0, completed.returncode, completed.stderr)
            assembled = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(2, len(assembled["patients"]))
            self.assertEqual(2, assembled["provenance"]["regenerated_patients"])
            self.assertEqual(
                {"alpha", "beta"},
                {row["cohort"] for row in assembled["patients"]},
            )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import json
import tempfile
import unittest
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path

from core import summarize_agent
from tools import deterministic_mngs_max_scorer, mngs_common


class FilmArrayDefaultTests(unittest.TestCase):
    def test_summary_collects_only_standard_filmarray_agent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            patient = Path(temporary) / "NGS_patient_1_json"
            agent_dir = patient / "agent_outputs"
            agent_dir.mkdir(parents=True)
            base = "NGS_patient_1"
            for suffix in (
                "cbc_other_lab_agent",
                "image_agent",
                "filmarray_gmtest_agent",
                "culture_agent",
            ):
                (agent_dir / f"{base}_{suffix}.json").write_text(
                    json.dumps({"source": suffix}),
                    encoding="utf-8",
                )

            sections = summarize_agent.collect_agent_outputs(
                patient,
                include_underlying=False,
            )

            paths = [path.name for _, path, _ in sections]
            self.assertIn(f"{base}_filmarray_gmtest_agent.json", paths)
            self.assertTrue(all("no_filmarray" not in path for path in paths))
            self.assertEqual(
                summarize_agent.derive_output_path(patient, include_underlying=False).name,
                f"{base}_final_summary_with_filmarray.json",
            )

    def test_auto_summary_lookup_never_falls_back_to_no_filmarray(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            patient = Path(temporary) / "NGS_patient_1_json"
            summary_dir = patient / "summary_outputs"
            summary_dir.mkdir(parents=True)
            no_filmarray = summary_dir / "NGS_patient_1_final_summary_no_filmarray.json"
            no_filmarray.write_text("{}", encoding="utf-8")

            self.assertIsNone(mngs_common.find_final_summary_file(patient, summary_mode="auto"))

            full = summary_dir / "NGS_patient_1_final_summary_with_filmarray.json"
            full.write_text("{}", encoding="utf-8")
            self.assertEqual(
                mngs_common.find_final_summary_file(patient, summary_mode="auto"),
                full,
            )

    def test_v20_cli_rejects_no_filmarray_mode(self) -> None:
        with redirect_stderr(StringIO()):
            with self.assertRaises(SystemExit):
                deterministic_mngs_max_scorer.parse_args(
                    ["patient-root", "--summary-mode", "no_filmarray"]
                )


if __name__ == "__main__":
    unittest.main()

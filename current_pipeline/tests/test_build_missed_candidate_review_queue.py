from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools import build_missed_candidate_review_queue as queue


class ReviewQueueSummaryPathTests(unittest.TestCase):
    def test_explicit_shadow_summary_never_falls_back_to_baseline(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            patient_dir = Path(temporary) / "NGS_patient_14_json"
            (patient_dir / "summary_outputs").mkdir(parents=True)
            baseline = queue.output_path_for(
                patient_dir, "final_summary_with_filmarray_deterministic"
            )
            baseline.write_text("{}", encoding="utf-8")

            self.assertIsNone(
                queue.final_summary_path_for(
                    patient_dir, "final_summary_evidence_v2_shadow_20260917"
                )
            )

            shadow = queue.output_path_for(
                patient_dir, "final_summary_evidence_v2_shadow_20260917"
            )
            shadow.write_text("{}", encoding="utf-8")
            self.assertEqual(
                shadow,
                queue.final_summary_path_for(
                    patient_dir, "final_summary_evidence_v2_shadow_20260917"
                ),
            )


if __name__ == "__main__":
    unittest.main()

import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from tools.evaluate_external_kh_frozen_release import run


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class EvaluateExternalKhFrozenReleaseTest(unittest.TestCase):
    def test_explicit_negative_predictions_count_as_false_positives(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            answers = root / "answers.csv"
            with answers.open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["patient_id", "answers", "label_status"],
                )
                writer.writeheader()
                writer.writerow({
                    "patient_id": 1,
                    "answers": "Staphylococcus aureus",
                    "label_status": "positive",
                })
                writer.writerow({
                    "patient_id": 2,
                    "answers": "-",
                    "label_status": "explicit_no_pathogen",
                })
            manifest = root / "manifest.json"
            manifest.write_text(json.dumps({
                "scope": "post-freeze evaluation only",
                "freeze_passed_before_unblinding": True,
                "answers_loaded_into_decision_pipeline": False,
                "outputs": {"full_answers_sha256": sha256_file(answers)},
                "source_contracts": {"decision_freeze_sha256": "freeze-hash"},
            }), encoding="utf-8")
            shadow = root / "shadow"
            shadow.mkdir()
            with (shadow / "possible_pathogen_decisions.csv").open(
                "w", encoding="utf-8-sig", newline=""
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "patient_id", "organism_name", "possible_reporting_role",
                    "selected_for_complete_report",
                ])
                writer.writeheader()
                writer.writerow({
                    "patient_id": 1,
                    "organism_name": "Staphylococcus aureus",
                    "possible_reporting_role": "strict_picked_report",
                    "selected_for_complete_report": "True",
                })
                writer.writerow({
                    "patient_id": 2,
                    "organism_name": "Candida albicans",
                    "possible_reporting_role": "strict_picked_report",
                    "selected_for_complete_report": "True",
                })

            result = run(shadow, answers, manifest, root / "output")
            primary = result["all_labeled_patients_including_explicit_no_pathogen"][
                "strict_picked"
            ]
            legacy = result["legacy_positive_patients_only"]["strict_picked"]
            all_possible = result[
                "all_labeled_patients_including_explicit_no_pathogen"
            ]["picked_plus_all_possible"]
            self.assertEqual(primary["true_positive"], 1)
            self.assertEqual(primary["false_positive"], 1)
            self.assertEqual(primary["false_negative"], 0)
            self.assertEqual(primary["precision"], 0.5)
            self.assertEqual(primary["recall"], 1.0)
            self.assertEqual(
                primary["explicit_no_pathogen_patients_with_predictions"], [2]
            )
            self.assertEqual(legacy["precision"], 1.0)
            self.assertEqual(all_possible, primary)


if __name__ == "__main__":
    unittest.main()

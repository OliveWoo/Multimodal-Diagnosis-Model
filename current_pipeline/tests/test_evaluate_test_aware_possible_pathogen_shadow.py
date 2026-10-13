import csv
import json
import tempfile
import unittest
from pathlib import Path

from tools.evaluate_test_aware_possible_pathogen_shadow import run


TEST_TMP = Path(__file__).parents[1] / ".test_tmp"
TEST_TMP.mkdir(exist_ok=True)


FIELDS = [
    "patient_id",
    "organism_name",
    "taxonomy_family",
    "clinical_decision",
    "possible_reporting_role",
    "selected_for_complete_report",
    "possible_route",
]


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class StatusAwareEvaluationTests(unittest.TestCase):
    def test_explicit_negative_and_strata_are_counted(self) -> None:
        with tempfile.TemporaryDirectory(dir=TEST_TMP) as raw:
            root = Path(raw)
            shadow = root / "shadow"
            answers = root / "answers.csv"
            output = root / "evaluation"
            write_csv(
                shadow / "possible_pathogen_decisions.csv",
                FIELDS,
                [
                    {
                        "patient_id": "1",
                        "organism_name": "Klebsiella pneumoniae",
                        "taxonomy_family": "typical_respiratory_pathogen",
                        "clinical_decision": "picked_shadow",
                        "possible_reporting_role": "strict_picked_report",
                        "selected_for_complete_report": "True",
                        "possible_route": "strict",
                    },
                    {
                        "patient_id": "2",
                        "organism_name": "Candida albicans",
                        "taxonomy_family": "candida_or_yeast",
                        "clinical_decision": "review_context_needed",
                        "possible_reporting_role": "fallback_best_available_possible_pathogen",
                        "selected_for_complete_report": "True",
                        "possible_route": "fallback",
                    },
                ],
            )
            write_csv(
                answers,
                ["patient_id", "answers", "cohort", "label_status"],
                [
                    {
                        "patient_id": "1",
                        "answers": "Klebsiella pneumoniae",
                        "cohort": "KH33",
                        "label_status": "positive",
                    },
                    {
                        "patient_id": "2",
                        "answers": "",
                        "cohort": "KH33",
                        "label_status": "explicit_no_pathogen",
                    },
                    {
                        "patient_id": "3",
                        "answers": "Escherichia coli",
                        "cohort": "two_hospital_A11",
                        "label_status": "positive",
                    },
                ],
            )

            result = run(shadow, answers, output)

            complete = result["metrics"]["complete_report"]
            self.assertEqual(complete["true_positive"], 1)
            self.assertEqual(complete["false_positive"], 1)
            self.assertEqual(complete["false_negative"], 1)
            self.assertEqual(complete["precision"], 0.5)
            self.assertEqual(complete["recall"], 0.5)
            self.assertEqual(
                complete["explicit_no_pathogen_patients_with_predictions"], [2]
            )
            self.assertIn("without_B11", result["stratified_metrics"])
            with (output / "missed_answers.csv").open(
                encoding="utf-8-sig"
            ) as handle:
                missed = list(csv.DictReader(handle))
            self.assertTrue(
                any(
                    row["endpoint"] == "complete_report"
                    and row["patient_id"] == "3"
                    and row["answer"] == "Escherichia coli"
                    for row in missed
                )
            )
            with (output / "explicit_no_pathogen_predictions.csv").open(
                encoding="utf-8-sig"
            ) as handle:
                audit = list(csv.DictReader(handle))
            self.assertEqual(audit[0]["patient_id"], "2")
            self.assertEqual(
                json.loads((output / "metrics.json").read_text(encoding="utf-8"))[
                    "schema_version"
                ],
                "test_aware_posthoc_evaluation.v2",
            )


if __name__ == "__main__":
    unittest.main()

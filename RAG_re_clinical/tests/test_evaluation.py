from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from rag_re_clinical.evaluation import (  # noqa: E402
    BASELINE_ARM,
    FIXED_ARM_NAMES,
    EvaluationError,
    evaluate_directory,
    render_markdown,
    write_markdown,
)


def make_candidate(
    name: str,
    *,
    baseline: bool = False,
    overrides: dict[str, bool] | None = None,
) -> dict:
    decisions: dict[str, bool] = {
        arm: baseline for arm in FIXED_ARM_NAMES
    }
    if not baseline:
        decisions = {arm: False for arm in FIXED_ARM_NAMES}
    decisions[BASELINE_ARM] = baseline
    decisions.update(overrides or {})
    return {
        "organism_name": name,
        "baseline_selected": baseline,
        "review_tier": None if baseline else "review_context_needed",
        "disposition": "baseline" if baseline else "eligible_non_e0",
        "modules": {},
        "arms": decisions,
    }


def make_artifact(
    patient_id: str,
    candidates: list[dict],
    *,
    manual: dict[str, list[str]] | None = None,
) -> dict:
    manual = manual or {}
    arms = {}
    for arm in FIXED_ARM_NAMES:
        predicted = [
            candidate["organism_name"]
            for candidate in candidates
            if candidate["arms"][arm] is True
        ]
        manual_names = list(manual.get(arm, []))
        arms[arm] = {
            "predicted_pathogens": predicted,
            "manual_review_pathogens": manual_names,
            "predicted_count": len(predicted),
            "manual_review_count": len(manual_names),
        }
    return {
        "schema_version": "rag_re_clinical.output.v1",
        "patient_id": str(patient_id),
        "run_status": "complete",
        "candidates": candidates,
        "arms": arms,
    }


def make_gold(cases: list[tuple[str, list[str]]]) -> dict:
    return {
        "schema_version": "rag_re.gold.v1",
        "label_semantics": "infection_source_positive_list",
        "complete_candidate_negative_labels": True,
        "aliases": {
            "Human alphaherpesvirus 1": "HSV-1",
            "Herpes simplex virus type 1": "HSV-1",
        },
        "cases": [
            {
                "patient_id": str(patient_id),
                "pathogens": pathogens,
                "uncertain_pathogens": [],
            }
            for patient_id, pathogens in cases
        ],
    }


class EvaluationTests(unittest.TestCase):
    def write_json(self, path: Path, payload: object) -> None:
        path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")

    def build_main_fixture(self, root: Path) -> tuple[Path, Path]:
        outputs = root / "outputs"
        outputs.mkdir()
        gold_path = root / "gold.json"
        self.write_json(
            gold_path,
            make_gold(
                [
                    (
                        "1",
                        ["Human alphaherpesvirus 1", "Klebsiella pneumoniae"],
                    ),
                    ("2", ["Klebsiella pneumoniae"]),
                    ("3", []),
                    ("4", ["Acinetobacter baumannii"]),
                    ("5", []),
                    ("6", []),
                ]
            ),
        )

        patient_1 = make_artifact(
            "1",
            [
                make_candidate("HSV-1", baseline=True),
                make_candidate(
                    "Klebsiella pneumoniae",
                    overrides={"CL2_CONVERGENT": True},
                ),
                make_candidate(
                    "Candida albicans"
                ),
                make_candidate("Acinetobacter nosocomialis"),
            ],
            manual={arm: ["Candida albicans"] for arm in FIXED_ARM_NAMES},
        )
        patient_2 = make_artifact(
            "2",
            [
                make_candidate(
                    "Klebsiella variicola",
                    overrides={"CL2_CONVERGENT": True},
                )
            ],
        )
        patient_3 = make_artifact(
            "3",
            [
                make_candidate(
                    "Candida albicans", overrides={"CL2_CONVERGENT": True}
                )
            ],
        )
        patient_6 = make_artifact("6", [])
        for patient_id, payload in (
            ("1", patient_1),
            ("2", patient_2),
            ("3", patient_3),
            ("6", patient_6),
        ):
            self.write_json(outputs / f"patient_{patient_id}.json", payload)
        self.write_json(
            outputs / "manifest.json",
            {
                "schema_version": "rag_re_clinical.batch_manifest.v1",
                "artifact_count": 4,
            },
        )
        return outputs, gold_path

    def test_exploratory_rescue_may_leave_a_primary_manual_queue(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs = root / "outputs"
            outputs.mkdir()
            gold = root / "gold.json"
            candidate = make_candidate(
                "Review species", overrides={"CL3_BALANCED": True}
            )
            artifact = make_artifact(
                "1",
                [candidate],
                manual={
                    arm: ["Review species"]
                    for arm in FIXED_ARM_NAMES
                    if arm != "CL3_BALANCED"
                },
            )
            self.write_json(outputs / "patient_1.json", artifact)
            self.write_json(gold, make_gold([("1", ["Review species"])]))
            report = evaluate_directory(outputs, gold)

        primary = report["scopes"]["all_labeled_exact_frozen_synonyms"]["arms"]
        self.assertEqual(primary["CL2_CONVERGENT"]["manual_review_count"], 1)
        self.assertEqual(primary["CL3_BALANCED"]["manual_review_count"], 0)
        self.assertEqual(primary["CL3_BALANCED"]["tp"], 1)

    def test_primary_sensitivity_rescue_missingness_and_review_counts(self):
        with tempfile.TemporaryDirectory() as temp:
            outputs, gold = self.build_main_fixture(Path(temp))
            report = evaluate_directory(outputs, gold)

        self.assertEqual(report["arm_names"], list(FIXED_ARM_NAMES))
        self.assertEqual(report["primary_arm"], "CL2_CONVERGENT")
        self.assertEqual(report["artifacts"]["batch_manifest_count"], 1)
        self.assertEqual(
            report["artifacts"]["missing_labeled_patient_ids"], ["4", "5"]
        )

        primary = report["scopes"]["all_labeled_exact_frozen_synonyms"]
        cl2 = primary["arms"]["CL2_CONVERGENT"]
        self.assertEqual((cl2["tp"], cl2["fp"], cl2["fn"]), (2, 2, 2))
        self.assertEqual(cl2["precision"], 0.5)
        self.assertEqual(cl2["recall"], 0.5)
        self.assertEqual(cl2["f0_5"], 0.5)
        self.assertEqual(cl2["f1"], 0.5)
        self.assertEqual(
            cl2["rescue_from_baseline"],
            {
                "added_tp": 1,
                "added_fp": 2,
                "added_predictions": 3,
                "ppv": 1 / 3,
            },
        )
        self.assertEqual(cl2["manual_review_count"], 1)
        self.assertEqual(cl2["manual_review_case_count"], 1)
        self.assertEqual(cl2["unknown_count"], 0)
        self.assertEqual(cl2["unknown_case_count"], 0)
        self.assertEqual(cl2["no_path_cases"], 3)
        self.assertEqual(cl2["no_path_correct_negative_cases"], 1)
        self.assertEqual(cl2["no_path_predicted_positive_cases"], 1)
        self.assertEqual(cl2["no_path_missing_artifact_cases"], 1)
        self.assertAlmostEqual(cl2["no_path_specificity"], 1 / 3)
        self.assertEqual(cl2["case_counts"]["artifact_missing"], 2)

        sensitivity = report["scopes"]["answer_positive_only_genus_relaxed"]
        cl2_sensitivity = sensitivity["arms"]["CL2_CONVERGENT"]
        self.assertEqual(
            (cl2_sensitivity["tp"], cl2_sensitivity["fp"], cl2_sensitivity["fn"]),
            (3, 0, 1),
        )
        self.assertEqual(cl2_sensitivity["precision"], 1.0)
        self.assertEqual(cl2_sensitivity["recall"], 0.75)
        self.assertAlmostEqual(cl2_sensitivity["f0_5"], 0.9375)
        self.assertAlmostEqual(cl2_sensitivity["f1"], 6 / 7)
        self.assertEqual(
            cl2_sensitivity["rescue_from_baseline"],
            {
                "added_tp": 2,
                "added_fp": 0,
                "added_predictions": 2,
                "ppv": 1.0,
            },
        )

    def test_genus_relaxed_matching_is_one_to_one(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs = root / "outputs"
            outputs.mkdir()
            gold = root / "gold.json"
            self.write_json(
                gold,
                make_gold(
                    [
                        (
                            "1",
                            ["Klebsiella pneumoniae", "Klebsiella oxytoca"],
                        )
                    ]
                ),
            )
            artifact = make_artifact(
                "1",
                [
                    make_candidate(
                        "Klebsiella variicola",
                        overrides={"CL2_CONVERGENT": True},
                    )
                ],
            )
            self.write_json(outputs / "patient_1.json", artifact)
            report = evaluate_directory(outputs, gold)

        primary = report["scopes"]["all_labeled_exact_frozen_synonyms"]["arms"]
        self.assertEqual(
            (primary["CL2_CONVERGENT"]["tp"], primary["CL2_CONVERGENT"]["fn"]),
            (0, 2),
        )
        sensitivity = report["scopes"]["answer_positive_only_genus_relaxed"][
            "arms"
        ]["CL2_CONVERGENT"]
        self.assertEqual((sensitivity["tp"], sensitivity["fp"], sensitivity["fn"]), (1, 0, 1))

    def test_frozen_cmv_alias_reconciles_root_and_candidate_names(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs = root / "outputs"
            outputs.mkdir()
            gold = root / "gold.json"
            self.write_json(gold, make_gold([("1", ["Human cytomegalovirus"])]))
            candidate = make_candidate("CMV", baseline=True)
            candidate["canonical_organism_name"] = "Human cytomegalovirus"
            self.write_json(
                outputs / "patient_1.json", make_artifact("1", [candidate])
            )
            report = evaluate_directory(outputs, gold)

        primary = report["scopes"]["all_labeled_exact_frozen_synonyms"]["arms"]
        for arm in FIXED_ARM_NAMES:
            self.assertEqual(
                (primary[arm]["tp"], primary[arm]["fp"], primary[arm]["fn"]),
                (1, 0, 0),
            )

    def test_empty_output_directory_keeps_end_to_end_false_negatives(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs = root / "outputs"
            outputs.mkdir()
            gold = root / "gold.json"
            self.write_json(
                gold,
                make_gold([("1", ["Klebsiella pneumoniae"]), ("2", [])]),
            )
            report = evaluate_directory(outputs, gold)

        metrics = report["scopes"]["all_labeled_exact_frozen_synonyms"]["arms"]
        for arm in FIXED_ARM_NAMES:
            self.assertEqual((metrics[arm]["tp"], metrics[arm]["fp"], metrics[arm]["fn"]), (0, 0, 1))
            self.assertEqual(metrics[arm]["no_path_specificity"], 0.0)
            self.assertEqual(metrics[arm]["no_path_missing_artifact_cases"], 1)

    def test_custom_arm_subset_automatically_includes_baseline(self):
        with tempfile.TemporaryDirectory() as temp:
            outputs, gold = self.build_main_fixture(Path(temp))
            report = evaluate_directory(outputs, gold, arm_names=["CL2_CONVERGENT"])
        self.assertEqual(report["arm_names"], [BASELINE_ARM, "CL2_CONVERGENT"])

    def test_rejects_prediction_outside_candidates(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs, gold = self.build_main_fixture(root)
            path = outputs / "patient_1.json"
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["arms"]["CL2_CONVERGENT"]["predicted_pathogens"].append(
                "Outside organism"
            )
            payload["arms"]["CL2_CONVERGENT"]["predicted_count"] += 1
            self.write_json(path, payload)
            with self.assertRaisesRegex(EvaluationError, "outside candidates"):
                evaluate_directory(outputs, gold)

    def test_rejects_e0_subset_violation(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs, gold = self.build_main_fixture(root)
            path = outputs / "patient_1.json"
            payload = json.loads(path.read_text(encoding="utf-8"))
            hsv = payload["candidates"][0]
            hsv["arms"]["CL2_CONVERGENT"] = False
            predicted = payload["arms"]["CL2_CONVERGENT"]["predicted_pathogens"]
            predicted.remove("HSV-1")
            payload["arms"]["CL2_CONVERGENT"]["predicted_count"] -= 1
            self.write_json(path, payload)
            with self.assertRaisesRegex(EvaluationError, "E0 immutability"):
                evaluate_directory(outputs, gold)

    def test_rejects_missing_arm_and_root_candidate_disagreement(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs, gold = self.build_main_fixture(root)
            path = outputs / "patient_1.json"
            original = json.loads(path.read_text(encoding="utf-8"))

            missing = copy.deepcopy(original)
            del missing["arms"]["CL4_TIERED_EXPLORATORY"]
            self.write_json(path, missing)
            with self.assertRaisesRegex(EvaluationError, "missing root arm"):
                evaluate_directory(outputs, gold)

            disagreeing = copy.deepcopy(original)
            disagreeing["candidates"][1]["arms"]["CL1_DIRECT"] = True
            self.write_json(path, disagreeing)
            with self.assertRaisesRegex(EvaluationError, "root predictions disagree"):
                evaluate_directory(outputs, gold)

    def test_unknown_schema_is_rejected_but_batch_manifest_is_skipped(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs, gold = self.build_main_fixture(root)
            self.write_json(
                outputs / "unknown.json",
                {"schema_version": "not.a.manifest", "patient_id": "999"},
            )
            with self.assertRaisesRegex(EvaluationError, "unsupported schema_version"):
                evaluate_directory(outputs, gold)

    def test_duplicate_patient_artifacts_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs, gold = self.build_main_fixture(root)
            payload = json.loads(
                (outputs / "patient_1.json").read_text(encoding="utf-8")
            )
            self.write_json(outputs / "patient_1_copy.json", payload)
            with self.assertRaisesRegex(EvaluationError, "Duplicate clinical artifact"):
                evaluate_directory(outputs, gold)

    def test_markdown_render_and_write(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            outputs, gold = self.build_main_fixture(root)
            report = evaluate_directory(outputs, gold)
            rendered = render_markdown(report)
            self.assertIn("CL2_CONVERGENT", rendered)
            self.assertIn("Missing patient IDs: 4, 5", rendered)
            markdown_path = root / "reports" / "evaluation.md"
            returned = write_markdown(report, markdown_path)
            self.assertEqual(returned, markdown_path)
            self.assertEqual(markdown_path.read_text(encoding="utf-8"), rendered)


if __name__ == "__main__":
    unittest.main()

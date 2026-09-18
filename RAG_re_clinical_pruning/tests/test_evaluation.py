from __future__ import annotations

import copy
import hashlib
import json
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from rag_re_clinical_pruning.evaluation import (  # noqa: E402
    BASELINE_ARM,
    CUMULATIVE_ARMS,
    EvaluationError,
    FIXED_ARM_NAMES,
    PRIMARY_ARM,
    evaluate_directory,
    render_markdown,
    write_json,
    write_markdown,
)
from evaluate_pruning import main as evaluate_cli_main  # noqa: E402


CANDIDATES = (
    ("Klebsiella pneumoniae", True),
    ("Candida albicans", True),
    ("Pseudomonas aeruginosa", False),
    ("Human cytomegalovirus", False),
)


def _default_state(baseline: bool) -> str:
    return "KEEP" if baseline else "NO_RESCUE"


def _make_artifact(patient_id: str = "p1", *, unsafe: bool = False) -> dict:
    states = {
        arm: {name: _default_state(baseline) for name, baseline in CANDIDATES}
        for arm in FIXED_ARM_NAMES
    }
    if unsafe:
        for arm in CUMULATIVE_ARMS[1:]:
            states[arm]["Klebsiella pneumoniae"] = "PRUNE"
            states[arm]["Candida albicans"] = "PRUNE"
    else:
        states["L1_WEAK_RAW"]["Candida albicans"] = "PRUNE"
        for arm in CUMULATIVE_ARMS[2:]:
            states[arm]["Candida albicans"] = "KEEP"
        for arm in CUMULATIVE_ARMS[3:]:
            states[arm]["Human cytomegalovirus"] = "MANUAL"
        for arm in CUMULATIVE_ARMS[5:]:
            states[arm]["Pseudomonas aeruginosa"] = "RESCUE"

    candidates = []
    for name, baseline in CANDIDATES:
        candidates.append(
            {
                "organism_name": name,
                "canonical_organism_name": name,
                "baseline_selected": baseline,
                "arm_states": {
                    arm: {"state": states[arm][name], "rule_ledger": []}
                    for arm in FIXED_ARM_NAMES
                },
            }
        )

    arms = {}
    for arm in FIXED_ARM_NAMES:
        state_by_name = states[arm]
        auto = [name for name, _ in CANDIDATES if state_by_name[name] in {"KEEP", "RESCUE"}]
        manual = [name for name, _ in CANDIDATES if state_by_name[name] == "MANUAL"]
        review = auto + manual
        pruned = [name for name, _ in CANDIDATES if state_by_name[name] == "PRUNE"]
        rescued = [name for name, _ in CANDIDATES if state_by_name[name] == "RESCUE"]
        no_rescue = [name for name, _ in CANDIDATES if state_by_name[name] == "NO_RESCUE"]
        state_counts = dict(Counter(state_by_name.values()))
        arms[arm] = {
            "auto_positive_pathogens": auto,
            "auto_positive_count": len(auto),
            "review_inclusive_pathogens": review,
            "review_inclusive_count": len(review),
            "manual_review_pathogens": manual,
            "manual_review_count": len(manual),
            "pruned_pathogens": pruned,
            "pruned_count": len(pruned),
            "rescued_pathogens": rescued,
            "rescued_count": len(rescued),
            "no_rescue_pathogens": no_rescue,
            "no_rescue_count": len(no_rescue),
            "state_counts": state_counts,
            "rule_sequence": [],
            "family": "standalone" if arm.startswith("S") else "cumulative",
        }
    return {
        "schema_version": "rag_re_clinical_pruning.output.v1",
        "patient_id": patient_id,
        "run_status": "complete",
        "candidates": candidates,
        "arms": arms,
    }


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _gold(*cases: dict, aliases: dict | None = None) -> dict:
    return {
        "schema_version": "clinical_gold.v1",
        "aliases": aliases or {},
        "cases": list(cases),
    }


def _run(tmp_path: Path, artifact: dict, gold: dict, **kwargs):
    outputs = tmp_path / "outputs"
    _write_json(outputs / "p1.json", artifact)
    gold_path = tmp_path / "gold.json"
    _write_json(gold_path, gold)
    return evaluate_directory(outputs, gold_path, **kwargs)


def test_metrics_attribution_deltas_and_safety(tmp_path: Path) -> None:
    artifact = _make_artifact()
    report = _run(
        tmp_path,
        artifact,
        _gold(
            {
                "patient_id": "p1",
                "pathogens": [
                    "Klebsiella pneumoniae",
                    "Pseudomonas aeruginosa",
                    "CMV",
                ],
            }
        ),
    )
    primary = report["scopes"]["all_labeled_exact_frozen_synonyms"]["arms"]
    row = primary[PRIMARY_ARM]
    assert row["auto"]["tp"] == 2
    assert row["auto"]["fp"] == 1
    assert row["auto"]["fn"] == 1
    assert row["review_inclusive"]["tp"] == 3
    assert row["review_inclusive"]["fp"] == 1
    assert row["attribution"]["manual"] == {
        "count": 1,
        "tp": 1,
        "fp": 0,
        "ppv": 1.0,
    }
    assert row["attribution"]["added_non_e0_rescue"]["tp"] == 1
    assert row["attribution"]["added_non_e0_rescue"]["fp"] == 0
    assert row["delta_vs_e0"]["auto"]["tp"] == 1
    assert row["delta_vs_e0"]["auto"]["fp"] == 0
    assert row["delta_vs_previous_cumulative"]["reference_arm"] == "L4_PLUS_GUARD_REVIEW"
    assert row["delta_vs_previous_cumulative"]["auto"]["tp"] == 1
    assert row["attribution"]["stage_added_non_e0_rescue"]["tp"] == 1
    assert row["recall_safety"]["status"] == "PASS"
    assert report["recall_safety_summary"][PRIMARY_ARM]["status"] == "PASS"

    l1 = primary["L1_WEAK_RAW"]
    assert l1["attribution"]["explicit_pruned"]["tp"] == 0
    assert l1["attribution"]["explicit_pruned"]["fp"] == 1
    l2 = primary["L2_DIRECT_RESTORE"]
    assert l2["attribution"]["restored_from_prune"]["tp"] == 0
    assert l2["attribution"]["restored_from_prune"]["fp"] == 1
    assert l2["attribution"]["stage_restored_from_prune"]["fp"] == 1


def test_primary_exact_and_positive_only_genus_relaxed_are_separate(tmp_path: Path) -> None:
    report = _run(
        tmp_path,
        _make_artifact(),
        _gold(
            {
                "patient_id": "p1",
                "pathogens": ["Klebsiella pneumoniae", "Pseudomonas putida"],
            }
        ),
    )
    exact = report["scopes"]["all_labeled_exact_frozen_synonyms"]["arms"][PRIMARY_ARM]
    relaxed = report["scopes"]["answer_positive_only_genus_relaxed"]["arms"][PRIMARY_ARM]
    assert exact["auto"]["tp"] == 1
    assert exact["auto"]["fp"] == 2
    assert relaxed["auto"]["tp"] == 2
    assert relaxed["auto"]["fp"] == 1


def test_patient_level_one_to_one_matching(tmp_path: Path) -> None:
    artifact = _make_artifact()
    report = _run(
        tmp_path,
        artifact,
        _gold(
            {
                "patient_id": "p1",
                "pathogens": ["Pseudomonas putida", "Pseudomonas fluorescens"],
            }
        ),
    )
    relaxed = report["scopes"]["answer_positive_only_genus_relaxed"]["arms"][PRIMARY_ARM]
    assert relaxed["auto"]["tp"] == 1
    assert relaxed["auto"]["fn"] == 1


def test_missing_artifacts_fail_and_remain_in_denominators(tmp_path: Path) -> None:
    outputs = tmp_path / "outputs"
    _write_json(outputs / "p1.json", _make_artifact())
    gold_path = tmp_path / "gold.json"
    _write_json(
        gold_path,
        _gold(
            {"patient_id": "p1", "pathogens": ["Klebsiella pneumoniae"]},
            {"patient_id": "missing-positive", "pathogens": ["Escherichia coli"]},
            {"patient_id": "missing-no-path", "pathogens": []},
        ),
    )
    report = evaluate_directory(outputs, gold_path)
    row = report["scopes"]["all_labeled_exact_frozen_synonyms"]["arms"][BASELINE_ARM]
    assert report["artifacts"]["missing_labeled_patient_ids"] == [
        "missing-no-path",
        "missing-positive",
    ]
    assert row["auto"]["tp"] == 1
    assert row["auto"]["fn"] == 1
    assert row["auto"]["no_path_cases"] == 1
    assert row["auto"]["no_path_correct_negative_cases"] == 0
    assert row["auto"]["no_path_missing_artifact_cases"] == 1
    assert row["auto"]["no_path_specificity"] == 0.0
    assert row["case_counts"]["included"] == 3
    assert row["case_counts"]["artifact_missing"] == 2


def test_recall_safety_unsafe_and_does_not_mutate_artifact(tmp_path: Path) -> None:
    artifact = _make_artifact(unsafe=True)
    before = copy.deepcopy(artifact)
    report = _run(
        tmp_path,
        artifact,
        _gold(
            {
                "patient_id": "p1",
                "pathogens": ["Klebsiella pneumoniae", "Candida albicans"],
            }
        ),
    )
    row = report["scopes"]["all_labeled_exact_frozen_synonyms"]["arms"]["L1_WEAK_RAW"]
    assert row["attribution"]["explicit_pruned"]["tp"] == 2
    assert row["recall_safety"]["status"] == "UNSAFE"
    assert row["recall_safety"]["criteria"]["auto_absolute_recall_drop"]["pass"] is False
    assert row["recall_safety"]["criteria"]["explicitly_pruned_tp"]["pass"] is False
    assert row["recall_safety"]["evaluation_only"] is True
    assert artifact == before


def test_uncertain_gold_is_neutral_not_positive_or_false_positive(tmp_path: Path) -> None:
    report = _run(
        tmp_path,
        _make_artifact(),
        _gold(
            {
                "patient_id": "p1",
                "pathogens": ["Klebsiella pneumoniae"],
                "uncertain_pathogens": ["Candida albicans"],
            }
        ),
    )
    baseline = report["scopes"]["all_labeled_exact_frozen_synonyms"]["arms"][BASELINE_ARM]
    assert baseline["auto"]["tp"] == 1
    assert baseline["auto"]["fp"] == 0


def test_manifest_skipped_and_unknown_json_rejected(tmp_path: Path) -> None:
    outputs = tmp_path / "outputs"
    _write_json(outputs / "p1.json", _make_artifact())
    _write_json(
        outputs / "manifest.json",
        {"schema_version": "rag_re_clinical_pruning.batch_manifest.v1"},
    )
    gold_path = tmp_path / "gold.json"
    _write_json(gold_path, _gold({"patient_id": "p1", "pathogens": []}))
    report = evaluate_directory(outputs, gold_path)
    assert report["artifacts"]["batch_manifest_count"] == 1
    _write_json(outputs / "junk.json", {"schema_version": "unexpected.v1"})
    with unittest.TestCase().assertRaisesRegex(EvaluationError, "unsupported schema_version"):
        evaluate_directory(outputs, gold_path)


def test_schema_and_provenance_violations_rejected(tmp_path: Path) -> None:
    scenarios = [
        (
            lambda artifact: (
                artifact["arms"][PRIMARY_ARM]["auto_positive_pathogens"].append(
                    "Outside organism"
                ),
                artifact["arms"][PRIMARY_ARM].update(
                    {
                        "auto_positive_count": artifact["arms"][PRIMARY_ARM][
                            "auto_positive_count"
                        ]
                        + 1
                    }
                ),
            ),
            "outside candidates",
        ),
        (
            lambda artifact: artifact["candidates"][0]["arm_states"][PRIMARY_ARM].update(
                {"state": "RESCUE"}
            ),
            "invalid for baseline_selected=True",
        ),
        (
            lambda artifact: artifact["candidates"][2]["arm_states"][PRIMARY_ARM].update(
                {"state": "PRUNE"}
            ),
            "invalid for baseline_selected=False",
        ),
        (
            lambda artifact: artifact["arms"].pop(PRIMARY_ARM),
            f"arm {PRIMARY_ARM} must be an object",
        ),
        (
            lambda artifact: artifact["arms"][PRIMARY_ARM].update(
                {"auto_positive_count": 999}
            ),
            "auto_positive_count does not match",
        ),
        (
            lambda artifact: artifact["arms"][PRIMARY_ARM].update(
                {"family": "wrong"}
            ),
            "family must be",
        ),
    ]
    for index, (mutator, message) in enumerate(scenarios):
        artifact = _make_artifact()
        mutator(artifact)
        with unittest.TestCase().assertRaisesRegex(EvaluationError, message):
            _run(
                tmp_path / str(index),
                artifact,
                _gold({"patient_id": "p1", "pathogens": []}),
            )


def test_root_lists_must_agree_with_candidate_states(tmp_path: Path) -> None:
    artifact = _make_artifact()
    arm = artifact["arms"][PRIMARY_ARM]
    arm["rescued_pathogens"] = []
    arm["rescued_count"] = 0
    with unittest.TestCase().assertRaisesRegex(EvaluationError, "rescued_pathogens disagrees"):
        _run(tmp_path, artifact, _gold({"patient_id": "p1", "pathogens": []}))


def test_duplicate_patient_artifacts_rejected(tmp_path: Path) -> None:
    outputs = tmp_path / "outputs"
    _write_json(outputs / "one.json", _make_artifact())
    _write_json(outputs / "two.json", _make_artifact())
    gold_path = tmp_path / "gold.json"
    _write_json(gold_path, _gold({"patient_id": "p1", "pathogens": []}))
    with unittest.TestCase().assertRaisesRegex(EvaluationError, "Duplicate pruning artifact"):
        evaluate_directory(outputs, gold_path)


def test_subset_request_adds_provenance_dependencies(tmp_path: Path) -> None:
    report = _run(
        tmp_path,
        _make_artifact(),
        _gold({"patient_id": "p1", "pathogens": []}),
        arm_names=[PRIMARY_ARM],
    )
    assert report["arm_names"] == [
        BASELINE_ARM,
        *CUMULATIVE_ARMS[: CUMULATIVE_ARMS.index(PRIMARY_ARM) + 1],
    ]


def test_protocol_markdown_and_json_writers(tmp_path: Path) -> None:
    protocol = PACKAGE_ROOT / "config" / "evaluation_protocol_v1.json"
    report = _run(
        tmp_path,
        _make_artifact(),
        _gold({"patient_id": "p1", "pathogens": ["Klebsiella pneumoniae"]}),
        protocol_path=protocol,
    )
    markdown = render_markdown(report)
    assert "Cumulative step deltas" in markdown
    assert "Pruned TP/FP" in markdown
    assert "evaluation-only" in markdown
    md_path = write_markdown(report, tmp_path / "report.md")
    json_path = write_json(report, tmp_path / "report.json")
    assert md_path.read_text(encoding="utf-8") == markdown
    loaded = json.loads(json_path.read_text(encoding="utf-8"))
    assert loaded["schema_version"] == "rag_re_clinical_pruning.evaluation.v1"


def test_input_artifact_hash_is_stable_across_evaluation(tmp_path: Path) -> None:
    outputs = tmp_path / "outputs"
    artifact_path = outputs / "p1.json"
    _write_json(artifact_path, _make_artifact())
    gold_path = tmp_path / "gold.json"
    _write_json(gold_path, _gold({"patient_id": "p1", "pathogens": []}))
    before = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
    report = evaluate_directory(outputs, gold_path)
    after = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
    assert before == after
    assert report["artifacts"]["sha256_by_patient"]["p1"] == before


def test_evaluation_cli_writes_both_reports(tmp_path: Path) -> None:
    outputs = tmp_path / "outputs"
    _write_json(outputs / "p1.json", _make_artifact())
    gold_path = tmp_path / "gold.json"
    _write_json(
        gold_path,
        _gold({"patient_id": "p1", "pathogens": ["Klebsiella pneumoniae"]}),
    )
    json_out = tmp_path / "reports" / "evaluation.json"
    markdown_out = tmp_path / "reports" / "evaluation.md"
    status = evaluate_cli_main(
        [
            "--outputs",
            str(outputs),
            "--gold",
            str(gold_path),
            "--protocol",
            str(PACKAGE_ROOT / "config" / "evaluation_protocol_v1.json"),
            "--json-out",
            str(json_out),
            "--markdown-out",
            str(markdown_out),
            "--quiet",
        ]
    )
    assert status == 0
    assert json.loads(json_out.read_text(encoding="utf-8"))["schema_version"] == (
        "rag_re_clinical_pruning.evaluation.v1"
    )
    assert markdown_out.read_text(encoding="utf-8").startswith(
        "# RAG_re clinical pruning evaluation"
    )


def test_generation_modules_do_not_import_evaluation_or_gold(tmp_path: Path) -> None:  # noqa: ARG001
    package = PACKAGE_ROOT / "rag_re_clinical_pruning"
    for filename in ("engine.py", "rules.py", "source_adapter.py", "states.py"):
        source = (package / filename).read_text(encoding="utf-8").lower()
        assert "import evaluation" not in source
        assert "from .evaluation" not in source
        assert "from rag_re_clinical_pruning.evaluation" not in source


def load_tests(loader, tests, pattern):  # noqa: ARG001
    """Run the dependency-free function tests with an isolated temp directory."""

    suite = unittest.TestSuite()
    for name, function in sorted(globals().items()):
        if not name.startswith("test_") or not callable(function):
            continue

        def run_test(function=function) -> None:
            with tempfile.TemporaryDirectory() as directory:
                function(Path(directory))

        suite.addTest(unittest.FunctionTestCase(run_test, description=name))
    return suite

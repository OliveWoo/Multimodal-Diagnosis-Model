from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

from rag_re_clinical_pruning_abc.engine import (
    _run_projected_patient,
    arm_definitions,
    load_config,
    run_patient,
)
from rag_re_clinical_pruning_abc.source_adapter import (
    SourceValidationError,
    load_and_join_patient,
    load_pruning_manifest_anchors,
    sha256_file,
)
from rag_re_clinical_pruning_abc.states import (
    ARM_ORDER,
    predicate_values,
    tri_and,
    tri_or,
    tri_two_of_three,
)


PROJECT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT / "config" / "post_l5_abc_policy_v1.json"


def candidate(
    name: str,
    *,
    baseline: bool,
    l5: str,
    d: str = "D_PASS",
    corr=(True, True, True),
    raw=(True, True, True),
) -> dict:
    raw_a, raw_b, raw_c = raw
    raw_a_audit = (
        {"positive": raw_a, "status": "ok", "support_count": 3 if raw_a else 2, "judgeable_count": 10, "support_ratio": 0.3 if raw_a else 0.2}
        if raw_a is not None
        else {"positive": None, "status": "not_requested", "support_count": 0, "judgeable_count": 0, "support_ratio": None}
    )
    raw_b_audit = {
        "positive": raw_b,
        "status": "ok",
        "site_aligned": True if raw_b else None,
        "trigger_count": 1 if raw_b else 0,
    }
    raw_c_audit = {
        "positive": raw_c,
        "status": "ok" if raw_c is not None else "insufficient_input",
        "site_aligned": True if raw_c is not None else None,
        "support_record_count": 1 if raw_c is not False else 0,
    }
    corr_c = corr[2]
    return {
        "organism_name": name,
        "canonical_organism_name": name,
        "baseline_selected": baseline,
        "review_tier": None if baseline else "review_context_needed",
        "l5_state": l5,
        "context": {
            "D_state": d,
            "F_state": "F_UNKNOWN",
            "clinical_C_grade": "C2" if corr_c is True else "CNEG" if corr_c is False else "C0",
        },
        "signals": {
            "CORR": dict(zip(("A", "B", "C"), corr)),
            "RAW": dict(zip(("A", "B", "C"), raw)),
        },
        "signal_audit": {"RAW_A": raw_a_audit, "RAW_B": raw_b_audit, "RAW_C": raw_c_audit},
    }


def joined(candidates: list[dict]) -> dict:
    e0 = sum(item["baseline_selected"] for item in candidates)
    return {
        "patient_id": "fixture",
        "sources": {
            "pruning": {
                "path": "pruning.json",
                "sha256": "a",
                "schema_version": "rag_re_clinical_pruning.output.v1",
                "decision_sha256": "d",
                "batch_manifest_path": "manifest.json",
                "batch_manifest_sha256": "m",
            },
            "clinical": {
                "path": "clinical.json",
                "sha256": "b",
                "schema_version": "rag_re_clinical.output.v1",
            },
            "raw": {
                "path": "raw.json",
                "sha256": "c",
                "schema_version": "rag_re.output.v1",
                "pipeline_fingerprint": "pipe",
                "config_hash": "cfg",
                "candidate_pool_sha256": "pool",
            },
        },
        "counts": {
            "candidate_count": len(candidates),
            "e0_count": e0,
            "non_e0_count": len(candidates) - e0,
        },
        "candidates": candidates,
    }


class TriStateTests(unittest.TestCase):
    def test_and_or_truth_tables(self):
        self.assertIs(tri_and((True, True)), True)
        self.assertIs(tri_and((True, None)), None)
        self.assertIs(tri_and((False, None)), False)
        self.assertIs(tri_or((False, False)), False)
        self.assertIs(tri_or((False, None)), None)
        self.assertIs(tri_or((True, None)), True)

    def test_two_of_three_truth_table(self):
        cases = {
            (True, True, None): True,
            (True, False, None): None,
            (False, False, None): False,
            (True, None, None): None,
            (False, None, None): None,
            (False, False, True): False,
        }
        for values, expected in cases.items():
            self.assertIs(tri_two_of_three(*values), expected)

    def test_predicates_use_and_not_or(self):
        values = predicate_values(True, False, None)
        self.assertIs(values["AB"], False)
        self.assertIs(values["AC"], None)
        self.assertIs(values["ABC"], False)
        self.assertIs(values["OR"], True)


class EngineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config, cls.meta = load_config(CONFIG_PATH)

    def test_frozen_arm_count_and_order(self):
        self.assertEqual(len(ARM_ORDER), 79)
        self.assertEqual(tuple(arm_definitions()), ARM_ORDER)
        self.assertEqual(ARM_ORDER[0], "L5_FROZEN")

    def test_transition_scopes_and_d_gate(self):
        rows = [
            candidate("rsc-pass", baseline=False, l5="NO_RESCUE", d="D_PASS"),
            candidate("rsc-guard", baseline=False, l5="NO_RESCUE", d="D_GUARDED"),
            candidate("rsc-block", baseline=False, l5="NO_RESCUE", d="D_BLOCK"),
            candidate("prune-pass", baseline=True, l5="PRUNE", d="D_PASS"),
            candidate("prune-guard", baseline=True, l5="PRUNE", d="D_UNKNOWN"),
            candidate("manual-pass", baseline=True, l5="MANUAL", d="D_PASS"),
            candidate("keep", baseline=True, l5="KEEP", d="D_PASS"),
        ]
        output = _run_projected_patient(joined(rows), self.config, config_meta=self.meta)
        by_name = {row["canonical_organism_name"]: row for row in output["candidates"]}
        self.assertEqual(by_name["rsc-pass"]["arm_states"]["S_CORR_RSC_A"]["state"], "RESCUE")
        self.assertEqual(by_name["rsc-guard"]["arm_states"]["S_CORR_RSC_A"]["state"], "MANUAL")
        self.assertEqual(by_name["rsc-block"]["arm_states"]["S_CORR_RSC_A"]["state"], "NO_RESCUE")
        self.assertEqual(by_name["prune-pass"]["arm_states"]["S_CORR_P_RST_A"]["state"], "KEEP")
        self.assertEqual(by_name["prune-guard"]["arm_states"]["S_CORR_P_RST_A"]["state"], "MANUAL")
        self.assertEqual(by_name["manual-pass"]["arm_states"]["S_CORR_M_RST_A"]["state"], "KEEP")
        self.assertEqual(by_name["keep"]["arm_states"]["S_CORR_RSC_A"]["state"], "KEEP")

    def test_unknown_never_auto_transitions(self):
        row = candidate(
            "unknown",
            baseline=False,
            l5="NO_RESCUE",
            corr=(None, True, False),
            raw=(None, True, True),
        )
        output = _run_projected_patient(joined([row]), self.config)
        state = output["candidates"][0]["arm_states"]["S_CORR_RSC_A"]
        self.assertEqual(state["state"], "NO_RESCUE")
        self.assertEqual(state["rule_ledger"][-1]["action"], "ABSTAIN")

    def test_ungated_is_diagnostic_and_can_cross_block(self):
        row = candidate("blocked", baseline=False, l5="NO_RESCUE", d="D_BLOCK")
        output = _run_projected_patient(joined([row]), self.config)
        item = output["candidates"][0]
        self.assertEqual(item["arm_states"]["S_CORR_RSC_A"]["state"], "NO_RESCUE")
        self.assertEqual(item["arm_states"]["X_CORR_RSC_UNGATED_A"]["state"], "RESCUE")
        self.assertEqual(output["arms"]["X_CORR_RSC_UNGATED_A"]["eligible_primary"], False)

    def test_cumulative_second_stage_equals_all_scope_two_of_three(self):
        rows = [
            candidate("r", baseline=False, l5="NO_RESCUE", corr=(True, True, None)),
            candidate("p", baseline=True, l5="PRUNE", corr=(True, True, None)),
            candidate("m", baseline=True, l5="MANUAL", corr=(True, True, None)),
        ]
        output = _run_projected_patient(joined(rows), self.config)
        for row in output["candidates"]:
            scoped = (
                "S_CORR_RSC_2OF3"
                if not row["baseline_selected"]
                else "S_CORR_P_RST_2OF3"
                if row["l5_state"] == "PRUNE"
                else "S_CORR_M_RST_2OF3"
            )
            self.assertEqual(
                row["arm_states"]["K_CORR_2_2OF3"]["state"],
                row["arm_states"][scoped]["state"],
            )

    def test_cumulative_ledger_attributes_only_first_effective_transition(self):
        row = candidate("r", baseline=False, l5="NO_RESCUE", corr=(True, True, True))
        output = _run_projected_patient(joined([row]), self.config)
        ledger = output["candidates"][0]["arm_states"]["K_CORR_3_OR_STRESS"]["rule_ledger"]
        self.assertEqual([entry["action"] for entry in ledger[1:]], ["RESCUE", "NO_CHANGE", "NO_CHANGE"])
        self.assertEqual(
            [entry["reason_code"] for entry in ledger[2:]],
            ["ALREADY_TRANSITIONED_BY_EARLIER_STAGE"] * 2,
        )

    def test_l5_auto_is_never_removed_by_post_l5_arms(self):
        rows = [
            candidate("keep", baseline=True, l5="KEEP", corr=(False, False, False), raw=(False, False, False)),
            candidate("rescue", baseline=False, l5="RESCUE", corr=(False, False, False), raw=(False, False, False)),
        ]
        output = _run_projected_patient(joined(rows), self.config)
        for row in output["candidates"]:
            for arm in ARM_ORDER:
                self.assertIn(row["arm_states"][arm]["state"], {"KEEP", "RESCUE"})

    def test_public_api_rejects_gold_like_extra(self):
        envelope = joined([candidate("x", baseline=False, l5="NO_RESCUE")])
        envelope["candidates"][0]["gold"] = True
        with self.assertRaises(SourceValidationError):
            run_patient(envelope, self.config, config_meta=self.meta)

    def test_public_api_rejects_integer_bool(self):
        envelope = joined([candidate("x", baseline=False, l5="NO_RESCUE")])
        envelope["candidates"][0]["signals"]["CORR"]["B"] = 1
        with self.assertRaises(SourceValidationError):
            run_patient(envelope, self.config, config_meta=self.meta)

    def test_public_api_rejects_config_meta_extra(self):
        envelope = joined([candidate("x", baseline=False, l5="NO_RESCUE")])
        meta = dict(self.meta)
        meta["answer"] = "leak"
        with self.assertRaises(SourceValidationError):
            run_patient(envelope, self.config, config_meta=meta)

    def test_decision_is_deterministic_excluding_timestamp(self):
        envelope = joined([candidate("x", baseline=False, l5="NO_RESCUE")])
        one = _run_projected_patient(envelope, self.config, generated_at_utc="one")
        two = _run_projected_patient(envelope, self.config, generated_at_utc="two")
        self.assertEqual(one["decision_sha256"], two["decision_sha256"])


class SourceAdapterTests(unittest.TestCase):
    def _write_fixture(self, directory: Path, *, raw_extra: dict | None = None):
        raw_candidate = {
            "organism_name": "Fixture organism",
            "canonical_organism_name": "Fixture organism",
            "baseline_selected": False,
            "modules": {
                "A": {
                    "positive": True,
                    "status": "ok",
                    "support_count": 3,
                    "judgeable_count": 10,
                    "support_ratio": 0.3,
                },
                "B": {"positive": True, "status": "ok", "site_aligned": True, "triggers": ["x"]},
                "C": {
                    "positive": True,
                    "status": "ok",
                    "site_aligned": True,
                    "supports": [{"module": "culture"}],
                },
            },
        }
        if raw_extra:
            raw_candidate.update(raw_extra)
        raw = {
            "schema_version": "rag_re.output.v1",
            "run_status": "complete",
            "patient_id": "1",
            "pipeline_fingerprint": "pipe",
            "config_hash": "cfg",
            "input_meta": {"candidate_pool_sha256": "pool"},
            "candidates": [raw_candidate],
        }
        raw_path = directory / "raw.json"
        raw_path.write_text(json.dumps(raw), encoding="utf-8")
        clinical = {
            "schema_version": "rag_re_clinical.output.v1",
            "run_status": "complete",
            "patient_id": "1",
            "provenance": {"frozen_rag_re_sha256": sha256_file(raw_path)},
            "candidates": [
                {
                    "organism_name": "Fixture organism",
                    "canonical_organism_name": "Fixture organism",
                    "baseline_selected": False,
                    "modules": {
                        "D": {"state": "D_PASS"},
                        "F": {"state": "F_UNKNOWN"},
                        "E": {"inputs": {"A": True}},
                        "B_STRICT": {"positive": True},
                        "C": {"grade": "C2"},
                    },
                }
            ],
        }
        clinical_path = directory / "clinical.json"
        clinical_path.write_text(json.dumps(clinical), encoding="utf-8")
        pruning = {
            "schema_version": "rag_re_clinical_pruning.output.v1",
            "run_status": "complete",
            "patient_id": "1",
            "decision_sha256": "d" * 64,
            "source": {"sha256": sha256_file(clinical_path)},
            "counts": {"candidate_count": 1, "e0_count": 0, "non_e0_count": 1},
            "candidates": [
                {
                    "organism_name": "Fixture organism",
                    "canonical_organism_name": "Fixture organism",
                    "baseline_selected": False,
                    "review_tier": "review_context_needed",
                    "input_features": {
                        "D": {"state": "D_PASS"},
                        "F": {"state": "F_UNKNOWN"},
                        "E": {"A_positive": True},
                        "B_STRICT": {"positive": True},
                        "C": {"grade": "C2"},
                    },
                    "arm_states": {
                        "L5_PLUS_DIRECT_CONVERGENT_RESCUE": {"state": "NO_RESCUE"}
                    },
                }
            ],
        }
        pruning_path = directory / "pruning.json"
        pruning_path.write_text(json.dumps(pruning), encoding="utf-8")
        manifest = {
            "schema_version": "rag_re_clinical_pruning.batch_manifest.v1",
            "output_schema_version": "rag_re_clinical_pruning.output.v1",
            "run_status": "complete",
            "primary_arm": "L5_PLUS_DIRECT_CONVERGENT_RESCUE",
            "patient_count": 1,
            "artifacts": [
                {
                    "patient_id": "1",
                    "path": str(pruning_path.resolve()),
                    "sha256": sha256_file(pruning_path),
                    "decision_sha256": "d" * 64,
                    "source_path": str(clinical_path.resolve()),
                    "source_sha256": sha256_file(clinical_path),
                    "candidate_count": 1,
                    "e0_count": 0,
                    "non_e0_count": 1,
                }
            ],
        }
        (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return pruning_path, clinical_path, raw_path

    def _anchor(self, directory: Path, paths):
        anchors, _ = load_pruning_manifest_anchors(
            directory, {"1": paths[0]}, {"1": paths[1]}
        )
        return anchors["1"]

    def _refresh_chain(self, directory: Path, paths):
        clinical = json.loads(paths[1].read_text(encoding="utf-8"))
        clinical["provenance"]["frozen_rag_re_sha256"] = sha256_file(paths[2])
        paths[1].write_text(json.dumps(clinical), encoding="utf-8")
        pruning = json.loads(paths[0].read_text(encoding="utf-8"))
        pruning["source"]["sha256"] = sha256_file(paths[1])
        paths[0].write_text(json.dumps(pruning), encoding="utf-8")
        manifest_path = directory / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        entry = manifest["artifacts"][0]
        entry["sha256"] = sha256_file(paths[0])
        entry["source_sha256"] = sha256_file(paths[1])
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        return self._anchor(directory, paths)

    def test_exact_join_and_projection(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = self._write_fixture(Path(directory))
            anchor = self._anchor(Path(directory), paths)
            result = load_and_join_patient(*paths, anchor=anchor)
            self.assertEqual(result["patient_id"], "1")
            self.assertEqual(result["candidates"][0]["signals"]["RAW"], {"A": True, "B": True, "C": True})

    def test_raw_sha_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = self._write_fixture(Path(directory))
            anchor = self._anchor(Path(directory), paths)
            raw = json.loads(paths[2].read_text(encoding="utf-8"))
            raw["pipeline_fingerprint"] = "tampered"
            paths[2].write_text(json.dumps(raw), encoding="utf-8")
            with self.assertRaises(SourceValidationError):
                load_and_join_patient(*paths, anchor=anchor)

    def test_manifest_rejects_pruning_l5_tamper(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = self._write_fixture(Path(directory))
            anchor = self._anchor(Path(directory), paths)
            pruning = json.loads(paths[0].read_text(encoding="utf-8"))
            pruning["candidates"][0]["arm_states"]["L5_PLUS_DIRECT_CONVERGENT_RESCUE"]["state"] = "RESCUE"
            paths[0].write_text(json.dumps(pruning), encoding="utf-8")
            with self.assertRaises(SourceValidationError):
                load_and_join_patient(*paths, anchor=anchor)

    def test_exact_name_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = self._write_fixture(Path(directory))
            raw = json.loads(paths[2].read_text(encoding="utf-8"))
            raw["candidates"][0]["canonical_organism_name"] = "Other organism"
            paths[2].write_text(json.dumps(raw), encoding="utf-8")
            clinical = json.loads(paths[1].read_text(encoding="utf-8"))
            clinical["provenance"]["frozen_rag_re_sha256"] = sha256_file(paths[2])
            paths[1].write_text(json.dumps(clinical), encoding="utf-8")
            pruning = json.loads(paths[0].read_text(encoding="utf-8"))
            pruning["source"]["sha256"] = sha256_file(paths[1])
            paths[0].write_text(json.dumps(pruning), encoding="utf-8")
            anchor = self._refresh_chain(Path(directory), paths)
            with self.assertRaises(SourceValidationError):
                load_and_join_patient(*paths, anchor=anchor)

    def test_raw_gold_extra_is_not_projected(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = self._write_fixture(Path(directory), raw_extra={"gold": True, "answer": "leak"})
            anchor = self._anchor(Path(directory), paths)
            result = load_and_join_patient(*paths, anchor=anchor)
            candidate_result = result["candidates"][0]
            self.assertNotIn("gold", candidate_result)
            self.assertNotIn("answer", candidate_result)

    def test_public_engine_rebinds_projection_to_locked_files(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = self._write_fixture(Path(directory))
            anchor = self._anchor(Path(directory), paths)
            projection = load_and_join_patient(*paths, anchor=anchor)
            config, meta = load_config(CONFIG_PATH)
            self.assertEqual(run_patient(projection, config, config_meta=meta)["run_status"], "complete")
            tampered = copy.deepcopy(projection)
            tampered["candidates"][0]["signals"]["CORR"]["B"] = False
            with self.assertRaises(SourceValidationError):
                run_patient(tampered, config, config_meta=meta)


if __name__ == "__main__":
    unittest.main()

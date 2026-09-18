from __future__ import annotations

import inspect
import json
import tempfile
import unittest
from pathlib import Path

from rag_re_clinical_pruning.cli import build_parser
from rag_re_clinical_pruning.engine import (
    EXPECTED_ARMS,
    OUTPUT_SCHEMA_VERSION,
    PRIMARY_ARM,
    PolicyContractError,
    load_policy_config,
    run_patient,
)
from rag_re_clinical_pruning.source_adapter import (
    SourceContractError,
    validate_source_artifact,
)
from rag_re_clinical_pruning.states import (
    RuleVote,
    State,
    StateContractError,
    apply_vote,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "pruning_policy_v1.json"


def candidate(
    name: str,
    *,
    e0: bool,
    d: str = "D_PASS",
    b: bool | None = False,
    absolute: bool | None = False,
    relative: bool | None = False,
    c: str = "C0",
    tier: str | None = None,
    a: bool | None = None,
) -> dict:
    return {
        "organism_name": name,
        "canonical_organism_name": name,
        "baseline_selected": e0,
        "review_tier": tier,
        "source_features": {
            "D": {"state": d},
            "B_STRICT": {
                "positive": b,
                "absolute_strength_axis": absolute,
                "relative_strength_axis": relative,
            },
            "C": {"grade": c},
            "E": {"A_positive": a},
            "F": {"state": "F_UNKNOWN"},
        },
    }


def projected(candidates: list[dict]) -> dict:
    e0_count = sum(row["baseline_selected"] for row in candidates)
    return {
        "patient_id": "synthetic",
        "source": {
            "path": "synthetic.json",
            "sha256": "a" * 64,
            "sha256_basis": "canonical_projected_inputs_v1",
            "schema_version": "rag_re_clinical.output.v1",
            "pipeline_version": "test",
            "artifact_sha256_basis": "test",
        },
        "candidates": candidates,
        "counts": {
            "candidate_count": len(candidates),
            "e0_count": e0_count,
            "non_e0_count": len(candidates) - e0_count,
        },
    }


def state(artifact: dict, organism: str, arm: str) -> str:
    row = next(row for row in artifact["candidates"] if row["organism_name"] == organism)
    return row["arm_states"][arm]["state"]


class TestFrozenPolicy(unittest.TestCase):
    def test_policy_has_exact_frozen_arm_order_and_primary(self) -> None:
        policy = load_policy_config(CONFIG)
        actual = [
            (
                arm["name"],
                arm["family"],
                arm["analysis_role"],
                tuple(arm["rule_sequence"]),
            )
            for arm in policy.arms
        ]
        self.assertEqual(actual, list(EXPECTED_ARMS))
        self.assertEqual(policy.raw["primary_arm"], PRIMARY_ARM)
        self.assertNotIn(
            "D_BLOCK_DIAGNOSTIC",
            {
                rule
                for arm in policy.arms
                if arm["family"] == "cumulative"
                for rule in arm["rule_sequence"]
            },
        )

    def test_policy_drift_is_rejected(self) -> None:
        raw = json.loads(CONFIG.read_text(encoding="utf-8"))
        raw["arms"][-1]["rule_sequence"].reverse()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad.json"
            path.write_text(json.dumps(raw), encoding="utf-8")
            with self.assertRaises(PolicyContractError):
                load_policy_config(path)


class TestStateContract(unittest.TestCase):
    def test_unknown_never_becomes_negative(self) -> None:
        vote = RuleVote("X", None, None, "UNKNOWN", {})
        state_e0, _ = apply_vote(State.KEEP, vote, baseline_selected=True)
        state_non_e0, _ = apply_vote(State.NO_RESCUE, vote, baseline_selected=False)
        state_pruned, _ = apply_vote(State.PRUNE, vote, baseline_selected=True)
        self.assertEqual(state_e0, State.KEEP)
        self.assertEqual(state_non_e0, State.MANUAL)
        self.assertEqual(state_pruned, State.MANUAL)

    def test_state_domains_cannot_cross(self) -> None:
        with self.assertRaises(StateContractError):
            apply_vote(
                State.KEEP,
                RuleVote("BAD", True, State.RESCUE, "BAD", {}),
                baseline_selected=True,
            )
        with self.assertRaises(StateContractError):
            apply_vote(
                State.NO_RESCUE,
                RuleVote("BAD", True, State.PRUNE, "BAD", {}),
                baseline_selected=False,
            )


class TestRulesAndArms(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        rows = [
            candidate("weak_c0", e0=True, c="C0"),
            candidate("weak_c1", e0=True, c="C1"),
            candidate("weak_c2", e0=True, c="C2"),
            candidate("weak_c3", e0=True, c="C3"),
            candidate("guard_c2", e0=True, d="D_GUARDED", b=True, c="C2"),
            candidate("guard_c3", e0=True, d="D_GUARDED", b=True, c="C3"),
            candidate("dblock", e0=True, d="D_BLOCK", b=True, c="C0"),
            candidate("weak_axis_unknown", e0=True, relative=None, c="C0"),
            candidate("direct", e0=False, d="D_GUARDED", c="C3"),
            candidate("direct_block", e0=False, d="D_BLOCK", c="C3"),
            candidate("convergent", e0=False, b=True, c="C2"),
            candidate("convergent_unknown", e0=False, b=None, c="C2"),
            candidate(
                "high_strict",
                e0=False,
                d="D_GUARDED",
                b=True,
                absolute=True,
                c="C0",
                tier="review_high_priority",
            ),
            candidate(
                "high_any",
                e0=False,
                d="D_GUARDED",
                b=False,
                absolute=True,
                c="C0",
                tier="review_high_priority",
            ),
            candidate(
                "high_any_unknown",
                e0=False,
                d="D_GUARDED",
                b=False,
                absolute=None,
                relative=False,
                c="C0",
                tier="review_high_priority",
            ),
            candidate(
                "context_a",
                e0=False,
                b=True,
                c="C1",
                tier="review_context_needed",
                a=True,
            ),
            candidate(
                "context_a_unknown",
                e0=False,
                b=True,
                c="C1",
                tier="review_context_needed",
                a=None,
            ),
        ]
        cls.artifact = run_patient(
            projected(rows),
            load_policy_config(CONFIG),
            generated_at_utc="2026-08-11T00:00:00Z",
        )

    def test_weak_raw_prune_and_evidence_restore_ladder(self) -> None:
        self.assertEqual(state(self.artifact, "weak_c0", "L1_WEAK_RAW"), "PRUNE")
        self.assertEqual(state(self.artifact, "weak_c0", "L3_CONTEXT_REVIEW"), "PRUNE")
        self.assertEqual(state(self.artifact, "weak_c1", "L1_WEAK_RAW"), "PRUNE")
        self.assertEqual(state(self.artifact, "weak_c1", "L3_CONTEXT_REVIEW"), "MANUAL")
        for name in ("weak_c2", "weak_c3"):
            self.assertEqual(state(self.artifact, name, "L1_WEAK_RAW"), "PRUNE")
            self.assertEqual(state(self.artifact, name, "L2_DIRECT_RESTORE"), "KEEP")

    def test_guard_and_dblock_are_distinct(self) -> None:
        self.assertEqual(state(self.artifact, "guard_c2", "L4_PLUS_GUARD_REVIEW"), "MANUAL")
        self.assertEqual(state(self.artifact, "guard_c3", "L4_PLUS_GUARD_REVIEW"), "KEEP")
        self.assertEqual(state(self.artifact, "dblock", "S_D_BLOCK_DIAGNOSTIC"), "MANUAL")
        self.assertEqual(state(self.artifact, "dblock", PRIMARY_ARM), "KEEP")

    def test_unknown_weak_axis_does_not_prune_e0(self) -> None:
        self.assertEqual(state(self.artifact, "weak_axis_unknown", "S_Q_WEAK_RAW"), "KEEP")
        ledger = next(
            row for row in self.artifact["candidates"] if row["organism_name"] == "weak_axis_unknown"
        )["arm_states"]["S_Q_WEAK_RAW"]["rule_ledger"]
        self.assertIsNone(ledger[-1]["matched"])
        self.assertEqual(ledger[-1]["action"], "NO_CHANGE")

    def test_direct_and_convergent_rescue(self) -> None:
        self.assertEqual(state(self.artifact, "direct", PRIMARY_ARM), "RESCUE")
        self.assertEqual(state(self.artifact, "direct_block", PRIMARY_ARM), "NO_RESCUE")
        self.assertEqual(state(self.artifact, "convergent", PRIMARY_ARM), "RESCUE")
        self.assertEqual(state(self.artifact, "convergent_unknown", PRIMARY_ARM), "MANUAL")

    def test_exploratory_rescues_enter_only_their_later_arms(self) -> None:
        self.assertEqual(state(self.artifact, "high_strict", PRIMARY_ARM), "NO_RESCUE")
        self.assertEqual(state(self.artifact, "high_strict", "L6_PLUS_HIGH_STRICT"), "RESCUE")
        self.assertEqual(state(self.artifact, "high_any", "L6_PLUS_HIGH_STRICT"), "NO_RESCUE")
        self.assertEqual(state(self.artifact, "high_any", "L7_PLUS_HIGH_ANY"), "RESCUE")
        self.assertEqual(state(self.artifact, "high_any_unknown", "L7_PLUS_HIGH_ANY"), "MANUAL")
        self.assertEqual(state(self.artifact, "context_a", "L7_PLUS_HIGH_ANY"), "NO_RESCUE")
        self.assertEqual(state(self.artifact, "context_a", "L8_PLUS_CONTEXT_A"), "RESCUE")
        self.assertEqual(state(self.artifact, "context_a_unknown", "L8_PLUS_CONTEXT_A"), "MANUAL")

    def test_output_schema_ledgers_and_root_lists_are_complete(self) -> None:
        artifact = self.artifact
        arm_names = [row[0] for row in EXPECTED_ARMS]
        self.assertEqual(artifact["schema_version"], OUTPUT_SCHEMA_VERSION)
        self.assertEqual(list(artifact["arms"]), arm_names)
        self.assertFalse(artifact["generation_contract"]["reads_gold"])
        for row in artifact["candidates"]:
            self.assertEqual(list(row["arm_states"]), arm_names)
            for arm in load_policy_config(CONFIG).arms:
                result = row["arm_states"][arm["name"]]
                self.assertEqual(
                    len(result["rule_ledger"]), len(arm["rule_sequence"]) + 1
                )
        for arm in artifact["arms"].values():
            self.assertEqual(arm["auto_positive_count"], len(arm["auto_positive_pathogens"]))
            self.assertEqual(
                arm["review_inclusive_count"], len(arm["review_inclusive_pathogens"])
            )

    def test_e0_baseline_is_exact_in_both_baseline_arms(self) -> None:
        expected = [
            row["organism_name"]
            for row in self.artifact["candidates"]
            if row["baseline_selected"]
        ]
        self.assertEqual(self.artifact["arms"]["S0_E0"]["auto_positive_pathogens"], expected)
        self.assertEqual(self.artifact["arms"]["L0_E0"]["auto_positive_pathogens"], expected)


class TestPredictionGoldIsolation(unittest.TestCase):
    def _raw_source(self) -> dict:
        return {
            "schema_version": "rag_re_clinical.output.v1",
            "pipeline_version": "test",
            "run_status": "complete",
            "patient_id": "x",
            "answer": "GOLD_MARKER",
            "candidates": [
                {
                    "organism_name": "Organism x",
                    "canonical_organism_name": "Organism x",
                    "baseline_selected": True,
                    "review_tier": None,
                    "gold_label": "GOLD_MARKER",
                    "rationale": "GOLD_MARKER",
                    "modules": {
                        "D": {"state": "D_PASS", "reason": "GOLD_MARKER"},
                        "B_STRICT": {
                            "positive": False,
                            "absolute_strength_axis": False,
                            "relative_strength_axis": False,
                            "reason": "GOLD_MARKER",
                        },
                        "C": {"grade": "C0", "records": ["GOLD_MARKER"]},
                        "E": {"inputs": {"A": None}, "reason": "GOLD_MARKER"},
                        "F": {"state": "F_UNKNOWN", "reason": "GOLD_MARKER"},
                    },
                }
            ],
            "arms": {
                "CL0_BASELINE": {
                    "predicted_pathogens": ["Organism x"],
                    "predicted_count": 1,
                }
            },
        }

    def test_adapter_projects_only_rule_inputs(self) -> None:
        projection = validate_source_artifact(self._raw_source())
        encoded = json.dumps(projection, ensure_ascii=False)
        self.assertNotIn("GOLD_MARKER", encoded)
        self.assertEqual(
            set(projection["candidates"][0]["source_features"]),
            {"D", "B_STRICT", "C", "E", "F"},
        )
        self.assertEqual(projection["source"]["path"], "<memory>")
        self.assertEqual(
            projection["source"]["sha256_basis"],
            "canonical_projected_inputs_v1",
        )
        run_patient(projection, load_policy_config(CONFIG))
        changed_only_gold = self._raw_source()
        changed_only_gold["answer"] = "DIFFERENT_GOLD"
        changed_only_gold["candidates"][0]["gold_label"] = "DIFFERENT_GOLD"
        self.assertEqual(
            projection["source"]["sha256"],
            validate_source_artifact(changed_only_gold)["source"]["sha256"],
        )

    def test_generation_api_has_no_gold_or_safety_inputs(self) -> None:
        parameters = set(inspect.signature(run_patient).parameters)
        forbidden = {"gold", "answer", "labels", "recall_safety", "configured_recall"}
        self.assertTrue(parameters.isdisjoint(forbidden))
        option_dests = {
            action.dest for action in build_parser()._actions  # noqa: SLF001 - contract test
        }
        self.assertTrue(option_dests.isdisjoint(forbidden))

    def test_cl0_mismatch_is_rejected(self) -> None:
        raw = self._raw_source()
        raw["arms"]["CL0_BASELINE"]["predicted_pathogens"] = []
        raw["arms"]["CL0_BASELINE"]["predicted_count"] = 0
        with self.assertRaises(SourceContractError):
            validate_source_artifact(raw)

    def test_schema_drift_in_review_tier_or_e_inputs_is_rejected(self) -> None:
        raw = self._raw_source()
        raw["candidates"][0]["review_tier"] = "new_unfrozen_tier"
        with self.assertRaises(SourceContractError):
            validate_source_artifact(raw)
        raw = self._raw_source()
        raw["candidates"][0]["modules"]["E"]["inputs"] = "malformed"
        with self.assertRaises(SourceContractError):
            validate_source_artifact(raw)

    def test_public_engine_revalidates_projected_feature_types_and_enums(self) -> None:
        policy = load_policy_config(CONFIG)
        mutations = (
            ("D", "state", "D_WAT"),
            ("C", "grade", "C9"),
            ("F", "state", "F_WAT"),
            ("B_STRICT", "positive", 1),
            ("B_STRICT", "absolute_strength_axis", "false"),
        )
        for module, field, invalid in mutations:
            with self.subTest(module=module, field=field, invalid=invalid):
                envelope = projected([candidate("bad", e0=False, b=True, c="C2")])
                envelope["candidates"][0]["source_features"][module][field] = invalid
                with self.assertRaises(ValueError):
                    run_patient(envelope, policy)

    def test_public_engine_rejects_extra_answer_like_projection_keys(self) -> None:
        envelope = projected([candidate("bad", e0=True)])
        envelope["candidates"][0]["gold_label"] = True
        with self.assertRaises(ValueError):
            run_patient(envelope, load_policy_config(CONFIG))


if __name__ == "__main__":
    unittest.main()

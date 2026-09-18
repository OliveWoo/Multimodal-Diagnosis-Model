from __future__ import annotations

import copy
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from .rules import RULE_EVALUATORS, evaluate_rule
from .source_adapter import (
    ALLOWED_C_GRADES,
    ALLOWED_D_STATES,
    ALLOWED_F_STATES,
    ALLOWED_REVIEW_TIERS,
    load_json,
    sha256_file,
)
from .states import (
    AUTO_POSITIVE_STATES,
    REVIEW_INCLUSIVE_STATES,
    State,
    apply_vote,
    initial_state,
    validate_state,
)


PathLike = str | Path
JSONDict = dict[str, Any]

CONFIG_SCHEMA_VERSION = "rag_re_clinical_pruning.config.v1"
OUTPUT_SCHEMA_VERSION = "rag_re_clinical_pruning.output.v1"
BATCH_MANIFEST_SCHEMA_VERSION = "rag_re_clinical_pruning.batch_manifest.v1"
PIPELINE_VERSION = "0.1.0"

PRIMARY_ARM = "L5_PLUS_DIRECT_CONVERGENT_RESCUE"

EXPECTED_ARMS: tuple[tuple[str, str, str, tuple[str, ...]], ...] = (
    ("S0_E0", "standalone", "baseline", ()),
    ("S_Q_WEAK_RAW", "standalone", "confirmatory", ("Q_WEAK_RAW",)),
    ("S_Q_GUARD", "standalone", "confirmatory", ("Q_GUARD",)),
    (
        "S_D_BLOCK_DIAGNOSTIC",
        "standalone",
        "diagnostic",
        ("D_BLOCK_DIAGNOSTIC",),
    ),
    ("S_R_DIRECT", "standalone", "confirmatory", ("R_DIRECT",)),
    ("S_R_CONVERGENT", "standalone", "confirmatory", ("R_CONVERGENT",)),
    ("S_R_HIGH_STRICT", "standalone", "exploratory", ("R_HIGH_STRICT",)),
    ("S_R_HIGH_ANY", "standalone", "exploratory", ("R_HIGH_ANY",)),
    ("S_R_CONTEXT_A", "standalone", "exploratory", ("R_CONTEXT_A",)),
    ("L0_E0", "cumulative", "baseline", ()),
    ("L1_WEAK_RAW", "cumulative", "confirmatory", ("Q_WEAK_RAW",)),
    (
        "L2_DIRECT_RESTORE",
        "cumulative",
        "confirmatory",
        ("Q_WEAK_RAW", "Q_DIRECT_RESTORE"),
    ),
    (
        "L3_CONTEXT_REVIEW",
        "cumulative",
        "confirmatory",
        ("Q_WEAK_RAW", "Q_DIRECT_RESTORE", "Q_CONTEXT_REVIEW"),
    ),
    (
        "L4_PLUS_GUARD_REVIEW",
        "cumulative",
        "confirmatory",
        ("Q_WEAK_RAW", "Q_DIRECT_RESTORE", "Q_CONTEXT_REVIEW", "Q_GUARD"),
    ),
    (
        PRIMARY_ARM,
        "cumulative",
        "primary",
        (
            "Q_WEAK_RAW",
            "Q_DIRECT_RESTORE",
            "Q_CONTEXT_REVIEW",
            "Q_GUARD",
            "R_DIRECT",
            "R_CONVERGENT",
        ),
    ),
    (
        "L6_PLUS_HIGH_STRICT",
        "cumulative",
        "exploratory",
        (
            "Q_WEAK_RAW",
            "Q_DIRECT_RESTORE",
            "Q_CONTEXT_REVIEW",
            "Q_GUARD",
            "R_DIRECT",
            "R_CONVERGENT",
            "R_HIGH_STRICT",
        ),
    ),
    (
        "L7_PLUS_HIGH_ANY",
        "cumulative",
        "exploratory",
        (
            "Q_WEAK_RAW",
            "Q_DIRECT_RESTORE",
            "Q_CONTEXT_REVIEW",
            "Q_GUARD",
            "R_DIRECT",
            "R_CONVERGENT",
            "R_HIGH_STRICT",
            "R_HIGH_ANY",
        ),
    ),
    (
        "L8_PLUS_CONTEXT_A",
        "cumulative",
        "exploratory",
        (
            "Q_WEAK_RAW",
            "Q_DIRECT_RESTORE",
            "Q_CONTEXT_REVIEW",
            "Q_GUARD",
            "R_DIRECT",
            "R_CONVERGENT",
            "R_HIGH_STRICT",
            "R_HIGH_ANY",
            "R_CONTEXT_A",
        ),
    ),
)


class PolicyContractError(ValueError):
    """Raised when configuration drifts from the frozen pruning experiment."""


@dataclass(frozen=True)
class PolicyConfig:
    path: str
    sha256: str
    raw: JSONDict

    @property
    def arms(self) -> list[JSONDict]:
        return self.raw["arms"]


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _validate_policy(raw: Mapping[str, Any]) -> None:
    if raw.get("schema_version") != CONFIG_SCHEMA_VERSION:
        raise PolicyContractError(
            f"Expected {CONFIG_SCHEMA_VERSION}, got {raw.get('schema_version')!r}"
        )
    if raw.get("policy_version") != "clinical_pruning_v1":
        raise PolicyContractError("policy_version must be clinical_pruning_v1")
    source_contract = raw.get("source_contract")
    if not isinstance(source_contract, Mapping):
        raise PolicyContractError("source_contract must be an object")
    if source_contract.get("schema_version") != "rag_re_clinical.output.v1":
        raise PolicyContractError("Policy must consume rag_re_clinical.output.v1")
    if source_contract.get("run_status") != "complete":
        raise PolicyContractError("Policy source run_status must be complete")
    if source_contract.get("generation_reads_gold") is not False:
        raise PolicyContractError("generation_reads_gold must be false")
    if source_contract.get("candidate_universe") != "locked_clinical_artifact_candidates":
        raise PolicyContractError("candidate_universe must remain locked")
    prediction_views = raw.get("prediction_views")
    if not isinstance(prediction_views, Mapping):
        raise PolicyContractError("prediction_views must be an object")
    if prediction_views.get("auto_positive_states") != ["KEEP", "RESCUE"]:
        raise PolicyContractError("auto_positive_states must be KEEP + RESCUE")
    if prediction_views.get("review_inclusive_states") != [
        "KEEP",
        "RESCUE",
        "MANUAL",
    ]:
        raise PolicyContractError(
            "review_inclusive_states must be KEEP + RESCUE + MANUAL"
        )
    if raw.get("primary_arm") != PRIMARY_ARM:
        raise PolicyContractError(f"primary_arm must be {PRIMARY_ARM}")

    rules = raw.get("rules")
    if not isinstance(rules, Mapping) or set(rules) != set(RULE_EVALUATORS):
        raise PolicyContractError("Policy rule IDs do not match the frozen implementation")

    arms = raw.get("arms")
    if not isinstance(arms, list):
        raise PolicyContractError("arms must be a list")
    actual: list[tuple[str, str, str, tuple[str, ...]]] = []
    for index, arm in enumerate(arms):
        if not isinstance(arm, Mapping):
            raise PolicyContractError(f"arms[{index}] must be an object")
        sequence = arm.get("rule_sequence")
        if not isinstance(sequence, list) or not all(
            isinstance(rule, str) for rule in sequence
        ):
            raise PolicyContractError(f"arms[{index}].rule_sequence must be strings")
        actual.append(
            (
                arm.get("name"),
                arm.get("family"),
                arm.get("analysis_role"),
                tuple(sequence),
            )
        )
        should_be_primary = arm.get("name") == PRIMARY_ARM
        if arm.get("eligible_primary") is not should_be_primary:
            raise PolicyContractError(
                f"arms[{index}].eligible_primary inconsistent with primary arm"
            )
    if tuple(actual) != EXPECTED_ARMS:
        raise PolicyContractError("Arm order or frozen rule sequences have changed")

    cumulative_rules = {
        rule
        for name, family, _role, sequence in EXPECTED_ARMS
        if family == "cumulative"
        for rule in sequence
    }
    if "D_BLOCK_DIAGNOSTIC" in cumulative_rules:
        raise PolicyContractError("D_BLOCK_DIAGNOSTIC cannot enter cumulative arms")


def load_policy_config(path: PathLike) -> PolicyConfig:
    source = Path(path).resolve()
    raw = load_json(source)
    _validate_policy(raw)
    return PolicyConfig(path=str(source), sha256=sha256_file(source), raw=copy.deepcopy(raw))


def _initial_ledger(baseline_selected: bool) -> JSONDict:
    state = initial_state(baseline_selected)
    return {
        "step": 0,
        "rule_id": "INITIAL_STATE",
        "matched": True,
        "from_state": None,
        "to_state": state.value,
        "action": "INITIALIZE",
        "reason_code": "E0_BASELINE_KEEP" if baseline_selected else "NON_E0_NO_RESCUE",
        "facts": {"baseline_selected": baseline_selected},
    }


def _run_candidate_arm(candidate: Mapping[str, Any], sequence: list[str]) -> JSONDict:
    baseline = candidate.get("baseline_selected")
    if not isinstance(baseline, bool):
        raise ValueError("Projected candidate baseline_selected must be boolean")
    state = initial_state(baseline)
    ledger = [_initial_ledger(baseline)]
    for step, rule_id in enumerate(sequence, start=1):
        vote = evaluate_rule(rule_id, candidate, state)
        state, event = apply_vote(state, vote, baseline_selected=baseline)
        event["step"] = step
        ledger.append(event)
    validate_state(state, baseline_selected=baseline)
    return {"state": state.value, "rule_ledger": ledger}


def _arm_summary(
    candidates: list[Mapping[str, Any]], arm_name: str, arm: Mapping[str, Any]
) -> JSONDict:
    buckets: dict[State, list[str]] = {state: [] for state in State}
    auto_positive: list[str] = []
    review_inclusive: list[str] = []
    for candidate in candidates:
        name = candidate["organism_name"]
        state = State(candidate["arm_states"][arm_name]["state"])
        buckets[state].append(name)
        if state in AUTO_POSITIVE_STATES:
            auto_positive.append(name)
        if state in REVIEW_INCLUSIVE_STATES:
            review_inclusive.append(name)

    state_counts = {state.value: len(buckets[state]) for state in State}
    summary: JSONDict = {
        "family": arm["family"],
        "analysis_role": arm["analysis_role"],
        "eligible_primary": arm["eligible_primary"],
        "rule_sequence": list(arm["rule_sequence"]),
        "state_counts": state_counts,
        "auto_positive_pathogens": auto_positive,
        "auto_positive_count": len(auto_positive),
        "review_inclusive_pathogens": review_inclusive,
        "review_inclusive_count": len(review_inclusive),
        "manual_review_pathogens": buckets[State.MANUAL],
        "manual_review_count": len(buckets[State.MANUAL]),
        "pruned_pathogens": buckets[State.PRUNE],
        "pruned_count": len(buckets[State.PRUNE]),
        "rescued_pathogens": buckets[State.RESCUE],
        "rescued_count": len(buckets[State.RESCUE]),
        "no_rescue_pathogens": buckets[State.NO_RESCUE],
        "no_rescue_count": len(buckets[State.NO_RESCUE]),
    }
    return summary


def _validate_projected_source(projected: Mapping[str, Any]) -> None:
    def require_exact_keys(
        value: Mapping[str, Any], expected: set[str], *, where: str
    ) -> None:
        actual = set(value)
        if actual != expected:
            missing = sorted(expected - actual)
            extra = sorted(str(key) for key in actual - expected)
            raise ValueError(f"{where} keys mismatch; missing={missing}, extra={extra}")

    def require_mapping(value: Any, *, where: str) -> Mapping[str, Any]:
        if not isinstance(value, Mapping):
            raise ValueError(f"{where} must be an object")
        return value

    def require_text(value: Any, *, where: str) -> str:
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{where} must be a non-empty string")
        return " ".join(value.strip().split())

    def require_bool_or_none(value: Any, *, where: str) -> bool | None:
        if value is True or value is False or value is None:
            return value
        raise ValueError(f"{where} must be true, false, or null")

    require_exact_keys(
        projected,
        {"patient_id", "source", "candidates", "counts"},
        where="projected source",
    )
    patient_id = projected.get("patient_id")
    if not isinstance(patient_id, str) or not patient_id.strip():
        raise ValueError("Projected source requires a patient_id")
    source = require_mapping(projected.get("source"), where="projected source.source")
    require_exact_keys(
        source,
        {
            "path",
            "sha256",
            "sha256_basis",
            "schema_version",
            "pipeline_version",
            "artifact_sha256_basis",
        },
        where="projected source.source",
    )
    if source.get("schema_version") != "rag_re_clinical.output.v1":
        raise ValueError("Projected source schema must be rag_re_clinical.output.v1")
    require_text(source.get("path"), where="projected source.source.path")
    source_hash = source.get("sha256")
    if not isinstance(source_hash, str) or re.fullmatch(r"[0-9a-fA-F]{64}", source_hash) is None:
        raise ValueError("projected source.source.sha256 must be a SHA-256 hex digest")
    if source.get("sha256_basis") not in {
        "source_file_bytes_v1",
        "canonical_projected_inputs_v1",
    }:
        raise ValueError("projected source.source.sha256_basis is unsupported")
    candidates = projected.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise ValueError("Projected source requires candidates")
    counts = projected.get("counts")
    if not isinstance(counts, Mapping):
        raise ValueError("Projected source requires counts")
    require_exact_keys(
        counts,
        {"candidate_count", "e0_count", "non_e0_count"},
        where="projected source.counts",
    )
    for count_name in ("candidate_count", "e0_count", "non_e0_count"):
        value = counts.get(count_name)
        if type(value) is not int or value < 0:
            raise ValueError(f"projected source.counts.{count_name} must be a non-negative integer")

    seen_names: set[str] = set()
    e0_count = 0
    for index, raw_candidate in enumerate(candidates):
        where = f"projected source.candidates[{index}]"
        row = require_mapping(raw_candidate, where=where)
        require_exact_keys(
            row,
            {
                "organism_name",
                "canonical_organism_name",
                "baseline_selected",
                "review_tier",
                "source_features",
            },
            where=where,
        )
        name = require_text(row.get("organism_name"), where=f"{where}.organism_name")
        require_text(
            row.get("canonical_organism_name"),
            where=f"{where}.canonical_organism_name",
        )
        name_key = name.casefold()
        if name_key in seen_names:
            raise ValueError(f"{where}.organism_name is duplicated")
        seen_names.add(name_key)
        baseline = row.get("baseline_selected")
        if type(baseline) is not bool:
            raise ValueError(f"{where}.baseline_selected must be boolean")
        e0_count += int(baseline)
        tier = row.get("review_tier")
        if tier not in ALLOWED_REVIEW_TIERS:
            raise ValueError(f"{where}.review_tier is unsupported: {tier!r}")

        features = require_mapping(row.get("source_features"), where=f"{where}.source_features")
        require_exact_keys(
            features,
            {"D", "B_STRICT", "C", "E", "F"},
            where=f"{where}.source_features",
        )
        d = require_mapping(features.get("D"), where=f"{where}.source_features.D")
        require_exact_keys(d, {"state"}, where=f"{where}.source_features.D")
        if d.get("state") not in ALLOWED_D_STATES:
            raise ValueError(f"{where}.source_features.D.state is unsupported")

        b = require_mapping(
            features.get("B_STRICT"), where=f"{where}.source_features.B_STRICT"
        )
        require_exact_keys(
            b,
            {"positive", "absolute_strength_axis", "relative_strength_axis"},
            where=f"{where}.source_features.B_STRICT",
        )
        for field in ("positive", "absolute_strength_axis", "relative_strength_axis"):
            require_bool_or_none(
                b.get(field), where=f"{where}.source_features.B_STRICT.{field}"
            )

        c = require_mapping(features.get("C"), where=f"{where}.source_features.C")
        require_exact_keys(c, {"grade"}, where=f"{where}.source_features.C")
        if c.get("grade") not in ALLOWED_C_GRADES:
            raise ValueError(f"{where}.source_features.C.grade is unsupported")

        e = require_mapping(features.get("E"), where=f"{where}.source_features.E")
        require_exact_keys(e, {"A_positive"}, where=f"{where}.source_features.E")
        require_bool_or_none(
            e.get("A_positive"), where=f"{where}.source_features.E.A_positive"
        )

        f = require_mapping(features.get("F"), where=f"{where}.source_features.F")
        require_exact_keys(f, {"state"}, where=f"{where}.source_features.F")
        if f.get("state") not in ALLOWED_F_STATES:
            raise ValueError(f"{where}.source_features.F.state is unsupported")

    if counts.get("candidate_count") != len(candidates):
        raise ValueError("Projected candidate_count mismatch")
    if counts.get("e0_count") != e0_count:
        raise ValueError("Projected e0_count mismatch")
    if counts.get("non_e0_count") != len(candidates) - e0_count:
        raise ValueError("Projected non_e0_count mismatch")


def _validate_output(artifact: Mapping[str, Any], policy: PolicyConfig) -> None:
    if artifact.get("schema_version") != OUTPUT_SCHEMA_VERSION:
        raise AssertionError("Output schema mismatch")
    if artifact.get("run_status") != "complete":
        raise AssertionError("Output must be complete")
    arm_names = [arm[0] for arm in EXPECTED_ARMS]
    arms = artifact.get("arms")
    candidates = artifact.get("candidates")
    if not isinstance(arms, Mapping) or list(arms) != arm_names:
        raise AssertionError("Output arms are missing or out of order")
    if not isinstance(candidates, list):
        raise AssertionError("Output candidates must be a list")

    for candidate in candidates:
        baseline = candidate["baseline_selected"]
        arm_states = candidate.get("arm_states")
        if not isinstance(arm_states, Mapping) or list(arm_states) != arm_names:
            raise AssertionError("Every candidate must have every arm state")
        for arm in policy.arms:
            result = arm_states[arm["name"]]
            state = State(result["state"])
            validate_state(state, baseline_selected=baseline)
            ledger = result.get("rule_ledger")
            if not isinstance(ledger, list) or len(ledger) != len(arm["rule_sequence"]) + 1:
                raise AssertionError("Rule ledger length does not match rule sequence")
            if [event["rule_id"] for event in ledger[1:]] != arm["rule_sequence"]:
                raise AssertionError("Rule ledger order does not match rule sequence")

    baseline_names = [
        row["organism_name"] for row in candidates if row["baseline_selected"]
    ]
    for baseline_arm in ("S0_E0", "L0_E0"):
        summary = arms[baseline_arm]
        if summary["auto_positive_pathogens"] != baseline_names:
            raise AssertionError(f"{baseline_arm} must exactly reproduce E0")
        if summary["review_inclusive_pathogens"] != baseline_names:
            raise AssertionError(f"{baseline_arm} review view must exactly reproduce E0")

    for arm_name, summary in arms.items():
        expected = _arm_summary(candidates, arm_name, summary)
        for key in (
            "state_counts",
            "auto_positive_pathogens",
            "auto_positive_count",
            "review_inclusive_pathogens",
            "review_inclusive_count",
            "manual_review_pathogens",
            "manual_review_count",
            "pruned_pathogens",
            "pruned_count",
            "rescued_pathogens",
            "rescued_count",
            "no_rescue_pathogens",
            "no_rescue_count",
        ):
            if summary.get(key) != expected[key]:
                raise AssertionError(f"{arm_name}.{key} is inconsistent")

    if artifact.get("primary_arm") != PRIMARY_ARM:
        raise AssertionError("Primary arm mismatch")
    generation_contract = artifact.get("generation_contract")
    if not isinstance(generation_contract, Mapping):
        raise AssertionError("Missing generation contract")
    if generation_contract.get("reads_gold") is not False:
        raise AssertionError("Generation contract must explicitly prohibit gold")


def run_patient(
    projected: Mapping[str, Any],
    policy: PolicyConfig,
    *,
    generated_at_utc: str | None = None,
) -> JSONDict:
    """Generate pruning states from the allowlisted locked-source projection.

    This function deliberately has no gold/answer/label parameter.  Recall
    safety belongs to the separate evaluator and cannot alter generation.
    """

    _validate_projected_source(projected)
    _validate_policy(policy.raw)
    generated_at = generated_at_utc or _utc_now()

    candidates: list[JSONDict] = []
    for source_candidate in projected["candidates"]:
        arm_states: JSONDict = {}
        for arm in policy.arms:
            arm_states[arm["name"]] = _run_candidate_arm(
                source_candidate, arm["rule_sequence"]
            )
        candidates.append(
            {
                "organism_name": source_candidate["organism_name"],
                "canonical_organism_name": source_candidate[
                    "canonical_organism_name"
                ],
                "baseline_selected": source_candidate["baseline_selected"],
                "review_tier": source_candidate.get("review_tier"),
                "input_features": copy.deepcopy(source_candidate["source_features"]),
                "arm_states": arm_states,
            }
        )

    arms: JSONDict = {}
    for arm in policy.arms:
        arms[arm["name"]] = _arm_summary(candidates, arm["name"], arm)

    artifact: JSONDict = {
        "schema_version": OUTPUT_SCHEMA_VERSION,
        "pipeline_version": PIPELINE_VERSION,
        "policy_version": policy.raw["policy_version"],
        "study_status": policy.raw.get("study_status"),
        "run_status": "complete",
        "generated_at_utc": generated_at,
        "patient_id": projected["patient_id"],
        "source": copy.deepcopy(projected["source"]),
        "config": {
            "path": policy.path,
            "sha256": policy.sha256,
            "schema_version": policy.raw["schema_version"],
            "policy_version": policy.raw["policy_version"],
        },
        "generation_contract": {
            "source_schema": "rag_re_clinical.output.v1",
            "candidate_universe": "locked_clinical_artifact_candidates",
            "reads_gold": False,
            "unknown_is_negative": False,
            "configured_recall_safety_used_for_generation": False,
        },
        "primary_arm": PRIMARY_ARM,
        "counts": copy.deepcopy(projected["counts"]),
        "candidates": candidates,
        "arms": arms,
    }
    decision_basis = {
        "schema_version": artifact["schema_version"],
        "policy_version": artifact["policy_version"],
        "patient_id": artifact["patient_id"],
        "source_sha256": artifact["source"].get("sha256"),
        "config_sha256": artifact["config"]["sha256"],
        "primary_arm": artifact["primary_arm"],
        "candidates": artifact["candidates"],
        "arms": artifact["arms"],
    }
    artifact["decision_sha256_basis"] = "canonical_json_decision_payload_v1"
    artifact["decision_sha256"] = _canonical_sha256(decision_basis)
    _validate_output(artifact, policy)
    return artifact


__all__ = [
    "BATCH_MANIFEST_SCHEMA_VERSION",
    "CONFIG_SCHEMA_VERSION",
    "EXPECTED_ARMS",
    "OUTPUT_SCHEMA_VERSION",
    "PIPELINE_VERSION",
    "PRIMARY_ARM",
    "PolicyConfig",
    "PolicyContractError",
    "load_policy_config",
    "run_patient",
]

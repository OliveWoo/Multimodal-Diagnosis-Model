"""Deterministic post-L5 A/B/C factorial decision engine."""

from __future__ import annotations

import copy
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .source_adapter import (
    CLINICAL_SCHEMA,
    PRUNING_MANIFEST_SCHEMA,
    PRUNING_SCHEMA,
    RAW_SCHEMA,
    SourceValidationError,
    canonical_json_sha256,
    revalidate_joined_envelope,
    sha256_file,
)
from .states import (
    ARM_ORDER,
    BASE_ARM,
    CUMULATIVE_PREDICATES,
    D_STATES,
    E0_STATES,
    FAMILIES,
    FUTURE_HOLDOUT_CANDIDATE_ARM,
    NON_E0_STATES,
    PREDICATES,
    SCOPES,
    cumulative_arm,
    predicate_values,
    standalone_arm,
    ungated_arm,
)

OUTPUT_SCHEMA = "rag_re_clinical_pruning_abc.output.v1"
CONFIG_SCHEMA = "rag_re_clinical_pruning_abc.config.v1"
PIPELINE_VERSION = "0.1.0"


def load_config(path: str | Path) -> tuple[dict[str, Any], dict[str, Any]]:
    source = Path(path)
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise SourceValidationError(f"Cannot read config {source}: {exc}") from exc
    validate_config(payload)
    return payload, {
        "path": str(source.resolve()),
        "sha256": sha256_file(source),
        "schema_version": CONFIG_SCHEMA,
        "policy_version": payload["policy_version"],
    }


def validate_config(config: Any) -> None:
    if not isinstance(config, dict):
        raise SourceValidationError("Config must be an object")
    if config.get("schema_version") != CONFIG_SCHEMA:
        raise SourceValidationError(f"Config schema must be {CONFIG_SCHEMA}")
    if config.get("families") != list(FAMILIES):
        raise SourceValidationError("Config families/order differs from frozen contract")
    if config.get("scopes") != list(SCOPES):
        raise SourceValidationError("Config scopes/order differs from frozen contract")
    if config.get("predicates") != list(PREDICATES):
        raise SourceValidationError("Config predicates/order differs from frozen contract")
    if config.get("cumulative_order") != list(CUMULATIVE_PREDICATES):
        raise SourceValidationError("Config cumulative order differs from frozen contract")
    if config.get("primary_arm") != BASE_ARM:
        raise SourceValidationError("Current development primary must remain frozen L5")
    if config.get("future_holdout_candidate_arm") != FUTURE_HOLDOUT_CANDIDATE_ARM:
        raise SourceValidationError("Future holdout candidate arm differs from frozen contract")
    if config.get("new_arm_primary_eligibility") is not False:
        raise SourceValidationError("All new arms must be ineligible as current primary")
    if config.get("include_ungated_non_e0_stress") is not True:
        raise SourceValidationError("Frozen experiment includes ungated stress arms")
    if config.get("generation_reads_gold") is not False:
        raise SourceValidationError("Generation must not read gold")


def _validate_projection(joined: Any) -> None:
    if not isinstance(joined, dict):
        raise SourceValidationError("Joined envelope must be an object")
    allowed_root = {"patient_id", "sources", "counts", "candidates"}
    if set(joined) != allowed_root:
        raise SourceValidationError("Joined envelope keys differ from the label-blind projection")
    if not isinstance(joined["patient_id"], str) or not joined["patient_id"].strip():
        raise SourceValidationError("Joined patient_id must be a nonempty string")
    if not isinstance(joined["sources"], dict) or set(joined["sources"]) != {"pruning", "clinical", "raw"}:
        raise SourceValidationError("Joined sources must contain pruning, clinical, and raw")
    source_contracts = {
        "pruning": {
            "path",
            "sha256",
            "schema_version",
            "decision_sha256",
            "batch_manifest_path",
            "batch_manifest_sha256",
        },
        "clinical": {"path", "sha256", "schema_version"},
        "raw": {
            "path",
            "sha256",
            "schema_version",
            "pipeline_fingerprint",
            "config_hash",
            "candidate_pool_sha256",
        },
    }
    expected_schemas = {
        "pruning": PRUNING_SCHEMA,
        "clinical": CLINICAL_SCHEMA,
        "raw": RAW_SCHEMA,
    }
    for layer, required in source_contracts.items():
        value = joined["sources"].get(layer)
        if not isinstance(value, dict) or set(value) != required:
            raise SourceValidationError(f"Joined {layer} source differs from strict allowlist")
        if value.get("schema_version") != expected_schemas[layer]:
            raise SourceValidationError(f"Joined {layer} source schema is invalid")
        for key in ("path", "sha256"):
            if not isinstance(value.get(key), str) or not value[key]:
                raise SourceValidationError(f"Joined {layer}.{key} is invalid")
    candidates = joined["candidates"]
    if not isinstance(candidates, list) or not candidates:
        raise SourceValidationError("Joined candidates must be a nonempty list")
    allowed_candidate = {
        "organism_name",
        "canonical_organism_name",
        "baseline_selected",
        "review_tier",
        "l5_state",
        "context",
        "signals",
        "signal_audit",
    }
    seen: set[str] = set()
    for index, candidate in enumerate(candidates):
        if not isinstance(candidate, dict) or set(candidate) != allowed_candidate:
            raise SourceValidationError(f"Candidate {index} differs from projected allowlist")
        name = candidate.get("canonical_organism_name")
        if not isinstance(name, str) or not name.strip() or name.casefold() in seen:
            raise SourceValidationError(f"Candidate {index} name is empty or duplicated")
        seen.add(name.casefold())
        organism_name = candidate.get("organism_name")
        if not isinstance(organism_name, str) or not organism_name.strip():
            raise SourceValidationError(f"Candidate {index} organism_name is invalid")
        baseline = candidate.get("baseline_selected")
        if type(baseline) is not bool:
            raise SourceValidationError(f"Candidate {index} baseline_selected must be boolean")
        state = candidate.get("l5_state")
        if state not in (E0_STATES if baseline else NON_E0_STATES):
            raise SourceValidationError(f"Candidate {index} L5 state is invalid")
        review_tier = candidate.get("review_tier")
        allowed_tiers = {None} if baseline else {"review_high_priority", "review_context_needed"}
        if review_tier not in allowed_tiers:
            raise SourceValidationError(f"Candidate {index} review tier is invalid for its baseline stratum")
        context = candidate.get("context")
        if not isinstance(context, dict) or set(context) != {"D_state", "F_state", "clinical_C_grade"}:
            raise SourceValidationError(f"Candidate {index} context is invalid")
        if context["D_state"] not in D_STATES:
            raise SourceValidationError(f"Candidate {index} D state is invalid")
        if context["F_state"] not in {"F_COHERENT", "F_MISMATCH", "F_UNKNOWN"}:
            raise SourceValidationError(f"Candidate {index} F state is invalid")
        if context["clinical_C_grade"] not in {"C0", "C1", "C2", "C3", "CNEG"}:
            raise SourceValidationError(f"Candidate {index} clinical C grade is invalid")
        signals = candidate.get("signals")
        if not isinstance(signals, dict) or set(signals) != set(FAMILIES):
            raise SourceValidationError(f"Candidate {index} signal families are invalid")
        for family in FAMILIES:
            values = signals[family]
            if not isinstance(values, dict) or set(values) != {"A", "B", "C"}:
                raise SourceValidationError(f"Candidate {index} {family} signals are invalid")
            for key, value in values.items():
                if value is not None and type(value) is not bool:
                    raise SourceValidationError(f"Candidate {index} {family}.{key} must be tri-state")
        expected_corr_c = (
            True
            if context["clinical_C_grade"] in {"C2", "C3"}
            else False
            if context["clinical_C_grade"] == "CNEG"
            else None
        )
        if signals["CORR"]["C"] is not expected_corr_c:
            raise SourceValidationError(f"Candidate {index} CORR.C conflicts with clinical C grade")
        audit = candidate.get("signal_audit")
        if not isinstance(audit, dict) or set(audit) != {"RAW_A", "RAW_B", "RAW_C"}:
            raise SourceValidationError(f"Candidate {index} raw signal audit is invalid")
        audit_contracts = {
            "RAW_A": {"positive", "status", "support_count", "judgeable_count", "support_ratio"},
            "RAW_B": {"positive", "status", "site_aligned", "trigger_count"},
            "RAW_C": {"positive", "status", "site_aligned", "support_record_count"},
        }
        for module, keys in audit_contracts.items():
            if not isinstance(audit[module], dict) or set(audit[module]) != keys:
                raise SourceValidationError(f"Candidate {index} {module} audit differs from allowlist")
        if signals["CORR"]["A"] is not signals["RAW"]["A"]:
            raise SourceValidationError(f"Candidate {index} CORR.A must equal frozen RAW.A")
        for key, module in (("A", "RAW_A"), ("B", "RAW_B"), ("C", "RAW_C")):
            if audit[module]["positive"] is not signals["RAW"][key]:
                raise SourceValidationError(f"Candidate {index} RAW.{key} conflicts with its audit")
    counts = joined["counts"]
    if not isinstance(counts, dict) or set(counts) != {"candidate_count", "e0_count", "non_e0_count"}:
        raise SourceValidationError("Joined counts are invalid")
    e0 = sum(candidate["baseline_selected"] for candidate in candidates)
    expected = {
        "candidate_count": len(candidates),
        "e0_count": e0,
        "non_e0_count": len(candidates) - e0,
    }
    if counts != expected:
        raise SourceValidationError("Joined counts do not reconcile")


def arm_definitions() -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {
        BASE_ARM: {
            "kind": "baseline",
            "family": None,
            "scope": "ALL",
            "predicate": None,
            "gate": "FROZEN_L5",
            "analysis_role": "primary_baseline",
            "eligible_primary": True,
            "future_holdout_candidate": False,
            "rule_sequence": [],
        }
    }
    for family in FAMILIES:
        for scope in SCOPES:
            for predicate in PREDICATES:
                name = standalone_arm(family, scope, predicate)
                result[name] = {
                    "kind": "standalone",
                    "family": family,
                    "scope": scope,
                    "predicate": predicate,
                    "gate": "D_AWARE",
                    "analysis_role": "exploratory_ablation",
                    "eligible_primary": False,
                    "future_holdout_candidate": False,
                    "rule_sequence": [predicate],
                }
    for family in FAMILIES:
        for predicate in PREDICATES:
            name = ungated_arm(family, predicate)
            result[name] = {
                "kind": "standalone",
                "family": family,
                "scope": "RSC",
                "predicate": predicate,
                "gate": "NONE_DIAGNOSTIC",
                "analysis_role": "unsafe_diagnostic_stress",
                "eligible_primary": False,
                "future_holdout_candidate": False,
                "rule_sequence": [predicate],
            }
    for family in FAMILIES:
        for index, predicate in enumerate(CUMULATIVE_PREDICATES, start=1):
            name = cumulative_arm(family, index, predicate)
            sequence = list(CUMULATIVE_PREDICATES[:index])
            result[name] = {
                "kind": "cumulative",
                "family": family,
                "scope": "ALL_THREE_SCOPES",
                "predicate": predicate,
                "gate": "D_AWARE",
                "analysis_role": "future_holdout_candidate" if name == FUTURE_HOLDOUT_CANDIDATE_ARM else "exploratory_ladder",
                "eligible_primary": False,
                "future_holdout_candidate": name == FUTURE_HOLDOUT_CANDIDATE_ARM,
                "rule_sequence": sequence,
            }
    if tuple(result) != ARM_ORDER:
        raise RuntimeError("Arm definition order drifted from frozen ARM_ORDER")
    return result


def _initial_state(candidate: dict[str, Any]) -> tuple[str, list[dict[str, Any]]]:
    state = candidate["l5_state"]
    return state, [
        {
            "step": 0,
            "rule_id": "FROZEN_L5",
            "scope": "BASE",
            "predicate": None,
            "predicate_value": None,
            "D_state": candidate["context"]["D_state"],
            "from_state": None,
            "to_state": state,
            "action": "INITIALIZE",
            "reason_code": "COPY_FROZEN_L5_STATE",
        }
    ]


def _scope_for_candidate(candidate: dict[str, Any]) -> str | None:
    if not candidate["baseline_selected"] and candidate["l5_state"] == "NO_RESCUE":
        return "RSC"
    if candidate["baseline_selected"] and candidate["l5_state"] == "PRUNE":
        return "P_RST"
    if candidate["baseline_selected"] and candidate["l5_state"] == "MANUAL":
        return "M_RST"
    return None


def _apply(
    candidate: dict[str, Any],
    state: str,
    *,
    family: str,
    predicate: str,
    requested_scope: str,
    ungated: bool,
    step: int,
) -> tuple[str, dict[str, Any]]:
    actual_scope = _scope_for_candidate(candidate)
    value = predicate_values(**{key.lower(): candidate["signals"][family][key] for key in ("A", "B", "C")})[predicate]
    d_state = candidate["context"]["D_state"]
    entry = {
        "step": step,
        "rule_id": f"{family}_{requested_scope}_{predicate}" + ("_UNGATED" if ungated else ""),
        "scope": requested_scope,
        "predicate": predicate,
        "predicate_value": value,
        "D_state": d_state,
        "from_state": state,
        "to_state": state,
        "action": "NO_CHANGE",
        "reason_code": "",
    }
    if actual_scope != requested_scope:
        entry["reason_code"] = "OUTSIDE_TRANSITION_SCOPE"
        return state, entry
    if state != candidate["l5_state"]:
        entry["reason_code"] = "ALREADY_TRANSITIONED_BY_EARLIER_STAGE"
        return state, entry
    if value is None:
        entry["action"] = "ABSTAIN"
        entry["reason_code"] = "PREDICATE_UNKNOWN"
        return state, entry
    if value is False:
        entry["reason_code"] = "PREDICATE_FALSE"
        return state, entry
    if ungated:
        if requested_scope != "RSC":
            raise RuntimeError("Ungated diagnostics are restricted to non-E0 rescue")
        new_state = "RESCUE"
        entry.update(to_state=new_state, action="RESCUE", reason_code="UNGATED_TRUE_DIAGNOSTIC")
        return new_state, entry
    if d_state == "D_PASS":
        new_state = "RESCUE" if requested_scope == "RSC" else "KEEP"
        action = "RESCUE" if requested_scope == "RSC" else "RESTORE_AUTO"
        entry.update(to_state=new_state, action=action, reason_code="D_PASS_AND_PREDICATE_TRUE")
        return new_state, entry
    if d_state in {"D_GUARDED", "D_UNKNOWN"}:
        new_state = "MANUAL"
        action = "FLAG_MANUAL" if state != "MANUAL" else "NO_CHANGE"
        entry.update(to_state=new_state, action=action, reason_code=f"{d_state}_TRUE_ROUTES_MANUAL")
        return new_state, entry
    entry["reason_code"] = "D_BLOCK_PREVENTS_TRANSITION"
    return state, entry


def _standalone_state(
    candidate: dict[str, Any], family: str, scope: str, predicate: str, *, ungated: bool = False
) -> dict[str, Any]:
    state, ledger = _initial_state(candidate)
    state, entry = _apply(
        candidate,
        state,
        family=family,
        predicate=predicate,
        requested_scope=scope,
        ungated=ungated,
        step=1,
    )
    ledger.append(entry)
    return {"state": state, "rule_ledger": ledger}


def _cumulative_state(candidate: dict[str, Any], family: str, sequence: list[str]) -> dict[str, Any]:
    state, ledger = _initial_state(candidate)
    scope = _scope_for_candidate(candidate)
    if scope is None:
        for step, predicate in enumerate(sequence, start=1):
            _, entry = _apply(
                candidate,
                state,
                family=family,
                predicate=predicate,
                requested_scope="RSC",
                ungated=False,
                step=step,
            )
            entry["scope"] = "ALL_THREE_SCOPES"
            entry["reason_code"] = "OUTSIDE_ALL_TRANSITION_SCOPES"
            ledger.append(entry)
        return {"state": state, "rule_ledger": ledger}
    for step, predicate in enumerate(sequence, start=1):
        state, entry = _apply(
            candidate,
            state,
            family=family,
            predicate=predicate,
            requested_scope=scope,
            ungated=False,
            step=step,
        )
        entry["scope"] = "ALL_THREE_SCOPES"
        ledger.append(entry)
    return {"state": state, "rule_ledger": ledger}


def _candidate_output(candidate: dict[str, Any], definitions: dict[str, dict[str, Any]]) -> dict[str, Any]:
    predicates = {
        family: predicate_values(**{key.lower(): candidate["signals"][family][key] for key in ("A", "B", "C")})
        for family in FAMILIES
    }
    arm_states: dict[str, dict[str, Any]] = {}
    initial, ledger = _initial_state(candidate)
    arm_states[BASE_ARM] = {"state": initial, "rule_ledger": ledger}
    for family in FAMILIES:
        for scope in SCOPES:
            for predicate in PREDICATES:
                arm_states[standalone_arm(family, scope, predicate)] = _standalone_state(
                    candidate, family, scope, predicate
                )
    for family in FAMILIES:
        for predicate in PREDICATES:
            arm_states[ungated_arm(family, predicate)] = _standalone_state(
                candidate, family, "RSC", predicate, ungated=True
            )
    for family in FAMILIES:
        for index, predicate in enumerate(CUMULATIVE_PREDICATES, start=1):
            sequence = list(CUMULATIVE_PREDICATES[:index])
            arm_states[cumulative_arm(family, index, predicate)] = _cumulative_state(
                candidate, family, sequence
            )
    if tuple(arm_states) != tuple(definitions):
        raise RuntimeError("Candidate arm order differs from definitions")
    return {
        "organism_name": candidate["organism_name"],
        "canonical_organism_name": candidate["canonical_organism_name"],
        "baseline_selected": candidate["baseline_selected"],
        "review_tier": candidate["review_tier"],
        "l5_state": candidate["l5_state"],
        "context": copy.deepcopy(candidate["context"]),
        "signals": copy.deepcopy(candidate["signals"]),
        "predicates": predicates,
        "signal_audit": copy.deepcopy(candidate["signal_audit"]),
        "arm_states": arm_states,
    }


def _arm_summary(
    name: str,
    definition: dict[str, Any],
    candidates: list[dict[str, Any]],
) -> dict[str, Any]:
    states = [candidate["arm_states"][name]["state"] for candidate in candidates]
    state_counts = dict(sorted(Counter(states).items()))

    def names_for(allowed: set[str]) -> list[str]:
        return [
            candidate["canonical_organism_name"]
            for candidate, state in zip(candidates, states)
            if state in allowed
        ]

    auto = names_for({"KEEP", "RESCUE"})
    review = names_for({"KEEP", "RESCUE", "MANUAL"})
    manual = names_for({"MANUAL"})
    pruned = names_for({"PRUNE"})
    rescued = names_for({"RESCUE"})
    no_rescue = names_for({"NO_RESCUE"})
    added_non_e0 = []
    restored_prune_auto = []
    restored_prune_review = []
    promoted_manual_auto = []
    changed = []
    for candidate, state in zip(candidates, states):
        original = candidate["l5_state"]
        canonical = candidate["canonical_organism_name"]
        if state != original:
            changed.append(canonical)
        if not candidate["baseline_selected"] and original == "NO_RESCUE" and state == "RESCUE":
            added_non_e0.append(canonical)
        if candidate["baseline_selected"] and original == "PRUNE" and state == "KEEP":
            restored_prune_auto.append(canonical)
        if candidate["baseline_selected"] and original == "PRUNE" and state in {"KEEP", "MANUAL"}:
            restored_prune_review.append(canonical)
        if candidate["baseline_selected"] and original == "MANUAL" and state == "KEEP":
            promoted_manual_auto.append(canonical)
    summary = copy.deepcopy(definition)
    summary.update(
        {
            "state_counts": state_counts,
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
            "transition_summary": {
                "changed_from_l5_pathogens": changed,
                "changed_from_l5_count": len(changed),
                "added_non_e0_auto_pathogens": added_non_e0,
                "added_non_e0_auto_count": len(added_non_e0),
                "restored_e0_prune_auto_pathogens": restored_prune_auto,
                "restored_e0_prune_auto_count": len(restored_prune_auto),
                "restored_e0_prune_review_pathogens": restored_prune_review,
                "restored_e0_prune_review_count": len(restored_prune_review),
                "promoted_e0_manual_auto_pathogens": promoted_manual_auto,
                "promoted_e0_manual_auto_count": len(promoted_manual_auto),
            },
        }
    )
    return summary


def _run_projected_patient(
    joined: dict[str, Any],
    config: dict[str, Any],
    *,
    config_meta: dict[str, Any] | None = None,
    generated_at_utc: str | None = None,
) -> dict[str, Any]:
    validate_config(config)
    _validate_projection(joined)
    definitions = arm_definitions()
    candidates = [_candidate_output(candidate, definitions) for candidate in joined["candidates"]]
    arms = {
        name: _arm_summary(name, definition, candidates)
        for name, definition in definitions.items()
    }
    if generated_at_utc is None:
        generated_at_utc = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    config_meta = copy.deepcopy(config_meta) if config_meta is not None else {
        "path": "<memory>",
        "sha256": canonical_json_sha256(config),
        "schema_version": CONFIG_SCHEMA,
        "policy_version": config["policy_version"],
    }
    result: dict[str, Any] = {
        "schema_version": OUTPUT_SCHEMA,
        "pipeline_version": PIPELINE_VERSION,
        "policy_version": config["policy_version"],
        "study_status": config["study_status"],
        "run_status": "complete",
        "generated_at_utc": generated_at_utc,
        "patient_id": joined["patient_id"],
        "sources": copy.deepcopy(joined["sources"]),
        "config": config_meta,
        "generation_contract": {
            "reads_gold": False,
            "candidate_join": "patient_scoped_exact_canonical_name_and_order",
            "unknown_is_negative": False,
            "frozen_l5_is_primary": True,
            "new_arms_eligible_primary": False,
            "raw_evidence_retrieval_reused": True,
            "new_api_or_llm_calls": False,
        },
        "primary_arm": BASE_ARM,
        "future_holdout_candidate_arm": FUTURE_HOLDOUT_CANDIDATE_ARM,
        "arm_order": list(ARM_ORDER),
        "counts": copy.deepcopy(joined["counts"]),
        "candidates": candidates,
        "arms": arms,
    }
    decision_basis = {
        "schema_version": result["schema_version"],
        "pipeline_version": result["pipeline_version"],
        "policy_version": result["policy_version"],
        "patient_id": result["patient_id"],
        "source_hashes": {key: value.get("sha256") for key, value in result["sources"].items()},
        "config_sha256": config_meta.get("sha256"),
        "arm_order": result["arm_order"],
        "candidates": candidates,
        "arms": arms,
    }
    result["decision_sha256_basis"] = "canonical_projected_inputs_config_and_decisions_v1"
    result["decision_sha256"] = canonical_json_sha256(decision_basis)
    return result


def _validate_config_meta(config: dict[str, Any], config_meta: Any) -> None:
    required = {"path", "sha256", "schema_version", "policy_version"}
    if not isinstance(config_meta, dict) or set(config_meta) != required:
        raise SourceValidationError("Public run_patient requires strict file-backed config metadata")
    if (
        config_meta.get("schema_version") != CONFIG_SCHEMA
        or config_meta.get("policy_version") != config.get("policy_version")
    ):
        raise SourceValidationError("Config metadata schema/policy mismatch")
    path = Path(config_meta.get("path", ""))
    if not path.is_file() or sha256_file(path) != config_meta.get("sha256"):
        raise SourceValidationError("Config metadata path/SHA is stale or invalid")
    try:
        file_config = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise SourceValidationError(f"Cannot revalidate config metadata: {exc}") from exc
    if canonical_json_sha256(file_config) != canonical_json_sha256(config):
        raise SourceValidationError("Supplied config differs from its anchored file")


def run_patient(
    joined: dict[str, Any],
    config: dict[str, Any],
    *,
    config_meta: dict[str, Any],
    generated_at_utc: str | None = None,
) -> dict[str, Any]:
    """Run only after re-binding the projection to locked files and manifest."""

    validate_config(config)
    _validate_config_meta(config, config_meta)
    _validate_projection(joined)
    revalidate_joined_envelope(joined)
    return _run_projected_patient(
        joined,
        config,
        config_meta=config_meta,
        generated_at_utc=generated_at_utc,
    )

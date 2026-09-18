"""Gold-isolated evaluation for the post-L5 corrected/raw A/B/C experiment."""

from __future__ import annotations

import hashlib
import json
import math
import random
import re
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .states import (
    ARM_ORDER,
    BASE_ARM,
    CUMULATIVE_PREDICATES,
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


OUTPUT_SCHEMA_VERSION = "rag_re_clinical_pruning_abc.output.v1"
MANIFEST_SCHEMA_VERSION = "rag_re_clinical_pruning_abc.batch_manifest.v1"
REPORT_SCHEMA_VERSION = "rag_re_clinical_pruning_abc.evaluation.v1"
PROTOCOL_SCHEMA_VERSION = "rag_re_clinical_pruning_abc.evaluation_protocol.v1"
EVALUATOR_VERSION = "0.1.0"
PRIMARY_SCOPE = "all_labeled_exact_frozen_synonyms"
SENSITIVITY_SCOPE = "answer_positive_only_genus_relaxed"
E0_REFERENCE = "E0_ORIGINAL"

FROZEN_SYNONYMS = {
    "cmv": "human cytomegalovirus",
    "hsv": "hsv-1",
    "herpes simplex virus 1": "hsv-1",
    "herpes simplex virus type 1": "hsv-1",
    "human alphaherpesvirus 1": "hsv-1",
    "human herpesvirus 1": "hsv-1",
    "pjp": "pneumocystis jirovecii",
    "pneumocystis jiroveci": "pneumocystis jirovecii",
    "candida albican": "candida albicans",
    "crkp": "klebsiella pneumoniae",
}

GENUS_STOPWORDS = {
    "human",
    "herpes",
    "herpesvirus",
    "influenza",
    "parainfluenza",
    "respiratory",
    "virus",
    "viral",
}

TRANSITION_LIST_FIELDS = (
    "changed_from_l5_pathogens",
    "added_non_e0_auto_pathogens",
    "restored_e0_prune_auto_pathogens",
    "restored_e0_prune_review_pathogens",
    "promoted_e0_manual_auto_pathogens",
)


class EvaluationError(ValueError):
    """Raised when evaluation inputs violate the frozen contract."""


@dataclass(frozen=True)
class CandidateView:
    name: str
    baseline_selected: bool
    l5_state: str
    signals: Mapping[str, Mapping[str, bool | None]]
    predicates: Mapping[str, Mapping[str, bool | None]]
    states: Mapping[str, str]


@dataclass(frozen=True)
class ArmView:
    metadata: Mapping[str, Any]
    auto: tuple[str, ...]
    review: tuple[str, ...]
    manual: tuple[str, ...]
    transitions: Mapping[str, tuple[str, ...]]


@dataclass(frozen=True)
class ArtifactView:
    patient_id: str
    path: str
    candidates: tuple[CandidateView, ...]
    arms: Mapping[str, ArmView]


def _arm_metadata() -> dict[str, dict[str, Any]]:
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
                result[standalone_arm(family, scope, predicate)] = {
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
            result[ungated_arm(family, predicate)] = {
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
            result[name] = {
                "kind": "cumulative",
                "family": family,
                "scope": "ALL_THREE_SCOPES",
                "predicate": predicate,
                "gate": "D_AWARE",
                "analysis_role": (
                    "future_holdout_candidate"
                    if name == FUTURE_HOLDOUT_CANDIDATE_ARM
                    else "exploratory_ladder"
                ),
                "eligible_primary": False,
                "future_holdout_candidate": name == FUTURE_HOLDOUT_CANDIDATE_ARM,
                "rule_sequence": list(CUMULATIVE_PREDICATES[:index]),
            }
    if tuple(result) != ARM_ORDER:
        raise RuntimeError("Evaluator arm metadata differs from ARM_ORDER")
    return result


ARM_METADATA = _arm_metadata()


def _normalize(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).strip().lower()
    return re.sub(r"\s+", " ", text)


def _require_name(value: Any, *, where: str) -> str:
    value = _normalize(value)
    if not value:
        raise EvaluationError(f"{where} must be a non-empty pathogen name")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_batch_manifest(
    payload: Any,
    *,
    path: Path,
    artifacts: Mapping[str, ArtifactView],
    hashes: Mapping[str, str],
) -> dict[str, Any]:
    if not isinstance(payload, dict) or payload.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        raise EvaluationError("Unsupported ABC batch manifest")
    if payload.get("run_status") != "complete":
        raise EvaluationError("ABC batch manifest must have run_status=complete")
    if payload.get("output_schema_version") != OUTPUT_SCHEMA_VERSION:
        raise EvaluationError("ABC batch manifest output schema mismatch")
    if payload.get("primary_arm") != BASE_ARM or payload.get("new_arm_primary_eligibility") is not False:
        raise EvaluationError("ABC batch manifest primary contract drifted")
    if payload.get("arm_order") != list(ARM_ORDER):
        raise EvaluationError("ABC batch manifest arm order drifted")
    raw_entries = payload.get("artifacts")
    if not isinstance(raw_entries, list) or payload.get("patient_count") != len(raw_entries):
        raise EvaluationError("ABC batch manifest artifact/patient counts are invalid")
    seen: set[str] = set()
    decision_hashes: list[str] = []
    candidate_total = e0_total = non_e0_total = 0
    for index, entry in enumerate(raw_entries):
        if not isinstance(entry, dict):
            raise EvaluationError(f"ABC manifest artifact[{index}] must be an object")
        patient_id = str(entry.get("patient_id") or "").strip()
        if not patient_id or patient_id in seen or patient_id not in artifacts:
            raise EvaluationError(f"ABC manifest patient {patient_id!r} is missing, duplicated, or stale")
        seen.add(patient_id)
        artifact_path = Path(artifacts[patient_id].path)
        recorded_path = entry.get("path")
        if not isinstance(recorded_path, str) or Path(recorded_path).resolve() != artifact_path.resolve():
            raise EvaluationError(f"ABC manifest path mismatch for patient {patient_id}")
        if entry.get("sha256") != hashes[patient_id] or _sha256(artifact_path) != hashes[patient_id]:
            raise EvaluationError(f"ABC manifest file SHA mismatch for patient {patient_id}")
        raw_artifact = _load_json(artifact_path)
        decision_sha = raw_artifact.get("decision_sha256")
        if not isinstance(decision_sha, str) or entry.get("decision_sha256") != decision_sha:
            raise EvaluationError(f"ABC manifest decision SHA mismatch for patient {patient_id}")
        counts = raw_artifact.get("counts")
        if not isinstance(counts, dict):
            raise EvaluationError(f"ABC artifact counts missing for patient {patient_id}")
        if entry.get("candidate_count") != counts.get("candidate_count") or entry.get("e0_count") != counts.get("e0_count"):
            raise EvaluationError(f"ABC manifest candidate counts mismatch for patient {patient_id}")
        decision_hashes.append(decision_sha)
        candidate_total += counts.get("candidate_count", 0)
        e0_total += counts.get("e0_count", 0)
        non_e0_total += counts.get("non_e0_count", 0)
    if seen != set(artifacts):
        raise EvaluationError("ABC manifest and artifact patient sets differ")
    expected_counts = {
        "candidate_count": candidate_total,
        "e0_count": e0_total,
        "non_e0_count": non_e0_total,
    }
    if payload.get("counts") != expected_counts:
        raise EvaluationError("ABC batch manifest aggregate counts do not reconcile")
    if payload.get("decision_set_sha256") != _canonical_sha256(decision_hashes):
        raise EvaluationError("ABC batch manifest decision set SHA mismatch")
    return {
        "path": str(path.resolve()),
        "sha256": _sha256(path),
        "validated": True,
        "patient_count": len(seen),
        "decision_set_sha256": payload.get("decision_set_sha256"),
    }


def _load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise EvaluationError(f"Could not read JSON from {path}: {exc}") from exc


def _build_aliases(gold_aliases: Any) -> dict[str, str]:
    aliases = {
        _require_name(key, where="frozen alias key"): _require_name(
            value, where="frozen alias value"
        )
        for key, value in FROZEN_SYNONYMS.items()
    }
    if gold_aliases is None:
        gold_aliases = {}
    if not isinstance(gold_aliases, dict):
        raise EvaluationError("gold.aliases must be an object")
    for raw_key, raw_value in gold_aliases.items():
        key = _require_name(raw_key, where="gold alias key")
        value = _require_name(raw_value, where=f"gold alias {raw_key!r}")
        if key in aliases and aliases[key] != value:
            raise EvaluationError(f"Gold alias {key!r} conflicts with frozen synonyms")
        aliases[key] = value
    for start in list(aliases):
        current = start
        visited: set[str] = set()
        while current in aliases:
            if current in visited:
                raise EvaluationError(f"Alias cycle detected at {start!r}")
            visited.add(current)
            current = aliases[current]
        aliases[start] = current
    return aliases


def _canonical(value: Any, aliases: Mapping[str, str]) -> str:
    current = _require_name(value, where="pathogen")
    visited: set[str] = set()
    while current in aliases:
        if current in visited:
            raise EvaluationError(f"Alias cycle for {value!r}")
        visited.add(current)
        current = aliases[current]
    return current


def _canonical_list(
    raw: Any,
    aliases: Mapping[str, str],
    *,
    where: str,
    reject_duplicates: bool = True,
) -> tuple[str, ...]:
    if not isinstance(raw, list):
        raise EvaluationError(f"{where} must be a list")
    values = tuple(_canonical(value, aliases) for value in raw)
    if reject_duplicates and len(set(values)) != len(values):
        raise EvaluationError(f"{where} contains duplicate or alias-equivalent names")
    return values if reject_duplicates else tuple(dict.fromkeys(values))


def _strict_tri(value: Any, *, where: str) -> bool | None:
    if value is None or type(value) is bool:
        return value
    raise EvaluationError(f"{where} must be true, false, or null")


def _genus(value: str) -> str | None:
    match = re.match(r"^([a-z][a-z-]+)(?:\s+)([a-z][a-z0-9.-]+)", value)
    if not match:
        return None
    token = match.group(1)
    return None if token in GENUS_STOPWORDS else token


def _names_match(left: str, right: str, *, genus_relaxed: bool) -> bool:
    if left == right:
        return True
    genus = _genus(left)
    return bool(genus_relaxed and genus and genus == _genus(right))


def _maximum_match_count(
    truth: Sequence[str], predictions: Sequence[str], *, genus_relaxed: bool
) -> int:
    edges = [
        [
            index
            for index, predicted in enumerate(predictions)
            if _names_match(gold, predicted, genus_relaxed=genus_relaxed)
        ]
        for gold in truth
    ]
    assigned: dict[int, int] = {}

    def augment(gold_index: int, visited: set[int]) -> bool:
        for prediction_index in edges[gold_index]:
            if prediction_index in visited:
                continue
            visited.add(prediction_index)
            if prediction_index not in assigned or augment(
                assigned[prediction_index], visited
            ):
                assigned[prediction_index] = gold_index
                return True
        return False

    return sum(augment(index, set()) for index in range(len(truth)))


def _safe_div(numerator: int | float, denominator: int | float) -> float | None:
    return numerator / denominator if denominator else None


def _fbeta(precision: float | None, recall: float | None, beta: float) -> float | None:
    if precision is None or recall is None:
        return None
    if precision == 0 and recall == 0:
        return 0.0
    beta2 = beta * beta
    denominator = beta2 * precision + recall
    return (1 + beta2) * precision * recall / denominator if denominator else 0.0


def _scope_for_candidate(candidate: CandidateView) -> str | None:
    if not candidate.baseline_selected and candidate.l5_state == "NO_RESCUE":
        return "RSC"
    if candidate.baseline_selected and candidate.l5_state == "PRUNE":
        return "P_RST"
    if candidate.baseline_selected and candidate.l5_state == "MANUAL":
        return "M_RST"
    return None


def _validate_count(raw: Mapping[str, Any], prefix: str, length: int, *, where: str) -> None:
    value = raw.get(f"{prefix}_count")
    if type(value) is not int or value != length:
        raise EvaluationError(f"{where}.{prefix}_count does not match its list")


def _expected_transition_sets(
    candidates: Sequence[CandidateView], arm: str
) -> dict[str, tuple[str, ...]]:
    changed: list[str] = []
    added: list[str] = []
    prune_auto: list[str] = []
    prune_review: list[str] = []
    manual_auto: list[str] = []
    for candidate in candidates:
        state = candidate.states[arm]
        original = candidate.l5_state
        if state != original:
            changed.append(candidate.name)
        if not candidate.baseline_selected and original == "NO_RESCUE" and state == "RESCUE":
            added.append(candidate.name)
        if candidate.baseline_selected and original == "PRUNE" and state == "KEEP":
            prune_auto.append(candidate.name)
        if candidate.baseline_selected and original == "PRUNE" and state in {"KEEP", "MANUAL"}:
            prune_review.append(candidate.name)
        if candidate.baseline_selected and original == "MANUAL" and state == "KEEP":
            manual_auto.append(candidate.name)
    return {
        "changed_from_l5_pathogens": tuple(changed),
        "added_non_e0_auto_pathogens": tuple(added),
        "restored_e0_prune_auto_pathogens": tuple(prune_auto),
        "restored_e0_prune_review_pathogens": tuple(prune_review),
        "promoted_e0_manual_auto_pathogens": tuple(manual_auto),
    }


def _validate_artifact(
    payload: Any,
    *,
    path: Path,
    aliases: Mapping[str, str],
) -> ArtifactView:
    if not isinstance(payload, dict):
        raise EvaluationError(f"ABC artifact {path} must be an object")
    if payload.get("schema_version") != OUTPUT_SCHEMA_VERSION:
        raise EvaluationError(
            f"ABC artifact {path} has unsupported schema_version {payload.get('schema_version')!r}"
        )
    patient_id = str(payload.get("patient_id") or "").strip()
    if not patient_id:
        raise EvaluationError(f"ABC artifact {path} is missing patient_id")
    if payload.get("run_status") != "complete":
        raise EvaluationError(f"Patient {patient_id} run_status must be complete")
    if payload.get("primary_arm") != BASE_ARM:
        raise EvaluationError(f"Patient {patient_id} primary_arm must be {BASE_ARM}")
    if payload.get("future_holdout_candidate_arm") != FUTURE_HOLDOUT_CANDIDATE_ARM:
        raise EvaluationError(f"Patient {patient_id} future holdout arm drifted")
    if payload.get("arm_order") != list(ARM_ORDER):
        raise EvaluationError(f"Patient {patient_id} arm_order differs from frozen contract")
    generation = payload.get("generation_contract")
    if not isinstance(generation, dict) or generation.get("reads_gold") is not False:
        raise EvaluationError(f"Patient {patient_id} generation contract must prohibit gold")

    raw_candidates = payload.get("candidates")
    if not isinstance(raw_candidates, list) or not raw_candidates:
        raise EvaluationError(f"Patient {patient_id} candidates must be nonempty")
    candidates: list[CandidateView] = []
    seen: set[str] = set()
    for index, raw_candidate in enumerate(raw_candidates):
        where = f"patient {patient_id} candidate[{index}]"
        if not isinstance(raw_candidate, dict):
            raise EvaluationError(f"{where} must be an object")
        name = _canonical(
            raw_candidate.get("canonical_organism_name")
            or raw_candidate.get("organism_name"),
            aliases,
        )
        if name in seen:
            raise EvaluationError(f"{where} duplicates alias-equivalent candidate {name!r}")
        seen.add(name)
        baseline = raw_candidate.get("baseline_selected")
        if type(baseline) is not bool:
            raise EvaluationError(f"{where}.baseline_selected must be boolean")
        l5_state = str(raw_candidate.get("l5_state") or "")
        if l5_state not in (E0_STATES if baseline else NON_E0_STATES):
            raise EvaluationError(f"{where}.l5_state is invalid")
        raw_signals = raw_candidate.get("signals")
        raw_predicates = raw_candidate.get("predicates")
        if not isinstance(raw_signals, dict) or set(raw_signals) != set(FAMILIES):
            raise EvaluationError(f"{where}.signals has invalid families")
        if not isinstance(raw_predicates, dict) or set(raw_predicates) != set(FAMILIES):
            raise EvaluationError(f"{where}.predicates has invalid families")
        signals: dict[str, dict[str, bool | None]] = {}
        predicates: dict[str, dict[str, bool | None]] = {}
        for family in FAMILIES:
            source = raw_signals.get(family)
            reported = raw_predicates.get(family)
            if not isinstance(source, dict) or set(source) != {"A", "B", "C"}:
                raise EvaluationError(f"{where}.signals.{family} is invalid")
            if not isinstance(reported, dict) or set(reported) != set(PREDICATES):
                raise EvaluationError(f"{where}.predicates.{family} is invalid")
            values = {
                key: _strict_tri(source[key], where=f"{where}.signals.{family}.{key}")
                for key in ("A", "B", "C")
            }
            expected = predicate_values(
                a=values["A"], b=values["B"], c=values["C"]
            )
            reported_values = {
                key: _strict_tri(
                    reported[key], where=f"{where}.predicates.{family}.{key}"
                )
                for key in PREDICATES
            }
            if reported_values != expected:
                raise EvaluationError(f"{where}.predicates.{family} disagrees with signals")
            signals[family] = values
            predicates[family] = reported_values
        raw_states = raw_candidate.get("arm_states")
        if not isinstance(raw_states, dict) or tuple(raw_states) != ARM_ORDER:
            raise EvaluationError(f"{where}.arm_states differs from ARM_ORDER")
        states: dict[str, str] = {}
        for arm in ARM_ORDER:
            raw_state = raw_states[arm]
            if not isinstance(raw_state, dict):
                raise EvaluationError(f"{where}.arm_states[{arm}] must be an object")
            state = str(raw_state.get("state") or "")
            if state not in (E0_STATES if baseline else NON_E0_STATES):
                raise EvaluationError(f"{where}.arm_states[{arm}].state is invalid")
            ledger = raw_state.get("rule_ledger")
            expected_length = 1 + len(ARM_METADATA[arm]["rule_sequence"])
            if not isinstance(ledger, list) or len(ledger) != expected_length:
                raise EvaluationError(f"{where}.arm_states[{arm}].rule_ledger length is invalid")
            if not ledger or not isinstance(ledger[0], dict) or ledger[0].get("rule_id") != "FROZEN_L5":
                raise EvaluationError(f"{where}.arm_states[{arm}] lacks FROZEN_L5 ledger anchor")
            if ledger[0].get("to_state") != l5_state:
                raise EvaluationError(f"{where}.arm_states[{arm}] ledger L5 state mismatch")
            if any(not isinstance(event, dict) for event in ledger) or ledger[-1].get("to_state") != state:
                raise EvaluationError(f"{where}.arm_states[{arm}] final ledger state mismatch")
            states[arm] = state
        if states[BASE_ARM] != l5_state:
            raise EvaluationError(f"{where} L5_FROZEN must equal l5_state")
        candidates.append(
            CandidateView(
                name=name,
                baseline_selected=baseline,
                l5_state=l5_state,
                signals=signals,
                predicates=predicates,
                states=states,
            )
        )

    counts = payload.get("counts")
    expected_counts = {
        "candidate_count": len(candidates),
        "e0_count": sum(candidate.baseline_selected for candidate in candidates),
        "non_e0_count": sum(not candidate.baseline_selected for candidate in candidates),
    }
    if counts != expected_counts:
        raise EvaluationError(f"Patient {patient_id} root counts do not reconcile")
    candidate_set = {candidate.name for candidate in candidates}
    raw_arms = payload.get("arms")
    if not isinstance(raw_arms, dict) or tuple(raw_arms) != ARM_ORDER:
        raise EvaluationError(f"Patient {patient_id} root arms differ from ARM_ORDER")
    arms: dict[str, ArmView] = {}
    for arm in ARM_ORDER:
        where = f"patient {patient_id} arm {arm}"
        raw_arm = raw_arms[arm]
        if not isinstance(raw_arm, dict):
            raise EvaluationError(f"{where} must be an object")
        expected_meta = ARM_METADATA[arm]
        for key, value in expected_meta.items():
            if raw_arm.get(key) != value:
                raise EvaluationError(f"{where}.{key} differs from frozen definition")
        list_specs = {
            "auto_positive": {"KEEP", "RESCUE"},
            "review_inclusive": {"KEEP", "RESCUE", "MANUAL"},
            "manual_review": {"MANUAL"},
            "pruned": {"PRUNE"},
            "rescued": {"RESCUE"},
            "no_rescue": {"NO_RESCUE"},
        }
        parsed: dict[str, tuple[str, ...]] = {}
        for prefix, allowed in list_specs.items():
            values = _canonical_list(
                raw_arm.get(f"{prefix}_pathogens"),
                aliases,
                where=f"{where}.{prefix}_pathogens",
            )
            _validate_count(raw_arm, prefix, len(values), where=where)
            if set(values) - candidate_set:
                raise EvaluationError(f"{where}.{prefix}_pathogens leaves candidate universe")
            expected = tuple(
                candidate.name
                for candidate in candidates
                if candidate.states[arm] in allowed
            )
            if values != expected:
                raise EvaluationError(f"{where}.{prefix}_pathogens disagrees with candidate states")
            parsed[prefix] = values
        state_counts = raw_arm.get("state_counts")
        expected_state_counts = dict(
            sorted(Counter(candidate.states[arm] for candidate in candidates).items())
        )
        if state_counts != expected_state_counts:
            raise EvaluationError(f"{where}.state_counts disagrees with candidates")
        raw_transition = raw_arm.get("transition_summary")
        if not isinstance(raw_transition, dict):
            raise EvaluationError(f"{where}.transition_summary must be an object")
        expected_transition = _expected_transition_sets(candidates, arm)
        transitions: dict[str, tuple[str, ...]] = {}
        for field in TRANSITION_LIST_FIELDS:
            values = _canonical_list(
                raw_transition.get(field),
                aliases,
                where=f"{where}.transition_summary.{field}",
            )
            prefix = field.removesuffix("_pathogens")
            _validate_count(raw_transition, prefix, len(values), where=f"{where}.transition_summary")
            if values != expected_transition[field]:
                raise EvaluationError(f"{where}.transition_summary.{field} is inconsistent")
            transitions[field] = values
        arms[arm] = ArmView(
            metadata=expected_meta,
            auto=parsed["auto_positive"],
            review=parsed["review_inclusive"],
            manual=parsed["manual_review"],
            transitions=transitions,
        )
    return ArtifactView(
        patient_id=patient_id,
        path=str(path.resolve()),
        candidates=tuple(candidates),
        arms=arms,
    )


def _load_gold(path: Path) -> tuple[dict[str, Any], dict[str, str]]:
    payload = _load_json(path)
    if not isinstance(payload, dict):
        raise EvaluationError("Gold must be an object")
    cases = payload.get("cases")
    if not isinstance(cases, list) or not cases:
        raise EvaluationError("gold.cases must be a non-empty list")
    aliases = _build_aliases(payload.get("aliases"))
    normalized: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise EvaluationError(f"gold.cases[{index}] must be an object")
        patient_id = str(case.get("patient_id") or "").strip()
        if not patient_id or patient_id in seen:
            raise EvaluationError(f"Gold patient_id missing or duplicated at index {index}")
        seen.add(patient_id)
        truth = _canonical_list(
            case.get("pathogens"),
            aliases,
            where=f"gold patient {patient_id}.pathogens",
            reject_duplicates=False,
        )
        uncertain = _canonical_list(
            case.get("uncertain_pathogens", []),
            aliases,
            where=f"gold patient {patient_id}.uncertain_pathogens",
            reject_duplicates=False,
        )
        if set(truth) & set(uncertain):
            raise EvaluationError(f"Gold patient {patient_id} positive/uncertain overlap")
        normalized.append(
            {"patient_id": patient_id, "truth": truth, "uncertain": uncertain}
        )
    return (
        {
            "schema_version": payload.get("schema_version"),
            "label_semantics": payload.get("label_semantics"),
            "complete_candidate_negative_labels": payload.get(
                "complete_candidate_negative_labels"
            ),
            "cases": normalized,
            "path": str(path.resolve()),
            "sha256": _sha256(path),
        },
        aliases,
    )


def _effect_added(
    truth: Sequence[str],
    base: Sequence[str],
    added: Sequence[str],
    *,
    genus_relaxed: bool,
) -> dict[str, Any]:
    base_values = tuple(dict.fromkeys(base))
    base_set = set(base_values)
    added_values = tuple(value for value in dict.fromkeys(added) if value not in base_set)
    before = _maximum_match_count(truth, base_values, genus_relaxed=genus_relaxed)
    after = _maximum_match_count(
        truth, base_values + added_values, genus_relaxed=genus_relaxed
    )
    tp = after - before
    return {
        "count": len(added_values),
        "tp": tp,
        "fp": len(added_values) - tp,
        "ppv": _safe_div(tp, len(added_values)),
        "pathogens": list(added_values),
    }


def _effect_removed(
    truth: Sequence[str],
    full: Sequence[str],
    removed: Sequence[str],
    *,
    genus_relaxed: bool,
) -> dict[str, Any]:
    full_values = tuple(dict.fromkeys(full))
    removed_set = set(removed)
    removed_values = tuple(value for value in full_values if value in removed_set)
    remaining = tuple(value for value in full_values if value not in removed_set)
    before = _maximum_match_count(truth, full_values, genus_relaxed=genus_relaxed)
    after = _maximum_match_count(truth, remaining, genus_relaxed=genus_relaxed)
    tp = before - after
    return {
        "count": len(removed_values),
        "tp": tp,
        "fp": len(removed_values) - tp,
        "ppv": _safe_div(tp, len(removed_values)),
        "pathogens": list(removed_values),
    }


def _prediction_counts(
    truth: Sequence[str], predictions: Sequence[str], *, genus_relaxed: bool
) -> dict[str, Any]:
    matched = _maximum_match_count(truth, predictions, genus_relaxed=genus_relaxed)
    return {
        "tp": matched,
        "fp": len(predictions) - matched,
        "fn": len(truth) - matched,
        "predicted_count": len(predictions),
        "predicted_pathogens": list(predictions),
    }


def _sum_effects(rows: Sequence[Mapping[str, Any]], key: str) -> dict[str, Any]:
    count = sum(row["attribution"][key]["count"] for row in rows)
    tp = sum(row["attribution"][key]["tp"] for row in rows)
    fp = sum(row["attribution"][key]["fp"] for row in rows)
    return {"count": count, "tp": tp, "fp": fp, "ppv": _safe_div(tp, count)}


def _aggregate_metrics(rows: Sequence[Mapping[str, Any]], kind: str) -> dict[str, Any]:
    tp = sum(row[kind]["tp"] for row in rows)
    fp = sum(row[kind]["fp"] for row in rows)
    fn = sum(row[kind]["fn"] for row in rows)
    precision = _safe_div(tp, tp + fp)
    recall = _safe_div(tp, tp + fn)
    no_path_total = sum(row["no_path_case"] for row in rows)
    no_path_correct = sum(
        row["no_path_case"]
        and row["artifact_present"]
        and row[kind]["predicted_count"] == 0
        for row in rows
    )
    no_path_positive = sum(
        row["no_path_case"]
        and row["artifact_present"]
        and row[kind]["predicted_count"] > 0
        for row in rows
    )
    no_path_missing = sum(
        row["no_path_case"] and not row["artifact_present"] for row in rows
    )
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "predicted_count": tp + fp,
        "truth_count": tp + fn,
        "precision": precision,
        "recall": recall,
        "f0_5": _fbeta(precision, recall, 0.5),
        "f1": _fbeta(precision, recall, 1.0),
        "no_path_specificity": _safe_div(no_path_correct, no_path_total),
        "no_path_cases": no_path_total,
        "no_path_correct_negative_cases": no_path_correct,
        "no_path_predicted_positive_cases": no_path_positive,
        "no_path_missing_artifact_cases": no_path_missing,
    }


def _metric_delta(current: Mapping[str, Any], reference: Mapping[str, Any]) -> dict[str, Any]:
    keys = (
        "tp",
        "fp",
        "fn",
        "predicted_count",
        "precision",
        "recall",
        "f0_5",
        "f1",
        "no_path_specificity",
    )
    return {
        key: (
            None
            if current.get(key) is None or reference.get(key) is None
            else current[key] - reference[key]
        )
        for key in keys
    }


def _coverage_for_artifact(artifact: ArtifactView, arm: str) -> dict[str, Any]:
    meta = ARM_METADATA[arm]
    family = meta["family"]
    if arm == BASE_ARM:
        return {
            "candidate_count": len(artifact.candidates),
            "eligible_count": 0,
            "predicates": {},
            "changed_count": 0,
        }
    eligible = [
        candidate
        for candidate in artifact.candidates
        if (
            _scope_for_candidate(candidate) in SCOPES
            if meta["scope"] == "ALL_THREE_SCOPES"
            else _scope_for_candidate(candidate) == meta["scope"]
        )
    ]
    coverage: dict[str, dict[str, int]] = {}
    for predicate in meta["rule_sequence"]:
        counts = Counter(
            "true"
            if candidate.predicates[family][predicate] is True
            else "false"
            if candidate.predicates[family][predicate] is False
            else "unknown"
            for candidate in eligible
        )
        coverage[predicate] = {
            "true": counts["true"],
            "false": counts["false"],
            "unknown": counts["unknown"],
        }
    return {
        "candidate_count": len(artifact.candidates),
        "eligible_count": len(eligible),
        "predicates": coverage,
        "changed_count": len(artifact.arms[arm].transitions["changed_from_l5_pathogens"]),
    }


def _sum_coverage(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    candidate_count = sum(row["predicate_coverage"]["candidate_count"] for row in rows)
    eligible_count = sum(row["predicate_coverage"]["eligible_count"] for row in rows)
    changed_count = sum(row["predicate_coverage"]["changed_count"] for row in rows)
    predicates: dict[str, dict[str, int]] = {}
    for row in rows:
        for predicate, counts in row["predicate_coverage"]["predicates"].items():
            target = predicates.setdefault(
                predicate, {"true": 0, "false": 0, "unknown": 0}
            )
            for key in target:
                target[key] += counts[key]
    return {
        "candidate_count": candidate_count,
        "eligible_count": eligible_count,
        "predicates": predicates,
        "changed_count": changed_count,
    }


def _evaluate_scope(
    *,
    cases: Sequence[Mapping[str, Any]],
    artifacts: Mapping[str, ArtifactView],
    positive_cases_only: bool,
    genus_relaxed: bool,
) -> dict[str, Any]:
    scoped_cases = [case for case in cases if not positive_cases_only or case["truth"]]
    arm_rows: dict[str, list[dict[str, Any]]] = {arm: [] for arm in ARM_ORDER}
    e0_rows: list[dict[str, Any]] = []
    for case in scoped_cases:
        patient_id = case["patient_id"]
        truth = tuple(case["truth"])
        uncertain = set(case["uncertain"])
        artifact = artifacts.get(patient_id)

        def filtered(values: Iterable[str]) -> tuple[str, ...]:
            return tuple(value for value in values if value not in uncertain)

        if artifact is None:
            e0_predictions: tuple[str, ...] = ()
        else:
            e0_predictions = filtered(
                candidate.name
                for candidate in artifact.candidates
                if candidate.baseline_selected
            )
        e0_counts = _prediction_counts(
            truth, e0_predictions, genus_relaxed=genus_relaxed
        )
        e0_rows.append(
            {
                "patient_id": patient_id,
                "artifact_present": artifact is not None,
                "no_path_case": not truth,
                "auto": e0_counts,
                "review_inclusive": dict(e0_counts),
            }
        )
        for arm in ARM_ORDER:
            if artifact is None:
                auto = review = l5_auto = l5_review = ()
                transitions = {field: () for field in TRANSITION_LIST_FIELDS}
                coverage = {
                    "candidate_count": 0,
                    "eligible_count": 0,
                    "predicates": {},
                    "changed_count": 0,
                }
                manual_count = 0
            else:
                view = artifact.arms[arm]
                auto = filtered(view.auto)
                review = filtered(view.review)
                l5_auto = filtered(artifact.arms[BASE_ARM].auto)
                l5_review = filtered(artifact.arms[BASE_ARM].review)
                transitions = {
                    field: filtered(view.transitions[field])
                    for field in TRANSITION_LIST_FIELDS
                }
                coverage = _coverage_for_artifact(artifact, arm)
                manual_count = len(view.manual)
            auto_counts = _prediction_counts(truth, auto, genus_relaxed=genus_relaxed)
            review_counts = _prediction_counts(
                truth, review, genus_relaxed=genus_relaxed
            )
            added_vs_l5_auto = tuple(value for value in auto if value not in set(l5_auto))
            removed_vs_l5_auto = tuple(value for value in l5_auto if value not in set(auto))
            added_vs_l5_review = tuple(
                value for value in review if value not in set(l5_review)
            )
            removed_vs_l5_review = tuple(
                value for value in l5_review if value not in set(review)
            )
            attribution = {
                "added_vs_l5_auto": _effect_added(
                    truth, l5_auto, added_vs_l5_auto, genus_relaxed=genus_relaxed
                ),
                "removed_vs_l5_auto": _effect_removed(
                    truth, l5_auto, removed_vs_l5_auto, genus_relaxed=genus_relaxed
                ),
                "added_vs_l5_review": _effect_added(
                    truth,
                    l5_review,
                    added_vs_l5_review,
                    genus_relaxed=genus_relaxed,
                ),
                "removed_vs_l5_review": _effect_removed(
                    truth,
                    l5_review,
                    removed_vs_l5_review,
                    genus_relaxed=genus_relaxed,
                ),
                "added_non_e0_auto": _effect_added(
                    truth,
                    tuple(
                        value
                        for value in auto
                        if value
                        not in set(transitions["added_non_e0_auto_pathogens"])
                    ),
                    transitions["added_non_e0_auto_pathogens"],
                    genus_relaxed=genus_relaxed,
                ),
                "restored_e0_prune_auto": _effect_added(
                    truth,
                    tuple(
                        value
                        for value in auto
                        if value
                        not in set(transitions["restored_e0_prune_auto_pathogens"])
                    ),
                    transitions["restored_e0_prune_auto_pathogens"],
                    genus_relaxed=genus_relaxed,
                ),
                "restored_e0_prune_review": _effect_added(
                    truth,
                    tuple(
                        value
                        for value in review
                        if value
                        not in set(transitions["restored_e0_prune_review_pathogens"])
                    ),
                    transitions["restored_e0_prune_review_pathogens"],
                    genus_relaxed=genus_relaxed,
                ),
                "promoted_e0_manual_auto": _effect_added(
                    truth,
                    tuple(
                        value
                        for value in auto
                        if value
                        not in set(transitions["promoted_e0_manual_auto_pathogens"])
                    ),
                    transitions["promoted_e0_manual_auto_pathogens"],
                    genus_relaxed=genus_relaxed,
                ),
            }
            arm_rows[arm].append(
                {
                    "patient_id": patient_id,
                    "artifact_present": artifact is not None,
                    "truth_count": len(truth),
                    "truth_pathogens": list(truth),
                    "no_path_case": not truth,
                    "auto": auto_counts,
                    "review_inclusive": review_counts,
                    "manual_review_count": manual_count,
                    "transitions": {key: list(value) for key, value in transitions.items()},
                    "attribution": attribution,
                    "predicate_coverage": coverage,
                }
            )
    e0_metrics = {
        "auto": _aggregate_metrics(e0_rows, "auto"),
        "review_inclusive": _aggregate_metrics(e0_rows, "review_inclusive"),
        "cases": e0_rows,
    }
    arms: dict[str, Any] = {}
    attribution_keys = (
        "added_vs_l5_auto",
        "removed_vs_l5_auto",
        "added_vs_l5_review",
        "removed_vs_l5_review",
        "added_non_e0_auto",
        "restored_e0_prune_auto",
        "restored_e0_prune_review",
        "promoted_e0_manual_auto",
    )
    for arm, rows in arm_rows.items():
        arms[arm] = {
            **ARM_METADATA[arm],
            "auto": _aggregate_metrics(rows, "auto"),
            "review_inclusive": _aggregate_metrics(rows, "review_inclusive"),
            "attribution": {key: _sum_effects(rows, key) for key in attribution_keys},
            "predicate_coverage": _sum_coverage(rows),
            "case_counts": {
                "included": len(rows),
                "artifact_present": sum(row["artifact_present"] for row in rows),
                "artifact_missing": sum(not row["artifact_present"] for row in rows),
                "answer_positive": sum(not row["no_path_case"] for row in rows),
                "no_path": sum(row["no_path_case"] for row in rows),
                "with_auto_positive": sum(row["auto"]["predicted_count"] > 0 for row in rows),
                "with_manual_review": sum(row["manual_review_count"] > 0 for row in rows),
            },
            "cases": rows,
        }
    l5 = arms[BASE_ARM]
    for arm in ARM_ORDER:
        arms[arm]["delta_vs_l5"] = {
            "reference_arm": BASE_ARM,
            "auto": _metric_delta(arms[arm]["auto"], l5["auto"]),
            "review_inclusive": _metric_delta(
                arms[arm]["review_inclusive"], l5["review_inclusive"]
            ),
        }
        arms[arm]["delta_vs_e0"] = {
            "reference_arm": E0_REFERENCE,
            "auto": _metric_delta(arms[arm]["auto"], e0_metrics["auto"]),
            "review_inclusive": _metric_delta(
                arms[arm]["review_inclusive"], e0_metrics["review_inclusive"]
            ),
        }
        arms[arm]["cumulative_marginal"] = None
    for family in FAMILIES:
        previous = BASE_ARM
        for index, predicate in enumerate(CUMULATIVE_PREDICATES, start=1):
            arm = cumulative_arm(family, index, predicate)
            current_rows = arms[arm]["cases"]
            previous_rows = arms[previous]["cases"]
            stage_rows: list[dict[str, Any]] = []
            for current, prior in zip(current_rows, previous_rows):
                truth = current["truth_pathogens"]
                current_auto = tuple(current["auto"]["predicted_pathogens"])
                prior_auto = tuple(prior["auto"]["predicted_pathogens"])
                current_review = tuple(
                    current["review_inclusive"]["predicted_pathogens"]
                )
                prior_review = tuple(prior["review_inclusive"]["predicted_pathogens"])
                added_auto = tuple(value for value in current_auto if value not in set(prior_auto))
                removed_auto = tuple(value for value in prior_auto if value not in set(current_auto))
                added_review = tuple(
                    value for value in current_review if value not in set(prior_review)
                )
                stage_rows.append(
                    {
                        "attribution": {
                            "added_auto": _effect_added(
                                truth,
                                prior_auto,
                                added_auto,
                                genus_relaxed=genus_relaxed,
                            ),
                            "removed_auto": _effect_removed(
                                truth,
                                prior_auto,
                                removed_auto,
                                genus_relaxed=genus_relaxed,
                            ),
                            "added_review": _effect_added(
                                truth,
                                prior_review,
                                added_review,
                                genus_relaxed=genus_relaxed,
                            ),
                        }
                    }
                )
            arms[arm]["cumulative_marginal"] = {
                "previous_arm": previous,
                "auto_metric_delta": _metric_delta(
                    arms[arm]["auto"], arms[previous]["auto"]
                ),
                "review_metric_delta": _metric_delta(
                    arms[arm]["review_inclusive"],
                    arms[previous]["review_inclusive"],
                ),
                "added_auto": _sum_effects(stage_rows, "added_auto"),
                "removed_auto": _sum_effects(stage_rows, "removed_auto"),
                "added_review": _sum_effects(stage_rows, "added_review"),
            }
            previous = arm
    return {
        "case_count": len(scoped_cases),
        "positive_cases_only": positive_cases_only,
        "genus_relaxed": genus_relaxed,
        "references": {E0_REFERENCE: e0_metrics},
        "arms": arms,
    }


def _arm_groups() -> dict[str, Any]:
    grouped: dict[str, Any] = {
        "baseline": [BASE_ARM],
        "standalone": {},
        "ungated_diagnostic": {},
        "cumulative": {},
    }
    for family in FAMILIES:
        grouped["standalone"][family] = {
            scope: [standalone_arm(family, scope, predicate) for predicate in PREDICATES]
            for scope in SCOPES
        }
        grouped["ungated_diagnostic"][family] = [
            ungated_arm(family, predicate) for predicate in PREDICATES
        ]
        grouped["cumulative"][family] = [
            cumulative_arm(family, index, predicate)
            for index, predicate in enumerate(CUMULATIVE_PREDICATES, start=1)
        ]
    return grouped


def _pareto_points(arms: Mapping[str, Any], arm_names: Sequence[str]) -> list[dict[str, Any]]:
    signatures: dict[tuple[int, int, int], list[str]] = defaultdict(list)
    for arm in arm_names:
        metrics = arms[arm]["auto"]
        signatures[(metrics["tp"], metrics["fp"], metrics["fn"])].append(arm)
    points: list[dict[str, Any]] = []
    for signature, equivalents in signatures.items():
        metrics = arms[equivalents[0]]["auto"]
        dominated = False
        for other_signature, other_equivalents in signatures.items():
            if other_signature == signature:
                continue
            other = arms[other_equivalents[0]]["auto"]
            if (
                other["precision"] is not None
                and other["recall"] is not None
                and metrics["precision"] is not None
                and metrics["recall"] is not None
                and other["precision"] >= metrics["precision"]
                and other["recall"] >= metrics["recall"]
                and (
                    other["precision"] > metrics["precision"]
                    or other["recall"] > metrics["recall"]
                )
            ):
                dominated = True
                break
        if not dominated:
            points.append(
                {
                    "representative_arm": equivalents[0],
                    "equivalent_arms": equivalents,
                    "tp": signature[0],
                    "fp": signature[1],
                    "fn": signature[2],
                    "precision": metrics["precision"],
                    "recall": metrics["recall"],
                    "f0_5": metrics["f0_5"],
                }
            )
    return sorted(
        points,
        key=lambda row: (
            -(row["precision"] if row["precision"] is not None else -1),
            -(row["recall"] if row["recall"] is not None else -1),
        ),
    )


def _pareto_report(scope: Mapping[str, Any]) -> dict[str, Any]:
    arms = scope["arms"]
    d_aware = [
        arm for arm in ARM_ORDER if ARM_METADATA[arm]["gate"] != "NONE_DIAGNOSTIC"
    ]
    groups: dict[str, list[str]] = {
        "all_arms": list(ARM_ORDER),
        "d_aware_only": d_aware,
    }
    for family in FAMILIES:
        groups[f"family_{family}"] = [
            arm for arm in ARM_ORDER if ARM_METADATA[arm]["family"] == family
        ]
        for transition_scope in SCOPES:
            groups[f"{family}_{transition_scope}_d_aware"] = [
                arm
                for arm in ARM_ORDER
                if ARM_METADATA[arm]["family"] == family
                and ARM_METADATA[arm]["scope"] == transition_scope
                and ARM_METADATA[arm]["gate"] == "D_AWARE"
            ]
    return {
        name: _pareto_points(arms, arm_names)
        for name, arm_names in groups.items()
        if arm_names
    }


def _percentile(values: Sequence[float], percentile: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * percentile
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def _bootstrap_scope(
    scope: Mapping[str, Any],
    *,
    targets: Sequence[str],
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    if replicates <= 0:
        return {"replicates": 0, "seed": seed, "arms": {}}
    base_rows = scope["arms"][BASE_ARM]["cases"]
    n = len(base_rows)
    rng = random.Random(seed)
    weights: list[list[int]] = []
    for _ in range(replicates):
        row = [0] * n
        for _ in range(n):
            row[rng.randrange(n)] += 1
        weights.append(row)
    output: dict[str, Any] = {}
    base_counts = [
        (row["auto"]["tp"], row["auto"]["fp"], row["auto"]["fn"])
        for row in base_rows
    ]
    for arm in targets:
        arm_rows = scope["arms"][arm]["cases"]
        arm_counts = [
            (row["auto"]["tp"], row["auto"]["fp"], row["auto"]["fn"])
            for row in arm_rows
        ]
        precision_deltas: list[float] = []
        recall_deltas: list[float] = []
        for weight in weights:
            btp = sum(w * value[0] for w, value in zip(weight, base_counts))
            bfp = sum(w * value[1] for w, value in zip(weight, base_counts))
            bfn = sum(w * value[2] for w, value in zip(weight, base_counts))
            atp = sum(w * value[0] for w, value in zip(weight, arm_counts))
            afp = sum(w * value[1] for w, value in zip(weight, arm_counts))
            afn = sum(w * value[2] for w, value in zip(weight, arm_counts))
            bp = _safe_div(btp, btp + bfp)
            br = _safe_div(btp, btp + bfn)
            ap = _safe_div(atp, atp + afp)
            ar = _safe_div(atp, atp + afn)
            if bp is not None and ap is not None:
                precision_deltas.append(ap - bp)
            if br is not None and ar is not None:
                recall_deltas.append(ar - br)
        point = scope["arms"][arm]["delta_vs_l5"]["auto"]
        output[arm] = {
            "point_precision_delta": point["precision"],
            "precision_delta_percentile_95_ci": [
                _percentile(precision_deltas, 0.025),
                _percentile(precision_deltas, 0.975),
            ],
            "point_recall_delta": point["recall"],
            "recall_delta_percentile_95_ci": [
                _percentile(recall_deltas, 0.025),
                _percentile(recall_deltas, 0.975),
            ],
            "probability_precision_delta_gt_zero": _safe_div(
                sum(value > 0 for value in precision_deltas), len(precision_deltas)
            ),
            "probability_recall_delta_ge_zero": _safe_div(
                sum(value >= 0 for value in recall_deltas), len(recall_deltas)
            ),
        }
    return {"replicates": replicates, "seed": seed, "arms": output}


def _cluster_sensitivity(scope: Mapping[str, Any], minimum_tp: int) -> dict[str, Any]:
    result: dict[str, Any] = {}
    l5_metrics = scope["arms"][BASE_ARM]["auto"]
    for arm in ARM_ORDER:
        arm_data = scope["arms"][arm]
        added_tp = arm_data["attribution"]["added_non_e0_auto"]["tp"]
        if added_tp < minimum_tp:
            continue
        clusters: dict[str, dict[str, int]] = defaultdict(
            lambda: {"candidate_count": 0, "essential_tp_count": 0}
        )
        for row in arm_data["cases"]:
            full = tuple(row["auto"]["predicted_pathogens"])
            truth = row["truth_pathogens"]
            transition = row["transitions"]["added_non_e0_auto_pathogens"]
            full_match = _maximum_match_count(
                truth, full, genus_relaxed=scope["genus_relaxed"]
            )
            for organism in transition:
                clusters[organism]["candidate_count"] += 1
                without = tuple(value for value in full if value != organism)
                lost = full_match - _maximum_match_count(
                    truth, without, genus_relaxed=scope["genus_relaxed"]
                )
                clusters[organism]["essential_tp_count"] += max(0, lost)
        detail: list[dict[str, Any]] = []
        for organism, counts in sorted(clusters.items()):
            tp = fp = fn = 0
            removed_count = 0
            for row in arm_data["cases"]:
                transition = set(row["transitions"]["added_non_e0_auto_pathogens"])
                predictions = tuple(
                    value
                    for value in row["auto"]["predicted_pathogens"]
                    if not (value == organism and value in transition)
                )
                removed_count += row["auto"]["predicted_count"] - len(predictions)
                matched = _maximum_match_count(
                    row["truth_pathogens"],
                    predictions,
                    genus_relaxed=scope["genus_relaxed"],
                )
                tp += matched
                fp += len(predictions) - matched
                fn += len(row["truth_pathogens"]) - matched
            precision = _safe_div(tp, tp + fp)
            recall = _safe_div(tp, tp + fn)
            metrics = {
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": precision,
                "recall": recall,
                "f0_5": _fbeta(precision, recall, 0.5),
                "f1": _fbeta(precision, recall, 1.0),
            }
            detail.append(
                {
                    "organism": organism,
                    **counts,
                    "removed_prediction_count": removed_count,
                    "lost_tp": arm_data["auto"]["tp"] - tp,
                    "removed_fp": arm_data["auto"]["fp"] - fp,
                    "metrics_after_cluster_removal": metrics,
                    "delta_vs_l5_after_cluster_removal": _metric_delta(
                        metrics, l5_metrics
                    ),
                }
            )
        max_cluster_tp = max(
            (row["essential_tp_count"] for row in detail), default=0
        )
        result[arm] = {
            "added_non_e0_tp": added_tp,
            "added_non_e0_fp": arm_data["attribution"]["added_non_e0_auto"]["fp"],
            "distinct_added_organism_count": len(clusters),
            "maximum_single_organism_tp_share": _safe_div(max_cluster_tp, added_tp),
            "clusters": detail,
        }
    return result


def _load_protocol(
    path: str | Path | None,
    *,
    bootstrap_replicates: int | None,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    settings = {
        "bootstrap": {
            "replicates": 0,
            "seed": 20260811,
            "target_policy": "pareto_and_named",
            "named_arms": [BASE_ARM, FUTURE_HOLDOUT_CANDIDATE_ARM],
        },
        "organism_cluster": {"minimum_added_tp_for_detail": 1},
    }
    meta = None
    if path is not None:
        protocol_path = Path(path)
        payload = _load_json(protocol_path)
        if not isinstance(payload, dict) or payload.get("schema_version") != PROTOCOL_SCHEMA_VERSION:
            raise EvaluationError("Unsupported ABC evaluation protocol schema")
        if payload.get("primary_anchor") != BASE_ARM:
            raise EvaluationError("Protocol primary_anchor differs from frozen L5")
        raw_bootstrap = payload.get("bootstrap")
        raw_cluster = payload.get("organism_cluster")
        if not isinstance(raw_bootstrap, dict) or not isinstance(raw_cluster, dict):
            raise EvaluationError("Protocol bootstrap/organism_cluster must be objects")
        settings = {
            "bootstrap": dict(raw_bootstrap),
            "organism_cluster": dict(raw_cluster),
        }
        meta = {"path": str(protocol_path.resolve()), "sha256": _sha256(protocol_path)}
    if bootstrap_replicates is not None:
        settings["bootstrap"]["replicates"] = bootstrap_replicates
    replicates = settings["bootstrap"].get("replicates")
    seed = settings["bootstrap"].get("seed")
    if type(replicates) is not int or replicates < 0 or type(seed) is not int:
        raise EvaluationError("Bootstrap replicates/seed must be integers")
    named = settings["bootstrap"].get("named_arms", [])
    if not isinstance(named, list) or any(arm not in ARM_ORDER for arm in named):
        raise EvaluationError("Protocol bootstrap named_arms contains unsupported arms")
    minimum = settings["organism_cluster"].get("minimum_added_tp_for_detail", 1)
    if type(minimum) is not int or minimum < 1:
        raise EvaluationError("minimum_added_tp_for_detail must be a positive integer")
    return settings, meta


def evaluate_directory(
    outputs_dir: str | Path,
    gold_path: str | Path,
    protocol_path: str | Path | None = None,
    *,
    bootstrap_replicates: int | None = None,
) -> dict[str, Any]:
    outputs_path = Path(outputs_dir)
    gold_file = Path(gold_path)
    if not outputs_path.is_dir():
        raise EvaluationError(f"outputs_dir is not a directory: {outputs_path}")
    if not gold_file.is_file():
        raise EvaluationError(f"gold_path is not a file: {gold_file}")
    settings, protocol_meta = _load_protocol(
        protocol_path, bootstrap_replicates=bootstrap_replicates
    )
    gold, aliases = _load_gold(gold_file)
    artifacts: dict[str, ArtifactView] = {}
    hashes: dict[str, str] = {}
    manifests: list[tuple[Path, dict[str, Any]]] = []
    for path in sorted(outputs_path.glob("*.json")):
        payload = _load_json(path)
        if isinstance(payload, dict) and payload.get("schema_version") == MANIFEST_SCHEMA_VERSION:
            manifests.append((path, payload))
            continue
        artifact = _validate_artifact(payload, path=path, aliases=aliases)
        if artifact.patient_id in artifacts:
            raise EvaluationError(f"Duplicate artifact for patient {artifact.patient_id}")
        artifacts[artifact.patient_id] = artifact
        hashes[artifact.patient_id] = _sha256(path)
    if len(manifests) > 1:
        raise EvaluationError("Output directory contains multiple ABC batch manifests")
    manifest_validation = (
        _validate_batch_manifest(
            manifests[0][1], path=manifests[0][0], artifacts=artifacts, hashes=hashes
        )
        if manifests
        else None
    )
    cases = gold["cases"]
    labeled_ids = {case["patient_id"] for case in cases}
    artifact_ids = set(artifacts)
    missing = sorted(labeled_ids - artifact_ids)
    extra = sorted(artifact_ids - labeled_ids)
    scopes = {
        PRIMARY_SCOPE: _evaluate_scope(
            cases=cases,
            artifacts=artifacts,
            positive_cases_only=False,
            genus_relaxed=False,
        ),
        SENSITIVITY_SCOPE: _evaluate_scope(
            cases=cases,
            artifacts=artifacts,
            positive_cases_only=True,
            genus_relaxed=True,
        ),
    }
    pareto = {name: _pareto_report(scope) for name, scope in scopes.items()}
    bootstrap: dict[str, Any] = {}
    named = set(settings["bootstrap"].get("named_arms", []))
    for index, (scope_name, scope) in enumerate(scopes.items()):
        targets = set(named)
        for point in pareto[scope_name]["d_aware_only"]:
            targets.update(point["equivalent_arms"])
        bootstrap[scope_name] = _bootstrap_scope(
            scope,
            targets=[arm for arm in ARM_ORDER if arm in targets],
            replicates=settings["bootstrap"]["replicates"],
            seed=settings["bootstrap"]["seed"] + index,
        )
    cluster = {
        name: _cluster_sensitivity(
            scope,
            settings["organism_cluster"]["minimum_added_tp_for_detail"],
        )
        for name, scope in scopes.items()
    }
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "evaluator_version": EVALUATOR_VERSION,
        "evaluator_sha256": _sha256(Path(__file__)),
        "primary_scope": PRIMARY_SCOPE,
        "sensitivity_scope": SENSITIVITY_SCOPE,
        "primary_anchor": BASE_ARM,
        "future_holdout_candidate_arm": FUTURE_HOLDOUT_CANDIDATE_ARM,
        "arm_order": list(ARM_ORDER),
        "arm_count": len(ARM_ORDER),
        "arm_groups": _arm_groups(),
        "protocol": {"source": protocol_meta, **settings},
        "gold": {key: value for key, value in gold.items() if key != "cases"}
        | {"labeled_case_count": len(cases)},
        "artifacts": {
            "outputs_dir": str(outputs_path.resolve()),
            "artifact_count": len(artifacts),
            "labeled_artifact_count": len(artifact_ids & labeled_ids),
            "missing_labeled_artifact_count": len(missing),
            "missing_labeled_patient_ids": missing,
            "extra_unlabeled_artifact_count": len(extra),
            "extra_unlabeled_patient_ids": extra,
            "batch_manifest_count": len(manifests),
            "batch_manifest_paths": [str(path.resolve()) for path, _ in manifests],
            "batch_manifest_validation": manifest_validation,
            "sha256_by_patient": hashes,
        },
        "scopes": scopes,
        "pareto": pareto,
        "paired_patient_bootstrap": bootstrap,
        "organism_cluster_sensitivity": cluster,
        "notes": [
            "Primary anchor remains frozen L5; every new arm is exploratory on this development cohort.",
            "Primary scope uses all labeled patients, exact species, frozen synonyms, and patient-scoped one-to-one matching.",
            "Sensitivity scope uses answer-positive patients and genus-relaxed one-to-one matching.",
            "Manual candidates are not auto positives; review-inclusive metrics are reported separately.",
            "Missing artifacts remain failures in denominators and missing no-path artifacts are not credited as true negatives.",
            "Predicate unknown remains unknown and never becomes a negative vote.",
            "Bootstrap intervals are within-development-sample stability summaries, not external validation.",
        ],
    }


def _percent(value: Any) -> str:
    if value is None or not isinstance(value, (int, float)) or not math.isfinite(value):
        return "NA"
    return f"{value * 100:.1f}%"


def render_markdown(report: Mapping[str, Any]) -> str:
    if report.get("schema_version") != REPORT_SCHEMA_VERSION:
        raise EvaluationError("Unsupported ABC evaluation report schema")
    lines = [
        "# Post-L5 corrected/raw A/B/C evaluation",
        "",
        f"- Labeled cases: {report['gold']['labeled_case_count']}",
        f"- Complete artifacts: {report['artifacts']['artifact_count']}",
        f"- Missing labeled artifacts retained as failures: {report['artifacts']['missing_labeled_artifact_count']}",
        f"- Arms: {report['arm_count']}",
        f"- Primary anchor: `{report['primary_anchor']}`",
        "- All new arms are exploratory; no post-hoc arm replaces frozen L5.",
    ]
    for scope_name in (report["primary_scope"], report["sensitivity_scope"]):
        scope = report["scopes"][scope_name]
        lines.extend(
            [
                "",
                f"## {scope_name}",
                "",
                f"Cases: {scope['case_count']}",
                "",
                "### Complete 79-arm metrics",
                "",
                "| Arm | Family/scope/gate | Auto TP/FP/FN | Auto P/R/F0.5/F1 | Review TP/FP/FN | Added nonE0 TP/FP | Restored prune TP/FP | Promoted manual TP/FP | Delta P/R vs L5 | No-path specificity |",
                "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for arm in report["arm_order"]:
            row = scope["arms"][arm]
            auto = row["auto"]
            review = row["review_inclusive"]
            attr = row["attribution"]
            delta = row["delta_vs_l5"]["auto"]
            descriptor = f"{row['family'] or '-'}/{row['scope']}/{row['gate']}"
            lines.append(
                f"| {arm} | {descriptor} | {auto['tp']}/{auto['fp']}/{auto['fn']} | "
                f"{_percent(auto['precision'])}/{_percent(auto['recall'])}/{_percent(auto['f0_5'])}/{_percent(auto['f1'])} | "
                f"{review['tp']}/{review['fp']}/{review['fn']} | "
                f"{attr['added_non_e0_auto']['tp']}/{attr['added_non_e0_auto']['fp']} | "
                f"{attr['restored_e0_prune_auto']['tp']}/{attr['restored_e0_prune_auto']['fp']} | "
                f"{attr['promoted_e0_manual_auto']['tp']}/{attr['promoted_e0_manual_auto']['fp']} | "
                f"{_percent(delta['precision'])}/{_percent(delta['recall'])} | "
                f"{_percent(auto['no_path_specificity'])} |"
            )
        lines.extend(["", "### Cumulative marginal transitions", ""])
        lines.extend(
            [
                "| Arm | Previous | Added auto TP/FP | Delta P/R | Added review TP/FP |",
                "|---|---|---:|---:|---:|",
            ]
        )
        for family in FAMILIES:
            for index, predicate in enumerate(CUMULATIVE_PREDICATES, start=1):
                arm = cumulative_arm(family, index, predicate)
                marginal = scope["arms"][arm]["cumulative_marginal"]
                lines.append(
                    f"| {arm} | {marginal['previous_arm']} | "
                    f"{marginal['added_auto']['tp']}/{marginal['added_auto']['fp']} | "
                    f"{_percent(marginal['auto_metric_delta']['precision'])}/{_percent(marginal['auto_metric_delta']['recall'])} | "
                    f"{marginal['added_review']['tp']}/{marginal['added_review']['fp']} |"
                )
        lines.extend(["", "### D-aware Pareto frontier", ""])
        lines.extend(
            [
                "| Representative arm | Equivalent arms | TP/FP/FN | Precision | Recall | F0.5 |",
                "|---|---|---:|---:|---:|---:|",
            ]
        )
        for point in report["pareto"][scope_name]["d_aware_only"]:
            lines.append(
                f"| {point['representative_arm']} | {len(point['equivalent_arms'])} | "
                f"{point['tp']}/{point['fp']}/{point['fn']} | {_percent(point['precision'])} | "
                f"{_percent(point['recall'])} | {_percent(point['f0_5'])} |"
            )
        bootstrap = report["paired_patient_bootstrap"][scope_name]
        lines.extend(
            [
                "",
                f"### Paired patient bootstrap ({bootstrap['replicates']} replicates)",
                "",
                "| Arm | Delta precision (95% CI) | Delta recall (95% CI) |",
                "|---|---:|---:|",
            ]
        )
        for arm, result in bootstrap["arms"].items():
            pci = result["precision_delta_percentile_95_ci"]
            rci = result["recall_delta_percentile_95_ci"]
            lines.append(
                f"| {arm} | {_percent(result['point_precision_delta'])} "
                f"({_percent(pci[0])} to {_percent(pci[1])}) | "
                f"{_percent(result['point_recall_delta'])} "
                f"({_percent(rci[0])} to {_percent(rci[1])}) |"
            )
        clusters = report["organism_cluster_sensitivity"][scope_name]
        lines.extend(["", "### Organism-cluster sensitivity", ""])
        if not clusters:
            lines.append("No arm added a non-E0 true positive in this scope.")
        else:
            lines.extend(
                [
                    "| Arm | Added TP/FP | Organisms | Maximum single-organism TP share |",
                    "|---|---:|---:|---:|",
                ]
            )
            for arm, result in clusters.items():
                lines.append(
                    f"| {arm} | {result['added_non_e0_tp']}/{result['added_non_e0_fp']} | "
                    f"{result['distinct_added_organism_count']} | "
                    f"{_percent(result['maximum_single_organism_tp_share'])} |"
                )
    lines.extend(["", "## Interpretation notes", ""])
    lines.extend(f"- {note}" for note in report.get("notes", []))
    return "\n".join(lines) + "\n"


def write_json(report: Mapping[str, Any], path: str | Path) -> Path:
    if report.get("schema_version") != REPORT_SCHEMA_VERSION:
        raise EvaluationError("Unsupported ABC evaluation report schema")
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return output


def write_markdown(report: Mapping[str, Any], path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_markdown(report), encoding="utf-8")
    return output


__all__ = [
    "E0_REFERENCE",
    "EvaluationError",
    "MANIFEST_SCHEMA_VERSION",
    "OUTPUT_SCHEMA_VERSION",
    "PRIMARY_SCOPE",
    "REPORT_SCHEMA_VERSION",
    "SENSITIVITY_SCOPE",
    "evaluate_directory",
    "render_markdown",
    "write_json",
    "write_markdown",
]

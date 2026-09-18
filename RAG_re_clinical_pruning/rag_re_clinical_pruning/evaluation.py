from __future__ import annotations

import hashlib
import json
import math
import re
import unicodedata
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


OUTPUT_SCHEMA_VERSION = "rag_re_clinical_pruning.output.v1"
MANIFEST_SCHEMA_VERSION = "rag_re_clinical_pruning.batch_manifest.v1"
REPORT_SCHEMA_VERSION = "rag_re_clinical_pruning.evaluation.v1"
BASELINE_ARM = "S0_E0"
PRIMARY_ARM = "L5_PLUS_DIRECT_CONVERGENT_RESCUE"

STANDALONE_ARMS = (
    BASELINE_ARM,
    "S_Q_WEAK_RAW",
    "S_Q_GUARD",
    "S_D_BLOCK_DIAGNOSTIC",
    "S_R_DIRECT",
    "S_R_CONVERGENT",
    "S_R_HIGH_STRICT",
    "S_R_HIGH_ANY",
    "S_R_CONTEXT_A",
)
CUMULATIVE_ARMS = (
    "L0_E0",
    "L1_WEAK_RAW",
    "L2_DIRECT_RESTORE",
    "L3_CONTEXT_REVIEW",
    "L4_PLUS_GUARD_REVIEW",
    PRIMARY_ARM,
    "L6_PLUS_HIGH_STRICT",
    "L7_PLUS_HIGH_ANY",
    "L8_PLUS_CONTEXT_A",
)
FIXED_ARM_NAMES = STANDALONE_ARMS + CUMULATIVE_ARMS
E0_EQUIVALENT_ARMS = {BASELINE_ARM, "L0_E0"}

DEFAULT_SAFETY_THRESHOLDS = {
    "auto_max_absolute_recall_drop": 0.10,
    "review_inclusive_max_absolute_recall_drop": 0.035,
    "max_explicitly_pruned_tp": 1,
}

# Frozen before outcomes are scored. Gold may add non-conflicting, versioned aliases.
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

E0_STATES = {"KEEP", "MANUAL", "PRUNE"}
NON_E0_STATES = {"RESCUE", "MANUAL", "NO_RESCUE"}
ALL_STATES = E0_STATES | NON_E0_STATES


class EvaluationError(ValueError):
    """Raised when artifacts, protocol, or gold violate the frozen contract."""


@dataclass(frozen=True)
class CandidateView:
    name: str
    baseline_selected: bool
    states: Mapping[str, str]


@dataclass(frozen=True)
class ArmView:
    auto: tuple[str, ...]
    review_inclusive: tuple[str, ...]
    manual: tuple[str, ...]
    pruned: tuple[str, ...]
    rescued: tuple[str, ...]
    no_rescue: tuple[str, ...]
    family: str


@dataclass(frozen=True)
class ArtifactView:
    patient_id: str
    path: str
    candidates: tuple[CandidateView, ...]
    arms: Mapping[str, ArmView]


def _normalize(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).strip().lower()
    return re.sub(r"\s+", " ", text)


def _require_name(value: Any, *, where: str) -> str:
    normalized = _normalize(value)
    if not normalized:
        raise EvaluationError(f"{where} must be a non-empty pathogen name")
    return normalized


def _strict_bool(value: Any, *, where: str) -> bool:
    if value is True:
        return True
    if value is False:
        return False
    raise EvaluationError(f"{where} must be a boolean")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
        existing = aliases.get(key)
        if existing is not None and existing != value:
            raise EvaluationError(
                f"Gold alias {raw_key!r} conflicts with frozen mapping "
                f"{existing!r} != {value!r}"
            )
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
            raise EvaluationError(f"Alias cycle encountered for {value!r}")
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


def _genus(value: str) -> str | None:
    match = re.match(r"^([a-z][a-z-]+)(?:\s+)([a-z][a-z0-9.-]+)", value)
    if not match:
        return None
    token = match.group(1)
    return None if token in GENUS_STOPWORDS else token


def _names_match(left: str, right: str, *, genus_relaxed: bool) -> bool:
    if left == right:
        return True
    left_genus = _genus(left)
    return bool(genus_relaxed and left_genus and left_genus == _genus(right))


def _maximum_match_count(
    truth: Sequence[str], predictions: Sequence[str], *, genus_relaxed: bool
) -> int:
    """Maximum one-to-one matching, scoped to one patient."""

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
    beta_squared = beta * beta
    denominator = beta_squared * precision + recall
    return (
        (1 + beta_squared) * precision * recall / denominator
        if denominator
        else 0.0
    )


def _resolve_arm_names(arm_names: Sequence[str] | None) -> tuple[str, ...]:
    if arm_names is None:
        return FIXED_ARM_NAMES
    if isinstance(arm_names, (str, bytes)):
        raise EvaluationError("arm_names must be a sequence, not a string")
    selected = tuple(str(value).strip() for value in arm_names)
    if not selected or any(not value for value in selected):
        raise EvaluationError("arm_names must contain non-empty names")
    if len(set(selected)) != len(selected):
        raise EvaluationError("arm_names contains duplicates")
    unsupported = sorted(set(selected) - set(FIXED_ARM_NAMES))
    if unsupported:
        raise EvaluationError(f"Unsupported pruning arm(s): {unsupported}")
    # Provenance and deltas require the frozen E0 plus all earlier cumulative arms.
    required = {BASELINE_ARM}
    cumulative_indices = [CUMULATIVE_ARMS.index(arm) for arm in selected if arm in CUMULATIVE_ARMS]
    if cumulative_indices:
        required.update(CUMULATIVE_ARMS[: max(cumulative_indices) + 1])
    ordered = tuple(arm for arm in FIXED_ARM_NAMES if arm in set(selected) | required)
    return ordered


def _expected_family(arm: str) -> str:
    return "standalone" if arm in STANDALONE_ARMS else "cumulative"


def _validate_count(raw_arm: Mapping[str, Any], field: str, length: int, *, where: str) -> None:
    count_field = field.removesuffix("_pathogens") + "_count"
    value = raw_arm.get(count_field)
    if type(value) is not int or value != length:
        raise EvaluationError(f"{where}.{count_field} does not match {field}")


def _validate_state_counts(
    raw: Any, expected: Counter[str], *, where: str
) -> None:
    if not isinstance(raw, dict):
        raise EvaluationError(f"{where}.state_counts must be an object")
    normalized: dict[str, int] = {}
    for key, value in raw.items():
        state = str(key).strip().upper()
        if state not in ALL_STATES or type(value) is not int or value < 0:
            raise EvaluationError(f"{where}.state_counts has invalid entry {key!r}: {value!r}")
        normalized[state] = normalized.get(state, 0) + value
    if {k: v for k, v in normalized.items() if v} != dict(expected):
        raise EvaluationError(f"{where}.state_counts disagrees with candidate arm_states")


def _validate_artifact(
    payload: Any,
    *,
    path: Path,
    aliases: Mapping[str, str],
    arm_names: Sequence[str],
) -> ArtifactView:
    if not isinstance(payload, dict):
        raise EvaluationError(f"Pruning artifact {path} must be an object")
    if payload.get("schema_version") != OUTPUT_SCHEMA_VERSION:
        raise EvaluationError(
            f"Pruning artifact {path} has unsupported schema_version "
            f"{payload.get('schema_version')!r}"
        )
    patient_id = str(payload.get("patient_id") or "").strip()
    if not patient_id:
        raise EvaluationError(f"Pruning artifact {path} is missing patient_id")
    if payload.get("run_status") != "complete":
        raise EvaluationError(
            f"Pruning artifact for patient {patient_id} must have run_status='complete'"
        )

    raw_candidates = payload.get("candidates")
    if not isinstance(raw_candidates, list):
        raise EvaluationError(f"Patient {patient_id} candidates must be a list")
    candidates: list[CandidateView] = []
    seen_names: set[str] = set()
    for index, raw_candidate in enumerate(raw_candidates):
        where = f"patient {patient_id} candidate[{index}]"
        if not isinstance(raw_candidate, dict):
            raise EvaluationError(f"{where} must be an object")
        name = _canonical(
            raw_candidate.get("canonical_organism_name")
            or raw_candidate.get("organism_name"),
            aliases,
        )
        if name in seen_names:
            raise EvaluationError(
                f"Patient {patient_id} has duplicate or alias-equivalent candidate {name!r}"
            )
        seen_names.add(name)
        baseline = _strict_bool(
            raw_candidate.get("baseline_selected"),
            where=f"{where}.baseline_selected",
        )
        raw_states = raw_candidate.get("arm_states")
        if not isinstance(raw_states, dict):
            raise EvaluationError(f"{where}.arm_states must be an object")
        states: dict[str, str] = {}
        for arm in arm_names:
            if arm not in raw_states:
                raise EvaluationError(f"{where}.arm_states is missing {arm}")
            raw_state = raw_states[arm]
            if not isinstance(raw_state, dict):
                raise EvaluationError(f"{where}.arm_states[{arm}] must be an object")
            state = str(raw_state.get("state") or "").strip().upper()
            valid_states = E0_STATES if baseline else NON_E0_STATES
            if state not in valid_states:
                raise EvaluationError(
                    f"{where}.arm_states[{arm}].state={state!r} invalid for "
                    f"baseline_selected={baseline}"
                )
            if not isinstance(raw_state.get("rule_ledger"), list):
                raise EvaluationError(f"{where}.arm_states[{arm}].rule_ledger must be a list")
            states[arm] = state
        candidates.append(
            CandidateView(name=name, baseline_selected=baseline, states=states)
        )

    candidate_set = {candidate.name for candidate in candidates}
    baseline_set = {candidate.name for candidate in candidates if candidate.baseline_selected}
    raw_arms = payload.get("arms")
    if not isinstance(raw_arms, dict):
        raise EvaluationError(f"Patient {patient_id} arms must be an object")
    arms: dict[str, ArmView] = {}
    list_fields = (
        "auto_positive_pathogens",
        "review_inclusive_pathogens",
        "manual_review_pathogens",
        "pruned_pathogens",
        "rescued_pathogens",
        "no_rescue_pathogens",
    )
    for arm in arm_names:
        where = f"patient {patient_id} arm {arm}"
        raw_arm = raw_arms.get(arm)
        if not isinstance(raw_arm, dict):
            raise EvaluationError(f"{where} must be an object")
        lists = {
            field: _canonical_list(raw_arm.get(field), aliases, where=f"{where}.{field}")
            for field in list_fields
        }
        for field, values in lists.items():
            _validate_count(raw_arm, field, len(values), where=where)
            outside = set(values) - candidate_set
            if outside:
                raise EvaluationError(
                    f"{where}.{field} contains pathogens outside candidates: {sorted(outside)}"
                )
        by_state: dict[str, set[str]] = {state: set() for state in ALL_STATES}
        state_counter: Counter[str] = Counter()
        for candidate in candidates:
            state = candidate.states[arm]
            by_state[state].add(candidate.name)
            state_counter[state] += 1
        expected_auto = by_state["KEEP"] | by_state["RESCUE"]
        expected_manual = by_state["MANUAL"]
        expected_review = expected_auto | expected_manual
        expected_pruned = by_state["PRUNE"]
        expected_rescued = by_state["RESCUE"]
        expected_no_rescue = by_state["NO_RESCUE"]
        comparisons = {
            "auto_positive_pathogens": expected_auto,
            "review_inclusive_pathogens": expected_review,
            "manual_review_pathogens": expected_manual,
            "pruned_pathogens": expected_pruned,
            "rescued_pathogens": expected_rescued,
            "no_rescue_pathogens": expected_no_rescue,
        }
        for field, expected in comparisons.items():
            if set(lists[field]) != expected:
                raise EvaluationError(
                    f"{where}.{field} disagrees with candidate arm_states"
                )
        _validate_state_counts(raw_arm.get("state_counts"), state_counter, where=where)
        if not isinstance(raw_arm.get("rule_sequence"), list):
            raise EvaluationError(f"{where}.rule_sequence must be a list")
        family = str(raw_arm.get("family") or "").strip().lower()
        if family != _expected_family(arm):
            raise EvaluationError(
                f"{where}.family must be {_expected_family(arm)!r}, got {family!r}"
            )
        arms[arm] = ArmView(
            auto=lists["auto_positive_pathogens"],
            review_inclusive=lists["review_inclusive_pathogens"],
            manual=lists["manual_review_pathogens"],
            pruned=lists["pruned_pathogens"],
            rescued=lists["rescued_pathogens"],
            no_rescue=lists["no_rescue_pathogens"],
            family=family,
        )

    for arm in E0_EQUIVALENT_ARMS & set(arm_names):
        if set(arms[arm].auto) != baseline_set:
            raise EvaluationError(f"Patient {patient_id} {arm} must exactly equal E0")
        if arms[arm].manual or arms[arm].pruned or arms[arm].rescued:
            raise EvaluationError(f"Patient {patient_id} {arm} must not prune/review/rescue")
    return ArtifactView(
        patient_id=patient_id,
        path=str(path.resolve()),
        candidates=tuple(candidates),
        arms=arms,
    )


def _load_gold(path: Path) -> tuple[dict[str, Any], dict[str, str]]:
    payload = _load_json(path)
    if not isinstance(payload, dict):
        raise EvaluationError("Gold file must contain an object")
    cases = payload.get("cases")
    if not isinstance(cases, list) or not cases:
        raise EvaluationError("gold.cases must be a non-empty list")
    aliases = _build_aliases(payload.get("aliases"))
    normalized_cases: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise EvaluationError(f"gold.cases[{index}] must be an object")
        patient_id = str(case.get("patient_id") or "").strip()
        if not patient_id:
            raise EvaluationError(f"gold.cases[{index}] is missing patient_id")
        if patient_id in seen_ids:
            raise EvaluationError(f"Gold contains duplicate patient_id {patient_id}")
        seen_ids.add(patient_id)
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
        overlap = set(truth) & set(uncertain)
        if overlap:
            raise EvaluationError(
                f"Gold patient {patient_id} has positive/uncertain overlap: {sorted(overlap)}"
            )
        normalized_cases.append(
            {"patient_id": patient_id, "truth": truth, "uncertain": uncertain}
        )
    return (
        {
            "schema_version": payload.get("schema_version"),
            "label_semantics": payload.get("label_semantics"),
            "complete_candidate_negative_labels": payload.get(
                "complete_candidate_negative_labels"
            ),
            "cases": normalized_cases,
            "path": str(path.resolve()),
            "sha256": _sha256(path),
        },
        aliases,
    )


def _metric_counts(
    truth: Sequence[str], predictions: Sequence[str], *, genus_relaxed: bool
) -> dict[str, int]:
    matched = _maximum_match_count(truth, predictions, genus_relaxed=genus_relaxed)
    return {"tp": matched, "fp": len(predictions) - matched, "fn": len(truth) - matched}


def _effect_added(
    truth: Sequence[str],
    base: Sequence[str],
    added: Sequence[str],
    *,
    genus_relaxed: bool,
) -> dict[str, Any]:
    base_unique = tuple(dict.fromkeys(base))
    added_unique = tuple(value for value in dict.fromkeys(added) if value not in set(base_unique))
    before = _maximum_match_count(truth, base_unique, genus_relaxed=genus_relaxed)
    after = _maximum_match_count(
        truth, base_unique + added_unique, genus_relaxed=genus_relaxed
    )
    tp = after - before
    return {
        "count": len(added_unique),
        "tp": tp,
        "fp": len(added_unique) - tp,
        "ppv": _safe_div(tp, len(added_unique)),
        "pathogens": list(added_unique),
    }


def _effect_removed(
    truth: Sequence[str],
    full: Sequence[str],
    removed: Sequence[str],
    *,
    genus_relaxed: bool,
) -> dict[str, Any]:
    removed_set = set(removed)
    full_unique = tuple(dict.fromkeys(full))
    removed_unique = tuple(value for value in full_unique if value in removed_set)
    remaining = tuple(value for value in full_unique if value not in removed_set)
    before = _maximum_match_count(truth, full_unique, genus_relaxed=genus_relaxed)
    after = _maximum_match_count(truth, remaining, genus_relaxed=genus_relaxed)
    tp = before - after
    return {
        "count": len(removed_unique),
        "tp": tp,
        "fp": len(removed_unique) - tp,
        "ppv": _safe_div(tp, len(removed_unique)),
        "pathogens": list(removed_unique),
    }


def _sum_effects(rows: Sequence[Mapping[str, Any]], key: str) -> dict[str, Any]:
    count = sum(int(row["attribution"][key]["count"]) for row in rows)
    tp = sum(int(row["attribution"][key]["tp"]) for row in rows)
    fp = sum(int(row["attribution"][key]["fp"]) for row in rows)
    return {"count": count, "tp": tp, "fp": fp, "ppv": _safe_div(tp, count)}


def _aggregate_prediction_metrics(
    rows: Sequence[Mapping[str, Any]], prediction_kind: str
) -> dict[str, Any]:
    tp = sum(int(row[prediction_kind]["tp"]) for row in rows)
    fp = sum(int(row[prediction_kind]["fp"]) for row in rows)
    fn = sum(int(row[prediction_kind]["fn"]) for row in rows)
    precision = _safe_div(tp, tp + fp)
    recall = _safe_div(tp, tp + fn)
    no_path_total = sum(bool(row["no_path_case"]) for row in rows)
    no_path_correct = sum(
        bool(row["no_path_case"] and row["artifact_present"] and row[prediction_kind]["predicted_count"] == 0)
        for row in rows
    )
    no_path_positive = sum(
        bool(row["no_path_case"] and row["artifact_present"] and row[prediction_kind]["predicted_count"] > 0)
        for row in rows
    )
    no_path_missing = sum(
        bool(row["no_path_case"] and not row["artifact_present"]) for row in rows
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
    output: dict[str, Any] = {}
    for key in keys:
        left = current.get(key)
        right = reference.get(key)
        output[key] = None if left is None or right is None else left - right
    return output


def _candidate_names_by_transition(
    artifact: ArtifactView,
    arm: str,
    previous_arm: str | None,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    prior_cumulative = CUMULATIVE_ARMS[: CUMULATIVE_ARMS.index(arm)] if arm in CUMULATIVE_ARMS else ()
    restored_ever: list[str] = []
    stage_restored: list[str] = []
    stage_non_e0: list[str] = []
    for candidate in artifact.candidates:
        current = candidate.states[arm]
        if candidate.baseline_selected:
            if current == "KEEP" and any(candidate.states[past] == "PRUNE" for past in prior_cumulative):
                restored_ever.append(candidate.name)
            if previous_arm is not None and candidate.states[previous_arm] == "PRUNE" and current == "KEEP":
                stage_restored.append(candidate.name)
        elif previous_arm is not None and candidate.states[previous_arm] != "RESCUE" and current == "RESCUE":
            stage_non_e0.append(candidate.name)
    return tuple(restored_ever), tuple(stage_restored), tuple(stage_non_e0)


def _evaluate_scope(
    *,
    cases: Sequence[Mapping[str, Any]],
    artifacts: Mapping[str, ArtifactView],
    arm_names: Sequence[str],
    positive_cases_only: bool,
    genus_relaxed: bool,
    safety_thresholds: Mapping[str, float | int],
) -> dict[str, Any]:
    scoped_cases = [case for case in cases if not positive_cases_only or case["truth"]]
    arm_rows: dict[str, list[dict[str, Any]]] = {arm: [] for arm in arm_names}
    previous_cumulative: dict[str, str | None] = {}
    for index, arm in enumerate(CUMULATIVE_ARMS):
        previous_cumulative[arm] = CUMULATIVE_ARMS[index - 1] if index else None

    for case in scoped_cases:
        patient_id = str(case["patient_id"])
        truth = tuple(case["truth"])
        uncertain = set(case["uncertain"])
        artifact = artifacts.get(patient_id)
        for arm in arm_names:
            previous_arm = previous_cumulative.get(arm)
            if artifact is None:
                auto: tuple[str, ...] = ()
                review: tuple[str, ...] = ()
                manual: tuple[str, ...] = ()
                e0: tuple[str, ...] = ()
                pruned: tuple[str, ...] = ()
                rescued: tuple[str, ...] = ()
                restored_ever: tuple[str, ...] = ()
                stage_restored: tuple[str, ...] = ()
                stage_non_e0: tuple[str, ...] = ()
                candidate_count = 0
            else:
                def neutral_filtered(values: Iterable[str]) -> tuple[str, ...]:
                    return tuple(value for value in values if value not in uncertain)

                view = artifact.arms[arm]
                auto = neutral_filtered(view.auto)
                review = neutral_filtered(view.review_inclusive)
                manual = neutral_filtered(view.manual)
                e0 = neutral_filtered(artifact.arms[BASELINE_ARM].auto)
                pruned = neutral_filtered(view.pruned)
                rescued = neutral_filtered(view.rescued)
                restored_ever, stage_restored, stage_non_e0 = _candidate_names_by_transition(
                    artifact, arm, previous_arm
                )
                restored_ever = neutral_filtered(restored_ever)
                stage_restored = neutral_filtered(stage_restored)
                stage_non_e0 = neutral_filtered(stage_non_e0)
                candidate_count = len(artifact.candidates)

            auto_counts = _metric_counts(truth, auto, genus_relaxed=genus_relaxed)
            review_counts = _metric_counts(truth, review, genus_relaxed=genus_relaxed)
            explicit_pruned = _effect_removed(
                truth, e0, pruned, genus_relaxed=genus_relaxed
            )
            manual_effect = _effect_added(
                truth, auto, manual, genus_relaxed=genus_relaxed
            )
            restored_effect = _effect_added(
                truth,
                tuple(value for value in auto if value not in set(restored_ever)),
                restored_ever,
                genus_relaxed=genus_relaxed,
            )
            rescued_effect = _effect_added(
                truth,
                tuple(value for value in auto if value not in set(rescued)),
                rescued,
                genus_relaxed=genus_relaxed,
            )
            stage_restored_effect = _effect_added(
                truth,
                tuple(value for value in auto if value not in set(stage_restored)),
                stage_restored,
                genus_relaxed=genus_relaxed,
            )
            stage_rescue_effect = _effect_added(
                truth,
                tuple(value for value in auto if value not in set(stage_non_e0)),
                stage_non_e0,
                genus_relaxed=genus_relaxed,
            )
            arm_rows[arm].append(
                {
                    "patient_id": patient_id,
                    "artifact_present": artifact is not None,
                    "truth_count": len(truth),
                    "candidate_count": candidate_count,
                    "no_path_case": not truth,
                    "auto": auto_counts | {
                        "predicted_count": len(auto),
                        "predicted_pathogens": list(auto),
                    },
                    "review_inclusive": review_counts | {
                        "predicted_count": len(review),
                        "predicted_pathogens": list(review),
                    },
                    "manual_review_pathogens": list(manual),
                    "attribution": {
                        "explicit_pruned": explicit_pruned,
                        "manual": manual_effect,
                        "restored_from_prune": restored_effect,
                        "added_non_e0_rescue": rescued_effect,
                        "stage_restored_from_prune": stage_restored_effect,
                        "stage_added_non_e0_rescue": stage_rescue_effect,
                    },
                }
            )

    output_arms: dict[str, Any] = {}
    for arm, rows in arm_rows.items():
        auto_metrics = _aggregate_prediction_metrics(rows, "auto")
        review_metrics = _aggregate_prediction_metrics(rows, "review_inclusive")
        output_arms[arm] = {
            "family": _expected_family(arm),
            "auto": auto_metrics,
            "review_inclusive": review_metrics,
            "attribution": {
                key: _sum_effects(rows, key)
                for key in (
                    "explicit_pruned",
                    "manual",
                    "restored_from_prune",
                    "added_non_e0_rescue",
                    "stage_restored_from_prune",
                    "stage_added_non_e0_rescue",
                )
            },
            "case_counts": {
                "included": len(rows),
                "artifact_present": sum(bool(row["artifact_present"]) for row in rows),
                "artifact_missing": sum(not row["artifact_present"] for row in rows),
                "answer_positive": sum(not row["no_path_case"] for row in rows),
                "no_path": sum(bool(row["no_path_case"]) for row in rows),
                "with_auto_positive": sum(row["auto"]["predicted_count"] > 0 for row in rows),
                "with_manual_review": sum(bool(row["manual_review_pathogens"]) for row in rows),
            },
            "cases": rows,
        }

    e0 = output_arms[BASELINE_ARM]
    for arm in arm_names:
        metrics = output_arms[arm]
        metrics["delta_vs_e0"] = {
            "reference_arm": BASELINE_ARM,
            "auto": _metric_delta(metrics["auto"], e0["auto"]),
            "review_inclusive": _metric_delta(
                metrics["review_inclusive"], e0["review_inclusive"]
            ),
        }
        previous = previous_cumulative.get(arm)
        if previous is not None and previous in output_arms:
            metrics["delta_vs_previous_cumulative"] = {
                "reference_arm": previous,
                "auto": _metric_delta(metrics["auto"], output_arms[previous]["auto"]),
                "review_inclusive": _metric_delta(
                    metrics["review_inclusive"],
                    output_arms[previous]["review_inclusive"],
                ),
            }
        else:
            metrics["delta_vs_previous_cumulative"] = None

        auto_recall = metrics["auto"]["recall"]
        review_recall = metrics["review_inclusive"]["recall"]
        e0_auto_recall = e0["auto"]["recall"]
        e0_review_recall = e0["review_inclusive"]["recall"]
        auto_drop = (
            None
            if auto_recall is None or e0_auto_recall is None
            else max(0.0, e0_auto_recall - auto_recall)
        )
        review_drop = (
            None
            if review_recall is None or e0_review_recall is None
            else max(0.0, e0_review_recall - review_recall)
        )
        pruned_tp = metrics["attribution"]["explicit_pruned"]["tp"]
        criteria = {
            "auto_absolute_recall_drop": {
                "value": auto_drop,
                "maximum": safety_thresholds["auto_max_absolute_recall_drop"],
                "pass": auto_drop is not None
                and auto_drop <= float(safety_thresholds["auto_max_absolute_recall_drop"]) + 1e-12,
            },
            "review_inclusive_absolute_recall_drop": {
                "value": review_drop,
                "maximum": safety_thresholds[
                    "review_inclusive_max_absolute_recall_drop"
                ],
                "pass": review_drop is not None
                and review_drop
                <= float(
                    safety_thresholds["review_inclusive_max_absolute_recall_drop"]
                )
                + 1e-12,
            },
            "explicitly_pruned_tp": {
                "value": pruned_tp,
                "maximum": safety_thresholds["max_explicitly_pruned_tp"],
                "pass": pruned_tp <= int(safety_thresholds["max_explicitly_pruned_tp"]),
            },
        }
        metrics["recall_safety"] = {
            "status": "PASS" if all(item["pass"] for item in criteria.values()) else "UNSAFE",
            "evaluation_only": True,
            "criteria": criteria,
        }
    return {
        "case_count": len(scoped_cases),
        "positive_cases_only": positive_cases_only,
        "genus_relaxed": genus_relaxed,
        "arms": output_arms,
    }


def _load_safety_thresholds(protocol_path: str | Path | None) -> tuple[dict[str, float | int], dict[str, Any] | None]:
    thresholds: dict[str, float | int] = dict(DEFAULT_SAFETY_THRESHOLDS)
    protocol_meta: dict[str, Any] | None = None
    if protocol_path is None:
        return thresholds, protocol_meta
    path = Path(protocol_path)
    payload = _load_json(path)
    if not isinstance(payload, dict) or payload.get("schema_version") != "rag_re_clinical_pruning.evaluation_protocol.v1":
        raise EvaluationError("Unsupported pruning evaluation protocol schema")
    if tuple(payload.get("arm_order", [])) != FIXED_ARM_NAMES:
        raise EvaluationError("Protocol arm_order differs from frozen evaluator arms")
    raw = payload.get("recall_safety")
    if not isinstance(raw, dict):
        raise EvaluationError("Protocol recall_safety must be an object")
    for key, default in DEFAULT_SAFETY_THRESHOLDS.items():
        value = raw.get(key)
        if isinstance(default, int):
            if type(value) is not int or value < 0:
                raise EvaluationError(f"Protocol {key} must be a non-negative integer")
        elif not isinstance(value, (int, float)) or isinstance(value, bool) or not 0 <= float(value) <= 1:
            raise EvaluationError(f"Protocol {key} must be between 0 and 1")
        thresholds[key] = value
    protocol_meta = {"path": str(path.resolve()), "sha256": _sha256(path)}
    return thresholds, protocol_meta


def evaluate_directory(
    outputs_dir: str | Path,
    gold_path: str | Path,
    arm_names: Sequence[str] | None = None,
    protocol_path: str | Path | None = None,
) -> dict[str, Any]:
    """Evaluate locked pruning artifacts without changing their decisions.

    Missing labeled artifacts stay in every denominator as empty predictions.
    Thus positive labels become false negatives, while missing no-path cases are
    conservatively not credited as correct negatives.
    """

    outputs_path = Path(outputs_dir)
    gold_file = Path(gold_path)
    if not outputs_path.is_dir():
        raise EvaluationError(f"outputs_dir is not a directory: {outputs_path}")
    if not gold_file.is_file():
        raise EvaluationError(f"gold_path is not a file: {gold_file}")
    selected_arms = _resolve_arm_names(arm_names)
    thresholds, protocol_meta = _load_safety_thresholds(protocol_path)
    gold, aliases = _load_gold(gold_file)

    artifacts: dict[str, ArtifactView] = {}
    artifact_hashes: dict[str, str] = {}
    manifests: list[str] = []
    for path in sorted(outputs_path.glob("*.json")):
        payload = _load_json(path)
        if isinstance(payload, dict) and payload.get("schema_version") == MANIFEST_SCHEMA_VERSION:
            manifests.append(str(path.resolve()))
            continue
        artifact = _validate_artifact(
            payload, path=path, aliases=aliases, arm_names=selected_arms
        )
        if artifact.patient_id in artifacts:
            raise EvaluationError(
                f"Duplicate pruning artifact for patient {artifact.patient_id}"
            )
        artifacts[artifact.patient_id] = artifact
        artifact_hashes[artifact.patient_id] = _sha256(path)

    cases = gold["cases"]
    labeled_ids = {str(case["patient_id"]) for case in cases}
    artifact_ids = set(artifacts)
    missing_ids = sorted(labeled_ids - artifact_ids)
    extra_ids = sorted(artifact_ids - labeled_ids)
    scope_specs = {
        "all_labeled_exact_frozen_synonyms": (False, False),
        "answer_positive_only_genus_relaxed": (True, True),
    }
    scopes = {
        scope_name: _evaluate_scope(
            cases=cases,
            artifacts=artifacts,
            arm_names=selected_arms,
            positive_cases_only=positive_only,
            genus_relaxed=genus_relaxed,
            safety_thresholds=thresholds,
        )
        for scope_name, (positive_only, genus_relaxed) in scope_specs.items()
    }
    safety_summary = {
        arm: {
            "status": (
                "PASS"
                if all(scopes[scope]["arms"][arm]["recall_safety"]["status"] == "PASS" for scope in scopes)
                else "UNSAFE"
            ),
            "evaluation_only": True,
            "scope_statuses": {
                scope: scopes[scope]["arms"][arm]["recall_safety"]["status"]
                for scope in scopes
            },
        }
        for arm in selected_arms
    }
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "primary_scope": "all_labeled_exact_frozen_synonyms",
        "sensitivity_scope": "answer_positive_only_genus_relaxed",
        "primary_arm": PRIMARY_ARM if PRIMARY_ARM in selected_arms else None,
        "arm_names": list(selected_arms),
        "protocol": {
            "source": protocol_meta,
            "recall_safety": thresholds,
            "safety_is_evaluation_only": True,
        },
        "gold": {key: value for key, value in gold.items() if key != "cases"}
        | {"labeled_case_count": len(cases)},
        "artifacts": {
            "outputs_dir": str(outputs_path.resolve()),
            "artifact_count": len(artifacts),
            "labeled_artifact_count": len(artifact_ids & labeled_ids),
            "missing_labeled_artifact_count": len(missing_ids),
            "missing_labeled_patient_ids": missing_ids,
            "extra_unlabeled_artifact_count": len(extra_ids),
            "extra_unlabeled_patient_ids": extra_ids,
            "batch_manifest_count": len(manifests),
            "batch_manifest_paths": manifests,
            "sha256_by_patient": artifact_hashes,
        },
        "scopes": scopes,
        "recall_safety_summary": safety_summary,
        "notes": [
            "Primary: all labeled patients, exact normalized species, frozen synonyms, one-to-one patient matching.",
            "Sensitivity: answer-positive patients only, one-to-one genus-relaxed matching.",
            "Review-inclusive predictions are auto positives plus the manual queue; manual is never an auto positive.",
            "Uncertain gold labels are neutral: matching predictions are excluded from TP and FP.",
            "Missing artifacts remain in denominators; missing positive cases become FN and missing no-path cases are not TN.",
            "Recall-safety PASS/UNSAFE is reporting only and never changes predictions.",
            "Attribution uses counterfactual one-to-one matching, not name-wise gold leakage.",
        ],
    }


def _percent(value: Any) -> str:
    if value is None or not isinstance(value, (int, float)) or not math.isfinite(value):
        return "NA"
    return f"{value * 100:.1f}%"


def render_markdown(report: Mapping[str, Any]) -> str:
    if report.get("schema_version") != REPORT_SCHEMA_VERSION:
        raise EvaluationError("Unsupported pruning evaluation report schema")
    lines = [
        "# RAG_re clinical pruning evaluation",
        "",
        f"- Labeled cases: {report['gold']['labeled_case_count']}",
        f"- Artifacts: {report['artifacts']['artifact_count']}",
        f"- Missing labeled artifacts (retained as failures): {report['artifacts']['missing_labeled_artifact_count']}",
        f"- Primary arm: {report.get('primary_arm') or 'not evaluated'}",
        "- Recall safety: evaluation-only; it does not alter any prediction.",
    ]
    missing = report["artifacts"].get("missing_labeled_patient_ids") or []
    if missing:
        lines.append("- Missing patient IDs: " + ", ".join(map(str, missing)))
    for scope_name in (report["primary_scope"], report["sensitivity_scope"]):
        scope = report["scopes"][scope_name]
        lines.extend(
            [
                "",
                f"## {scope_name}",
                "",
                f"Cases: {scope['case_count']}",
                "",
                "| Arm | Auto TP/FP/FN | Auto P/R/F0.5 | Review TP/FP/FN | Review P/R | Pruned TP/FP | Manual TP/FP | Restored TP/FP | nonE0 TP/FP | Delta P/Delta R vs E0 | Safety |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for arm in report["arm_names"]:
            row = scope["arms"][arm]
            auto = row["auto"]
            review = row["review_inclusive"]
            attr = row["attribution"]
            delta = row["delta_vs_e0"]["auto"]
            lines.append(
                f"| {arm} | {auto['tp']}/{auto['fp']}/{auto['fn']} | "
                f"{_percent(auto['precision'])}/{_percent(auto['recall'])}/{_percent(auto['f0_5'])} | "
                f"{review['tp']}/{review['fp']}/{review['fn']} | "
                f"{_percent(review['precision'])}/{_percent(review['recall'])} | "
                f"{attr['explicit_pruned']['tp']}/{attr['explicit_pruned']['fp']} | "
                f"{attr['manual']['tp']}/{attr['manual']['fp']} | "
                f"{attr['restored_from_prune']['tp']}/{attr['restored_from_prune']['fp']} | "
                f"{attr['added_non_e0_rescue']['tp']}/{attr['added_non_e0_rescue']['fp']} | "
                f"{_percent(delta['precision'])}/{_percent(delta['recall'])} | "
                f"{row['recall_safety']['status']} |"
            )
        lines.extend(
            [
                "",
                "### Cumulative step deltas",
                "",
                "| Arm | Previous | Auto Delta TP | Auto Delta FP | Auto Delta P | Auto Delta R | Stage restored TP/FP | Stage nonE0 TP/FP |",
                "|---|---|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for arm in CUMULATIVE_ARMS:
            if arm not in scope["arms"]:
                continue
            row = scope["arms"][arm]
            previous = row["delta_vs_previous_cumulative"]
            if previous is None:
                lines.append(f"| {arm} | -- | 0 | 0 | 0.0% | 0.0% | 0/0 | 0/0 |")
                continue
            delta = previous["auto"]
            attr = row["attribution"]
            lines.append(
                f"| {arm} | {previous['reference_arm']} | {delta['tp']} | {delta['fp']} | "
                f"{_percent(delta['precision'])} | {_percent(delta['recall'])} | "
                f"{attr['stage_restored_from_prune']['tp']}/{attr['stage_restored_from_prune']['fp']} | "
                f"{attr['stage_added_non_e0_rescue']['tp']}/{attr['stage_added_non_e0_rescue']['fp']} |"
            )
        primary = scope["arms"].get(report.get("primary_arm"))
        if primary:
            lines.extend(
                [
                    "",
                    f"Primary-arm no-path specificity (auto/review): {_percent(primary['auto']['no_path_specificity'])}/{_percent(primary['review_inclusive']['no_path_specificity'])}.",
                ]
            )
    lines.extend(["", "## Interpretation notes", ""])
    lines.extend(f"- {note}" for note in report.get("notes", []))
    return "\n".join(lines) + "\n"


def write_json(report: Mapping[str, Any], path: str | Path) -> Path:
    if report.get("schema_version") != REPORT_SCHEMA_VERSION:
        raise EvaluationError("Unsupported pruning evaluation report schema")
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return output_path


def write_markdown(report: Mapping[str, Any], path: str | Path) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_markdown(report), encoding="utf-8")
    return output_path


__all__ = [
    "BASELINE_ARM",
    "CUMULATIVE_ARMS",
    "DEFAULT_SAFETY_THRESHOLDS",
    "EvaluationError",
    "FIXED_ARM_NAMES",
    "PRIMARY_ARM",
    "STANDALONE_ARMS",
    "evaluate_directory",
    "render_markdown",
    "write_json",
    "write_markdown",
]

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

from .states import RuleVote, State


Candidate = Mapping[str, Any]
RuleEvaluator = Callable[[Candidate, State], RuleVote]


def _features(candidate: Candidate) -> Mapping[str, Any]:
    features = candidate.get("source_features")
    if not isinstance(features, Mapping):
        raise ValueError("candidate.source_features must be an object")
    return features


def _module(candidate: Candidate, name: str) -> Mapping[str, Any]:
    value = _features(candidate).get(name)
    if not isinstance(value, Mapping):
        raise ValueError(f"candidate.source_features.{name} must be an object")
    return value


def _is_e0(candidate: Candidate) -> bool:
    value = candidate.get("baseline_selected")
    if not isinstance(value, bool):
        raise ValueError("candidate.baseline_selected must be boolean")
    return value


def _review_tier(candidate: Candidate) -> str | None:
    value = candidate.get("review_tier")
    return value if isinstance(value, str) and value else None


def _and_nullable(*values: bool | None) -> bool | None:
    """Kleene AND: false dominates; otherwise unknown propagates."""

    if any(value is False for value in values):
        return False
    if any(value is None for value in values):
        return None
    return True


def _or_nullable(*values: bool | None) -> bool | None:
    """Kleene OR: true dominates; otherwise unknown propagates."""

    if any(value is True for value in values):
        return True
    if any(value is None for value in values):
        return None
    return False


def q_weak_raw(candidate: Candidate, current: State) -> RuleVote:
    d_state = _module(candidate, "D").get("state")
    b = _module(candidate, "B_STRICT")
    absolute = b.get("absolute_strength_axis")
    relative = b.get("relative_strength_axis")
    applicable = _is_e0(candidate) and d_state == "D_PASS"
    if not applicable:
        matched: bool | None = False
        reason = "NOT_E0_D_PASS"
    else:
        # The rule is specifically "both axes are false".  A missing axis is
        # indeterminate and can never be converted into PRUNE.
        matched = _and_nullable(
            False if absolute is True else (True if absolute is False else None),
            False if relative is True else (True if relative is False else None),
        )
        reason = (
            "E0_D_PASS_BOTH_RAW_AXES_FALSE"
            if matched is True
            else "RAW_AXIS_PRESENT"
            if matched is False
            else "RAW_AXIS_UNKNOWN"
        )
    return RuleVote(
        rule_id="Q_WEAK_RAW",
        matched=matched,
        target_state=State.PRUNE if matched is True else None,
        reason_code=reason,
        facts={
            "baseline_selected": _is_e0(candidate),
            "D_state": d_state,
            "B_absolute_strength_axis": absolute,
            "B_relative_strength_axis": relative,
            "current_state": current.value,
        },
    )


def q_direct_restore(candidate: Candidate, current: State) -> RuleVote:
    grade = _module(candidate, "C").get("grade")
    matched = _is_e0(candidate) and current == State.PRUNE and grade in {"C2", "C3"}
    return RuleVote(
        rule_id="Q_DIRECT_RESTORE",
        matched=matched,
        target_state=State.KEEP if matched else None,
        reason_code="DIRECT_EVIDENCE_RESTORE" if matched else "RESTORE_NOT_APPLICABLE",
        facts={
            "baseline_selected": _is_e0(candidate),
            "current_state": current.value,
            "C_grade": grade,
        },
    )


def q_context_review(candidate: Candidate, current: State) -> RuleVote:
    grade = _module(candidate, "C").get("grade")
    matched = _is_e0(candidate) and current == State.PRUNE and grade == "C1"
    return RuleVote(
        rule_id="Q_CONTEXT_REVIEW",
        matched=matched,
        target_state=State.MANUAL if matched else None,
        reason_code="CONTEXT_EVIDENCE_REVIEW" if matched else "REVIEW_NOT_APPLICABLE",
        facts={
            "baseline_selected": _is_e0(candidate),
            "current_state": current.value,
            "C_grade": grade,
        },
    )


def q_guard(candidate: Candidate, current: State) -> RuleVote:
    d_state = _module(candidate, "D").get("state")
    grade = _module(candidate, "C").get("grade")
    matched = _is_e0(candidate) and d_state == "D_GUARDED" and grade != "C3"
    return RuleVote(
        rule_id="Q_GUARD",
        matched=matched,
        target_state=State.MANUAL if matched else None,
        reason_code="GUARDED_NON_C3_REVIEW" if matched else "GUARD_NOT_APPLICABLE",
        facts={
            "baseline_selected": _is_e0(candidate),
            "current_state": current.value,
            "D_state": d_state,
            "C_grade": grade,
        },
    )


def d_block_diagnostic(candidate: Candidate, current: State) -> RuleVote:
    d_state = _module(candidate, "D").get("state")
    matched = _is_e0(candidate) and d_state == "D_BLOCK"
    return RuleVote(
        rule_id="D_BLOCK_DIAGNOSTIC",
        matched=matched,
        target_state=State.MANUAL if matched else None,
        reason_code="D_BLOCK_DIAGNOSTIC_REVIEW" if matched else "D_BLOCK_NOT_APPLICABLE",
        facts={
            "baseline_selected": _is_e0(candidate),
            "current_state": current.value,
            "D_state": d_state,
        },
    )


def r_direct(candidate: Candidate, current: State) -> RuleVote:
    d_state = _module(candidate, "D").get("state")
    grade = _module(candidate, "C").get("grade")
    matched = (not _is_e0(candidate)) and grade == "C3" and d_state != "D_BLOCK"
    return RuleVote(
        rule_id="R_DIRECT",
        matched=matched,
        target_state=State.RESCUE if matched else None,
        reason_code="NON_E0_C3_NONBLOCK_RESCUE" if matched else "DIRECT_RESCUE_NOT_APPLICABLE",
        facts={
            "baseline_selected": _is_e0(candidate),
            "current_state": current.value,
            "D_state": d_state,
            "C_grade": grade,
        },
    )


def r_convergent(candidate: Candidate, current: State) -> RuleVote:
    d_state = _module(candidate, "D").get("state")
    b_positive = _module(candidate, "B_STRICT").get("positive")
    grade = _module(candidate, "C").get("grade")
    categorical_gate = (not _is_e0(candidate)) and d_state == "D_PASS" and grade == "C2"
    matched = _and_nullable(True, b_positive) if categorical_gate else False
    return RuleVote(
        rule_id="R_CONVERGENT",
        matched=matched,
        target_state=State.RESCUE if matched is True else None,
        reason_code=(
            "NON_E0_PASS_BSTRICT_C2_RESCUE"
            if matched is True
            else "B_STRICT_UNKNOWN"
            if matched is None
            else "CONVERGENT_RESCUE_NOT_APPLICABLE"
        ),
        facts={
            "baseline_selected": _is_e0(candidate),
            "current_state": current.value,
            "D_state": d_state,
            "B_STRICT_positive": b_positive,
            "C_grade": grade,
        },
    )


def r_high_strict(candidate: Candidate, current: State) -> RuleVote:
    d_state = _module(candidate, "D").get("state")
    b_positive = _module(candidate, "B_STRICT").get("positive")
    grade = _module(candidate, "C").get("grade")
    tier = _review_tier(candidate)
    categorical_gate = (
        (not _is_e0(candidate))
        and tier == "review_high_priority"
        and grade == "C0"
        and d_state != "D_BLOCK"
    )
    matched = _and_nullable(True, b_positive) if categorical_gate else False
    return RuleVote(
        rule_id="R_HIGH_STRICT",
        matched=matched,
        target_state=State.RESCUE if matched is True else None,
        reason_code=(
            "HIGH_PRIORITY_STRICT_RESCUE"
            if matched is True
            else "B_STRICT_UNKNOWN"
            if matched is None
            else "HIGH_STRICT_NOT_APPLICABLE"
        ),
        facts={
            "baseline_selected": _is_e0(candidate),
            "current_state": current.value,
            "review_tier": tier,
            "D_state": d_state,
            "B_STRICT_positive": b_positive,
            "C_grade": grade,
        },
    )


def r_high_any(candidate: Candidate, current: State) -> RuleVote:
    d_state = _module(candidate, "D").get("state")
    b = _module(candidate, "B_STRICT")
    absolute = b.get("absolute_strength_axis")
    relative = b.get("relative_strength_axis")
    grade = _module(candidate, "C").get("grade")
    tier = _review_tier(candidate)
    categorical_gate = (
        (not _is_e0(candidate))
        and tier == "review_high_priority"
        and grade == "C0"
        and d_state != "D_BLOCK"
    )
    matched = _or_nullable(absolute, relative) if categorical_gate else False
    return RuleVote(
        rule_id="R_HIGH_ANY",
        matched=matched,
        target_state=State.RESCUE if matched is True else None,
        reason_code=(
            "HIGH_PRIORITY_ANY_RAW_AXIS_RESCUE"
            if matched is True
            else "RAW_AXIS_UNKNOWN"
            if matched is None
            else "HIGH_ANY_NOT_APPLICABLE"
        ),
        facts={
            "baseline_selected": _is_e0(candidate),
            "current_state": current.value,
            "review_tier": tier,
            "D_state": d_state,
            "B_absolute_strength_axis": absolute,
            "B_relative_strength_axis": relative,
            "C_grade": grade,
        },
    )


def r_context_a(candidate: Candidate, current: State) -> RuleVote:
    d_state = _module(candidate, "D").get("state")
    b_positive = _module(candidate, "B_STRICT").get("positive")
    a_positive = _module(candidate, "E").get("A_positive")
    tier = _review_tier(candidate)
    categorical_gate = (
        (not _is_e0(candidate))
        and tier == "review_context_needed"
        and d_state == "D_PASS"
    )
    matched = _and_nullable(a_positive, b_positive) if categorical_gate else False
    return RuleVote(
        rule_id="R_CONTEXT_A",
        matched=matched,
        target_state=State.RESCUE if matched is True else None,
        reason_code=(
            "CONTEXT_A_STRICT_PASS_RESCUE"
            if matched is True
            else "CONTEXT_A_OR_B_UNKNOWN"
            if matched is None
            else "CONTEXT_A_NOT_APPLICABLE"
        ),
        facts={
            "baseline_selected": _is_e0(candidate),
            "current_state": current.value,
            "review_tier": tier,
            "D_state": d_state,
            "E_A_positive": a_positive,
            "B_STRICT_positive": b_positive,
        },
    )


RULE_EVALUATORS: dict[str, RuleEvaluator] = {
    "Q_WEAK_RAW": q_weak_raw,
    "Q_DIRECT_RESTORE": q_direct_restore,
    "Q_CONTEXT_REVIEW": q_context_review,
    "Q_GUARD": q_guard,
    "D_BLOCK_DIAGNOSTIC": d_block_diagnostic,
    "R_DIRECT": r_direct,
    "R_CONVERGENT": r_convergent,
    "R_HIGH_STRICT": r_high_strict,
    "R_HIGH_ANY": r_high_any,
    "R_CONTEXT_A": r_context_a,
}


def evaluate_rule(rule_id: str, candidate: Candidate, current: State) -> RuleVote:
    try:
        evaluator = RULE_EVALUATORS[rule_id]
    except KeyError as exc:
        raise ValueError(f"Unknown pruning rule: {rule_id}") from exc
    return evaluator(candidate, current)


__all__ = ["RULE_EVALUATORS", "evaluate_rule"]

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any


class State(str, Enum):
    KEEP = "KEEP"
    MANUAL = "MANUAL"
    PRUNE = "PRUNE"
    RESCUE = "RESCUE"
    NO_RESCUE = "NO_RESCUE"


E0_STATES = frozenset({State.KEEP, State.MANUAL, State.PRUNE})
NON_E0_STATES = frozenset({State.RESCUE, State.MANUAL, State.NO_RESCUE})
AUTO_POSITIVE_STATES = frozenset({State.KEEP, State.RESCUE})
REVIEW_INCLUSIVE_STATES = frozenset(
    {State.KEEP, State.RESCUE, State.MANUAL}
)


class StateContractError(ValueError):
    """Raised when a state or transition violates the frozen state domains."""


def initial_state(baseline_selected: bool) -> State:
    return State.KEEP if baseline_selected else State.NO_RESCUE


def validate_state(state: State, *, baseline_selected: bool) -> None:
    allowed = E0_STATES if baseline_selected else NON_E0_STATES
    if state not in allowed:
        scope = "E0" if baseline_selected else "non-E0"
        raise StateContractError(f"{state.value} is invalid for {scope}")


def transition_action(previous: State, current: State) -> str:
    if previous == current:
        return "NO_CHANGE"
    actions = {
        (State.KEEP, State.PRUNE): "PRUNE",
        (State.KEEP, State.MANUAL): "FLAG_MANUAL",
        (State.PRUNE, State.KEEP): "RESTORE_KEEP",
        (State.PRUNE, State.MANUAL): "RESTORE_MANUAL",
        (State.MANUAL, State.KEEP): "CONFIRM_KEEP",
        (State.NO_RESCUE, State.RESCUE): "RESCUE",
        (State.NO_RESCUE, State.MANUAL): "FLAG_MANUAL",
        (State.MANUAL, State.RESCUE): "RESCUE_FROM_MANUAL",
    }
    try:
        return actions[(previous, current)]
    except KeyError as exc:
        raise StateContractError(
            f"Unsupported transition {previous.value}->{current.value}"
        ) from exc


@dataclass(frozen=True)
class RuleVote:
    rule_id: str
    matched: bool | None
    target_state: State | None
    reason_code: str
    facts: dict[str, Any]

    def __post_init__(self) -> None:
        if self.matched is not True and self.matched is not False and self.matched is not None:
            raise StateContractError(
                f"Rule {self.rule_id} matched must be true, false, or null"
            )
        if self.matched is True and self.target_state is None:
            raise StateContractError(
                f"Matched rule {self.rule_id} must provide a target state"
            )
        if self.matched is not True and self.target_state is not None:
            raise StateContractError(
                f"Unmatched/unknown rule {self.rule_id} cannot provide a target state"
            )


def apply_vote(
    previous: State,
    vote: RuleVote,
    *,
    baseline_selected: bool,
) -> tuple[State, dict[str, Any]]:
    if vote.matched is True:
        assert vote.target_state is not None
        current = vote.target_state
    elif vote.matched is None:
        # An indeterminate pruning predicate must never become PRUNE.  A
        # potentially applicable non-E0 rescue predicate is reviewable rather
        # than silently becoming NO_RESCUE.
        current = (
            State.MANUAL
            if previous in {State.NO_RESCUE, State.PRUNE}
            else previous
        )
    else:
        current = previous
    validate_state(current, baseline_selected=baseline_selected)
    action = transition_action(previous, current)
    return current, {
        "rule_id": vote.rule_id,
        "matched": vote.matched,
        "from_state": previous.value,
        "to_state": current.value,
        "action": action,
        "reason_code": vote.reason_code,
        "facts": vote.facts,
    }

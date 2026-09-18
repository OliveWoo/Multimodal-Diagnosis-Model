"""Frozen arm names and three-valued logic for the post-L5 experiment."""

from __future__ import annotations

from typing import Iterable, Optional

Tri = Optional[bool]

FAMILIES = ("CORR", "RAW")
SCOPES = ("RSC", "P_RST", "M_RST")
PREDICATES = ("A", "B", "C", "AB", "AC", "BC", "ABC", "2OF3", "OR")
CUMULATIVE_PREDICATES = ("ABC", "2OF3", "OR")
BASE_ARM = "L5_FROZEN"
FUTURE_HOLDOUT_CANDIDATE_ARM = "K_CORR_2_2OF3"

E0_STATES = {"KEEP", "MANUAL", "PRUNE"}
NON_E0_STATES = {"RESCUE", "MANUAL", "NO_RESCUE"}
D_STATES = {"D_PASS", "D_GUARDED", "D_BLOCK", "D_UNKNOWN"}
C_GRADES = {"C0", "C1", "C2", "C3", "CNEG"}


def tri_and(values: Iterable[Tri]) -> Tri:
    items = tuple(values)
    if any(value is False for value in items):
        return False
    if all(value is True for value in items):
        return True
    return None


def tri_or(values: Iterable[Tri]) -> Tri:
    items = tuple(values)
    if any(value is True for value in items):
        return True
    if all(value is False for value in items):
        return False
    return None


def tri_two_of_three(a: Tri, b: Tri, c: Tri) -> Tri:
    values = (a, b, c)
    true_count = sum(value is True for value in values)
    unknown_count = sum(value is None for value in values)
    if true_count >= 2:
        return True
    if true_count + unknown_count < 2:
        return False
    return None


def predicate_values(a: Tri, b: Tri, c: Tri) -> dict[str, Tri]:
    return {
        "A": a,
        "B": b,
        "C": c,
        "AB": tri_and((a, b)),
        "AC": tri_and((a, c)),
        "BC": tri_and((b, c)),
        "ABC": tri_and((a, b, c)),
        "2OF3": tri_two_of_three(a, b, c),
        "OR": tri_or((a, b, c)),
    }


def standalone_arm(family: str, scope: str, predicate: str) -> str:
    return f"S_{family}_{scope}_{predicate}"


def ungated_arm(family: str, predicate: str) -> str:
    return f"X_{family}_RSC_UNGATED_{predicate}"


def cumulative_arm(family: str, index: int, predicate: str) -> str:
    suffix = "OR_STRESS" if predicate == "OR" else predicate
    return f"K_{family}_{index}_{suffix}"


def arm_order() -> tuple[str, ...]:
    arms: list[str] = [BASE_ARM]
    for family in FAMILIES:
        for scope in SCOPES:
            for predicate in PREDICATES:
                arms.append(standalone_arm(family, scope, predicate))
    for family in FAMILIES:
        for predicate in PREDICATES:
            arms.append(ungated_arm(family, predicate))
    for family in FAMILIES:
        for index, predicate in enumerate(CUMULATIVE_PREDICATES, start=1):
            arms.append(cumulative_arm(family, index, predicate))
    return tuple(arms)


ARM_ORDER = arm_order()


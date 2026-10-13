"""Shared calendar-date event alignment for clinical evidence.

The primary retrospective episode window includes the two calendar dates
before and after the mNGS collection date. Exact hour differences remain
available for audit and sensitivity analyses, but do not control inclusion.
"""

from __future__ import annotations

from datetime import datetime
from typing import Iterable


DEFAULT_CALENDAR_DAYS_EACH_SIDE = 2
PRIMARY_WINDOW_MODE = "inclusive_calendar_dates"


def signed_calendar_day_delta(observed: datetime, index: datetime) -> int:
    return (observed.date() - index.date()).days


def within_calendar_day_window(
    observed: datetime,
    index: datetime,
    *,
    days_each_side: int = DEFAULT_CALENDAR_DAYS_EACH_SIDE,
) -> bool:
    return abs(signed_calendar_day_delta(observed, index)) <= days_each_side


def nearest_calendar_event(
    observed: datetime,
    events: Iterable[datetime],
    *,
    days_each_side: int = DEFAULT_CALENDAR_DAYS_EACH_SIDE,
) -> dict[str, object] | None:
    candidates = []
    for event in events:
        calendar_delta = signed_calendar_day_delta(observed, event)
        hour_delta = (observed - event).total_seconds() / 3600.0
        candidates.append((abs(calendar_delta), abs(hour_delta), calendar_delta, hour_delta, event))
    if not candidates:
        return None
    _, _, calendar_delta, hour_delta, event = min(candidates)
    return {
        "event": event,
        "within_window": abs(calendar_delta) <= days_each_side,
        "signed_calendar_day_delta": calendar_delta,
        "signed_delta_hours": hour_delta,
        "window_mode": PRIMARY_WINDOW_MODE,
        "calendar_days_each_side": days_each_side,
    }

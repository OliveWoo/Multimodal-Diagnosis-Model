from __future__ import annotations

import re
from datetime import date, datetime, timedelta
from typing import Any, Iterable, Mapping, Sequence


DEFAULT_WINDOW_BEFORE_DAYS = 2
DEFAULT_WINDOW_AFTER_DAYS = 2

ENTRY_DATE_KEYS = ("collected_time", "reported_time")
MNGS_DATE_KEY = "collected_time"

DATE_PATTERN = re.compile(r"(?P<year>20\d{2}|\d{2})[/-](?P<month>\d{1,2})[/-](?P<day>\d{1,2})")
COMPACT_DATE_PATTERN = re.compile(r"(?P<year>20\d{2})(?P<month>\d{2})(?P<day>\d{2})")
SPECIMEN_PREFIX_DATE_PATTERN = re.compile(r"^(?P<year>\d{2})(?P<month>\d{2})(?P<day>\d{2})")
SPECIMEN_MONTH_NAME_PATTERN = re.compile(r"^(?P<year>\d{2})(?P<month>[A-Za-z]{3})(?P<day>\d{2})")
EMPTY_DATE_VALUES = {"", "-", "nan", "none", "nat"}
MONTH_NAME_TO_NUMBER = {
    "jan": 1,
    "feb": 2,
    "mar": 3,
    "apr": 4,
    "may": 5,
    "jun": 6,
    "jul": 7,
    "aug": 8,
    "sep": 9,
    "oct": 10,
    "nov": 11,
    "dec": 12,
}


def parse_date(value: Any) -> date | None:
    if value is None:
        return None
    if isinstance(value, datetime):
        return value.date()
    if isinstance(value, date):
        return value

    to_pydatetime = getattr(value, "to_pydatetime", None)
    if callable(to_pydatetime):
        try:
            return to_pydatetime().date()
        except (TypeError, ValueError):
            return None

    text = str(value).strip()
    if text.lower() in EMPTY_DATE_VALUES:
        return None

    match = DATE_PATTERN.search(text)
    if not match:
        match = COMPACT_DATE_PATTERN.search(text)
    if not match:
        return None

    year = int(match.group("year"))
    if year < 100:
        year += 2000
    try:
        return date(year, int(match.group("month")), int(match.group("day")))
    except ValueError:
        return None


def mngs_collected_dates(mngs_entries: Iterable[Mapping[str, Any]] | None) -> list[date]:
    dates: list[date] = []
    for entry in mngs_entries or []:
        if not isinstance(entry, Mapping):
            continue
        parsed = parse_date(entry.get(MNGS_DATE_KEY))
        if parsed is not None:
            dates.append(parsed)
    return dates


def cutoff_entries_from_dates(values: Iterable[Any] | None) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    seen: set[date] = set()
    for value in values or []:
        parsed = parse_date(value)
        if parsed is None or parsed in seen:
            continue
        seen.add(parsed)
        entries.append({MNGS_DATE_KEY: parsed.isoformat()})
    return entries


def infer_date_from_specimen_id(value: Any) -> date | None:
    text = str(value or "").strip()
    if not text:
        return None

    match = SPECIMEN_PREFIX_DATE_PATTERN.match(text)
    if match:
        year = 2000 + int(match.group("year"))
        try:
            return date(year, int(match.group("month")), int(match.group("day")))
        except ValueError:
            pass

    match = SPECIMEN_MONTH_NAME_PATTERN.match(text)
    if not match:
        return None
    month = MONTH_NAME_TO_NUMBER.get(match.group("month").lower())
    if month is None:
        return None
    year = 2000 + int(match.group("year"))
    try:
        return date(year, month, int(match.group("day")))
    except ValueError:
        return None


def cutoff_entries_from_manifest(manifest: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not isinstance(manifest, Mapping):
        return []

    manifest_dates = manifest.get("mngs_collected_dates")
    if not manifest_dates:
        manifest_dates = (manifest.get("identifiers") or {}).get("mngs_collected_dates")
    if isinstance(manifest_dates, list):
        entries = cutoff_entries_from_dates(manifest_dates)
        if entries:
            return entries

    specimen_id = (manifest.get("identifiers") or {}).get("specimen_id")
    inferred = infer_date_from_specimen_id(specimen_id)
    if inferred is None:
        return []
    return [{MNGS_DATE_KEY: inferred.isoformat()}]


def entry_observation_date(entry: Mapping[str, Any]) -> date | None:
    for key in ENTRY_DATE_KEYS:
        parsed = parse_date(entry.get(key))
        if parsed is not None:
            return parsed
    return None


def date_in_mngs_window(
    observed_date: date,
    cutoff_dates: Sequence[date],
    *,
    window_before_days: int = DEFAULT_WINDOW_BEFORE_DAYS,
    window_after_days: int = DEFAULT_WINDOW_AFTER_DAYS,
) -> bool:
    for cutoff_date in cutoff_dates:
        start_date = cutoff_date - timedelta(days=window_before_days)
        end_date = cutoff_date + timedelta(days=window_after_days)
        if start_date <= observed_date <= end_date:
            return True
    return False


def filter_entries_by_mngs_window(
    entries: Iterable[Mapping[str, Any]] | None,
    mngs_entries: Iterable[Mapping[str, Any]] | None,
    *,
    window_before_days: int = DEFAULT_WINDOW_BEFORE_DAYS,
    window_after_days: int = DEFAULT_WINDOW_AFTER_DAYS,
) -> list[Mapping[str, Any]]:
    entry_list = list(entries or [])
    cutoff_dates = mngs_collected_dates(mngs_entries)
    if not cutoff_dates:
        return entry_list

    filtered: list[Mapping[str, Any]] = []
    for entry in entry_list:
        if not isinstance(entry, Mapping):
            continue
        observed_date = entry_observation_date(entry)
        if observed_date is None:
            continue
        if date_in_mngs_window(
            observed_date,
            cutoff_dates,
            window_before_days=window_before_days,
            window_after_days=window_after_days,
        ):
            filtered.append(entry)
    return filtered

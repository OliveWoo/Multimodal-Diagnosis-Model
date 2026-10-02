from __future__ import annotations

import json
from datetime import date, datetime, time
from pathlib import Path
from typing import Any

import pandas as pd

from utils import sanitize_filename

SHEET_NAME = "gm test"


def _ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _normalize_value(value: Any) -> str:
    if value is None:
        return ""
    dt_value: datetime | None = None
    if isinstance(value, pd.Timestamp):
        if pd.isna(value):
            return ""
        dt_value = value.to_pydatetime()
    elif isinstance(value, datetime):
        dt_value = value
    elif isinstance(value, date):
        return value.strftime("%Y-%m-%d")
    elif isinstance(value, time):
        if value.second == 0 and value.microsecond == 0:
            return value.strftime("%H:%M")
        return value.strftime("%H:%M:%S")
    elif hasattr(value, "dtype") and "datetime64" in str(value.dtype):
        timestamp = pd.Timestamp(value)
        if pd.isna(timestamp):
            return ""
        dt_value = timestamp.to_pydatetime()

    if dt_value is not None:
        if (
            dt_value.hour == 0
            and dt_value.minute == 0
            and dt_value.second == 0
            and dt_value.microsecond == 0
        ):
            return dt_value.strftime("%Y-%m-%d")
        if dt_value.second == 0 and dt_value.microsecond == 0:
            return dt_value.strftime("%Y-%m-%d %H:%M")
        return dt_value.strftime("%Y-%m-%d %H:%M:%S")

    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return str(value)


def _sheet_to_rows(frame: pd.DataFrame) -> list[list[str]]:
    if frame.empty:
        return []
    rows: list[list[str]] = []
    for row in frame.itertuples(index=False, name=None):
        rows.append([_normalize_value(cell) for cell in row])
    return rows


def parse(excel_path: Path | str, output_dir: Path | str) -> Path:
    excel_path = Path(excel_path)
    output_dir = _ensure_dir(Path(output_dir))

    frame = pd.read_excel(excel_path, sheet_name=SHEET_NAME, header=None, dtype=object)
    raw_data = _sheet_to_rows(frame)

    payload = {
        "source_file": excel_path.name,
        "sheet": SHEET_NAME,
        "raw_data": raw_data,
    }

    sanitized_workbook = sanitize_filename(excel_path.stem)
    sanitized_sheet = sanitize_filename(SHEET_NAME)
    raw_filename = f"{sanitized_workbook}_{sanitized_sheet}_Raw.json"
    output_path = output_dir / raw_filename
    with output_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)

    return output_path


__all__ = ["SHEET_NAME", "parse"]

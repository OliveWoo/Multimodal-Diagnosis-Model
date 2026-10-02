from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from utils import sanitize_filename

SHEET_NAME = "admission diagnosis"


def _ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _normalize_value(value: Any) -> str:
    if value is None:
        return ""
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

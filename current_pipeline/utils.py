import re
from datetime import date, datetime
from pathlib import Path

import pandas as pd

_OVERWRITE_MODE = "ask"


def stringify_key(value):
    if isinstance(value, (datetime, date)):
        return value.strftime("%Y-%m-%d")
    return str(value)


def to_python_scalar(value):
    if pd.isna(value):
        return None
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (datetime, date)):
        if isinstance(value, date) and not isinstance(value, datetime):
            value = datetime(value.year, value.month, value.day)
        return value.strftime("%Y-%m-%dT%H:%M:%S")
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return value


def split_name_unit(raw_label):
    match = re.match(r"^(?P<name>.*?)(?:\s*\((?P<unit>[^()]+)\))?$", raw_label.strip())
    if not match:
        return raw_label.strip(), None
    name = match.group("name").strip()
    unit = match.group("unit")
    if unit:
        unit = unit.strip()
    return name, unit

def sorted_labels(labels):
    def sort_key(label):
        if isinstance(label, (datetime, date)):
            return (0, datetime(label.year, label.month, label.day))
        text = str(label)
        for fmt in ("%Y-%m-%d", "%Y/%m/%d", "%Y%m%d"):
            try:
                return (0, datetime.strptime(text, fmt))
            except ValueError:
                continue
        return (1, text)

    return sorted(labels, key=sort_key)


def sanitize_filename(name: str) -> str:
    sanitized = re.sub(r"\s+", "_", name.strip())
    return sanitized or "_"


def confirm_overwrite(path: Path, prompt: str | None = None) -> bool:
    """Prompt whether to overwrite existing files, supporting global decisions."""
    global _OVERWRITE_MODE  # noqa: PLW0603

    if not path.exists():
        return True

    if _OVERWRITE_MODE == "overwrite_all":
        return True
    if _OVERWRITE_MODE == "skip_all":
        return False

    message = prompt or f"檔案 {path} 已存在，要覆蓋嗎？"
    while True:
        choice = input(f"{message} [y/n/all/skip]: ").strip().lower()
        if choice in {"y", "yes"}:
            return True
        if choice in {"n", "no"}:
            return False
        if choice == "all":
            _OVERWRITE_MODE = "overwrite_all"
            return True
        if choice == "skip":
            _OVERWRITE_MODE = "skip_all"
            return False
        print("請輸入 y / n / all / skip。")

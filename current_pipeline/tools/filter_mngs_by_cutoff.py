from __future__ import annotations

"""
根據 cutoff JSON（patient_id -> collected_time）過濾 mNGS grouped 資料。

規則：
- 對每個病人，只保留 collected_time「日期」<= 該病人 cutoff 日期的檢體（同一天保留，時間不重要）。
- 病人沒有 cutoff 時，不會刪除該病人的任何檢體。

預設：
- 來源：outputs/specimen_lookup_summary_mNGS_grouped_with_type.json
- cutoff：outputs/patient_cutoffs_from_mngs.json
- 輸出：outputs/specimen_lookup_summary_mNGS_grouped_with_type_filtered.json
"""

import argparse
import json
from datetime import datetime, date
from pathlib import Path
from typing import Any, Mapping, Sequence

DEFAULT_SOURCE = Path("outputs") / "specimen_lookup_summary_mNGS_grouped_with_type.json"
DEFAULT_CUTOFF = Path("outputs") / "patient_cutoffs_from_mngs.json"
DEFAULT_OUTPUT = Path("outputs") / "specimen_lookup_summary_mNGS_grouped_with_type_filtered.json"

# 支援的日期/時間格式（取日期部分，不含時區）
DATETIME_FORMATS = (
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%d %H:%M",
    "%Y-%m-%d",
    "%Y/%m/%d %H:%M:%S",
    "%Y/%m/%d %H:%M",
    "%Y/%m/%d",
)


def _parse_date(value: Any) -> date | None:
    text = str(value or "").strip()
    if not text:
        return None
    for fmt in DATETIME_FORMATS:
        try:
            return datetime.strptime(text, fmt).date()
        except ValueError:
            continue
    return None


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def filter_by_cutoff(
    source: Mapping[str, list[Mapping[str, Any]]],
    cutoff_map: Mapping[str, str],
) -> dict[str, list[Mapping[str, Any]]]:
    result: dict[str, list[Mapping[str, Any]]] = {}
    for patient_id, entries in source.items():
        if not isinstance(entries, list):
            raise TypeError(f"病人 {patient_id} 的資料不是 list，而是 {type(entries)}")

        cutoff_date = _parse_date(cutoff_map.get(patient_id))
        if cutoff_date is None:
            # 無 cutoff，全部保留
            result[str(patient_id)] = entries
            continue

        kept: list[Mapping[str, Any]] = []
        for entry in entries:
            if not isinstance(entry, Mapping):
                continue
            collected_date = _parse_date(entry.get("collected_time"))
            if collected_date is None:
                # 無法解析日期，保留以免誤刪
                kept.append(entry)
                continue
            if collected_date <= cutoff_date:
                kept.append(entry)

        if kept:
            result[str(patient_id)] = kept
    return result


def write_json(payload: Any, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="依病人 cutoff 日期過濾 mNGS grouped 資料（日期<=cutoff）。")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE, help=f"來源 grouped JSON (default: {DEFAULT_SOURCE})")
    parser.add_argument("--cutoff", type=Path, default=DEFAULT_CUTOFF, help=f"cutoff JSON (default: {DEFAULT_CUTOFF})")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help=f"輸出檔案路徑 (default: {DEFAULT_OUTPUT})")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    source = _load_json(args.source)
    cutoff_map = _load_json(args.cutoff)

    if not isinstance(source, Mapping):
        raise TypeError(f"來源 JSON 應為 dict，實際為 {type(source)}")
    if not isinstance(cutoff_map, Mapping):
        raise TypeError(f"cutoff JSON 應為 dict，實際為 {type(cutoff_map)}")

    filtered = filter_by_cutoff(source, cutoff_map)
    write_json(filtered, args.output)

    kept_patients = len(filtered)
    kept_records = sum(len(v) for v in filtered.values())
    print(f"完成過濾，保留 {kept_patients} 位病人、{kept_records} 筆檢體，輸出至 {args.output}")


if __name__ == "__main__":
    main()

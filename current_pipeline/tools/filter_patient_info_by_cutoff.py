from __future__ import annotations

"""
Filter patient_info JSON files by patient cutoff dates.

Rules:
- Only `collected_time` is used as the timestamp field.
- Entries without a parseable `collected_time` are dropped.
- If a patient has a cutoff date, keep entries within [cutoff - window_before, cutoff + window_after].
- If a patient has no cutoff date, keep entries only by timestamp validity (no date-window filtering).
"""

import argparse
import json
import re
import shutil
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import Any, Mapping, Sequence

DEFAULT_PATIENT_ROOT = Path(r"C:\Users\User\Desktop\patient_info\data(all)")
DEFAULT_CUTOFF = Path("outputs") / "patient_cutoffs_from_mngs.json"

DATETIME_FORMATS = (
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%d %H:%M",
    "%Y-%m-%d",
    "%Y/%m/%d %H:%M:%S",
    "%Y/%m/%d %H:%M",
    "%Y/%m/%d",
)

TIME_KEY = "collected_time"


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


def _pick_entry_date(entry: Mapping[str, Any]) -> date | None:
    return _parse_date(entry.get(TIME_KEY))


def _load_cutoff_map(path: Path) -> dict[str, date]:
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    if not isinstance(raw, Mapping):
        raise TypeError(f"cutoff JSON should be a dict, got {type(raw)}")

    parsed: dict[str, date] = {}
    for pid, val in raw.items():
        dt = _parse_date(val)
        if dt:
            parsed[str(pid)] = dt
    return parsed


def _extract_patient_id(folder: Path) -> str:
    match = re.search(r"NGS_patient_(.+?)_json", folder.name)
    if match:
        return match.group(1)
    return folder.name


def _ensure_output_root(patient_root: Path, output_root: Path | None) -> Path:
    if output_root:
        return output_root
    return patient_root.parent / f"{patient_root.name}_filtered"


def _filter_entries(
    entries: list[Any],
    cutoff_date: date | None,
    *,
    window_before: int,
    window_after: int,
) -> list[Any]:
    kept: list[Any] = []
    if cutoff_date is None:
        start_date = None
        end_date = None
    else:
        start_date = cutoff_date - timedelta(days=window_before)
        end_date = cutoff_date + timedelta(days=window_after)

    for entry in entries:
        if not isinstance(entry, Mapping):
            continue

        entry_date = _pick_entry_date(entry)
        if entry_date is None:
            continue

        if start_date is None or end_date is None or start_date <= entry_date <= end_date:
            kept.append(entry)
    return kept


def _process_file(
    src: Path,
    dst: Path,
    cutoff_date: date | None,
    *,
    window_before: int,
    window_after: int,
) -> tuple[int, int]:
    with src.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if isinstance(data, list):
        filtered = _filter_entries(data, cutoff_date, window_before=window_before, window_after=window_after)
        original_count = len(data)
        kept_count = len(filtered)
    else:
        filtered = data
        original_count = 0
        kept_count = 0

    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("w", encoding="utf-8") as handle:
        json.dump(filtered, handle, ensure_ascii=False, indent=2)

    return original_count, kept_count


def filter_patient_folder(
    patient_folder: Path,
    output_folder: Path,
    cutoff_date: date | None,
    *,
    window_before: int,
    window_after: int,
) -> tuple[int, int]:
    total_original = 0
    total_kept = 0
    output_folder.mkdir(parents=True, exist_ok=True)

    for src_file in patient_folder.glob("*.json"):
        dst_file = output_folder / src_file.name
        orig, kept = _process_file(
            src_file,
            dst_file,
            cutoff_date,
            window_before=window_before,
            window_after=window_after,
        )
        total_original += orig
        total_kept += kept

    for src_file in patient_folder.iterdir():
        if src_file.suffix.lower() != ".json":
            dst_file = output_folder / src_file.name
            if src_file.is_file():
                shutil.copy2(src_file, dst_file)
            elif src_file.is_dir():
                shutil.copytree(src_file, dst_file, dirs_exist_ok=True)

    return total_original, total_kept


def filter_all_patients(
    patient_root: Path,
    cutoff_map: Mapping[str, date],
    output_root: Path,
    *,
    window_before: int,
    window_after: int,
) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    total_patients = 0
    processed_patients = 0
    total_records = 0
    total_kept = 0

    for folder in patient_root.iterdir():
        if not folder.is_dir():
            continue

        patient_id = _extract_patient_id(folder)
        total_patients += 1
        cutoff_date = cutoff_map.get(patient_id)

        dst_folder = output_root / folder.name
        orig, kept = filter_patient_folder(
            folder,
            dst_folder,
            cutoff_date,
            window_before=window_before,
            window_after=window_after,
        )

        total_records += orig
        total_kept += kept
        processed_patients += 1

        if cutoff_date:
            print(f"[{patient_id}] cutoff={cutoff_date} original={orig} kept={kept}")
        else:
            print(
                f"[{patient_id}] no cutoff, removed entries without valid collected_time: "
                f"original={orig}, kept={kept}"
            )

    print(
        f"Done: processed {processed_patients}/{total_patients} patients, "
        f"records {total_records} -> kept {total_kept}. Output: {output_root}"
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter patient_info JSON files by patient cutoff date. "
            "Only entries with valid collected_time are kept."
        )
    )
    parser.add_argument(
        "--patient-root",
        type=Path,
        default=DEFAULT_PATIENT_ROOT,
        help=f"Patient folder root containing NGS_patient_*_json folders (default: {DEFAULT_PATIENT_ROOT})",
    )
    parser.add_argument(
        "--cutoff",
        type=Path,
        default=DEFAULT_CUTOFF,
        help=f"Cutoff mapping JSON path (default: {DEFAULT_CUTOFF})",
    )
    parser.add_argument(
        "--window-before-days",
        type=int,
        default=2,
        help="Keep entries from this many days before cutoff (default: 2).",
    )
    parser.add_argument(
        "--window-after-days",
        type=int,
        default=2,
        help="Keep entries up to this many days after cutoff (default: 2).",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help="Output root. If omitted, defaults to <patient-root>_filtered.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    cutoff_map = _load_cutoff_map(args.cutoff)
    output_root = _ensure_output_root(args.patient_root, args.output_root)
    filter_all_patients(
        args.patient_root,
        cutoff_map,
        output_root,
        window_before=args.window_before_days,
        window_after=args.window_after_days,
    )


if __name__ == "__main__":
    main()

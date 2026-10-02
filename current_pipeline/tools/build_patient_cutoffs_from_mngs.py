from __future__ import annotations

"""
Generate per-patient cutoff timestamps using the mNGS grouped JSON.

The script reads outputs/specimen_lookup_summary_mNGS_grouped_with_type.json,
extracts collected_time for each patient, and writes a simple mapping:
{
  "<patient_id>": "YYYY-MM-DD HH:MM:SS"
}
By default it uses the earliest collected_time per patient (or the latest when
--strategy latest is provided).
"""

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, Sequence

DEFAULT_SOURCE = Path("outputs") / "specimen_lookup_summary_mNGS_grouped_with_type.json"
DEFAULT_OUTPUT = Path("outputs") / "patient_cutoffs_from_mngs.json"

DATETIME_FORMATS = (
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%d %H:%M",
    "%Y/%m/%d %H:%M:%S",
    "%Y/%m/%d %H:%M",
)


def _parse_datetime(value: Any) -> datetime | None:
    text = str(value or "").strip()
    if not text:
        return None
    for fmt in DATETIME_FORMATS:
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None


def _load_grouped(path: Path) -> Mapping[str, list[Mapping[str, Any]]]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, Mapping):
        raise TypeError(f"Expected grouped JSON to be a mapping, got {type(data)}")
    return data  # type: ignore[return-value]


def build_cutoffs(
    grouped_payload: Mapping[str, list[Mapping[str, Any]]],
    *,
    strategy: str = "earliest",
) -> dict[str, str]:
    use_latest = strategy == "latest"
    cutoffs: dict[str, str] = {}

    for patient_id, entries in grouped_payload.items():
        if not isinstance(entries, list):
            raise TypeError(f"Expected a list of entries for patient {patient_id}, got {type(entries)}")

        timestamps: list[datetime] = []
        for entry in entries:
            if not isinstance(entry, Mapping):
                continue
            dt = _parse_datetime(entry.get("collected_time"))
            if dt:
                timestamps.append(dt)

        if not timestamps:
            continue

        chosen = max(timestamps) if use_latest else min(timestamps)
        cutoffs[str(patient_id)] = chosen.strftime("%Y-%m-%d %H:%M:%S")

    return cutoffs


def write_cutoffs(cutoffs: Mapping[str, str], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(cutoffs, handle, ensure_ascii=False, indent=2)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build patient cutoffs from specimen_lookup_summary_mNGS_grouped_with_type.json.",
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=DEFAULT_SOURCE,
        help=f"Path to specimen_lookup_summary_mNGS_grouped_with_type.json (default: {DEFAULT_SOURCE})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Where to write the cutoff mapping JSON (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--strategy",
        choices=["earliest", "latest"],
        default="earliest",
        help="Pick earliest or latest collected_time per patient (default: earliest).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    payload = _load_grouped(args.source)
    cutoffs = build_cutoffs(payload, strategy=args.strategy)
    write_cutoffs(cutoffs, args.output)

    print(f"Wrote {len(cutoffs)} patient cutoffs to {args.output}")
    sample = list(cutoffs.items())[:5]
    if sample:
        print("Sample:", sample)


if __name__ == "__main__":
    main()

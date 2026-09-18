#!/usr/bin/env python3
"""Build final patient JSON folders for the mNGS-window workflow.

For each patient and each original test-type JSON file, keep the union of:

1. All records whose selected clinical time falls inside the mNGS anchor window.
2. All-time microbiology records that match candidate microbes found in that
   window.

No metadata files are written to the output folder. Output files keep the
original JSON shape as much as possible and only remove non-matching records.
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any

from build_candidate_window_evidence import (
    canonical_key,
    choose_record_time,
    gm_related_to_candidate,
    names_match,
    record_type_from_path,
)


PATIENT_DIR_RE = re.compile(r"NGS_patient_(\d+)_json$")
MICROBE_RECORD_TYPES = {
    "culture",
    "filmarray",
    "molecular_microbiology",
    "mngs_grouped",
    "mNGS_grouped",
    "mngs_name_reads",
    "all_RK_NTC_microbes",
    "gm_test",
}
NON_TEST_RECORD_TYPES = {"source_mapping"}
ALWAYS_KEEP_RECORD_TYPES = {"admission_diagnosis", "underlying"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Keep all records inside the mNGS anchor window, plus all-time "
            "candidate-matching microbiology evidence."
        )
    )
    parser.add_argument(
        "--input-root",
        required=True,
        type=Path,
        help="Root containing original NGS_patient_*_json directories.",
    )
    parser.add_argument(
        "--candidate-root",
        required=True,
        type=Path,
        help="Root containing mNGS-anchor candidate evidence bundles.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        help="Output root for filtered NGS_patient_*_json directories.",
    )
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output.")
    parser.add_argument(
        "--write-empty-files",
        action="store_true",
        help="Write [] for original test-type files with no retained records.",
    )
    return parser.parse_args()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def write_json(path: Path, data: Any, pretty: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        if pretty:
            json.dump(data, f, ensure_ascii=False, indent=2)
            f.write("\n")
        else:
            json.dump(data, f, ensure_ascii=False, separators=(",", ":"))
            f.write("\n")


def patient_id_from_dir(path: Path) -> int | None:
    match = PATIENT_DIR_RE.match(path.name)
    return int(match.group(1)) if match else None


def candidate_file(candidate_root: Path, patient_id: int) -> Path | None:
    matches = sorted(
        candidate_root.glob(f"NGS_patient_{patient_id}_candidate_window_*_evidence.json")
    )
    return matches[0] if matches else None


def parse_iso_time(value: Any) -> datetime | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        return None


def candidate_names(bundle: dict[str, Any]) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    for candidate in bundle.get("candidates", []):
        values = [candidate.get("name"), *candidate.get("observed_names", [])]
        for value in values:
            if not isinstance(value, str) or not value.strip():
                continue
            key = canonical_key(value)
            if key in seen:
                continue
            seen.add(key)
            names.append(value)
    return names


def inside_window(record: dict[str, Any], start: datetime | None, end: datetime | None) -> bool:
    if start is None or end is None:
        return False
    record_time, _ = choose_record_time(record)
    return record_time is not None and start <= record_time <= end


def matches_any_candidate(name: Any, candidates: list[str]) -> bool:
    if not isinstance(name, str) or not name.strip():
        return False
    return any(names_match(name, candidate) for candidate in candidates)


def candidate_micro_filter(
    record_type: str, record: dict[str, Any], candidates: list[str]
) -> dict[str, Any] | None:
    if record_type == "culture":
        return record if matches_any_candidate(record.get("organism"), candidates) else None

    if record_type in {"filmarray", "molecular_microbiology"}:
        return record if matches_any_candidate(record.get("target"), candidates) else None

    if record_type == "gm_test":
        return record if any(gm_related_to_candidate(candidate) for candidate in candidates) else None

    if record_type == "mngs_name_reads":
        return record if matches_any_candidate(record.get("name"), candidates) else None

    if record_type in {"mngs_grouped", "mNGS_grouped", "all_RK_NTC_microbes"}:
        pathogens = record.get("pathogens")
        if not isinstance(pathogens, dict):
            return None
        filtered_pathogens: dict[str, list[dict[str, Any]]] = {}
        kept_count = 0
        for group, items in pathogens.items():
            kept_items = []
            if isinstance(items, list):
                for item in items:
                    if isinstance(item, dict) and matches_any_candidate(item.get("name"), candidates):
                        kept_items.append(item)
            filtered_pathogens[group] = kept_items
            kept_count += len(kept_items)
        if kept_count == 0:
            return None
        copied = dict(record)
        copied["pathogens"] = filtered_pathogens
        return copied

    return None


def keep_record(
    record_type: str,
    record: dict[str, Any],
    candidates: list[str],
    start: datetime | None,
    end: datetime | None,
) -> dict[str, Any] | None:
    if record_type in ALWAYS_KEEP_RECORD_TYPES:
        return record
    if inside_window(record, start, end):
        return record
    if record_type in MICROBE_RECORD_TYPES:
        return candidate_micro_filter(record_type, record, candidates)
    return None


def filter_json_data(
    record_type: str,
    data: Any,
    candidates: list[str],
    start: datetime | None,
    end: datetime | None,
) -> Any | None:
    if isinstance(data, list):
        kept = [
            filtered
            for item in data
            if isinstance(item, dict)
            for filtered in [keep_record(record_type, item, candidates, start, end)]
            if filtered is not None
        ]
        return kept if kept else None

    if isinstance(data, dict):
        return keep_record(record_type, data, candidates, start, end)

    return None


def main() -> int:
    args = parse_args()
    input_root = args.input_root.resolve()
    candidate_root = args.candidate_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    patients = 0
    patient_dirs_written = 0
    files_written = 0
    records_written = 0

    for patient_dir in sorted(input_root.iterdir(), key=lambda p: p.name):
        if not patient_dir.is_dir():
            continue
        patient_id = patient_id_from_dir(patient_dir)
        if patient_id is None:
            continue
        cfile = candidate_file(candidate_root, patient_id)
        if cfile is None:
            continue
        bundle = load_json(cfile)
        start = parse_iso_time(bundle.get("candidate_window", {}).get("start"))
        end = parse_iso_time(bundle.get("candidate_window", {}).get("end"))
        candidates = candidate_names(bundle)
        patient_file_count = 0
        patients += 1

        for source_file in sorted(patient_dir.glob("*.json")):
            record_type = record_type_from_path(source_file, patient_id)
            if record_type in NON_TEST_RECORD_TYPES:
                continue
            filtered = filter_json_data(
                record_type, load_json(source_file), candidates, start, end
            )
            if filtered is None:
                if args.write_empty_files:
                    write_json(output_root / patient_dir.name / source_file.name, [], args.pretty)
                    files_written += 1
                    patient_file_count += 1
                continue
            count = len(filtered) if isinstance(filtered, list) else 1
            if count == 0:
                continue
            write_json(output_root / patient_dir.name / source_file.name, filtered, args.pretty)
            files_written += 1
            patient_file_count += 1
            records_written += count

        if patient_file_count:
            patient_dirs_written += 1

    print(
        "patients={patients} patient_dirs_written={dirs} files_written={files} records_written={records} output={output}".format(
            patients=patients,
            dirs=patient_dirs_written,
            files=files_written,
            records=records_written,
            output=output_root,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

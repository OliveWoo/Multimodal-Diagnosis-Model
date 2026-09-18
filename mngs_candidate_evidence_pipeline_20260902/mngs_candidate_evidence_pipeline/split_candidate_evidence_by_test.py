#!/usr/bin/env python3
"""Write candidate-matching evidence as separate files by test type.

This consumes candidate-window evidence bundles and the original patient JSON
directories. Output files preserve the original JSON shape as much as possible:
only non-matching records are removed, and no metadata wrapper is added.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from build_candidate_window_evidence import (
    canonical_key,
    extract_mngs_pathogens,
    gm_related_to_candidate,
    names_match,
    normalize_name,
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split candidate evidence into separate test-type JSON files."
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
        help="Root containing NGS_patient_*_candidate_window_*_evidence.json files.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        help="Output root for filtered NGS_patient_*_json directories.",
    )
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON output.")
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


def candidate_names(bundle: dict[str, Any]) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    for candidate in bundle.get("candidates", []):
        for value in [candidate.get("name"), *candidate.get("observed_names", [])]:
            if not isinstance(value, str) or not value.strip():
                continue
            key = canonical_key(value)
            if key in seen:
                continue
            seen.add(key)
            names.append(value)
    return names


def matches_any_candidate(name: Any, candidates: list[str]) -> bool:
    if not isinstance(name, str) or not name.strip():
        return False
    return any(names_match(name, candidate) for candidate in candidates)


def filter_record(record_type: str, record: dict[str, Any], candidates: list[str]) -> dict[str, Any] | None:
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


def filter_json_data(record_type: str, data: Any, candidates: list[str]) -> Any | None:
    if isinstance(data, list):
        kept = [
            filtered
            for item in data
            if isinstance(item, dict)
            for filtered in [filter_record(record_type, item, candidates)]
            if filtered is not None
        ]
        return kept if kept else None

    if isinstance(data, dict):
        filtered = filter_record(record_type, data, candidates)
        return filtered

    return None


def main() -> int:
    args = parse_args()
    input_root = args.input_root.resolve()
    candidate_root = args.candidate_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    patients = 0
    files_written = 0
    records_written = 0
    skipped_no_candidate = 0

    for patient_dir in sorted(input_root.iterdir(), key=lambda p: p.name):
        if not patient_dir.is_dir():
            continue
        patient_id = patient_id_from_dir(patient_dir)
        if patient_id is None:
            continue
        cfile = candidate_file(candidate_root, patient_id)
        if cfile is None:
            skipped_no_candidate += 1
            continue
        bundle = load_json(cfile)
        candidates = candidate_names(bundle)
        if not candidates:
            skipped_no_candidate += 1
            continue

        patients += 1
        patient_output_dir = output_root / patient_dir.name
        for source_file in sorted(patient_dir.glob("*.json")):
            record_type = record_type_from_path(source_file, patient_id)
            if record_type not in MICROBE_RECORD_TYPES:
                continue
            filtered = filter_json_data(record_type, load_json(source_file), candidates)
            if filtered is None:
                continue
            if isinstance(filtered, list):
                count = len(filtered)
            else:
                count = 1
            if count == 0:
                continue
            write_json(patient_output_dir / source_file.name, filtered, args.pretty)
            files_written += 1
            records_written += count

    print(
        "patients={patients} skipped_no_candidate={skipped} files_written={files} records_written={records} output={output}".format(
            patients=patients,
            skipped=skipped_no_candidate,
            files=files_written,
            records=records_written,
            output=output_root,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import pathogen_normalization as pathogen_names  # noqa: E402
from tools.mngs_common import deduplicate_ranked_mngs_payload  # noqa: E402


DEFAULT_EXCEL_MAPPING = Path("excel_species_list_from_date_mapping.json")
DEFAULT_PATIENT_ROOT = Path("outputs") / "patient_info_filtered"
DEFAULT_REPORT = Path("outputs") / "chosen_ranked_mngs_from_excel_mapping_report.json"
SOURCE_SUFFIX = "_mNGS_ranked_candidates.json"
CHOSEN_SUFFIX = "_mNGS_ranked_candidates_chosen.json"

ALIASES = pathogen_names.alias_expansions()


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build per-patient *_mNGS_ranked_candidates_chosen.json files by "
            "filtering ranked mNGS candidates to species listed in the Excel/date mapping."
        )
    )
    parser.add_argument("--excel-mapping", type=Path, default=DEFAULT_EXCEL_MAPPING)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--no-overwrite", action="store_true")
    return parser.parse_args(argv)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict JSON from {path}, got {type(payload).__name__}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def normalize_name(value: Any) -> str:
    return pathogen_names.canonical_key(value)


def expanded_name_keys(value: Any) -> set[str]:
    return pathogen_names.expanded_name_keys(value)


def patient_id_from_dir(path: Path) -> str:
    name = path.name.removesuffix("_json")
    return name.removeprefix("NGS_patient_")


def ranked_path_for_patient(patient_dir: Path, patient_id: str) -> Path:
    base = f"NGS_patient_{patient_id}"
    return patient_dir / f"{base}{SOURCE_SUFFIX}"


def chosen_path_for_patient(patient_dir: Path, patient_id: str) -> Path:
    base = f"NGS_patient_{patient_id}"
    return patient_dir / f"{base}{CHOSEN_SUFFIX}"


def collect_excel_species(mapping_payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    patients = mapping_payload.get("patients")
    if not isinstance(patients, dict):
        raise KeyError("excel mapping JSON missing patients dict")

    result: dict[str, dict[str, Any]] = {}
    for patient_id, patient_payload in patients.items():
        if not isinstance(patient_payload, dict):
            continue
        organisms = patient_payload.get("organisms")
        if not isinstance(organisms, list):
            organisms = []
        names: list[str] = []
        normalized_keys: set[str] = set()
        for organism in organisms:
            if not isinstance(organism, dict):
                continue
            name = str(organism.get("name") or "").strip()
            if not name:
                continue
            names.append(name)
            normalized_keys.update(expanded_name_keys(name))
        result[str(patient_id)] = {
            "source_names": sorted(set(names), key=lambda item: item.lower()),
            "normalized_keys": normalized_keys,
            "specimen_ids": patient_payload.get("specimen_ids", []),
        }
    return result


def filter_ranked_payload(
    ranked_payload: dict[str, Any],
    *,
    patient_id: str,
    excel_entry: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    wanted_keys: set[str] = set(excel_entry.get("normalized_keys") or set())
    matched_keys: set[str] = set()
    matched_names: set[str] = set()
    input_count = 0
    output_count = 0

    filtered_records: list[dict[str, Any]] = []
    for record in ranked_payload.get("records") or []:
        if not isinstance(record, dict):
            continue
        filtered_candidates: list[dict[str, Any]] = []
        for candidate in record.get("candidates") or []:
            if not isinstance(candidate, dict):
                continue
            input_count += 1
            organism_name = candidate.get("organism_name") or candidate.get("name")
            candidate_keys = expanded_name_keys(organism_name)
            if candidate_keys & wanted_keys:
                filtered_candidates.append(candidate)
                matched_keys.update(candidate_keys & wanted_keys)
                matched_names.add(str(organism_name or "").strip())
        if filtered_candidates:
            new_record = dict(record)
            new_record["candidates"] = filtered_candidates
            new_record["status"] = "parsed"
            filtered_records.append(new_record)
            output_count += len(filtered_candidates)

    source_names = set(excel_entry.get("source_names") or [])
    unmatched_excel_names = [
        name
        for name in sorted(source_names, key=lambda item: item.lower())
        if not (expanded_name_keys(name) & matched_keys)
    ]

    metadata = dict(ranked_payload.get("ranking_metadata") or {})
    metadata["chosen_filter"] = {
        "enabled": True,
        "source": "excel_species_list_from_date_mapping",
        "excel_species_count": len(source_names),
        "matched_ranked_species_count": len(matched_names),
        "input_candidate_rows": input_count,
        "output_candidate_rows": output_count,
        "unmatched_excel_species": unmatched_excel_names,
    }

    chosen_payload = deduplicate_ranked_mngs_payload({
        "patient_id": str(ranked_payload.get("patient_id") or patient_id),
        "records": filtered_records,
        "ranking_metadata": metadata,
    })
    dedup = (chosen_payload.get("ranking_metadata") or {}).get("deduplication") or {}
    report_entry = {
        "patient_id": patient_id,
        "excel_species_count": len(source_names),
        "matched_ranked_species_count": len(matched_names),
        "input_candidate_rows": input_count,
        "output_candidate_rows": output_count,
        "output_candidate_rows_after_dedup": dedup.get("output_candidate_rows", output_count),
        "matched_ranked_species": sorted(matched_names, key=lambda item: item.lower()),
        "unmatched_excel_species": unmatched_excel_names,
    }
    return chosen_payload, report_entry


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    mapping_payload = load_json(args.excel_mapping)
    excel_by_patient = collect_excel_species(mapping_payload)

    report: dict[str, Any] = {
        "excel_mapping": str(args.excel_mapping),
        "patient_root": str(args.patient_root),
        "written": 0,
        "skipped_existing": 0,
        "skipped_missing_excel_patient": 0,
        "skipped_missing_ranked": 0,
        "patients": [],
    }

    patient_dirs = sorted(path for path in args.patient_root.glob("NGS_patient_*_json") if path.is_dir())
    for patient_dir in patient_dirs:
        patient_id = patient_id_from_dir(patient_dir)
        excel_entry = excel_by_patient.get(patient_id)
        if excel_entry is None:
            report["skipped_missing_excel_patient"] += 1
            continue
        source_path = ranked_path_for_patient(patient_dir, patient_id)
        if not source_path.exists():
            report["skipped_missing_ranked"] += 1
            continue
        destination = chosen_path_for_patient(patient_dir, patient_id)
        if destination.exists() and args.no_overwrite:
            report["skipped_existing"] += 1
            continue
        ranked_payload = load_json(source_path)
        chosen_payload, report_entry = filter_ranked_payload(
            ranked_payload,
            patient_id=patient_id,
            excel_entry=excel_entry,
        )
        write_json(destination, chosen_payload)
        report_entry["source_path"] = str(source_path)
        report_entry["output_path"] = str(destination)
        report["patients"].append(report_entry)
        report["written"] += 1

    write_json(args.report, report)
    print(f"written={report['written']}")
    print(f"skipped_missing_excel_patient={report['skipped_missing_excel_patient']}")
    print(f"skipped_missing_ranked={report['skipped_missing_ranked']}")
    print(f"report={args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

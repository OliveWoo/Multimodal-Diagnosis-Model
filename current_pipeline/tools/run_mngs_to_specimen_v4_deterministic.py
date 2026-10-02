from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

from tools.run_mngs_to_specimen_agent import (
    GROUP_KEYS,
    _sortable_patient_key,
    augment_curated_output,
    build_specimen_payload,
    derive_specimen_entries_from_source,
    get_metadata_entries,
    load_json,
    write_json,
)


DEFAULT_SOURCE = Path("outputs") / "mngs_candidate_microbes_rich_ranked.json"
DEFAULT_METADATA = Path("outputs") / "mngs_candidate_microbes_rich_metadata.json"
DEFAULT_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs_rich_v4_deterministic.json"


def empty_grouped_pathogens() -> dict[str, list[dict[str, Any]]]:
    return {group_name: [] for group_name in GROUP_KEYS}


def build_empty_curated(llm_input: dict[str, Any]) -> dict[str, Any]:
    return {
        "specimen_code": llm_input.get("specimen_code", ""),
        "collected_time": llm_input.get("collected_time", ""),
        "specimen_site": llm_input.get("specimen_site", ""),
        "pathogens_by_likelihood": {
            "high": empty_grouped_pathogens(),
            "medium": empty_grouped_pathogens(),
            "low_colonizer": empty_grouped_pathogens(),
        },
        "pathogens": empty_grouped_pathogens(),
    }


def run_deterministic(
    source_json: Path,
    metadata_json: Path,
    output_json: Path,
    patient_id: str | None = None,
) -> dict[str, Any]:
    source_payload = load_json(source_json)
    patients = source_payload.get("patients", {})
    if not isinstance(patients, dict):
        raise TypeError("Expected source payload to contain a dict at 'patients'")

    metadata_payload: dict[str, Any] = {}
    if metadata_json.exists():
        metadata_payload = load_json(metadata_json)

    patient_ids = [str(patient_id)] if patient_id else sorted(patients.keys(), key=_sortable_patient_key)
    results: dict[str, list[dict[str, Any]]] = {}

    for current_patient_id in patient_ids:
        patient_payload = patients.get(current_patient_id)
        if not isinstance(patient_payload, dict):
            raise KeyError(f"Patient id not found or not a dict: {current_patient_id}")

        metadata_entries = get_metadata_entries(metadata_payload, current_patient_id)
        if not metadata_entries:
            metadata_entries = derive_specimen_entries_from_source(patient_payload)
        if not metadata_entries:
            results[current_patient_id] = []
            continue

        patient_results: list[dict[str, Any]] = []
        for metadata_entry in metadata_entries:
            llm_input = build_specimen_payload(
                current_patient_id,
                patient_payload,
                metadata_entry,
                metadata_payload,
            )
            curated = augment_curated_output(
                build_empty_curated(llm_input),
                llm_input.get("candidates", []),
                llm_input.get("related_specimens", []),
            )
            curated["_generation_mode"] = "deterministic_v4_no_llm"
            patient_results.append(curated)
            print(f"patient={current_patient_id} specimen={metadata_entry.get('specimen_code', '')} completed")

        results[current_patient_id] = patient_results

    write_json(output_json, results)
    return results


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build mNGS-to-specimen tiered outputs using the deterministic v4 "
            "post-processing rules only. This does not call the LLM."
        )
    )
    parser.add_argument("--source-json", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--metadata-json", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--patient-id")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    results = run_deterministic(
        args.source_json,
        args.metadata_json,
        args.output,
        patient_id=args.patient_id,
    )
    specimen_count = sum(len(specimens) for specimens in results.values())
    print(f"Wrote deterministic v4 results to {args.output}")
    print(f"patient_count={len(results)}")
    print(f"specimen_count={specimen_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

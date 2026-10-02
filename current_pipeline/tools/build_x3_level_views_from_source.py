from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from tools.run_mngs_to_specimen_agent import augment_curated_output, normalize_organism_name

DEFAULT_INPUT = (
    Path("Specimen") / "specimen_lookup_summary_mNGS_grouped_with_type_filtered_code_ntc.json"
)
DEFAULT_OUTPUT_FULL = Path("outputs") / "mngs_to_specimen_outputs_x3_from_source.json"
DEFAULT_OUTPUT_LEVEL = Path("outputs") / "mngs_to_specimen_level_views_x3_from_source.json"

GROUP_KEYS = ("bacterial", "viral", "fungal", "others")
GROUP_TO_CATEGORY = {
    "bacterial": "1.Bac",
    "fungal": "2.Fungi",
    "viral": "3.Virus",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build x3-style tiered outputs directly from source specimen JSON "
            "(without using mngs_to_specimen_outputs.json)."
        )
    )
    parser.add_argument("--input-json", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-full-json", type=Path, default=DEFAULT_OUTPUT_FULL)
    parser.add_argument("--output-level-json", type=Path, default=DEFAULT_OUTPUT_LEVEL)
    return parser.parse_args()


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict payload from {path}, got {type(payload)}")
    return payload


def _to_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _to_reads(value: Any) -> int | float:
    numeric = _to_float(value)
    return int(numeric) if float(numeric).is_integer() else numeric


def _normalize_group(pathogens: Any) -> dict[str, list[dict[str, Any]]]:
    normalized = {key: [] for key in GROUP_KEYS}
    if not isinstance(pathogens, dict):
        return normalized
    for group_name in GROUP_KEYS:
        values = pathogens.get(group_name)
        if isinstance(values, list):
            normalized[group_name] = [item for item in values if isinstance(item, dict)]
    return normalized


def _build_candidates_from_pathogens(pathogens: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for group_name, source_category in GROUP_TO_CATEGORY.items():
        for item in pathogens.get(group_name, []):
            name = str(item.get("name", "")).strip()
            if not name:
                continue
            candidates.append(
                {
                    "organism_name": name,
                    "sec_hit": _to_reads(item.get("reads", 0)),
                    "source_category": source_category,
                }
            )
    return candidates


def _build_related_specimens(
    all_entries: list[dict[str, Any]],
    current_specimen_code: str,
) -> list[dict[str, Any]]:
    related: list[dict[str, Any]] = []
    for entry in all_entries:
        if not isinstance(entry, dict):
            continue
        specimen_code = str(entry.get("specimen_code", "")).strip()
        if not specimen_code or specimen_code == current_specimen_code:
            continue
        pathogens = _normalize_group(entry.get("pathogens"))
        related.append(
            {
                "specimen_code": specimen_code,
                "collected_time": entry.get("collected_time", ""),
                "specimen_site": entry.get("specimen_site", ""),
                "related_candidates": _build_candidates_from_pathogens(pathogens),
            }
        )
    return related


def _merge_groups(*groups: dict[str, list[dict[str, Any]]]) -> dict[str, list[dict[str, Any]]]:
    merged = {key: [] for key in GROUP_KEYS}
    seen = {key: set() for key in GROUP_KEYS}
    for group in groups:
        normalized_group = _normalize_group(group)
        for group_name in GROUP_KEYS:
            for item in normalized_group[group_name]:
                name = str(item.get("name", "")).strip()
                if not name:
                    continue
                dedupe_key = normalize_organism_name(name)
                if not dedupe_key or dedupe_key in seen[group_name]:
                    continue
                seen[group_name].add(dedupe_key)
                merged[group_name].append({"name": name, "reads": _to_reads(item.get("reads", 0))})
    return merged


def _build_level_views_payload(full_payload: dict[str, Any]) -> dict[str, Any]:
    by_patient: dict[str, Any] = {}
    for patient_id, entries in full_payload.items():
        if not isinstance(entries, list):
            continue
        output_entries: list[dict[str, Any]] = []
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            tiered = entry.get("pathogens_by_likelihood", {})
            high = _normalize_group(tiered.get("high"))
            medium = _normalize_group(tiered.get("medium"))
            low = _normalize_group(tiered.get("low_colonizer", tiered.get("low", {})))
            output_entries.append(
                {
                    "specimen_code": entry.get("specimen_code", ""),
                    "collected_time": entry.get("collected_time", ""),
                    "specimen_site": entry.get("specimen_site", ""),
                    "level_views": {
                        "high": high,
                        "high_medium": _merge_groups(high, medium),
                        "all_levels": _merge_groups(high, medium, low),
                    },
                }
            )
        by_patient[str(patient_id)] = output_entries
    return by_patient


def main() -> None:
    args = parse_args()
    source_payload = _load_json(args.input_json)

    full_output: dict[str, Any] = {}
    for patient_id, entries in source_payload.items():
        if not isinstance(entries, list):
            continue

        patient_outputs: list[dict[str, Any]] = []
        normalized_entries = [entry for entry in entries if isinstance(entry, dict)]
        for entry in normalized_entries:
            specimen_code = str(entry.get("specimen_code", "")).strip()
            pathogens = _normalize_group(entry.get("pathogens"))
            curated = {
                "specimen_code": specimen_code,
                "collected_time": entry.get("collected_time", ""),
                "specimen_site": entry.get("specimen_site", ""),
                "pathogens": pathogens,
            }
            candidates = _build_candidates_from_pathogens(pathogens)
            related = _build_related_specimens(normalized_entries, specimen_code)
            patient_outputs.append(augment_curated_output(curated, candidates, related))
        full_output[str(patient_id)] = patient_outputs

    level_output = _build_level_views_payload(full_output)

    args.output_full_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_full_json.open("w", encoding="utf-8") as handle:
        json.dump(full_output, handle, ensure_ascii=False, indent=2)

    args.output_level_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_level_json.open("w", encoding="utf-8") as handle:
        json.dump(level_output, handle, ensure_ascii=False, indent=2)

    specimen_count = sum(len(v) for v in full_output.values() if isinstance(v, list))
    print(
        f"Wrote full x3 output to {args.output_full_json} "
        f"(patients={len(full_output)}, specimens={specimen_count})"
    )
    print(f"Wrote x3 level views to {args.output_level_json}")


if __name__ == "__main__":
    main()

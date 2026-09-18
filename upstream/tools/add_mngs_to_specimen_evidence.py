from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any, Sequence

from tools.run_mngs_to_specimen_agent import CATEGORY_TO_GROUP, GROUP_KEYS, normalize_organism_name


DEFAULT_INPUT = Path("outputs") / "mngs_to_specimen_outputs_rich_all_rk_ntc_v4.json"
DEFAULT_RANKED = Path("outputs") / "mngs_candidate_microbes_rich_all_rk_ntc_ranked.json"
DEFAULT_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs_rich_all_rk_ntc_opt_with_evidence.json"
DEFAULT_REPORT = Path("outputs") / "mngs_to_specimen_outputs_rich_all_rk_ntc_opt_with_evidence_report.json"


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def to_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def evidence_sort_key(candidate: dict[str, Any]) -> tuple[int, float, float, str]:
    ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
    rank_priority = ranking.get("rank_priority")
    try:
        rank_value = int(rank_priority)
    except (TypeError, ValueError):
        rank_value = 999
    return (
        rank_value,
        -to_float(candidate.get("sec_hit", candidate.get("reads", 0))),
        -to_float(ranking.get("final_score", 0)),
        str(candidate.get("organism_name", "")),
    )


def build_evidence_index(
    ranked_payload: dict[str, Any],
) -> dict[tuple[str, str, str, str], dict[str, Any]]:
    index: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    patients = ranked_payload.get("patients", {})
    if not isinstance(patients, dict):
        return index

    for patient_id, patient_payload in patients.items():
        if not isinstance(patient_payload, dict):
            continue
        for record in patient_payload.get("records", []):
            if not isinstance(record, dict):
                continue
            specimen_code = str(record.get("specimen_code") or "").strip()
            if not specimen_code:
                continue
            for candidate in record.get("candidates", []):
                if not isinstance(candidate, dict):
                    continue
                source_category = str(candidate.get("source_category", ""))
                group_name = CATEGORY_TO_GROUP.get(source_category)
                normalized = normalize_organism_name(candidate.get("organism_name"))
                if not group_name or not normalized:
                    continue
                key = (str(patient_id), specimen_code, group_name, normalized)
                current = index.get(key)
                if current is None or evidence_sort_key(candidate) < evidence_sort_key(current):
                    index[key] = candidate
    return index


def evidence_payload(candidate: dict[str, Any]) -> dict[str, Any]:
    ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
    selection_reasons = list(candidate.get("selection_reasons") or [])
    return {
        "source_name": candidate.get("organism_name", ""),
        "source_category": candidate.get("source_category", ""),
        "source_reads": candidate.get("reads", candidate.get("sec_hit", 0)),
        "sec_hit": candidate.get("sec_hit", candidate.get("reads", 0)),
        "rk_ntc": selection_reasons,
        "selection_reasons": selection_reasons,
        "rank_priority": ranking.get("rank_priority"),
        "rank_rule": ranking.get("rank_rule"),
        "reads_percentile": ranking.get("reads_percentile"),
        "possibility_level": ranking.get("possibility_level"),
        "code_score": ranking.get("code_score"),
        "final_score": ranking.get("final_score"),
    }


def enrich_item(
    item: dict[str, Any],
    *,
    patient_id: str,
    specimen_code: str,
    group_name: str,
    evidence_index: dict[tuple[str, str, str, str], dict[str, Any]],
) -> tuple[dict[str, Any], bool]:
    enriched = dict(item)
    normalized = normalize_organism_name(item.get("name"))
    candidate = evidence_index.get((patient_id, specimen_code, group_name, normalized))
    if candidate is None:
        enriched["evidence_matched"] = False
        return enriched, False
    enriched.update(evidence_payload(candidate))
    enriched["evidence_matched"] = True
    return enriched, True


def enrich_group(
    group_payload: Any,
    *,
    patient_id: str,
    specimen_code: str,
    evidence_index: dict[tuple[str, str, str, str], dict[str, Any]],
) -> tuple[dict[str, list[dict[str, Any]]], int, int]:
    enriched_group: dict[str, list[dict[str, Any]]] = {group_name: [] for group_name in GROUP_KEYS}
    matched = 0
    total = 0
    if not isinstance(group_payload, dict):
        return enriched_group, matched, total
    for group_name in GROUP_KEYS:
        values = group_payload.get(group_name)
        if not isinstance(values, list):
            continue
        for item in values:
            if not isinstance(item, dict):
                continue
            total += 1
            enriched_item, did_match = enrich_item(
                item,
                patient_id=patient_id,
                specimen_code=specimen_code,
                group_name=group_name,
                evidence_index=evidence_index,
            )
            if did_match:
                matched += 1
            enriched_group[group_name].append(enriched_item)
    return enriched_group, matched, total


def enrich_outputs(
    output_payload: dict[str, Any],
    ranked_payload: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    evidence_index = build_evidence_index(ranked_payload)
    enriched_payload = copy.deepcopy(output_payload)
    total_items = 0
    matched_items = 0
    unmatched: list[dict[str, str]] = []

    for patient_id, specimen_outputs in enriched_payload.items():
        if not isinstance(specimen_outputs, list):
            continue
        for specimen in specimen_outputs:
            if not isinstance(specimen, dict):
                continue
            specimen_code = str(specimen.get("specimen_code") or "").strip()
            levels = specimen.get("pathogens_by_likelihood")
            if isinstance(levels, dict):
                for tier_name in ("high", "medium", "low_colonizer"):
                    enriched_group, matched, total = enrich_group(
                        levels.get(tier_name),
                        patient_id=str(patient_id),
                        specimen_code=specimen_code,
                        evidence_index=evidence_index,
                    )
                    levels[tier_name] = enriched_group
                    matched_items += matched
                    total_items += total
                    for group_name, values in enriched_group.items():
                        for item in values:
                            if not item.get("evidence_matched"):
                                unmatched.append(
                                    {
                                        "patient_id": str(patient_id),
                                        "specimen_code": specimen_code,
                                        "tier": tier_name,
                                        "group": group_name,
                                        "name": str(item.get("name", "")),
                                    }
                                )

            pathogens, matched, total = enrich_group(
                specimen.get("pathogens"),
                patient_id=str(patient_id),
                specimen_code=specimen_code,
                evidence_index=evidence_index,
            )
            specimen["pathogens"] = pathogens
            matched_items += matched
            total_items += total
            for group_name, values in pathogens.items():
                for item in values:
                    if not item.get("evidence_matched"):
                        unmatched.append(
                            {
                                "patient_id": str(patient_id),
                                "specimen_code": specimen_code,
                                "tier": "pathogens",
                                "group": group_name,
                                "name": str(item.get("name", "")),
                            }
                        )

    report = {
        "output_items_seen": total_items,
        "evidence_matched_items": matched_items,
        "evidence_unmatched_items": total_items - matched_items,
        "unmatched": unmatched,
    }
    return enriched_payload, report


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add ranked mNGS evidence fields to mNGS-to-specimen tiered output JSON."
    )
    parser.add_argument("--input-json", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--ranked-json", type=Path, default=DEFAULT_RANKED)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-json", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_payload = load_json(args.input_json)
    ranked_payload = load_json(args.ranked_json)
    enriched, report = enrich_outputs(output_payload, ranked_payload)
    write_json(args.output_json, enriched)
    report.update(
        {
            "input_json": str(args.input_json),
            "ranked_json": str(args.ranked_json),
            "output_json": str(args.output_json),
        }
    )
    write_json(args.report_json, report)
    print(f"Wrote enriched output to {args.output_json}")
    print(f"Wrote evidence report to {args.report_json}")
    print(
        f"evidence matched {report['evidence_matched_items']}/"
        f"{report['output_items_seen']} items"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any, Sequence


DEFAULT_INPUT = Path("mngs_candidate_microbes_for_llm_kmuh40_tsg_vghtpe_read_count_filtered.json")
DEFAULT_OUTPUT = Path("mngs_candidate_microbes_for_llm_kmuh40_tsg_vghtpe_read_count_filtered_ranked.json")

CODE_SCORES = {
    1: 1.0,
    2: 0.85,
    3: 0.65,
    4: 0.45,
    5: 0.30,
}

RANK_RULES = {
    1: "Code_NTC=00 + Code_RK_NTC=RK00",
    2: "Code_RK_NTC=RK8 + Code_NTC=00",
    3: "Code_RK_NTC=RK00",
    4: "Code_RK_NTC=RK8",
    5: "Code_NTC=00",
}

SOURCE_CATEGORY_BY_TYPE = {
    "bacteria": "1.Bac",
    "bacterial": "1.Bac",
    "fungi": "2.Fungi",
    "fungal": "2.Fungi",
    "virus": "3.Virus",
    "viral": "3.Virus",
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rank mNGS candidate microbes using the same corrected priority and reads "
            "percentile rules as mngs_candidate_microbes_ranked.json."
        )
    )
    parser.add_argument("--input-json", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT)
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


def to_float(value: Any) -> float:
    if value is None:
        return 0.0
    if isinstance(value, bool):
        return float(int(value))
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def normalize_reason(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    return text.upper().replace(" ", "")


def rank_priority(selection_reasons: Any) -> int:
    normalized = {normalize_reason(reason) for reason in selection_reasons or []}
    has_ntc00 = "CODE_NTC=00" in normalized or "CODENTC=00" in normalized
    has_rk00 = "CODE_RK_NTC=RK00" in normalized or "CODERKNTC=RK00" in normalized
    has_rk8 = "CODE_RK_NTC=RK8" in normalized or "CODE_RK_NTC=RK08" in normalized
    has_rk8 = has_rk8 or "CODERKNTC=RK8" in normalized or "CODERKNTC=RK08" in normalized

    if has_ntc00 and has_rk00:
        return 1
    if has_rk8 and has_ntc00:
        return 2
    if has_rk00:
        return 3
    if has_rk8:
        return 4
    if has_ntc00:
        return 5
    return 999


def possibility_level(final_score: float) -> str:
    if final_score >= 0.70:
        return "high"
    if final_score >= 0.45:
        return "medium"
    return "low"


def read_value(item: dict[str, Any]) -> float:
    return to_float(item.get("reads", item.get("sec_hit", 0)))


def build_read_percentiles(payload: dict[str, Any]) -> dict[float, float]:
    values: list[float] = []
    for items in payload.values():
        if not isinstance(items, list):
            continue
        for item in items:
            if isinstance(item, dict):
                values.append(math.log10(read_value(item) + 1.0))

    if not values:
        return {}

    sorted_values = sorted(values)
    total = len(sorted_values)
    percentiles: dict[float, float] = {}
    index = 0
    for value in sorted_values:
        while index < total and sorted_values[index] <= value:
            index += 1
        percentiles[value] = round(index / total, 4)
    return percentiles


def ranked_candidate(item: dict[str, Any], reads_percentiles: dict[float, float]) -> dict[str, Any]:
    reads = read_value(item)
    reads_out: int | float = int(reads) if reads.is_integer() else reads
    priority = rank_priority(item.get("selection_reasons"))
    code_score = CODE_SCORES.get(priority, 0.0)
    reads_log = math.log10(reads + 1.0)
    reads_percentile = reads_percentiles.get(reads_log, 0.0)
    final_score = round((0.60 * code_score) + (0.40 * reads_percentile), 4)
    organism_name = str(item.get("organism_name") or item.get("name") or "").strip()
    type_value = str(item.get("type") or "").strip()

    candidate = {
        "organism_name": organism_name,
        "sec_hit": reads_out,
        "selection_reasons": list(item.get("selection_reasons") or []),
        "source_category": SOURCE_CATEGORY_BY_TYPE.get(type_value.lower(), type_value or "Unknown"),
        "reads": reads_out,
        "type": type_value,
        "ranking": {
            "rank_priority": priority,
            "rank_rule": RANK_RULES.get(priority, "unranked"),
            "code_score": code_score,
            "reads_percentile": reads_percentile,
            "final_score": final_score,
            "possibility_level": possibility_level(final_score),
        },
    }
    return candidate


def sort_candidate(candidate: dict[str, Any]) -> tuple[int, float, str]:
    ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
    priority = ranking.get("rank_priority")
    try:
        priority_value = int(priority)
    except (TypeError, ValueError):
        priority_value = 999
    return priority_value, -to_float(candidate.get("reads")), str(candidate.get("organism_name", "")).lower()


def build_ranked_payload(payload: dict[str, Any]) -> dict[str, Any]:
    reads_percentiles = build_read_percentiles(payload)
    patients: dict[str, Any] = {}
    summary_counts: Counter[str] = Counter()
    rank_counts: Counter[str] = Counter()
    total_candidates = 0

    for patient_id in sorted(payload):
        items = payload.get(patient_id)
        if not isinstance(items, list):
            continue

        candidates = [
            ranked_candidate(item, reads_percentiles)
            for item in items
            if isinstance(item, dict)
        ]
        candidates.sort(key=sort_candidate)

        for candidate in candidates:
            ranking = candidate["ranking"]
            summary_counts[str(ranking["possibility_level"])] += 1
            rank_counts[str(ranking["rank_priority"])] += 1
        total_candidates += len(candidates)

        patients[str(patient_id)] = {
            "records": [
                {
                    "specimen_code": str(patient_id),
                    "chip_ids": [],
                    "seq_id": None,
                    "status": "parsed" if candidates else "parsed_no_match",
                    "candidates": candidates,
                }
            ]
        }

    return {
        "patients": patients,
        "ranking_metadata": {
            "version": "v2_corrected_priority",
            "notes": [
                "Corrected professor priority is used as the primary signal.",
                "Priority: 1) NTC00+RK00, 2) RK8+NTC00, 3) RK00, 4) RK8, 5) NTC00",
                "reads (sec_hit) is transformed by log10(reads+1) and converted to empirical percentile as secondary signal.",
                "final_score = 0.60 * code_score + 0.40 * reads_percentile",
                "high: final_score >= 0.70; medium: 0.45 <= final_score < 0.70; low: final_score < 0.45",
            ],
            "summary_counts": dict(sorted(summary_counts.items())),
            "rank_counts": dict(sorted(rank_counts.items(), key=lambda item: int(item[0]) if item[0].isdigit() else 999)),
            "patient_count": len(patients),
            "candidate_count": total_candidates,
        },
    }


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    payload = load_json(args.input_json)
    ranked = build_ranked_payload(payload)
    write_json(args.output_json, ranked)
    metadata = ranked["ranking_metadata"]
    print(f"Wrote {args.output_json}")
    print(f"patients={metadata['patient_count']} candidates={metadata['candidate_count']}")
    print(f"summary_counts={metadata['summary_counts']}")
    print(f"rank_counts={metadata['rank_counts']}")


if __name__ == "__main__":
    main()

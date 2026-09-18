from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

from tools import compare_mngs_to_specimen_with_excel as cmp

DEFAULT_INPUT = Path("outputs") / "mngs_to_specimen_outputs.json"
DEFAULT_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs_optimized.json"
DEFAULT_REPORT = Path("outputs") / "mngs_to_specimen_vs_excel_report_optimized.json"
DEFAULT_CONFIG = Path("outputs") / "mngs_low_colonizer_optimization.json"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Greedy optimizer for low_colonizer output. Maximizes precision while "
            "enforcing recall >= min_recall."
        )
    )
    parser.add_argument("--input-json", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--specimen-excel", type=Path, default=None)
    parser.add_argument("--min-recall", type=float, default=0.95)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-json", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--config-json", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--scope",
        choices=("by_specimen", "by_patient"),
        default="by_specimen",
        help="Which recall/precision scope is constrained and optimized.",
    )
    parser.add_argument(
        "--view",
        choices=("high", "high_medium", "all_levels"),
        default="all_levels",
        help="Which likelihood view to optimize.",
    )
    parser.add_argument(
        "--min-output-only-count",
        type=int,
        default=1,
        help="Only consider names with at least this many output_only occurrences.",
    )
    parser.add_argument(
        "--allow-remove-matched",
        action="store_true",
        help=(
            "Allow removing names that currently have matched occurrences. "
            "Default behavior only removes names with matched_count == 0."
        ),
    )
    return parser.parse_args(argv)


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict payload from {path}, got {type(payload)}")
    return payload


def _extract_names_from_groups(pathogens_group: Any) -> set[str]:
    names: set[str] = set()
    if not isinstance(pathogens_group, dict):
        return names
    for values in pathogens_group.values():
        if not isinstance(values, list):
            continue
        for item in values:
            if not isinstance(item, dict):
                continue
            normalized = cmp._normalize_name(item.get("name"))
            if normalized:
                names.add(normalized)
    return names


def _build_level_sets(
    payload: dict[str, Any],
) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, set[str]]]:
    high_by_specimen: dict[str, set[str]] = defaultdict(set)
    medium_by_specimen: dict[str, set[str]] = defaultdict(set)
    low_by_specimen: dict[str, set[str]] = defaultdict(set)
    for entries in payload.values():
        if not isinstance(entries, list):
            continue
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            specimen_code = cmp._normalize_identifier(entry.get("specimen_code", ""))
            if not specimen_code:
                continue
            tiered = entry.get("pathogens_by_likelihood")
            if isinstance(tiered, dict):
                high_names = _extract_names_from_groups(tiered.get("high", {}))
                medium_names = _extract_names_from_groups(tiered.get("medium", {}))
                low_names = _extract_names_from_groups(
                    tiered.get("low_colonizer", tiered.get("low", {}))
                )
            else:
                legacy = _extract_names_from_groups(entry.get("pathogens", {}))
                high_names = set(legacy)
                medium_names = set()
                low_names = set()
            high_by_specimen[specimen_code].update(high_names)
            medium_by_specimen[specimen_code].update(medium_names)
            low_by_specimen[specimen_code].update(low_names)
    return dict(high_by_specimen), dict(medium_by_specimen), dict(low_by_specimen)


def _compose_view_sets(
    high_by_specimen: dict[str, set[str]],
    medium_by_specimen: dict[str, set[str]],
    low_by_specimen: dict[str, set[str]],
    view: str,
) -> dict[str, set[str]]:
    composed: dict[str, set[str]] = {}
    specimen_codes = set(high_by_specimen) | set(medium_by_specimen) | set(low_by_specimen)
    for specimen_code in specimen_codes:
        high_names = high_by_specimen.get(specimen_code, set())
        medium_names = medium_by_specimen.get(specimen_code, set())
        low_names = low_by_specimen.get(specimen_code, set())
        if view == "high":
            composed[specimen_code] = set(high_names)
        elif view == "high_medium":
            composed[specimen_code] = set(high_names) | set(medium_names)
        else:
            composed[specimen_code] = set(high_names) | set(medium_names) | set(low_names)
    return composed


def _aggregate_by_patient(
    by_specimen: dict[str, set[str]],
    payload: dict[str, Any],
) -> dict[str, set[str]]:
    by_patient: dict[str, set[str]] = defaultdict(set)
    for patient_id, entries in payload.items():
        if not isinstance(entries, list):
            continue
        patient_key = cmp._normalize_identifier(patient_id)
        if not patient_key:
            continue
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            specimen_code = cmp._normalize_identifier(entry.get("specimen_code", ""))
            if specimen_code and specimen_code in by_specimen:
                by_patient[patient_key].update(by_specimen[specimen_code])
    return dict(by_patient)


def _compute_metrics(lhs: dict[str, set[str]], rhs: dict[str, set[str]]) -> dict[str, Any]:
    shared = sorted(set(lhs) & set(rhs))
    total_match = 0
    total_lhs = 0
    total_rhs = 0
    details: dict[str, dict[str, int]] = {}
    for key in shared:
        left = lhs[key]
        right = rhs[key]
        matched = len(left & right)
        total_match += matched
        total_lhs += len(left)
        total_rhs += len(right)
        details[key] = {
            "matched_count": matched,
            "output_count": len(left),
            "excel_count": len(right),
            "output_only_count": len(left - right),
            "excel_only_count": len(right - left),
        }
    precision = total_match / total_lhs if total_lhs else 0.0
    recall = total_match / total_rhs if total_rhs else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "shared_keys": len(shared),
        "total_output_species": total_lhs,
        "total_excel_species": total_rhs,
        "total_matched_species": total_match,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "details": details,
    }


def _build_name_stats(details: dict[str, Any]) -> dict[str, dict[str, int]]:
    stats: dict[str, dict[str, int]] = defaultdict(lambda: {"output_only": 0, "matched": 0})
    for key, key_detail in details.items():
        output_set = set(key_detail.get("output_set", []))
        excel_set = set(key_detail.get("excel_set", []))
        for name in output_set - excel_set:
            stats[name]["output_only"] += 1
        for name in output_set & excel_set:
            stats[name]["matched"] += 1
    return dict(stats)


def _build_details_with_sets(lhs: dict[str, set[str]], rhs: dict[str, set[str]]) -> dict[str, Any]:
    shared = sorted(set(lhs) & set(rhs))
    details: dict[str, Any] = {}
    for key in shared:
        left = set(lhs[key])
        right = set(rhs[key])
        details[key] = {
            "output_set": sorted(left),
            "excel_set": sorted(right),
        }
    return details


def _optimize_low_colonizer(
    *,
    low_by_specimen: dict[str, set[str]],
    fixed_view_by_specimen: dict[str, set[str]],
    payload: dict[str, Any],
    scope: str,
    excel_by_specimen: dict[str, set[str]],
    excel_by_patient: dict[str, set[str]],
    min_recall: float,
    min_output_only_count: int,
    allow_remove_matched: bool,
) -> tuple[dict[str, set[str]], list[str], dict[str, Any], dict[str, Any]]:
    def score_scope(by_specimen_output: dict[str, set[str]]) -> dict[str, Any]:
        if scope == "by_patient":
            by_patient_output = _aggregate_by_patient(by_specimen_output, payload)
            return _compute_metrics(by_patient_output, excel_by_patient)
        return _compute_metrics(by_specimen_output, excel_by_specimen)

    current_low = {k: set(v) for k, v in low_by_specimen.items()}
    current_output = {
        specimen_code: set(fixed_view_by_specimen.get(specimen_code, set())) | set(current_low.get(specimen_code, set()))
        for specimen_code in set(fixed_view_by_specimen) | set(current_low)
    }
    baseline_metrics = score_scope(current_output)
    initial_metrics = dict(baseline_metrics)

    detail_sets = _build_details_with_sets(current_output, excel_by_specimen)
    stats = _build_name_stats(detail_sets)
    candidates = []
    for name, stat in stats.items():
        output_only_count = stat["output_only"]
        matched_count = stat["matched"]
        if output_only_count < min_output_only_count:
            continue
        if not allow_remove_matched and matched_count > 0:
            continue
        candidates.append((name, output_only_count, matched_count))
    candidates.sort(key=lambda item: (item[2], -item[1], item[0]))

    removed: list[str] = []
    current_precision = baseline_metrics["precision"]
    current_recall = baseline_metrics["recall"]

    name_to_specs: dict[str, set[str]] = defaultdict(set)
    for specimen_code, names in current_low.items():
        for name in names:
            name_to_specs[name].add(specimen_code)

    for name, _, _ in candidates:
        target_specs = name_to_specs.get(name, set())
        if not target_specs:
            continue

        trial_low = {k: set(v) for k, v in current_low.items()}
        for specimen_code in target_specs:
            if specimen_code in trial_low:
                trial_low[specimen_code].discard(name)

        trial_output = {
            specimen_code: set(fixed_view_by_specimen.get(specimen_code, set())) | set(trial_low.get(specimen_code, set()))
            for specimen_code in set(fixed_view_by_specimen) | set(trial_low)
        }
        trial_metrics = score_scope(trial_output)
        trial_recall = trial_metrics["recall"]
        trial_precision = trial_metrics["precision"]

        if trial_recall < min_recall:
            continue
        if trial_precision > current_precision or (
            abs(trial_precision - current_precision) < 1e-12
            and trial_metrics["total_output_species"] < baseline_metrics["total_output_species"]
        ):
            current_low = trial_low
            current_output = trial_output
            current_precision = trial_precision
            current_recall = trial_recall
            baseline_metrics = trial_metrics
            removed.append(name)

    final_metrics = score_scope(current_output)
    return current_low, removed, initial_metrics, final_metrics


def _apply_removed_names_to_payload(payload: dict[str, Any], removed_names: set[str]) -> dict[str, Any]:
    if not removed_names:
        return payload
    updated: dict[str, Any] = {}
    for patient_id, entries in payload.items():
        if not isinstance(entries, list):
            updated[patient_id] = entries
            continue
        new_entries: list[Any] = []
        for entry in entries:
            if not isinstance(entry, dict):
                new_entries.append(entry)
                continue
            new_entry = dict(entry)
            tiered = new_entry.get("pathogens_by_likelihood")
            if isinstance(tiered, dict):
                low = tiered.get("low_colonizer")
                if isinstance(low, dict):
                    new_low: dict[str, Any] = {}
                    for group_name, values in low.items():
                        if not isinstance(values, list):
                            new_low[group_name] = values
                            continue
                        filtered = []
                        for item in values:
                            if not isinstance(item, dict):
                                continue
                            normalized = cmp._normalize_name(item.get("name"))
                            if normalized in removed_names:
                                continue
                            filtered.append(item)
                        new_low[group_name] = filtered
                    new_tiered = dict(tiered)
                    new_tiered["low_colonizer"] = new_low
                    new_entry["pathogens_by_likelihood"] = new_tiered
            new_entries.append(new_entry)
        updated[patient_id] = new_entries
    return updated


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    specimen_excel = args.specimen_excel or cmp._resolve_default_excel()
    payload = _load_json(args.input_json)
    excel_by_specimen, excel_by_patient = cmp._load_excel_lookup(specimen_excel)

    high_by_specimen, medium_by_specimen, low_by_specimen = _build_level_sets(payload)
    if args.view == "all_levels":
        fixed_view = _compose_view_sets(high_by_specimen, medium_by_specimen, {}, "high_medium")
    elif args.view == "high_medium":
        fixed_view = _compose_view_sets(high_by_specimen, {}, {}, "high")
    else:
        fixed_view = _compose_view_sets({}, {}, {}, "high")

    _optimized_low, removed_names, before_metrics, after_metrics = _optimize_low_colonizer(
        low_by_specimen=low_by_specimen if args.view == "all_levels" else {},
        fixed_view_by_specimen=fixed_view,
        payload=payload,
        scope=args.scope,
        excel_by_specimen=excel_by_specimen,
        excel_by_patient=excel_by_patient,
        min_recall=args.min_recall,
        min_output_only_count=args.min_output_only_count,
        allow_remove_matched=args.allow_remove_matched,
    )

    updated_payload = _apply_removed_names_to_payload(payload, set(removed_names))
    _write_json(args.output_json, updated_payload)

    output_by_specimen, output_by_patient = cmp._load_output_lookup_by_level(args.output_json)
    by_specimen_report = cmp._build_comparison_by_level(output_by_specimen, excel_by_specimen)
    by_patient_report = cmp._build_comparison_by_level(output_by_patient, excel_by_patient)
    report = {
        "input_json": str(args.output_json),
        "specimen_excel": str(specimen_excel),
        "likelihood_views": list(cmp.LIKELIHOOD_KEYS),
        "recommended_recall_view": "all_levels",
        "summary": cmp._build_summary(by_specimen_report, by_patient_report),
        "by_specimen": by_specimen_report,
        "by_patient": by_patient_report,
    }
    _write_json(args.report_json, report)

    optimization = {
        "source_input_json": str(args.input_json),
        "optimized_output_json": str(args.output_json),
        "report_json": str(args.report_json),
        "specimen_excel": str(specimen_excel),
        "scope": args.scope,
        "view": args.view,
        "min_recall": args.min_recall,
        "min_output_only_count": args.min_output_only_count,
        "allow_remove_matched": args.allow_remove_matched,
        "before": before_metrics,
        "after": after_metrics,
        "removed_name_count": len(removed_names),
        "removed_names": sorted(removed_names),
    }
    _write_json(args.config_json, optimization)

    print(
        "Optimization complete. "
        f"removed={len(removed_names)} "
        f"precision {before_metrics['precision']:.4f}->{after_metrics['precision']:.4f} "
        f"recall {before_metrics['recall']:.4f}->{after_metrics['recall']:.4f}"
    )
    print(f"Wrote optimized output to {args.output_json}")
    print(f"Wrote optimized report to {args.report_json}")
    print(f"Wrote optimization summary to {args.config_json}")


if __name__ == "__main__":
    main()

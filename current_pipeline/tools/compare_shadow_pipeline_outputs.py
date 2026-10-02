from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Iterable

from . import pathogen_normalization


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def _patient_id(patient_dir: Path) -> int:
    try:
        return int(patient_dir.name.removeprefix("NGS_patient_").removesuffix("_json"))
    except ValueError:
        return 10**9


def _patient_dirs(root: Path) -> list[Path]:
    return sorted(
        (path for path in root.glob("NGS_patient_*_json") if path.is_dir()),
        key=_patient_id,
    )


def _output_path(patient_dir: Path, suffix: str) -> Path:
    patient_name = patient_dir.name.removesuffix("_json")
    return patient_dir / "summary_outputs" / f"{patient_name}_{suffix}.json"


def _name(item: dict[str, Any]) -> str:
    return str(
        item.get("organism_name")
        or item.get("pathogen_name")
        or item.get("name")
        or ""
    ).strip()


def _key(item: dict[str, Any]) -> str:
    raw_name = _name(item)
    return pathogen_normalization.canonical_key(raw_name) or raw_name.casefold()


def _index(items: Iterable[Any]) -> dict[str, dict[str, Any]]:
    indexed: dict[str, dict[str, Any]] = {}
    for raw_item in items:
        if not isinstance(raw_item, dict):
            continue
        key = _key(raw_item)
        if key:
            indexed[key] = raw_item
    return indexed


def _picked(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    summary = payload.get("best_available_summary")
    if not isinstance(summary, dict):
        summary = {}
    return _index(summary.get("picked_pathogens") or [])


def _candidates(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return _index(payload.get("pathogen_candidates") or [])


def _queue(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return _index(payload.get("review_queue") or [])


def _display(item: dict[str, Any] | None, fallback: str) -> str:
    if item is None:
        return fallback
    return _name(item) or fallback


def _level(item: dict[str, Any]) -> str:
    value = item.get("integrated_causative_level")
    if value is None:
        value = item.get("basis_level")
    if value is None:
        value = item.get("level")
    return "" if value is None else str(value)


def _compare_maps(
    before: dict[str, dict[str, Any]],
    after: dict[str, dict[str, Any]],
) -> tuple[list[str], list[str]]:
    before_only = [
        _display(before[key], key) for key in sorted(set(before) - set(after))
    ]
    after_only = [
        _display(after[key], key) for key in sorted(set(after) - set(before))
    ]
    return before_only, after_only


def _compare_levels(
    before: dict[str, dict[str, Any]],
    after: dict[str, dict[str, Any]],
) -> list[dict[str, str]]:
    changes: list[dict[str, str]] = []
    for key in sorted(set(before) & set(after)):
        before_level = _level(before[key])
        after_level = _level(after[key])
        if before_level != after_level:
            changes.append(
                {
                    "organism": _display(after.get(key) or before.get(key), key),
                    "before": before_level,
                    "after": after_level,
                }
            )
    return changes


def compare_patient(
    patient_dir: Path,
    baseline_scorer_suffix: str,
    candidate_scorer_suffix: str,
    baseline_queue_suffix: str,
    candidate_queue_suffix: str,
) -> dict[str, Any]:
    paths = {
        "baseline_scorer": _output_path(patient_dir, baseline_scorer_suffix),
        "candidate_scorer": _output_path(patient_dir, candidate_scorer_suffix),
        "baseline_queue": _output_path(patient_dir, baseline_queue_suffix),
        "candidate_queue": _output_path(patient_dir, candidate_queue_suffix),
    }
    missing = [label for label, path in paths.items() if not path.is_file()]
    row: dict[str, Any] = {
        "patient": patient_dir.name,
        "missing_outputs": missing,
        "picked_removed": [],
        "picked_added": [],
        "candidate_removed": [],
        "candidate_added": [],
        "level_changes": [],
        "queue_removed": [],
        "queue_added": [],
    }
    if missing:
        return row

    baseline_scorer = _load_json(paths["baseline_scorer"])
    candidate_scorer = _load_json(paths["candidate_scorer"])
    baseline_queue = _load_json(paths["baseline_queue"])
    candidate_queue = _load_json(paths["candidate_queue"])

    row["picked_removed"], row["picked_added"] = _compare_maps(
        _picked(baseline_scorer), _picked(candidate_scorer)
    )
    row["candidate_removed"], row["candidate_added"] = _compare_maps(
        _candidates(baseline_scorer), _candidates(candidate_scorer)
    )
    row["level_changes"] = _compare_levels(
        _candidates(baseline_scorer), _candidates(candidate_scorer)
    )
    row["queue_removed"], row["queue_added"] = _compare_maps(
        _queue(baseline_queue), _queue(candidate_queue)
    )
    return row


def _join_names(values: list[str]) -> str:
    return " | ".join(values)


def _level_text(values: list[dict[str, str]]) -> str:
    return " | ".join(
        f"{value['organism']}: {value['before']} -> {value['after']}"
        for value in values
    )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "patient",
        "missing_outputs",
        "picked_removed",
        "picked_added",
        "candidate_removed",
        "candidate_added",
        "level_changes",
        "queue_removed",
        "queue_added",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "patient": row["patient"],
                    "missing_outputs": _join_names(row["missing_outputs"]),
                    "picked_removed": _join_names(row["picked_removed"]),
                    "picked_added": _join_names(row["picked_added"]),
                    "candidate_removed": _join_names(row["candidate_removed"]),
                    "candidate_added": _join_names(row["candidate_added"]),
                    "level_changes": _level_text(row["level_changes"]),
                    "queue_removed": _join_names(row["queue_removed"]),
                    "queue_added": _join_names(row["queue_added"]),
                }
            )


def _write_markdown(path: Path, report: dict[str, Any]) -> None:
    summary = report["summary"]
    lines = [
        "# Shadow Pipeline Comparison",
        "",
        "This report is answer-blind. It compares pipeline outputs only.",
        "",
        "## Summary",
        "",
        f"- Patients: {summary['patients']}",
        f"- Patients with missing outputs: {summary['patients_with_missing_outputs']}",
        f"- Patients with Picked membership changes: {summary['patients_with_picked_changes']}",
        f"- Patients with scorer candidate membership changes: {summary['patients_with_candidate_changes']}",
        f"- Patients with candidate level changes: {summary['patients_with_level_changes']}",
        f"- Patients with queue membership changes: {summary['patients_with_queue_changes']}",
        "",
        "## Differences",
        "",
    ]
    changed_rows = [
        row
        for row in report["patients"]
        if row["missing_outputs"]
        or row["picked_removed"]
        or row["picked_added"]
        or row["candidate_removed"]
        or row["candidate_added"]
        or row["level_changes"]
        or row["queue_removed"]
        or row["queue_added"]
    ]
    if not changed_rows:
        lines.append("No differences found.")
    for row in changed_rows:
        lines.extend([f"### {row['patient']}", ""])
        if row["missing_outputs"]:
            lines.append(f"- Missing: {_join_names(row['missing_outputs'])}")
        if row["picked_removed"]:
            lines.append(f"- Picked removed: {_join_names(row['picked_removed'])}")
        if row["picked_added"]:
            lines.append(f"- Picked added: {_join_names(row['picked_added'])}")
        if row["candidate_removed"]:
            lines.append(f"- Scorer candidates removed: {_join_names(row['candidate_removed'])}")
        if row["candidate_added"]:
            lines.append(f"- Scorer candidates added: {_join_names(row['candidate_added'])}")
        if row["level_changes"]:
            lines.append(f"- Level changes: {_level_text(row['level_changes'])}")
        if row["queue_removed"]:
            lines.append(f"- Queue removed: {_join_names(row['queue_removed'])}")
        if row["queue_added"]:
            lines.append(f"- Queue added: {_join_names(row['queue_added'])}")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare frozen and shadow deterministic scorer/review queue outputs "
            "without reading benchmark answers."
        )
    )
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--baseline-scorer-suffix", required=True)
    parser.add_argument("--candidate-scorer-suffix", required=True)
    parser.add_argument("--baseline-queue-suffix", required=True)
    parser.add_argument("--candidate-queue-suffix", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    rows = [
        compare_patient(
            patient_dir,
            args.baseline_scorer_suffix,
            args.candidate_scorer_suffix,
            args.baseline_queue_suffix,
            args.candidate_queue_suffix,
        )
        for patient_dir in _patient_dirs(args.patient_root)
    ]
    summary = {
        "patients": len(rows),
        "patients_with_missing_outputs": sum(bool(row["missing_outputs"]) for row in rows),
        "patients_with_picked_changes": sum(
            bool(row["picked_removed"] or row["picked_added"]) for row in rows
        ),
        "patients_with_candidate_changes": sum(
            bool(row["candidate_removed"] or row["candidate_added"]) for row in rows
        ),
        "patients_with_level_changes": sum(bool(row["level_changes"]) for row in rows),
        "patients_with_queue_changes": sum(
            bool(row["queue_removed"] or row["queue_added"]) for row in rows
        ),
    }
    report = {
        "metadata": {
            "patient_root": str(args.patient_root.resolve()),
            "baseline_scorer_suffix": args.baseline_scorer_suffix,
            "candidate_scorer_suffix": args.candidate_scorer_suffix,
            "baseline_queue_suffix": args.baseline_queue_suffix,
            "candidate_queue_suffix": args.candidate_queue_suffix,
            "answer_blind": True,
        },
        "summary": summary,
        "patients": rows,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "shadow_pipeline_comparison.json"
    csv_path = args.output_dir / "shadow_pipeline_comparison.csv"
    markdown_path = args.output_dir / "shadow_pipeline_comparison.md"
    json_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    _write_csv(csv_path, rows)
    _write_markdown(markdown_path, report)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"json={json_path}")
    print(f"csv={csv_path}")
    print(f"markdown={markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

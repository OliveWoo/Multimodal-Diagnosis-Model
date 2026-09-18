from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Callable


REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import pathogen_normalization as pathogen_names  # noqa: E402


DEFAULT_ANSWER_CSV = REPO_ROOT / "outputs" / "runs" / "legacy_reports" / "filtered" / (
    "20260715_20260728_pulmonary_opt"
) / "kmuh_filtered_latest_pulmonary_opt_answer_comparison.csv"
DEFAULT_PATIENT_ROOT = REPO_ROOT / "outputs" / "patient_info_KH_0728_2Days"
DEFAULT_DATABASE_CSV = REPO_ROOT / "outputs" / "patient_database_latest" / "patient_database_index.csv"
DEFAULT_OVERRIDES = REPO_ROOT / "rules" / "kh_benchmark_answer_overrides_20260914.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "outputs" / "runs" / "2026-09-14_KH_answer_revision_metrics"

SPLIT_RE = re.compile(r"[;\n\r\uFF1B\u3001]+")
RESISTANCE_TARGETS = {
    "ctxm",
    "ndm",
    "mecacandmrej",
    "imp",
    "oxa48like",
    "kpc",
    "vim",
}
WINDOW = timedelta(days=2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply versioned KH answer revisions and recalculate strict benchmark metrics."
    )
    parser.add_argument("--answer-csv", type=Path, default=DEFAULT_ANSWER_CSV)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--patient-database-csv", type=Path, default=DEFAULT_DATABASE_CSV)
    parser.add_argument("--overrides-json", type=Path, default=DEFAULT_OVERRIDES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_datetime(value: Any) -> datetime | None:
    text = str(value or "").strip()
    if not text:
        return None
    for pattern in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return datetime.strptime(text, pattern)
        except ValueError:
            continue
    return None


def mngs_reference_dates(patient_id: int, patient_root: Path) -> list[datetime]:
    path = patient_root / f"NGS_patient_{patient_id}_json" / (
        f"NGS_patient_{patient_id}_all_RK_NTC_microbes.json"
    )
    if not path.exists():
        return []
    payload = read_json(path)
    records = payload if isinstance(payload, list) else [payload]
    return [
        parsed
        for record in records
        if isinstance(record, dict)
        for parsed in [parse_datetime(record.get("collected_time"))]
        if parsed is not None
    ]


def within_mngs_window(value: Any, reference_dates: list[datetime]) -> bool:
    event_time = parse_datetime(value)
    if event_time is None or not reference_dates:
        return False
    return any(abs(event_time - reference) <= WINDOW for reference in reference_dates)


def split_names(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text or text in {"-", "NA"}:
        return []
    return [part.strip() for part in SPLIT_RE.split(text) if part.strip() and part.strip() != "-"]


def canonical(value: Any) -> str:
    return pathogen_names.canonical_key(pathogen_names.clean_display_text(value))


def unique_names(values: list[str]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        key = canonical(value)
        if not key or key.isdigit() or key in seen:
            continue
        seen.add(key)
        result.append(value.strip())
    return result


def load_original_answers(path: Path) -> tuple[list[dict[str, str]], dict[int, list[str]]]:
    rows = read_csv(path)
    answers: dict[int, list[str]] = {}
    for row in rows:
        patient_id = int(str(row.get("patient_id") or "0"))
        answers[patient_id] = split_names(row.get("answers") or row.get("answer"))
    return rows, answers


def load_database_rows(path: Path) -> dict[int, dict[str, str]]:
    rows = read_csv(path)
    return {
        int(row["patient_id"]): row
        for row in rows
        if row.get("dataset") == "KH" and str(row.get("patient_id") or "").isdigit()
    }


def load_overrides(path: Path) -> tuple[dict[int, list[str]], dict[int, list[str]], dict[str, Any]]:
    payload = read_json(path)
    names: dict[int, list[str]] = {}
    display: dict[int, list[str]] = {}
    for patient_id, items in payload.get("patients", {}).items():
        pid = int(patient_id)
        names[pid] = [str(item["organism"]).strip() for item in items]
        display[pid] = [f'{item["organism"]} ({item["reads"]})' for item in items]
    return names, display, payload


def apply_overrides(
    original: dict[int, list[str]], overrides: dict[int, list[str]]
) -> dict[int, list[str]]:
    revised = {patient_id: list(names) for patient_id, names in original.items()}
    for patient_id, names in overrides.items():
        revised[patient_id] = list(names)
    return revised


def match_type(output: str, answer: str, *, genus_relaxed: bool) -> str:
    match_type = pathogen_names.name_match_type(output, answer)
    if match_type in {"exact_or_alias", "approved_group_member_match"}:
        return match_type
    if genus_relaxed and match_type == "genus_relaxed":
        return match_type
    return ""


def pairwise_match(
    outputs: list[str], answers: list[str], *, genus_relaxed: bool
) -> dict[str, Any]:
    remaining = set(range(len(answers)))
    matched: list[dict[str, str]] = []
    output_only: list[str] = []
    for output in outputs:
        best_index: int | None = None
        best_type = ""
        for index in sorted(remaining):
            candidate_type = match_type(output, answers[index], genus_relaxed=genus_relaxed)
            if not candidate_type:
                continue
            if candidate_type == "exact_or_alias":
                best_index = index
                best_type = candidate_type
                break
            if best_index is None or (
                candidate_type == "approved_group_member_match" and best_type == "genus_relaxed"
            ):
                best_index = index
                best_type = candidate_type
        if best_index is None:
            output_only.append(output)
            continue
        remaining.remove(best_index)
        matched.append(
            {"output": output, "answer": answers[best_index], "match_type": best_type}
        )
    return {
        "matched": matched,
        "output_only": output_only,
        "answer_only": [answers[index] for index in sorted(remaining)],
    }


def model_outputs(patient_id: int, database_rows: dict[int, dict[str, str]]) -> list[str]:
    row = database_rows.get(patient_id, {})
    return unique_names(split_names(row.get("final_pathogens")))


def mngs_outputs(patient_id: int, patient_root: Path) -> list[str]:
    path = patient_root / f"NGS_patient_{patient_id}_json" / (
        f"NGS_patient_{patient_id}_all_RK_NTC_microbes.json"
    )
    if not path.exists():
        return []
    payload = read_json(path)
    names: list[str] = []
    records = payload if isinstance(payload, list) else [payload]
    for record in records:
        if not isinstance(record, dict):
            continue
        for items in (record.get("pathogens") or {}).values():
            if not isinstance(items, list):
                continue
            names.extend(
                str(item.get("name") or "").strip()
                for item in items
                if isinstance(item, dict) and item.get("name")
            )
    return unique_names(names)


def filmarray_outputs(patient_id: int, patient_root: Path) -> list[str]:
    return hospital_module_outputs(patient_id, patient_root, "filmarray_gmtest")


def culture_outputs(patient_id: int, patient_root: Path) -> list[str]:
    return hospital_module_outputs(patient_id, patient_root, "culture")


def hospital_module_outputs(patient_id: int, patient_root: Path, module: str) -> list[str]:
    path = patient_root / f"NGS_patient_{patient_id}_json" / "summary_outputs" / (
        f"NGS_patient_{patient_id}_final_summary_with_filmarray_deterministic.json"
    )
    if not path.exists():
        return []
    payload = read_json(path)
    names = [
        str(item.get("organism_name") or "").strip()
        for item in payload.get("hospital_organism_evidence", [])
        if isinstance(item, dict)
        and module in (item.get("evidence_modules") or [])
        and item.get("organism_name")
    ]
    return unique_names(names)


def f1_score(precision: float, recall: float) -> float:
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0


def evaluate(
    answers_by_patient: dict[int, list[str]],
    output_loader: Callable[[int], list[str]],
    *,
    genus_relaxed: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    total_answers = 0
    total_outputs = 0
    total_matches = 0
    details: list[dict[str, Any]] = []
    for patient_id in sorted(answers_by_patient):
        answers = answers_by_patient[patient_id]
        if not answers:
            continue
        outputs = output_loader(patient_id)
        comparison = pairwise_match(outputs, answers, genus_relaxed=genus_relaxed)
        matches = len(comparison["matched"])
        total_answers += len(answers)
        total_outputs += len(outputs)
        total_matches += matches
        details.append(
            {
                "patient_id": patient_id,
                "answers": answers,
                "outputs": outputs,
                **comparison,
            }
        )
    precision = total_matches / total_outputs if total_outputs else 0.0
    recall = total_matches / total_answers if total_answers else 0.0
    summary = {
        "matches": total_matches,
        "outputs": total_outputs,
        "answers": total_answers,
        "precision": precision,
        "recall": recall,
        "f1": f1_score(precision, recall),
        "precision_ratio": f"{total_matches}/{total_outputs}",
        "recall_ratio": f"{total_matches}/{total_answers}",
    }
    return summary, details


def write_revised_answers(
    path: Path,
    original_rows: list[dict[str, str]],
    revised: dict[int, list[str]],
    override_display: dict[int, list[str]],
    database_rows: dict[int, dict[str, str]],
) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["patient_id", "answers", "answers_with_reads", "revision_status"],
        )
        writer.writeheader()
        for source_row in original_rows:
            patient_id = int(source_row["patient_id"])
            names = revised.get(patient_id, [])
            if patient_id in override_display:
                display = override_display[patient_id]
                status = "revised_20260914"
            else:
                display = split_names(database_rows.get(patient_id, {}).get("hospital_answer"))
                status = "unchanged"
            writer.writerow(
                {
                    "patient_id": patient_id,
                    "answers": "; ".join(names) if names else "-",
                    "answers_with_reads": "\n".join(display),
                    "revision_status": status,
                }
            )


def main() -> None:
    args = parse_args()
    original_rows, original_answers = load_original_answers(args.answer_csv)
    database_rows = load_database_rows(args.patient_database_csv)
    overrides, override_display, override_payload = load_overrides(args.overrides_json)
    revised_answers = apply_overrides(original_answers, overrides)

    loaders: dict[str, Callable[[int], list[str]]] = {
        "model_r5": lambda patient_id: model_outputs(patient_id, database_rows),
        "mngs": lambda patient_id: mngs_outputs(patient_id, args.patient_root),
        "filmarray": lambda patient_id: filmarray_outputs(patient_id, args.patient_root),
        "culture": lambda patient_id: culture_outputs(patient_id, args.patient_root),
    }
    original_metrics: dict[str, dict[str, Any]] = {"strict": {}, "genus_relaxed": {}}
    revised_metrics: dict[str, dict[str, Any]] = {"strict": {}, "genus_relaxed": {}}
    revised_details: dict[str, dict[str, Any]] = {"strict": {}, "genus_relaxed": {}}
    for mode, allow_genus in (("strict", False), ("genus_relaxed", True)):
        for method, loader in loaders.items():
            original_metrics[mode][method], _ = evaluate(
                original_answers, loader, genus_relaxed=allow_genus
            )
            revised_metrics[mode][method], revised_details[mode][method] = evaluate(
                revised_answers, loader, genus_relaxed=allow_genus
            )

    expected_original = {
        "strict": {
            "model_r5": (52, 76, 58),
            "mngs": (53, 674, 58),
            "filmarray": (20, 34, 58),
            "culture": (23, 52, 58),
        },
        "genus_relaxed": {
            "model_r5": (53, 76, 58),
            "mngs": (54, 674, 58),
            "filmarray": (23, 34, 58),
            "culture": (26, 52, 58),
        },
    }
    validation: dict[str, dict[str, Any]] = {"strict": {}, "genus_relaxed": {}}
    for mode, expected_by_method in expected_original.items():
        for method, expected in expected_by_method.items():
            actual = original_metrics[mode][method]
            actual_tuple = (actual["matches"], actual["outputs"], actual["answers"])
            validation[mode][method] = {
                "expected_matches_outputs_answers": list(expected),
                "actual_matches_outputs_answers": list(actual_tuple),
                "passed": actual_tuple == expected,
            }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    revised_answer_csv = args.output_dir / "kh_answers_clinical_revision_20260914.csv"
    write_revised_answers(
        revised_answer_csv,
        original_rows,
        revised_answers,
        override_display,
        database_rows,
    )

    report = {
        "schema_version": "kh_benchmark_metrics.v1",
        "matching": {
            "strict": "exact/alias plus approved group-member; no generic same-genus matches",
            "genus_relaxed": "strict matches plus generic same-genus matches",
        },
        "source_answer_csv": str(args.answer_csv.resolve()),
        "revised_answer_csv": str(revised_answer_csv.resolve()),
        "overrides_json": str(args.overrides_json.resolve()),
        "revision_id": override_payload.get("revision_id"),
        "original_baseline_validation": validation,
        "original_metrics": original_metrics,
        "revised_metrics": revised_metrics,
        "revised_details": revised_details,
    }
    report_path = args.output_dir / "kh_revised_benchmark_metrics.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    summary_path = args.output_dir / "RESULTS.md"
    lines = [
        "# KH revised benchmark metrics",
        "",
        f"Revision: `{override_payload.get('revision_id')}`",
        "",
        "Both strict and genus-relaxed metrics are reported. The prior presentation's FilmArray and culture values used genus-relaxed matching.",
        "",
        "## Strict",
        "",
        "| Method | Precision | Recall | F1 |",
        "|---|---:|---:|---:|",
    ]
    for method in ("model_r5", "mngs", "filmarray", "culture"):
        metric = revised_metrics["strict"][method]
        lines.append(
            f'| {method} | {metric["precision_ratio"]} = {metric["precision"]:.3f} | '
            f'{metric["recall_ratio"]} = {metric["recall"]:.3f} | {metric["f1"]:.3f} |'
        )
    lines.extend([
        "",
        "## Genus relaxed",
        "",
        "| Method | Precision | Recall | F1 |",
        "|---|---:|---:|---:|",
    ])
    for method in ("model_r5", "mngs", "filmarray", "culture"):
        metric = revised_metrics["genus_relaxed"][method]
        lines.append(
            f'| {method} | {metric["precision_ratio"]} = {metric["precision"]:.3f} | '
            f'{metric["recall_ratio"]} = {metric["recall"]:.3f} | {metric["f1"]:.3f} |'
        )
    lines.extend(["", "## Original baseline validation", ""])
    for mode, by_method in validation.items():
        for method, item in by_method.items():
            lines.append(f'- `{mode}/{method}`: {"PASS" if item["passed"] else "FAIL"}')
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps({"validation": validation, "revised_metrics": revised_metrics}, indent=2))
    print(revised_answer_csv)
    print(report_path)
    print(summary_path)


if __name__ == "__main__":
    main()

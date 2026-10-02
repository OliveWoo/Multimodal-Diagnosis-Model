from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import mngs_big_agent as mngs_agent  # noqa: E402


DEFAULT_INPUT = Path("outputs") / "patient_info_filtered"
DEFAULT_OUTPUT_DIR = Path("outputs") / "mngs_max_variability"
DEFAULT_MODEL = "gpt-5"
DEFAULT_REPEATS = 3
DEFAULT_SUMMARY_MODE = "full"
DEFAULT_IGNORED_FIELD_NAMES = {
    "cross_module_reasoning",
    "decision_trace",
    "integrated_reasoning",
    "integration_trace",
    "why_picked",
    "data_gaps",
    "evidence",
    "priority_flags",
    "caution_flags",
    "applied_rules",
    "best_record_rule",
    "guardrail_rule",
    "level_rule",
    "support_modules",
    "filmarray_gmtest_note",
}
CAUTION_FLAG_ORDER = {
    "Possible_colonization": 1,
    "Possible_reactivation": 2,
    "Evidence_limited": 3,
    "Host_expansion_applied": 4,
    "Protected_low_read_pathogen": 5,
    "Not_recommended_as_sole_treatment_basis": 6,
}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the same mNGS_max prompt multiple times and compare material outputs."
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        type=Path,
        default=[DEFAULT_INPUT],
        help=f"Patient folders or parent folders (default: {DEFAULT_INPUT}).",
    )
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS)
    parser.add_argument(
        "--repeat-workers",
        type=int,
        default=1,
        help=(
            "Number of repeat runs to execute concurrently per patient. "
            "Use 1 for serial execution; use 3 with --repeats 3 to run all repeats in parallel."
        ),
    )
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--patients",
        nargs="*",
        help="Specific patient folder names or ids, for example NGS_patient_10_json or 10.",
    )
    parser.add_argument(
        "--summary-mode",
        choices=["full", "no_filmarray"],
        default=DEFAULT_SUMMARY_MODE,
    )
    parser.add_argument(
        "--ranked-mngs",
        type=Path,
        default=mngs_agent.DEFAULT_RANKED_MNGS_PATH,
        help=f"Ranked mNGS JSON (default: {mngs_agent.DEFAULT_RANKED_MNGS_PATH}).",
    )
    parser.add_argument("--diff-limit", type=int, default=200)
    parser.add_argument(
        "--strict-all-fields",
        action="store_true",
        help="Treat reasoning/data_gaps/evidence text as material differences.",
    )
    args = parser.parse_args()

    if args.repeats < 2:
        parser.error("--repeats must be >= 2")
    if args.repeat_workers < 1:
        parser.error("--repeat-workers must be >= 1")

    patient_dirs = select_patients(
        args.inputs,
        patients=args.patients,
        limit=args.limit,
        summary_mode=args.summary_mode,
    )
    run_id = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_dir = args.output_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    ignored = set() if args.strict_all_fields else DEFAULT_IGNORED_FIELD_NAMES
    report: dict[str, Any] = {
        "check_name": "mngs_max_prompt_variability",
        "run_id": run_id,
        "model": args.model,
        "summary_mode": args.summary_mode,
        "repeats": args.repeats,
        "repeat_workers": args.repeat_workers,
        "api_key_fingerprint": api_key_fingerprint(),
        "patient_count": len(patient_dirs),
        "patients_with_all_runs_equal": 0,
        "patients_with_variability": 0,
        "patients_with_ignored_only_variability": 0,
        "patients_with_errors": 0,
        "material_comparison": {
            "ignored_field_names": sorted(ignored),
            "pass_fail_basis": "all_fields" if args.strict_all_fields else "material_fields_only",
        },
        "results": [],
    }

    prompt_path = mngs_agent.SUMMARY_MODE_TO_PROMPT_PATH[args.summary_mode]
    template = mngs_agent.load_prompt(prompt_path)
    for patient_dir in patient_dirs:
        result = run_patient(
            patient_dir,
            template=template,
            prompt_path=prompt_path,
            ranked_mngs_path=args.ranked_mngs,
            repeats=args.repeats,
            repeat_workers=args.repeat_workers,
            model=args.model,
            summary_mode=args.summary_mode,
            run_dir=run_dir,
            diff_limit=args.diff_limit,
            ignored_field_names=ignored,
        )
        report["results"].append(result)
        if result["error_count"]:
            report["patients_with_errors"] += 1
        elif result["all_runs_equal"]:
            report["patients_with_all_runs_equal"] += 1
            if result.get("ignored_changed_paths"):
                report["patients_with_ignored_only_variability"] += 1
        else:
            report["patients_with_variability"] += 1

    report_path = run_dir / "variability_report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"run_id={run_id}")
    print(f"patients={report['patient_count']}")
    print(f"repeats={args.repeats}")
    print(f"repeat_workers={args.repeat_workers}")
    print(f"summary_mode={args.summary_mode}")
    print(f"api_key={report['api_key_fingerprint']}")
    print(f"all_equal={report['patients_with_all_runs_equal']}")
    print(f"variable={report['patients_with_variability']}")
    print(f"ignored_only_variable={report['patients_with_ignored_only_variability']}")
    print(f"errors={report['patients_with_errors']}")
    print(f"report={report_path}")
    return 1 if report["patients_with_variability"] or report["patients_with_errors"] else 0


def select_patients(
    inputs: list[Path],
    *,
    patients: list[str] | None,
    limit: int | None,
    summary_mode: str,
) -> list[Path]:
    patient_dirs = mngs_agent.collect_patient_directories(inputs, summary_mode=summary_mode)
    if patients:
        requested = {normalize_patient_selector(value) for value in patients}
        patient_dirs = [
            path
            for path in patient_dirs
            if normalize_patient_selector(path.name) in requested
            or normalize_patient_selector(mngs_agent.extract_patient_identifier(path)) in requested
        ]
    if limit is not None:
        patient_dirs = patient_dirs[:limit]
    return patient_dirs


def run_patient(
    patient_dir: Path,
    *,
    template: str,
    prompt_path: Path,
    ranked_mngs_path: Path | None,
    repeats: int,
    repeat_workers: int,
    model: str,
    summary_mode: str,
    run_dir: Path,
    diff_limit: int,
    ignored_field_names: set[str],
) -> dict[str, Any]:
    patient_id = mngs_agent.extract_patient_identifier(patient_dir)
    patient_run_dir = run_dir / patient_id
    patient_run_dir.mkdir(parents=True, exist_ok=True)

    result: dict[str, Any] = {
        "patient_id": patient_id,
        "patient_dir": str(patient_dir),
        "summary_mode": summary_mode,
        "prompt_path": str(prompt_path),
        "prompt_sha256": "",
        "ranked_mngs_included": False,
        "all_runs_equal": False,
        "all_runs_exact_equal": False,
        "run_count": repeats,
        "repeat_workers": min(repeat_workers, repeats),
        "parsed_count": 0,
        "error_count": 0,
        "run_hashes": [],
        "changed_paths": [],
        "material_changed_paths": [],
        "ignored_changed_paths": [],
        "errors": [],
    }

    final_summary_path = mngs_agent.find_final_summary_file(patient_dir, summary_mode=summary_mode)
    mngs_path = mngs_agent.find_mngs_name_reads_file(patient_dir)
    if final_summary_path is None:
        result["errors"].append("final_summary_missing")
        result["error_count"] = 1
        return result
    if mngs_path is None:
        result["errors"].append("mngs_name_reads_missing")
        result["error_count"] = 1
        return result

    final_summary = mngs_agent.load_json(final_summary_path)
    mngs_grouped = mngs_agent.load_json(mngs_path)
    ranked_mngs = mngs_agent.load_ranked_mngs_for_patient(ranked_mngs_path, patient_dir)
    result["ranked_mngs_included"] = ranked_mngs is not None
    prompt = mngs_agent.build_prompt(
        template,
        final_summary=final_summary,
        mngs_grouped=mngs_grouped,
        final_path=final_summary_path,
        mngs_path=mngs_path,
        ranked_mngs=ranked_mngs,
        ranked_path=ranked_mngs_path,
    )
    result["prompt_sha256"] = sha256_text(prompt)

    parsed_runs: list[dict[str, Any]] = []
    repeat_results = execute_repeats(
        prompt,
        model=model,
        repeats=repeats,
        repeat_workers=repeat_workers,
        patient_run_dir=patient_run_dir,
        ignored_field_names=ignored_field_names,
        ranked_mngs=ranked_mngs,
        mngs_grouped=mngs_grouped,
        final_summary=final_summary,
    )
    for repeat_result in repeat_results:
        if error := repeat_result.get("error"):
            result["errors"].append(str(error))
            result["error_count"] += 1
            continue
        parsed = repeat_result["parsed"]
        parsed_runs.append(parsed)
        result["parsed_count"] += 1
        result["run_hashes"].append(repeat_result["run_hash"])

    if len(parsed_runs) >= 2:
        baseline = parsed_runs[0]
        material_baseline = normalize_for_material_compare(
            strip_ignored_fields(baseline, ignored_field_names)
        )
        material_changed = set()
        ignored_changed = set()
        exact_changed = set()
        for parsed in parsed_runs[1:]:
            material_parsed = normalize_for_material_compare(
                strip_ignored_fields(parsed, ignored_field_names)
            )
            material, ignored = diff_paths(
                material_baseline,
                material_parsed,
                "$",
                diff_limit,
                ignored_field_names=set(),
            )
            material_changed.update(material)
            _, ignored = diff_paths(
                baseline,
                parsed,
                "$",
                diff_limit,
                ignored_field_names=ignored_field_names,
            )
            ignored_changed.update(ignored)
            exact_changed.update(diff_paths_all(baseline, parsed, "$", diff_limit))
            if len(material_changed) >= diff_limit:
                break
        result["material_changed_paths"] = sorted(material_changed)[:diff_limit]
        result["ignored_changed_paths"] = sorted(ignored_changed)[:diff_limit]
        result["changed_paths"] = result["material_changed_paths"]
        result["all_runs_equal"] = (
            not result["material_changed_paths"] and result["error_count"] == 0
        )
        result["all_runs_exact_equal"] = not exact_changed and result["error_count"] == 0

    return result


def execute_repeats(
    prompt: str,
    *,
    model: str,
    repeats: int,
    repeat_workers: int,
    patient_run_dir: Path,
    ignored_field_names: set[str],
    ranked_mngs: dict[str, Any] | None,
    mngs_grouped: dict[str, Any],
    final_summary: dict[str, Any],
) -> list[dict[str, Any]]:
    workers = min(repeat_workers, repeats)
    indexes = list(range(1, repeats + 1))
    if workers == 1:
        return [
            execute_single_repeat(
                index,
                prompt=prompt,
                model=model,
                patient_run_dir=patient_run_dir,
                ignored_field_names=ignored_field_names,
                ranked_mngs=ranked_mngs,
                mngs_grouped=mngs_grouped,
                final_summary=final_summary,
            )
            for index in indexes
        ]

    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                execute_single_repeat,
                index,
                prompt=prompt,
                model=model,
                patient_run_dir=patient_run_dir,
                ignored_field_names=ignored_field_names,
                ranked_mngs=ranked_mngs,
                mngs_grouped=mngs_grouped,
                final_summary=final_summary,
            ): index
            for index in indexes
        }
        for future in as_completed(futures):
            index = futures[future]
            try:
                results.append(future.result())
            except Exception as exc:  # pragma: no cover - operational CLI behavior
                results.append(
                    {
                        "index": index,
                        "run": f"run_{index:02d}",
                        "error": f"run_{index:02d}: unexpected_worker_failed: {exc}",
                    }
                )
    return sorted(results, key=lambda item: item["index"])


def execute_single_repeat(
    index: int,
    *,
    prompt: str,
    model: str,
    patient_run_dir: Path,
    ignored_field_names: set[str],
    ranked_mngs: dict[str, Any] | None,
    mngs_grouped: dict[str, Any],
    final_summary: dict[str, Any],
) -> dict[str, Any]:
    run_name = f"run_{index:02d}"
    raw_path = patient_run_dir / f"{run_name}.raw.txt"
    json_path = patient_run_dir / f"{run_name}.json"
    try:
        raw, usage = mngs_agent.send_to_llm(prompt, model=model)
    except Exception as exc:  # pragma: no cover - operational CLI behavior
        return {
            "index": index,
            "run": run_name,
            "error": f"{run_name}: llm_request_failed: {exc}",
        }

    raw_path.write_text(raw, encoding="utf-8")
    parsed = mngs_agent.validate_json(raw)
    if parsed is None:
        return {
            "index": index,
            "run": run_name,
            "error": f"{run_name}: json_parse_failed",
        }
    parsed = mngs_agent.postprocess_mngs_max_output(
        parsed,
        ranked_mngs=ranked_mngs,
        mngs_grouped=mngs_grouped,
        final_summary=final_summary,
    )

    canonical = canonical_json(parsed)
    material_value = normalize_for_material_compare(
        strip_ignored_fields(parsed, ignored_field_names)
    )
    json_path.write_text(json.dumps(parsed, ensure_ascii=False, indent=2), encoding="utf-8")
    return {
        "index": index,
        "run": run_name,
        "parsed": parsed,
        "run_hash": {
            "run": run_name,
            "sha256": sha256_text(canonical),
            "material_sha256": sha256_text(canonical_json(material_value)),
            "usage": usage,
            "path": str(json_path),
        },
    }


def diff_paths(
    left: Any,
    right: Any,
    path: str,
    limit: int,
    *,
    ignored_field_names: set[str],
) -> tuple[list[str], list[str]]:
    material_diffs: list[str] = []
    ignored_diffs: list[str] = []

    def add_material(value: str) -> None:
        if len(material_diffs) < limit:
            material_diffs.append(value)

    def add_ignored(value: str) -> None:
        if len(ignored_diffs) < limit:
            ignored_diffs.append(value)

    def walk(a: Any, b: Any, current: str) -> None:
        if len(material_diffs) >= limit:
            return
        if type(a) is not type(b):
            add_material(current)
            return
        if isinstance(a, dict):
            keys = set(a) | set(b)
            for key in sorted(keys):
                next_path = f"{current}.{key}"
                if key in ignored_field_names:
                    left_value = a.get(key, _MISSING)
                    right_value = b.get(key, _MISSING)
                    if left_value != right_value:
                        add_ignored(next_path)
                    continue
                if key not in a or key not in b:
                    add_material(next_path)
                else:
                    walk(a[key], b[key], next_path)
                if len(material_diffs) >= limit:
                    return
            return
        if isinstance(a, list):
            if len(a) != len(b):
                add_material(f"{current}.length")
            for index, (left_item, right_item) in enumerate(zip(a, b)):
                walk(left_item, right_item, f"{current}[{index}]")
                if len(material_diffs) >= limit:
                    return
            return
        if a != b:
            add_material(current)

    walk(left, right, path)
    return material_diffs, ignored_diffs


def diff_paths_all(left: Any, right: Any, path: str, limit: int) -> list[str]:
    diffs: list[str] = []

    def add(value: str) -> None:
        if len(diffs) < limit:
            diffs.append(value)

    def walk(a: Any, b: Any, current: str) -> None:
        if len(diffs) >= limit:
            return
        if type(a) is not type(b):
            add(current)
            return
        if isinstance(a, dict):
            keys = set(a) | set(b)
            for key in sorted(keys):
                next_path = f"{current}.{key}"
                if key not in a or key not in b:
                    add(next_path)
                else:
                    walk(a[key], b[key], next_path)
                if len(diffs) >= limit:
                    return
            return
        if isinstance(a, list):
            if len(a) != len(b):
                add(f"{current}.length")
            for index, (left_item, right_item) in enumerate(zip(a, b)):
                walk(left_item, right_item, f"{current}[{index}]")
                if len(diffs) >= limit:
                    return
            return
        if a != b:
            add(current)

    walk(left, right, path)
    return diffs


def strip_ignored_fields(value: Any, ignored_field_names: set[str]) -> Any:
    if isinstance(value, dict):
        return {
            key: strip_ignored_fields(item, ignored_field_names)
            for key, item in value.items()
            if key not in ignored_field_names
        }
    if isinstance(value, list):
        return [strip_ignored_fields(item, ignored_field_names) for item in value]
    return value


def normalize_for_material_compare(value: Any) -> Any:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return normalize_number(value)
    if isinstance(value, str):
        return normalize_numeric_string(value)
    if isinstance(value, dict):
        normalized = {
            key: normalize_for_material_compare(item)
            for key, item in value.items()
        }
        if is_lower_priority_candidate(normalized):
            normalized.pop("module_support_summary", None)
        for key in (
            "pathogen_candidates",
            "excluded_candidates",
            "picked_pathogens",
            "not_picked_pathogens",
        ):
            items = normalized.get(key)
            if isinstance(items, list):
                normalized[key] = sorted(items, key=material_candidate_sort_key)
        flags = normalized.get("caution_flags")
        if isinstance(flags, list):
            normalized["caution_flags"] = sorted(
                flags,
                key=lambda item: (CAUTION_FLAG_ORDER.get(str(item), 999), str(item)),
            )
        return normalized
    if isinstance(value, list):
        return [normalize_for_material_compare(item) for item in value]
    return value


def normalize_number(value: int | float) -> int | float:
    if isinstance(value, float):
        if not math.isfinite(value):
            return value
        if value.is_integer():
            return int(value)
    return value


def normalize_numeric_string(value: str) -> str | int | float:
    text = value.strip()
    if not text:
        return value
    try:
        number = float(text)
    except ValueError:
        return value
    if not math.isfinite(number):
        return value
    return normalize_number(number)


def is_lower_priority_candidate(value: dict[str, Any]) -> bool:
    if "module_support_summary" not in value:
        return False
    rank = level_rank(
        value.get("integrated_causative_level")
        or value.get("observed_level")
        or value.get("basis_level")
    )
    return rank >= 4


def material_candidate_sort_key(item: Any) -> tuple[Any, ...]:
    if not isinstance(item, dict):
        return (999, "", canonical_json(item))
    return (
        level_rank(
            item.get("integrated_causative_level")
            or item.get("observed_best_level")
            or item.get("causative_level")
            or item.get("level")
        ),
        normalize_sort_text(
            item.get("organism_name")
            or item.get("name")
            or item.get("pathogen_name")
        ),
        normalize_sort_text(item.get("classification") or item.get("type")),
        normalize_sort_text(item.get("exclusion_reason_code") or item.get("reason")),
        canonical_json(item),
    )


def level_rank(value: Any) -> int:
    text = str(value or "").strip().lower()
    for index in range(1, 6):
        if text == f"level {index}" or text == f"level{index}":
            return index
    if text in {"not_available", "not available"}:
        return 98
    return 99


def normalize_sort_text(value: Any) -> str:
    return "".join(ch for ch in str(value or "").lower() if ch.isalnum())


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def normalize_patient_selector(value: str) -> str:
    value = value.strip().removesuffix("_json")
    return value.removeprefix("NGS_patient_")


def api_key_fingerprint() -> dict[str, Any]:
    value = os.environ.get("OPENAI_API_KEY") or ""
    if not value:
        return {"is_set": False, "prefix": "", "suffix": "", "length": 0}
    return {
        "is_set": True,
        "prefix": value[:12],
        "suffix": value[-4:],
        "length": len(value),
    }


class _Missing:
    pass


_MISSING = _Missing()


if __name__ == "__main__":
    sys.exit(main())

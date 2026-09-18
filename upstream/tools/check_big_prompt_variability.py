from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from core.summarize_agent import (  # noqa: E402
    AgentOutputMissingError,
    build_prompt,
    collect_agent_outputs,
    collect_patient_directories,
    extract_patient_identifier,
    load_big_prompt,
    send_to_llm,
    validate_json,
)


DEFAULT_INPUT = Path("outputs") / "patient_info_filtered"
DEFAULT_OUTPUT_DIR = Path("outputs") / "big_prompt_variability"
DEFAULT_MODEL = "gpt-5"
DEFAULT_REPEATS = 3
DEFAULT_IGNORED_FIELD_NAMES = {
    "cross_module_reasoning",
    "decision_trace",
    "integrated_reasoning",
    "integration_trace",
    "reasoning",
    "summary_reasoning",
    "data_gaps",
    "priority_flags",
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
        description=(
            "Run the same big_prompt with-filmarray input multiple times and compare outputs."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        type=Path,
        default=[DEFAULT_INPUT],
        help=f"Patient folders or parent folders (default: {DEFAULT_INPUT}).",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=DEFAULT_REPEATS,
        help=f"Number of LLM runs per patient (default: {DEFAULT_REPEATS}).",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Model to use (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Where to store run outputs and report (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Only test the first N matching patient folders.",
    )
    parser.add_argument(
        "--patients",
        nargs="*",
        help=(
            "Specific patient folder names or patient ids to test, "
            "for example NGS_patient_10_json or NGS_patient_10."
        ),
    )
    parser.add_argument(
        "--diff-limit",
        type=int,
        default=200,
        help="Maximum diff paths stored per patient (default: 200).",
    )
    parser.add_argument(
        "--strict-all-fields",
        action="store_true",
        help=(
            "Treat narrative fields such as reasoning/data_gaps as material "
            "differences. By default they are reported separately and ignored "
            "for pass/fail."
        ),
    )
    args = parser.parse_args()

    if args.repeats < 2:
        parser.error("--repeats must be >= 2")

    patient_dirs = select_patients(args.inputs, patients=args.patients, limit=args.limit)
    run_id = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_dir = args.output_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    report = {
        "check_name": "big_prompt_with_filmarray_variability",
        "run_id": run_id,
        "model": args.model,
        "repeats": args.repeats,
        "api_key_fingerprint": api_key_fingerprint(),
        "patient_count": len(patient_dirs),
        "patients_with_all_runs_equal": 0,
        "patients_with_variability": 0,
        "patients_with_ignored_only_variability": 0,
        "patients_with_errors": 0,
        "material_comparison": {
            "ignored_field_names": (
                [] if args.strict_all_fields else sorted(DEFAULT_IGNORED_FIELD_NAMES)
            ),
            "pass_fail_basis": (
                "all_fields" if args.strict_all_fields else "material_fields_only"
            ),
        },
        "results": [],
    }

    template = load_big_prompt()
    for patient_dir in patient_dirs:
        result = run_patient(
            patient_dir,
            template=template,
            repeats=args.repeats,
            model=args.model,
            run_dir=run_dir,
            diff_limit=args.diff_limit,
            ignored_field_names=(
                set() if args.strict_all_fields else DEFAULT_IGNORED_FIELD_NAMES
            ),
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
) -> list[Path]:
    patient_dirs = collect_patient_directories(inputs)
    if patients:
        requested = {normalize_patient_selector(value) for value in patients}
        patient_dirs = [
            path
            for path in patient_dirs
            if normalize_patient_selector(path.name) in requested
            or normalize_patient_selector(path.name.removesuffix("_json")) in requested
        ]
    if limit is not None:
        patient_dirs = patient_dirs[:limit]
    return patient_dirs


def run_patient(
    patient_dir: Path,
    *,
    template: str,
    repeats: int,
    model: str,
    run_dir: Path,
    diff_limit: int,
    ignored_field_names: set[str],
) -> dict[str, Any]:
    patient_id = extract_patient_identifier(patient_dir)
    patient_run_dir = run_dir / patient_id
    patient_run_dir.mkdir(parents=True, exist_ok=True)

    result: dict[str, Any] = {
        "patient_id": patient_id,
        "patient_dir": str(patient_dir),
        "prompt_sha256": "",
        "all_runs_equal": False,
        "all_runs_exact_equal": False,
        "run_count": repeats,
        "parsed_count": 0,
        "error_count": 0,
        "run_hashes": [],
        "changed_paths": [],
        "material_changed_paths": [],
        "ignored_changed_paths": [],
        "errors": [],
    }

    try:
        sections = collect_agent_outputs(
            patient_dir,
            include_underlying=False,
        )
    except AgentOutputMissingError as exc:
        result["errors"].append(str(exc))
        result["error_count"] = 1
        return result

    prompt = build_prompt(template, sections)
    result["prompt_sha256"] = sha256_text(prompt)

    parsed_runs: list[dict[str, Any]] = []
    for index in range(1, repeats + 1):
        run_name = f"run_{index:02d}"
        try:
            raw = send_to_llm(prompt, model=model)
        except Exception as exc:  # pragma: no cover - operational CLI behavior
            result["errors"].append(f"{run_name}: llm_request_failed: {exc}")
            result["error_count"] += 1
            continue

        (patient_run_dir / f"{run_name}.raw.txt").write_text(raw, encoding="utf-8")
        parsed = validate_json(raw)
        if parsed is None:
            result["errors"].append(f"{run_name}: json_parse_failed")
            result["error_count"] += 1
            continue

        parsed_runs.append(parsed)
        result["parsed_count"] += 1
        canonical = canonical_json(parsed)
        material_value = normalize_for_material_compare(
            strip_ignored_fields(parsed, ignored_field_names)
        )
        result["run_hashes"].append(
            {
                "run": run_name,
                "sha256": sha256_text(canonical),
                "material_sha256": sha256_text(canonical_json(material_value)),
                "path": str(patient_run_dir / f"{run_name}.json"),
            }
        )
        (patient_run_dir / f"{run_name}.json").write_text(
            json.dumps(parsed, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

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
    if isinstance(value, dict):
        normalized = {
            key: normalize_for_material_compare(item)
            for key, item in value.items()
        }
        for key in ("pathogen_candidates", "excluded_candidates"):
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


def material_candidate_sort_key(item: Any) -> tuple[Any, ...]:
    if not isinstance(item, dict):
        return (999, "", canonical_json(item))
    return (
        level_rank(
            item.get("integrated_causative_level")
            or item.get("observed_best_level")
            or item.get("causative_level")
        ),
        normalize_sort_text(item.get("organism_name")),
        normalize_sort_text(item.get("classification")),
        normalize_sort_text(item.get("exclusion_reason_code")),
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
    return value.strip().removesuffix("_json")


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

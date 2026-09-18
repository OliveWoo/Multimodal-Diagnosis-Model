from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any, Sequence

from tools import mngs_common as mngs


DEFAULT_PATIENT_ROOT = Path("outputs") / "patient_info_filtered"
DEFAULT_OUTPUT_JSON = Path("outputs") / "final_summary_llm_vs_deterministic_diff.json"
DEFAULT_OUTPUT_CSV = Path("outputs") / "final_summary_llm_vs_deterministic_diff.csv"

LLM_SUFFIX = "final_summary_with_filmarray"
DETERMINISTIC_SUFFIX = "final_summary_with_filmarray_deterministic"

MATERIAL_TOP_LEVEL_FIELDS = [
    "final_infection_likelihood",
    "dominant_source",
    "dominant_pathogen_type",
    "support_direction",
    "data_quality_tier",
]

HOST_FIELDS = [
    "host_vulnerability_tier",
    "expanded_candidate_policy",
    "opportunistic_coverage_level",
]

NAME_SANITIZER = re.compile(r"[^a-z0-9]+")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare LLM final_summary_with_filmarray with deterministic summary outputs."
    )
    parser.add_argument("patient_root", nargs="?", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--patients", nargs="*", help="Optional patient IDs.")
    parser.add_argument("--llm-suffix", default=LLM_SUFFIX)
    parser.add_argument("--deterministic-suffix", default=DETERMINISTIC_SUFFIX)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT_JSON)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT_CSV)
    return parser.parse_args(argv)


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def patient_id_from_dir(patient_dir: Path) -> str:
    parts = patient_dir.name.split("_")
    if len(parts) >= 3 and parts[0] == "NGS" and parts[1] == "patient":
        return parts[2]
    raise ValueError(f"Cannot parse patient ID from {patient_dir}")


def patient_sort_key(patient_dir: Path) -> tuple[int, str]:
    patient_id = patient_id_from_dir(patient_dir)
    return (int(patient_id), patient_id) if patient_id.isdigit() else (999999, patient_id)


def patient_dirs(patient_root: Path, patients: Sequence[str] | None) -> list[Path]:
    wanted = {
        str(item).removeprefix("NGS_patient_").removesuffix("_json")
        for item in patients or []
    }
    dirs = sorted(
        (path for path in patient_root.glob("NGS_patient_*_json") if path.is_dir()),
        key=patient_sort_key,
    )
    if wanted:
        dirs = [path for path in dirs if patient_id_from_dir(path) in wanted]
    return dirs


def summary_path(patient_dir: Path, suffix: str) -> Path:
    patient_id = patient_id_from_dir(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{patient_id}_{suffix}.json"


def normalize_scalar(value: Any) -> Any:
    if isinstance(value, str):
        text = value.strip()
        lowered = text.lower()
        if lowered in {"true", "false"}:
            return lowered == "true"
        return text
    return value


def normalize_field_value(field: str, value: Any) -> Any:
    normalized = normalize_scalar(value)
    if field == "dominant_source":
        lowered = str(normalized or "").strip().lower()
        if lowered in {"blood", "bloodstream", "systemic"}:
            return "Systemic"
    return normalized


def canonical_name(value: Any) -> str:
    text = str(value or "").strip().lower()
    if not text:
        return ""
    return mngs.normalize_organism_name(text)


def candidate_map(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    for item in (payload.get("pathogen_candidates", []) or []) + (
        payload.get("hospital_organism_evidence", []) or []
    ):
        if not isinstance(item, dict):
            continue
        key = canonical_name(item.get("organism_name"))
        if key:
            output[key] = item
    return output


def display_name_map(*maps: dict[str, dict[str, Any]]) -> dict[str, str]:
    output: dict[str, str] = {}
    for mapping in maps:
        for key, item in mapping.items():
            output.setdefault(key, str(item.get("organism_name") or key))
    return output


def field_mismatches(llm_payload: dict[str, Any], det_payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    mismatches: dict[str, dict[str, Any]] = {}
    for field in MATERIAL_TOP_LEVEL_FIELDS:
        llm_value = normalize_field_value(field, llm_payload.get(field))
        det_value = normalize_field_value(field, det_payload.get(field))
        if llm_value != det_value:
            mismatches[field] = {"llm": llm_value, "deterministic": det_value}

    llm_host = llm_payload.get("host_context") or {}
    det_host = det_payload.get("host_context") or {}
    for field in HOST_FIELDS:
        llm_value = normalize_scalar(llm_host.get(field))
        det_value = normalize_scalar(det_host.get(field))
        if llm_value != det_value:
            mismatches[f"host_context.{field}"] = {
                "llm": llm_value,
                "deterministic": det_value,
            }
    return mismatches


def module_availability_mismatches(
    llm_payload: dict[str, Any], det_payload: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    llm_modules = llm_payload.get("module_availability") or {}
    det_modules = det_payload.get("module_availability") or {}
    keys = sorted(set(llm_modules) | set(det_modules))
    mismatches: dict[str, dict[str, Any]] = {}
    for key in keys:
        llm_value = normalize_scalar(llm_modules.get(key))
        det_value = normalize_scalar(det_modules.get(key))
        if llm_value != det_value:
            mismatches[key] = {"llm": llm_value, "deterministic": det_value}
    return mismatches


def candidate_level_mismatches(
    llm_candidates: dict[str, dict[str, Any]],
    det_candidates: dict[str, dict[str, Any]],
    names: dict[str, str],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for key in sorted(set(llm_candidates) & set(det_candidates), key=lambda item: names.get(item, item).lower()):
        llm_item = llm_candidates[key]
        det_item = det_candidates[key]
        llm_level = normalize_scalar(
            llm_item.get("integrated_causative_level", llm_item.get("best_hospital_level"))
        )
        det_level = normalize_scalar(
            det_item.get("integrated_causative_level", det_item.get("best_hospital_level"))
        )
        llm_class = normalize_scalar(llm_item.get("classification"))
        det_class = normalize_scalar(det_item.get("classification"))
        if llm_level != det_level or llm_class != det_class:
            output.append(
                {
                    "organism_name": names.get(key, key),
                    "llm_level": llm_level,
                    "deterministic_level": det_level,
                    "llm_classification": llm_class,
                    "deterministic_classification": det_class,
                }
            )
    return output


def candidate_detail(item: dict[str, Any] | None) -> str:
    if not item:
        return ""
    modules = item.get("evidence_modules") or []
    if isinstance(modules, list):
        module_text = "+".join(str(module) for module in modules)
    else:
        module_text = str(modules)
    return (
        f"{item.get('organism_name', '')}"
        f"|level={item.get('integrated_causative_level', item.get('best_hospital_level', ''))}"
        f"|class={item.get('classification', '')}"
        f"|modules={module_text}"
    )


def compare_patient(patient_dir: Path, llm_suffix: str, deterministic_suffix: str) -> dict[str, Any]:
    patient_id = patient_id_from_dir(patient_dir)
    llm_path = summary_path(patient_dir, llm_suffix)
    det_path = summary_path(patient_dir, deterministic_suffix)
    if not llm_path.exists() or not det_path.exists():
        return {
            "patient_id": patient_id,
            "status": "missing_file",
            "llm_path": str(llm_path),
            "deterministic_path": str(det_path),
            "missing": [
                name
                for name, path in (("llm", llm_path), ("deterministic", det_path))
                if not path.exists()
            ],
        }

    llm_payload = load_json(llm_path)
    det_payload = load_json(det_path)
    llm_candidates = candidate_map(llm_payload)
    det_candidates = candidate_map(det_payload)
    names = display_name_map(llm_candidates, det_candidates)
    llm_keys = set(llm_candidates)
    det_keys = set(det_candidates)
    both = sorted(llm_keys & det_keys, key=lambda key: names.get(key, key).lower())
    llm_only = sorted(llm_keys - det_keys, key=lambda key: names.get(key, key).lower())
    det_only = sorted(det_keys - llm_keys, key=lambda key: names.get(key, key).lower())

    top_level_mismatches = field_mismatches(llm_payload, det_payload)
    availability_mismatches = module_availability_mismatches(llm_payload, det_payload)
    level_mismatches = candidate_level_mismatches(llm_candidates, det_candidates, names)

    return {
        "patient_id": patient_id,
        "status": "ok",
        "llm_path": str(llm_path),
        "deterministic_path": str(det_path),
        "top_level_mismatches": top_level_mismatches,
        "module_availability_mismatches": availability_mismatches,
        "llm_candidate_count": len(llm_keys),
        "deterministic_candidate_count": len(det_keys),
        "candidate_overlap_count": len(both),
        "llm_only_count": len(llm_only),
        "deterministic_only_count": len(det_only),
        "candidate_level_mismatch_count": len(level_mismatches),
        "both_candidates": [names.get(key, key) for key in both],
        "llm_only_candidates": [candidate_detail(llm_candidates.get(key)) for key in llm_only],
        "deterministic_only_candidates": [candidate_detail(det_candidates.get(key)) for key in det_only],
        "candidate_level_mismatches": level_mismatches,
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    rows = [
        compare_patient(patient_dir, args.llm_suffix, args.deterministic_suffix)
        for patient_dir in patient_dirs(args.patient_root, args.patients)
    ]
    ok_rows = [row for row in rows if row["status"] == "ok"]
    summary = {
        "patients_scanned": len(rows),
        "patients_compared": len(ok_rows),
        "patients_missing_files": sum(1 for row in rows if row["status"] != "ok"),
        "patients_with_top_level_mismatch": sum(1 for row in ok_rows if row["top_level_mismatches"]),
        "patients_with_module_availability_mismatch": sum(
            1 for row in ok_rows if row["module_availability_mismatches"]
        ),
        "patients_with_candidate_set_difference": sum(
            1 for row in ok_rows if row["llm_only_count"] or row["deterministic_only_count"]
        ),
        "patients_with_candidate_level_mismatch": sum(
            1 for row in ok_rows if row["candidate_level_mismatch_count"]
        ),
        "total_llm_candidates": sum(row["llm_candidate_count"] for row in ok_rows),
        "total_deterministic_candidates": sum(row["deterministic_candidate_count"] for row in ok_rows),
        "total_candidate_overlap": sum(row["candidate_overlap_count"] for row in ok_rows),
        "total_llm_only_candidates": sum(row["llm_only_count"] for row in ok_rows),
        "total_deterministic_only_candidates": sum(row["deterministic_only_count"] for row in ok_rows),
        "total_candidate_level_mismatches": sum(
            row["candidate_level_mismatch_count"] for row in ok_rows
        ),
    }
    report = {
        "patient_root": str(args.patient_root),
        "llm_suffix": args.llm_suffix,
        "deterministic_suffix": args.deterministic_suffix,
        "comparison_note": (
            "For deterministic v2 summaries, organism comparison uses hospital_organism_evidence "
            "because summary no longer decides pathogen_candidates/excluded_candidates."
        ),
        "summary": summary,
        "rows": rows,
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "patient_id",
                "status",
                "top_level_mismatches",
                "module_availability_mismatches",
                "llm_candidate_count",
                "deterministic_candidate_count",
                "candidate_overlap_count",
                "llm_only_candidates",
                "deterministic_only_candidates",
                "candidate_level_mismatches",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "patient_id": row.get("patient_id", ""),
                    "status": row.get("status", ""),
                    "top_level_mismatches": json.dumps(
                        row.get("top_level_mismatches", {}), ensure_ascii=False
                    ),
                    "module_availability_mismatches": json.dumps(
                        row.get("module_availability_mismatches", {}), ensure_ascii=False
                    ),
                    "llm_candidate_count": row.get("llm_candidate_count", ""),
                    "deterministic_candidate_count": row.get("deterministic_candidate_count", ""),
                    "candidate_overlap_count": row.get("candidate_overlap_count", ""),
                    "llm_only_candidates": "; ".join(row.get("llm_only_candidates", [])),
                    "deterministic_only_candidates": "; ".join(
                        row.get("deterministic_only_candidates", [])
                    ),
                    "candidate_level_mismatches": json.dumps(
                        row.get("candidate_level_mismatches", []), ensure_ascii=False
                    ),
                }
            )

    print(f"Wrote JSON report: {args.output_json}")
    print(f"Wrote CSV report: {args.output_csv}")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

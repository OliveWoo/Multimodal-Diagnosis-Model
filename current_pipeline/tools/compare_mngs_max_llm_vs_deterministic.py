from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any, Sequence

from tools import mngs_big_agent as mngs


DEFAULT_PATIENT_ROOT = Path("outputs") / "patient_info_rich_normalized"
DEFAULT_OUTPUT_JSON = Path("outputs") / "mngs_max_llm_vs_deterministic_diff.json"
DEFAULT_OUTPUT_CSV = Path("outputs") / "mngs_max_llm_vs_deterministic_diff.csv"

COMPARISONS = {
    "full_ranked": {
        "llm_suffix": "mNGS_max_agent_full",
        "deterministic_suffix": "mNGS_max_deterministic_full_ranked",
    },
    "opt_chosen": {
        "llm_suffix": "mNGS_max_agent_all_rk_ntc_opt_chosen_full",
        "deterministic_suffix": "mNGS_max_deterministic_opt_chosen_full",
    },
}

NAME_SANITIZER = re.compile(r"[^a-z0-9]+")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare old LLM mNGS max picked pathogens with deterministic scorer outputs."
    )
    parser.add_argument("patient_root", nargs="?", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--patients", nargs="*", help="Optional patient IDs.")
    parser.add_argument(
        "--comparison",
        choices=["all", *COMPARISONS.keys()],
        default="all",
        help="Which pair to compare.",
    )
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


def patient_dirs(patient_root: Path, patients: Sequence[str] | None) -> list[Path]:
    wanted = {str(item).removeprefix("NGS_patient_").removesuffix("_json") for item in patients or []}
    dirs = sorted(
        (path for path in patient_root.glob("NGS_patient_*_json") if path.is_dir()),
        key=lambda path: int(patient_id_from_dir(path)),
    )
    if wanted:
        dirs = [path for path in dirs if patient_id_from_dir(path) in wanted]
    return dirs


def canonical_name(value: Any) -> str:
    text = str(value or "").strip().lower()
    if not text:
        return ""
    return mngs.normalize_organism_name(text)


def summary_path(patient_dir: Path, suffix: str) -> Path:
    patient_id = patient_id_from_dir(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{patient_id}_{suffix}.json"


def find_summary_path(patient_dir: Path, suffix: str) -> Path:
    path = summary_path(patient_dir, suffix)
    if path.exists():
        return path
    archived = path.parent / "archive_llm_baseline" / path.name
    if archived.exists():
        return archived
    return path


def picked_map(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    picked = payload.get("best_available_summary", {}).get("picked_pathogens", [])
    output: dict[str, dict[str, Any]] = {}
    for item in picked:
        if not isinstance(item, dict):
            continue
        key = canonical_name(item.get("organism_name"))
        if key:
            output[key] = item
    return output


def candidate_map(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    for item in payload.get("pathogen_candidates", []) or []:
        if not isinstance(item, dict):
            continue
        key = canonical_name(item.get("organism_name"))
        if key:
            output[key] = item
    return output


def excluded_map(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    output: dict[str, dict[str, Any]] = {}
    for item in payload.get("excluded_candidates", []) or []:
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


def deterministic_non_pick_reason(key: str, det_payload: dict[str, Any]) -> str:
    candidates = candidate_map(det_payload)
    excluded = excluded_map(det_payload)
    candidate = candidates.get(key)
    excluded_candidate = excluded.get(key)
    if candidate is None:
        return "not_in_deterministic_candidate_pool"
    level = candidate.get("integrated_causative_level", "Unknown")
    signal = candidate.get("mngs_signal_tier", "Unknown")
    guardrail = (candidate.get("key_evidence") or {}).get("guardrail_rule") or ""
    if candidate.get("is_likely_colonizer_or_background") is True:
        return f"deterministic_guardrail_background_or_colonizer; level={level}; signal={signal}; guardrail={guardrail}"
    if excluded_candidate:
        code = excluded_candidate.get("exclusion_reason_code", "excluded")
        return f"deterministic_excluded; code={code}; level={level}; signal={signal}; guardrail={guardrail}"
    return f"deterministic_not_picked; level={level}; signal={signal}; guardrail={guardrail}"


def picked_detail(item: dict[str, Any] | None) -> str:
    if not item:
        return ""
    return (
        f"{item.get('organism_name', '')}"
        f"|level={item.get('basis_level', item.get('integrated_causative_level', ''))}"
        f"|signal={item.get('mngs_signal_tier', '')}"
        f"|rank={item.get('rank_priority', '')}"
        f"|reads={item.get('reads', '')}"
    )


def compare_pair(patient_dir: Path, label: str, llm_suffix: str, deterministic_suffix: str) -> dict[str, Any]:
    llm_path = find_summary_path(patient_dir, llm_suffix)
    det_path = find_summary_path(patient_dir, deterministic_suffix)
    patient_id = patient_id_from_dir(patient_dir)
    if not llm_path.exists() or not det_path.exists():
        return {
            "patient_id": patient_id,
            "comparison": label,
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
    llm_picked = picked_map(llm_payload)
    det_picked = picked_map(det_payload)
    names = display_name_map(llm_picked, det_picked, candidate_map(det_payload))

    llm_keys = set(llm_picked)
    det_keys = set(det_picked)
    both = sorted(llm_keys & det_keys, key=lambda key: names.get(key, key).lower())
    llm_only = sorted(llm_keys - det_keys, key=lambda key: names.get(key, key).lower())
    det_only = sorted(det_keys - llm_keys, key=lambda key: names.get(key, key).lower())

    return {
        "patient_id": patient_id,
        "comparison": label,
        "status": "ok",
        "llm_path": str(llm_path),
        "deterministic_path": str(det_path),
        "llm_count": len(llm_keys),
        "deterministic_count": len(det_keys),
        "both_count": len(both),
        "llm_only_count": len(llm_only),
        "deterministic_only_count": len(det_only),
        "both": [names.get(key, key) for key in both],
        "llm_only": [
            {
                "organism_name": names.get(key, key),
                "llm_detail": picked_detail(llm_picked.get(key)),
                "deterministic_non_pick_reason": deterministic_non_pick_reason(key, det_payload),
            }
            for key in llm_only
        ],
        "deterministic_only": [
            {
                "organism_name": names.get(key, key),
                "deterministic_detail": picked_detail(det_picked.get(key)),
            }
            for key in det_only
        ],
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    selected = COMPARISONS if args.comparison == "all" else {args.comparison: COMPARISONS[args.comparison]}
    rows: list[dict[str, Any]] = []
    for patient_dir in patient_dirs(args.patient_root, args.patients):
        for label, config in selected.items():
            rows.append(
                compare_pair(
                    patient_dir,
                    label,
                    config["llm_suffix"],
                    config["deterministic_suffix"],
                )
            )

    summary: dict[str, Any] = {}
    for label in selected:
        ok_rows = [row for row in rows if row["comparison"] == label and row["status"] == "ok"]
        summary[label] = {
            "patients_compared": len(ok_rows),
            "patients_missing_files": sum(
                1 for row in rows if row["comparison"] == label and row["status"] != "ok"
            ),
            "total_llm_picked": sum(row["llm_count"] for row in ok_rows),
            "total_deterministic_picked": sum(row["deterministic_count"] for row in ok_rows),
            "total_overlap": sum(row["both_count"] for row in ok_rows),
            "total_llm_only": sum(row["llm_only_count"] for row in ok_rows),
            "total_deterministic_only": sum(row["deterministic_only_count"] for row in ok_rows),
        }

    report = {
        "patient_root": str(args.patient_root),
        "comparisons": selected,
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
                "comparison",
                "status",
                "llm_count",
                "deterministic_count",
                "both",
                "llm_only",
                "llm_only_reasons",
                "deterministic_only",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "patient_id": row.get("patient_id", ""),
                    "comparison": row.get("comparison", ""),
                    "status": row.get("status", ""),
                    "llm_count": row.get("llm_count", ""),
                    "deterministic_count": row.get("deterministic_count", ""),
                    "both": "; ".join(row.get("both", [])),
                    "llm_only": "; ".join(item["organism_name"] for item in row.get("llm_only", [])),
                    "llm_only_reasons": "; ".join(
                        f"{item['organism_name']} => {item['deterministic_non_pick_reason']}"
                        for item in row.get("llm_only", [])
                    ),
                    "deterministic_only": "; ".join(
                        item["organism_name"] for item in row.get("deterministic_only", [])
                    ),
                }
            )

    print(f"Wrote JSON report: {args.output_json}")
    print(f"Wrote CSV report: {args.output_csv}")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Extract a frozen, traceable pathogen-selection package from upstream outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PATIENT_ID_RE = re.compile(r"NGS_patient_(\d+)_", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract deterministic_max.best_available_summary.picked_pathogens "
            "and the corresponding upstream candidate list."
        )
    )
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--hospital", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--pipeline-version", default="upstream_frozen")
    parser.add_argument("--pretty", action="store_true")
    return parser.parse_args()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, value: Any, pretty: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        if pretty:
            json.dump(value, handle, ensure_ascii=False, indent=2)
        else:
            json.dump(value, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")


def prepare_output_root(path: Path) -> Path:
    resolved = path.resolve()
    if resolved.exists() and not resolved.is_dir():
        raise SystemExit(f"output root is not a directory: {resolved}")
    if resolved.is_dir() and any(resolved.iterdir()):
        raise SystemExit(f"output root must be new or empty: {resolved}")
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def canonical_key(name: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", compact(name).casefold())


def string_list(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    return [text for item in value if (text := compact(item))]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def patient_id_from(path: Path, data: dict[str, Any]) -> str:
    value = compact(data.get("patient_id"))
    if value:
        return value
    match = PATIENT_ID_RE.search(path.name)
    if not match:
        raise ValueError("patient_id is absent from both JSON and filename")
    return match.group(1)


def picked_entry(item: dict[str, Any], rank: int) -> dict[str, Any]:
    name = compact(item.get("organism_name"))
    if not name:
        raise ValueError(f"picked pathogen at rank {rank} has no organism_name")
    return {
        "original_name": name,
        "canonical_name": name,
        "canonical_key": compact(item.get("canonical_key")) or canonical_key(name),
        "selection_status": "picked",
        "rank": rank,
        "classification": item.get("classification"),
        "picked_role": item.get("picked_role"),
        "evidence_source": item.get("evidence_source"),
        "basis_level": item.get("basis_level"),
        "selection_reasons": string_list(item.get("why_picked")),
        "caution_flags": string_list(item.get("caution_flags")),
        "upstream_features": {
            key: item.get(key)
            for key in (
                "mngs_signal_tier",
                "rank_priority",
                "reads",
                "support_modules",
                "hospital_evidence_detail",
            )
            if key in item
        },
    }


def exclusion_reasons(item: dict[str, Any]) -> list[str]:
    reasons: list[str] = []
    exclusion = item.get("formal_pick_exclusion_rule")
    if isinstance(exclusion, dict):
        reason = compact(exclusion.get("reason"))
        if reason:
            reasons.append(reason)
    return reasons


def candidate_entry(
    item: dict[str, Any], rank: int, picked_by_key: dict[str, dict[str, Any]]
) -> dict[str, Any] | None:
    name = compact(item.get("organism_name"))
    if not name:
        return None
    key = compact(item.get("canonical_key")) or canonical_key(name)
    picked = picked_by_key.get(key)
    selected_rank = picked.get("rank") if picked else None
    selected_reasons = picked.get("selection_reasons", []) if picked else []
    return {
        "original_name": name,
        "canonical_name": picked.get("canonical_name", name) if picked else name,
        "canonical_key": key,
        "selection_status": "picked" if picked else "not_picked",
        "candidate_rank": rank,
        "selected_rank": selected_rank,
        "classification": item.get("classification"),
        "picked_role": picked.get("picked_role") if picked else None,
        "evidence_source": picked.get("evidence_source") if picked else None,
        "basis_level": picked.get("basis_level") if picked else None,
        "consideration_reasons": string_list(item.get("integrated_reasoning")),
        "decision_reasons": selected_reasons if picked else exclusion_reasons(item),
        "decision_reason_status": (
            "available"
            if (selected_reasons if picked else exclusion_reasons(item))
            else "not_provided_by_upstream"
        ),
        "upstream_features": {
            key_name: item.get(key_name)
            for key_name in (
                "integrated_causative_level",
                "mngs_signal_tier",
                "rank_priority",
                "rank_rule",
                "reads",
                "reads_tier",
                "reads_percentile",
                "dominance_tier",
                "specimen_alignment",
                "specimen_class",
                "source_category",
                "module_support_summary",
                "key_evidence",
                "possibility_level",
                "formal_pick_exclusion_rule",
            )
            if key_name in item
        },
    }


def extract_one(
    path: Path, hospital: str, dataset: str, pipeline_version: str
) -> dict[str, Any]:
    data = load_json(path)
    if not isinstance(data, dict):
        raise ValueError("top-level JSON must be an object")
    patient_id = patient_id_from(path, data)
    deterministic = data.get("deterministic_max")
    if not isinstance(deterministic, dict):
        raise ValueError("deterministic_max is missing")
    summary = deterministic.get("best_available_summary")
    if not isinstance(summary, dict):
        raise ValueError("deterministic_max.best_available_summary is missing")

    picked_items = summary.get("picked_pathogens")
    if not isinstance(picked_items, list):
        raise ValueError("picked_pathogens is not an array")
    selected = [picked_entry(item, rank) for rank, item in enumerate(picked_items, 1)]
    selected_keys = {item["canonical_key"] for item in selected}
    if len(selected_keys) != len(selected):
        raise ValueError("picked_pathogens contains duplicate canonical keys")
    picked_by_key = {item["canonical_key"]: item for item in selected}

    candidate_items = deterministic.get("pathogen_candidates")
    if not isinstance(candidate_items, list):
        candidate_items = []
    candidates = [
        candidate
        for rank, item in enumerate(candidate_items, 1)
        if isinstance(item, dict)
        and (candidate := candidate_entry(item, rank, picked_by_key)) is not None
    ]
    candidate_keys = {item["canonical_key"] for item in candidates}
    for selected_item in selected:
        if selected_item["canonical_key"] not in candidate_keys:
            candidates.append(
                {
                    "original_name": selected_item["original_name"],
                    "canonical_name": selected_item["canonical_name"],
                    "canonical_key": selected_item["canonical_key"],
                    "selection_status": "picked",
                    "candidate_rank": None,
                    "selected_rank": selected_item["rank"],
                    "classification": selected_item.get("classification"),
                    "consideration_reasons": [],
                    "decision_reasons": selected_item["selection_reasons"],
                    "decision_reason_status": (
                        "available"
                        if selected_item["selection_reasons"]
                        else "not_provided_by_upstream"
                    ),
                    "upstream_features": selected_item.get("upstream_features", {}),
                }
            )

    return {
        "schema_version": "selected_pathogens.v1",
        "patient_id": patient_id,
        "hospital": hospital,
        "dataset": dataset,
        "selection_source": "deterministic_max.best_available_summary.picked_pathogens",
        "pipeline_version": pipeline_version,
        "source": {
            "file": str(path.resolve()),
            "sha256": sha256_file(path),
            "report_type": data.get("report_type"),
            "rule_version": deterministic.get("rule_version"),
        },
        "selected_pathogens": selected,
        "candidate_pathogens": candidates,
        "summary": {
            "selection_mode": summary.get("selection_mode"),
            "selected_count": len(selected),
            "upstream_candidate_count": len(candidates),
        },
    }


def main() -> int:
    args = parse_args()
    input_root = args.input_root.resolve()
    if not input_root.is_dir():
        raise SystemExit(f"input root is not a directory: {input_root}")
    output_root = prepare_output_root(args.output_root)
    outputs: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []

    for path in sorted(input_root.glob("*.json")):
        try:
            package = extract_one(path, args.hospital, args.dataset, args.pipeline_version)
            patient_id = package["patient_id"]
            output_file = output_root / f"NGS_patient_{patient_id}_selected_pathogens.json"
            write_json(output_file, package, args.pretty)
            outputs.append(
                {
                    "patient_id": patient_id,
                    "hospital": args.hospital,
                    "selected_count": package["summary"]["selected_count"],
                    "candidate_count": package["summary"]["upstream_candidate_count"],
                    "selected_names": [
                        item["original_name"] for item in package["selected_pathogens"]
                    ],
                    "output_file": str(output_file),
                }
            )
        except Exception as exc:
            failures.append({"source_file": str(path), "error": f"{type(exc).__name__}: {exc}"})

    manifest = {
        "schema_version": "selected_pathogens.manifest.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_root": str(input_root),
        "output_root": str(output_root),
        "hospital": args.hospital,
        "dataset": args.dataset,
        "pipeline_version": args.pipeline_version,
        "summary": {
            "patients_processed": len(outputs),
            "patients_failed": len(failures),
            "selected_pathogens": sum(item["selected_count"] for item in outputs),
            "upstream_candidates": sum(item["candidate_count"] for item in outputs),
        },
        "patients": outputs,
        "failures": failures,
    }
    write_json(output_root / "selected_pathogens_manifest.json", manifest, True)
    with (output_root / "selected_pathogens_summary.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "hospital",
                "patient_id",
                "selected_count",
                "selected_names",
                "candidate_count",
                "output_file",
            ),
        )
        writer.writeheader()
        for item in outputs:
            row = dict(item)
            row["selected_names"] = "; ".join(item["selected_names"])
            writer.writerow(row)

    print(
        f"patients={len(outputs)} failed={len(failures)} "
        f"selected={manifest['summary']['selected_pathogens']} output={output_root}"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

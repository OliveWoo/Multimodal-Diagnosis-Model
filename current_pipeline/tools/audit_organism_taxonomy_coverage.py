"""Audit central taxonomy coverage without using benchmark answers."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

from tools import organism_taxonomy_classifier as taxonomy


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def patient_id(patient_dir: Path) -> str:
    return patient_dir.name.removeprefix("NGS_patient_").removesuffix("_json")


def add_observation(
    observations: dict[str, dict[str, Any]],
    *,
    name: Any,
    biological_class: Any,
    patient: str,
    reads: Any = 0,
) -> None:
    display_name = str(name or "").strip()
    if not display_name:
        return
    profile = taxonomy.classify_organism(display_name, biological_class=biological_class)
    key = str(profile.get("canonical_key") or display_name.lower())
    item = observations.setdefault(
        key,
        {
            "organism_name": display_name,
            "profile": profile,
            "patient_ids": set(),
            "observation_count": 0,
            "max_reads": 0,
            "total_reads": 0,
        },
    )
    item["patient_ids"].add(patient)
    item["observation_count"] += 1
    try:
        numeric_reads = max(0, int(float(reads or 0)))
        item["max_reads"] = max(item["max_reads"], numeric_reads)
        item["total_reads"] += numeric_reads
    except (TypeError, ValueError):
        pass


def audit_roots(roots: list[Path]) -> dict[str, Any]:
    observations: dict[str, dict[str, Any]] = {}
    read_errors: list[dict[str, str]] = []
    patient_count = 0
    ranked_file_count = 0
    for root in roots:
        for patient_dir in sorted(root.glob("NGS_patient_*_json")):
            if not patient_dir.is_dir():
                continue
            patient_count += 1
            pid = patient_id(patient_dir)
            ranked_path = patient_dir / f"NGS_patient_{pid}_mNGS_ranked_candidates.json"
            if ranked_path.is_file():
                ranked_file_count += 1
                try:
                    payload = read_json(ranked_path)
                except Exception as exc:  # pragma: no cover - CLI diagnostics
                    read_errors.append({"path": str(ranked_path), "error": str(exc)})
                    continue
                for record in payload.get("records") or []:
                    if not isinstance(record, dict):
                        continue
                    for candidate in record.get("candidates") or []:
                        if not isinstance(candidate, dict):
                            continue
                        add_observation(
                            observations,
                            name=candidate.get("organism_name") or candidate.get("name"),
                            biological_class=candidate.get("source_category"),
                            patient=pid,
                            reads=candidate.get("reads", candidate.get("sec_hit", 0)),
                        )

    mapping_counts = Counter()
    family_counts = Counter()
    class_counts = Counter()
    conflicts: list[dict[str, Any]] = []
    unmapped: list[dict[str, Any]] = []
    serialized: list[dict[str, Any]] = []
    total_observations = 0
    mapped_observations = 0
    total_reads = 0
    mapped_reads = 0
    for item in observations.values():
        profile = item["profile"]
        is_mapped = profile.get("mapping_status") != "unmapped"
        mapping_counts[str(profile.get("mapping_status"))] += 1
        family_counts[str(profile.get("primary_rule_family"))] += 1
        class_counts[str(profile.get("biological_class"))] += 1
        total_observations += item["observation_count"]
        total_reads += item["total_reads"]
        if is_mapped:
            mapped_observations += item["observation_count"]
            mapped_reads += item["total_reads"]
        row = {
            "organism_name": item["organism_name"],
            "canonical_key": profile.get("canonical_key"),
            "primary_rule_family": profile.get("primary_rule_family"),
            "mapping_status": profile.get("mapping_status"),
            "classification_confidence": profile.get("classification_confidence"),
            "biological_class": profile.get("biological_class"),
            "biological_class_conflict": profile.get("biological_class_conflict"),
            "patient_count": len(item["patient_ids"]),
            "patient_ids": sorted(
                item["patient_ids"],
                key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value),
            ),
            "observation_count": item["observation_count"],
            "max_reads": item["max_reads"],
            "total_reads": item["total_reads"],
        }
        serialized.append(row)
        if profile.get("biological_class_conflict"):
            conflicts.append(row)
        if profile.get("mapping_status") == "unmapped":
            unmapped.append(row)

    serialized.sort(key=lambda row: (-row["patient_count"], -row["max_reads"], row["organism_name"]))
    unmapped.sort(key=lambda row: (-row["max_reads"], -row["patient_count"], row["organism_name"]))
    conflicts.sort(key=lambda row: (-row["max_reads"], row["organism_name"]))
    total = len(serialized)
    mapped = total - len(unmapped)
    return {
        "audit_version": "organism_taxonomy_coverage_v1",
        "answer_blind": True,
        "roots": [str(root) for root in roots],
        "patient_directory_count": patient_count,
        "ranked_file_count": ranked_file_count,
        "unique_organism_count": total,
        "mapped_unique_organism_count": mapped,
        "mapped_unique_organism_fraction": round(mapped / total, 4) if total else 0,
        "candidate_observation_count": total_observations,
        "mapped_candidate_observation_count": mapped_observations,
        "mapped_candidate_observation_fraction": (
            round(mapped_observations / total_observations, 4) if total_observations else 0
        ),
        "candidate_read_count": total_reads,
        "mapped_candidate_read_count": mapped_reads,
        "mapped_candidate_read_fraction": round(mapped_reads / total_reads, 4) if total_reads else 0,
        "mapping_status_counts": dict(mapping_counts),
        "primary_rule_family_counts": dict(family_counts),
        "biological_class_counts": dict(class_counts),
        "biological_class_conflicts": conflicts,
        "unmapped_organisms": unmapped,
        "all_organisms": serialized,
        "read_errors": read_errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit answer-blind central organism taxonomy coverage.")
    parser.add_argument("roots", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit_roots(args.roots)
    text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8-sig")
        print(f"wrote {args.output}")
        print(
            f"mapped={report['mapped_unique_organism_count']}/{report['unique_organism_count']} "
            f"fraction={report['mapped_unique_organism_fraction']}"
        )
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

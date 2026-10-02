"""Create an answer-blind, organism-level queue for taxonomy review."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from tools.organism_taxonomy_classifier import classify_organism


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _ncbi_lookup(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    return {row["input_name"].casefold(): row for row in _read_csv(path)}


def build_queue(
    candidates_csv: Path,
    output_dir: Path,
    *,
    ncbi_resolution_csv: Path | None = None,
    include_inherited: bool = False,
) -> dict[str, Any]:
    """Whitelist organism fields only; never serialize patient-level input rows."""
    rows = _read_csv(candidates_csv)
    if not rows or "organism_name" not in rows[0]:
        raise ValueError("Input must contain an organism_name column")
    ncbi = _ncbi_lookup(ncbi_resolution_csv)
    unique: dict[str, dict[str, Any]] = {}
    for row in rows:
        name = row["organism_name"].strip()
        if not name:
            continue
        supplied_class = (
            row.get("category")
            or row.get("biological_class")
            or row.get("source_category")
            or None
        )
        profile = classify_organism(name, biological_class=supplied_class)
        # Compact taxonomy output already contains the effective reviewed
        # mapping.  Prefer it when present so the LLM queue is exactly the
        # unresolved remainder of that release rather than a second, slightly
        # different classifier pass.
        mapping_status = row.get("taxonomy_mapping_status") or profile["mapping_status"]
        if mapping_status == "exact_species" or (mapping_status != "unmapped" and not include_inherited):
            continue
        key = profile["canonical_key"]
        if not key:
            continue
        category = str(
            supplied_class
            or profile.get("supplied_biological_class")
            or profile.get("biological_class")
            or "unknown"
        )
        existing = unique.get(key)
        if existing is not None:
            if category not in existing["input_classes"]:
                existing["input_classes"].append(category)
            if name not in existing["reported_names"]:
                existing["reported_names"].append(name)
            continue
        unique[key] = {
            "canonical_key": key,
            "organism_name": name,
            "reported_names": [name],
            "input_classes": [category],
            "classification_status": mapping_status,
            "current_rule_source": profile["rule_source"],
            "current_family": row.get("taxonomy_family") or profile["primary_rule_family"],
            "review_status": "pending_identity_and_literature_review",
            "ncbi": {},
        }

    for item in unique.values():
        record = ncbi.get(item["organism_name"].casefold())
        if record:
            item["ncbi"] = {
                "resolution_status": record.get("resolution_status", ""),
                "taxid": record.get("ncbi_taxid", ""),
                "scientific_name": record.get("ncbi_scientific_name", ""),
                "rank": record.get("ncbi_rank", ""),
                "division": record.get("ncbi_division", ""),
                "lineage": record.get("ncbi_lineage", ""),
            }

    queue = sorted(unique.values(), key=lambda item: item["canonical_key"])
    output_dir.mkdir(parents=True, exist_ok=True)
    queue_path = output_dir / "organism_review_queue.jsonl"
    with queue_path.open("w", encoding="utf-8") as handle:
        for item in queue:
            handle.write(json.dumps(item, ensure_ascii=False) + "\n")
    ncbi_input_path = output_dir / "ncbi_lookup_input.csv"
    with ncbi_input_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["organism_name", "category"])
        writer.writeheader()
        for item in queue:
            writer.writerow({"organism_name": item["organism_name"], "category": "|".join(item["input_classes"])})
    summary = {
        "schema_version": "organism_taxonomy_review_queue_v1",
        "answer_blind": True,
        "patient_level_fields_exported": False,
        "input_file": str(candidates_csv.resolve()),
        "input_sha256": _source_hash(candidates_csv),
        "input_columns": list(rows[0]),
        "ncbi_resolution_file": str(ncbi_resolution_csv.resolve()) if ncbi_resolution_csv else None,
        "ncbi_resolution_sha256": _source_hash(ncbi_resolution_csv) if ncbi_resolution_csv else None,
        "unique_review_organisms": len(queue),
        "unique_unmapped_organisms": sum(item["classification_status"] == "unmapped" for item in queue),
        "include_inherited": include_inherited,
        "queue_path": str(queue_path.resolve()),
        "ncbi_lookup_input": str(ncbi_input_path.resolve()),
        "safety_note": "A queue item is not a clinical pathogen or an approved taxonomy rule.",
    }
    (output_dir / "queue_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build an organism-only taxonomy review queue from any candidate CSV.")
    parser.add_argument("candidates_csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ncbi-resolution-csv", type=Path)
    parser.add_argument("--include-inherited", action="store_true", help="Also audit genus/key-prefix inherited species; exact profiles remain excluded")
    args = parser.parse_args()
    print(json.dumps(build_queue(args.candidates_csv, args.output_dir, ncbi_resolution_csv=args.ncbi_resolution_csv, include_inherited=args.include_inherited), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

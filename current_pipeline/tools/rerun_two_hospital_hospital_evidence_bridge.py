"""Rebuild two-hospital microbiology evidence without running the mNGS scorer.

This command deliberately consumes only raw Culture, FilmArray, GM, and
targeted-molecular JSON files.  It preserves source timestamps for later
plus/minus-48-hour alignment, but it does not require or score mNGS data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools import build_deterministic_summary as deterministic_summary
from tools import normalized_agent_fallback


SCHEMA_VERSION = "two_hospital_hospital_evidence_bridge.v1"


def read_json(path: Path | None) -> Any:
    if path is None:
        return None
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def patient_id(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)", path.name)
    if not match:
        raise ValueError(f"Cannot parse patient id from {path.name}")
    return int(match.group(1))


def source_path(patient_dir: Path, suffix: str) -> Path | None:
    path = patient_dir / f"NGS_patient_{patient_id(patient_dir)}_{suffix}.json"
    return path if path.is_file() else None


def source_contract(path: Path | None) -> dict[str, Any]:
    return {
        "status": "available" if path is not None else "missing",
        "path": str(path.resolve()) if path is not None else None,
        "sha256": sha256_file(path) if path is not None else None,
    }


def evidence_times(item: dict[str, Any]) -> tuple[list[str], list[str]]:
    collected: list[str] = []
    reported: list[str] = []
    for rows in (item.get("module_evidence") or {}).values():
        for row in rows if isinstance(rows, list) else []:
            if not isinstance(row, dict):
                continue
            collected_time = str(row.get("collected_time") or "").strip()
            reported_time = str(row.get("reported_time") or "").strip()
            if collected_time and collected_time != "Unknown":
                collected.append(collected_time)
            if reported_time and reported_time != "Unknown":
                reported.append(reported_time)
    return list(dict.fromkeys(collected)), list(dict.fromkeys(reported))


def build_patient_packet(patient_dir: Path) -> dict[str, Any]:
    paths = {
        "culture": source_path(patient_dir, "culture"),
        "filmarray": source_path(patient_dir, "filmarray"),
        "gm": source_path(patient_dir, "gm_test"),
        "targeted_molecular": source_path(patient_dir, "molecular_microbiology"),
    }
    raw = {name: read_json(path) for name, path in paths.items()}
    culture = (
        normalized_agent_fallback.culture_agent_from_source(raw["culture"])
        if raw["culture"] is not None
        else None
    )
    filmarray_gm = (
        normalized_agent_fallback.filmarray_gm_agent_from_sources(
            raw["filmarray"], raw["gm"]
        )
        if raw["filmarray"] is not None or raw["gm"] is not None
        else None
    )
    molecular = (
        normalized_agent_fallback.molecular_microbiology_agent_from_source(
            raw["targeted_molecular"]
        )
        if raw["targeted_molecular"] is not None
        else None
    )
    evidence = deterministic_summary.collect_hospital_organism_evidence(
        culture=culture,
        filmarray=filmarray_gm,
        image=None,
        molecular=molecular,
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "patient_id": patient_id(patient_dir),
        "mngs_scoring_performed": False,
        "event_alignment_status": "pending_mngs_anchor",
        "source_contracts": {
            name: source_contract(path) for name, path in paths.items()
        },
        "hospital_organism_evidence": evidence,
        "constraints": [
            "Only Culture, FilmArray, GM, and targeted molecular sources were normalized.",
            "No mNGS candidate generation, ranking, scorer, promotion, or reporting stage was run.",
            "Collected/reported timestamps are preserved; plus/minus-48-hour alignment waits for the fixed mNGS collection anchor.",
        ],
    }


def run(patient_root: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    packet_root = output_root / "patient_packets"
    packet_root.mkdir(parents=True, exist_ok=True)
    flat_rows: list[dict[str, Any]] = []
    patients = sorted(
        (path for path in patient_root.glob("NGS_patient_*_json") if path.is_dir()),
        key=patient_id,
    )
    for patient_dir in patients:
        packet = build_patient_packet(patient_dir)
        pid = int(packet["patient_id"])
        packet_path = packet_root / f"NGS_patient_{pid}_hospital_evidence_bridge.json"
        packet_path.write_text(
            json.dumps(packet, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        for item in packet["hospital_organism_evidence"]:
            collected, reported = evidence_times(item)
            flat_rows.append({
                "patient_id": pid,
                "organism_name": item.get("organism_name"),
                "classification": item.get("classification"),
                "best_hospital_level": item.get("best_hospital_level"),
                "evidence_modules": ";".join(item.get("evidence_modules") or []),
                "collected_times": ";".join(collected),
                "reported_times": ";".join(reported),
                "event_alignment_status": "pending_mngs_anchor",
            })
    csv_path = output_root / "hospital_organism_evidence.csv"
    fields = [
        "patient_id", "organism_name", "classification", "best_hospital_level",
        "evidence_modules", "collected_times", "reported_times",
        "event_alignment_status",
    ]
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat_rows)
    module_counts: Counter[str] = Counter()
    for row in flat_rows:
        module_counts.update(filter(None, row["evidence_modules"].split(";")))
    summary = {
        "schema_version": SCHEMA_VERSION,
        "patient_count": len(patients),
        "patients_with_hospital_evidence": len(
            {int(row["patient_id"]) for row in flat_rows}
        ),
        "hospital_organism_evidence_count": len(flat_rows),
        "module_evidence_counts": dict(sorted(module_counts.items())),
        "mngs_scoring_performed": False,
        "event_alignment_status": "pending_mngs_anchor",
        "flat_evidence_csv": str(csv_path.resolve()),
    }
    (output_root / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("output_root", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.patient_root, args.output_root), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

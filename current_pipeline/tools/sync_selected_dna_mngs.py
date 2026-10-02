"""Audit and optionally install selected DNA mNGS records for the KH cohort."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PATIENT_ROOT = REPO_ROOT / "outputs" / "patient_info_KH_0728_2Days"
DEFAULT_SOURCE_ROOT = DEFAULT_PATIENT_ROOT / "33_patients_selected_DNA_all_RK_NTC_microbes"
DEFAULT_RUN_ROOT = REPO_ROOT / "outputs" / "runs" / "2026-09-18_KH_0728_selected_DNA_sync"
SOURCE_RE = re.compile(r"^NGS_patient_(\d+)_DNA_[^/\\]+_all_RK_NTC_microbes\.json$")
CATEGORIES = ("bacterial", "viral", "fungal", "others")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_records(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, list) or not payload or any(not isinstance(record, dict) for record in payload):
        raise ValueError(f"Expected nonempty list of mNGS records: {path}")
    for record in payload:
        if not isinstance(record.get("pathogens"), dict):
            raise ValueError(f"Missing pathogens object: {path}")
        for category in CATEGORIES:
            if not isinstance(record["pathogens"].get(category), list):
                raise ValueError(f"Missing {category} list: {path}")
    return payload


def specimen_identity(records: list[dict[str, Any]]) -> list[tuple[Any, Any, Any]]:
    return [
        (record.get("specimen_code"), record.get("collected_time"), record.get("specimen_site"))
        for record in records
    ]


def organisms(records: list[dict[str, Any]]) -> dict[str, list[Any]]:
    result: dict[str, list[Any]] = {}
    for record in records:
        for category in CATEGORIES:
            for item in record["pathogens"][category]:
                if not isinstance(item, dict) or not item.get("name"):
                    raise ValueError(f"Malformed {category} organism: {item!r}")
                result.setdefault(str(item["name"]), []).append(item.get("reads"))
    return result


def build_plan(source_root: Path, patient_root: Path, expected_count: int) -> list[dict[str, Any]]:
    source_files = sorted(source_root.glob("NGS_patient_*_DNA_*_all_RK_NTC_microbes.json"))
    patient_dirs = sorted(patient_root.glob("NGS_patient_*_json"))
    if len(source_files) != expected_count or len(patient_dirs) != expected_count:
        raise ValueError(
            f"Expected {expected_count} source files and patient directories; "
            f"found {len(source_files)} and {len(patient_dirs)}"
        )

    plan: list[dict[str, Any]] = []
    seen: set[int] = set()
    for source in source_files:
        match = SOURCE_RE.fullmatch(source.name)
        if match is None:
            raise ValueError(f"Unexpected source filename: {source}")
        patient_id = int(match.group(1))
        if patient_id in seen:
            raise ValueError(f"Multiple DNA sources for patient {patient_id}")
        seen.add(patient_id)
        target = patient_root / f"NGS_patient_{patient_id}_json" / f"NGS_patient_{patient_id}_all_RK_NTC_microbes.json"
        if not target.is_file():
            raise FileNotFoundError(target)
        source_records = read_records(source)
        current_records = read_records(target)
        if specimen_identity(source_records) != specimen_identity(current_records):
            raise ValueError(f"Specimen code/time/site differs for patient {patient_id}; manual review required")
        source_orgs = organisms(source_records)
        current_orgs = organisms(current_records)
        source_bytes = source.read_bytes()
        current_bytes = target.read_bytes()
        plan.append(
            {
                "patient_id": patient_id,
                "source": str(source.resolve()),
                "target": str(target.resolve()),
                "specimens": specimen_identity(source_records),
                "same_json": source_records == current_records,
                "source_sha256": sha256(source_bytes),
                "current_sha256": sha256(current_bytes),
                "source_organism_count": sum(map(len, source_orgs.values())),
                "current_organism_count": sum(map(len, current_orgs.values())),
                "source_only": sorted(source_orgs.keys() - current_orgs.keys()),
                "current_only": sorted(current_orgs.keys() - source_orgs.keys()),
                "changed_reads": {
                    name: {"source": source_orgs[name], "current": current_orgs[name]}
                    for name in sorted(source_orgs.keys() & current_orgs.keys())
                    if source_orgs[name] != current_orgs[name]
                },
            }
        )
    expected_ids = {int(re.search(r"NGS_patient_(\d+)_json$", path.name).group(1)) for path in patient_dirs}
    if seen != expected_ids:
        raise ValueError(f"Patient sets differ: missing sources={sorted(expected_ids - seen)}; extra={sorted(seen - expected_ids)}")
    return sorted(plan, key=lambda row: row["patient_id"])


def install(plan: list[dict[str, Any]], backup_root: Path) -> None:
    changed = [row for row in plan if not row["same_json"]]
    backup_root.mkdir(parents=True, exist_ok=True)
    for row in changed:
        target = Path(row["target"])
        backup = backup_root / target.name
        if backup.exists() and sha256(backup.read_bytes()) != row["current_sha256"]:
            raise FileExistsError(f"Backup already exists with different contents: {backup}")
    for row in changed:
        target = Path(row["target"])
        backup = backup_root / target.name
        if not backup.exists():
            shutil.copy2(target, backup)
        if sha256(backup.read_bytes()) != row["current_sha256"]:
            raise RuntimeError(f"Backup verification failed: {backup}")
    for row in changed:
        source = Path(row["source"])
        target = Path(row["target"])
        temporary = target.with_name(target.name + ".selected_dna_tmp")
        try:
            shutil.copy2(source, temporary)
            if sha256(temporary.read_bytes()) != row["source_sha256"]:
                raise RuntimeError(f"Copy verification failed: {temporary}")
            os.replace(temporary, target)
        finally:
            temporary.unlink(missing_ok=True)
        if sha256(target.read_bytes()) != row["source_sha256"]:
            raise RuntimeError(f"Installed file hash differs from source: {target}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--expected-count", type=int, default=33)
    parser.add_argument("--apply", action="store_true", help="Back up and replace only semantically different canonical mNGS files")
    args = parser.parse_args()

    plan = build_plan(args.source_root, args.patient_root, args.expected_count)
    changed = [row for row in plan if not row["same_json"]]
    if args.apply:
        install(plan, args.run_root / "original_mngs_backup")
    args.run_root.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": "selected_dna_mngs_sync.v1",
        "applied": args.apply,
        "patient_count": len(plan),
        "different_patient_count": len(changed),
        "backup_root": str((args.run_root / "original_mngs_backup").resolve()) if args.apply else None,
        "patients": plan,
    }
    report_path = args.run_root / ("applied_report.json" if args.apply else "dry_run_report.json")
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"patients={len(plan)} different={len(changed)} applied={args.apply}")
    print(report_path)


if __name__ == "__main__":
    main()

"""Build a 33-patient pipeline input from the enriched workbook shadow.

Only patients with a current per-test DNA/RNA mNGS bundle are forwarded.  The
remaining staging patients are written to an explicit pending list and are not
represented by empty or negative mNGS inputs.  Source files are copied and
renamed to the ``NGS_patient_*`` contract used by the current KH pipeline.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from collections import Counter
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "outputs" / "extracted_workbooks_enriched_shadow_v4_52mngs"
DEFAULT_OUTPUT = ROOT / "outputs" / "extracted_workbooks_pipeline_input_v3_51patients"
PATIENT_RE = re.compile(r"^extracted_patient_(\d+)_json$")

CLINICAL_SUFFIXES = (
    "CBC",
    "other_lab",
    "culture",
    "culture_unlinked_ast",
    "filmarray",
    "gm_test",
    "molecular_microbiology",
    "image",
    "underlying",
    "admission_diagnosis",
)
CORE_CLINICAL_SUFFIXES = tuple(
    suffix for suffix in CLINICAL_SUFFIXES if suffix != "culture_unlinked_ast"
)
MNGS_SUFFIXES = (
    "all_mNGS_DNA_RNA_RK_NTC_by_species",
    "mNGS_ranked_candidates",
)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_empty_payload(path: Path) -> bool:
    try:
        payload = read_json(path)
    except (OSError, json.JSONDecodeError):
        return True
    return payload in (None, [], {})


def copy_traced(source: Path, destination: Path) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    source_hash = sha256_file(source)
    destination_hash = sha256_file(destination)
    if source_hash != destination_hash:
        raise RuntimeError(f"Copy hash mismatch: {source} -> {destination}")
    return {
        "source": str(source.resolve()),
        "destination": str(destination.resolve()),
        "sha256": source_hash,
        "bytes": source.stat().st_size,
    }


def mngs_counts(path: Path) -> dict[str, int]:
    payload = read_json(path)
    cases = payload.get("cases") or [] if isinstance(payload, dict) else []
    tests = sum(len(case.get("tests") or []) for case in cases if isinstance(case, dict))
    return {"case_count": len(cases), "test_count": tests}


def build_patient(source_dir: Path, patients_root: Path, flat_mngs_root: Path) -> dict[str, Any]:
    match = PATIENT_RE.fullmatch(source_dir.name)
    if match is None:
        raise ValueError(f"Unexpected patient directory: {source_dir}")
    patient_id = int(match.group(1))
    source_prefix = f"extracted_patient_{patient_id}"
    target_prefix = f"NGS_patient_{patient_id}"
    multi_assay = source_dir / f"{source_prefix}_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
    ranked = source_dir / f"{source_prefix}_mNGS_ranked_candidates.json"
    enrichment_manifest_path = source_dir / "enrichment_manifest.json"
    enrichment_manifest = (
        read_json(enrichment_manifest_path) if enrichment_manifest_path.is_file() else {}
    )

    clinical_status: dict[str, str] = {}
    for suffix in CLINICAL_SUFFIXES:
        path = source_dir / f"{source_prefix}_{suffix}.json"
        if not path.is_file():
            clinical_status[suffix] = "file_missing"
        elif is_empty_payload(path):
            clinical_status[suffix] = "present_empty_unknown_not_negative"
        else:
            clinical_status[suffix] = "present_nonempty"

    if not multi_assay.is_file():
        return {
            "patient_id": patient_id,
            "status": "mngs_pending_not_forwarded",
            "source_patient_dir": str(source_dir.resolve()),
            "clinical_status": clinical_status,
        }
    destination = patients_root / f"{target_prefix}_json"
    destination.mkdir(parents=True, exist_ok=False)
    copied: list[dict[str, Any]] = []
    mapped_files: dict[str, list[str]] = {}
    for suffix in CLINICAL_SUFFIXES + MNGS_SUFFIXES:
        source_path = source_dir / f"{source_prefix}_{suffix}.json"
        if not source_path.is_file():
            continue
        target_name = f"{target_prefix}_{suffix}.json"
        copied.append(copy_traced(source_path, destination / target_name))
        mapped_files[suffix] = [target_name]

    flat_name = f"{target_prefix}_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
    flat_copy = copy_traced(multi_assay, flat_mngs_root / flat_name)
    copied.append(flat_copy)
    (destination / "agent_outputs").mkdir()
    (destination / "summary_outputs").mkdir()

    counts = mngs_counts(multi_assay)
    manifest = {
        "schema_version": "extracted_workbooks_pipeline_input.v3",
        "patient_id": patient_id,
        "status": "ready_for_small_agents_and_multi_assay_inventory",
        "source_patient_dir": str(source_dir.resolve()),
        "source_enrichment_manifest": str(enrichment_manifest_path.resolve()),
        "same_number_linkage_audit": enrichment_manifest.get("same_number_linkage_audit"),
        "clinical_status": clinical_status,
        "mngs": {
            "status": "current_per_test_bundle_present",
            **counts,
            "ranked_candidates_status": (
                "legacy_ranked_present_for_context_only"
                if ranked.is_file()
                else "pending_deterministic_derivation_from_multi_assay_compact_packet"
            ),
            "empty_legacy_mngs_grouped_intentionally_omitted": True,
        },
        "copied_files": copied,
        "not_copied": [
            "legacy_raw_sources",
            "old agent_outputs or summary_outputs",
            "empty legacy mNGS_grouped.json",
        ],
    }
    write_json(destination / f"{target_prefix}_pipeline_input_manifest.json", manifest)
    write_json(destination / "source_mapping.json", mapped_files)
    return manifest


def run(source_root: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output root: {output_root}")
    patients_root = output_root / "patients"
    flat_mngs_root = output_root / "multi_assay_source"
    patients_root.mkdir(parents=True, exist_ok=True)
    flat_mngs_root.mkdir(parents=True, exist_ok=True)

    source_dirs = sorted(
        (
            path for path in source_root.iterdir()
            if path.is_dir() and PATIENT_RE.fullmatch(path.name)
        ),
        key=lambda path: int(PATIENT_RE.fullmatch(path.name).group(1)),
    )
    records = [build_patient(path, patients_root, flat_mngs_root) for path in source_dirs]
    ready = [record for record in records if record["status"].startswith("ready_")]
    pending = [record for record in records if record["status"] == "mngs_pending_not_forwarded"]
    missing_modules: Counter[str] = Counter()
    for record in ready:
        missing_modules.update(
            suffix for suffix, status in record["clinical_status"].items()
            if suffix in CORE_CLINICAL_SUFFIXES and status != "present_nonempty"
        )

    pending_payload = {
        "schema_version": "mngs_pending_patients.v1",
        "interpretation": "Clinical data may be present, but these patients are excluded until a current per-test mNGS bundle is supplied.",
        "patient_count": len(pending),
        "patients": pending,
    }
    write_json(output_root / "mngs_pending_patients.json", pending_payload)
    summary = {
        "schema_version": "extracted_workbooks_pipeline_input_summary.v3",
        "source_root": str(source_root.resolve()),
        "output_root": str(output_root.resolve()),
        "source_patient_count": len(records),
        "forwarded_patient_count": len(ready),
        "forwarded_patient_ids": [record["patient_id"] for record in ready],
        "pending_mngs_patient_count": len(pending),
        "pending_mngs_patient_ids": [record["patient_id"] for record in pending],
        "multi_assay_case_count": sum(record["mngs"]["case_count"] for record in ready),
        "multi_assay_test_count": sum(record["mngs"]["test_count"] for record in ready),
        "clinical_missing_or_empty_counts_in_forwarded_cohort": dict(sorted(missing_modules.items())),
        "patient_root_for_small_agents": str(patients_root.resolve()),
        "source_root_for_inventory_builder": str(flat_mngs_root.resolve()),
        "ready_for_small_agents": len(ready) == len(records),
        "ready_for_multi_assay_inventory": len(ready) == len(records),
        "ready_for_integrated_release": False,
        "next_required_steps": [
            "build multi-assay inventory from multi_assay_source",
            "run the five clinical small agents on patients",
            "build the current deterministic clinical summary contract",
            "run integrated decision stages in a new shadow output",
        ],
    }
    write_json(output_root / "pipeline_input_summary.json", summary)
    (output_root / "README.md").write_text(
        "# Extracted workbooks pipeline input v3 (51 patients with clinical staging)\n\n"
        "This is a shadow input. It does not modify the source staging or prior KH outputs.\n\n"
        f"- Forwarded patients: {len(ready)}\n"
        f"- mNGS pending: {len(pending)}\n"
        "- Zero-read rows remain source QC records and are not positive mNGS evidence.\n"
        "- The old empty `mNGS_grouped.json` files are intentionally omitted.\n"
        "- Empty clinical modules mean unavailable/unknown, not negative.\n",
        encoding="utf-8",
    )
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    print(json.dumps(run(args.source_root, args.output_root), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

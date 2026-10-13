"""Build an answer-blind external input root containing no clinical evidence.

The generated patient directories intentionally contain only an empty hospital
summary.  The multi-assay inventory remains the sole event and organism source
when this root is passed to ``run_external_integrated_release_v2``.  Missing
clinical evidence is represented as unavailable/unknown, never as negative.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


SUMMARY_SUFFIX = "mngs_only_no_clinical_evidence_v1"


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def patient_id_from_inventory(path: Path) -> int:
    prefix = "NGS_patient_"
    suffix = "_multi_assay_mngs_evidence.json"
    if not path.name.startswith(prefix) or not path.name.endswith(suffix):
        raise ValueError(f"Unexpected inventory filename: {path.name}")
    return int(path.name[len(prefix):-len(suffix)])


def empty_summary(patient_id: int) -> dict[str, Any]:
    return {
        "schema_version": "mngs_only_empty_clinical_summary.v1",
        "patient_id": patient_id,
        "answer_blind": True,
        "decision_input_mode": "mngs_only",
        "host_context": {
            "host_vulnerability_tier": "Unknown",
            "opportunistic_coverage_level": "Unknown",
            "expanded_candidate_policy": False,
            "opportunistic_host_support": False,
            "opportunistic_host_review_context": False,
            "immunocompromise_tier": "Unknown",
            "acute_instability_tier": "Unknown",
        },
        "hospital_organism_evidence": [],
        "pathogen_candidates": [],
        "excluded_candidates": [],
        "priority_flags": [],
        "module_summaries": {
            "culture": {"availability": "not_provided_to_mngs_only_run"},
            "filmarray_gmtest": {"availability": "not_provided_to_mngs_only_run"},
            "molecular_microbiology": {"availability": "not_provided_to_mngs_only_run"},
            "image": {"availability": "not_provided_to_mngs_only_run"},
            "cbc_other_lab": {
                "availability": "not_provided_to_mngs_only_run",
                "host_vulnerability_tier": "Unknown",
                "opportunistic_coverage_level": "Unknown",
                "expanded_candidate_policy": False,
                "opportunistic_host_support": False,
            },
        },
        "source_files": {},
        "data_gaps": [
            "culture_not_provided",
            "filmarray_gm_not_provided",
            "targeted_molecular_not_provided",
            "imaging_not_provided",
            "cbc_other_lab_not_provided",
            "history_and_medication_not_provided",
        ],
        "constraints": [
            "This summary is an intentional mNGS-only isolation control.",
            "Unavailable clinical evidence is unknown and is not negative evidence.",
            "No benchmark answer was loaded while creating this input root.",
        ],
    }


def run(inventory_root: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty directory: {output_root}")
    inventory_paths = sorted(
        (inventory_root / "patients").glob(
            "NGS_patient_*_multi_assay_mngs_evidence.json"
        ),
        key=patient_id_from_inventory,
    )
    if not inventory_paths:
        raise FileNotFoundError(f"No patient inventory files under {inventory_root}")

    # The integrated runner expects a phenotype root even when phenotype data
    # is intentionally absent.  An empty directory preserves that absence as
    # explicit unknown evidence rather than blocking an mNGS-only run.
    phenotype_root = output_root / "phenotype_unavailable"
    phenotype_root.mkdir(parents=True, exist_ok=True)

    patient_ids: list[int] = []
    for inventory_path in inventory_paths:
        patient_id = patient_id_from_inventory(inventory_path)
        inventory = read_json(inventory_path)
        if int(inventory.get("patient_id")) != patient_id:
            raise ValueError(f"Inventory patient mismatch: {inventory_path}")
        patient_ids.append(patient_id)
        summary_path = (
            output_root
            / f"NGS_patient_{patient_id}_json"
            / "summary_outputs"
            / f"NGS_patient_{patient_id}_{SUMMARY_SUFFIX}.json"
        )
        write_json(summary_path, empty_summary(patient_id))

    manifest = {
        "schema_version": "mngs_only_external_input.v1",
        "answer_blind": True,
        "patient_count": len(patient_ids),
        "patient_ids": patient_ids,
        "inventory_root": str(inventory_root.resolve()),
        "output_root": str(output_root.resolve()),
        "summary_suffix": SUMMARY_SUFFIX,
        "phenotype_root": str(phenotype_root.resolve()),
        "clinical_files_copied": 0,
        "hospital_evidence_items": 0,
        "constraints": [
            "Only patient identity and empty clinical summaries are generated here.",
            "mNGS cases, tests, reads, QC, rank, specimen, and taxonomy are read from the separate frozen inventory.",
            "No culture, FilmArray, GM, targeted PCR, imaging, CBC, history, medication, phenotype, or benchmark answer is included.",
        ],
    }
    write_json(output_root / "mngs_only_input_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory_root", type=Path)
    parser.add_argument("output_root", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.inventory_root, args.output_root), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

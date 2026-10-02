"""Audit whether the two-hospital cohort can validate the current KH pipeline.

The audit is intentionally conservative: it distinguishes candidate identity
transportability from full per-test scorer compatibility. Historical outcome
metrics are not recomputed when required molecule/QC/normalization fields are
missing.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from tools.organism_taxonomy_classifier import classify_organism


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PATIENT_ROOT = ROOT / "outputs" / "patient_info_two_hospitals_0611_rich_moderate"
DEFAULT_MANIFEST = (
    ROOT
    / "outputs"
    / "runs"
    / "2026-09-21_pipeline_method_alignment"
    / "two_hospitals_answered_rich_moderate_20_20260826.manifest.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "outputs"
    / "runs"
    / "2026-09-30_two_hospital_current_contract_compatibility_v1"
)
PATIENT_RE = re.compile(r"NGS_patient_(\d+)", re.IGNORECASE)


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _as_float(value: Any) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _infer_molecule(source_file: Any) -> str:
    text = str(source_file or "").upper()
    has_dna = "DNA" in text
    has_rna = "RNA" in text
    if has_dna and not has_rna:
        return "DNA"
    if has_rna and not has_dna:
        return "RNA"
    return "UNKNOWN"


def _biological_class(source_category: Any) -> str:
    text = str(source_category or "").lower()
    if "bac" in text:
        return "bacterium"
    if "fung" in text:
        return "fungus"
    if "virus" in text:
        return "virus"
    if "paras" in text:
        return "parasite"
    return "unknown"


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: Iterable[str]) -> None:
    fieldnames = list(fields)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _ranked_metadata(patient_dir: Path, patient_id: str) -> tuple[set[str], set[str]]:
    path = patient_dir / f"NGS_patient_{patient_id}_mNGS_ranked_candidates.json"
    if not path.exists():
        return set(), set()
    payload = _read_json(path)
    records = payload.get("records", []) if isinstance(payload, dict) else []
    times = {str(row.get("collected_time") or "").strip() for row in records}
    sites = {str(row.get("specimen_site") or "").strip() for row in records}
    return {value for value in times if value}, {value for value in sites if value}


def audit(patient_root: Path, manifest_path: Path) -> tuple[
    dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]
]:
    manifest = _read_json(manifest_path)
    answered = {str(value) for value in manifest.get("patient_selection", [])}
    patient_dirs: dict[str, Path] = {}
    for path in patient_root.glob("NGS_patient_*_json"):
        match = PATIENT_RE.search(path.name)
        if match and path.is_dir():
            patient_dirs[match.group(1)] = path

    patient_rows: list[dict[str, Any]] = []
    unique_candidates: dict[tuple[str, str], dict[str, Any]] = {}
    for patient_id in sorted(patient_dirs, key=int):
        patient_dir = patient_dirs[patient_id]
        raw_path = patient_dir / f"NGS_patient_{patient_id}_raw_RK_NTC_rows.json"
        raw_rows = _read_json(raw_path) if raw_path.exists() else []
        raw_rows = raw_rows if isinstance(raw_rows, list) else []
        seq_ids = {str(row.get("seq_id") or "").strip() for row in raw_rows}
        seq_ids.discard("")
        specimen_ids = {str(row.get("specimen_id") or "").strip() for row in raw_rows}
        specimen_ids.discard("")
        molecule_by_seq: dict[str, set[str]] = defaultdict(set)
        for row in raw_rows:
            seq_id = str(row.get("seq_id") or "").strip()
            molecule_by_seq[seq_id].add(_infer_molecule(row.get("source_file")))
            organism = str(row.get("organism_name") or "").strip()
            category = str(row.get("source_category") or "").strip()
            if organism:
                key = (organism.lower(), category.lower())
                source_row = row.get("row", {}) if isinstance(row.get("row"), dict) else {}
                profile = classify_organism(
                    organism, biological_class=_biological_class(category)
                )
                if key not in unique_candidates:
                    unique_candidates[key] = {
                        "organism_name": organism,
                        "source_category": category,
                        "source_biological_class": _biological_class(category),
                        "source_taxid": str(source_row.get("taxid") or "").strip(),
                        "canonical_key": profile["canonical_key"],
                        "central_taxid": profile.get("taxid") or "",
                        "mapping_status": profile["mapping_status"],
                        "primary_rule_family": profile["primary_rule_family"],
                        "biological_class_conflict": profile["biological_class_conflict"],
                        "patient_count": 0,
                        "row_count": 0,
                        "_patients": set(),
                    }
                target = unique_candidates[key]
                target["row_count"] += 1
                target["_patients"].add(patient_id)

        molecule_known_tests = sum(
            bool(values) and values <= {"DNA", "RNA"} for values in molecule_by_seq.values()
        )
        molecule_unknown_tests = len(molecule_by_seq) - molecule_known_tests
        source_taxid_rows = sum(
            bool(str((row.get("row") or {}).get("taxid") or "").strip())
            for row in raw_rows
            if isinstance(row.get("row"), dict)
        )
        reads_rows = sum(
            _as_float((row.get("row") or {}).get("Sec.hit")) is not None
            for row in raw_rows
            if isinstance(row.get("row"), dict)
        )
        candidate_qc_rows = sum(
            any(
                str((row.get("row") or {}).get(field) or "").strip()
                for field in ("Code_NTC", "Code_RK_NTC", "NTC_ratio", "NTC_ratio.adj")
            )
            for row in raw_rows
            if isinstance(row.get("row"), dict)
        )
        collected_times, specimen_sites = _ranked_metadata(patient_dir, patient_id)
        patient_rows.append(
            {
                "patient_id": patient_id,
                "in_answered_20_manifest": patient_id in answered,
                "raw_file_exists": raw_path.exists(),
                "candidate_row_count": len(raw_rows),
                "distinct_test_count": len(seq_ids),
                "distinct_specimen_count": len(specimen_ids),
                "molecule_known_test_count": molecule_known_tests,
                "molecule_unknown_test_count": molecule_unknown_tests,
                "has_both_dna_and_rna": {
                    molecule
                    for values in molecule_by_seq.values()
                    for molecule in values
                }
                >= {"DNA", "RNA"},
                "candidate_reads_row_count": reads_rows,
                "source_taxid_row_count": source_taxid_rows,
                "candidate_qc_row_count": candidate_qc_rows,
                "has_structured_input_reads": False,
                "has_structured_nonhost_denominator": False,
                "has_structured_protocol_condition": False,
                "collected_time_available": bool(collected_times),
                "specimen_site_available": bool(specimen_sites),
                "current_full_contract_compatible": bool(raw_rows)
                and molecule_unknown_tests == 0
                and False,
                "compatibility_reason": (
                    "missing_structured_input_reads_nonhost_denominator_protocol"
                    + ("_and_molecule_type" if molecule_unknown_tests else "")
                ),
            }
        )

    taxonomy_rows: list[dict[str, Any]] = []
    for row in unique_candidates.values():
        row["patient_count"] = len(row.pop("_patients"))
        taxonomy_rows.append(row)
    taxonomy_rows.sort(key=lambda row: (row["mapping_status"], row["organism_name"].lower()))

    gates = [
        {
            "gate": "patient_identifier",
            "required_for": "full_pipeline",
            "available": all(row["patient_id"] for row in patient_rows),
            "decision": "pass",
            "detail": "Directory and manifest patient identifiers are available.",
        },
        {
            "gate": "per_test_seq_id",
            "required_for": "test_aware_scorer",
            "available": all(row["distinct_test_count"] > 0 for row in patient_rows if row["candidate_row_count"]),
            "decision": "pass",
            "detail": "Raw candidate rows retain seq_id; empty-candidate P86 remains explicit.",
        },
        {
            "gate": "structured_dna_rna_identity",
            "required_for": "cross_molecule_support",
            "available": all(row["molecule_unknown_test_count"] == 0 for row in patient_rows),
            "decision": "fail",
            "detail": "Molecule can only be inferred from some source paths and is unknown for many tests.",
        },
        {
            "gate": "candidate_reads",
            "required_for": "analytical_strength",
            "available": all(
                row["candidate_reads_row_count"] == row["candidate_row_count"]
                for row in patient_rows
            ),
            "decision": "pass",
            "detail": "Sec.hit is retained for all candidate rows.",
        },
        {
            "gate": "candidate_taxid",
            "required_for": "identity_audit",
            "available": all(
                row["source_taxid_row_count"] == row["candidate_row_count"]
                for row in patient_rows
            ),
            "decision": "pass",
            "detail": "Source taxid is retained in raw candidate rows.",
        },
        {
            "gate": "candidate_level_ntc_codes",
            "required_for": "legacy_selection_replay",
            "available": all(
                row["candidate_qc_row_count"] == row["candidate_row_count"]
                for row in patient_rows
            ),
            "decision": "pass",
            "detail": "Candidate-level NTC/RK fields exist, but do not replace test-level QC denominators.",
        },
        {
            "gate": "input_reads_and_nonhost_denominator",
            "required_for": "rpm_normalization",
            "available": False,
            "decision": "fail",
            "detail": "No structured total-input/non-host denominator is preserved.",
        },
        {
            "gate": "structured_protocol_condition",
            "required_for": "technical_repeat_interpretation",
            "available": False,
            "decision": "fail",
            "detail": "Adapter/centrifugation/rerun condition is not represented as a structured field.",
        },
        {
            "gate": "collection_time_and_specimen_site",
            "required_for": "clinical_timeline",
            "available": all(
                row["collected_time_available"] and row["specimen_site_available"]
                for row in patient_rows
            ),
            "decision": "pass",
            "detail": (
                "Flattened ranked artifacts retain date/site for most patients, but P76, "
                "P86, and P102 are incomplete; P76 is part of the answered-20 manifest."
            ),
        },
        {
            "gate": "current_contract_ground_truth_freeze",
            "required_for": "held_out_metrics",
            "available": False,
            "decision": "fail",
            "detail": "Historical answered manifest exists, but no untouched current-contract validation freeze exists.",
        },
    ]
    for gate in gates:
        gate["decision"] = "pass" if gate["available"] else "fail"

    counts = {
        "patient_directories": len(patient_rows),
        "answered_manifest_patients": len(answered),
        "answered_manifest_distinct_tests": sum(
            row["distinct_test_count"]
            for row in patient_rows
            if row["in_answered_20_manifest"]
        ),
        "answered_manifest_molecule_known_tests": sum(
            row["molecule_known_test_count"]
            for row in patient_rows
            if row["in_answered_20_manifest"]
        ),
        "answered_manifest_molecule_unknown_tests": sum(
            row["molecule_unknown_test_count"]
            for row in patient_rows
            if row["in_answered_20_manifest"]
        ),
        "raw_candidate_rows": sum(row["candidate_row_count"] for row in patient_rows),
        "distinct_tests": sum(row["distinct_test_count"] for row in patient_rows),
        "molecule_known_tests": sum(row["molecule_known_test_count"] for row in patient_rows),
        "molecule_unknown_tests": sum(row["molecule_unknown_test_count"] for row in patient_rows),
        "patients_with_both_inferred_dna_rna": sum(row["has_both_dna_and_rna"] for row in patient_rows),
        "unique_candidate_names": len(taxonomy_rows),
        "taxonomy_mapping_status": dict(Counter(row["mapping_status"] for row in taxonomy_rows)),
        "taxonomy_family": dict(Counter(row["primary_rule_family"] for row in taxonomy_rows)),
        "taxonomy_biological_class_conflicts": sum(
            row["biological_class_conflict"] for row in taxonomy_rows
        ),
        "source_taxid_present_unique": sum(bool(row["source_taxid"]) for row in taxonomy_rows),
        "central_taxid_present_unique": sum(bool(row["central_taxid"]) for row in taxonomy_rows),
        "full_contract_compatible_patients": sum(
            row["current_full_contract_compatible"] for row in patient_rows
        ),
    }
    summary = {
        "schema_version": "external_validation_compatibility_audit.v1",
        "answer_blind_for_pipeline_rules": True,
        "counts": counts,
        "gates": gates,
        "full_pipeline_external_validation_compatible": all(
            gate["available"] for gate in gates
        ),
        "partial_validation_ready": {
            "candidate_identity_and_taxonomy_transportability": True,
            "legacy_candidate_selection_replay": True,
            "current_test_aware_scorer": False,
            "current_cross_molecule_logic": False,
            "current_rpm_guardrail": False,
            "current_held_out_performance_metrics": False,
        },
        "historical_metric_use": (
            "Historical two-hospital metrics remain descriptive for the old artifact chain and "
            "must not be compared as if they were produced by the current KH per-test contract."
        ),
        "minimum_data_request": [
            "Per seq_id DNA/RNA molecule type (not inferred from file paths).",
            "Per test total input reads and preferably non-host reads for RPM denominators.",
            "Structured protocol/technical-repeat fields (rerun, re-extraction, adapter, fraction).",
            "A frozen untouched patient/answer subset evaluated only after the current policy is locked.",
        ],
    }
    return summary, patient_rows, taxonomy_rows, gates


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary, patient_rows, taxonomy_rows, gates = audit(args.patient_root, args.manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "compatibility_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    _write_csv(
        args.output_dir / "patient_contract_compatibility.csv",
        patient_rows,
        [
            "patient_id",
            "in_answered_20_manifest",
            "raw_file_exists",
            "candidate_row_count",
            "distinct_test_count",
            "distinct_specimen_count",
            "molecule_known_test_count",
            "molecule_unknown_test_count",
            "has_both_dna_and_rna",
            "candidate_reads_row_count",
            "source_taxid_row_count",
            "candidate_qc_row_count",
            "has_structured_input_reads",
            "has_structured_nonhost_denominator",
            "has_structured_protocol_condition",
            "collected_time_available",
            "specimen_site_available",
            "current_full_contract_compatible",
            "compatibility_reason",
        ],
    )
    _write_csv(
        args.output_dir / "candidate_taxonomy_transportability.csv",
        taxonomy_rows,
        [
            "organism_name",
            "source_category",
            "source_biological_class",
            "source_taxid",
            "canonical_key",
            "central_taxid",
            "mapping_status",
            "primary_rule_family",
            "biological_class_conflict",
            "patient_count",
            "row_count",
        ],
    )
    _write_csv(
        args.output_dir / "validation_gate_matrix.csv",
        gates,
        ["gate", "required_for", "available", "decision", "detail"],
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

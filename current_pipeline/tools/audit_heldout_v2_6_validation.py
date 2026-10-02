"""Audit the answer-blind 18-patient held-out v2.6 shadow.

This audit never reads benchmark answers.  It verifies the per-test source
contract, the frozen development/held-out partition, and candidate continuity
from compact entry through the integrated report.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RELEASE = (
    ROOT
    / "outputs/runs/2026-09-30_51patients_integrated_release_v2_6_heldout_freeze_v1"
)
DEFAULT_INVENTORY = (
    ROOT / "outputs/runs/2026-09-30_51patients_multi_assay_inventory_no_zero_v1"
)
DEFAULT_OUTPUT = (
    ROOT / "outputs/runs/2026-09-30_18patients_heldout_v2_6_validation_audit_v1"
)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def candidate_key(row: dict[str, str]) -> tuple[str, str]:
    return (
        str(row.get("patient_id") or "").strip(),
        str(row.get("organism_name") or "").strip().casefold(),
    )


def audit(release_root: Path, inventory_root: Path) -> tuple[
    dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]
]:
    freeze = read_json(release_root / "answer_blind_freeze.json")
    partition = freeze["cohort_partition"]
    heldout_ids = {str(value) for value in partition["heldout_patient_ids"]}
    development_ids = {
        str(value) for value in partition["development_overlap_patient_ids"]
    }
    test_rows = [
        row
        for row in read_csv(inventory_root / "test_qc.csv")
        if row["patient_id"] in heldout_ids
    ]
    observation_rows = [
        row
        for row in read_csv(inventory_root / "positive_observations.csv")
        if row["patient_id"] in heldout_ids
    ]
    union_rows = [
        row
        for row in read_csv(inventory_root / "candidate_union.csv")
        if row["patient_id"] in heldout_ids
    ]
    compact_rows = [
        row
        for row in read_csv(release_root / "01_compact_taxonomy/all_decisions.csv")
        if row["patient_id"] in heldout_ids
    ]
    scorer_rows = [
        row
        for row in read_csv(
            release_root / "02_analytical_scorer/patient_organism_decisions.csv"
        )
        if row["patient_id"] in heldout_ids
    ]
    integrated_rows = [
        row
        for row in read_csv(release_root / "integrated_decisions.csv")
        if row["patient_id"] in heldout_ids
    ]
    report_rows = [
        row
        for row in read_csv(
            release_root / "10_reporting_v3_route_a/complete_report.csv"
        )
        if row["patient_id"] in heldout_ids
    ]

    test_keys = [
        (row["patient_id"], row["case_review_id"], row["seq_id"])
        for row in test_rows
    ]
    required_test_fields = (
        "patient_id",
        "specimen_code",
        "collected_time",
        "specimen_site",
        "seq_id",
        "nucleic_type",
        "qc_status",
    )
    missing_required = Counter()
    for row in test_rows:
        for field in required_test_fields:
            if not str(row.get(field) or "").strip():
                missing_required[field] += 1
    missing_denominator_rows = [
        row for row in test_rows if not str(row.get("input_reads") or "").strip()
    ]
    missing_denominator_without_marker = [
        row
        for row in missing_denominator_rows
        if not str(row.get("normalization_missingness") or "").strip()
    ]
    invalid_molecule = [
        row for row in test_rows if row.get("nucleic_type") not in {"DNA", "RNA"}
    ]
    nonpositive_observations = [
        row
        for row in observation_rows
        if float(row.get("reads") or 0) <= 0
    ]
    invalid_selected = [
        row
        for row in observation_rows
        if str(row.get("selected") or "").strip().lower()
        not in {"true", "false"}
    ]
    protocol_annotated = [
        row
        for row in test_rows
        if str(row.get("condition") or "").strip() not in {"", "-"}
    ]

    compact_scorer_keys = {
        candidate_key(row)
        for row in compact_rows
        if row.get("disposition") == "scorer_entry"
    }
    scorer_forward_keys = {
        candidate_key(row)
        for row in scorer_rows
        if as_bool(row.get("forward_to_clinical_scorer"))
    }
    scorer_keys = {candidate_key(row) for row in scorer_rows}
    integrated_keys = {candidate_key(row) for row in integrated_rows}
    missing_compact_to_scorer = sorted(compact_scorer_keys - scorer_keys)
    scorer_not_forwarded = sorted(scorer_keys - scorer_forward_keys)
    missing_scorer_to_integrated = sorted(scorer_forward_keys - integrated_keys)
    integrated_hospital_only = sorted(integrated_keys - scorer_forward_keys)
    disappearance_rows = [
        {"patient_id": patient, "organism_name_key": organism, "stage": stage}
        for stage, keys in (
            ("compact_scorer_entry_missing_from_patient_scorer", missing_compact_to_scorer),
            ("scorer_forward_missing_from_integrated", missing_scorer_to_integrated),
            ("integrated_hospital_only_or_noncompact", integrated_hospital_only),
        )
        for patient, organism in keys
    ]

    patient_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in test_rows:
        patient_counts[row["patient_id"]]["tests"] += 1
        if not str(row.get("input_reads") or "").strip():
            patient_counts[row["patient_id"]]["tests_without_input_reads"] += 1
    for row in compact_rows:
        patient_counts[row["patient_id"]][f"compact_{row['disposition']}"] += 1
    for row in integrated_rows:
        patient_counts[row["patient_id"]][
            f"tier_{row['final_reporting_tier']}"
        ] += 1
    for row in report_rows:
        patient_counts[row["patient_id"]]["complete_report"] += 1
    patient_rows = []
    for patient_id in sorted(heldout_ids, key=int):
        counts = patient_counts[patient_id]
        patient_rows.append({
            "patient_id": patient_id,
            "tests": counts["tests"],
            "tests_without_input_reads": counts["tests_without_input_reads"],
            "compact_scorer_entry": counts["compact_scorer_entry"],
            "compact_clinical_review": counts["compact_clinical_review"],
            "compact_qc_only": counts["compact_qc_only"],
            "Picked": counts["tier_Picked"],
            "Possible": counts["tier_Possible"],
            "Fallback-Possible": counts["tier_Fallback-Possible"],
            "Context": counts["tier_Context"],
            "complete_report": counts["complete_report"],
        })

    gates = [
        {
            "gate": "frozen_development_heldout_partition",
            "passed": not bool(development_ids & heldout_ids)
            and len(heldout_ids) == 18,
            "detail": "Held-out IDs are excluded from the frozen 33-patient development set.",
        },
        {
            "gate": "answer_blind_release_freeze",
            "passed": bool(freeze.get("passed"))
            and not bool(freeze.get("answers_loaded")),
            "detail": "Decisions and held-out outputs were frozen before any answer file was loaded.",
        },
        {
            "gate": "per_test_identity_and_timing",
            "passed": len(test_keys) == len(set(test_keys))
            and not missing_required,
            "detail": "Patient/case/seq_id keys, specimen, collection time, molecule, and QC are required.",
        },
        {
            "gate": "dna_rna_molecule_contract",
            "passed": not invalid_molecule,
            "detail": "Every held-out test has an explicit DNA or RNA molecule type.",
        },
        {
            "gate": "positive_read_contract",
            "passed": not nonpositive_observations,
            "detail": "Positive observations require reads > 0.",
        },
        {
            "gate": "lab_selection_and_qc_status",
            "passed": not invalid_selected
            and all(str(row.get("test_qc_status") or "").strip() for row in observation_rows),
            "detail": "Lab selection and per-test QC are explicit for every positive observation.",
        },
        {
            "gate": "missing_denominator_is_explicit_unknown",
            "passed": not missing_denominator_without_marker,
            "detail": "Missing input reads disable RPM for that test and are never interpreted as zero.",
        },
        {
            "gate": "candidate_continuity",
            "passed": not missing_compact_to_scorer and not missing_scorer_to_integrated,
            "detail": "Every compact scorer entry reaches the patient scorer and every forwarded scorer candidate reaches integrated output.",
        },
        {
            "gate": "heldout_patient_output_coverage",
            "passed": {row["patient_id"] for row in integrated_rows} == heldout_ids
            and {row["patient_id"] for row in report_rows} == heldout_ids,
            "detail": "All 18 held-out patients have integrated and complete-report output.",
        },
    ]
    core_pass = all(gate["passed"] for gate in gates)
    summary = {
        "schema_version": "heldout_v2_6_validation_audit.v1",
        "answer_blind": True,
        "release_root": str(release_root.resolve()),
        "release_summary_sha256": sha256_file(release_root / "summary.json"),
        "inventory_root": str(inventory_root.resolve()),
        "inventory_summary_sha256": sha256_file(
            inventory_root / "inventory_summary.json"
        ),
        "cohort": {
            "development_overlap_count": len(development_ids),
            "heldout_patient_count": len(heldout_ids),
            "heldout_patient_ids": sorted(map(int, heldout_ids)),
            "sets_disjoint": not bool(development_ids & heldout_ids),
        },
        "source_contract": {
            "heldout_test_count": len(test_rows),
            "nucleic_type_counts": dict(Counter(row["nucleic_type"] for row in test_rows)),
            "test_qc_status_counts": dict(Counter(row["qc_status"] for row in test_rows)),
            "positive_observation_count": len(observation_rows),
            "candidate_union_count": len(union_rows),
            "tests_with_input_reads": len(test_rows) - len(missing_denominator_rows),
            "tests_without_input_reads": len(missing_denominator_rows),
            "missing_denominator_without_explicit_marker": len(
                missing_denominator_without_marker
            ),
            "protocol_condition_annotated_test_count": len(protocol_annotated),
            "technical_repeat_biological_independence_known": False,
        },
        "output": {
            "compact_candidate_count": len(compact_rows),
            "scorer_patient_organism_count": len(scorer_rows),
            "scorer_not_forwarded_count": len(scorer_not_forwarded),
            "integrated_candidate_count": len(integrated_rows),
            "complete_report_candidate_count": len(report_rows),
            "final_tier_counts": dict(
                Counter(row["final_reporting_tier"] for row in integrated_rows)
            ),
            "reporting_role_counts": dict(
                Counter(row["reporting_role"] for row in integrated_rows)
            ),
            "missing_compact_to_scorer": len(missing_compact_to_scorer),
            "missing_scorer_to_integrated": len(missing_scorer_to_integrated),
            "integrated_hospital_only_or_noncompact": len(integrated_hospital_only),
        },
        "gates": gates,
        "all_core_gates_pass": core_pass,
        "validation_scope": {
            "operational_transportability_ready": core_pass,
            "rpm_guardrail_evaluable_on_all_tests": not missing_denominator_rows,
            "technical_repeat_independence_validation_ready": False,
            "precision_recall_f1_ready": False,
            "reason_metrics_pending": (
                "No independent physician-answer file exists for the 18 held-out patients."
            ),
        },
    }
    return summary, gates, patient_rows, disappearance_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-root", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--inventory-root", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary, gates, patients, disappearance = audit(
        args.release_root, args.inventory_root
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json(args.output_dir / "summary.json", summary)
    write_csv(
        args.output_dir / "validation_gates.csv",
        gates,
        ["gate", "passed", "detail"],
    )
    write_csv(
        args.output_dir / "heldout_patient_output_counts.csv",
        patients,
        [
            "patient_id",
            "tests",
            "tests_without_input_reads",
            "compact_scorer_entry",
            "compact_clinical_review",
            "compact_qc_only",
            "Picked",
            "Possible",
            "Fallback-Possible",
            "Context",
            "complete_report",
        ],
    )
    write_csv(
        args.output_dir / "candidate_disappearance_audit.csv",
        disappearance,
        ["patient_id", "organism_name_key", "stage"],
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

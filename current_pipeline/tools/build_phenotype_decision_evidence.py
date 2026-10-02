"""Build a lossless, queryable shadow ledger from imported phenotype workbooks.

The ledger retains every normalized worksheet row. Indexes are pointers to rows,
not clinical decisions. Nothing here changes pathogen tiers or benchmark answers.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools.import_phenotype_context import (
    EXPECTED_PHENOTYPES,
    PHENOTYPE_GROUPS,
    REQUIRED_SHEETS,
    file_sha256,
    phenotype_model_role,
)


SCHEMA_VERSION = "phenotype_decision_evidence_v1.0"


def key(value: Any) -> str:
    return str(value).strip() if value is not None and str(value).strip() else "__MISSING__"


def build_ledger(full: dict[str, Any], context: dict[str, Any], full_hash: str) -> dict[str, Any]:
    if full["patient_id"] != context["patient_id"]:
        raise ValueError("patient_id differs between full evidence and compact context")
    if full["source_workbook_sha256"] != context["source"]["workbook_sha256"]:
        raise ValueError("workbook SHA-256 differs between full evidence and compact context")
    if set(full["sheets"]) != set(REQUIRED_SHEETS):
        raise ValueError("full evidence does not contain exactly the eight expected sheets")

    records: dict[str, list[dict[str, Any]]] = {}
    index: dict[str, dict[str, list[str]]] = {
        name: defaultdict(list)
        for name in (
            "phenotype_group", "phenotype", "phenotype_role", "phenotype_qc",
            "fact_type", "timeline_event_type", "organism", "specimen",
            "lab_test", "lab_date", "temporal_relation", "source_type",
            "source_date", "evidence_page",
        )
    }
    seen_ids: set[str] = set()

    for sheet in REQUIRED_SHEETS:
        source_rows = full["sheets"][sheet]
        if not isinstance(source_rows, list):
            raise ValueError(f"{sheet} is not a list")
        records[sheet] = []
        for row in source_rows:
            source_row = row.get("_source_row")
            if not isinstance(source_row, int) or source_row < 2:
                raise ValueError(f"{sheet} has an invalid source row")
            evidence_id = f"{sheet}:{source_row}"
            if evidence_id in seen_ids:
                raise ValueError(f"duplicate evidence ID: {evidence_id}")
            seen_ids.add(evidence_id)
            records[sheet].append({"id": evidence_id, "row": row})

            if sheet == "Phenotypes_Long":
                phenotype = key(row.get("phenotype"))
                group = next(
                    (name for name, members in PHENOTYPE_GROUPS.items() if phenotype in members),
                    "unrecognized",
                )
                index["phenotype_group"][group].append(evidence_id)
                index["phenotype"][phenotype].append(evidence_id)
                index["phenotype_role"][phenotype_model_role(row)].append(evidence_id)
                index["phenotype_qc"][key(row.get("phenotype_qc_status"))].append(evidence_id)
            elif sheet == "Clinical_Facts":
                index["fact_type"][key(row.get("fact_type"))].append(evidence_id)
            elif sheet == "Clinical_Timeline":
                index["timeline_event_type"][key(row.get("event_type"))].append(evidence_id)
            elif sheet == "Microbiology":
                index["organism"][key(row.get("organism"))].append(evidence_id)
                index["specimen"][key(row.get("specimen"))].append(evidence_id)
            elif sheet == "Lab_Results":
                index["lab_test"][key(row.get("test_name"))].append(evidence_id)
                index["lab_date"][key(row.get("date"))].append(evidence_id)

            if row.get("temporal_relation") is not None:
                index["temporal_relation"][key(row["temporal_relation"])].append(evidence_id)
            if row.get("source_type") is not None:
                index["source_type"][key(row["source_type"])].append(evidence_id)
            if row.get("source_date") is not None:
                index["source_date"][key(row["source_date"])].append(evidence_id)
            if row.get("evidence_page") is not None:
                index["evidence_page"][key(row["evidence_page"])].append(evidence_id)

    phenotype_names = [item["row"].get("phenotype") for item in records["Phenotypes_Long"]]
    missing_phenotypes = sorted(set(EXPECTED_PHENOTYPES) - set(phenotype_names))
    duplicate_phenotypes = sorted(name for name, count in Counter(phenotype_names).items() if count > 1)
    if missing_phenotypes or duplicate_phenotypes:
        raise ValueError(
            f"incomplete phenotype contract: missing={missing_phenotypes}, duplicate={duplicate_phenotypes}"
        )

    source_counts = {sheet: len(full["sheets"][sheet]) for sheet in REQUIRED_SHEETS}
    output_counts = {sheet: len(records[sheet]) for sheet in REQUIRED_SHEETS}
    if source_counts != output_counts or sum(source_counts.values()) != len(seen_ids):
        raise ValueError("source rows were not preserved one-to-one")

    review_ids = [item["id"] for item in records["Review_Only"]]
    review_ids.extend(index["phenotype_role"].get("manual_review", []))
    return {
        "schema_version": SCHEMA_VERSION,
        "patient_id": full["patient_id"],
        "patient_number": context["patient_number"],
        "source": {
            "workbook_file": full["source_workbook"],
            "workbook_sha256": full["source_workbook_sha256"],
            "full_evidence_sha256": full_hash,
            "patient_qc_status": context["source"].get("patient_qc_status"),
            "system_issues": context["source"].get("system_issues", []),
            "source_pdf_available": False,
            "patient_identity_link_to_raw_verified": False,
        },
        "records": records,
        "index": {name: dict(sorted(buckets.items())) for name, buckets in index.items()},
        "review_queue_ids": list(dict.fromkeys(review_ids)),
        "integrity": {
            "source_row_counts": source_counts,
            "output_row_counts": output_counts,
            "total_rows": len(seen_ids),
            "all_normalized_rows_preserved": True,
        },
        "use_constraints": [
            "Shadow evidence only; do not automatically select or re-tier a pathogen.",
            "The source PDF was not supplied; workbook-to-PDF completeness is unverified.",
            "Raw patient identity linkage and index mNGS time require separate verification.",
            "Phenotype NO means not established, not proof of absence.",
            "AFTER_ONSET evidence is not a baseline etiologic prior.",
            "Lab rows may include parser artifacts and must be validated before clinical interpretation.",
            "Treatment response is not a pathogen-specific causal link without timing and competing-therapy review.",
            "Benchmark answers are not read by this builder; source narratives may still contain diagnoses.",
        ],
    }


def build_all(input_root: Path, output_root: Path, workbook_root: Path | None = None) -> dict[str, Any]:
    full_paths = sorted(input_root.glob("patient_*_phenotype_evidence_full.json"))
    if not full_paths:
        raise ValueError(f"no phenotype full-evidence files in {input_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    patients: list[dict[str, Any]] = []
    for full_path in full_paths:
        context_path = full_path.with_name(
            full_path.name.replace("_phenotype_evidence_full.json", "_phenotype_context.json")
        )
        if not context_path.is_file():
            raise ValueError(f"missing context: {context_path}")
        full = json.loads(full_path.read_text(encoding="utf-8"))
        context = json.loads(context_path.read_text(encoding="utf-8"))
        ledger = build_ledger(full, context, file_sha256(full_path))
        if workbook_root is not None:
            workbook_path = workbook_root / ledger["source"]["workbook_file"]
            if not workbook_path.is_file():
                raise ValueError(f"missing workbook: {workbook_path}")
            if file_sha256(workbook_path) != ledger["source"]["workbook_sha256"]:
                raise ValueError(f"workbook changed since import: {workbook_path}")
            ledger["source"]["workbook_sha256_verified"] = True
        else:
            ledger["source"]["workbook_sha256_verified"] = False
        output_path = output_root / full_path.name.replace(
            "_phenotype_evidence_full.json", "_phenotype_decision_evidence_v1.json"
        )
        output_path.write_text(json.dumps(ledger, ensure_ascii=False, indent=2), encoding="utf-8")
        patients.append({
            "patient_id": ledger["patient_id"],
            "file": output_path.name,
            "row_counts": ledger["integrity"]["output_row_counts"],
            "total_rows": ledger["integrity"]["total_rows"],
            "review_queue_count": len(ledger["review_queue_ids"]),
            "patient_qc_status": ledger["source"]["patient_qc_status"],
            "sha256": file_sha256(output_path),
        })
    totals = {
        sheet: sum(patient["row_counts"][sheet] for patient in patients)
        for sheet in REQUIRED_SHEETS
    }
    audit = {
        "schema_version": SCHEMA_VERSION,
        "patient_count": len(patients),
        "row_counts": totals,
        "total_rows": sum(totals.values()),
        "workbook_sha256_verified": workbook_root is not None,
        "patients": patients,
    }
    (output_root / "decision_evidence_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workbook-root", type=Path)
    args = parser.parse_args()
    audit = build_all(args.input_root, args.output, args.workbook_root)
    print(json.dumps({key: value for key, value in audit.items() if key != "patients"}, indent=2))


if __name__ == "__main__":
    main()

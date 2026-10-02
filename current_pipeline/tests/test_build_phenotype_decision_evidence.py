from __future__ import annotations

import unittest

from tools.build_phenotype_decision_evidence import build_ledger
from tools.import_phenotype_context import EXPECTED_PHENOTYPES, REQUIRED_SHEETS


def fixture() -> tuple[dict, dict]:
    sheets = {name: [] for name in REQUIRED_SHEETS}
    sheets["Phenotypes_Wide"] = [{"_source_row": 2, "patient_qc_status": "SYSTEM_REVIEW_REQUIRED"}]
    sheets["Phenotypes_Long"] = [
        {
            "_source_row": number + 2,
            "phenotype": name,
            "status": "YES" if number == 0 else "NO",
            "confidence": "HIGH" if number == 0 else "LOW",
            "temporal_relation": "PRE_PNEUMONIA" if number == 0 else "UNCLEAR",
            "phenotype_qc_status": "PASS",
        }
        for number, name in enumerate(EXPECTED_PHENOTYPES)
    ]
    sheets["Review_Only"] = [{"_source_row": 2, "review_reason": "source needs review"}]
    sheets["Microbiology"] = [{"_source_row": 2, "organism": "Example species", "specimen": "BAL"}]
    sheets["Clinical_Timeline"] = [{"_source_row": 2, "event_type": "MEDICATION", "date": "2026-01-02"}]
    sheets["Lab_Results"] = [
        {"_source_row": 2, "test_name": "WBC", "value": "4.2", "date": "2026-01-03"},
        {"_source_row": 3, "test_name": "Unrecognized assay", "value": "text", "date": None},
    ]
    sheets["Evidence"] = [{"_source_row": 2, "evidence": "original evidence", "evidence_page": 7}]
    sheets["Clinical_Facts"] = [
        {"_source_row": 2, "fact_type": "TREATMENT", "temporal_relation": "AFTER_ONSET"}
    ]
    full = {
        "patient_id": "patient 1",
        "source_workbook": "patient 1_phenotype.xlsx",
        "source_workbook_sha256": "a" * 64,
        "sheets": sheets,
    }
    context = {
        "patient_id": "patient 1",
        "patient_number": 1,
        "source": {
            "workbook_sha256": "a" * 64,
            "patient_qc_status": "SYSTEM_REVIEW_REQUIRED",
            "system_issues": ["example"],
        },
    }
    return full, context


class PhenotypeDecisionEvidenceTests(unittest.TestCase):
    def test_preserves_all_rows_and_indexes_without_dropping_unknown_labs(self) -> None:
        full, context = fixture()
        ledger = build_ledger(full, context, "b" * 64)
        for sheet, rows in full["sheets"].items():
            self.assertEqual(rows, [item["row"] for item in ledger["records"][sheet]])
        self.assertEqual(["Lab_Results:3"], ledger["index"]["lab_test"]["Unrecognized assay"])
        self.assertEqual(["Lab_Results:3"], ledger["index"]["lab_date"]["__MISSING__"])
        self.assertIn("Review_Only:2", ledger["review_queue_ids"])
        self.assertEqual(["Clinical_Facts:2"], ledger["index"]["temporal_relation"]["AFTER_ONSET"])
        self.assertEqual("SYSTEM_REVIEW_REQUIRED", ledger["source"]["patient_qc_status"])

    def test_rejects_mismatched_source_hash(self) -> None:
        full, context = fixture()
        context["source"]["workbook_sha256"] = "c" * 64
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            build_ledger(full, context, "b" * 64)

    def test_rejects_duplicate_row_id(self) -> None:
        full, context = fixture()
        full["sheets"]["Lab_Results"][1]["_source_row"] = 2
        with self.assertRaisesRegex(ValueError, "duplicate evidence ID"):
            build_ledger(full, context, "b" * 64)

    def test_rejects_missing_phenotype(self) -> None:
        full, context = fixture()
        full["sheets"]["Phenotypes_Long"].pop()
        with self.assertRaisesRegex(ValueError, "incomplete phenotype contract"):
            build_ledger(full, context, "b" * 64)


if __name__ == "__main__":
    unittest.main()

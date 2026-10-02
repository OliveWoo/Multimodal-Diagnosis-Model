from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools import candidate_phenotype_context as phenotype_context
from tools import import_phenotype_context as importer


def phenotype_record(
    index: int,
    name: str,
    *,
    status: str = "NO",
    confidence: str = "HIGH",
    temporal_relation: str = "PRE_PNEUMONIA",
    qc: str = "PASS",
) -> dict[str, object]:
    return {
        "id": f"Phenotypes_Long:{index}",
        "row": {
            "phenotype": name,
            "status": status,
            "confidence": confidence,
            "temporal_relation": temporal_relation,
            "phenotype_qc_status": qc,
            "reason": f"reason for {name}",
            "evidence": f"evidence for {name}",
            "source_type": "PROGRESS_NOTE",
            "source_date": "2024-01-01",
        },
    }


def synthetic_ledger(patient_number: int = 9) -> dict[str, object]:
    records = [
        phenotype_record(index, name)
        for index, name in enumerate(importer.EXPECTED_PHENOTYPES, start=1)
    ]
    for item in records:
        row = item["row"]
        if row["phenotype"] == "IMMUNOCOMPROMISED":
            row["status"] = "YES"
        if row["phenotype"] == "FUNGAL_INFECTION_RISK":
            row["status"] = "YES"
            row["temporal_relation"] = "AFTER_ONSET"
    return {
        "schema_version": "phenotype_decision_evidence_v1.0",
        "patient_id": f"patient {patient_number}",
        "source": {
            "workbook_file": f"patient {patient_number}_phenotype.xlsx",
            "workbook_sha256": "abc123",
            "patient_qc_status": "PASS",
        },
        "records": {
            "Phenotypes_Long": records,
            "Evidence": [],
            "Review_Only": [],
        },
    }


class CandidatePhenotypeContextTests(unittest.TestCase):
    def write_ledger(self, root: Path, patient_number: int = 9) -> None:
        path = root / f"patient_{patient_number}_phenotype_decision_evidence_v1.json"
        path.write_text(
            json.dumps(synthetic_ledger(patient_number), ensure_ascii=False),
            encoding="utf-8",
        )

    def test_all_22_are_retained_but_only_family_relevant_rows_are_expanded(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.write_ledger(root)
            shadow = phenotype_context.build_patient_shadow(
                root,
                patient_id="9",
                candidates=[{"organism_name": "Pneumocystis jirovecii"}],
                linkage_mode="same_number_shadow",
            )

        self.assertEqual(22, len(shadow["all_phenotype_overview"]))
        candidate = shadow["candidate_contexts"][0]
        self.assertEqual(
            {"IMMUNOCOMPROMISED", "FUNGAL_INFECTION_RISK"},
            {item["phenotype"] for item in candidate["expanded_phenotypes"]},
        )

    def test_no_is_not_negative_and_after_onset_is_not_an_etiologic_prior(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.write_ledger(root)
            shadow = phenotype_context.build_patient_shadow(
                root,
                patient_id="9",
                candidates=[{"organism_name": "Pneumocystis jirovecii"}],
                linkage_mode="same_number_shadow",
            )

        expanded = {
            item["phenotype"]: item for item in shadow["candidate_contexts"][0]["expanded_phenotypes"]
        }
        self.assertEqual("etiologic_prior", expanded["IMMUNOCOMPROMISED"]["model_role"])
        self.assertTrue(expanded["IMMUNOCOMPROMISED"]["eligible_etiologic_prior"])
        self.assertEqual("outcome_or_context", expanded["FUNGAL_INFECTION_RISK"]["model_role"])
        self.assertFalse(expanded["FUNGAL_INFECTION_RISK"]["eligible_etiologic_prior"])
        no_rows = [item for item in shadow["all_phenotype_overview"] if item["status"] == "NO"]
        self.assertTrue(no_rows)
        self.assertTrue(all(item["model_role"] == "not_established" for item in no_rows))

    def test_same_number_link_is_shadow_only(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            self.write_ledger(root)
            shadow = phenotype_context.build_patient_shadow(
                root,
                patient_id="9",
                candidates=[],
                linkage_mode="same_number_shadow",
            )

        self.assertEqual("provisional_same_number_link", shadow["linkage"]["status"])
        self.assertTrue(shadow["linkage"]["shadow_reasoning_allowed"])
        self.assertFalse(shadow["linkage"]["production_link_verified"])

    def test_missing_ledger_is_explicit_and_does_not_invent_context(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            shadow = phenotype_context.build_patient_shadow(
                Path(temp_dir),
                patient_id="14",
                candidates=[{"organism_name": "Pneumocystis jirovecii"}],
                linkage_mode="same_number_shadow",
            )

        self.assertEqual("phenotype_file_missing", shadow["linkage"]["status"])
        self.assertFalse(shadow["linkage"]["shadow_reasoning_allowed"])
        self.assertEqual([], shadow["all_phenotype_overview"])
        self.assertEqual(
            "phenotype_file_missing",
            shadow["candidate_contexts"][0]["summary"]["unavailable_reason"],
        )


if __name__ == "__main__":
    unittest.main()

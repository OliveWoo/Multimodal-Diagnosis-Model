from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from core import llm_parser
from tools import source_provenance
from tools.audit_source_to_normalized_contract import (
    audit_pair,
    infer_text_issue_stage,
    text_integrity_summary,
)


class SourceProvenanceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.raw = {
            "source_file": "NGS patient 1.xlsx",
            "sheet": "culture",
            "raw_data": [
                ["2024/01/02 08:30", "BALF", "Candida albicans", ">10^5 CFU/mL"],
                ["20240103", "Blood", "No growth", ""],
            ],
        }

    def test_provenance_has_stable_hash_and_anchors(self) -> None:
        first = source_provenance.build_source_provenance(self.raw)
        second = source_provenance.build_source_provenance(self.raw)
        self.assertEqual(first["raw_payload_sha256"], second["raw_payload_sha256"])
        self.assertEqual(2, first["raw_row_count"])
        self.assertEqual(7, first["raw_nonempty_cell_count"])
        self.assertEqual(2, len(first["critical_anchor_index"]["dates"]))
        self.assertEqual(1, len(first["critical_anchor_index"]["quantities"]))
        self.assertEqual(
            {"balf", "blood"},
            {item["canonical"] for item in first["critical_anchor_index"]["specimens"]},
        )

    def test_save_result_attaches_raw_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            raw_path = root / "NGS_patient_1_culture_Raw.json"
            raw_path.write_text("{}", encoding="utf-8")
            output_path = llm_parser.save_result(
                {"records": []}, raw_path, root / "normalized", self.raw
            )
            payload = llm_parser.load_json(output_path)
            provenance = payload["_source_provenance"]
            self.assertEqual("culture", provenance["sheet"])
            self.assertEqual("none_provenance_only", provenance["tier_effect"])

    def test_pair_audit_finds_missing_date_without_changing_data(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            raw_path = root / "NGS_patient_1_culture_Raw.json"
            normalized_path = root / "NGS_patient_1_culture.json"
            raw_path.write_text(json.dumps(self.raw), encoding="utf-8")
            normalized_path.write_text(
                json.dumps(
                    [{"specimen_type": "BALF", "quantity": ">10^5 CFU/mL"}]
                ),
                encoding="utf-8",
            )
            result = audit_pair(raw_path, normalized_path)
            self.assertEqual("paired", result["status"])
            self.assertEqual(2, result["anchor_summary"]["dates"]["missing_count"])
            self.assertEqual(1, result["anchor_summary"]["quantities"]["retained_count"])

    def test_text_integrity_detects_private_use_and_replacement_chars(self) -> None:
        result = text_integrity_summary(
            {"raw_data": [["clean"], ["broken \ue123 text \ufffd"]]}
        )
        self.assertTrue(result["has_suspicious_text"])
        self.assertEqual(1, result["suspicious_node_count"])
        self.assertEqual(1, result["issue_counts"]["private_use_character"])
        self.assertEqual(1, result["issue_counts"]["unicode_replacement_character"])

    def test_text_integrity_does_not_flag_normal_chinese(self) -> None:
        result = text_integrity_summary({"text": "下呼吸道檢體培養未檢出病原菌"})
        self.assertFalse(result["has_suspicious_text"])

    def test_text_issue_stage_separates_raw_and_normalized(self) -> None:
        clean = {"has_suspicious_text": False}
        suspicious = {"has_suspicious_text": True}
        self.assertEqual(
            "introduced_or_exposed_during_normalization",
            infer_text_issue_stage(clean, suspicious),
        )
        self.assertEqual(
            "present_in_raw_and_normalized",
            infer_text_issue_stage(suspicious, suspicious),
        )


if __name__ == "__main__":
    unittest.main()

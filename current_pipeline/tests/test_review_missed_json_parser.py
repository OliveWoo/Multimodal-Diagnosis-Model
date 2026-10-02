from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tools import review_missed_mngs_candidates as review


class ReviewMissedJsonParserTests(unittest.TestCase):
    def test_compact_phenotype_event_view_preserves_evidence_without_repeated_hashes(self) -> None:
        value = {
            "schema_version": "test.v1",
            "patient_id": 22,
            "organism_name": "Stenotrophomonas maltophilia",
            "canonical_key": "stenotrophomonasmaltophilia",
            "phenotype_status": "provisional_link",
            "source": {
                "source_workbook": "patient 22_phenotype.xlsx",
                "source_workbook_sha256": "secret-hash",
                "source_pdf_available": False,
                "patient_qc_status": "SYSTEM_REVIEW_REQUIRED",
            },
            "event_views": [{
                "index_event": {
                    "event_id": "P22:mngs:1",
                    "specimen_site": "BALF",
                    "source_ref": {
                        "file": r"D:\full\path\NGS_patient_22_mNGS.json",
                        "sha256": "index-hash",
                        "json_path": "records[0]",
                    },
                },
                "linkage": {
                    "patient_identity_assessment": "local_strong_same_patient_support",
                    "duplicate_export_cluster": [],
                    "production_link_verified": False,
                },
                "phenotype_context": [{
                    "source_id": "Phenotypes_Long:14",
                    "phenotype": "HAP_PHENOTYPE",
                    "status": "YES",
                    "confidence": "HIGH",
                    "reason": "Hospital-acquired timing",
                    "evidence": "Pneumonia in both lungs.",
                    "event_alignment": {
                        "source_id": "Phenotypes_Long:14",
                        "source_ref": {
                            "ledger_file": r"D:\full\path\patient_22_ledger.json",
                            "ledger_sha256": "ledger-hash",
                            "sheet": "Phenotypes_Long",
                            "source_row": 14,
                        },
                        "source_text": "Pneumonia in both lungs.",
                    },
                }, {
                    "source_id": "Phenotypes_Long:6",
                    "phenotype": "RECENT_HOSPITALIZATION_90D",
                    "status": "NO",
                    "confidence": "LOW",
                    "model_role": "not_established",
                    "event_alignment": {
                        "source_id": "Phenotypes_Long:6",
                        "source_ref": {
                            "ledger_file": r"D:\full\path\patient_22_ledger.json",
                            "ledger_sha256": "ledger-hash",
                            "sheet": "Phenotypes_Long",
                            "source_row": 6,
                        },
                    },
                }],
                "dated_clinical_facts": [{
                    "source_id": "Clinical_Facts:8",
                    "description": "Bilateral pleural effusion",
                    "source_text": "Pneumonia in both lungs. Bilateral pleural effusion.",
                    "source_ref": {
                        "ledger_file": r"D:\full\path\patient_22_ledger.json",
                        "ledger_sha256": "ledger-hash",
                        "sheet": "Clinical_Facts",
                        "source_row": 8,
                    },
                }],
                "exact_organism_workbook_microbiology": [],
                "dated_clinical_timeline": [],
                "lab_parser_rows_in_packet": 193,
                "lab_values_supplied_to_model": False,
                "model_use": "provisional_shadow_review_only_event_unverified",
            }],
            "constraints": ["No change to deterministic Picked is permitted."],
        }

        compact = review.compact_phenotype_event_view(value)
        encoded = json.dumps(compact, ensure_ascii=False)

        self.assertIn("Pneumonia in both lungs.", encoded)
        self.assertIn("Bilateral pleural effusion", encoded)
        self.assertIn("HAP_PHENOTYPE", encoded)
        self.assertNotIn("RECENT_HOSPITALIZATION_90D", encoded)
        self.assertEqual(1, compact["event_views"][0]["unsupported_low_confidence_no_omitted"])
        self.assertIn("patient_22_ledger.json", encoded)
        self.assertNotIn("D:\\\\full\\\\path", encoded)
        self.assertNotIn("secret-hash", encoded)
        self.assertNotIn("ledger-hash", encoded)
        self.assertNotIn("index-hash", encoded)

    def test_codex_cli_backend_uses_login_without_api_key(self) -> None:
        payload = {
            "patient_id": "P1",
            "review_purpose": "safety_review_only",
            "overall_assessment_zh": "測試",
            "review_high_priority": [],
            "review_context_needed": [],
            "review_low_specificity": [],
            "omitted_with_reason": [],
            "oral_aspiration_flora_pattern": {
                "enabled": True,
                "pattern_detected": False,
                "representative_organisms": [],
                "grouped_non_representative_organisms": [],
                "interpretation_zh": "",
            },
            "llm_limitations": [],
        }

        def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
            self.assertNotIn("OPENAI_API_KEY", kwargs["env"])
            self.assertIn("--output-schema", command)
            result_path = Path(command[command.index("--output-last-message") + 1])
            result_path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
            return SimpleNamespace(returncode=0, stderr="", stdout="")

        with tempfile.TemporaryDirectory() as tmp, patch.object(
            review, "find_codex_executable", return_value="codex.exe"
        ), patch.object(review.subprocess, "run", side_effect=fake_run), patch.dict(
            os.environ, {"OPENAI_API_KEY": "invalid-key"}
        ):
            result_path = Path(tmp) / "result.tmp"
            text, usage = review.send_to_codex_cli(
                "prompt",
                model="gpt-5.6-luna",
                reasoning_effort="medium",
                result_path=result_path,
            )

        self.assertEqual(payload, json.loads(text))
        self.assertIsNone(usage)
        self.assertFalse(result_path.exists())
        self.assertFalse(result_path.with_suffix(".schema.json").exists())

    def test_codex_cli_recovery_mode_omits_output_schema(self) -> None:
        payload = {"patient_id": "P22", "review_purpose": "safety_review_only"}

        def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
            self.assertNotIn("--output-schema", command)
            result_path = Path(command[command.index("--output-last-message") + 1])
            result_path.write_text(json.dumps(payload), encoding="utf-8")
            return SimpleNamespace(returncode=0, stderr="", stdout="")

        with tempfile.TemporaryDirectory() as tmp, patch.object(
            review, "find_codex_executable", return_value="codex.exe"
        ), patch.object(review.subprocess, "run", side_effect=fake_run):
            result_path = Path(tmp) / "result.tmp"
            text, usage = review.send_to_codex_cli(
                "prompt",
                model="gpt-6-luna",
                reasoning_effort="low",
                result_path=result_path,
                use_output_schema=False,
            )

        self.assertEqual(payload, json.loads(text))
        self.assertIsNone(usage)
        self.assertFalse(result_path.exists())
        self.assertFalse(result_path.with_suffix(".schema.json").exists())

    def test_repairs_trailing_commas_outside_strings(self) -> None:
        payload = review.parse_json_response(
            '{"items": [1, 2,], "nested": {"value": 3,}, '
            '"literal": "keep comma, }",}'
        )

        self.assertEqual([1, 2], payload["items"])
        self.assertEqual({"value": 3}, payload["nested"])
        self.assertEqual("keep comma, }", payload["literal"])

    def test_repairs_unterminated_container_with_trailing_comma(self) -> None:
        payload = review.parse_json_response('{"items": [1, 2],')

        self.assertEqual([1, 2], payload["items"])

    def test_d3_hospital_gnb_is_retained_as_context(self) -> None:
        tier, reason = review.direct_review_tier_for_item(
            {"organism_name": "Klebsiella pneumoniae", "confidence": "moderate"},
            source_section="review_high_priority",
            evidence={
                "organism_name": "Klebsiella pneumoniae",
                "classification": "Bacterial",
                "reads": 703,
                "reads_tier": "R2_medium",
                "reads_percentile": 1.0,
                "dominance_tier": "D3_dominant",
                "specimen_class": "S2_lower_respiratory",
            },
            payload={"hospital_side_summary": {}, "mngs_to_specimen": {}},
            picked_keys=set(),
            picked_genera=set(),
        )

        self.assertEqual("review_context_needed", tier)
        self.assertIn("D3-dominant", reason)


if __name__ == "__main__":
    unittest.main()

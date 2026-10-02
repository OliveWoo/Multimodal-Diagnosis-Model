from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from tools.auto_decide_organism_taxonomy import (
    _assert_answer_blind,
    apply_staged_overlay,
    clinical_prompt_for,
    ensure_codex_cli_authenticated,
    evaluate_reviews,
    identity_prompt_for,
    run_decisions,
    validate_identity_review,
)


def _item() -> dict:
    return {
        "organism_name": "Novelibacter testii",
        "canonical_key": "novelibactertestii",
        "reported_names": ["Novelibacter testii"],
        "input_classes": ["bacterium"],
        "ncbi": {
            "resolution_status": "exact_scientific_name",
            "taxid": "12345",
            "scientific_name": "Novelibacter testii",
            "rank": "species",
            "division": "Bacteria",
            "lineage": "cellular organisms; Bacteria; Novelibacter",
        },
    }


def _identity() -> dict:
    return {
        "decision": "verified",
        "canonical_name": "Novelibacter testii",
        "taxid": 12345,
        "taxonomic_rank": "species",
        "biological_class": "bacterium",
        "name_relation": "exact",
        "confidence": "high",
        "rationale": "All supplied identity fields agree.",
        "conflict_notes": [],
    }


def _clinical() -> dict:
    return {
        "decision": "classify",
        "primary_rule_family": "environmental_low_specificity",
        "secondary_rule_families": [],
        "clinical_traits": ["environmental_association", "direct_support_preferred"],
        "supporting_pmids": ["42"],
        "counterevidence_pmids": [],
        "evidence_strength": "sufficient",
        "confidence": "high",
        "human_infection_evidence": "Rare human isolation is documented.",
        "specimen_relevance": "Direct confirmation is preferred.",
        "colonization_or_contamination_caveat": "Environmental contamination remains possible.",
        "rationale": "The supplied source supports an environmental low-specificity default.",
        "uncertainty_reasons": [],
    }


def _literature() -> dict:
    return {"articles": [{"pmid": "42", "title": "Clinical report", "abstract": "Human isolation."}]}


class AutoTaxonomyDecisionTests(unittest.TestCase):
    def test_prompts_are_answer_blind(self) -> None:
        identity_messages = identity_prompt_for(_item())
        _assert_answer_blind(identity_messages)
        clinical_messages = clinical_prompt_for(
            _item(),
            {
                "canonical_name": "Novelibacter testii",
                "taxid": 12345,
                "taxonomic_rank": "species",
                "biological_class": "bacterium",
            },
            _literature(),
        )
        _assert_answer_blind(clinical_messages)
        text = json.dumps(identity_messages + clinical_messages).lower()
        self.assertNotIn("patient_id", text)
        self.assertNotIn("benchmark_answer", text)

    def test_identity_conflict_forces_provisional(self) -> None:
        review = _identity()
        review["taxid"] = 999
        status, reasons, identity = validate_identity_review(_item(), review)
        self.assertEqual("provisional_unmapped", status)
        self.assertIn("LLM_TAXID_CONFLICT", reasons)
        self.assertIsNone(identity)

    def test_high_confidence_supported_review_is_auto_approved(self) -> None:
        with patch(
            "tools.auto_decide_organism_taxonomy.classify_organism",
            return_value={"mapping_status": "unmapped"},
        ):
            verdict, profile = evaluate_reviews(
                _item(),
                _literature(),
                _identity(),
                _clinical(),
                identity_model="identity-model",
                clinical_model="clinical-model",
                decision_input_sha256="abc",
            )
        self.assertEqual("auto_approved", verdict["status"])
        self.assertEqual("environmental_low_specificity", profile["primary_rule_family"])
        self.assertEqual(12345, profile["taxid"])

    def test_unavailable_pmid_forces_provisional(self) -> None:
        clinical = _clinical()
        clinical["supporting_pmids"] = ["999"]
        with patch(
            "tools.auto_decide_organism_taxonomy.classify_organism",
            return_value={"mapping_status": "unmapped"},
        ):
            verdict, profile = evaluate_reviews(
                _item(),
                _literature(),
                _identity(),
                clinical,
                identity_model="identity-model",
                clinical_model="clinical-model",
                decision_input_sha256="abc",
            )
        self.assertEqual("provisional_unmapped", verdict["status"])
        self.assertIn("PMID_NOT_IN_RETRIEVED_SOURCE_PACK", verdict["validator_reasons"])
        self.assertIsNone(profile)

    def test_class_family_conflict_forces_provisional(self) -> None:
        clinical = _clinical()
        clinical["primary_rule_family"] = "candida_or_yeast"
        with patch(
            "tools.auto_decide_organism_taxonomy.classify_organism",
            return_value={"mapping_status": "unmapped"},
        ):
            verdict, _ = evaluate_reviews(
                _item(),
                _literature(),
                _identity(),
                clinical,
                identity_model="identity-model",
                clinical_model="clinical-model",
                decision_input_sha256="abc",
            )
        self.assertIn("BIOLOGICAL_CLASS_FAMILY_CONFLICT", verdict["validator_reasons"])

    def test_preview_calls_no_model_and_skips_existing_exact(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.jsonl"
            queue.write_text(json.dumps(_item()) + "\n", encoding="utf-8")
            literature = root / "literature.json"
            literature.write_text(json.dumps({"novelibactertestii": _literature()}), encoding="utf-8")
            with patch(
                "tools.auto_decide_organism_taxonomy.classify_organism",
                return_value={"mapping_status": "unmapped"},
            ):
                result = run_decisions(
                    queue,
                    literature,
                    root / "out",
                    identity_model="identity-model",
                    clinical_model="clinical-model",
                    backend="codex-cli",
                    execute=False,
                )
            self.assertEqual(0, result["model_calls"])
            self.assertEqual(1, result["eligible"])
            self.assertTrue((root / "out" / "prompt_preview.json").exists())

    def test_only_key_limits_preview_selection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = _item()
            second = {**_item(), "organism_name": "Second bacter testii", "canonical_key": "secondbactertestii"}
            queue = root / "queue.jsonl"
            queue.write_text("\n".join((json.dumps(first), json.dumps(second))) + "\n", encoding="utf-8")
            literature = root / "literature.json"
            literature.write_text("{}", encoding="utf-8")
            with patch(
                "tools.auto_decide_organism_taxonomy.classify_organism",
                return_value={"mapping_status": "unmapped"},
            ):
                result = run_decisions(
                    queue,
                    literature,
                    root / "out",
                    identity_model="identity-model",
                    clinical_model="clinical-model",
                    backend="codex-cli",
                    execute=False,
                    only_keys=("secondbactertestii",),
                )
            self.assertEqual(2, result["eligible"])
            self.assertEqual(1, result["selected"])
            preview = json.loads((root / "out" / "prompt_preview.json").read_text(encoding="utf-8"))
            self.assertEqual("secondbactertestii", preview[0]["canonical_key"])

    def test_execute_timeout_is_recorded_as_partial_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.jsonl"
            queue.write_text(json.dumps(_item()) + "\n", encoding="utf-8")
            literature = root / "literature.json"
            literature.write_text(json.dumps({"novelibactertestii": _literature()}), encoding="utf-8")
            with (
                patch(
                    "tools.auto_decide_organism_taxonomy.classify_organism",
                    return_value={"mapping_status": "unmapped"},
                ),
                patch(
                    "tools.auto_decide_organism_taxonomy._call_codex_cli",
                    side_effect=TimeoutError("timed out"),
                ),
                patch("tools.auto_decide_organism_taxonomy.ensure_codex_cli_authenticated"),
            ):
                result = run_decisions(
                    queue,
                    literature,
                    root / "out",
                    identity_model="identity-model",
                    clinical_model="clinical-model",
                    backend="codex-cli",
                    execute=True,
                    call_timeout_seconds=30,
                )
            self.assertEqual("partial_failure", result["completion_status"])
            self.assertEqual(1, result["model_invocation_attempts_this_run"])
            self.assertEqual(0, result["model_calls_this_run"])
            failure = json.loads(
                (root / "out" / "failures" / "novelibactertestii.json").read_text(encoding="utf-8")
            )
            self.assertEqual("identity_review", failure["stage"])

    def test_codex_cli_authentication_fails_fast_with_login_command(self) -> None:
        status = Mock(returncode=0, stdout="Not logged in\n", stderr="")
        with (
            patch("tools.auto_decide_organism_taxonomy.find_codex_executable", return_value="C:/codex.exe"),
            patch("tools.auto_decide_organism_taxonomy.subprocess.run", return_value=status),
            patch.dict("tools.auto_decide_organism_taxonomy.os.environ", {}, clear=True),
        ):
            with self.assertRaisesRegex(RuntimeError, "login --device-auth"):
                ensure_codex_cli_authenticated()

    def test_apply_requires_matching_safe_audit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            staged = root / "staged.json"
            staged.write_text(
                json.dumps({"schema_version": "organism_taxonomy_auto_v1", "exact_profiles": {}}),
                encoding="utf-8",
            )
            audit = root / "audit.json"
            audit.write_text(
                json.dumps({"safe_to_consider_apply": False, "staged_overlay": str(staged.resolve())}),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "did not mark"):
                apply_staged_overlay(staged, audit)

    def test_safe_audit_can_apply_to_separate_auto_overlay(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            staged = root / "staged.json"
            staged.write_text(
                json.dumps(
                    {
                        "schema_version": "organism_taxonomy_auto_v1",
                        "exact_profiles": {
                            "novelibactertestii": {
                                "rule_id": "TAX-AUTO-NOVELIBACTERTESTII",
                                "taxid": 12345,
                                "taxonomic_rank": "species",
                                "biological_class": "bacterium",
                                "primary_rule_family": "environmental_low_specificity",
                                "clinical_traits": ["environmental_association"],
                            }
                        },
                    }
                ),
                encoding="utf-8",
            )
            audit = root / "audit.json"
            audit.write_text(
                json.dumps({"safe_to_consider_apply": True, "staged_overlay": str(staged.resolve())}),
                encoding="utf-8",
            )
            auto_overlay = root / "organism_taxonomy_auto_v1.json"
            with patch("tools.auto_decide_organism_taxonomy.AUTO_OVERLAY", auto_overlay), patch(
                "tools.auto_decide_organism_taxonomy.classify_organism",
                return_value={"mapping_status": "unmapped"},
            ):
                result = apply_staged_overlay(staged, audit)
            self.assertEqual("applied_after_safe_shadow_audit", result["mode"])
            applied = json.loads(auto_overlay.read_text(encoding="utf-8"))
            self.assertIn("novelibactertestii", applied["exact_profiles"])


if __name__ == "__main__":
    unittest.main()

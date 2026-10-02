import csv
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.approve_organism_taxonomy_reviews import compile_approvals, validate_approval
from tools.build_organism_taxonomy_review_queue import build_queue
from tools.draft_organism_taxonomy_reviews import (
    _safe_api_error,
    find_codex_executable,
    make_draft,
    make_draft_codex_cli,
    prompt_for,
    run,
)
from tools.retrieve_organism_taxonomy_literature import retrieve


class TestOrganismTaxonomyReviewWorkflow(unittest.TestCase):
    def test_finds_codex_desktop_cli_when_not_on_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            executable = Path(tmp) / "OpenAI" / "Codex" / "bin" / "build-id" / "codex.exe"
            executable.parent.mkdir(parents=True)
            executable.touch()
            with patch("tools.draft_organism_taxonomy_reviews.shutil.which", return_value=None), patch.dict(
                "os.environ", {"LOCALAPPDATA": tmp}
            ):
                self.assertEqual(find_codex_executable(), str(executable))

    def test_queue_deduplicates_and_excludes_patient_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "candidates.csv"
            with source.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["patient_id", "organism_name", "category", "max_reads_single_test", "benchmark_answer"])
                writer.writeheader()
                writer.writerow({"patient_id": "SECRET-P1", "organism_name": "Novelibacter testii", "category": "bacterial", "max_reads_single_test": 12345, "benchmark_answer": "SECRET-ANSWER"})
                writer.writerow({"patient_id": "SECRET-P2", "organism_name": "Novelibacter testii", "category": "bacterial", "max_reads_single_test": 67890, "benchmark_answer": "SECRET-ANSWER"})
                writer.writerow({"patient_id": "SECRET-P3", "organism_name": "Staphylococcus aureus", "category": "bacterial", "max_reads_single_test": 99, "benchmark_answer": "SECRET-ANSWER"})
            summary = build_queue(source, root / "out")
            self.assertEqual(summary["unique_unmapped_organisms"], 1)
            queue_text = (root / "out" / "organism_review_queue.jsonl").read_text(encoding="utf-8")
            for forbidden in ("SECRET-P1", "SECRET-P2", "SECRET-ANSWER", "12345", "67890", "patient_id", "reads"):
                self.assertNotIn(forbidden, queue_text)
            item = json.loads(queue_text)
            self.assertEqual(item["canonical_key"], "novelibactertestii")
            prompt_text = json.dumps(prompt_for(item))
            self.assertNotIn("SECRET", prompt_text)

    def test_preview_does_not_call_api(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.jsonl"
            queue.write_text(json.dumps({"organism_name": "Novelibacter testii", "canonical_key": "novelibactertestii", "ncbi": {}}) + "\n", encoding="utf-8")
            result = run(queue, root / "preview", model="unused", limit=1, execute=False)
            self.assertEqual(result["api_calls"], 0)
            self.assertTrue((root / "preview" / "prompt_preview.json").exists())

    def test_malformed_api_key_is_rejected_before_client_creation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.jsonl"
            queue.write_text(json.dumps({"organism_name": "Novelibacter testii", "canonical_key": "novelibactertestii", "ncbi": {}}) + "\n", encoding="utf-8")
            with patch.dict(os.environ, {"OPENAI_API_KEY": '$env:OPENAI_API_KEY = "sk-invalid"'}):
                with self.assertRaisesRegex(RuntimeError, "invalid local format"):
                    run(queue, root / "drafts", model="unused", limit=1, execute=True)

    def test_prompted_key_ignores_bad_environment_key(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.jsonl"
            queue.write_text(json.dumps({"organism_name": "Novelibacter testii", "canonical_key": "novelibactertestii", "ncbi": {}}) + "\n", encoding="utf-8")
            with patch.dict(os.environ, {"OPENAI_API_KEY": "$env:OPENAI_API_KEY = bad"}), patch(
                "tools.draft_organism_taxonomy_reviews.getpass.getpass", return_value="still-not-a-key"
            ):
                with self.assertRaisesRegex(RuntimeError, "invalid local format"):
                    run(queue, root / "drafts", model="unused", limit=1, execute=True, prompt_api_key=True)

    def test_inherited_species_can_be_reviewed_separately(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "candidates.csv"
            with source.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["organism_name", "category"])
                writer.writeheader()
                writer.writerow({"organism_name": "Acinetobacter novelis", "category": "bacterial"})
            self.assertEqual(build_queue(source, root / "default")["unique_unmapped_organisms"], 0)
            summary = build_queue(source, root / "inherited", include_inherited=True)
            self.assertEqual(summary["unique_review_organisms"], 1)
            self.assertEqual(summary["unique_unmapped_organisms"], 0)
            item = json.loads((root / "inherited" / "organism_review_queue.jsonl").read_text(encoding="utf-8"))
            self.assertEqual(item["classification_status"], "genus_inherited")

    def test_queue_accepts_compact_taxonomy_output_without_category_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "all_decisions.csv"
            with source.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "patient_id",
                        "organism_name",
                        "taxonomy_mapping_status",
                        "taxonomy_family",
                    ],
                )
                writer.writeheader()
                writer.writerow({
                    "patient_id": "SECRET-P1",
                    "organism_name": "Novelibacter testii",
                    "taxonomy_mapping_status": "unmapped",
                    "taxonomy_family": "unmapped_or_uncertain",
                })
            summary = build_queue(source, root / "out")
            self.assertEqual(summary["unique_unmapped_organisms"], 1)
            item = json.loads(
                (root / "out" / "organism_review_queue.jsonl").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(item["input_classes"], ["unknown"])
            self.assertNotIn("SECRET-P1", json.dumps(item))

    def test_pubmed_pack_is_organism_only_and_resumable(self):
        class FakePubMed:
            def __init__(self, **_kwargs):
                pass

            def search(self, _query, _count):
                return 1, ["42"]

            def fetch(self, _pmids):
                return [{"pmid": "42", "title": "Novelibacter testii clinical report", "abstract": "Novelibacter testii was studied.", "year": 2020, "publication_types": ["Case Reports"], "url": "https://pubmed.ncbi.nlm.nih.gov/42/"}]

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.jsonl"
            queue.write_text(json.dumps({"organism_name": "Novelibacter testii", "canonical_key": "novelibactertestii", "ncbi": {}, "patient_id": "SECRET-P1"}) + "\n", encoding="utf-8")
            output = root / "literature.json"
            with patch("tools.retrieve_organism_taxonomy_literature.PubMedClient", FakePubMed):
                first = retrieve(queue, output, cache_dir=root / "cache")
                second = retrieve(queue, output, cache_dir=root / "cache")
            self.assertEqual(first["retrieved_this_run"], 1)
            self.assertEqual(second["retrieved_this_run"], 0)
            self.assertNotIn("SECRET-P1", output.read_text(encoding="utf-8"))

    def test_model_cannot_cite_unavailable_pmid(self):
        class FakeResponse:
            status = "completed"
            output_text = json.dumps({
                "identity_assessment": "plausible", "family_hypothesis": "environmental_low_specificity",
                "clinical_traits_hypothesis": [], "human_infection_evidence": "", "specimen_relevance": "",
                "contamination_or_colonization_caveat": "", "verification_questions": [],
                "literature_search_queries": [], "supporting_pmids": ["999"], "confidence": "low",
            })

        class FakeClient:
            class responses:
                @staticmethod
                def create(**_kwargs):
                    return FakeResponse()

        item = {"organism_name": "Novelibacter testii", "canonical_key": "novelibactertestii", "ncbi": {}}
        with self.assertRaisesRegex(ValueError, "unavailable PMID"):
            make_draft(item, "test-model", FakeClient(), {"articles": [{"pmid": "42"}]})

    def test_prompt_limits_pubmed_abstracts(self):
        item = {"organism_name": "Novelibacter testii", "canonical_key": "novelibactertestii", "ncbi": {}}
        literature = {"articles": [{"pmid": str(index), "title": "title", "abstract": "x" * 5000} for index in range(8)]}
        payload = json.loads(prompt_for(item, literature)[1]["content"])
        self.assertEqual(len(payload["pubmed_search_hits"]), 5)
        self.assertEqual(len(payload["pubmed_search_hits"][0]["abstract"]), 2400)

    def test_api_error_is_sanitized(self):
        class FakeError(Exception):
            status_code = 401
            body = {"error": {"message": "Incorrect API key provided: sk-secret. You can inspect keys.", "code": "invalid_api_key", "param": None}}
            request_id = "req_test"

        result = _safe_api_error(FakeError(), "novelibactertestii", "test-model")
        self.assertNotIn("sk-secret", json.dumps(result))
        self.assertEqual(result["error_code"], "invalid_api_key")

    def test_codex_cli_backend_validates_structured_result(self):
        item = {"organism_name": "Novelibacter testii", "canonical_key": "novelibactertestii", "ncbi": {}}
        proposal = {
            "identity_assessment": "not_resolved",
            "family_hypothesis": "unmapped_or_uncertain",
            "clinical_traits_hypothesis": [],
            "human_infection_evidence": "No supplied evidence.",
            "specimen_relevance": "Unknown.",
            "contamination_or_colonization_caveat": "Unknown.",
            "verification_questions": ["Verify identity."],
            "literature_search_queries": ["Novelibacter testii infection"],
            "supporting_pmids": [],
            "confidence": "low",
        }

        def fake_run(command, **_kwargs):
            output_path = Path(command[command.index("--output-last-message") + 1])
            output_path.write_text(json.dumps(proposal), encoding="utf-8")
            return type("Completed", (), {"returncode": 0, "stdout": "", "stderr": ""})()

        with tempfile.TemporaryDirectory() as tmp, patch(
            "tools.draft_organism_taxonomy_reviews.shutil.which", return_value="codex.exe"
        ), patch("tools.draft_organism_taxonomy_reviews.subprocess.run", side_effect=fake_run):
            draft = make_draft_codex_cli(item, "gpt-test", None, Path(tmp))
        self.assertEqual(draft["backend"], "codex_cli_chatgpt_login")
        self.assertEqual(draft["proposal"]["family_hypothesis"], "unmapped_or_uncertain")

    def test_llm_draft_is_not_an_approval(self):
        row = {
            "organism_name": "Novelibacter testii",
            "canonical_key": "novelibactertestii",
            "decision": "llm_draft_unverified",
            "primary_rule_family": "environmental_low_specificity",
        }
        with self.assertRaises(ValueError):
            validate_approval(row)

    def test_proteus_style_review_preserves_nonpulmonary_caveat(self):
        row = {
            "organism_name": "Proteus mirabilis",
            "canonical_key": "proteusmirabilis",
            "decision": "approved",
            "evidence_verified": True,
            "manual_review": {
                "reviewer": "Reviewer A",
                "reviewed_at": "2026-09-24",
                "rationale": "Identity, pneumonia literature, and urinary-source caveat verified.",
            },
            "taxid": 584,
            "taxonomic_rank": "species",
            "biological_class": "bacterium",
            "primary_rule_family": "hospital_or_nonfermenter_gnb",
            "clinical_traits": ["pulmonary_pathogen_possible", "urinary_source_prone"],
            "source_references": [
                {"role": "identity", "url": "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=584"},
                {"role": "clinical", "url": "https://pubmed.ncbi.nlm.nih.gov/20502932/"},
            ],
        }
        with patch(
            "tools.approve_organism_taxonomy_reviews.classify_organism",
            return_value={"mapping_status": "unmapped"},
        ):
            key, profile = validate_approval(row)
        self.assertEqual("proteusmirabilis", key)
        self.assertEqual("hospital_or_nonfermenter_gnb", profile["primary_rule_family"])
        self.assertIn("urinary_source_prone", profile["clinical_traits"])

    def test_approval_needs_human_and_clinical_reference(self):
        row = {
            "organism_name": "Novelibacter testii",
            "canonical_key": "novelibactertestii",
            "decision": "approved",
            "evidence_verified": True,
            "manual_review": {"reviewer": "Reviewer A", "reviewed_at": "2026-09-22", "rationale": "Verified identity and clinical ecology."},
            "taxid": 123,
            "taxonomic_rank": "species",
            "biological_class": "bacterium",
            "primary_rule_family": "environmental_low_specificity",
            "clinical_traits": ["direct_support_preferred"],
            "source_references": [{"role": "identity", "url": "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=123"}],
        }
        with self.assertRaises(ValueError):
            validate_approval(row)
        row["source_references"].append({"role": "clinical", "url": "https://pubmed.ncbi.nlm.nih.gov/123/"})
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "approved.json"
            source.write_text(json.dumps([row]), encoding="utf-8")
            output = root / "staged.json"
            summary = compile_approvals(source, output)
            self.assertEqual(summary["mode"], "staged")
            profile = json.loads(output.read_text(encoding="utf-8"))["exact_profiles"]["novelibactertestii"]
            self.assertEqual(profile["primary_rule_family"], "environmental_low_specificity")


if __name__ == "__main__":
    unittest.main()

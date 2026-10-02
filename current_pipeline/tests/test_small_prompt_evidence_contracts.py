from __future__ import annotations

import unittest
from pathlib import Path


PROMPT_ROOT = Path(__file__).resolve().parents[1] / "agents" / "prompts"


class SmallPromptEvidenceContractTests(unittest.TestCase):
    def prompt_text(self, name: str) -> str:
        return (PROMPT_ROOT / name).read_text(encoding="utf-8")

    def assert_prompt_contains(self, name: str, required: tuple[str, ...]) -> None:
        text = self.prompt_text(name)
        for token in required:
            with self.subTest(prompt=name, token=token):
                self.assertIn(token, text)

    def test_host_prompt_preserves_medications_labs_and_missingness(self) -> None:
        self.assert_prompt_contains(
            "CBC_Lab_underlying_prompt.txt",
            (
                "AdmissionDiagnosis",
                '"structured_raw_evidence"',
                '"host_risk_medications"',
                '"cbc_data_status"',
                '"lab_observations"',
                '"source_text"',
                '"unit"',
                "alc_potentially_derivable",
            ),
        )

    def test_culture_prompt_preserves_each_record(self) -> None:
        self.assert_prompt_contains(
            "culture_agent_prompt.txt",
            (
                '"culture_observations"',
                '"observation_id"',
                '"evidence_observation_ids"',
                '"raw_quantity"',
                '"quantity_unit"',
                '"collected_time"',
                '"reported_time"',
                '"source_text"',
            ),
        )

    def test_assay_prompt_preserves_numeric_and_negative_records(self) -> None:
        self.assert_prompt_contains(
            "filmarray_GMtest_prompt.txt",
            (
                '"assay_observations"',
                '"evidence_observation_ids"',
                '"numeric_value"',
                '"reference_or_cutoff"',
                '"semiquant_bin"',
                "equivocal",
                "invalid",
            ),
        )

    def test_image_prompt_preserves_studies_and_source_sentences(self) -> None:
        self.assert_prompt_contains(
            "image_agent_prompt.txt",
            (
                '"imaging_studies"',
                '"study_id"',
                '"raw_findings"',
                '"raw_impression"',
                '"source_sentence"',
                '"diagnostic_mentions"',
                '"microbiologic_confirmation": false',
            ),
        )

    def test_molecular_prompt_preserves_quantitative_and_negative_records(self) -> None:
        self.assert_prompt_contains(
            "molecular_microbiology_prompt.txt",
            (
                '"molecular_observations"',
                '"evidence_observation_ids"',
                '"ct_value"',
                '"viral_load_or_copies"',
                '"reference_or_cutoff"',
                "not_detected",
                "indeterminate",
            ),
        )


if __name__ == "__main__":
    unittest.main()

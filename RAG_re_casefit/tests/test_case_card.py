import copy
import unittest

from rag_re_casefit.case_card import build_case_card, prompt_case_card, validate_case_card
from rag_re_casefit.cli import _articles


class CaseCardTests(unittest.TestCase):
    def fixture(self):
        return {
            "patient_id": "14",
            "gold": ["Candida tropicalis"],
            "deterministic_max": {
                "dominant_source": "Respiratory",
                "host_context": {"host_vulnerability_tier": "V3"},
            },
        }, {
            "organism_name": "Candida tropicalis",
            "specimen_class": "S2_lower_respiratory",
            "reads": 999,
            "key_evidence": {"specimen_type": "Lower BAL"},
            "integrated_reasoning": "gold-like upstream conclusion",
        }

    def test_sparse_current_data_builds_without_guessing_host(self):
        root, candidate = self.fixture()
        card = build_case_card(root, candidate)
        self.assertEqual(card["target_syndrome"], "lower_respiratory_infection")
        self.assertEqual(card["clinical_site"], "lower_respiratory")
        self.assertEqual(card["specimen_type"], "bal_or_balf")
        self.assertEqual(card["host_factors"], [])
        self.assertIn("raw_host_factors", card["missing_fields"])

    def test_prompt_projection_excludes_identifiers_gold_and_upstream_signals(self):
        root, candidate = self.fixture()
        card = build_case_card(root, candidate)
        prompt = prompt_case_card(card)
        blob = repr(prompt).casefold()
        self.assertNotIn("patient_id", prompt)
        self.assertNotIn("14", blob)
        self.assertNotIn("gold", blob)
        self.assertNotIn("reads", blob)
        self.assertNotIn("v3", blob)
        self.assertNotIn("candida", blob)

    def test_extra_key_is_rejected(self):
        root, candidate = self.fixture()
        card = build_case_card(root, candidate)
        bad = copy.deepcopy(card)
        bad["answer"] = "Candida tropicalis"
        with self.assertRaises(ValueError):
            validate_case_card(bad)

    def test_can_reuse_frozen_articles_from_rag_output(self):
        value = {
            "candidates": [{
                "organism_name": "CMV",
                "canonical_organism_name": "Human cytomegalovirus",
                "literature_evidence": {
                    "retrieval": {"articles": [{"pmid": "1", "abstract": "x"}]}
                },
            }]
        }
        self.assertEqual(_articles(value, "Human cytomegalovirus")[0]["pmid"], "1")


if __name__ == "__main__":
    unittest.main()

import unittest

from rag_re_casefit.rules import aggregate_judgments, apply_traceability_gate, derive_article_verdict
from rag_re_casefit.evaluation import canonical


def extraction(**overrides):
    value = {
        "organism_scope": "exact_species",
        "human_clinical_context": "causal_infection",
        "site_alignment": "match",
        "syndrome_alignment": "match",
        "specimen_alignment": "unknown",
        "host_alignment": "unknown",
        "phenotype_alignment": "unknown",
        "temporal_alignment": "unknown",
    }
    value.update(overrides)
    return value


class RuleTests(unittest.TestCase):
    def test_covid_and_sars_cov_2_are_frozen_synonyms(self):
        self.assertEqual(canonical("COVID-19"), canonical("SARS-CoV-2"))

    def test_strong_requires_an_extra_patient_specific_match(self):
        result = derive_article_verdict(extraction(specimen_alignment="match"))
        self.assertEqual(result["verdict"], "strong_match")

    def test_hard_gates_without_extra_detail_are_partial(self):
        result = derive_article_verdict(extraction())
        self.assertEqual(result["verdict"], "partial_match")

    def test_unknown_hard_gate_is_insufficient_not_negative(self):
        result = derive_article_verdict(extraction(site_alignment="unknown"))
        self.assertEqual(result["verdict"], "insufficient")
        self.assertFalse(result["counterevidence"])

    def test_colonization_at_same_site_is_counterevidence(self):
        result = derive_article_verdict(extraction(human_clinical_context="colonization_or_contamination"))
        self.assertEqual(result["verdict"], "mismatch")
        self.assertTrue(result["counterevidence"])
        self.assertFalse(result["retrieval_mismatch"])

    def test_detection_only_is_neutral_non_support_not_counterevidence(self):
        result = derive_article_verdict(extraction(human_clinical_context="detection_only"))
        self.assertEqual(result["verdict"], "mismatch")
        self.assertFalse(result["counterevidence"])
        self.assertFalse(result["retrieval_mismatch"])

    def test_unrelated_partial_site_reactivation_is_not_counterevidence(self):
        result = derive_article_verdict(extraction(
            human_clinical_context="reactivation_or_bystander",
            site_alignment="partial",
            syndrome_alignment="partial",
        ))
        self.assertEqual(result["verdict"], "mismatch")
        self.assertFalse(result["counterevidence"])

    def test_support_without_valid_span_becomes_insufficient(self):
        decision = derive_article_verdict(extraction(specimen_alignment="match"))
        result = apply_traceability_gate(
            decision, support_span_valid=False, counter_span_valid=False
        )
        self.assertEqual(result["verdict"], "insufficient")
        self.assertEqual(result["traceability_status"], "invalid_required_span")

    def test_counter_without_valid_span_becomes_neutral(self):
        decision = derive_article_verdict(extraction(
            human_clinical_context="colonization_or_contamination"
        ))
        result = apply_traceability_gate(
            decision, support_span_valid=True, counter_span_valid=False
        )
        self.assertFalse(result["counterevidence"])

    def test_wrong_species_is_retrieval_mismatch_not_counterevidence(self):
        result = derive_article_verdict(extraction(organism_scope="genus_only"))
        self.assertEqual(result["verdict"], "mismatch")
        self.assertFalse(result["counterevidence"])
        self.assertTrue(result["retrieval_mismatch"])

    def test_aggregate_high_never_auto_rescues(self):
        rows = [
            {"pmid": str(index), "verdict": "strong_match" if index < 2 else "partial_match", "counterevidence": False, "retrieval_mismatch": False}
            for index in range(5)
        ]
        result = aggregate_judgments(rows)
        self.assertEqual(result["casefit_evidence_tier"], "high")
        self.assertFalse(result["auto_rescue"])

    def test_support_plus_counterevidence_is_conflicted(self):
        rows = [
            {"pmid": str(index), "verdict": "partial_match", "counterevidence": False, "retrieval_mismatch": False}
            for index in range(4)
        ]
        rows.append({"pmid": "5", "verdict": "mismatch", "counterevidence": True, "retrieval_mismatch": False})
        result = aggregate_judgments(rows)
        self.assertEqual(result["casefit_evidence_tier"], "conflicted")

    def test_duplicate_pmids_are_rejected(self):
        rows = [
            {"pmid": "1", "verdict": "partial_match", "counterevidence": False, "retrieval_mismatch": False},
            {"pmid": "1", "verdict": "partial_match", "counterevidence": False, "retrieval_mismatch": False},
        ]
        with self.assertRaises(ValueError):
            aggregate_judgments(rows)


if __name__ == "__main__":
    unittest.main()

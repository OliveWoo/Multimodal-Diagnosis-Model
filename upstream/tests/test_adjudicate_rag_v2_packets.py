from __future__ import annotations

import unittest

from tools import adjudicate_rag_v2_packets as adjudicate


def packet() -> dict:
    return {
        "case_id": "P30_Bacteroides_fragilis",
        "literature_evidence": {
            "knowledge_card_key": "adult_pulmonary:bacteroidesfragilis",
            "articles": [{"pmid": "123", "doi": "10.1/example"}],
        },
    }


def result() -> dict:
    return {
        "schema_version": "rag_adjudication_v2.0",
        "case_id": "P30_Bacteroides_fragilis",
        "knowledge_card": {
            "knowledge_card_key": "adult_pulmonary:bacteroidesfragilis",
            "pulmonary_pathogenicity": "recognized_but_uncommon",
            "respiratory_colonization_risk": "moderate",
            "typical_supporting_evidence": ["lung abscess"],
        },
        "retrieval_summary": {
            "search_completed": True,
            "sources": [{
                "citation_id": "PMID:123",
                "title": "Example",
                "year": 2020,
                "pmid_or_doi": "123",
                "study_type": "case_series",
                "quality": "moderate",
                "case_similarity": "moderate",
                "direction": "supports_pathogenicity",
                "key_point": "Example",
            }],
            "weighted_evidence": {
                "supports_pathogenicity": 1,
                "supports_colonization_or_background": 0,
                "uncertainty": 1,
            },
        },
        "patient_fit": {
            "supporting_local_evidence": ["lower respiratory"],
            "contradicting_local_evidence": ["D0"],
            "missing_key_evidence": ["culture"],
            "competing_explanation_strength": "moderate",
        },
        "decision": {
            "recommended_action": "context_only",
            "clinician_visibility": "show_context",
            "confidence": "low",
            "one_sentence_reason": "Evidence is incomplete.",
        },
        "audit": {"human_review_required": True, "limitations": ["Sparse literature"]},
    }


class AdjudicateRagV2Tests(unittest.TestCase):
    def test_valid_result_passes(self) -> None:
        adjudicate.validate_adjudication(result(), packet())

    def test_invented_citation_fails(self) -> None:
        value = result()
        value["retrieval_summary"]["sources"][0]["pmid_or_doi"] = "999"
        with self.assertRaises(ValueError):
            adjudicate.validate_adjudication(value, packet())

    def test_prefixed_composite_citation_passes(self) -> None:
        value = result()
        value["retrieval_summary"]["sources"][0]["pmid_or_doi"] = "PMID:123; DOI:10.1/example"
        adjudicate.validate_adjudication(value, packet())

    def test_case_id_mismatch_fails(self) -> None:
        value = result()
        value["case_id"] = "P99_X"
        with self.assertRaises(ValueError):
            adjudicate.validate_adjudication(value, packet())


if __name__ == "__main__":
    unittest.main()

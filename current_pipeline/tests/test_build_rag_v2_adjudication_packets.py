from __future__ import annotations

import unittest

from tools import build_rag_v2_adjudication_packets as packets


class RagV2AdjudicationPacketTests(unittest.TestCase):
    def test_packet_keeps_patient_and_literature_evidence_separate(self) -> None:
        case = {
            "schema_version": "rag_production_case_v2.0",
            "case_id": "P14_Candida_tropicalis",
            "dataset_label": "test",
            "patient_id": "14",
            "organism": {"display_name": "Candida tropicalis"},
            "current_model_state": {"source_tier": "review_high_priority"},
            "local_evidence": {"mngs_evidence": {"reads": 100}},
            "rag_task": {
                "knowledge_card_key": "adult_pulmonary:candidatropicalis",
                "fixed_questions": ["Is it pathogenic?"],
                "required_output_schema": "rag_adjudication_v2.0",
            },
            "provenance": {"merged_suffix": "x"},
        }
        retrieval = {
            "schema_version": "rag_pubmed_retrieval_v2.1",
            "knowledge_card_key": "adult_pulmonary:candidatropicalis",
            "retrieved_at": "2026-08-13T00:00:00Z",
            "retrieval_evidence_volume": "limited",
            "selected_article_count": 1,
            "selected_track_counts": {"pulmonary_pathogenicity": 1},
            "missing_evidence_tracks": [],
            "retrieval_errors": [],
            "articles": [{"pmid": "1", "evidence_direction": "unassessed"}],
            "adjudication_note": "not a conclusion",
        }
        packet = packets.make_packet(case, retrieval)
        self.assertEqual(packet["local_evidence"]["mngs_evidence"]["reads"], 100)
        self.assertEqual(packet["literature_evidence"]["articles"][0]["pmid"], "1")
        self.assertFalse(packet["provenance"]["answer_source_read"])

    def test_forbidden_answer_field_is_rejected(self) -> None:
        case = {
            "case_id": "P1_X",
            "organism": {"display_name": "X"},
            "local_evidence": {"answer_hit": True},
            "rag_task": {"knowledge_card_key": "k"},
        }
        retrieval = {"knowledge_card_key": "k", "articles": []}
        with self.assertRaises(ValueError):
            packets.make_packet(case, retrieval)


if __name__ == "__main__":
    unittest.main()

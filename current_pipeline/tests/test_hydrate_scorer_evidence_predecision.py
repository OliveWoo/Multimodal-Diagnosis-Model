import copy
import json
import tempfile
import unittest
from pathlib import Path

from tools.hydrate_scorer_evidence_predecision import (
    decision_signature,
    graph_edges,
    sync_history_consumption_ledger,
)


class PredecisionEvidenceGraphTests(unittest.TestCase):
    def test_graph_contains_mngs_and_hospital_source_edges(self):
        candidate = {
            "organism_name": "Candida albicans",
            "analytical_profile": {
                "per_test_signals": [{
                    "seq_id": "DNA01",
                    "reads": 42,
                    "rpm_total": 7.5,
                    "rank_in_retained_universe": 2,
                    "qc_status": "complete",
                }]
            },
            "hospital_profile": {
                "hospital_evidence_detail": {
                    "source_observations": [{
                        "observation_id": "CUL-001",
                        "source_record_key": "row-1",
                        "specimen_type": "BALF",
                        "collected_time": "2026-01-01 10:00",
                        "raw_result": "positive",
                    }]
                }
            },
        }
        edges = graph_edges(candidate)
        self.assertEqual(2, len(edges))
        self.assertEqual(
            {"mngs_per_test_signal", "hospital_source_observation"},
            {edge["evidence_type"] for edge in edges},
        )

    def test_provenance_fields_do_not_change_decision_signature(self):
        payload = {
            "clinical_scorer_forward": [{
                "organism_name": "Pneumocystis jirovecii",
                "decision": "review_high_priority",
                "analytic_route": "Priority",
                "formal_pick_allowed": False,
            }]
        }
        hydrated = copy.deepcopy(payload)
        hydrated["clinical_scorer_forward"][0]["hospital_evidence_detail"] = {
            "source_observations": []
        }
        self.assertEqual(decision_signature(payload), decision_signature(hydrated))

    def test_decision_change_is_detectable(self):
        payload = {
            "clinical_scorer_forward": [{
                "organism_name": "Example organism",
                "decision": "review_context_needed",
            }]
        }
        changed = copy.deepcopy(payload)
        changed["clinical_scorer_forward"][0]["decision"] = "picked_shadow"
        self.assertNotEqual(decision_signature(payload), decision_signature(changed))

    def test_history_consumption_is_written_to_shared_graph_and_ledger(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            graph_root = root / "graph"
            graph_dir = graph_root / "evidence_graphs"
            graph_dir.mkdir(parents=True)
            history_dir = root / "history" / "patient_outputs"
            history_dir.mkdir(parents=True)
            graph_path = graph_dir / "P3_candidate_evidence_edges.json"
            graph_path.write_text(
                json.dumps({"patient_id": 3, "edges": []}), encoding="utf-8"
            )
            history_path = (
                history_dir / "NGS_patient_3_history_adjusted_route_shadow.json"
            )
            history_path.write_text(json.dumps({
                "candidates": [{
                    "organism_name": "Candida albicans",
                    "history_policy_family": "reactivation",
                    "history_rule_ids": ["H-1"],
                    "history_evidence_ids": ["P3:HISTCLAIM:abc"],
                    "consumed_evidence_ids": ["P3:HISTCLAIM:abc"],
                    "history_evidence_consumption": [{
                        "evidence_id": "P3:HISTCLAIM:abc",
                        "downstream_reuse": "provenance_only_not_independent_support",
                    }],
                }]
            }), encoding="utf-8")
            (graph_root / "summary.json").write_text(
                json.dumps({"patient_count": 1}), encoding="utf-8"
            )

            result = sync_history_consumption_ledger(
                graph_root, root / "history"
            )

            graph = json.loads(graph_path.read_text(encoding="utf-8"))
            self.assertEqual(1, result["ledger_entry_count"])
            self.assertEqual("history_claim", graph["edges"][0]["evidence_type"])
            self.assertEqual(
                "history_adjusted_route",
                graph["evidence_consumption_ledger"][0]["consumed_at_stage"],
            )
            self.assertTrue((graph_root / "evidence_consumption_ledger.csv").is_file())


if __name__ == "__main__":
    unittest.main()

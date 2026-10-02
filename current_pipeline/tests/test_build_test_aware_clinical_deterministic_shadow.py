import json
import unittest
from pathlib import Path

from tools.build_test_aware_clinical_deterministic_shadow import (
    DEFAULT_POLICY,
    build_candidate,
    decide,
    generic_identity,
)
from tools.organism_taxonomy_classifier import classify_organism


class TestClinicalDeterministicShadow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.policy = json.loads(Path(DEFAULT_POLICY).read_text(encoding="utf-8"))

    def test_pjp_exact_image_differential_promotes_to_high_not_picked(self):
        candidate = {
            "patient_id": "14", "organism_name": "Pneumocystis jirovecii",
            "category": "fungal", "taxonomy_family": "high_consequence_opportunistic",
            "decision": "review_context_needed", "integrated_level": "Level 4",
            "case_review_id": "CASE-A", "host_support": False,
            "analytical_profile": {"best_rank_in_retained_universe": 4,
                                   "reproducibility_axis": True},
        }
        timeline = {
            "patient_id": "14",
            "index_events": [{"event_id": "P14:mngs:1", "case_review_id": "CASE-A"}],
            "clinical_evidence": [{
                "evidence_id": "P14:E1", "evidence_role": "imaging_context",
                "source_text": "Differential diagnosis: PJP.",
                "diagnostic_mentions": [{"concept_name": "Pneumocystis jirovecii",
                                         "assertion": "differential"}],
                "episode_evidence_status": {"P14:mngs:1": "episode_observation"},
                "available_at_index": {"P14:mngs:1": "post_index_not_available"},
            }],
        }
        phenotype = {"phenotype_status": "phenotype_file_missing", "events": [{
            "index_event": timeline["index_events"][0], "candidate_packets": [{
                "canonical_key": "pneumocystisjirovecii", "phenotype_context": [],
                "exact_organism_workbook_microbiology": [],
            }],
        }]}
        result = build_candidate(candidate, timeline, phenotype, self.policy)
        self.assertEqual(result["clinical_decision"], "review_high_priority")
        self.assertFalse(result["formal_pick_allowed"])
        self.assertEqual(len(result["image_evidence"]), 1)
        self.assertFalse(result["image_evidence"][0]["microbiologic_confirmation"])

    def test_weak_image_signal_stays_context(self):
        taxonomy = classify_organism("Pneumocystis jirovecii")
        candidate = {"organism_name": "Pneumocystis jirovecii",
                     "taxonomy_family": "high_consequence_opportunistic",
                     "decision": "review_context_needed",
                     "analytical_profile": {"best_rank_in_retained_universe": 20,
                                            "reproducibility_axis": False}}
        decision, _, _ = decide(candidate, taxonomy, [{"assertion": "differential"}],
                                {"routed_phenotype_context": []}, self.policy)
        self.assertEqual(decision, "review_context_needed")

    def test_generic_hospital_pick_is_demoted_to_high(self):
        name = "Acinetobacter calcoaceticus-baumannii complex"
        taxonomy = classify_organism(name)
        candidate = {"organism_name": name, "evidence_source": "hospital_only",
                     "decision": "picked_shadow"}
        self.assertTrue(generic_identity(name, taxonomy, self.policy))
        decision, _, _ = decide(candidate, taxonomy, [],
                                {"routed_phenotype_context": []}, self.policy)
        self.assertEqual(decision, "review_high_priority")

    def test_phenotype_alone_does_not_change_context_tier(self):
        taxonomy = classify_organism("Candida tropicalis")
        candidate = {"organism_name": "Candida tropicalis",
                     "taxonomy_family": "candida_or_yeast",
                     "decision": "review_context_needed",
                     "analytical_profile": {"best_rank_in_retained_universe": 1,
                                            "reproducibility_axis": True}}
        decision, rules, _ = decide(candidate, taxonomy, [],
                                    {"routed_phenotype_context": [{"phenotype": "ASPIRATION"}]},
                                    self.policy)
        self.assertEqual(decision, "review_context_needed")
        self.assertIn("CLIN-PHENO1-PROVISIONAL-CONTEXT-NO-TIER-EFFECT", rules)

    def test_exact_existing_pick_is_preserved(self):
        taxonomy = classify_organism("Pseudomonas aeruginosa")
        candidate = {"organism_name": "Pseudomonas aeruginosa",
                     "decision": "picked_shadow"}
        decision, _, _ = decide(candidate, taxonomy, [],
                                {"routed_phenotype_context": []}, self.policy)
        self.assertEqual(decision, "picked_shadow")

    def test_history_priority_route_enters_high_without_direct_pick(self):
        candidate = {
            "patient_id": "30",
            "organism_name": "Bacteroides fragilis",
            "category": "bacterial",
            "taxonomy_family": "oral_aspiration_or_anaerobe",
            "decision": "review_context_needed",
            "history_routed_decision": "review_high_priority",
            "analytic_route": "Context",
            "history_adjusted_route": "Priority",
            "history_route_changed": True,
            "history_policy_family": "aspiration",
            "history_rule_ids": ["HIST-ASP3-CONTEXT-TO-PRIORITY"],
            "history_evidence_ids": ["P30:HISTCLAIM:1"],
            "consumed_evidence_ids": ["P30:HISTCLAIM:1"],
            "history_evidence_consumption": [{
                "evidence_id": "P30:HISTCLAIM:1",
                "consumed_at_stage": "history_adjusted_route",
                "downstream_reuse": "provenance_only_not_independent_support",
            }],
            "future_picked_gate_review_required": True,
            "future_picked_gate_review_status": "pending_independent_evidence_audit",
            "future_picked_gate_review_requirements": ["independent disease evidence"],
            "analytical_profile": {
                "best_rank_in_retained_universe": 3,
                "selected_positive_test_count": 3,
                "reproducibility_axis": True,
                "cross_molecule_selected": True,
            },
        }
        timeline = {"patient_id": "30", "index_events": [], "clinical_evidence": []}
        phenotype = {"phenotype_status": "available", "events": []}
        result = build_candidate(candidate, timeline, phenotype, self.policy)
        self.assertEqual(result["incoming_decision"], "review_context_needed")
        self.assertEqual(result["history_routed_decision"], "review_high_priority")
        self.assertEqual(result["clinical_decision"], "review_high_priority")
        self.assertFalse(result["formal_pick_allowed"])
        self.assertEqual(result["consumed_evidence_ids"], ["P30:HISTCLAIM:1"])
        self.assertEqual(result["history_policy_family"], "aspiration")
        self.assertTrue(result["future_picked_gate_review_required"])
        self.assertEqual(
            result["future_picked_gate_review_status"],
            "pending_independent_evidence_audit",
        )
        self.assertIn("CLIN-HIST1-HISTORY-ROUTE-TO-SCORER", result["rule_ids"])


if __name__ == "__main__":
    unittest.main()

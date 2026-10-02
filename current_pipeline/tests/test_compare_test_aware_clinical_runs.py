from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.compare_test_aware_clinical_runs import compare_runs


class CompareTestAwareClinicalRunsTests(unittest.TestCase):
    def write_run(self, root: Path, candidates: list[dict], *, history: bool) -> None:
        clinical = root / "clinical"
        promotion = root / "promotion"
        (promotion / "patient_outputs").mkdir(parents=True)
        clinical.mkdir(parents=True)
        (clinical / "summary.json").write_text(
            json.dumps(
                {
                    "policy_sha256": "clinical-sha",
                    "source_scorer_root": "frozen-scorer",
                    "timeline_root": "frozen-timeline",
                    "phenotype_packet_root": "frozen-phenotype",
                    "history_route_root": "history" if history else None,
                }
            ),
            encoding="utf-8",
        )
        (promotion / "summary.json").write_text(
            json.dumps(
                {
                    "policy_sha256": "promotion-sha",
                    "source_clinical_shadow_root": str(clinical.resolve()),
                }
            ),
            encoding="utf-8",
        )
        (promotion / "patient_outputs" / "patient_1.json").write_text(
            json.dumps(
                {
                    "patient_id": "1",
                    "all_forwarded_candidates": candidates,
                }
            ),
            encoding="utf-8",
        )

    def test_reports_isolated_history_changes(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            control_candidate = {
                "organism_name": "Bacteroides fragilis",
                "organism_key": "bacteroidesfragilis",
                "clinical_decision": "review_context_needed",
            }
            history_candidate = {
                **control_candidate,
                "clinical_decision": "review_high_priority",
                "analytic_route": "Context",
                "history_adjusted_route": "Priority",
                "history_route_changed": True,
                "history_rule_ids": ["HIST-ASP3-CONTEXT-TO-PRIORITY"],
                "history_evidence_ids": ["H1"],
                "consumed_evidence_ids": ["H1"],
                "promotion_gate": {
                    "evaluated": True,
                    "outcome": "retained_high_priority",
                    "route": "no_promotion",
                    "blockers": ["no family-specific route"],
                    "axes": {"host_support_eligible_for_promotion": False},
                },
            }
            added_candidate = {
                "organism_name": "Dialister invisus",
                "organism_key": "dialisterinvisus",
                "clinical_decision": "review_context_needed",
                "analytic_route": "Audit",
                "history_adjusted_route": "Context",
                "history_route_changed": True,
                "history_evidence_ids": ["H2"],
                "consumed_evidence_ids": ["H2"],
                "promotion_gate": {"evaluated": False, "outcome": "not_applicable"},
            }
            self.write_run(root / "control", [control_candidate], history=False)
            self.write_run(
                root / "candidate",
                [history_candidate, added_candidate],
                history=True,
            )

            report = compare_runs(
                root / "control" / "promotion",
                root / "candidate" / "promotion",
            )

            self.assertTrue(report["policy_and_input_checks"]["all_equal"])
            self.assertEqual(report["summary"]["candidate_added_count"], 1)
            self.assertEqual(report["summary"]["existing_candidate_tier_change_count"], 1)
            self.assertEqual(report["summary"]["history_route_changed_count"], 2)
            self.assertEqual(report["summary"]["incremental_route_effect_change_count"], 2)
            self.assertEqual(report["summary"]["history_changed_ending_picked_count"], 0)
            self.assertEqual(report["summary"]["double_count_guard_violation_count"], 0)
            self.assertEqual(len(report["changes"]), 2)

    def test_detects_consumed_history_reuse(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            candidate = {
                "organism_name": "Example organism",
                "organism_key": "exampleorganism",
                "clinical_decision": "review_high_priority",
                "history_route_changed": True,
                "consumed_evidence_ids": ["H1"],
                "promotion_gate": {
                    "evaluated": True,
                    "axes": {"host_support_eligible_for_promotion": True},
                },
            }
            self.write_run(root / "control", [], history=False)
            self.write_run(root / "candidate", [candidate], history=True)

            report = compare_runs(
                root / "control" / "promotion",
                root / "candidate" / "promotion",
            )

            self.assertEqual(report["summary"]["double_count_guard_violation_count"], 1)

    def test_omits_unchanged_history_effect_between_two_history_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            shared = {
                "organism_name": "Shared organism",
                "organism_key": "sharedorganism",
                "clinical_decision": "review_context_needed",
                "analytic_route": "Audit",
                "history_adjusted_route": "Context",
                "history_route_changed": True,
                "history_policy_family": "aspiration",
                "history_rule_ids": ["HIST-ASP2-AUDIT-TO-CONTEXT"],
                "consumed_evidence_ids": ["H1"],
            }
            added_effect = {
                "organism_name": "New organism",
                "organism_key": "neworganism",
                "clinical_decision": "review_high_priority",
                "analytic_route": "Context",
                "history_adjusted_route": "Priority",
                "history_route_changed": True,
                "history_policy_family": "opportunistic_host_risk",
                "history_rule_ids": ["HIST-OPP3-CONTEXT-TO-PRIORITY"],
                "consumed_evidence_ids": ["H2"],
            }
            before_new = {
                **added_effect,
                "clinical_decision": "review_high_priority",
                "history_adjusted_route": "Context",
                "history_route_changed": False,
                "history_rule_ids": [],
                "consumed_evidence_ids": [],
            }
            self.write_run(root / "control", [shared, before_new], history=True)
            self.write_run(root / "candidate", [shared, added_effect], history=True)

            report = compare_runs(
                root / "control" / "promotion",
                root / "candidate" / "promotion",
            )

            self.assertEqual(report["summary"]["history_route_changed_count"], 2)
            self.assertEqual(report["summary"]["incremental_route_effect_change_count"], 1)
            self.assertEqual(len(report["changes"]), 1)
            self.assertEqual(report["changes"][0]["organism_name"], "New organism")


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import copy
import unittest
from unittest.mock import patch

from rag_re_clinical.engine import ARM_DEFINITIONS, run_patient


def config_fixture():
    return {
        "schema_version": "rag_re_clinical.config.v1",
        "module_a": {"strong_min_support": 3, "strong_min_judgeable": 5},
        "decision": {
            "primary_arm": "CL2_CONVERGENT",
            "baseline_immutable": True,
        },
    }


def frozen_fixture():
    def row(name, baseline, tier=None):
        return {
            "organism_name": name,
            "canonical_organism_name": name,
            "baseline_selected": baseline,
            "review_tier": tier,
            "candidate_provenance": "fixture",
            "source_context": {
                "site_code": "lower_respiratory",
                "specimen_site_code": "lower_respiratory",
                "specimen_alignment": "Aligned",
            },
            "rule_input": {
                "organism_name": name,
                "specimen_alignment": "Aligned",
                "specimen_class": "S2_lower_respiratory",
            },
            "modules": {
                "A": {
                    "positive": False,
                    "support_count": 0,
                    "judgeable_count": 10,
                }
            },
        }

    return {
        "schema_version": "rag_re.output.v1",
        "run_status": "complete",
        "patient_id": "7",
        "pipeline_fingerprint": "fixture",
        "input_meta": {"candidate_pool_sha256": "fixture"},
        "candidates": [
            row("Baseline species", True),
            row("Rescue species", False, "review_high_priority"),
        ],
    }


def merge_fixture():
    return {
        "patient_id": "7",
        "deterministic_max": {
            "pathogen_candidates": [
                {"organism_name": "Baseline species"},
                {"organism_name": "Rescue species"},
            ],
            "best_available_summary": {
                "picked_pathogens": [{"organism_name": "Baseline species"}]
            },
        },
        "llm_missed_candidate_review": {
            "review_high_priority": [
                {"organism_name": "Rescue species", "evidence_snapshot": {}}
            ],
            "review_context_needed": [],
        },
    }


def module_result(e_state):
    return {
        "D": {"state": "D_PASS", "status": "ok"},
        "B_STRICT": {"positive": True, "status": "partial"},
        "C": {"grade": "C3" if e_state == "E3_DIRECT_HIGH" else "C0", "status": "ok"},
        "E": {"state": e_state, "status": "ok"},
        "F": {"state": "F_UNKNOWN", "status": "insufficient_input"},
    }


class EngineTests(unittest.TestCase):
    def test_every_arm_retains_e0_even_when_its_clinical_modules_block(self):
        frozen = frozen_fixture()

        def modules(_merged, candidate, _config):
            if candidate["baseline_selected"]:
                return {
                    "D": {"state": "D_BLOCK", "status": "ok"},
                    "B_STRICT": {"positive": False, "status": "ok"},
                    "C": {"grade": "CNEG", "status": "ok"},
                    "E": {"state": "E_BLOCKED", "status": "ok"},
                    "F": {"state": "F_UNKNOWN", "status": "insufficient_input"},
                }
            return module_result("E_NO_CONVERGENCE")

        with patch("rag_re_clinical.engine.evaluate_modules", side_effect=modules):
            result = run_patient(frozen, merge_fixture(), config_fixture())
        for arm in ARM_DEFINITIONS:
            self.assertIn("Baseline species", result["arms"][arm]["predicted_pathogens"])
            self.assertTrue(result["candidates"][0]["arms"][arm])
        self.assertEqual(
            result["arms"]["CL0_BASELINE"]["predicted_pathogens"],
            ["Baseline species"],
        )

    def test_c3_rescue_is_added_without_changing_baseline(self):
        frozen = frozen_fixture()

        def modules(_merged, candidate, _config):
            return module_result(
                "E_NO_CONVERGENCE"
                if candidate["baseline_selected"]
                else "E3_DIRECT_HIGH"
            )

        with patch("rag_re_clinical.engine.evaluate_modules", side_effect=modules):
            result = run_patient(frozen, merge_fixture(), config_fixture())
        self.assertEqual(result["arms"]["CL0_BASELINE"]["predicted_count"], 1)
        for arm in (
            "CL1_DIRECT",
            "CL2_CONVERGENT",
            "CL3_BALANCED",
            "CL4_TIERED_EXPLORATORY",
        ):
            self.assertEqual(
                result["arms"][arm]["predicted_pathogens"],
                ["Baseline species", "Rescue species"],
            )
        self.assertEqual(result["candidates"][1]["disposition"], "AUTO_RESCUE")

    def test_disposition_uses_the_configured_primary_arm(self):
        frozen = frozen_fixture()
        cfg = config_fixture()
        cfg["decision"]["primary_arm"] = "CL3_BALANCED"

        def modules(_merged, candidate, _config):
            return module_result(
                "E_NO_CONVERGENCE"
                if candidate["baseline_selected"]
                else "E1_BIOLOGIC_PLUS_SIGNAL"
            )

        with patch("rag_re_clinical.engine.evaluate_modules", side_effect=modules):
            result = run_patient(frozen, merge_fixture(), cfg)
        rescue = result["candidates"][1]
        self.assertFalse(rescue["arms"]["CL2_CONVERGENT"])
        self.assertTrue(rescue["arms"]["CL3_BALANCED"])
        self.assertEqual(rescue["disposition"], "AUTO_RESCUE")

    def test_injected_answer_fields_do_not_change_modules_or_predictions(self):
        frozen = frozen_fixture()
        merge = merge_fixture()
        clean = run_patient(frozen, merge, config_fixture())

        tainted_frozen = copy.deepcopy(frozen)
        tainted_merge = copy.deepcopy(merge)
        tainted_frozen["gold"] = {"pathogens": ["Rescue species"]}
        tainted_frozen["candidates"][1]["answer"] = True
        tainted_merge["clinical_final_diagnosis"] = ["Rescue species"]
        tainted_merge["deterministic_max"]["pathogen_candidates"][1]["gold"] = True
        tainted = run_patient(tainted_frozen, tainted_merge, config_fixture())

        self.assertEqual(
            [row["modules"] for row in clean["candidates"]],
            [row["modules"] for row in tainted["candidates"]],
        )
        self.assertEqual(clean["arms"], tainted["arms"])

    def test_label_blind_patient_context_reaches_clinical_modules(self):
        frozen = frozen_fixture()
        merge = merge_fixture()
        merge["deterministic_max"]["host_context"] = {
            "host_vulnerability_tier": "V2"
        }
        observed = []

        def modules(merged, candidate, _config):
            observed.append(merged.get("patient_context"))
            return module_result("E_NO_CONVERGENCE")

        with patch("rag_re_clinical.engine.evaluate_modules", side_effect=modules):
            run_patient(frozen, merge, config_fixture())
        self.assertTrue(observed)
        self.assertTrue(
            all(
                context["deterministic_patient_context"]["host_context"]
                ["host_vulnerability_tier"]
                == "V2"
                for context in observed
            )
        )


if __name__ == "__main__":
    unittest.main()

import json
import tempfile
import unittest
from pathlib import Path

from tools.apply_test_aware_history_route_shadow import (
    apply_aspiration_cluster_guardrail,
    apply_hap_mdr_cluster_guardrail,
    apply_opportunistic_cluster_guardrail,
    adjust_candidate,
    aspiration_history_profile,
    colonization_history_profile,
    hap_mdr_history_profile,
    opportunistic_history_profile,
    phenotype_disagreement,
    reactivation_history_profile,
    read_json,
)


POLICY = read_json(Path("rules/test_aware_history_route_v1.json"))
POLICY_V2 = read_json(Path("rules/test_aware_history_route_v2.json"))
POLICY_V3 = read_json(Path("rules/test_aware_history_route_v3.json"))
POLICY_V4 = read_json(Path("rules/test_aware_history_route_v4.json"))
POLICY_V4B = read_json(Path("rules/test_aware_history_route_v4b.json"))
POLICY_V4C = read_json(Path("rules/test_aware_history_route_v4c_combined.json"))


def history_records():
    return [
        {
            "evidence_id": "P30:HIST:admission:2",
            "text": "Acute respiratory failure due to aspiration pneumonia",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[2].diagnosis"},
        },
        {
            "evidence_id": "P30:HIST:underlying:3",
            "text": "Hypopharyngeal cancer causing severe upper esophagus stenosis",
            "source_kind": "underlying_condition",
            "source_ref": {"file": "underlying.json", "json_path": "[0].underlying[3]"},
        },
    ]


def candidate(family="oral_aspiration_or_anaerobe", decision="review_context_needed"):
    return {
        "organism_name": "Bacteroides fragilis",
        "taxonomy_family": family,
        "decision": decision,
        "analytical_profile": {
            "best_rank_in_retained_universe": 3,
            "selected_positive_test_count": 3,
            "reproducibility_axis": True,
            "cross_molecule_selected": True,
        },
    }


def colonization_profile(
    organism="Candida albicans",
    category="explicit_candidate_colonization",
    specimen_scope="respiratory_nonsterile",
    direct=True,
):
    record = {
        "evidence_id": "P50:COLONIZATION:1",
        "text": f"{organism} colonization in BALF",
        "categories": [category],
        "organism_name": organism if not direct else None,
        "organism_key": None if direct else organism.lower().replace(" ", ""),
        "specimen_type": "BALF" if specimen_scope == "respiratory_nonsterile" else "blood",
        "specimen_category": None,
        "specimen_scope": specimen_scope,
        "source_kind": (
            "underlying_condition" if direct else "provisional_culture_agent_label"
        ),
        "source_refs": [],
        "direct_source_explicit": direct,
    }
    return colonization_history_profile([record], POLICY_V4B)


def colonization_candidate(decision="review_context_needed"):
    item = candidate(family="other", decision=decision)
    item["organism_name"] = "Candida albicans"
    item["specimen_context"] = "lower_respiratory"
    item["analytical_profile"].update({
        "best_rank_in_retained_universe": 15,
        "selected_positive_test_count": 1,
        "reproducibility_axis": False,
        "cross_molecule_selected": False,
    })
    return item


class HistoryRouteTests(unittest.TestCase):
    def test_p30_like_aspiration_profile_moves_context_to_priority(self):
        profile = aspiration_history_profile(history_records(), POLICY)
        result = adjust_candidate(candidate(), profile, POLICY)
        self.assertEqual(result["analytic_route"], "Context")
        self.assertEqual(result["history_adjusted_route"], "Priority")
        self.assertEqual(result["history_routed_decision"], "review_high_priority")
        self.assertEqual(result["downstream_route"]["scorer_review_tier"], "High")
        self.assertEqual(result["downstream_route"]["projected_final_tier"], "High")
        self.assertFalse(result["downstream_route"]["picked_allowed_from_history_alone"])
        self.assertTrue(result["downstream_route"]["picked_promotion_blockers"])
        self.assertTrue(result["consumed_evidence_ids"])
        self.assertTrue(all(
            row["downstream_reuse"] == "provenance_only_not_independent_support"
            for row in result["history_evidence_consumption"]
        ))

    def test_aspiration_history_does_not_move_unrelated_family(self):
        profile = aspiration_history_profile(history_records(), POLICY)
        result = adjust_candidate(
            candidate(family="hospital_or_nonfermenter_gnb"), profile, POLICY
        )
        self.assertEqual(result["analytic_route"], result["history_adjusted_route"])
        self.assertFalse(result["consumed_evidence_ids"])

    def test_prior_event_without_durable_risk_does_not_reach_priority(self):
        profile = aspiration_history_profile(history_records()[:1], POLICY)
        result = adjust_candidate(candidate(), profile, POLICY)
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertIn("without durable anatomic risk", " ".join(result["history_adjustment_blockers"]))

    def test_hold_moves_only_one_level(self):
        profile = aspiration_history_profile(history_records(), POLICY)
        result = adjust_candidate(candidate(decision="hold"), profile, POLICY)
        self.assertEqual(result["analytic_route"], "Hold")
        self.assertEqual(result["history_adjusted_route"], "Audit")

    def test_low_confidence_no_is_not_negative_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = {
                "all_phenotype_overview": [{
                    "phenotype": "ASPIRATION_RISK",
                    "status": "NO",
                    "confidence": "LOW",
                    "temporal_relation": "UNCLEAR",
                    "evidence": None,
                }]
            }
            path = root / "patient_30_candidate_phenotype_shadow_v1.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            result = phenotype_disagreement(30, root)
            self.assertEqual(result["interpretation"], "not_established_not_negative_evidence")

    def test_audit_to_context_uses_two_strict_cluster_representatives(self):
        records = [{
            "evidence_id": "P18:HIST:admission:0",
            "text": "Acute respiratory failure, favor pneumonia or choking related",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0].diagnosis"},
        }]
        profile = aspiration_history_profile(records, POLICY)
        values = []
        for name, rank in (
            ("Dialister invisus", 6),
            ("Porphyromonas catoniae", 8),
            ("Prevotella nigrescens", 10),
        ):
            item = candidate(decision="review_low_specificity")
            item["organism_name"] = name
            item["analytical_profile"].update({
                "best_rank_in_retained_universe": rank,
                "selected_positive_test_count": 1,
                "reproducibility_axis": False,
                "cross_molecule_selected": False,
            })
            values.append(adjust_candidate(item, profile, POLICY))
        audit = apply_aspiration_cluster_guardrail(values, POLICY)
        self.assertEqual(
            audit["retained_representatives"],
            ["Dialister invisus", "Porphyromonas catoniae"],
        )
        self.assertEqual(audit["suppressed_to_audit"], ["Prevotella nigrescens"])
        self.assertEqual(
            [item["history_adjusted_route"] for item in values],
            ["Context", "Context", "Audit"],
        )

    def test_broad_oral_organism_is_not_promoted_from_audit(self):
        records = [{
            "evidence_id": "P18:HIST:admission:0",
            "text": "Acute respiratory failure, favor pneumonia or choking related",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0].diagnosis"},
        }]
        profile = aspiration_history_profile(records, POLICY)
        item = candidate(decision="review_low_specificity")
        item["organism_name"] = "Neisseria sicca"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 2,
            "selected_positive_test_count": 1,
            "reproducibility_axis": False,
            "cross_molecule_selected": False,
        })
        result = adjust_candidate(item, profile, POLICY)
        self.assertEqual(result["history_adjusted_route"], "Audit")
        self.assertIn("not a strict aspiration-anaerobe", " ".join(result["history_adjustment_blockers"]))

    def test_high_confidence_opportunistic_history_moves_reproducible_context(self):
        records = [
            {
                "evidence_id": "P28:HIST:admission:0",
                "text": "Acute leukemia, with pancytopenia",
                "source_kind": "admission_diagnosis",
                "source_ref": {"file": "admission.json", "json_path": "[0].diagnosis"},
            }
        ]
        opportunistic = opportunistic_history_profile(records, POLICY_V2)
        item = candidate(
            family="mold_or_opportunistic_fungus",
            decision="review_context_needed",
        )
        item["organism_name"] = "Aspergillus fumigatus"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 5,
            "selected_positive_test_count": 2,
            "reproducibility_axis": True,
            "cross_molecule_selected": True,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V2),
            POLICY_V2,
            opportunistic_history=opportunistic,
        )
        self.assertTrue(opportunistic["has_high_confidence_immunosuppression"])
        self.assertEqual(result["history_adjusted_route"], "Priority")
        self.assertIn("HIST-OPP3-CONTEXT-TO-PRIORITY", result["history_rule_ids"])
        self.assertTrue(result["consumed_evidence_ids"])
        self.assertFalse(result["downstream_route"]["picked_allowed_from_history_alone"])

    def test_nonspecific_host_susceptibility_does_not_change_route(self):
        records = [
            {
                "evidence_id": "P12:HIST:underlying:0",
                "text": "Diabetes mellitus and end-stage renal disease",
                "source_kind": "underlying_condition",
                "source_ref": {"file": "underlying.json", "json_path": "[0]"},
            }
        ]
        opportunistic = opportunistic_history_profile(records, POLICY_V2)
        item = candidate(family="mold_or_opportunistic_fungus")
        item["organism_name"] = "Aspergillus versicolor"
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V2),
            POLICY_V2,
            opportunistic_history=opportunistic,
        )
        self.assertFalse(opportunistic["qualifies_for_bounded_adjustment"])
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertFalse(result["consumed_evidence_ids"])

    def test_strong_host_risk_does_not_rescue_rank_23_signal(self):
        records = [
            {
                "evidence_id": "P33:HIST:underlying:0",
                "text": "B-cell lymphoma with marrow involvement, stage IV",
                "source_kind": "underlying_condition",
                "source_ref": {"file": "underlying.json", "json_path": "[0]"},
            }
        ]
        opportunistic = opportunistic_history_profile(records, POLICY_V2)
        item = candidate(
            family="high_consequence_opportunistic",
            decision="review_context_needed",
        )
        item["organism_name"] = "Pneumocystis jirovecii"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 23,
            "selected_positive_test_count": 3,
            "reproducibility_axis": True,
            "cross_molecule_selected": True,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V2),
            POLICY_V2,
            opportunistic_history=opportunistic,
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertIn("rank at most ten", " ".join(result["history_adjustment_blockers"]))
        self.assertFalse(result["consumed_evidence_ids"])

    def test_opportunistic_single_test_audit_signal_stays_audit(self):
        records = [
            {
                "evidence_id": "P24:HIST:underlying:0",
                "text": "Rituximab chemotherapy",
                "source_kind": "underlying_condition",
                "source_ref": {"file": "underlying.json", "json_path": "[0]"},
            }
        ]
        opportunistic = opportunistic_history_profile(records, POLICY_V2)
        item = candidate(
            family="high_consequence_opportunistic",
            decision="review_low_specificity",
        )
        item["organism_name"] = "Pneumocystis jirovecii"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 8,
            "selected_positive_test_count": 1,
            "reproducibility_axis": False,
            "cross_molecule_selected": False,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V2),
            POLICY_V2,
            opportunistic_history=opportunistic,
        )
        self.assertEqual(result["history_adjusted_route"], "Audit")
        self.assertFalse(result["consumed_evidence_ids"])

    def test_opportunistic_cluster_is_capped_at_two_representatives(self):
        records = [
            {
                "evidence_id": "P28:HIST:admission:0",
                "text": "Acute leukemia, with pancytopenia",
                "source_kind": "admission_diagnosis",
                "source_ref": {"file": "admission.json", "json_path": "[0]"},
            }
        ]
        opportunistic = opportunistic_history_profile(records, POLICY_V2)
        values = []
        for name, rank in (
            ("Aspergillus fumigatus", 2),
            ("Aspergillus terreus", 3),
            ("Trichosporon asahii", 4),
        ):
            item = candidate(family="mold_or_opportunistic_fungus")
            item["organism_name"] = name
            item["analytical_profile"].update({
                "best_rank_in_retained_universe": rank,
                "selected_positive_test_count": 2,
                "reproducibility_axis": True,
                "cross_molecule_selected": True,
            })
            values.append(adjust_candidate(
                item,
                aspiration_history_profile([], POLICY_V2),
                POLICY_V2,
                opportunistic_history=opportunistic,
            ))
        audit = apply_opportunistic_cluster_guardrail(values, POLICY_V2)
        self.assertEqual(len(audit["retained_representatives"]), 2)
        self.assertEqual(audit["suppressed_to_analytic_route"], ["Trichosporon asahii"])
        self.assertEqual(
            [item["history_adjusted_route"] for item in values],
            ["Priority", "Priority", "Context"],
        )

    def test_hap_setting_moves_strong_rank_five_context_to_priority(self):
        records = [{
            "evidence_id": "P35:HIST:admission:0",
            "text": "Sepsis, focus on HAP (hospital-acquired pneumonia)",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0]"},
        }]
        hap_mdr = hap_mdr_history_profile(records, POLICY_V3)
        item = candidate(
            family="hospital_or_nonfermenter_gnb",
            decision="review_context_needed",
        )
        item["organism_name"] = "Stenotrophomonas maltophilia"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 5,
            "selected_positive_test_count": 4,
            "reproducibility_axis": True,
            "cross_molecule_selected": True,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V3),
            POLICY_V3,
            opportunistic_history=opportunistic_history_profile([], POLICY_V3),
            hap_mdr_history=hap_mdr,
        )
        self.assertEqual(result["history_adjusted_route"], "Priority")
        self.assertIn("HIST-HAP3-SETTING-CONTEXT-TO-PRIORITY", result["history_rule_ids"])
        self.assertTrue(result["consumed_evidence_ids"])

    def test_hap_setting_does_not_move_rank_seven_context(self):
        records = [{
            "evidence_id": "P32:HIST:underlying:0",
            "text": "status post tracheostomy",
            "source_kind": "underlying_condition",
            "source_ref": {"file": "underlying.json", "json_path": "[0]"},
        }]
        hap_mdr = hap_mdr_history_profile(records, POLICY_V3)
        item = candidate(family="hospital_or_nonfermenter_gnb")
        item["organism_name"] = "Stenotrophomonas maltophilia"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 7,
            "selected_positive_test_count": 2,
            "reproducibility_axis": True,
            "cross_molecule_selected": False,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V3),
            POLICY_V3,
            opportunistic_history=opportunistic_history_profile([], POLICY_V3),
            hap_mdr_history=hap_mdr,
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertFalse(result["consumed_evidence_ids"])

    def test_candidate_matched_crab_uses_exact_mdr_rank_ten_route(self):
        records = [{
            "evidence_id": "P19:HIST:admission:0",
            "text": "Prior UTI due to CRAB (carbapenem-resistant Acinetobacter baumannii)",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0]"},
        }]
        hap_mdr = hap_mdr_history_profile(records, POLICY_V3)
        item = candidate(family="hospital_or_nonfermenter_gnb")
        item["organism_name"] = "Acinetobacter baumannii"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 8,
            "selected_positive_test_count": 2,
            "reproducibility_axis": True,
            "cross_molecule_selected": False,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V3),
            POLICY_V3,
            opportunistic_history=opportunistic_history_profile([], POLICY_V3),
            hap_mdr_history=hap_mdr,
        )
        self.assertEqual(result["history_adjusted_route"], "Priority")
        self.assertIn("HIST-MDR3-EXACT-MATCH-CONTEXT-TO-PRIORITY", result["history_rule_ids"])

    def test_crab_does_not_apply_to_unmatched_pseudomonas(self):
        records = [{
            "evidence_id": "P19:HIST:admission:0",
            "text": "Prior UTI due to CRAB (carbapenem-resistant Acinetobacter baumannii)",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0]"},
        }]
        hap_mdr = hap_mdr_history_profile(records, POLICY_V3)
        item = candidate(family="hospital_or_nonfermenter_gnb")
        item["organism_name"] = "Pseudomonas aeruginosa"
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V3),
            POLICY_V3,
            opportunistic_history=opportunistic_history_profile([], POLICY_V3),
            hap_mdr_history=hap_mdr,
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertFalse(result["history_evidence_ids"])

    def test_hap_cluster_is_capped_at_two_representatives(self):
        records = [{
            "evidence_id": "P35:HIST:admission:0",
            "text": "Sepsis, focus on HAP (hospital-acquired pneumonia)",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0]"},
        }]
        hap_mdr = hap_mdr_history_profile(records, POLICY_V3)
        values = []
        for name, rank in (
            ("Pseudomonas aeruginosa", 2),
            ("Acinetobacter baumannii", 3),
            ("Stenotrophomonas maltophilia", 4),
        ):
            item = candidate(family="hospital_or_nonfermenter_gnb")
            item["organism_name"] = name
            item["analytical_profile"].update({
                "best_rank_in_retained_universe": rank,
                "selected_positive_test_count": 2,
                "reproducibility_axis": True,
                "cross_molecule_selected": True,
            })
            values.append(adjust_candidate(
                item,
                aspiration_history_profile([], POLICY_V3),
                POLICY_V3,
                opportunistic_history=opportunistic_history_profile([], POLICY_V3),
                hap_mdr_history=hap_mdr,
            ))
        audit = apply_hap_mdr_cluster_guardrail(values, POLICY_V3)
        self.assertEqual(len(audit["retained_representatives"]), 2)
        self.assertEqual(
            audit["suppressed_to_analytic_route"],
            ["Stenotrophomonas maltophilia"],
        )

    def test_prior_hsv_history_is_context_only(self):
        records = [{
            "evidence_id": "P21:HIST:underlying:0",
            "text": "Suspect herpes simplex infection history",
            "source_kind": "underlying_condition",
            "source_ref": {"file": "underlying.json", "json_path": "[0]"},
        }]
        reactivation = reactivation_history_profile(records, POLICY_V4)
        item = candidate(family="herpes_or_reactivation_virus")
        item["organism_name"] = "Human alphaherpesvirus 1"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 1,
            "selected_positive_test_count": 4,
            "reproducibility_axis": True,
            "cross_molecule_selected": True,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V4),
            POLICY_V4,
            reactivation_history=reactivation,
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertIn("context only", " ".join(result["history_adjustment_blockers"]))
        self.assertTrue(result["history_evidence_ids"])
        self.assertFalse(result["consumed_evidence_ids"])

    def test_resolved_reactivation_demotes_only_weak_context_signal(self):
        records = [{
            "evidence_id": "P40:HIST:admission:0",
            "text": "Resolved HSV reactivation with asymptomatic HSV shedding and no end-organ disease",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0]"},
        }]
        reactivation = reactivation_history_profile(records, POLICY_V4)
        item = candidate(family="herpes_or_reactivation_virus")
        item["organism_name"] = "Human alphaherpesvirus 1"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 15,
            "selected_positive_test_count": 1,
            "reproducibility_axis": False,
            "cross_molecule_selected": False,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V4),
            POLICY_V4,
            reactivation_history=reactivation,
        )
        self.assertEqual(result["analytic_route"], "Context")
        self.assertEqual(result["history_adjusted_route"], "Audit")
        self.assertIn(
            "HIST-REACT1-RESOLVED-WEAK-CONTEXT-TO-AUDIT",
            result["history_rule_ids"],
        )
        self.assertTrue(result["consumed_evidence_ids"])

    def test_resolved_reactivation_does_not_demote_strong_signal(self):
        records = [{
            "evidence_id": "P40:HIST:admission:0",
            "text": "Resolved CMV reactivation without end-organ disease",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0]"},
        }]
        reactivation = reactivation_history_profile(records, POLICY_V4)
        item = candidate(family="herpes_or_reactivation_virus")
        item["organism_name"] = "Human betaherpesvirus 5"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 2,
            "selected_positive_test_count": 3,
            "reproducibility_axis": True,
            "cross_molecule_selected": True,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V4),
            POLICY_V4,
            reactivation_history=reactivation,
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertFalse(result["consumed_evidence_ids"])

    def test_convergent_current_reactivation_moves_context_to_priority_and_flags_review(self):
        records = [{
            "evidence_id": "P40:HIST:admission:0",
            "text": "HSV reactivation with rising quantitative PCR viral load and pneumonia",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0]"},
        }]
        reactivation = reactivation_history_profile(records, POLICY_V4)
        item = candidate(family="herpes_or_reactivation_virus")
        item["organism_name"] = "Human alphaherpesvirus 1"
        item["analytical_profile"].update({
            "best_rank_in_retained_universe": 1,
            "selected_positive_test_count": 3,
            "reproducibility_axis": True,
            "cross_molecule_selected": True,
        })
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V4),
            POLICY_V4,
            reactivation_history=reactivation,
        )
        self.assertEqual(result["history_adjusted_route"], "Priority")
        self.assertIn(
            "HIST-REACT2-CONVERGENT-CONTEXT-TO-PRIORITY",
            result["history_rule_ids"],
        )
        self.assertTrue(result["future_picked_gate_review_required"])
        self.assertEqual(
            result["future_picked_gate_review_status"],
            "pending_independent_evidence_audit",
        )
        self.assertTrue(result["future_picked_gate_review_requirements"])
        self.assertFalse(result["downstream_route"]["picked_allowed_from_history_alone"])

    def test_reactivation_history_must_match_candidate(self):
        records = [{
            "evidence_id": "P40:HIST:admission:0",
            "text": "CMV reactivation with rising quantitative PCR viral load and pneumonia",
            "source_kind": "admission_diagnosis",
            "source_ref": {"file": "admission.json", "json_path": "[0]"},
        }]
        reactivation = reactivation_history_profile(records, POLICY_V4)
        item = candidate(family="herpes_or_reactivation_virus")
        item["organism_name"] = "Human alphaherpesvirus 1"
        result = adjust_candidate(
            item,
            aspiration_history_profile([], POLICY_V4),
            POLICY_V4,
            reactivation_history=reactivation,
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertFalse(result["history_evidence_ids"])

    def test_v4b_overlay_inherits_v4a_reactivation_rules(self):
        self.assertEqual(
            POLICY_V4B["reactivation_family"], POLICY_V4["reactivation_family"]
        )
        self.assertIn("policy_composition", POLICY_V4B)
        self.assertGreater(
            len(POLICY_V4B["constraints"]), len(POLICY_V4["constraints"])
        )

    def test_v4c_combines_all_history_modules_without_behavior_change(self):
        self.assertEqual(POLICY_V4C["policy_role"], "frozen_combined_history_route_shadow")
        self.assertEqual(POLICY_V4C["reactivation_family"], POLICY_V4B["reactivation_family"])
        self.assertEqual(
            POLICY_V4C["colonization_history_patterns"],
            POLICY_V4B["colonization_history_patterns"],
        )
        self.assertEqual(
            POLICY_V4C["merge_contract"]["behavior_change_from_v4b"], False
        )
        self.assertEqual(len(POLICY_V4C["combined_modules"]), 5)
        self.assertGreaterEqual(len(POLICY_V4C["policy_composition_chain"]), 2)

    def test_overlay_base_hash_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "base.json").write_text(
                json.dumps({"policy_id": "base", "constraints": []}),
                encoding="utf-8",
            )
            (root / "overlay.json").write_text(json.dumps({
                "base_policy": "base.json",
                "base_policy_sha256": "0" * 64,
                "policy_id": "overlay",
            }), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                read_json(root / "overlay.json")

    def test_structural_airway_context_alone_does_not_demote(self):
        profile = colonization_history_profile([{
            "evidence_id": "P50:HIST:1",
            "text": "COPD with tracheostomy",
            "categories": ["structural_airway_context_only"],
            "organism_name": None,
            "organism_key": None,
            "specimen_scope": "unknown",
            "source_kind": "underlying_condition",
            "source_refs": [],
            "direct_source_explicit": True,
        }], POLICY_V4B)
        result = adjust_candidate(
            colonization_candidate(), aspiration_history_profile([], POLICY_V4B),
            POLICY_V4B, colonization_history=profile,
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertFalse(
            result["colonization_interpretation"]["has_candidate_match"]
        )

    def test_mismatched_colonization_organism_does_not_demote(self):
        result = adjust_candidate(
            colonization_candidate(), aspiration_history_profile([], POLICY_V4B),
            POLICY_V4B,
            colonization_history=colonization_profile("Pseudomonas aeruginosa"),
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertFalse(result["colonization_interpretation"]["has_candidate_match"])

    def test_explicit_same_compartment_weak_colonization_moves_context_to_audit(self):
        result = adjust_candidate(
            colonization_candidate(), aspiration_history_profile([], POLICY_V4B),
            POLICY_V4B, colonization_history=colonization_profile(),
        )
        self.assertEqual(result["history_adjusted_route"], "Audit")
        self.assertIn(
            "HIST-COL2-EXPLICIT-WEAK-CONTEXT-TO-AUDIT",
            result["history_rule_ids"],
        )
        self.assertEqual(result["history_policy_family"], "colonization_interpretation")
        self.assertTrue(result["consumed_evidence_ids"])

    def test_explicit_same_compartment_very_weak_colonization_moves_audit_to_hold(self):
        result = adjust_candidate(
            colonization_candidate("review_low_specificity"),
            aspiration_history_profile([], POLICY_V4B), POLICY_V4B,
            colonization_history=colonization_profile(),
        )
        self.assertEqual(result["history_adjusted_route"], "Hold")
        self.assertIn(
            "HIST-COL3-EXPLICIT-VERY-WEAK-AUDIT-TO-HOLD",
            result["history_rule_ids"],
        )

    def test_explicit_current_not_infection_moves_priority_to_context(self):
        profile = colonization_profile(category="explicit_current_not_infection")
        profile["matched_evidence"][0]["text"] = (
            "Candida albicans in BALF is not infection"
        )
        item = colonization_candidate("review_high_priority")
        result = adjust_candidate(
            item, aspiration_history_profile([], POLICY_V4B), POLICY_V4B,
            colonization_history=profile,
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertIn(
            "HIST-COL1-CURRENT-NOT-INFECTION-PRIORITY-TO-CONTEXT",
            result["history_rule_ids"],
        )

    def test_exact_hospital_support_blocks_colonization_demotion(self):
        item = colonization_candidate()
        item["hospital_profile"] = {"exact_or_alias_evidence": True}
        result = adjust_candidate(
            item, aspiration_history_profile([], POLICY_V4B), POLICY_V4B,
            colonization_history=colonization_profile(),
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertIn("culture/PCR", " ".join(result["history_adjustment_blockers"]))

    def test_provisional_sterile_contaminant_label_cannot_demote(self):
        profile = colonization_profile(
            category="provisional_model_contaminant_label",
            specimen_scope="sterile_or_systemic",
            direct=False,
        )
        result = adjust_candidate(
            colonization_candidate(), aspiration_history_profile([], POLICY_V4B),
            POLICY_V4B, colonization_history=profile,
        )
        self.assertEqual(result["history_adjusted_route"], "Context")
        self.assertIn(
            "sterile or systemic", " ".join(result["history_adjustment_blockers"])
        )


if __name__ == "__main__":
    unittest.main()

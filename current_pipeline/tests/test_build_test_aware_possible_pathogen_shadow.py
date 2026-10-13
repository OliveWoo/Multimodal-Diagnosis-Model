import copy
import json
import unittest
from pathlib import Path

from tools.build_test_aware_possible_pathogen_shadow import (
    ANALYTICAL_DIRECT_ROLE,
    ASPERGILLUS_GROUP_ROLE,
    DEFAULT_POLICY,
    FALLBACK_ROLE,
    LOW_SPECIFICITY_ROLE,
    MNGS_ONLY_ROLE,
    NOT_SELECTED_ROLE,
    OPPORTUNISTIC_FUNGUS_ROLE,
    PJP_ROLE,
    REACTIVATION_ROLE,
    RESPIRATORY_VIRUS_ROLE,
    STRICT_ROLE,
    SYSTEMIC_ROLE,
    annotate_candidate,
    apply_best_available_fallback,
    apply_graceful_degradation,
    decision_evidence_profile,
    load_policy,
)


GRACEFUL_POLICY = Path(
    "rules/test_aware_unified_reporting_v4_graceful_degradation.json"
)
FAMILY_AWARE_POLICY = Path(
    "rules/test_aware_unified_reporting_v5_family_aware_fallback.json"
)
EVIDENCE_CONVERGENCE_POLICY = Path(
    "rules/test_aware_unified_reporting_v6_evidence_convergence.json"
)
RESPIRATORY_VIRUS_POLICY = Path(
    "rules/test_aware_unified_reporting_v7_respiratory_virus_possible.json"
)
ABSTENTION_POLICY = Path(
    "rules/test_aware_unified_reporting_v8_abstention.json"
)
HOST_FUNGUS_POLICY = Path(
    "rules/test_aware_unified_reporting_v9_host_context_opportunistic_fungus.json"
)


def candidate(
    *, family="hospital_or_nonfermenter_gnb", decision="review_high_priority",
    rank=2, count=2, reproducible=True, cross=False, direct_level=None,
    direct_rows=None, context="lower_respiratory",
):
    role = (
        "existing_picked_pathogen" if decision == "picked_shadow"
        else "high_priority_review_candidate" if decision == "review_high_priority"
        else "context_review_candidate"
    )
    return {
        "patient_id": "1",
        "organism_name": "Pseudomonas aeruginosa",
        "clinical_decision": decision,
        "clinical_level": "Level 2" if decision == "picked_shadow" else "Level 3",
        "formal_pick_allowed": decision == "picked_shadow",
        "selection_role": role,
        "taxonomy_profile": {
            "primary_rule_family": family,
            "taxonomic_rank": "species",
            "mapping_status": "exact_species",
        },
        "source_candidate": {"specimen_context": context},
        "promotion_gate": {
            "axes": {
                "taxonomy_family": family,
                "precise_identity": True,
                "specimen_context": context,
                "best_rank_in_retained_universe": rank,
                "selected_positive_test_count": count,
                "reproducibility_axis": reproducible,
                "cross_molecule_selected": cross,
                "direct_hospital_level": direct_level,
                "direct_evidence_timing_profile": {
                    "rows": direct_rows or [],
                },
            }
        },
    }


class TestPossiblePathogenShadow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.policy = json.loads(Path(DEFAULT_POLICY).read_text(encoding="utf-8"))
        cls.graceful_policy = json.loads(
            GRACEFUL_POLICY.read_text(encoding="utf-8")
        )
        cls.family_aware_policy = json.loads(
            FAMILY_AWARE_POLICY.read_text(encoding="utf-8")
        )
        cls.evidence_convergence_policy = json.loads(
            EVIDENCE_CONVERGENCE_POLICY.read_text(encoding="utf-8")
        )
        cls.respiratory_virus_policy = json.loads(
            RESPIRATORY_VIRUS_POLICY.read_text(encoding="utf-8")
        )
        cls.abstention_policy, _ = load_policy(ABSTENTION_POLICY)
        cls.host_fungus_policy, _ = load_policy(HOST_FUNGUS_POLICY)

    @staticmethod
    def add_recorded_immunosuppressants(item, count=2):
        medications = ["prednisolone", "methotrexate"][:count]
        item["host_evidence"] = [
            {
                "evidence_id": f"HOST-{index}",
                "evidence_role": "host_medication",
                "label": label,
                "assertion": (
                    "recorded_pre_sample_exposure_not_administration_verified"
                ),
                "observed_at": f"2023-{10 + index:02d}",
                "source_ref": {
                    "source_claim_group": f"admission:{index}"
                },
            }
            for index, label in enumerate(medications, start=1)
        ]
        return item

    @staticmethod
    def add_mngs_signal(
        item, *, rank=1, priority=None, qc_status="partial_evaluable",
        reads=100,
    ):
        item["analytical_profile"] = {
            "per_test_signals": [{
                "seq_id": "TEST-1",
                "reads": reads,
                "rank_in_retained_universe": rank,
                "analytical_rank_priority": priority,
                "rpm_total": None,
                "qc_status": qc_status,
            }]
        }
        return item

    def test_strict_pick_is_preserved_without_mutating_clinical_fields(self):
        source = candidate(decision="picked_shadow")
        result = annotate_candidate(source, self.policy)
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual("existing_picked_pathogen", result["selection_role"])
        self.assertEqual(STRICT_ROLE, result["possible_reporting_role"])
        self.assertTrue(result["selected_for_complete_report"])

    def test_v9_pjp_uses_two_recorded_pre_sample_medications_for_possible_only(self):
        item = self.add_recorded_immunosuppressants(
            candidate(
                family="high_consequence_opportunistic",
                decision="review_context_needed",
                rank=2,
                count=2,
                reproducible=True,
                cross=False,
            )
        )
        item["organism_name"] = "Pneumocystis jirovecii"

        result = annotate_candidate(item, self.host_fungus_policy)

        self.assertEqual(PJP_ROLE, result["possible_reporting_role"])
        self.assertEqual("Possible", result["final_reporting_tier"])
        self.assertFalse(result["formal_pick_allowed"])
        host = result["possible_pathogen_gate"]["axes"][
            "recorded_host_review_context"
        ]
        self.assertTrue(host["available"])
        self.assertFalse(host["administration_verified"])
        self.assertFalse(host["missing_dose_is_negative"])

    def test_v9_single_unverified_medication_is_not_enough_for_pjp(self):
        item = self.add_recorded_immunosuppressants(
            candidate(
                family="high_consequence_opportunistic",
                decision="review_context_needed",
                rank=2,
                count=2,
                reproducible=True,
                cross=False,
            ),
            count=1,
        )
        item["organism_name"] = "Pneumocystis jirovecii"

        result = annotate_candidate(item, self.host_fungus_policy)

        self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_v9_lower_respiratory_vzv_reproducibility_plus_host_is_possible(self):
        item = self.add_recorded_immunosuppressants(
            candidate(
                family="herpes_or_reactivation_virus",
                decision="review_context_needed",
                rank=3,
                count=2,
                reproducible=True,
                cross=False,
                context="lower_respiratory",
            )
        )
        item["organism_name"] = "Human alphaherpesvirus 3"

        result = annotate_candidate(item, self.host_fungus_policy)

        self.assertEqual(REACTIVATION_ROLE, result["possible_reporting_role"])
        self.assertEqual(
            "reactivation_family_reproducible_signal_with_recorded_host_context",
            result["possible_pathogen_gate"]["route"],
        )
        self.assertFalse(result["formal_pick_allowed"])

    def test_v9_provisional_phenotype_alone_does_not_promote_respiratory_herpes(self):
        item = candidate(
            family="herpes_or_reactivation_virus",
            decision="review_context_needed",
            rank=1,
            count=3,
            reproducible=True,
            cross=False,
            context="lower_respiratory",
        )
        item["organism_name"] = "Human alphaherpesvirus 1"
        item["phenotype_evidence"] = {
            "production_link_verified": False,
            "routed_phenotype_context": [{
                "source_id": "PHENO-1",
                "phenotype": "IMMUNOCOMPROMISED",
                "status": "YES",
                "confidence": "HIGH",
                "temporal_relation": "PRE_PNEUMONIA",
                "phenotype_qc_status": "PASS",
            }],
        }

        result = annotate_candidate(item, self.host_fungus_policy)

        self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_v9_exact_opportunistic_fungus_with_direct_convergence_is_possible(self):
        item = candidate(
            family="mold_or_opportunistic_fungus",
            decision="review_context_needed",
            rank=1,
            count=4,
            reproducible=True,
            cross=True,
            direct_level=3,
            direct_rows=[{
                "observation_id": "CUL-1",
                "positive": True,
                "timing": "within_event_window",
                "specimen_type": "BALF",
                "target": "Trichosporon asahii",
            }],
        )
        item["organism_name"] = "Trichosporon asahii"

        result = annotate_candidate(item, self.host_fungus_policy)

        self.assertEqual(
            OPPORTUNISTIC_FUNGUS_ROLE, result["possible_reporting_role"]
        )
        self.assertEqual("Possible", result["final_reporting_tier"])
        self.assertFalse(result["formal_pick_allowed"])

    def test_v9_opportunistic_fungus_without_direct_evidence_stays_context(self):
        item = candidate(
            family="mold_or_opportunistic_fungus",
            decision="review_context_needed",
            rank=1,
            count=4,
            reproducible=True,
            cross=True,
            direct_level=3,
        )
        item["organism_name"] = "Trichosporon asahii"

        result = annotate_candidate(item, self.host_fungus_policy)

        self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_reproducible_top_three_typical_bacterium_is_possible(self):
        result = annotate_candidate(candidate(), self.policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])
        self.assertEqual(ANALYTICAL_DIRECT_ROLE, result["possible_reporting_role"])
        self.assertFalse(result["formal_pick_allowed"])

    def test_event_aligned_direct_positive_without_level_is_possible(self):
        rows = [{
            "observation_id": "ASSAY-1",
            "positive": True,
            "timing": "within_event_window",
            "specimen_type": "BAL",
        }]
        result = annotate_candidate(
            candidate(rank=8, count=3, direct_rows=rows), self.policy
        )
        self.assertEqual(ANALYTICAL_DIRECT_ROLE, result["possible_reporting_role"])
        self.assertIn("generic hospital level is unavailable", " ".join(
            result["possible_pathogen_gate"]["cautions"]
        ))

    def test_level_three_direct_culture_is_possible(self):
        rows = [{
            "observation_id": "CUL-1",
            "positive": True,
            "timing": "within_event_window",
            "specimen_type": "Lower BAL",
        }]
        result = annotate_candidate(
            candidate(
                family="environmental_low_specificity", rank=1, count=3,
                cross=True, direct_level=3, direct_rows=rows,
            ),
            self.policy,
        )
        self.assertEqual(ANALYTICAL_DIRECT_ROLE, result["possible_reporting_role"])

    def test_context_top_rank_cross_molecule_low_specificity_is_possible(self):
        result = annotate_candidate(
            candidate(
                family="skin_airway_colonizer_prone",
                decision="review_context_needed", rank=1, count=3, cross=True,
            ),
            self.policy,
        )
        self.assertEqual(LOW_SPECIFICITY_ROLE, result["possible_reporting_role"])
        self.assertEqual("review_context_needed", result["clinical_decision"])

    def test_high_top_rank_cross_molecule_low_specificity_stays_possible(self):
        result = annotate_candidate(
            candidate(
                family="skin_airway_colonizer_prone",
                decision="review_high_priority", rank=1, count=3, cross=True,
            ),
            self.policy,
        )
        self.assertEqual(LOW_SPECIFICITY_ROLE, result["possible_reporting_role"])
        self.assertEqual("review_high_priority", result["clinical_decision"])

    def test_insufficient_context_signal_stays_unselected(self):
        result = annotate_candidate(
            candidate(
                family="skin_airway_colonizer_prone",
                decision="review_context_needed", rank=2, count=2, cross=False,
            ),
            self.policy,
        )
        self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_special_families_cannot_use_generic_route(self):
        for family in self.policy["special_family_guardrails"]:
            with self.subTest(family=family):
                item = candidate(family=family, rank=1, count=6, cross=True)
                result = annotate_candidate(item, self.policy)
                self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_v6_cmv_cross_molecule_requires_independent_clinical_axis(self):
        item = candidate(
            family="herpes_or_reactivation_virus",
            decision="review_context_needed",
            rank=2,
            count=2,
            reproducible=True,
            cross=True,
        )
        item["organism_name"] = "Human betaherpesvirus 5"
        without_support = annotate_candidate(
            item, self.evidence_convergence_policy
        )
        self.assertEqual(NOT_SELECTED_ROLE, without_support["possible_reporting_role"])

        item["promotion_gate"]["axes"]["direct_evidence_timing_profile"]["rows"] = [{
            "observation_id": "PCR-1",
            "positive": True,
            "timing": "within_event_window",
            "specimen_type": "BAL",
            "target": "CMV PCR",
        }]
        with_support = annotate_candidate(item, self.evidence_convergence_policy)
        self.assertEqual(
            "reactivation_or_shedding_possible_pathogen",
            with_support["possible_reporting_role"],
        )

    def test_v6_cmv_cross_molecule_host_risk_alone_is_not_clinical_support(self):
        item = candidate(
            family="herpes_or_reactivation_virus",
            decision="review_context_needed",
            rank=2,
            count=3,
            reproducible=True,
            cross=True,
        )
        item["organism_name"] = "Human betaherpesvirus 5"
        item["promotion_gate"]["axes"]["host_support"] = True
        result = annotate_candidate(item, self.evidence_convergence_policy)
        self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_v6_systemic_rank_one_vzv_with_host_risk_is_possible(self):
        item = candidate(
            family="herpes_or_reactivation_virus",
            decision="review_context_needed",
            rank=1,
            count=1,
            reproducible=False,
            cross=False,
            context="blood_or_systemic",
        )
        item["organism_name"] = "Human alphaherpesvirus 3"
        item["promotion_gate"]["axes"]["host_support"] = True
        result = annotate_candidate(item, self.evidence_convergence_policy)
        self.assertEqual(
            "reactivation_or_shedding_possible_pathogen",
            result["possible_reporting_role"],
        )
        self.assertEqual(
            "reactivation_family_top_rank_with_recorded_host_context",
            result["possible_pathogen_gate"]["route"],
        )

    def test_v6_lower_respiratory_rank_one_vzv_host_risk_stays_context(self):
        item = candidate(
            family="herpes_or_reactivation_virus",
            decision="review_context_needed",
            rank=1,
            count=1,
            reproducible=False,
            cross=False,
            context="lower_respiratory",
        )
        item["organism_name"] = "Human alphaherpesvirus 3"
        item["promotion_gate"]["axes"]["host_support"] = True
        result = annotate_candidate(item, self.evidence_convergence_policy)
        self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_v6_non_cmv_reactivation_route_preserves_six_test_threshold(self):
        item = candidate(
            family="herpes_or_reactivation_virus",
            decision="review_context_needed",
            rank=1,
            count=6,
            reproducible=True,
            cross=True,
            context="lower_respiratory",
        )
        item["organism_name"] = "Human alphaherpesvirus 1"
        result = annotate_candidate(item, self.evidence_convergence_policy)
        self.assertEqual(
            "reactivation_or_shedding_possible_pathogen",
            result["possible_reporting_role"],
        )

        item["promotion_gate"]["axes"]["selected_positive_test_count"] = 5
        result = annotate_candidate(item, self.evidence_convergence_policy)
        self.assertEqual(NOT_SELECTED_ROLE, result["possible_reporting_role"])

    def test_v6_single_top_rank_cmv_with_direct_support_is_possible(self):
        item = candidate(
            family="herpes_or_reactivation_virus",
            decision="review_context_needed",
            rank=1,
            count=1,
            reproducible=False,
            cross=False,
        )
        item["organism_name"] = "Human betaherpesvirus 5"
        item["promotion_gate"]["axes"]["direct_evidence_timing_profile"]["rows"] = [{
            "observation_id": "PCR-1",
            "positive": True,
            "timing": "within_event_window",
            "specimen_type": "Blood",
            "target": "CMV PCR",
        }]
        result = annotate_candidate(item, self.evidence_convergence_policy)
        self.assertEqual(
            "reactivation_or_shedding_possible_pathogen",
            result["possible_reporting_role"],
        )

    def test_v6_aspergillus_group_mycology_and_host_risk_is_possible(self):
        item = candidate(
            family="mold_or_opportunistic_fungus",
            decision="review_high_priority",
            rank=None,
            count=0,
            reproducible=False,
            cross=False,
            direct_level=3,
            context="sterile_or_systemic",
        )
        item["organism_name"] = "Aspergillus spp."
        item["candidate_source"] = "hospital_only"
        item["taxonomy_profile"]["taxonomic_rank"] = "genus_or_group_label"
        item["taxonomy_profile"]["mapping_status"] = "genus_inherited"
        item["promotion_gate"]["axes"]["precise_identity"] = False
        item["promotion_gate"]["axes"]["direct_evidence_timing_profile"]["rows"] = [{
            "observation_id": "GM-1",
            "positive": True,
            "timing": "outside_event_window",
            "collected_time": "2024-11-19",
            "specimen_type": "Blood",
            "target": "Aspergillus galactomannan Ag",
        }]
        item["linked_events"] = [{"collected_time": "2024-11-21 22:00"}]
        item["host_evidence"] = [{
            "evidence_id": "HOST-1",
            "label": "Active malignancy",
            "assertion": "recorded_history",
        }]
        result = annotate_candidate(item, self.evidence_convergence_policy)
        self.assertEqual(ASPERGILLUS_GROUP_ROLE, result["possible_reporting_role"])
        evidence = result["possible_pathogen_gate"]["axes"][
            "aspergillus_group_evidence"
        ]
        self.assertEqual("species_unresolved", evidence["species_resolution"])
        self.assertEqual(
            "date_only_interval_overlaps_window",
            evidence["mycology_rows"][0]["event_alignment_mode"],
        )

    def test_v6_klebsiella_identity_does_not_imply_carbapenem_resistance(self):
        item = candidate()
        item["organism_name"] = "Klebsiella pneumoniae"
        result = annotate_candidate(item, self.evidence_convergence_policy)
        self.assertEqual("unknown", result["resistance_interpretation"]["status"])
        self.assertTrue(
            result["resistance_interpretation"][
                "identity_and_resistance_are_separate"
            ]
        )

    def test_repeated_blood_observations_trigger_systemic_review_with_caution(self):
        rows = [
            {
                "observation_id": "BC-1", "positive": True,
                "timing": "within_event_window", "specimen_type": "Blood",
            },
            {
                "observation_id": "BC-2", "positive": True,
                "timing": "within_event_window", "specimen_type": "Blood",
            },
        ]
        source = candidate(
            family="skin_airway_colonizer_prone",
            decision="review_context_needed", rank=1, count=3, cross=True,
            direct_rows=rows,
        )
        snapshot = copy.deepcopy(source)
        result = annotate_candidate(source, self.policy)
        self.assertEqual(SYSTEMIC_ROLE, result["possible_reporting_role"])
        self.assertFalse(
            result["possible_pathogen_gate"]["axes"]["repeated_sterile_independence_verified"]
        )
        self.assertEqual(snapshot, source)

    def test_mngs_only_strict_pick_is_capped_at_possible(self):
        source = self.add_mngs_signal(candidate(decision="picked_shadow"))
        annotated = annotate_candidate(source, self.graceful_policy)
        profile = decision_evidence_profile([annotated])
        result = apply_graceful_degradation(
            [annotated], self.graceful_policy, evidence_profile=profile
        )[0]
        self.assertEqual("mngs_only", profile["decision_evidence_mode"])
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual(MNGS_ONLY_ROLE, result["possible_reporting_role"])
        self.assertEqual("Possible", result["final_reporting_tier"])
        self.assertTrue(result["selected_for_complete_report"])

    def test_partial_qc_fallback_uses_within_test_rank(self):
        worse = self.add_mngs_signal(
            candidate(decision="hold"), rank=2, priority=1
        )
        worse["organism_name"] = "Pseudomonas aeruginosa"
        better = self.add_mngs_signal(
            candidate(decision="hold"), rank=1, priority=5
        )
        better["organism_name"] = "Klebsiella pneumoniae"
        annotated = [
            annotate_candidate(item, self.graceful_policy)
            for item in (worse, better)
        ]
        profile = decision_evidence_profile(annotated)
        annotated = apply_graceful_degradation(
            annotated, self.graceful_policy, evidence_profile=profile
        )
        result = apply_best_available_fallback(
            annotated, self.graceful_policy, evidence_profile=profile
        )
        selected = [item for item in result if item["selected_for_complete_report"]]
        self.assertEqual(1, len(selected))
        self.assertEqual("Klebsiella pneumoniae", selected[0]["organism_name"])
        self.assertEqual(FALLBACK_ROLE, selected[0]["possible_reporting_role"])
        features = selected[0]["possible_pathogen_gate"]["axes"][
            "fallback_features"
        ]
        self.assertEqual(1, features["best_rankable_rank"])
        self.assertEqual(1, features["rankable_positive_test_count"])

    def test_hospital_only_positive_can_supply_low_confidence_fallback(self):
        source = candidate(decision="hold", count=0, reproducible=False)
        source["candidate_source"] = "hospital_only"
        annotated = annotate_candidate(source, self.graceful_policy)
        profile = decision_evidence_profile([annotated])
        result = apply_best_available_fallback(
            [annotated], self.graceful_policy, evidence_profile=profile
        )[0]
        self.assertEqual("hospital_only", profile["decision_evidence_mode"])
        self.assertEqual(FALLBACK_ROLE, result["possible_reporting_role"])
        self.assertEqual("Fallback-Possible", result["final_reporting_tier"])

    def test_family_aware_fallback_prefers_high_consequence_over_reactivation(self):
        pjp = self.add_mngs_signal(
            candidate(
                family="high_consequence_opportunistic",
                decision="review_context_needed",
                rank=2,
                count=1,
                reproducible=False,
            ),
            rank=2,
        )
        pjp["organism_name"] = "Pneumocystis jirovecii"
        cmv = self.add_mngs_signal(
            candidate(
                family="herpes_or_reactivation_virus",
                decision="review_context_needed",
                rank=2,
                count=1,
                reproducible=False,
            ),
            rank=2,
        )
        cmv["organism_name"] = "Human betaherpesvirus 5"
        annotated = [
            annotate_candidate(item, self.family_aware_policy)
            for item in (cmv, pjp)
        ]
        profile = decision_evidence_profile(annotated)
        result = apply_best_available_fallback(
            annotated, self.family_aware_policy, evidence_profile=profile
        )
        selected = [item for item in result if item["selected_for_complete_report"]]
        self.assertEqual(["Pneumocystis jirovecii"], [item["organism_name"] for item in selected])

    def test_family_priority_never_overrides_clinical_decision_tier(self):
        context_pjp = self.add_mngs_signal(
            candidate(
                family="high_consequence_opportunistic",
                decision="review_context_needed",
                rank=1,
            )
        )
        high_virus = self.add_mngs_signal(
            candidate(
                family="herpes_or_reactivation_virus",
                decision="review_high_priority",
                rank=8,
            ),
            rank=8,
        )
        annotated = [
            annotate_candidate(item, self.family_aware_policy)
            for item in (context_pjp, high_virus)
        ]
        profile = decision_evidence_profile(annotated)
        result = apply_best_available_fallback(
            annotated, self.family_aware_policy, evidence_profile=profile
        )
        selected = [item for item in result if item["selected_for_complete_report"]]
        self.assertEqual("herpes_or_reactivation_virus", selected[0]["possible_pathogen_gate"]["axes"]["taxonomy_family"])

    def test_curated_exact_respiratory_group_label_can_enter_fallback(self):
        virus = self.add_mngs_signal(
            candidate(
                family="other_respiratory_virus",
                decision="review_context_needed",
                rank=2,
                count=1,
                reproducible=False,
            ),
            rank=2,
        )
        virus["organism_name"] = "Severe acute respiratory syndrome-related coronavirus"
        virus["taxonomy_profile"].update({
            "taxonomic_rank": "group_or_species_label",
            "mapping_status": "exact_species",
            "rule_source": "curated_exact_species",
        })
        virus["promotion_gate"]["axes"]["precise_identity"] = False
        bacterium = self.add_mngs_signal(
            candidate(
                family="typical_respiratory_pathogen",
                decision="review_context_needed",
                rank=3,
                count=1,
                reproducible=False,
            ),
            rank=3,
        )
        bacterium["organism_name"] = "Klebsiella pneumoniae"
        annotated = [
            annotate_candidate(item, self.family_aware_policy)
            for item in (bacterium, virus)
        ]
        profile = decision_evidence_profile(annotated)
        result = apply_best_available_fallback(
            annotated, self.family_aware_policy, evidence_profile=profile
        )
        selected = [item for item in result if item["selected_for_complete_report"]]
        self.assertEqual(
            "Severe acute respiratory syndrome-related coronavirus",
            selected[0]["organism_name"],
        )
        self.assertTrue(
            selected[0]["possible_pathogen_gate"]["axes"]["fallback_features"][
                "identity_exception_applied"
            ]
        )

    def test_curated_respiratory_group_label_rna_can_enter_possible(self):
        virus = self.add_mngs_signal(
            candidate(
                family="other_respiratory_virus",
                decision="review_high_priority",
                rank=2,
                count=1,
                reproducible=False,
                context="lower_respiratory",
            ),
            rank=2,
            qc_status="evaluable",
            reads=4050,
        )
        virus["organism_name"] = (
            "Severe acute respiratory syndrome-related coronavirus"
        )
        virus["taxonomy_profile"].update({
            "taxonomic_rank": "group_or_species_label",
            "mapping_status": "exact_species",
            "rule_source": "curated_exact_species",
        })
        virus["promotion_gate"]["axes"]["precise_identity"] = False
        virus["analytical_profile"]["per_test_signals"][0][
            "nucleic_type"
        ] = "RNA"
        result = annotate_candidate(virus, self.respiratory_virus_policy)
        self.assertTrue(result["selected_for_complete_report"])
        self.assertEqual(
            RESPIRATORY_VIRUS_ROLE, result["possible_reporting_role"]
        )
        self.assertEqual(
            "curated_respiratory_virus_group_label_lower_respiratory_rna",
            result["possible_pathogen_gate"]["route"],
        )
        self.assertTrue(
            result["possible_pathogen_gate"]["axes"]
            ["respiratory_virus_identity"]["species_confirmation_required"]
        )

    def test_curated_respiratory_group_label_dna_only_stays_context(self):
        virus = self.add_mngs_signal(
            candidate(
                family="other_respiratory_virus",
                decision="review_high_priority",
                rank=2,
                count=1,
                reproducible=False,
                context="lower_respiratory",
            ),
            rank=2,
            qc_status="evaluable",
            reads=4050,
        )
        virus["organism_name"] = (
            "Severe acute respiratory syndrome-related coronavirus"
        )
        virus["taxonomy_profile"].update({
            "taxonomic_rank": "group_or_species_label",
            "mapping_status": "exact_species",
            "rule_source": "curated_exact_species",
        })
        virus["promotion_gate"]["axes"]["precise_identity"] = False
        virus["analytical_profile"]["per_test_signals"][0][
            "nucleic_type"
        ] = "DNA"
        result = annotate_candidate(virus, self.respiratory_virus_policy)
        self.assertFalse(result["selected_for_complete_report"])

    def test_prefix_inferred_respiratory_species_stays_context(self):
        virus = self.add_mngs_signal(
            candidate(
                family="other_respiratory_virus",
                decision="review_high_priority",
                rank=1,
                count=1,
                reproducible=False,
                context="lower_respiratory",
            ),
            rank=1,
            qc_status="evaluable",
            reads=2000,
        )
        virus["organism_name"] = "rhinovirus A89"
        virus["taxonomy_profile"].update({
            "taxonomic_rank": "species",
            "mapping_status": "key_prefix_inferred",
            "rule_source": "prefix_fallback",
        })
        virus["analytical_profile"]["per_test_signals"][0][
            "nucleic_type"
        ] = "RNA"
        result = annotate_candidate(virus, self.respiratory_virus_policy)
        self.assertFalse(result["selected_for_complete_report"])

    def test_v8_context_rank_below_floor_abstains(self):
        item = self.add_mngs_signal(
            candidate(
                family="typical_respiratory_pathogen",
                decision="review_context_needed",
                rank=8,
                count=1,
                reproducible=False,
            ),
            rank=8,
            qc_status="partial_evaluable",
        )
        annotated = [annotate_candidate(item, self.abstention_policy)]
        result = apply_best_available_fallback(
            annotated,
            self.abstention_policy,
            evidence_profile=decision_evidence_profile(annotated),
        )
        self.assertFalse(result[0]["selected_for_complete_report"])
        eligibility = result[0]["possible_pathogen_gate"]["axes"][
            "fallback_eligibility"
        ]
        self.assertFalse(eligibility["eligible"])
        self.assertIn("fallback_rank_below_evidence_floor", eligibility["failures"])

    def test_v8_high_rank_within_floor_can_fallback(self):
        item = self.add_mngs_signal(
            candidate(
                family="typical_respiratory_pathogen",
                decision="review_high_priority",
                rank=8,
                count=1,
                reproducible=False,
            ),
            rank=8,
            qc_status="partial_evaluable",
        )
        annotated = [annotate_candidate(item, self.abstention_policy)]
        result = apply_best_available_fallback(
            annotated,
            self.abstention_policy,
            evidence_profile=decision_evidence_profile(annotated),
        )
        self.assertTrue(result[0]["selected_for_complete_report"])
        self.assertEqual(FALLBACK_ROLE, result[0]["possible_reporting_role"])


if __name__ == "__main__":
    unittest.main()

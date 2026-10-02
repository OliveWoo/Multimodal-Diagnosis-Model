import json
import unittest
from pathlib import Path

from tools.promote_test_aware_clinical_high_shadow import (
    DEFAULT_POLICY,
    build_patient_promotion_context,
    candida_evidence_profile,
    promote_candidate,
    promotion_decision,
)
from tools.promote_test_aware_clinical_high_shadow_v5 import (
    promote_candidate as promote_candidate_v5,
)
from tools.promote_test_aware_clinical_high_shadow_v6 import (
    promote_candidate as promote_candidate_v6,
)
from tools.promote_test_aware_clinical_high_shadow_v7 import (
    promote_candidate as promote_candidate_v7,
)
from tools.promote_test_aware_clinical_high_shadow_v8 import (
    build_patient_promotion_context as build_patient_promotion_context_v8,
    promote_candidate as promote_candidate_v8,
)


def candidate(
    name="Pseudomonas aeruginosa",
    family="hospital_or_nonfermenter_gnb",
    *,
    decision="review_high_priority",
    rank=2,
    count=2,
    reproducible=True,
    cross=False,
    direct=None,
    host=False,
    image=False,
    source="mngs_positive_reads",
    context="lower_respiratory",
):
    item = {
        "patient_id": "1",
        "organism_name": name,
        "candidate_source": source,
        "incoming_decision": decision,
        "clinical_decision": decision,
        "clinical_level": "Level 3" if decision == "review_high_priority" else "Level 4",
        "formal_pick_allowed": decision == "picked_shadow",
        "rule_ids": [],
        "reasons": [],
        "taxonomy_profile": {
            "primary_rule_family": family,
            "taxonomic_rank": "species",
            "mapping_status": "exact_species",
        },
        "analytical_profile": {
            "best_rank_in_retained_universe": rank,
            "selected_positive_test_count": count,
            "reproducibility_axis": reproducible,
            "cross_molecule_selected": cross,
        } if source == "mngs_positive_reads" else None,
        "hospital_profile": {
            "best_direct_level": f"Level {direct}" if direct else None,
        } if source == "mngs_positive_reads" else None,
        "hospital_evidence_detail": {
            "best_hospital_level": f"Level {direct}" if direct else None,
        } if source == "hospital_only" else None,
        "host_support_flag_from_upstream": host,
        "host_evidence": ([{
            "evidence_id": "HOST-001",
            "evidence_role": "host_medication",
            "label": "verified immunosuppressive exposure",
            "assertion": "verified_active_exposure",
            "episode_evidence_status": "host_exposure_verified",
        }] if host else []),
        "image_evidence": ([{"assertion": "differential"}] if image else []),
        "source_candidate": {
            "specimen_context": context,
            "collected_time": "2024-11-14 09:30",
        },
        "linked_events": [],
    }
    return item


class TestClinicalHighPromotion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.policy = json.loads(Path(DEFAULT_POLICY).read_text(encoding="utf-8"))

    def test_reproducible_top_three_respiratory_bacterium_promotes(self):
        result = promote_candidate(candidate(), self.policy)
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual("bacterial_reproducible_top3", result["promotion_gate"]["route"])
        self.assertEqual("possible_concurrent_pathogen", result["selection_role"])
        self.assertTrue(result["selected_for_report"])

    def test_existing_pick_keeps_distinct_reporting_role(self):
        result = promote_candidate(candidate(decision="picked_shadow"), self.policy)
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual("existing_picked_pathogen", result["selection_role"])
        self.assertTrue(result["selected_for_report"])

    def test_conservative_policy_keeps_analytical_only_bacterium_high(self):
        policy = json.loads(json.dumps(self.policy))
        policy["allow_analytical_only_bacterial_promotion"] = False
        result = promote_candidate(candidate(), policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])

    def test_conservative_policy_allows_direct_supported_bacterium(self):
        policy = json.loads(json.dumps(self.policy))
        policy["allow_analytical_only_bacterial_promotion"] = False
        result = promote_candidate(candidate(direct=1), policy)
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual("bacterial_direct_support", result["promotion_gate"]["route"])

    def test_targeted_assay_bridge_requires_exact_event_aligned_cross_molecule_top3(self):
        policy = json.loads(json.dumps(self.policy))
        policy["allow_analytical_only_bacterial_promotion"] = False
        policy["require_event_aligned_direct_evidence"] = True
        policy["targeted_respiratory_assay_bridge"] = {
            "enabled": True,
            "allowed_families": ["typical_respiratory_pathogen"],
            "allowed_mapping_statuses": ["exact_species"],
            "maximum_rank": 3,
            "minimum_positive_tests": 3,
        }
        item = candidate(
            "Haemophilus influenzae", "typical_respiratory_pathogen",
            rank=2, count=3, reproducible=True, cross=True,
        )
        item["linked_events"] = [{"collected_time": "2024-11-14 09:30"}]
        item["hospital_profile"] = {
            "best_direct_level": None,
            "exact_or_alias_evidence": True,
            "hospital_evidence_detail": {
                "source_observations": [{
                    "observation_id": "ASSAY-HINF",
                    "detection_status": "detected",
                    "raw_result": "Detected; 10^4 copy/mL",
                    "reported_time": "2024-11-14 10:00",
                    "specimen_type": "BAL",
                    "target": "Haemophilus influenza",
                }],
            },
        }
        result = promote_candidate_v5(item, policy)
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual(
            "bacterial_targeted_respiratory_assay_bridge",
            result["promotion_gate"]["route"],
        )

        item["analytical_profile"]["best_rank_in_retained_universe"] = 4
        retained = promote_candidate_v5(item, policy)
        self.assertEqual("review_high_priority", retained["clinical_decision"])

    def test_v6_analytically_compelling_route_requires_two_tests_at_rpm_five(self):
        policy = json.loads(
            Path(
                "rules/test_aware_clinical_promotion_v6_analytically_compelling_rpm5.json"
            ).read_text(encoding="utf-8")
        )
        item = candidate(
            "Stenotrophomonas maltophilia",
            "hospital_or_nonfermenter_gnb",
            rank=2,
            count=3,
            reproducible=True,
            cross=True,
        )
        item["analytical_profile"]["per_test_signals"] = [
            {
                "seq_id": "DNA-1", "nucleic_type": "DNA",
                "qc_status": "partial_evaluable",
                "normalization_available": True, "rpm_total": 5.1,
            },
            {
                "seq_id": "DNA-2", "nucleic_type": "DNA",
                "qc_status": "partial_evaluable",
                "normalization_available": True, "rpm_total": 9.7,
            },
            {
                "seq_id": "RNA-1", "nucleic_type": "RNA",
                "qc_status": "partial_evaluable",
                "normalization_available": True, "rpm_total": 1.2,
            },
        ]
        promoted = promote_candidate_v6(item, policy)
        self.assertEqual("picked_shadow", promoted["clinical_decision"])
        self.assertEqual(
            "bacterial_analytically_compelling_cross_molecule",
            promoted["promotion_gate"]["route"],
        )

        item["analytical_profile"]["per_test_signals"][1]["rpm_total"] = 4.9
        retained = promote_candidate_v6(item, policy)
        self.assertEqual("review_high_priority", retained["clinical_decision"])

    def test_v6_analytical_route_rejects_genus_inherited_and_explicit_colonization(self):
        policy = json.loads(
            Path(
                "rules/test_aware_clinical_promotion_v6_analytically_compelling_rpm5.json"
            ).read_text(encoding="utf-8")
        )
        item = candidate(
            "Klebsiella pneumoniae", "typical_respiratory_pathogen",
            rank=2, count=3, reproducible=True, cross=True,
        )
        item["analytical_profile"]["per_test_signals"] = [
            {
                "seq_id": "DNA-1", "nucleic_type": "DNA",
                "qc_status": "evaluable", "normalization_available": True,
                "rpm_total": 20.0,
            },
            {
                "seq_id": "DNA-2", "nucleic_type": "DNA",
                "qc_status": "evaluable", "normalization_available": True,
                "rpm_total": 15.0,
            },
            {
                "seq_id": "RNA-1", "nucleic_type": "RNA",
                "qc_status": "evaluable", "normalization_available": True,
                "rpm_total": 10.0,
            },
        ]
        item["taxonomy_profile"]["mapping_status"] = "genus_inherited"
        retained = promote_candidate_v6(item, policy)
        self.assertEqual("review_high_priority", retained["clinical_decision"])

        item["taxonomy_profile"]["mapping_status"] = "exact_species"
        item["colonization_interpretation"] = {
            "has_explicit_colonization": True,
            "has_explicit_current_not_infection": False,
        }
        retained = promote_candidate_v6(item, policy)
        self.assertEqual("review_high_priority", retained["clinical_decision"])

    @staticmethod
    def add_cross_site_culture(item, *, organism="Corynebacterium striatum"):
        item["hospital_profile"] = {
            "exact_or_alias_evidence": True,
            "best_direct_level": None,
            "hospital_evidence_detail": {
                "module_evidence": {
                    "culture": [{
                        "specimen_type": "Blood",
                        "specimen_category": "Sterile_Site",
                        "growth_purity": "isolated_pure",
                        "evidence_observation_ids": ["CUL-EXACT"],
                    }],
                },
                "source_observations": [{
                    "observation_id": "CUL-EXACT",
                    "organism_name": organism,
                    "specimen_type": "Blood",
                    "specimen_category": "Sterile_Site",
                    "collected_time": "2024-11-13 21:37",
                    "raw_result": "Isolated",
                    "quantitation_status": "positive_detected_quantity_unknown",
                }],
            },
        }
        return item

    def test_v7_cross_site_convergence_can_promote_context_colonizer(self):
        policy = json.loads(
            Path(
                "rules/test_aware_clinical_promotion_v7_cross_site_convergence.json"
            ).read_text(encoding="utf-8")
        )
        item = candidate(
            "Corynebacterium striatum",
            "skin_airway_colonizer_prone",
            decision="review_context_needed",
            rank=1,
            count=3,
            reproducible=True,
            cross=True,
        )
        item["analytical_profile"]["per_test_signals"] = [
            {
                "seq_id": "DNA-1", "nucleic_type": "DNA",
                "qc_status": "partial_evaluable",
                "normalization_available": True, "rpm_total": 100.0,
            },
            {
                "seq_id": "DNA-2", "nucleic_type": "DNA",
                "qc_status": "partial_evaluable",
                "normalization_available": False, "rpm_total": None,
            },
            {
                "seq_id": "RNA-1", "nucleic_type": "RNA",
                "qc_status": "partial_evaluable",
                "normalization_available": False, "rpm_total": None,
            },
        ]
        item["source_candidate"]["colonization_interpretation"] = {
            "has_provisional_contaminant_label": True,
            "has_explicit_colonization": False,
            "has_explicit_current_not_infection": False,
        }
        self.add_cross_site_culture(item)
        promoted = promote_candidate_v7(item, policy)
        self.assertEqual("picked_shadow", promoted["clinical_decision"])
        self.assertEqual(
            "colonizer_prone_cross_site_convergence",
            promoted["promotion_gate"]["route"],
        )
        self.assertIn(
            "PROMO-CS5-PROVISIONAL-CONTAMINANT-LABEL-RETAINED-AS-CAUTION",
            promoted["rule_ids"],
        )

    def test_v7_cross_site_convergence_requires_exact_culture_and_no_explicit_counterevidence(self):
        policy = json.loads(
            Path(
                "rules/test_aware_clinical_promotion_v7_cross_site_convergence.json"
            ).read_text(encoding="utf-8")
        )
        item = candidate(
            "Corynebacterium striatum",
            "skin_airway_colonizer_prone",
            decision="review_context_needed",
            rank=1,
            count=3,
            reproducible=True,
            cross=True,
        )
        item["analytical_profile"]["per_test_signals"] = [{
            "seq_id": "DNA-1", "nucleic_type": "DNA",
            "qc_status": "evaluable", "normalization_available": True,
            "rpm_total": 100.0,
        }]
        self.add_cross_site_culture(item, organism="Corynebacterium jeikeium")
        retained = promote_candidate_v7(item, policy)
        self.assertEqual("review_context_needed", retained["clinical_decision"])

        self.add_cross_site_culture(item)
        item["source_candidate"]["colonization_interpretation"] = {
            "has_explicit_colonization": True,
            "has_explicit_current_not_infection": False,
        }
        retained = promote_candidate_v7(item, policy)
        self.assertEqual("review_context_needed", retained["clinical_decision"])

    def test_v8_exact_hpiv3_represents_event_aligned_generic_panel_target(self):
        policy = json.loads(
            Path(
                "rules/test_aware_clinical_promotion_v8_respiratory_group_convergence.json"
            ).read_text(encoding="utf-8")
        )
        exact = candidate(
            "Human respirovirus 3",
            "other_respiratory_virus",
            rank=2,
            count=1,
            reproducible=False,
        )
        group = candidate(
            "Parainfluenza Virus",
            "other_respiratory_virus",
            source="hospital_only",
            direct=2,
        )
        group["hospital_evidence_detail"]["source_observations"] = [{
            "observation_id": "ASSAY-024",
            "test_type": "filmarray",
            "specimen_type": "BAL",
            "reported_time": "2024-11-14 10:00",
            "target": "Parainfluenza Virus",
            "raw_result": "Detected",
            "detection_status": "detected",
        }]
        context = build_patient_promotion_context_v8([exact, group], policy)
        promoted = promote_candidate_v8(exact, policy, context)
        consolidated = promote_candidate_v8(group, policy, context)
        self.assertEqual("picked_shadow", promoted["clinical_decision"])
        self.assertEqual(
            "respiratory_virus_exact_mngs_plus_group_assay",
            promoted["promotion_gate"]["route"],
        )
        self.assertEqual(
            "review_context_needed", consolidated["clinical_decision"]
        )
        self.assertEqual(
            "broad_group_assay_consolidated_under_exact_mngs_member",
            consolidated["promotion_gate"]["route"],
        )

    def test_v8_group_link_does_not_cross_event_window(self):
        policy = json.loads(
            Path(
                "rules/test_aware_clinical_promotion_v8_respiratory_group_convergence.json"
            ).read_text(encoding="utf-8")
        )
        exact = candidate(
            "Human respirovirus 3", "other_respiratory_virus", rank=2, count=1
        )
        group = candidate(
            "Parainfluenza Virus",
            "other_respiratory_virus",
            source="hospital_only",
            direct=2,
        )
        group["hospital_evidence_detail"]["source_observations"] = [{
            "observation_id": "ASSAY-024",
            "test_type": "filmarray",
            "specimen_type": "BAL",
            "reported_time": "2024-10-01 10:00",
            "target": "Parainfluenza Virus",
            "raw_result": "Detected",
            "detection_status": "detected",
        }]
        context = build_patient_promotion_context_v8([exact, group], policy)
        self.assertEqual([], context["respiratory_virus_group_convergence_links"])
        retained = promote_candidate_v8(exact, policy, context)
        self.assertEqual("review_high_priority", retained["clinical_decision"])

    def test_single_test_top_one_without_direct_support_stays_high(self):
        result = promote_candidate(candidate(rank=1, count=1, reproducible=False), self.policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])

    def test_candida_in_respiratory_specimen_stays_high(self):
        result = promote_candidate(
            candidate("Candida tropicalis", "candida_or_yeast", direct=1), self.policy
        )
        self.assertEqual("review_high_priority", result["clinical_decision"])
        self.assertIn("no event-aligned positive", result["promotion_gate"]["blockers"][0])

    @staticmethod
    def add_candida_observation(
        item, *, organism="Candida dubliniensis", specimen="Blood",
        category="Sterile_Site", collected="2024-11-13 21:37", result="Isolated",
    ):
        detail = {
            "best_hospital_level": "Level 2",
            "source_observations": [{
                "observation_id": "CUL-001",
                "organism_name": organism,
                "specimen_type": specimen,
                "specimen_category": category,
                "collected_time": collected,
                "raw_result": result,
                "quantitation_status": (
                    "positive_detected_quantity_unknown" if result == "Isolated" else ""
                ),
            }],
        }
        if item["candidate_source"] == "hospital_only":
            item["hospital_evidence_detail"] = detail
        else:
            item["hospital_profile"]["hospital_evidence_detail"] = detail
        return item

    def test_exact_candida_blood_positive_within_event_promotes(self):
        item = self.add_candida_observation(
            candidate("Candida dubliniensis", "candida_or_yeast", direct=2)
        )
        context = build_patient_promotion_context([item], self.policy)
        result = promote_candidate(item, self.policy, context)
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual(
            "candida_exact_species_invasive_event", result["promotion_gate"]["route"]
        )

    def test_negative_candida_blood_record_does_not_promote(self):
        item = self.add_candida_observation(
            candidate("Candida dubliniensis", "candida_or_yeast", direct=2),
            result="No growth",
        )
        result = promote_candidate(item, self.policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])
        self.assertEqual(
            "negative",
            result["promotion_gate"]["axes"]["candida_invasive_evidence"]["rows"][0]["assertion"],
        )

    def test_missing_candida_timing_is_not_negative_or_removed(self):
        item = self.add_candida_observation(
            candidate("Candida dubliniensis", "candida_or_yeast", direct=2), collected=None
        )
        profile = candida_evidence_profile(item, self.policy)
        result = promote_candidate(item, self.policy)
        self.assertEqual(1, profile["positive_invasive_timing_unknown_count"])
        self.assertFalse(profile["missing_evidence_is_negative"])
        self.assertEqual("review_high_priority", result["clinical_decision"])

    def test_respiratory_candida_does_not_use_invasive_route(self):
        item = self.add_candida_observation(
            candidate("Candida albicans", "candida_or_yeast", direct=2),
            organism="Candida albicans", specimen="Lower BAL",
            category="Lower_Respiratory",
        )
        result = promote_candidate(item, self.policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])

    def test_generic_yeast_can_represent_unresolved_invasive_event(self):
        item = candidate(
            "Yeast/Candida evidence", "unmapped_or_uncertain",
            source="hospital_only", direct=2,
        )
        item["taxonomy_profile"]["taxonomic_rank"] = "genus_or_group_label"
        item = self.add_candida_observation(item, organism="Yeast")
        context = build_patient_promotion_context([item], self.policy)
        result = promote_candidate(item, self.policy, context)
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual(
            "candida_group_invasive_event_species_unresolved",
            result["promotion_gate"]["route"],
        )

    def test_generic_yeast_stays_supporting_when_exact_species_exists(self):
        exact = self.add_candida_observation(
            candidate("Candida dubliniensis", "candida_or_yeast", direct=2)
        )
        generic = candidate(
            "Yeast/Candida evidence", "unmapped_or_uncertain",
            source="hospital_only", direct=2,
        )
        generic["taxonomy_profile"]["taxonomic_rank"] = "genus_or_group_label"
        generic = self.add_candida_observation(generic, organism="Yeast")
        context = build_patient_promotion_context([exact, generic], self.policy)
        result = promote_candidate(generic, self.policy, context)
        self.assertEqual("review_high_priority", result["clinical_decision"])
        self.assertIn("exact Candida species already represents", result["promotion_gate"]["blockers"][0])

    def test_pending_bacterial_identification_is_not_a_candida_group(self):
        item = candidate(
            "G(+) Bacilli, identification to follow.", "unmapped_or_uncertain",
            source="hospital_only", direct=2,
        )
        item["taxonomy_profile"]["taxonomic_rank"] = "genus_or_group_label"
        item = self.add_candida_observation(
            item, organism="G(+) Bacilli, identification to follow."
        )
        context = build_patient_promotion_context([item], self.policy)
        result = promote_candidate(item, self.policy, context)
        self.assertEqual("review_high_priority", result["clinical_decision"])
        self.assertIsNone(
            result["promotion_gate"]["axes"].get("candida_invasive_evidence")
        )

    def test_mold_needs_direct_l1_l2_and_reproducibility(self):
        promoted = promote_candidate(
            candidate("Aspergillus flavus", "mold_or_opportunistic_fungus", rank=5, direct=1),
            self.policy,
        )
        retained = promote_candidate(
            candidate("Aspergillus fumigatus", "mold_or_opportunistic_fungus",
                      rank=5, count=1, reproducible=False, direct=1),
            self.policy,
        )
        self.assertEqual("picked_shadow", promoted["clinical_decision"])
        self.assertEqual("review_high_priority", retained["clinical_decision"])

    def test_pjp_host_route_requires_reproducible_top_three(self):
        result = promote_candidate(
            candidate("Pneumocystis jirovecii", "high_consequence_opportunistic", host=True),
            self.policy,
        )
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual("pjp_host_risk_triangulation", result["promotion_gate"]["route"])

    def test_pjp_broad_host_flag_without_verified_evidence_stays_high(self):
        item = candidate(
            "Pneumocystis jirovecii", "high_consequence_opportunistic", host=True
        )
        item["host_evidence"] = [{
            "evidence_id": "HOST-UNVERIFIED",
            "evidence_role": "host_medication",
            "label": "prednisolone",
            "assertion": "recorded_exposure_not_administration_verified",
        }]
        result = promote_candidate(item, self.policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])
        self.assertFalse(
            result["promotion_gate"]["axes"]["host_support_eligible_for_promotion"]
        )

    def test_consumed_history_cannot_be_reused_as_pjp_host_axis(self):
        item = candidate(
            "Pneumocystis jirovecii", "high_consequence_opportunistic", host=True
        )
        item["consumed_evidence_ids"] = ["P1:HISTCLAIM:host"]
        result = promote_candidate(item, self.policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])
        self.assertFalse(
            result["promotion_gate"]["axes"]["host_support_eligible_for_promotion"]
        )
        self.assertIn("PROMO-S3-CONSUMED-HISTORY-NOT-REUSED", result["rule_ids"])

    def test_pjp_image_route_requires_cross_molecule_reproducibility(self):
        promoted = promote_candidate(
            candidate("Pneumocystis jirovecii", "high_consequence_opportunistic",
                      rank=4, count=3, reproducible=True, cross=True, image=True),
            self.policy,
        )
        retained = promote_candidate(
            candidate("Pneumocystis jirovecii", "high_consequence_opportunistic",
                      rank=4, count=2, reproducible=True, cross=False, image=True),
            self.policy,
        )
        self.assertEqual("picked_shadow", promoted["clinical_decision"])
        self.assertEqual("review_high_priority", retained["clinical_decision"])

    def test_generic_complex_cannot_promote(self):
        item = candidate("Acinetobacter calcoaceticus-baumannii complex")
        item["taxonomy_profile"]["taxonomic_rank"] = "genus_or_group_label"
        result = promote_candidate(item, self.policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])
        self.assertFalse(result["promotion_gate"]["axes"]["precise_identity"])

    def test_direct_respiratory_assay_target_can_promote(self):
        item = candidate(
            "respiratory syncytial virus", "other_respiratory_virus",
            source="hospital_only", direct=2,
        )
        result = promote_candidate(item, self.policy)
        self.assertEqual("picked_shadow", result["clinical_decision"])
        self.assertEqual("respiratory_assay_target", result["promotion_gate"]["route"])

    def test_event_aligned_policy_rejects_direct_evidence_with_unknown_time(self):
        policy = json.loads(json.dumps(self.policy))
        policy["require_event_aligned_direct_evidence"] = True
        item = candidate(
            "respiratory syncytial virus", "other_respiratory_virus",
            source="hospital_only", direct=2,
        )
        item["hospital_evidence_detail"]["source_observations"] = [{
            "observation_id": "ASSAY-UNKNOWN-TIME",
            "detection_status": "detected",
            "raw_result": "Detected",
            "specimen_type": "BAL",
        }]
        result = promote_candidate(item, policy)
        self.assertEqual("review_high_priority", result["clinical_decision"])

    def test_event_aligned_policy_accepts_report_time_within_window(self):
        policy = json.loads(json.dumps(self.policy))
        policy["require_event_aligned_direct_evidence"] = True
        item = candidate(
            "respiratory syncytial virus", "other_respiratory_virus",
            source="hospital_only", direct=2,
        )
        item["linked_events"] = [{"collected_time": "2024-11-14 09:30"}]
        item["hospital_evidence_detail"]["source_observations"] = [{
            "observation_id": "ASSAY-TIMED",
            "detection_status": "detected",
            "raw_result": "Detected",
            "reported_time": "2024-11-14 10:00",
            "specimen_type": "BAL",
        }]
        result = promote_candidate(item, policy)
        self.assertEqual("picked_shadow", result["clinical_decision"])

    def test_non_high_candidate_is_not_reclassified(self):
        item = candidate(decision="review_context_needed")
        decision, gate, _, _ = promotion_decision(item, self.policy)
        self.assertEqual("review_context_needed", decision)
        self.assertFalse(gate["evaluated"])
        result = promote_candidate(item, self.policy)
        self.assertEqual("context_review_candidate", result["selection_role"])
        self.assertFalse(result["selected_for_report"])


if __name__ == "__main__":
    unittest.main()

# All fixtures below are synthetic, not patient exports.
import copy
import json
import tempfile
import unittest
from pathlib import Path

from tools.bind_viral_evidence_shadow import (
    POLICY_PATH, align, bind_candidate, claim_from_fact, result_status,
    row_facts, target_matches, current_activity_assertion, deduplicate_claims,
    reconcile_frozen_observations,
    run,
)
from tools.viral_evidence_axes_shadow import DEFAULT_POLICY


class ViralEvidenceBindingTests(unittest.TestCase):
    def setUp(self):
        self.policy = json.loads(POLICY_PATH.read_text(encoding="utf-8"))
        self.axes = json.loads(DEFAULT_POLICY.read_text(encoding="utf-8"))
        self.event = {"event_id": "P1:mngs:1", "case_review_id": "CASE-A",
                      "collected_time": "2024-01-10 13:00", "specimen_site": "blood"}

    def candidate(self):
        return {"patient_id": "1", "organism_name": "Human betaherpesvirus 5",
                "taxonomy_profile": {"mapping_status": "exact_species"},
                "linked_events": [self.event], "source_candidate": {"case_review_id": "CASE-A"},
                "analytical_profile": {"per_test_signals": [
                    {"reads": 12, "qc_status": "evaluable", "nucleic_type": "DNA", "seq_id": "S1"}]},
                "final_reporting_tier": "Possible", "selected_for_complete_report": True}

    def fact(self, **overrides):
        return {"evidence_id": "clinical-1", "patient_id": "1", "module": "molecular_microbiology",
                "fact_kind": "molecular", "target_matches": {"humanbetaherpesvirus5": "exact_or_alias"},
                "collected_time": "2024-01-11", "specimen": "Blood", "result_status": "positive",
                "result": "Detected", "assay_method": "PCR", "blockers": [],
                "temporal_usable": True, "source_refs": [{"file": "source.json", "sha256": "hash", "json_path": "[0]"}],
                **overrides}

    def bind(self, facts, candidate=None):
        return bind_candidate(candidate or self.candidate(), facts, run="A11",
            provenance={"file": "frozen.json", "sha256": "hash", "json_path": "[0]"},
            policy=self.policy, axes_policy=self.axes)

    def test_cmv_vl_is_cmv_and_never_hsv(self):
        self.assertEqual(target_matches("CMV-VL", self.policy), {"humanbetaherpesvirus5": "exact_or_alias"})

    def test_generic_hsv_does_not_identify_type_one(self):
        matches = target_matches("HSV PCR", self.policy)
        self.assertEqual(matches["humanalphaherpesvirus1"], "group_only")

    def test_combined_hsv_panel_does_not_identify_either_type(self):
        matches = target_matches("HSV-1/2 PCR", self.policy)
        self.assertEqual(set(matches.values()), {"group_only"})

    def test_mixed_positive_negative_result_requires_review(self):
        self.assertEqual(result_status("Detected / Not Detected"), "conflicting_or_composite_result")

    def test_negative_is_not_positive(self):
        self.assertEqual(result_status("Target Not Detected"), "negative")
        result = self.bind([self.fact(result_status="negative", result="Target Not Detected")])
        self.assertEqual(result["axes"]["detection"], "D1")
        self.assertEqual(len(result["claims"]), 1)
        self.assertIn("result_is_negative", result["evidence_links"][0]["blockers"])

    def test_calendar_window_and_report_time_do_not_replace_collection(self):
        self.assertEqual(align("2024-01-12 23:59", "2024-01-10 00:00")["status"], "within_window")
        result = self.bind([self.fact(collected_time="2024-01-01", reported_time="2024-01-11")])
        self.assertEqual(len(result["claims"]), 1)
        self.assertEqual(result["evidence_links"][0]["event_alignment"]["status"], "outside_window")

    def test_report_only_time_stays_unknown(self):
        result = self.bind([self.fact(collected_time=None, reported_time="2024-01-11")])
        self.assertEqual(result["axes"]["detection"], "D1")
        self.assertIn("unknown", result["evidence_links"][0]["blockers"])

    def test_unresolved_source_never_auto_released(self):
        fact = self.fact(blockers=["source_record_requires_review"])
        result = self.bind([fact])
        self.assertEqual(result["axes"]["detection"], "D1")
        self.assertFalse(result["evidence_links"][0]["accepted_claim_type"])

    def test_serology_does_not_prove_activity(self):
        result = self.bind([self.fact(fact_kind="serology", result="Reactive")])
        self.assertEqual((result["axes"]["detection"], result["axes"]["activity"]), ("D1", "A0"))

    def test_blood_pcr_supports_detection_not_pneumonia(self):
        result = self.bind([self.fact()])
        self.assertEqual((result["axes"]["detection"], result["axes"]["organ_support"], result["axes"]["activity"]),
                         ("D2", "O0", "A0"))

    def test_unknown_method_does_not_invent_independent_assay(self):
        result = self.bind([self.fact(assay_method=None)])
        self.assertEqual(result["axes"]["detection"], "D1")
        self.assertEqual(result["claims"][-1]["independence_class"], "not_established")

    def test_different_compartment_can_corroborate_detection_only(self):
        result = self.bind([self.fact(assay_method=None, specimen="Sputum (Endo)")])
        self.assertEqual((result["axes"]["detection"], result["axes"]["organ_support"]), ("D2", "O0"))

    def test_technical_branches_count_as_one_biological_evidence_group(self):
        c = self.candidate()
        c["analytical_profile"]["per_test_signals"] *= 8
        result = self.bind([], c)
        self.assertEqual(result["axes"]["detection"], "D1")
        self.assertEqual(len(result["claims"]), 1)
        self.assertEqual(result["axes"]["unique_supporting_evidence_count"], 1)

    def test_pulmonary_context_does_not_create_organ_disease(self):
        fact = self.fact(fact_kind="index_symptom", module="mngs_case_symptom", target_matches={},
                         pulmonary_context_mentioned=True, source_text="Pneumonia", case_review_id="CASE-A")
        result = self.bind([fact])
        self.assertEqual((result["axes"]["organ_support"], result["axes"]["shadow_tier"]), ("O1", "Context"))

    def test_other_case_context_is_not_pooled(self):
        fact = self.fact(fact_kind="index_symptom", target_matches={}, pulmonary_context_mentioned=True,
                         source_text="Pneumonia", case_review_id="CASE-B")
        result = self.bind([fact])
        self.assertEqual(result["axes"]["organ_support"], "O0")
        self.assertIn("different_mngs_case", result["evidence_links"][0]["blockers"])

    def test_numeric_measurement_keeps_activity_unknown(self):
        result = self.bind([self.fact(result_status="numeric_measurement", result="702000", unit="IU/mL")])
        self.assertEqual(result["axes"]["activity"], "A0")

    def test_literal_unit_recovered_without_creating_assay_method(self):
        row = {"target": "CMV-VL", "test": "Non-panel Molecular Test", "result": "123",
               "source": {"raw_result": "(IU/mL)123"}}
        facts = row_facts("SYNTHETIC", "901", "molecular_microbiology", row, 0, POLICY_PATH, self.policy)
        self.assertEqual(facts[0]["unit"], "IU/mL")
        self.assertIsNone(facts[0]["assay_method"])

    def test_separate_numeric_value_is_preserved_without_unit_inference(self):
        row = {"target": "CMV", "test": "Non-panel Molecular Test", "result": "Detected", "value": "777"}
        fact = row_facts("SYNTHETIC", "902", "molecular_microbiology", row, 0, POLICY_PATH, self.policy)[0]
        self.assertEqual(fact["numeric_value"], 777)
        self.assertEqual(fact["numeric_source_field"], "value")
        self.assertIsNone(fact["unit"])
        self.assertEqual(fact["result_status"], "positive")

    def test_missing_specimen_is_not_filled_with_mngs_specimen(self):
        result = self.bind([self.fact(specimen=None)])
        self.assertEqual(result["axes"]["detection"], "D1")
        self.assertIn("specimen_missing", result["evidence_links"][0]["blockers"])

    def test_current_production_fields_are_not_mutated(self):
        c = self.candidate()
        original = copy.deepcopy(c)
        result = self.bind([self.fact()], c)
        self.assertEqual(c, original)
        self.assertEqual(result["current_tier"], "Possible")
        self.assertFalse(result["production_changed"])

    def test_dated_explicit_reactivation_is_an_assertion_not_proven_disease(self):
        fact = self.fact(fact_kind="narrative", module="admission_diagnosis", source_text="Current CMV reactivation")
        result = self.bind([fact])
        self.assertEqual((result["axes"]["activity"], result["axes"]["organ_support"]), ("A1", "O0"))
        self.assertFalse(result["claims"][-1]["human_clinical_adjudication"])
        self.assertEqual(result["claims"][-1]["specimen"], "not_applicable_clinical_context")

    def test_suspected_or_other_virus_reactivation_does_not_leak(self):
        for statement in ["Possible CMV reactivation", "No CMV reactivation", "CMV infection, HSV-1 reactivation"]:
            self.assertFalse(current_activity_assertion(statement, "humanbetaherpesvirus5", self.policy))

    def test_deduplication_keeps_all_pointers_without_counting_copies(self):
        claim = self.bind([self.fact()])["claims"][-1]
        duplicate = copy.deepcopy(claim)
        duplicate["source_refs"] = [{"file": "copy.json", "json_path": "[4]", "sha256": "copyhash"}]
        deduped = deduplicate_claims([claim, duplicate])
        self.assertEqual(len(deduped), 1)
        self.assertEqual(len(deduped[0]["source_refs"]), 2)

    def test_preexisting_hospital_observation_is_reconciled(self):
        c = self.candidate()
        c["hospital_profile"] = {"hospital_evidence_detail": {"organism_name": "CMV", "source_observations": [
            {"observation_id": "MOL-001", "target": "CMV", "collected_time": "2024-01-11", "specimen_type": "Blood", "raw_result": "Detected"}]}}
        result = reconcile_frozen_observations(c, [self.fact()], self.policy)
        self.assertEqual((result["requested"], result["resolved"]), (1, 1))

    def test_qpcr_numeric_detection_does_not_prove_replication(self):
        result = self.bind([self.fact(result_status="numeric_measurement", result="120", numeric_value=120,
                                      target="CMV-VL", unit="IU/mL", assay_method="qPCR")])
        self.assertEqual((result["axes"]["detection"], result["axes"]["activity"], result["axes"]["organ_support"]),
                         ("D2", "A0", "O0"))

    def test_below_limit_numeric_value_does_not_become_positive(self):
        result = self.bind([self.fact(result_status="below_limit", result="<34.5", numeric_value=34.5,
                                      numeric_comparator="<", target="CMV-VL", unit="IU/mL", assay_method="qPCR")])
        self.assertEqual(len(result["claims"]), 1)
        self.assertEqual(result["axes"]["detection"], "D1")

    def test_complete_binding_run_preserves_candidates_and_source_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            input_root, output_root = base / "input", base / "output"
            patients, inventory = base / "patients", base / "inventory"
            folder = patients / "NGS_patient_1_json"
            folder.mkdir(parents=True)
            (input_root / "review").mkdir(parents=True)
            (input_root / "A11").mkdir()
            (inventory / "patients").mkdir(parents=True)
            def save(path, value):
                path.write_text(json.dumps(value), encoding="utf-8")
            candidate = self.candidate()
            save(input_root / "review/all_affected_evidence_packets.json", [
                {"cohort": "A11", "patient_id": "1", "organism_name": candidate["organism_name"],
                 "current_candidate_evidence_unchanged": candidate}])
            save(input_root / "A11/run_manifest.json", {"input": {"patient_root": str(patients), "inventory_root": str(inventory)}})
            source = folder / "NGS_patient_1_molecular_microbiology.json"
            save(source, [{"test": "CMV PCR", "target": "CMV", "sample": "Blood",
                           "collected_time": "2024-01-11", "result": "Detected"}])
            original = source.read_bytes()
            save(inventory / "patients/NGS_patient_1_multi_assay_mngs_evidence.json", {"cases": [
                {"case_review_id": "CASE-A", "collected_time": "2024-01-10 13:00", "specimen_site": "blood",
                 "case_review": {"symptom_description": "Pneumonia", "clinical_result": "MUST_NOT_CONSUME"}}]})
            summary = run(input_root, output_root)
            self.assertEqual((summary["candidate_count"], summary["candidate_loss"]), (1, 0))
            self.assertFalse(summary["benchmark_answers_loaded"])
            self.assertFalse(summary["production_decisions_changed"])
            self.assertEqual(source.read_bytes(), original)
            self.assertNotIn("MUST_NOT_CONSUME", (output_root / "source_facts.json").read_text(encoding="utf-8"))
            packet = json.loads((output_root / "bound_candidate_packets.json").read_text(encoding="utf-8"))[0]
            self.assertEqual((packet["axes"]["detection"], packet["axes"]["organ_support"]), ("D2", "O1"))
            with self.assertRaises(FileExistsError):
                run(input_root, output_root)


if __name__ == "__main__":
    unittest.main()

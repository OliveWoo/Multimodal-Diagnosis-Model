# All fixtures below are synthetic, not patient exports.
import copy
import unittest

from tools.viral_evidence_axes_shadow import interpret


class ViralEvidenceAxesTest(unittest.TestCase):
    def candidate(self):
        return {"organism_name": "Human betaherpesvirus 5",
                "taxonomy_profile": {"mapping_status": "exact_species"},
                "analytical_profile": {"per_test_signals": [
                    {"reads": 10, "rpm_total": 5, "rank_in_retained_universe": 1,
                     "nucleic_type": "DNA", "qc_status": "evaluable"}]}}

    def claim(self, kind, **extra):
        return {"evidence_id": "pcr-1", "organism_key": "humanbetaherpesvirus5",
                "episode_id": "P1:E1", "specimen": "blood", "collected_time": "2024-01-01",
                "source_file": "source.json", "source_row": "1", "review_complete": True,
                "event_aligned": True, "claim_type": kind, **extra}

    def test_technical_repeat_molecule_reads_and_rank_cannot_create_disease(self):
        c = self.candidate()
        before = interpret(c, episode_id="P1:E1")
        for _ in range(12):
            row = copy.deepcopy(c["analytical_profile"]["per_test_signals"][0])
            row.update(reads=1000000, nucleic_type="RNA", rpm_total=9999)
            c["analytical_profile"]["per_test_signals"].append(row)
        after = interpret(c, episode_id="P1:E1")
        self.assertEqual([(x["detection"], x["activity"], x["organ_support"], x["shadow_tier"])
                          for x in [before, after]], [("D1", "A0", "O0", "Context")] * 2)

    def test_blood_activity_does_not_prove_pneumonia(self):
        c = self.claim("validated_activity_interpretation")
        result = interpret(self.candidate(), episode_id="P1:E1", claims=[c])
        self.assertEqual((result["activity"], result["organ_support"], result["shadow_tier"]),
                         ("A2", "O0", "Context"))

    def test_proven_lung_claim_is_not_blocked_by_single_dna_test(self):
        c = self.claim("reviewed_proven_pulmonary_disease", organ="lung",
                       specimen="lung tissue", specimen_context="lower_respiratory_or_lung_tissue",
                       attribution_basis="lung_specific_pathology",
                       clinical_interpretation="source-reviewed tissue IHC and compatible disease")
        self.assertEqual(interpret(self.candidate(), episode_id="P1:E1", claims=[c])["shadow_tier"], "Picked")

    def test_wrong_episode_or_blood_organ_claim_is_rejected(self):
        c = self.claim("reviewed_proven_pulmonary_disease", episode_id="P1:E2", organ="lung")
        result = interpret(self.candidate(), episode_id="P1:E1", claims=[c])
        self.assertEqual(result["shadow_tier"], "Context")
        self.assertTrue(result["rejected_claims"])

    def test_missing_confirmation_remains_unknown(self):
        result = interpret(self.candidate(), episode_id="P1:E1")
        self.assertFalse(result["clinical_absence_inferred_from_missing_claims"])
        self.assertEqual(result["activity"], "A0")

    def test_same_evidence_is_not_two_independent_supports(self):
        claims = [self.claim("detection", independence_class="independent_clinical_assay"),
                  self.claim("validated_activity_interpretation")]
        result = interpret(self.candidate(), episode_id="P1:E1", claims=claims)
        self.assertEqual(result["unique_supporting_evidence_count"], 1)
        self.assertEqual(result["organ_support"], "O0")

    def test_one_hospital_detection_alone_is_not_corroborated(self):
        c = self.candidate()
        c["analytical_profile"]["per_test_signals"] = []
        claim = self.claim("detection", independence_class="independent_clinical_assay")
        result = interpret(c, episode_id="P1:E1", claims=[claim])
        self.assertEqual(result["detection"], "D1")

    def test_reviewed_virus_specific_attribution_does_not_require_lung_sampling(self):
        c = self.claim("reviewed_probable_pulmonary_disease", organ="lung", specimen="lesion",
                       attribution_basis="virus_specific_clinical_definition",
                       clinical_interpretation="source-reviewed virus-specific syndrome and lung attribution")
        result = interpret(self.candidate(), episode_id="P1:E1", claims=[c])
        self.assertEqual(result["shadow_tier"], "Possible")


if __name__ == "__main__":
    unittest.main()

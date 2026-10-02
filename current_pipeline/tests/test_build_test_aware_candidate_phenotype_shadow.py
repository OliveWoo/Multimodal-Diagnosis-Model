import json
import tempfile
import unittest
from pathlib import Path

from tools.build_test_aware_candidate_phenotype_shadow import enrich_contexts, run


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


class TestTestAwareCandidatePhenotypeShadow(unittest.TestCase):
    def test_enriches_mngs_and_hospital_only_event_scope(self):
        shadow = {"candidate_contexts": [
            {"organism_name": "Pneumocystis jirovecii", "canonical_key": "pneumocystisjirovecii"},
            {"organism_name": "Staphylococcus aureus", "canonical_key": "staphylococcusaureus"},
        ]}
        candidates = [
            {"organism_name": "Pneumocystis jirovecii", "organism_key": "pneumocystisjirovecii",
             "decision": "review_context_needed", "integrated_level": "Level 4",
             "case_decision_refs": [{"case_review_id": "CASE-B"}]},
            {"organism_name": "Staphylococcus aureus", "organism_key": "staphylococcusaureus",
             "decision": "review_high_priority", "integrated_level": "Level 2",
             "evidence_source": "hospital_only"},
        ]
        result = enrich_contexts(shadow, candidates)
        pjp, aureus = result["candidate_contexts"]
        self.assertEqual(pjp["candidate_event_refs"], ["CASE-B"])
        self.assertEqual(pjp["event_scope"], "listed_multi_assay_cases")
        self.assertEqual(aureus["candidate_event_refs"], [])
        self.assertEqual(aureus["event_scope"], "patient_episode_all_events")
        self.assertEqual(aureus["phenotype_deterministic_effect"], "none_provisional_link")

    def test_enrich_uses_central_alias_normalization(self):
        shadow = {"candidate_contexts": [{
            "organism_name": "Nakaseomyces glabratus",
            "canonical_key": "candidaglabrata",
        }]}
        candidates = [{
            "organism_name": "Nakaseomyces glabratus",
            "organism_key": "nakaseomycesglabratus",
            "decision": "review_context_needed",
            "case_review_id": "CASE-1",
        }]
        result = enrich_contexts(shadow, candidates)
        self.assertEqual(result["candidate_contexts"][0]["candidate_event_refs"], ["CASE-1"])

    def test_run_preserves_forwarded_candidate_count_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            scorer = root / "scorer"
            phenotype = root / "phenotype"
            output = root / "output"
            write(scorer / "patient_outputs" / "NGS_patient_8_test_aware_deterministic_shadow.json", {
                "patient_id": "8",
                "clinical_scorer_forward": [{
                    "organism_name": "Staphylococcus haemolyticus",
                    "organism_key": "staphylococcushaemolyticus",
                    "decision": "review_context_needed",
                    "integrated_level": "Level 4",
                    "case_decision_refs": [{"case_review_id": "CASE-8"}],
                }],
            })
            summary = run(scorer, phenotype, output)
            self.assertEqual(summary["counts"]["candidates"], 1)
            saved = json.loads((output / "patient_8_candidate_phenotype_shadow_v1.json").read_text())
            self.assertEqual(saved["candidate_contexts"][0]["candidate_event_refs"], ["CASE-8"])
            self.assertEqual(saved["linkage"]["status"], "phenotype_file_missing")
            with self.assertRaises(FileExistsError):
                run(scorer, phenotype, output)


if __name__ == "__main__":
    unittest.main()

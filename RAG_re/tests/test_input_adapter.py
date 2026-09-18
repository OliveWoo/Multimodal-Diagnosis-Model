import json
import tempfile
import unittest
from pathlib import Path

from rag_re.input_adapter import (
    candidate_literature_eligible,
    candidate_name,
    load_bundle,
    patient_id_from,
    select_candidate_pool,
    source_context,
)


class InputAdapterTests(unittest.TestCase):
    def test_nested_wrapper_patient_id_is_preserved(self):
        payload = {
            "result": {
                "patient_id": "nested-7",
                "deterministic_max": {"pathogen_candidates": []},
            }
        }
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "unknown.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            root, bundle = load_bundle(path)
        self.assertEqual(patient_id_from(path, root, bundle), "nested-7")

    def test_ambiguous_sibling_wrappers_are_rejected(self):
        payload = {
            "agent_bundle": {"pathogen_candidates": []},
            "infection_bundle": {"pathogen_candidates": []},
        }
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "ambiguous.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "ambiguous infection bundles"):
                load_bundle(path)

    def test_direct_plus_nested_bundle_and_conflicting_ids_are_rejected(self):
        payload = {
            "pathogen_candidates": [],
            "agent_bundle": {"pathogen_candidates": []},
        }
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "ambiguous.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "ambiguous infection bundles"):
                load_bundle(path)

        root = {"patient_id": "root", "agent_bundle": {"pathogen_candidates": []}}
        bundle = root["agent_bundle"]
        bundle["patient_id"] = "bundle"
        with self.assertRaisesRegex(ValueError, "conflicting patient/case IDs"):
            patient_id_from(Path("unknown.json"), root, bundle)

    def test_clinical_site_is_separate_from_candidate_specimen_site(self):
        context = source_context(
            {"dominant_source": "Respiratory"},
            {"specimen_class": "S1_sterile_systemic"},
        )
        self.assertEqual(context["site_code"], "lower_respiratory")
        self.assertEqual(context["specimen_site_code"], "bloodstream")
        self.assertEqual(context["clinical_site_source"], "dominant_source")

    def test_formal_merge_uses_only_picked_and_required_review_tiers(self):
        picked = {"organism_name": "Picked species", "basis_level": "Level 2"}
        deterministic = {
            "best_available_summary": {"picked_pathogens": [picked]},
            "pathogen_candidates": [
                {"organism_name": "Picked species", "integrated_causative_level": "Level 2"},
                {"organism_name": "High species", "integrated_causative_level": "Level 4"},
                {"organism_name": "Low species", "integrated_causative_level": "Level 5"},
            ],
        }
        root = {
            "deterministic_max": deterministic,
            "llm_missed_candidate_review": {
                "review_high_priority": [{"organism_name": "High species"}],
                "review_context_needed": [
                    {
                        "organism_name": "Snapshot-only species",
                        "evidence_snapshot": {
                            "rank_priority": 1,
                            "reads_percentile": 0.9,
                            "specimen_class": "S2_lower_respiratory",
                            "support_modules": ["host"],
                        },
                    }
                ],
                "review_low_specificity": [{"organism_name": "Low species"}],
                "omitted_with_reason": [{"organism_name": "Omitted species"}],
            },
        }
        selected, meta = select_candidate_pool(
            root,
            deterministic,
            {"merge_review_tiers": ["review_high_priority", "review_context_needed"]},
        )
        self.assertEqual(
            [candidate_name(row) for row in selected],
            ["Picked species", "High species", "Snapshot-only species"],
        )
        self.assertFalse(candidate_literature_eligible(selected[0]))
        self.assertTrue(candidate_literature_eligible(selected[1]))
        self.assertTrue(candidate_literature_eligible(selected[2]))
        self.assertEqual(meta["synthetic_review_candidate_count"], 1)
        self.assertEqual(meta["excluded_review_tier_counts"]["review_low_specificity"], 1)
        self.assertEqual(meta["literature_eligible_candidate_count"], 2)

        full_a, full_a_meta = select_candidate_pool(
            root,
            deterministic,
            {
                "merge_review_tiers": ["review_high_priority", "review_context_needed"],
                "literature_on_baseline_picked": True,
            },
        )
        self.assertTrue(all(candidate_literature_eligible(row) for row in full_a))
        self.assertEqual(full_a_meta["literature_candidate_scope"], "all_selected_candidates")


if __name__ == "__main__":
    unittest.main()

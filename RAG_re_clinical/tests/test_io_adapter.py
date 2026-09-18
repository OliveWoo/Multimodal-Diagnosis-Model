import copy
import json
import tempfile
import unittest
from pathlib import Path

from rag_re_clinical.io_adapter import (
    AmbiguousCandidateMatchError,
    DuplicateCandidateKeyError,
    PatientMismatchError,
    discover_jsons,
    index_merge_files,
    load_json,
    merge_candidate_context,
    patient_id,
    sha256_file,
)


def merge_fixture(patient="7"):
    exact = {
        "organism_name": "Exact species",
        "classification": "Bacterial",
        "reads": 99,
        "key_evidence": {
            "hospital_evidence_detail": {
                "culture": [
                    {
                        "specimen_type": "Blood",
                        "growth_purity": "isolated_pure",
                    }
                ]
            }
        },
        "is_likely_colonizer_or_background": False,
    }
    reviewed = {
        "organism_name": "Reviewed species",
        "classification": "Viral",
        "reads": 11,
        "key_evidence": {"support_modules": ["host"]},
    }
    return {
        "patient_id": patient,
        "report_type": "formal merge",
        "deterministic_max": {
            "rule_version": "fixture",
            "host_context": {"host_vulnerability_tier": "V2"},
            "best_available_summary": {
                "selection_mode": "formal",
                "picked_pathogens": [
                    {
                        "organism_name": "Exact species",
                        "picked_role": "Primary",
                        "caution_flags": ["keep-separated"],
                    }
                ],
            },
            "pathogen_candidates": [exact, reviewed],
            "excluded_candidates": [],
        },
        "llm_missed_candidate_review": {
            "review_high_priority": [
                {
                    "organism_name": "Reviewed species",
                    "confidence": "moderate",
                    "evidence_snapshot": {"reads": 11, "reads_tier": "R1_low"},
                }
            ],
            "review_context_needed": [
                {
                    "organism_name": "Synthetic species",
                    "confidence": "low",
                    "evidence_snapshot": {
                        "rank_priority": "2",
                        "reads": 17,
                        "reads_tier": "R2_medium",
                        "specimen_class": "S2_lower_respiratory",
                        "support_modules": ["host"],
                    },
                }
            ],
            "review_low_specificity": [],
            "review_output_policy": "fixture-policy",
        },
        "final_name_reconciliation": {"enabled": True, "merged_items": []},
    }


def frozen_fixture(patient="7"):
    return {
        "schema_version": "rag_re.output.v1",
        "patient_id": patient,
        "candidates": [
            {
                "organism_name": "Exact species",
                "canonical_organism_name": "Exact species",
                "baseline_selected": True,
            },
            {
                "organism_name": "Reviewed species",
                "canonical_organism_name": "Reviewed species",
                "baseline_selected": False,
            },
            {
                "organism_name": "Synthetic species",
                "canonical_organism_name": "Synthetic species",
                "baseline_selected": False,
            },
        ],
    }


class IoAdapterTests(unittest.TestCase):
    def test_layers_exact_and_synthetic_candidates_without_mutating_inputs(self):
        merge = merge_fixture()
        frozen = frozen_fixture()
        original_merge = copy.deepcopy(merge)
        result = merge_candidate_context(frozen, merge)

        self.assertEqual(result["patient_id"], "7")
        self.assertEqual(
            result["counts"],
            {
                "candidate_count": 3,
                "exact_species": 2,
                "synthetic_unmatched": 1,
                "ambiguous": 0,
            },
        )
        exact, reviewed, synthetic = result["candidates"]
        self.assertEqual(exact["match_scope"], "exact_species")
        self.assertEqual(exact["deterministic_candidate"]["reads"], 99)
        self.assertEqual(exact["picked_metadata"]["picked_role"], "Primary")
        self.assertEqual(
            exact["deterministic_candidate"]["key_evidence"]
            ["hospital_evidence_detail"]["culture"][0]["growth_purity"],
            "isolated_pure",
        )
        self.assertEqual(reviewed["review_metadata"]["tier"], "review_high_priority")
        self.assertEqual(reviewed["effective_candidate"]["reads"], 11)
        self.assertEqual(synthetic["match_scope"], "synthetic_unmatched")
        self.assertIsNone(synthetic["deterministic_candidate"])
        self.assertEqual(
            synthetic["synthetic_fallback"]["candidate_fields"]["reads"], 17
        )
        self.assertEqual(
            synthetic["synthetic_fallback"]["evidence_snapshot"]["support_modules"],
            ["host"],
        )
        self.assertEqual(merge, original_merge)

    def test_matching_is_alnum_exact_but_not_genus_or_fuzzy(self):
        frozen = frozen_fixture()
        frozen["candidates"] = [
            {"organism_name": "Exact-species"},
            {"organism_name": "Reviewed"},
        ]
        result = merge_candidate_context(frozen, merge_fixture())
        self.assertEqual(result["candidates"][0]["match_scope"], "exact_species")
        self.assertEqual(result["candidates"][1]["match_scope"], "synthetic_unmatched")
        self.assertIsNone(result["candidates"][1]["review_metadata"])

    def test_duplicate_and_ambiguous_keys_are_rejected(self):
        duplicate = merge_fixture()
        duplicate["deterministic_max"]["pathogen_candidates"].append(
            {"organism_name": "Exact-species"}
        )
        with self.assertRaises(DuplicateCandidateKeyError):
            merge_candidate_context(frozen_fixture(), duplicate)

        ambiguous = merge_fixture()
        ambiguous["llm_missed_candidate_review"]["review_context_needed"].append(
            {"organism_name": "Reviewed species", "evidence_snapshot": {}}
        )
        with self.assertRaises(AmbiguousCandidateMatchError) as caught:
            merge_candidate_context(frozen_fixture(), ambiguous)
        self.assertEqual(caught.exception.context["match_scope"], "ambiguous")

    def test_patient_mismatch_and_duplicate_merge_files_are_rejected(self):
        with self.assertRaises(PatientMismatchError):
            merge_candidate_context(frozen_fixture("1"), merge_fixture("2"))

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for name in ("one.json", "two.json"):
                (root / name).write_text(json.dumps(merge_fixture()), encoding="utf-8")
            with self.assertRaises(PatientMismatchError):
                index_merge_files(root)

    def test_file_helpers_support_bom_and_deterministic_discovery(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            second = root / "b.json"
            first = root / "a.json"
            first.write_text("\ufeff" + json.dumps(merge_fixture("1")), encoding="utf-8")
            second.write_text(json.dumps(merge_fixture("2")), encoding="utf-8")
            paths = discover_jsons(root)
            self.assertEqual([path.name for path in paths], ["a.json", "b.json"])
            self.assertEqual(load_json(first)["patient_id"], "1")
            self.assertEqual(patient_id(load_json(first), first), "1")
            self.assertEqual(len(sha256_file(first)), 64)
            self.assertEqual(set(index_merge_files(root)), {"1", "2"})


if __name__ == "__main__":
    unittest.main()

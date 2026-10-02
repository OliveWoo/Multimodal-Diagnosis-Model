import csv
import json
import tempfile
import unittest
from pathlib import Path

from tools.build_independent_multi_assay_screening_shadow import build_shadow, route_candidate


POLICY = json.loads(Path("rules/independent_multi_assay_screening_v1.json").read_text(encoding="utf-8"))


def observation(*, seq_id="S1", reads=10, selected=True, molecule="DNA", rank=1):
    return {
        "seq_id": seq_id,
        "reads": reads,
        "selected": selected,
        "detected": True,
        "test_eligible_for_union": True,
        "availability_status": "evaluable",
        "nucleic_type": molecule,
        "rank_by_reads_in_retained_universe": rank,
        "test_qc_status": "partial_evaluable",
    }


def candidate(*observations):
    return {
        "organism_name": "Example organism",
        "organism_key": "exampleorganism",
        "category": "bacterial",
        "observations": list(observations),
    }


def taxonomy(family="typical_respiratory_pathogen", mapping="exact_species", biological_class="bacterium"):
    return {
        "primary_rule_family": family,
        "mapping_status": mapping,
        "biological_class": biological_class,
    }


class TestIndependentRouting(unittest.TestCase):
    def test_selected_positive_priority_enters_scorer(self):
        result = route_candidate(candidate(observation()), taxonomy(), POLICY)
        self.assertEqual(result["disposition"], "scorer_entry")

    def test_filtered_only_and_zero_reads_do_not_enter(self):
        result = route_candidate(
            candidate(observation(reads=0), observation(seq_id="S2", selected=False)), taxonomy(), POLICY
        )
        self.assertEqual(result["disposition"], "qc_only")
        self.assertEqual(result["selected_positive_test_count"], 0)

    def test_unknown_selected_positive_requires_taxonomy_review_even_if_low_rank(self):
        result = route_candidate(
            candidate(observation(rank=17)), taxonomy("unmapped_or_uncertain", "unmapped"), POLICY
        )
        self.assertEqual(result["disposition"], "clinical_review")
        self.assertIn("TAXONOMY_REVIEW_REQUIRED", result["rule_ids"])

    def test_rna_only_nonviral_is_review_not_scorer(self):
        result = route_candidate(
            candidate(observation(molecule="RNA")), taxonomy(), POLICY
        )
        self.assertEqual(result["disposition"], "clinical_review")

    def test_cross_test_is_technical_repeat_not_independent_confirmation(self):
        result = route_candidate(
            candidate(observation(rank=8), observation(seq_id="S2", molecule="RNA", rank=9)),
            taxonomy(), POLICY,
        )
        self.assertEqual(result["disposition"], "scorer_entry")
        self.assertTrue(result["technical_repeat"])
        self.assertTrue(result["cross_molecule_selected"])
        self.assertIn("SAME_CASE", "|".join(result["rule_ids"]))

    def test_context_family_does_not_enter_scorer(self):
        result = route_candidate(candidate(observation()), taxonomy("candida_or_yeast"), POLICY)
        self.assertEqual(result["disposition"], "clinical_review")

    def test_cases_stay_separate_and_baseline_only_audit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory = root / "inventory"
            patients = inventory / "patients"
            patients.mkdir(parents=True)
            selected = candidate(observation())
            selected["organism_name"] = "Streptococcus pneumoniae"
            selected["organism_key"] = "streptococcuspneumoniae"
            weak = candidate(observation(rank=10))
            weak["organism_name"] = selected["organism_name"]
            weak["organism_key"] = selected["organism_key"]
            document = {
                "patient_id": "1",
                "cases": [
                    {
                        "case_review_id": case_id,
                        "tests": [{"seq_id": "S1", "qc": {"status": "partial_evaluable"}}],
                        "candidate_union": {"selected_union": [item], "detected_but_filtered_union": []},
                    }
                    for case_id, item in (("A", selected), ("B", weak))
                ],
            }
            (patients / "NGS_patient_1_multi_assay_mngs_evidence.json").write_text(
                json.dumps(document), encoding="utf-8"
            )
            (inventory / "selected_dna_comparison.csv").write_text(
                "patient_id,organism_key,comparison_status\n"
                "1,streptococcuspneumoniae,retained_from_selected_dna\n", encoding="utf-8"
            )
            output = root / "shadow"
            summary = build_shadow(inventory, output, Path("rules/independent_multi_assay_screening_v1.json"))
            self.assertEqual(summary["case_organism_count"], 2)
            with (output / "all_decisions.csv").open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual({row["case_review_id"] for row in rows}, {"A", "B"})
            self.assertEqual({row["disposition"] for row in rows}, {"scorer_entry", "clinical_review"})
            self.assertTrue(all(row["old_selected_dna_status_audit_only"] == "retained_from_selected_dna" for row in rows))


if __name__ == "__main__":
    unittest.main()

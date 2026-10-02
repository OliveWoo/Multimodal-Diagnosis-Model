import json
import unittest
from pathlib import Path

from tools.build_test_aware_mngs_screening_shadow import route_candidate


POLICY = json.loads(Path("rules/multi_assay_mngs_screening_policy.json").read_text(encoding="utf-8"))


def candidate(**overrides):
    value = {
        "union_status": "selected_any_test",
        "category": "bacterial",
        "detected_test_count": 1,
        "dna_detected_test_count": 1,
        "rna_detected_test_count": 0,
        "cross_molecule_support": False,
        "best_analytical_rank_priority": 5,
        "max_reads_single_test": 10,
    }
    value.update(overrides)
    return value


def taxonomy(family, biological_class="bacterial", needs_review=False):
    return {
        "primary_rule_family": family,
        "biological_class": biological_class,
        "needs_literature_review": needs_review,
    }


class TestRouting(unittest.TestCase):
    def test_baseline_is_always_retained(self):
        action, rules, _ = route_candidate(
            candidate(union_status="filtered_only"),
            taxonomy("environmental_low_specificity"),
            "retained_from_selected_dna",
            POLICY,
        )
        self.assertEqual(action, "retain_baseline")
        self.assertIn("BASELINE_001", rules)

    def test_cross_molecule_priority_family_is_forwarded(self):
        action, _, _ = route_candidate(
            candidate(
                detected_test_count=2,
                dna_detected_test_count=1,
                rna_detected_test_count=1,
                cross_molecule_support=True,
            ),
            taxonomy("typical_respiratory_pathogen"),
            "new_from_multi_assay_union",
            POLICY,
        )
        self.assertEqual(action, "forward_priority_review")

    def test_cross_molecule_does_not_double_count_for_mold_priority(self):
        action, _, _ = route_candidate(
            candidate(
                category="fungal",
                detected_test_count=2,
                dna_detected_test_count=1,
                rna_detected_test_count=1,
                cross_molecule_support=True,
                best_analytical_rank_priority=5,
            ),
            taxonomy("mold_or_opportunistic_fungus", biological_class="fungus"),
            "new_from_multi_assay_union",
            POLICY,
        )
        self.assertEqual(action, "forward_context_review")

    def test_rna_only_nonviral_is_context_not_priority(self):
        action, _, _ = route_candidate(
            candidate(
                detected_test_count=1,
                dna_detected_test_count=0,
                rna_detected_test_count=1,
                best_analytical_rank_priority=1,
            ),
            taxonomy("hospital_or_nonfermenter_gnb"),
            "new_from_multi_assay_union",
            POLICY,
        )
        self.assertEqual(action, "forward_context_review")

    def test_strong_context_family_is_only_context(self):
        action, _, _ = route_candidate(
            candidate(detected_test_count=2, best_analytical_rank_priority=1),
            taxonomy("candida_or_yeast", biological_class="fungal"),
            "new_from_multi_assay_union",
            POLICY,
        )
        self.assertEqual(action, "forward_context_review")

    def test_weak_unmapped_is_held(self):
        action, _, _ = route_candidate(
            candidate(max_reads_single_test=1),
            taxonomy("unmapped_or_uncertain", needs_review=True),
            "new_from_multi_assay_union",
            POLICY,
        )
        self.assertEqual(action, "hold_not_forwarded")

    def test_filtered_only_is_never_forwarded(self):
        action, _, _ = route_candidate(
            candidate(union_status="filtered_only", best_analytical_rank_priority=1),
            taxonomy("high_consequence_opportunistic"),
            "new_from_multi_assay_union",
            POLICY,
        )
        self.assertEqual(action, "hold_not_forwarded")


if __name__ == "__main__":
    unittest.main()

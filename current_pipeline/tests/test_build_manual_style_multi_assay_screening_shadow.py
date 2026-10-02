import json
import unittest
from pathlib import Path

from tests.test_build_independent_multi_assay_screening_shadow import candidate, observation, taxonomy
from tools.build_manual_style_multi_assay_screening_shadow import route_candidate


POLICY = json.loads(Path("rules/manual_style_multi_assay_entry_v2.json").read_text(encoding="utf-8"))


class TestManualStyleEntry(unittest.TestCase):
    def test_positive_context_enters_low_tier_even_without_top_rank(self):
        result = route_candidate(
            candidate(observation(rank=3)), taxonomy("candida_or_yeast"), POLICY
        )
        self.assertEqual((result["disposition"], result["screening_tier"]), ("scorer_entry", "low_colonizer"))

    def test_repeated_background_signal_enters_low_tier(self):
        result = route_candidate(
            candidate(observation(rank=4), observation(seq_id="S2", rank=5)),
            taxonomy("gi_urinary_or_nonpulmonary_prone"), POLICY,
        )
        self.assertEqual((result["disposition"], result["screening_tier"]), ("scorer_entry", "low_colonizer"))

    def test_rna_only_priority_is_not_removed(self):
        result = route_candidate(
            candidate(observation(molecule="RNA", rank=3)), taxonomy(), POLICY
        )
        self.assertEqual((result["disposition"], result["screening_tier"]), ("scorer_entry", "low_colonizer"))

    def test_unmapped_still_requires_taxonomy_review(self):
        result = route_candidate(
            candidate(observation()), taxonomy("unmapped_or_uncertain", "unmapped"), POLICY
        )
        self.assertEqual(result["disposition"], "clinical_review")

    def test_filtered_only_and_zero_reads_remain_qc(self):
        result = route_candidate(
            candidate(observation(reads=0), observation(seq_id="S2", selected=False)), taxonomy(), POLICY
        )
        self.assertEqual(result["disposition"], "qc_only")

    def test_commensal_virome_remains_qc_even_when_repeated(self):
        result = route_candidate(
            candidate(observation(), observation(seq_id="S2")),
            taxonomy("commensal_virome_or_endogenous_element"), POLICY,
        )
        self.assertEqual(result["disposition"], "qc_only")


if __name__ == "__main__":
    unittest.main()

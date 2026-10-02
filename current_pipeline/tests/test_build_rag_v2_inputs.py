from __future__ import annotations

import unittest

from tools import build_rag_v2_inputs as rag_v2


class RagV2InputTests(unittest.TestCase):
    def test_same_genus_is_not_same_organism(self) -> None:
        summary = {
            "hospital_organism_evidence": [
                {"organism_name": "Candida albicans", "evidence_modules": ["culture"]}
            ]
        }
        result = rag_v2.hospital_evidence_for_name(summary, "Candida tropicalis")
        self.assertEqual(result["same_organism_or_alias"], [])
        self.assertEqual(
            [item["organism_name"] for item in result["same_genus_or_related"]],
            ["Candida albicans"],
        )

    def test_forbidden_benchmark_key_is_detected(self) -> None:
        payload = {"case_id": "P14_Candida_tropicalis", "nested": {"benchmark_answer": "x"}}
        self.assertEqual(rag_v2.find_forbidden_keys(payload), ["$.nested.benchmark_answer"])

    def test_allowed_organism_name_is_not_value_scanned(self) -> None:
        payload = {
            "case_id": "P14_Candida_tropicalis",
            "organism": {"display_name": "Candida tropicalis"},
        }
        self.assertEqual(rag_v2.find_forbidden_keys(payload), [])

    def test_unterminated_or_relaxed_genus_is_not_needed_for_exact_match(self) -> None:
        summary = {
            "hospital_organism_evidence": [
                {"organism_name": "SARS-CoV-2", "evidence_modules": ["molecular_microbiology"]}
            ]
        }
        result = rag_v2.hospital_evidence_for_name(summary, "COVID-19")
        self.assertEqual(len(result["same_organism_or_alias"]), 1)

    def test_negative_evidence_delegates_to_central_profile(self) -> None:
        result = rag_v2.negative_evidence(
            {
                "dominance_tier": "D0_not_top",
                "reads_tier": "R1_low",
                "support_modules": ["host"],
                "non_host_support_modules": [],
                "specimen_class": "S2_lower_respiratory",
            },
            {"same_organism_or_alias": [], "same_genus_or_related": []},
            [],
            organism_name="Candida albicans",
        )
        codes = [item["code"] for item in result]
        self.assertIn("CANDIDA_NO_STRONG_INVASIVE_EVIDENCE", codes)


if __name__ == "__main__":
    unittest.main()

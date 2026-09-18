from __future__ import annotations

import unittest

from tools import build_deterministic_summary as summary


def label(name: str, level: str = "Level 2") -> dict[str, object]:
    return {
        "organism_name": name,
        "classification": "Bacterial",
        "causative_level": level,
        "key_evidence": {"specimen": "BALF"},
    }


class ApprovedGroupMemberSummaryTests(unittest.TestCase):
    def test_group_and_species_are_linked_but_not_collapsed(self) -> None:
        evidence = summary.collect_hospital_organism_evidence(
            culture={
                "organism_labels": [
                    label("Acinetobacter baumannii"),
                    label("acineto calc baumannii complex", "Level 1"),
                ]
            },
            filmarray=None,
            image=None,
            molecular=None,
        )

        self.assertEqual(2, len(evidence))
        by_name = {item["organism_name"]: item for item in evidence}
        member = by_name["Acinetobacter baumannii"]
        group = by_name["acineto calc baumannii complex"]
        self.assertTrue(group["group_level_support_only"])
        self.assertEqual(
            "approved_group_member_match",
            group["approved_detected_members"][0]["relationship"],
        )
        self.assertTrue(member["approved_group_support"][0]["not_exact_species_support"])
        self.assertEqual("Level 2", member["best_hospital_level"])
        self.assertEqual("Level 1", group["best_hospital_level"])

    def test_same_genus_species_is_not_group_member_support(self) -> None:
        evidence = summary.collect_hospital_organism_evidence(
            culture={
                "organism_labels": [
                    label("Acinetobacter pittii"),
                    label("Acinetobacter baumannii complex"),
                ]
            },
            filmarray=None,
            image=None,
            molecular=None,
        )

        by_name = {item["organism_name"]: item for item in evidence}
        self.assertNotIn("approved_group_support", by_name["Acinetobacter pittii"])
        self.assertNotIn("approved_detected_members", by_name["Acinetobacter baumannii complex"])


if __name__ == "__main__":
    unittest.main()

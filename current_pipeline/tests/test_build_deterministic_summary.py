from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools import build_deterministic_summary as summary


def label(name: str, level: str = "Level 2") -> dict[str, object]:
    return {
        "organism_name": name,
        "classification": "Bacterial",
        "causative_level": level,
        "key_evidence": {"specimen": "BALF"},
    }


class ApprovedGroupMemberSummaryTests(unittest.TestCase):
    def test_raw_molecular_source_uses_deterministic_fallback_when_agent_missing(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            patient_dir = Path(temp_dir) / "NGS_patient_999_json"
            patient_dir.mkdir()
            source_path = patient_dir / "NGS_patient_999_molecular_microbiology.json"
            source_path.write_text(
                json.dumps(
                    [
                        {
                            "test": "Non-panel Molecular Test",
                            "sample": "BALF",
                            "target": "Pneumocystis jirovecii",
                            "result": "Detected",
                            "reported_time": "2026-09-01",
                        }
                    ]
                ),
                encoding="utf-8",
            )

            payload = summary.build_summary(patient_dir)

            self.assertEqual(
                "normalized_molecular_microbiology_source",
                payload["evidence_preservation"]["deterministic_source_fallbacks"][
                    "molecular_microbiology"
                ],
            )
            molecular = payload["module_summaries"]["molecular_microbiology"]
            self.assertEqual("available", payload["module_availability"]["molecular_microbiology"])
            self.assertEqual("Likely", molecular["infection_likelihood"])
            pjp = next(
                item
                for item in payload["hospital_organism_evidence"]
                if item["organism_name"] == "Pneumocystis jirovecii"
            )
            self.assertEqual("Level 2", pjp["best_hospital_level"])
            self.assertEqual(
                str(source_path), payload["source_files"]["raw_molecular_microbiology"]
            )

    def test_raw_molecular_source_replaces_legacy_agent_without_observation_chain(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            patient_dir = Path(temp_dir) / "NGS_patient_998_json"
            agent_dir = patient_dir / "agent_outputs"
            agent_dir.mkdir(parents=True)
            (patient_dir / "NGS_patient_998_molecular_microbiology.json").write_text(
                json.dumps(
                    [
                        {
                            "sample": "Blood",
                            "target": "EBV",
                            "result": "Detected",
                        },
                        {
                            "sample": "Blood",
                            "target": "CMV",
                            "result": "Not Detected",
                        },
                    ]
                ),
                encoding="utf-8",
            )
            (agent_dir / "NGS_patient_998_molecular_microbiology_agent.json").write_text(
                json.dumps(
                    {
                        "rule_version": "legacy_v1",
                        "organism_labels": [
                            {
                                "organism_name": "CMV",
                                "classification": "Viral",
                                "causative_level": "Level 2",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            payload = summary.build_summary(patient_dir)

            self.assertEqual(
                "normalized_molecular_microbiology_source_complete_overlay",
                payload["evidence_preservation"]["deterministic_source_fallbacks"][
                    "molecular_microbiology"
                ],
            )
            molecular_names = {
                item["organism_name"]
                for item in payload["hospital_organism_evidence"]
                if "molecular_microbiology" in item.get("evidence_modules", [])
            }
            self.assertEqual({"EBV"}, molecular_names)

    def test_acute_severity_alone_is_not_opportunistic_host_support(self) -> None:
        context = summary.host_context(
            {
                "host_state": {
                    "immunocompromise_tier": "H1",
                    "acute_instability_tier": "A3",
                    "host_vulnerability_tier": "V3",
                    "expanded_candidate_policy": True,
                    "opportunistic_coverage_level": "O2",
                },
                "structured_host_evidence": {
                    "medication_profile": {
                        "host_risk_evidence_strength": "intermediate"
                    }
                },
            }
        )
        self.assertFalse(context["opportunistic_host_support"])
        self.assertTrue(context["opportunistic_host_review_context"])
        self.assertIn("severity_not_opportunistic_support", context["opportunistic_host_support_basis"][0])

    def test_provisional_h2_explicit_false_is_not_reenabled_by_tier(self) -> None:
        context = summary.host_context(
            {
                "host_state": {
                    "immunocompromise_tier": "H2",
                    "acute_instability_tier": "A2",
                    "host_vulnerability_tier": "V3",
                    "expanded_candidate_policy": True,
                    "opportunistic_coverage_level": "O2",
                    "opportunistic_host_support": False,
                    "opportunistic_host_support_status": "insufficient",
                },
                "structured_host_evidence": {
                    "medication_profile": {
                        "host_risk_evidence_strength": "intermediate"
                    }
                },
            }
        )
        self.assertFalse(context["opportunistic_host_support"])
        self.assertTrue(context["opportunistic_host_review_context"])
        self.assertEqual(
            "insufficient", context["opportunistic_host_support_status"]
        )

    def test_shadow_agent_suffix_tag_is_explicit_and_does_not_change_baseline(self) -> None:
        patient_dir = Path("D:/example/NGS_patient_14_json")
        self.assertEqual(
            "cbc_other_lab_agent",
            summary.tagged_agent_suffix("cbc_other_lab_agent", None),
        )
        self.assertEqual(
            "cbc_other_lab_agent_evidence_v2_shadow_20260917",
            summary.tagged_agent_suffix(
                "cbc_other_lab_agent", "evidence_v2_shadow_20260917"
            ),
        )
        self.assertEqual(
            "NGS_patient_14_cbc_other_lab_agent_evidence_v2_shadow_20260917.json",
            summary.agent_path(
                patient_dir,
                summary.tagged_agent_suffix(
                    "cbc_other_lab_agent", "evidence_v2_shadow_20260917"
                ),
            ).name,
        )

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

    def test_parainfluenza_panel_group_supports_but_does_not_replace_exact_hpiv3(self) -> None:
        evidence = summary.collect_hospital_organism_evidence(
            culture=None,
            filmarray={
                "organism_labels": [
                    {
                        **label("Parainfluenza Virus"),
                        "classification": "Viral",
                    }
                ]
            },
            image=None,
            molecular=None,
        )
        # The panel result remains broad.  The approved relationship is used
        # only when an exact mNGS member is evaluated downstream.
        self.assertEqual(1, len(evidence))
        self.assertEqual("Parainfluenza Virus", evidence[0]["organism_name"])

    def test_image_respiratory_virus_mention_is_context_only_evidence(self) -> None:
        evidence = summary.collect_hospital_organism_evidence(
            culture=None,
            filmarray=None,
            molecular=None,
            image={
                "diagnostic_mentions": [
                    {
                        "concept_name": "SARS-CoV-2",
                        "concept_type": "organism_or_disease",
                        "biological_class": "virus",
                        "assertion": "mentioned",
                        "evidence_role": "imaging_context_only",
                        "microbiologic_confirmation": False,
                        "source_sentence": "Covid-19 pneumonia.",
                    }
                ]
            },
        )

        self.assertEqual(1, len(evidence))
        item = evidence[0]
        self.assertEqual("SARS-CoV-2", item["organism_name"])
        self.assertEqual("Level 4", item["best_hospital_level"])
        row = item["module_evidence"]["image"][0]
        self.assertTrue(row["context_candidate_eligible"])
        self.assertFalse(row["pulmonary_causality_support"])

    def test_image_bacterial_mention_is_not_promoted_to_hospital_evidence(self) -> None:
        evidence = summary.collect_hospital_organism_evidence(
            culture=None,
            filmarray=None,
            molecular=None,
            image={
                "diagnostic_mentions": [
                    {
                        "concept_name": "Klebsiella pneumoniae",
                        "concept_type": "organism_or_disease",
                        "biological_class": "bacteria",
                        "assertion": "mentioned",
                        "microbiologic_confirmation": False,
                    }
                ]
            },
        )

        self.assertEqual([], evidence)

    def test_culture_observation_count_is_preserved_for_downstream_rules(self) -> None:
        culture_label = label("Trichosporon asahii")
        culture_label["classification"] = "Fungal"
        culture_label["evidence_observation_ids"] = ["CUL-001", "CUL-002", "CUL-003"]
        culture_label["key_evidence"] = {
            "specimen_type": "Lower BAL",
            "specimen_category": "Lower_Respiratory",
            "quantity_tier": "Q3",
            "quantitation_status": "reported",
        }

        evidence = summary.collect_hospital_organism_evidence(
            culture={"organism_labels": [culture_label]},
            filmarray=None,
            image=None,
            molecular=None,
        )

        culture_row = evidence[0]["module_evidence"]["culture"][0]
        self.assertEqual(3, culture_row["observation_count"])
        self.assertEqual(
            ["CUL-001", "CUL-002", "CUL-003"],
            culture_row["evidence_observation_ids"],
        )

    def test_culture_observation_count_is_preserved_for_downstream_rules(self) -> None:
        culture_label = label("Trichosporon asahii")
        culture_label["classification"] = "Fungal"
        culture_label["evidence_observation_ids"] = ["CUL-001", "CUL-002", "CUL-003"]
        culture_label["key_evidence"] = {
            "specimen_type": "Lower BAL",
            "specimen_category": "Lower_Respiratory",
            "quantity_tier": "Q3",
            "quantitation_status": "reported",
        }

        evidence = summary.collect_hospital_organism_evidence(
            culture={"organism_labels": [culture_label]},
            filmarray=None,
            image=None,
            molecular=None,
        )

        culture_row = evidence[0]["module_evidence"]["culture"][0]
        self.assertEqual(3, culture_row["observation_count"])
        self.assertEqual(
            ["CUL-001", "CUL-002", "CUL-003"],
            culture_row["evidence_observation_ids"],
        )


if __name__ == "__main__":
    unittest.main()

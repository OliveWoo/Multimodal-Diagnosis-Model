import csv
import json
import tempfile
import unittest
from pathlib import Path

from tools.build_mandatory_primary_reporting_shadow import run


TEST_TMP = Path(__file__).parents[1] / ".test_tmp"
TEST_TMP.mkdir(exist_ok=True)


REPORT_FIELDS = [
    "patient_id",
    "organism_name",
    "taxonomy_family",
    "clinical_decision",
    "original_selection_role",
    "possible_reporting_role",
    "selected_for_complete_report",
    "possible_outcome",
    "possible_route",
    "possible_blockers",
    "possible_cautions",
    "specimen_context",
    "precise_identity",
    "best_rank",
    "selected_positive_test_count",
    "reproducibility_axis",
    "cross_molecule_selected",
    "direct_hospital_level",
    "event_aligned_direct_positive_count",
    "repeated_sterile_positive_count",
    "identity_resolution",
    "taxonomy_rule_source",
    "resistance_target",
    "resistance_status",
    "final_reporting_tier",
    "fallback_rank",
    "fallback_primary",
    "final_reporting_confidence",
    "patient_decision_evidence_mode",
    "candidate_decision_evidence_mode",
    "positive_or_contextual_sources",
]


def analytical_candidate(name: str, family: str, rank: int, reads: int) -> dict:
    return {
        "patient_id": "2",
        "organism_name": name,
        "taxonomy_family": family,
        "specimen_context": "blood_or_systemic",
        "decision": "review_low_specificity",
        "analytical_profile": {
            "reproducibility_axis": False,
            "cross_molecule_selected": False,
            "per_test_signals": [
                {
                    "reads": reads,
                    "rank_in_retained_universe": rank,
                    "qc_status": "evaluable",
                }
            ],
        },
        "hospital_profile": {"best_direct_level": None},
        "source_entry": {
            "taxonomy_family": family,
            "specimen_context": "blood_or_systemic",
            "taxonomy_mapping_status": "genus_inherited",
            "taxonomy_identity_resolution": "species_name_with_genus_inherited_routing",
            "taxonomy_rule_source": "curated_genus_profile",
        },
    }


class MandatoryPrimaryTests(unittest.TestCase):
    def test_generic_family_ranking_adds_fallback_without_promoting_certainty(self) -> None:
        with tempfile.TemporaryDirectory(dir=TEST_TMP) as raw:
            root = Path(raw)
            release = root / "release"
            report_dir = release / "10_reporting_v3_route_a"
            analytical_dir = release / "02_analytical_scorer" / "patient_outputs"
            report_dir.mkdir(parents=True)
            analytical_dir.mkdir(parents=True)
            with (report_dir / "possible_pathogen_decisions.csv").open(
                "w", encoding="utf-8-sig", newline=""
            ) as handle:
                writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS)
                writer.writeheader()
                writer.writerow({
                    "patient_id": "1",
                    "organism_name": "Escherichia coli",
                    "taxonomy_family": "typical_respiratory_pathogen",
                    "clinical_decision": "picked_shadow",
                    "possible_reporting_role": "strict_picked_report",
                    "selected_for_complete_report": "True",
                    "specimen_context": "blood_or_systemic",
                    "best_rank": "1",
                    "final_reporting_tier": "Picked",
                })
            (analytical_dir / "NGS_patient_1_test_aware_deterministic_shadow.json").write_text(
                json.dumps({"patient_organism_decisions": []}), encoding="utf-8"
            )
            (analytical_dir / "NGS_patient_2_test_aware_deterministic_shadow.json").write_text(
                json.dumps({
                    "patient_organism_decisions": [
                        analytical_candidate(
                            "Staphylococcus epidermidis",
                            "skin_airway_colonizer_prone",
                            1,
                            100,
                        ),
                        analytical_candidate(
                            "Klebsiella pneumoniae",
                            "typical_respiratory_pathogen",
                            6,
                            18,
                        ),
                    ]
                }),
                encoding="utf-8",
            )
            policy = Path(__file__).parents[1] / "rules" / "mandatory_primary_reporting_shadow_v1.json"
            output = root / "output"

            summary = run(release, output, policy)

            self.assertEqual(summary["primary_candidate_count"], 2)
            self.assertEqual(summary["mandatory_fallback_primary_count"], 1)
            with (output / "primary_candidates.csv").open(
                encoding="utf-8-sig"
            ) as handle:
                rows = list(csv.DictReader(handle))
            patient_two = next(row for row in rows if row["patient_id"] == "2")
            self.assertEqual(patient_two["organism_name"], "Klebsiella pneumoniae")
            self.assertEqual(patient_two["final_reporting_tier"], "Fallback-BestAvailable")
            self.assertEqual(patient_two["clinical_decision"], "review_low_specificity")
            self.assertEqual(patient_two["clinical_certainty_changed"], "False")


if __name__ == "__main__":
    unittest.main()

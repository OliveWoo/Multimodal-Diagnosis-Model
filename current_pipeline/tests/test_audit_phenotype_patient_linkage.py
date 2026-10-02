from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
import json

from tools import audit_phenotype_patient_linkage as linkage


class PhenotypePatientLinkageAuditTests(unittest.TestCase):
    def test_lab_normalization_matches_equivalent_values(self) -> None:
        self.assertEqual(
            linkage.lab_feature("2024/02/22 10:55", "WBC (x1000/uL)", "7.60"),
            linkage.lab_feature("2024-02-22", "WBC", 7.6),
        )

    def test_patient_number_does_not_affect_pair_score(self) -> None:
        shared = frozenset({"2024-01-01|wbc|8"})
        empty = frozenset()
        kh = linkage.PatientFeatures(9, "kh", shared, empty, empty, empty)
        phenotype = linkage.PatientFeatures(42, "phenotype", shared, empty, empty, empty)
        result = linkage.compare_pair(
            kh,
            phenotype,
            {
                "labs": {"2024-01-01|wbc|8": 1},
                "microbiology": {},
                "diagnoses": {},
                "events": {},
            },
            2,
        )
        self.assertGreater(result["score"], 0)
        self.assertEqual(42, result["phenotype_patient_number"])

    def test_strong_status_requires_same_number_to_be_best(self) -> None:
        same = {
            "phenotype_patient_number": 9,
            "score": 40,
            "evidence": {
                "exact_dated_lab_match_count": 10,
                "phenotype_unique_lab_match_count": 8,
            },
        }
        other = {
            "phenotype_patient_number": 10,
            "score": 20,
            "evidence": {"exact_dated_lab_match_count": 2},
        }
        self.assertEqual(
            "local_strong_same_patient_support",
            linkage.classify_link(same=same, best=same, second=other),
        )
        self.assertEqual(
            "local_ambiguous_best_match_is_different_patient",
            linkage.classify_link(same=same, best=other, second=other),
        )

    def test_sync_report_with_numeric_patient_ids_is_loaded(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "report.json"
            path.write_text(
                json.dumps({"patients": [{"patient_id": 9, "specimens": [
                    ["sample-9", "2024-02-22 10:55:00", "Blood"]
                ]}]}),
                encoding="utf-8",
            )
            self.assertEqual(
                "sample-9", linkage.load_mngs_metadata(path)[9]["specimens"][0]["specimen_id"]
            )

    def test_identical_labs_form_duplicate_export_cluster_not_identity_conflict(self) -> None:
        labs = frozenset(f"2024-01-01|wbc|{number}" for number in range(20))
        empty = frozenset()
        phenotype = [
            linkage.PatientFeatures(number, str(number), labs, empty, empty, empty)
            for number in (9, 11)
        ]
        clusters, lookup = linkage.duplicate_export_clusters(phenotype)
        self.assertEqual([[9, 11]], clusters)
        same = {"phenotype_patient_number": 9, "score": 20, "evidence": {}}
        other = {"phenotype_patient_number": 11, "score": 21, "evidence": {}}
        self.assertEqual(
            "local_same_patient_duplicate_export_cluster",
            linkage.classify_link(
                same=same,
                best=other,
                second=other,
                same_identity_cluster=lookup[9],
            ),
        )

    def test_blood_index_does_not_use_lower_respiratory_date_as_episode_match(self) -> None:
        match = {"evidence": {
            "explicit_index_mngs_event_match_count": 0,
            "same_date_lower_respiratory_event_match_count": 1,
        }}
        self.assertEqual(
            (0, 0),
            linkage.episode_anchor_strength(match, index_is_lower_respiratory=False),
        )
        self.assertEqual(
            "index_infection_episode_unresolved",
            linkage.classify_episode(
                match,
                phenotype_available=True,
                index_is_lower_respiratory=False,
                ambiguous=False,
            ),
        )
        self.assertEqual(
            "multiple_exports_share_index_date_anchor",
            linkage.classify_episode(
                None,
                phenotype_available=True,
                index_is_lower_respiratory=True,
                ambiguous=True,
            ),
        )


if __name__ == "__main__":
    unittest.main()

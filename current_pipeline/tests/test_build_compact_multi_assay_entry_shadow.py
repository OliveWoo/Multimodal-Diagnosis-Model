import csv
import json
import tempfile
import unittest
from pathlib import Path

from tests.test_build_independent_multi_assay_screening_shadow import candidate, observation
from tools.build_compact_multi_assay_entry_shadow import build_shadow, compact_candidate


POLICY = json.loads(Path("rules/manual_style_multi_assay_entry_v2.json").read_text(encoding="utf-8"))


def make_test(seq_id, *, issues=None):
    return {
        "seq_id": seq_id,
        "qc": {"status": "partial_evaluable" if issues else "evaluable", "issues": issues or [], "warnings": []},
    }


class TestCompactEntryShadow(unittest.TestCase):
    def test_missing_normalization_denominator_keeps_positive_organism_reads(self):
        item = candidate(observation(seq_id="DNA1", reads=35, rank=4))
        item["organism_name"] = "Streptococcus pneumoniae"
        item["organism_key"] = "streptococcuspneumoniae"
        case = {"case_review_id": "A", "specimen_site": "BALF"}
        profile = compact_candidate(item, case, {"DNA1": make_test("DNA1", issues=["invalid_or_missing_input_reads"])}, POLICY)
        self.assertEqual(profile["disposition"], "scorer_entry")
        self.assertEqual(profile["selected_test_signals"][0]["reads"], 35)
        self.assertFalse(profile["selected_test_signals"][0]["normalization_available"])
        self.assertIn("RPM_UNAVAILABLE_USE_WITHIN_TEST_SIGNAL_ONLY", profile["rule_ids"])

    def test_two_tests_remain_separate_and_filtered_is_traceable(self):
        item = candidate(
            observation(seq_id="DNA1", reads=40, molecule="DNA", rank=2),
            observation(seq_id="RNA1", reads=7, molecule="RNA", rank=3),
            observation(seq_id="DNA2", reads=5, selected=False, molecule="DNA", rank=4),
            observation(seq_id="RNA2", reads=0, selected=False, molecule="RNA", rank=5),
        )
        item["organism_name"] = "Streptococcus pneumoniae"
        item["organism_key"] = "streptococcuspneumoniae"
        tests = {seq_id: make_test(seq_id) for seq_id in ("DNA1", "RNA1", "DNA2", "RNA2")}
        profile = compact_candidate(item, {"case_review_id": "A", "specimen_site": "Blood"}, tests, POLICY)
        self.assertEqual([signal["reads"] for signal in profile["selected_test_signals"]], [40, 7])
        self.assertEqual(profile["filtered_positive_test_refs"], ["DNA2"])
        self.assertEqual(profile["specimen_context"], "blood_or_systemic")
        self.assertIn("BLOOD_SOURCE_REQUIRES_SITE_INTERPRETATION", profile["rule_ids"])
        self.assertNotIn("reads_sum_across_tests", profile)

    def test_patient_packet_excludes_qc_only_and_preserves_case_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inventory = root / "inventory"
            patient_dir = inventory / "patients"
            patient_dir.mkdir(parents=True)
            clinical = candidate(observation(seq_id="S1", rank=2))
            clinical["organism_name"] = "Streptococcus pneumoniae"
            clinical["organism_key"] = "streptococcuspneumoniae"
            qc = candidate(observation(seq_id="S1", reads=0))
            document = {
                "patient_id": "9",
                "cases": [
                    {
                        "case_review_id": case_id,
                        "specimen_site": site,
                        "tests": [make_test("S1")],
                        "candidate_union": {"selected_union": [clinical], "detected_but_filtered_union": [qc]},
                    }
                    for case_id, site in (("A", "BALF"), ("B", "Blood"))
                ],
            }
            (patient_dir / "NGS_patient_9_multi_assay_mngs_evidence.json").write_text(json.dumps(document), encoding="utf-8")
            out = root / "out"
            summary = build_shadow(inventory, out)
            self.assertEqual(summary["case_organism_count"], 4)
            packet = json.loads((out / "patient_packets" / "NGS_patient_9_compact_entry_shadow.json").read_text(encoding="utf-8"))
            self.assertEqual([case["case_review_id"] for case in packet["cases"]], ["A", "B"])
            self.assertTrue(all(len(case["scorer_entries"]) == 1 for case in packet["cases"]))
            self.assertTrue(all(len(case["clinical_review"]) == 0 for case in packet["cases"]))
            self.assertTrue(all(len(case["qc_audit_refs"]) == 1 for case in packet["cases"]))
            with (out / "all_decisions.csv").open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([row["disposition"] for row in rows].count("qc_only"), 2)
            with self.assertRaises(FileExistsError):
                build_shadow(inventory, out)


if __name__ == "__main__":
    unittest.main()

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools import run_kh_integrated_release_v2 as release


class IntegratedReleaseV2Test(unittest.TestCase):
    def test_manifest_taxonomy_components_cover_all_runtime_layers(self):
        components = release._taxonomy_components()

        self.assertEqual(
            [item["role"] for item in components],
            [
                "central_taxonomy_base",
                "taxonomy_gap_resolution_overlay",
                "auto_taxonomy_overlay",
                "reviewed_taxonomy_overlay",
            ],
        )
        self.assertTrue(all(item["sha256"] for item in components))

    def test_current_source_summary_freeze_is_run_local_and_source_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            patient_root = root / "patients"
            output_root = root / "stage00"
            patient_dirs = [
                patient_root / f"NGS_patient_{patient_id}_json"
                for patient_id in range(1, 34)
            ]
            for patient_id, patient_dir in enumerate(patient_dirs, start=1):
                summary_dir = patient_dir / "summary_outputs"
                summary_dir.mkdir(parents=True)
                frozen = {
                    "host_context": {
                        "immunocompromise_tier": "H1",
                        "opportunistic_host_support": False,
                    },
                    "module_summaries": {"culture": {"preserved": True}},
                    "hospital_organism_evidence": [{"organism_name": "Preserved"}],
                    "evidence_preservation": {},
                    "source_files": {"culture": "frozen-culture.json"},
                }
                frozen_path = summary_dir / (
                    f"NGS_patient_{patient_id}_"
                    + release.FROZEN_HOSPITAL_SUMMARY_SUFFIX
                    + ".json"
                )
                frozen_path.write_text(json.dumps(frozen), encoding="utf-8")

            def fake_build_summary(patient_dir: Path):
                patient_id = int(patient_dir.name.split("_")[2])
                return {
                    "host_context": {
                        "immunocompromise_tier": "H2" if patient_id == 16 else "H1",
                        "opportunistic_host_support": False,
                    },
                    "module_summaries": {"cbc_other_lab": {"current": True}},
                    "evidence_preservation": {
                        "structured_host_evidence": {"current": True}
                    },
                    "source_files": {"raw_cbc": f"p{patient_id}-cbc.json"},
                }

            with (
                patch.object(
                    release.deterministic_summary,
                    "collect_patient_dirs",
                    return_value=patient_dirs,
                ),
                patch.object(
                    release.deterministic_summary,
                    "build_summary",
                    side_effect=fake_build_summary,
                ),
            ):
                summary = release.build_current_source_summaries(
                    patient_root, output_root
                )

            self.assertEqual(summary["patient_count"], 33)
            self.assertEqual(summary["immunocompromise_tier_counts"], {"H1": 32, "H2": 1})
            self.assertEqual(
                summary["opportunistic_host_support_counts"],
                {"true": 0, "false": 33},
            )
            self.assertEqual(summary["host_context_change_count"], 1)
            p16 = output_root / (
                "NGS_patient_16_" + release.CURRENT_SUMMARY_SUFFIX + ".json"
            )
            self.assertTrue(p16.is_file())
            self.assertEqual(
                json.loads(p16.read_text(encoding="utf-8"))["host_context"][
                    "immunocompromise_tier"
                ],
                "H2",
            )
            p16_payload = json.loads(p16.read_text(encoding="utf-8"))
            self.assertEqual(
                p16_payload["hospital_organism_evidence"],
                [{"organism_name": "Preserved"}],
            )
            self.assertTrue(p16_payload["module_summaries"]["culture"]["preserved"])
            self.assertTrue(
                p16_payload["release_input_overlay"][
                    "host_context_rebuilt_from_current_sources"
                ]
            )

    def test_current_source_summary_freeze_rejects_incomplete_cohort(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with patch.object(
                release.deterministic_summary,
                "collect_patient_dirs",
                return_value=[root / "NGS_patient_1_json"],
            ):
                with self.assertRaisesRegex(ValueError, "Expected 33"):
                    release.build_current_source_summaries(root, root / "output")


if __name__ == "__main__":
    unittest.main()

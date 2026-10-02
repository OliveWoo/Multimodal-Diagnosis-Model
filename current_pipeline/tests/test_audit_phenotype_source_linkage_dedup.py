import tempfile
import unittest
from pathlib import Path

from tools.audit_phenotype_source_linkage_dedup import (
    _cluster_map,
    _confirmed_episode_map,
    _patient_number,
    _sha256,
)


class PhenotypeSourceLinkageHelpersTest(unittest.TestCase):
    def test_patient_number_and_cluster_map(self):
        self.assertEqual(_patient_number("patient 33_phenotype.xlsx"), 33)
        self.assertEqual(_patient_number("patient_9_phenotype_context.json"), 9)
        self.assertIsNone(_patient_number("unrelated.xlsx"))
        mapping = _cluster_map([{"value": [9, 11], "Count": 2}])
        self.assertEqual(mapping[9], "semantic_cluster_1:9-11")
        self.assertEqual(mapping[11], "semantic_cluster_1:9-11")
        self.assertEqual(_cluster_map([[25, 26]])[26], "semantic_cluster_1:25-26")

    def test_sha256_is_content_based(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            first = Path(temp_dir) / "first.bin"
            second = Path(temp_dir) / "second.bin"
            first.write_bytes(b"same")
            second.write_bytes(b"same")
            self.assertEqual(_sha256(first), _sha256(second))

    def test_confirmed_episode_policy_keeps_same_patient_episodes_separate(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            policy = Path(temp_dir) / "episodes.json"
            policy.write_text(
                '{"clusters":[{"members":[11,9],'
                '"patient_identity":"same_patient_confirmed",'
                '"episode_relationship":"distinct_pneumonia_episodes_confirmed",'
                '"action":"keep_separate"}]}',
                encoding="utf-8",
            )
            by_member, by_cluster = _confirmed_episode_map(policy)
            self.assertEqual("same_patient_confirmed", by_member[9]["patient_identity"])
            self.assertEqual(
                "distinct_pneumonia_episodes_confirmed",
                by_cluster[(9, 11)]["episode_relationship"],
            )


if __name__ == "__main__":
    unittest.main()

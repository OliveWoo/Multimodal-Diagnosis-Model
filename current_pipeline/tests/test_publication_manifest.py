"""Verify the public packaging did not rewrite the decision algorithms."""
import ast
import hashlib
import json
import unittest
from pathlib import Path

from tools.validate_frozen_policy_bundle import validate

ROOT = Path(__file__).resolve().parents[1]


class PublicationManifestTests(unittest.TestCase):
    def test_core_decision_function_asts_match_the_source_logic(self):
        manifest = json.loads((ROOT / "publication_manifest_20261020.json").read_text(encoding="utf-8"))
        self.assertGreaterEqual(len(manifest["unchanged_function_contracts"]), 10)
        for contract in manifest["unchanged_function_contracts"]:
            with self.subTest(function=contract["function"]):
                module = ast.parse((ROOT / contract["file"]).read_text(encoding="utf-8"))
                node = next(n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == contract["function"])
                digest = hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest()
                self.assertEqual(digest, contract["ast_sha256"])

    def test_all_published_rule_hashes_match(self):
        manifest = json.loads((ROOT / "publication_manifest_20261020.json").read_text(encoding="utf-8"))
        for item in manifest["rule_files"]:
            with self.subTest(path=item["path"]):
                self.assertEqual(hashlib.sha256((ROOT / item["path"]).read_bytes()).hexdigest(), item["sha256"])

    def test_historical_bundle_and_new_bundle_both_validate(self):
        self.assertTrue(validate()["passed"])
        self.assertTrue(validate(ROOT / "rules/current_pipeline_policy_bundle_integrated_cmv_20261020.json")["passed"])
        manifest = json.loads((ROOT / "publication_manifest_20261020.json").read_text(encoding="utf-8"))
        original = ROOT / "rules/current_pipeline_policy_bundle_v6_v9_calendar_v7_frozen_20261013.json"
        self.assertEqual(hashlib.sha256(original.read_bytes()).hexdigest(), manifest["historical_20261013_bundle_sha256"])


if __name__ == "__main__":
    unittest.main()

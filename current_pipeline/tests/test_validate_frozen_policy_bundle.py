from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest

from tools.validate_frozen_policy_bundle import DEFAULT_BUNDLE, validate


class FrozenPolicyBundleTest(unittest.TestCase):
    def test_current_bundle_passes(self) -> None:
        result = validate(DEFAULT_BUNDLE)
        self.assertTrue(result["passed"])
        self.assertEqual(result["component_count"], 5)

    def test_hash_drift_fails(self) -> None:
        source = json.loads(DEFAULT_BUNDLE.read_text(encoding="utf-8"))
        source["components"][0]["sha256"] = "0" * 64
        with TemporaryDirectory() as temp:
            path = Path(temp) / "bundle.json"
            path.write_text(json.dumps(source), encoding="utf-8")
            with self.assertRaises(ValueError):
                validate(path)


if __name__ == "__main__":
    unittest.main()

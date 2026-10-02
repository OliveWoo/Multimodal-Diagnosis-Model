import unittest

from tools.audit_external_validation_compatibility import (
    _as_float,
    _biological_class,
    _infer_molecule,
)


class ExternalValidationCompatibilityHelpersTest(unittest.TestCase):
    def test_molecule_inference_is_conservative(self):
        self.assertEqual(_infer_molecule("folder/pipeline_DNA/file.tsv"), "DNA")
        self.assertEqual(_infer_molecule("folder/pipeline_RNA/file.tsv"), "RNA")
        self.assertEqual(_infer_molecule("folder/file.tsv"), "UNKNOWN")
        self.assertEqual(_infer_molecule("DNA_RNA/file.tsv"), "UNKNOWN")

    def test_source_class_and_numeric_parser(self):
        self.assertEqual(_biological_class("1.Bac"), "bacterium")
        self.assertEqual(_biological_class("2.Fungi"), "fungus")
        self.assertEqual(_biological_class("3.Virus"), "virus")
        self.assertEqual(_as_float("12"), 12.0)
        self.assertIsNone(_as_float(""))


if __name__ == "__main__":
    unittest.main()

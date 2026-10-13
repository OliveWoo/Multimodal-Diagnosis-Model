import json
import tempfile
import unittest
from pathlib import Path

from tools.run_historical_ober_r5_replay import load_cohort_specs, trigger_window


class HistoricalOberReplayTests(unittest.TestCase):
    def test_cohort_paths_are_relative_to_config_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / 'cohorts.json'
            config.write_text(json.dumps({'cohorts': [{'label': 'synthetic', 'packets': 'packets.jsonl', 'manifest': 'manifest.json'}]}), encoding='utf-8')
            specs = load_cohort_specs(config)
            self.assertEqual(specs[0].packets, root / 'packets.jsonl')
            self.assertEqual(specs[0].manifest, root / 'manifest.json')

    def test_duplicate_or_path_like_cohort_label_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / 'cohorts.json'
            for label in ('../outside', 'synthetic'):
                entries = [{'label': label, 'packets': 'p.jsonl', 'manifest': 'm.json'}] * 2
                config.write_text(json.dumps({'cohorts': entries}), encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'Invalid or duplicate'):
                    load_cohort_specs(config)

    def test_deterministic_span_is_bounded_and_traceable(self):
        abstract = 'The synthetic study describes infection at the stated site with limited available evidence.'
        excerpt, valid = trigger_window('', abstract, ['infection'])
        self.assertTrue(valid)
        self.assertIn(excerpt, abstract)
        self.assertLessEqual(len(excerpt.split()), 20)

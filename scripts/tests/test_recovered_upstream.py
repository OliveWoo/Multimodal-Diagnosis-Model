"""Synthetic version and replay integrity checks; no clinical records or gold."""
from pathlib import Path
import contextlib, copy, hashlib, io, json, sys, tempfile, unittest

SCRIPTS = Path(__file__).resolve().parents[1]
ROOT = SCRIPTS.parent
sys.path[:0] = [str(SCRIPTS), str(SCRIPTS/'offline_guard'), str(ROOT/'upstream')]
import sitecustomize
import replay_v19_upstream as replay
from tools import deterministic_mngs_max_scorer as current
from tools import deterministic_mngs_max_scorer_v19_recovered as historical


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


class RecoveredUpstreamTests(unittest.TestCase):
    def test_general_high_signal_hospital_rule_distinguishes_versions(self):
        folder = ROOT/'examples/upstream_v19'
        ranked = read(folder/'ranked.synthetic.json')
        summary = read(folder/'summary.synthetic.json')
        args = dict(patient_dir=Path('NGS_patient_900001_json'), ranked_mngs=ranked, final_summary=summary)
        v19 = historical.score_payload(**copy.deepcopy(args))['pathogen_candidates'][0]
        v20 = current.score_payload(**copy.deepcopy(args))['pathogen_candidates'][0]
        self.assertEqual((v19['mngs_signal_tier'], v19['integrated_causative_level']), ('M2_moderate', 'Level 2'))
        self.assertEqual((v20['mngs_signal_tier'], v20['integrated_causative_level']), ('M1_strong', 'Level 1'))
        self.assertNotIn('taxonomy_profile', v19)
        self.assertIn('taxonomy_profile', v20)

    def test_synthetic_replay_and_refusal_to_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'replay'
            with contextlib.redirect_stdout(io.StringIO()):
                status = replay.replay(ROOT/'examples/upstream_v19/manifest.synthetic.json', output)
            self.assertEqual(status, 0)
            before = (output/'receipt.json').read_bytes()
            with self.assertRaisesRegex(ValueError, 'already exists'):
                replay.replay(ROOT/'examples/upstream_v19/manifest.synthetic.json', output)
            self.assertEqual(before, (output/'receipt.json').read_bytes())

    def test_changed_input_hash_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'input.json'
            original = b'{"records": []}'
            path.write_bytes(original+b' ')
            with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
                replay.read_reference({'path':path.name, 'sha256':hashlib.sha256(original).hexdigest()}, path.parent)

    def test_missing_input_hash_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'SHA256'):
            replay.read_reference({'path':'absent.json'}, Path('.'))

    def test_patient_mismatch_is_refused(self):
        with self.assertRaisesRegex(ValueError, 'patient ID mismatch'):
            replay.ranked_for_patient({'records':[], 'patient_id':'900002'}, '900001', historical.mngs)

    def test_changed_rules_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'docs').mkdir()
            (root/'rules.json').write_text('{}', encoding='utf-8')
            profile = {'profile_id':'recovered-v19-20260918', 'files':[{'path':'rules.json', 'sha256':'0'*64}]}
            (root/replay.PROFILE).write_text(json.dumps(profile), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'Code/rules changed'):
                replay.verify_profile(root)

    def test_comparison_preserves_level_order_reason_and_version(self):
        a = {'rule_version':'v19', 'items':[{'level':1, 'reason':'x'}, {'level':2}], 'source_files':{'file':'old'}}
        b = copy.deepcopy(a)
        b['source_files']['file'] = 'relocated'
        self.assertEqual(replay.sans_source_paths(a), replay.sans_source_paths(b))
        for key, value in [('rule_version','v20'), ('items', list(reversed(a['items']))), ('items',[{'level':1,'reason':'changed'}, {'level':2}])]:
            changed = copy.deepcopy(a);changed[key] = value
            self.assertNotEqual(replay.sans_source_paths(a), replay.sans_source_paths(changed))

    def test_duplicate_case_refused_before_output_creation(self):
        folder = ROOT/'examples/upstream_v19'
        manifest = read(folder/'manifest.synthetic.json')
        for label in ['ranked','summary','expected_scorer']:
            manifest['cases'][0][label]['path'] = str(folder/manifest['cases'][0][label]['path'])
        manifest['cases'].append(copy.deepcopy(manifest['cases'][0]))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'manifest.json';path.write_text(json.dumps(manifest),encoding='utf-8')
            output = Path(tmp)/'output'
            with self.assertRaisesRegex(ValueError, 'duplicate case'):
                replay.replay(path, output)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()

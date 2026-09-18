"""Synthetic checks for incorrect identities, frozen image bindings and routes."""
from pathlib import Path
from types import SimpleNamespace
import copy
import tempfile
import unittest
import sys
import subprocess
import os

from pipeline import (Inputs, FrozenVision, compare_expected, identity_rows,
                      radiology_route, read, run, safe_name, semantic_sha, sha,
                      validate_identity, write, SHARED)
sys.path.insert(0, str(SHARED))


class FrontendGuards(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def rows(self):
        return [dict(patient_key='SYN_A', case_code='S001', specimen_id='SYN_001', hospital_id='900001'),
                dict(patient_key='SYN_B', case_code='S002', specimen_id='SYN_002', hospital_id='900001')]

    def test_same_person_different_specimens_are_not_duplicate(self):
        validate_identity(self.rows())

    def test_patient_order_is_locked(self):
        rows = self.rows()
        with self.assertRaisesRegex(ValueError, 'identity/order'):
            validate_identity(list(reversed(rows)), semantic_sha(identity_rows(rows)))

    def test_duplicate_specimen_identity_is_rejected(self):
        rows = self.rows(); rows.append(dict(rows[0], patient_key='DIFFERENT_FILENAME'))
        with self.assertRaisesRegex(ValueError, 'Duplicate specimen'):
            validate_identity(rows)

    def test_duplicate_patient_folder_is_rejected(self):
        rows = self.rows(); rows[1]['patient_key'] = rows[0]['patient_key']
        with self.assertRaisesRegex(ValueError, 'Duplicate patient key'):
            validate_identity(rows)

    def test_changed_file_is_rejected(self):
        p = self.root / 'input.json'; write(p, [1]); ref = {'path': p.name, 'sha256': sha(p)}
        write(p, [2])
        with self.assertRaisesRegex(ValueError, 'SHA256'):
            Inputs(self.root).file(ref)

    def test_change_during_run_is_rejected(self):
        p = self.root / 'input.json'; write(p, [1]); inputs = Inputs(self.root)
        inputs.file({'path': p.name, 'sha256': sha(p)}); write(p, [2])
        with self.assertRaisesRegex(ValueError, 'changed during'):
            inputs.unchanged()

    def test_existing_output_is_never_overwritten(self):
        with self.assertRaisesRegex(ValueError, 'already exists'):
            run(self.root / 'missing.json', self.root)

    def test_output_name_cannot_escape(self):
        for name in ['../patient.json', 'C:\\patient.json', '.', '..', '/patient.json']:
            with self.assertRaises(ValueError): safe_name(name)

    def test_image_routing_is_restored_on_failure(self):
        original = lambda payload: ['radiology']
        linker = SimpleNamespace(extract_radiology_entries=original)
        with self.assertRaises(RuntimeError):
            with radiology_route(linker, 'labs_only'):
                self.assertEqual(linker.extract_radiology_entries({}), [])
                raise RuntimeError('stop')
        self.assertIs(linker.extract_radiology_entries, original)

    def test_full_routing_preserves_image_content(self):
        linker = SimpleNamespace(extract_radiology_entries=lambda payload: ['radiology'])
        with radiology_route(linker, 'radiology_and_labs'):
            self.assertEqual(linker.extract_radiology_entries({}), ['radiology'])

    def test_unknown_routing_is_rejected(self):
        with self.assertRaises(ValueError):
            with radiology_route(SimpleNamespace(), 'auto'): pass

    def frozen(self):
        parsed = self.root / 'saved/slides/slide_001/parsed.json'
        image = parsed.parent / 'images/raw_01.png'; image.parent.mkdir(parents=True)
        image.write_bytes(b'synthetic-image-bytes')
        value = dict(source_file='synthetic.pptx', slide_index=1, images=[dict(image_index=1,
                     raw_path='old/machine/images/raw_01.png', vision=dict(status='success',
                     raw_text=['Synthetic text'], findings=[], confidence=1.0, model_used='synthetic-fixture'))])
        write(parsed, value)
        return parsed, image, value

    def test_wrong_ppt_cache_is_rejected(self):
        parsed, _, _ = self.frozen()
        with self.assertRaisesRegex(ValueError, 'filename mismatch'):
            FrozenVision(Path('another.pptx'), [{'path': str(parsed), 'sha256': sha(parsed)}], Inputs(self.root))

    def test_failed_model_cache_is_not_accepted(self):
        parsed, _, value = self.frozen(); value['images'][0]['vision']['status'] = 'error'; write(parsed, value)
        with self.assertRaisesRegex(ValueError, 'Unsuccessful'):
            FrozenVision(Path('synthetic.pptx'), [{'path': str(parsed), 'sha256': sha(parsed)}], Inputs(self.root))

    def test_duplicate_frozen_slides_are_rejected(self):
        parsed, _, _ = self.frozen(); ref = {'path': str(parsed), 'sha256': sha(parsed)}
        with self.assertRaisesRegex(ValueError, 'Duplicate/invalid'):
            FrozenVision(Path('synthetic.pptx'), [ref, ref], Inputs(self.root))

    def test_changed_image_bytes_cannot_reuse_cached_response(self):
        parsed, _, _ = self.frozen()
        frozen = FrozenVision(Path('synthetic.pptx'), [{'path': str(parsed), 'sha256': sha(parsed)}], Inputs(self.root))
        image = self.root / 'new/slides/slide_001/images/raw_01.png'
        image.parent.mkdir(parents=True); image.write_bytes(b'changed image')
        with self.assertRaisesRegex(ValueError, 'bytes differ'):
            frozen.extract(image)

    def test_missing_image_response_fails_without_network(self):
        frozen = FrozenVision(Path('synthetic.pptx'), [], Inputs(self.root))
        with self.assertRaisesRegex(ValueError, 'No frozen response'):
            frozen.extract(self.root / 'new/slides/slide_001/images/raw_01.png')

    def test_exact_image_binding_returns_saved_text(self):
        parsed, image, _ = self.frozen()
        frozen = FrozenVision(Path('synthetic.pptx'), [{'path': str(parsed), 'sha256': sha(parsed)}], Inputs(self.root))
        self.assertEqual(frozen.extract(image).raw_text, ['Synthetic text'])
        with self.assertRaisesRegex(ValueError, 'reused'):
            frozen.extract(image)

    def test_semantic_comparison_does_not_ignore_clinical_changes(self):
        expected = self.root / 'expected.json'; write(expected, [{'value': 1}])
        out = self.root / 'new'; write(out / 'result.json', [{'value': 2}])
        with self.assertRaisesRegex(ValueError, 'differs'):
            compare_expected([{'output': 'result.json', 'reference': {'path': str(expected), 'sha256': sha(expected)}}], Inputs(self.root), out)
        self.assertFalse(read(out / 'comparison.private.json')[0]['equal'])

    def test_comparison_cannot_read_outside_output(self):
        out = self.root / 'new'; out.mkdir()
        with self.assertRaisesRegex(ValueError, 'escapes'):
            compare_expected([{'output': '../outside.json'}], Inputs(self.root), out)

    def test_json_formatting_is_not_a_content_change(self):
        self.assertEqual(semantic_sha({'b': 2, 'a': 1}), semantic_sha({'a': 1, 'b': 2}))
        self.assertNotEqual(semantic_sha([1, 2]), semantic_sha([2, 1]))

    def test_public_excel_ppt_example_reaches_standardized_content(self):
        repo = SHARED.parents[1]
        output = self.root / 'run'
        env = dict(os.environ, PYTHONIOENCODING='utf-8', PYTHONUTF8='1')
        result = subprocess.run([sys.executable, '-B', str(repo / 'frontend/pipeline.py'),
            '--manifest', str(repo / 'examples/frontend/manifest.synthetic.json'), '--output', str(output)],
            cwd=self.root, env=env, capture_output=True, text=True, encoding='utf-8', timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt = read(output / 'receipt.private.json')
        self.assertEqual((receipt['cases'], receipt['slides'], receipt['unmatched_slides'], receipt['fresh_model_calls']), (1, 1, 0, 0))
        patient = output / 'standardized/NGS_patient_1_json'
        mngs = read(patient / 'NGS_patient_1_mNGS_grouped.json')
        self.assertEqual(mngs[0]['pathogens']['bacterial'], [{'name': 'Synthetic bacterium', 'reads': 120}])
        radiology = read(patient / 'NGS_patient_1_image.json')
        self.assertEqual(len(radiology), 1)
        self.assertIn('fictional record', str(radiology[0]['findings']))


if __name__ == '__main__':
    unittest.main()

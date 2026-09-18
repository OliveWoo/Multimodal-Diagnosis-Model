"""Cohort-level identity and real CLI tests using wholly synthetic JSON."""
import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from cohorts import PROFILES, ROOT, run, validate_binding, validate_anchor, check_filename
from pipeline import read, write


class CohortIsolation(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.examples = ROOT / 'examples/frontend/cohorts'
        self.kh = read(self.examples / 'KH.synthetic.json')
        self.two = read(self.examples / 'two_hospital.synthetic.json')
        self.index = {r['legacy_patient_id'].replace('NGS_patient_', 'P'): r
                      for r in read(self.examples / 'source_index.synthetic.json')}

    def test_kh_local_id_can_differ_from_global_ordinal(self):
        row = validate_binding('KH', self.kh['cases'][0], self.index)
        self.assertEqual(row['legacy_patient_id'], 'NGS_patient_8')
        self.assertEqual(self.kh['cases'][0]['patient_id'], 'P3')

    def test_kh_same_number_join_is_not_assumed(self):
        case = copy.deepcopy(self.kh['cases'][0]); case['master_patient_id'] = 'P3'
        with self.assertRaisesRegex(ValueError, 'missing from pinned index'):
            validate_binding('KH', case, self.index)

    def test_kh_case_code_cannot_follow_wrong_clinical_id(self):
        case = copy.deepcopy(self.kh['cases'][0]); case['patient_id'] = 'P8'
        with self.assertRaisesRegex(ValueError, 'K case code'):
            validate_binding('KH', case, self.index)

    def test_two_hospital_keeps_global_id(self):
        validate_binding('two_hospital', self.two['cases'][0], self.index)
        case = copy.deepcopy(self.two['cases'][0]); case['patient_id'] = 'P1'
        with self.assertRaisesRegex(ValueError, 'must not be renumbered'):
            validate_binding('two_hospital', case, self.index)

    def test_kh_case_cannot_enter_two_hospital(self):
        with self.assertRaisesRegex(ValueError, 'Hospital group'):
            validate_binding('two_hospital', self.kh['cases'][0], self.index)

    def test_wrong_specimen_rejected_even_if_patient_number_matches(self):
        case = copy.deepcopy(self.two['cases'][0]); case['specimen_code'] = 'OTHER'
        with self.assertRaisesRegex(ValueError, 'Specimen differs'):
            validate_binding('two_hospital', case, self.index)

    def test_flat_rk_ntc_is_not_silently_used_as_grouped_mngs(self):
        with self.assertRaisesRegex(ValueError, 'Flat RK/NTC'):
            validate_anchor('two_hospital', self.two['cases'][0], [{'name':'Synthetic bacterium','rk_ntc':[]}])

    def test_empty_mngs_requires_explicit_status(self):
        case = copy.deepcopy(self.two['cases'][0])
        with self.assertRaisesRegex(ValueError, 'Empty specimen'):
            validate_anchor('two_hospital', case, [])
        case['anchor_status'] = 'empty_explicit'
        self.assertIn('index_only', validate_anchor('two_hospital', case, []))

    def test_kh_cannot_omit_specimen_anchor(self):
        case = copy.deepcopy(self.kh['cases'][0]); case['anchor_status'] = 'empty_explicit'
        with self.assertRaises(ValueError): validate_anchor('KH', case, [])

    def test_model_free_specimen_check(self):
        case = self.kh['cases'][0]
        with self.assertRaisesRegex(ValueError, 'does not match'):
            validate_anchor('KH', case, [{'specimen_code':'SYN-WRONG','pathogens':{}}])

    def test_numbered_duplicate_filename_is_explicitly_supported(self):
        check_filename(Path('NGS_patient_54_CBC (1).json'), 'P54', 'CBC')
        with self.assertRaisesRegex(ValueError, 'different local patient'):
            check_filename(Path('NGS_patient_8_CBC.json'), 'P3', 'CBC')

    def test_cross_cohort_manifest_fails_before_writing(self):
        out = self.root / 'wrong'
        with self.assertRaisesRegex(ValueError, 'other hospital group'):
            run('two_hospital', self.examples / 'KH.synthetic.json', out)
        self.assertFalse(out.exists())

    def test_existing_output_is_never_overwritten(self):
        with self.assertRaisesRegex(ValueError, 'already exists'):
            run('KH', self.examples / 'KH.synthetic.json', self.root)

    def cli(self, cohort, entry):
        output = self.root / cohort
        env = dict(os.environ, PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
        p = subprocess.run([sys.executable, '-B', str(ROOT/f'frontend/{entry}/run.py'), '--manifest',
                            str(self.examples/f'{cohort}.synthetic.json'), '--output', str(output)],
                           cwd=self.root, env=env, capture_output=True, text=True, encoding='utf-8', timeout=30)
        self.assertEqual(p.returncode, 0, p.stderr)
        receipt = read(output/'receipt.private.json')
        self.assertFalse(receipt['full_raw_to_ober_reproduced'])
        self.assertEqual(receipt['fresh_model_calls'], 0)
        return output

    def test_kh_cli_preserves_local_name_and_global_mapping(self):
        output = self.cli('KH', 'kh')
        self.assertTrue((output/'KH/patient_info/NGS_patient_3_json/NGS_patient_3_CBC.json').exists())
        self.assertFalse((output/'KH/patient_info/NGS_patient_8_json').exists())
        mapping = read(output/'KH/identity_mapping.private.json')[0]
        self.assertEqual((mapping['patient_id'], mapping['master_patient_id']), ('P3','P8'))

    def test_two_hospital_cli_preserves_global_name(self):
        output = self.cli('two_hospital', 'two_hospital')
        self.assertTrue((output/'two_hospital/patient_info/NGS_patient_54_json/NGS_patient_54_CBC.json').exists())
        self.assertFalse((output/'two_hospital/patient_info/NGS_patient_1_json').exists())


if __name__ == '__main__': unittest.main()

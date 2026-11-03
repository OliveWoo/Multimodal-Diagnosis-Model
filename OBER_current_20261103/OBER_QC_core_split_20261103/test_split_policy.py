import copy
import unittest

from split_policy import CORE_FIELDS, core_profile, split_consensus


def record(qualified=False):
    value = dict(scope='exact_species', human_any=True, qualified=qualified,
                 basis='direct' if qualified else 'none', nonhuman_only=False,
                 noncausal=False, secondary_other=False, status='resolved',
                 internal_issues=[], source_summaries={}, signature=[qualified, False],
                 human_ids=['S2'], qualifying_ids=['S1'], identity_ids=['S1'])
    return value


def pending():
    return dict(status='unresolved_qc', qualified=None, basis='uncertain')


class SplitTests(unittest.TestCase):
    def test_recovers_only_negative_core_and_keeps_details(self):
        left, right = record(), record()
        right['noncausal'] = True
        result = split_consensus(left, right, pending())
        self.assertIs(result['qualified'], False)
        self.assertEqual(result['status'], 'qc_core_consensus_negative')
        self.assertTrue(result['detail_QC_pending'])
        self.assertIsNone(result['noncausal'])
        self.assertEqual(result['detail_QC_pair'][1]['noncausal'], True)

    def test_does_not_resolve_positive_auxiliary_disagreement(self):
        left, right = record(True), record(True)
        right['noncausal'] = True
        self.assertEqual(split_consensus(left, right, pending())['status'], 'unresolved_qc')

    def test_does_not_resolve_any_core_disagreement(self):
        for key in CORE_FIELDS:
            with self.subTest(key=key):
                left, right = record(), record()
                right[key] = not left[key] if isinstance(left[key], bool) else 'different'
                self.assertEqual(split_consensus(left, right, pending())['status'], 'unresolved_qc')

    def test_does_not_resolve_single_unresolved_or_contradiction(self):
        for change in [dict(status='unresolved_qc'), dict(internal_issues=['contradiction'])]:
            left, right = record(), record()
            right.update(change)
            self.assertEqual(split_consensus(left, right, pending())['status'], 'unresolved_qc')

    def test_preserves_existing_positive(self):
        left, right = record(True), record(True)
        baseline = copy.deepcopy(left)
        baseline['status'] = 'qc_consensus'
        result = split_consensus(left, right, baseline)
        for key in baseline:
            self.assertEqual(result[key], baseline[key])
        self.assertFalse(result['core_recovered_from_detail_disagreement'])

    def test_no_input_mutation_and_symmetric_core(self):
        left, right = record(), record()
        right['noncausal'] = True
        right['human_ids'] = ['S3']
        originals = copy.deepcopy([left, right])
        a = split_consensus(left, right, pending())
        b = split_consensus(right, left, pending())
        self.assertEqual([left, right], originals)
        self.assertEqual(core_profile(a), core_profile(b))
        self.assertEqual(a['human_ids'], ['S2', 'S3'])

    def test_unknown_detail_is_not_converted_to_absent(self):
        left, right = record(), record()
        right['noncausal'] = None
        result = split_consensus(left, right, pending())
        self.assertIsNone(result['noncausal'])
        self.assertTrue(result['detail_QC_pending'])


if __name__ == '__main__':
    unittest.main()

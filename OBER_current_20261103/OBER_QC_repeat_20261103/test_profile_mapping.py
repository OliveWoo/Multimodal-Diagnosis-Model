import unittest
import repeat_entry as e


class MappingTests(unittest.TestCase):
    def test_actual_normalized_scope_is_read(self):
        r={'status':'qc_consensus','scope':'exact_species','human_any':True,'qualified':True,'basis':'direct',
            'nonhuman_only':False,'noncausal':False,'secondary_other':False}
        self.assertEqual(e.profile(r)[1],'exact_species')
    def test_scope_only_change_detected(self):
        a={'status':'qc_consensus','scope':'genus_only','human_any':False,'qualified':False,'basis':'none'}
        b=a|{'scope':'other'};self.assertNotEqual(e.profile(a),e.profile(b))
    def test_pending_remains_unknown_not_resolved_identity(self):
        self.assertEqual(e.profile({'status':'unresolved_qc','scope':'unclear'}),['unresolved_qc'])


if __name__=='__main__':unittest.main()

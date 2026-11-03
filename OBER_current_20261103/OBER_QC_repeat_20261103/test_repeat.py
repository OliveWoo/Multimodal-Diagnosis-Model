import unittest
import copy
import compare_runs as c
import repeat_driver as d


class RepeatTests(unittest.TestCase):
    def test_runs_have_distinct_paths(self):
        self.assertNotEqual(d.ROUNDS['A'],d.ROUNDS['B']);self.assertNotEqual(d.PARENT,d.ROUNDS['A']);self.assertNotEqual(d.PARENT,d.ROUNDS['B'])
    def test_parent_path_not_mutable_context(self):
        before=d.PARENT;old=d.s.RUN
        try:d.bind('B');self.assertEqual(d.PARENT,before)
        finally:d.s.RUN=old;d.s.c.RUN=old
    def test_unknown_is_not_resolved_agreement(self):
        self.assertEqual(c.profile({'status':'unresolved_qc'}),['unresolved_qc'])
        resolved={'status':'qc_consensus',**dict(zip(c.SIGNATURE_FIELDS,['exact_species',True,True,'direct',False,False,False]))}
        self.assertNotEqual(c.profile(resolved),['unresolved_qc'])
    def test_wording_not_semantic_difference(self):
        r={'status':'qc_consensus',**dict(zip(c.SIGNATURE_FIELDS,['exact_species',True,True,'direct',False,False,False]))}
        other=copy.deepcopy(r);other['source_summaries']={'x':'Different wording.'};self.assertEqual(c.profile(r),c.profile(other))
    def test_changed_basis_detected(self):
        a=['exact_species',True,True,'direct',False,False,False];b=copy.deepcopy(a);b[3]='disease_related'
        self.assertEqual(c.changes(a,b),['basis'])
    def test_pair_sets_not_global_organism_sets(self):
        a=[{'case_id':'P1','organism':'X','decisions':{'k2':True}},{'case_id':'P2','organism':'X','decisions':{'k2':False}}]
        b=copy.deepcopy(a);b[0]['decisions']['k2']=False;b[1]['decisions']['k2']=True
        self.assertNotEqual(c.selected(a,'k2'),c.selected(b,'k2'));self.assertEqual(len(c.selected(a,'k2')^c.selected(b,'k2')),2)


if __name__=='__main__':unittest.main()

import copy
import unittest
import qc_common as s
import qc_rules as r
import qc_review as q


class QCTests(unittest.TestCase):
    def setUp(self):
        self.p={'article_key':'test','organism':'Example microbe','documented_clinical_names_not_scoring_aliases':[],
            'publication_types':['Review'],'article':{'pmid':'1','title':'Title','abstract':'Example microbe causes pneumonia in human patients. Disease D has pneumonia as a complication.',
                'sentences':[{'id':'S1','start':0,'end':5,'text':'Title'},
                    {'id':'S2','start':6,'end':60,'text':'Example microbe causes pneumonia in human patients.'},
                    {'id':'S3','start':61,'end':120,'text':'Disease D has pneumonia as a complication.'}]}}
        check=lambda status,ids:{'status':status,'sentence_ids':ids,'summary':'Synthetic test.'}
        self.q={'article_key':'test','organism_scope':'exact_species','identity_sentence_ids':['S2'],
            'study_design':'clinical_review','study_design_sentence_ids':[],
            'human_any':check('supported',['S2']),
            'direct_pulmonary':check('unclear',[])|{'level':'unclear'},
            'disease_related':{'status':'supported','chain':'pathogen_disease_lung','agent_disease_sentence_ids':['S2'],
                'disease_lung_sentence_ids':['S3'],'summary':'Synthetic disease-associated manifestation.'},
            'noncausal':check('absent',[]),'secondary_other':check('absent',[])}

    def test_complication_not_automatically_zero(self):
        x=s.normalize(self.q,self.p);self.assertTrue(x['qualified']);self.assertEqual(x['basis'],'disease_related')
    def test_unknown_direct_mechanism_not_veto(self):self.assertEqual(s.normalize(self.q,self.p)['status'],'resolved')
    def test_explicit_secondary_not_original_direct_support(self):
        self.q['disease_related']['chain']='secondary_other_agent';self.q['direct_pulmonary'].update(status='absent',level='none')
        self.assertFalse(s.normalize(self.q,self.p)['qualified'])
    def test_unknown_chain_retained_QC(self):
        self.q['disease_related']['chain']='unclear';self.q['direct_pulmonary'].update(status='absent',level='none')
        self.assertEqual(s.normalize(self.q,self.p)['status'],'unresolved_qc')
    def test_two_QC_disagree_stays_unresolved(self):
        a=s.normalize(self.q,self.p);different=copy.deepcopy(a);different['signature'][3]='direct'
        self.assertEqual(s.consensus(a,different)['status'],'unresolved_qc')
    def test_two_QC_agree_accept_without_requiring_same_wording(self):
        a=s.normalize(self.q,self.p);b=copy.deepcopy(a);b['source_summaries']['human_any']='Another concise wording.'
        self.assertEqual(s.consensus(a,b)['status'],'qc_consensus')
    def test_internal_contradiction_QC_not_silent_accept(self):
        self.q['human_any']['status']='absent';self.q['human_any']['sentence_ids']=[]
        self.assertEqual(s.normalize(self.q,self.p)['status'],'unresolved_qc')
    def test_case_only_not_mature_etiology(self):
        self.q['disease_related'].update(status='absent',chain='not_applicable')
        self.q['direct_pulmonary'].update(status='supported',level='case_only',sentence_ids=['S2'])
        self.assertFalse(s.normalize(self.q,self.p)['qualified'])
    def test_no_cross_paper_source(self):
        self.q['human_any']['sentence_ids']=['OTHER']
        with self.assertRaises(ValueError):q.validate(self.q,self.p)
    def test_title_only_clinical_claim_rejected(self):
        self.q['human_any']['sentence_ids']=['S1']
        with self.assertRaises(ValueError):q.validate(self.q,self.p)
    def test_both_links_required_for_complication(self):
        self.q['disease_related']['agent_disease_sentence_ids']=[]
        with self.assertRaises(ValueError):q.validate(self.q,self.p)
    def test_valid_source_ID_schema(self):q.validate(self.q,self.p)
    def test_metadata_design_locked(self):
        self.q['study_design']='single_case'
        with self.assertRaises(ValueError):q.validate(self.q,self.p)
    def test_rule_positive_explicit_human_cause_sentence(self):
        self.p['article']['abstract']='Example microbe causes pneumonia in human patients.'
        self.p['article']['sentences']=self.p['article']['sentences'][:2]
        self.assertEqual(r.route(self.p)['route'],'rule')
    def test_complication_routed_QC_not_dismissed(self):self.assertEqual(r.route(self.p)['route'],'QC')
    def test_no_keyword_is_QC_not_negative(self):
        self.p['article']['abstract']='A general biological description.';self.p['article']['sentences']=self.p['article']['sentences'][:1]
        self.assertEqual(r.route(self.p)['route'],'QC')
    def test_nonhuman_is_QC_not_auto_human_support(self):
        self.p['article']['abstract']='Example microbe causes pneumonia in human patients. An in vitro experiment.'
        self.assertEqual(r.route(self.p)['route'],'QC')
    def test_negation_not_auto_positive(self):
        self.p['article']['abstract']='Example microbe does not cause pneumonia in human patients.'
        self.p['article']['sentences'][1]['text']=self.p['article']['abstract']
        self.assertEqual(r.route(self.p)['route'],'QC')
    def test_no_abstract_unavailable(self):
        self.p['article']['abstract']='';self.assertEqual(r.route(self.p)['route'],'unavailable')
    def test_rule_deterministic(self):self.assertEqual(r.route(self.p),r.route(copy.deepcopy(self.p)))
    def test_normalization_not_modify_raw_response(self):
        before=copy.deepcopy(self.q);s.normalize(self.q,self.p);self.assertEqual(before,self.q)
    def test_meta_not_proof_of_causality(self):
        self.q['direct_pulmonary'].update(status='absent',level='none');self.q['disease_related'].update(status='absent',chain='not_applicable')
        self.assertFalse(s.normalize(self.q,self.p)['qualified'])


if __name__=='__main__':unittest.main()

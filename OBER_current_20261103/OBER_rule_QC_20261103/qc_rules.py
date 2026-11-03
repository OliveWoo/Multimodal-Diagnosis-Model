"""Positive-only high-precision routing, never absence-of-keyword medical rejection."""
import re
import qc_common as s

LUNG=re.compile(r'\b(?:bronchopneumonia|pneumonia|pneumonitis|(?:lower\s+respiratory(?:\s+tract)?|pulmonary|lung)\s+infections?)\b',re.I)
HUMAN=re.compile(r'\b(?:humans?|patients?|people|adults?|children|infants?)\b',re.I)
COMPLEX=re.compile(r'\b(?:complications?|secondary|superinfection|colonization|colonisation|reactivation|bystander|contamination|controversial)\b',re.I)
NEGATED=re.compile(r'\b(?:not|no|never|uncertain|unclear|might|whether|possibly|potentially|association|associated|detected|detection|isolated)\b',re.I)
NONHUMAN=re.compile(r'\b(?:mice|mouse|murine|rats?|canine|dogs?|cats?|veterinary|in\s+vitro|animal\s+model)\b',re.I)


def names(p):
    values=[p['organism']]+p.get('documented_clinical_names_not_scoring_aliases',[])
    words=p['organism'].split()
    if len(words)==2 and all(w.isalpha() for w in words):values += [words[0][0]+'. '+words[1],words[0][0]+'.'+words[1]]
    return sorted(set(values))


def route(p):
    a=p['article']
    if not a.get('abstract'):return {'route':'unavailable','reasons':['no_abstract'],'rule_record':None}
    reasons=[]
    design=s.metadata_design(p.get('publication_types',[]))
    if design is None:reasons.append('design_not_resolved_by_publication_metadata')
    text=a['abstract']
    if COMPLEX.search(text):reasons.append('mixed_or_complication_or_secondary_context')
    if NONHUMAN.search(text):reasons.append('human_nonhuman_scope_needs_review')
    matches=[]
    for sentence in a['sentences']:
        if sentence['id'] not in s.body_ids(p):continue
        t=sentence['text']
        if not LUNG.search(t) or not HUMAN.search(t) or NEGATED.search(t):continue
        for name in names(p):
            n=r'(?<!\w)'+re.escape(name).replace(r'\ ',r'\s+')+r'(?!\w)'
            forward=n+r'[^.;]{0,140}\b(?:causes?|causing|caused)\b[^.;]{0,120}'+LUNG.pattern
            backward=LUNG.pattern+r'[^.;]{0,120}\b(?:caused\s+by|due\s+to|attributed\s+to)\b[^.;]{0,100}'+n
            if re.search(forward,t,re.I) or re.search(backward,t,re.I):matches.append(sentence['id']);break
    if not matches:reasons.append('no_unambiguous_literal_human_pulmonary_cause_sentence')
    if reasons:return {'route':'QC','reasons':reasons,'rule_record':None}
    first=list(dict.fromkeys(matches));level='case_only' if design=='single_case' else 'recognized'
    claim=lambda status,ids,summary:{'status':status,'sentence_ids':ids,'summary':summary}
    q={'article_key':p['article_key'],'organism_scope':'exact_species','identity_sentence_ids':first,'study_design':design,'study_design_sentence_ids':[],
        'human_any':claim('supported',first,'規則定位到同句人類／患者及候選致病陳述。'),
        'direct_pulmonary':claim('supported',first,'規則定位到明確候選—致病連接詞—肺炎陳述。')|{'level':level},
        'disease_related':{'status':'absent','chain':'not_applicable','agent_disease_sentence_ids':[],
            'disease_lung_sentence_ids':[],'summary':'未由規則建立另外的疾病—併發症鏈。'},
        'noncausal':claim('absent',[],'自動規則未遇到須另外解析的非因果詞。'),
        'secondary_other':claim('absent',[],'自動規則未遇到繼發感染詞。')}
    record=s.normalize(q,p);record['status']='accepted_rule';record['raw_rule_output']=q
    return {'route':'rule','reasons':['literal_same_sentence_with_metadata_and_no_flagged_ambiguity'],'rule_record':record}

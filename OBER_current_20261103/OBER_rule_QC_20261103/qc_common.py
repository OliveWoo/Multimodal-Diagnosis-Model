"""Isolated rule/QC experiment; no mutation of earlier projects/artifacts."""
import copy
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'OBER_article_only_isolation_20261102'))
import iso_common as old
c=old.c;a=old.a;judge=old.judge;olda1=old.olda1
RUN=HERE/'runs/qc_v1';PARENT=old.RUN
c.RUN=RUN;c.WORKERS=20
load,save,require,digest,write_new=old.load,old.save,old.require,old.digest,old.write_new
obj,enum,ids,string=old.obj,old.enum,old.ids,old.string
MODEL=old.MODEL
ASPECTS=['human_any','direct_pulmonary','disease_related','noncausal','secondary_other']


def freeze(stage,paths):
    save(RUN/(stage+'_freeze.json'),{'files':{str(p.resolve()):digest(p.read_bytes()) for p in sorted(set(paths))},
        'historically_viewed_development_data':True,'independent_test':False})


def verify(stage):
    for p,h in load(RUN/(stage+'_freeze.json'))['files'].items():require(digest(Path(p).read_bytes())==h,'Changed frozen artifact '+p)


def response_files(folder):return old.response_files(folder)


def metadata_design(types):
    values=set(types)
    if values&{'Practice Guideline','Guideline','Consensus Statement'}:return 'guideline_consensus'
    if values&{'Systematic Review','Meta-Analysis','Scoping Review'}:return 'systematic_review'
    if 'Review' in values:return 'clinical_review'
    if 'Case Reports' in values:return 'single_case'
    return None


def body_ids(p):
    article=p['article'];boundary=len(article['title'])+1
    return [x['id'] for x in article['sentences'] if x['end']>boundary]


def cited(p,values):
    require(values and set(values)<=set(body_ids(p)),'Positive clinical claim needs actual abstract sentence IDs')


def snippet(p,values):
    by={x['id']:x for x in p['article']['sentences']}
    for i in values:
        if i in body_ids(p):return ' '.join(by[i]['text'].split()[:20])
    return ''


def or3(*values):return True if any(x is True for x in values) else None if any(x is None for x in values) else False
def and3(*values):return False if any(x is False for x in values) else None if any(x is None for x in values) else True
def state(value):return True if value=='supported' else False if value=='absent' else None


def normalize(q,p):
    exact=True if q['organism_scope']=='exact_species' else None if q['organism_scope']=='unclear' else False
    human=state(q['human_any']['status']);direct=state(q['direct_pulmonary']['status'])
    related_state=state(q['disease_related']['status']);chain=q['disease_related']['chain']
    related=and3(related_state,True if chain=='pathogen_disease_lung' else None if chain=='unclear' else False)
    mature=q['study_design'] in judge.MATURE_DESIGNS if q['study_design']!='unclear' else None
    level=q['direct_pulmonary']['level']
    direct_mature=and3(direct,True if level in ['recognized','observational'] else None if level=='unclear' else False,mature)
    complication_mature=and3(related,mature)
    qualified=and3(exact,or3(direct_mature,complication_mature))
    basis='direct' if direct_mature is True else 'disease_related' if complication_mature is True else 'uncertain' if qualified is None else 'none'
    issues=[]
    if (direct is True or related is True) and human is False:issues.append('pulmonary_claim_but_human_causation_absent')
    if direct is True and level=='none':issues.append('direct_attribution_marked_supported_but_level_none')
    if related_state is True and chain=='not_applicable':issues.append('disease_relation_supported_but_chain_not_applicable')
    if issues:qualified=None;human=None;basis='uncertain'
    # Uncertain direct mechanism does not veto an explicit disease-related chain.
    relevant=list(dict.fromkeys(q['direct_pulmonary']['sentence_ids']+q['disease_related']['agent_disease_sentence_ids']+q['disease_related']['disease_lung_sentence_ids']))
    human_ids=list(dict.fromkeys(q['human_any']['sentence_ids']+relevant))
    nonhuman=q['study_design']=='laboratory_nonhuman'
    signature=[q['organism_scope'],human,qualified,basis,nonhuman,state(q['noncausal']['status']),state(q['secondary_other']['status'])]
    unresolved=bool(issues or human is None or qualified is None or exact is None)
    return {'scope':q['organism_scope'],'human_any':human,'qualified':qualified,'basis':basis,
        'mature_design':mature,'direct_mature':direct_mature,'disease_related_complete':related,
        'nonhuman_only':nonhuman,'design_ids':q.get('study_design_sentence_ids',[]),
        'human_ids':human_ids,'qualifying_ids':relevant,'identity_ids':q['identity_sentence_ids'],
        'noncausal':state(q['noncausal']['status']),'secondary_other':state(q['secondary_other']['status']),
        'signature':signature,'status':'unresolved_qc' if unresolved else 'resolved',
        'internal_issues':issues,'open_nondecision_fields':[k for k in ['noncausal','secondary_other'] if state(q[k]['status']) is None],
        'source_summaries':{k:q[k]['summary'] for k in ASPECTS}}


def unresolved(reason):
    return {'scope':'unclear','human_any':None,'qualified':None,'basis':'uncertain','mature_design':None,
        'direct_mature':None,'disease_related_complete':None,'human_ids':[],'qualifying_ids':[],'identity_ids':[],
        'noncausal':None,'secondary_other':None,'signature':None,'status':'unresolved_qc',
        'internal_issues':[reason],'source_summaries':{}}


def consensus(left,right):
    if left['signature']!=right['signature'] or left['status']!='resolved' or right['status']!='resolved':
        out=unresolved('paired_QC_disagreement_or_unresolved_semantics')
        out['paired_signatures']=[left['signature'],right['signature']];return out
    out=copy.deepcopy(left)
    for key in ['human_ids','qualifying_ids','identity_ids']:out[key]=list(dict.fromkeys(left[key]+right[key]))
    out['status']='qc_consensus';return out

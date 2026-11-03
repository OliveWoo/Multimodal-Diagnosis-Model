"""Unique module names and immutable artifacts for an article-only replay."""
import copy
import json
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'OBER_A1_S_k_20261031'))
import a1s_common as a
import a1_articles as olda1
import run_resilient_judge as transport
c=a.c;judge=a.source.judge
RUN=HERE/'runs/isolation_v1'
OLD_A1=ROOT/'OBER_A1_S_k_20261031/runs/a1s_alias_v2'
OLD_JUDGE=a.source.PARENT/'judgements'
MODEL=c.MODEL
c.RUN=RUN;c.WORKERS=24;c.api=transport.body_resilient_api
load,save,require,digest,write_new=c.load,c.save,c.require,c.digest,c.write_new
obj,enum,ids,string=c.obj,c.enum,c.ids,c.string


def array(item):return {'type':'array','items':item}


def freeze(stage,paths):
    save(RUN/(stage+'_freeze.json'),{'files':{str(p.resolve()):digest(p.read_bytes()) for p in sorted(set(paths))},
        'historically_viewed_development_data':True,'not_independent_test':True})


def verify(stage):
    for p,h in load(RUN/(stage+'_freeze.json'))['files'].items():require(digest(Path(p).read_bytes())==h,'Frozen artifact changed: '+p)


def response_files(folder):
    return [folder/'result.json']+list(folder.glob('*_request.json'))+list(folder.glob('*_raw_response.json'))


def nonpaper_pass(q):
    d=judge.selection(q,'broad')
    return all(d[k] for k in ['clinical_gate','case_bridge_gate','revision_gate','counter_gate'])


def quoted(text,article,source_ids):
    compact=lambda x:' '.join(x.split())
    sentences={x['id']:x for x in article['sentences']}
    require(source_ids and set(source_ids)<=set(sentences),'Unknown/missing article source IDs')
    require(judge.has_abstract_citation(article,source_ids),'Title-only claim cannot establish clinical role')
    fragment=compact(text)
    require(fragment and len(fragment.split())<=20,'Article excerpt must have 1..20 words')
    require(fragment.casefold() in compact(article['title']+' '+article['abstract']).casefold(),'Excerpt absent from article')
    require(fragment.casefold() in compact(' '.join(sentences[i]['text'] for i in source_ids)).casefold(),'Excerpt absent from cited sentences')


def projections(q):
    attributed=[x for x in q['human_attributions'] if x['organism_scope']=='exact_species']
    pulmonary=[x for x in attributed if x['lower_respiratory']]
    other=[x for x in q['other_contexts'] if x['organism_scope']=='exact_species']
    pulmonary_other=[x for x in other if x['lower_respiratory']]
    order={'recognized':0,'observational':1,'case_only':2}
    attributed=sorted(attributed,key=lambda x:(not x['lower_respiratory'],order[x['level']],x['claim_id']))
    pulmonary=sorted(pulmonary,key=lambda x:(order[x['level']],x['claim_id']))
    context='unclear';chosen=None
    if attributed:context='causal_infection';chosen=attributed[0]
    else:
        for kind,target in [('colonization','colonization_or_contamination'),('reactivation','reactivation_or_bystander'),
                            ('detection','detection_only'),('nonhuman','nonhuman')]:
            matching=[x for x in other if x['kind']==kind]
            if matching:context=target;chosen=matching[0];break
    identity={'exact_species':'aligned_resolution','other':'other','genus_only':'uncertain','unclear':'uncertain'}[q['organism_scope']]
    if pulmonary:
        claim={'recognized':'recognized_etiology','observational':'observational_etiology','case_only':'causal_case_only'}[pulmonary[0]['level']]
        stance='mixed' if any(x['kind'] in ['colonization','reactivation','uncertainty'] for x in pulmonary_other) else 'support'
        causal= pulmonary[0]
    else:
        claim='detection_only' if any(x['kind'] in ['detection','colonization','reactivation'] for x in pulmonary_other) else 'no_causal_claim'
        stance='counter' if any(x['kind'] in ['colonization','reactivation'] for x in pulmonary_other) else 'neutral'
        causal=None
    evidence=list(dict.fromkeys(q['scope_sentence_ids']+[sid for x in pulmonary+pulmonary_other for sid in x['sentence_ids']]))
    support=next(iter(pulmonary or attributed),None)
    counter=next((x for x in pulmonary_other if x['kind'] in ['colonization','reactivation']),None)
    return {'G_T':{'identity_match':identity,'study_design':q['study_design'],'clinical_claim':claim,'stance':stance,
                   'sentence_ids':evidence},
        'A1':{'organism_scope':q['organism_scope'],'human_clinical_context':context,
              'scope_sentence_ids':q['scope_sentence_ids'],'context_sentence_ids':chosen['sentence_ids'] if chosen else [],
              'support_span':support['excerpt'] if support else '', 'counter_span':counter['excerpt'] if counter else ''},
        'mixed_contexts_retained':bool(pulmonary and pulmonary_other)}


def substitute_judge(original,ledger_map):
    result=copy.deepcopy(original)
    for pmid,assess in result['article_assessments'].items():
        if pmid not in ledger_map:continue
        projected=projections(ledger_map[pmid])['G_T']
        for key in ['identity_match','study_design','clinical_claim','stance']:assess[key]=projected[key]
        assess['sentence_ids']=list(dict.fromkeys(assess['sentence_ids']+projected['sentence_ids']))
    for key in original:
        if key!='article_assessments':require(result[key]==original[key],'Changed frozen clinical/revision judgement')
    for pmid in original['article_assessments']:
        for key in ['applicability','question_ids','reason']:
            require(result['article_assessments'][pmid][key]==original['article_assessments'][pmid][key],'Changed frozen article applicability')
    return result


def substitute_a1(original,q):
    result=copy.deepcopy(original);p=projections(q)['A1']
    for key in ['organism_scope','human_clinical_context']:result[key]=p[key]
    for field,key in [('organism_scope','scope_sentence_ids'),('human_clinical_context','context_sentence_ids')]:
        result['field_evidence'][field]={'article_sentence_ids':p[key],'patient_fact_ids':[]}
    result['support_span']=p['support_span'];result['counter_span']=p['counter_span']
    for key in a.legacy.FIELDS[2:]:
        require(result[key]==original[key] and result['field_evidence'][key]==original['field_evidence'][key],'Changed six fixed case alignment fields')
    return result

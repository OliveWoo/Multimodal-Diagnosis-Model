"""One existing paper per request; frozen legacy A1 truth table and k sweep."""
import argparse
import collections
import copy
import json
import re
from pathlib import Path
import a1s_common as a
from run_resilient_judge import body_resilient_api

PROMPT=a.legacy.PROMPT
IDENTIFIERS=re.compile(r'開單醫師|報告醫師|判讀醫師|病歷號碼?|身分證|患者姓名|病患姓名|病人姓名|出生日期|\bDOB\b',re.I)


def make_view(packet,index_type):
    facts=copy.deepcopy([f for f in packet['patient_facts'] if f['kind'] in ['index','background','imaging']])
    for f in facts:
        a.require(not IDENTIFIERS.search(f['text']),'Identifier field in A1 input; transmission blocked at fact '+f['id'])
    facts += [{'id':'P_SITE','kind':'protocol_target','temporal_role':'protocol','value':'lower_respiratory',
               'caveat':'Study comparison target, not confirmation of the patient disease.'},
              {'id':'P_SYNDROME','kind':'protocol_target','temporal_role':'protocol','value':'lower_respiratory_infection',
               'caveat':'Study comparison target, not patient-level etiologic adjudication.'}]
    return {'candidate_token':packet['candidate_token'],'clinical_site':'lower_respiratory',
        'target_syndrome':'lower_respiratory_infection','index_specimen_type':index_type,
        'time_semantics':'+/-2 calendar dates for primary collection. Earlier is history, timing_unknown cannot establish contemporaneous phenotype or onset.',
        'absence_semantics':'Unprovided, unperformed and unknown are not negative tests. Do not infer unrecorded host or treatment details.',
        'facts':facts,'withheld':['mNGS reads/RPM/rank','culture/PCR/FilmArray results',
            'upstream selection or reasons','G/T reflections and S pass','reference answers']}


def schema(packet):
    s=copy.deepcopy(a.legacy.SCHEMA)
    s['properties']['candidate_token']=a.c.enum([packet['case_card']['candidate_token']])
    s['properties']['pmid']=a.c.enum([packet['article']['pmid']])
    fields=s['properties']['field_evidence']['properties']
    for key in a.legacy.FIELDS:
        fields[key]['properties']['article_sentence_ids']=a.c.ids([x['id'] for x in packet['article']['sentences']])
        fields[key]['properties']['patient_fact_ids']=a.c.ids([x['id'] for x in packet['case_card']['facts']])
    return s


def validate(q,p):
    a.require(set(q)==set(a.legacy.SCHEMA['properties']),'Extra/missing A1 fields')
    a.require(q['candidate_token']==p['case_card']['candidate_token'] and q['pmid']==p['article']['pmid'],'Wrong candidate/paper')
    a.legacy.derive_article_verdict({k:q[k] for k in a.legacy.FIELDS})
    checked=a.legacy.validate_extraction(q,p['case_card'],p['article'])
    # Source-ID/format errors may be corrected once; a genuine unknown or an
    # unavailable support excerpt is retained, never retried to get more support.
    a.require(not checked['field_trace_failures'],'A1 field source errors: '+','.join(checked['field_trace_failures']))


def prepare():
    a.source.verify('prediction');a.source.verify('source_guard')
    a.source.parent_verify('judge_inputs');a.source.parent_verify('judge_outputs')
    a.require(not (a.RUN/'plan.json').exists(),'A1/S already prepared')
    source=a.load(a.source.RUN/'predictions_without_answers.json')
    guard={r['token']:r for r in a.load(a.source.RUN/'source_guard_predictions_without_answers.json')}
    rows=[];jobs=[];files=[];stats=collections.Counter()
    original_base=a.ROOT/'OBER_case_lookup_20261031/runs/lookup_v1/cases'
    for r in source:
        active=guard[r['token']]['decisions']['M_C0_Q0_S1_NONE']
        row={k:r[k] for k in ['case_id','group','token','organism','native_flags']}|{
            'S_pass':active,'S_raw_pass':r['features']['T_broad'] and r['features']['G_broad'],
            'source_deferred':guard[r['token']]['deferred_due_to_known_later_input'],
            'folder':f"candidates/{r['case_id']}/{r['token']}"}
        if not active:rows.append(row);continue
        packets={arm:a.load(a.source.PARENT/'judgements'/r['case_id']/r['token']/arm/'input.json') for arm in ['G','T']}
        a.require(packets['G']['patient_facts']==packets['T']['patient_facts'],'G/T clinical facts differ')
        base=a.load(original_base/r['case_id']/'full_source_facts.json')
        view=make_view(packets['T'],base['index_specimen_type'])
        folder=a.RUN/row['folder'];a.save(folder/'a1_patient_view.json',view);files.append(folder/'a1_patient_view.json')
        byid={};arms={}
        for arm in ['T','G']:
            arms[arm]=[x['pmid'] for x in packets[arm]['papers']]
            a.require(len(arms[arm])<=10 and len(set(arms[arm]))==len(arms[arm]),'Duplicate/too many PMIDs')
            for paper in packets[arm]['papers']:
                item={k:paper[k] for k in ['pmid','title','abstract','sentences']}
                if item['pmid'] in byid:
                    a.require(all(item[k]==byid[item['pmid']][k] for k in ['title','abstract']),
                        'Same PMID has inconsistent article text; cannot share A1')
                else:byid[item['pmid']]=item
        row['arm_pmids']=arms;row['unique_pmid_count']=len(byid)
        row['repeat_all_articles']=a.source.native.passes(r['native_flags'],'both')
        a.save(folder/'arm_pmids.json',arms);a.save(folder/'papers.json',list(byid.values()))
        files += [folder/'arm_pmids.json',folder/'papers.json']
        for pmid,paper in byid.items():
            packet={'organism':r['organism'],'case_card':view,'article':paper}
            for repeat in range(2 if row['repeat_all_articles'] else 1):
                dest=folder/'articles'/pmid/f'r{repeat}'
                a.save(dest/'input.json',packet);a.save(dest/'schema.json',schema(packet))
                files += [dest/'input.json',dest/'schema.json']
                jobs.append({'case_id':r['case_id'],'token':r['token'],'pmid':pmid,'replicate':repeat,
                    'folder':dest.relative_to(a.RUN).as_posix(),'has_abstract':bool(paper['abstract'])})
                stats['unique_main_articles' if repeat==0 else 'repeat_articles']+=1
                stats['main_API_calls' if repeat==0 else 'repeat_API_calls']+=bool(paper['abstract'])
        rows.append(row)
    a.save(a.RUN/'plan.json',{'candidates':rows,'jobs':jobs,'model':a.MODEL,'reasoning':'medium',
        'primary':a.PRIMARY,'rules':a.configs(),'scope_cases':49,'scope_context':319,
        'conditional_A1_on_S_only':True,'A1_standalone_full_pool_not_available':True,
        'primary_replicate':0,'repeat_selection':'All unique papers of every S-pass native-BOTH candidate; fixed once, not by labels.',
        'no_gold_loaded_during_prepare':True,'historical_labels_seen':True,'no_new_case_lookup':True,
        'legacy_A1_truth_table_unchanged':True,'stats':dict(stats),
        'parent_module_prediction_sha256':a.digest((a.source.RUN/'prediction_freeze.json').read_bytes()),
        'parent_guard_sha256':a.digest((a.source.RUN/'source_guard_freeze.json').read_bytes())})
    a.write_new(a.RUN/'prompt_snapshot.txt',PROMPT.encode())
    a.freeze('inputs',files+[a.RUN/'plan.json',a.RUN/'prompt_snapshot.txt',Path(__file__),a.HERE/'a1s_common.py',
        a.HERE/'PROTOCOL_zh.md',Path(a.legacy.__file__),a.legacy.SOURCE_PROMPT,
        a.ROOT/'RAG_re_casefit/rag_re_casefit/rules.py',a.ROOT/'OBER_reason_reflection_20261031/reason_common.py',
        a.ROOT/'OBER_reason_reflection_20261031/run_resilient_judge.py'])
    print({'S_candidates':sum(r['S_pass'] for r in rows),'all_context':len(rows),'settings':len(a.configs()),**dict(stats)},flush=True)


def run():
    a.verify('inputs')
    if (a.RUN/'outputs_freeze.json').exists():a.verify('outputs');print('A1/S complete; no new calls.');return
    a.c.RUN=a.RUN;a.c.WORKERS=24;a.c.api=body_resilient_api
    plan=a.load(a.RUN/'plan.json')
    def worker(job,env):
        folder=a.RUN/job['folder'];p=a.load(folder/'input.json')
        if not job['has_abstract']:
            a.save(folder/'unavailable.json',a.legacy.no_abstract(p['case_card'],p['article']));return
        a.c.infer(folder,PROMPT,p,a.load(folder/'schema.json'),validate,env,3500)
        print({'A1':job['case_id'],'pmid':job['pmid'],'repeat':job['replicate']},flush=True)
    a.c.dispatch(plan['jobs'],worker,'a1s')
    files=[]
    for job in plan['jobs']:
        folder=a.RUN/job['folder'];p=a.load(folder/'input.json')
        if job['has_abstract']:
            q=a.c.validate_saved(folder,validate);checked=a.legacy.validate_extraction(q,p['case_card'],p['article'])
            files += [folder/'result.json']+list(folder.glob('*_request.json'))+list(folder.glob('*_raw_response.json'))
        else:checked=a.load(folder/'unavailable.json');files.append(folder/'unavailable.json')
        a.save(folder/'grading.json',checked);files.append(folder/'grading.json')
    features=[]
    for row in plan['candidates']:
        counts=None;repeat_counts=None
        if row['S_pass']:
            folder=a.RUN/row['folder'];counts={};repeat_counts={}
            for arm in ['G','T']:
                judgments=[dict(a.load(folder/'articles'/pmid/'r0/grading.json'),pmid=pmid) for pmid in row['arm_pmids'][arm]]
                counts[arm]=a.legacy.aggregate_judgments(judgments)
                if row['repeat_all_articles']:
                    repeated=[dict(a.load(folder/'articles'/pmid/'r1/grading.json'),pmid=pmid) for pmid in row['arm_pmids'][arm]]
                    repeat_counts[arm]=a.legacy.aggregate_judgments(repeated)
            a.save(folder/'A1_counts.json',{'main':counts,'repeat':repeat_counts or None});files.append(folder/'A1_counts.json')
        features.append(row|{'counts':counts,'repeat_counts':repeat_counts or None,
            'decisions':{s['id']:a.admit(row,counts,s) for s in plan['rules']},
            'repeat_decisions':{s['id']:a.admit(row,repeat_counts,s) for s in plan['rules'] if s['signal']=='both'} if repeat_counts else None})
    a.save(a.RUN/'predictions_without_answers.json',features);files.append(a.RUN/'predictions_without_answers.json')
    a.freeze('outputs',files)
    print({'all_k_choices_frozen':True,'A1_candidates':sum(r['S_pass'] for r in features),'reference_read':False},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','run'])
    {'prepare':prepare,'run':run}[p.parse_args().command]()

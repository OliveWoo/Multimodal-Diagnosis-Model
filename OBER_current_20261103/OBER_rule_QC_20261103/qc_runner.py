"""Frozen rules -> two QC calls only where routed -> frozen legacy case replay."""
import argparse
import copy
import json
from pathlib import Path
import qc_common as s
import qc_rules as rules
import qc_review as review

STRATEGIES=['RULE_ONLY','QC_R0','QC_R1','QC_CONSENSUS','QC_DIRECT_ONLY']


def prepare():
    if (s.RUN/'inputs_freeze.json').exists():s.verify('inputs');print('Inputs already frozen.');return
    for stage in ['inputs','ledger_outputs','source_snippets','predictions']:
        for path,h in s.load(s.PARENT/(stage+'_freeze.json'))['files'].items():
            s.require(s.digest(Path(path).read_bytes())==h,'Changed parent artifact '+path)
    source=s.load(s.PARENT/'plan.json');rows=copy.deepcopy(source['potential'])
    metadata={};files=[]
    for row in rows:
        for arm in ['G','T']:
            f=s.old.OLD_JUDGE/row['case_id']/row['token']/arm;p=s.load(f/'input.json')
            s.c.validate_saved(f,s.judge.validate)
            files += [f/'input.json',f/'schema.json']+s.response_files(f)
            for a in p['papers']:
                key=(row['organism'],a['pmid']);types=a.get('publication_types',[])
                if key in metadata:s.require(set(types)==set(metadata[key]),'Publication metadata differs across lanes')
                metadata[key]=types
        f=s.RUN/'candidates'/row['case_id']/row['token']
        for name in ['a1_patient_view.json','papers.json']:
            s.save(f/name,s.load(s.PARENT/row['folder']/name));files.append(f/name)
        row['folder']=f.relative_to(s.RUN).as_posix()
    articles=[];jobs=[]
    for j in source['ledger_jobs']:
        if j['replicate']!=0:continue
        parent=s.PARENT/j['folder'];p=s.load(parent/'input.json')
        p['publication_types']=metadata[p['organism'],p['article']['pmid']]
        f=s.RUN/'articles'/p['article_key'];route=rules.route(p)
        s.save(f/'source.json',p);s.save(f/'routing.json',route);files += [f/'source.json',f/'routing.json']
        articles.append({'key':p['article_key'],'pmid':p['article']['pmid'],'organism':p['organism'],
            'folder':f.relative_to(s.RUN).as_posix(),'route':route['route'],'reasons':route['reasons']})
        if route['route']=='QC':
            for rep in range(2):
                dest=f/f'r{rep}';s.save(dest/'input.json',p);s.save(dest/'schema.json',review.schema(p))
                jobs.append({'token':p['article_key'],'replicate':rep,'folder':dest.relative_to(s.RUN).as_posix()})
                files += [dest/'input.json',dest/'schema.json']
    s.require(len(articles)==58 and len(rows)==5,'Changed frozen pilot scope')
    plan={'tracked':source['tracked'],'potential':rows,'articles':articles,'QC_jobs':jobs,'strategies':STRATEGIES,
        'primary':'QC_CONSENSUS_F0_match_k2','QC_is_evidence_review_not_sequencing_quality':True,
        'independent_test':False,'historical_labels_seen':True,'conditional_replay':True,
        'model':s.MODEL,'reasoning':'medium','no_new_retrieval_or_patient_facts':True}
    s.save(s.RUN/'plan.json',plan)
    deps=list(s.HERE.glob('*.py'))+[s.HERE/'PROTOCOL_zh.md',s.HERE/'README.md']
    deps += [Path(m.__file__) for m in [s.old,s.old.olda1,s.old.judge,s.c,s.a.legacy]]
    deps += [Path(s.a.legacy.SOURCE_PROMPT),Path(s.old.transport.__file__)]
    deps += [s.ROOT/'RAG_re_casefit/rag_re_casefit/rules.py',s.PARENT/'plan.json',s.PARENT/'alignment_jobs.json',s.PARENT/'predictions_without_answers.json']
    s.freeze('inputs',files+deps+[s.RUN/'plan.json'])
    print({'articles':58,'rule_accepted':sum(x['route']=='rule' for x in articles),
        'QC_articles':sum(x['route']=='QC' for x in articles),'planned_public_QC_calls':len(jobs),'answers_not_read':True},flush=True)


def run_QC(smoke=False):
    s.verify('inputs')
    if (s.RUN/'QC_outputs_freeze.json').exists():s.verify('QC_outputs');print('QC outputs cached.');return
    jobs=s.load(s.RUN/'plan.json')['QC_jobs'];chosen=jobs[:2] if smoke else jobs
    def worker(j,env):
        f=s.RUN/j['folder'];s.c.infer(f,review.PROMPT,s.load(f/'input.json'),s.load(f/'schema.json'),review.validate,env,4500)
        s.c.validate_saved(f,review.validate);print({'QC':j['token'],'replicate':j['replicate']},flush=True)
    if chosen:s.c.dispatch(chosen,worker,'QC')
    if smoke:print('Two preplanned source-ID QC outputs validated; no scoring.');return
    files=[]
    for j in jobs:s.c.validate_saved(s.RUN/j['folder'],review.validate);files+=s.response_files(s.RUN/j['folder'])
    s.freeze('QC_outputs',files)


def derive_articles():
    s.verify('inputs');s.verify('QC_outputs')
    if (s.RUN/'article_decisions_freeze.json').exists():s.verify('article_decisions');return
    records={};paths=[]
    for a in s.load(s.RUN/'plan.json')['articles']:
        f=s.RUN/a['folder'];p=s.load(f/'source.json');route=s.load(f/'routing.json')
        if a['route']=='rule':values={k:copy.deepcopy(route['rule_record']) for k in STRATEGIES};norms=None
        elif a['route']=='unavailable':values={k:s.unresolved('unavailable_abstract') for k in STRATEGIES};norms=None
        else:
            norms=[s.normalize(s.c.validate_saved(f/f'r{rep}',review.validate),p) for rep in range(2)]
            values={'RULE_ONLY':s.unresolved('QC_not_reviewed_in_rule_only_control'),
                'QC_R0':norms[0],'QC_R1':norms[1],'QC_CONSENSUS':s.consensus(*norms)}
            values['QC_DIRECT_ONLY']=copy.deepcopy(values['QC_CONSENSUS'])
        # Direct-only is an explicit diagnostic comparator, not a relabeling of
        # disease-related evidence as absent biological capability.
        if values['QC_DIRECT_ONLY']['basis']=='disease_related':
            values['QC_DIRECT_ONLY']['qualified']=False;values['QC_DIRECT_ONLY']['basis']='withheld_disease_related_by_comparator'
        out={**a,'normalized_pair':norms,'strategies':values}
        s.save(f/'article_decisions.json',out);paths.append(f/'article_decisions.json');records[a['key']]=out
    s.save(s.RUN/'article_decisions_without_answers.json',records);paths.append(s.RUN/'article_decisions_without_answers.json')
    s.freeze('article_decisions',paths)


def sources(row):
    return {a:s.c.validate_saved(s.old.OLD_JUDGE/row['case_id']/row['token']/a,s.judge.validate) for a in ['G','T']}


def article_gate(row,strategy):
    records=s.load(s.RUN/'article_decisions_without_answers.json');qs=sources(row);result={}
    for arm,q in qs.items():
        old=s.judge.selection(q,'broad');cited={x['pmid'] for x in q['revision']['paper_refs']}
        qualified=[];pending=[]
        for pmid,assessment in q['article_assessments'].items():
            if pmid not in row['article_keys'] or pmid not in cited or assessment['applicability'] not in ['aligned','transportable']:continue
            item=records[row['article_keys'][pmid]]['strategies'][strategy]
            if item['status']=='unresolved_qc' or item['qualified'] is None:pending.append(pmid)
            elif item['qualified'] is True and item['qualifying_ids']:qualified.append({'pmid':pmid,'basis':item['basis']})
        fixed={k:old[k] for k in ['clinical_gate','case_bridge_gate','revision_gate','counter_gate']}
        result[arm]={**fixed,'article_support_gate':bool(qualified),'pass':all(fixed.values()) and bool(qualified),
            'qualified_articles':qualified,'pending_QC_articles':pending}
    return result


def prepare_A1():
    s.verify('article_decisions')
    if (s.RUN/'A1_inputs_freeze.json').exists():s.verify('A1_inputs');return
    parent_jobs=s.load(s.PARENT/'alignment_jobs.json');bytoken={}
    for j in parent_jobs:bytoken.setdefault(j['token'],[]).append(j)
    jobs=[];files=[];casegates=[]
    for row in s.load(s.RUN/'plan.json')['potential']:
        gates={strategy:article_gate(row,strategy) for strategy in STRATEGIES}
        needed=row['old_S'] or any(all(x['pass'] for x in z.values()) for z in gates.values())
        casegates.append({'token':row['token'],'needed':needed,'gates':gates})
        if not needed:continue
        if row['token'] in bytoken:
            for j in bytoken[row['token']]:
                job=copy.deepcopy(j);job['reused']=True;jobs.append(job);f=Path(j['folder'])
                files += [f/'input.json',f/'schema.json']+ (s.response_files(f) if j['has_abstract'] else [f/'unavailable.json'])
            continue
        f=s.RUN/row['folder'];view=s.load(f/'a1_patient_view.json')
        for paper in s.load(f/'papers.json'):
            for rep in range(2):
                dest=f/'legacy_A1'/paper['pmid']/f'r{rep}';p={'organism':row['organism'],'case_card':view,'article':paper}
                aliases=s.a.legacy.CLINICAL_NAMES.get(row['organism'],[])
                if aliases:p['documented_clinical_names_not_scoring_aliases']=aliases
                s.save(dest/'input.json',p);s.save(dest/'schema.json',s.olda1.schema(p))
                jobs.append({k:row[k] for k in ['token','case_id','organism']}|{'pmid':paper['pmid'],'replicate':rep,
                    'folder':str(dest.resolve()),'reused':False,'has_abstract':bool(paper['abstract'])})
                files += [dest/'input.json',dest/'schema.json']
    s.save(s.RUN/'A1_jobs.json',jobs);s.save(s.RUN/'case_gates_without_answers.json',casegates)
    s.freeze('A1_inputs',files+[s.RUN/'A1_jobs.json',s.RUN/'case_gates_without_answers.json'])
    print({'new_legacy_A1_calls':sum(not x['reused'] and x['has_abstract'] for x in jobs),
        'reused_A1_outputs':sum(x['reused'] and x['has_abstract'] for x in jobs)},flush=True)


def run_A1():
    s.verify('A1_inputs')
    if (s.RUN/'A1_outputs_freeze.json').exists():s.verify('A1_outputs');return
    jobs=s.load(s.RUN/'A1_jobs.json');new=[j for j in jobs if not j['reused']]
    def worker(j,env):
        f=Path(j['folder']);p=s.load(f/'input.json')
        if not j['has_abstract']:s.save(f/'unavailable.json',s.a.legacy.no_abstract(p['case_card'],p['article']));return
        s.c.infer(f,s.olda1.PROMPT,p,s.load(f/'schema.json'),s.olda1.validate,env,3500)
        print({'A1':j['case_id'],'pmid':j['pmid'],'replicate':j['replicate']},flush=True)
    if new:s.c.dispatch(new,worker,'A1')
    files=[]
    for j in jobs:
        f=Path(j['folder'])
        if j['has_abstract']:s.c.validate_saved(f,s.olda1.validate);files+=s.response_files(f)
        else:files.append(f/'unavailable.json')
    s.freeze('A1_outputs',files)


def a1_records(row,strategy,rep):
    jobs=[j for j in s.load(s.RUN/'A1_jobs.json') if j['token']==row['token'] and j['replicate']==rep]
    if not jobs:return None
    articles=s.load(s.RUN/'article_decisions_without_answers.json');records={}
    for j in jobs:
        f=Path(j['folder']);p=s.load(f/'input.json')
        if not j['has_abstract']:record=s.a.legacy.no_abstract(p['case_card'],p['article']);pending=False
        else:
            original=s.c.validate_saved(f,s.olda1.validate);q=copy.deepcopy(original)
            item=articles[row['article_keys'][j['pmid']]]['strategies'][strategy]
            pending=item['status']=='unresolved_qc'
            scope='unclear' if pending else item['scope'];human=None if pending else item['human_any']
            nonhuman=item.get('nonhuman_only',False) and not pending
            q['organism_scope']=scope;q['human_clinical_context']='nonhuman' if nonhuman else 'causal_infection' if human is True else 'unclear'
            q['field_evidence']['organism_scope']={'article_sentence_ids':[] if scope=='unclear' else item['identity_ids'],'patient_fact_ids':[]}
            q['field_evidence']['human_clinical_context']={'article_sentence_ids':item.get('design_ids',[]) if nonhuman else item['human_ids'] if human is True else [],'patient_fact_ids':[]}
            packet=s.load(s.RUN/'articles'/row['article_keys'][j['pmid']]/'source.json')
            q['support_span']=s.snippet(packet,item['qualifying_ids'] or item['human_ids']) if human is True else ''
            q['counter_span']=''
            for key in s.a.legacy.FIELDS[2:]:
                s.require(q[key]==original[key] and q['field_evidence'][key]==original['field_evidence'][key],'Changed frozen A1 case fields')
            record=s.a.legacy.validate_extraction(q,p['case_card'],p['article'])
            s.require(not record['field_trace_failures'],'Invalid source substitution in derived A1')
        records[j['pmid']]=dict(record,pmid=j['pmid'],evidence_review_qc_pending=pending,derived_not_raw_LLM_response=True)
    return records


def configs():
    return [{'id':f'{strategy}_F{rep}_{count}_k{k}','strategy':strategy,'fit':rep,'count':count,'k':k}
        for strategy in STRATEGIES for rep in range(2) for count in ['match','strong'] for k in range(11)]


def predict():
    for name in ['inputs','QC_outputs','article_decisions','A1_inputs','A1_outputs']:s.verify(name)
    if (s.RUN/'predictions_freeze.json').exists():s.verify('predictions');return
    plan=s.load(s.RUN/'plan.json');active={r['token']:r for r in plan['potential']};rulespec=configs();features=[];files=[]
    for row in plan['tracked']:
        out={k:row[k] for k in ['case_id','group','token','organism','status']};out['decisions']={r['id']:False for r in rulespec}
        if row['token'] in active:
            r=active[row['token']];variants={};out['folder']=r['folder']
            for strategy in STRATEGIES:
                gates=article_gate(r,strategy);passed=all(g['pass'] for g in gates.values())
                for rep in range(2):
                    prefix=f'{strategy}_F{rep}';records=a1_records(r,strategy,rep)
                    counts={a:s.a.legacy.aggregate_judgments([records[pmid] for pmid in r['arm_pmids'][a]]) for a in ['G','T']} if records else None
                    s.require(not passed or counts is not None,'Passing source gate needs actual A1 extraction')
                    pending={a:[pmid for pmid in r['arm_pmids'][a] if records and records[pmid]['evidence_review_qc_pending']] for a in ['G','T']}
                    variants[prefix]={'S_pass':passed,'arm_gates':gates,'counts':counts,'pending_A1_QC_articles':pending,
                        'needs_source_review':any(g['pending_QC_articles'] for g in gates.values()) or any(pending.values())}
                    if records:
                        dest=s.RUN/r['folder']/'replay'/prefix/'A1_derived_records.json';s.save(dest,records);files.append(dest)
                    for rule in rulespec:
                        if rule['strategy']!=strategy or rule['fit']!=rep:continue
                        count=0 if counts is None else min(counts[a]['strong_match_count']+(counts[a]['partial_match_count'] if rule['count']=='match' else 0) for a in ['G','T'])
                        out['decisions'][rule['id']]=bool(passed and (rule['k']==0 or count>=rule['k']))
            out['variants']=variants
        features.append(out)
    s.save(s.RUN/'predictions_without_answers.json',features);s.save(s.RUN/'rules.json',rulespec)
    s.freeze('predictions',files+[s.RUN/'predictions_without_answers.json',s.RUN/'rules.json'])
    print({'primary_additions':sum(r['decisions'][plan['primary']] for r in features),'all_source_QC_and_predictions_frozen':True,'answers_not_read':True},flush=True)


def pipeline():
    prepare();run_QC();derive_articles();prepare_A1();run_A1();predict()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','smoke','pipeline','predict'])
    {'prepare':prepare,'smoke':lambda:run_QC(True),'pipeline':pipeline,'predict':predict}[p.parse_args().stage]()

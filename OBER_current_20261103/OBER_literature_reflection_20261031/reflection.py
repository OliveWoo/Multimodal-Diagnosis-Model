"""Paired source-only / literature reflection. No organism-specific selectors."""
import argparse
import collections
import concurrent.futures
import copy
import csv
import io
import json
import random
import time
from pathlib import Path
import data_and_pubmed as d

HERE=d.HERE;ROOT=d.ROOT;RUN=d.RUN
load,save,require,digest,write_new=d.load,d.save,d.require,d.digest,d.write_new
MODEL=d.clinical.MODEL
PROMPT=(HERE/'reflection_prompt.txt').read_text(encoding='utf-8')
POLICIES=['strict','balanced','broad']


def obj(fields):return d.clinical.obj(fields)
def enum(values):return {'type':'string','enum':values}
def strings():return {'type':'array','items':{'type':'string'}}
def string():return {'type':'string'}


REF=obj({'pmid':string(),'sentence_ids':strings()})
CHECK=obj({'status':enum(['supported','plausible','unknown','contradicted']),
    'patient_fact_ids':strings(),'reason':string()})
ITEM=obj({'candidate_token':string(),
    'detection':CHECK,'episode_site':CHECK,'pulmonary_event':CHECK,
    'premise_audits':{'type':'array','items':obj({'premise_id':string(),
        'verdict':enum(['justified','reason_to_reconsider','not_established','unclear']),
        'patient_fact_ids':strings(),'paper_refs':{'type':'array','items':REF},'reason':string()})},
    'article_assessments':{'type':'array','items':obj({'pmid':string(),
        'organism_scope':enum(['exact_species','group_only','other','unclear']),
        'human_context':enum(['causal_infection','colonization_contamination','detection_reactivation','nonhuman','unclear']),
        'role':enum(['support','counter','neutral','irrelevant','insufficient']),
        'patient_applicability':enum(['case_relevant','generic_only','mismatch','unknown']),
        'patient_fact_ids':strings(),'sentence_ids':strings(),'reason':string()})},
    'case_attribution':obj({'level':enum(['supported','plausible','generic_only','conflicting','unknown']),
        'patient_fact_ids':strings(),'paper_refs':{'type':'array','items':REF},
        'discriminating_link':string(),'limitations':string()}),
    'reflection':obj({'verdict':enum(['case_specific_revision','generic_possibility','upstream_retained','uncertain','conflicted']),
        'premise_ids':strings(),'patient_fact_ids':strings(),'paper_refs':{'type':'array','items':REF},
        'what_changed':string(),'remaining_uncertainty':string(),'strongest_counterargument':string()})})
SCHEMA=obj({'results':{'type':'array','items':ITEM}})


def literature(name):return load(RUN/'literature'/digest(name.encode())[:16]/'pubmed.json')


def prepare():
    require(not (RUN/'plan.json').exists(),'Prepared run exists; do not overwrite')
    source=load(RUN/'source_plan.json');require((RUN/'retrieval_summary.json').exists(),'Complete PubMed retrieval first')
    batches=[];ordered=sorted(source['candidates'],key=lambda x:(x['case_id'],x['token']))
    grouped=collections.defaultdict(list)
    for row in ordered:grouped[row['case_id']].append(row)
    for cid,rows in grouped.items():
        for start in range(0,len(rows),3):
            block=rows[start:start+3];shared={};items=[]
            for row in block:
                packet=load(RUN/row['folder']/'packet.json')
                lit=literature(row['organism'])
                for f in packet['facts']:
                    if f['kind']!='mngs':
                        require(f['id'] not in shared or shared[f['id']]==f,'Conflicting same-case source ID')
                        shared[f['id']]=f
                items.append({'candidate_token':row['token'],'organism':row['organism'],
                    'upstream_premises':packet['upstream_premises'],
                    'candidate_measurements':[f for f in packet['facts'] if f['kind']=='mngs'],
                    'source_limits':packet['source_limits'],'papers':lit['articles']})
            ident=f'{cid}_{start//3+1:02}'
            for arm in ['N','L']:
                input_items=copy.deepcopy(items)
                if arm=='N':
                    for item in input_items:item['papers']=[]
                else:input_items=[item for item in input_items if any(a['abstract'] for a in item['papers'])]
                payload={'shared_patient_facts':list(shared.values()),'candidates':input_items,
                    'comparison_contract':'Each candidate independently; no ranking. Missing is unknown, not negative.'}
                require(len(d.clinical.sources.old.encoded(payload))<800000,'Batch too large; no silent truncation')
                folder=RUN/'batches'/ident/arm
                save(folder/'input.json',payload)
                batches.append({'batch_id':ident,'arm':arm,'case_id':cid,'folder':folder.relative_to(RUN).as_posix(),
                    'tokens':[i['candidate_token'] for i in input_items], 'skip':not bool(input_items)})
    deps=[Path(__file__),HERE/'data_and_pubmed.py',HERE/'reflection_prompt.txt',HERE/'PROTOCOL_zh.md',
        ROOT/'OBER_clinical_framework_20261031/experiment.py',ROOT/'OBER_clinical_framework_20261031/clinical_rules.py',
        ROOT/'OBER_v7_review_20261030/api_probe.py',ROOT/'OBER_v7_review_20261030/inspect_release.py']
    usable={row['token']:bool(literature(row['organism'])['abstract_count']) for row in ordered}
    save(RUN/'plan.json',{'schema':'ober.literature_reflection.v1','objective':'micro F1','primary':'L_balanced',
        'policies':POLICIES,'model':MODEL,'reasoning':'medium','max_workers':6,'batch_candidate_limit':3,
        'technical_repairs_per_batch':1,'independent_test':False,'historical_labels_previously_seen':True,
        'reference_read_during_prepare':False,'organism_specific_selection_rules':False,
        'source_only_and_literature_separate_requests':True,'candidates':ordered,'batches':batches,
        'usable_literature':usable,'paired_comparison_scope':'Same candidates with at least one retrieved abstract; N_all also reported.',
        'dependencies':{p.relative_to(ROOT).as_posix():digest(p.read_bytes()) for p in deps},
        'prompt_sha256':digest(PROMPT.encode()),'schema_sha256':digest(d.clinical.sources.old.encoded(SCHEMA))})
    write_new(RUN/'program_snapshot.py',Path(__file__).read_bytes())
    save(RUN/'input_freeze.json',{'reference_read':False,'files':{
        p.relative_to(RUN).as_posix():digest(p.read_bytes()) for p in RUN.rglob('*') if p.is_file()}})
    verify_inputs()
    print(json.dumps({'candidates':len(ordered),'batches':len(batches),'planned_model_calls':sum(not b['skip'] for b in batches),
        'usable_candidate_literature':sum(usable.values()),'paired_scoring_settings':6,'additional_N_all_settings':3}),flush=True)


def verify_inputs():
    plan=load(RUN/'plan.json')
    for n,h in load(RUN/'input_freeze.json')['files'].items():require(digest((RUN/n).read_bytes())==h,'Frozen input changed: '+n)
    for n,h in plan['dependencies'].items():require(digest((ROOT/n).read_bytes())==h,'Dependency changed: '+n)
    require(plan['source_only_and_literature_separate_requests'],'Invalid paired design')
    return plan


def validate(output,payload):
    require(set(output)=={'results'},'Unexpected top-level output')
    items={i['candidate_token']:i for i in payload['candidates']}
    require(len(output['results'])==len(items) and {x['candidate_token'] for x in output['results']}==set(items),'Missing/duplicate/wrong candidate result')
    shared={f['id']:f for f in payload['shared_patient_facts']}
    for result in output['results']:
        require(set(result)==set(ITEM['properties']),'Unexpected reflection fields')
        item=items[result['candidate_token']]
        facts={**shared,**{f['id']:f for f in item['candidate_measurements']}}
        papers={a['pmid']:a for a in item['papers']}
        sentences={s['id'] for a in item['papers'] for s in a['sentences']}
        premises={p['id'] for p in item['upstream_premises']}
        def factids(ids,required=False):
            require(set(ids)<=set(facts),'Unknown/wrong-candidate patient fact ID')
            require(not required or bool(ids),'Known claim lacks patient source')
        def refs(refs):
            for ref in refs:
                require(ref['pmid'] in papers,'Unknown/wrong-candidate PMID')
                require(ref['sentence_ids'] and all(s in sentences and s.startswith(ref['pmid']+':') for s in ref['sentence_ids']),'Unknown/wrong-paper sentence ID')
        for name in ['detection','episode_site','pulmonary_event']:
            check=result[name];factids(check['patient_fact_ids'],check['status']!='unknown')
            if name=='detection' and check['status']=='supported':
                require(any(facts[i]['kind'] in ['mngs','microbiology'] for i in check['patient_fact_ids']),'Detection cites no assay source')
        require(len(result['premise_audits'])==len(premises) and {p['premise_id'] for p in result['premise_audits']}==premises,'Missing upstream premise audit')
        for audit in result['premise_audits']:factids(audit['patient_fact_ids']);refs(audit['paper_refs'])
        require(len(result['article_assessments'])==len(papers) and {a['pmid'] for a in result['article_assessments']}==set(papers),'Missing/duplicate paper assessment')
        for a in result['article_assessments']:
            factids(a['patient_fact_ids'])
            require(all(s in sentences and s.startswith(a['pmid']+':') for s in a['sentence_ids']),'Invalid article sentence ID')
            if not papers[a['pmid']]['abstract']:require(a['role']=='insufficient','Title-only article cannot be supporting')
            if a['role'] in ['support','counter']:require(bool(a['sentence_ids']),'Paper stance lacks sentence source')
            if a['role']=='support' and a['patient_applicability']=='case_relevant':factids(a['patient_fact_ids'],True)
        ca=result['case_attribution'];rf=result['reflection']
        factids(ca['patient_fact_ids'],ca['level'] in ['supported','plausible']);refs(ca['paper_refs'])
        factids(rf['patient_fact_ids'],rf['verdict']=='case_specific_revision');refs(rf['paper_refs'])
        require(set(rf['premise_ids'])<=premises,'Unknown reflection premise ID')
        require(rf['verdict']!='case_specific_revision' or bool(rf['premise_ids']),'Revision lacks target premise')
    return output['results']


def selection(result,arm,policy):
    det=result['detection'];site=result['episode_site'];event=result['pulmonary_event']
    ca=result['case_attribution'];rf=result['reflection']
    specific=ca['level'] in ['supported','plausible'] and bool(ca['patient_fact_ids']) and bool(ca['discriminating_link'].strip())
    clinical=det['status']=='supported' and site['status'] in ['supported','plausible'] and event['status'] in ['supported','plausible']
    revision=rf['verdict']=='case_specific_revision' or policy=='broad' and rf['verdict']=='uncertain'
    if policy=='strict':clinical=all(x['status']=='supported' for x in [det,site,event]) and ca['level']=='supported'
    cited={x['pmid'] for x in rf['paper_refs']}
    qualifying=[a['pmid'] for a in result['article_assessments'] if a['role']=='support'
        and a['organism_scope']=='exact_species' and a['human_context']=='causal_infection'
        and a['patient_applicability']=='case_relevant' and a['patient_fact_ids'] and a['sentence_ids'] and a['pmid'] in cited]
    paper_ok=arm=='N' or bool(qualifying)
    return {'add':bool(clinical and specific and revision and paper_ok),'clinical_gate':clinical,
        'specific_case_link':specific,'revision_basis':revision,'literature_used':bool(qualifying),
        'supporting_pmids_used':qualifying,'literature_requirement_met':paper_ok,
        'not_a_clinical_confirmation':True}


def api(folder,name,request,env):
    sha=digest(d.clinical.sources.old.encoded(request));path=folder/(name+'_raw_response.json')
    save(folder/(name+'_request.json'),request)
    if path.exists():
        record=load(path);require(record['request_sha256']==sha,'Stale request cache');return record['response']
    base=env.get('OPENAI_BASE_URL','https://api.openai.com/v1').rstrip('/')
    require(base=='https://api.openai.com/v1','Unexpected API endpoint')
    started=time.monotonic();response=json.loads(d.clinical.sources.old.http(base+'/responses',request,env['OPENAI_API_KEY'],timeout=300))
    save(path,{'request_sha256':sha,'elapsed_seconds':time.monotonic()-started,'response':response})
    return response


def run_batch(row,env):
    folder=RUN/row['folder'];payload=load(folder/'input.json')
    if row['skip']:
        save(folder/'result.json',{'results':[],'status':'no_usable_literature','not_negative_diagnosis':True});return
    if (folder/'result.json').exists():return
    errors=[]
    for replicate in range(2):
        instructions=PROMPT
        if errors:instructions+='\nCorrect only these schema/source-ID problems, preserving the same evidentiary task:\n'+json.dumps(errors)
        req={'model':MODEL,'reasoning':{'effort':'medium'},'store':False,'max_output_tokens':14000,
            'instructions':instructions,'input':json.dumps(payload,ensure_ascii=False),
            'text':{'format':{'type':'json_schema','name':'source_grounded_reflection','strict':True,'schema':SCHEMA}}}
        raw=api(folder,f'reflection_r{replicate}',req,env)
        require(raw.get('status')=='completed','Incomplete response; not a negative reflection')
        text=''.join(c.get('text','') for i in raw.get('output',[]) for c in i.get('content',[]) if c.get('type')=='output_text')
        output=json.loads(text)
        try:results=validate(output,payload)
        except (ValueError,KeyError,TypeError) as exc:
            errors=[str(exc)];save(folder/f'validation_r{replicate}.json',{'valid':False,'errors':errors})
            if replicate==1:raise
            continue
        save(folder/f'validation_r{replicate}.json',{'valid':True})
        save(folder/'result.json',{'results':results,'status':'completed','replicate_used':replicate})
        print(json.dumps({'done':row['batch_id'],'arm':row['arm'],'candidates':len(results),'repair':replicate}),flush=True);return


def run(smoke=False):
    plan=verify_inputs()
    if (RUN/'output_freeze.json').exists():verify_outputs();print('Complete frozen run verified; no API calls.');return
    env=d.clinical.environment();require(bool(env.get('OPENAI_API_KEY')),'Missing credential')
    todo=plan['batches']
    if smoke:
        first=next(b['batch_id'] for b in todo if b['arm']=='L' and not b['skip'])
        todo=[b for b in todo if b['batch_id']==first]
    errors=[]
    def work(row):
        try:run_batch(row,env)
        except Exception as exc:
            errors.append({'batch':row['batch_id'],'arm':row['arm'],'type':type(exc).__name__,'message':str(exc)[:200]})
            print(json.dumps({'failed':row['batch_id'],'arm':row['arm'],'type':type(exc).__name__,'message':str(exc)[:140]}),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:list(pool.map(work,todo))
    require(not errors,'Technical failures; do not evaluate: '+json.dumps(errors))
    if smoke:print('Paired smoke completed; no reference read.',flush=True);return
    by_token={r['token']:{} for r in plan['candidates']}
    for batch in plan['batches']:
        for result in load(RUN/batch['folder']/'result.json')['results']:
            require(batch['arm'] not in by_token[result['candidate_token']],'Duplicate reflection')
            by_token[result['candidate_token']][batch['arm']]=result
    decisions=[]
    for row in plan['candidates']:
        item={**row,'decisions':{},'literature_available':plan['usable_literature'][row['token']]}
        for arm in ['N','L']:
            if arm not in by_token[row['token']]:
                require(arm=='L' and not item['literature_available'],'Missing model reflection')
                for p in POLICIES:item['decisions'][arm+'_'+p]={'add':False,'status':'literature_unavailable_not_negative'}
                continue
            result=by_token[row['token']][arm]
            save(RUN/row['folder']/(arm+'_reflection.json'),result)
            for p in POLICIES:item['decisions'][arm+'_'+p]=selection(result,arm,p)
        decisions.append(item)
    save(RUN/'decisions_without_answers.json',decisions)
    files={p.relative_to(RUN).as_posix():digest(p.read_bytes()) for p in RUN.rglob('*') if p.is_file()
        and not p.name in ['input_freeze.json','output_freeze.json']}
    save(RUN/'output_freeze.json',{'reference_read':False,'plan_sha256':digest((RUN/'plan.json').read_bytes()),'files':files})
    print('All source-only and literature reflections/decisions frozen before evaluation.',flush=True)


def verify_outputs():
    plan=verify_inputs();freeze=load(RUN/'output_freeze.json')
    require(digest((RUN/'plan.json').read_bytes())==freeze['plan_sha256'],'Changed plan')
    for n,h in freeze['files'].items():require(digest((RUN/n).read_bytes())==h,'Changed frozen output '+n)
    return plan


def csv_save(path,rows):
    out=io.StringIO(newline='');w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write_new(path,('\ufeff'+out.getvalue()).encode('utf-8'))


def evaluate():
    plan=verify_outputs()
    from inspect_release import diagnostic_name as canon,POSITIVE_TIERS
    truths={x['case_id']:{canon(n) for n in x['reference']} for x in
        load(ROOT/'OBER_v7_review_20261030/cohort49_review.json')['cases']}
    ds=load(RUN/'decisions_without_answers.json')
    source_candidates=load(d.SOURCE/'candidates_not_model_input.json')
    manifest=load(ROOT/'OBER_v7_review_20261030/cohort49_manifest.json')['cases']
    baseline={p['case_id']:{canon(c['organism']) for c in source_candidates if c['case_id']==p['case_id'] and c['baseline_tier'] in POSITIVE_TIERS} for p in manifest}
    settings=['BASE']+['N_all_'+p for p in POLICIES]+['N_paired_'+p for p in POLICIES]+['L_'+p for p in POLICIES]
    rows=[];per=[];adds=[];predictions={}
    for setting in settings:
        for p in manifest:
            cid=p['case_id'];extra=[]
            if setting!='BASE':
                key='N_'+setting.split('_')[-1] if setting.startswith('N_') else setting
                extra=[x for x in ds if x['case_id']==cid and x['decisions'][key]['add']
                    and (not setting.startswith('N_paired') or x['literature_available'])]
            predicted=baseline[cid]|{canon(x['organism']) for x in extra};truth=truths[cid]
            predictions[(setting,cid)]=predicted
            per.append({'setting':setting,'case_id':cid,'group':p['analysis_group'],
                **d.clinical.metrics(len(predicted&truth),len(predicted-truth),len(truth-predicted))})
            for x in extra:adds.append({'setting':setting,'case_id':cid,'group':p['analysis_group'],
                'organism':x['organism'],'matches_reference':canon(x['organism']) in truth,
                'token':x['token'],'literature_available':x['literature_available']})
        for group in ['ALL49','KH33','KH5','A11']:
            local=[r for r in per if r['setting']==setting and (group=='ALL49' or r['group']==group)]
            additions=[a for a in adds if a['setting']==setting and (group=='ALL49' or a['group']==group)]
            rows.append({'setting':setting,'group':group,'cases':len(local),
                **d.clinical.metrics(*(sum(a[k] for a in local) for k in ['TP','FP','FN'])),
                'added_TP':sum(a['matches_reference'] for a in additions),'added_FP':sum(not a['matches_reference'] for a in additions)})
    require(next(r for r in rows if r['setting']=='BASE' and r['group']=='ALL49')['TP']==73,'Changed baseline contract')
    rng=random.Random(20261101);index={(r['setting'],r['case_id']):r for r in per};deltas={'L_vs_BASE':[],'L_vs_N_paired':[]}
    cids=list(truths)
    for _ in range(5000):
        sample=rng.choices(cids,k=49);score={}
        for s in ['BASE','N_paired_balanced','L_balanced']:
            score[s]=d.clinical.metrics(*(sum(index[(s,c)][k] for c in sample) for k in ['TP','FP','FN']))['F1']
        deltas['L_vs_BASE'].append(score['L_balanced']-score['BASE']);deltas['L_vs_N_paired'].append(score['L_balanced']-score['N_paired_balanced'])
    intervals={k:[sorted(v)[125],sorted(v)[4874]] for k,v in deltas.items()}
    usage=collections.Counter();repairs=0
    for batch in plan['batches']:
        for q in (RUN/batch['folder']).glob('reflection_r*_raw_response.json'):
            u=load(q)['response'].get('usage') or {};usage.update({k:u.get(k,0) for k in ['input_tokens','output_tokens','total_tokens']});usage['api_calls']+=1
        repairs+=(RUN/batch['folder']/'reflection_r1_request.json').exists()
    transitions=[]
    for x in ds:
        n=x['decisions']['N_balanced']['add'];l=x['decisions']['L_balanced']['add']
        if x['literature_available']:
            transitions.append({'case_id':x['case_id'],'organism':x['organism'],'N_add':n,'L_add':l,
                'matches_reference':canon(x['organism']) in truths[x['case_id']]})
    report={'all_metrics':rows,'primary':'L_balanced','bootstrap_delta_F1':intervals,
        'usable_literature_candidates':sum(plan['usable_literature'].values()),'candidate_count':319,
        'usage':dict(usage),'technical_repair_batches':repairs,'independent_test':False,
        'N_L_comparison_note':'Same eligible source inputs; stochastic model variation remains a possible contributor.'}
    save(RUN/'evaluation.json',report);csv_save(RUN/'metrics.csv',rows);csv_save(RUN/'case_metrics.csv',per)
    if adds:csv_save(RUN/'additions.csv',adds)
    csv_save(RUN/'paired_transitions.csv',transitions)
    for p in manifest:
        cid=p['case_id'];lines=['# '+cid+'：文獻反思','', '以下為研究性possible-source審查，不是臨床確診。','']
        for x in [x for x in ds if x['case_id']==cid]:
            folder=RUN/x['folder'];lines+=['## '+x['organism'],'']
            packet=load(folder/'packet.json');facts={f['id']:f for f in packet['facts']}
            lines+=['### 原始未選理由']+[f"- {u['id']}：{u['text']}" for u in packet['upstream_premises']]
            for arm,title in [('N','不提供文獻'),('L','提供文獻反思')]:
                lines+=['','### '+title]
                if not (folder/(arm+'_reflection.json')).exists():lines+=['文獻摘要不足，保留待審；不是感染陰性。'];continue
                a=load(folder/(arm+'_reflection.json'));rf=a['reflection']
                lines += [f"反思：{rf['verdict']}；主設定補入：{x['decisions'][arm+'_balanced']['add']}",
                    '修正依據：'+rf['what_changed'],'仍有疑慮：'+rf['remaining_uncertainty'],
                    '最強反對理由：'+rf['strongest_counterargument'],'病例連結：'+a['case_attribution']['discriminating_link']]
                for ref in rf['paper_refs']:lines.append('- 論文：[PMID '+ref['pmid']+'](https://pubmed.ncbi.nlm.nih.gov/'+ref['pmid']+'/)，句ID：'+','.join(ref['sentence_ids']))
                for fid in rf['patient_fact_ids']:
                    lines.append('- 病例來源 '+fid+'（'+facts[fid]['kind']+'／'+facts[fid]['temporal_role']+'）')
                if arm=='L':
                    lines+=['','|PMID|文獻角色|病例適用性|判讀|','|---|---|---|---|']
                    for a1 in a['article_assessments']:lines.append(f"|[{a1['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{a1['pmid']}/)|{a1['role']}|{a1['patient_applicability']}|{a1['reason'].replace('|','/')}|")
            lines+=['','反思原始JSON、病例fact文字及來源定位見此候選子資料夾。','']
        write_new(RUN/'patients'/cid/'reflection_zh.md','\n'.join(lines).encode('utf-8'))
    lines=['# OBER逐病原文獻反思：配對實驗','',
        '49例、319個Context全部納入來源與檢索；論文不足者保留未知。原97筆正式陽性不刪。',
        f"具至少一篇摘要的候選：{sum(plan['usable_literature'].values())}/319。每菌獨立最多10篇，非每患者共用10篇。",'',
        '|設定|TP|FP|FN|Precision|Recall|F1|新增TP/FP|','|---|---:|---:|---:|---:|---:|---:|---|']
    for r0 in rows:
        if r0['group']=='ALL49':lines.append(f"|{r0['setting']}|{r0['TP']}|{r0['FP']}|{r0['FN']}|{r0['precision']:.2%}|{r0['recall']:.2%}|{r0['F1']:.2%}|{r0['added_TP']}/{r0['added_FP']}|")
    lines+=['','主分析是L_balanced；N_paired使用相同具摘要候選範圍，N_all另列所有候選來源反思，不把覆蓋率差異冒稱論文效果。',
        '沒有逐菌專屬選入規則、reads/RPM閾值、DNA/RNA硬門檻或文獻票數k；只有strict/balanced/broad三組全域條件。',
        '比較是來源反思與文獻反思，並非禁止LLM判讀病例；反思必須有來源、原理由查核、逐篇適用性與反證。',
        '資料已反覆用於開發，不是未見test；bootstrap無法修正開發偏差，N/L也可能有抽樣變異。',
        '書籍／手冊條目若出現在PubMed結果，保留記錄類型，不能冒稱原始人類病例研究；本輪只讀title/abstract。', '',
        '## 配對bootstrap ΔF1（病例層級）','',json.dumps(intervals,ensure_ascii=False),'',
        '## API用量','',json.dumps(dict(usage),ensure_ascii=False),'token是保存用量，不是美元帳單。', '',
        '每病例的人讀報告見patients/<case_id>/reflection_zh.md；每候選完整N/L JSON見cases/<case_id>/<token>/。']
    write_new(RUN/'REPORT_zh.md','\n'.join(lines).encode('utf-8'))
    print(json.dumps({'ALL49':[r for r in rows if r['group']=='ALL49'],'bootstrap':intervals,'usage':dict(usage)},ensure_ascii=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','verify','smoke','run','evaluate'])
    {'prepare':prepare,'verify':verify_inputs,'smoke':lambda:run(True),'run':run,'evaluate':evaluate}[p.parse_args().command]()

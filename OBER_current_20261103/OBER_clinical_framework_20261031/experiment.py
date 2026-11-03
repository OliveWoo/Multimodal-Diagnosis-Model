"""Full49 clinical-fact extraction -> frozen deterministic rescue routes -> evaluation."""
import argparse
import collections
import concurrent.futures
import csv
import io
import json
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT/'OBER_v7_review_20261030'))
import information_pilot as sources
import audit_A1_B_causes as raw_audit
from api_probe import environment
from inspect_release import load, require, digest, write_new, POSITIVE_TIERS
import clinical_rules as rules

RUN = HERE/'runs/full49_v1_transportfix'
MODEL = sources.old.MODEL
PROMPT = (HERE/'extraction_prompt.txt').read_text(encoding='utf-8')


def save(path, value):
    write_new(path, sources.old.encoded(value))


def enum(values): return {'type': 'string', 'enum': values}
def obj(fields): return sources.old.obj(fields)
def text(): return {'type': 'string'}


EVIDENCE = obj({'fact_id': text(), 'quote': text()})
SCHEMA = obj({'case_token': text(),
    'features': {'type': 'array', 'items': obj({'feature': enum(rules.FEATURES),
        'status': enum(['present', 'absent', 'uncertain']),
        'evidence': {'type': 'array', 'items': EVIDENCE}, 'note': text()})},
    'observations': {'type': 'array', 'items': obj({'fact_id': text(), 'organism_as_written': text(),
        'organism_scope': enum(['species', 'genus', 'group', 'unspecified']),
        'method': enum(rules.METHODS), 'compartment': enum(rules.COMPARTMENTS),
        'detection': enum(['positive', 'negative', 'equivocal', 'not_result']),
        'quote': text(), 'quantity_text': text(), 'invasive_or_cytopathic_evidence': {'type': 'boolean'}})},
    'limitations': {'type': 'array', 'items': text()}})


def flatten(value, prefix=''):
    if isinstance(value, dict):
        return '\n'.join(flatten(v, prefix+k+'.') for k,v in value.items()
            if k not in ['record_id', 'source_rows', 'source_cells', 'review_required']
            and not (k=='interpretation' and v=='May overlap laboratory-source rows; not an independent confirmation or ground truth'))
    if isinstance(value, list):
        return '\n'.join(flatten(v, prefix) for v in value)
    return prefix.rstrip('.')+': '+str(value) if value is not None else ''


def normalized(value): return ' '.join(value.split())


def extractable(base):
    facts, excluded = [], []
    images = sorted((f for f in base['facts'] if f['kind']=='imaging' and f['temporal_role']!='primary_window'),
        key=lambda f: abs(f['recorded_day_relative_to_index']) if f['recorded_day_relative_to_index'] is not None else 999999)
    image_ids = {f['id'] for f in images[:3]}
    for f in base['facts']:
        v = f['value']; reason = None
        if f['kind']=='microbiology':
            if v.get('record_type')=='AST' or v.get('usable_for_result_inference') is False:
                reason = 'AST_or_no_usable_result'
            elif f['temporal_role']=='earlier_context': reason = 'earlier_test_not_primary_confirmation'
        if f['kind']=='imaging' and f['temporal_role']!='primary_window' and f['id'] not in image_ids:
            reason = 'unresolved_image_context_cap_3'
        if reason: excluded.append({'id':f['id'],'reason':reason}); continue
        facts.append({'id':f['id'],'kind':f['kind'],'temporal_role':f['temporal_role'],
            'recorded_day_relative_to_index':f['recorded_day_relative_to_index'],
            'time_caveat':f['time_caveat'],'text':flatten(v)})
    return facts, excluded


def compartment(value):
    k = rules.key(value)
    if any(x in k for x in ['balf','bronchoalveolar','lowerbal','sputum','endotracheal','trachealaspirate']) or k=='bal':
        return 'lower_respiratory'
    if any(x in k for x in ['blood','plasma','serum']): return 'blood'
    if any(x in k for x in ['nasal','nasopharyn','throat']): return 'upper_respiratory'
    return 'unknown'


def prepare():
    require(not (RUN/'plan.json').exists(), 'Already prepared; use verify/run, never overwrite')
    manifest = load(ROOT/'OBER_v7_review_20261030/cohort49_manifest.json')['cases']
    all_rows = sources.v4.all_rows()
    audit = raw_audit.audit_b()
    require(audit['counts']['raw_branch_identity_bound']==946, 'Unbound raw branch')
    by_branch = {(b['case_id'],b['organism'],b['seq_id']):b for b in audit['branches']}
    patients=[]; source_hashes={}; candidates=[]
    for rec in manifest:
        cid=rec['case_id']; base, provenance=sources.build_case(rec)
        facts, excluded=extractable(base)
        token=digest(cid.encode())[:14]
        folder=RUN/'cases'/cid
        packet={'case_token':token,'index_specimen_type':base['index_specimen_type'],
            'time_semantics':'D0 index collection; +/-2 calendar dates; timing_unknown is never primary confirmation.',
            'absence_semantics':'Missing/unprovided is unknown, not negative.', 'facts':facts}
        require(len(sources.old.encoded(packet))<400000,'Input too large; no silent truncation')
        require(not any(sources.BLOCK.search(f['text']) for f in facts),'Unsafe source block')
        save(folder/'input.json',packet); save(folder/'source_provenance.json',provenance)
        save(folder/'fact_selection.json',{'provided':len(facts),'available':len(base['facts']),'excluded':excluded})
        source_hashes.update(provenance['source_files'])
        raw=load(sources.raw_file(rec))['files']['input/mngs_raw.json']
        category={rules.key(o['name']):o['category'] for o in raw['organisms']}
        tests={t['seq_id']:t for c in raw['cases'] for t in c['tests']}
        linked={e['case_review_id']:e for e in rec['release_linked_events']}
        index=sources.old.date_of(base['index_date'])
        for row in [r for r in all_rows if r['record']['case_id']==cid]:
            c=row['candidate']; name=c['organism_name']; branches=[]
            for signal in c['analytical_evidence'].get('per_test_signals') or []:
                a=by_branch[(cid,name,signal['seq_id'])]; binding=a['raw_event_binding']
                event=linked.get(binding['raw_case_review_id'])
                require(event is not None,'Raw branch has no explicit release event')
                day=sources.old.date_of(event['collected_time'])
                require(day is not None,'Missing linked collection date')
                role='primary_window' if abs((day-index).days)<=2 else 'outside_primary_window'
                raw_site=tests[signal['seq_id']].get('parameters',{}).get('specimen_type')
                site=compartment(raw_site)
                fallback_site=compartment(event['specimen_site'])
                require(site=='unknown' or fallback_site=='unknown' or site==fallback_site,'Conflicting raw specimen')
                branches.append({'seq_id':signal['seq_id'],'compartment':site if site!='unknown' else fallback_site,
                    'specimen_raw':raw_site,'release_specimen':event['specimen_site'], 'temporal_role':role,
                    'reads':a['raw_reads'],'code_priority':a['code_priority_candidate'],
                    'percentile_all_positive':a['all_available_positive_universe']['reads_percentile_candidate'],
                    'percentile_denominator':a['all_available_positive_universe']['denominator'],
                    'rpm':a.get('rpm_recomputed'),'qc_status':signal['qc_status'],
                    'provenance':{'raw_case_review_id':binding['raw_case_review_id'],
                        'independence':'technical branch, not independent clinical confirmation'}})
            candidate={'case_id':cid,'organism':name,'organism_key':c['organism_key'],
                'baseline_tier':c['reporting']['tier'],'category':category.get(rules.key(name),'unknown'), 'branches':branches}
            candidates.append(candidate)
        patients.append({'case_id':cid,'group':rec['analysis_group'],'case_token':token,
            'folder':folder.relative_to(RUN).as_posix(),'fact_count':len(facts)})
    save(RUN/'candidates_not_model_input.json',candidates)
    save(RUN/'raw_component_audit.json',audit['counts'])
    deps=[Path(__file__),HERE/'clinical_rules.py',HERE/'extraction_prompt.txt',HERE/'PROTOCOL_zh.md',
        Path(sources.__file__),Path(sources.old.__file__),Path(sources.v4.__file__),Path(raw_audit.__file__),
        ROOT/'OBER_v7_review_20261030/inspect_release.py',ROOT/'OBER_v7_review_20261030/api_probe.py',
        ROOT/'OBER_v7_review_20261030/cohort49_manifest.json',
        ROOT/'LLM_test/data/patient_info_standardized/legacy_patient_index.json']
    save(RUN/'plan.json',{'schema':'ober.clinical_routes.v1','objective':'micro F1','independent_test':False,
        'development_labels_previously_seen':True,'answers_read_during_this_prepare':False,
        'case_count':49,'candidate_count':416,'context_count':319,'baseline_count':97,
        'model':MODEL,'reasoning':'medium','max_workers':3,'technical_repair_limit_per_case':1,
        'primary':'C3_primary','rules':rules.specifications(),'patients':patients,
        'sources':source_hashes,'dependencies':{p.relative_to(ROOT).as_posix():digest(p.read_bytes()) for p in deps},
        'new_top10_pubmed_retrieval':False,'prompt_sha256':digest(PROMPT.encode())})
    write_new(RUN/'framework_program_snapshot.py',Path(__file__).read_bytes())
    write_new(RUN/'clinical_rules_snapshot.py',(HERE/'clinical_rules.py').read_bytes())
    save(RUN/'input_freeze.json',{'answers_read':False,'files':{
        p.relative_to(RUN).as_posix():digest(p.read_bytes()) for p in RUN.rglob('*') if p.is_file()}})
    verify_inputs()
    print(json.dumps({'prepared_cases':49,'candidates':len(candidates),'arms':len(rules.specifications()),
        'context':sum(c['baseline_tier']=='Context' for c in candidates),'model':MODEL}),flush=True)


def verify_inputs():
    plan=load(RUN/'plan.json')
    for n,h in load(RUN/'input_freeze.json')['files'].items(): require(digest((RUN/n).read_bytes())==h,'Changed input: '+n)
    for n,h in {**plan['dependencies'],**plan['sources']}.items(): require(digest((ROOT/n).read_bytes())==h,'Changed source: '+n)
    return plan


def validate(output, packet):
    require(set(output)==set(SCHEMA['properties']),'Unexpected extraction fields')
    require(output['case_token']==packet['case_token'],'Wrong case token')
    require(len(output['features'])==len(rules.FEATURES) and
        {f['feature'] for f in output['features']}==set(rules.FEATURES),'Missing/duplicate fixed feature')
    facts={f['id']:f for f in packet['facts']}
    def citation(e, kind=None):
        require(e['fact_id'] in facts,'Unknown fact ID: '+e['fact_id'])
        f=facts[e['fact_id']]
        require(bool(e['quote'].strip()) and normalized(e['quote']) in normalized(f['text']),
            'Quote not exact in '+e['fact_id']+': '+e['quote'][:60])
        require(len(e['quote'])<=400,'Quote exceeds short citation limit')
        require(kind is None or f['kind'] in kind,'Wrong evidence kind for '+e['fact_id'])
    for f in output['features']:
        require(f['status'] in ['present','absent','uncertain'],'Wrong feature status')
        require(f['status']=='uncertain' or f['evidence'],'Known feature lacks evidence: '+f['feature'])
        for e in f['evidence']: citation(e,['background','imaging','index'])
    for o in output['observations']:
        citation(o,['microbiology'])
        require(o['method'] in rules.METHODS and o['compartment'] in rules.COMPARTMENTS,'Invalid observation category')
        require(not o['invasive_or_cytopathic_evidence'] or o['method']=='microscopy_histology','PCR/culture is not histology')
    return {**output,'fact_roles':{i:f['temporal_role'] for i,f in facts.items()},'traceability_valid':True}


def run_case(row, env):
    folder=RUN/row['folder'];packet=load(folder/'input.json')
    errors=[]
    for replicate in range(2):
        instructions=PROMPT
        if errors: instructions+='\nFix these source/schema validation failures only; do not infer new facts:\n'+json.dumps(errors,ensure_ascii=False)
        req={'model':MODEL,'reasoning':{'effort':'medium'},'store':False,'max_output_tokens':12000,
            'instructions':instructions,'input':json.dumps(packet,ensure_ascii=False),
            'text':{'format':{'type':'json_schema','name':'clinical_source_ledger','strict':True,'schema':SCHEMA}}}
        raw=api_call(folder,f'extraction_r{replicate}',req,env)
        require(raw.get('status')=='completed','Incomplete API response; not a negative result')
        content=''.join(c.get('text','') for item in raw.get('output',[]) for c in item.get('content',[]) if c.get('type')=='output_text')
        extracted=json.loads(content)
        try:
            ledger=validate(extracted,packet)
        except (ValueError,KeyError,TypeError) as exc:
            errors=[str(exc)];save(folder/f'extraction_r{replicate}_validation.json',{'valid':False,'errors':errors})
            if replicate==1: raise
            continue
        save(folder/f'extraction_r{replicate}_validation.json',{'valid':True})
        save(folder/'ledger.json',ledger)
        print(json.dumps({'done':row['case_id'],'features':len(ledger['features']),
            'observations':len(ledger['observations']),'technical_repairs':replicate}),flush=True)
        return


def api_call(folder, name, request, env):
    # Own transport: imported historical recovery adapters monkeypatch their
    # module globals and must not choose a cache path for this isolated run.
    sha=digest(sources.old.encoded(request));path=folder/(name+'_raw_response.json')
    save(folder/(name+'_request.json'),request)
    if path.exists():
        record=load(path);require(record['request_sha256']==sha,'Stale API cache')
        return record['response']
    base=env.get('OPENAI_BASE_URL','https://api.openai.com/v1').rstrip('/')
    require(base=='https://api.openai.com/v1','Unexpected endpoint')
    started=time.monotonic()
    raw=json.loads(sources.old.http(base+'/responses',request,env['OPENAI_API_KEY'],timeout=240))
    save(path,{'request_sha256':sha,'elapsed_seconds':time.monotonic()-started,'response':raw})
    return raw


def verify_outputs():
    plan=verify_inputs();freeze=load(RUN/'output_freeze.json')
    require(freeze['plan_sha256']==digest((RUN/'plan.json').read_bytes()),'Changed plan')
    for n,h in freeze['files'].items(): require(digest((RUN/n).read_bytes())==h,'Changed output: '+n)
    return plan


def run(smoke=False):
    plan=verify_inputs()
    if (RUN/'output_freeze.json').exists(): verify_outputs();print('Frozen complete run verified. No API calls.');return
    env=environment();require(bool(env.get('OPENAI_API_KEY')),'Missing credential')
    errors=[]
    def work(row):
        try: run_case(row,env)
        except Exception as exc:
            errors.append({'case':row['case_id'],'type':type(exc).__name__,'message':str(exc)[:250]})
            print(json.dumps({'failed':row['case_id'],'type':type(exc).__name__,'message':str(exc)[:150]}),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(work,plan['patients'][:2] if smoke else plan['patients']))
    require(not errors,'Technical failures; no evaluation: '+json.dumps(errors))
    if smoke: print('Two-case smoke complete; no labels read.',flush=True);return
    candidates=load(RUN/'candidates_not_model_input.json');decisions=[]
    for patient in plan['patients']:
        ledger=load(RUN/patient['folder']/'ledger.json')
        for c in [c for c in candidates if c['case_id']==patient['case_id']]:
            decisions.append({'case_id':c['case_id'],'organism':c['organism'],'baseline_tier':c['baseline_tier'],
                'arms':{s['id']:rules.decide(c,ledger,s) for s in plan['rules']}})
    save(RUN/'decisions_without_answers.json',decisions)
    files={p.relative_to(RUN).as_posix():digest(p.read_bytes()) for p in (RUN/'cases').rglob('*') if p.is_file()}
    files['decisions_without_answers.json']=digest((RUN/'decisions_without_answers.json').read_bytes())
    save(RUN/'output_freeze.json',{'answers_read':False,'plan_sha256':digest((RUN/'plan.json').read_bytes()),'files':files})
    print('49 cases, 416 candidates, 52 rule arms frozen BEFORE evaluation.',flush=True)


def metrics(tp,fp,fn):
    p=tp/(tp+fp) if tp+fp else 0;r=tp/(tp+fn) if tp+fn else 0
    return {'TP':tp,'FP':fp,'FN':fn,'precision':p,'recall':r,'F1':2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0}


def csv_save(path,rows):
    s=io.StringIO(newline='');w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write_new(path,('\ufeff'+s.getvalue()).encode('utf-8'))


def evaluate():
    plan=verify_outputs()
    # Only this command reads the physician answers. No scoring aliases are used
    # by extraction or the clinical observation join.
    from inspect_release import diagnostic_name as canon
    truth={r['case_id']:{canon(n) for n in r['reference']} for r in
        load(ROOT/'OBER_v7_review_20261030/cohort49_review.json')['cases']}
    decisions=load(RUN/'decisions_without_answers.json')
    base={p['case_id']:{canon(d['organism']) for d in decisions if d['case_id']==p['case_id']
        and d['baseline_tier'] in POSITIVE_TIERS} for p in plan['patients']}
    rows=[];case_rows=[];additions=[];predictions={}
    for spec in plan['rules']:
        arm=spec['id'];predictions[arm]={}
        for p in plan['patients']:
            cid=p['case_id'];extra=[d for d in decisions if d['case_id']==cid and
                d['baseline_tier']=='Context' and d['arms'][arm]['add']]
            predicted=base[cid]|{canon(d['organism']) for d in extra};predictions[arm][cid]=predicted
            t=truth[cid];m=metrics(len(predicted&t),len(predicted-t),len(t-predicted))
            case_rows.append({'arm':arm,'case_id':cid,'group':p['group'],**m})
            for d in extra: additions.append({'arm':arm,'case_id':cid,'group':p['group'],'organism':d['organism'],
                'matches_reference':canon(d['organism']) in t,'routes':'|'.join(d['arms'][arm]['routes'])})
        for group in ['ALL49','KH33','KH5','A11']:
            cases=[r for r in case_rows if r['arm']==arm and (group=='ALL49' or r['group']==group)]
            m=metrics(*(sum(r[k] for r in cases) for k in ['TP','FP','FN']))
            rows.append({'arm':arm,'group':group,'cases':len(cases),**m,
                'added_TP':sum(a['matches_reference'] for a in additions if a['arm']==arm and (group=='ALL49' or a['group']==group)),
                'added_FP':sum(not a['matches_reference'] for a in additions if a['arm']==arm and (group=='ALL49' or a['group']==group))})
    baseline=next(r for r in rows if r['arm']=='C0_baseline' and r['group']=='ALL49')
    require((baseline['TP'],baseline['FP'],baseline['FN'])==(73,24,14),'Baseline scoring contract changed')
    primary=next(r for r in rows if r['arm']==plan['primary'] and r['group']=='ALL49')
    grid=sorted([r for r in rows if r['arm'].startswith('G_') and r['group']=='ALL49'],key=lambda r:(-r['F1'],r['arm']))
    rng=random.Random(20261031);delta=[];patients=[p['case_id'] for p in plan['patients']]
    per={(r['arm'],r['case_id']):r for r in case_rows}
    for _ in range(5000):
        sample=rng.choices(patients,k=len(patients));scores=[]
        for arm in ['C0_baseline',plan['primary']]:
            m=metrics(*(sum(per[(arm,c)][k] for c in sample) for k in ['TP','FP','FN']));scores.append(m['F1'])
        delta.append(scores[1]-scores[0])
    delta.sort();ci=[delta[125],delta[4874]]
    use=collections.Counter();validity=collections.Counter()
    for p in plan['patients']:
        f=RUN/p['folder'];ledger=load(f/'ledger.json')
        validity.update({'valid_cases':1,'observations':len(ledger['observations'])})
        for q in f.glob('extraction_r*_raw_response.json'):
            resp=load(q)['response'];usage=resp.get('usage') or {};use.update({k:usage.get(k,0) for k in ['input_tokens','output_tokens','total_tokens']});use['api_calls']+=1
        md=['# '+p['case_id']+' — clinical evidence routes','', '研究用；非臨床確診。既有名單保持不變。', '', '## 原始事實抽取']
        md += [f"- {x['feature']}: {x['status']}；{x['note']}；來源 "+','.join(e['fact_id'] for e in x['evidence']) for x in ledger['features']]
        md += ['', '## 主分析 C3 的候選判定', '', '|病原|原始層級|是否補入|路徑|未滿足／警示|', '|---|---|---|---|---|']
        for d in [d for d in decisions if d['case_id']==p['case_id']]:
            q=d['arms'][plan['primary']]
            md.append(f"|{d['organism']}|{d['baseline_tier']}|{q['add'] if d['baseline_tier']=='Context' else '原樣保留'}|{','.join(q['routes'])}|{','.join(q['blocking_or_missing'])}|")
        write_new(f/'REPORT_zh.md','\n'.join(md).encode('utf-8'))
    result={'primary':primary,'baseline':baseline,'exploratory_best':grid[0],
        'primary_delta_F1_bootstrap_95_percentile':ci,'bootstrap_unit':'case/episode; not proven independent natural persons; development only',
        'all_rows':rows,'usage':dict(use),'validity':dict(validity),'independent_test':False}
    save(RUN/'evaluation.json',result);csv_save(RUN/'metrics_all_52_arms.csv',rows)
    csv_save(RUN/'case_metrics.csv',case_rows)
    if additions: csv_save(RUN/'all_additions.csv',additions)
    head=['# 全49例 OBER Clinical Evidence Routes 實跑結果','',
        '52臂（4個層遞＋48組預先凍結敏感度）；49例、416候選、319 Context均完成。',
        '不覆寫舊版；LLM僅抽事實，沒有診斷／是否選入欄位。主分析C3在答案揭露前固定。', '',
        '|臂|TP|FP|FN|Precision|Recall|F1|新增TP/FP|','|---|---:|---:|---:|---:|---:|---:|---|']
    selected=[r for r in rows if r['group']=='ALL49' and r['arm'].startswith('C')]+grid[:5]
    for r in selected: head.append(f"|{r['arm']}|{r['TP']}|{r['FP']}|{r['FN']}|{r['precision']:.2%}|{r['recall']:.2%}|{r['F1']:.2%}|{r['added_TP']}/{r['added_FP']}|")
    head += ['',f"主分析 ΔF1={primary['F1']-baseline['F1']:+.4f}；病例層級 bootstrap 95%區間 [{ci[0]:+.4f}, {ci[1]:+.4f}]。", '',
        '網格最高值是探索結果，不是獨立test表現；病例資料已反覆用於開發。區間不消除參數探索偏差，同一自然人不同事件也可能相關。',
        '主分析若沒有改善必須如實呈現，不可用網格最佳取代預先設定的C3。不同Candida species仍分開。', '',
        '評分沿用專案既有回溯名稱契約，包含部分群層級／同義名稱映射；不是完全純exact species。臨床檢驗配對使用另一份更嚴格別名表。', '',
        '## 跨cohort主分析','', '|群組|TP|FP|FN|Precision|Recall|F1|','|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        if r['arm']==plan['primary']: head.append(f"|{r['group']}|{r['TP']}|{r['FP']}|{r['FN']}|{r['precision']:.2%}|{r['recall']:.2%}|{r['F1']:.2%}|")
    head += ['', '## 主分析新增項目', '']
    head += [f"- {a['case_id']} / {a['organism']}：{'TP' if a['matches_reference'] else 'FP'}；{a['routes']}" for a in additions if a['arm']==plan['primary']]
    head += ['', '## API與資料品質', '',json.dumps({'usage':dict(use),'validity':dict(validity)},ensure_ascii=False),
        'token是本地保存用量，不是美元帳單；包含必要的技術引文修復。引用能定位不保證醫學解讀正確。',
        '沒有新的PubMed top-10/A1判讀；外部依據為固定臨床指引。所有補入是研究性possible，不宣稱已確診。', '',
        '詳見PROTOCOL_zh.md、plan.json、input/output_freeze.json、每病例input/ledger/source_provenance及REPORT_zh.md。']
    write_new(RUN/'REPORT_zh.md','\n'.join(head).encode('utf-8'))
    print(json.dumps({'baseline':baseline,'primary':primary,'exploratory_best':grid[0],'usage':dict(use)},ensure_ascii=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','verify','smoke','run','evaluate'])
    command=parser.parse_args().command
    {'prepare':prepare,'verify':verify_inputs,'smoke':lambda:run(True),'run':run,'evaluate':evaluate}[command]()

"""Verify completed reruns, evaluate with unchanged scoring, compare old/new/OBER."""
import json
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
import direct_raw_runner as r
import evaluate_direct_raw as ev
from input_identity import sha256_file
from prepare_corrected_reruns import RUNS,EXPECTED_CHANGED

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
REPORT=BASE/'reports/idfixed_20260906'
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,obj): p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def pct(n): return f'{n*100:.2f}%'

def validate_runs():
    _,specs=r.load_manifest(r.DEFAULT_MANIFEST)
    audits={s.patient_id:r.build_patient_payload(s)[1] for s in specs}
    prompt_hash=r.sha256_text(r.DEFAULT_PROMPT.read_text(encoding='utf-8-sig').strip())
    verifications=[]
    for model,(old_name,new_name) in RUNS.items():
        new_dir=r.DEFAULT_OUTPUT_ROOT/new_name
        summary=read(new_dir/'run_summary.json')
        expected_ids={s.patient_id for s in specs}
        if summary.get('failures') or set(summary['success_patient_ids'])!=expected_ids:
            raise ValueError(f'{model}: run is incomplete or has failures')
        if {p.stem.removesuffix('_direct_raw') for p in new_dir.glob('P*_direct_raw.json')}!=expected_ids:
            raise ValueError(f'{model}: output cohort differs from manifest')
        for name in ['direct_raw_runner.py','input_identity.py','prompt_v1.txt','input_mapping_v2.json','cohort_manifest.json']:
            if sha256_file(new_dir/'source_snapshots'/name)!=sha256_file(BASE/name):
                raise ValueError(f'{model}: current source differs from frozen snapshot: {name}')
        reuse=read(new_dir/'reuse_provenance.json')
        if set(reuse['must_rerun'])!=EXPECTED_CHANGED: raise ValueError('Rerun scope changed')
        reused={row['patient_id']:row for row in reuse['reused']}
        if len(reused)!=19: raise ValueError('Incorrect reuse count')
        for spec in specs:
            result_file=new_dir/f'{spec.patient_id}_direct_raw.json'
            result=read(result_file)
            r.validate_cached_output(result,spec,audits[spec.patient_id],prompt_hash,f'gpt-5.6-{model}','medium')
            old_file=r.DEFAULT_OUTPUT_ROOT/old_name/f'{spec.patient_id}_direct_raw.json'
            if spec.patient_id in reused:
                frozen=reused[spec.patient_id]
                if sha256_file(old_file)!=frozen['source_sha256'] or sha256_file(result_file)!=frozen['source_sha256']:
                    raise ValueError('Original/reused output bytes changed')
            else:
                if result.get('input_identity')!=audits[spec.patient_id]['input_identity']:
                    raise ValueError('New output does not attest to current input mapping')
                if result['payload_sha256']==read(old_file)['payload_sha256']:
                    raise ValueError('Affected input unexpectedly unchanged')
            verifications.append({'model':model,'patient_id':spec.patient_id,'reused':spec.patient_id in reused,'output_sha256':sha256_file(result_file),'payload_sha256':result['payload_sha256']})
    return verifications

def evaluate(run_dir,out_prefix):
    args=[sys.executable,str(BASE/'evaluate_direct_raw.py'),'--run-dir',str(run_dir),'--out-json',str(REPORT/(out_prefix+'.json')),'--out-md',str(REPORT/(out_prefix+'.md'))]
    done=subprocess.run(args,capture_output=True,text=True,encoding='utf-8')
    if done.returncode: raise RuntimeError(done.stderr or done.stdout)
    return read(REPORT/(out_prefix+'.json'))

def main():
    checks=validate_runs()  # No incomplete cohort can be scored silently.
    REPORT.mkdir(parents=True,exist_ok=True)
    compared={}
    for model,(old,new) in RUNS.items():
        before=evaluate(r.DEFAULT_OUTPUT_ROOT/old,f'{model}_before')
        after=evaluate(r.DEFAULT_OUTPUT_ROOT/new,f'{model}_after')
        changes=[]
        prior={x['patient_id']:x for x in before['per_patient']}
        for row in after['per_patient']:
            old_row=prior[row['patient_id']]
            if row['gold']!=old_row['gold']: raise ValueError('Gold changed')
            if row['cohort'].startswith('two_hospital') and row!=old_row: raise ValueError('Unchanged two-hospital result changed')
            if 'P'+row['patient_id'] in EXPECTED_CHANGED:
                changes.append({'patient_id':row['patient_id'],'gold':row['gold'],'old_predicted':old_row['predicted'],'new_predicted':row['predicted'],'old_tp':len(old_row['tp']),'new_tp':len(row['tp']),'old_fp':len(old_row['fp']),'new_fp':len(row['fp']),'old_fn':len(old_row['fn']),'new_fn':len(row['fn'])})
        compared[model]={'before':before['metrics'],'after':after['metrics'],'affected_22_before':ev.summarize([x for x in before['per_patient'] if 'P'+x['patient_id'] in EXPECTED_CHANGED]),'affected_22_after':ev.summarize([x for x in after['per_patient'] if 'P'+x['patient_id'] in EXPECTED_CHANGED]),'changed_patients':changes}
        model_checks=[x for x in checks if x['model']==model]
        new_ids=[x['patient_id'] for x in model_checks if not x['reused']]
        reused_ids=[x['patient_id'] for x in model_checks if x['reused']]
        new_outputs=[read(r.DEFAULT_OUTPUT_ROOT/new/f'{pid}_direct_raw.json') for pid in new_ids]
        if len(new_ids)!=22 or len(reused_ids)!=19: raise ValueError('Unexpected batch totals')
        write(r.DEFAULT_OUTPUT_ROOT/new/'batch_completion.json',{
            'schema_version':'direct_raw_benchmark.corrected_batch_completion.v1',
            'verified_at_utc':datetime.now(timezone.utc).isoformat(),
            'status':'complete', 'model':f'gpt-5.6-{model}',
            'verified_patient_count':len(model_checks),
            'new_response_count':len(new_ids),'new_response_patient_ids':new_ids,
            'legacy_reused_count':len(reused_ids),'legacy_reused_patient_ids':reused_ids,
            'new_response_usage':{k:sum(x.get('usage',{}).get(k,0) for x in new_outputs) for k in ['input_tokens','output_tokens','total_tokens']},
            'models_returned':sorted({x['model_returned'] for x in new_outputs}),
            'note':'Batch-wide accounting, including the Luna P14 pilot. run_summary.json counts only the most recent invocation; validated reuse within this corrected batch is not legacy reuse. Token totals cover saved successful new responses only, not possible failed requests.',
        })
        # Store a convenient current evaluation beside its patient results.
        write(r.DEFAULT_OUTPUT_ROOT/new/'evaluation_exact.json',after)
        (r.DEFAULT_OUTPUT_ROOT/new/'evaluation_exact_zh.md').write_text(ev.render_markdown(after),encoding='utf-8')
    gold_files=[Path(p) for p in after['gold_sources']]
    aliases=ev.build_aliases([read(p) for p in gold_files])
    ober_rows=[]
    for row in after['per_patient']:
        cohort='KH' if row['cohort'].startswith('kmuh') else 'two_hospital'
        source=ROOT/f"OBER_patient_delivery/outputs/20260905_r5_final_019feb4d/patients/{cohort}/P{row['patient_id']}/evidence/final_decision.json"
        decision=read(source)
        predicted={ev.canonical(x['organism'],aliases) for x in decision['selected']}
        gold=set(row['gold']); tp=gold&predicted; fp=predicted-gold; fn=gold-predicted
        precision=len(tp)/len(predicted) if predicted else 0; recall=len(tp)/len(gold) if gold else 0
        ober_rows.append({'patient_id':row['patient_id'],'cohort':row['cohort'],'gold':sorted(gold),'predicted':sorted(predicted),'tp':sorted(tp),'fp':sorted(fp),'fn':sorted(fn),'precision':precision,'recall':recall,'f1':ev.f1(precision,recall),'source_sha256':sha256_file(source)})
    ober={c:ev.summarize([x for x in ober_rows if x['cohort']==c]) for c in ['kmuh_answer_positive_30','two_hospital_patient_quality_A_11']}
    ober['all_41']=ev.summarize(ober_rows)
    result={'schema_version':'direct_raw_benchmark.idfix_comparison.v1','matching':'unchanged_exact_species_plus_frozen_synonyms','models':compared,'OBER_R5_same_41':ober,'OBER_per_patient':ober_rows,'verification':checks,'source_hashes':{'mapping':sha256_file(BASE/'input_mapping_v2.json'),'prompt':r.sha256_text(r.DEFAULT_PROMPT.read_text(encoding='utf-8-sig').strip()),'evaluator':sha256_file(BASE/'evaluate_direct_raw.py'),'gold':{p.name:sha256_file(p) for p in gold_files}},'notes':['22 patients rerun per model; 19 byte-identical verified prior responses reused.','Old results used mismatched mNGS on 22 KH cases; old/new is an input-correction comparison, not a model upgrade.','OBER was not rerun; existing final decisions rescored on identical 41 patients and Gold.','Broad/non-species labels and frozen synonyms unchanged; this is not an untouched clinical test set.','KH P33 age discrepancy and P32 missing demographics remain; clinical records were not modified.']}
    write(REPORT/'comparison.json',result)
    lines=['# Direct-Raw 編號修正後重跑結果','','只修正病例／mNGS 配對；同一份提示詞、medium reasoning、4000 output-token cap、同一評分程式及 Gold。每模型重跑 22 位、驗證沿用 19 位。','','## 合併 41 位：修正前後','','| 模型 | Precision 前 → 後 | Recall 前 → 後 | F1 前 → 後 | TP 前 → 後 | FP 前 → 後 | FN 前 → 後 |','|---|---:|---:|---:|---:|---:|---:|']
    for model,value in compared.items():
        b=value['before']['all_41']; a=value['after']['all_41']
        lines.append(f"| {model.capitalize()} | {pct(b['precision'])} → {pct(a['precision'])} | {pct(b['recall'])} → {pct(a['recall'])} | {pct(b['f1'])} → {pct(a['f1'])} | {b['tp']} → {a['tp']} | {b['fp']} → {a['fp']} | {b['fn']} → {a['fn']} |")
    for cohort,label in [('kmuh_answer_positive_30','高醫 30 位'),('two_hospital_patient_quality_A_11','Two-hospital A 11 位'),('all_41','合併 41 位')]:
        lines+=['',f'## 修正後：{label}','','| 方法 | Gold | 預測 | TP | FP | FN | Precision | Recall | F1 |','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for name,m in [(model.capitalize(),v['after'][cohort]) for model,v in compared.items()]+[('OBER R5（原決策，未重跑）',ober[cohort])]:
            lines.append(f"| {name} | {m['gold_pair_count']} | {m['predicted_pair_count']} | {m['tp']} | {m['fp']} | {m['fn']} | {pct(m['precision'])} | {pct(m['recall'])} | {pct(m['f1'])} |")
    lines+=['','## 結果來源','','- 所有 123 份回答均核對目前 payload、prompt、模型、reasoning 與病人 ID。','- 沿用的 57 份回答與舊檔 byte-for-byte 相同；新 66 份有明確 input_identity。','- 每批 batch_completion.json 統計完整批次的 22 新回答與 19 舊回答沿用；Luna P14 先行試跑也計入 22 份。run_summary.json 僅統計最後一次命令。','- OBER 使用原本 R5 決策，未修改門檻、選菌或理由。','- 修正前數字只供錯配影響分析，不是有效的同輸入模型比較。','- 此為既有 answer-positive／人工篩選資料的內部比較，不能當作外部獨立測試證明。','- 分數是與 recorded clinical-answer labels 的一致性，不代表臨床確診；不自動修正 P33 年齡等既有資料爭議。','','## 每模型逐病人明細','','- [Luna 修正後](luna_after.md)','- [Terra 修正後](terra_after.md)','- [Sol 修正後](sol_after.md)','- [機器可讀比較與逐病人變化](comparison.json)','']
    (REPORT/'comparison_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps({'new_metrics':{k:v['after'] for k,v in compared.items()},'OBER_R5_same_41':ober,'report':str(REPORT/'comparison_zh.md')},ensure_ascii=False,indent=2))

if __name__=='__main__': main()

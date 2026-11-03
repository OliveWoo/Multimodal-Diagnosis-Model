"""Post-hoc large-module coupling. No species/case-specific selection rules."""
import argparse
import collections
import math
from pathlib import Path
import reason_common as c

NRUN=c.HERE/'runs/native_signal_v1'
POLICIES=['absolute','relative','either','both']


def flags(analytical,burden):
    signals=analytical.get('per_test_signals') or []
    accepted=set(burden.get('accepted_qc_statuses') or [])
    threshold=burden.get('per_test_rpm_threshold')
    valid=[]
    for s in signals:
        value=s.get('rpm_total')
        if isinstance(value,(int,float)) and math.isfinite(value) and value>=0 and s.get('qc_status') in accepted:
            valid.append(value)
    absolute=None
    if signals and isinstance(threshold,(int,float)) and math.isfinite(threshold):
        if any(x>=threshold for x in valid) and burden.get('all_selected_tests_qc_accepted') is True:absolute=True
        elif len(valid)==len(signals):absolute=False
    band=analytical.get('rank_band')
    relative=True if band in ['top_1','top_3'] else False if band in ['top_10','below_top_10'] else None
    return {'absolute':absolute,'relative':relative,'native_rpm_threshold':threshold,
        'native_rank_band':band,'positive_branch_count':len(signals),'valid_rpm_branch_count':len(valid)}


def passes(f,policy):
    a=f['absolute'] is True;b=f['relative'] is True
    return {'absolute':a,'relative':b,'either':a or b,'both':a and b}[policy]


def prepare():
    c.verify('judge_inputs');c.verify('judge_outputs')
    c.require(not (NRUN/'plan.json').exists(),'Native module already prepared')
    rows=c.load(c.RUN/'source_plan.json')['candidates']
    upstream={(x['record']['case_id'],x['candidate']['organism_name']):x for x in c.old.d.clinical.sources.v4.all_rows()}
    components={(x['case_id'],x['organism']):x for x in c.load(c.old.d.SOURCE/'candidates_not_model_input.json')}
    feature={};profiles=[]
    for row in rows:
        key=(row['case_id'],row['organism']);u=upstream[key];d=u['detail'];ap=d.get('analytical_profile') or {}
        burden=((d.get('promotion_gate') or {}).get('axes') or {}).get('analytically_compelling_burden_profile') or {}
        ids={s['seq_id'] for s in ap.get('per_test_signals') or []}
        retained={s['seq_id'] for s in components[key]['branches']}
        c.require(ids==retained,'Native analytical profile uses different episode branches')
        f=flags(ap,burden);feature[row['token']]=f
        profiles.append({k:row[k] for k in ['case_id','group','token','organism']}|f|{
            'aligned_seq_ids':sorted(ids),'native_burden':burden,'source_result_file':u['record']['patient_result_file']})
    predictions=[]
    for row in c.load(c.RUN/'decisions_without_answers.json'):
        f=feature[row['token']];broad=row['decisions']['broad']['add']
        predictions.append({k:row[k] for k in ['case_id','group','token','organism','arm']}|{
            'original_broad_add':broad,'add':{p:bool(broad and passes(f,p)) for p in POLICIES}})
    c.save(NRUN/'native_features_not_model_input.json',profiles)
    c.save(NRUN/'predictions_without_answers.json',predictions)
    c.save(NRUN/'plan.json',{'posthoc_after_primary_results_seen':True,'independent_test':False,
        'primary_exploratory_direction':'T_either','policies':POLICIES,'no_new_model_calls':True,
        'source_judgement_freeze_sha256':c.digest((c.RUN/'judge_outputs_freeze.json').read_bytes()),
        'script_sha256':c.digest(Path(__file__).read_bytes()),
        'protocol_sha256':c.digest((c.HERE/'NATIVE_SIGNAL_PROTOCOL_zh.md').read_bytes()),
        'reference_read_during_prepare':False,'all_native_seq_ids_match_retained_episode':True})
    c.save(NRUN/'prediction_freeze.json',{'files':{p.name:c.digest(p.read_bytes()) for p in NRUN.iterdir() if p.is_file()}})
    print({'frozen_candidates':len(profiles),'policies':POLICIES,'no_new_API_calls':True,
        'absolute_flags':dict(collections.Counter(str(x['absolute']) for x in profiles)),
        'relative_flags':dict(collections.Counter(str(x['relative']) for x in profiles))})


def evaluate():
    c.verify('judge_outputs');plan=c.load(NRUN/'plan.json')
    c.require(plan['script_sha256']==c.digest(Path(__file__).read_bytes()),'Changed native module code')
    c.require(plan['protocol_sha256']==c.digest((c.HERE/'NATIVE_SIGNAL_PROTOCOL_zh.md').read_bytes()),'Changed native module protocol')
    c.require(plan['source_judgement_freeze_sha256']==c.digest((c.RUN/'judge_outputs_freeze.json').read_bytes()),'Changed source judgement freeze')
    for n,h in c.load(NRUN/'prediction_freeze.json')['files'].items():c.require(c.digest((NRUN/n).read_bytes())==h,'Changed native predictions')
    from inspect_release import diagnostic_name as canon,POSITIVE_TIERS
    refs=c.load(c.ROOT/'OBER_v7_review_20261030/cohort49_review.json')['cases']
    truth={p['case_id']:{canon(n) for n in p['reference']} for p in refs}
    source=c.load(c.old.d.SOURCE/'candidates_not_model_input.json')
    base={cid:{canon(x['organism']) for x in source if x['case_id']==cid and x['baseline_tier'] in POSITIVE_TIERS} for cid in truth}
    rows=c.load(NRUN/'predictions_without_answers.json');groups={x['case_id']:x['group'] for x in rows}
    settings=['BASE']+[a+'_'+p for a in ['G','T'] for p in ['broad']+POLICIES]
    per={};pred={}
    def score(setting,cids):return c.old.d.clinical.metrics(*(sum(per[setting,c][k] for c in cids) for k in ['TP','FP','FN']))
    for setting in settings:
        for cid,t in truth.items():
            extra=set()
            if setting!='BASE':
                arm,policy=setting.split('_')
                extra={canon(x['organism']) for x in rows if x['case_id']==cid and x['arm']==arm
                    and (x['original_broad_add'] if policy=='broad' else x['add'][policy])}
            p=base[cid]|extra;pred[setting,cid]=p
            per[setting,cid]=c.old.d.clinical.metrics(len(p&t),len(p-t),len(t-p))
    metrics=[]
    for setting in settings:
        for group in ['ALL49','KH33','KH5','A11']:
            cids=[cid for cid in truth if group=='ALL49' or groups[cid]==group];m=score(setting,cids);b=score('BASE',cids)
            metrics.append({'setting':setting,'group':group,**m,'net_added_TP':m['TP']-b['TP'],'net_added_FP':m['FP']-b['FP']})
    choices=['BASE','T_both','T_either','T_absolute','T_relative'];folds=[];cv={}
    for held in sorted(truth):
        train=[cid for cid in truth if cid!=held]
        def key(setting):
            added=sum(len(pred[setting,cid]-base[cid]) for cid in train)
            return score(setting,train)['F1'],-added,-choices.index(setting)
        selected=max(choices,key=key);cv[held]=per[selected,held]
        folds.append({'held_case':held,'selected_setting':selected,'training_F1':score(selected,train)['F1']})
    cvscore=c.old.d.clinical.metrics(*(sum(x[k] for x in cv.values()) for k in ['TP','FP','FN']))
    c.old.csv_save(NRUN/'metrics.csv',metrics);c.old.csv_save(NRUN/'episode_cv.csv',folds)
    result={'metrics':metrics,'leave_one_episode_out':cvscore,
        'selected_counts':dict(collections.Counter(x['selected_setting'] for x in folds)),
        'posthoc_development_not_independent_test':True}
    c.save(NRUN/'evaluation.json',result)
    lines=['# 文獻反思＋原生分析訊號：事後大模組對照','',
        '此為主實驗計分後的探索，不取代T_balanced主結果。沒有新API呼叫，也沒有按病人／菌種訂規則。',
        '每個原生seq_id集合均已和時間窗對齊來源核對；不使用病原family-specific eligible欄位。','',
        '|設定|TP|FP|FN|Precision|Recall|F1|淨增TP／FP|','|---|---:|---:|---:|---:|---:|---:|---:|']
    for x in metrics:
        if x['group']=='ALL49':lines.append(f"|{x['setting']}|{x['TP']}|{x['FP']}|{x['FN']}|{x['precision']:.2%}|{x['recall']:.2%}|{x['F1']:.2%}|{x['net_added_TP']}／{x['net_added_FP']}|")
    lines+=['','absolute沿用上游每測試RPM旗標及QC；relative沿用top_1/top_3排名區間；either為OR，both為AND。',
        '上游RPM值不被宣稱是臨床致病閾值。未知不等於陰性，但無法由對應軸提供正向補入支持。',
        'DNA/RNA多分支不被計為多次獨立臨床確認；跨測試reads不加總。','',
        '## Leave-one-episode-out穩定性檢查','',
        f"Precision {cvscore['precision']:.2%}；Recall {cvscore['recall']:.2%}；F1 {cvscore['F1']:.2%}。",
        '每折只用其他48例選BASE或T的四種組合；選擇分布：'+str(result['selected_counts']),
        '沒有證明episode皆屬不同自然人，這不是獨立test，也不能消除先前反覆開發的偏差。',
        '全49例的最高分只能作探索；所有四組完整呈現，不把較好的組合事後換成已驗證主假說。','',
        '完整設定與cohort表見runs/native_signal_v1/metrics.csv；旗標來源見native_features_not_model_input.json。']
    c.write_new(c.HERE/'NATIVE_SIGNAL_RESULTS_zh.md','\n'.join(lines).encode('utf-8'))
    print({'ALL49':[x for x in metrics if x['group']=='ALL49'],'episode_cv':cvscore,'selected_counts':result['selected_counts']})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','evaluate'])
    {'prepare':prepare,'evaluate':evaluate}[p.parse_args().command]()

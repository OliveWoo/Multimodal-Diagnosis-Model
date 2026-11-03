"""All49 development scoring with QC coverage and unresolveds kept visible."""
import collections
import json
from pathlib import Path
import qc_common as s
import qc_runner as runner
import qc_rules as rules
import qc_review as review


def usage():
    totals=collections.Counter();models=collections.Counter();seen=set()
    for p in s.RUN.rglob('*_raw_response.json'):
        r=s.load(p)['response'];rid=r.get('id')
        if rid in seen:continue
        seen.add(rid);totals['saved_responses']+=1;models[r.get('model','unknown')]+=1
        for k in ['input_tokens','output_tokens','total_tokens']:totals[k]+=(r.get('usage') or {}).get(k,0)
    old_ids=set()
    for j in s.load(s.RUN/'A1_jobs.json'):
        if j['reused'] and j['has_abstract']:
            for p in Path(j['folder']).glob('*_raw_response.json'):old_ids.add(s.load(p)['response'].get('id'))
    out={'new_saved_response_usage':dict(totals),'new_actual_models':dict(models),'reused_old_A1_saved_responses_not_new_charges':len(old_ids),
        'new_transport_incidents':len(list(s.RUN.rglob('transport_error_*.json')))+len(list(s.RUN.rglob('http_body_recovery_*.json'))),
        'not_invoice':True,'interrupted_response_billing_may_be_unrecorded':True}
    s.save(s.RUN/'api_usage.json',out);return out


def main():
    for stage in ['inputs','QC_outputs','article_decisions','A1_inputs','A1_outputs','predictions']:s.verify(stage)
    plan=s.load(s.RUN/'plan.json');features=s.load(s.RUN/'predictions_without_answers.json');specs=runner.configs()
    s.require(s.load(s.RUN/'rules.json')==specs,'Changed rules')
    decisions=s.load(s.RUN/'article_decisions_without_answers.json');articles=[];repeats=[]
    for a in plan['articles']:
        f=s.RUN/a['folder'];p=s.load(f/'source.json');s.require(rules.route(p)==s.load(f/'routing.json'),'Nondeterministic rule routing')
        out=decisions[a['key']]
        if a['route']=='QC':
            normalized=[s.normalize(s.c.validate_saved(f/f'r{i}',review.validate),p) for i in range(2)]
            s.require(out['normalized_pair']==normalized,'Changed normalized QC response')
            s.require(out['strategies']['QC_CONSENSUS']==s.consensus(*normalized),'Changed QC consensus')
            repeats.append({'key':a['key'],'pmid':a['pmid'],'organism':a['organism'],
                'same_decision_relevant_signature':normalized[0]['signature']==normalized[1]['signature'],
                'r0_status':normalized[0]['status'],'r1_status':normalized[1]['status'],
                'r0_basis':normalized[0]['basis'],'r1_basis':normalized[1]['basis'],
                'consensus_status':out['strategies']['QC_CONSENSUS']['status']})
        c=out['strategies']['QC_CONSENSUS'];articles.append({**a,'primary_status':c['status'],'primary_basis':c['basis'],
            'qualified':c['qualified'],'human_any':c['human_any'],'rule_reasons':'|'.join(a['reasons']),
            'unresolved_reason':'|'.join(c['internal_issues'])})
    for r in features:
        for spec in specs:
            v=r.get('variants',{}).get(f"{spec['strategy']}_F{spec['fit']}");expected=False
            if v and v['S_pass']:
                counts=v['counts'];s.require(counts is not None,'Missing A1 must not be zero')
                score=min(counts[a]['strong_match_count']+(counts[a]['partial_match_count'] if spec['count']=='match' else 0) for a in ['G','T'])
                expected=spec['k']==0 or score>=spec['k']
            s.require(r['decisions'][spec['id']]==expected,'Changed final decision')
    # First reference read in this runner, after fixed predictions; historical
    # labels have been seen in earlier development and are not a new test.
    from inspect_release import diagnostic_name as canon,POSITIVE_TIERS
    ref=s.ROOT/'OBER_v7_review_20261030/cohort49_review.json';refs=s.load(ref)['cases']
    truth={r['case_id']:{canon(n) for n in r['reference']} for r in refs};groups={r['case_id']:r['analysis_group'] for r in refs}
    clinical=s.load(s.c.old.d.SOURCE/'candidates_not_model_input.json')
    base={cid:{canon(r['organism']) for r in clinical if r['case_id']==cid and r['baseline_tier'] in POSITIVE_TIERS} for cid in truth}
    s.require(len(truth)==49 and sum(map(len,truth.values()))==87 and sum(map(len,base.values()))==97,'Changed scoring denominator')
    parent=s.load(s.PARENT/'predictions_without_answers.json')
    controls={'BASE':[],
        'OLD_A1_S':[r for r in parent if r['decisions']['OLD_F0_match_k2']],
        'PREVIOUS_SHARED_R0':[r for r in parent if r['decisions']['SHARED_L0_F0_match_k2']],
        'PREVIOUS_SHARED_R1':[r for r in parent if r['decisions']['SHARED_L1_F0_match_k2']]}
    choices=controls|{x['id']:[r for r in features if r['decisions'][x['id']]] for x in specs}
    metrics=[];per={};added=[];metric=s.c.old.d.clinical.metrics
    for name,picks in choices.items():
        for cid,t in truth.items():
            p=base[cid]|{canon(r['organism']) for r in picks if r['case_id']==cid}
            s.require(base[cid]<=p,'Original formal prediction pruned');per[name,cid]=metric(len(p&t),len(p-t),len(t-p))
        for g in ['ALL49','KH33','KH5','A11']:
            cids=[cid for cid in truth if g=='ALL49' or groups[cid]==g]
            m=metric(*(sum(per[name,cid][k] for cid in cids) for k in ['TP','FP','FN']))
            b=metric(*(sum(per['BASE',cid][k] for cid in cids) for k in ['TP','FP','FN']))
            metrics.append({'setting':name,'group':g,'cases':len(cids),**m,'net_added_TP':m['TP']-b['TP'],'net_added_FP':m['FP']-b['FP']})
        for r in picks:added.append({'setting':name,'case_id':r['case_id'],'organism':r['organism'],
            'reference_match':canon(r['organism']) in truth[r['case_id']],'new_evaluation_pair':canon(r['organism']) not in base[r['case_id']]})
    idx={(x['setting'],x['group']):x for x in metrics}
    s.require((idx['OLD_A1_S','ALL49']['TP'],idx['OLD_A1_S','ALL49']['FP'],idx['OLD_A1_S','ALL49']['FN'])==(76,24,11),'Old control changed')
    dist=[]
    for r in features:
        for prefix,v in r.get('variants',{}).items():
            c=v['counts'];dist.append({k:r[k] for k in ['case_id','group','organism','token']}|{'variant':prefix,'S_pass':v['S_pass'],
                'reference_match':canon(r['organism']) in truth[r['case_id']],'needs_source_review':v['needs_source_review'],
                'G_match':None if c is None else c['G']['strong_match_count']+c['G']['partial_match_count'],
                'T_match':None if c is None else c['T']['strong_match_count']+c['T']['partial_match_count'],
                'G_pending_QC_pmids':'|'.join(v['arm_gates']['G']['pending_QC_articles']),
                'T_pending_QC_pmids':'|'.join(v['arm_gates']['T']['pending_QC_articles'])})
    coverage={'article_count':58,'rule_accepted':sum(a['route']=='rule' for a in articles),
        'QC_routed':sum(a['route']=='QC' for a in articles),'primary_unresolved_QC':sum(a['primary_status']=='unresolved_qc' for a in articles),
        'QC_pair_decision_signature_disagreements':sum(not x['same_decision_relevant_signature'] for x in repeats),
        'QC_pairs':len(repeats),'primary_qualified_articles':sum(a['qualified'] is True and a['primary_status']!='unresolved_qc' for a in articles),
        'QC_is_not_sequencing_quality':True}
    result={'primary':idx[plan['primary'],'ALL49'],'controls':{k:idx[k,'ALL49'] for k in controls},'metrics':metrics,'coverage':coverage,
        'reference_sha256':s.digest(ref.read_bytes()),'independent_test':False,'fixed_case_replay':True,
        'original_formal_pairs_preserved':97,'same_model_QC_consensus_not_independent_experts':True,
        'unresolved_not_removed_from_evaluation_denominator':True,'QC_repeat':repeats}
    s.save(s.RUN/'evaluation.json',result);s.c.old.csv_save(s.RUN/'metrics.csv',metrics)
    s.c.old.csv_save(s.RUN/'candidate_distribution.csv',dist);s.c.old.csv_save(s.RUN/'article_QC_distribution.csv',articles)
    s.c.old.csv_save(s.RUN/'QC_repeat.csv',repeats);s.c.old.csv_save(s.RUN/'additions.csv',added);use=usage()
    patients(features,clinical,truth,POSITIVE_TIERS)
    md=['# 規則＋文獻QC：實驗結果','',
        '固定舊病例判斷的49例開發replay，87答案配對，原97筆正式結果不刪。不是完整新版部署效能或獨立test。QC是文獻語義待查核，不是定序品質。','',
        '## 先看規則與QC覆蓋','', '```json',json.dumps(coverage,ensure_ascii=False,indent=2),'```','',
        '規則只自動接受可定位明確正向句；找不到關鍵詞不判陰性。QC兩次不同或關鍵語義未知留待查，不追加第三次抽到較有利結果。','',
        '## 固定主設定與預定對照','', '|設定|TP|FP|FN|Precision|Recall|F1|新增TP／FP|','|---|---:|---:|---:|---:|---:|---:|---:|']
    names=['BASE','OLD_A1_S','PREVIOUS_SHARED_R0','RULE_ONLY_F0_match_k2','QC_R0_F0_match_k2','QC_R1_F0_match_k2',plan['primary'],'QC_DIRECT_ONLY_F0_match_k2','QC_CONSENSUS_F1_match_k2']
    for name in names:
        x=idx[name,'ALL49'];md.append(f"|{name}|{x['TP']}|{x['FP']}|{x['FN']}|{x['precision']:.2%}|{x['recall']:.2%}|{x['F1']:.2%}|{x['net_added_TP']}／{x['net_added_FP']}|")
    md += ['', 'F0/F1是固定原A1六個病例匹配欄位的兩份既有／按原prompt抽取結果；QC_R0/R1是本輪兩次文獻語義判讀。舊G/T病例審查未重抽。','',
        'direct-only對照只撤下疾病相關支持的文獻資格，不改其人類致病角色與病例欄位。因此差異不能全說成QC降低隨機性；其中也包含文獻支持定義的修正。','',
        '## 主設定的k敏感度（主k2不換）','', '|k|TP／FP／FN|Precision|Recall|F1|','|---:|---|---:|---:|---:|']
    for k in range(11):
        x=idx[f'QC_CONSENSUS_F0_match_k{k}','ALL49'];md.append(f"|{k}|{x['TP']}／{x['FP']}／{x['FN']}|{x['precision']:.2%}|{x['recall']:.2%}|{x['F1']:.2%}|")
    md += ['', '## 主設定新增','', '|病例|病原|與固定答案配對|','|---|---|---|']
    for x in added:
        if x['setting']==plan['primary']:md.append(f"|{x['case_id']}|{x['organism']}|{'TP' if x['reference_match'] else 'FP'}|")
    md += ['', '## 限制與判讀','',
        '- 兩次同模型可以同錯；來源ID驗證、固定規則與快取均不等於醫學語義驗證。未證明全流程變異下降。',
        '- source QC的signature只涵蓋決策相關的類別，不要求措辭或引用句順序完全相同；部分不影響已成立支持鏈的機制未知仍保留。',
        '- 存在未解QC的病例仍留在49例分母。尚未補入會反映在FN，但這是系統的行動結果，不是醫學上判定病原不存在。',
        '- 疾病相關肺部併發症支持possible-source情境，不偽裝成已證明直接肺部病因；明確其他病原的繼發感染不自動變成原候選直接病因。',
        '- 文獻相關修正僅能影響5個固定病例閘門容許的候選；這個小範圍不能突破舊版已達到的正確補入上限。',
        '- 沒有病例／病原特例，舊版不覆寫，沒有因F1更換主設定或剔除答案。','',
        '## 新呼叫用量','', '```json',json.dumps(use,ensure_ascii=False,indent=2),'```','',
        '逐例資料在runs/qc_v1/patients/；來源QC明細在article_QC_distribution.csv，全部設定在metrics.csv。']
    s.write_new(s.HERE/'RESULTS_zh.md','\n'.join(md).encode('utf-8'))
    s.freeze('evaluation',[s.RUN/'evaluation.json',s.RUN/'metrics.csv',s.RUN/'api_usage.json',s.HERE/'RESULTS_zh.md'])
    print(json.dumps({'primary':result['primary'],'coverage':coverage,'single_r0':idx['QC_R0_F0_match_k2','ALL49'],
        'single_r1':idx['QC_R1_F0_match_k2','ALL49'],'direct_only':idx['QC_DIRECT_ONLY_F0_match_k2','ALL49'],'usage':use},ensure_ascii=False),flush=True)


def patients(features,clinical,truth,positive):
    article_records=s.load(s.RUN/'article_decisions_without_answers.json');plan=s.load(s.RUN/'plan.json');rowmap={r['token']:r for r in plan['potential']}
    for cid in sorted(truth):
        case=[r for r in features if r['case_id']==cid];original=[r['organism'] for r in clinical if r['case_id']==cid and r['baseline_tier'] in positive]
        added=[r['organism'] for r in case if r['decisions'][plan['primary']]];dest=s.RUN/'patients'/cid
        pending=[r['organism'] for r in case if r.get('variants',{}).get('QC_CONSENSUS_F0',{}).get('needs_source_review')]
        s.save(dest/'research_selection.json',{'case_id':cid,'setting':plan['primary'],'original_formal_preserved':original,
            'research_added':added,'unresolved_source_review_candidates':pending,'research_replay_only':True})
        md=[f'# {cid}：規則＋文獻QC','', '文獻QC不是檢體品質；未解QC不等於沒有感染。這是固定病例判斷的研究replay。','',
            '- 原正式清單：'+('、'.join(original) or '無'),'- 主設定補入：'+('、'.join(added) or '無'),
            '- 仍含待查文獻的候選：'+('、'.join(pending) or '無'),'', '|候選|流程狀態|主設定補入|','|---|---|---|']
        for r in case:md.append(f"|{r['organism']}|{r['status']}|{r['decisions'][plan['primary']]}|")
        for r in case:
            if not r.get('variants'):continue
            v=r['variants']['QC_CONSENSUS_F0'];md += ['', '## '+r['organism'],'',f"S={v['S_pass']}；仍含待查來源={v['needs_source_review']}。"]
            for a in ['G','T']:
                gate=v['arm_gates'][a];md += ['', a+'合格來源：'+json.dumps(gate['qualified_articles'],ensure_ascii=False),
                    a+'待查來源：'+('、'.join(gate['pending_QC_articles']) or '無')]
            row=rowmap[r['token']]
            for pmid,key in row['article_keys'].items():
                record=article_records[key];x=record['strategies']['QC_CONSENSUS'];p=s.load(s.RUN/record['folder']/'source.json')
                md += ['',f"### [{pmid}：{p['article']['title']}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)", '',
                    '分流：'+record['route']+'；狀態：'+x['status']+'；支持型別：'+x['basis'],
                    '','規則分流原因：'+'、'.join(record['reasons'])]
                if x['status']=='unresolved_qc':md += ['', '未解QC：'+'、'.join(x['internal_issues'])]
                for name,summary in x['source_summaries'].items():md.append('- '+name+'：'+summary)
                ids=list(dict.fromkeys(x['qualifying_ids']+x['human_ids']))
                byid={z['id']:z['text'] for z in p['article']['sentences']}
                if ids:md += ['', '來源原句（不是模型生成引文）：']+['- ['+i+'] '+byid[i] for i in ids]
        s.write_new(dest/'REPORT_zh.md','\n'.join(md).encode('utf-8'))


if __name__=='__main__':main()

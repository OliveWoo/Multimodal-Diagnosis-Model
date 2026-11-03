"""Compare genuinely new model responses, semantics, QC states and selections."""
import collections
import json
from pathlib import Path
import repeat_driver as d

s=d.s
SIGNATURE_FIELDS=['organism_scope','human_any','qualified','basis','nonhuman_only','noncausal','secondary_other']


def profile(record):
    if record['status']=='unresolved_qc':return ['unresolved_qc']
    return [record['status']]+[record.get(k) for k in SIGNATURE_FIELDS]


def changes(left,right):return [name for i,name in enumerate(SIGNATURE_FIELDS) if left[i]!=right[i]]


def selected(rows,setting):return {(r['case_id'],r['organism']) for r in rows if r['decisions'][setting]}


def actual(folder):
    saved=s.load(folder/'result.json');tag=saved['selected_response_tag']
    first=s.load(folder/'r0_request.json');chosen=s.load(folder/(tag+'_request.json'))
    raw=s.load(folder/(tag+'_raw_response.json'))['response']
    s.c.validate_saved(folder,d.review.validate)
    return {'first_request_sha256':s.digest(s.c.old.d.clinical.sources.old.encoded(first)),
        'selected_request_sha256':s.digest(s.c.old.d.clinical.sources.old.encoded(chosen)),
        'response_id':raw['id'],'selected_response_tag':tag,'model':raw.get('model'),
        'created_at':raw.get('created_at'),'system_fingerprint':raw.get('system_fingerprint'),
        'output':saved['output']}


def usage_summary():
    summary={};totals=collections.Counter()
    for label,root in [('round_A',d.ROUNDS['A']),('round_B',d.ROUNDS['B']),('common_new_A1',d.COMMON)]:
        ct=collections.Counter();models=collections.Counter();ids=set()
        for p in root.rglob('*_raw_response.json'):
            r=s.load(p)['response'];rid=r.get('id')
            if rid in ids:continue
            ids.add(rid);ct['saved_responses']+=1;models[r.get('model','unknown')]+=1
            u=r.get('usage') or {}
            for k in ['input_tokens','output_tokens','total_tokens']:ct[k]+=u.get(k,0)
            ct['prompt_cached_input_tokens']+=(u.get('input_tokens_details') or {}).get('cached_tokens',0)
        summary[label]={'usage':dict(ct),'actual_models':dict(models)}
        totals.update(ct)
    summary['all_new_saved_responses']=dict(totals)
    summary['not_invoice']=True;summary['prompt_cache_is_not_local_model_answer_reuse']=True
    return summary


def compare():
    d.verify_file_manifest(d.BUNDLE/'setup_freeze.json');d.verify_file_manifest(d.BUNDLE/'both_predictions_freeze.json')
    d.verify_file_manifest(d.PARENT/'QC_outputs_freeze.json')
    plans={label:s.load(folder/'plan.json') for label,folder in d.ROUNDS.items()}
    articles={label:s.load(folder/'article_decisions_without_answers.json') for label,folder in d.ROUNDS.items()}
    predictions={label:s.load(folder/'predictions_without_answers.json') for label,folder in d.ROUNDS.items()}
    evaluations={label:s.load(folder/'evaluation.json') for label,folder in d.ROUNDS.items()}
    for label,folder in d.ROUNDS.items():
        for stage in ['inputs','QC_outputs','article_decisions','A1_inputs','A1_outputs','predictions','evaluation']:d.verify_file_manifest(folder/(stage+'_freeze.json'))
    s.require(plans['A']['QC_jobs']==plans['B']['QC_jobs'],'Two blocks use different QC slot lists')
    audits=[];slot_rows=[];models={label:collections.Counter() for label in ['A','B']};new_ids={label:set() for label in ['A','B']}
    parent_ids=set()
    for p in d.PARENT.glob('articles/*/r*/*_raw_response.json'):parent_ids.add(s.load(p)['response'].get('id'))
    field_counts=collections.Counter();all_original_requests_same=True;selected_requests_same=True
    for job in plans['A']['QC_jobs']:
        actuals={label:actual(folder/job['folder']) for label,folder in d.ROUNDS.items()}
        old_first=s.load(d.PARENT/job['folder']/'r0_request.json')
        old_hash=s.digest(s.c.old.d.clinical.sources.old.encoded(old_first))
        a,b=actuals['A'],actuals['B']
        identical=a['first_request_sha256']==b['first_request_sha256']==old_hash
        s.require(identical,'Initial model requests changed across blocks or versus frozen baseline')
        all_original_requests_same=all_original_requests_same and identical
        selected_requests_same=selected_requests_same and a['selected_request_sha256']==b['selected_request_sha256']
        for label,value in actuals.items():
            s.require(value['response_id'] not in new_ids[label],'Duplicate selected Response ID within one new block')
            new_ids[label].add(value['response_id']);models[label][value['model']]+=1
        key=job['token'];i=job['replicate'];aa=articles['A'][key]['normalized_pair'][i];bb=articles['B'][key]['normalized_pair'][i]
        fields=changes(aa['signature'],bb['signature']);field_counts.update(fields)
        slot_rows.append({'article_key':key,'pmid':articles['A'][key]['pmid'],'organism':articles['A'][key]['organism'],
            'slot':i,'same_raw_JSON':a['output']==b['output'],'same_decision_signature':aa['signature']==bb['signature'],
            'same_normalized_status':aa['status']==bb['status'],'changed_fields':'|'.join(fields),
            'A_signature':json.dumps(aa['signature'],ensure_ascii=False),'B_signature':json.dumps(bb['signature'],ensure_ascii=False)})
        audits.append({'article_key':key,'slot':i,'identical_initial_request_to_baseline':identical,
            'A':{k:v for k,v in a.items() if k!='output'},'B':{k:v for k,v in b.items() if k!='output'}})
    s.require(not new_ids['A']&new_ids['B'],'New blocks share selected Response IDs')
    s.require(not (new_ids['A']|new_ids['B'])&parent_ids,'A new block reused a baseline Response ID')
    all_raw_ids={}
    for label,root in d.ROUNDS.items():
        all_raw_ids[label]={s.load(p)['response'].get('id') for p in root.glob('articles/*/r*/*_raw_response.json')}
    s.require(not all_raw_ids['A']&all_raw_ids['B'] and not (all_raw_ids['A']|all_raw_ids['B'])&parent_ids,'Raw API responses reused across blocks')
    article_rows=[]
    for key,a in articles['A'].items():
        b=articles['B'][key];x=a['strategies']['QC_CONSENSUS'];y=b['strategies']['QC_CONSENSUS']
        xa=x['status']=='unresolved_qc';ya=y['status']=='unresolved_qc';both_resolved=not xa and not ya
        article_rows.append({'article_key':key,'pmid':a['pmid'],'organism':a['organism'],'route':a['route'],
            'A_status':x['status'],'B_status':y['status'],'QC_status_flip':xa!=ya,'both_pending_QC':xa and ya,
            'both_resolved':both_resolved,'same_resolved_signature':x['signature']==y['signature'] if both_resolved else None,
            'final_article_profile_changed':profile(x)!=profile(y),'A_qualified':x['qualified'],'B_qualified':y['qualified'],
            'A_basis':x['basis'],'B_basis':y['basis']})
    setting=plans['A']['primary'];s.require(setting==plans['B']['primary'],'Different primary setting')
    pick={label:selected(rows,setting) for label,rows in predictions.items()}
    case_rows=[]
    by={label:{r['token']:r for r in rows} for label,rows in predictions.items()}
    for row in plans['A']['potential']:
        a=by['A'][row['token']];b=by['B'][row['token']];va=a['variants']['QC_CONSENSUS_F0'];vb=b['variants']['QC_CONSENSUS_F0']
        count=lambda v,arm:None if v['counts'] is None else v['counts'][arm]['strong_match_count']+v['counts'][arm]['partial_match_count']
        case_rows.append({'case_id':row['case_id'],'organism':row['organism'],'A_selected':a['decisions'][setting],
            'B_selected':b['decisions'][setting],'selection_changed':a['decisions'][setting]!=b['decisions'][setting],
            'A_S':va['S_pass'],'B_S':vb['S_pass'],'A_G_match':count(va,'G'),'A_T_match':count(va,'T'),
            'B_G_match':count(vb,'G'),'B_T_match':count(vb,'T'),
            'A_pending_QC':sum(len(x) for x in va['pending_A1_QC_articles'].values()),
            'B_pending_QC':sum(len(x) for x in vb['pending_A1_QC_articles'].values())})
    curves=[]
    for k in range(11):
        name=f'QC_CONSENSUS_F0_match_k{k}';picks={label:selected(rows,name) for label,rows in predictions.items()}
        scores={label:next(x for x in evaluations[label]['metrics'] if x['setting']==name and x['group']=='ALL49') for label in ['A','B']}
        curves.append({'k':k,'A_F1':scores['A']['F1'],'B_F1':scores['B']['F1'],
            'same_selected_pairs':picks['A']==picks['B'],'changed_case_organism_pairs':len(picks['A']^picks['B'])})
    tally={'new_selected_QC_responses_A':len(new_ids['A']),'new_selected_QC_responses_B':len(new_ids['B']),
        'identical_initial_request_bodies':all_original_requests_same,'identical_selected_request_bodies':selected_requests_same,
        'response_ID_sets_disjoint_from_each_other_and_baseline':True,
        'cross_block_single_slot_pairs':len(slot_rows),'cross_block_raw_JSON_differences':sum(not x['same_raw_JSON'] for x in slot_rows),
        'cross_block_decision_signature_differences':sum(not x['same_decision_signature'] for x in slot_rows),
        'cross_block_signature_field_changes':dict(field_counts),
        'final_article_profiles_changed':sum(x['final_article_profile_changed'] for x in article_rows),
        'QC_status_flips':sum(x['QC_status_flip'] for x in article_rows),'both_pending_QC':sum(x['both_pending_QC'] for x in article_rows),
        'both_resolved_articles':sum(x['both_resolved'] for x in article_rows),
        'both_resolved_but_signature_changed':sum(x['same_resolved_signature'] is False for x in article_rows),
        'same_primary_selected_pairs':pick['A']==pick['B'],'primary_selected_pair_symmetric_difference':len(pick['A']^pick['B'])}
    api=usage_summary()
    result={'scope':'two fresh QC consensus blocks with fixed G/T and A1; not full OB ER reruns',
        'tally':tally,'round_A':{'primary':evaluations['A']['primary'],'coverage':evaluations['A']['coverage'],'selected':sorted(pick['A'])},
        'round_B':{'primary':evaluations['B']['primary'],'coverage':evaluations['B']['coverage'],'selected':sorted(pick['B'])},
        'only_A':sorted(pick['A']-pick['B']),'only_B':sorted(pick['B']-pick['A']),
        'usage':api,'actual_models':{k:dict(v) for k,v in models.items()},'independent_test':False,'full_pipeline_stability_not_established':True}
    s.save(d.BUNDLE/'comparison.json',result);s.save(d.BUNDLE/'fresh_API_audit.json',audits)
    s.c.old.csv_save(d.BUNDLE/'single_QC_slot_changes.csv',slot_rows);s.c.old.csv_save(d.BUNDLE/'consensus_article_changes.csv',article_rows)
    s.c.old.csv_save(d.BUNDLE/'candidate_changes.csv',case_rows);s.c.old.csv_save(d.BUNDLE/'k_curve_comparison.csv',curves)
    md=['# 兩個全新QC輪次：重跑比較','',
        '固定同一批病例判斷、A1與文獻；A、B各自重新執行57篇×2次QC，形成兩組共識。不是上次兩個答案，也不是只重讀快取。這不是全部OBER從raw data端到端重跑。','',
        '## 兩輪主設定k2','', '|輪次|TP|FP|FN|Precision|Recall|F1|輪內兩次QC分歧|未解QC|','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for label in ['A','B']:
        x=evaluations[label]['primary'];c=evaluations[label]['coverage']
        md.append(f"|{label}|{x['TP']}|{x['FP']}|{x['FN']}|{x['precision']:.2%}|{x['recall']:.2%}|{x['F1']:.2%}|{c['QC_pair_decision_signature_disagreements']}／{c['QC_pairs']}|{c['primary_unresolved_QC']}|")
    md += ['', '## 跨輪非確定性','',
        f"- 對應單次QC slot共有{len(slot_rows)}對：原始JSON不同{tally['cross_block_raw_JSON_differences']}對，決策相關signature不同{tally['cross_block_decision_signature_differences']}對。措辭不同與分類不同分開計算。",
        f"- 58篇最終文章結果有{tally['final_article_profiles_changed']}篇改變；其中接受／QC狀態翻轉{tally['QC_status_flips']}篇，兩輪都接受但分類不同{tally['both_resolved_but_signature_changed']}篇。",
        f"- 兩輪都保留QC有{tally['both_pending_QC']}篇；這只表示兩次都未解，不視為醫學判讀已一致。",
        f"- 主設定最終病例—病原清單相同：{tally['same_primary_selected_pairs']}；集合對稱差{tally['primary_selected_pair_symmetric_difference']}筆。F1與清單有分開查核。",'',
        '## 有沒有真的新呼叫','',
        f"- A/B各{len(new_ids['A'])}／{len(new_ids['B'])}個新選用Response ID，彼此及舊版本完全互斥。",
        f"- 每個初始API request body都與另一輪及原凍結版本相同：{all_original_requests_same}。最終選用request body相同：{selected_requests_same}（若有格式修補／容量追加會另記）。",
        '- 它們可能使用平台prompt cache，不能因此說重用了模型答案；本地沒有複製QC result/raw_response。實際IDs、request hash、model、時間與fingerprint（若有）在fresh_API_audit.json。', '',
        '## 主設定補入比較','', '|病例|病原|A|B|G/T match A|G/T match B|','|---|---|---|---|---|---|']
    for r in case_rows:md.append(f"|{r['case_id']}|{r['organism']}|{r['A_selected']}|{r['B_selected']}|{r['A_G_match']}／{r['A_T_match']}|{r['B_G_match']}／{r['B_T_match']}|")
    md += ['', '## k敏感度，不更換主k2','', '|k|A F1|B F1|清單相同|','|---:|---:|---:|---|']
    for r in curves:md.append(f"|{r['k']}|{r['A_F1']:.2%}|{r['B_F1']:.2%}|{r['same_selected_pairs']}|")
    md += ['', '## 解讀界線','',
        '- 相同配置出現分類差異，是本次觀察到的操作非確定性；不能把相同總分解讀成模型完全沒有變動。',
        '- QC共識不是把模型變成確定性；它可能把不一致轉成待查，但兩輪對哪些文章待查也可能變動。',
        '- G/T、A1、檢索與上游選菌固定，所以沒有證明這些層穩定，也沒有估計完整OBER端到端的隨機性。',
        '- 同一模型名稱未必代表可確認相同的後端快照；保留實際model/fingerprint資訊，無法觀察的後端差異不冒稱已排除。',
        '- 只測兩輪、49例反覆開發資料，不足證明長期穩定機率或臨床正確性，並非獨立test。',
        '- 沒有重抽第三輪、挑較好輪次、改規則或依答案換k。兩輪無答案選擇全凍結後才評分。','',
        '## 實際新用量','', '```json',json.dumps(api,ensure_ascii=False,indent=2),'```']
    path=d.HERE/'COMPARISON_zh.md';s.write_new(path,'\n'.join(md).encode('utf-8'))
    d.freeze_bundle('comparison',[Path(__file__),path,d.BUNDLE/'comparison.json',d.BUNDLE/'fresh_API_audit.json',d.BUNDLE/'single_QC_slot_changes.csv',d.BUNDLE/'consensus_article_changes.csv',d.BUNDLE/'candidate_changes.csv',d.BUNDLE/'k_curve_comparison.csv'])
    print(json.dumps({'tally':tally,'A':evaluations['A']['primary'],'B':evaluations['B']['primary'],'usage':api['all_new_saved_responses']},ensure_ascii=False),flush=True)


if __name__=='__main__':compare()

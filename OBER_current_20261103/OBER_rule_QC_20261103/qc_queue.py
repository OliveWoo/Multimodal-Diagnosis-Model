"""Human-readable unresolved article queue. No labels or new model calls."""
from pathlib import Path
import qc_common as s

FIELDS=['物種解析度','人類病因支持','文獻資格','支持型別','純非人類研究','非因果情境','其他病原繼發感染']
NAMES={'human_any':'人類疾病病因','direct_pulmonary':'直接肺部病因','disease_related':'疾病相關肺部併發症',
       'noncausal':'定植／再活化等非因果情境','secondary_other':'其他病原的繼發感染'}


def main():
    s.verify('article_decisions');s.verify('predictions')
    records=s.load(s.RUN/'article_decisions_without_answers.json');pred=s.load(s.RUN/'predictions_without_answers.json')
    plan=s.load(s.RUN/'plan.json');rows=plan['potential'];queue=[]
    for key,item in sorted(records.items()):
        c=item['strategies']['QC_CONSENSUS']
        if c['status']!='unresolved_qc':continue
        affected=[];gate_impact=False
        for row in rows:
            if key not in row['article_keys'].values():continue
            p=next(x for x in pred if x['token']==row['token']);variant=p['variants']['QC_CONSENSUS_F0']
            lanes=[a for a in ['G','T'] if item['pmid'] in row['arm_pmids'][a]]
            blocks=[a for a in lanes if not variant['arm_gates'][a]['article_support_gate'] and item['pmid'] in variant['arm_gates'][a]['pending_QC_articles']]
            gate_impact=gate_impact or bool(blocks)
            affected.append({'case_id':row['case_id'],'organism':row['organism'],'lanes':lanes,
                'could_affect_current_unmet_article_gate':blocks,'current_primary_added':p['decisions'][plan['primary']]})
        f=s.RUN/item['folder'];raw=[s.load(f/f'r{i}/result.json')['output'] for i in range(2)]
        sig=c.get('paired_signatures',[])
        changed=[name for i,name in enumerate(FIELDS) if len(sig)==2 and sig[0][i]!=sig[1][i]]
        queue.append({'key':key,'pmid':item['pmid'],'organism':item['organism'],
            'priority':'可能影響目前未過的文獻閘門' if gate_impact else '其他文章語義待查，不一定改變目前選擇',
            'reason':'關鍵分類不一致：'+'、'.join(changed) if changed else '兩次皆保留關鍵語義未知',
            'affected_cases':affected,'source_snapshot':str((f/'source.json').resolve()),'paired_raw_QC':raw})
    queue.sort(key=lambda x:(not x['priority'].startswith('可能影響'),x['key']))
    s.save(s.RUN/'human_QC_queue.json',{'no_answers_used_to_rank':True,'third_LLM_call_not_requested':True,'articles':queue})
    md=['# 文獻QC待查清單','',
        '本清單只整理已凍結的模型分歧／不明，不改分類、選擇或分數，不新增API。排序只看是否可能影響目前未過的文獻閘門，不看醫院答案。QC不是「病原不致病」。','',
        '共9篇：8篇關鍵分類不一致、1篇兩次都留未知。沒有第三次抽到一致。應核對原文／來源，必要時由人工裁決；同模型一致也可能同錯。','',
        '每篇建議確認：是否有同解析度的人類致病主張？肺部支持是直接歸因、疾病併發症、僅檢出，還是其他病原繼發感染？原文能支持到哪一層？不要看病例答案決定文章分類。']
    for x in queue:
        md += ['', f"## [{x['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{x['pmid']}/) — {x['organism']}", '',
            '- 優先情況：'+x['priority'],'- 待查原因：'+x['reason'],
            '- 涉及病例：'+'；'.join(a['case_id']+'（'+','.join(a['lanes'])+'）' for a in x['affected_cases']),
            '- 原文快照：'+x['source_snapshot']]
        for name,label in NAMES.items():
            md += ['', '### '+label]
            for i,q in enumerate(x['paired_raw_QC']):
                part=q[name];sourceids=part.get('sentence_ids',part.get('agent_disease_sentence_ids',[])+part.get('disease_lung_sentence_ids',[]))
                details='；level='+part['level'] if 'level' in part else '；chain='+part['chain'] if 'chain' in part else ''
                md += ['', f"- r{i}：{part['status']}{details}；{part['summary']}",
                    '  - 原句ID：'+('、'.join(sourceids) or '無支持句／該次未作正向主張')]
    path=s.HERE/'QC_REVIEW_zh.md';s.write_new(path,'\n'.join(md).encode('utf-8'))
    s.freeze('human_QC_queue',[Path(__file__),path,s.RUN/'human_QC_queue.json'])
    print({'unresolved_articles':len(queue),'potential_article_gate_impact':sum(x['priority'].startswith('可能影響') for x in queue),
        'no_new_API_or_reference_access':True},flush=True)


if __name__=='__main__':main()

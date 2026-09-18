import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';

export const readJson = async p => JSON.parse((await fs.readFile(p, 'utf8')).replace(/^\uFEFF/, ''));
export const sha256 = data => crypto.createHash('sha256').update(data).digest('hex');
export const jsonText = value => JSON.stringify(value, null, 2) + '\n';
export const key = name => String(name).trim().toLowerCase().replace(/\s+/g, ' ');
export const str = value => value == null ? '未提供' : typeof value === 'string' ? value : JSON.stringify(value);
export const arr = x => Array.isArray(x) ? x : [];
export const stageLabel = stage => stage === 'workflow_final' ? '工作流最終結果' : '上游 picked（非 OBER 最終結果）';

export function assertSelection(expected, explanations) {
  const encode = items => items.map(x => [x.organism, x.rank]);
  for (const [label, items] of [['decision', expected], ['rationale', explanations]]) {
    if (new Set(items.map(x => x.organism)).size !== items.length) throw Error(`${label}: duplicate organism`);
    if (items.some(x => typeof x.organism !== 'string' || !Number.isInteger(x.rank) || x.rank < 1)) throw Error(`${label}: invalid organism/rank`);
    if (new Set(items.map(x => x.rank)).size !== items.length) throw Error(`${label}: duplicate rank`);
  }
  if (JSON.stringify(encode([...expected].sort((a,b)=>a.rank-b.rank))) !== JSON.stringify(encode([...explanations].sort((a,b)=>a.rank-b.rank)))) {
    throw Error(`Final selection and reason-generator output differ. Expected ${JSON.stringify(encode(expected))}; got ${JSON.stringify(encode(explanations))}. Regenerate reasons from the correct frozen final decision; no silent fallback.`);
  }
}

export function verifyPatient(id, ...records) {
  if (!/^\d+$/.test(String(id))) throw Error('Invalid patient identifier');
  for (const record of records) if (String(record.patient_id ?? record['病例識別']?.['病人編號']) !== String(id)) throw Error('Cross-patient input mismatch');
}

const verdictZh = value => ({strong_match:'強匹配',partial_match:'部分匹配',insufficient:'資訊不足',retrieval_mismatch:'檢索不符',mismatch:'檢索／病例不匹配',counterevidence:'反向證據',neutral_non_support:'中性／不支持'}[value] ?? value ?? '未提供');

export async function loadCasefits(root, id, selectedPackage) {
  if (!root) return {entries: [], rag: null};
  const manifest = await readJson(path.join(root, 'manifest.json'));
  const entries = [];
  for (const row of manifest.outputs.filter(x => String(x.patient_id) === String(id))) {
    const file = path.resolve(root, row.output_file);
    if (path.dirname(file) !== path.resolve(root)) throw Error('Unsafe casefit path');
    const bytes = await fs.readFile(file);
    if (row.output_sha256 && sha256(bytes) !== row.output_sha256) throw Error(`Casefit file hash mismatch: ${file}`);
    const data = JSON.parse(bytes.toString('utf8').replace(/^\uFEFF/, ''));
    verifyPatient(id, data);
    if (data.provenance.latest_merge_sha256 !== selectedPackage.source.sha256) throw Error(`Casefit/selection source versions differ: P${id}`);
    entries.push({file, data});
  }
  let rag = null;
  if (entries.length) {
    const reference = entries[0].data.provenance;
    const file = path.join(manifest.frozen_rag_dir, reference.frozen_rag_file);
    const bytes = await fs.readFile(file);
    if (sha256(bytes) !== reference.frozen_rag_sha256) throw Error(`Frozen RAG hash mismatch: P${id}`);
    const data = JSON.parse(bytes.toString('utf8').replace(/^\uFEFF/, ''));
    verifyPatient(id, data);
    // A1 can explicitly reuse articles retrieved by an older RAG run while
    // judging a newer case card. Bibliography is reusable, old patient modules
    // are not. The A1 latest-merge hash above still must match exactly.
    rag = {file, data, sameInput:data.input_meta.input_sha256 === selectedPackage.source.sha256};
  }
  return {entries, rag};
}

export function makeRecord({cohort, result, chain, technical, selection, casefits, decision, stage, stageNote}) {
  const id = String(result.patient_id);
  verifyPatient(id, result, chain, technical, selection);
  const gen = result.generated_rationale;
  verifyPatient(id, gen);
  const expected = decision ? decision.selected : arr(selection.selected_pathogens).map(x=>({organism:x.canonical_name, rank:x.rank}));
  if (stage === 'workflow_final' && (!decision || decision.stage !== 'workflow_final')) throw Error('workflow_final requires explicit final_decision.json for every patient');
  if (decision) verifyPatient(id, decision);
  assertSelection(expected, gen.selected_pathogen_explanations);
  if (gen.no_selected_pathogen !== (expected.length === 0)) throw Error('No-pick flag mismatch');
  const sourceCandidates = arr(chain['步驟2至5_逐候選病原判讀']);
  const chosen = sourceCandidates.filter(x=>x['最終處置']==='選入');
  if (!decision) assertSelection(expected, chosen.map(x=>({organism:x['病原'],rank:x['最終排序']})));
  const candidates = sourceCandidates.map((c,i)=>({
    '候選ID':`C${String(i+1).padStart(3,'0')}`, '病原':c['病原'], '候選來源':'上游候選清單',
    '原上游處置':c['最終處置'], '最終處置':expected.some(x=>x.organism===c['病原'])?'選入':'未選入',
    '最終排序':expected.find(x=>x.organism===c['病原'])?.rank ?? null,
    '為什麼考慮':arr(c['為什麼會被考慮']), '為什麼選入或排除':arr(c['為什麼選入或排除']),
    '理由來源':'既有上游規則／判讀紀錄（非本次新生成）', '關鍵判斷資訊':c['關鍵判斷資訊'] ?? {},
    '限制':arr(c['重要限制']), '病例證據':arr(c['病例內支持或反對證據']).map((e,j)=>({
      '證據ID':`C${String(i+1).padStart(3,'0')}-E${String(j+1).padStart(3,'0')}`,
      '摘要':e['證據摘要'], '既有解讀':e['如何解讀'], '可核對欄位':e['可核對欄位'] ?? {},
      '原始定位':e['原始位置'], '原始檔案':null, '定位狀態':'待核對',
    })), '外部證據':[], '外部查核狀態':'本次輸入未提供查核紀錄（不表示沒有文獻）',
    'OBER模組紀錄':null,
  }));
  const findCandidate = name => candidates.find(x=>key(x['病原'])===key(name));
  const warnings=[];
  for (const {file,data} of casefits.entries) {
    let c=findCandidate(data.candidate_organism);
    if (!c) {
      c={'候選ID':`C${String(candidates.length+1).padStart(3,'0')}`,'病原':data.candidate_organism,'候選來源':'OBER review 候選（原上游清單未列）','原上游處置':'未提供','最終處置':'未選入','最終排序':null,'為什麼考慮':[`由 ${data.review_tier} 送文獻判讀。`],'為什麼選入或排除':['未在本批次凍結選菌清單；未提供與該決策階段一致的排除理由，不臆測。'],'理由來源':'候選來源與凍結清單狀態；缺少決策理由','關鍵判斷資訊':{},'限制':[],'病例證據':[],'外部證據':[],'OBER模組紀錄':null};
      candidates.push(c);
      warnings.push(`${c['病原']} 來自 OBER review；未用同屬、縮寫或推測同義詞自動合併。`);
    }
    const ragCandidate=casefits.rag?.data.candidates.find(x=>key(x.organism_name)===key(data.candidate_organism) || key(x.canonical_organism_name)===key(data.canonical_organism));
    c['OBER模組紀錄']=casefits.rag?.sameInput ? ragCandidate?.modules ?? null : null;
    if(casefits.rag && !casefits.rag.sameInput) warnings.push('A1 對本版病例判讀，但重用較早 RAG 檢索的文獻；只串接其文獻資料，不混入旧版病人模組結論。');
    c['A1彙總']=data.A1_aggregate;
    c['外部查核狀態']='已有保存的 A1 判讀；與本病人致病性的判定不同';
    c['外部證據適用範圍']='歷史 OBER 文獻查核，作為查核附件；不宣稱是當時上游 picked 決策或簡易理由使用的證據。';
    const articles=ragCandidate?.literature_evidence?.retrieval?.articles ?? [];
    c['外部證據']=data.A1_judgments.map((j,index)=>{
      const article=articles.find(a=>String(a.pmid)===String(j.pmid));
      const normalize=s=>String(s??'').replace(/\s+/g,' ').trim();
      const spanMatch=Boolean(j.support_span && article?.abstract && normalize(article.abstract).includes(normalize(j.support_span)));
      const counterMatch=Boolean(j.counter_span && article?.abstract && normalize(article.abstract).includes(normalize(j.counter_span)));
      return {'證據ID':`${c['候選ID']}-L${String(index+1).padStart(3,'0')}`,'來源類型':'PubMed 文獻／病例報告','PMID':String(j.pmid),'標題':article?.title ?? '本地判讀紀錄未附標題','URL':j.article_url,'年份':article?.year??null,'判讀結論':verdictZh(j.verdict),'判讀理由原文':j.rationale,'支持原文':j.support_span??'','反向原文':j.counter_span??'','上游支持原文驗證':j.support_span_valid??null,'本次支持原文逐字核對':spanMatch,'本次反向原文逐字核對':counterMatch,'核對範圍':'只對已保存摘要核對文字；未重新連網驗證來源或醫學有效性','狀態':j.status,'判讀快照':`evidence/casefit/${path.basename(file)}`,'快照定位':`/A1_judgments/${index}`,'摘要快照':casefits.rag?'evidence/frozen_rag.json':null};
    });
  }
  if (decision) {
    if (!Array.isArray(decision.candidates)) throw Error('Final decision requires all candidate decisions/reasons');
    for (const item of arr(decision.candidates)) {
      const c=findCandidate(item.organism);
      if (!c) throw Error(`Final decision contains unknown candidate: ${item.organism}`);
      c['最終處置']=expected.some(x=>x.organism===item.organism)?'選入':'未選入';
      c['最終排序']=expected.find(x=>x.organism===item.organism)?.rank??null;
      c['為什麼選入或排除']=arr(item.reasons);
      c['理由來源']='明確提供的工作流最終決策紀錄';
    }
    for (const item of expected) if (!findCandidate(item.organism)) throw Error(`Final selection missing candidate evidence: ${item.organism}`);
    if (candidates.some(c=>!decision.candidates.some(x=>x.organism===c['病原']))) throw Error('Final decisions do not cover all candidates');
  }
  const llmCandidates=[...chosen,...sourceCandidates.filter(x=>x['最終處置']!=='選入').slice(0,5)];
  const llmMap={};
  llmCandidates.forEach((c,i)=>arr(c['病例內支持或反對證據']).forEach((e,j)=>{
    const cId=findCandidate(c['病原'])['候選ID'];
    llmMap[`C${String(i+1).padStart(2,'0')}-E${String(j+1).padStart(2,'0')}`]=`${cId}-E${String(j+1).padStart(3,'0')}`;
  }));
  const allModelRefs=[...arr(gen.patient_level_evidence_ids),...gen.selected_pathogen_explanations.flatMap(x=>[...arr(x.evidence_ids),...arr(x.limitations).flatMap(l=>arr(l.evidence_ids))])];
  if (allModelRefs.some(e=>!llmMap[e])) throw Error('Unknown reason-generator evidence ID');
  for (const x of gen.selected_pathogen_explanations) {
    const c=findCandidate(x.organism);
    if (!c) throw Error(`Rationale candidate missing: ${x.organism}`);
    const own=new Set(c['病例證據'].map(e=>e['證據ID']));
    if ([...arr(x.evidence_ids),...arr(x.limitations).flatMap(l=>arr(l.evidence_ids))].some(e=>!own.has(llmMap[e]))) throw Error('Cross-candidate rationale citation');
    c['理由產生器說明']={
      '原因':x.concise_reason,
      '關鍵證據摘要':x.key_evidence_summary,
      '限制':arr(x.limitations).map(l=>({...l,'對應證據ID':arr(l.evidence_ids).map(e=>llmMap[e])})),
      '需要醫師複核':Boolean(x.requires_clinician_review),
      '待醫師確認':x.clinician_review_question,
      '原始引用ID':arr(x.evidence_ids),
      '對應證據ID':arr(x.evidence_ids).map(e=>llmMap[e]),
    };
  }
  const missing=candidates.filter(c=>c['最終處置']==='選入'&&!c['病例證據'].length).map(c=>c['病原']);
  if(missing.length) warnings.push(`缺少可定位病例證據：${missing.join('、')}。不代表醫院沒有原始資料。`);
  if(!result.test_identifiers.length) warnings.push('檢驗編號未提供。');
  if(result.hospital.includes('待補')) warnings.push('兩院的逐病人院別尚待補。');
  if(stage!=='workflow_final') warnings.unshift('這是上游 picked 交付示範，不是已合併 OBER 的最終結果。');
  const patientData={...(chain['步驟1_病人資料摘要']??{})};
  // The original directory is needed by the exporter to archive evidence, but
  // it must not leak into the user-facing JSON report.
  delete patientData['原始病例資料夾'];
  const decisionLabel=stage==='workflow_final'&&decision?.status==='computational_final_pending_rationale_and_clinical_review'
    ? 'OBER R5 計算最終結果（待臨床複核）' : stageLabel(stage);
  return {
    '格式版本':'ober.patient_record.v1','資料組':cohort,'醫院':result.hospital,'病人編號':id,'檢驗編號':result.test_identifiers,
    '決策階段':decisionLabel,'決策階段代碼':stage,'決策來源':decision?.source ?? technical.selection_source,'批次說明':stageNote,
    '決策驗證狀態':decision?.validation ?? null,
    '病人資料':patientData,'候選診斷':candidates,
    '最終鑑別診斷排序':expected.sort((a,b)=>a.rank-b.rank).map(x=>({'排序':x.rank,'病原':x.organism,'原因':gen.selected_pathogen_explanations.find(y=>y.organism===x.organism).concise_reason,'原因來源':'理由產生器（既有 LLM 回覆）'})),
    '病例總結':gen.patient_result_summary,'病例總結證據ID':arr(gen.patient_level_evidence_ids).map(e=>llmMap[e]),
    '無入選病原':gen.no_selected_pathogen,'資料缺口與提醒':warnings,'來源清單':[],
    '理由產生器':{'模型':result.model_returned,'Response ID':result.response_id,'提示詞SHA256':result.prompt_sha256,'輸入SHA256':result.input_sha256,'來源快照':'evidence/llm_rationale.json','说明':'僅入選菌／無選菌總結是 LLM 生成；未選候選理由沿用既有決策紀錄。'},
    '查核狀態':'待人工臨床複核；檔案可追溯不等於醫學正確率已驗證',
  };
}

export function simpleView(record) {
  return {'格式版本':'ober.simple.v1','醫院':record['醫院'],'病人編號':record['病人編號'],'檢驗編號':record['檢驗編號'],'決策階段':record['決策階段'],'選擇菌種與原因':record['最終鑑別診斷排序'],'無入選病原':record['無入選病原'],'無入選時說明':record['無入選病原']?record['病例總結']:null,'原因來源':{'模型':record['理由產生器']['模型'],'說明':'LLM 僅將凍結選擇整理成短理由，不改變病原選擇；Response ID 與雜湊請見可追溯版。'},'提醒':record['資料缺口與提醒'],'詳細查核':'traceable.markdown'};
}

export function simpleMarkdown(simple) {
  const lines=[`# P${simple['病人編號']}｜簡易版`,'',`醫院：${simple['醫院']}`,`病人編號：P${simple['病人編號']}`,`檢驗編號：${simple['檢驗編號'].join('、')||'未提供'}`,`決策階段：${simple['決策階段']}`,''];
  if(simple['無入選病原']) lines.push('## 無入選病原','',simple['無入選時說明'],'');
  for(const x of simple['選擇菌種與原因']) lines.push(`## ${x['排序']}. ${x['病原']}`,'',x['原因'],'');
  lines.push('理由產生器只解釋已凍結的病原選擇，不重新診斷或改變排序。','');
  if(simple['提醒'].length) lines.push('## 資料提醒','',...simple['提醒'].map(x=>'- '+x),'');
  lines.push('[查看可追溯版](traceable.markdown)');
  return lines.join('\n')+'\n';
}

export function traceMarkdown(r) {
  const lines=[`# P${r['病人編號']}｜可追溯性報告`,'',`醫院：${r['醫院']}；檢驗編號：${r['檢驗編號'].join('、')||'未提供'}`,'',`決策階段：${r['決策階段']}`,'','病人資料 → 候選診斷 → 考慮／排除理由 → 病例與文獻證據 → 來源核對 → 最終排序','','這是可檢查的決策與證據紀錄，不是模型內部思考逐字稿。','','## 1. 病人資料',''];
  for(const [k,v] of Object.entries(r['病人資料'])) if(k!=='原始病例資料夾') lines.push(`- ${k}：${str(v)}`);
  lines.push('','## 2. 候選診斷','',...r['候選診斷'].map(c=>`- ${c['候選ID']} ${c['病原']}｜${c['最終處置']}｜${c['候選來源']}`),'','## 3. 為什麼考慮／選入／排除每個候選','');
  for(const c of r['候選診斷']) {
    lines.push(`### ${c['候選ID']} ${c['病原']}（${c['最終處置']}）`,'',`理由來源：${c['理由來源']}`,'','為什麼考慮：','',...(c['為什麼考慮'].length?c['為什麼考慮']:['原流程未提供；不臆測。']).map(x=>'- '+x),'','為什麼選入／排除：','',...(c['為什麼選入或排除'].length?c['為什麼選入或排除']:['原流程未提供；不臆測。']).map(x=>'- '+x),'');
    if(c['理由產生器說明']) {
      const llm=c['理由產生器說明'];
      lines.push('理由產生器說明：','',llm['原因'],'',`關鍵證據摘要：${llm['關鍵證據摘要']||'未提供'}`,`引用：${llm['對應證據ID'].join('、')||'無可定位證據；見資料提醒'}`,'');
      if(llm['限制'].length) {
        lines.push('理由產生器指出的限制：','');
        for(const item of llm['限制']) lines.push(`- ${item.text}${item['對應證據ID'].length?`（證據：${item['對應證據ID'].join('、')}）`:''}`);
        lines.push('');
      }
      lines.push(`需要醫師複核：${llm['需要醫師複核']?'是':'否'}`);
      if(llm['待醫師確認']) lines.push(`待醫師確認：${llm['待醫師確認']}`);
      lines.push('');
    }
    if(c['限制'].length) lines.push('既有資料限制：','',...c['限制'].map(x=>'- '+x),'');
    lines.push(`病例證據：${c['病例證據'].map(x=>x['證據ID']).join('、')||'未提供'}`,`外部查核：${c['外部查核狀態']}`,'');
  }
  lines.push('## 4. 支持／反對判斷的病例、文獻與來源','','外部文獻表示曾有相似病例或相關證據，不等於此病人已確定致病。歷史 OBER 查核未必是本次簡易理由當時使用的依據。','');
  for(const c of r['候選診斷']) {
    lines.push(`### ${c['候選ID']} ${c['病原']}`,'');
    for(const e of c['病例證據']) lines.push(`#### ${e['證據ID']}｜病例內證據`,'',e['摘要']??'未提供',`可核對欄位：${str(e['可核對欄位'])}`,`原始定位：${e['原始定位']}`,e['原始檔案']?`[開啟來源檔案](<${e['原始檔案']}>)`:'原始檔案未找到',`定位狀態：${e['定位狀態']}`,'');
    if(!c['外部證據'].length) lines.push(`外部文獻：${c['外部查核狀態']}`,'');
    for(const e of c['外部證據']) {
      lines.push(`#### ${e['證據ID']}｜${e['標題']}`,'',`PMID：${e['PMID']}｜${e['判讀結論']}｜狀態：${e['狀態']}`,`[PubMed 原頁](<${e['URL']}>)`,`既有 A1 判讀理由（保留原文）：${e['判讀理由原文']}`,'');
      if(e['支持原文']) lines.push(`支持原文：${e['支持原文']}`,`保存摘要逐字核對：${e['本次支持原文逐字核對']?'符合（僅文字核對）':'未確認，勿視為有效引文'}`,'');
      if(e['反向原文']) lines.push(`反向原文：${e['反向原文']}`,'');
      lines.push(`[查看原始判讀紀錄](<${e['判讀快照']}>)；JSON 定位：${e['快照定位']}`,'');
    }
  }
  lines.push('## 5. 可回查的來源與限制','',...r['來源清單'].map(s=>`- [${s['用途']}](<${s['封存檔案']}>)｜SHA256：${s['SHA256']}`),'',...r['資料缺口與提醒'].map(x=>'- '+x),'','## 6. 最終鑑別診斷排序','',`排序來源階段：${r['決策階段']}；不替原流程重排。`,'');
  if(r['無入選病原']) lines.push('無入選病原。',r['病例總結'],'');
  for(const x of r['最終鑑別診斷排序']) lines.push(`${x['排序']}. ${x['病原']}：${x['原因']}`,'');
  lines.push(r['查核狀態'],'','[返回簡易版](simple.markdown)');
  return lines.join('\n')+'\n';
}

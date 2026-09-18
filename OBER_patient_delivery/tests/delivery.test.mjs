import test from 'node:test';
import assert from 'node:assert/strict';
import {assertSelection,verifyPatient,makeRecord,simpleView,simpleMarkdown,traceMarkdown} from '../delivery_model.mjs';

function fixture(selected=true) {
  const organism='Candida albicans';
  const explanation={organism,rank:1,concise_reason:'測試理由，不是醫學結論。',key_evidence_summary:'測試證據。',evidence_ids:['C01-E01'],limitations:[{text:'測試限制。',evidence_ids:['C01-E01']}],requires_clinician_review:true,clinician_review_question:'請核對。'};
  return {cohort:'test',stage:'upstream_picked',stageNote:'test',decision:null,casefits:{entries:[],rag:null},
    result:{patient_id:'1',hospital:'測試院',test_identifiers:['001A'],generated_rationale:{patient_id:'1',no_selected_pathogen:!selected,patient_result_summary:'病例總結。',patient_level_evidence_ids:['C01-E01'],selected_pathogen_explanations:selected?[explanation]:[]},model_returned:'frozen-test',response_id:'test',prompt_sha256:'x',input_sha256:'y'},
    chain:{'病例識別':{'病人編號':'1'},'步驟1_病人資料摘要':{'基本資料':{'年齡':'60'},'原始病例資料夾':'C:/private/source'},'步驟2至5_逐候選病原判讀':[{'病原':organism,'最終處置':selected?'選入':'未選入','最終排序':selected?1:null,'為什麼會被考慮':['測試考慮理由'],'為什麼選入或排除':['測試決策理由'],'重要限制':[],'病例內支持或反對證據':[{'證據摘要':'合成測試紀錄','如何解讀':'不作醫學推論','原始位置':'NGS_patient_1_culture.json#0','可核對欄位':{'檢體':'BALF'}}]}]},
    technical:{patient_id:'1',selection_source:'test-picked'},selection:{patient_id:'1',selected_pathogens:selected?[{canonical_name:organism,rank:1}]:[]}};
}
test('species and rank cannot silently change',()=>{
  assert.throws(()=>assertSelection([{organism:'Candida albicans',rank:1}],[{organism:'Candida tropicalis',rank:1}]),/differ/);
  assert.throws(()=>assertSelection([{organism:'A',rank:1}],[{organism:'A',rank:2}]),/differ/);
});
test('duplicate predictions rejected',()=>assert.throws(()=>assertSelection([],[{organism:'A',rank:1},{organism:'A',rank:2}]),/duplicate/));
test('patient identity mismatch rejected',()=>assert.throws(()=>verifyPatient('1',{patient_id:'2'}),/Cross-patient/));
test('two views preserve the identical decision and LLM text',()=>{
  const r=makeRecord(fixture()),s=simpleView(r);
  assert.deepEqual(s['選擇菌種與原因'],r['最終鑑別診斷排序']);
  assert.match(simpleMarkdown(s),/測試理由，不是醫學結論/);
  assert.match(traceMarkdown(r),/測試理由，不是醫學結論/);
  assert.match(traceMarkdown(r),/關鍵證據摘要：測試證據/);
  assert.match(traceMarkdown(r),/測試限制。/);
  assert.match(traceMarkdown(r),/需要醫師複核：是/);
  assert.match(traceMarkdown(r),/待醫師確認：請核對。/);
  assert.equal(r['候選診斷'][0]['理由產生器說明']['對應證據ID'][0],'C001-E001');
  assert.equal(r['候選診斷'][0]['理由產生器說明']['限制'][0]['對應證據ID'][0],'C001-E001');
  assert.equal(s['原因來源']['Response ID'],undefined);
  assert.equal(r['病人資料']['原始病例資料夾'],undefined);
});
test('no-pick patient is retained in both views',()=>{
  const r=makeRecord(fixture(false));
  assert.equal(r['無入選病原'],true);assert.equal(r['最終鑑別診斷排序'].length,0);
  assert.match(simpleMarkdown(simpleView(r)),/無入選病原/);assert.match(traceMarkdown(r),/無入選病原/);
});
test('workflow-final claim requires an explicit final decision',()=>{
  const f=fixture();f.stage='workflow_final';assert.throws(()=>makeRecord(f),/requires explicit/);
});
test('unknown references rejected',()=>{
  const f=fixture();f.result.generated_rationale.selected_pathogen_explanations[0].evidence_ids=['C99-E01'];assert.throws(()=>makeRecord(f),/Unknown/);
});
test('explicit final decision works and covers candidates',()=>{
  const f=fixture();f.stage='workflow_final';f.decision={patient_id:'1',stage:'workflow_final',source:'synthetic_final_test',selected:[{organism:'Candida albicans',rank:1}],candidates:[{organism:'Candida albicans',reasons:['凍結最終規則']}]};
  const r=makeRecord(f);assert.equal(r['決策階段代碼'],'workflow_final');assert.equal(r['候選診斷'][0]['為什麼選入或排除'][0],'凍結最終規則');
  f.decision.candidates=[];assert.throws(()=>makeRecord(f),/cover all/);
});

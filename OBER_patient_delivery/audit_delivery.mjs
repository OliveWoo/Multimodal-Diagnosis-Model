import fs from 'node:fs/promises';
import path from 'node:path';
import {readJson, sha256} from './delivery_model.mjs';

function argsOf(argv) {
  if(argv.length!==2||argv[0]!=='--output') throw Error('Usage: node audit_delivery.mjs --output completed_batch');
  return {output:path.resolve(argv[1])};
}

const same=(a,b)=>JSON.stringify(a)===JSON.stringify(b);

async function main() {
  const {output}=argsOf(process.argv.slice(2));
  const manifest=await readJson(path.join(output,'manifest.json'));
  if(manifest.status!=='complete') throw Error('Manifest is not complete');
  const finalListJson=await fs.readFile(path.join(output,manifest.final_selection_list.json),'utf8');
  const finalListMarkdown=await fs.readFile(path.join(output,manifest.final_selection_list.markdown),'utf8');
  if(sha256(finalListJson)!==manifest.final_selection_list.json_sha256||sha256(finalListMarkdown)!==manifest.final_selection_list.markdown_sha256) throw Error('Final-selection list hash mismatch');
  const finalList=JSON.parse(finalListJson);
  const totals={patients:0,selected:0,candidates:0,literature:0,rationale_blocks:0,limitations:0,review_questions:0,rescued:0};
  const rescueRows=[];
  for(const row of manifest.patients) {
    const patientDir=path.join(output,row.patient_dir);
    const files={};
    for(const name of ['simple.json','simple.markdown','traceable.json','traceable.markdown']) {
      const bytes=await fs.readFile(path.join(patientDir,name));
      if(sha256(bytes)!==row.files[name]) throw Error(`Primary file hash mismatch: ${row.patient_dir}/${name}`);
      files[name]=bytes.toString('utf8');
    }
    if(/[A-Za-z]:\\Users\\/i.test(files['simple.json']+files['traceable.json'])||files['traceable.json'].includes('"原始路徑"')) {
      throw Error(`User-facing JSON exposes an absolute local path: ${row.patient_dir}`);
    }
    const simple=JSON.parse(files['simple.json']), trace=JSON.parse(files['traceable.json']);
    if(trace['決策階段代碼']!=='workflow_final'||!trace['決策階段'].includes('待臨床複核')) throw Error(`Wrong final-stage label: ${row.patient_dir}`);
    if(!same(simple['選擇菌種與原因'],trace['最終鑑別診斷排序'])) throw Error(`Simple/trace decision mismatch: ${row.patient_dir}`);
    if(simple['原因來源']['Response ID']||simple['原因來源']['提示詞SHA256']) throw Error(`Simple view contains technical IDs: ${row.patient_dir}`);
    for(const source of trace['來源清單']) {
      if(path.isAbsolute(source['封存檔案'])||source['原始路徑']) throw Error(`Unsafe source entry: ${row.patient_dir}`);
    }
    const md=files['traceable.markdown'];
    for(const candidate of trace['候選診斷']) {
      const llm=candidate['理由產生器說明'];
      if(!llm) continue;
      totals.rationale_blocks++;
      for(const required of [llm['原因'],llm['關鍵證據摘要']]) if(required&&!md.includes(required)) throw Error(`Rationale text omitted: ${row.patient_dir}/${candidate['病原']}`);
      for(const limitation of llm['限制']) {
        totals.limitations++;
        if(!md.includes(limitation.text)) throw Error(`Limitation omitted: ${row.patient_dir}/${candidate['病原']}`);
      }
      if(llm['待醫師確認']) {
        totals.review_questions++;
        if(!md.includes(llm['待醫師確認'])) throw Error(`Review question omitted: ${row.patient_dir}/${candidate['病原']}`);
      }
      if(!md.includes(`需要醫師複核：${llm['需要醫師複核']?'是':'否'}`)) throw Error(`Review flag omitted: ${row.patient_dir}/${candidate['病原']}`);
    }
    const decision=await readJson(path.join(patientDir,'evidence','final_decision.json'));
    if(decision.validation?.gold_answer_used_in_this_run!==false) throw Error(`Gold-use validation flag missing: ${row.patient_dir}`);
    for(const selected of decision.selected) if(selected.selection_origin==='OBER_R5_rescue') {
      totals.rescued++;
      rescueRows.push(`${row.cohort}/P${row.patient_id} ${selected.organism}`);
    }
    totals.patients++;
    totals.selected+=trace['最終鑑別診斷排序'].length;
    totals.candidates+=trace['候選診斷'].length;
    totals.literature+=trace['候選診斷'].reduce((sum,c)=>sum+c['外部證據'].length,0);
  }
  for(const workbook of ['simple_summary.xlsx','traceable_summary.xlsx']) {
    const scan=await fs.readFile(path.join(output,'_verification',workbook+'.formula_scan.json'),'utf8');
    if(!scan.includes('matched 0 entries')) throw Error(`Workbook formula scan is not clean: ${workbook}`);
  }
  for(const field of ['patients','selected','candidates','literature']) if(totals[field]!==manifest.summary[field]) throw Error(`Manifest reconciliation failed: ${field}`);
  if(finalList.summary.patients!==totals.patients||finalList.summary.selected!==totals.selected) throw Error('Final-selection list totals do not reconcile');
  if(totals.rationale_blocks!==totals.selected) throw Error('Not every selected pathogen has a rationale block');
  console.log(JSON.stringify({status:'verified',totals,rescues:rescueRows},null,2));
}

main().catch(error=>{console.error(error.stack);process.exitCode=1;});

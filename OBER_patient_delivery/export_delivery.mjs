import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {readJson, sha256, jsonText, arr, verifyPatient, loadCasefits, makeRecord, simpleView, simpleMarkdown, traceMarkdown} from './delivery_model.mjs';

const here=path.dirname(fileURLToPath(import.meta.url));
function argsOf(argv) {
  const args={};
  for(let i=0;i<argv.length;i+=2) {if(!argv[i].startsWith('--')||!argv[i+1]) throw Error('Usage: node export_delivery.mjs --config config.json --output new_folder'); args[argv[i].slice(2)]=argv[i+1];}
  if(!args.config||!args.output) throw Error('Required: --config and --output');
  return args;
}
export async function prepareBatch(configPath) {
  const config=await readJson(configPath), base=path.dirname(path.resolve(configPath));
  if(!['upstream_picked','workflow_final'].includes(config.decision_stage)) throw Error('Declare decision_stage: upstream_picked or workflow_final');
  const resolve=p=>path.resolve(base,p);
  const accepted=await readJson(resolve(config.accepted_rationales));
  if(accepted.status!=='structurally_validated') throw Error('Input rationale collection is not structurally validated');
  const batches=[], seen=new Set();
  for(const originalEntry of accepted.patients) {
    let entry=originalEntry;
    const override=config.rationale_overrides?.[`${entry.cohort}/${entry.result.patient_id}`];
    if(override) {
      const resultFile=resolve(override),replacement=await readJson(resultFile);
      verifyPatient(entry.result.patient_id,replacement);
      entry={...entry,result:replacement,output_file:resultFile};
    }
    const result=entry.result, cohort=entry.cohort, id=String(result.patient_id);
    if(!/^[A-Za-z0-9_-]+$/.test(cohort)||!/^\d+$/.test(id)) throw Error('Unsafe cohort/patient path');
    if(seen.has(`${cohort}/${id}`)) throw Error('Duplicate patient identity'); seen.add(`${cohort}/${id}`);
    const cc=config.cohorts[cohort]; if(!cc) throw Error(`Missing cohort config: ${cohort}`);
    const resultFile=path.resolve(entry.output_file);
    if(JSON.stringify(result)!==JSON.stringify(await readJson(resultFile))) throw Error('Accepted rationale differs from original file');
    const chainFile=path.resolve(result.input_source),chain=await readJson(chainFile);
    const technicalFile=path.resolve(chain['完整技術追溯']['檔案']),technical=await readJson(technicalFile);
    const selectedFile=path.join(resolve(cc.selected_root),`NGS_patient_${id}_selected_pathogens.json`),selection=await readJson(selectedFile);
    verifyPatient(id,result,chain,technical,selection);
    if(selection.source.sha256!==sha256(await fs.readFile(selection.source.file))) throw Error(`Changed original selection source: ${cohort}/P${id}`);
    const casefits=await loadCasefits(cc.casefit_root?resolve(cc.casefit_root):null,id,selection);
    let decision=null,decisionFile=null;
    if(cc.final_decision_root) {decisionFile=path.join(resolve(cc.final_decision_root),`P${id}_final_decision.json`);decision=await readJson(decisionFile);}
    const record=makeRecord({cohort,result,chain,technical,selection,casefits,decision,stage:config.decision_stage,stageNote:config.selection_stage_note??''});
    batches.push({record,chain,technical,selection,casefits,files:{resultFile,chainFile,technicalFile,selectedFile,decisionFile}});
  }
  return {config,batches};
}

async function snapshot(source,destination,relative,label,record) {
  const bytes=await fs.readFile(source);
  await fs.mkdir(path.dirname(destination),{recursive:true}); await fs.writeFile(destination,bytes);
  record['來源清單'].push({'用途':label,'原始來源檔名':path.basename(source),'封存檔案':relative,'SHA256':sha256(bytes)});
}

export async function writePatients(batch, output) {
  const summaries=[];
  for(const item of batch.batches) {
    const r=item.record, relative=`patients/${r['資料組']}/P${r['病人編號']}`,patientDir=path.join(output,relative);
    await fs.mkdir(patientDir,{recursive:true});
    for(const [field,filename,label] of [['resultFile','llm_rationale.json','理由產生器原始回覆'],['chainFile','input_chain.json','輸入病例與候選'],['technicalFile','technical_trace.json','上游完整判讀紀錄'],['selectedFile','selected_pathogens.json','凍結選菌與來源']]) {
      await snapshot(item.files[field],path.join(patientDir,'evidence',filename),`evidence/${filename}`,label,r);
    }
    if(item.files.decisionFile) await snapshot(item.files.decisionFile,path.join(patientDir,'evidence/final_decision.json'),'evidence/final_decision.json','工作流最終決策',r);
    for(const entry of item.casefits.entries) await snapshot(entry.file,path.join(patientDir,'evidence/casefit',path.basename(entry.file)),`evidence/casefit/${path.basename(entry.file)}`,'A1 文獻判讀快照',r);
    if(item.casefits.rag) await snapshot(item.casefits.rag.file,path.join(patientDir,'evidence/frozen_rag.json'),'evidence/frozen_rag.json','保存的檢索摘要與 RAG 模組',r);
    const rawRoot=path.resolve(item.chain['步驟1_病人資料摘要']['原始病例資料夾']);
    const rawFiles=new Map();
    for(const c of r['候選診斷']) for(const e of c['病例證據']) {
      const file=String(e['原始定位']??'').split('#')[0];
      if(!file || path.basename(file)!==file || !file.startsWith(`NGS_patient_${r['病人編號']}_`)) {e['定位狀態']='來源定位格式不完整／病人編號不符，未複製';continue;}
      if(!rawFiles.has(file)) {
        const raw=path.join(rawRoot,file);
        try {await snapshot(raw,path.join(patientDir,'evidence/raw',file),`evidence/raw/${file}`,'原始病例檢驗',r);rawFiles.set(file,true);}
        catch(err) {if(err.code!=='ENOENT') throw err;rawFiles.set(file,false);}
      }
      if(rawFiles.get(file)) {e['原始檔案']=`evidence/raw/${file}`;e['定位狀態']='原始檔案已封存；索引及菌種／數值對應仍需人工核對';}
      else {e['定位狀態']='本地原始檔案缺失；保留上游摘錄';r['資料缺口與提醒'].push(`未找到原始檔案：${file}`);}
    }
    r['資料缺口與提醒']=[...new Set(r['資料缺口與提醒'])];
    // Only four primary documents; evidence is kept in a separate audit subfolder.
    const simple=simpleView(r);
    const files={'simple.json':jsonText(simple),'simple.markdown':simpleMarkdown(simple),'traceable.json':jsonText(r),'traceable.markdown':traceMarkdown(r)};
    for(const [name,text] of Object.entries(files)) await fs.writeFile(path.join(patientDir,name),text,'utf8');
    const reread=await readJson(path.join(patientDir,'simple.json'));
    if(JSON.stringify(reread['選擇菌種與原因'])!==JSON.stringify(r['最終鑑別診斷排序'])) throw Error('Simple/trace output mismatch');
    if(r['候選診斷'].some(c=>!files['traceable.markdown'].includes(c['病原']))) throw Error('Candidate omitted from Markdown');
    if(r['最終鑑別診斷排序'].some(c=>!files['simple.markdown'].includes(c['原因']))) throw Error('Reason text omitted from Markdown');
    summaries.push({cohort:r['資料組'],patient_id:r['病人編號'],patient_dir:relative,selected:r['最終鑑別診斷排序'].length,candidates:r['候選診斷'].length,literature:r['候選診斷'].reduce((n,c)=>n+c['外部證據'].length,0),files:Object.fromEntries(Object.entries(files).map(([n,t])=>[n,sha256(t)]))});
  }
  return summaries;
}

async function main() {
  const args=argsOf(process.argv.slice(2)), output=path.resolve(args.output);
  try {if((await fs.readdir(output)).length) throw Error('Output must be new or empty. Existing results are never overwritten.');} catch(err) {if(err.code!=='ENOENT') throw err;}
  // All source, identity and selection checks happen before writing patient outputs.
  const batch=await prepareBatch(path.resolve(args.config));
  // Fail before output if the spreadsheet runtime cannot be loaded.
  const {buildWorkbooks}=await import('./workbooks.mjs');
  await fs.mkdir(output,{recursive:true});
  await fs.writeFile(path.join(output,'RUNNING.json'),jsonText({status:'in_progress',started_at:new Date().toISOString()}));
  const summaries=await writePatients(batch,output);
  const xlsx=await buildWorkbooks(batch.batches.map(x=>x.record),output);
  const records=batch.batches.map(x=>x.record);
  const finalList={
    schema_version:'ober.final_selection_list.v1',decision_stage:batch.config.decision_stage,
    decision_label:records[0]?.['決策階段']??'未提供',batch_note:batch.config.selection_stage_note??'',
    summary:{patients:records.length,selected:records.reduce((n,r)=>n+r['最終鑑別診斷排序'].length,0)},
    patients:records.map(r=>({cohort:r['資料組'],hospital:r['醫院'],patient_id:r['病人編號'],test_identifiers:r['檢驗編號'],selected:r['最終鑑別診斷排序']})),
  };
  const finalListJson=jsonText(finalList);
  const finalListMarkdown=['# OBER R5 最終選擇清單','',`決策階段：${finalList.decision_label}`,'',finalList.batch_note,'',`共 ${finalList.summary.patients} 位病人、${finalList.summary.selected} 筆病人×病原選擇。`,'','| 資料組 | 病人 | 醫院 | 檢驗編號 | 最終選擇 |','|---|---:|---|---|---|',...finalList.patients.map(p=>`| ${p.cohort} | P${p.patient_id} | ${p.hospital} | ${p.test_identifiers.join('、')||'未提供'} | ${p.selected.map(x=>`${x['排序']}. ${x['病原']}`).join('<br>')||'無入選病原'} |`),'','每個病原的理由與限制請見病人資料夾中的 `simple.markdown` 與 `traceable.markdown`。',''].join('\n');
  await fs.writeFile(path.join(output,'final_selection_list.json'),finalListJson,'utf8');
  await fs.writeFile(path.join(output,'final_selection_list.markdown'),finalListMarkdown,'utf8');
  const manifest={schema_version:'ober.delivery_manifest.v1',status:'complete',created_at:new Date().toISOString(),decision_stage:batch.config.decision_stage,config_path:path.resolve(args.config),config_sha256:sha256(await fs.readFile(args.config)),offline_export:true,patients:summaries,summary:{patients:summaries.length,selected:summaries.reduce((n,p)=>n+p.selected,0),candidates:summaries.reduce((n,p)=>n+p.candidates,0),literature:summaries.reduce((n,p)=>n+p.literature,0)},final_selection_list:{json:'final_selection_list.json',json_sha256:sha256(finalListJson),markdown:'final_selection_list.markdown',markdown_sha256:sha256(finalListMarkdown)},workbooks:xlsx};
  await fs.writeFile(path.join(output,'manifest.json'),jsonText(manifest));
  const index=['# 病人交付包','',`資料階段：${batch.config.decision_stage}。${batch.config.selection_stage_note??''}`,'',`本批次 ${summaries.length} 位；${manifest.summary.selected} 筆入選理由。`,'','- [最終選擇清單（Markdown）](final_selection_list.markdown)','- [最終選擇清單（JSON）](final_selection_list.json)','- [簡易版總表](simple_summary.xlsx)','- [可追溯版總表](traceable_summary.xlsx)','- [程式與使用說明](../../README.md)','','## 每位病人','',...summaries.map(p=>`- ${p.cohort} P${p.patient_id}：[簡易版](<${p.patient_dir}/simple.markdown>)｜[可追溯版](<${p.patient_dir}/traceable.markdown>)`),''];
  await fs.writeFile(path.join(output,'README.md'),index.join('\n'),'utf8');
  // Remove only the exact marker created by this run; no source files are deleted.
  await fs.unlink(path.join(output,'RUNNING.json'));
  console.log(JSON.stringify({status:'complete',output,summary:manifest.summary},null,2));
}
if(process.argv[1] && path.resolve(process.argv[1])===fileURLToPath(import.meta.url)) main().catch(err=>{console.error(err.stack);process.exitCode=1;});

import fs from 'node:fs/promises';
import path from 'node:path';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';
import {jsonText, str} from './delivery_model.mjs';

const colors={navy:'#17324D',teal:'#0F766E',pale:'#E8F3F1',ink:'#243746',light:'#F4F7FA',amber:'#FFF1CC'};
const patientKey=r=>`${r['資料組']}/P${r['病人編號']}`;
const rel=r=>`patients/${r['資料組']}/P${r['病人編號']}`;
const literal=v=>typeof v==='string'&&v.startsWith('=')?"'"+v:v;
const formulaString=s=>String(s).replace(/"/g,'""');
const hyperlink=(target,label)=>`=HYPERLINK("${formulaString(target)}","${formulaString(label)}")`;
const letters=n=>{let s='';for(let x=n;x>0;x=Math.floor((x-1)/26))s=String.fromCharCode(65+(x-1)%26)+s;return s;};
let tableSequence=0;

function dataSheet(wb,name,headers,rows,widths={}) {
  const s=wb.worksheets.add(name), count=Math.max(rows.length,1),cols=headers.length;
  s.showGridLines=false;
  s.getRangeByIndexes(0,0,count+1,cols).format={font:{name:'Microsoft JhengHei',size:11,color:colors.ink},verticalAlignment:'top',wrapText:true};
  s.getRangeByIndexes(0,0,1,cols).values=[headers];
  if(rows.length)s.getRangeByIndexes(1,0,rows.length,cols).values=rows.map(row=>row.map(literal));
  s.getRangeByIndexes(0,0,1,cols).format={fill:colors.navy,font:{bold:true,color:'#FFFFFF'},rowHeight:34,wrapText:true};
  for(let col=0;col<cols;col++)s.getRangeByIndexes(0,col,count+1,1).format.columnWidth=widths[col]??20;
  rows.forEach((row,i)=>{
    const estimated=Math.max(1,...row.map((v,c)=>String(v??'').split('\n').reduce((n,line)=>n+Math.max(1,Math.ceil([...line].reduce((sum,ch)=>sum+(ch.charCodeAt(0)>255?2:1),0)/(widths[c]??20))),0)));
    s.getRangeByIndexes(i+1,0,1,cols).format.rowHeight=Math.min(409,Math.max(35,estimated*15+12));
  });
  s.freezePanes.freezeRows(1);
  if(rows.length){const t=s.tables.add(`A1:${letters(cols)}${rows.length+1}`,true,`DeliveryData${++tableSequence}`);t.style='TableStyleMedium2';t.showFilterButton=true;}
  return s;
}

function formulaColumn(s,col,formulas) {if(formulas.length)s.getRangeByIndexes(1,col,formulas.length,1).formulas=formulas.map(x=>[x]);}

function instructions(wb,title,records,metricRefs) {
  const s=wb.worksheets.add('使用說明');s.showGridLines=false;
  s.getRange('A1:F1').merge();s.getRange('A1').values=[[title]];
  s.getRange('A1:F1').format={fill:colors.navy,font:{name:'Microsoft JhengHei',size:18,bold:true,color:'#FFFFFF'},rowHeight:42};
  s.getRange('A2:F19').format={font:{name:'Microsoft JhengHei',size:11,color:colors.ink},wrapText:true,verticalAlignment:'top'};
  s.getRange('A2:A19').format.columnWidth=25;s.getRange('B2:B19').format.columnWidth=34;s.getRange('C2:F19').format.columnWidth=16;
  const notes=[
    ['決策資料階段',records[0]?.['決策階段']??'未提供'],
    ['病人數',null],['入選病原筆數',null],
    ['本表用途','由同一份病人紀錄生成；總覽可篩選，明細可回查。'],
    ['計數單位','病例＝資料組＋病人編號；菌數是病人×病原，不是不同菌種數。'],
    ['臨床驗證狀態','待人工複核。引用存在或原文逐字相符，不等於診斷正確。'],
    ['欄位缺漏','未提供≠陰性；空的最終排序表示未選入，不自行加排序。'],
    ['原因來源','簡易理由為 LLM 對凍結選擇的說明；候選考慮／排除理由沿用明確工作流決策紀錄。'],
    ['文獻來源','只使用已保存 A1／PubMed 紀錄；未新增網路檢索。非所有候選皆有文獻。'],
    ['檔案連結','請保留 Excel 與 patients 資料夾的相對位置，整個批次一起搬移。'],
    ['顏色說明','深藍＝欄名；淡黃＝資料／來源缺口，非醫療風險分級。'],
    ['編輯注意','本表是報告快照；若更改選菌請回原工作流重跑，以保持 JSON／Markdown 一致。'],
  ];
  notes.forEach((row,i)=>{const r=i+3;s.getRange(`A${r}`).values=[[row[0]]];s.getRange(`B${r}:F${r}`).merge();if(row[1]!==null)s.getRange(`B${r}`).values=[[row[1]]];s.getRange(`A${r}:F${r}`).format.rowHeight=[3,4].includes(i)?48:40;});
  s.getRange('B4').formulas=[[metricRefs.patients]];s.getRange('B5').formulas=[[metricRefs.selected]];
  s.getRange('B4:F5').format={fill:colors.pale,font:{size:16,bold:true,color:colors.teal},numberFormat:'0'};
  s.getRange('B3:F3').format.fill=records.some(r=>r['決策階段代碼']!=='workflow_final')?colors.amber:colors.pale;
  return s;
}

async function verifyAndSave(wb,file,output,expected,sheetRanges) {
  const notes=wb.worksheets.getItem('使用說明');
  const vals=notes.getRange('B4:B5').values;
  if(Number(vals[0][0])!==expected.patients||Number(vals[1][0])!==expected.selected) throw Error(`Workbook count reconciliation failed: ${file}; ${JSON.stringify(vals)}`);
  const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A',options:{useRegex:true,maxResults:30},maxChars:3000});
  // The count checks above also force evaluation; keep the full compact scan for audit.
  await fs.mkdir(path.join(output,'_verification'),{recursive:true});
  await fs.writeFile(path.join(output,'_verification',file+'.formula_scan.json'),typeof errors.ndjson==='string'?errors.ndjson:jsonText(errors));
  for(const [sheetName,range] of Object.entries(sheetRanges)) {
    const preview=await wb.render({sheetName,range,scale:1,format:'png'});
    await fs.writeFile(path.join(output,'_verification',`${file}_${sheetName}.png`),new Uint8Array(await preview.arrayBuffer()));
  }
  const exported=await SpreadsheetFile.exportXlsx(wb);await exported.save(path.join(output,file));
  console.log(`Workbook verified: ${file} (${expected.patients} patients / ${expected.selected} selected)`);
  return {file,counts:expected,formula_scan:`_verification/${file}.formula_scan.json`};
}

export async function buildWorkbooks(records,output) {
  const expected={patients:records.length,selected:records.reduce((n,r)=>n+r['最終鑑別診斷排序'].length,0)};
  const simple=Workbook.create(),details=[];
  for(const r of records) {
    const choices=r['最終鑑別診斷排序'];
    if(!choices.length)details.push([patientKey(r),r['醫院'],`P${r['病人編號']}`,r['檢驗編號'].join('; ')||'未提供',null,'無入選病原',r['病例總結'],0,`${rel(r)}/simple.markdown`]);
    for(const x of choices)details.push([patientKey(r),r['醫院'],`P${r['病人編號']}`,r['檢驗編號'].join('; ')||'未提供',x['排序'],x['病原'],x['原因'],1,`${rel(r)}/simple.markdown`]);
  }
  const ds=dataSheet(simple,'選菌與理由',['病例鍵','醫院','病人編號','檢驗編號','排序','選擇菌種','原因（理由產生器）','入選筆數','個案檔案'],details,{0:23,1:26,2:12,3:23,4:8,5:34,6:92,7:12,8:38});
  const simpleRows=records.map(r=>[patientKey(r),r['醫院'],`P${r['病人編號']}`,r['檢驗編號'].join('; ')||'未提供',null,r['最終鑑別診斷排序'].map(x=>x['病原']).join('\n')||'無入選病原',r['病例總結'],`${rel(r)}/simple.markdown`,`${rel(r)}/simple.json`,r['決策階段'],r['資料缺口與提醒'].join('\n')]);
  const ps=dataSheet(simple,'病人總覽',['病例鍵','醫院','病人編號','檢驗編號','選菌數','選擇菌種','病例總結（理由產生器）','閱讀簡易版','JSON','決策階段','資料提醒'],simpleRows,{0:23,1:26,2:12,3:23,4:10,5:42,6:88,7:18,8:18,9:35,10:70});
  formulaColumn(ps,4,records.map((r,i)=>`=SUMIF('選菌與理由'!$A$2:$A$${details.length+1},A${i+2},'選菌與理由'!$H$2:$H$${details.length+1})`));
  formulaColumn(ps,7,records.map(r=>hyperlink(`${rel(r)}/simple.markdown`,'開啟簡易版')));formulaColumn(ps,8,records.map(r=>hyperlink(`${rel(r)}/simple.json`,'開啟 JSON')));
  if(records.length)ps.getRange(`K2:K${records.length+1}`).conditionalFormats.add('notContainsBlanks',{format:{fill:colors.amber}});
  instructions(simple,'簡易版｜病人選菌與理由',records,{patients:`=COUNTA('病人總覽'!A2:A${records.length+1})`,selected:`=SUM('選菌與理由'!H2:H${details.length+1})`});
  const first=await verifyAndSave(simple,'simple_summary.xlsx',output,expected,{'使用說明':'A1:F15','病人總覽':'A1:F5','選菌與理由':'C1:G4'});

  const trace=Workbook.create(),candidateRows=[],evidenceRows=[],literatureRows=[];
  for(const r of records)for(const c of r['候選診斷']) {
    candidateRows.push([patientKey(r),c['候選ID'],c['病原'],c['最終處置'],c['最終排序'],c['為什麼考慮'].join('\n'),c['為什麼選入或排除'].join('\n'),c['理由來源'],c['理由產生器說明']?.['原因']??'未另請理由產生器生成',c['理由產生器說明']?.['需要醫師複核']?'是':c['理由產生器說明']?'否':'未產生',c['理由產生器說明']?.['待醫師確認']??'',c['病例證據'].map(e=>e['證據ID']).join('; '),c['外部查核狀態'],c['最終處置']==='選入'?1:0,`${rel(r)}/traceable.markdown`]);
    for(const e of c['病例證據']) {const f=e['可核對欄位'];evidenceRows.push([patientKey(r),c['病原'],e['證據ID'],f['時間']??'未提供',f['檢體']??'未提供',f['結果']??'未提供',f['reads']??null,e['摘要']??'',e['原始定位'],e['原始檔案']?`${rel(r)}/${e['原始檔案']}`:'未找到',e['定位狀態']]);}
    for(const e of c['外部證據'])literatureRows.push([patientKey(r),c['病原'],e['證據ID'],e['PMID'],e['標題'],e['判讀結論'],e['判讀理由原文'],e['支持原文'],e['本次支持原文逐字核對']?'符合（僅文字）':'未確認',e['URL'],`${rel(r)}/${e['判讀快照']}`,e['快照定位'],e['狀態']]);
  }
  dataSheet(trace,'候選判讀',['病例鍵','候選ID','候選病原','最終處置','最終排序','為什麼考慮','為什麼選入或排除','理由來源','理由產生器說明','需要醫師複核','待醫師確認','病例證據ID','外部查核狀態','入選筆數','可追溯檔案'],candidateRows,{0:23,1:12,2:35,3:12,4:12,5:85,6:85,7:40,8:88,9:18,10:75,11:45,12:45,13:12,14:38});
  dataSheet(trace,'病例證據',['病例鍵','病原','證據ID','時間','檢體','結果','reads（如有）','證據摘要','原始定位','來源檔案','定位狀態'],evidenceRows,{0:23,1:36,2:20,3:23,4:22,5:20,6:16,7:90,8:65,9:75,10:65});
  dataSheet(trace,'文獻證據',['病例鍵','病原','證據ID','PMID','標題','A1判讀','既有判讀理由原文','支持原文','保存摘要逐字核對','來源URL','判讀快照','快照JSON定位','狀態'],literatureRows,{0:23,1:36,2:20,3:16,4:100,5:18,6:105,7:105,8:23,9:56,10:65,11:30,12:18});
  const traceRows=records.map(r=>[patientKey(r),r['醫院'],`P${r['病人編號']}`,r['檢驗編號'].join('; ')||'未提供',null,null,null,null,r['最終鑑別診斷排序'].map(x=>`${x['排序']}. ${x['病原']}`).join('\n')||'無入選病原',`${rel(r)}/traceable.markdown`,`${rel(r)}/traceable.json`,r['資料缺口與提醒'].join('\n'),r['決策階段']]);
  const ts=dataSheet(trace,'病人總覽',['病例鍵','醫院','病人編號','檢驗編號','候選數','選菌數','病例證據筆數','文獻判讀筆數','最終排序','閱讀可追溯版','JSON','資料提醒','決策階段'],traceRows,{0:23,1:26,2:12,3:23,4:12,5:10,6:16,7:16,8:45,9:22,10:18,11:85,12:40});
  const end=n=>Math.max(n+1,2);
  formulaColumn(ts,4,records.map((r,i)=>`=COUNTIF('候選判讀'!$A$2:$A$${end(candidateRows.length)},A${i+2})`));
  formulaColumn(ts,5,records.map((r,i)=>`=SUMIF('候選判讀'!$A$2:$A$${end(candidateRows.length)},A${i+2},'候選判讀'!$N$2:$N$${end(candidateRows.length)})`));
  formulaColumn(ts,6,records.map((r,i)=>`=COUNTIF('病例證據'!$A$2:$A$${end(evidenceRows.length)},A${i+2})`));
  formulaColumn(ts,7,records.map((r,i)=>`=COUNTIF('文獻證據'!$A$2:$A$${end(literatureRows.length)},A${i+2})`));
  formulaColumn(ts,9,records.map(r=>hyperlink(`${rel(r)}/traceable.markdown`,'開啟可追溯版')));formulaColumn(ts,10,records.map(r=>hyperlink(`${rel(r)}/traceable.json`,'開啟 JSON')));
  if(records.length)ts.getRange(`L2:L${records.length+1}`).conditionalFormats.add('notContainsBlanks',{format:{fill:colors.amber}});
  instructions(trace,'可追溯版｜候選、判讀與證據',records,{patients:`=COUNTA('病人總覽'!A2:A${records.length+1})`,selected:`=SUM('候選判讀'!N2:N${candidateRows.length+1})`});
  const second=await verifyAndSave(trace,'traceable_summary.xlsx',output,expected,{'使用說明':'A1:F15','病人總覽':'A1:H5','候選判讀':'B1:K3','病例證據':'A1:G5','文獻證據':'B1:F3'});
  return [first,second];
}

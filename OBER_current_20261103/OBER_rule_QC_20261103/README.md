# Rule-first evidence QC

QC是文獻語義待查核，不是mNGS檢體品質。研究可能感染源，不等於臨床確診。固定58篇公開摘要、5個可能受文獻閘門影響的候選；49例全部計分，原97筆正式結果不刪。主要規則與限制見PROTOCOL_zh.md。

在workspace根目錄執行：

```powershell
$oberPython = 'C:/Users/chung/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $oberPython -X utf8 -m unittest discover -s OBER_rule_QC_20261103 -p test_*.py
& $oberPython -X utf8 OBER_rule_QC_20261103/qc_runner.py prepare
& $oberPython -X utf8 OBER_rule_QC_20261103/qc_runner.py smoke
& $oberPython -X utf8 OBER_rule_QC_20261103/qc_runner.py pipeline
& $oberPython -X utf8 OBER_rule_QC_20261103/qc_eval.py
```

成功API結果直接命中快取，不重抽。程式、輸入、規則、schema均先凍結，更動須另建版本。沿用既有環境金鑰、官方OpenAI endpoint、gpt-5.6-luna medium，沒有存金鑰。只對QC文章固定兩次判讀；若新開啟候選才依原A1補缺少的病例匹配，不傳醫院答案。引用僅用實際句ID，展示節錄由程式擷取。格式／來源錯誤最多修補一次，不按臨床支持或F1重試；拒答、不完整、技術失敗不是臨床陰性。

`runs/qc_v1/`包含：

- articles/<key>/source.json、routing.json：原文、出版類型、規則／QC原因。
- articles/<key>/r0|r1/：QC實際請求／回應。
- article_decisions_without_answers.json：規則／單次QC／兩次一致／direct-only對照，保留未解QC。
- predictions_without_answers.json：所有Context候選與k曲線，先凍結再評分。
- candidate_distribution.csv、article_QC_distribution.csv、QC_repeat.csv、metrics.csv：分布、待查核、變動及P/R/F1。
- patients/<case>/REPORT_zh.md與research_selection.json：每位病人的研究清單、待查核事項與證據定位。
- evaluation.json、api_usage.json：完整評分與實際新用量／舊回應重用分開。

主設定QC_CONSENSUS_F0_match_k2。疾病相關肺部併發症不偽裝成直接病因；direct-only是預定診斷對照，避免把定義變更的效果誤稱成QC減少隨機性。遇分歧保持待查，不排除病人分母。此為固定病例判斷的開發replay，不是完整新版OBER部署效能或獨立test。

兩次同模型也可能同錯。快取可重現既存結果，不保證重新呼叫模型完全相同。來源與schema通過不等於醫學語義正確。API格式依[OpenAI Docs](https://developers.openai.com/api/docs/guides/structured-outputs?api-mode=responses)。

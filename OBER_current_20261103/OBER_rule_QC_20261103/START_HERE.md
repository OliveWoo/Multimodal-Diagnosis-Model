# 規則＋文獻QC：已完成版本

先看RESULTS_zh.md（數據）及FINDINGS_zh.md（解讀）；九篇未解來源見QC_REVIEW_zh.md。QC指文獻語義待查核，不是定序品質。

```powershell
$oberPython = 'C:/Users/chung/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $oberPython -X utf8 -m unittest discover -s OBER_rule_QC_20261103 -p test_*.py
& $oberPython -X utf8 OBER_rule_QC_20261103/qc_runner.py pipeline
& $oberPython -X utf8 OBER_rule_QC_20261103/qc_eval.py
& $oberPython -X utf8 OBER_rule_QC_20261103/qc_queue.py
```

已有結果會命中快取，不新增API；若要改規則、輸入或模型，另建版本，不覆寫凍結紀錄。沿用官方API與既有環境金鑰，程式未存金鑰。此次114份新回應均為公開摘要QC、沒有格式重試；110份舊A1回應重用，沒有新增病例API呼叫。用量不是帳單。

`runs/qc_v1/`包含全部原文與分流、實際API回應、文章QC決策、49份病人報告、研究清單、評分與用量；`human_QC_queue.json`是可機讀的九篇待查佇列。舊版不覆寫。

這是49例固定病例判斷的開發replay，不是獨立test或完整新版OBER效能。F1恢復81.28%但沒有超過舊版；未證明全流程隨機性已下降。詳細程序見README.md與PROTOCOL_zh.md。

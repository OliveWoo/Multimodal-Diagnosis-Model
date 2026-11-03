# 最終執行入口：兩個fresh QC共識輪次

```powershell
$oberPython = 'C:/Users/chung/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $oberPython -X utf8 -m unittest discover -s OBER_QC_repeat_20261103 -p test_*.py
& $oberPython -X utf8 OBER_QC_repeat_20261103/repeat_entry.py pipeline
& $oberPython -X utf8 OBER_QC_repeat_20261103/repeat_entry.py evaluate
```

每輪57篇各判兩次，A/B共228個新QC呼叫。G/T與A1為明確固定控制，不是從raw data端到端重跑整套OBER。原正式97筆不改，49例仍全部計分。

兩輪輸出在runs/repeat_v1/round_A與round_B；共用A1控制在common_controls。比較結果是COMPARISON_zh.md，以及comparison.json、fresh_API_audit.json、single_QC_slot_changes.csv、consensus_article_changes.csv、candidate_changes.csv、k_curve_comparison.csv。

入口包含新API之前完成的報告欄位映射修正，見REPORT_MAPPING_zh.md；不改模型推論或選擇。所有原凍結程式與舊版都保留。

本次成功呼叫可斷點接續；重讀同一份本次快取不算額外fresh run。若要第三輪，另建新實驗；不可因結果不理想重抽、換k或覆寫。兩輪同分不等於文章分類無變動，也不是長期穩定或醫學正確性的證明。

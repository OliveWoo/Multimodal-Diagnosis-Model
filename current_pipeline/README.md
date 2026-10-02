# Current deterministic integrated pipeline（2026-10-01）

此目錄保存目前新版 research candidate 的可執行程式、規則與對應離線測試。它與 repository 根目錄的 legacy deterministic v20＋OBER／R5 流程並存，不會覆寫舊版。

## 主入口

```powershell
cd current_pipeline
python -m tools.run_kh_integrated_release_v2 `
  --output-root C:\private_runs\current_F `
  --inventory-root C:\private_data\per_test_inventory `
  --patient-root C:\private_data\patient_summaries `
  --phenotype-root C:\private_data\phenotype_evidence `
  --linkage-audit C:\private_data\phenotype_linkage_audit.json `
  --answers C:\private_data\answers.csv `
  --reference-root C:\private_runs\reference
```

決策階段會先完成並寫入 answer-blind freeze；`--answers` 只供 freeze 後的 post-hoc evaluation 使用。若要做前瞻運行，應將 decision run 與後續 evaluation 分成不同權限／不同批次。

## 執行階段

1. current hospital-source summary；
2. per-test DNA／RNA compact entry＋identity／taxonomy；
3. analytical deterministic scorer；
4. index-event clinical timeline；
5. candidate phenotype 與 event packets；
6. history v4C route；
7. clinical deterministic scorer；
8. hospital evidence rehydration；
9. promotion v8；
10. reporting v3 Route A＋互斥 fallback。

完整流程圖：

- [舊版 deterministic＋OBER／R5](../docs/KH_OLD_COMPLETE_PICKED_OBER_PIPELINE_20261001_zh.md)
- [新版 current F deterministic 流程](../docs/KH_NEW_COMPLETE_DETERMINISTIC_PIPELINE_20261001_zh.md)
- [菌名、taxid 與 taxonomy 詳細流程](../docs/KH_ORGANISM_TAXONOMY_DETAILED_FLOW_V2_6_20261001_zh.md)

## 目錄

- `tools/`：目前正式 runner 的完整 local import closure，加上 taxonomy review／approval 工具。
- `rules/`：production 與研究比較規則；不含 benchmark answer override、held-out answer override 或人工逐病例 queue。
- `tests/`：integrated release、scorer、history、promotion、reporting、taxonomy 與 evidence preservation 測試。

## 測試

```powershell
cd current_pipeline
python -m unittest discover -s tests -p "test_*.py"
```

本次公開前在乾淨 clone 中執行結果：**426 tests passed，5 skipped，0 failures**。部分 integration tests 會建立合成暫存檔；不需要真實病例資料。若測試明確尋找 private frozen artifact，應在公開環境標記為 unavailable，而不是把病例資料加入 repository。

`PUBLIC_SOURCE_MANIFEST.csv` 保存本目錄每個公開檔案的相對路徑、SHA-256 與大小，供教授或其他研究者核對實際檢視的版本。

## 資料治理界線

此公開目錄刻意排除：

- 病人 JSON、Excel、PDF、PPT 與原始院端資料；
- benchmark／held-out 答案與逐病例人工 override；
- `outputs/`、run artifacts、模型 prompts／raw responses；
- `.env`、API key、cache、logs 與本機絕對路徑設定。

本模型是研究工具，輸出是待臨床複核的候選與 audit evidence，不是自動臨床診斷。

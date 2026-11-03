# 核心結論／細節QC分離：保守離線小實驗

2026-11-03：使用者認為目前結果可接受，要求只稍微測試，不覆寫目前版本。

## 事先固定的規則

目前基準：`QC_CONSENSUS_F0_match_k2`，49例／87答案，P76.00%、R87.36%、F1 81.28%。保留此基準，不自動採用小實驗。

本實驗策略：`QC_SPLIT_NEGATIVE`。

1. 原本兩次完整共識成立的文章，保留原結論。
2. 原本未解QC的文章，只有在兩次判讀各自resolved，且以下核心五欄完全相同、`qualified`明確為false時，才保留「核心一致不符合支持條件」。
   - 病原範圍`scope`、人類致病情境`human_any`、支持資格`qualified`、支持型別`basis`、純非人類研究`nonhuman_only`。
3. 這些文章的`noncausal`／`secondary_other`若不一致或未知，標記獨立的`detail_QC_pending`，保存兩次原判讀，不硬選其中一次。
4. 正向支持有分歧、核心五欄有分歧、單次未解或內部矛盾，仍保留原本未解QC。**不允許靠拆分把文章升成正向支持。**
5. 同樣規則套用所有文章，不為病例／病原加例外。不改提示詞或聲稱修好了LLM的語義分類。

注意：不符合本版「成熟文獻支持閘門」，不等於否定醫學致病性，也不一定代表A1病例匹配必須為零。來源核心恢復後可能改變A1的衍生人類／病原欄位與match計數，所以仍需真的重播下游，不能推定選菌一定不變。

## 範圍與評估

- 使用`OBER_QC_repeat_20261103/runs/repeat_v1/round_A`和`round_B`已保存的228份回應。**沒有新的模型重跑、沒有新增API呼叫。**
- 同一批58個病原—文章配對；G/T、A1原始病例匹配、文獻與原97筆正式選擇不變。
- 核心／細節判讀與兩輪下游預測先凍結，再讀答案計分。
- 固定主設定F0、match、k2；k0–10僅敏感度，不挑最佳k。
- 比較核心未解QC、細節QC、兩輪核心狀態差異、match計數、最後病例—病原清單與P/R/F1。
- 這是已反覆觀察的49例開發資料之條件式重播，不是新的隨機性實驗或獨立test。保留細節QC，不把待查件數變少當作醫學準確率提高。

依[OpenAI Docs評估原則](https://developers.openai.com/api/docs/guides/evaluation-best-practices)，先固定評估條件並保留紀錄；本實驗不使用已退場／退場中的線上Evals服務，只在本地計分。

## 執行

從專案根目錄：

```powershell
$oberPython = 'C:/Users/chung/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $oberPython -X utf8 -m unittest discover -s OBER_QC_core_split_20261103 -p 'test_*.py'
& $oberPython -X utf8 OBER_QC_core_split_20261103/offline_replay.py predict
& $oberPython -X utf8 OBER_QC_core_split_20261103/offline_replay.py evaluate
```

程式沒有API執行入口，推論／dispatch函式也會被封鎖。結果另存`runs/offline_v1`；舊目錄僅讀。

輸出：`RESULTS_zh.md`、`runs/offline_v1/comparison.json`、`article_comparison.csv`、`metrics.csv`、`k_curve.csv`、`candidate_changes.csv`及各輪衍生預測／A1紀錄。每一輪的原始模型判讀仍留在原位置，可回查。

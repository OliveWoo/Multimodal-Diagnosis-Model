# 標準化 JSON 到最終交付的單一編排入口

根目錄 `run_full_pipeline.py` 是整張主流程圖的**編排器**。它不重寫菌種規則，也不在 main 裡自行決定候選；它只做三件事：

1. 驗證 cohort、病人 ID、設定檔與每一步預期 input。
2. 依固定順序呼叫現有程式，上一階段 output 直接成為下一階段 input。
3. 記錄每個 command、工作目錄、完成狀態與產物 SHA256。

編排順序如下：

```text
standardized patient JSON
  → run_deterministic_v20.py
  → build_missed_candidate_review_queue.py
  → review_missed_mngs_candidates.py                 [LLM]
  → merge_deterministic_max_with_missed_review.py
  → extract_selected_pathogens.py
  ├→ evidence run_pipeline.py
  │   → generate_llm_rationales.py                   [LLM]
  │   → audit_llm_rationale_batch.py
  └→ RAG_re batch                                    [PubMed/LLM]
      → run_casefit_batch.py                         [LLM; missing frozen evidence may fallback]
  → generate_r5_final_decisions.py                   [loads clinical rules/config]
  → prepare_r5_rationale_inputs.py
  → generate_llm_rationales.py                       [LLM; R5-aligned second pass]
  → assemble_r5_rationales.py
  → export_delivery.mjs
  → audit_delivery.mjs
```

主程式內只有兩種 handoff adapter：把 per-patient merged JSON 複製成下一套工具需要的 flat directory，以及生成 OBER delivery config。複製時會核對 SHA256；這兩步都不改菌種內容。

## 使用方式

先複製並修改 [`examples/full_pipeline_config.template.json`](../examples/full_pipeline_config.template.json)。真實病例路徑、cohort 名稱、hospital、dataset、mNGS source kind 與模型設定都放在 config；不要寫死進程式。

先只檢查計畫。這不建立 output、不連網、不呼叫模型：

```powershell
python -B run_full_pipeline.py `
  --config examples/full_pipeline_config.template.json `
  --output-root local_outputs/full_run_001 `
  --plan
```

確認列出的 input、病人與 31 個 stage（範例兩 cohort）後，再以相同指令移除 `--plan` 執行。正式執行需要：

- `OPENAI_API_KEY` 在環境中，或 config 的 `env_file` 指向只存於本機的檔案；
- PubMed／模型所需網路；
- Python 套件與 Node.js；
- delivery XLSX 所需 `@oai/artifact-tool` runtime；
- 一個**尚不存在**的 output root。

任何 stage 失敗時，不會跳到後面；已完成資料保留在該次 output，並寫出 `run_failed.json`。完整成功才會產生 `run_manifest.json` 與 `12_delivery/`。為避免把舊結果與新結果混合，目前不提供就地覆寫或自動 resume。

## 重要界線

- `rag_re.skip_literature=false` 才是完整 A/B/C live 路徑；設為 `true` 是明確的無文獻分支，結果不等價。
- R5 finalizer 直接載入 clinical config 與 clinical rules；`RAG_re_clinical_pruning*` 是研究比較支線，不在每次主線中重跑。
- 理由是事後、證據受限的可讀說明，不可反向改選菌。流程先審核 upstream rationale，R5 後再產生一次對齊最終選擇的理由。
- 目前第二次 rationale 對每位病人重跑，換取簡單且明確的全批次一致性；因此要先用 `--plan` 核對病人數與預估成本。
- 程式完整接線不等於臨床驗證；最終輸出仍標示 pending clinical review。

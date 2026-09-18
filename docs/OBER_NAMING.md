# OBER 名稱與模組對照

對外研究系統名稱為 **OBER**。目前實作沿用一部分 `RAG_re` 目錄、Python package 與 schema 名稱；它們是現用 OBER 的實際依賴，不是另選回 earlier RAG 系統。

| 現有目錄 | 在 OBER 中的角色 | 此次處置 |
|---|---|---|
| `RAG_re` | 文獻檢索與判讀、候選輸入介面、基礎證據規則；casefit 會匯入 ArticleJudge、PubMedClient 和規則函式 | 保留必要原碼，內部名稱暫時不變 |
| `RAG_re_casefit` | A1 病例適用性與對應資料；R5 載入其 FROZEN_SYNONYMS | 保留 |
| `RAG_re_clinical` | B/C/D 臨床證據規則；R5 載入 evaluate_modules | 保留 |
| `OBER_patient_delivery` | R5 最終決策、理由對齊及交付 | 保留，為最終決策入口 |
| `RAG_re_clinical_pruning` / `RAG_re_clinical_pruning_abc` | 研究對照支線；不是每次 R5 的必要前置步驟 | 獨立保留作研究參考，不標為 R5 主線 |
| earlier `legacy_rag` | 舊 casecard / rag_engine_v3 系統，現用 OBER 未匯入它 | 已移出交付repo，留私人盤點備份 |

## 程式證據

- `OBER_patient_delivery/generate_r5_final_decisions.py` 的 `runtime_modules()` 明確載入 `RAG_re_casefit/rag_re_casefit/rescue_evaluation.py` 和 `RAG_re_clinical/rag_re_clinical/clinical_rules.py`。
- `RAG_re_casefit/rag_re_casefit/batch.py` 直接匯入 `rag_re.judge.ArticleJudge`、`rag_re.pubmed.PubMedClient` 與 `rag_re.rules.module_a_from_judgments`。
- 保存的 artifact 使用 `rag_re_casefit.batch_manifest.v1`、`rag_re_casefit.candidate_output.v1` 等 schema。不能把 schema 字串一起盲目更名，否則已驗證的歷史資料可能無法讀取。

此次調整只變更公開範圍及閱讀入口，沒有改演算法、prompt、threshold、Python import 或歷史 schema。若後續統一內部名稱，應另立純命名遷移版本，逐一修改依賴並以55例凍結結果验证，不把改名視為新方法。

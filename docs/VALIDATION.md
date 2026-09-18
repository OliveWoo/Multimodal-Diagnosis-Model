# 本次驗證與限制

> 2026-09-18更新：已還原並驗證v19 scorer，保存summary＋ranked起算55/55完整內容相同；請用 `scripts/replay_v19_upstream.py`。v20 計分規則保留，共用 model-free helper 已與退役的 LLM mNGS-max 分支解耦。raw→summary不是完整重現。方法見 [版本恢復與操作](UPSTREAM_RECOVERY.md)。

2026-09-17；研究原檔不變。測試與重播輸出都放在私人 audit 目錄或新建暫存目錄；沒有付費 API 請求或 GitHub 上傳。

| 驗證 | 結果 | 實際證明範圍 |
|---|---|---|
| 公開候選 Python suite | 265 tests passed，0 skipped | 上游76；RAG_re36；clinical47；casefit20；pruning47；abc20；delivery8；evidence4；Direct-Raw7 |
| Node delivery model | 8 tests passed | 資料交付模型，未涵蓋 artifact-tool XLSX 真正匯出 |
| synthetic handoff | passed | 合成 merged JSON → input adapter 正確排除 low tier → R5 保留1個picked並補入1個候選；A1訊號為明示 fixture |
| 私有 R5 整批重播 | 55/55完整決策 JSON bytes/hash 相同 | 114 picked＋9 rescued＝123 final，518候選；重播的是凍結上游與A1，不重呼叫LLM |
| 上游銜接 | 55/55 input hashes 相同；55/55在上游交付目錄找到完全相同副本 | 這批 OBER 輸入可追溯至另一位成員的交付內容 |
| A1 lineage | 62/62 artifact、latest merge、frozen RAG hashes 均相符 | 橫跨主線與其他研究批次；62不能當55例分母。15個artifact記錄過fallback，不保證文章內容完全保存 |
| Direct-Raw 評估重算 | 四模型41例，共164份保存prediction；metrics與per-patient均相同 | 使用identity-fixed run與指定修訂gold；不重新生成prediction |
| OBER同41例評估 | 41/41 predicted sets、gold sets吻合；TP69/FP33/FN7 | 僅對這一比較cohort與gold版本有效 |
| 私有Direct-Raw測試 | 23 passed | 原病例identity與答案修訂；包含私有fixture，沒有收入公開候選 |
| shared來源測試（隔離前） | 50 tests run，其中46 passed、4因缺原workbook路徑skipped | 含病例內容之7份測試後續全數隔離；公開候選未假稱保有這50項測試 |
| Table 1 原Excel解析 | 成功；11個JSON中8個與保存內容一致，3個不同 | underlying／diagnosis／mNGS差異需baseline provenance補足，不能宣称全流程重現 |
| Table 1 baseline追查 | 3個來源候選中2個合併後，10個clinical/mNGS section全部吻合 | 證明內容可重建；仍不能僅由內容判定當初用了哪個副本 |
| 原有合成evidence完整入口 | 成功產生8份JSON | run_pipeline的candidate／combined／split確定性流程，無模型呼叫 |

最初移除私有資料後，Direct-Raw兩個setUpClass因cohort manifest缺失報錯；已據此辨識並隔離真實資料測試，公開runner只執行可公開的7項。此類資料依賴沒有用假資料掩蓋。

目前未驗證：付費模型生成、live PubMed、從所有原始Excel/PPT重建歷史clinical資料、完整XLSX匯出、乾淨新機安裝、臨床有效性與公開權利。完整驗證log及逐病例比較保留在私人audit附件；本文件只提供不含逐病例內容的摘要。

## OBER 公開範圍修訂

舊 `legacy_rag` 已移出此交付repo。修改前以Python AST核對，未發現現用程式對該舊引擎的import。移出後重新執行離線測試、合成handoff及SOURCE_MANIFEST核對；保留程式內容不變。

## 2026-09-18 恢復入口驗證

新增8項輸入／規則hash、病例ID、重複case、覆寫拒絕、一般M1/M2版本分歧及比較嚴格性測試。現有與新增測試合計273個Python、8個Node，全部通過。v19 `.py` 入口重新驗證55/55完整scorer及55/55保存merge的選菌順序／候選等級。詳細數據見 [UPSTREAM_RECOVERY_VALIDATION.json](UPSTREAM_RECOVERY_VALIDATION.json)。

## 前端格式工具驗證（前一批）

新增21項前端測試，包含真實讀取合成Excel＋PPT的完整入口、圖片bytes／回覆錯配、缺回覆、重複病例、index改序、路由還原、拒絕覆寫與內容差異回傳失敗。全 repo 合計294個Python、8個Node測試通過，總計302項。

實資料透過最終入口完成36個run：33份舊Excel、1個master＋39份PPT、2個Table 1分支。master標準分項1,070/1,070相同，200頁PPT解析（僅排除raw_path）相同。比較範圍、232份可追既有Raw的空白列差異、32份無歷史Raw，以及Table 1 identity界線見 [FRONTEND.md](FRONTEND.md)。

## 高醫／兩院分線驗證

新增15項分線測試，涵蓋兩種namespace、KH local／global不同號、兩院不可重新編號、跨院誤接、檢體錯配、flat／grouped schema、空mNGS明確標記，以及兩個CLI實際執行。全repo共309個Python＋8個Node＝317項通過。

## 退役 LLM mNGS-max 後的驗證

移除 `mngs_big_agent`、舊 scheduler／batch wrapper、兩份 mNGS-max prompt、舊 LLM 比較／清理工具與其專屬 postprocess 測試後，共用 model-free helper 已抽至 `upstream/tools/mngs_common.py`。清理版全repo共308個Python＋8個Node＝316項通過；v19 synthetic replay 與預期完整內容相同，v20 synthetic scorer 清理前後輸出 SHA256 完全相同。歷史 source snapshot 仍保留於 `reference_versions/`。

其後再移除現行主線自動生成 no-FilmArray agent／summary 的平行分支，新增3項測試固定 with-FilmArray agent、summary 檔名、auto lookup 與 v20 CLI 契約。全repo共311個Python＋8個Node＝319項通過；既有 deterministic summary＋ranked 輸入下，v20 synthetic scorer 清理前後輸出 SHA256 仍完全相同。recovered v3/v19 與 `reference_versions/` 仍保留 no-FilmArray 歷史相容性。

## 標準化 JSON 到 v20 的單一入口

新增根目錄 `run_deterministic_v20.py`，串接 ranked mNGS、with-FilmArray deterministic summary 與 v20 scorer。4項新測試涵蓋 flat／per-patient 輸入、無 `agent_outputs/`、mNGS來源歧義拒絕及輸出覆寫拒絕；全repo共315個Python＋8個Node＝323項通過。另以公開合成兩院 `mNGS_grouped` 與 KH `all_RK_NTC` 各完成一個 root CLI run，皆產生 `final_results.json`、`run_manifest.json` 與每例完整 scorer JSON，全程沒有網路或 LLM 呼叫。此入口終點是 deterministic v20 可能菌種；OBER R5 仍需額外 casefit、文獻與 clinical config，不由標準化 JSON 自動推定。

實資料分線整理：KH33例297份、兩院48例480份，均与指定保存上游JSON內容相同；是讀取／驗證／複製凍結資料，不是LLM重新計算。master到兩院48例的432份clinical全部相同，42份grouped不同、6份皆空。可公開彙總見 [COHORT_VALIDATION.json](COHORT_VALIDATION.json)。

## 完整流程編排入口

新增根目錄 `run_full_pipeline.py` 與 `examples/full_pipeline_config.template.json`。公開合成 KH／two_hospital config 的 `--plan` 成功解析2個 cohort、2個病人與31個固定 stage，且確認 plan mode 不建立 output。新增測試涵蓋完整 stage 順序、merged handoff bytes/SHA256 不變，以及兩 cohort 病人編號重疊時，R5 rationale 仍以 `COHORT=PATH` 明確對應。`audit_llm_rationale_batch.py`、`prepare_r5_rationale_inputs.py` 與 `assemble_r5_rationales.py` 保留原介面並增加 generic cohort mapping。

完整離線套件共319個Python＋8個Node＝327項通過；另實際完成公開合成 KH 的 deterministic → review queue → dry review → merge handoff smoke。這次只驗證編排計畫、接線 adapter 與不需外部服務的階段；沒有呼叫付費模型、live PubMed，也沒有宣稱 XLSX runtime 已可在乾淨環境安裝。完整 run 的臨床與外部依賴界線見 [FULL_PIPELINE.md](FULL_PIPELINE.md)。

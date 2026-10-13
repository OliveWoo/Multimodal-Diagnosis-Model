# 2026-10-13 公開流程更新

本次更新公開程式、通用規則、合成測試及執行方式，不公開病例、答案、模型逐例回應或執行產物。原有根目錄完整流程保留；新版位於 `current_pipeline/`。此次公開與離線測試不代表重新完成臨床效能驗證。

## 這次同步的內容

- 小 agents：空 FilmArray／GM 回報 Unknown；空分子檢驗保留 data gap；影像按來源去重；培養初步報告與最終菌名連結。
- 多次 DNA／RNA：按病人、事件與 test 保留來源，不跨 test 加總 reads；技術重現性不等於獨立臨床證據。
- Taxonomy：自動偵測未映射菌名，經 NCBI、文獻、身分 reviewer、臨床家族 reviewer 與固定驗證器，產生本次 run 的 overlay；未完成分類複核不得進 scorer，不確定者維持 Provisional／Hold。
- 時間窗：同期直接證據使用採檢日前後兩個日曆日；窗外資料與來源稽核不因此消失，也不得直接冒充當次陽性支持。
- 舊版相容：代表 DNA 的選擇、病史路由、受限升級及原有標籤保持獨立；不把新版 Possible 混進舊版 Priority。
- OBER：公開歷史重播，以及 Luna 初審、疑難案例 Sol 複核、最後固定條件驗證的工具。兩階段複核和歷史 OBER 是不同實驗設定，必須明確選擇，不能混合結果。
- Reporting：新版 Picked／Possible／Fallback／Context 與複核優先序是不同概念；mandatory-primary shadow 只補最佳可用候選，不提高其臨床確定度。
- 評估：陽性與明確無病原病例分開處理；答案只在輸出凍結後由 evaluator 載入。已參與規則開發的資料不能宣稱為未解盲外部驗證。

## 新版執行順序

先安裝 `current_pipeline/requirements-dev.txt`、啟用虛擬環境，再切換至 `current_pipeline`。以下私人路徑需自行指定，不是倉庫內附送的病例。

### 1. 小 prompts

`SourceRoot` 必須是轉接器支援的標準化資料包，不是任意原始 Excel。轉接器可讀取 `patients/<cohort>/NGS_patient_*_json`，使用補強後的臨床輸入建立隔離工作區，不複製答案與舊衍生結果。

```powershell
.\scripts\run_small_agents.ps1 `
  -SourceRoot C:\private_data\standardized_bundle `
  -RunRoot C:\private_runs\small_agents_luna_high `
  -Model gpt-5.6-luna -ReasoningEffort high -MaxWorkers 4
```

可加 `-PreflightOnly` 只做單一病人 CBC／其他檢驗預檢。再次執行同一命令會續跑缺少模組；不同模型必須使用新的 RunRoot。程式只使用終端機既有的環境變數，不詢問、不列印 API key。

### 2. 建立 inventory 與摘要

```powershell
python -m tools.build_multi_assay_mngs_evidence `
  --source-root C:\private_data\multi_assay_mngs `
  --output-root C:\private_runs\inventory --expected-patients 10

python -m tools.build_deterministic_summary `
  C:\private_runs\small_agents_luna_high\patients `
  --output-suffix decision_input
```

`expected-patients` 必須改成實際人數。缺少候選或來源資料時先檢查 inventory 的 QC 與事件對照，不為了取得結果而自動填值。

### 3. 自動分類與 deterministic scorer

```powershell
.\scripts\run_taxonomy_complete_pipeline.ps1 `
  -InventoryRoot C:\private_runs\inventory `
  -TaxonomyRoot C:\private_runs\taxonomy_review `
  -PatientRoot C:\private_runs\small_agents_luna_high\patients `
  -SourceSummarySuffix decision_input `
  -UpstreamExtractionModel gpt-5.6-luna/high `
  -OutputRoot C:\private_runs\integrated_release `
  -Backend codex-cli -IdentityModel gpt-6-luna `
  -ClinicalModel gpt-6.1-sol -TaxonomyReasoningEffort high
```

使用 API backend 時改為 `-Backend openai-api`。實際能否呼叫取決於帳號模型權限及服務端支援，這份公開程式不保證所有帳號可用相同模型。

PhenotypeRoot 有資料時才提供。若只有 mNGS，省略 PatientRoot／SourceSummarySuffix，改加 `-MngsOnly`，使缺失臨床資料明確標為 unknown。`-TaxonomyOnly` 可只跑分類複核。TaxonomyRoot 若已屬於不同 inventory 或模型設定，程式會要求新的目錄，避免混用快取。

公開預設的凍結清單是 `current_v6_v9_calendar_v7_frozen_20261013`。它對應目前程式的日曆日時間窗，不代表舊的 10/8 精確小時版本；沒有沿用舊版凍結雜湊宣稱等價。

## 舊版與 OBER

完整原始流程仍使用根目錄 `run_full_pipeline.py` 和 `examples/full_pipeline_config.template.json`；可先加 `--plan` 檢查執行計畫。舊版相容 runner 使用 `--contract` 指定私人資料路徑與病人清單，不能直接套用開發資料的預設路徑。

歷史 OBER 重播先準備已建立的 answer-blind 病例包及其 manifest：

```powershell
python -m tools.run_historical_ober_r5_replay `
  --cohort-config C:\private_data\ober_cohorts.json `
  --output-dir C:\private_runs\historical_ober --prepare-only
```

私人設定的結構可參考 `../examples/historical_ober_cohorts.template.json`；相對路徑以設定檔位置為準。移除 `--prepare-only` 才會送出模型請求。這裡的 OBER 包含證據複核及固定接受條件，程式檔保留 R5 名稱是歷史相容性，不是另一套獨立流程。

新版兩階段 OBER 則先以 `tools.adjudicate_rag_v2_packets` 完成 Luna，經 `tools.select_ober_sol_escalations` 選取疑難案例，送 Sol 後再用 `tools.merge_ober_two_stage_adjudications` 合併；各工具的 `--help` 說明輸入與 schema。Sol 不能繞過 clinician-visibility boundary 或固定接受條件。分類用的 identity／clinical-family reviewers 與個案 OBER reviewers 不是同一階段。

## 驗證與公開邊界

```powershell
python -m pytest tests -q
python -m tools.validate_frozen_policy_bundle
```

`PUBLIC_SOURCE_MANIFEST.csv` 保存公開檔案位元組的 SHA-256。根目錄的 `SOURCE_MANIFEST.csv` 同時保留 `original_source_sha256`；公開檔案使用 Git 的 LF 換行，避免 Windows CRLF 讓乾淨 clone 的驗證失敗。範例使用的封裝 profile 已重新綁定，舊 profile 保存在 `reference_versions/publication_profiles_pre_20261013/`，不把封裝雜湊修正當成歷史規則重跑。

其他 reporting／升級規則屬於可選研究設定，不能把檔案都存在解讀為同一次執行全部使用。任何結果比較仍需一致的資料、答案版本、輸出範圍與 matching policy。

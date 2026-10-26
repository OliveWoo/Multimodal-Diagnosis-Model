# 最新整合流程：2026-10-20 邏輯，2026-10-26 公開同步

本次同步的是目前新版研究流程的程式、通用規則與合成測試。新版入口不包含 OBER，不公開病例、醫師答案、逐病人結果、來源補充內容或模型快取。根目錄舊流程及 OBER 工具保留為歷史功能，沒有混入這個新版入口。

## 更新內容

- 統一底層預設：analytical v6、history v4C、clinical v1、promotion v9、reporting v7。
- 接上病毒來源 facts／claims 與檢出、活動／再活化、器官支持三軸。
- CMV 使用來源可核對、同期、獨立呼吸道檢出的報告路徑。排名、RPM、DNA／RNA 組合及技術重測次數不直接成立 CMV 肺炎。
- 沒有方法資訊就保留 unknown，不推定 PCR／qPCR、不推定病毒活動；血液陽性不直接證明肺炎。
- 全病例重算 fallback：有一般 Picked／Possible 就不額外補 Fallback-Possible。Fallback 保留為完整陽性端點的一部分，但不提高疾病確定度。
- 多瘤病毒獨立於疱疹病毒 family，保留持續感染／脫落及器官支持限制。
- 加入已核准的 Arthrobacter woluwensis 身分映射，不因身分已確認而直接升級。
- 輸出凍結後才由獨立評估器讀答案。未匹配候選、資料缺失及明確無病原病例不能從統計中省略。

這些是研究整合與可追溯性修正，不代表新的臨床確診標準，也不宣稱本次發布提高 F1。

## 公開版與本機研究版的差別

本機原入口固定使用特定資料批次、分組名單、既有來源確認及答案檔案。這些不能放進公開倉庫。本次將私有內容外移為必要的 `--config`，通用判讀函式沿用研究版；沒有重新調整選菌門檻。

公開版另存帶有 `20261020` 的政策及 taxonomy 檔名。政策 ID 加上日期辨識，`derived_from_policy_id` 保留原 ID。這是版本與公開包裝差異，不是新的臨床門檻。10/13 已發布的凍結政策檔案及舊 bundle 保持原樣，不拿舊雜湊替新規則背書。

個別來源確認一律由私有設定傳入，沒有公開任何特定病例的確認內容。空模板只表示沒有提供補充，不表示方法獨立，也不能用它重現原本有補充的研究病例。

## 執行前準備

在 repository 的 `current_pipeline/` 安裝相依套件。上游臨床抽取、inventory 與 taxonomy completion 要先完成，詳見 [10/13 上游步驟](PIPELINE_UPDATE_20261013.md)。本次入口不會重新呼叫全部 LLM，也不能直接讀取任意原始 Excel 或未轉換資料夾。

請自行準備以下私有檔案，不要 commit：

1. 每組基底流程的 `run_manifest.json`，含已完成的 inventory、病人資料、phenotype 路徑與 summary suffix。
2. 完成驗證的 taxonomy completion manifest 及其 overlay。
3. 感染事件分離設定。沒有已確認群組時可明確提供 `{"clusters": []}`，不能猜測或合併事件。
4. 已核對的 source clarifications，沒有補充時明確提供 `{"clarifications": []}`。來源列、內容雜湊、物種及事件必須能唯一對回。
5. 相同資料、相同基底邏輯的凍結 reference 結果，供 baseline replay 比較。不同資料或版本不能強行繞過差異。
6. 原本開發分組的 `answer_blind_freeze.json`，供既有 runner 保存分組標記。這份檔案只提供開發病例 ID 集合，不是答案，也不代表其他資料現在仍是未參與開發的外部驗證。

### 私有 inference 設定示例

以下路徑只是佔位示例，沒有病例內容：

```json
{
  "schema_version": "public_integrated_run.v1",
  "source_manifests": {
    "COHORT_A": "C:/private_runs/cohort_a/run_manifest.json"
  },
  "taxonomy_completion_manifest": "C:/private_runs/taxonomy/taxonomy_completion_manifest.json",
  "episode_identity_file": "C:/private_data/episodes.private.json",
  "source_clarifications": "C:/private_data/source_clarifications.private.json",
  "reference_root": "C:/private_runs/frozen_reference",
  "development_freeze": "C:/private_runs/development/answer_blind_freeze.json"
}
```

設定內的路徑可相對於設定檔解析。來源 manifest 中的 inventory／patient／phenotype 路徑與 completion 中的 overlay 必須是存在的絕對路徑。每組 baseline 預期位於 `reference_root/<input-run>/patient_results/`。

私有來源及執行輸出可放在公開 checkout 之外。Taxonomy 來源註記會保存其可追溯路徑，不要求將臨床快照搬進 Git。

如需將一個輸入批次分成多個輸出組，可在私有設定加 `partitions`，格式為「輸入組：輸出組：病例 ID 清單」。入口要求名單不重疊、完整涵蓋 inventory；未分組時直接沿用輸入組名稱。公開程式沒有內建研究病例分組名單。

```powershell
cd current_pipeline
python -m tools.run_current_integrated_pipeline `
  --config C:/private_data/integrated_run.private.json `
  --output-root C:/private_runs/integrated_cmv_new_run
```

輸出目錄必須不存在。設定不完整、分類 overlay 雜湊錯誤、基底結果不一致、來源錯接或 fallback 衝突時停止，不憑空補證據。

## 凍結後評估

另外建立只供評估使用的私有設定：

```json
{
  "answers": {
    "COHORT_A": "C:/private_answers/cohort_a.csv"
  }
}
```

```powershell
python -m tools.evaluate_current_integrated_cmv `
  --run-root C:/private_runs/integrated_cmv_new_run `
  --answer-config C:/private_answers/evaluation.private.json
```

評估器先檢查凍結檔案及 inference 設定，再讀答案。答案組別必須與本次輸出相符。計分沿用既有 alias／group 契約，不新增答案導向匹配。

正式陽性為 Picked＋Possible＋Fallback-Possible，Context 不列陽性；明確無病原病例的額外預測也計 FP。統計單位是病例－病原配對，不是病人 accuracy。

## 實際順序與輸出

1. 驗證私有設定與來源，凍結程式、政策及設定雜湊。
2. 重新執行 compact entry、analytical、timeline／phenotype／history、clinical、rehydration、promotion 及 reporting 基底。
3. 與指定 reference 比較所有候選及基底入選名單。
4. 從本輪來源接回全部相關病毒候選的 facts／claims，核對來源 pointer 與內容雜湊。
5. 對 CMV 套用證據路徑，其餘病毒保留註記但不自動套 CMV 疾病規則。
6. 從全病例候選重算互斥 fallback，產生一致的逐病人 JSON 及 CSV。
7. 凍結結果，另行載入答案計分。

每組 `patient_results/` 保存最終名單與全部候選，`reporting/complete_report.csv` 保存正式陽性。病毒來源在 `02_viral_evidence_binding/`，`answer_blind_freeze.json` 保存輸出雜湊，`postfreeze_evaluation/` 另存統計與逐候選差異。

來源快照、設定、報告與結果仍含臨床內容，只能存放在授權的本機私有位置，不應提交到公開倉庫。

## 方法界線

- mNGS 數值處理是程式分支，與五個臨床語意 agents 分開。
- 多 test 保留逐筆向量，不跨 test 加總 reads；技術重測不是獨立臨床證據。
- 同期窗口是採檢日前後兩個日曆日，不是精確到小時的 ±48h；報告時間不能代替採檢時間。
- 缺少資料維持 unknown，不當陰性；mNGS reads／RPM 不等於 qPCR 病毒量。
- 系統 Picked／Possible／Fallback 是研究報告層級，不等同臨床 proven／probable。
- 通過合成測試與來源稽核不能證明臨床安全或排除過擬合。已參與開發的資料不因重跑而成為盲性外部驗證。

## 離線測試

```powershell
python -m unittest discover -s tests
python -m tools.validate_frozen_policy_bundle
python -m tools.validate_frozen_policy_bundle --bundle rules/current_pipeline_policy_bundle_integrated_cmv_20261020.json
```

第二個命令仍驗證歷史 10/13 bundle，第三個驗證這次的新 bundle。新版本政策、taxonomy 與相同判讀函式的來源雜湊另見 `current_pipeline/publication_manifest_20261020.json`。無病例、答案、API key 或外部模型服務仍可跑合成測試；不表示沒有私有資料也能重現研究效能。

2026-10-26 公開前核對：unittest 共 611 項，5 項因既定外部／私有條件跳過，其餘通過。新舊政策 bundle 均通過雜湊驗證。另以授權本機私有輸入完成整合重播與凍結後評估，候選、入選名單及各組 TP／FP／FN 與原凍結結果一致。這些臨床輸入與重播產物沒有加入公開倉庫，也不構成新的獨立外部驗證。

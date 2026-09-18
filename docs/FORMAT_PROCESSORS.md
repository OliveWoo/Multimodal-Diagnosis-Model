# 共用格式工具：Excel／PPT／Table 1

> 正式整理層級是[高醫／兩院兩條主線](FRONTEND.md)。107筆master包含K53、T19、V35；33份Excel路徑不是33位高醫病人。

本次已把三種實際存在的輸入格式接成 `frontend/pipeline.py`，以原始 Excel／PPT 實測，保留原研究程式。入口在獨立程序載入 `reference_versions/shared_upstream`，不會把它與 `upstream/tools` 混成同一個 Python package。所有輸出寫入全新目錄。

**本頁只證明格式層的抽取／重播；不代表高醫與兩院兩條主線已從原始資料完整重建。** master＋PPT 分支可離線重現保存的標準化資料；其中圖片文字使用保存的模型回覆。沒有重新呼叫 vision、formatter、clinical agent 或 OBER 的付費模型，也沒有把 107 個檢體病例認定為論文的 55 例 cohort。

## 可直接執行的公開範例

在 repo 根目錄，以 Python 3.13 與 `requirements-audit.txt` 的依賴執行：

```powershell
python -B frontend/pipeline.py --manifest examples/frontend/manifest.synthetic.json --output local_outputs/frontend_demo_001
```

輸入是完全虛構的四分頁 Excel 與一頁原生文字 PPT。輸出包含 1 個檢體病例、1 份影像文字紀錄、mNGS grouped records、10 份標準分項、`source_mapping.json` 與 `legacy_patient_index.json`。不是臨床範例，也沒有生成文獻或 final decision。程式拒絕覆寫既有輸出；重跑請換目錄名。

## 共用格式工具與子步驟

下列程式除新入口外，都位於 `reference_versions/shared_upstream/`。

| 分支 | 實際程式與順序 | 輸入／輸出與銜接 |
|---|---|---|
| 跨病例 master | `tools/split_excel_master_by_patient.py` → PPT 匯入 → `tools/export_legacy_patient_json.py` | `基本資料`、`微生物診斷資料`、`mNGS`、`Read count` → 以 case code／specimen ID／hospital ID 組成的病例鍵 → `NGS_patient_<ordinal>_json`。ordinal 取決於 patient index 順序，不能跨 cohort 直接配號。 |
| PPT 補入 | `ppt_patient_parser/cli.py` → `tools/link_ppt_slides_to_excel_patients.py` → `tools/import_ppt_labs_to_patient_exports.py` | PPT 原生文字與逐圖片抽取 → `slides/slide_*/parsed.json`、`debug.json`、`images/` → 病例 bundle、image、CBC、other_lab。以 mNGS 時間窗篩選，保留各來源。 |
| 分人多分頁 Excel | `core/raw_parser.py`、`parsers/*.py` → `core/llm_parser.py` 的 prompt 組裝 | 每份 Excel 拆出 8 份 Raw JSON。新入口依 Raw 的 `sheet` 在 shared 的 `prompts/` 精確找提示詞，產生待執行請求；不呼叫 `send_to_llm`。只有明確提供相同 Raw 內容的凍結回覆時才重播 normalized JSON。 |
| 單頁 Table 1 | `tools/parse_table1_workbook.py` 與 `tools/parse_table1_standardized.py` | 同一工作表分別產生 Raw 分段與確定性標準化結果；不是一定要先 LLM formatter。可明確指定 baseline 補保留 mNGS 等欄位，baseline 與原檔都不修改。 |

**更正文檔錯誤：shared linker 讀的是 `parsed.json`，不是先前報告寫的 `slide.json`。** 成員交付的 `upstream/ppt_patient_parser/cli.py` 只寫原生 `slide_text.txt` 與抽出圖片，没有這套 parsed／vision／linker 銜接。因此本入口使用 shared 分支；不是依檔名或修改時間選版本。

`reference_versions` 的名稱不代表此部分已廢棄；它是本次驗證過的前端依賴。逐檔原始來源與 hash 仍在 [SOURCE_MANIFEST.csv](../SOURCE_MANIFEST.csv)，新增入口／固定清單在 [FRONTEND_CODE_PROFILE.json](FRONTEND_CODE_PROFILE.json)。

## 實資料驗證結果

| 檢查 | 結果與比較範圍 |
|---|---|
| Office 檔盤點 | 93 個 Excel／PPT 路徑、80 組不同內容 hash；逐檔路徑與用途留在私人 registry。包含前端原料、重複副本及評估／輔助工作簿，不把所有 Excel 都當原始病例輸入。 |
| master identity | 三份不同 bytes 的 master 均產生同樣 107 個檢體病例及相同 index 順序。shared master 與成員 `Specimen/` 副本的抽取臨床內容相同；成員根目錄版本的部分 culture／FilmArray／GM／molecular 內容不同。 |
| 完整 PPT | 40 個來源路徑中 39 份內容不重複；200 頁、212 張圖片。所有圖片 bytes 與對應保存回覆的圖片相符；200/200 parsed JSON 相同，僅排除輸出位置 `raw_path`。 |
| master＋PPT → standardized | 新入口重算 107 × 10 = **1,070/1,070 份分項 JSON 完全相同**，沒有排除臨床欄位或排序；比較為 JSON 內容，不宣稱序列化 bytes 相同。 |
| 下游端保存副本 | 與 OBER 專案 `LLM_test/data/patient_info_standardized` 另存副本也是 1,070/1,070 相同。證明介面／內容相接，不證明之後所有 LLM agent、ranked、merge 都重新生成。 |
| 舊格式 Excel | 33 個來源路徑（25組不同內容、27個檔名病例標籤），全部拆出 264 份 Raw JSON。230 份與既有 Raw 完全相同；另 2 份只有空白列差異，去掉全空白列後內容相同；32 份找不到保存 Raw 對照，不以別人的結果填補。已產生 264 份 formatter 請求，API 呼叫數為 0。 |
| Table 1 | 無 baseline 可產生 10 份分項與 mapping；使用候選 baseline 時，10 份分項＋mapping 共 11/11 保存 JSON 相同。這只支持內容相容；尚無歷史 manifest 證明 baseline 病人身分與當時選用來源。 |

### 為何 PPT 有兩種匯入方式

保存的 38 份批次 PPT 走 `radiology_and_labs`。另一份獨立 PPT 有保存的 19 頁完整圖片回覆，頁面也出現在檢驗匯入 bundle，但不出現在保存的 image 分項來源，因此重播 manifest 對它指定 `labs_only`。

如果把這份獨立 PPT 也完整匯入影像，會在 14 個病例新增 35 筆 image，僅有 1,056/1,070 分項相同。按來源指定 `labs_only` 後才達到 1,070/1,070。這是由 bundle、image provenance 與實際內容比較推得的重播設定，**不是找回了當年的完整命令列紀錄**。程式沒有依病人 ID 刪資料，也沒有讀取 gold 配答案。

### 未匹配與多檢體對照

200 頁中 28 頁沒有匹配病例，15 頁會依相同 hospital ID 連到多個檢體病例；共 187 個頁面到病例連結、54 個病例有 bundle。未匹配頁中有 11 頁能解析出 patient ID，但仍不能據此假定一定属于這批 master。新入口逐頁保留 keys、reason、來源 hash 與 route；重播 manifest 固定整份 binding 的 hash，避免靜默改配。

這是研究原有解析與對照行為的重現，不能把未匹配頁視為已納入，也不能把同病人的不同檢體當成重複刪除。lab import 與原保存一致：寫入 CBC 784 筆、other_lab 2,283 筆；時間窗排除資訊留在 receipt。

## Manifest 使用方式

所有檔案參照都是 `{ "path": "...", "sha256": "..." }`。相對路徑以 manifest 所在目錄為基準；可用 PowerShell `Get-FileHash -Algorithm SHA256 -LiteralPath <檔案>` 取得 hash。不得只填同名檔或依最新日期自動挑資料。

共同欄位：

- `schema_version: 1`、`mode`、`run_kind`、`code_profile_sha256`。
- `run_kind: "extract"` 表示新抽取；`"replay"` 必須帶明確的保存輸出 `expected`，不能用「成功執行」代替重現判定。
- `code_profile_sha256` 是 repo 的 `docs/FRONTEND_CODE_PROFILE.json` SHA256；它固定前端入口、原研究模組與提示詞。修改程式需另建 profile／run，不沿用舊重播證明。
- `expected` 是陣列，每項含新輸出內的相對 `output` 路徑與保存檔的 `reference`。程式會逐檔比較 JSON；差異回傳非零 exit code。

| mode | 專用欄位 |
|---|---|
| `master_ppt` | `master`；`presentations[]` 的 `source`、`route`、`frozen_slides[]`。重播另需 `expected_identity_sha256`（index 中 case code／specimen ID／hospital ID 的有序清單）與 `expected_binding_sha256`（有序逐頁對照清單）。公開範例展示可執行的最小格式。 |
| `legacy_excel` | `workbook`。可选 `frozen_formatter[]`：每項為 `raw_name`、`raw_semantic_sha256`、`response`、`output_name`、可選 `prompt`；重播 normalized 時另需 `identity_namespace`，且必須完整涵蓋該工作簿所有 Raw。 |
| `table1` | `workbook`、`patient_id`、`identity_namespace`、可選 `sheet`（預設 `Table 1`）。可選 `baseline` 是10個標準分項 suffix 到檔案參照的映射，並須填 `baseline_identity_evidence`。文字理由不是病人身分證明；若未證實，只能標記為相容性研究。 |

圖片回覆必須使用完整保存的 `parsed.json`，其 `images/` sibling 中須有保存的原圖片。程式檢查 PPT SHA、來源檔名、slide／image index、實際抽出圖片 SHA、保存回覆狀態及 parsed 內容。任一圖片沒有回覆、回覆失敗、圖片不符、slide 被原 parser 跳過，都使該次執行失敗；不會轉用付費 API，也不把空 OCR 當成功。

`identity_sha256`／`binding_sha256` 使用 `frontend.pipeline.semantic_sha` 的 canonical JSON 定義。建立新資料集應先以 `extract` 審閱私人 index／bindings，再固定供後續重播；這是人工確認來源的步驟，不能只從當次輸出自動鎖定便當成歷史正確性證據。

## 輸出與上游銜接

| 新輸出 | 用途 |
|---|---|
| `patient_exports/` | master 拆分、PPT 匯入的隔離工作目錄；原始 patient exports 不動。 |
| `ppt/` | 每份來源單獨輸出 parsed、debug、原圖片與抽取文字。 |
| `standardized/NGS_patient_*_json/` | 上游 clinical agent 所需的 section JSON、`mNGS_grouped`、`source_mapping`。 |
| `standardized/legacy_patient_index.json` | ordinal 與 case／specimen／hospital 三元 identity 的對照，須隨私人資料保存。 |
| `raw/`、`formatter_requests.private.json` | 舊 Excel 或 Table 1 的抽取材料；待模型處理請求只在舊多分頁分支產生。請求不是模型回覆。 |
| `ppt_bindings.private.json`、`ppt_batches.private.json` | 各來源 route、逐頁對照、未匹配、fanout 與匯入摘要。 |
| `comparison.private.json`、`receipt.private.json` | 比對結果、程式與所有已讀輸入 hash、模型呼叫數與執行狀態；失敗的部分輸出保留供追查，不標為成功。 |

之後仍須依 [WORKFLOW.md](WORKFLOW.md) 進行 cutoff／agent／summary、mNGS ranked、v19 scorer／review／merge，再交 OBER。不能直接把 standardized 當作 OBER merged input。既有 55 例的 main replay 仍使用 [UPSTREAM_RECOVERY.md](UPSTREAM_RECOVERY.md) 明確指定的凍結 summary／ranked，沒有被這次新抽取輸出取代。

### 特別保留的版本界線

- 三份 master 不能只看檔名選最新。shared master 與 `Specimen/` 版抽取相同，但 bytes 不同；沒有證據可唯一判定當年究竟讀哪一份。主要重播固定已具保存對照的 shared 檔。
- `TSGH_v4.xlsx` 缺 master parser 需要的 `mNGS`／`Read count` 工作表，屬另一種來源／參考表。其19組檢體／病人配對在master中都能找到，但基本／微生物表的部分共同欄位不同，mNGS共同欄位只有9/19列完全相同。新入口會拒絕當成相同master，未擅自改欄位或混入不同內容。
- 舊 Excel 的 shared prompts 已配對可產生請求；成員上游根目錄缺 prompts 的問題由明確選用 shared frontend 解開。但當年的模型版本、請求／回覆與正式論文 run 仍不因目前提示詞存在而獲證明。
- Table 1 baseline 必須保留「內容相容但身分來源待確認」標記，不能只因同叫 patient 1 就跨資料集補資料。

## 公開與私人資料

repo 只新增入口、測試、文件、程式／prompt hash 與人工合成 Excel／PPT。實際 Excel、PPT、截圖、parsed、模型回覆、病例對照、formatter requests、逐例 receipt／差異都留私人盤點目錄。`.gitignore` 只放行兩份具名合成 Office 範例，並忽略 `*.private.json`。

格式測試通過不代表兩條cohort歷史鏈或整個研究已可以公開發布：論文正式 run、共同作者授權與下游已知交付問題仍列於 [PUBLICATION.md](PUBLICATION.md)。

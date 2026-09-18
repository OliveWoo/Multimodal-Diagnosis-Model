# 上游版本恢復與重播入口（2026-09-18）

已處理先前「現況 v20 與保存 v19 結果不一致」的 scorer 版本問題。新增可閱讀的 v19 還原程式與固定輸入重播入口，**從保存的 summary＋ranked mNGS 重算，55/55 完整 scorer 內容相同**，包括版本、候選、選菌順序、等級、規則與理由；僅排除來源路徑容器。

這支持此入口重現這批保存的 scorer 結果。**不代表從原始 Excel/PPT 重新呼叫模型後，已重現整篇論文，也不代表臨床正確性已獲驗證。**

## 使用哪個程式

| 用途 | 程式 | 狀態與保留理由 |
|---|---|---|
| 重播已核對的 v19 結果 | `scripts/replay_v19_upstream.py` | 建議入口；明確指定 manifest 與全新輸出目錄，先核對資料與程式／規則 hash；不自動找同名替代檔 |
| v19 計分邏輯 | `upstream/tools/deterministic_mngs_max_scorer_v19_recovered.py` | 從提供的 CPython 3.13 快取還原；計分函式保留，之後僅將 model-free 共用 helper import 改指向 `mngs_common`；不是找回原始文字／commit |
| 現況 v20 計分邏輯 | `upstream/tools/deterministic_mngs_max_scorer.py` | 後續版本的計分規則保留；共用 helper import 已與退役的 LLM agent 解耦，不能代替 v19 重播入口 |
| 快取中的 summary v3 邏輯 | `upstream/tools/build_deterministic_summary_v3_recovered.py` | 50 個函式及整個模組相符；供歷史差異研究，尚不作為 55 例完整 raw→summary 重建入口 |
| 現況 summary v4 | `upstream/tools/build_deterministic_summary.py` | evidence preservation 邏輯保留；共用菌名 helper 改由 `mngs_common` 提供 |
| 快取對照工具 | `scripts/verify_recovered_cache.py` | 以指定快取 SHA256、bytecode magic 先確認來源，再比對編譯內容；不執行快取 |

歷史 byte-identical 快照仍保留於 `reference_versions/`。現行 `upstream/` 的 deterministic runtime 已做移除 LLM mNGS-max 的結構重整，`SOURCE_MANIFEST.csv` 會區分原樣快照與 refactored runtime；兩份快取還原檔另記於 [RECOVERED_SOURCE_MANIFEST.csv](../RECOVERED_SOURCE_MANIFEST.csv)，不冒充原始歷史 source 文字。

## 為什麼這次有依據

1. 在交付的 `tools/__pycache__` 找到內部為 v19 的 scorer 快取；`.bak` 中找到的相關備份內容則是 v7。沒有用檔名或日期決定版本。
2. 另外只讀檢查 14,306 個 Git 物件，對其中 2,384 個符合大小範圍的 blob 搜尋計分／summary 原碼；找到的 13 份候選均不等同這份快取。這不是宣稱每個 Git blob 都做全文比對。
3. 按指令差異還原原碼，並比較 bytecode、常數（含巢狀函式）、名稱／區域變數／closure、flags、參數、stack size、exception table；只排除檔名與除錯行號資訊。scorer 130/130、summary 50/50 函式及整個模組一致。
4. scorer 的直接研究依賴已固定；原本內嵌於 `mngs_big_agent` 的 model-free I/O、菌名與去重 helper 現已原樣抽至 `mngs_common`，v19/v20 scorer 均不再依賴已退役的 LLM mNGS-max 程式。外部 JSON 規則沒有歷史 hash；本次固定現存版本後，仍須以 55 例重播確認結果相符。
5. 使用 `.py` 還原版與新入口再次重算，55/55 完整 scorer 內容相同；55/55 的選菌順序、候選集合及等級與保存 merge 中的 deterministic 區塊相接。程式不讀 gold，不按病例 ID 配答案。

公開的 [來源證據](UPSTREAM_RECOVERY_PROVENANCE.json) 與 [程式／規則固定清單](UPSTREAM_RECOVERY_PROFILE.json) 包含 hash，沒有原始病例內容或 `.pyc`。所有逐病例輸入、差異、輸出與快取仍在私人盤點資料中。

## v19 與 v20 差在哪裡

典型肺炎細菌有下呼吸道 culture／FilmArray 支持時，v19 的相應分支直接設為 `M2_moderate`；v20 在已滿足 `m1_strong` 條件時保留 `M1_strong`，另加入中央 taxonomy 欄位。這會改變 Level，繼而影響以 Level 排序的 picked 清單。

這批資料以現況 v20 重算：選菌集合 55/55 相同，但 15 例各有 1 個候選等級不同，5 例 picked 順序不同。v19 還原版解決的是**歷史結果重播的版本配對**，不是對 v20 規則好壞作臨床裁決。兩版都保留，入口明確分開。

## 公開的合成範例

在 repo 根目錄執行：

```powershell
python -B scripts/replay_v19_upstream.py --manifest examples/upstream_v19/manifest.synthetic.json --output local_outputs/v19_demo_001
python -B scripts/demo_offline.py --output local_outputs/ober_demo_001.json
python -B scripts/check_offline.py --output local_outputs/check_001
python -B scripts/verify_sources.py
```

第一個指令示範 summary＋ranked→v19 scorer，第二個示範另外的合成 merge→OBER adapter→R5；**它們是兩個階段範例，沒有宣稱第一個輸出已經產生 LLM review、A1 證據或第二個範例的輸入。** 輸出目錄若存在，重播指令會拒絕覆寫。

合成病例 900001 是人工構造，病原名稱採用實際規則涵蓋的菌名以測試一般條件，不來自病人，也不是臨床範例建議。expected scorer fixture 由已核對快取產生，用來示範介面；歷史重現證據來自私人 55 例，兩者分開。

## 私人 55 例重播

本地另外建立 `private_upstream_replay_20260918`，內有 manifest、原始來源 hash 對照、172 份去重後的 JSON 內容複本（7,539,456 bytes）。這些複本固定的是**2026-09-18 核對的現存資料**，不補造當年未記錄的 hash。保持私人，不放進公開 repo。

```powershell
python -B scripts/replay_v19_upstream.py --manifest <私人快照目錄>/manifest.json --output <全新的私人重播目錄>
```

每例明確指定 patient ID、ranked、summary、expected scorer、frozen merge 的路徑及 SHA256；先驗證全部輸入，再進行計算。hash 不符、病例不符、重複 case key 或輸出已存在均拒絕；內容比對不符則輸出結果與 receipt 並回傳非零狀態。receipt 記錄每個輸入 hash、程式／規則 profile hash、比較範圍與結果。

重播只重新計算 scorer。保存 merge 的 advisory review／taxonomy 仍是凍結輸入；此指令核對其 deterministic 銜接，不假裝重做 LLM review 或重建完全相同的 merge JSON。**不要直接拿重新序列化的 scorer／merge 取代 OBER 已凍結的來源檔**，因為下游會檢查原始檔案 bytes/hash。

## 更前面的 summary 仍有何限制

快取 v3 summary builder 從目前 patient directory 重建：29/55 完整非路徑內容一致；23 例只有新增 `D-SUM-01B`（GM 證據合併說明）這一行；另 3 例的院內證據不同，其中 1 例亦有 FilmArray/GM probable_source 差異。把新建 summary 接入還原 scorer，52/55 完整 scorer 相同。

因此主要重播入口明確使用已保存、已核對 hash 的 summary；它沒有將目前 agent 資料重新生成的 summary 偷換為論文輸入。3 例仍需查當時 agent／GM 輸入快照及執行設定；既有 prompt 配對、原始資料 identity／baseline、付費模型生成與最終 XLSX 匯出環境等缺口仍見 [資料介面](DATA_CONTRACTS.md) 與 [發布前待辦](PUBLICATION.md)。

這 55 例與指定論文圖表是否完全相同的 cohort/run，仍需論文實際採用的 experiment registry 證明。授權、真實資料公開權限亦尚未確認。

# GitHub 發布準備

候選資料夾採用程式／規則／prompt allowlist，沒有複製真實病例檔、API回覆、研究報告、封存包、金鑰、node_modules或Git歷史。13份含實際或不能确认合成的病例內容之程式／測試，另存私人附件。程式本身亦可能藏有病例值，因此副檔名不是隱私判定依據。

本次 token-pattern 掃描是輔助檢查，不是授權／匿名化證書。只考慮此資料夾的明列檔案；父目錄的詳細盤點、hash表、來源路徑、逐病例log、PDF文字／預覽與quarantine資料都屬私人附件。

## 發布前需要完成

| 項目 | 目前狀態 | 要補上的具體資料 |
|---|---|---|
| 論文版本 | 待確認 | 教授主表／圖對應哪個run；55例R5／41例比較／新單例重跑各自用途；experiment registry |
| 高醫／兩院前端 | 兩條主線目錄、保存輸入整理入口及各自驗證已補上 | KH297、兩院480份保存分項比對通過；缺失Excel、formatter歷史、兩院54→48選取理由與grouped生成recipe仍待補證據 |
| 流程可攜 | 主線離線測試通過 | 新環境安裝驗證、消除腳本硬編碼本機預設、補明確CLI範例與config驗證 |
| 文獻／理由追溯 | 仍有已知缺口 | fallback文章快照、定位evidence、修正export欄位與rationale hash；修正版另立版本 |
| 交付XLSX | 未完成可攜安裝 | artifact-tool取得與授權，或另立經驗證的可公開匯出器 |
| 原始碼權利 | 未確認 | 共同作者／各修改者／所屬機構授權；選定LICENSE、作者與citation |
| 病例資料 | 留在私人目錄 | 資料提供方／研究倫理與發布許可、完整去識別審查；以合成樣例作公開quickstart |
| 文獻內容 | 未公開附帶 | PubMed metadata、摘要、全文／embedding cache各自權利與再散布範圍 |

## 不上傳的類型

- `env.env`、任何API token、key、私有憑證；原資料中已發現明文token型態，應由持有人撤銷／更換，勿把值貼進報告或Git。若曾入Git，需另查歷史；本次未改寫歷史。
- 真實病例Excel/PPT/PDF（具名合成frontend範例除外）、原始與normalized JSON、檢體／醫院代號對照、日期、出生資料、醫師簽名／執照字樣、模型輸入與回覆、gold、逐病例報告。
- 一次性病例override、內嵌實際病例的測試／gold斷言。可在後續另寫完全合成的測試，不改名便當作已匿名化。
- 重複ZIP/RAR/7z/TAR、cache、venv、node_modules、工作表與圖片大量衍生檔。研究保存可以留私人封存，但不整包放進程式repo。
- 未確認授權的第三方程式、文章全文／摘要副本、影像與報表；在權利清楚前不得自行標為MIT或其他開源授權。

## 建議補的文檔／範例

本次已補流程、介面、版本、環境、離線驗證、公開資料政策、逐檔來源／依賴表、空環境變數範例及合成handoff。已新增master＋PPT的完全合成Office檔及前端入口；後續可補其他格式的合成工作簿、A1文章fixture、portable delivery config、模型與PubMed凍結manifest規格、可重算論文彙总表的腳本、CI與正式CITATION。

不新增虛構作者、授權、DOI或論文結果。此階段不初始化遠端、不推送，也不把「測試通過」當成所有權與病例發布許可。

## 上游恢復版新增檔案

v19還原碼已驗證對應這批保存scorer結果；不因此解除原作者授權、論文run指定或原始病例公開權限。新增示範僅為合成資料。`private_upstream_replay_20260918`、`.pyc`與逐病例receipt仍是私人附件，不能加入公開repo。兩份還原碼的來源和hash在 `RECOVERED_SOURCE_MANIFEST.csv`，目前沒有宣稱原始歷史source文字或commit已找回。

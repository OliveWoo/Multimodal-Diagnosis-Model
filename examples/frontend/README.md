# 完全合成的 Excel＋PPT 範例

這兩份 Office 檔是本次從零建立的虛構資料，沒有複製病人內容或截圖。它們示範現有 parser 的輸入結構，不是病原／診斷教學資料。

```powershell
python -B frontend/pipeline.py --manifest examples/frontend/manifest.synthetic.json --output local_outputs/frontend_demo_001
```

從 repo 根目錄執行。會得到 1 個檢體病例、1 頁 PPT 的影像文字紀錄與 10 類標準化 JSON；圖片數為 0，沒有模型呼叫。圖片凍結回覆的錯配拒絕另外由合成單元測試及私人 212 張實圖核對。

`synthetic_master.xlsx` 是四張 parser 必需的工作表；`synthetic.pptx` 使用可讀取的原生文字。範例中的 `900001`、`S001`、`SYN-20250102` 與 `Synthetic bacterium` 均為人造測試標記。

產物不得直接當作 OBER 最終輸入；還需研究上游的 summary、ranked、review／merge。詳細格式與實資料比較範圍見 [前端文件](../../docs/FRONTEND.md)。

## 本範例的範圍修正

這是格式工具測試，不是高醫或兩院完整cohort。兩個資料來源主線的独立合成範例位於 [cohorts](cohorts/README.md)；正式入口說明見 [前端總覽](../../docs/FRONTEND.md)。

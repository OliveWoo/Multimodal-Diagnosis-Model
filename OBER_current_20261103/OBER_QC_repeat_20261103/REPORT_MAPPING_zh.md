# 比較欄位映射修正（新API執行前）

原比較函式以organism_scope讀取已正規化記錄，但normalize實際欄位名為scope。repeat_entry.py將報告層的scope正確對應，避免兩輪都讀到None而漏計解析度變動；既有signature比較本來就包含解析度。

修正在新A/B任何QC結果產生前固定、測試並另凍結。沒有改模型、prompt、schema、request、病例、規則、門檻、原始回應或選擇。原凍結程式不覆寫。後續以repeat_entry.py pipeline／evaluate為入口。

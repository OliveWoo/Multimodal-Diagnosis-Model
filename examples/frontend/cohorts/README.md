# 高醫／兩院的獨立合成範例

全部JSON由人工建立，沒有病人資料。兩個manifest共用人工source index，但各有獨立輸入與輸出namespace：KH臨床P3對global P8，兩院P54保持global P54。這是為了實際驗證不同ID規則，不是研究病例編號。

從repo根目錄執行 `python -B frontend/kh/run.py --manifest examples/frontend/cohorts/KH.synthetic.json --output local_outputs/kh_demo_001`，或 `python -B frontend/two_hospital/run.py --manifest examples/frontend/cohorts/two_hospital.synthetic.json --output local_outputs/two_demo_001`。

入口從保存JSON整理patient_info並驗證identity，不會呼叫模型。真實歷史內容對照另在私人33例／48例驗證中完成。完整規格見 [前端總覽](../../../docs/FRONTEND.md)。

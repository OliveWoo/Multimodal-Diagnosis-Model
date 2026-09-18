# 前端先按資料來源分成兩條主線

| 主線 | 專用入口 | 編號規則 |
|---|---|---|
| [高醫 KH](kh/README.md) | `frontend/kh/run.py` | 臨床 Pn 對應 Knnn；global master ordinal 必須另外對照，不能直接使用同一個 n。 |
| [兩院：北榮、三總](two_hospital/README.md) | `frontend/two_hospital/run.py` | 保留三院 master index 的 global ID；只接受 T／V case，不從 P1 重新編號。 |

兩個入口目前負責 **既有臨床／mNGS 輸入的分線整理、identity 檢查和上游保存內容比對**。它們不重新生成 LLM 回覆、不冒充最原始 Excel 到 OBER 的完整重算。

```powershell
python -B frontend/kh/run.py --manifest examples/frontend/cohorts/KH.synthetic.json --output local_outputs/kh_demo_001
python -B frontend/two_hospital/run.py --manifest examples/frontend/cohorts/two_hospital.synthetic.json --output local_outputs/two_demo_001
```

兩份示範都是完全合成 JSON。高醫示範刻意讓臨床 P3 對應 global P8，以驗證程式不做相同流水號連接；兩院保留 P54。它們不來自真實病例。

`pipeline.py` 是共用的 **格式工具**：分人 Excel、master＋PPT、Table 1。107筆三院master的格式重播，不代表高醫與兩院兩條主線都已完成歷史追溯。Table 1 是補充格式，未確認來源身分時不歸入任何一院。

流程、分母、版本差異及尚缺材料見 [前端總覽](../docs/FRONTEND.md)；低階格式操作見 [格式工具](../docs/FORMAT_PROCESSORS.md)。

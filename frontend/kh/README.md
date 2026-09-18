# 高醫（KH）最前端

本線的臨床病人編號是 `kmuh_case_code`：`Pn` 對應 `Knnn`。mNGS 可以來自共用三院master，但必須依 case code＋檢體代碼找到正確 global ID，不能拿同名 `NGS_patient_n` 直接配。

## 流程與選用依據

1. 分人多分頁 Excel 用 `reference_versions/shared_upstream/core/raw_parser.py` 和 `parsers/` 拆出8類Raw；`core/llm_parser.py` 與 shared `prompts/` 組裝 formatter 請求。新格式入口 `frontend/pipeline.py --manifest ...` 不呼叫模型。
2. 現存高醫交付來源是成員資料根下 `0728_KH_filtter_data/KH_rawData_2Days`。33個病例的8類臨床JSON及 `all_RK_NTC_microbes`，與 `outputs/patient_info_KH_0728_2Days` 共 **297/297 分項內容相同**。
3. `frontend/kh/run.py` 固定上述来源的私人 manifest，驗證病例／檢體與檔案hash，再複製到新的 `KH/patient_info/`。它保留KH臨床編號，另存global映射。
4. 後續接 cutoff、clinical agents／summary、KH chosen ranked mNGS、scorer／review／merge，再進OBER。當前55例交付中的KH為30例；33個來源病例不等於30例交付或論文分母。

0726的2日版本有17例culture、4例FilmArray、2例GM不同；0728的3日版本亦有不同時間範圍。選用0728的2日資料是依實際上游內容吻合，沒有用日期較新當作選版理由。

## 執行入口

```powershell
python -B frontend/kh/run.py --manifest examples/frontend/cohorts/KH.synthetic.json --output local_outputs/kh_demo_001
```

實資料使用私人manifest。必需欄位包括 `cohort: "KH"`、`clinical_id_namespace: "kmuh_case_code"`、固定source index、各例的臨床ID、global ID、K case code、檢體代碼與9個來源JSON參照。完整schema以合成manifest為例；來源與比較檔均須有SHA256。

## 明確的缺口

- 先前的33份Excel是33個路徑、25組不同內容、27個檔名病例標籤；和本線33個來源病例不能混為一談。
- 以檔名病例標籤核對，現有解壓縮檔案中只有24/33有Excel候選，9/33未找到同標籤候選。另查已建立的壓縮檔清單共15個Excel成員，也沒有這9個標籤的檔名命中；不排除另名檔案。24份標籤能對上也不等於已證明當年每次formatter輸入與回覆。
- Raw→formatter→0728資料的完整歷史鏈尚未重建；新入口從明確保存的臨床／mNGS JSON起算。
- 兩個K case code在master出現不只一列；本次以保存mNGS中的檢體代碼排除歧義，33/33都找到唯一列，未使用第一筆同名候選。

輸出與來源資料不得上傳公開repo；公開資料夾只有程式、文件和合成範例。

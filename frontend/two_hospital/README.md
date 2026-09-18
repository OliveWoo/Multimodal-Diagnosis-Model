# 兩院（北榮、三總）最前端

本線使用 `legacy_global`：保留共用master index的病例ID，只接受T／V case。篩選成兩院後不能重新從P1編號，也不能把K病例納入。

## 流程與選用依據

1. 共用三院master Excel經 `split_excel_master_by_patient.py` 拆分；shared PPT parser產生 `parsed.json`，再用linker與lab importer補影像、CBC、other_lab。
2. `export_legacy_patient_json.py` 將完整master輸出成global編號。107個檢體病例包含K53、T19、V35；T＋V的54筆只是master候選範圍。
3. 成員交付 `兩院資料_0611` 實際有48個病例，附README明示剔除K。本次逐例核對global index；48例的9類臨床資料與保存上游完全相同。master重播選取同48例，也得到 **432/432份臨床JSON相同**。
4. mNGS必須另外保留版本：交付的 `mNGS_grouped` 與上游48/48份相同，但與簡單master抽取的grouped有42份不同；另6份皆為空。不能拿master抽取替換研究使用的grouped。
5. `frontend/two_hospital/run.py` 固定48例的9類臨床＋grouped來源，整理到新的 `two_hospital/patient_info/`。之後接兩院自己的ranked／summary／scorer／review／merge，再進OBER。

來源檔和目錄的` (1)`只在明確的單檔hash指定後由新輸出名稱去除；原檔保留。沒有使用舊 `prepare_patient_batch.py --overwrite`，避免碰到其中的目錄刪除與同名覆寫行為。

## 執行入口

```powershell
python -B frontend/two_hospital/run.py --manifest examples/frontend/cohorts/two_hospital.synthetic.json --output local_outputs/two_demo_001
```

實資料manifest必須寫 `cohort: "two_hospital"`、`clinical_id_namespace: "legacy_global"`、固定source index、每例T／V case與檢體，以及10個來源JSON參照。輸出保持原global ID；沒有重排病人編號。

## 不可混用的範圍與格式

- 54：master中的T／V候選；48：交付兩院來源；27：保存的rich/moderate root；25：目前55例OBER交付中的兩院；11：41例Direct-Raw比較中的兩院A組。這些集合各有明確清單。
- 25例交付的實際summary來源是rich/moderate root的20例加完整root的5例，**不能認為25例全部來自27例root**。
- 兩院來源中的 `all_RK_NTC_microbes` 是flat候選格式，上游同名檔多為grouped specimen格式，47/48份內容不同。入口不把來源的同名檔直接當作上游grouped；原flat資料仍留私人原始材料供追查。
- 6例grouped為空，入口只驗證其global index身分並明確標記 `empty_mngs_anchor_identity_from_index_only`，不捏造病原或稱為已驗證mNGS檢體內容。
- 54→48的完整歷史排除理由，以及從更原始實驗室資料生成42份較完整grouped的recipe，尚缺足夠證據。48→27／25／11只能依保存清單，不可依答案自行補選。

目前完成的是兩院臨床前端內容對照及保存輸入整理，不是重新產生全部實驗室mNGS、模型回覆與OBER結果。

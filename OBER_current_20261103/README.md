# OBER：目前保留的文獻反思＋QC版本

版本日期：2026-11-03。這是**程式碼、提示詞、合成測試與研究方法快照**，不是病人資料包，也不是已完成獨立臨床驗證的診斷產品。

OBER（Organism-Based Evidence-guided Rescue）保留上游正式選擇，逐一反思未選Context病原是否應補入。LLM抽取及審查來源，程式依固定條件決定補入；不直接要求LLM重新排名全部病原。

## 先看哪些檔案

- [WORKFLOW_zh.md](WORKFLOW_zh.md)：目前整體流程與各模組。
- [CURRENT_STATUS.md](CURRENT_STATUS.md)：目前保留的基準、分數與未採用實驗。
- [PRIVATE_DATA.md](PRIVATE_DATA.md)：私有資料界線及重播限制。
- [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json)：逐檔SHA256；來源程式及提示詞是本地研究版本的原樣複本。

## 三個版本不要混淆

|目錄|角色|目前是否採用|
|---|---|---|
|`OBER_rule_QC_20261103`|規則分流＋每篇兩次文獻QC；主設定`QC_CONSENSUS_F0_match_k2`|目前保留的研究基準|
|`OBER_QC_repeat_20261103`|兩輪新QC呼叫的驗證程式；入口是`repeat_entry.py`|檢查模組變動，不另立新診斷規則|
|`OBER_QC_core_split_20261103`|核心負向結論與細節QC分離的保守離線試驗|**不是預設版本**|

其他日期目錄是目前程式仍引用的來源、檢索、病例反思、訊號、A1或基礎工具依賴，不代表每個歷史實驗都被納入主規則，也不是建議逐一從頭跑完。

## 不需病人資料的測試

Python 3.10以上。先進入本目錄；以下合成單元測試只用Python標準函式庫，不呼叫API：

```powershell
python -X utf8 -m unittest discover -s OBER_rule_QC_20261103 -p 'test_*.py'
python -X utf8 -m unittest discover -s OBER_QC_repeat_20261103 -p 'test_*.py'
python -X utf8 -m unittest discover -s OBER_QC_core_split_20261103 -p 'test_*.py'
```

此快照不包含原始`runs/`、病歷、醫院答案、逐例報告、檢驗編號、API請求／回應或憑證。不能只clone本目錄就重現49例的數字。

## 在原本私有工作區重播

以下是研究入口說明，**不是公開快照可直接執行的無資料範例**：

```powershell
# 目前規則＋QC基準；需要既有來源及前段凍結產物
python -X utf8 OBER_rule_QC_20261103/qc_runner.py pipeline
python -X utf8 OBER_rule_QC_20261103/qc_eval.py

# 兩輪QC驗證：使用含報告映射修正的入口
python -X utf8 OBER_QC_repeat_20261103/repeat_entry.py pipeline
python -X utf8 OBER_QC_repeat_20261103/repeat_entry.py evaluate

# 已保存A/B輸出的保守離線試驗；不呼叫API
python -X utf8 OBER_QC_core_split_20261103/offline_replay.py predict
python -X utf8 OBER_QC_core_split_20261103/offline_replay.py evaluate
```

完整重播須在原私有工作區、保留原凍結路徑與檔案雜湊的條件下執行。部分歷史檔案含本地路徑，原文件中的bundled Python路徑不是跨電腦安裝指引。這次發布沒有為了可攜性改寫已凍結研究程式，也未重新整理成適用所有新資料的單一端到端CLI。

新API執行使用環境中的`OPENAI_API_KEY`，不要把金鑰寫入程式或commit。不得因分數不理想重抽、覆寫舊結果或事後更換主k。要新增全新重跑應另建實驗版本，並先取得資料外送授權。

## 結果與限制

固定49例開發資料、87個病例—答案配對，原正式97筆保留。兩輪QC的主設定均為TP76、FP24、FN11，Precision **76.00%**、Recall **87.36%**、F1 **81.28%**。

這些是固定G/T與A1病例判讀後的條件式重播結果。中間文獻分類仍有變動；相同最終分數不等於沒有非確定性。非獨立test，也不能當作完整OBER重新執行、臨床正確性或部署效能的證明。

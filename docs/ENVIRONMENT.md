# 環境與操作

> 2026-09-18更新：已還原並驗證v19 scorer，保存summary＋ranked起算55/55完整內容相同；請用 `scripts/replay_v19_upstream.py`。v20 計分規則保留，共用 model-free helper 已與退役的 LLM mNGS-max 分支解耦。raw→summary不是完整重現。方法見 [版本恢復與操作](UPSTREAM_RECOVERY.md)。
## 本次已驗證的環境

- Windows，Python 3.13.7。
- Python 直接套件版本在根目錄 `requirements-audit.txt`；是目前機器的離線驗證快照，沒有冒充原論文 lockfile，也未以全新 virtualenv 安裝驗證。
- Node 版本見 `docs/environment_snapshot.json`。Node unit tests 使用內建 test runner；不需要載入 XLSX 匯出依賴。
- 每個來源原有的 `requirements.txt` 原樣保留；reference `chung_test` 的舊 OpenAI/Pydantic/FastAPI pin 與現行主線不同，應使用獨立環境，勿合併安裝。

建立環境時可從以下指令開始；本次沒有執行套件下載／安裝：

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements-audit.txt
.\.venv\Scripts\python -B scripts/check_offline.py
```

Linux/macOS 虛擬環境使用 `.venv/bin/python`。目前跨平台未實測；現有研究工具中有 Windows 路徑預設值。

## 不需 key 的例子

```powershell
python -B run_deterministic_v20.py --input-root examples/frontend/cohorts/two_hospital --output-root local_outputs/deterministic_v20_demo_001
python -B scripts/demo_offline.py --output local_outputs/demo_001.json
python -B scripts/check_offline.py --output local_outputs/check_001
```

`run_deterministic_v20.py` 是標準化病人 JSON 到 v20 可能菌種的正式離線入口。它會複製輸入到全新 run 目錄、建立 ranked mNGS 與 with-FilmArray deterministic summary，再執行 v20 scorer；不修改來源目錄。主要結果在 `final_results.json`，完整來源與產物 hash 在 `run_manifest.json`。輸出目錄若已存在會直接拒絕執行。

Node 不在 PATH 時，第二條加 `--node <node executable>`。公開套件已排除含病例資料的測試；private cohort 的 locked integration 仍明確排除。

獨立候選證據工具的合成範例（從該 package 目錄執行）：

```powershell
python -B run_pipeline.py --input-root examples/input --mngs-root examples/mngs --output-root run_output/demo_001 --anchor-source mngs --candidate-window-days 2 --window-mode date --mode both --pretty
```

此入口只作確定性資料整理；不包含 `generate_llm_rationales.py`。請用新的 output root。

## 需外部服務的流程

formatter、clinical agent、mNGS_to_specimen_agent、missed review、OBER ArticleJudge／CaseFitJudge、理由生成、Direct-Raw inference 使用 OpenAI；PubMed 檢索需要網路；私人封存的舊 KB 路徑可能需要 embedding。`.env.example` 僅列空變數，不會自動載入或觸發請求。各來源的環境載入行為不同，應透過各自 CLI/config 明確設定，不引用原 `env.env`。

本次未執行上述付費流程。不能以「重播」命名就認定離線，特別是 A1 batch 存在 fallback。

需跑完整主線時，根目錄 `run_full_pipeline.py` 是 paid-inference 編排入口；它與離線 `run_deterministic_v20.py` 分開。先執行下列 `--plan`，只驗證設定並列出每個 stage，不建立 output 或呼叫外部服務：

```powershell
python -B run_full_pipeline.py --config examples/full_pipeline_config.template.json --output-root local_outputs/full_run_001 --plan
```

移除 `--plan` 才會執行 missed review、live RAG／PubMed、casefit、兩輪 rationale、R5 與 delivery。入口拒絕既有 output root，逐步更新 `RUNNING.json`，失敗寫 `run_failed.json`，完整成功寫 `run_manifest.json`。設定與限制見 [FULL_PIPELINE.md](FULL_PIPELINE.md)。

## XLSX 與選配 API

`OBER_patient_delivery/workbooks.mjs` 使用 `@oai/artifact-tool`。原材料沒有可攜的公開安裝與授權說明，本次沒有 vendor 本機 runtime，也沒有製造不實的 package-lock。正式 repo 必須提供可取得的正確依賴，或另外做經 golden-file 驗證的匯出器替換；替換版需明確標記，不能稱為歷史原碼。

舊 FastAPI 入口屬參考分支，需獨立確認 uvicorn／multipart 等依賴，不列入主線已驗證環境。

## 新前端離線執行

前端共用根目錄 `requirements-audit.txt` 的 Python 依賴。執行 `python -B frontend/pipeline.py --manifest examples/frontend/manifest.synthetic.json --output local_outputs/frontend_demo_001` 即可讀取合成 Office 檔並轉成標準 JSON；使用者不需要 Node 或 artifact-tool 來執行前端。合成 Office 檔由 artifact-tool 建立並經輸出／畫面检查，build 環境不是前端讀取需求。

實資料請使用私人 manifest 與全新私人輸出目錄；不要直接執行會預設寫回原資料旁的舊 importer。程式一開始套用網路封鎖並移除 API key，所用研究 modules／prompts 由 `docs/FRONTEND_CODE_PROFILE.json` 固定。此防誤觸機制不是用來執行惡意程式的隔離沙箱。

## 分線保存輸入整理入口

先依資料來源使用`frontend/kh/run.py`或`frontend/two_hospital/run.py`，公開合成manifest位於`examples/frontend/cohorts/`。操作見 [前端總覽](FRONTEND.md)。兩者共用Python環境、沒有付費API呼叫；Office格式工具仍是`frontend/pipeline.py`。新分線入口使用`COHORT_CODE_PROFILE.json`，原格式入口仍使用其原`FRONTEND_CODE_PROFILE.json`，不破壞先前固定的格式重播。

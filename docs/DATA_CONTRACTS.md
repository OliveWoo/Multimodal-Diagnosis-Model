# 資料介面、銜接與缺口

> 前端已改按高醫KH與兩院分線。KH保存輸入297份、兩院保存輸入480份與對應上游相同；原107筆master重播屬格式層驗證，不代表兩條原始資料鏈均已重建。入口、身分規則與尚缺材料見 [FRONTEND.md](FRONTEND.md)。

> 2026-09-18更新：已還原並驗證v19 scorer，保存summary＋ranked起算55/55完整內容相同；請用 `scripts/replay_v19_upstream.py`。v20 計分規則保留，共用 model-free helper 已與退役的 LLM mNGS-max 分支解耦。raw→summary不是完整重現。方法見 [版本恢復與操作](UPSTREAM_RECOVERY.md)。
> 新增上游實測證據：保存v19與現況v20非完全相同；55例選菌集合相同、5例排序不同、15例候選等級不同。詳見 [上游版本核對](UPSTREAM_VERSION_STATUS.md)。下文OBER凍結重播結果不代表上游程式已核定。

| 接口 | 必需欄位／檔案 | 檢查與限制 |
|---|---|---|
| Excel/PPT → patient JSON | 各 section、patient index、specimen code、來源檔與時間 | Excel 格式不一致，需選對 parser。跨醫院的流水號不可直接對接；shared exporter 的 index 不能丟。PPT 可含影像內姓名／病歷號。 |
| Raw → normalized | `*_Raw.json`、formatter prompt、模型／API 參數 | 新前端固定 shared 分支，8 類 Raw 均可精確配對該分支提示詞並產生待模型處理請求。成員 upstream 缺 prompts，沒有把其他分支提示詞偷放回其中；當年正式模型／prompt 仍待歷史證據。 |
| clinical → summary | standard section → agent outputs 或 normalized fallback → final summary | 使用 LLM summary 或 deterministic summary 會改變 lineage；新 evidence-preservation 程式不代表歷史 artifact 原本就有 provenance。 |
| mNGS → chosen ranked | grouped／mapping、reads、分類、檢體與時間 | 從全 raw mNGS、filtered mNGS、chosen mNGS 得到的是不同候選全集。Direct-Raw baseline 必須保留其定義的 raw 範圍。 |
| 上游 → OBER | `patient_id`, `deterministic_max.pathogen_candidates`, `deterministic_max.best_available_summary.picked_pathogens`, `llm_missed_candidate_review` | input adapter 可接不同 wrapper；多個 candidate bundle 或相衝 ID 會拒絕。正式 merge 僅 high/context 作文獻補漏池；不能把 low/omitted 全部混入。 |
| frozen RAG → A1 | candidate organism、article metadata/abstract、latest merge、source hashes | 缺候選／文章會觸發網路與 LLM fallback；部分現存 A1 artifact 未內嵌 fallback 文獻快照，hash 相符不等於文獻可完整重建。 |
| A1/B/C → R5 | A1 strong/partial/match，B site／absolute／relative，C grade，picked | counts 必須對同一 candidate、同一病例；R5 source provenance 與 artifact hash 需固定。rule 與 clinical config 要一起保存。 |
| R5 → rationale → export | final decision、selected evidence、casefit、accepted rationale、delivery config | `upstream_picked` 與 `workflow_final` 兩種 config 不可互換；新版9個 rescued selection 需有相符理由與 evidence。 |
| prediction → evaluation | frozen cohort、namespace identity mapping、prediction hash、gold版本、aliases／匹配規則 | 相同 P number 不保證相同病人。分母55交付、33KH研究、41比較不可混用。答案修訂必須另立版本，不回填成原模型輸出。 |

## 已確認的缺口

1. **formatter 執行界線**：已補 shared 分支的精確 prompt 配對及264份離線請求；新入口不呼叫付費模型。歷史模型／prompt 版本與缺失回覆不能由目前配對反推。
2. **PPT 分支已隔離接通**：新入口固定 shared parser → parsed.json → linker／lab importer，已驗證200頁。成員 upstream 的文字／图片 parser 仍保留原碼，不混接；新 adapter 明確控制每份來源的影像／檢驗匯入方式。
3. **Table 1 不是全部來源**：實測單份提供的 Excel，11個輸出 JSON 中8個與保存內容相同；不同者是 underlying、admission_diagnosis、mNGS_grouped。進一步以3個已提供的 baseline 候選在記憶體中合併，其中2個可得到與保存結果完全相同的10個 clinical/mNGS section。已確認可重建的內容銜接，但缺少唯一指認歷史 baseline 路徑的執行 manifest。
4. **歷史絕對路徑**：有 `D:\CSIE_PROJECT`、舊資料集根目錄等；必須依 run config 與執行 cwd 解釋。部分 `../RAG_re/...` 是相對 component cwd，不能一律用 manifest 所在目錄判定缺檔。
5. **原始輸入追溯不足**：歷史 Raw/normalized 沒有完整 source hashes；部分資料集只交付 normalized，沒有相配 Raw／原 workbook。現行 source-provenance 修正無法反推遺失原檔。
6. **文獻可追溯性不足**：A1 fallback 沒有完整保存文章正文欄位；交付器從舊 frozen RAG 取文章，可能丟失新增文獻標題／摘要。
7. **交付欄位已知問題**：`delivery_model.mjs` 的「原上游處置」取自已更新的 final disposition，rescued 可能被誤標為原上游選中；quote 對照來源／大小寫規則亦與 A1 不一致。本階段保留原碼並標記，未把修正版冒充既有結果版本。
8. **理由與證據**：原交付 audit 記錄部分 selected 缺可定位病例 evidence、理由模型只收到 A1 counts 等問題；不能因交付格式驗證通過就稱理由均經醫師驗證。R5 input preparation 的內容變更與 input hash 更新也需正式檢查／修正。
9. **XLSX 匯出環境**：`workbooks.mjs` 依賴 `@oai/artifact-tool`，既有本機 Node runtime junction 不是可攜安裝規格；尚未提供可公開復現的套件取得方式。本次只跑 Node 資料模型測試，沒有重新匯出全部 XLSX。
10. **完整歷史環境**：原多數 requirements 未鎖版，LLM／PubMed 結果具時點差異。只有本次離線驗證環境快照，不是論文執行環境證明。

## 路徑規則

- OBER 子專案維持本 repo 的兄弟目錄；它們會透過 parent 路徑載入其他 package。
- upstream 與 shared upstream 分開在各自目錄執行，避免同名 `tools`／`core` import 混用。
- CLI 所有輸入／輸出使用明確參數，避免研究腳本的舊預設路徑。新的輸出用獨立目錄；不要在原病例目錄執行會就地更新的 importer／cleanup。
- 私有重播資料放 repo 外，或 `private_data/`（已忽略）。病例 manifest 依原程式定義解析相對路徑，不以搬動 JSON 方式強行修復。
- 完整 reproduction manifest 至少包含：case namespace、cohort、原始輸入 SHA、normalized SHA、upstream merged SHA、RAG/A1 artifact SHA、prompt/config/code SHA、model/參數、final-decision SHA、gold/alias/evaluation-policy SHA。

## 高醫／兩院前端身分契約

KH使用`kmuh_case_code`，clinical Pn必須對Knnn，並以specimen連到master global ID。兩院使用`legacy_global`，只接受T／V case且保留global ID。兩者manifest和輸出目錄分離。

KH保存來源33例包含8類clinical＋structured all_RK_NTC；兩院48例包含9類clinical＋mNGS_grouped。兩院來源的flat all_RK_NTC不能因檔名相同就接到需要grouped的上游；6例grouped為空會明確標記，未補造數據。297／480是保存內容比對，1,070是三院master格式重播，皆不是完整raw→OBER的等價證明。見 [FRONTEND.md](FRONTEND.md)。

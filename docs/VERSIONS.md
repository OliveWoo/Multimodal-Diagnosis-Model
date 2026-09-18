# 版本、結果與論文狀態

> 2026-09-18更新：已還原並驗證v19 scorer，保存summary＋ranked起算55/55完整內容相同；請用 `scripts/replay_v19_upstream.py`。v20 計分規則保留，共用 model-free helper 已與退役的 LLM mNGS-max 分支解耦。raw→summary不是完整重現。方法見 [版本恢復與操作](UPSTREAM_RECOVERY.md)。

## 判定規則

以內容 SHA256、import／行為差異、artifact manifest、source hash、prompt/config hash 與離線重播建立關係。日期與 `latest` 只作定位線索。沒有完整證據者標為「待確認」，不擅自提升成論文正式版本。

| 材料 | 本次能確認的角色 | 不能據此推論的事 |
|---|---|---|
| `upstream/` 現況原檔（不含另名還原檔） | 另一位成員交付的目前工作樹；含 evidence preservation／provenance 及多條後续開發工具 | 不等於 Git HEAD，也不等於早期研究結果的生成程式 |
| `reference_versions/shared_upstream/` | 原共同專案另一個完整來源；包含 master Excel、PPT／vision、Table 1 匯入分支 | 不能因路徑較舊就丟棄，也不能直接整套取代現行上游 |
| earlier RAG（私人盤點備份） | 舊 RAG 實作及其差異分析保留於私人材料；已從這份 OBER 交付repo移出 | 不等於 RAG_re／R5 的前後相容版本 |
| RAG_re／clinical／casefit | 本機現況程式、prompt/config；保存 artifacts 可追 source hashes | hash 只能證明檔案一致，不能證明每篇 fallback 文獻完整保存或模型重呼叫會相同 |
| 保存的55例 R5 | 55份 complete decision JSON 本次重播全部相同；114 upstream picked＋9 rescued＝123 final selected | 尚未被使用者指定為論文主結果；開發選規則已看過答案，非 untouched external test |
| Direct-Raw identity-fixed 41例 | 四模型164份保存 prediction 重算分數與 per-patient 完全相同；24份 py/txt/json run snapshot 與現況同名來源 byte-identical | 不等於最早同名41例，也不是與 OBER 同模型同成本的控制變因 ablation |
| 後续單例重跑 | 有後續兩病例單獨重跑 artifact，可獨立追溯 | 不可自動覆蓋55例整批或41例比較表 |

## 實際差異

- upstream vs shared 的指定主程式範圍：41份相同、18份不同、74份只在現行 upstream、16份只在 shared。18份差異包括 llm_parser、summary／agent、mapping、PPT parser 等，不只換檔名。
- upstream 內的 RAG 副本 vs 使用者 OBER 目錄的 RAG：35份相同、3份不同、使用者目錄多3份。差異主檔是 `rag_engine_v3.py`、`run_casecard_cli.py` 與 folder picker。使用者版保留 explicit selected/picked bundle 資訊，避免 fallback 取代已指定的選菌。
- 現行 upstream formatter prompt 目錄缺失；shared 9份中7份可在 upstream Git HEAD 找到，原始 bytes 因換行不同，正規化文字相同；另2份無該 HEAD 對應。此事只支持「可回收的 prompt 候選」，不支持「已配對論文版」。
- 現行上游9月 evidence-preservation／shadow 修正與8月正式包不同。上游文件記錄 picked／queue 未變但部分 level 改變，不能只比 picked count 就當演算法版本相同。
- Gold 有修訂：OBER／Direct-Raw 保存的修訂包含4病例變更；另一位成員的上游修訂文件列3病例。需確認來源、匹配政策與涵蓋 cohort，不能只因總 gold pair 數相同就混用。

## 論文 freeze 所需登錄

每張論文 table／figure 指向一個 `experiment_id`，附：cohort、input identity mapping、code commit或每檔SHA、prompt/configSHA、model/參數、frozen prediction、gold版本、評估腳本／alias政策、metrics JSON。推理輸出和後續 answer revision 分開保存。新增重跑使用新 experiment_id，不就地替換。

正式論文主表／run 尚待確認。本 repo 不建立代表正式論文的 tag、CITATION 作者名或 DOI；權利人與教授確認後再填。

## 另名還原版本

新增scorer v19與summary v3還原檔、快取SHA256及驗證方法見 [UPSTREAM_RECOVERY.md](UPSTREAM_RECOVERY.md)。v19已通過55例scorer重播，但還原檔不是原始source文字／commit；summary重建仍有3例院內證據差異。

## 前端版本選擇的新增證據

新增入口明確使用 shared parser／linker／exporter；歷史來源快照保留於 `reference_versions/`，現行 `upstream/` 已做 deterministic runtime 結構重整。shared master與成員Specimen副本在抽取內容上相同，成員根目錄master存在分項內容差異；沒有以檔名或日期挑選。來源層級的PPT影像／lab路由重現1,070/1,070份標準分項。Table 1 baseline只有內容相容證據，未認定唯一歷史來源；107個檢體病例也未自動升格為論文55例。詳見 [FRONTEND.md](FRONTEND.md)。

## 依高醫／兩院修正前端版本定位

高醫0728兩日資料與保存上游297/297相同，0726兩日及0728三日版本有分項差異，不能互換。兩院48例clinical＋grouped與保存上游480/480相同；簡單master中的42份nonempty grouped與此不同，來源flat all_RK_NTC也與上游同名grouped不同，均保留版本界線。新的cohort staging只驗證保存輸入，不聲稱找回所有歷史生成recipe；兩個主線目錄、identity規則與差異見 [FRONTEND.md](FRONTEND.md)。

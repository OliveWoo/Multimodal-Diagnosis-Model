# 上游版本核對狀態（2026-09-18 更新）

**已找到並還原可重播保存結果的 v19 scorer。** 新入口為 `scripts/replay_v19_upstream.py`，使用 `upstream/tools/deterministic_mngs_max_scorer_v19_recovered.py`。`deterministic_mngs_max_scorer.py` 仍是 v20；兩版計分規則分開保留，但共用 model-free helper 已抽至 `mngs_common.py`，兩者不要混用。

| 驗證 | 結果 |
|---|---|
| 還原 scorer 對原快取 | 130/130 函式及整個模組編譯邏輯一致 |
| 還原 summary 對原快取 | 50/50 函式及整個模組編譯邏輯一致 |
| 保存 summary＋ranked → 還原 v19 scorer | 55/55 完整非路徑內容相同，包含排序、等級、理由及版本 |
| v19 結果 → 保存 merge 的 deterministic 銜接 | 55/55 選菌順序、候選集合及等級相同；不是全 merge JSON 重新生成 |
| 現況 v20 → 同批輸入 | 選菌集合55/55；排序50/55；候選等級40/55相同 |
| 目前 agent 資料 → 還原 summary v3 | 29/55 完整非路徑內容相同；23例僅GM說明行不同，3例院內證據不同 |
| 重建 summary → 還原 scorer | 52/55 完整非路徑內容相同 |

本次排除的來源容器只有 `source_files`、`merged_source_files`，不忽略版本、taxonomy、排序或理由。保存的 summary＋ranked 是主重播入口；不把尚未完全重現的 summary 重建鏈當成已完成。

已建立獨立私人輸入快照、公開合成範例、固定程式／規則的 hash 清單與拒絕錯配的重播入口。原始研究檔完全保留，沒有 API 流程或 GitHub 推送。

完整方法、指令、差異原因及限制見 [UPSTREAM_RECOVERY.md](UPSTREAM_RECOVERY.md)。原先 v20 比較保留於 [UPSTREAM_REPLAY_SUMMARY.json](UPSTREAM_REPLAY_SUMMARY.json)；新增 v19 驗證彙總見 [UPSTREAM_RECOVERY_VALIDATION.json](UPSTREAM_RECOVERY_VALIDATION.json)。

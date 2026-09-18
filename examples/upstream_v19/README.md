# 合成 v19 重播範例

人工構造的病例 900001，使用實際規則涵蓋的菌名測試 M1/M2 分支；不含真實病人資料。expected scorer 由已核對的 v19 快取計算，不是獨立臨床 ground truth。

在 repo 根目錄執行：

```powershell
python -B scripts/replay_v19_upstream.py --manifest examples/upstream_v19/manifest.synthetic.json --output local_outputs/v19_demo_001
```

程式拒絕覆寫已存在的輸出。這個範例只涵蓋保存 summary＋ranked→scorer，沒有呼叫模型或產生 OBER 的文獻證據。

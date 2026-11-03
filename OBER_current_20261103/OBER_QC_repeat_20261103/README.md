# 兩個真正新QC輪次

本輪228個預定新QC呼叫；不是重用上一輪的兩個答案。每輪57篇×2次判讀，形成一組共識；共A、B兩輪。固定G/T與A1病例判斷作控制，所以此實驗只測文獻QC層，不是從原始病歷跑兩次全部OBER。完整定義見PROTOCOL_zh.md。

```powershell
$oberPython = 'C:/Users/chung/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $oberPython -X utf8 -m unittest discover -s OBER_QC_repeat_20261103 -p test_*.py
& $oberPython -X utf8 OBER_QC_repeat_20261103/repeat_driver.py prepare
& $oberPython -X utf8 OBER_QC_repeat_20261103/repeat_driver.py pipeline
& $oberPython -X utf8 OBER_QC_repeat_20261103/repeat_driver.py evaluate
```

已成功的本次呼叫可接續，但不會複製舊QC輸出來製造一致。新實驗需要另建版本目錄；不可將再次讀快取當成第三個fresh run。沿用原官方API憑證載入器，不在程式存金鑰。

輸出在`runs/repeat_v1/round_A|round_B/`，各有實際API請求／回應、文章共識、無答案選擇、49例報告及分數。common_controls/保存兩輪共用的固定A1來源；若新QC開啟原未評候選，最多補28個原版A1呼叫，這部分只做一次供兩輪共用並另記成本，不混稱為QC重複。

`COMPARISON_zh.md`及`comparison.json`／CSV比較：每輪內分歧、跨輪共識／QC變動、最終病原清單、P/R/F1與完整k曲線。相同初始request hash、互斥Response ID與時間戳用來查核未使用舊答案。平台prompt cache不等於模型輸出快取。

主設定保持QC_CONSENSUS_F0_match_k2，兩輪全部無答案選擇凍結後才計分。兩輪同分不代表分類全相同，也不能保證未來每次相同。已反覆開發的49例，非獨立test；更不等於醫學正確性驗證。

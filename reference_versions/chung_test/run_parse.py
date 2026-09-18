# run_parse.py
# 用法：
#   python run_parse.py ipf_case_01.txt ipf_case_01.text_json.json
# -----------------------------------------------

import json
import sys
from pathlib import Path
from text_parser import parse_text_to_textjson

def main():
    if len(sys.argv) < 3:
        print("使用方式：python run_parse.py <輸入文字檔> <輸出JSON>")
        print("例如：python run_parse.py ipf_case_01.txt ipf_case_01.text_json.json")
        sys.exit(1)

    src = Path(sys.argv[1])
    dst = Path(sys.argv[2])

    if not src.exists():
        print(f"❌ 找不到檔案：{src}")
        sys.exit(1)

    print(f"📖 讀取檔案：{src}")
    text = src.read_text(encoding="utf-8", errors="ignore")

    # 印出前幾行檢查文字格式
    print("\n== 前15行內容（供除錯） ==")
    for i, ln in enumerate(text.splitlines()[:15], 1):
        print(f"{i:02d} | {repr(ln)}")

    print("\n🚀 解析中，請稍候...\n")

    try:
        out = parse_text_to_textjson(text, source=str(src))
    except Exception as e:
        print("❌ 解析時發生錯誤：", e)
        raise

    # 輸出結果
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"✅ 已輸出結果至：{dst.resolve()}")

    # 顯示摘要
    lab_batches = len(out.get("sections", {}).get("laboratory", {}).get("batches", []))
    imaging_flags = out.get("sections", {}).get("imaging", {}).get("flags", [])
    print(f"\n🧪 抽取到 {lab_batches} 組實驗資料批次。")
    print(f"🩻 影像標籤：{imaging_flags}")

if __name__ == "__main__":
    main()

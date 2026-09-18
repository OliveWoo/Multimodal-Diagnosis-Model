import os, json
from typing import Tuple
from openai import OpenAI
from dotenv import load_dotenv
from pydantic import ValidationError
from schemas import CaseCard, empty_casecard
from metrics import METER

load_dotenv()
client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
MODEL = os.environ.get("CASECARD_MODEL", "gpt-4o")

SYS_PROMPT = """你是臨床資料抽取器（Clinical IE）。
請將輸入的各段文本轉為嚴格 JSON 的 Case Card：
- 產出 summary.key_labs / imaging_flags / suspected_dx
- laboratory → batches[{date, items[item,unit,value,reference,notes,source_span,confidence]}]
- pft → rows[{date, FVC_L, FVC_pct, FEV1_L, FEV1_pct, TLC_L, TLC_pct, DLCO_pct, source_span}]
- imaging 保留 text 並抽 flags（UIP, reticulation, honeycombing, traction bronchiectasis, subpleural, basal, cardiomegaly）
- pathology 保留 text
規則：
- 日期：YYYY/MM 或 YYYY/MM/DD
- titer（例 1:40 speckled）整段放 value
- 雙分析物合併行（Albumin / Total protein ... 4.4 / 8.1 3.50-5.50 / 6.4-8.3）拆成兩筆
- 註記(# 之後) → notes；reference 僅放參考區間或陰陽性
- 為每個關鍵項加上 source_span 與 confidence(0-1)
僅輸出 JSON，勿加任何說明。
"""

USER_TEMPLATE = """[Laboratory]
{LAB}

[Imaging]
{IMG}

[PFT]
{PFT}

[Pathology]
{PAT}
"""

def _pack_sections(text_json: dict) -> Tuple[str, str, str, str]:
    sec = text_json.get("sections", {})
    lab = sec.get("laboratory", {}).get("text", "")
    img = sec.get("imaging", {}).get("text", "")
    pft = sec.get("pft", {}).get("text", "")
    pat = sec.get("pathology", {}).get("text", "")
    return lab, img, pft, pat

def parse_textjson_to_casecard(text_json: dict) -> CaseCard:
    lab, img, pft, pat = _pack_sections(text_json)
    user = USER_TEMPLATE.format(LAB=lab, IMG=img, PFT=pft, PAT=pat)

    last_err = None
    last_usage = None

    for attempt in range(2):  # 0,1
        resp = client.chat.completions.create(
            model=MODEL,
            temperature=0,
            response_format={"type": "json_object"},
            messages=[
                {"role": "system", "content": SYS_PROMPT},
                {"role": "user", "content": user if not last_err else (user + f"\n\n⚠️上次驗證錯誤：{last_err}\n請修正後輸出嚴格 JSON。")}
            ],
        )
        usage = getattr(resp, "usage", None)
        last_usage = {
            "model": MODEL,
            "input_tokens": getattr(usage, "prompt_tokens", None),
            "output_tokens": getattr(usage, "completion_tokens", None),
            "total_tokens": getattr(usage, "total_tokens", None),
            "attempt": attempt + 1
        }
        # 累計用量（就算驗證失敗也要計）
        METER.add(MODEL, last_usage["input_tokens"] or 0, last_usage["output_tokens"] or 0)

        content = resp.choices[0].message.content
        try:
            data = json.loads(content)
            cc = CaseCard.model_validate(data)
            # 寫 usage 到 meta
            if not isinstance(cc.meta, dict):
                cc.meta = {}
            cc.meta.setdefault("usage", last_usage)
            return cc
        except (json.JSONDecodeError, ValidationError) as e:
            last_err = str(e)
            continue

    # 回退：回傳空殼與錯誤訊息
    cc = empty_casecard()
    cc.meta["errors"].append({"stage": "casecard_parse", "message": last_err})
    if last_usage:
        cc.meta["usage"] = last_usage
    return cc

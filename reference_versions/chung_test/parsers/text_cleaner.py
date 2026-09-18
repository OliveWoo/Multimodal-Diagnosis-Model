

# parsers/text_cleaner.py
from __future__ import annotations
import json
import os
from typing import Any, Dict

from openai import OpenAI
from schemas import TextJSON
from metrics import METER

_MODEL = os.getenv("CASECARD_MODEL", "gpt-4o-mini")

client = OpenAI()
def preclean(s: str) -> str:
    """軟清洗：統一符號、去掉隱藏字元"""
    return (
        s.replace("•", "- ")
         .replace("・", "- ")
         .replace("\u200b", "")
         .replace("\ufeff", "")
         .replace("\u00a0", " ")
    )
def parse_raw_to_textjson(raw: str, source: str = "") -> TextJSON:
    # ✅ 新增這行：先做清洗
    raw = preclean(raw)

SYSTEM_PROMPT = (
  "You are a precise clinical text parser. Extract text for four sections: "
  "laboratory, imaging, pft, pathology. Do NOT leave any field empty; if a section "
  "truly does not exist, write 'none'. Respond ONLY as JSON with keys: "
  '{"laboratory": "...", "imaging": "...", "pft": "...", "pathology": "..."}'
)

def _make_prompt(raw: str) -> str:
    return (
        "請閱讀以下臨床病例原文，抽出四個區塊的文字摘要（非結構化也可以）：\n"
        "1) laboratory（檢驗相關數據與日期）\n"
        "2) imaging（放射影像/CT 描述與結論）\n"
        "3) pft（肺功能檢查的數值/趨勢）\n"
        "4) pathology（病理/細胞學；若未施作請敘明）\n\n"
        "輸出格式必須是 JSON：\n"
        '{ "laboratory": "...", "imaging": "...", "pft": "...", "pathology": "..." }\n\n'
        "以下是病例原文：\n"
        f"{raw}"
    )

def parse_raw_to_textjson(raw: str, source: str = "") -> TextJSON:
    """
    用 LLM 做輕清洗 → 直接產出 text_json 的四個段落（皆為純文字）。
    """
    user_prompt = _make_prompt(raw)

    completion = client.chat.completions.create(
        model=_MODEL,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_prompt},
        ],
        # 要求回傳就是 JSON，減少「空字串」或格式跑掉的機率
        response_format={"type": "json_object"},
    )

    usage = getattr(completion, "usage", None)
    usage_obj: Dict[str, Any] = {}
    if usage:
        usage_obj = {
            "model": _MODEL,
            "input_tokens": getattr(usage, "prompt_tokens", None) or getattr(usage, "input_tokens", None),
            "output_tokens": getattr(usage, "completion_tokens", None) or getattr(usage, "output_tokens", None),
            "total_tokens": getattr(usage, "total_tokens", None),
        }
        # 計入全域 token 統計
        try:
            METER.add(model=_MODEL,
                      input_tokens=usage_obj.get("input_tokens") or 0,
                      output_tokens=usage_obj.get("output_tokens") or 0)
        except Exception:
            pass

    # 解析 LLM 的 JSON 內容
    content = completion.choices[0].message.content
    try:
        parsed = json.loads(content or "{}")
    except json.JSONDecodeError:
        parsed = {}

    # 保底：確保四個欄位都存在且是字串
    lab = str(parsed.get("laboratory", "") or "")
    img = str(parsed.get("imaging", "") or "")
    pft = str(parsed.get("pft", "") or "")
    pat = str(parsed.get("pathology", "") or "")

    # 組成 TextJSON（回到你的 response_model）
    return TextJSON(
        meta={"source": source or "", "usage": usage_obj},
        sections={
            "laboratory": {"text": lab},
            "imaging": {"text": img},
            "pft": {"text": pft},
            "pathology": {"text": pat},
        },
    )

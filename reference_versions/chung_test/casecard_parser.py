# parsers/casecard_parser.py
from __future__ import annotations
import os, json
from typing import Any, Dict
from openai import OpenAI
from metrics import METER
from schemas import CaseCard, CaseCardSections, CaseCardSummary

client = OpenAI()
_CASE_MODEL = os.getenv("CASECARD_MODEL_CASE", "gpt-4o-mini")

SYSTEM_PROMPT = (
    "You are a clinical information structurer. Given a text_json (sections.laboratory/imaging/pft/pathology), "
    "extract a structured CaseCard. Return JSON ONLY, matching the exact schema. "
    "If data is missing, output empty arrays/strings but DO NOT omit keys."
)

SCHEMA_SPEC = """
Return JSON exactly with this shape:

{
  "summary": {
    "key_labs": {
      "Creatinine": {"value": "2.27 mg/dL", "date": "2014/02", "confidence": 0.9}
    },
    "imaging_flags": ["reticulation","honeycombing","traction bronchiectasis","subpleural","basal","UIP"],
    "suspected_dx": ["Usual Interstitial Pneumonia (UIP)"]
  },
  "sections": {
    "laboratory": { "batches": [
      { "date": "YYYY/MM or YYYY/MM/DD", "items": [
        {"item":"Creatinine","unit":"mg/dL","value":2.27,"reference":"0.70-1.20","notes":"","source_span":"","confidence":0.9}
      ]}
    ]},
    "pft": { "rows": [
      {"date":"2017/10","FVC_L":1.82,"FVC_pct":57,"FEV1_L":1.50,"FEV1_pct":68,"TLC_L":2.91,"TLC_pct":56,"DLCO_pct":43}
    ]},
    "imaging": { "text":"free text summary", "flags":["UIP","basal","subpleural"] },
    "pathology": { "text":"free text summary" }
  },
  "meta": { "errors": [] }
}
"""

def _make_prompt_from_textjson(tx: Dict[str, Any]) -> str:
    lab = tx.get("sections", {}).get("laboratory", {}).get("text", "") or ""
    img = tx.get("sections", {}).get("imaging", {}).get("text", "") or ""
    pft = tx.get("sections", {}).get("pft", {}).get("text", "") or ""
    pat = tx.get("sections", {}).get("pathology", {}).get("text", "") or ""
    return (
        "Sections extracted from the case:\n\n"
        f"[laboratory]\n{lab}\n\n"
        f"[imaging]\n{img}\n\n"
        f"[pft]\n{pft}\n\n"
        f"[pathology]\n{pat}\n\n"
        "Build the CaseCard. " + SCHEMA_SPEC
    )

def parse_textjson_to_casecard(textjson: Dict[str, Any]) -> CaseCard:
    prompt = _make_prompt_from_textjson(textjson)

    completion = client.chat.completions.create(
        model=_CASE_MODEL,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ],
        response_format={"type": "json_object"},
        temperature=0.1,
    )

    # usage & metrics
    usage = getattr(completion, "usage", None)
    usage_obj = {}
    if usage:
        usage_obj = {
            "model": _CASE_MODEL,
            "input_tokens": getattr(usage, "prompt_tokens", None) or getattr(usage, "input_tokens", None),
            "output_tokens": getattr(usage, "completion_tokens", None) or getattr(usage, "output_tokens", None),
            "total_tokens": getattr(usage, "total_tokens", None),
        }
        try:
            METER.add(model=_CASE_MODEL,
                      input_tokens=usage_obj.get("input_tokens") or 0,
                      output_tokens=usage_obj.get("output_tokens") or 0)
        except Exception:
            pass

    # 解析 LLM JSON
    content = completion.choices[0].message.content or "{}"
    try:
        data = json.loads(content)
    except Exception:
        data = {}

    # ---------- 保底填補（避免 key 缺漏） ----------
    data.setdefault("summary", {})
    data["summary"].setdefault("key_labs", {})
    data["summary"].setdefault("imaging_flags", [])
    data["summary"].setdefault("suspected_dx", [])

    data.setdefault("sections", {})
    sec = data["sections"]
    sec.setdefault("laboratory", {"batches": []})
    sec.setdefault("pft", {"rows": []})
    sec.setdefault("imaging", {"text": "", "flags": []})
    sec.setdefault("pathology", {"text": ""})

    data.setdefault("meta", {})
    data["meta"].setdefault("errors", [])

    # 把原始文字補回 imaging/pathology（萬一 LLM 沒產）
    src = textjson.get("sections", {})
    if not sec["imaging"].get("text"):
        sec["imaging"]["text"] = (src.get("imaging", {}).get("text") or "").strip()
    if not sec["pathology"].get("text"):
        sec["pathology"]["text"] = (src.get("pathology", {}).get("text") or "").strip()
    if not sec["imaging"].get("flags"):
        sec["imaging"]["flags"] = list(data["summary"].get("imaging_flags", []))

    # meta.source & usage
    data["meta"]["source"] = (textjson.get("meta", {}) or {}).get("source", "")
    data["meta"]["usage"] = usage_obj

    # 交給 Pydantic 做驗證，型別不對會自動報錯成 422/500
    return CaseCard.model_validate(data)

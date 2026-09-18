# parsers/casecard_parser_v2.py
from __future__ import annotations
import os, json, re
from typing import Any, Dict
from openai import OpenAI
from metrics import METER
from schemas_v2 import CaseCardV2

client = OpenAI()
_CASE_MODEL = os.getenv("CASECARD_MODEL_CASE", "gpt-4o-mini")

SYSTEM_PROMPT_V2 = """
你是臨床資料結構化小幫手。請將輸入的四段文字（laboratory/imaging/pft/pathology）整理成 CaseCardV2。
規則：
- 僅輸出單一 JSON，不能有額外文字。
- 不臆測：沒有明確資訊就留空或省略。
- Outline 是「RAG 索引」：包含 suspected_dx、imaging_flags、key_labs（帶日期/信心）、pft_trend 的句子。
- Sections 保留完整內容的結構化表示：laboratory 以日期分批、pft 逐列、imaging/pathology 保留原句與 flags。
- key_labs 只放最關鍵 3–8 個項目（例如：Creatinine、ANA、Anti-SSA/Ro 等）。
- 日期採原文格式（例如 2014/12），不要自動補日。
- 僅允許以下鍵：
  { "outline":{ "case_uid","version","patient_id","encounter_date","age","sex","chief_complaint","suspected_dx","imaging_flags","key_labs","pft_trend","notes","sources" },
    "sections":{ "laboratory","pft","imaging","pathology" },
    "meta":{ "source","usage" } }
- JSON 中不要出現註解、不要多餘文字。
"""

USER_PROMPT_TEMPLATE = """以下是四段原始文本，請輸出符合 CaseCardV2 的 JSON：

[laboratory]
{LAB_TEXT}

[imaging]
{IMG_TEXT}

[pft]
{PFT_TEXT}

[pathology]
{PAT_TEXT}

請將「大綱」放在 outline，「完整內容」放在 sections；其中：
- sections.laboratory：以日期分組的 batches 陣列（每批包含 items[]，item/unit/value/reference）
- sections.pft：逐列日期與數值（FVC_L/FVC_pct/FEV1_L/FEV1_pct/TLC_L/TLC_pct/DLCO_pct）
- sections.imaging：text 為原句彙整、flags 為關鍵影像語彙（cardiomegaly、reticulation、honeycombing、traction bronchiectasis、subpleural、basal、UIP）
- sections.pathology：text 保留原句
- outline.key_labs：挑出最關鍵實驗數值並附日期（如：Creatinine、ANA、Anti-SSA/Ro）
- outline.suspected_dx：若影像像 UIP，放 "Usual Interstitial Pneumonia (UIP)"
- outline.pft_trend：一句話形容趨勢（例：FVC 與 DLCO 逐步下降）

只輸出 JSON，且僅允許上面規定的鍵。
"""

def _make_prompt_v2(tx: Dict[str, Any]) -> str:
    lab = tx.get("sections", {}).get("laboratory", {}).get("text", "") or ""
    img = tx.get("sections", {}).get("imaging", {}).get("text", "") or ""
    pft = tx.get("sections", {}).get("pft", {}).get("text", "") or ""
    pat = tx.get("sections", {}).get("pathology", {}).get("text", "") or ""
    return USER_PROMPT_TEMPLATE.format(
        LAB_TEXT=lab.strip(),
        IMG_TEXT=img.strip(),
        PFT_TEXT=pft.strip(),
        PAT_TEXT=pat.strip()
    )

def _safe_json_loads(s: str) -> Dict[str, Any]:
    # 去除 BOM/前後雜字，只保留第一個 { ... } 區塊
    s = s.strip()
    # 若模型不小心多了前置文字，抓第一個 JSON object
    m = re.search(r'\{.*\}\s*$', s, flags=re.S)
    if m:
        s = m.group(0)
    return json.loads(s)

def _normalize_sections(data: Dict[str, Any]) -> None:
    """把模型可能偏離的結構拉回既定格式。就地修改 data。"""
    sections = data.setdefault("sections", {})

    # laboratory：強制成 {"batches": [...]}
    lab = sections.get("laboratory", {})
    if isinstance(lab, list):
        sections["laboratory"] = {"batches": lab}
    elif isinstance(lab, dict) and "batches" not in lab:
        # 有些模型可能輸出 {"laboratory": [...]} 或 {"laboratory": {"rows":[...]}}
        if "rows" in lab and isinstance(lab["rows"], list):
            sections["laboratory"] = {"batches": lab["rows"]}
        else:
            sections["laboratory"] = {"batches": []}

    # pft：強制為 list
    pft = sections.get("pft", [])
    if not isinstance(pft, list):
        sections["pft"] = []

    # imaging/pathology：確保是 dict，並補預設鍵
    img = sections.get("imaging", {})
    if not isinstance(img, dict):
        img = {}
    img.setdefault("text", "")
    img.setdefault("flags", [])
    sections["imaging"] = img

    pat = sections.get("pathology", {})
    if not isinstance(pat, dict):
        pat = {}
    pat.setdefault("text", "")
    sections["pathology"] = pat

def parse_textjson_to_casecard_v2(textjson: Dict[str, Any]) -> CaseCardV2:
    prompt = _make_prompt_v2(textjson)
    completion = client.chat.completions.create(
        model=_CASE_MODEL,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT_V2},
            {"role": "user", "content": prompt},
        ],
        response_format={"type": "json_object"},
        temperature=0.2,
    )

    # token usage
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
            METER.add(
                model=_CASE_MODEL,
                input_tokens=usage_obj.get("input_tokens") or 0,
                output_tokens=usage_obj.get("output_tokens") or 0
            )
        except Exception:
            pass

    raw = completion.choices[0].message.content or "{}"
    data = _safe_json_loads(raw)

    # 補 meta
    data.setdefault("meta", {})
    data["meta"]["source"] = (textjson.get("meta", {}) or {}).get("source", "")
    data["meta"]["usage"] = usage_obj

    # 正規化 sections 結構
    _normalize_sections(data)

    # 強制版本與來源
    data.setdefault("outline", {})
    data["outline"].setdefault("version", 2)
    if "sources" not in data["outline"]:
        data["outline"]["sources"] = [data["meta"]["source"]] if data["meta"].get("source") else []

    # 驗證＆回傳
    return CaseCardV2.model_validate(data)

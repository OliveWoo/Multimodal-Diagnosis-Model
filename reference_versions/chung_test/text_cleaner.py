import os, re, json
from openai import OpenAI
from pydantic import ValidationError
from dotenv import load_dotenv
from schemas import TextJSON, empty_textjson
from metrics import METER

load_dotenv()
client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
MODEL = os.environ.get("CLEANER_MODEL", "gpt-4o-mini")

SYS_PROMPT = """你是一個臨床文本清洗器（Clinical text cleaner）。
任務：
1) 將輸入原文切成四段：Laboratory、Imaging、PFT、Pathology。
2) 去除頁碼、圖說(Figure)、頁尾註解、連續空白與奇怪符號（保留醫療內容）。
3) 僅輸出 JSON：{"meta":{},"sections":{"laboratory":{"text":""},"imaging":{"text":""},"pft":{"text":""},"pathology":{"text":""}}}
4) 不要虛構（No hallucination）。若缺某段，該段 text 留空字串。
"""

USER_TEMPLATE = """原始文本如下（可能沒有換行或被擠在一行）：
請輸出嚴格 JSON（僅 JSON，不要多餘文字），並確保四段 keys 都存在。
"""

# ========== 預清洗函式 ==========
def _pre_normalize(s: str) -> str:
    s = s.replace("\u00A0", " ").replace("\u200B", "").replace("\uFF0F", "/").replace("\u3000", " ")
    s = re.sub(r"[ \t]+", " ", s)
    # 提示模型斷行 cue：在日期與常見標題前插入換行
    s = re.sub(r"(?=\d{4}[./\-]\d{2}(?:[./\-]\d{2})?)", "\n", s)
    s = re.sub(r"(?i)\s+(?=Radiological imaging studies|Serial pulmonary function tests|Cytologic and pathologic studies)\b", "\n", s)
    return s.strip()

# ========== 主函式 ==========
def parse_raw_to_textjson(raw_text: str, source: str = "") -> TextJSON:
    raw = _pre_normalize(raw_text)

    resp = client.chat.completions.create(
        model=MODEL,
        temperature=0,
        response_format={"type": "json_object"},
        messages=[
            {"role": "system", "content": SYS_PROMPT},
            {"role": "user", "content": USER_TEMPLATE.format(RAW=raw)}
        ],
    )

    content = resp.choices[0].message.content
    usage = getattr(resp, "usage", None)
    usage_dict = {
        "model": MODEL,
        "input_tokens": getattr(usage, "prompt_tokens", None),
        "output_tokens": getattr(usage, "completion_tokens", None),
        "total_tokens": getattr(usage, "total_tokens", None),
    }

    data = json.loads(content)
    # 保底：補齊四段 key
    tj = empty_textjson(source=source).model_dump()
    tj["sections"].update(data.get("sections", {}))
    tj["meta"]["usage"] = usage_dict

    # 全域累加
    METER.add(MODEL, usage_dict["input_tokens"] or 0, usage_dict["output_tokens"] or 0)

    try:
        return TextJSON.model_validate(tj)
    except ValidationError:
        fallback = empty_textjson(source=source)
        fallback.sections["laboratory"]["text"] = raw
        fallback.meta["usage"] = usage_dict
        return fallback

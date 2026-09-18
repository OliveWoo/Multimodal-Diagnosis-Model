# text_parser.py（把以下段落覆蓋到你的檔案：helpers + regex 定義 + _parse_lab_line + _parse_labs）

import re
from typing import List, Dict, Any
# --- 清理 reference：去除腳註井號與多餘語句 ---
def _clean_ref(s: str | None) -> str | None:
    if not s:
        return None
    # 去掉前導井號與冒號空白
    s = s.strip()
    s = s.lstrip("#").strip()
    # 若出現「# :」或「 #:」等腳註起點，截斷其後雜訊
    s = re.split(r"\s+#\s*:|\s+#:|\s+#\s", s, maxsplit=1)[0].strip()
    return s or None

# --- 專門處理「Albumin / Total protein (unit) valA / valB refA / refB」一行 ---
_dual_line = re.compile(
    r"^\s*(Albumin)\s*/\s*(Total protein)\s*\(\s*([^)]+)\s*\)\s*"
    r"([0-9]+(?:\.[0-9]+)?)\s*/\s*([0-9]+(?:\.[0-9]+)?)\s+"
    r"([0-9.]+(?:\s*[-–]\s*[0-9.]+)?)\s*/\s*([0-9.]+(?:\s*[-–]\s*[0-9.]+)?)\s*$",
    re.I
)

def _expand_dual_analyte(line: str):
    """
    若符合 Albumin / Total protein 這種合併一行的格式，回傳兩筆 item；
    否則回傳 None。
    """
    m = _dual_line.match(line)
    if not m:
        return None
    a1, a2, unit, v1, v2, r1, r2 = m.groups()
    return [
        {"item": a1, "unit": unit, "value": _to_num(v1), "reference": _clean_ref(r1)},
        {"item": a2, "unit": unit, "value": _to_num(v2), "reference": _clean_ref(r2)},
    ]

# ---- 文字正規化：處理 NBSP、全形符號、零寬空白等 ----
def _normalize_text(s: str) -> str:
    import re
    # 常見怪字元修正
    s = s.replace("\u00A0", " ")   # NBSP
    s = s.replace("\u200B", "")    # zero-width space
    s = s.replace("\uFF0F", "/")   # 全形／ -> /
    s = s.replace("\u3000", " ")   # 全形空白 -> space
    s = re.sub(r"[ \t]+", " ", s)  # 壓縮多重空白
    s = re.sub(r"[ \t]+$", "", s, flags=re.M)

    # 若幾乎都是一行，啟用「自動斷行」
    if len(s.splitlines()) <= 2:
        # 1) 在每個日期前面斷行
        s = re.sub(r"(?=\d{4}[./\-]\d{2}(?:[./\-]\d{2})?)", "\n", s)

        # 2) 在各大段標題前斷行
        for hdr in [
            r"Radiological imaging studies",
            r"Serial pulmonary function tests",
            r"Cytologic and pathologic studies",
            r"Laboratory panels", r"Laboratory", r"Labs"
        ]:
            s = re.sub(fr"(?i)(?={hdr})", "\n", s)

        # 3) 在常見檢驗名出現處斷行（避免同一日期的多項目黏在一行）
        analytes = [
            r"Neutrophils", r"Eosinophils", r"Basophils", r"Monocytes", r"Lymphocytes",
            r"Hemoglobin", r"Platelets", r"Creatinine",
            r"Alanine aminotransferase", r"\[ALT\]", r"ALT",
            r"Albumin", r"Total protein",
            r"Creatinine kinase", r"\[CK\]", r"CK",
            r"Glycated hemoglobin", r"\[HbA1c\]", r"HbA1c",
            r"Anti\-nuclear antibodies", r"\[ANA\]", r"ANA",
            r"Rheumatoid factor", r"\[RF\]", r"RF",
            r"Complement 3", r"\[C3\]", r"C3",
            r"Complement 4", r"\[C4\]", r"C4",
            r"Anti\-SSA/Ro antibodies", r"Anti\-SSB/La antibodies",
            r"Anti\-Scl\-70 antibodies", r"Anti\-Jo\-1 antibodies",
            r"Hepatitis C antibody", r"Hepatitis B surface antigen"
        ]
        pat = re.compile(r"(?<!\n)\s+(?=(" + r"|".join(analytes) + r")\b)")
        s = pat.sub("\n", s)

    return s

# ---- 段落標題（可擴充同義詞）----
SEC_PATTERNS = [
    ("laboratory", re.compile(r"^\s*(Laboratory panels|Laboratory|Labs)\b.*?$", re.I|re.M)),
    ("imaging",    re.compile(r"^\s*(Radiological imaging studies|Imaging|Chest CT)\b.*?$", re.I|re.M)),
    ("pathology",  re.compile(r"^\s*(Cytologic and pathologic studies|Pathology)\b.*?$", re.I|re.M)),
    ("pft",        re.compile(r"^\s*(Serial pulmonary function tests|PFT)\b.*?$", re.I|re.M)),
]

# ---- 日期（放寬：允許 / - .）----
DATE_LINE = re.compile(r"^\s*(\d{4}[./\-]\d{2}(?:[./\-]\d{2})?)\b")

# ---- 檢驗行 regex ----
# 例：WBC ( 10^3/uL ) 9.8 3.4-9.1   （參考值最後一段變成可選）
LAB_WITH_UNIT = re.compile(
    r"^\s*([A-Za-z0-9\-\+\.,/ \[\]]+?)\s*\(\s*([^)]+?)\s*\)\s*([<>]=?\s*[\d.:]+|\d+(?:\.\d+)?)\s*([^\s].*?)?\s*$"
)
# 例：Anti-nuclear antibodies, [ANA] 1:40 speckled < 1:40
LAB_LOOSE = re.compile(
    r"^\s*([A-Za-z0-9\-\+\.,/ \[\]]+?)\s+("
    r"\d+:\d+(?:\s+[A-Za-z]+)?|"                # ← titer: 1:40 / 1:40 speckled
    r"[<>]=?\s*[\d.:/]+|negative|positive|"
    r"Normal|Abnormal|normal|abnormal|"
    r"\d+(?:\.\d+)?"
    r")\s*(.*)$"
)


IMAGING_FLAGS = [
    "UIP","reticulation","honeycombing","traction bronchiectasis",
    "subpleural","basal","cardiomegaly","pleural effusion",
    "consolidation","ground-glass","atelectasis","bronchiectasis"
]

def _split_sections(text: str) -> Dict[str, str]:
    blocks = {}
    indices = []
    for key, pat in SEC_PATTERNS:
        m = pat.search(text)
        if m:
            indices.append((m.start(), key))
    indices.sort()
    for i, (pos, key) in enumerate(indices):
        end = indices[i+1][0] if i+1 < len(indices) else len(text)
        blocks[key] = text[pos:end].strip()
    return blocks

def _to_num(s: str):
    ss = s.strip()
    # 含比較符號/比例/比值就保持字串
    if ss.startswith((">","<","=")) or ":" in ss or "/" in ss:
        return ss
    try:
        return float(ss) if "." in ss else int(ss)
    except Exception:
        return ss

# ---- 新增：單行解析器（安全、避免 m1/m2 未定義）----
def _parse_lab_line(line: str):
    # 先判斷是否為雙分析物合併一行
    dual = _expand_dual_analyte(line)
    if dual:
        return dual  # 注意：這裡回傳「list」，待上層展開

    m = LAB_WITH_UNIT.match(line)
    if m:
        item, unit, value, ref = m.groups()
        ref = _clean_ref(ref)
        return {"item": item.strip(), "unit": unit.strip(), "value": _to_num(value), "reference": ref}

    m = LAB_LOOSE.match(line)
    if m:
        item, value, ref = m.groups()
        ref = _clean_ref(ref)
        return {"item": item.strip(), "unit": None, "value": value.strip(), "reference": ref}

    return None


# ---- 修正版：_parse_labs（使用 _parse_lab_line，無 m1/m2 區域變數風險）----
def _parse_labs(block: str) -> Dict[str, Any]:
    lines = [l for l in block.splitlines()]
    # 去掉首行標題
    if lines and lines[0].strip():
        lines = lines[1:]
    batches = []
    cur_date = None
    cur_items: List[Dict[str, Any]] = []

    def flush():
        nonlocal cur_date, cur_items, batches
        if cur_date and cur_items:
            batches.append({"date": cur_date, "items": cur_items})
        cur_date, cur_items = None, []

    for raw_ln in lines:
        ln = raw_ln.strip()
        if not ln:
            continue

        # 新日期一行
        md = DATE_LINE.match(ln)
        if md:
            flush()
            cur_date = md.group(1)
            tail = ln[md.end():].strip()
            if tail:
                rec = _parse_lab_line(tail)
                if rec:
                    cur_items.append(rec)
            continue

        # 一般 lab 行
        rec = _parse_lab_line(ln)
        if rec:
            cur_items.append(rec)
            continue

        # 兜不上的行 → 試著接到上一個 item's reference（處理換行的參考值/單位）
        if cur_items:
            last = cur_items[-1]
            if last.get("reference"):
                last["reference"] = (str(last["reference"]) + " " + ln).strip()
            else:
                last["reference"] = ln

    flush()
    return {"batches": batches}

# 入口：記得先正規化文字！
def parse_text_to_textjson(text: str, source: str|None=None) -> Dict[str, Any]:
    text = _normalize_text(text)
    blocks = _split_sections(text)
    out: Dict[str, Any] = {"meta":{"source": source or ""},"sections":{}}
    if "laboratory" in blocks:
        out["sections"]["laboratory"] = _parse_labs(blocks["laboratory"])
    if "imaging" in blocks:
        # 省略：你原本的 _parse_imaging 可沿用
        pass
    if "pft" in blocks:
        # 省略：你原本的 _parse_pft 可沿用
        pass
    if "pathology" in blocks:
        out["sections"]["pathology"] = {"text": blocks["pathology"].strip()}
    return out

# summarize_level_changes.py
# 用途：
# 讀取 max_agent.json（infection bundle）+ out.json（RAG 輸出）
# 1) 找出哪些菌種 Level 有調整（升/降/維持）
# 2) 依「Best Available Summary」規則輸出：
#    - 若存在 Level 1/2：只輸出 Level 1–2 病原（1–2 個）
#    - 若最高僅 Level 3：輸出 1 個 + 明確標示證據不足與替代解釋
#    - 若全部 Level 4–5：輸出 1 個 + 明確說明目前無可信致病源/不得單獨投藥
# 3) 依指定格式輸出（上半部：1/2/3 只列菌名；下半部：補種類/主次/簡潔原因）

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


LEVEL_ORDER = {"Level 1": 1, "Level 2": 2, "Level 3": 3, "Level 4": 4, "Level 5": 5}


def load_json(p: Path) -> Dict[str, Any]:
    with p.open("r", encoding="utf-8") as f:
        obj = json.load(f)
    if not isinstance(obj, dict):
        raise ValueError(f"{p} is not a JSON object")
    return obj


def as_level_str(x: Any) -> Optional[str]:
    if x is None:
        return None
    s = str(x).strip()
    if not s:
        return None
    m = re.search(r"(Level\s*[1-5])", s, flags=re.IGNORECASE)
    if not m:
        return None
    lv = m.group(1)
    lv = re.sub(r"\s+", " ", lv, flags=re.IGNORECASE).strip()
    lv = lv[0].upper() + lv[1:]  # "level 3" -> "Level 3"
    return lv


def level_num(level_str: Optional[str]) -> Optional[int]:
    if not level_str:
        return None
    ls = as_level_str(level_str)
    if not ls:
        return None
    return LEVEL_ORDER.get(ls)


def first_sentence_cn(text: str) -> str:
    t = (text or "").strip()
    if not t:
        return ""
    for sep in ["。", "；", "\n"]:
        if sep in t:
            out = t.split(sep, 1)[0].strip()
            return out + ("。" if sep == "。" and not out.endswith("。") else "")
    return t[:120].strip() + ("…" if len(t) > 120 else "")


def extract_dominant_source(bundle: Dict[str, Any]) -> str:
    dom = bundle.get("dominant_source")
    if isinstance(dom, str) and dom.strip():
        return dom.strip()
    return "未提供（Unknown）"


def get_bundle_candidates(bundle: Dict[str, Any]) -> List[Dict[str, Any]]:
    cand = bundle.get("pathogen_candidates", [])
    return cand if isinstance(cand, list) else []


def extract_rag_level_reassessment(out_obj: Dict[str, Any]) -> List[Dict[str, Any]]:
    rag = out_obj.get("rag_answer", out_obj)
    if not isinstance(rag, dict):
        return []
    for key in ["level_reassessment_conclusion", "level_reassessment_conclusions"]:
        v = rag.get(key)
        if v is None:
            continue
        if isinstance(v, list):
            return [x for x in v if isinstance(x, dict)]
        if isinstance(v, dict):
            return [v]
    return []


def normalize_decision(s: Any) -> str:
    t = str(s or "").strip()
    if "升" in t:
        return "升級"
    if "降" in t:
        return "降級"
    if "維" in t:
        return "維持"
    return t or "維持"


def decision_rank(decision: Optional[str]) -> int:
    d = normalize_decision(decision)
    if d == "升級":
        return 0
    if d == "維持":
        return 1
    return 2  # 降級


def primary_or_secondary(recommended_level: Optional[str], decision: Optional[str]) -> str:
    n = level_num(recommended_level)
    if n is not None and n <= 2:
        return "主要"
    return "次要"


def rank_key(
    *,
    decision: str,
    orig_level: Optional[str],
    rec_level: Optional[str],
    bundle_conf: float,
) -> Tuple[int, int, int, float]:
    """
    一般排序（用於選候選時做 tie-break）：
    1) recommended_level 越低越優先（1 最重要）
    2) 升級 > 維持 > 降級
    3) original_level 越低越優先
    4) bundle_conf 越高越優先
    """
    rec_n = level_num(rec_level) or 9
    orig_n = level_num(orig_level) or 9
    return (rec_n, decision_rank(decision), orig_n, -bundle_conf)


def best_available_selection(rows: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], str]:
    """
    依新規則選 Best Available Summary（BAS）
    回傳： (selected_rows, bas_mode)
      bas_mode ∈ {"L12", "L3", "L45", "EMPTY"}
    """
    if not rows:
        return ([], "EMPTY")

    # 先找出各 level 群
    l12 = [r for r in rows if level_num(r.get("recommended_level")) in (1, 2)]
    l3 = [r for r in rows if level_num(r.get("recommended_level")) == 3]
    l45 = [r for r in rows if level_num(r.get("recommended_level")) in (4, 5)]

    if l12:
        l12.sort(key=lambda r: rank_key(
            decision=r.get("decision", "維持"),
            orig_level=r.get("original_level"),
            rec_level=r.get("recommended_level"),
            bundle_conf=float(r.get("bundle_conf") or 0.0),
        ))
        return (l12[:3], "L12")

    if l3:
        l3.sort(key=lambda r: rank_key(
            decision=r.get("decision", "維持"),
            orig_level=r.get("original_level"),
            rec_level=r.get("recommended_level"),
            bundle_conf=float(r.get("bundle_conf") or 0.0),
        ))
        return ([l3[0]], "L3")

    if l45:
        l45.sort(key=lambda r: rank_key(
            decision=r.get("decision", "維持"),
            orig_level=r.get("original_level"),
            rec_level=r.get("recommended_level"),
            bundle_conf=float(r.get("bundle_conf") or 0.0),
        ))
        return ([l45[0]], "L45")

    return ([], "EMPTY")


def add_bas_disclaimer(reason: str, bas_mode: str) -> str:
    """
    依 BAS 規則補上必要的「證據不足/替代解釋/不得單獨投藥」提示
    """
    base = (reason or "").strip()

    if bas_mode == "L3":
        extra = "【證據不足提醒】目前最高僅 Level 3，可能替代解釋包含：定殖／再活化／污染；不建議僅憑此單一證據擴大投藥，需以無菌部位/跨模組/連續趨勢再驗證。"
        return (base + " " + extra).strip() if base else extra

    if bas_mode == "L45":
        extra = "【重要提醒】目前僅 Level 4–5，代表「目前無可信致病源」。此結果不得作為單獨投藥依據，應優先追求更高品質證據（無菌部位培養/標的 PCR/抗原、影像對位、重複檢體）。"
        return (base + " " + extra).strip() if base else extra

    return base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max_agent", required=True, help="Path to max_agent infection bundle JSON")
    ap.add_argument("--out", required=True, help="Path to out.json (RAG output JSON)")
    ap.add_argument("--patient_id", default="", help="PatientID to print (optional)")
    ap.add_argument("--save", default="", help="Optional path to save output text")
    args = ap.parse_args()

    max_agent_path = Path(args.max_agent)
    out_path = Path(args.out)

    bundle = load_json(max_agent_path)
    out_obj = load_json(out_path)

    dom_source = extract_dominant_source(bundle)
    candidates = get_bundle_candidates(bundle)

    # bundle 端 organism -> info
    bundle_map: Dict[str, Dict[str, Any]] = {}
    for c in candidates:
        if not isinstance(c, dict):
            continue
        name = str(c.get("organism_name") or "").strip()
        if not name:
            continue
        bundle_map[name] = {
            "classification": str(c.get("classification") or "").strip(),
            "original_level": as_level_str(c.get("integrated_causative_level")),
            "bundle_conf": float(c.get("integrated_causative_confidence") or 0.0),
        }

    reassess = extract_rag_level_reassessment(out_obj)

    # rows：整合 RAG reassess + bundle info
    rows: List[Dict[str, Any]] = []
    if reassess:
        for r in reassess:
            pathogen = str(r.get("pathogen") or "").strip()
            if not pathogen:
                continue

            binfo = bundle_map.get(pathogen, {})
            classification = str(binfo.get("classification") or "").strip()
            bundle_conf = float(binfo.get("bundle_conf") or 0.0)

            original_level = as_level_str(r.get("original_level")) or binfo.get("original_level")
            recommended_level = as_level_str(r.get("recommended_level")) or original_level
            decision = normalize_decision(r.get("decision"))

            rationale = str(r.get("rationale") or "").strip()
            refs = r.get("supporting_references")
            ref_list: List[str] = []
            if isinstance(refs, list):
                ref_list = [str(x).strip() for x in refs if str(x).strip()]
            elif isinstance(refs, str) and refs.strip():
                ref_list = [refs.strip()]

            concise = first_sentence_cn(rationale)
            if ref_list:
                concise = (concise + " " if concise else "") + f"（{ref_list[0]}）"

            rows.append(
                {
                    "pathogen": pathogen,
                    "classification": classification,
                    "original_level": original_level,
                    "recommended_level": recommended_level,
                    "decision": decision,
                    "bundle_conf": bundle_conf,
                    "reason": concise or "未提供調整理由。",
                }
            )
    else:
        # 若 out.json 沒有 level_reassessment_conclusion：退回 bundle
        for name, info in bundle_map.items():
            rows.append(
                {
                    "pathogen": name,
                    "classification": info.get("classification", ""),
                    "original_level": info.get("original_level"),
                    "recommended_level": info.get("original_level"),
                    "decision": "維持",
                    "bundle_conf": float(info.get("bundle_conf") or 0.0),
                    "reason": "未提供 level_reassessment_conclusion，暫以 bundle 綜合信心作為 BestAvailable 候選。",
                }
            )

    # ✅ 依新 BAS 規則選出輸出名單
    selected, bas_mode = best_available_selection(rows)

    # 組輸出
    lines: List[str] = []
    pid = str(args.patient_id or "").strip()
    lines.append(f"PatientID: {pid}" if pid else "PatientID:")
    lines.append(f"感染源: {dom_source}")

    # 上半部：只列菌名（依你最新要求：1.(菌名)）
    if selected:
        for i, r in enumerate(selected, 1):
            lines.append(f"{i}.({r['pathogen']})")
    else:
        lines.append("（無可輸出病原：缺少 pathogen_candidates 或資料不完整）")

    lines.append("-" * 123)

    # 下半部 details
    for i, r in enumerate(selected, 1):
        mainsec = primary_or_secondary(r.get("recommended_level"), r.get("decision"))
        reason = add_bas_disclaimer(str(r.get("reason") or ""), bas_mode)

        lines.append(f"{i}.")
        lines.append(f"    菌種名稱: {r['pathogen']}")
        lines.append(f"    菌種種類: {r['classification'] or 'Unknown'}")
        lines.append(f"    主次關係:(主要/次要) {mainsec}")
        lines.append(f"    原因:(簡潔) {reason}")
        lines.append("")

    out_text = "\n".join(lines).rstrip() + "\n"
    print(out_text)

    if args.save:
        save_path = Path(args.save)
        save_path.write_text(out_text, encoding="utf-8")
        print(f"✅ 已另存輸出：{save_path}")


if __name__ == "__main__":
    main()

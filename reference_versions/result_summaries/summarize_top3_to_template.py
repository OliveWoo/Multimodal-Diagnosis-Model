import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

LEVEL_TO_INT = {"Level 1": 1, "Level 2": 2, "Level 3": 3, "Level 4": 4, "Level 5": 5}

def load_json(p: Path) -> Dict[str, Any]:
    with p.open("r", encoding="utf-8") as f:
        obj = json.load(f)
    if not isinstance(obj, dict):
        raise ValueError(f"JSON root must be object: {p}")
    return obj

def norm_level(x: Any) -> Optional[str]:
    if not isinstance(x, str):
        return None
    s = x.strip()
    if s.startswith("Level") and " " not in s:
        s = "Level " + s.replace("Level", "").strip()
    return s if s in LEVEL_TO_INT else None

def decision_weight(decision: str) -> int:
    if "升級" in decision:
        return 3
    if "降級" in decision:
        return 2
    if "維持" in decision:
        return 1
    return 0

def level_delta(orig: Optional[str], rec: Optional[str]) -> int:
    if orig not in LEVEL_TO_INT or rec not in LEVEL_TO_INT:
        return 0
    return abs(LEVEL_TO_INT[orig] - LEVEL_TO_INT[rec])

def extract_patient_id_from_filename(path: Path) -> str:
    # 例：NGS_patient_3_mNGS_max_agent.json -> 3
    m = re.search(r"patient[_-](\d+)", path.name, flags=re.IGNORECASE)
    if m:
        return m.group(1)
    return path.stem

def short_reason(rationale: str, refs: List[str], max_len: int = 60) -> str:
    r = (rationale or "").replace("\n", " ").strip()
    if len(r) > max_len:
        r = r[:max_len].rstrip() + "…"
    if refs:
        # 只列前 2 個 PMID
        r += "（" + "；".join(refs[:2]) + "）"
    return r

def build_bundle_map(max_agent: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    organism_name -> {level, confidence, classification}
    """
    out: Dict[str, Dict[str, Any]] = {}
    cand = max_agent.get("pathogen_candidates", [])
    if not isinstance(cand, list):
        return out
    for p in cand:
        if not isinstance(p, dict):
            continue
        name = str(p.get("organism_name") or "").strip()
        if not name:
            continue
        lvl = norm_level(p.get("integrated_causative_level"))
        conf = p.get("integrated_causative_confidence")
        conf_val = float(conf) if isinstance(conf, (int, float)) else None
        cls = str(p.get("classification") or "").strip()  # Bacterial/Fungal/Viral
        out[name] = {"level": lvl, "confidence": conf_val, "classification": cls}
    return out

def get_reassessments(out_json: Dict[str, Any]) -> List[Dict[str, Any]]:
    rag = out_json.get("rag_answer", {})
    if not isinstance(rag, dict):
        return []
    lr = rag.get("level_reassessment_conclusion")
    if lr is None:
        return []
    if isinstance(lr, dict):
        lr_list = [lr]
    elif isinstance(lr, list):
        lr_list = lr
    else:
        return []

    cleaned: List[Dict[str, Any]] = []
    for x in lr_list:
        if not isinstance(x, dict):
            continue
        pathogen = str(x.get("pathogen") or "").strip()
        orig = norm_level(x.get("original_level"))
        rec = norm_level(x.get("recommended_level"))
        decision = str(x.get("decision") or "").strip()
        rationale = str(x.get("rationale") or "").strip()
        refs = x.get("supporting_references", [])
        if not isinstance(refs, list):
            refs = []
        refs = [str(r).strip() for r in refs if str(r).strip()]
        if pathogen:
            cleaned.append({
                "pathogen": pathogen,
                "original_level": orig,
                "recommended_level": rec,
                "decision": decision,
                "rationale": rationale,
                "supporting_references": refs,
            })
    return cleaned

def main():
    ap = argparse.ArgumentParser(description="Compare max_agent vs out; print Top 1/2/3 in your template.")
    ap.add_argument("--max_agent", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top_k", type=int, default=3)
    ap.add_argument("--reason_len", type=int, default=60)
    ap.add_argument("--save", default=None, help="Optional: save as .txt")
    args = ap.parse_args()

    max_agent_path = Path(args.max_agent)
    out_path = Path(args.out)

    max_agent = load_json(max_agent_path)
    out_json = load_json(out_path)

    patient_id = extract_patient_id_from_filename(max_agent_path)

    # 感染源（你 bundle 通常會有 dominant_source / dominant_pathogen_type）
    dominant_source = str(max_agent.get("dominant_source") or "").strip()
    dominant_type = str(max_agent.get("dominant_pathogen_type") or "").strip()
    infection_source_line = " / ".join([x for x in [dominant_source, dominant_type] if x]) or "（未提供）"

    bundle_map = build_bundle_map(max_agent)
    reassess = get_reassessments(out_json)

    # 只保留「真的有調整」：升級/降級，或 level 不同
    changed: List[Dict[str, Any]] = []
    for x in reassess:
        dec = x.get("decision", "")
        orig = x.get("original_level")
        rec = x.get("recommended_level")
        is_changed = ("維持" not in dec) or (orig != rec)
        if not is_changed:
            continue

        info = bundle_map.get(x["pathogen"], {})
        changed.append({
            **x,
            "bundle_level": info.get("level"),
            "bundle_confidence": info.get("confidence"),
            "classification": info.get("classification") or "（未知）",
        })

    # 排序（你可以改規則）
    changed.sort(
        key=lambda z: (
            decision_weight(str(z.get("decision", ""))),
            level_delta(z.get("original_level"), z.get("recommended_level")),
            (z.get("bundle_confidence") if isinstance(z.get("bundle_confidence"), (int, float)) else -1.0),
        ),
        reverse=True,
    )

    top = changed[: max(0, args.top_k)]

    # 產出你要的格式
    lines: List[str] = []
    lines.append(f"PatientID: {patient_id}")
    lines.append(f"感染源: {infection_source_line}")
    lines.append("1.")
    lines.append("2.")
    lines.append("3.")
    lines.append("-" * 115)

    for i in range(3):
        idx = i + 1
        lines.append(f"{idx}.")
        if i < len(top):
            t = top[i]
            name = t["pathogen"]
            cls = t.get("classification", "（未知）")
            role = "主要" if i == 0 else "次要"
            reason = short_reason(t.get("rationale", ""), t.get("supporting_references", []), max_len=args.reason_len)
            lines.append(f"    菌種名稱: {name}")
            lines.append(f"    菌種種類: {cls}")
            lines.append(f"    主次關係:(主要/次要) {role}")
            lines.append(f"    原因:(簡潔) {reason}")
        else:
            lines.append(f"    菌種名稱:")
            lines.append(f"    菌種種類:")
            lines.append(f"    主次關係:(主要/次要)")
            lines.append(f"    原因:(簡潔)")
        lines.append("")

    output_text = "\n".join(lines)
    print(output_text)

    if args.save:
        save_path = Path(args.save)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        save_path.write_text(output_text, encoding="utf-8")

if __name__ == "__main__":
    main()

"""Build a concise Traditional Chinese report from deterministic max picks and LLM missed-candidate review."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any

DEFAULT_MAX_SUFFIX = "mNGS_max_deterministic_resp_commensal_dominance_guardrail_opt_chosen_full"
DEFAULT_REVIEW_SUFFIX = "mNGS_missed_candidate_review"

CLASSIFICATION_ZH = {
    "Bacterial": "細菌",
    "Viral": "病毒",
    "Fungal": "真菌",
    "Parasite": "寄生蟲",
}
ROLE_ZH = {"Primary": "主要候選病原", "Secondary": "次要候選病原"}
CONFIDENCE_ZH = {"high": "高", "moderate": "中等", "low": "低"}
SIGNAL_ZH = {
    "M1_strong": "M1 強訊號",
    "M2_moderate": "M2 中等訊號",
    "M3_protected": "M3 受保護或需保留觀察",
    "M4_weak": "M4 弱訊號",
    "M5_background": "M5 背景訊號",
}
READS_TIER_ZH = {
    "R0_trace": "R0 極低量訊號",
    "R1_low": "R1 低量訊號",
    "R2_medium": "R2 中等訊號",
    "R3_high": "R3 高訊號",
    "R4_very_high": "R4 極高訊號",
}
DOMINANCE_ZH = {
    "D0_not_top": "D0 非優勢菌",
    "D1_low": "D1 低度優勢",
    "D2_moderate": "D2 中度優勢",
    "D3_dominant": "D3 明顯優勢菌",
}
SPECIMEN_CLASS_ZH = {
    "S1_sterile_systemic": "S1 無菌或系統性檢體",
    "S2_lower_respiratory": "S2 下呼吸道檢體",
    "S3_unknown_or_low_value": "S3 未知或判讀價值較低的檢體",
}
CAUTION_ZH = {
    "Evidence_limited": "證據仍有限，需結合臨床情境判讀",
    "Not_recommended_as_sole_treatment_basis": "不建議僅依此結果決定治療",
}
SUPPORT_ZH = {
    "mNGS": "mNGS",
    "culture": "培養",
    "filmarray_gmtest": "FilmArray 或 GM test",
    "image": "影像",
    "host": "宿主風險",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("patient_root", type=Path, help="Folder containing NGS_patient_<id>_json folders")
    parser.add_argument("--output", type=Path, required=True, help="Combined output JSON")
    parser.add_argument("--max-suffix", default=DEFAULT_MAX_SUFFIX)
    parser.add_argument("--review-suffix", default=DEFAULT_REVIEW_SUFFIX)
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def patient_sort_key(patient_dir: Path) -> tuple[int, str]:
    try:
        return int(patient_dir.name.split("_")[2]), patient_dir.name
    except (IndexError, ValueError):
        return 10**9, patient_dir.name


def translated_code(code: Any, mapping: dict[str, str]) -> str:
    value = str(code or "-")
    return f"{value}（{mapping.get(value, '未定義')}）" if value != "-" else "-"


def module_support_zh(support: dict[str, Any]) -> list[str]:
    result = []
    for key, label in SUPPORT_ZH.items():
        value = support.get(key)
        if value and value != "Not_available":
            result.append(f"{label}：{value}")
    return result


def find_candidate(max_data: dict[str, Any], organism_name: str) -> dict[str, Any]:
    for candidate in max_data.get("pathogen_candidates", []):
        if candidate.get("organism_name") == organism_name:
            return candidate
    return {}


def build_pick_explanation(pick: dict[str, Any], candidate: dict[str, Any]) -> str:
    rank = pick.get("rank_priority", candidate.get("rank_priority", "-"))
    reads = pick.get("reads", candidate.get("reads", "-"))
    reads_tier = translated_code(candidate.get("reads_tier"), READS_TIER_ZH)
    percentile = candidate.get("reads_percentile", "-")
    signal = translated_code(pick.get("mngs_signal_tier"), SIGNAL_ZH)
    level = pick.get("basis_level", candidate.get("integrated_causative_level", "-"))
    dominance = translated_code(candidate.get("dominance_tier"), DOMINANCE_ZH)
    specimen = translated_code(candidate.get("specimen_class"), SPECIMEN_CLASS_ZH)
    support = module_support_zh(candidate.get("module_support_summary", {}))
    support_text = "、".join(support) if support else "無額外模組支持"
    return (
        f"程式依固定規則選入：rank_priority={rank}，reads={reads}，reads_tier={reads_tier}，"
        f"reads_percentile={percentile}，mNGS_signal_tier={signal}，整合判斷為 {level}。"
        f"優勢程度為 {dominance}；檢體分類為 {specimen}；支持資訊：{support_text}。"
    )


def readable_pick(pick: dict[str, Any], max_data: dict[str, Any]) -> dict[str, Any]:
    organism_name = str(pick.get("organism_name", ""))
    candidate = find_candidate(max_data, organism_name)
    cautions = [CAUTION_ZH.get(flag, flag) for flag in pick.get("caution_flags", [])]
    return {
        "菌種": organism_name,
        "分類": CLASSIFICATION_ZH.get(str(pick.get("classification", "")), pick.get("classification", "-")),
        "角色": ROLE_ZH.get(str(pick.get("picked_role", "")), pick.get("picked_role", "-")),
        "整合層級": pick.get("basis_level", candidate.get("integrated_causative_level", "-")),
        "mNGS訊號等級": translated_code(pick.get("mngs_signal_tier"), SIGNAL_ZH),
        "rank_priority": pick.get("rank_priority", candidate.get("rank_priority", "-")),
        "reads": pick.get("reads", candidate.get("reads", "-")),
        "reads_tier": translated_code(candidate.get("reads_tier"), READS_TIER_ZH),
        "reads_percentile": candidate.get("reads_percentile", "-"),
        "dominance_tier": translated_code(candidate.get("dominance_tier"), DOMINANCE_ZH),
        "檢體分類": translated_code(candidate.get("specimen_class"), SPECIMEN_CLASS_ZH),
        "支持資訊": module_support_zh(candidate.get("module_support_summary", {})),
        "選入原因": build_pick_explanation(pick, candidate),
        "注意事項": cautions,
    }


def readable_missed(item: dict[str, Any]) -> dict[str, Any]:
    confidence = str(item.get("confidence", ""))
    return {
        "菌種": item.get("organism_name", ""),
        "人工覆核優先度": f"{confidence}（{CONFIDENCE_ZH.get(confidence, '未定義')}）" if confidence else "-",
        "LLM提出原因": item.get("rationale_zh", ""),
        "若忽略的風險": item.get("risk_if_ignored", ""),
        "未自動升級原因": item.get("why_not_auto_upgrade", ""),
        "建議下一步確認": item.get("suggested_next_check", ""),
    }


def rag_visible_review_items(review_data: dict[str, Any]) -> list[dict[str, Any]]:
    high = [item for item in review_data.get("review_high_priority") or [] if isinstance(item, dict)]
    context = [item for item in review_data.get("review_context_needed") or [] if isinstance(item, dict)]
    if high or context:
        return high + context
    return [item for item in review_data.get("possible_missed_pathogens") or [] if isinstance(item, dict)]


def build_patient(patient_dir: Path, max_suffix: str, review_suffix: str) -> dict[str, Any] | None:
    patient_id = patient_dir.name.split("_")[2]
    summary_dir = patient_dir / "summary_outputs"
    max_path = summary_dir / f"NGS_patient_{patient_id}_{max_suffix}.json"
    review_path = summary_dir / f"NGS_patient_{patient_id}_{review_suffix}.json"
    if not max_path.exists() or not review_path.exists():
        return None
    max_data = read_json(max_path)
    review_data = read_json(review_path)
    summary = max_data.get("best_available_summary", {})
    picked = [readable_pick(item, max_data) for item in summary.get("picked_pathogens", [])]
    missed = [readable_missed(item) for item in rag_visible_review_items(review_data)]
    return {
        "patient_id": patient_id,
        "正式max已選入菌種數": len(picked),
        "正式max已選入菌種": picked,
        "LLM建議人工覆核的可能漏掉菌種數": len(missed),
        "LLM建議人工覆核的可能漏掉菌種": missed,
    }


def main() -> int:
    args = parse_args()
    patient_dirs = sorted(args.patient_root.glob("NGS_patient_*_json"), key=patient_sort_key)
    patients = []
    missing = []
    for patient_dir in patient_dirs:
        patient = build_patient(patient_dir, args.max_suffix, args.review_suffix)
        if patient is None:
            missing.append(patient_dir.name)
        else:
            patients.append(patient)
    output = {
        "報告用途": "正式 deterministic max 選入結果與 LLM missed-candidate safety review 整併閱讀版",
        "重要說明": [
            "正式 max 已選入菌種為 deterministic scorer 的正式輸出。",
            "LLM RAG-visible review 項目僅供安全覆核，不會自動加入正式結果。",
            "本檔刻意省略 review_low_specificity、omitted_with_reason 與其他未列入 RAG-visible review 的菌種。",
        ],
        "產生時間": datetime.now().isoformat(timespec="seconds"),
        "病人數": len(patients),
        "正式max已選入菌種總數": sum(p["正式max已選入菌種數"] for p in patients),
        "LLM建議人工覆核的可能漏掉菌種總數": sum(p["LLM建議人工覆核的可能漏掉菌種數"] for p in patients),
        "缺少必要輸入檔案的病人資料夾": missing,
        "病人結果": patients,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8-sig")
    print(f"output={args.output}")
    print(f"patients={len(patients)}")
    print(f"formal_picked_total={output['正式max已選入菌種總數']}")
    print(f"possible_missed_total={output['LLM建議人工覆核的可能漏掉菌種總數']}")
    print(f"missing_inputs={len(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Build a human-reviewable rule matrix for every organism in KH release F."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools.evaluate_multi_assay_candidate_entry import load_answer_rows
from tools.recalculate_kh_benchmark_metrics import match_type


DEFAULT_RELEASE = Path("outputs/runs/2026-10-01_KH_ablation_F_current_full_v1")
DEFAULT_OUTPUT = Path("outputs/runs/2026-10-01_KH_organism_rule_review_v1")
DEFAULT_DOC = Path("docs/workflow/KH_CURRENT_ORGANISM_RULE_REVIEW_20261001_zh.md")

TIER_ORDER = {"Picked": 0, "Possible": 1, "Fallback-Possible": 2, "Context": 3, "Absent": 4}
PRIORITY_ORDER = {
    "P0_policy_conflict": 0,
    "P1_manual_review": 1,
    "P2_verify_evidence_chain": 2,
    "P3_context_only": 3,
}

FAMILY_RULES = {
    "typical_respiratory_pathogen": {
        "zh": "典型呼吸道病原",
        "analytical": "下呼吸道 top-1 且可重現可直接 Picked；top-3 先進 High。Level 1/2 exact hospital evidence 可直接 Picked。",
        "reporting": "High 若 exact、rank≤3、≥2 positive tests 且可重現，可進 Possible；直接證據 route 可放寬到 rank≤10。",
        "review": "確認 technical repeat 未被當成多個獨立臨床證據，並核對 culture/PCR 是否在 ±48 小時。",
    },
    "hospital_or_nonfermenter_gnb": {
        "zh": "院內／非發酵革蘭陰性菌",
        "analytical": "與典型呼吸道菌相同；top-1＋可重現可 Picked，top-3 進 High。",
        "reporting": "可走一般 bacterial Possible、targeted-assay bridge 或 analytically-compelling Picked。",
        "review": "特別檢查定植、環境來源及 HAP/MDR 病史是否只是背景，而非重複計分。",
    },
    "other_respiratory_virus": {
        "zh": "其他呼吸道病毒",
        "analytical": "top-3 mNGS 先進 High；exact Level 1/2 直接檢驗可 Picked。",
        "reporting": "可由呼吸道 direct assay Picked；Human respirovirus 3 可用 exact mNGS member＋group assay bridge。",
        "review": "group assay 只能支持群組，不能冒充 exact subtype；目前 v8 whitelist 僅 Human respirovirus 3。",
    },
    "high_consequence_opportunistic": {
        "zh": "高後果伺機性感染",
        "analytical": "analytical-only 不能自動 Picked；需宿主風險、影像或病原特異支持。",
        "reporting": "PJP 有獨立的 Picked／Possible 規則；其他菌需 family-specific gate。",
        "review": "宿主風險只提供背景，不能單獨證明菌種；缺 PCR／BDG 不當作陰性。",
    },
    "mold_or_opportunistic_fungus": {
        "zh": "黴菌／伺機性真菌",
        "analytical": "top-3＋可重現＋宿主支持進 High；Level 1/2 direct evidence＋宿主支持才可在 analytical stage Picked。",
        "reporting": "promotion 要求 exact identity、相容檢體、可重現與 ±48 小時 direct Level 1/2 evidence。",
        "review": "確認 direct evidence 是同一菌、同一 episode，並區分 airway detection 與侵襲性感染。",
    },
    "candida_or_yeast": {
        "zh": "Candida／酵母菌",
        "analytical": "屬 context-prone family；即使 direct L1/2 也先進 High，不由 mNGS 單獨 Picked。",
        "reporting": "Picked 需 exact species 與 ±48 小時血液／無菌部位／組織陽性；group label 不轉移 species identity。",
        "review": "肺部 Candida 訊號不等同 Candida pneumonia；Fallback-Possible 只能視為低信心 best available。",
    },
    "herpes_or_reactivation_virus": {
        "zh": "疱疹／再活化病毒",
        "analytical": "不做 analytical-only Picked；保留 Context／High 供再活化解讀。",
        "reporting": "cross-molecule、rank≤3 且通常需 ≥6 tests，可進 Possible；或宿主支持的單一 rank-1 訊號進 Possible。",
        "review": "沒有 reactivation wording 是未知而非陰性；要升 Picked 仍需 viral load、組織或器官侵犯證據。",
    },
    "oral_aspiration_or_anaerobe": {
        "zh": "口腔／吸入性肺炎／厭氧菌",
        "analytical": "預設不由分析訊號自動 Picked；強訊號留 Context。",
        "reporting": "v4C 將相符 aspiration history 的 Context 升 Priority；exact、rank≤3、≥3 tests、cross-molecule 且 QC 完整才進 Possible。",
        "review": "病史只能升一級且不能重複當菌種證據；Picked 仍需要獨立同菌或侵襲性證據。",
    },
    "skin_airway_colonizer_prone": {
        "zh": "皮膚／氣道定植傾向菌",
        "analytical": "預設不 analytical-only Picked；強或重現訊號進 Context。",
        "reporting": "rank-1＋≥3 tests＋cross-molecule 可進 Possible；若 exact sterile-site pure culture ±48h，可走 cross-site Picked。",
        "review": "需逐案確認污染／定植與真正感染；technical repetition 不能取代獨立採血。",
    },
    "environmental_low_specificity": {
        "zh": "環境型／低特異性菌",
        "analytical": "通常只在 Context/Audit；Route A 要 exact、下呼吸道、至少 2 technical branches、每支 rank-1、RPM≥10，最多 High。",
        "reporting": "Route A 只能 Possible；direct event-aligned evidence 或無其他候選時 fallback 才可能進完整報告。",
        "review": "最容易 overfit；要確認 normalization、technical branch 關係及污染可能性。",
    },
    "gi_urinary_or_nonpulmonary_prone": {
        "zh": "腸胃／泌尿／非肺部來源傾向菌",
        "analytical": "低特異性處理，不能自動 Picked；足夠強時保留 Context。",
        "reporting": "rank-1＋≥3 tests＋cross-molecule，或 rank≤4＋≥6 tests 的強 cross-molecule 訊號，可進 Possible。",
        "review": "確認不是其他感染部位、轉位或非肺部事件；缺 culture/PCR 不直接視為陰性。",
    },
    "unmapped_or_uncertain": {
        "zh": "未確認／generic identity",
        "analytical": "不能 Picked；generic group／spp／identification-to-follow 只能留 audit。",
        "reporting": "family-specific guardrail 阻止進完整陽性報告。",
        "review": "依使用者既定決定，應 audit-visible 但不進病人 Context；目前 7 筆仍標 Context，需確認並修正呈現層。",
    },
}

ROUTE_ZH = {
    "existing_clinical_picked": "已通過 Picked，reporting 保留",
    "event_aligned_direct_detection_or_level3_culture": "事件對齊 direct detection／Level 3 culture → Possible",
    "reproducible_top_ranked_typical_or_hospital_bacterium": "典型／院內菌 top-ranked＋可重現 → Possible",
    "reproducible_cross_molecule_top_ranked_low_specificity_signal": "低特異性菌 rank-1＋cross-molecule 重現 → Possible",
    "analytically_compelling_low_specificity_repeat": "低特異性 technical repeat Route A → Possible",
    "pneumocystis_reproducible_top_rank_with_recorded_host_context": "PJP top-ranked 重現＋宿主背景 → Possible",
    "reactivation_family_cross_molecule_reproducible_signal": "再活化病毒強 cross-molecule → Possible",
    "reactivation_family_top_rank_with_recorded_host_context": "再活化病毒 rank-1＋宿主背景 → Possible",
    "strong_repeated_cross_molecule_signal_despite_source_caution": "來源可疑但極強 repeated DNA/RNA → Possible",
    "aspiration_history_priority_with_complete_cross_molecule_signal": "aspiration 病史 Priority＋完整 cross-molecule → Possible",
    "no_picked_or_possible_best_available_fallback": "全病人無一般 Picked/Possible 時的 best-available fallback",
    "family_specific_guardrail_required": "特殊 family guardrail；保留 Context",
    "possible_pathogen_threshold_not_met": "未達 Possible 門檻；保留 Context",
    "analytical_only_low_normalized_burden_context": "所有 per-test RPM 皆低於 5；降回 Context",
    "analytical_only_inherited_taxonomy_without_host_context": "非 exact taxonomy 且缺宿主背景；降回 Context",
    "bacterial_analytically_compelling_cross_molecule": "exact bacteria：top-3、≥3 tests、cross-molecule、≥2 tests RPM≥5 → Picked",
    "bacterial_targeted_respiratory_assay_bridge": "exact bacteria＋事件對齊 targeted assay＋完整 cross-molecule → Picked",
    "colonizer_prone_cross_site_convergence": "定植傾向菌的下呼吸道 mNGS＋無菌部位 pure culture convergence → Picked",
    "respiratory_virus_exact_mngs_plus_group_assay": "exact respiratory-virus mNGS member＋群組 assay → Picked",
    "respiratory_assay_target": "醫院端呼吸道 assay target → Picked",
    "candida_exact_species_invasive_event": "exact Candida＋事件對齊侵襲性證據 → Picked",
    "opportunistic_fungus_direct_support": "伺機性真菌＋direct Level 1/2＋重現 → Picked",
    "pjp_image_cross_molecule_triangulation": "PJP 影像 context＋cross-molecule 重現 → Picked",
    "preserve_existing_tier": "不需 promotion，保留原層級",
    "no_promotion": "未通過 promotion，保留 High／Context",
    "broad_group_assay_consolidated_under_exact_mngs_member": "廣義 group label 收在 exact mNGS member 下，不重複報告",
}


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise TypeError(path)
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["organism_name"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).casefold() == "true"


def as_int(value: Any) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def benchmark_status(
    patient: int,
    name: str,
    selected: bool,
    answers: dict[int, list[str]],
) -> str:
    gold = answers.get(patient, [])
    if not gold:
        return "unlabeled_patient"
    matched = any(match_type(name, item, genus_relaxed=False) for item in gold)
    if selected:
        return "reported_TP" if matched else "reported_FP"
    return "answer_not_reported" if matched else "context_not_scored_as_FP"


def review_priority(row: dict[str, Any]) -> tuple[str, str]:
    family = row["taxonomy_family"]
    tier = row["final_reporting_tier"]
    mapping = row["taxonomy_mapping_status"]
    if family == "unmapped_or_uncertain":
        return "P0_policy_conflict", "未確認／generic identity 目前仍標 Context；應確認是否改成純 Audit。"
    if tier == "Fallback-Possible":
        return "P1_manual_review", "Fallback 是低信心 best available，必須人工確認。"
    if tier in {"Picked", "Possible"} and (
        family in {
            "environmental_low_specificity",
            "skin_airway_colonizer_prone",
            "gi_urinary_or_nonpulmonary_prone",
            "herpes_or_reactivation_virus",
            "oral_aspiration_or_anaerobe",
            "candida_or_yeast",
            "high_consequence_opportunistic",
        }
        or mapping != "exact_species"
    ):
        return "P1_manual_review", "特殊／低特異性 family 或非 exact mapping 已進報告，需逐案確認 gate 是否合理。"
    if tier == "Picked" or tier == "Possible":
        return "P2_verify_evidence_chain", "已進報告；確認 per-test、時序與 direct evidence provenance。"
    return "P3_context_only", "目前只在 Context，不影響陽性報告；確認分類與保留理由即可。"


def candidate_rows(release_root: Path, answers: dict[int, list[str]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    paths = sorted((release_root / "patient_results").glob("patient_*_integrated_decision.json"))
    for path in paths:
        payload = read_json(path)
        patient = int(payload["patient_id"])
        for item in payload.get("all_candidates") or []:
            taxonomy = item.get("taxonomy") or {}
            analytical = item.get("analytical_evidence") or {}
            clinical = item.get("clinical_evidence") or {}
            decision = item.get("decision_path") or {}
            promotion = item.get("promotion") or {}
            reporting = item.get("reporting") or {}
            audit = item.get("audit") or {}
            selected = as_bool(reporting.get("selected_for_complete_report"))
            row = {
                "patient_id": patient,
                "organism_name": item.get("organism_name"),
                "candidate_source": item.get("candidate_source"),
                "taxid": taxonomy.get("taxid") if taxonomy.get("taxid") is not None else "",
                "taxonomic_rank": taxonomy.get("rank") or "",
                "taxonomy_family": taxonomy.get("family") or "unmapped_or_uncertain",
                "taxonomy_mapping_status": taxonomy.get("mapping_status") or "",
                "classification_confidence": taxonomy.get("classification_confidence") or "",
                "taxonomy_rule_ids": "|".join(taxonomy.get("matched_rule_ids") or []),
                "best_rank": analytical.get("best_rank") if analytical.get("best_rank") is not None else "",
                "positive_test_count": analytical.get("selected_positive_test_count") or 0,
                "reproducibility_axis": bool(analytical.get("reproducibility_axis")),
                "cross_molecule": bool(analytical.get("cross_molecule_selected")),
                "specimen_context": clinical.get("specimen_context") or "",
                "direct_hospital_level": clinical.get("direct_hospital_level") or "",
                "event_aligned_direct_positive_count": clinical.get("event_aligned_direct_positive_count") or 0,
                "host_support": bool(clinical.get("host_support")),
                "precise_identity": bool(clinical.get("precise_identity")),
                "analytical_route": decision.get("analytical_route") or "",
                "history_adjusted_route": decision.get("history_adjusted_route") or "",
                "pre_promotion_decision": decision.get("pre_promotion_clinical_decision") or "",
                "post_promotion_decision": decision.get("post_promotion_clinical_decision") or "",
                "promotion_outcome": promotion.get("outcome") or "",
                "promotion_route": promotion.get("route") or "",
                "promotion_blockers": "|".join(promotion.get("blockers") or []),
                "final_reporting_tier": decision.get("final_reporting_tier") or reporting.get("tier") or "Context",
                "selected_for_complete_report": selected,
                "reporting_role": reporting.get("role") or "",
                "reporting_route": reporting.get("route") or "",
                "reporting_blockers": "|".join(reporting.get("blockers") or []),
                "reporting_cautions": "|".join(reporting.get("cautions") or []),
                "history_rule_ids": "|".join(audit.get("history_rule_ids") or []),
                "decision_rule_ids": "|".join(audit.get("decision_rule_ids") or []),
                "decision_reasons": "|".join(audit.get("reasons") or []),
                "benchmark_status": benchmark_status(patient, str(item.get("organism_name") or ""), selected, answers),
            }
            priority, reason = review_priority(row)
            row["review_priority"] = priority
            row["review_reason_zh"] = reason
            output.append(row)
    return sorted(output, key=lambda row: (row["patient_id"], str(row["organism_name"])))


def organism_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["organism_name"])].append(row)
    output: list[dict[str, Any]] = []
    for name, items in sorted(grouped.items(), key=lambda pair: pair[0].casefold()):
        tiers = Counter(str(item["final_reporting_tier"]) for item in items)
        statuses = Counter(str(item["benchmark_status"]) for item in items)
        priorities = sorted({str(item["review_priority"]) for item in items}, key=PRIORITY_ORDER.get)
        best_priority = priorities[0]
        priority_reason = next(item["review_reason_zh"] for item in items if item["review_priority"] == best_priority)
        ranks = [as_int(item["best_rank"]) for item in items]
        ranks = [value for value in ranks if value is not None]
        routes = sorted({str(item["reporting_route"]) for item in items if item["reporting_route"]})
        promotion_routes = sorted({str(item["promotion_route"]) for item in items if item["promotion_route"] not in {"", "preserve_existing_tier", "no_promotion"}})
        families = sorted({str(item["taxonomy_family"]) for item in items})
        family = families[0] if len(families) == 1 else "|".join(families)
        family_rule = FAMILY_RULES.get(families[0], {}) if len(families) == 1 else {}
        route_descriptions = [ROUTE_ZH.get(route, route) for route in promotion_routes + routes]
        output.append({
            "organism_name": name,
            "taxid": "|".join(sorted({str(item["taxid"]) for item in items if str(item["taxid"])})),
            "taxonomy_family": family,
            "family_zh": family_rule.get("zh", "需確認"),
            "mapping_status": "|".join(sorted({str(item["taxonomy_mapping_status"]) for item in items})),
            "taxonomy_rule_ids": "|".join(sorted({rule for item in items for rule in str(item["taxonomy_rule_ids"]).split("|") if rule})),
            "patient_count": len({int(item["patient_id"]) for item in items}),
            "patients": "|".join(str(value) for value in sorted({int(item["patient_id"]) for item in items})),
            "occurrence_count": len(items),
            "tier_counts": "|".join(f"{tier}:{tiers[tier]}" for tier in sorted(tiers, key=TIER_ORDER.get)),
            "best_tier": min(tiers, key=TIER_ORDER.get),
            "reported_count": sum(as_bool(item["selected_for_complete_report"]) for item in items),
            "benchmark_status_counts": "|".join(f"{key}:{value}" for key, value in sorted(statuses.items())),
            "best_rank": min(ranks) if ranks else "",
            "maximum_positive_test_count": max(int(item["positive_test_count"]) for item in items),
            "any_cross_molecule": any(as_bool(item["cross_molecule"]) for item in items),
            "direct_hospital_levels": "|".join(sorted({str(item["direct_hospital_level"]) for item in items if str(item["direct_hospital_level"])})),
            "actual_promotion_routes": "|".join(promotion_routes),
            "actual_reporting_routes": "|".join(routes),
            "actual_rule_summary_zh": "；".join(dict.fromkeys(route_descriptions)),
            "family_rule_zh": family_rule.get("analytical", "") + " " + family_rule.get("reporting", ""),
            "review_priority": best_priority,
            "review_reason_zh": priority_reason,
            "review_question_zh": family_rule.get("review", "請確認 taxonomy 與 clinical role。"),
        })
    return output


def markdown_table(rows: list[dict[str, Any]], fields: list[tuple[str, str]]) -> str:
    lines = ["| " + " | ".join(title for _, title in fields) + " |", "|" + "|".join("---" for _ in fields) + "|"]
    for row in rows:
        values = []
        for key, _ in fields:
            text = str(row.get(key, "")).replace("|", "／").replace("\n", " ")
            values.append(text)
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def build_document(rows: list[dict[str, Any]], organisms: list[dict[str, Any]]) -> str:
    family_counts = Counter(row["taxonomy_family"] for row in rows)
    family_reported = Counter(row["taxonomy_family"] for row in rows if row["selected_for_complete_report"])
    family_lines = []
    for family, config in FAMILY_RULES.items():
        family_lines.append({
            "family": family,
            "中文": config["zh"],
            "候選數": family_counts.get(family, 0),
            "報告數": family_reported.get(family, 0),
            "analytical／Picked": config["analytical"],
            "Possible／限制": config["reporting"],
            "你要確認": config["review"],
        })

    non_strict = [row for row in rows if row["final_reporting_tier"] in {"Possible", "Fallback-Possible"}]
    promoted = [row for row in rows if row["final_reporting_tier"] == "Picked" and row["promotion_outcome"] == "promoted_to_picked_shadow"]
    conflicts = [row for row in rows if row["review_priority"] == "P0_policy_conflict"]

    appendices = []
    for family in FAMILY_RULES:
        selected = [row for row in organisms if row["taxonomy_family"] == family]
        if not selected:
            continue
        appendices.append(
            f"### {FAMILY_RULES[family]['zh']}（`{family}`）\n\n"
            + markdown_table(selected, [
                ("organism_name", "菌名"),
                ("taxid", "taxid"),
                ("mapping_status", "mapping"),
                ("patients", "病人"),
                ("tier_counts", "目前結果"),
                ("actual_rule_summary_zh", "實際套用路徑"),
                ("review_priority", "覆核優先級"),
            ])
        )

    return f"""# KH 主流程：每個菌種的規則覆核表（2026-10-01）

## 1. 這份文件怎麼看

本文件依目前完整 F release 的 {len(rows)} 筆 patient-organism 候選，整理 {len(organisms)} 個實際菌名的 taxonomy、scorer、promotion 與 reporting 規則。它描述的是「目前程式實際怎麼判」，不是臨床 guideline，也不是對每個菌的醫學正確性做最終背書。

同一菌在不同病人可能因 rank、test 數、DNA/RNA、檢體、醫院端證據或病史不同而落在不同層級；完整 {len(rows)} 筆逐病人規則請看 `patient_organism_rule_matrix.csv`。

## 2. 全域共同規則

1. 只接受 reads>0 的 per-test 訊號；不同 test 的 reads／RPM 永不相加。
2. rank 1 為 top、rank≤3 為 strong、rank≤10 為 moderate。
3. ≥2 個 selected tests 或 DNA/RNA 都檢出，記為一個 reproducibility axis；technical repeat 不是獨立臨床證據。
4. 單一 test rank≤3，或可重現且 rank≤10，可救回 Context；更弱的訊號只留 Audit/Hold。
5. direct evidence 原則上需與 mNGS 在 ±48 小時內，且 exact／approved alias 才算同菌支持。
6. generic group、spp、complex、evidence、identification-to-follow 不能正式 Picked。
7. 病史每次最多升降一級、不能刪除候選、不能單靠病史 Picked，也不能在 promotion gate 重複計分。
8. 缺少 PCR、culture、viral load 或藥物劑量視為未知，不視為陰性。
9. Picked 是高信心子集；Possible 是 recall layer；Fallback-Possible 是沒有任何一般 Picked／Possible 時的低信心 best available；Context 不列為陽性報告。

## 3. Family-specific 規則

{markdown_table(family_lines, [('family','family'),('中文','中文'),('候選數','候選數'),('報告數','完整報告數'),('analytical／Picked','analytical／Picked 規則'),('Possible／限制','Possible／限制'),('你要確認','你要確認')])}

## 4. 目前最需要你確認的規則衝突

你先前已決定「taxonomy 未確認菌保留於 Audit，但不進病人 Context」。目前以下 {len(conflicts)} 筆 generic／unmapped hospital labels 雖未進完整陽性報告，仍以 `final_reporting_tier=Context` 保存：

{markdown_table(conflicts, [('patient_id','病人'),('organism_name','名稱'),('candidate_source','來源'),('post_promotion_decision','臨床層級'),('final_reporting_tier','目前 reporting tier'),('reporting_route','阻擋路徑')])}

建議修正方式是保留原始 observation 與 audit trail，但把 user-facing tier 改成 `Audit`，不參與 Context 清單、fallback 排序或報告數量；generic label 對應到 exact species 時只作 supporting evidence，不另報一個病原。

## 5. 不是原始 Picked、但目前進入完整報告的菌

以下為全部 Possible／Fallback-Possible，最值得逐項判斷規則是否合理：

{markdown_table(non_strict, [('patient_id','病人'),('organism_name','菌名'),('taxonomy_family','family'),('taxonomy_mapping_status','mapping'),('best_rank','best rank'),('positive_test_count','tests'),('cross_molecule','DNA/RNA'),('direct_hospital_level','direct level'),('final_reporting_tier','結果'),('reporting_route','實際規則'),('benchmark_status','benchmark')])}

特別提醒：

- P6 `Acinetobacter ursingii` 與 P11 `Candida albicans` 是 Fallback-Possible，不是一般 Possible，更不是 Picked。
- P15 `Brucella intermedia` 是 taxonomy v2.6 後新增的 report-level FP，應優先確認其 family 是否不應直接視為 hospital GNB。
- 多個 Corynebacterium 與 Staphylococcus 候選是 low-specificity Possible；其規則依賴 rank、test 數及 cross-molecule，而非宣稱已證實感染。
- P30 `Bacteroides fragilis` 是目前唯一由 v4C aspiration route 新增的 TP。

## 6. High／Context 經 promotion 升成 Picked 的菌

{markdown_table(promoted, [('patient_id','病人'),('organism_name','菌名'),('taxonomy_family','family'),('best_rank','rank'),('positive_test_count','tests'),('cross_molecule','DNA/RNA'),('direct_hospital_level','direct level'),('promotion_route','promotion 規則'),('benchmark_status','benchmark')])}

這 15 筆需特別確認的是：P10 Corynebacterium cross-site convergence、P10 Human respirovirus 3 的 group-assay bridge、P14 PJP 的 image＋cross-molecule、三個 analytical-only bacterial promotion，以及 Candida／真菌 direct evidence 的 ±48 小時與 specimen identity。

## 7. 物種專屬或窄版例外

| 對象 | 現行規則 | 上限／限制 | 建議確認 |
|---|---|---|---|
| Pneumocystis jirovecii | Possible：rank≤3、≥2 tests、可重現＋宿主背景；Picked：可重現且宿主 triage，或影像 context＋cross-molecule | 缺 PCR／BDG 不扣分；宿主背景單獨不能 Picked | P9 Possible 與 P14 Picked 是否符合預期 |
| Human respirovirus 3 | exact mNGS member rank≤3＋±48h lower-respiratory group assay Level 1/2 | whitelist 目前只含此 species；group assay 不當 exact subtype | 是否要擴充為通用 approved group/member 表，而非單物種 whitelist |
| Candida／Yeast | exact Candida species＋±48h blood/sterile/tissue positive 可 Picked | airway evidence不等於 Candida pneumonia；group label 不取代 species | P11 Candida albicans fallback 是否應保留 |
| Herpes／reactivation virus | 強 cross-molecule 或 rank-1＋宿主背景可 Possible | 無 reactivation 標示不當陰性；Picked 仍需器官侵犯／定量證據 | P9 VZV 單一 test Possible 是否過寬 |
| Corynebacterium／skin-airway colonizer | top-ranked repeated cross-molecule 可 Possible；exact sterile-site culture convergence 可 Picked | technical repeats 非獨立 specimen；污染標記仍保留 caution | P10 Picked、P13/17/29/32/37/38 Possible 逐案確認 |
| analytical-only bacteria | exact、top-3、≥3 tests、cross-molecule、全部 normalized/QC 合格、至少 2 tests RPM≥5、無 colonization 反證可 Picked | 不接受 genus-inherited；RPM 不跨 test 相加 | 是否仍接受 RPM 5 作研究門檻 |
| aspiration anaerobe | v4C Priority＋exact、rank≤3、≥3 tests、cross-molecule、全部 QC evaluable可 Possible | 病史不得重複計分，不能直接 Picked | P30 B. fragilis 是目前唯一通過者 |
| Fallback | 病人完全沒有一般 Picked／Possible 時，依 tier、evaluable rank、cross-molecule、重現性、檢體、舊 rank 及 taxonomy penalty 排序 | 最多 2 隻；第二名必須 core rank 完全相同 | 是否希望 fallback 永遠至少選 1 隻，或允許「無可報菌」 |

## 8. 每個實際菌種的目前規則與結果

以下 `實際套用路徑` 是這 33 位資料中真正觸發的規則；沒有觸發 Picked／Possible 的菌通常顯示 family guardrail 或 threshold not met。

{chr(10).join(appendices)}

## 9. 建議覆核順序

1. 先決定 7 個 unmapped/generic label 是否由 Context 改成純 Audit。
2. 逐一確認 25 筆 Possible／Fallback-Possible，優先 P15 Brucella、P6 Acinetobacter、P11 Candida、P9 VZV。
3. 確認 15 筆 promotion→Picked，特別是 analytical-only bacteria、P10 Corynebacterium 與 Human respirovirus 3。
4. 再檢查各 family 的 Context-only 菌；這些目前不影響陽性報告，可以分批覆核。
5. 規則確認後凍結 policy；不要再依這 33 位的 F1 新增物種特例，留待 held-out／external validation。

## 10. 對應檔案

- 一菌一列：`outputs/runs/2026-10-01_KH_organism_rule_review_v1/organism_rule_matrix.csv`
- 一病人一菌一列：`outputs/runs/2026-10-01_KH_organism_rule_review_v1/patient_organism_rule_matrix.csv`
- 需優先確認：`outputs/runs/2026-10-01_KH_organism_rule_review_v1/high_priority_review.csv`
- 本文件只說明目前程式規則；臨床合理性仍需教授／臨床端確認。
"""


def run(release_root: Path, output_root: Path, doc_path: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    summary = read_json(release_root / "summary.json")
    answers_path = Path(summary["reporting_evaluation"]["answer_path"])
    answers = load_answer_rows(answers_path)
    rows = candidate_rows(release_root, answers)
    organisms = organism_rows(rows)
    if len(rows) != summary["candidate_count"]:
        raise ValueError(f"Candidate count mismatch: {len(rows)} != {summary['candidate_count']}")

    output_root.mkdir(parents=True)
    write_csv(output_root / "patient_organism_rule_matrix.csv", rows)
    write_csv(output_root / "organism_rule_matrix.csv", organisms)
    high_priority = [row for row in rows if row["review_priority"] in {"P0_policy_conflict", "P1_manual_review"}]
    write_csv(output_root / "high_priority_review.csv", high_priority)

    document = build_document(rows, organisms)
    doc_path.parent.mkdir(parents=True, exist_ok=True)
    doc_path.write_text(document, encoding="utf-8")

    manifest = {
        "schema_version": "kh_organism_rule_review.v1",
        "status": "complete_frozen_from_release_F",
        "source_release": str(release_root.resolve()),
        "source_manifest_sha256": sha256_file(release_root / "run_manifest.json"),
        "patient_count": summary["patient_count"],
        "candidate_count": len(rows),
        "unique_organism_count": len(organisms),
        "taxonomy_family_count": len({row["taxonomy_family"] for row in rows}),
        "final_tier_counts": dict(Counter(row["final_reporting_tier"] for row in rows)),
        "review_priority_counts": dict(Counter(row["review_priority"] for row in rows)),
        "unmapped_generic_context_count": sum(row["review_priority"] == "P0_policy_conflict" for row in rows),
        "limitations": [
            "This describes deterministic behavior in the fixed 33-patient development cohort.",
            "Review priority flags structural/policy consistency risk, not clinical truth.",
            "Clinical validity still requires domain review and held-out/external validation.",
        ],
        "outputs": {
            "organism_rule_matrix": str((output_root / "organism_rule_matrix.csv").resolve()),
            "organism_rule_matrix_sha256": sha256_file(output_root / "organism_rule_matrix.csv"),
            "patient_organism_rule_matrix": str((output_root / "patient_organism_rule_matrix.csv").resolve()),
            "patient_organism_rule_matrix_sha256": sha256_file(output_root / "patient_organism_rule_matrix.csv"),
            "high_priority_review": str((output_root / "high_priority_review.csv").resolve()),
            "high_priority_review_sha256": sha256_file(output_root / "high_priority_review.csv"),
            "report": str(doc_path.resolve()),
            "report_sha256": sha256_file(doc_path),
        },
    }
    write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--doc", type=Path, default=DEFAULT_DOC)
    args = parser.parse_args()
    print(json.dumps(run(args.release, args.output, args.doc), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

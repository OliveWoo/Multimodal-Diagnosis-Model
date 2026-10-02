from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


TIER_ORDER = {"Picked": 0, "Possible": 1, "Fallback-Possible": 2, "Context": 3, "Audit": 4}

ROUTE_ZH = {
    "existing_clinical_picked": "既有 deterministic scorer 已達 Picked",
    "reproducible_cross_molecule_top_ranked_low_specificity_signal": "低特異性菌：rank 1、多 test、DNA/RNA 一致",
    "reproducible_top_ranked_typical_or_hospital_bacterium": "典型／院內菌：top-ranked 且可重現",
    "event_aligned_direct_detection_or_level3_culture": "±48 小時院端同菌支持或 Level 3 culture",
    "strong_repeated_cross_molecule_signal_despite_source_caution": "來源有疑慮，但訊號強且跨 DNA/RNA 重現",
    "no_picked_or_possible_best_available_fallback": "無一般 Picked／Possible 時的 best-available fallback",
    "analytically_compelling_low_specificity_repeat": "低特異性菌 Route A technical-repeat 規則",
    "aspiration_history_priority_with_complete_cross_molecule_signal": "吸入風險＋完整跨 DNA/RNA 訊號",
    "reactivation_family_cross_molecule_reproducible_signal": "再活化病毒：強跨 DNA/RNA 重現",
    "reactivation_family_top_rank_with_recorded_host_context": "再活化病毒：rank 1＋宿主背景",
    "pneumocystis_reproducible_top_rank_with_recorded_host_context": "PJP：top-ranked 可重現＋宿主背景",
    "bacterial_targeted_respiratory_assay_bridge": "exact bacteria＋院端 targeted respiratory assay",
    "bacterial_analytically_compelling_cross_molecule": "exact bacteria＋top 3＋≥3 tests＋跨 DNA/RNA＋RPM gate",
    "colonizer_prone_cross_site_convergence": "定植傾向菌＋相容無菌部位／跨部位同菌證據",
    "respiratory_virus_exact_mngs_plus_group_assay": "exact mNGS virus member＋院端 group assay",
    "respiratory_assay_target": "院端呼吸道 assay target",
    "pjp_image_cross_molecule_triangulation": "PJP：影像 context＋跨 DNA/RNA",
    "opportunistic_fungus_direct_support": "伺機性真菌＋±48 小時 direct support",
    "candida_exact_species_invasive_event": "exact Candida＋血液／無菌部位／組織事件",
}

FAMILY_CAUTION = {
    "candida_or_yeast": "肺部 Candida／yeast 訊號可能反映定植；需區分侵襲性事件。",
    "herpes_or_reactivation_virus": "可能為再活化或 shedding；缺少 viral load／器官侵犯證據不等於陰性。",
    "skin_airway_colonizer_prone": "具有皮膚污染或氣道定植可能，需結合獨立檢體與臨床事件。",
    "environmental_low_specificity": "環境型／低特異性菌；technical repeat 不是獨立生物學證據。",
    "gi_urinary_or_nonpulmonary_prone": "可能來自肺外感染或轉位；需確認本次肺部事件相容性。",
    "oral_aspiration_or_anaerobe": "與吸入相關但病史不是菌種特異證據。",
    "high_consequence_opportunistic": "宿主風險只提供背景，不能單獨證實病原。",
    "mold_or_opportunistic_fungus": "需區分 airway detection、定植與侵襲性感染。",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["empty"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: Any) -> bool:
    return str(value).strip().casefold() in {"true", "1", "yes"}


def as_int(value: Any) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def pct(value: float | None) -> str:
    return "—" if value is None else f"{value:.3f}"


def md_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    lines = ["| " + " | ".join(title for _, title in columns) + " |", "|" + "|".join("---" for _ in columns) + "|"]
    for row in rows:
        values = []
        for key, _ in columns:
            value = str(row.get(key, "")).replace("|", "／").replace("\n", " ")
            values.append(value)
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def effective_route(row: dict[str, str]) -> tuple[str, str]:
    if row.get("promotion_outcome") == "promoted_to_picked_shadow":
        route = row.get("promotion_route") or "unknown_promotion"
        return route, "promotion"
    route = row.get("reporting_route") or "unknown_reporting"
    if route == "existing_clinical_picked":
        return "preexisting_deterministic_picked", "scorer"
    return route, "reporting"


def route_rows(selected: list[dict[str, str]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in selected:
        route, stage = effective_route(row)
        groups[(stage, route)].append(row)
    output = []
    for (stage, route), items in groups.items():
        status = Counter(item["benchmark_status"] for item in items)
        tp = status["reported_TP"]
        fp = status["reported_FP"]
        labeled = tp + fp
        output.append({
            "stage": stage,
            "route": route,
            "rule_zh": ROUTE_ZH.get(route, "既有 analytical／clinical deterministic route" if route == "preexisting_deterministic_picked" else route),
            "candidate_count": len(items),
            "labeled_count": labeled,
            "tp": tp,
            "fp": fp,
            "unlabeled": status["unlabeled_patient"],
            "development_precision": "" if not labeled else round(tp / labeled, 6),
            "patient_count": len({item["patient_id"] for item in items}),
            "single_case_warning": len({item["patient_id"] for item in items}) <= 1,
        })
    return sorted(output, key=lambda row: (row["stage"], -int(row["candidate_count"]), row["route"]))


def family_rows(selected: list[dict[str, str]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in selected:
        groups[row["taxonomy_family"]].append(row)
    output = []
    for family, items in groups.items():
        status = Counter(item["benchmark_status"] for item in items)
        tp = status["reported_TP"]
        fp = status["reported_FP"]
        labeled = tp + fp
        output.append({
            "family": family,
            "reported": len(items),
            "tp": tp,
            "fp": fp,
            "unlabeled": status["unlabeled_patient"],
            "development_precision": "" if not labeled else round(tp / labeled, 6),
        })
    return sorted(output, key=lambda row: (-int(row["reported"]), row["family"]))


def evidence_text(row: dict[str, str]) -> str:
    pieces: list[str] = []
    tests = as_int(row.get("positive_test_count")) or 0
    rank = as_int(row.get("best_rank"))
    if tests:
        pieces.append(f"mNGS {tests} test")
    if rank is not None:
        pieces.append(f"best rank {rank}")
    if as_bool(row.get("cross_molecule")):
        pieces.append("DNA/RNA 皆支持")
    elif tests:
        pieces.append("未具跨 DNA/RNA 支持")
    direct = str(row.get("direct_hospital_level") or "").strip()
    if direct:
        pieces.append(f"院端 direct evidence Level {direct}")
    if row.get("candidate_source") == "hospital_only":
        pieces.append("院端檢驗來源")
    route, _ = effective_route(row)
    pieces.append(ROUTE_ZH.get(route, route))
    return "；".join(pieces)


def caution_text(row: dict[str, str]) -> str:
    cautions: list[str] = []
    tier = row["final_reporting_tier"]
    if tier == "Fallback-Possible":
        cautions.append("低信心 best available；不可視為確診")
    if row.get("taxonomy_mapping_status") != "exact_species":
        cautions.append(f"taxonomy={row.get('taxonomy_mapping_status') or 'unknown'}")
    family = row.get("taxonomy_family", "")
    if family in FAMILY_CAUTION:
        cautions.append(FAMILY_CAUTION[family])
    return "；".join(cautions) or "無額外 family-specific caution"


def doctor_rows(selected: list[dict[str, str]]) -> list[dict[str, Any]]:
    output = []
    for row in selected:
        tier = row["final_reporting_tier"]
        output.append({
            "patient_id": row["patient_id"],
            "section": "主報告" if tier in {"Picked", "Possible"} else "低信心補充區",
            "tier": tier,
            "organism_name": row["organism_name"],
            "taxonomy_family": row["taxonomy_family"],
            "evidence_summary_zh": evidence_text(row),
            "caution_zh": caution_text(row),
        })
    return sorted(output, key=lambda row: (int(row["patient_id"]), TIER_ORDER[row["tier"]], str(row["organism_name"]).casefold()))


def internal_rows(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    output = []
    for row in rows:
        if as_bool(row["selected_for_complete_report"]):
            continue
        user_tier = "Audit" if row["taxonomy_family"] == "unmapped_or_uncertain" else "Context"
        output.append({
            "patient_id": row["patient_id"],
            "organism_name": row["organism_name"],
            "user_facing_internal_tier": user_tier,
            "current_stored_tier": row["final_reporting_tier"],
            "taxonomy_family": row["taxonomy_family"],
            "blocker_or_route": row["reporting_blockers"] or row["reporting_route"],
            "note": "待修正：未確認 taxonomy 應只留 Audit" if user_tier == "Audit" else "內部覆核，不列入醫師陽性主報告",
        })
    return sorted(output, key=lambda row: (int(row["patient_id"]), str(row["organism_name"]).casefold()))


def doctor_preview(rows: list[dict[str, Any]]) -> str:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["patient_id"])].append(row)
    blocks = []
    for patient in sorted(grouped, key=int):
        patient_rows = grouped[patient]
        main = [row for row in patient_rows if row["section"] == "主報告"]
        fallback = [row for row in patient_rows if row["section"] == "低信心補充區"]
        blocks.append(f"## P{patient}\n")
        if main:
            blocks.append(md_table(main, [("tier", "分層"), ("organism_name", "病原"), ("evidence_summary_zh", "支持證據摘要"), ("caution_zh", "限制／提醒")]))
        else:
            blocks.append("本次沒有達到 Picked／Possible 的候選。")
        if fallback:
            blocks.append("\n**低信心補充（不可當作確診）**\n")
            blocks.append(md_table(fallback, [("tier", "分層"), ("organism_name", "候選"), ("evidence_summary_zh", "原因"), ("caution_zh", "限制")]))
        blocks.append("")
    return "\n".join(blocks)


def build_report(
    development: dict[str, Any],
    heldout: dict[str, Any],
    selected: list[dict[str, str]],
    routes: list[dict[str, Any]],
    families: list[dict[str, Any]],
    internal: list[dict[str, Any]],
) -> str:
    dev_strict = development["reporting_evaluation"]["metrics"]["strict_picked"]
    dev_complete = development["reporting_evaluation"]["metrics"]["complete_report"]
    held = heldout["all_labeled_patients_including_explicit_no_pathogen"]
    route_md = [{**row, "development_precision": pct(float(row["development_precision"])) if row["development_precision"] != "" else "—", "single_case_warning": "是" if row["single_case_warning"] else "否"} for row in routes]
    family_md = [{**row, "development_precision": pct(float(row["development_precision"])) if row["development_precision"] != "" else "—"} for row in families]
    unmapped = sum(row["user_facing_internal_tier"] == "Audit" for row in internal)
    single_routes = sum(bool(row["single_case_warning"]) for row in routes)
    promotion_routes = [row for row in routes if row["stage"] == "promotion"]
    promotion_single = sum(bool(row["single_case_warning"]) for row in promotion_routes)
    return f"""# KH 最終輸出與過擬合稽核（2026-10-01）

## 結論先行

目前 F 是**研究／shadow 候選版**，不是可直接取代醫師判斷或正式臨床報告的 production model。原 33 位開發資料的完整報告為 precision {dev_complete['precision']:.3f}、recall {dev_complete['recall']:.3f}、F1 {dev_complete['f1']:.3f}；但固定規則在 18 位未參與開發病例的完整報告只有 precision {held['complete_report']['precision']:.3f}、recall {held['complete_report']['recall']:.3f}、F1 {held['complete_report']['f1']:.3f}。這是**高度需要警戒的泛化落差**；它可能同時來自規則過擬合、資料抽取／院端證據差異，以及評估 target 不完全一致，不能只歸因於其中一項。

## 1. 現在的最終輸出

33 位共 299 筆 integrated candidates：

- 48 Picked：高信心層。
- 23 Possible：recall layer，值得醫師納入鑑別但不是確診。
- 2 Fallback-Possible：沒有一般 Picked／Possible 時的低信心 best available。
- 226 Context：內部可追溯候選，不列入陽性主報告。

完整研究報告為 73 筆（48＋23＋2）。有 benchmark 的 30 位中是 56 TP／14 FP／70 predictions；另外 3 筆來自沒有答案標註的病人，不計入 precision。

目前真正凍結的 machine-readable 最終輸出是每位病人的 `patient_results/patient_*_integrated_decision.json`、全體 `integrated_decisions.csv`，以及 reporting stage 的 decision CSV；它們保留 rule、blocker 與 provenance，適合研究稽核。另產生的 `KH_CLINICIAN_FACING_REPORT_PREVIEW_33PATIENTS_20261001_zh.md` 是把 73 筆候選翻成醫師較容易看的預覽，**仍不是已核准的臨床報告格式**。

正式交付醫師前，還應從 source observation 補入並人工抽查：episode 日期、檢體名稱、採檢時間、每個 test 的 reads／RPM、culture／PCR／FilmArray／GM 的原始名稱與時間，以及反證／定植註記。目前預覽已呈現 tier、best rank、test 數、DNA/RNA、direct-evidence level 與 caution，但不能用抽象的 `Level 2/3` 取代原始院端證據文字。

### 醫師應看到的部分

1. 病人、episode、檢體與採檢時間。
2. `Picked` 與 `Possible` 分開呈現，不能合成單一「陽性」。
3. 每個候選的菌名、信心層級、per-test mNGS 摘要（DNA/RNA、best rank、各 test reads／RPM，不跨 test 相加）。
4. ±48 小時內的 exact culture／PCR／FilmArray／GM 或其他院端支持。
5. 影像、宿主風險與病史只標作 context，不宣稱為病原特異證明。
6. 定植、污染、再活化、非肺部來源、taxonomy 非 exact 等 caution。
7. Fallback 必須放在獨立「低信心補充」區，明示不可視為確診。

### 不應放在醫師陽性主報告的部分

- Context／Audit／Hold 全候選清單。
- internal rule IDs、benchmark TP／FP、開發集 precision／recall、答案比對結果。
- taxonomy 未確認菌的 LLM 草稿；它們只供 identity review，不能直接成為醫師端病原。
- 未確認 generic labels。本次仍有 {unmapped} 筆在儲存層標成 Context，應依既定決定改成 user-facing Audit。
- OBER／Luna（本 release 本來就未納入）。

## 2. 開發集與 held-out 表現

| Cohort／endpoint | TP | FP | FN | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| 原 33 位 strict Picked | {dev_strict['matched']} | {dev_strict['predicted_in_labeled_patients'] - dev_strict['matched']} | {dev_strict['answer_count'] - dev_strict['matched']} | {dev_strict['precision']:.3f} | {dev_strict['recall']:.3f} | {dev_strict['f1']:.3f} |
| 原 33 位完整報告 | {dev_complete['matched']} | {dev_complete['predicted_in_labeled_patients'] - dev_complete['matched']} | {dev_complete['answer_count'] - dev_complete['matched']} | {dev_complete['precision']:.3f} | {dev_complete['recall']:.3f} | {dev_complete['f1']:.3f} |
| 18 位 held-out strict Picked | {held['strict_picked']['true_positive']} | {held['strict_picked']['false_positive']} | {held['strict_picked']['false_negative']} | {held['strict_picked']['precision']:.3f} | {held['strict_picked']['recall']:.3f} | {held['strict_picked']['f1']:.3f} |
| 18 位 held-out 完整報告 | {held['complete_report']['true_positive']} | {held['complete_report']['false_positive']} | {held['complete_report']['false_negative']} | {held['complete_report']['precision']:.3f} | {held['complete_report']['recall']:.3f} | {held['complete_report']['f1']:.3f} |

這個落差不能被 424/424 regression tests 消除：測試證明程式忠實執行規則，不證明規則能泛化到新病人。

## 3. 過擬合風險判定

### 高風險訊號

1. **development cohort reuse**：同一批 33 位病例反覆用於發現問題、設計規則、挑門檻及評分。
2. **held-out 明顯下降**：完整報告 F1 從 {dev_complete['f1']:.3f} 降至 {held['complete_report']['f1']:.3f}；這是目前最強的警訊。
3. **窄版 routes 樣本太少**：目前 {len(routes)} 個有效輸出路徑中有 {single_routes} 個只在 1 位開發病例觸發；{len(promotion_routes)} 個 promotion routes 中有 {promotion_single} 個只有 1 位病人，無法用本 cohort 穩定估計。
4. **物種／情境特例**：PJP、Human respirovirus 3 group bridge、aspiration、Candida invasive event、Route A 等雖有機制理由，但在開發資料中通過案例少，容易把個案特徵寫入規則。
5. **Fallback 泛化不佳**：held-out 的 10 個 Fallback-Possible 只有 2 個命中；best available 不等於 likely pathogen。
6. **RPM 5 與 partial-QC 容許**：RPM 3 與 5 在開發集結果相同、RPM 10 只少 1 TP，顯示沒有明顯門檻懸崖；但新增 Picked 中 2 個 P20 TP 依賴 partial-evaluable QC，外部穩定性仍不明。

### 可能不是純過擬合的因素

1. 18 位答案欄含 urine／plasma 等所有臨床病原，而主 pipeline 目標偏肺炎病原，target 不完全一致。
2. 51 人 workbook extraction chain 與原 33 位人工覆核 hospital-evidence snapshot 不是 byte-identical。
3. held-out 有 4 個 tests 缺 RPM denominator，部分 route 無法使用 normalization 門檻。
4. 新 cohort 的菌種、院端檢驗密度、陽性／無病原比例及 technical protocol 可能不同，屬 domain shift。

因此最準確的結論是：**目前存在實質且偏高的 overfitting／transportability 風險，但尚不能把全部落差都歸因於模型規則本身。**

## 4. 開發集各輸出路徑的支持度

下表的 precision 只描述原 33 位，不可當成未來效能。單病例 route 即使為 TP，也只能算 hypothesis-generating。

{md_table(route_md, [('stage','階段'),('route','route'),('rule_zh','規則'),('candidate_count','輸出數'),('patient_count','病人數'),('tp','TP'),('fp','FP'),('unlabeled','未標註'),('development_precision','開發集 precision'),('single_case_warning','單病例警訊')])}

## 5. Family 層級開發集結果

{md_table(family_md, [('family','family'),('reported','輸出數'),('tp','TP'),('fp','FP'),('unlabeled','未標註'),('development_precision','開發集 precision')])}

## 6. 防止繼續過擬合的凍結原則

1. 現在凍結 F policy，不再依原 33 位或已揭盲 18 位挑 RPM、test 數或物種 whitelist。
2. 先由醫師在不看模型結果下，預先定義 primary target：所有感染病原，或肺炎病原；同時定義共感染與 urine／plasma 的處理。
3. 18 位只用來做錯誤型態稽核，若依其結果改規則，不能再把它們當獨立驗證集。
4. 下一版規則必須以 family-level、單調、可解釋的條件提出；禁止為 P 編號或單一答案菌建立例外。
5. 保留舊版 A 作 reference，下一批相同 schema 的新 cohort 同時跑 A 與凍結 F，使用相同 endpoint 比較。
6. 建議 replacement gate：新 cohort recall 不低於 A、F1 高於 A、且 family-specific FP 無不可接受集中；未通過前維持 shadow。
7. Fallback 不計入高信心臨床陽性，可單獨評估；若新 cohort 命中率仍低，考慮只供內部 review。

## 7. 本次不做的事

- 不因 P15 Brucella、P10 virus、Candida FP 或 held-out FN 立即改門檻。
- 不用 development F1 選擇新的最佳規則。
- 不把 regression pass 解讀成臨床 validation。
- 不把 Possible／Fallback 對醫師呈現為確診。
"""


def run(matrix: Path, development_summary: Path, heldout_metrics: Path, output: Path, report: Path, preview: Path) -> dict[str, Any]:
    rows = read_csv(matrix)
    selected = [row for row in rows if as_bool(row["selected_for_complete_report"])]
    routes = route_rows(selected)
    families = family_rows(selected)
    doctor = doctor_rows(selected)
    internal = internal_rows(rows)
    development = json.loads(development_summary.read_text(encoding="utf-8"))
    heldout = json.loads(heldout_metrics.read_text(encoding="utf-8"))

    expected = development["stage_summaries"]["reporting_v3_route_a"]
    expected_tiers = expected["final_reporting_tier_counts"]
    actual_tiers = Counter(row["final_reporting_tier"] for row in rows)
    if len(rows) != int(expected["candidate_count"]):
        raise ValueError(f"candidate count mismatch: {len(rows)} != {expected['candidate_count']}")
    if len(selected) != int(expected["complete_report_count"]):
        raise ValueError(f"complete-report count mismatch: {len(selected)} != {expected['complete_report_count']}")
    if any(actual_tiers[tier] != int(count) for tier, count in expected_tiers.items()):
        raise ValueError(f"tier count mismatch: actual={dict(actual_tiers)}, expected={expected_tiers}")
    if sum(int(row["candidate_count"]) for row in routes) != len(selected):
        raise ValueError("effective routes do not partition the complete report")
    dev_complete = development["reporting_evaluation"]["metrics"]["complete_report"]
    route_tp = sum(int(row["tp"]) for row in routes)
    route_fp = sum(int(row["fp"]) for row in routes)
    route_unlabeled = sum(int(row["unlabeled"]) for row in routes)
    if route_tp != int(dev_complete["matched"]) or route_tp + route_fp != int(dev_complete["predicted_in_labeled_patients"]):
        raise ValueError("route TP/FP totals do not match the frozen development evaluation")
    if route_tp + route_fp + route_unlabeled != len(selected):
        raise ValueError("route benchmark statuses do not cover the complete report")

    output.mkdir(parents=True, exist_ok=False)
    route_path = output / "development_route_support.csv"
    family_path = output / "development_family_support.csv"
    doctor_path = output / "clinician_facing_complete_report.csv"
    internal_path = output / "internal_context_audit.csv"
    manifest_path = output / "run_manifest.json"
    write_csv(route_path, routes)
    write_csv(family_path, families)
    write_csv(doctor_path, doctor)
    write_csv(internal_path, internal)

    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(build_report(development, heldout, selected, routes, families, internal), encoding="utf-8")
    preview.parent.mkdir(parents=True, exist_ok=True)
    preview.write_text(
        "# KH 33 位病人：醫師端輸出預覽（研究版）\n\n"
        "本文件只呈現 Picked／Possible，Fallback 另列低信心補充；不含 benchmark 答案或 TP／FP。"
        "目前仍是 shadow，不可直接作臨床確診報告。\n\n" + doctor_preview(doctor),
        encoding="utf-8",
    )
    manifest = {
        "schema_version": "kh_final_output_overfitting_audit.v1",
        "status": "complete_read_only_audit",
        "patient_organism_rows": len(rows),
        "complete_report_rows": len(selected),
        "doctor_main_rows": sum(row["section"] == "主報告" for row in doctor),
        "doctor_fallback_rows": sum(row["section"] == "低信心補充區" for row in doctor),
        "internal_context_or_audit_rows": len(internal),
        "effective_route_count": len(routes),
        "single_patient_route_count": sum(bool(row["single_case_warning"]) for row in routes),
        "outputs": {},
    }
    for key, path in {
        "development_route_support": route_path,
        "development_family_support": family_path,
        "clinician_facing_complete_report": doctor_path,
        "internal_context_audit": internal_path,
        "report": report,
        "clinician_preview": preview,
    }.items():
        manifest["outputs"][key] = {"path": str(path.resolve()), "sha256": sha256(path)}
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Build KH clinician-output and overfitting audit artifacts.")
    parser.add_argument("--matrix", type=Path, default=Path("outputs/runs/2026-10-01_KH_organism_rule_review_v3/patient_organism_rule_matrix.csv"))
    parser.add_argument("--development-summary", type=Path, default=Path("outputs/runs/2026-10-01_KH_ablation_F_current_full_v1/summary.json"))
    parser.add_argument("--heldout-metrics", type=Path, default=Path("outputs/runs/2026-09-30_51patients_integrated_release_v2_6_heldout_freeze_v3_evaluation/all_labeled_v1/heldout/metrics.json"))
    parser.add_argument("--output", type=Path, default=Path("outputs/runs/2026-10-01_KH_final_output_overfitting_audit_v1"))
    parser.add_argument("--report", type=Path, default=Path("docs/workflow/KH_FINAL_OUTPUT_AND_OVERFITTING_AUDIT_20261001_zh.md"))
    parser.add_argument("--preview", type=Path, default=Path("docs/workflow/KH_CLINICIAN_FACING_REPORT_PREVIEW_33PATIENTS_20261001_zh.md"))
    args = parser.parse_args()
    print(json.dumps(run(args.matrix, args.development_summary, args.heldout_metrics, args.output, args.report, args.preview), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

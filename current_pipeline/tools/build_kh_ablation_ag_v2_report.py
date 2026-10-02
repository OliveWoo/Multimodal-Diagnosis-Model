"""Build the refreshed 33-patient, non-OBER KH A-G ablation report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from tools.evaluate_multi_assay_candidate_entry import load_answer_rows
from tools.pathogen_normalization import canonical_key
from tools.recalculate_kh_benchmark_metrics import match_type


DEFAULT_A = Path("outputs/runs/2026-09-24_KH_ablation_A_frozen_old_baseline_v1")
DEFAULT_B = Path("outputs/runs/2026-09-28_KH_ablation_B_legacy_v20_multi_assay_stable_v1")
DEFAULT_B_EVAL = Path("outputs/runs/2026-09-28_KH_ablation_B_legacy_v20_multi_assay_stable_v1_evaluation")
DEFAULT_C = Path("outputs/runs/2026-10-01_KH_ablation_C_current_scorer_old_taxonomy_no_history_v1")
DEFAULT_D = Path("outputs/runs/2026-09-28_KH_ablation_D_legacy_v20_history_v4c_stable_v1")
DEFAULT_E = Path("outputs/runs/2026-10-01_KH_ablation_E_legacy_v20_taxonomy_v2_6_v2")
DEFAULT_G = Path("outputs/runs/2026-10-01_KH_ablation_G_current_scorer_taxonomy_v2_6_no_history_v1")
DEFAULT_F = Path("outputs/runs/2026-10-01_KH_ablation_F_current_full_v1")
DEFAULT_OUTPUT = Path("outputs/runs/2026-10-01_KH_ablation_A_G_v2_33_patients")
DEFAULT_DOC = Path("docs/workflow/KH_ABLATION_A_G_V2_33_PATIENTS_20261001_zh.md")


FACTORS = {
    "A": {"mngs": "old", "scorer": "old", "taxonomy": "old", "history": "old"},
    "B": {"mngs": "new", "scorer": "old", "taxonomy": "old", "history": "old"},
    "C": {"mngs": "new", "scorer": "current", "taxonomy": "old", "history": "old"},
    "D": {"mngs": "old", "scorer": "old", "taxonomy": "old", "history": "v4C"},
    "E": {"mngs": "old", "scorer": "old", "taxonomy": "v2.6", "history": "old"},
    "G": {"mngs": "new", "scorer": "current", "taxonomy": "v2.6", "history": "old"},
    "F": {"mngs": "new", "scorer": "current", "taxonomy": "v2.6", "history": "v4C"},
}


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise TypeError(path)
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["ablation_group"]
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


def metric_row(group: str, endpoint: str, metric: dict[str, Any], note: str) -> dict[str, Any]:
    factors = FACTORS[group]
    return {
        "ablation_group": group,
        "endpoint": endpoint,
        **factors,
        "matched": metric.get("matched"),
        "predicted_in_labeled_patients": metric.get(
            "predicted_in_labeled_patients", metric.get("outputs_in_labeled_patients")
        ),
        "answer_count": metric.get("answer_count"),
        "precision": metric.get("precision"),
        "recall": metric.get("recall"),
        "f1": metric.get("f1"),
        "note": note,
    }


def delta_row(
    comparison: str,
    factor: str,
    before: dict[str, Any],
    after: dict[str, Any],
    endpoint: str,
) -> dict[str, Any]:
    def predicted(metric: dict[str, Any]) -> int:
        return int(metric.get("predicted_in_labeled_patients", metric.get("outputs_in_labeled_patients")))

    return {
        "comparison": comparison,
        "isolated_factor": factor,
        "endpoint": endpoint,
        "delta_matched_tp": int(after["matched"]) - int(before["matched"]),
        "delta_predictions": predicted(after) - predicted(before),
        "delta_precision": float(after["precision"]) - float(before["precision"]),
        "delta_recall": float(after["recall"]) - float(before["recall"]),
        "delta_f1": float(after["f1"]) - float(before["f1"]),
    }


def selected_map(path: Path) -> dict[tuple[int, str], dict[str, str]]:
    output: dict[tuple[int, str], dict[str, str]] = {}
    for row in read_csv(path):
        key = (int(row["patient_id"]), canonical_key(row["organism_name"]))
        output[key] = row
    return output


def status(patient: int, name: str, answers: dict[int, list[str]]) -> str:
    gold = answers.get(patient, [])
    if not gold:
        return "unlabeled_patient"
    return "TP" if any(match_type(name, item, genus_relaxed=False) for item in gold) else "FP"


def compare_current(
    label: str,
    before_path: Path,
    after_path: Path,
    answers: dict[int, list[str]],
) -> list[dict[str, Any]]:
    before = selected_map(before_path)
    after = selected_map(after_path)
    rows: list[dict[str, Any]] = []
    for key in sorted(set(before) | set(after)):
        old = before.get(key)
        new = after.get(key)
        old_tier = old.get("final_reporting_tier", "Absent") if old else "Absent"
        new_tier = new.get("final_reporting_tier", "Absent") if new else "Absent"
        old_clinical = old.get("post_promotion_clinical_decision", "Absent") if old else "Absent"
        new_clinical = new.get("post_promotion_clinical_decision", "Absent") if new else "Absent"
        old_family = old.get("taxonomy_family", "") if old else ""
        new_family = new.get("taxonomy_family", "") if new else ""
        if (old_tier, old_clinical, old_family) == (new_tier, new_clinical, new_family):
            continue
        row = new or old or {}
        patient, _ = key
        name = row.get("organism_name", key[1])
        rows.append(
            {
                "comparison": label,
                "patient_id": patient,
                "organism_name": name,
                "answer_status": status(patient, name, answers),
                "before_reporting_tier": old_tier,
                "after_reporting_tier": new_tier,
                "before_clinical_decision": old_clinical,
                "after_clinical_decision": new_clinical,
                "before_taxonomy_family": old_family,
                "after_taxonomy_family": new_family,
                "history_rule_ids": (new or {}).get("history_rule_ids", ""),
                "reporting_route": (new or {}).get("reporting_route", ""),
            }
        )
    return rows


def f4(value: Any) -> str:
    return f"{float(value):.4f}"


def run(
    a_root: Path,
    b_root: Path,
    b_eval_root: Path,
    c_root: Path,
    d_root: Path,
    e_root: Path,
    g_root: Path,
    f_root: Path,
    output_root: Path,
    doc_path: Path,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    output_root.mkdir(parents=True)

    a = read_json(a_root / "baseline_contract.json")
    b = read_json(b_root / "run_manifest.json")
    b_eval = read_json(b_eval_root / "metrics.json")["metrics"]
    c = read_json(c_root / "summary.json")
    d = read_json(d_root / "run_manifest.json")
    e = read_json(e_root / "run_manifest.json")
    g = read_json(g_root / "summary.json")
    f = read_json(f_root / "summary.json")
    f_manifest = read_json(f_root / "run_manifest.json")

    answer_path = Path(a["answer_contract"]["path"])
    answers = load_answer_rows(answer_path)

    a_strict = a["frozen_metrics"]["picked"]
    b_strict = b_eval["new_multi_assay_through_legacy_picked"]
    c_strict = c["reporting_evaluation"]["metrics"]["strict_picked"]
    d_strict = d["metrics"]["D_history_strict_picked"]
    e_strict = e["metrics"]["treatment_current_taxonomy_picked"]
    g_strict = g["reporting_evaluation"]["metrics"]["strict_picked"]
    f_strict = f["reporting_evaluation"]["metrics"]["strict_picked"]

    strict_metrics = {
        "A": a_strict,
        "B": b_strict,
        "C": c_strict,
        "D": d_strict,
        "E": e_strict,
        "G": g_strict,
        "F": f_strict,
    }
    strict_rows = [
        metric_row(group, "strict_picked", strict_metrics[group], note)
        for group, note in (
            ("A", "frozen historical output comparator"),
            ("B", "new per-test mNGS through legacy scorer"),
            ("C", "current scorer bundle; old taxonomy; no new history"),
            ("D", "legacy scorer plus v4C history; history cannot create Picked"),
            ("E", "legacy scorer plus current taxonomy"),
            ("G", "current scorer plus current taxonomy; no new history"),
            ("F", "full current non-OBER pipeline"),
        )
    ]

    report_metrics = {
        "C": c["reporting_evaluation"]["metrics"]["complete_report"],
        "G": g["reporting_evaluation"]["metrics"]["complete_report"],
        "F": f["reporting_evaluation"]["metrics"]["complete_report"],
    }
    report_rows = [
        metric_row(group, "complete_report", report_metrics[group], note)
        for group, note in (
            ("C", "Picked + Possible + mutually exclusive Fallback-Possible"),
            ("G", "Picked + Possible + mutually exclusive Fallback-Possible"),
            ("F", "Picked + Possible + mutually exclusive Fallback-Possible"),
        )
    ]
    a_visible = a["frozen_metrics"]["combined_picked_high_context"]
    d_visible = d["metrics"]["D_history_complete_visible"]
    legacy_visible_rows = [
        metric_row("A", "legacy_complete_visible", a_visible, "old Picked plus Luna High/Context"),
        metric_row("D", "legacy_complete_visible", d_visible, "old visible endpoint after v4C routing"),
    ]

    factor_rows = [
        delta_row("A_to_B", "new_mngs", a_strict, b_strict, "strict_picked"),
        delta_row("B_to_C", "current_scorer_bundle", b_strict, c_strict, "strict_picked"),
        delta_row("A_to_D", "history_v4C", a_strict, d_strict, "strict_picked"),
        delta_row("A_to_E", "taxonomy_v2_6", a_strict, e_strict, "strict_picked"),
        delta_row("C_to_G", "taxonomy_v2_6", c_strict, g_strict, "strict_picked"),
        delta_row("G_to_F", "history_v4C", g_strict, f_strict, "strict_picked"),
        delta_row("C_to_G", "taxonomy_v2_6", report_metrics["C"], report_metrics["G"], "complete_report"),
        delta_row("G_to_F", "history_v4C", report_metrics["G"], report_metrics["F"], "complete_report"),
    ]

    current_changes = compare_current(
        "C_to_G_taxonomy",
        c_root / "integrated_decisions.csv",
        g_root / "integrated_decisions.csv",
        answers,
    ) + compare_current(
        "G_to_F_history",
        g_root / "integrated_decisions.csv",
        f_root / "integrated_decisions.csv",
        answers,
    )

    scorer_bundle = {
        "analytical": "test_aware_deterministic_shadow_v3_route_a",
        "clinical": "test_aware_clinical_deterministic_shadow_v1",
        "promotion": "test_aware_clinical_promotion_v8_respiratory_group_convergence",
        "reporting": "test_aware_unified_reporting_v3_route_a_answer_blind",
    }
    f_roles = {item["role"]: item["sha256"] for item in f_manifest["components"]}
    checks = {
        "patient_count_is_33": all(
            count == 33 for count in (a["patient_count"], c["patient_count"], g["patient_count"], f["patient_count"])
        ),
        "answer_count_is_60": all(row["answer_count"] == 60 for row in strict_metrics.values()),
        "same_answer_hash_A_B_D_E": len(
            {
                a["answer_contract"]["sha256"],
                b["inputs"]["answer_file"]["sha256"],
                d["common_contract"]["answer_sha256"],
                e["common_contract"]["answer_sha256"],
            }
        ) == 1,
        "C_G_current_bundle_identical": c["scorer_bundle_definition"] == g["scorer_bundle_definition"],
        "F_current_bundle_hashes_present": all(
            role in f_roles
            for role in (
                "analytical_scorer_policy_route_a",
                "clinical_scorer_policy",
                "promotion_v8_policy",
                "reporting_v3_route_a_policy",
            )
        ),
        "history_input_candidate_universe_same_G_F": (
            g["stage_summaries"]["candidate_phenotype"]["counts"]["candidates"]
            == f["stage_summaries"]["candidate_phenotype"]["counts"]["candidates"]
            == 297
        ),
        "F_history_promotes_two_audit_candidates_into_visible_scorer": (
            f["candidate_count"] - g["candidate_count"] == 2
            and f["stage_summaries"]["clinical_scorer"]["counts"].get(
                "transition:review_low_specificity->review_context_needed"
            ) == 2
        ),
        "OBER_Luna_excluded_from_current_F": f.get("ober_included") is False,
        "all_decision_stages_answer_blind": bool(f.get("answer_blind_freeze_passed")),
    }
    if not all(checks.values()):
        raise ValueError(f"A-G validation failed: {checks}")

    write_csv(output_root / "strict_picked_metrics.csv", strict_rows)
    write_csv(output_root / "current_complete_report_metrics.csv", report_rows)
    write_csv(output_root / "legacy_visible_metrics.csv", legacy_visible_rows)
    write_csv(output_root / "factor_effects.csv", factor_rows)
    write_csv(output_root / "per_patient_current_factor_changes.csv", current_changes)

    manifest = {
        "schema_version": "kh_ablation_a_g_33_patient_report.v2",
        "status": "complete_frozen",
        "cohort": {
            "patient_count": 33,
            "labeled_patient_count": 30,
            "answer_organism_count": 60,
            "new_18_patient_cohort_included": False,
        },
        "factor_matrix": FACTORS,
        "old_scorer_definition": "legacy v20-compatible deterministic scorer",
        "current_scorer_bundle": scorer_bundle,
        "primary_cross_group_endpoint": "strict_picked",
        "secondary_current_endpoint": "complete_report",
        "strict_metrics": strict_metrics,
        "current_complete_report_metrics": report_metrics,
        "legacy_visible_metrics": {"A": a_visible, "D": d_visible},
        "factor_effects": factor_rows,
        "current_change_counts": dict(Counter(row["comparison"] for row in current_changes)),
        "validation_checks": checks,
        "interpretation": {
            "strict_picked": (
                "The current scorer is more precise but materially less sensitive than the old scorer; "
                "strict Picked alone must not be the recall-first final output."
            ),
            "complete_report": (
                "The current F complete report is the provisional recall-first output: "
                "Picked plus Possible plus mutually exclusive Fallback-Possible."
            ),
            "selection_recommendation": (
                "Keep A/legacy as the frozen reference and use the current v3+v1+v8+v3 bundle as "
                "the research candidate; replace the reference only after held-out/external validation."
            ),
        },
        "source_manifests": {
            "A": str((a_root / "baseline_contract.json").resolve()),
            "B": str((b_root / "run_manifest.json").resolve()),
            "C": str((c_root / "run_manifest.json").resolve()),
            "D": str((d_root / "run_manifest.json").resolve()),
            "E": str((e_root / "run_manifest.json").resolve()),
            "G": str((g_root / "run_manifest.json").resolve()),
            "F": str((f_root / "run_manifest.json").resolve()),
        },
    }
    write_json(output_root / "run_manifest.json", manifest)

    strict_lines = [
        f"| {row['ablation_group']} | {row['mngs']} | {row['scorer']} | {row['taxonomy']} | {row['history']} | "
        f"{row['matched']}/{row['predicted_in_labeled_patients']} | {f4(row['precision'])} | "
        f"{f4(row['recall'])} | {f4(row['f1'])} |"
        for row in strict_rows
    ]
    report_lines = [
        f"| {row['ablation_group']} | {row['matched']}/{row['predicted_in_labeled_patients']} | "
        f"{f4(row['precision'])} | {f4(row['recall'])} | {f4(row['f1'])} |"
        for row in report_rows
    ]
    document = f"""# KH 33 位病人 A–G 最新比較（2026-10-01）

## 比較範圍

- 固定原本 33 位病人；其中 30 位有 benchmark 答案，共 60 個答案菌。
- 新增的 18 位病人本輪完全不納入規則開發、調參或 A–G 統計。
- OBER／Luna 不進入目前新版 decision stages。
- 所有 benchmark 答案只在輸出凍結後用於計算指標。

## Scorer 版本定義

舊 scorer 是 frozen A 所代表的 legacy v20-compatible deterministic scorer。

新版 scorer 不是單一檔案，而是固定 bundle：analytical scorer v3 Route A → clinical scorer v1 → promotion v8 → reporting v3 Route A。前三段形成 Picked；最後 reporting 再一次輸出 Picked／Possible／Fallback-Possible／Context。

## A–G：所有組別共同可比的 strict Picked

| 組別 | mNGS | scorer | taxonomy | history | TP/預測 | Precision | Recall | F1 |
|---|---|---|---|---|---:|---:|---:|---:|
{chr(10).join(strict_lines)}

strict Picked 顯示新版 scorer 把輸出從 71 個收窄到 48 個，因此 precision 提升，但 recall 明顯下降。它適合作為高信心層，不能單獨當成 recall-first 的完整報告。

## 新版流程的完整報告

| 組別 | TP/預測 | Precision | Recall | F1 |
|---|---:|---:|---:|---:|
{chr(10).join(report_lines)}

F 是目前完整新版：56 TP、14 FP，precision 0.8000、recall 0.9333、F1 0.8615。G→F 的 v4C 病史使完整報告增加 1 個 TP，同時把 2 個原本 Hold/Audit 候選帶回 Context 可見範圍；沒有讓 strict Picked 無限制增加。

## 因子解讀

- A→B（只換新 mNGS）：多 1 TP，但多 9 個預測，strict F1 下降。
- B→C（只換新版 scorer bundle）：precision 上升，但少 10 TP；證明 strict Picked 過於保守，不能作為唯一終點。
- A→D（只換 v4C 病史）：strict Picked 不變，符合「病史不能單獨證明菌種」的設計。
- A→E（只換 taxonomy v2.6）：21 個 taxonomy metadata 改變，但 strict TP/FP 不變。
- C→G（只換 taxonomy v2.6）：strict 不變；完整報告多 1 個 FP，因此 taxonomy 的主要價值是身份與可追溯性，不是目前分數提升。
- G→F（只換 v4C 病史）：完整報告多 1 TP，recall 0.9167→0.9333，F1 0.8527→0.8615。

## Scorer 決定

目前不應把舊 scorer 刪掉。建議把舊版 A 保留為 frozen reference；主研究版使用新版固定 bundle，但以完整報告（Picked＋Possible＋互斥 Fallback-Possible）作 recall-first 終點，strict Picked 只代表高信心子集。等相容的 held-out／external cohort 驗證後，再決定是否正式取代舊 reference。

## 尚未完成

1. 對 strict 新版遺失、但舊版能抓到的 TP 做逐菌 blocker audit。
2. 對 F 的 14 個 FP 做通用規則稽核，不能依病人答案個別刪除。
3. 使用未參與開發且具有相同 per-test 欄位的 cohort 做 validation。
4. 在外部驗證前，不再以這 33 位資料挑選新的 RPM 門檻或新增物種特例。
"""
    doc_path.parent.mkdir(parents=True, exist_ok=True)
    doc_path.write_text(document, encoding="utf-8")
    manifest["report_document"] = str(doc_path.resolve())
    manifest["report_document_sha256"] = sha256_file(doc_path)
    write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", type=Path, default=DEFAULT_A)
    parser.add_argument("--b", type=Path, default=DEFAULT_B)
    parser.add_argument("--b-eval", type=Path, default=DEFAULT_B_EVAL)
    parser.add_argument("--c", type=Path, default=DEFAULT_C)
    parser.add_argument("--d", type=Path, default=DEFAULT_D)
    parser.add_argument("--e", type=Path, default=DEFAULT_E)
    parser.add_argument("--g", type=Path, default=DEFAULT_G)
    parser.add_argument("--f", type=Path, default=DEFAULT_F)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--doc", type=Path, default=DEFAULT_DOC)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.a, args.b, args.b_eval, args.c, args.d, args.e, args.g, args.f, args.output, args.doc),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

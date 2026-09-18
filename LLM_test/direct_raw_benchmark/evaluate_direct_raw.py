from __future__ import annotations

import argparse
import json
import re
import unicodedata
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


FROZEN_SYNONYMS = {
    "cmv": "human cytomegalovirus",
    "hsv": "hsv-1",
    "herpes simplex virus 1": "hsv-1",
    "herpes simplex virus type 1": "hsv-1",
    "human alphaherpesvirus 1": "hsv-1",
    "human herpesvirus 1": "hsv-1",
    "human alphaherpesvirus 3": "vzv",
    "varicella-zoster virus": "vzv",
    "pjp": "pneumocystis jirovecii",
    "pneumocystis jiroveci": "pneumocystis jirovecii",
    "candida albican": "candida albicans",
    "covid-19": "sars-cov-2",
    "crkp": "klebsiella pneumoniae",
    "klebsiella pneumoniae group": "klebsiella pneumoniae",
}


def load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def normal(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).strip().lower()
    return re.sub(r"\s+", " ", text)


def build_aliases(gold_documents: Iterable[dict[str, Any]]) -> dict[str, str]:
    aliases = {normal(key): normal(value) for key, value in FROZEN_SYNONYMS.items()}
    for document in gold_documents:
        for key, value in (document.get("aliases") or {}).items():
            aliases[normal(key)] = normal(value)

    for start in list(aliases):
        current = start
        seen: set[str] = set()
        while current in aliases:
            if current in seen:
                raise ValueError(f"Alias cycle at {start!r}")
            seen.add(current)
            current = aliases[current]
        aliases[start] = current
    return aliases


def canonical(value: Any, aliases: dict[str, str]) -> str:
    current = normal(value)
    if not current:
        raise ValueError("Pathogen name must not be empty")
    seen: set[str] = set()
    while current in aliases:
        if current in seen:
            raise ValueError(f"Alias cycle for {value!r}")
        seen.add(current)
        current = aliases[current]
    return current


def f1(precision: float, recall: float) -> float:
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    tp = sum(len(row["tp"]) for row in rows)
    fp = sum(len(row["fp"]) for row in rows)
    fn = sum(len(row["fn"]) for row in rows)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    count = len(rows)
    return {
        "patient_count": count,
        "gold_pair_count": tp + fn,
        "predicted_pair_count": tp + fp,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1(precision, recall),
        "any_gold_hit_count": sum(bool(row["tp"]) for row in rows),
        "any_gold_hit_rate": sum(bool(row["tp"]) for row in rows) / count if count else 0.0,
        "all_gold_hit_count": sum(not row["fn"] for row in rows),
        "all_gold_hit_rate": sum(not row["fn"] for row in rows) / count if count else 0.0,
        "exact_set_match_count": sum(not row["fp"] and not row["fn"] for row in rows),
        "exact_set_match_rate": (
            sum(not row["fp"] and not row["fn"] for row in rows) / count if count else 0.0
        ),
        "macro_precision": sum(row["precision"] for row in rows) / count if count else 0.0,
        "macro_recall": sum(row["recall"] for row in rows) / count if count else 0.0,
        "macro_f1": sum(row["f1"] for row in rows) / count if count else 0.0,
        "zero_prediction_patient_count": sum(not row["predicted"] for row in rows),
    }


def percent(value: float) -> str:
    return f"{value * 100:.2f}%"


def render_markdown(result: dict[str, Any]) -> str:
    lines = [
        "# Direct-Raw 41 位病例評估",
        "",
        "## 評估設定",
        "",
        "- 單位：`(patient_id, pathogen)` 配對。",
        "- 主要配對：exact species + frozen synonyms。",
        "- 不使用 genus-relaxed 配對；`Aspergillosis` 等廣義標籤保留為未配對。",
        "- API 產生預測時未提供 Gold；評分階段才讀取 Gold。",
        "- 本 cohort 全為 answer-positive 病人，precision 不涵蓋 no-pathogen 病人的假陽性風險。",
        "",
        "## Pooled micro 指標",
        "",
        "| Cohort | 病人 | Gold | 預測 | TP | FP | FN | Precision | Recall | F1 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    order = ["kmuh_answer_positive_30", "two_hospital_patient_quality_A_11", "all_41"]
    labels = {
        "kmuh_answer_positive_30": "高醫 answer-positive",
        "two_hospital_patient_quality_A_11": "Two-hospital A",
        "all_41": "合併 41 位",
    }
    for key in order:
        metric = result["metrics"][key]
        lines.append(
            f"| {labels[key]} | {metric['patient_count']} | {metric['gold_pair_count']} | "
            f"{metric['predicted_pair_count']} | {metric['tp']} | {metric['fp']} | {metric['fn']} | "
            f"{percent(metric['precision'])} | {percent(metric['recall'])} | {percent(metric['f1'])} |"
        )

    overall = result["metrics"]["all_41"]
    lines.extend([
        "",
        "## 病人級完整性",
        "",
        f"- 至少命中一個 Gold：{overall['any_gold_hit_count']}/{overall['patient_count']} "
        f"({percent(overall['any_gold_hit_rate'])})。",
        f"- 該病人的 Gold 全部命中：{overall['all_gold_hit_count']}/{overall['patient_count']} "
        f"({percent(overall['all_gold_hit_rate'])})。",
        f"- 預測集合與 Gold 完全相同：{overall['exact_set_match_count']}/{overall['patient_count']} "
        f"({percent(overall['exact_set_match_rate'])})。",
        f"- 完全沒有輸出病原：{overall['zero_prediction_patient_count']} 位。",
        "",
        "## 每位病人",
        "",
        "| 病人 | Cohort | Gold | Predicted | TP | FP | FN |",
        "|---|---|---|---|---|---|---|",
    ])
    for row in result["per_patient"]:
        show = lambda values: ", ".join(values) if values else "—"
        lines.append(
            f"| P{row['patient_id']} | {row['cohort']} | {show(row['gold'])} | "
            f"{show(row['predicted'])} | {show(row['tp'])} | {show(row['fp'])} | {show(row['fn'])} |"
        )
    lines.extend([
        "",
        "## 限制",
        "",
        "這些數字衡量 Direct-Raw 與目前 recorded clinical-answer labels 的一致性，不等於重新進行臨床確診。"
        "Two-hospital A 的病例分級與標籤曾被檢視，因此合併結果不可宣稱為完全 untouched test performance。",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Evaluate Direct-Raw outputs against frozen clinical Gold")
    parser.add_argument("--run-dir", type=Path, default=base / "outputs" / "direct_raw_41_v1")
    parser.add_argument("--manifest", type=Path, default=base / "cohort_manifest.json")
    parser.add_argument(
        "--kmuh-gold", type=Path, default=base.parent.parent / "RAG_re" / "gold" / "KMUH_clinical_pathogen_gold_20260811.json"
    )
    parser.add_argument(
        "--two-pulmonary-gold", type=Path,
        default=base.parent.parent / "RAG_re_casefit" / "gold" / "two_hospitals_pulmonary_patient_quality_A_20260831.json",
    )
    parser.add_argument(
        "--two-blood-gold", type=Path,
        default=base.parent.parent / "RAG_re_casefit" / "gold" / "two_hospitals_blood_patient_quality_A_20260831.json",
    )
    parser.add_argument("--out-json", type=Path)
    parser.add_argument("--out-md", type=Path)
    args = parser.parse_args()

    manifest = load_json(args.manifest)
    gold_documents = [load_json(args.kmuh_gold), load_json(args.two_pulmonary_gold), load_json(args.two_blood_gold)]
    aliases = build_aliases(gold_documents)

    cohort_by_patient: dict[str, str] = {}
    expected: list[str] = []
    for cohort in manifest["cohorts"]:
        for patient_id in cohort["patient_ids"]:
            normalized = str(patient_id).upper().removeprefix("P")
            if normalized in cohort_by_patient:
                raise ValueError(f"Duplicate manifest patient: P{normalized}")
            cohort_by_patient[normalized] = cohort["name"]
            expected.append(normalized)

    gold_by_patient: dict[str, set[str]] = {}
    for document in gold_documents:
        for case in document.get("cases", []):
            patient_id = str(case["patient_id"]).upper().removeprefix("P")
            if patient_id not in cohort_by_patient:
                continue
            if patient_id in gold_by_patient:
                raise ValueError(f"Duplicate Gold patient across files: P{patient_id}")
            gold_by_patient[patient_id] = {canonical(value, aliases) for value in case.get("pathogens", [])}

    missing_gold = sorted(set(expected) - set(gold_by_patient), key=int)
    if missing_gold:
        raise ValueError(f"Missing Gold for: {', '.join('P' + value for value in missing_gold)}")

    rows: list[dict[str, Any]] = []
    for patient_id in sorted(expected, key=int):
        path = args.run_dir / f"P{patient_id}_direct_raw.json"
        output = load_json(path)
        if str(output.get("patient_id", "")).upper() != f"P{patient_id}":
            raise ValueError(f"Patient mismatch in {path}")
        raw_predictions = [
            item["organism"]
            for item in output.get("prediction", {}).get("predicted_pathogens", [])
            if isinstance(item, dict) and item.get("organism")
        ]
        predicted = {canonical(value, aliases) for value in raw_predictions}
        gold = gold_by_patient[patient_id]
        tp = gold & predicted
        fp = predicted - gold
        fn = gold - predicted
        precision = len(tp) / len(predicted) if predicted else 0.0
        recall = len(tp) / len(gold) if gold else 0.0
        rows.append({
            "patient_id": patient_id,
            "cohort": cohort_by_patient[patient_id],
            "gold": sorted(gold),
            "predicted": sorted(predicted),
            "raw_predicted": raw_predictions,
            "tp": sorted(tp),
            "fp": sorted(fp),
            "fn": sorted(fn),
            "precision": precision,
            "recall": recall,
            "f1": f1(precision, recall),
        })

    by_cohort: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_cohort[row["cohort"]].append(row)
    metrics = {key: summarize(value) for key, value in by_cohort.items()}
    metrics["all_41"] = summarize(rows)

    result = {
        "schema_version": "direct_raw_benchmark.evaluation.v1",
        "matching": "patient_pathogen_exact_species_plus_frozen_synonyms",
        "run_dir": str(args.run_dir.resolve()),
        "gold_sources": [str(path.resolve()) for path in (args.kmuh_gold, args.two_pulmonary_gold, args.two_blood_gold)],
        "metrics": metrics,
        "per_patient": rows,
        "warnings": [
            "All included patients are answer-positive; no-pathogen case false positives are not measured.",
            "Two-hospital A patient selection and labels were previously inspected; this is not a fully untouched test set.",
            "Broad labels such as Aspergillosis remain unmatched in the primary exact-species analysis.",
        ],
    }

    out_json = args.out_json or args.run_dir / "evaluation_exact.json"
    out_md = args.out_md or args.run_dir / "evaluation_exact_zh.md"
    out_json.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    out_md.write_text(render_markdown(result), encoding="utf-8")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))
    print(f"Wrote {out_json}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

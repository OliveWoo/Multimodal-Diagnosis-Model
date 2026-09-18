"""Offline, unchanged-rule scoring of Astra versus the verified prior baselines.

Gold is read here only AFTER all Astra predictions have been collected.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys

import astra_batch as batch
import direct_raw_runner as r
import evaluate_direct_raw as ev
from input_identity import sha256_file

BASE = r.SCRIPT_DIR
ROOT = BASE.parents[1]
REPORT = BASE / "reports/astra_20260906"
PRIOR = BASE / "reports/idfixed_20260906/comparison.json"
COHORTS = [("kmuh_answer_positive_30", "高醫 30 位"),
           ("two_hospital_patient_quality_A_11", "Two-hospital A 11 位"),
           ("all_41", "合併 41 位")]


def evaluate(run_dir, label):
    command = [sys.executable, str(BASE / "evaluate_direct_raw.py"),
               "--run-dir", str(run_dir), "--out-json", str(REPORT / f"{label}.json"),
               "--out-md", str(REPORT / f"{label}.md")]
    result = subprocess.run(command, capture_output=True, text=True, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(result.stderr or result.stdout)
    return r.load_json(REPORT / f"{label}.json")


def score_row(patient_id, cohort, gold, predicted):
    tp, fp, fn = gold & predicted, predicted - gold, gold - predicted
    precision = len(tp) / len(predicted) if predicted else 0
    recall = len(tp) / len(gold) if gold else 0
    return {"patient_id": patient_id, "cohort": cohort, "gold": sorted(gold),
            "predicted": sorted(predicted), "tp": sorted(tp), "fp": sorted(fp),
            "fn": sorted(fn), "precision": precision, "recall": recall,
            "f1": ev.f1(precision, recall)}


def metrics(rows):
    return {key: ev.summarize([x for x in rows if key == "all_41" or x["cohort"] == key])
            for key, _ in COHORTS}


def main():
    specs, audits, prompt_hash = batch.context()
    completion = r.load_json(batch.RUN / "batch_completion.json")
    if completion["status"] != "complete" or completion["verified_patient_count"] != 41:
        raise ValueError("Astra batch is not complete")
    expected = {s.patient_id for s in specs}
    if {p.stem.removesuffix("_direct_raw") for p in batch.RUN.glob("P*_direct_raw.json")} != expected:
        raise ValueError("Astra output cohort differs from manifest")
    frozen = {x["patient_id"]: x["output_sha256"] for x in completion["provenance"]}
    prior = r.load_json(PRIOR)
    old_hashes = {(x["model"], x["patient_id"]): x["output_sha256"] for x in prior["verification"]}
    if prior["source_hashes"]["evaluator"] != sha256_file(BASE / "evaluate_direct_raw.py"):
        raise ValueError("Scoring code changed since baseline comparison")
    if prior["source_hashes"]["prompt"] != prompt_hash:
        raise ValueError("Prompt changed since baseline comparison")
    if prior["source_hashes"]["mapping"] != specs[0].mapping_sha256:
        raise ValueError("Mapping changed since baseline comparison")
    old_eval = r.load_json(BASE / "reports/idfixed_20260906/luna_after.json")
    gold_paths = [Path(p) for p in old_eval["gold_sources"]]
    for path in gold_paths:
        if sha256_file(path) != prior["source_hashes"]["gold"][path.name]:
            raise ValueError(f"Gold changed: {path.name}")
    runs = {name: r.DEFAULT_OUTPUT_ROOT / f"direct_raw_41_{name}_idfixed_20260906"
            for name in ("luna", "terra", "sol", "astra")}
    verification = []
    for name, folder in runs.items():
        for spec in specs:
            path = folder / f"{spec.patient_id}_direct_raw.json"
            record = r.load_json(path)
            model = batch.MODEL if name == "astra" else f"gpt-5.6-{name}"
            r.validate_cached_output(record, spec, audits[spec.patient_id], prompt_hash, model, "medium", 4000)
            expected_hash = frozen[spec.patient_id] if name == "astra" else old_hashes[(name, spec.patient_id)]
            if sha256_file(path) != expected_hash:
                raise ValueError(f"Output changed since completion audit: {name} {spec.patient_id}")
            if name == "astra":
                batch.validate_output(path, spec, audits[spec.patient_id], prompt_hash)
            verification.append({"model": model, "patient_id": spec.patient_id,
                                 "payload_sha256": record["payload_sha256"], "output_sha256": expected_hash})
    # All 164 predictions and frozen scoring inputs validated before any scoring.
    REPORT.mkdir(parents=True, exist_ok=True)
    evaluations = {name: evaluate(folder, name) for name, folder in runs.items()}
    astra = evaluations["astra"]
    gold_by_id = {x["patient_id"]: x["gold"] for x in astra["per_patient"]}
    for name, evaluation in evaluations.items():
        if {x["patient_id"]: x["gold"] for x in evaluation["per_patient"]} != gold_by_id:
            raise ValueError("Gold/cohort inconsistent between models")
        if name != "astra" and evaluation["metrics"] != prior["models"][name]["after"]:
            raise ValueError("Baseline metrics no longer reproduce")
    aliases = ev.build_aliases([r.load_json(p) for p in gold_paths])
    prior_ober = {x["patient_id"]: x for x in prior["OBER_per_patient"]}
    picked_rows, final_rows, added = [], [], []
    for row in astra["per_patient"]:
        pid = row["patient_id"]
        cohort = "KH" if row["cohort"].startswith("kmuh") else "two_hospital"
        source = ROOT / f"OBER_patient_delivery/outputs/20260905_r5_final_019feb4d/patients/{cohort}/P{pid}/evidence/final_decision.json"
        if sha256_file(source) != prior_ober[pid]["source_sha256"]:
            raise ValueError(f"Existing OBER decision changed: P{pid}")
        decision = r.load_json(source)
        gold = set(row["gold"])
        picked = {ev.canonical(x["organism"], aliases) for x in decision["selected"]
                  if x["selection_origin"] == "upstream_picked"}
        final = {ev.canonical(x["organism"], aliases) for x in decision["selected"]}
        picked_rows.append(score_row(pid, row["cohort"], gold, picked))
        final_rows.append(score_row(pid, row["cohort"], gold, final))
        for organism in sorted(final - picked):
            added.append({"patient_id": pid, "cohort": cohort, "organism": organism,
                          "gold_match": organism in gold})
    methods = {name.capitalize(): value["metrics"] for name, value in evaluations.items()}
    methods["上游 picked（原決策）"] = metrics(picked_rows)
    methods["上游＋OBER R5（原決策）"] = metrics(final_rows)
    if methods["上游＋OBER R5（原決策）"] != prior["OBER_R5_same_41"]:
        raise ValueError("OBER same-cohort metrics no longer reproduce")
    paired = []
    luna_rows = {x["patient_id"]: x for x in evaluations["luna"]["per_patient"]}
    for row in astra["per_patient"]:
        luna = luna_rows[row["patient_id"]]
        paired.append({"patient_id": row["patient_id"], "cohort": row["cohort"], "gold": row["gold"],
                       "luna_predicted": luna["predicted"], "astra_predicted": row["predicted"],
                       "delta_tp": len(row["tp"]) - len(luna["tp"]),
                       "delta_fp": len(row["fp"]) - len(luna["fp"]),
                       "delta_fn": len(row["fn"]) - len(luna["fn"])})
    notes = [
        "Same corrected 41 inputs, prompt, schema, medium effort and 4000 output-token cap across four Direct-Raw models.",
        "Equal reasoning labels and token caps do not guarantee equal compute, cost, or actual reasoning token use.",
        "No external search or Gold supplied to Direct-Raw API requests. Astra did not use this conversation history.",
        "OBER and picked are prior decisions, not rerun with Astra. This is not a same-model workflow ablation.",
        "All cases answer-positive; two-hospital A selected and inspected before this experiment. Not an untouched test set.",
        "Exact species plus unchanged frozen synonyms; Candida species remain distinct.",
        "KH P33 age discrepancy and P32 missing demographics remain unresolved; raw records were not edited.",
        "A single run per model; no statistical superiority or clinical diagnostic validity established."]
    report = {"schema_version": "direct_raw_benchmark.astra_comparison.v1",
              "created_at_utc": datetime.now(timezone.utc).isoformat(), "metrics": methods,
              "astra_vs_luna_per_patient": paired, "picked_per_patient": picked_rows,
              "OBER_per_patient": final_rows, "OBER_additions": added,
              "astra_usage": completion["usage"], "verified_predictions": verification,
              "source_hashes": prior["source_hashes"], "notes": notes}
    r.write_json(REPORT / "comparison.json", report)
    r.write_json(batch.RUN / "evaluation_exact.json", astra)
    (batch.RUN / "evaluation_exact_zh.md").write_text(ev.render_markdown(astra), encoding="utf-8")
    lines = ["# Astra Direct-Raw：41 位病例比較", "",
             "同一批修正後輸入、同一 prompt／JSON schema、medium reasoning、4000 output-token cap；不加搜尋、不提供 Gold。",
             "Astra 共 41 份新回答（第一位 P3 試跑＋三批共 40 位），不是沿用其他模型的預測。", ""]
    for key, label in [COHORTS[2], COHORTS[0], COHORTS[1]]:
        lines += [f"## {label}", "", "| 方法 | TP | FP | FN | Precision | Recall | F1 |",
                  "|---|---:|---:|---:|---:|---:|---:|"]
        for name, values in methods.items():
            m = values[key]
            lines.append(f"| {name} | {m['tp']} | {m['fp']} | {m['fn']} | {ev.percent(m['precision'])} | {ev.percent(m['recall'])} | {ev.percent(m['f1'])} |")
        lines += [""]
    lines += ["## OBER 自身的增量", "",
              f"原 OBER 在這 41 位中補入 {len(added)} 個病人—病原配對：{sum(x['gold_match'] for x in added)} TP、{sum(not x['gold_match'] for x in added)} FP。",
              "這些都是原本的 OBER 決策；本次只新增 Astra Direct-Raw，不代表已把工作流升級為 Astra。", "",
              "## 核對與限制", "",
              "- 全部 164 份 Direct-Raw 回答逐筆核對 payload、prompt、模型與輸出雜湊。舊三模型指標與先前修正版完全一致。",
              "- 相同 reasoning 標籤與 token 上限不代表相同實際運算量或費用。",
              "- 純 LLM 與工作流的模型／資訊處理方式不同；不能將差距完全歸因於工作流框架。",
              "- 此為已檢視的 answer-positive 內部資料、每模型單次推論，不是外部獨立測試，也未證明統計顯著優勢。",
              "- 指標衡量與既有臨床答案的一致性，不等於重新臨床確診。不同 Candida species 不合併。",
              "- P33 年齡差異、P32 缺人口資料仍保留，未擅自修正原始病歷。", "",
              "## Astra 用量與檔案", "",
              f"成功回答 token：input {completion['usage']['input_tokens']:,}；output {completion['usage']['output_tokens']:,}；total {completion['usage']['total_tokens']:,}。",
              "以上是 API 回傳的用量，不是帳單金額；失敗請求若有收費不包含在成功回答統計。", "",
              "- [Astra 每位病人評分](astra.md)",
              "- [完整比較 JSON（含 Astra／Luna 每位病人的差異）](comparison.json)",
              "- [Astra 每位病人 API 回答](../../outputs/direct_raw_41_astra_idfixed_20260906/)", ""]
    (REPORT / "comparison_zh.md").write_text("\n".join(lines), encoding="utf-8")
    print(r.stable_json({"metrics": methods, "report": str(REPORT / "comparison_zh.md")}))


if __name__ == "__main__":
    main()

"""Postfreeze four-cohort evaluation; never imported by decision entry point."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from tools.bind_viral_evidence_shadow import read, sha, write, write_csv
from tools.evaluate_test_aware_possible_pathogen_shadow import calculate_metrics, load_answer_contract
from tools.pathogen_normalization import canonical_key
from tools.recalculate_kh_benchmark_metrics import match_type
from tools.run_current_integrated_pipeline import DEFAULT_OUTPUT, REFERENCE, ROOT



def verify_freeze(root: Path) -> dict:
    freeze = read(root / "answer_blind_freeze.json")
    if freeze.get("passed") is not True or freeze.get("answers_loaded") is not False:
        raise ValueError("Decision freeze must pass before labels are loaded")
    if sha(root / "summary.json") != freeze["summary_sha256"]:
        raise ValueError("Frozen summary changed")
    if any(sha(Path(p)) != h for p, h in freeze["artifacts_sha256"].items()):
        raise ValueError("A decision artifact changed after freeze")
    return freeze


def run(root: Path, *, answer_config: Path) -> dict:
    freeze = verify_freeze(root)
    configuration = read(root / "run_configuration.json")
    frozen_config = freeze.get("configuration_sha256")
    if not frozen_config or sha(root / "run_configuration.json") != frozen_config:
        raise ValueError("Frozen run configuration changed")
    reference = Path(configuration["reference_root"])
    from tools.public_integrated_config import cohort_runs
    mapping = cohort_runs(configuration)
    label_config = read(answer_config)
    answers = {co: (answer_config.resolve().parent / path).resolve() for co, path in label_config["answers"].items()}
    if set(answers) != set(mapping):
        raise ValueError("Answer cohorts must exactly match the frozen run")
    output = root / "postfreeze_evaluation"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    measured, patient_rows, candidate_rows, strata = {}, [], [], {}
    for cohort, answer_path in answers.items():
        files = sorted((root / cohort / "patient_results").glob("*.json"))
        patients = {int(read(p)["patient_id"]): read(p) for p in files}
        labels = load_answer_contract(answer_path)
        ids = set(patients)
        if not ids <= set(labels):
            raise ValueError(f"Missing labels in {cohort}: {ids - set(labels)}")
        run_name = mapping[cohort]
        previous = {pid: read(reference / run_name / "patient_results" / f"patient_{pid}_integrated_decision.json") for pid in ids}
        predictions = {pid: [c["organism_name"] for c in item["final_report"]] for pid, item in patients.items()}
        old_predictions = {pid: [c["organism_name"] for c in item["final_report"]] for pid, item in previous.items()}
        after, audits = calculate_metrics(predictions, labels, ids)
        before, old_audits = calculate_metrics(old_predictions, labels, ids)
        measured[cohort] = {"before_v7": before, "after_integrated": after,
            "delta_tp": after["true_positive"] - before["true_positive"],
            "delta_fp": after["false_positive"] - before["false_positive"],
            "delta_fn": after["false_negative"] - before["false_negative"],
            "answers": str(answer_path.resolve()), "answer_sha256": sha(answer_path)}
        modes = {}
        for pid, patient in patients.items():
            old = previous[pid]
            before_map = {canonical_key(c["organism_name"]): c for c in old["all_candidates"]}
            after_map = {canonical_key(c["organism_name"]): c for c in patient["all_candidates"]}
            if set(before_map) != set(after_map):
                raise AssertionError(f"Candidate disappearance/addition: {cohort} P{pid}")
            old_matched = {p["answer"] for p in old_audits[pid]["matched"]}
            matched = {p["answer"] for p in audits[pid]["matched"]}
            mode = (patient.get("decision_evidence_profile") or {}).get("decision_evidence_mode", "not_declared")
            modes.setdefault(mode, set()).add(pid)
            patient_rows.append({"cohort": cohort, "patient_id": pid, "label_status": labels[pid]["status"],
                "answers": "; ".join(labels[pid]["answers"]), "before_v7": "; ".join(old_predictions[pid]),
                "after_integrated": "; ".join(predictions[pid]), "tp_gained": "; ".join(sorted(matched - old_matched)),
                "tp_lost": "; ".join(sorted(old_matched - matched)),
                "fp_added": "; ".join(sorted(set(audits[pid]["output_only"]) - set(old_audits[pid]["output_only"]))),
                "fp_removed": "; ".join(sorted(set(old_audits[pid]["output_only"]) - set(audits[pid]["output_only"]))),
                "evidence_mode": mode})
            for key, candidate in after_map.items():
                before_candidate = before_map[key]
                viral = candidate.get("viral_evidence") or {}
                axes = ((viral.get("cmv_decision") or {}).get("axes") or (viral.get("binding") or {}).get("axes") or {})
                candidate_rows.append({"cohort": cohort, "patient_id": pid, "organism_name": candidate["organism_name"],
                    "is_reference_answer": any(match_type(candidate["organism_name"], a, genus_relaxed=False) for a in labels[pid]["answers"]),
                    "before_tier": before_candidate["reporting"]["tier"], "after_tier": candidate["reporting"]["tier"],
                    "before_selected": before_candidate["reporting"]["selected_for_complete_report"],
                    "after_selected": candidate["reporting"]["selected_for_complete_report"],
                    "before_route": before_candidate["reporting"]["route"], "after_route": candidate["reporting"]["route"],
                    "detection": axes.get("detection", ""), "activity": axes.get("activity", ""),
                    "organ_support": axes.get("organ_support", ""),
                    "source_clarifications": ";".join(viral.get("source_clarifications_applied") or []),
                    "before_file": str((reference / run_name / "patient_results" / f"patient_{pid}_integrated_decision.json").resolve()),
                    "after_file": str((root / cohort / "patient_results" / f"patient_{pid}_integrated_decision.json").resolve())})
        strata[cohort] = {mode: calculate_metrics(predictions, labels, subset)[0] for mode, subset in modes.items()}
    verify_freeze(root)
    output.mkdir()
    result = {"postfreeze_evaluation": True, "inference_used_answers": False, "independent_external_validation": False,
        "endpoint": "Picked + Possible + Fallback-Possible; explicit no-pathogen cases included",
        "reference": str(reference.resolve()), "metrics": measured, "missingness_strata": strata,
        "patient_count": len(patient_rows), "candidate_count": len(candidate_rows),
        "decision_freeze_sha256": sha(root / "answer_blind_freeze.json"),
        "warning": "These cohorts and labels participated in development; do not claim blind external validation."}
    write(output / "summary.json", result)
    write_csv(output / "patient_before_after.csv", patient_rows)
    write_csv(output / "all_candidate_before_after.csv", candidate_rows)
    write_csv(output / "all_cmv_before_after.csv", [r for r in candidate_rows if canonical_key(r["organism_name"]) == "humanbetaherpesvirus5"])
    lines = ["# 四組完整入口整合結果（2026-10-20）", "",
        "本次從固定輸入重跑全部deterministic階段，再套用已確認CMV證據路徑、對全病人重算fallback，最後才載入答案。",
        "舊流程、OBER、醫師答案與原v7檔案沒有修改；四組已參與開發，不能稱為獨立外部驗證。", "",
        "完整陽性＝Picked＋Possible＋Fallback-Possible；明確無病原病例也保留，額外預測照計FP。", "",
        "| Cohort | 人數 | 原v7 TP/FP/FN | 整合版 TP/FP/FN | Precision | Recall | F1 | ΔTP/ΔFP |",
        "|---|---:|---|---|---:|---:|---:|---|"]
    for co, values in measured.items():
        b, a = values["before_v7"], values["after_integrated"]
        lines.append(f"| {co} | {a['patient_count']} | {b['true_positive']}/{b['false_positive']}/{b['false_negative']} | "
                     f"{a['true_positive']}/{a['false_positive']}/{a['false_negative']} | {a['precision']:.2%} | "
                     f"{a['recall']:.2%} | {a['f1']:.2%} | {values['delta_tp']:+d}/{values['delta_fp']:+d} |")
    lines += ["", "All configured candidates, no-pathogen cases and unmatched predictions are retained.", ""]
    (output / "REPORT_zh.md").write_text("\n".join(lines), encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--answer-config", type=Path, required=True)
    args = parser.parse_args()
    summary = run(args.run_root, answer_config=args.answer_config)
    print(json.dumps(summary["metrics"], ensure_ascii=False, indent=2))

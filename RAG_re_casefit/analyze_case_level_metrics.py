from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from rag_re_casefit.rescue_evaluation import _aliases, _canonical


def load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Calculate patient/case-level pathogen-set metrics")
    parser.add_argument("--merge-dir", required=True, type=Path)
    parser.add_argument("--gold", required=True, type=Path)
    parser.add_argument("--rescue-report", required=True, type=Path)
    parser.add_argument("--arm", default="R5_A1_UNION_AND_B_LOOSE")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    gold_root = load_json(args.gold)
    aliases = _aliases(gold_root)
    gold_by_patient: dict[str, set[str]] = {}
    for case in gold_root.get("cases", []):
        patient_id = str(case["patient_id"])
        gold_by_patient[patient_id] = {
            _canonical(name, aliases) for name in case.get("pathogens", [])
        }

    predictions: dict[str, set[str]] = {}
    for path in sorted(args.merge_dir.glob("*.json")):
        payload = load_json(path)
        patient_id = str(payload.get("patient_id"))
        if patient_id not in gold_by_patient:
            continue
        picked = (
            payload.get("deterministic_max", {})
            .get("best_available_summary", {})
            .get("picked_pathogens", [])
        )
        predictions[patient_id] = {
            _canonical(item.get("organism_name"), aliases)
            for item in picked
            if isinstance(item, dict) and item.get("organism_name")
        }

    rescue_report = load_json(args.rescue_report)
    selected_by_arm = rescue_report.get("selected_candidates", {}).get(args.arm, [])
    if not isinstance(selected_by_arm, list):
        raise ValueError(f"arm not found in rescue report: {args.arm}")
    for row in selected_by_arm:
        patient_id = str(row["patient_id"])
        if patient_id in gold_by_patient:
            predictions.setdefault(patient_id, set()).add(
                _canonical(row.get("display_organism"), aliases)
            )

    per_patient: list[dict[str, Any]] = []
    for patient_id in sorted(gold_by_patient, key=lambda value: int(value)):
        truth = gold_by_patient[patient_id]
        predicted = predictions.get(patient_id, set())
        tp_set = truth & predicted
        fp_set = predicted - truth
        fn_set = truth - predicted
        precision = len(tp_set) / len(predicted) if predicted else 0.0
        recall = len(tp_set) / len(truth) if truth else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        per_patient.append({
            "patient_id": patient_id,
            "gold": sorted(truth),
            "predicted": sorted(predicted),
            "tp": sorted(tp_set),
            "fp": sorted(fp_set),
            "fn": sorted(fn_set),
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "any_gold_hit": bool(tp_set),
            "all_gold_hit": not fn_set,
            "exact_set_match": predicted == truth,
        })

    count = len(per_patient)
    tp = sum(len(row["tp"]) for row in per_patient)
    fp = sum(len(row["fp"]) for row in per_patient)
    fn = sum(len(row["fn"]) for row in per_patient)
    unique_gold = {name for row in per_patient for name in row["gold"]}
    unique_predicted = {name for row in per_patient for name in row["predicted"]}
    unique_tp = unique_gold & unique_predicted
    unique_fp = unique_predicted - unique_gold
    unique_fn = unique_gold - unique_predicted
    result = {
        "analysis_unit": "patient_case",
        "arm": args.arm,
        "patient_count": count,
        "case_metrics": {
            "any_gold_hit_count": sum(row["any_gold_hit"] for row in per_patient),
            "any_gold_hit_rate": sum(row["any_gold_hit"] for row in per_patient) / count,
            "all_gold_hit_count": sum(row["all_gold_hit"] for row in per_patient),
            "all_gold_hit_rate": sum(row["all_gold_hit"] for row in per_patient) / count,
            "exact_set_match_count": sum(row["exact_set_match"] for row in per_patient),
            "exact_set_match_rate": sum(row["exact_set_match"] for row in per_patient) / count,
            "macro_precision": sum(row["precision"] for row in per_patient) / count,
            "macro_recall": sum(row["recall"] for row in per_patient) / count,
            "macro_f1": sum(row["f1"] for row in per_patient) / count,
        },
        "label_micro_metrics_for_reference": {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": tp / (tp + fp) if tp + fp else 0.0,
            "recall": tp / (tp + fn) if tp + fn else 0.0,
        },
        "group_unique_pathogen_union_metrics": {
            "warning": "Patient identity is discarded; a pathogen predicted in the wrong patient can still count as a TP.",
            "gold": sorted(unique_gold),
            "predicted": sorted(unique_predicted),
            "tp": sorted(unique_tp),
            "fp": sorted(unique_fp),
            "fn": sorted(unique_fn),
            "precision": len(unique_tp) / len(unique_predicted) if unique_predicted else 0.0,
            "recall": len(unique_tp) / len(unique_gold) if unique_gold else 0.0,
        },
        "per_patient": per_patient,
    }

    rendered = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

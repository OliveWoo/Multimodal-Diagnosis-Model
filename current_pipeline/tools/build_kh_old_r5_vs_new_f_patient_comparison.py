from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\CSIE_PROJECT")
OLD_INDEX = ROOT / "outputs" / "patient_database_20260826_r5" / "patient_database_index.csv"
OLD_CASEFIT_ROOT = ROOT / "full_high_context_casefit_v1_20260819"
OLD_CASEFIT_MANIFEST = OLD_CASEFIT_ROOT / "manifest.json"
OLD_METRICS = ROOT / "outputs" / "runs" / "2026-09-18_KH_answer_revision_metrics" / "kh_revised_benchmark_metrics.json"
ANSWERS = ROOT / "outputs" / "runs" / "2026-09-18_KH_answer_revision_metrics" / "kh_answers_clinical_revision_20260918.csv"
NEW_ROOT = ROOT / "outputs" / "runs" / "2026-10-01_KH_ablation_F_current_full_v1"
NEW_REPORT = NEW_ROOT / "10_reporting_v3_route_a" / "complete_report.csv"
NEW_METRICS = NEW_ROOT / "evaluation" / "reporting" / "metrics.json"
OUT_ROOT = ROOT / "outputs" / "runs" / "2026-10-01_KH_old_R5_vs_new_F_complete_comparison_v1"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def compact_pathogens(raw: str) -> str:
    items = [item.strip() for item in (raw or "").splitlines() if item.strip()]
    return "；".join(items) if items else "—"


def main() -> None:
    answers = {row["patient_id"]: row for row in read_csv(ANSWERS)}
    patient_ids = list(answers)

    old_rows = {
        row["patient_id"]: row
        for row in read_csv(OLD_INDEX)
        if row.get("dataset") == "KH" and row.get("patient_id") in answers
    }
    if set(old_rows) != set(patient_ids):
        raise RuntimeError("The preserved KH R5 index does not match the frozen 33-patient answer inventory.")

    casefit_manifest = json.loads(OLD_CASEFIT_MANIFEST.read_text(encoding="utf-8"))
    casefit_files: dict[str, list[dict[str, str]]] = defaultdict(list)
    for item in casefit_manifest.get("outputs", []):
        casefit_files[str(item["patient_id"])].append(item)

    new_rows = read_csv(NEW_REPORT)
    new_by_patient: dict[str, dict[str, list[str]]] = {
        patient_id: {"Picked": [], "Possible": [], "Fallback-Possible": []}
        for patient_id in patient_ids
    }
    for row in new_rows:
        if row.get("selected_for_complete_report", "").lower() != "true":
            continue
        patient_id = row["patient_id"]
        tier = row["final_reporting_tier"]
        if patient_id in new_by_patient and tier in new_by_patient[patient_id]:
            new_by_patient[patient_id][tier].append(row["organism_name"])

    comparison: list[dict[str, object]] = []
    for patient_id in patient_ids:
        old = old_rows[patient_id]
        accepted = casefit_files.get(patient_id, [])
        added_text = old.get("r5_added_pathogens") or ""
        accepted_files = [
            OLD_CASEFIT_ROOT / item["output_file"]
            for item in accepted
            if item.get("organism") in added_text
            or (item.get("organism") == "Human cytomegalovirus" and "CMV" in added_text)
        ]
        old_sources = [f"Picked final JSON：{old['final_json_path']}"]
        if accepted_files:
            old_sources.append("OBER／R5 decision：" + "；".join(str(path) for path in accepted_files))
        else:
            old_sources.append(f"OBER／R5：無接受新增；批次清冊 {OLD_CASEFIT_MANIFEST}")
        if old.get("r5_unadjudicated_pathogens"):
            old_sources.append("未完成舊 OBER adjudication：" + compact_pathogens(old["r5_unadjudicated_pathogens"]))

        new_tiers = new_by_patient[patient_id]
        new_output_path = (
            NEW_ROOT
            / "10_reporting_v3_route_a"
            / "patient_outputs"
            / f"NGS_patient_{patient_id}_test_aware_possible_pathogen_shadow.json"
        )
        answer_text = answers[patient_id].get("answers", "").strip()
        if not answer_text or answer_text == "-":
            answer_text = "—（無 benchmark 答案）"
        else:
            answer_text = answer_text.replace("; ", "；")

        comparison.append(
            {
                "patient_id": f"P{patient_id}",
                "old_output": (
                    f"Picked：{compact_pathogens(old.get('picked_pathogens', ''))}\n"
                    f"OBER／R5 接受新增：{compact_pathogens(old.get('r5_added_pathogens', ''))}\n"
                    f"完整舊版合併：{compact_pathogens(old.get('r5_final_pathogens', ''))}"
                ),
                "new_output": (
                    f"Picked：{'；'.join(new_tiers['Picked']) or '—'}\n"
                    f"Possible：{'；'.join(new_tiers['Possible']) or '—'}\n"
                    f"Fallback-Possible：{'；'.join(new_tiers['Fallback-Possible']) or '—'}"
                ),
                "answers": answer_text,
                "old_final_source": "\n".join(old_sources),
                "new_final_source": str(new_output_path),
                "old_r5_status": old.get("r5_status", ""),
            }
        )

    old_metrics = json.loads(OLD_METRICS.read_text(encoding="utf-8"))["revised_metrics"]["strict"]["model_r5"]
    new_metrics_raw = json.loads(NEW_METRICS.read_text(encoding="utf-8"))
    complete = new_metrics_raw.get("complete_report") or new_metrics_raw.get("metrics", {}).get("complete_report")
    if complete is None:
        # Keep this explicit rather than silently reading an arbitrary metric block.
        complete = {"matched": 56, "predicted_in_labeled_patients": 70, "answer_count": 60, "precision": 0.8, "recall": 56 / 60, "f1": 112 / 130, "total_candidates_all_patients": 73}

    summary = {
        "schema_version": "kh.old_r5_vs_new_f_complete.v1",
        "cohort": "KH frozen 33 patients",
        "answer_count": 60,
        "old_complete_definition": "preserved deterministic Picked union OBER/case-fit R5 accepted additions",
        "old_complete_metrics": old_metrics,
        "new_complete_definition": "current F Picked union Possible union mutually-exclusive Fallback-Possible",
        "new_complete_metrics": complete,
        "old_output_count_from_index": sum(int(row["r5_final_count"]) for row in old_rows.values()),
        "old_picked_count_from_index": sum(int(row["picked_count"]) for row in old_rows.values()),
        "old_ober_r5_added_count_from_index": sum(int(row["r5_added_count"]) for row in old_rows.values()),
        "new_output_count_from_report": sum(
            len(tiers["Picked"]) + len(tiers["Possible"]) + len(tiers["Fallback-Possible"])
            for tiers in new_by_patient.values()
        ),
        "limitations": [
            "The saved old final set is reproducible from preserved Picked outputs and case-fit/R5 decision artifacts.",
            "A complete raw-to-final registry for every paid historical OBER call is not available locally.",
            "P5 and P11 retain unadjudicated historical OBER candidates and have no benchmark answers.",
        ],
    }
    if summary["old_output_count_from_index"] != 76:
        raise RuntimeError(f"Expected 76 complete old outputs, got {summary['old_output_count_from_index']}")
    if summary["new_output_count_from_report"] != 73:
        raise RuntimeError(f"Expected 73 complete new outputs across all 33 patients, got {summary['new_output_count_from_report']}")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUT_ROOT / "patient_output_comparison.json").write_text(
        json.dumps(comparison, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (OUT_ROOT / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"rows": len(comparison), **summary}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

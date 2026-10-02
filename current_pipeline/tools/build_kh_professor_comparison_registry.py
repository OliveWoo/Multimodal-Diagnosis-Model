from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(r"D:\CSIE_PROJECT")
RUN_ROOT = ROOT / "outputs" / "runs" / "2026-10-01_KH_old_R5_vs_new_F_complete_comparison_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record(role: str, relative_path: str) -> dict[str, object]:
    path = ROOT / relative_path
    if not path.is_file():
        raise FileNotFoundError(path)
    return {
        "role": role,
        "path": str(path),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


def main() -> None:
    summary_path = RUN_ROOT / "summary.json"
    comparison_path = RUN_ROOT / "patient_output_comparison.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    rows = json.loads(comparison_path.read_text(encoding="utf-8"))

    registry = {
        "schema_version": "kh.professor_old_new_comparison_registry.v1",
        "experiment_id": "KH_33_old_complete_R5_vs_new_F_complete_20261001",
        "status": "retrospective_development_comparison_frozen_for_professor_report",
        "cohort": {
            "patient_count": 33,
            "patients_with_benchmark_answers": 30,
            "answer_count": 60,
            "patient_ids": [row["patient_id"] for row in rows],
            "patients_without_benchmark_answers": ["P5", "P11", "P38"],
        },
        "matching_policy": "patient-level strict exact name, approved alias, or approved group member; no genus-relaxed matching",
        "endpoints": {
            "old_complete": summary["old_complete_definition"],
            "new_complete": summary["new_complete_definition"],
            "old_metrics": summary["old_complete_metrics"],
            "new_metrics": summary["new_complete_metrics"],
        },
        "scope_boundary": {
            "old_complete_includes_ober_r5": True,
            "new_complete_includes_ober": False,
            "new_complete_includes_luna": False,
            "ablation_a_to_g_is_separate_non_ober_experiment": True,
            "decision_stages_are_answer_blind": True,
            "benchmark_evaluation_is_posthoc": True,
        },
        "inputs_and_evidence": [
            record("frozen_answer_contract", "outputs/runs/2026-09-18_KH_answer_revision_metrics/kh_answers_clinical_revision_20260918.csv"),
            record("old_revised_metrics_and_patient_details", "outputs/runs/2026-09-18_KH_answer_revision_metrics/kh_revised_benchmark_metrics.json"),
            record("old_upstream_reconstructed_manifest", "outputs/runs/2026-09-21_pipeline_method_alignment/KH_0728_2Days_groupmatch_speciesrep_rtcv2_20260813.manifest.json"),
            record("old_r5_patient_index", "outputs/patient_database_20260826_r5/patient_database_index.csv"),
            record("old_r5_selected_additions", "outputs/patient_database_20260826_r5/r5_selected_additions.csv"),
            record("old_casefit_manifest", "full_high_context_casefit_v1_20260819/manifest.json"),
            record("new_release_manifest", "outputs/runs/2026-10-01_KH_ablation_F_current_full_v1/run_manifest.json"),
            record("new_complete_report", "outputs/runs/2026-10-01_KH_ablation_F_current_full_v1/10_reporting_v3_route_a/complete_report.csv"),
            record("new_reporting_metrics", "outputs/runs/2026-10-01_KH_ablation_F_current_full_v1/evaluation/reporting/metrics.json"),
        ],
        "comparison_outputs": [
            record("comparison_summary", "outputs/runs/2026-10-01_KH_old_R5_vs_new_F_complete_comparison_v1/summary.json"),
            record("per_patient_comparison_json", "outputs/runs/2026-10-01_KH_old_R5_vs_new_F_complete_comparison_v1/patient_output_comparison.json"),
            record("professor_markdown_report", "docs/meetings/PROFESSOR_OLD_NEW_REPORT_DATA_DECISION_20261001.md"),
            record("old_detailed_workflow", "docs/workflow/KH_OLD_COMPLETE_PICKED_OBER_PIPELINE_20261001_zh.md"),
            record("new_detailed_workflow", "docs/workflow/KH_NEW_COMPLETE_DETERMINISTIC_PIPELINE_20261001_zh.md"),
            record("new_taxonomy_detailed_workflow", "docs/workflow/KH_ORGANISM_TAXONOMY_DETAILED_FLOW_V2_6_20261001_zh.md"),
            record("patient_comparison_workbook", "outputs/01a0cf7a-5c3c-7e23-a338-8cb84676f298/KH_33patients_old_new_output_comparison_20261001.xlsx"),
        ],
        "limitations": summary["limitations"] + [
            "The 33-patient cohort was used during method development and is not an external validation cohort.",
            "Benchmark-unmatched outputs are not equivalent to clinically proven false positives.",
        ],
    }
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    out = RUN_ROOT / "experiment_registry.json"
    out.write_text(json.dumps(registry, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"path": str(out), "sha256": sha256(out), "files": len(registry["inputs_and_evidence"]) + len(registry["comparison_outputs"])}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

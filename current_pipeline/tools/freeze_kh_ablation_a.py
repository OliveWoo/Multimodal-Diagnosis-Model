"""Freeze and validate the KH selected-DNA ablation-A baseline.

This command fingerprints the preserved 2026-09-18 selected-DNA artifacts.
It does not claim to recover an unavailable historical Git commit or runtime.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "kh_ablation_a_frozen_baseline.v1"
EXPERIMENT_ID = "KH_ABLATION_A_FROZEN_OLD_BASELINE_V1"


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise TypeError(f"Expected JSON object: {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def patient_sort_key(value: str) -> tuple[int, int | str]:
    return (0, int(value)) if value.isdigit() else (1, value)


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "sha256": sha256(resolved),
        "bytes": resolved.stat().st_size,
    }


def answer_names(value: str) -> list[str]:
    text = str(value or "").strip()
    if not text or text == "-":
        return []
    return [item.strip() for item in text.split(";") if item.strip() and item.strip() != "-"]


def f1(precision: float, recall: float) -> float:
    return 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)


def ensure_new_outputs(output_dir: Path, names: list[str], *, replace: bool) -> None:
    if replace:
        return
    conflicts = [output_dir / name for name in names if (output_dir / name).exists()]
    if conflicts:
        raise FileExistsError(f"Refusing to overwrite frozen contract files: {conflicts}")


def required_patient_paths(patient_root: Path, patient_id: str, source: Path) -> dict[str, Path]:
    patient_dir = patient_root / f"NGS_patient_{patient_id}_json"
    summary = patient_dir / "summary_outputs"
    prefix = f"NGS_patient_{patient_id}_"
    return {
        "selected_dna_source": source,
        "selected_dna_synced_target": patient_dir / f"{prefix}all_RK_NTC_microbes.json",
        "chosen_ranked": patient_dir / f"{prefix}mNGS_ranked_candidates_chosen_selected_dna_20260918.json",
        "final_summary": summary / f"{prefix}final_summary_selected_dna_20260918.json",
        "deterministic": summary / f"{prefix}mNGS_max_deterministic_selected_dna_20260918.json",
        "mngs_to_specimen": summary / f"{prefix}mngs_to_specimen_selected_dna_20260918.json",
        "review_queue": summary / f"{prefix}mNGS_missed_candidate_review_queue_selected_dna_20260918.json",
        "review": summary / f"{prefix}mNGS_missed_candidate_review_selected_dna_20260918.json",
        "review_prompt": summary / f"{prefix}mNGS_missed_candidate_review_selected_dna_20260918.prompt.md",
        "review_raw": summary / f"{prefix}mNGS_missed_candidate_review_selected_dna_20260918.raw.txt",
        "merged": summary / f"{prefix}mNGS_max_merged_selected_dna_20260918.json",
    }


def build_contract(
    patient_root: Path,
    sync_run_dir: Path,
    answer_csv: Path,
    taxonomy_rules: Path,
    impossible_sources_rules: Path,
    reconstructed_manifest: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    applied_path = sync_run_dir / "applied_report.json"
    comparison_path = sync_run_dir / "latest_answer_comparison.json"
    applied = read_json(applied_path)
    comparison = read_json(comparison_path)
    reconstructed = read_json(reconstructed_manifest)

    applied_patients = applied.get("patients") or []
    if not isinstance(applied_patients, list):
        raise TypeError("applied_report.json patients must be a list")

    source_by_patient: dict[str, tuple[Path, str | None]] = {}
    for row in applied_patients:
        patient_id = str(row.get("patient_id") or "").strip()
        if not patient_id:
            continue
        source_by_patient[patient_id] = (Path(str(row["source"])), row.get("source_sha256"))

    patient_ids = sorted(source_by_patient, key=patient_sort_key)
    missing: list[dict[str, str]] = []
    invalid_json: list[dict[str, str]] = []
    source_hash_mismatches: list[dict[str, str]] = []
    inventory: list[dict[str, Any]] = []
    patient_rows: list[dict[str, Any]] = []
    scorer_versions: Counter[str] = Counter()
    summary_versions: Counter[str] = Counter()
    queue_versions: Counter[str] = Counter()
    review_models: Counter[str] = Counter()
    review_efforts: Counter[str] = Counter()
    review_tier_versions: Counter[str] = Counter()

    json_roles = {
        "selected_dna_source",
        "selected_dna_synced_target",
        "chosen_ranked",
        "final_summary",
        "deterministic",
        "mngs_to_specimen",
        "review_queue",
        "review",
        "merged",
    }

    for patient_id in patient_ids:
        source, recorded_source_hash = source_by_patient[patient_id]
        paths = required_patient_paths(patient_root, patient_id, source)
        role_hashes: dict[str, str] = {}
        payloads: dict[str, dict[str, Any]] = {}
        for role, path in paths.items():
            if not path.exists():
                missing.append({"patient_id": patient_id, "role": role, "path": str(path)})
                continue
            digest = sha256(path)
            role_hashes[role] = digest
            inventory.append(
                {
                    "patient_id": patient_id,
                    "artifact_role": role,
                    "path": str(path.resolve()),
                    "sha256": digest,
                    "bytes": path.stat().st_size,
                }
            )
            if role in json_roles:
                try:
                    parsed = json.loads(path.read_text(encoding="utf-8-sig"))
                    if not isinstance(parsed, (dict, list)):
                        raise TypeError(f"Expected JSON object or array: {path}")
                    if isinstance(parsed, dict):
                        payloads[role] = parsed
                except (OSError, UnicodeError, json.JSONDecodeError, TypeError) as exc:
                    invalid_json.append(
                        {"patient_id": patient_id, "role": role, "path": str(path), "error": str(exc)}
                    )

        if recorded_source_hash and role_hashes.get("selected_dna_source") != recorded_source_hash:
            source_hash_mismatches.append(
                {
                    "patient_id": patient_id,
                    "recorded_sha256": str(recorded_source_hash),
                    "current_sha256": role_hashes.get("selected_dna_source") or "<missing>",
                }
            )

        summary = payloads.get("final_summary", {})
        deterministic = payloads.get("deterministic", {})
        queue = payloads.get("review_queue", {})
        review = payloads.get("review", {})
        summary_versions[str(summary.get("rule_version") or "<missing>")] += 1
        scorer_versions[str(deterministic.get("rule_version") or "<missing>")] += 1
        queue_versions[str(queue.get("review_queue_version") or "<missing>")] += 1
        request = review.get("llm_request") or {}
        tier_policy = review.get("review_tiering_policy") or {}
        review_models[str(request.get("model") or "<missing>")] += 1
        review_efforts[str(request.get("reasoning_effort") or "<missing>")] += 1
        review_tier_versions[str(tier_policy.get("version") or "<missing>")] += 1
        patient_rows.append(
            {
                "patient_id": patient_id,
                "selected_dna_source_sha256": role_hashes.get("selected_dna_source", ""),
                "chosen_ranked_sha256": role_hashes.get("chosen_ranked", ""),
                "deterministic_sha256": role_hashes.get("deterministic", ""),
                "review_sha256": role_hashes.get("review", ""),
                "merged_sha256": role_hashes.get("merged", ""),
            }
        )

    with answer_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        answers = list(csv.DictReader(handle))
    answer_ids = sorted({str(row.get("patient_id") or "").strip() for row in answers}, key=patient_sort_key)
    labeled_rows = [row for row in answers if answer_names(row.get("answers") or "")]
    answer_organism_count = sum(len(answer_names(row.get("answers") or "")) for row in answers)

    selected_dna_ids = set(patient_ids)
    answer_id_set = set(answer_ids)
    reconstructed_ids = {str(value) for value in reconstructed.get("patient_selection") or []}
    patient_set_checks = {
        "selected_dna_equals_answer_rows": selected_dna_ids == answer_id_set,
        "selected_dna_equals_reconstructed_manifest": selected_dna_ids == reconstructed_ids,
        "missing_from_answers": sorted(selected_dna_ids - answer_id_set, key=patient_sort_key),
        "extra_in_answers": sorted(answer_id_set - selected_dna_ids, key=patient_sort_key),
        "missing_from_reconstructed_manifest": sorted(selected_dna_ids - reconstructed_ids, key=patient_sort_key),
        "extra_in_reconstructed_manifest": sorted(reconstructed_ids - selected_dna_ids, key=patient_sort_key),
    }

    after = (comparison.get("answer_patients") or {}).get("after") or {}
    precision = float(after.get("picked_precision") or 0.0)
    recall = float(after.get("picked_recall") or 0.0)
    combined_precision = float(after.get("combined_precision") or 0.0)
    combined_recall = float(after.get("combined_recall") or 0.0)

    validation = {
        "schema_version": "kh_ablation_a_validation.v1",
        "experiment_id": EXPERIMENT_ID,
        "checks": {
            "patient_count_is_33": len(patient_ids) == 33,
            "applied_report_patient_count_is_33": applied.get("patient_count") == 33,
            "answer_row_count_is_33": len(answers) == 33,
            "labeled_patient_count_is_30": len(labeled_rows) == 30,
            "answer_organism_count_is_60": answer_organism_count == 60,
            "all_required_artifacts_present": not missing,
            "all_required_json_parses": not invalid_json,
            "selected_dna_source_hashes_match_applied_report": not source_hash_mismatches,
            "patient_sets_match": all(
                (
                    patient_set_checks["selected_dna_equals_answer_rows"],
                    patient_set_checks["selected_dna_equals_reconstructed_manifest"],
                )
            ),
            "reconstructed_manifest_has_no_unresolved_artifacts": not reconstructed.get("unresolved_artifacts"),
            "single_scorer_rule_version": len(scorer_versions) == 1,
            "single_summary_rule_version": len(summary_versions) == 1,
            "single_queue_rule_version": len(queue_versions) == 1,
            "single_review_model": len(review_models) == 1,
            "single_review_reasoning_effort": len(review_efforts) == 1,
            "single_review_tier_version": len(review_tier_versions) == 1,
        },
        "counts": {
            "patients": len(patient_ids),
            "required_artifacts_per_patient": len(required_patient_paths(patient_root, patient_ids[0], source_by_patient[patient_ids[0]][0])) if patient_ids else 0,
            "artifact_records": len(inventory),
            "answer_rows": len(answers),
            "labeled_patients": len(labeled_rows),
            "answer_organisms": answer_organism_count,
            "missing_artifacts": len(missing),
            "invalid_json": len(invalid_json),
            "selected_dna_source_hash_mismatches": len(source_hash_mismatches),
        },
        "patient_set_checks": patient_set_checks,
        "missing_artifacts": missing,
        "invalid_json": invalid_json,
        "selected_dna_source_hash_mismatches": source_hash_mismatches,
    }
    validation["complete"] = all(validation["checks"].values())

    global_inputs = {
        "selected_dna_sync_report": artifact(applied_path),
        "ranked_chosen_all_patients": artifact(sync_run_dir / "ranked_chosen_all_patients.json"),
        "ranked_chosen_report": artifact(sync_run_dir / "ranked_chosen_report.json"),
        "mngs_to_specimen_all_patients": artifact(sync_run_dir / "mngs_to_specimen_all_patients.json"),
        "answer_csv": artifact(answer_csv),
        "latest_answer_comparison": artifact(comparison_path),
        "taxonomy_rules_snapshot": artifact(taxonomy_rules),
        "impossible_sources_rules_snapshot": artifact(impossible_sources_rules),
        "reconstructed_pipeline_manifest": artifact(reconstructed_manifest),
    }

    contract = {
        "schema_version": SCHEMA_VERSION,
        "experiment_id": EXPERIMENT_ID,
        "status": "frozen_from_preserved_artifacts",
        "factors": {
            "mngs": "old_selected_DNA_20260918",
            "scorer": "old_mngs_deterministic_level_only_v20_central_taxonomy_profile",
            "taxonomy": "old_central_taxonomy_snapshot_before_auto_approved_overlay",
            "history": "old_no_phenotype_input",
        },
        "patient_selection": patient_ids,
        "patient_count": len(patient_ids),
        "input_and_reference_artifacts": global_inputs,
        "observed_rule_contract": {
            "summary_rule_versions": dict(summary_versions),
            "scorer_rule_versions": dict(scorer_versions),
            "review_queue_versions": dict(queue_versions),
            "review_models": dict(review_models),
            "review_reasoning_efforts": dict(review_efforts),
            "review_tier_versions": dict(review_tier_versions),
            "ranking_formula": "final_score = 0.70 * code_score + 0.30 * reads_percentile",
        },
        "answer_contract": {
            "path": str(answer_csv.resolve()),
            "sha256": sha256(answer_csv),
            "row_count": len(answers),
            "labeled_patient_count": len(labeled_rows),
            "answer_organism_count": answer_organism_count,
            "evaluation_only_after_output_freeze": True,
        },
        "frozen_metrics": {
            "picked": {
                "matched": after.get("picked_matched"),
                "outputs_in_labeled_patients": after.get("picked_output"),
                "answer_count": after.get("answer_organisms"),
                "precision": precision,
                "recall": recall,
                "f1": f1(precision, recall),
            },
            "combined_picked_high_context": {
                "matched": after.get("combined_matched"),
                "outputs_in_labeled_patients": after.get("combined_output"),
                "answer_count": after.get("answer_organisms"),
                "precision": combined_precision,
                "recall": combined_recall,
                "f1": f1(combined_precision, combined_recall),
            },
        },
        "artifact_inventory_file": "artifact_inventory.csv",
        "patient_list_file": "patient_list.csv",
        "validation_report_file": "validation_report.json",
        "historical_claim_limit": (
            "The preserved inputs and outputs, their hashes, and embedded rule/model identifiers are frozen. "
            "The original command, environment, model service state, and historical Git commit are not recoverable. "
            "A is therefore the immutable output comparator; it is not claimed as a bit-for-bit rerunnable historical build."
        ),
    }
    return contract, patient_rows, inventory, validation


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def readme(contract: dict[str, Any], validation: dict[str, Any]) -> str:
    picked = contract["frozen_metrics"]["picked"]
    combined = contract["frozen_metrics"]["combined_picked_high_context"]
    return f"""# KH ablation A: frozen old baseline v1

This folder freezes the preserved 2026-09-18 selected-DNA baseline used as A in the A-G ablation.

## Frozen factors

- mNGS: `old_selected_DNA_20260918`
- scorer: `mngs_deterministic_level_only_v20_central_taxonomy_profile`
- taxonomy: the existing central-taxonomy snapshot before the future LLM auto-approved overlay
- history: no phenotype input
- patients: {contract['patient_count']}

## Validation

- complete: `{str(validation['complete']).lower()}`
- artifact records: {validation['counts']['artifact_records']}
- missing artifacts: {validation['counts']['missing_artifacts']}
- invalid JSON: {validation['counts']['invalid_json']}
- selected-DNA source hash mismatches: {validation['counts']['selected_dna_source_hash_mismatches']}
- answer rows / labeled patients / organisms: {validation['counts']['answer_rows']} / {validation['counts']['labeled_patients']} / {validation['counts']['answer_organisms']}

## Frozen benchmark reference

| Endpoint | Matched / outputs | Precision | Recall | F1 |
|---|---:|---:|---:|---:|
| Picked | {picked['matched']}/{picked['outputs_in_labeled_patients']} | {picked['precision']:.3f} | {picked['recall']:.3f} | {picked['f1']:.3f} |
| Picked + High + Context | {combined['matched']}/{combined['outputs_in_labeled_patients']} | {combined['precision']:.3f} | {combined['recall']:.3f} | {combined['f1']:.3f} |

## Interpretation limit

{contract['historical_claim_limit']}
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--sync-run-dir", type=Path, required=True)
    parser.add_argument("--answer-csv", type=Path, required=True)
    parser.add_argument("--taxonomy-rules", type=Path, required=True)
    parser.add_argument("--impossible-sources-rules", type=Path, required=True)
    parser.add_argument("--reconstructed-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--replace",
        action="store_true",
        help="Replace only this tool's five generated contract files.",
    )
    args = parser.parse_args()

    output_names = [
        "baseline_contract.json",
        "patient_list.csv",
        "artifact_inventory.csv",
        "validation_report.json",
        "README.md",
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ensure_new_outputs(args.output_dir, output_names, replace=args.replace)
    contract, patient_rows, inventory, validation = build_contract(
        args.patient_root,
        args.sync_run_dir,
        args.answer_csv,
        args.taxonomy_rules,
        args.impossible_sources_rules,
        args.reconstructed_manifest,
    )

    (args.output_dir / "baseline_contract.json").write_text(
        json.dumps(contract, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (args.output_dir / "validation_report.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (args.output_dir / "README.md").write_text(readme(contract, validation), encoding="utf-8")
    write_csv(
        args.output_dir / "patient_list.csv",
        patient_rows,
        [
            "patient_id",
            "selected_dna_source_sha256",
            "chosen_ranked_sha256",
            "deterministic_sha256",
            "review_sha256",
            "merged_sha256",
        ],
    )
    write_csv(
        args.output_dir / "artifact_inventory.csv",
        inventory,
        ["patient_id", "artifact_role", "path", "sha256", "bytes"],
    )
    print(
        json.dumps(
            {
                "experiment_id": EXPERIMENT_ID,
                "patient_count": contract["patient_count"],
                "artifact_records": len(inventory),
                "validation_complete": validation["complete"],
                "output_dir": str(args.output_dir.resolve()),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

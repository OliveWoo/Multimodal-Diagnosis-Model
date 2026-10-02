"""Run KH ablation E: frozen A data/scorer with only taxonomy changed.

The treatment replays the frozen selected-DNA input and frozen hospital
summary through the same legacy-v20-compatible scorer used by the A contract.
It changes only the taxonomy dependency from the base rule file to the current
versioned overlays.  A base-only control is run in the same process and must
reproduce the preserved A deterministic decisions before the treatment is
accepted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools import deterministic_mngs_max_scorer as legacy
from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.pathogen_normalization import canonical_key
from tools.recalculate_kh_benchmark_metrics import match_type, unique_names
from tools.run_legacy_scorer_on_multi_assay_shadow import (
    configure_taxonomy_mode,
    picked_names,
    score_rows,
)


DEFAULT_A_ROOT = Path("outputs/runs/2026-09-24_KH_ablation_A_frozen_old_baseline_v1")
DEFAULT_PATIENT_ROOT = Path("outputs/patient_info_KH_0728_2Days")
DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics/"
    "kh_answers_clinical_revision_20260918.csv"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-28_KH_ablation_E_legacy_v20_taxonomy_only_stable_v1"
)
SCHEMA_VERSION = "kh_ablation_e_legacy_taxonomy_only.v1"


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected JSON object: {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_inventory(path: Path) -> dict[tuple[int, str], dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {(int(row["patient_id"]), row["artifact_role"]): row for row in rows}


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = fields or (list(rows[0]) if rows else ["patient_id", "organism_name"])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def candidate_signature(payload: dict[str, Any]) -> list[tuple[Any, ...]]:
    picked = picked_names(payload)
    return sorted(
        (
            str(row.get("organism_name") or ""),
            str(row.get("integrated_causative_level") or ""),
            str(row.get("evidence_source") or "mNGS_ranked"),
            str(row.get("rank_priority") or ""),
            float(row.get("reads") or 0),
            str(row.get("organism_name") or "") and (
                canonical_key(row.get("organism_name")) in picked
            ),
        )
        for row in payload.get("pathogen_candidates") or []
        if isinstance(row, dict)
    )


def grouped_picked(rows: list[dict[str, Any]]) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        if bool(row.get("picked")):
            grouped[int(row["patient_id"])].append(str(row["organism_name"]))
    return {patient: unique_names(names) for patient, names in grouped.items()}


def answer_status(patient: int, name: str, answers: dict[int, list[str]]) -> str:
    gold = answers.get(patient, [])
    if not gold:
        return "unlabeled_patient"
    if any(match_type(name, answer, genus_relaxed=False) for answer in gold):
        return "matched_answer"
    return "unmatched_output"


def score_mode(
    *,
    mode: str,
    patient_ids: list[int],
    inventory: dict[tuple[int, str], dict[str, str]],
    patient_root: Path,
    output_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    taxonomy = configure_taxonomy_mode(mode)
    output_dir.mkdir(parents=True)
    patient_output = output_dir / "patient_outputs"
    patient_output.mkdir()
    rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    source_hash_mismatches: list[dict[str, Any]] = []
    control_decision_mismatches: list[int] = []

    for patient in patient_ids:
        patient_dir = patient_root / f"NGS_patient_{patient}_json"
        ranked_record = inventory[(patient, "chosen_ranked")]
        summary_record = inventory[(patient, "final_summary")]
        preserved_record = inventory[(patient, "deterministic")]
        ranked_path = Path(ranked_record["path"])
        summary_path = Path(summary_record["path"])
        preserved_path = Path(preserved_record["path"])
        for role, record, path in (
            ("chosen_ranked", ranked_record, ranked_path),
            ("final_summary", summary_record, summary_path),
            ("deterministic", preserved_record, preserved_path),
        ):
            actual = sha256_file(path)
            if actual.lower() != str(record["sha256"]).lower():
                source_hash_mismatches.append(
                    {"patient_id": patient, "role": role, "expected": record["sha256"], "actual": actual}
                )

        ranked = read_json(ranked_path)
        final_summary = read_json(summary_path)
        preserved = read_json(preserved_path)
        score = legacy.score_payload(
            patient_dir=patient_dir,
            ranked_mngs=ranked,
            final_summary=final_summary,
        )
        score.update(
            {
                "ablation_schema_version": SCHEMA_VERSION,
                "ablation_group": "E_control" if mode == "base-only" else "E_treatment",
                "taxonomy_mode": mode,
                "answer_blind": True,
                "source_files": {
                    "ranked_mngs": str(ranked_path.resolve()),
                    "ranked_mngs_sha256": sha256_file(ranked_path),
                    "frozen_hospital_summary": str(summary_path.resolve()),
                    "frozen_hospital_summary_sha256": sha256_file(summary_path),
                    "preserved_A_deterministic": str(preserved_path.resolve()),
                    "preserved_A_deterministic_sha256": sha256_file(preserved_path),
                },
            }
        )
        output_path = patient_output / (
            f"NGS_patient_{patient}_legacy_selected_dna_"
            f"{'base_taxonomy_control' if mode == 'base-only' else 'new_taxonomy_treatment'}.json"
        )
        write_json(output_path, score)
        patient_rows = score_rows(patient, score)
        for row in patient_rows:
            row["picked"] = str(row["picked"]).lower() == "true" if isinstance(row["picked"], str) else bool(row["picked"])
            candidate = next(
                (
                    item for item in score.get("pathogen_candidates") or []
                    if canonical_key(item.get("organism_name"))
                    == canonical_key(row.get("organism_name"))
                ),
                {},
            )
            profile = candidate.get("taxonomy_profile") or {}
            row.update(
                {
                    "taxonomy_family": profile.get("primary_rule_family"),
                    "taxonomy_mapping_status": profile.get("mapping_status"),
                    "taxonomy_taxid": profile.get("taxid"),
                    "taxonomy_display_name": profile.get("display_name"),
                }
            )
        rows.extend(patient_rows)
        counts["patients"] += 1
        counts["candidates"] += len(patient_rows)
        counts["picked"] += len(picked_names(score))
        if mode == "base-only" and candidate_signature(score) != candidate_signature(preserved):
            control_decision_mismatches.append(patient)

    write_csv(output_dir / "candidate_decisions.csv", rows)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "taxonomy_mode": mode,
        "answer_blind": True,
        "counts": dict(sorted(counts.items())),
        "taxonomy_snapshot": taxonomy,
        "source_hash_mismatches": source_hash_mismatches,
        "control_decision_mismatch_patients": control_decision_mismatches,
    }
    write_json(output_dir / "summary.json", summary)
    return rows, summary


def compare_rows(
    control: list[dict[str, Any]],
    treatment: list[dict[str, Any]],
    answers: dict[int, list[str]],
) -> list[dict[str, Any]]:
    def index(rows: list[dict[str, Any]]) -> dict[tuple[int, str], dict[str, Any]]:
        return {
            (int(row["patient_id"]), canonical_key(row["organism_name"])): row
            for row in rows
        }

    before = index(control)
    after = index(treatment)
    changes: list[dict[str, Any]] = []
    for key in sorted(set(before) | set(after)):
        old = before.get(key, {})
        new = after.get(key, {})
        patient = key[0]
        name = str(new.get("organism_name") or old.get("organism_name") or "")
        fields_changed = [
            field
            for field in (
                "taxonomy_family",
                "taxonomy_mapping_status",
                "taxonomy_taxid",
                "taxonomy_display_name",
                "integrated_level",
                "picked",
                "formal_pick_exclusion_rule",
            )
            if old.get(field) != new.get(field)
        ]
        if not fields_changed and old and new:
            continue
        old_picked = bool(old.get("picked"))
        new_picked = bool(new.get("picked"))
        if not old:
            change = "candidate_added"
        elif not new:
            change = "candidate_removed"
        elif old_picked != new_picked:
            change = "picked_gained" if new_picked else "picked_lost"
        elif old.get("integrated_level") != new.get("integrated_level"):
            change = "level_changed"
        else:
            change = "taxonomy_metadata_changed"
        changes.append(
            {
                "patient_id": patient,
                "organism_name": name,
                "change": change,
                "changed_fields": "|".join(fields_changed),
                "answer_status": answer_status(patient, name, answers),
                "before_picked": old_picked,
                "after_picked": new_picked,
                "before_level": old.get("integrated_level", ""),
                "after_level": new.get("integrated_level", ""),
                "before_family": old.get("taxonomy_family", ""),
                "after_family": new.get("taxonomy_family", ""),
                "before_mapping_status": old.get("taxonomy_mapping_status", ""),
                "after_mapping_status": new.get("taxonomy_mapping_status", ""),
                "before_taxid": old.get("taxonomy_taxid", ""),
                "after_taxid": new.get("taxonomy_taxid", ""),
            }
        )
    return changes


def run(a_root: Path, patient_root: Path, answers_path: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    contract_path = a_root / "baseline_contract.json"
    inventory_path = a_root / "artifact_inventory.csv"
    contract = read_json(contract_path)
    patient_ids = [int(value) for value in contract["patient_selection"]]
    inventory = read_inventory(inventory_path)
    answers = load_answer_rows(answers_path)
    output_root.mkdir(parents=True)

    control_rows, control_summary = score_mode(
        mode="base-only",
        patient_ids=patient_ids,
        inventory=inventory,
        patient_root=patient_root,
        output_dir=output_root / "control_base_taxonomy",
    )
    treatment_rows, treatment_summary = score_mode(
        mode="current-overlays",
        patient_ids=patient_ids,
        inventory=inventory,
        patient_root=patient_root,
        output_dir=output_root / "treatment_new_taxonomy",
    )
    if control_summary["source_hash_mismatches"]:
        raise ValueError("Frozen A source hash validation failed")
    if control_summary["control_decision_mismatch_patients"]:
        raise ValueError(
            "Base-taxonomy control did not reproduce A decisions: "
            f"{control_summary['control_decision_mismatch_patients']}"
        )

    control_picked = grouped_picked(control_rows)
    treatment_picked = grouped_picked(treatment_rows)
    changes = compare_rows(control_rows, treatment_rows, answers)
    evaluation_dir = output_root / "evaluation"
    write_csv(evaluation_dir / "per_patient_candidate_changes.csv", changes)
    metrics = {
        "control_base_taxonomy_picked": evaluate_set(control_picked, answers),
        "treatment_new_taxonomy_picked": evaluate_set(treatment_picked, answers),
    }
    evaluation = {
        "schema_version": "kh_ablation_e_evaluation.v1",
        "scope": "Post-hoc evaluation after answer-blind scorer outputs were written.",
        "metrics": metrics,
        "change_counts": dict(Counter(row["change"] for row in changes)),
        "changed_candidate_count": len(changes),
    }
    write_json(evaluation_dir / "metrics.json", evaluation)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "experiment_id": "KH_ABLATION_E_LEGACY_V20_TAXONOMY_ONLY_STABLE_V1",
        "status": "frozen_taxonomy_only_shadow",
        "answer_blind_generation": True,
        "benchmark_evaluation_is_posthoc": True,
        "factors": {
            "mngs": contract["factors"]["mngs"],
            "scorer": contract["factors"]["scorer"],
            "taxonomy": "current_versioned_gap_auto_reviewed_overlays",
            "history": contract["factors"]["history"],
        },
        "common_contract": {
            "patient_count": len(patient_ids),
            "patient_selection": [str(value) for value in patient_ids],
            "A_contract_path": str(contract_path.resolve()),
            "A_contract_sha256": sha256_file(contract_path),
            "answer_path": str(answers_path.resolve()),
            "answer_sha256": sha256_file(answers_path),
        },
        "control_validation": {
            "source_hash_mismatch_count": len(control_summary["source_hash_mismatches"]),
            "decision_mismatch_patient_count": len(control_summary["control_decision_mismatch_patients"]),
            "reproduces_preserved_A_deterministic_decisions": True,
        },
        "implementation": {
            "runner_path": str(Path(__file__).resolve()),
            "runner_sha256": sha256_file(Path(__file__)),
            "legacy_scorer_path": str(Path(legacy.__file__).resolve()),
            "legacy_scorer_sha256": sha256_file(Path(legacy.__file__)),
        },
        "metrics": metrics,
        "change_counts": evaluation["change_counts"],
        "taxonomy_snapshots": {
            "control": control_summary["taxonomy_snapshot"],
            "treatment": treatment_summary["taxonomy_snapshot"],
        },
        "outputs": {
            "control_candidate_sha256": sha256_file(output_root / "control_base_taxonomy/candidate_decisions.csv"),
            "treatment_candidate_sha256": sha256_file(output_root / "treatment_new_taxonomy/candidate_decisions.csv"),
            "change_table_sha256": sha256_file(evaluation_dir / "per_patient_candidate_changes.csv"),
            "metrics_sha256": sha256_file(evaluation_dir / "metrics.json"),
        },
        "limitations": [
            "This isolates taxonomy inside the frozen legacy-v20-compatible scorer; it does not use the developing new scorer.",
            "LLM review, new history, promotion v4, Possible, and OBER are excluded from E.",
        ],
    }
    write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-root", type=Path, default=DEFAULT_A_ROOT)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run(args.a_root, args.patient_root, args.answers, args.output), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

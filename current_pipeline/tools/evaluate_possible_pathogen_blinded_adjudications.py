"""Unblind frozen OBER decisions and evaluate possible-pathogen visibility."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools.build_index_event_clinical_timeline_shadow import sha256_file
from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.recalculate_kh_benchmark_metrics import match_type, unique_names


DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics/kh_answers_clinical_revision_20260918.csv"
)
SHOW_VISIBILITY = {"show_primary", "show_secondary", "show_context"}


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def grouped(rows: list[dict[str, Any]]) -> dict[int, list[str]]:
    output: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        output[int(row["patient_id"])].append(str(row["organism_name"]))
    return {patient: unique_names(names) for patient, names in output.items()}


def run(
    frozen_dir: Path, private_map_path: Path, shadow_root: Path,
    answer_path: Path, output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    freeze_manifest = read_json(frozen_dir / "freeze_manifest.json")
    if not freeze_manifest.get("ready_for_posthoc_unblinding"):
        raise ValueError("Adjudications have not passed the blind freeze gate")
    frozen_path = frozen_dir / "adjudications.frozen.jsonl"
    if sha256_file(frozen_path) != freeze_manifest.get("frozen_adjudications_sha256"):
        raise ValueError("Frozen adjudication hash mismatch")

    private_map = read_json(private_map_path)
    identity_by_id = {
        str(row["blind_case_id"]): row for row in private_map.get("rows") or []
    }
    records = read_jsonl(frozen_path)
    if set(identity_by_id) != {str(record.get("case_id")) for record in records}:
        raise ValueError("Frozen adjudications and private case map do not have identical case IDs")

    source_rows = read_csv(shadow_root / "possible_pathogen_decisions.csv")
    strict_rows = [row for row in source_rows if row["clinical_decision"] == "picked_shadow"]
    initial_possible_rows = [
        row for row in source_rows
        if row["possible_reporting_role"] not in {
            "strict_picked_report", "not_selected_for_complete_report"
        }
    ]
    adjudicated_rows = []
    visible_possible_rows = []
    for record in records:
        blind_id = str(record["case_id"])
        identity = identity_by_id[blind_id]
        adjudication = record["adjudication"]
        decision = adjudication["decision"]
        row = {
            "blind_case_id": blind_id,
            "patient_id": int(identity["patient_id"]),
            "organism_name": identity["organism_name"],
            "possible_reporting_role": identity["possible_reporting_role"],
            "recommended_action": decision["recommended_action"],
            "clinician_visibility": decision["clinician_visibility"],
            "confidence": decision["confidence"],
            "one_sentence_reason": decision["one_sentence_reason"],
            "human_review_required": adjudication["audit"]["human_review_required"],
        }
        adjudicated_rows.append(row)
        if decision["clinician_visibility"] in SHOW_VISIBILITY:
            visible_possible_rows.append(row)

    strict_group = grouped(strict_rows)
    initial_group = grouped(strict_rows + initial_possible_rows)
    ober_group = grouped(strict_rows + visible_possible_rows)
    answers = load_answer_rows(answer_path)
    metrics = {
        "strict_picked": evaluate_set(strict_group, answers),
        "initial_complete_report": evaluate_set(initial_group, answers),
        "strict_plus_ober_visible_possible": evaluate_set(ober_group, answers),
    }

    for row in adjudicated_rows:
        patient_answers = answers.get(row["patient_id"])
        if patient_answers is None:
            row["posthoc_benchmark_outcome"] = "unlabeled_patient"
            row["matched_answer"] = ""
            continue
        matched = next((answer for answer in patient_answers if match_type(
            row["organism_name"], answer, genus_relaxed=False
        )), None)
        row["posthoc_benchmark_outcome"] = (
            "matched_benchmark_answer" if matched else "unmatched_in_labeled_patient"
        )
        row["matched_answer"] = matched or ""

    output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(adjudicated_rows[0]) if adjudicated_rows else []
    with (output_dir / "unblinded_adjudication_audit.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(adjudicated_rows)
    result = {
        "schema_version": "possible_pathogen_blinded_adjudication_evaluation.v1",
        "posthoc_only": True,
        "freeze_manifest": str((frozen_dir / "freeze_manifest.json").resolve()),
        "frozen_adjudications_sha256": freeze_manifest["frozen_adjudications_sha256"],
        "adjudication_count": len(adjudicated_rows),
        "ober_visible_possible_count": len(visible_possible_rows),
        "decision_counts": dict(Counter(row["recommended_action"] for row in adjudicated_rows)),
        "visibility_counts": dict(Counter(row["clinician_visibility"] for row in adjudicated_rows)),
        "metrics": metrics,
        "possible_posthoc_outcome_by_visibility": {
            visibility: dict(Counter(
                row["posthoc_benchmark_outcome"] for row in adjudicated_rows
                if row["clinician_visibility"] == visibility
            ))
            for visibility in sorted({row["clinician_visibility"] for row in adjudicated_rows})
        },
        "held_out_validation": {
            "performed": False,
            "reason": "No second cohort with the same multi-assay per-test source contract and compatible adjudicated labels was supplied to this run.",
        },
    }
    (output_dir / "metrics.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("frozen_dir", type=Path)
    parser.add_argument("private_map", type=Path)
    parser.add_argument("shadow_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    args = parser.parse_args()
    print(json.dumps(run(
        args.frozen_dir, args.private_map, args.shadow_root,
        args.answers, args.output_dir,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

"""Build the common A-G contract report from frozen A, B, D, and E runs.

C, G, and F are deliberately left pending until a new scorer version is
explicitly frozen.  OBER is outside this primary ablation report.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools.evaluate_multi_assay_candidate_entry import load_answer_rows
from tools.recalculate_kh_benchmark_metrics import match_type


DEFAULT_CONTRACT = Path("rules/kh_ablation_common_contract_v1.json")
DEFAULT_A = Path("outputs/runs/2026-09-24_KH_ablation_A_frozen_old_baseline_v1")
DEFAULT_B = Path("outputs/runs/2026-09-28_KH_ablation_B_legacy_v20_multi_assay_stable_v1")
DEFAULT_B_EVAL = Path("outputs/runs/2026-09-28_KH_ablation_B_legacy_v20_multi_assay_stable_v1_evaluation")
DEFAULT_D = Path("outputs/runs/2026-09-28_KH_ablation_D_legacy_v20_history_v4c_stable_v1")
DEFAULT_E = Path("outputs/runs/2026-09-28_KH_ablation_E_legacy_v20_taxonomy_only_stable_v1")
DEFAULT_OUTPUT = Path("outputs/runs/2026-09-28_KH_ablation_A_B_D_E_common_v1")
SCHEMA_VERSION = "kh_ablation_common_report.v1"


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected object: {path}")
    return payload


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metric_row(
    group: str,
    endpoint: str,
    status: str,
    factors: dict[str, str],
    metric: dict[str, Any] | None,
    note: str = "",
) -> dict[str, Any]:
    metric = metric or {}
    return {
        "ablation_group": group,
        "endpoint": endpoint,
        "status": status,
        "mngs": factors["mngs"],
        "scorer": factors["scorer"],
        "taxonomy": factors["taxonomy"],
        "history": factors["history"],
        "matched": metric.get("matched", ""),
        "predicted_in_labeled_patients": metric.get(
            "predicted_in_labeled_patients", metric.get("outputs_in_labeled_patients", "")
        ),
        "answer_count": metric.get("answer_count", ""),
        "precision": metric.get("precision", ""),
        "recall": metric.get("recall", ""),
        "f1": metric.get("f1", ""),
        "note": note,
    }


def status_for(patient: int, organism: str, answers: dict[int, list[str]]) -> str:
    gold = answers.get(patient, [])
    if not gold:
        return "unlabeled_patient"
    return (
        "matched_answer"
        if any(match_type(organism, answer, genus_relaxed=False) for answer in gold)
        else "unmatched_output"
    )


def normalized_changes(
    b_eval: Path,
    d_root: Path,
    e_root: Path,
    answers: dict[int, list[str]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in read_csv(b_eval / "picked_changes_vs_selected_dna.csv"):
        patient = int(item["patient_id"])
        name = item["organism_name"]
        added = item["change"] == "added_by_multi_assay"
        rows.append(
            {
                "ablation_group": "B",
                "patient_id": patient,
                "organism_name": name,
                "endpoint": "strict_formal_picked",
                "before_state": "Not_Picked" if added else "Picked",
                "after_state": "Picked" if added else "Not_Picked",
                "change_type": "picked_gained" if added else "picked_lost",
                "answer_status": status_for(patient, name, answers),
                "rule_ids": "",
                "evidence_ids": "",
                "detail": item["change"],
            }
        )
    for item in read_csv(d_root / "evaluation/per_patient_candidate_changes.csv"):
        rows.append(
            {
                "ablation_group": "D",
                "patient_id": int(item["patient_id"]),
                "organism_name": item["organism_name"],
                "endpoint": "complete_visible_picked_high_context",
                "before_state": item["before_decision"],
                "after_state": item["after_decision"],
                "change_type": item["visibility_change"],
                "answer_status": item["answer_status"],
                "rule_ids": item["history_rule_ids"],
                "evidence_ids": item["consumed_evidence_ids"],
                "detail": f"{item['before_route']}->{item['after_route']}",
            }
        )
    for item in read_csv(e_root / "evaluation/per_patient_candidate_changes.csv"):
        rows.append(
            {
                "ablation_group": "E",
                "patient_id": int(item["patient_id"]),
                "organism_name": item["organism_name"],
                "endpoint": "taxonomy_metadata",
                "before_state": item["before_family"],
                "after_state": item["after_family"],
                "change_type": item["change"],
                "answer_status": item["answer_status"],
                "rule_ids": "",
                "evidence_ids": "",
                "detail": item["changed_fields"],
            }
        )
    return sorted(rows, key=lambda row: (row["ablation_group"], row["patient_id"], row["organism_name"]))


def patient_change_summary(changes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in changes:
        grouped[(row["ablation_group"], int(row["patient_id"]))].append(row)
    output: list[dict[str, Any]] = []
    for (group, patient), items in sorted(grouped.items()):
        def names(change: str, answer_status: str | None = None) -> str:
            selected = [
                row["organism_name"]
                for row in items
                if row["change_type"] == change
                and (answer_status is None or row["answer_status"] == answer_status)
            ]
            return "|".join(selected)

        output.append(
            {
                "ablation_group": group,
                "patient_id": patient,
                "tp_gained": names("picked_gained", "matched_answer")
                or names("report_added", "matched_answer"),
                "tp_lost": names("picked_lost", "matched_answer")
                or names("report_removed", "matched_answer"),
                "fp_added": names("picked_gained", "unmatched_output")
                or names("report_added", "unmatched_output"),
                "fp_removed": names("picked_lost", "unmatched_output")
                or names("report_removed", "unmatched_output"),
                "tier_or_route_changes": "|".join(
                    row["organism_name"]
                    for row in items
                    if row["change_type"] in {"report_tier_changed", "route_only_changed"}
                ),
                "taxonomy_metadata_changes": names("taxonomy_metadata_changed"),
            }
        )
    return output


def run(
    contract_path: Path,
    a_root: Path,
    b_root: Path,
    b_eval: Path,
    d_root: Path,
    e_root: Path,
    output_root: Path,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    contract = read_json(contract_path)
    a = read_json(a_root / "baseline_contract.json")
    b = read_json(b_root / "run_manifest.json")
    d = read_json(d_root / "run_manifest.json")
    e = read_json(e_root / "run_manifest.json")
    b_metrics = read_json(b_eval / "metrics.json")["metrics"]
    factors = contract["factor_matrix"]
    answer_path = Path(a["answer_contract"]["path"])
    answers = load_answer_rows(answer_path)

    common_checks = {
        "A_patient_count_is_33": a["patient_count"] == contract["common_controls"]["expected_patient_count"],
        "A_labeled_patient_count_is_30": a["answer_contract"]["labeled_patient_count"]
        == contract["common_controls"]["expected_labeled_patient_count"],
        "A_answer_count_is_60": a["answer_contract"]["answer_organism_count"]
        == contract["common_controls"]["expected_answer_organism_count"],
        "B_references_same_A_contract": b["inputs"]["frozen_A_contract"]["sha256"]
        == sha256_file(a_root / "baseline_contract.json"),
        "D_references_same_A_contract": d["common_contract"]["A_contract_sha256"]
        == sha256_file(a_root / "baseline_contract.json"),
        "E_references_same_A_contract": e["common_contract"]["A_contract_sha256"]
        == sha256_file(a_root / "baseline_contract.json"),
        "B_answer_hash_matches_A": b["inputs"]["answer_file"]["sha256"] == a["answer_contract"]["sha256"],
        "D_answer_hash_matches_A": d["common_contract"]["answer_sha256"] == a["answer_contract"]["sha256"],
        "E_answer_hash_matches_A": e["common_contract"]["answer_sha256"] == a["answer_contract"]["sha256"],
        "OBER_excluded": contract["common_controls"]["ober_in_primary_ablation"] is False,
    }
    if not all(common_checks.values()):
        raise ValueError(f"Common contract validation failed: {common_checks}")

    summary_rows = [
        metric_row(
            "A", "strict_formal_picked", "complete_frozen", factors["A"],
            a["frozen_metrics"]["picked"], "Preserved historical output comparator",
        ),
        metric_row(
            "A", "complete_visible_picked_high_context", "complete_frozen", factors["A"],
            a["frozen_metrics"]["combined_picked_high_context"], "Picked plus Luna High/Context",
        ),
        metric_row(
            "B", "strict_formal_picked", "complete_frozen", factors["B"],
            b_metrics["new_multi_assay_through_legacy_picked"], "New multi-assay mNGS only",
        ),
        metric_row(
            "D", "strict_formal_picked", "complete_frozen", factors["D"],
            d["metrics"]["D_history_strict_picked"], "History never creates Picked by itself",
        ),
        metric_row(
            "D", "complete_visible_picked_high_context", "complete_frozen", factors["D"],
            d["metrics"]["D_history_complete_visible"], "A-compatible v4C history routing",
        ),
        metric_row(
            "E", "strict_formal_picked", "complete_frozen", factors["E"],
            e["metrics"]["treatment_new_taxonomy_picked"], "New taxonomy only",
        ),
    ]
    for group in ("C", "G", "F"):
        summary_rows.append(
            metric_row(
                group,
                "strict_formal_picked",
                "pending_new_scorer_freeze",
                factors[group],
                None,
                "Do not run formally until the new scorer version and hash are frozen",
            )
        )

    changes = normalized_changes(b_eval, d_root, e_root, answers)
    required = set(contract["normalized_candidate_schema"]["required_fields"])
    schema_valid = all(required <= set(row) for row in changes)
    if not schema_valid:
        raise ValueError("Normalized per-patient change schema validation failed")
    patient_summary = patient_change_summary(changes)

    output_root.mkdir(parents=True)
    summary_fields = list(summary_rows[0])
    change_fields = [
        "ablation_group", "patient_id", "organism_name", "endpoint", "before_state",
        "after_state", "change_type", "answer_status", "rule_ids", "evidence_ids", "detail",
    ]
    patient_fields = [
        "ablation_group", "patient_id", "tp_gained", "tp_lost", "fp_added",
        "fp_removed", "tier_or_route_changes", "taxonomy_metadata_changes",
    ]
    write_csv(output_root / "ablation_summary.csv", summary_rows, summary_fields)
    write_csv(output_root / "per_patient_changes.csv", changes, change_fields)
    write_csv(output_root / "per_patient_change_summary.csv", patient_summary, patient_fields)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "report_id": "KH_ABLATION_A_B_D_E_COMMON_V1",
        "status": "A_B_D_E_complete_C_G_F_pending",
        "common_contract": {
            "path": str(contract_path.resolve()),
            "sha256": sha256_file(contract_path),
        },
        "common_checks": common_checks,
        "patient_selection": a["patient_selection"],
        "patient_count": a["patient_count"],
        "answer_contract": a["answer_contract"],
        "evaluator_contract": {
            "matching": contract["common_controls"]["matching_policy"],
            "evaluate_module": str(Path("tools/evaluate_multi_assay_candidate_entry.py").resolve()),
            "evaluate_module_sha256": sha256_file(Path("tools/evaluate_multi_assay_candidate_entry.py")),
            "matching_module": str(Path("tools/recalculate_kh_benchmark_metrics.py").resolve()),
            "matching_module_sha256": sha256_file(Path("tools/recalculate_kh_benchmark_metrics.py")),
        },
        "report_builder": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256_file(Path(__file__)),
        },
        "group_status": {
            "A": "complete_frozen",
            "B": "complete_frozen",
            "C": "pending_new_scorer_freeze",
            "D": "complete_frozen",
            "E": "complete_frozen",
            "G": "pending_new_scorer_freeze",
            "F": "pending_new_scorer_freeze",
        },
        "normalized_schema_valid": schema_valid,
        "change_counts": dict(Counter(row["ablation_group"] for row in changes)),
        "outputs": {
            "ablation_summary_sha256": sha256_file(output_root / "ablation_summary.csv"),
            "per_patient_changes_sha256": sha256_file(output_root / "per_patient_changes.csv"),
            "per_patient_change_summary_sha256": sha256_file(output_root / "per_patient_change_summary.csv"),
        },
        "ober": "excluded_from_primary_ablation",
    }
    write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--a-root", type=Path, default=DEFAULT_A)
    parser.add_argument("--b-root", type=Path, default=DEFAULT_B)
    parser.add_argument("--b-eval", type=Path, default=DEFAULT_B_EVAL)
    parser.add_argument("--d-root", type=Path, default=DEFAULT_D)
    parser.add_argument("--e-root", type=Path, default=DEFAULT_E)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.contract, args.a_root, args.b_root, args.b_eval, args.d_root, args.e_root, args.output),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

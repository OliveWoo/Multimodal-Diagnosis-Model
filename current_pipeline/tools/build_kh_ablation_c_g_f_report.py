"""Freeze C/G/F manifests and build the complete non-OBER A-G report.

The A-G primary endpoint stops after the clinical deterministic scorer.
Promotion v4 and Possible-pathogen v1 are retained as supplementary downstream
endpoints so their effect is visible without being confused with the scorer
ablation itself.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from tools.evaluate_multi_assay_candidate_entry import load_answer_rows
from tools.pathogen_normalization import canonical_key
from tools.recalculate_kh_benchmark_metrics import match_type


DEFAULT_CONTRACT = Path("rules/kh_ablation_common_contract_v1.json")
DEFAULT_RELEASE = Path("rules/kh_test_aware_release_v1.json")
DEFAULT_OLD_COMMON = Path("outputs/runs/2026-09-28_KH_ablation_A_B_D_E_common_v1")
DEFAULT_B = Path("outputs/runs/2026-09-28_KH_ablation_B_legacy_v20_multi_assay_stable_v1")
DEFAULT_C = Path("outputs/runs/2026-09-28_KH_ablation_C_new_scorer_v1")
DEFAULT_G = Path("outputs/runs/2026-09-28_KH_ablation_G_new_scorer_taxonomy_v1")
DEFAULT_F = Path("outputs/runs/2026-09-28_KH_ablation_F_full_non_ober_v1")
DEFAULT_OUTPUT = Path("outputs/runs/2026-09-28_KH_ablation_A_G_full_non_ober_v1")
SCHEMA_VERSION = "kh_ablation_full_non_ober_report.v1"


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected JSON object: {path}")
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
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_contract(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def validate_release(release_path: Path) -> dict[str, Any]:
    release = read_json(release_path)
    if release.get("status") != "frozen_shadow_for_ablation":
        raise ValueError(f"Release is not frozen for ablation: {release_path}")
    mismatches: list[str] = []
    for component in release.get("component_files") or []:
        path = Path(component["path"])
        if not path.is_file():
            mismatches.append(f"missing:{path}")
        elif sha256_file(path) != component["sha256"]:
            mismatches.append(f"hash:{path}")
    if mismatches:
        raise ValueError(f"Frozen release validation failed: {mismatches}")
    return release


def clinical_metrics(root: Path) -> dict[str, Any]:
    return read_json(root / "clinical_evaluation/metrics.json")["metrics"]


def promotion_metrics(root: Path) -> dict[str, Any]:
    return read_json(root / "promotion_v4_evaluation/metrics.json")["metrics"]


def possible_metrics(root: Path) -> dict[str, Any]:
    return read_json(root / "possible_pathogen_v1_evaluation/metrics.json")["metrics"]


def metric_row(
    group: str,
    endpoint: str,
    factors: dict[str, str],
    metric: dict[str, Any],
    note: str,
) -> dict[str, Any]:
    return {
        "ablation_group": group,
        "endpoint": endpoint,
        "status": "complete_frozen",
        "mngs": factors["mngs"],
        "scorer": factors["scorer"],
        "taxonomy": factors["taxonomy"],
        "history": factors["history"],
        "matched": metric.get("matched", ""),
        "predicted_in_labeled_patients": metric.get("predicted_in_labeled_patients", ""),
        "answer_count": metric.get("answer_count", ""),
        "precision": metric.get("precision", ""),
        "recall": metric.get("recall", ""),
        "f1": metric.get("f1", ""),
        "note": note,
    }


def answer_status(patient: int, organism: str, answers: dict[int, list[str]]) -> str:
    gold = answers.get(patient, [])
    if not gold:
        return "unlabeled_patient"
    return (
        "matched_answer"
        if any(match_type(organism, answer, genus_relaxed=False) for answer in gold)
        else "unmatched_output"
    )


def tier(decision: str) -> str:
    return {
        "picked_shadow": "Picked",
        "review_high_priority": "High",
        "review_context_needed": "Context",
        "review_low_specificity": "Hold",
    }.get(decision, "Absent")


def _best_rows(rows: Iterable[dict[str, str]]) -> dict[tuple[int, str], dict[str, str]]:
    order = {"Picked": 0, "High": 1, "Context": 2, "Audit": 3, "Hold": 4, "Absent": 5}
    output: dict[tuple[int, str], dict[str, str]] = {}
    for row in rows:
        patient = int(row["patient_id"])
        name = row["organism_name"]
        key = (patient, canonical_key(name))
        row = dict(row)
        row["_tier"] = tier(row.get("clinical_decision", ""))
        if key not in output or order[row["_tier"]] < order[output[key]["_tier"]]:
            output[key] = row
    return output


def _old_b_picks(path: Path) -> dict[tuple[int, str], dict[str, str]]:
    output: dict[tuple[int, str], dict[str, str]] = {}
    for row in read_csv(path):
        if str(row.get("picked", "")).casefold() != "true":
            continue
        patient = int(row["patient_id"])
        name = row["organism_name"]
        output[(patient, canonical_key(name))] = {**row, "_tier": "Picked"}
    return output


def compare_candidate_maps(
    group: str,
    before: dict[tuple[int, str], dict[str, str]],
    after: dict[tuple[int, str], dict[str, str]],
    answers: dict[int, list[str]],
    *,
    strict_only: bool = False,
) -> list[dict[str, Any]]:
    changes: list[dict[str, Any]] = []
    for key in sorted(set(before) | set(after)):
        old = before.get(key)
        new = after.get(key)
        old_tier = old.get("_tier", "Absent") if old else "Absent"
        new_tier = new.get("_tier", "Absent") if new else "Absent"
        old_picked = old_tier == "Picked"
        new_picked = new_tier == "Picked"
        if strict_only and old_picked == new_picked:
            continue
        if not strict_only and old_tier == new_tier:
            continue
        if old_picked != new_picked:
            change_type = "picked_gained" if new_picked else "picked_lost"
            endpoint = "strict_formal_picked"
        elif old is None:
            change_type = "candidate_added"
            endpoint = "picked_high_context"
        elif new is None:
            change_type = "candidate_removed"
            endpoint = "picked_high_context"
        else:
            change_type = "report_tier_changed"
            endpoint = "picked_high_context"
        source = new or old or {}
        patient = key[0]
        name = source.get("organism_name", key[1])
        changes.append(
            {
                "ablation_group": group,
                "patient_id": patient,
                "organism_name": name,
                "endpoint": endpoint,
                "before_state": old_tier,
                "after_state": new_tier,
                "change_type": change_type,
                "answer_status": answer_status(patient, name, answers),
                "rule_ids": source.get("rule_ids", ""),
                "evidence_ids": source.get("consumed_evidence_ids", ""),
                "detail": f"{old_tier}->{new_tier}",
            }
        )
    return changes


def patient_change_summary(changes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in changes:
        grouped[(row["ablation_group"], int(row["patient_id"]))].append(row)
    output: list[dict[str, Any]] = []
    for (group, patient), items in sorted(grouped.items()):
        def names(kind: str, status: str | None = None) -> str:
            return "|".join(
                row["organism_name"]
                for row in items
                if row["change_type"] == kind
                and (status is None or row["answer_status"] == status)
            )

        output.append(
            {
                "ablation_group": group,
                "patient_id": patient,
                "tp_gained": names("picked_gained", "matched_answer"),
                "tp_lost": names("picked_lost", "matched_answer"),
                "fp_added": names("picked_gained", "unmatched_output"),
                "fp_removed": names("picked_lost", "unmatched_output"),
                "tier_or_route_changes": "|".join(
                    row["organism_name"]
                    for row in items
                    if row["change_type"] in {
                        "report_tier_changed", "candidate_added", "candidate_removed"
                    }
                ),
                "taxonomy_metadata_changes": names("taxonomy_metadata_changed"),
            }
        )
    return output


def build_group_manifest(
    group: str,
    root: Path,
    factors: dict[str, str],
    release_path: Path,
    scorer_root: Path,
    source_compact_summary: Path,
) -> dict[str, Any]:
    scorer_summary_path = scorer_root / "summary.json"
    clinical_summary_path = root / "clinical/summary.json"
    clinical_metrics_path = root / "clinical_evaluation/metrics.json"
    promotion_summary_path = root / "promotion_v4/summary.json"
    promotion_metrics_path = root / "promotion_v4_evaluation/metrics.json"
    possible_summary_path = root / "possible_pathogen_v1/summary.json"
    possible_metrics_path = root / "possible_pathogen_v1_evaluation/metrics.json"
    scorer = read_json(scorer_summary_path)
    clinical = read_json(clinical_summary_path)
    metrics = read_json(clinical_metrics_path)
    promotion = read_json(promotion_summary_path)
    possible = read_json(possible_summary_path)

    checks = {
        "patient_count_is_33": scorer.get("patient_count") == 33
        and clinical.get("patient_count") == 33,
        "answer_blind_scorer": scorer.get("answer_blind") is True,
        "answer_blind_clinical": clinical.get("answer_blind") is True,
        "answer_count_is_60": metrics["metrics"]["clinical_picked_shadow"].get("answer_count") == 60,
        "compact_summary_hash_matches": scorer.get("source_compact_entry_summary_sha256")
        == sha256_file(source_compact_summary),
        "analytical_policy_is_v2": scorer.get("policy_id") == "test_aware_deterministic_shadow_v2",
        "clinical_policy_is_v1": clinical.get("policy_id")
        == "test_aware_clinical_deterministic_shadow_v1",
        "promotion_policy_is_v4": promotion.get("policy_id")
        == "test_aware_clinical_promotion_v4_independent_clinical_axis",
        "possible_policy_is_v1": possible.get("policy_id")
        == "test_aware_possible_pathogen_v1_answer_blind_generic",
    }
    if group in {"C", "G"}:
        checks["history_is_absent"] = clinical.get("history_route_root") is None
    else:
        history = read_json(root / "history_route/summary.json")
        checks["history_is_v4c"] = history.get("policy_id") == "test_aware_history_route_v4c_combined"
        checks["clinical_uses_history"] = clinical.get("history_route_root") is not None
    if not all(checks.values()):
        raise ValueError(f"Group {group} validation failed: {checks}")

    manifest = {
        "schema_version": "kh_ablation_new_scorer_group.v1",
        "group": group,
        "status": "complete_frozen",
        "factors": factors,
        "answer_blind_generation": True,
        "benchmark_evaluation_is_posthoc": True,
        "primary_endpoint_boundary": "clinical_deterministic_scorer_before_promotion",
        "release": file_contract(release_path),
        "frozen_inputs": {
            "compact_summary": file_contract(source_compact_summary),
            "scorer_summary": file_contract(scorer_summary_path),
            "clinical_summary": file_contract(clinical_summary_path),
        },
        "checks": checks,
        "primary_metrics": metrics["metrics"],
        "supplementary_metrics": {
            "promotion_v4": read_json(promotion_metrics_path)["metrics"],
            "possible_pathogen_v1": read_json(possible_metrics_path)["metrics"],
        },
        "outputs": {
            "clinical_metrics": file_contract(clinical_metrics_path),
            "promotion_summary": file_contract(promotion_summary_path),
            "promotion_metrics": file_contract(promotion_metrics_path),
            "possible_summary": file_contract(possible_summary_path),
            "possible_metrics": file_contract(possible_metrics_path),
        },
        "ober": "excluded",
    }
    if group == "F":
        manifest["frozen_inputs"]["history_summary"] = file_contract(
            root / "history_route/summary.json"
        )
        manifest["frozen_inputs"]["resolved_history_policy"] = file_contract(
            root / "history_route/resolved_policy_snapshot.json"
        )
    write_json(root / "run_manifest.json", manifest)
    return manifest


def run(
    contract_path: Path,
    release_path: Path,
    old_common_root: Path,
    b_root: Path,
    c_root: Path,
    g_root: Path,
    f_root: Path,
    output_root: Path,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    contract = read_json(contract_path)
    if contract["common_controls"].get("ober_in_primary_ablation") is not False:
        raise ValueError("Common contract does not exclude OBER")
    validate_release(release_path)
    factors = contract["factor_matrix"]

    base_compact_summary = Path(
        "outputs/runs/2026-09-23_KH_compact_multi_assay_entry_source_contract_v5/summary.json"
    )
    taxonomy_compact_summary = Path(
        "outputs/runs/2026-09-25_KH_compact_multi_assay_entry_taxonomy_auto_v1/summary.json"
    )
    manifests = {
        "C": build_group_manifest(
            "C", c_root, factors["C"], release_path, c_root / "scorer", base_compact_summary
        ),
        "G": build_group_manifest(
            "G", g_root, factors["G"], release_path, g_root / "scorer", taxonomy_compact_summary
        ),
        "F": build_group_manifest(
            "F", f_root, factors["F"], release_path, g_root / "scorer", taxonomy_compact_summary
        ),
    }

    old_rows = [
        row for row in read_csv(old_common_root / "ablation_summary.csv")
        if row["ablation_group"] not in {"C", "G", "F"}
    ]
    summary_rows: list[dict[str, Any]] = list(old_rows)
    for group, root in (("C", c_root), ("G", g_root), ("F", f_root)):
        cm = clinical_metrics(root)
        pm = promotion_metrics(root)
        xm = possible_metrics(root)
        summary_rows.extend(
            [
                metric_row(
                    group,
                    "strict_formal_picked",
                    factors[group],
                    cm["clinical_picked_shadow"],
                    "Primary A-G endpoint; before promotion v4",
                ),
                metric_row(
                    group,
                    "picked_plus_high",
                    factors[group],
                    cm["clinical_picked_plus_high_review"],
                    "Secondary review endpoint; before promotion v4",
                ),
                metric_row(
                    group,
                    "strict_after_promotion_v4",
                    factors[group],
                    pm["clinical_picked_shadow"],
                    "Supplementary; promotion is outside the factor matrix",
                ),
                metric_row(
                    group,
                    "complete_report_after_possible_v1",
                    factors[group],
                    xm["complete_report"],
                    "Supplementary reporting role; does not change formal Picked",
                ),
            ]
        )

    old_manifest = read_json(old_common_root / "run_manifest.json")
    answer_path = Path(old_manifest["answer_contract"]["path"])
    answers = load_answer_rows(answer_path)
    changes = read_csv(old_common_root / "per_patient_changes.csv")
    b_picks = _old_b_picks(b_root / "legacy_multi_assay_candidate_decisions.csv")
    c_rows = _best_rows(read_csv(c_root / "clinical/clinical_decisions.csv"))
    g_rows = _best_rows(read_csv(g_root / "clinical/clinical_decisions.csv"))
    f_rows = _best_rows(read_csv(f_root / "clinical/clinical_decisions.csv"))
    changes.extend(compare_candidate_maps("C", b_picks, c_rows, answers, strict_only=True))
    changes.extend(compare_candidate_maps("G", c_rows, g_rows, answers))
    changes.extend(compare_candidate_maps("F", g_rows, f_rows, answers))
    changes = sorted(
        changes,
        key=lambda row: (
            "ABCD EGF".replace(" ", "").find(row["ablation_group"]),
            int(row["patient_id"]),
            row["organism_name"],
        ),
    )
    patient_summary = patient_change_summary(changes)

    output_root.mkdir(parents=True)
    summary_fields = [
        "ablation_group", "endpoint", "status", "mngs", "scorer", "taxonomy",
        "history", "matched", "predicted_in_labeled_patients", "answer_count",
        "precision", "recall", "f1", "note",
    ]
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

    strict = {
        group: clinical_metrics(root)["clinical_picked_shadow"]
        for group, root in (("C", c_root), ("G", g_root), ("F", f_root))
    }
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "report_id": "KH_ABLATION_A_G_FULL_NON_OBER_V1",
        "status": "A_B_C_D_E_G_F_complete_frozen",
        "common_contract": file_contract(contract_path),
        "test_aware_release": file_contract(release_path),
        "patient_count": 33,
        "answer_contract": file_contract(answer_path),
        "group_status": {
            group: "complete_frozen" for group in ("A", "B", "C", "D", "E", "G", "F")
        },
        "new_group_manifests": {
            group: file_contract(root / "run_manifest.json")
            for group, root in (("C", c_root), ("G", g_root), ("F", f_root))
        },
        "primary_strict_metrics": strict,
        "change_counts": dict(Counter(row["ablation_group"] for row in changes)),
        "outputs": {
            "ablation_summary": file_contract(output_root / "ablation_summary.csv"),
            "per_patient_changes": file_contract(output_root / "per_patient_changes.csv"),
            "per_patient_change_summary": file_contract(
                output_root / "per_patient_change_summary.csv"
            ),
        },
        "methodology": {
            "primary_boundary": "clinical deterministic scorer before promotion v4",
            "promotion_v4": "supplementary downstream endpoint",
            "possible_pathogen_v1": "supplementary reporting endpoint",
            "answers": "posthoc only",
            "matching": contract["common_controls"]["matching_policy"],
            "ober": "excluded",
        },
    }
    write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
    parser.add_argument("--old-common", type=Path, default=DEFAULT_OLD_COMMON)
    parser.add_argument("--b-root", type=Path, default=DEFAULT_B)
    parser.add_argument("--c-root", type=Path, default=DEFAULT_C)
    parser.add_argument("--g-root", type=Path, default=DEFAULT_G)
    parser.add_argument("--f-root", type=Path, default=DEFAULT_F)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(
        json.dumps(
            run(
                args.contract,
                args.release,
                args.old_common,
                args.b_root,
                args.c_root,
                args.g_root,
                args.f_root,
                args.output,
            ),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

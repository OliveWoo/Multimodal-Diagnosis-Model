"""Run KH ablation D: frozen A pipeline plus the frozen v4C history route.

The preserved A pipeline predates the test-aware route interface.  This module
therefore creates an explicit compatibility bridge from A's deterministic and
Luna review tiers into the route schema, while preserving the frozen A
candidate names, old taxonomy profiles, selected-DNA analytical evidence, and
formal Picked decisions.  The only treatment is the v4C history route.

No new scorer, clinical scorer, promotion policy, Possible layer, or OBER step
is used in this experiment.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools.apply_test_aware_history_route_shadow import run as run_history_route
from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.pathogen_normalization import canonical_key
from tools.recalculate_kh_benchmark_metrics import match_type, pairwise_match, unique_names


DEFAULT_A_ROOT = Path("outputs/runs/2026-09-24_KH_ablation_A_frozen_old_baseline_v1")
DEFAULT_PATIENT_ROOT = Path("outputs/patient_info_KH_0728_2Days")
DEFAULT_POLICY = Path("rules/test_aware_history_route_v4c_combined.json")
DEFAULT_ANSWERS = Path(
    "outputs/runs/2026-09-18_KH_answer_revision_metrics/"
    "kh_answers_clinical_revision_20260918.csv"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-28_KH_ablation_D_legacy_v20_history_v4c_stable_v1"
)
SCHEMA_VERSION = "kh_ablation_d_legacy_history_bridge.v1"


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
    fieldnames = fields or (list(rows[0]) if rows else ["patient_id"])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def as_rank(value: Any) -> int | None:
    try:
        rank = int(float(str(value)))
    except (TypeError, ValueError):
        return None
    return rank if rank > 0 else None


def picked_keys(deterministic: dict[str, Any]) -> set[str]:
    return {
        canonical_key(item.get("organism_name"))
        for item in (deterministic.get("best_available_summary") or {}).get("picked_pathogens") or []
        if isinstance(item, dict) and item.get("organism_name")
    }


def review_tiers(review: dict[str, Any]) -> dict[str, tuple[str, dict[str, Any]]]:
    tiers: dict[str, tuple[str, dict[str, Any]]] = {}
    sections = (
        ("review_low_specificity", "review_low_specificity"),
        ("review_low_specificity_audit", "review_low_specificity"),
        ("review_context_needed", "review_context_needed"),
        ("review_high_priority", "review_high_priority"),
    )
    for section, decision in sections:
        for item in review.get(section) or []:
            if not isinstance(item, dict):
                continue
            key = canonical_key(item.get("organism_name"))
            if key:
                tiers[key] = (decision, item)
    return tiers


def specimen_context(candidate: dict[str, Any]) -> str:
    specimen = str(candidate.get("specimen_class") or "")
    if specimen.startswith("S1_"):
        return "sterile_or_systemic"
    if specimen.startswith("S2_"):
        return "lower_respiratory"
    return "unknown"


def bridge_candidate(
    candidate: dict[str, Any],
    *,
    decision: str,
    review_item: dict[str, Any] | None,
) -> dict[str, Any]:
    output = copy.deepcopy(candidate)
    if review_item:
        for key in ("taxonomy_profile", "classification"):
            if not output.get(key) and review_item.get(key):
                output[key] = copy.deepcopy(review_item[key])
    profile = output.get("taxonomy_profile") or {}
    rank = as_rank(output.get("rank_priority"))
    reads = float(output.get("reads") or 0)
    positive = bool(output.get("evidence_source") != "hospital_only" and reads > 0)
    hospital_detail = (
        output.get("hospital_evidence_detail")
        or (output.get("key_evidence") or {}).get("hospital_evidence_detail")
        or {}
    )
    output.update(
        {
            "decision": decision,
            "taxonomy_family": profile.get("primary_rule_family") or "unmapped_or_uncertain",
            "analytical_profile": {
                "best_rank_in_retained_universe": rank,
                "selected_positive_test_count": 1 if positive else 0,
                "reproducibility_axis": False,
                "cross_molecule_selected": False,
                "legacy_selected_dna_single_test": positive,
            },
            "specimen_context": specimen_context(output),
            "hospital_profile": {
                "exact_or_alias_evidence": bool(hospital_detail),
                "best_direct_level": hospital_detail.get("best_hospital_level"),
                "hospital_evidence_detail": hospital_detail,
            },
            "ablation_bridge": {
                "schema_version": SCHEMA_VERSION,
                "source_tier": decision,
                "single_selected_dna_test": True,
                "no_cross_molecule_or_reproducibility_axis": True,
                "old_taxonomy_profile_preserved": True,
            },
        }
    )
    return output


def build_bridge(
    *,
    a_root: Path,
    patient_ids: list[int],
    inventory: dict[tuple[int, str], dict[str, str]],
    output_dir: Path,
) -> dict[str, Any]:
    patient_output = output_dir / "patient_outputs"
    patient_output.mkdir(parents=True)
    counts: Counter[str] = Counter()
    source_hash_mismatches: list[dict[str, Any]] = []
    baseline_rows: list[dict[str, Any]] = []

    for patient in patient_ids:
        deterministic_record = inventory[(patient, "deterministic")]
        merged_record = inventory[(patient, "merged")]
        deterministic_path = Path(deterministic_record["path"])
        merged_path = Path(merged_record["path"])
        for role, record, path in (
            ("deterministic", deterministic_record, deterministic_path),
            ("merged", merged_record, merged_path),
        ):
            actual = sha256_file(path)
            if actual.lower() != str(record["sha256"]).lower():
                source_hash_mismatches.append(
                    {"patient_id": patient, "role": role, "expected": record["sha256"], "actual": actual}
                )
        deterministic = read_json(deterministic_path)
        merged = read_json(merged_path)
        review = merged.get("llm_missed_candidate_review") or {}
        tier_by_key = review_tiers(review)
        picked = picked_keys(deterministic)
        patient_candidates: list[dict[str, Any]] = []
        hospital_only: list[dict[str, Any]] = []
        seen: set[str] = set()
        for candidate in deterministic.get("pathogen_candidates") or []:
            if not isinstance(candidate, dict):
                continue
            key = canonical_key(candidate.get("organism_name"))
            if not key or key in seen:
                continue
            seen.add(key)
            tier, review_item = tier_by_key.get(key, ("review_low_specificity", {}))
            if key in picked:
                tier = "picked_shadow"
            bridged = bridge_candidate(candidate, decision=tier, review_item=review_item)
            target = hospital_only if candidate.get("evidence_source") == "hospital_only" else patient_candidates
            target.append(bridged)
            counts[f"baseline:{tier}"] += 1
            counts["candidates"] += 1
            baseline_rows.append(
                {
                    "patient_id": patient,
                    "organism_name": candidate.get("organism_name"),
                    "baseline_decision": tier,
                    "baseline_route": {
                        "picked_shadow": "Priority",
                        "review_high_priority": "Priority",
                        "review_context_needed": "Context",
                        "review_low_specificity": "Audit",
                    }.get(tier, "Hold"),
                    "taxonomy_family": bridged["taxonomy_family"],
                    "rank_priority": candidate.get("rank_priority"),
                    "reads": candidate.get("reads"),
                    "evidence_source": candidate.get("evidence_source", "mNGS_ranked"),
                }
            )
        # A small number of Luna-visible review items can originate from the
        # preserved review queue even when the deterministic scorer omitted
        # them from ``pathogen_candidates``.  They are part of A's frozen
        # combined endpoint and must remain in the D control universe.
        for key, (tier, review_item) in sorted(tier_by_key.items()):
            if key in seen:
                continue
            snapshot = copy.deepcopy(review_item.get("evidence_snapshot") or {})
            candidate = {
                **snapshot,
                "organism_name": review_item.get("organism_name"),
                "classification": review_item.get("classification")
                or snapshot.get("classification"),
                "taxonomy_profile": copy.deepcopy(review_item.get("taxonomy_profile") or {}),
                "evidence_source": snapshot.get("evidence_source") or "mNGS_ranked",
            }
            bridged = bridge_candidate(candidate, decision=tier, review_item=review_item)
            patient_candidates.append(bridged)
            seen.add(key)
            counts[f"baseline:{tier}"] += 1
            counts["candidates"] += 1
            baseline_rows.append(
                {
                    "patient_id": patient,
                    "organism_name": candidate.get("organism_name"),
                    "baseline_decision": tier,
                    "baseline_route": {
                        "review_high_priority": "Priority",
                        "review_context_needed": "Context",
                        "review_low_specificity": "Audit",
                    }.get(tier, "Hold"),
                    "taxonomy_family": bridged["taxonomy_family"],
                    "rank_priority": candidate.get("rank_priority"),
                    "reads": candidate.get("reads"),
                    "evidence_source": candidate.get("evidence_source"),
                }
            )
        payload = {
            "schema_version": SCHEMA_VERSION,
            "answer_blind": True,
            "patient_id": str(patient),
            "patient_organism_decisions": patient_candidates,
            "hospital_only_decisions": hospital_only,
            "source": {
                "frozen_A_contract": str((a_root / "baseline_contract.json").resolve()),
                "preserved_deterministic": str(deterministic_path.resolve()),
                "preserved_deterministic_sha256": sha256_file(deterministic_path),
                "preserved_merged": str(merged_path.resolve()),
                "preserved_merged_sha256": sha256_file(merged_path),
            },
            "constraints": [
                "Formal Picked and Luna review tiers are preserved from frozen A before history routing.",
                "The selected-DNA baseline contributes at most one positive mNGS test and no reproducibility axis.",
                "Embedded old taxonomy profiles are preserved; current taxonomy overlays are not loaded.",
            ],
        }
        write_json(
            patient_output / f"NGS_patient_{patient}_test_aware_deterministic_shadow.json",
            payload,
        )
        counts["patients"] += 1

    write_csv(output_dir / "baseline_candidate_decisions.csv", baseline_rows)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "answer_blind": True,
        "counts": dict(sorted(counts.items())),
        "source_hash_mismatches": source_hash_mismatches,
        "A_contract_sha256": sha256_file(a_root / "baseline_contract.json"),
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def group_decisions(rows: list[dict[str, Any]], allowed: set[str]) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = defaultdict(list)
    for row in rows:
        if str(row.get("decision") or row.get("baseline_decision")) in allowed:
            grouped[int(row["patient_id"])].append(str(row["organism_name"]))
    return {patient: unique_names(names) for patient, names in grouped.items()}


def answer_status(patient: int, name: str, answers: dict[int, list[str]]) -> str:
    gold = answers.get(patient, [])
    if not gold:
        return "unlabeled_patient"
    if any(match_type(name, answer, genus_relaxed=False) for answer in gold):
        return "matched_answer"
    return "unmatched_output"


def after_rows(history_root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted((history_root / "patient_outputs").glob("*.json")):
        payload = read_json(path)
        patient = int(payload["patient_id"])
        for item in payload.get("candidates") or []:
            rows.append(
                {
                    "patient_id": patient,
                    "organism_name": item.get("organism_name"),
                    "decision": item.get("history_routed_decision"),
                    "analytic_route": item.get("analytic_route"),
                    "history_adjusted_route": item.get("history_adjusted_route"),
                    "history_route_changed": bool(item.get("history_route_changed")),
                    "history_policy_family": item.get("history_policy_family"),
                    "history_rule_ids": "|".join(item.get("history_rule_ids") or []),
                    "history_evidence_ids": "|".join(item.get("history_evidence_ids") or []),
                    "consumed_evidence_ids": "|".join(item.get("consumed_evidence_ids") or []),
                    "blockers": "|".join(item.get("history_adjustment_blockers") or []),
                }
            )
    return rows


def per_patient_changes(
    baseline: list[dict[str, Any]],
    after: list[dict[str, Any]],
    answers: dict[int, list[str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    before_index = {
        (int(row["patient_id"]), canonical_key(row["organism_name"])): row for row in baseline
    }
    after_index = {
        (int(row["patient_id"]), canonical_key(row["organism_name"])): row for row in after
    }
    changes: list[dict[str, Any]] = []
    for key in sorted(set(before_index) | set(after_index)):
        old = before_index.get(key, {})
        new = after_index.get(key, {})
        old_decision = str(old.get("baseline_decision") or "hold")
        new_decision = str(new.get("decision") or "hold")
        old_route = str(old.get("baseline_route") or "Hold")
        new_route = str(new.get("history_adjusted_route") or "Hold")
        if old_decision == new_decision and old_route == new_route:
            continue
        patient = key[0]
        name = str(new.get("organism_name") or old.get("organism_name") or "")
        before_visible = old_decision in {"picked_shadow", "review_high_priority", "review_context_needed"}
        after_visible = new_decision in {"picked_shadow", "review_high_priority", "review_context_needed"}
        if not before_visible and after_visible:
            visibility = "report_added"
        elif before_visible and not after_visible:
            visibility = "report_removed"
        else:
            visibility = "report_tier_changed" if old_decision != new_decision else "route_only_changed"
        changes.append(
            {
                "patient_id": patient,
                "organism_name": name,
                "answer_status": answer_status(patient, name, answers),
                "visibility_change": visibility,
                "before_decision": old_decision,
                "after_decision": new_decision,
                "before_route": old_route,
                "after_route": new_route,
                "history_policy_family": new.get("history_policy_family", ""),
                "history_rule_ids": new.get("history_rule_ids", ""),
                "history_evidence_ids": new.get("history_evidence_ids", ""),
                "consumed_evidence_ids": new.get("consumed_evidence_ids", ""),
                "blockers": new.get("blockers", ""),
            }
        )

    baseline_visible = group_decisions(
        baseline, {"picked_shadow", "review_high_priority", "review_context_needed"}
    )
    after_visible = group_decisions(
        after, {"picked_shadow", "review_high_priority", "review_context_needed"}
    )
    summaries: list[dict[str, Any]] = []
    for patient in sorted(set(baseline_visible) | set(after_visible) | set(answers)):
        before_names = baseline_visible.get(patient, [])
        after_names = after_visible.get(patient, [])
        added = [name for name in after_names if not any(match_type(name, old, genus_relaxed=False) for old in before_names)]
        removed = [name for name in before_names if not any(match_type(name, new, genus_relaxed=False) for new in after_names)]
        gold = answers.get(patient, [])
        tp_gained = [name for name in added if any(match_type(name, answer, genus_relaxed=False) for answer in gold)]
        fp_added = [name for name in added if gold and name not in tp_gained]
        tp_lost = [name for name in removed if any(match_type(name, answer, genus_relaxed=False) for answer in gold)]
        fp_removed = [name for name in removed if gold and name not in tp_lost]
        summaries.append(
            {
                "patient_id": patient,
                "before_visible_count": len(before_names),
                "after_visible_count": len(after_names),
                "tp_gained": "|".join(tp_gained),
                "tp_lost": "|".join(tp_lost),
                "fp_added": "|".join(fp_added),
                "fp_removed": "|".join(fp_removed),
                "unlabeled_added": "|".join(added) if not gold else "",
                "unlabeled_removed": "|".join(removed) if not gold else "",
            }
        )
    return changes, summaries


def run(
    a_root: Path,
    patient_root: Path,
    policy_path: Path,
    answers_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    contract_path = a_root / "baseline_contract.json"
    contract = read_json(contract_path)
    patient_ids = [int(value) for value in contract["patient_selection"]]
    inventory = read_inventory(a_root / "artifact_inventory.csv")
    output_root.mkdir(parents=True)

    bridge_root = output_root / "legacy_A_route_bridge"
    bridge_summary = build_bridge(
        a_root=a_root,
        patient_ids=patient_ids,
        inventory=inventory,
        output_dir=bridge_root,
    )
    if bridge_summary["source_hash_mismatches"]:
        raise ValueError("Frozen A source hash validation failed")
    history_root = output_root / "history_route_v4c"
    history_summary = run_history_route(
        bridge_root,
        patient_root,
        history_root,
        policy_path=policy_path,
        phenotype_root=None,
        patients=set(patient_ids),
    )

    # Benchmark data are loaded only after the answer-blind bridge and history outputs exist.
    answers = load_answer_rows(answers_path)
    with (bridge_root / "baseline_candidate_decisions.csv").open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        baseline = list(csv.DictReader(handle))
    after = after_rows(history_root)
    evaluation_root = output_root / "evaluation"
    write_csv(evaluation_root / "final_candidate_decisions.csv", after)
    changes, summaries = per_patient_changes(baseline, after, answers)
    write_csv(evaluation_root / "per_patient_candidate_changes.csv", changes)
    write_csv(
        evaluation_root / "per_patient_summary.csv",
        summaries,
        [
            "patient_id", "before_visible_count", "after_visible_count", "tp_gained",
            "tp_lost", "fp_added", "fp_removed", "unlabeled_added", "unlabeled_removed",
        ],
    )

    strict = {"picked_shadow"}
    visible = {"picked_shadow", "review_high_priority", "review_context_needed"}
    metrics = {
        "A_bridge_strict_picked": evaluate_set(group_decisions(baseline, strict), answers),
        "D_history_strict_picked": evaluate_set(group_decisions(after, strict), answers),
        "A_bridge_complete_visible": evaluate_set(group_decisions(baseline, visible), answers),
        "D_history_complete_visible": evaluate_set(group_decisions(after, visible), answers),
    }
    expected = contract["frozen_metrics"]
    picked_metric = metrics["A_bridge_strict_picked"]
    visible_metric = metrics["A_bridge_complete_visible"]
    baseline_validation = {
        "strict_matches_frozen_A": (
            picked_metric["matched"] == expected["picked"]["matched"]
            and picked_metric["predicted_in_labeled_patients"]
            == expected["picked"]["outputs_in_labeled_patients"]
            and picked_metric["answer_count"] == expected["picked"]["answer_count"]
        ),
        "complete_visible_matches_frozen_A": (
            visible_metric["matched"] == expected["combined_picked_high_context"]["matched"]
            and visible_metric["predicted_in_labeled_patients"]
            == expected["combined_picked_high_context"]["outputs_in_labeled_patients"]
            and visible_metric["answer_count"]
            == expected["combined_picked_high_context"]["answer_count"]
        ),
    }
    if not all(baseline_validation.values()):
        raise ValueError(f"A compatibility bridge validation failed: {baseline_validation}")

    evaluation = {
        "schema_version": "kh_ablation_d_evaluation.v1",
        "scope": "Post-hoc evaluation; history generation was answer blind.",
        "metrics": metrics,
        "baseline_validation": baseline_validation,
        "route_change_count": len(changes),
        "route_change_counts": dict(Counter(row["visibility_change"] for row in changes)),
        "tp_gained": sum(bool(row["tp_gained"]) for row in summaries),
        "tp_lost": sum(bool(row["tp_lost"]) for row in summaries),
        "fp_added": sum(bool(row["fp_added"]) for row in summaries),
        "fp_removed": sum(bool(row["fp_removed"]) for row in summaries),
    }
    write_json(evaluation_root / "metrics.json", evaluation)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "experiment_id": "KH_ABLATION_D_LEGACY_V20_HISTORY_V4C_STABLE_V1",
        "status": "frozen_history_only_compatibility_shadow",
        "answer_blind_generation": True,
        "benchmark_evaluation_is_posthoc": True,
        "factors": {
            "mngs": contract["factors"]["mngs"],
            "scorer": contract["factors"]["scorer"],
            "taxonomy": contract["factors"]["taxonomy"],
            "history": "test_aware_history_route_v4c_combined",
        },
        "common_contract": {
            "patient_count": len(patient_ids),
            "patient_selection": [str(value) for value in patient_ids],
            "A_contract_path": str(contract_path.resolve()),
            "A_contract_sha256": sha256_file(contract_path),
            "answer_path": str(answers_path.resolve()),
            "answer_sha256": sha256_file(answers_path),
        },
        "history_policy": {
            "path": str(policy_path.resolve()),
            "sha256": sha256_file(policy_path),
            "policy_id": history_summary["policy_id"],
            "resolved_policy_sha256": history_summary["resolved_policy_sha256"],
        },
        "implementation": {
            "bridge_runner_path": str(Path(__file__).resolve()),
            "bridge_runner_sha256": sha256_file(Path(__file__)),
            "history_route_tool_path": str(
                Path("tools/apply_test_aware_history_route_shadow.py").resolve()
            ),
            "history_route_tool_sha256": sha256_file(
                Path("tools/apply_test_aware_history_route_shadow.py")
            ),
        },
        "baseline_validation": baseline_validation,
        "counts": history_summary["counts"],
        "metrics": metrics,
        "route_change_counts": evaluation["route_change_counts"],
        "outputs": {
            "bridge_summary_sha256": sha256_file(bridge_root / "summary.json"),
            "history_summary_sha256": sha256_file(history_root / "summary.json"),
            "candidate_changes_sha256": sha256_file(evaluation_root / "per_patient_candidate_changes.csv"),
            "patient_summary_sha256": sha256_file(evaluation_root / "per_patient_summary.csv"),
            "metrics_sha256": sha256_file(evaluation_root / "metrics.json"),
        },
        "limitations": [
            "A predates the test-aware route interface, so D uses an explicit compatibility bridge whose baseline metrics must reproduce frozen A.",
            "The old selected-DNA input provides one analytical test and cannot create DNA/RNA or repeated-test reproducibility.",
            "History may change reporting routes but never creates Picked by itself.",
            "New scorer, new taxonomy, promotion v4, Possible, and OBER are excluded from D.",
        ],
    }
    write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-root", type=Path, default=DEFAULT_A_ROOT)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.a_root, args.patient_root, args.policy, args.answers, args.output),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

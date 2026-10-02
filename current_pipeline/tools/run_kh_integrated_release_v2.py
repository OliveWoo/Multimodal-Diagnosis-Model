"""Build the clean 33-patient KH integrated release v2 from fixed evidence.

Decision stages are completed and hashed before benchmark answers are loaded.
The release rebuilds compact taxonomy, the analytical scorer, the 48-hour
timeline, phenotype packets, v4C history routing, clinical scoring, evidence
rehydration, promotion v8, unified reporting v3 Route A, and mutually
exclusive fallback.  OBER/Luna are deliberately outside this entry point.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from tools import apply_test_aware_history_route_shadow as history
from tools import build_compact_multi_assay_entry_shadow as compact
from tools import build_deterministic_summary as deterministic_summary
from tools import build_index_event_clinical_timeline_shadow as timeline
from tools import build_phenotype_event_packets as phenotype_packets
from tools import build_test_aware_candidate_phenotype_shadow as candidate_phenotype
from tools import build_test_aware_clinical_deterministic_shadow as clinical
from tools import build_test_aware_deterministic_shadow_scorer as scorer
from tools import build_test_aware_possible_pathogen_shadow as reporting
from tools import evaluate_test_aware_clinical_deterministic_shadow as eval_clinical
from tools import evaluate_test_aware_possible_pathogen_shadow as eval_reporting
from tools import promote_test_aware_clinical_high_shadow_v8 as promotion
from tools import rehydrate_clinical_hospital_evidence_shadow as rehydrate
from tools.pathogen_normalization import canonical_key
from tools.run_kh_integrated_decision_pipeline import build_integrated_results


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs/runs/2026-09-30_KH_integrated_release_v2_2"
DEFAULT_INVENTORY = ROOT / "outputs/runs/2026-09-23_KH_multi_assay_mngs_inventory_source_contract_v3"
DEFAULT_PATIENT_ROOT = ROOT / "outputs/patient_info_KH_0728_2Days"
DEFAULT_PHENOTYPE_ROOT = ROOT / "outputs/phenotype_decision_evidence_v1"
DEFAULT_LINKAGE_AUDIT = ROOT / "outputs/runs/2026-09-22_KH_phenotype_patient_linkage_audit_v1/phenotype_patient_linkage_audit.json"
DEFAULT_ANSWERS = ROOT / "outputs/runs/2026-09-18_KH_answer_revision_metrics/kh_answers_clinical_revision_20260918.csv"
DEFAULT_REFERENCE = ROOT / "outputs/runs/2026-09-29_KH_integrated_decision_pipeline_v1"
FROZEN_HOSPITAL_SUMMARY_SUFFIX = "final_summary_selected_dna_20260918"
CURRENT_SUMMARY_SUFFIX = "release_input_current_host"

COMPACT_POLICY = ROOT / "rules/manual_style_multi_assay_entry_v2.json"
ANALYTICAL_POLICY = ROOT / "rules/test_aware_deterministic_shadow_v3_route_a.json"
HISTORY_POLICY = ROOT / "rules/test_aware_history_route_v4c_combined.json"
CLINICAL_POLICY = ROOT / "rules/test_aware_clinical_deterministic_shadow_v1.json"
PROMOTION_POLICY = ROOT / "rules/test_aware_clinical_promotion_v8_respiratory_group_convergence.json"
REPORTING_POLICY = ROOT / "rules/test_aware_unified_reporting_v3_route_a.json"
IDENTITY_ALIAS_FILE = ROOT / "rules/pathogen_aliases.json"
TAXONOMY_BASE_FILE = ROOT / "rules/organism_taxonomy_rules.json"
TAXONOMY_GAP_FILE = ROOT / "rules/organism_taxonomy_gap_resolution_v1.json"
TAXONOMY_AUTO_FILE = Path(
    os.environ.get(
        "MNGS_TAXONOMY_AUTO_OVERLAY",
        ROOT / "rules/organism_taxonomy_auto_v1.json",
    )
)
TAXONOMY_REVIEWED_FILE = ROOT / "rules/organism_taxonomy_reviewed_v1.json"
EPISODE_IDENTITY_FILE = ROOT / "rules/phenotype_episode_identity_v1.json"
SCHEMA_VERSION = "kh_integrated_release.v2"


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else [
        "patient_id", "organism_name", "change_type", "before_tier", "after_tier"
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _candidate_key(row: dict[str, Any]) -> tuple[str, str]:
    return (
        str(row.get("patient_id") or ""),
        canonical_key(row.get("organism_name")),
    )


def _benchmark_index(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    if not path.is_file():
        return {}
    return {_candidate_key(row): row for row in read_csv(path)}


def build_before_after(
    reference_root: Path, current_root: Path, output_dir: Path
) -> dict[str, Any]:
    before_rows = read_csv(reference_root / "integrated_decisions.csv")
    after_rows = read_csv(current_root / "integrated_decisions.csv")
    before = {_candidate_key(row): row for row in before_rows}
    after = {_candidate_key(row): row for row in after_rows}
    before_benchmark = _benchmark_index(
        reference_root / "evaluation/reporting/prediction_audit.csv"
    )
    after_benchmark = _benchmark_index(
        current_root / "evaluation/reporting/prediction_audit.csv"
    )
    changes: list[dict[str, Any]] = []
    for key in sorted(set(before) | set(after)):
        old = before.get(key)
        new = after.get(key)
        old_selected = str((old or {}).get("selected_for_complete_report") or "").casefold() == "true"
        new_selected = str((new or {}).get("selected_for_complete_report") or "").casefold() == "true"
        old_tier = str((old or {}).get("final_reporting_tier") or "absent")
        new_tier = str((new or {}).get("final_reporting_tier") or "absent")
        old_route = str((old or {}).get("reporting_route") or "")
        new_route = str((new or {}).get("reporting_route") or "")
        old_clinical = str((old or {}).get("post_promotion_clinical_decision") or "")
        new_clinical = str((new or {}).get("post_promotion_clinical_decision") or "")
        if (old_selected, old_tier, old_route, old_clinical) == (
            new_selected, new_tier, new_route, new_clinical
        ):
            continue
        if old is None:
            change_type = "candidate_added"
        elif new is None:
            change_type = "candidate_removed"
        elif not old_selected and new_selected:
            change_type = "added_to_complete_report"
        elif old_selected and not new_selected:
            change_type = "removed_from_complete_report"
        else:
            change_type = "tier_or_route_changed"
        benchmark = after_benchmark.get(key) or before_benchmark.get(key) or {}
        changes.append({
            "patient_id": key[0],
            "organism_name": (new or old or {}).get("organism_name"),
            "change_type": change_type,
            "before_selected": old_selected,
            "after_selected": new_selected,
            "before_tier": old_tier,
            "after_tier": new_tier,
            "before_clinical_decision": old_clinical,
            "after_clinical_decision": new_clinical,
            "before_reporting_route": old_route,
            "after_reporting_route": new_route,
            "posthoc_outcome": benchmark.get("posthoc_outcome", ""),
            "matched_answer": benchmark.get("matched_answer", ""),
        })
    write_csv(output_dir / "per_patient_before_after.csv", changes)
    return {
        "reference_candidate_count": len(before),
        "current_candidate_count": len(after),
        "changed_candidate_count": len(changes),
        "change_type_counts": {
            value: sum(row["change_type"] == value for row in changes)
            for value in sorted({row["change_type"] for row in changes})
        },
        "changes_file": str((output_dir / "per_patient_before_after.csv").resolve()),
    }


def validate_episode_separation(patient_ids: set[int]) -> dict[str, Any]:
    policy = read_json(EPISODE_IDENTITY_FILE)
    invalid = [
        item for item in policy.get("clusters") or []
        if item.get("action") != "keep_separate"
        or item.get("episode_relationship") != "distinct_pneumonia_episodes_confirmed"
    ]
    cohort_members = sorted({
        int(member)
        for item in policy.get("clusters") or []
        for member in item.get("members") or []
        if int(member) in patient_ids
    })
    return {
        "passed": not invalid,
        "rule_file": str(EPISODE_IDENTITY_FILE.resolve()),
        "rule_sha256": sha256_file(EPISODE_IDENTITY_FILE),
        "cohort_episode_ids_kept_separate": cohort_members,
        "cross_episode_evidence_merge_performed": False,
        "invalid_rule_count": len(invalid),
    }


def _component(path: Path, role: str) -> dict[str, Any]:
    return {"role": role, "path": str(path.resolve()), "sha256": sha256_file(path)}


def _taxonomy_components() -> list[dict[str, Any]]:
    """Return every central-taxonomy layer needed to reproduce routing."""

    return [
        _component(TAXONOMY_BASE_FILE, "central_taxonomy_base"),
        _component(TAXONOMY_GAP_FILE, "taxonomy_gap_resolution_overlay"),
        _component(TAXONOMY_AUTO_FILE, "auto_taxonomy_overlay"),
        _component(TAXONOMY_REVIEWED_FILE, "reviewed_taxonomy_overlay"),
    ]


def build_current_source_summaries(
    patient_root: Path, output_dir: Path
) -> dict[str, Any]:
    """Freeze reviewed hospital evidence plus current host context in the release.

    Stage 4 reviewed the hospital-evidence semantics carried by the frozen
    selected-DNA summary.  Stage 5 subsequently improved CBC, underlying, and
    medication parsing.  Rebuilding every module with a newer generic summary
    builder would change already-reviewed hospital semantics, so this stage
    makes an explicit, narrow host-only overlay and records every changed host
    field.  Patient folders are never mutated.
    """
    patient_dirs = deterministic_summary.collect_patient_dirs([patient_root], None)
    if len(patient_dirs) != 33:
        raise ValueError(
            f"Expected 33 patient directories for current summaries, got {len(patient_dirs)}"
        )
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Current-summary output already exists: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    tier_counts: dict[str, int] = {}
    support_counts = {"true": 0, "false": 0}
    host_changes: list[dict[str, Any]] = []
    files: list[dict[str, str]] = []
    for patient_dir in patient_dirs:
        patient_name = deterministic_summary.extract_patient_identifier(patient_dir)
        patient_id = patient_name.removeprefix("NGS_patient_")
        frozen_path = (
            patient_dir
            / "summary_outputs"
            / f"NGS_patient_{patient_id}_{FROZEN_HOSPITAL_SUMMARY_SUFFIX}.json"
        )
        if not frozen_path.is_file():
            raise FileNotFoundError(frozen_path)
        frozen = read_json(frozen_path)
        current = deterministic_summary.build_summary(patient_dir)
        payload = copy.deepcopy(frozen)
        before_host = copy.deepcopy(payload.get("host_context") or {})
        after_host = copy.deepcopy(current.get("host_context") or {})
        payload["host_context"] = after_host
        payload.setdefault("module_summaries", {})["cbc_other_lab"] = copy.deepcopy(
            (current.get("module_summaries") or {}).get("cbc_other_lab") or {}
        )
        payload.setdefault("evidence_preservation", {})[
            "structured_host_evidence"
        ] = copy.deepcopy(
            (current.get("evidence_preservation") or {}).get(
                "structured_host_evidence"
            )
            or {}
        )
        current_sources = current.get("source_files") or {}
        for key in (
            "cbc_other_lab",
            "raw_cbc",
            "raw_other_lab",
            "raw_underlying",
            "raw_admission_diagnosis",
        ):
            if key in current_sources:
                payload.setdefault("source_files", {})[key] = current_sources[key]
        tracked_fields = (
            "immunocompromise_tier",
            "acute_instability_tier",
            "host_vulnerability_tier",
            "opportunistic_host_support",
            "opportunistic_host_support_status",
            "opportunistic_host_review_context",
        )
        changed_fields = {
            field: {"before": before_host.get(field), "after": after_host.get(field)}
            for field in tracked_fields
            if before_host.get(field) != after_host.get(field)
        }
        if changed_fields:
            host_changes.append({
                "patient_id": patient_id,
                "changed_fields": changed_fields,
            })
        payload["release_input_overlay"] = {
            "overlay_type": "current_deterministic_host_context_only",
            "frozen_hospital_summary": str(frozen_path.resolve()),
            "frozen_hospital_summary_sha256": sha256_file(frozen_path),
            "hospital_evidence_semantics_changed": False,
            "host_context_rebuilt_from_current_sources": True,
            "changed_host_fields": changed_fields,
        }
        destination = output_dir / (
            f"NGS_patient_{patient_id}_{CURRENT_SUMMARY_SUFFIX}.json"
        )
        write_json(destination, payload)
        host = payload.get("host_context") or {}
        tier = str(host.get("immunocompromise_tier") or "Unknown")
        tier_counts[tier] = tier_counts.get(tier, 0) + 1
        support_key = "true" if host.get("opportunistic_host_support") is True else "false"
        support_counts[support_key] += 1
        files.append({
            "patient_id": patient_id,
            "path": str(destination.resolve()),
            "sha256": sha256_file(destination),
        })
    summary = {
        "schema_version": "kh_current_source_summary_freeze.v1",
        "answer_blind": True,
        "patient_count": len(files),
        "summary_suffix": CURRENT_SUMMARY_SUFFIX,
        "source_patient_root": str(patient_root.resolve()),
        "frozen_hospital_summary_suffix": FROZEN_HOSPITAL_SUMMARY_SUFFIX,
        "immunocompromise_tier_counts": dict(sorted(tier_counts.items())),
        "opportunistic_host_support_counts": support_counts,
        "host_context_change_count": len(host_changes),
        "host_context_changes": host_changes,
        "patient_files": files,
        "constraints": [
            "Reviewed frozen hospital-evidence semantics are preserved.",
            "Only host context and its CBC/underlying/medication provenance are rebuilt from current sources.",
            "The release does not mutate patient summary_outputs directories.",
            "Host context cannot independently create an organism identity or Picked result.",
        ],
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def run(
    output_root: Path = DEFAULT_OUTPUT,
    *,
    inventory_root: Path = DEFAULT_INVENTORY,
    patient_root: Path = DEFAULT_PATIENT_ROOT,
    phenotype_root: Path = DEFAULT_PHENOTYPE_ROOT,
    linkage_audit: Path = DEFAULT_LINKAGE_AUDIT,
    answers: Path = DEFAULT_ANSWERS,
    reference_root: Path = DEFAULT_REFERENCE,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty release: {output_root}")
    required = [
        inventory_root / "inventory_summary.json", patient_root, phenotype_root,
        linkage_audit, answers, reference_root / "integrated_decisions.csv",
        COMPACT_POLICY, ANALYTICAL_POLICY, HISTORY_POLICY, CLINICAL_POLICY,
        PROMOTION_POLICY, REPORTING_POLICY, IDENTITY_ALIAS_FILE,
        TAXONOMY_REVIEWED_FILE, EPISODE_IDENTITY_FILE,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required release inputs: {missing}")
    output_root.mkdir(parents=True)

    stage = {name: output_root / name for name in (
        "00_current_source_summary", "01_compact_taxonomy", "02_analytical_scorer", "03_clinical_timeline",
        "04_candidate_phenotype", "05_phenotype_event_packets", "06_history_v4c",
        "07_clinical_scorer", "08_evidence_rehydrated", "09_promotion_v8",
        "10_reporting_v3_route_a",
    )}

    stage_summaries: dict[str, Any] = {}
    stage_summaries["current_source_summary"] = build_current_source_summaries(
        patient_root, stage["00_current_source_summary"]
    )
    stage_summaries["compact_taxonomy"] = compact.build_shadow(
        inventory_root, stage["01_compact_taxonomy"], COMPACT_POLICY
    )
    stage_summaries["analytical_scorer"] = scorer.build_shadow(
        stage["01_compact_taxonomy"], patient_root, stage["02_analytical_scorer"],
        summary_suffix=CURRENT_SUMMARY_SUFFIX,
        summary_root=stage["00_current_source_summary"],
        policy_path=ANALYTICAL_POLICY,
    )
    stage_summaries["clinical_timeline"] = timeline.build_shadow(
        patient_root, stage["03_clinical_timeline"],
        multi_assay_root=inventory_root / "patients",
    )
    stage_summaries["candidate_phenotype"] = candidate_phenotype.run(
        stage["02_analytical_scorer"], phenotype_root,
        stage["04_candidate_phenotype"], linkage_mode="same_number_shadow",
    )
    stage_summaries["phenotype_event_packets"] = phenotype_packets.run(
        stage["03_clinical_timeline"], phenotype_root, linkage_audit,
        stage["04_candidate_phenotype"], stage["05_phenotype_event_packets"],
    )
    stage_summaries["history_v4c"] = history.run(
        stage["02_analytical_scorer"], patient_root, stage["06_history_v4c"],
        policy_path=HISTORY_POLICY, phenotype_root=stage["04_candidate_phenotype"],
    )
    stage_summaries["clinical_scorer"] = clinical.run(
        stage["02_analytical_scorer"], stage["03_clinical_timeline"],
        stage["05_phenotype_event_packets"], stage["07_clinical_scorer"],
        policy_path=CLINICAL_POLICY, history_route_root=stage["06_history_v4c"],
    )
    stage_summaries["evidence_rehydrated"] = rehydrate.run(
        stage["07_clinical_scorer"], stage["02_analytical_scorer"],
        stage["08_evidence_rehydrated"],
    )
    stage_summaries["promotion_v8"] = promotion.run(
        stage["08_evidence_rehydrated"], stage["09_promotion_v8"],
        policy_path=PROMOTION_POLICY,
    )
    stage_summaries["reporting_v3_route_a"] = reporting.run(
        stage["09_promotion_v8"], stage["10_reporting_v3_route_a"],
        policy_path=REPORTING_POLICY,
    )
    stage_summaries["integrated_output"] = build_integrated_results(
        stage["10_reporting_v3_route_a"], output_root
    )

    patient_ids = {
        int(path.name.split("_")[2])
        for path in (stage["02_analytical_scorer"] / "patient_outputs").glob(
            "NGS_patient_*_test_aware_deterministic_shadow.json"
        )
    }
    patient_count_assertions = {
        name: summary.get("patient_count") == 33
        for name, summary in stage_summaries.items()
        if isinstance(summary, dict) and "patient_count" in summary
    }
    episode_validation = validate_episode_separation(patient_ids)
    answer_blind_freeze = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "answers_loaded": False,
        "patient_ids": sorted(patient_ids),
        "patient_count": len(patient_ids),
        "patient_count_assertions": patient_count_assertions,
        "episode_separation": episode_validation,
        "decision_artifacts": {
            "integrated_decisions_sha256": sha256_file(output_root / "integrated_decisions.csv"),
            "complete_report_sha256": sha256_file(stage["10_reporting_v3_route_a"] / "complete_report.csv"),
            "reporting_summary_sha256": sha256_file(stage["10_reporting_v3_route_a"] / "summary.json"),
        },
        "ober_included": False,
        "passed": (
            len(patient_ids) == 33
            and all(patient_count_assertions.values())
            and episode_validation["passed"]
        ),
    }
    write_json(output_root / "answer_blind_freeze.json", answer_blind_freeze)
    if not answer_blind_freeze["passed"]:
        raise RuntimeError("Answer-blind release completeness assertions failed")

    evaluation_root = output_root / "evaluation"
    clinical_evaluation = eval_clinical.run(
        stage["09_promotion_v8"], answers, evaluation_root / "clinical"
    )
    reporting_evaluation = eval_reporting.run(
        stage["10_reporting_v3_route_a"], answers, evaluation_root / "reporting"
    )
    comparison = build_before_after(reference_root, output_root, output_root / "comparison")
    reference_metrics = read_json(reference_root / "evaluation/reporting/metrics.json")
    benchmark_comparison = {
        "scope": "posthoc_only_after_answer_blind_freeze",
        "reference_root": str(reference_root.resolve()),
        "reference_complete_report": reference_metrics["metrics"]["complete_report"],
        "current_complete_report": reporting_evaluation["metrics"]["complete_report"],
        "before_after": comparison,
    }
    write_json(output_root / "comparison/benchmark_comparison.json", benchmark_comparison)

    components = [
        _component(COMPACT_POLICY, "compact_entry_policy"),
        _component(ANALYTICAL_POLICY, "analytical_scorer_policy_route_a"),
        _component(HISTORY_POLICY, "history_v4c_policy"),
        _component(CLINICAL_POLICY, "clinical_scorer_policy"),
        _component(PROMOTION_POLICY, "promotion_v8_policy"),
        _component(REPORTING_POLICY, "reporting_v3_route_a_policy"),
        _component(IDENTITY_ALIAS_FILE, "production_identity_aliases"),
        *_taxonomy_components(),
        _component(EPISODE_IDENTITY_FILE, "episode_separation_policy"),
    ]
    for module, role in (
        (compact, "compact_implementation"),
        (deterministic_summary, "current_source_summary_implementation"),
        (scorer, "analytical_scorer_implementation"),
        (timeline, "timeline_implementation"),
        (candidate_phenotype, "candidate_phenotype_implementation"),
        (phenotype_packets, "phenotype_packet_implementation"),
        (history, "history_router_implementation"),
        (clinical, "clinical_scorer_implementation"),
        (rehydrate, "evidence_rehydration_implementation"),
        (promotion, "promotion_v8_implementation"),
        (reporting, "reporting_implementation"),
    ):
        components.append(_component(Path(module.__file__), role))

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "release_id": output_root.name,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "retrospective_integrated_release_candidate",
        "answer_blind_decision_stages": True,
        "answer_blind_freeze_sha256": sha256_file(output_root / "answer_blind_freeze.json"),
        "execution_boundary": (
            "Rebuilt from the fixed 33-patient multi-assay per-test inventory and current "
            "hospital-source summaries through taxonomy, v4C history, clinical scoring, "
            "promotion v8, reporting v3 Route A, and mutually exclusive fallback."
        ),
        "input": {
            "inventory_root": str(inventory_root.resolve()),
            "inventory_summary_sha256": sha256_file(inventory_root / "inventory_summary.json"),
            "patient_root": str(patient_root.resolve()),
            "phenotype_root": str(phenotype_root.resolve()),
            "linkage_audit": str(linkage_audit.resolve()),
            "linkage_audit_sha256": sha256_file(linkage_audit),
        },
        "components": components,
        "stage_roots": {name: str(path.resolve()) for name, path in stage.items()},
        "patient_count": len(patient_ids),
        "candidate_count": stage_summaries["integrated_output"]["candidate_count"],
        "episode_separation": episode_validation,
        "identity_aliases_active": True,
        "route_a_possible_active": True,
        "route_b_fallback_active": True,
        "fallback_mutually_exclusive_recomputed": True,
        "ober_included": False,
        "evaluation_attached_posthoc": True,
        "constraints": [
            "This is retrospective development evaluation, not prospective or external validation.",
            "Phenotype links without source-holder verification remain contextual and cannot independently create Picked.",
            "The four confirmed same-patient export clusters remain separate pneumonia episodes.",
            "OBER/Luna are excluded from every decision stage in this release.",
        ],
    }
    write_json(output_root / "run_manifest.json", manifest)

    summary = {
        "schema_version": SCHEMA_VERSION,
        "release_id": manifest["release_id"],
        "answer_blind_freeze_passed": answer_blind_freeze["passed"],
        "patient_count": len(patient_ids),
        "candidate_count": manifest["candidate_count"],
        "stage_summaries": stage_summaries,
        "clinical_evaluation": clinical_evaluation,
        "reporting_evaluation": reporting_evaluation,
        "benchmark_comparison": benchmark_comparison,
        "ober_included": False,
    }
    write_json(output_root / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--inventory-root", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--phenotype-root", type=Path, default=DEFAULT_PHENOTYPE_ROOT)
    parser.add_argument("--linkage-audit", type=Path, default=DEFAULT_LINKAGE_AUDIT)
    parser.add_argument("--answers", type=Path, default=DEFAULT_ANSWERS)
    parser.add_argument("--reference-root", type=Path, default=DEFAULT_REFERENCE)
    args = parser.parse_args()
    print(json.dumps(run(
        args.output_root,
        inventory_root=args.inventory_root,
        patient_root=args.patient_root,
        phenotype_root=args.phenotype_root,
        linkage_audit=args.linkage_audit,
        answers=args.answers,
        reference_root=args.reference_root,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

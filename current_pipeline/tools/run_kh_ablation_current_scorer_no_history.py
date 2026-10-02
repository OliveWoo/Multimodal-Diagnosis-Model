"""Run current KH deterministic decision engine without the new history route.

This runner is used for the refreshed 33-patient ablation C and G groups:

* C: new multi-assay mNGS, current scorer bundle, old compact taxonomy, no history.
* G: new multi-assay mNGS, current scorer bundle, current v2.6 taxonomy, no history.

The scorer bundle is frozen as analytical scorer v3 Route A, clinical scorer v1,
promotion v8, and unified reporting v3 Route A.  Decisions are frozen before the
60-organism development answers are loaded for post-hoc evaluation.  OBER/Luna
are outside this runner.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

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
from tools import run_kh_integrated_release_v2 as integrated


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OLD_COMPACT = (
    ROOT
    / "outputs/runs/2026-09-23_KH_compact_multi_assay_entry_source_contract_v5"
)
DEFAULT_INVENTORY = integrated.DEFAULT_INVENTORY
DEFAULT_PATIENT_ROOT = integrated.DEFAULT_PATIENT_ROOT
DEFAULT_PHENOTYPE_ROOT = integrated.DEFAULT_PHENOTYPE_ROOT
DEFAULT_LINKAGE_AUDIT = integrated.DEFAULT_LINKAGE_AUDIT
DEFAULT_ANSWERS = integrated.DEFAULT_ANSWERS
SCHEMA_VERSION = "kh_ablation_current_scorer_no_history.v1"


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def component(path: Path, role: str) -> dict[str, Any]:
    return {
        "role": role,
        "path": str(path.resolve()),
        "sha256": integrated.sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def run(
    *,
    group: str,
    taxonomy_mode: str,
    output_root: Path,
    old_compact_root: Path = DEFAULT_OLD_COMPACT,
    inventory_root: Path = DEFAULT_INVENTORY,
    patient_root: Path = DEFAULT_PATIENT_ROOT,
    phenotype_root: Path = DEFAULT_PHENOTYPE_ROOT,
    linkage_audit: Path = DEFAULT_LINKAGE_AUDIT,
    answers: Path = DEFAULT_ANSWERS,
) -> dict[str, Any]:
    if group not in {"C", "G"}:
        raise ValueError("group must be C or G")
    expected_mode = "old" if group == "C" else "current_v2_6"
    if taxonomy_mode != expected_mode:
        raise ValueError(f"Group {group} requires taxonomy_mode={expected_mode}")
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")

    required = [
        inventory_root / "inventory_summary.json",
        patient_root,
        phenotype_root,
        linkage_audit,
        answers,
        integrated.COMPACT_POLICY,
        integrated.ANALYTICAL_POLICY,
        integrated.CLINICAL_POLICY,
        integrated.PROMOTION_POLICY,
        integrated.REPORTING_POLICY,
    ]
    if taxonomy_mode == "old":
        required.extend(
            [old_compact_root / "summary.json", old_compact_root / "patient_packets"]
        )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing ablation inputs: {missing}")

    output_root.mkdir(parents=True)
    stage = {
        name: output_root / name
        for name in (
            "00_current_source_summary",
            "01_compact_taxonomy",
            "02_analytical_scorer",
            "03_clinical_timeline",
            "04_candidate_phenotype",
            "05_phenotype_event_packets",
            "07_clinical_scorer",
            "08_evidence_rehydrated",
            "09_promotion_v8",
            "10_reporting_v3_route_a",
        )
    }

    summaries: dict[str, Any] = {}
    summaries["current_source_summary"] = integrated.build_current_source_summaries(
        patient_root, stage["00_current_source_summary"]
    )
    if taxonomy_mode == "current_v2_6":
        summaries["compact_taxonomy"] = integrated.compact.build_shadow(
            inventory_root,
            stage["01_compact_taxonomy"],
            integrated.COMPACT_POLICY,
        )
        compact_root = stage["01_compact_taxonomy"]
    else:
        compact_root = old_compact_root
        summaries["compact_taxonomy"] = integrated.read_json(
            old_compact_root / "summary.json"
        )

    summaries["analytical_scorer"] = scorer.build_shadow(
        compact_root,
        patient_root,
        stage["02_analytical_scorer"],
        summary_suffix=integrated.CURRENT_SUMMARY_SUFFIX,
        summary_root=stage["00_current_source_summary"],
        policy_path=integrated.ANALYTICAL_POLICY,
    )
    summaries["clinical_timeline"] = timeline.build_shadow(
        patient_root,
        stage["03_clinical_timeline"],
        multi_assay_root=inventory_root / "patients",
    )
    summaries["candidate_phenotype"] = candidate_phenotype.run(
        stage["02_analytical_scorer"],
        phenotype_root,
        stage["04_candidate_phenotype"],
        linkage_mode="same_number_shadow",
    )
    summaries["phenotype_event_packets"] = phenotype_packets.run(
        stage["03_clinical_timeline"],
        phenotype_root,
        linkage_audit,
        stage["04_candidate_phenotype"],
        stage["05_phenotype_event_packets"],
    )
    summaries["clinical_scorer"] = clinical.run(
        stage["02_analytical_scorer"],
        stage["03_clinical_timeline"],
        stage["05_phenotype_event_packets"],
        stage["07_clinical_scorer"],
        policy_path=integrated.CLINICAL_POLICY,
        history_route_root=None,
    )
    summaries["evidence_rehydrated"] = rehydrate.run(
        stage["07_clinical_scorer"],
        stage["02_analytical_scorer"],
        stage["08_evidence_rehydrated"],
    )
    summaries["promotion_v8"] = promotion.run(
        stage["08_evidence_rehydrated"],
        stage["09_promotion_v8"],
        policy_path=integrated.PROMOTION_POLICY,
    )
    summaries["reporting_v3_route_a"] = reporting.run(
        stage["09_promotion_v8"],
        stage["10_reporting_v3_route_a"],
        policy_path=integrated.REPORTING_POLICY,
    )
    summaries["integrated_output"] = integrated.build_integrated_results(
        stage["10_reporting_v3_route_a"], output_root
    )

    patient_ids = {
        int(path.name.split("_")[2])
        for path in (stage["02_analytical_scorer"] / "patient_outputs").glob(
            "NGS_patient_*_test_aware_deterministic_shadow.json"
        )
    }
    patient_assertions = {
        name: summary.get("patient_count") == 33
        for name, summary in summaries.items()
        if isinstance(summary, dict) and "patient_count" in summary
    }
    clinical_history_root = summaries["clinical_scorer"].get("history_route_root")
    freeze = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "schema_version": SCHEMA_VERSION,
        "ablation_group": group,
        "answers_loaded": False,
        "patient_ids": sorted(patient_ids),
        "patient_count": len(patient_ids),
        "patient_count_assertions": patient_assertions,
        "history_route_root": clinical_history_root,
        "history_disabled": clinical_history_root is None,
        "decision_artifacts": {
            "integrated_decisions_sha256": integrated.sha256_file(
                output_root / "integrated_decisions.csv"
            ),
            "complete_report_sha256": integrated.sha256_file(
                stage["10_reporting_v3_route_a"] / "complete_report.csv"
            ),
        },
        "passed": (
            len(patient_ids) == 33
            and all(patient_assertions.values())
            and clinical_history_root is None
        ),
    }
    write_json(output_root / "answer_blind_freeze.json", freeze)
    if not freeze["passed"]:
        raise RuntimeError(f"Group {group} answer-blind freeze failed")

    evaluation_root = output_root / "evaluation"
    clinical_evaluation = eval_clinical.run(
        stage["09_promotion_v8"], answers, evaluation_root / "clinical"
    )
    reporting_evaluation = eval_reporting.run(
        stage["10_reporting_v3_route_a"], answers, evaluation_root / "reporting"
    )

    components = [
        component(integrated.ANALYTICAL_POLICY, "analytical_scorer_v3_route_a"),
        component(integrated.CLINICAL_POLICY, "clinical_scorer_v1"),
        component(integrated.PROMOTION_POLICY, "promotion_v8"),
        component(integrated.REPORTING_POLICY, "reporting_v3_route_a"),
        component(compact_root / "summary.json", "compact_taxonomy_input"),
    ]
    if taxonomy_mode == "current_v2_6":
        components.extend(
            component(item, role)
            for item, role in (
                (integrated.TAXONOMY_BASE_FILE, "taxonomy_base"),
                (integrated.TAXONOMY_GAP_FILE, "taxonomy_gap_overlay"),
                (integrated.TAXONOMY_AUTO_FILE, "taxonomy_auto_overlay"),
                (integrated.TAXONOMY_REVIEWED_FILE, "taxonomy_reviewed_overlay"),
            )
        )

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "experiment_id": f"KH_ABLATION_{group}_CURRENT_SCORER_NO_HISTORY_V1",
        "status": "complete_frozen",
        "ablation_group": group,
        "factors": {
            "mngs": "new_multi_assay",
            "scorer": "current_deterministic_engine_v3_v1_v8",
            "taxonomy": taxonomy_mode,
            "history": "old_none",
        },
        "scorer_bundle_definition": {
            "analytical": summaries["analytical_scorer"].get("policy_id"),
            "clinical": summaries["clinical_scorer"].get("policy_id"),
            "promotion": summaries["promotion_v8"].get("policy_id"),
            "reporting": summaries["reporting_v3_route_a"].get("policy_id"),
        },
        "answer_blind_decision_stages": True,
        "answer_blind_freeze_sha256": integrated.sha256_file(
            output_root / "answer_blind_freeze.json"
        ),
        "components": components,
        "patient_count": len(patient_ids),
        "candidate_count": summaries["integrated_output"]["candidate_count"],
        "history_disabled": True,
        "ober_included": False,
        "evaluation_attached_posthoc": True,
        "limitations": [
            "This is the fixed 33-patient retrospective development benchmark.",
            "History v4C is deliberately disabled to isolate C/G.",
            "OBER/Luna are excluded from all decision stages.",
            "Complete-report metrics are supplementary to strict Picked and must not be called scorer-only metrics.",
        ],
    }
    write_json(output_root / "run_manifest.json", manifest)
    summary = {
        **manifest,
        "stage_summaries": summaries,
        "clinical_evaluation": clinical_evaluation,
        "reporting_evaluation": reporting_evaluation,
    }
    write_json(output_root / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", required=True, choices=("C", "G"))
    parser.add_argument(
        "--taxonomy-mode", required=True, choices=("old", "current_v2_6")
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--old-compact-root", type=Path, default=DEFAULT_OLD_COMPACT)
    args = parser.parse_args()
    print(
        json.dumps(
            run(
                group=args.group,
                taxonomy_mode=args.taxonomy_mode,
                output_root=args.output_root,
                old_compact_root=args.old_compact_root,
            ),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

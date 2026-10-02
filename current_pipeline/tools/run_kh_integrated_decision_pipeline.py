"""Run the frozen KH decision boundary through one auditable entry point.

This orchestrator deliberately starts from the frozen, history-integrated
clinical scorer output.  It does not recompute taxonomy, history extraction,
or the upstream analytical scorer.  It applies the cumulative v8 promotion
policy, builds the unified Picked/Possible/Fallback-Possible/Context report,
emits one compact patient result, and optionally runs post-hoc evaluation and
parity checks against previously frozen outputs.

Benchmark answers are never passed to either decision stage.  They are used
only after the answer-blind outputs have been written.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from tools import build_test_aware_possible_pathogen_shadow as reporting
from tools import evaluate_test_aware_clinical_deterministic_shadow as eval_clinical
from tools import evaluate_test_aware_possible_pathogen_shadow as eval_reporting
from tools import promote_test_aware_clinical_high_shadow_v8 as promotion_v8


DEFAULT_PROMOTION_POLICY = Path(
    "rules/test_aware_clinical_promotion_v6_analytically_compelling_rpm5.json"
)
DEFAULT_REPORTING_POLICY = Path("rules/test_aware_unified_reporting_v2.json")
SCHEMA_VERSION = "kh_integrated_decision_pipeline.v1"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().casefold() in {"1", "true", "yes"}


def _join(values: Iterable[Any]) -> str:
    return "|".join(str(value) for value in values if value not in (None, ""))


def _candidate_key(row: dict[str, Any], *, include_source: bool) -> tuple[str, ...]:
    key = (
        str(row.get("patient_id") or "").strip(),
        str(row.get("organism_name") or "").strip().casefold(),
    )
    if include_source:
        return key + (str(row.get("candidate_source") or "").strip().casefold(),)
    return key


def _index_unique(
    rows: list[dict[str, Any]], *, include_source: bool, label: str
) -> dict[tuple[str, ...], dict[str, Any]]:
    output: dict[tuple[str, ...], dict[str, Any]] = {}
    for row in rows:
        key = _candidate_key(row, include_source=include_source)
        if key in output:
            raise ValueError(f"Duplicate {label} candidate key: {key}")
        output[key] = row
    return output


def compare_candidate_views(
    current_rows: list[dict[str, Any]],
    reference_rows: list[dict[str, Any]],
    *,
    fields: tuple[str, ...],
    include_source: bool,
    stage: str,
) -> dict[str, Any]:
    """Compare candidate identity and semantic decision fields exactly."""
    current = _index_unique(
        current_rows, include_source=include_source, label=f"current {stage}"
    )
    reference = _index_unique(
        reference_rows, include_source=include_source, label=f"reference {stage}"
    )
    missing = sorted(set(reference) - set(current))
    extra = sorted(set(current) - set(reference))
    differences: list[dict[str, Any]] = []
    for key in sorted(set(current) & set(reference)):
        for field in fields:
            current_value = str(current[key].get(field) or "")
            reference_value = str(reference[key].get(field) or "")
            if current_value != reference_value:
                differences.append({
                    "stage": stage,
                    "candidate_key": " | ".join(key),
                    "field": field,
                    "reference_value": reference_value,
                    "current_value": current_value,
                })
    for key in missing:
        differences.append({
            "stage": stage,
            "candidate_key": " | ".join(key),
            "field": "__candidate__",
            "reference_value": "present",
            "current_value": "missing",
        })
    for key in extra:
        differences.append({
            "stage": stage,
            "candidate_key": " | ".join(key),
            "field": "__candidate__",
            "reference_value": "missing",
            "current_value": "present",
        })
    return {
        "stage": stage,
        "passed": not differences,
        "current_candidate_count": len(current),
        "reference_candidate_count": len(reference),
        "missing_candidate_count": len(missing),
        "extra_candidate_count": len(extra),
        "field_difference_count": sum(
            row["field"] != "__candidate__" for row in differences
        ),
        "differences": differences,
    }


def _compact_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    taxonomy = candidate.get("taxonomy_profile") or {}
    analytical = candidate.get("analytical_profile") or {}
    promotion = candidate.get("promotion_gate") or {}
    promotion_axes = promotion.get("axes") or {}
    possible = candidate.get("possible_pathogen_gate") or {}
    possible_axes = possible.get("axes") or {}
    history_route = candidate.get("history_adjusted_route")
    pre_promotion = candidate.get("pre_promotion_decision")
    final_tier = candidate.get("final_reporting_tier")
    return {
        "patient_id": str(candidate.get("patient_id") or ""),
        "organism_name": candidate.get("organism_name"),
        "organism_key": candidate.get("organism_key"),
        "candidate_source": candidate.get("candidate_source"),
        "candidate_event_refs": candidate.get("candidate_event_refs") or [],
        "decision_path": {
            "analytical_route": candidate.get("analytic_route"),
            "history_adjusted_route": history_route,
            "pre_promotion_clinical_decision": pre_promotion,
            "post_promotion_clinical_decision": candidate.get("clinical_decision"),
            "final_reporting_tier": final_tier,
        },
        "taxonomy": {
            "display_name": taxonomy.get("display_name"),
            "taxid": taxonomy.get("taxid"),
            "rank": taxonomy.get("taxonomic_rank"),
            "family": taxonomy.get("primary_rule_family"),
            "mapping_status": taxonomy.get("mapping_status"),
            "classification_confidence": taxonomy.get(
                "classification_confidence"
            ),
            "matched_rule_ids": taxonomy.get("matched_rule_ids") or [],
        },
        "analytical_evidence": {
            "best_rank": analytical.get("best_rank_in_retained_universe"),
            "selected_positive_test_count": analytical.get(
                "selected_positive_test_count"
            ),
            "reproducibility_axis": analytical.get("reproducibility_axis"),
            "cross_molecule_selected": analytical.get("cross_molecule_selected"),
            "per_test_signals": analytical.get("per_test_signals") or [],
        },
        "clinical_evidence": {
            "direct_hospital_level": promotion_axes.get("direct_hospital_level"),
            "event_aligned_direct_positive_count": (
                possible_axes.get("event_aligned_direct_positive_count")
            ),
            "host_support": promotion_axes.get("host_support"),
            "specimen_context": promotion_axes.get("specimen_context"),
            "precise_identity": promotion_axes.get("precise_identity"),
        },
        "promotion": {
            "outcome": promotion.get("outcome"),
            "route": promotion.get("route"),
            "blockers": promotion.get("blockers") or [],
        },
        "reporting": {
            "role": candidate.get("possible_reporting_role"),
            "selected_for_complete_report": _truthy(
                candidate.get("selected_for_complete_report")
            ),
            "tier": final_tier,
            "confidence": candidate.get("final_reporting_confidence"),
            "route": possible.get("route"),
            "outcome": possible.get("outcome"),
            "blockers": possible.get("blockers") or [],
            "cautions": possible.get("cautions") or [],
            "fallback_rank": possible_axes.get("fallback_rank"),
            "fallback_primary": possible_axes.get("fallback_primary"),
        },
        "audit": {
            "history_rule_ids": candidate.get("history_rule_ids") or [],
            "decision_rule_ids": candidate.get("rule_ids") or [],
            "reasons": candidate.get("reasons") or [],
            "history_adjustment_blockers": (
                candidate.get("history_adjustment_blockers") or []
            ),
            "consumed_history_evidence_ids": (
                candidate.get("consumed_evidence_ids") or []
            ),
        },
    }


def build_integrated_results(reporting_root: Path, output_root: Path) -> dict[str, Any]:
    """Emit one patient-facing result and one flat audit table."""
    source_paths = sorted(
        (reporting_root / "patient_outputs").glob(
            "NGS_patient_*_test_aware_possible_pathogen_shadow.json"
        ),
        key=lambda path: int(path.name.split("_")[2]),
    )
    if not source_paths:
        raise ValueError(f"No reporting patient outputs found under {reporting_root}")
    patient_dir = output_root / "patient_results"
    patient_dir.mkdir(parents=True, exist_ok=True)
    flat_rows: list[dict[str, Any]] = []
    tier_counts: Counter[str] = Counter()
    clinical_counts: Counter[str] = Counter()

    for source_path in source_paths:
        source = read_json(source_path)
        candidates = [
            _compact_candidate(candidate)
            for candidate in source.get("all_forwarded_candidates") or []
        ]
        for candidate in candidates:
            path = candidate["decision_path"]
            reporting_data = candidate["reporting"]
            taxonomy = candidate["taxonomy"]
            analytical = candidate["analytical_evidence"]
            clinical = candidate["clinical_evidence"]
            promotion = candidate["promotion"]
            audit = candidate["audit"]
            tier_counts[str(path.get("final_reporting_tier") or "missing")] += 1
            clinical_counts[str(
                path.get("post_promotion_clinical_decision") or "missing"
            )] += 1
            flat_rows.append({
                "patient_id": candidate["patient_id"],
                "organism_name": candidate["organism_name"],
                "candidate_source": candidate["candidate_source"],
                "taxonomy_family": taxonomy["family"],
                "taxonomy_mapping_status": taxonomy["mapping_status"],
                "analytical_route": path["analytical_route"],
                "history_adjusted_route": path["history_adjusted_route"],
                "pre_promotion_clinical_decision": path[
                    "pre_promotion_clinical_decision"
                ],
                "post_promotion_clinical_decision": path[
                    "post_promotion_clinical_decision"
                ],
                "promotion_outcome": promotion["outcome"],
                "promotion_route": promotion["route"],
                "final_reporting_tier": path["final_reporting_tier"],
                "reporting_role": reporting_data["role"],
                "selected_for_complete_report": reporting_data[
                    "selected_for_complete_report"
                ],
                "reporting_route": reporting_data["route"],
                "best_rank": analytical["best_rank"],
                "selected_positive_test_count": analytical[
                    "selected_positive_test_count"
                ],
                "reproducibility_axis": analytical["reproducibility_axis"],
                "cross_molecule_selected": analytical["cross_molecule_selected"],
                "specimen_context": clinical["specimen_context"],
                "direct_hospital_level": clinical["direct_hospital_level"],
                "history_rule_ids": _join(audit["history_rule_ids"]),
                "decision_rule_ids": _join(audit["decision_rule_ids"]),
                "promotion_blockers": _join(promotion["blockers"]),
                "reporting_blockers": _join(reporting_data["blockers"]),
                "reporting_cautions": _join(reporting_data["cautions"]),
                "reasons": _join(audit["reasons"]),
            })

        patient_id = str(source.get("patient_id") or "")
        write_json(
            patient_dir / f"patient_{patient_id}_integrated_decision.json",
            {
                "schema_version": SCHEMA_VERSION,
                "answer_blind": True,
                "patient_id": patient_id,
                "candidate_count": len(candidates),
                "final_report": [
                    item for item in candidates
                    if item["reporting"]["selected_for_complete_report"]
                ],
                "all_candidates": candidates,
                "source_reporting_file": str(source_path.resolve()),
                "source_reporting_sha256": sha256_file(source_path),
            },
        )

    csv_path = output_root / "integrated_decisions.csv"
    fields = list(flat_rows[0]) if flat_rows else []
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat_rows)
    return {
        "patient_count": len(source_paths),
        "candidate_count": len(flat_rows),
        "clinical_decision_counts": dict(sorted(clinical_counts.items())),
        "final_reporting_tier_counts": dict(sorted(tier_counts.items())),
        "integrated_decisions_file": str(csv_path.resolve()),
        "patient_results_root": str(patient_dir.resolve()),
    }


def _parity_report(
    decision_root: Path,
    reporting_root: Path,
    *,
    reference_promotion_root: Path | None,
    reference_reporting_root: Path | None,
) -> dict[str, Any]:
    if reference_promotion_root is None and reference_reporting_root is None:
        return {"status": "not_requested", "passed": None, "stages": []}
    if reference_promotion_root is None or reference_reporting_root is None:
        raise ValueError(
            "Both reference promotion and reporting roots are required for parity"
        )
    promotion_fields = (
        "incoming_decision",
        "pre_promotion_decision",
        "clinical_decision",
        "clinical_level",
        "promotion_outcome",
        "promotion_route",
        "promotion_blockers",
        "formal_pick_allowed",
        "selection_role",
        "selected_for_report",
        "rule_ids",
    )
    reporting_fields = (
        "clinical_decision",
        "original_selection_role",
        "possible_reporting_role",
        "selected_for_complete_report",
        "possible_outcome",
        "possible_route",
        "possible_blockers",
        "possible_cautions",
        "final_reporting_tier",
        "fallback_rank",
        "fallback_primary",
        "final_reporting_confidence",
    )
    stages = [
        compare_candidate_views(
            read_csv(decision_root / "clinical_decisions.csv"),
            read_csv(reference_promotion_root / "clinical_decisions.csv"),
            fields=promotion_fields,
            include_source=True,
            stage="promotion",
        ),
        compare_candidate_views(
            read_csv(reporting_root / "possible_pathogen_decisions.csv"),
            read_csv(reference_reporting_root / "possible_pathogen_decisions.csv"),
            fields=reporting_fields,
            include_source=False,
            stage="reporting",
        ),
    ]
    return {
        "status": "completed",
        "passed": all(stage["passed"] for stage in stages),
        "reference_promotion_root": str(reference_promotion_root.resolve()),
        "reference_reporting_root": str(reference_reporting_root.resolve()),
        "stages": stages,
    }


def _write_parity_artifacts(output_root: Path, parity: dict[str, Any]) -> None:
    write_json(output_root / "parity_report.json", parity)
    differences = [
        difference
        for stage in parity.get("stages") or []
        for difference in stage.get("differences") or []
    ]
    fields = [
        "stage", "candidate_key", "field", "reference_value", "current_value"
    ]
    with (output_root / "parity_differences.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(differences)
    lines = [
        "# Integrated pipeline parity check",
        "",
        f"- Status: `{parity.get('status')}`",
        f"- Passed: `{parity.get('passed')}`",
        f"- Difference rows: `{len(differences)}`",
        "",
    ]
    for stage in parity.get("stages") or []:
        lines.extend([
            f"## {stage['stage']}",
            "",
            f"- Passed: `{stage['passed']}`",
            f"- Current candidates: `{stage['current_candidate_count']}`",
            f"- Reference candidates: `{stage['reference_candidate_count']}`",
            f"- Missing / extra: `{stage['missing_candidate_count']}` / "
            f"`{stage['extra_candidate_count']}`",
            f"- Field differences: `{stage['field_difference_count']}`",
            "",
        ])
    (output_root / "parity_report.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def run(
    input_root: Path,
    output_root: Path,
    *,
    promotion_policy: Path = DEFAULT_PROMOTION_POLICY,
    reporting_policy: Path = DEFAULT_REPORTING_POLICY,
    patients: set[int] | None = None,
    rpm_threshold: float | None = None,
    qc_mode: str | None = None,
    answers: Path | None = None,
    reference_promotion_root: Path | None = None,
    reference_reporting_root: Path | None = None,
    require_parity: bool = False,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    decision_root = output_root / "decision_engine"
    reporting_root = output_root / "reporting"

    decision_summary = promotion_v8.run(
        input_root,
        decision_root,
        policy_path=promotion_policy,
        patients=patients,
        rpm_threshold=rpm_threshold,
        qc_mode=qc_mode,
    )
    reporting_summary = reporting.run(
        decision_root,
        reporting_root,
        policy_path=reporting_policy,
        patients=patients,
    )
    integrated_summary = build_integrated_results(reporting_root, output_root)

    evaluation: dict[str, Any] | None = None
    if answers is not None:
        evaluation_root = output_root / "evaluation"
        clinical_metrics = eval_clinical.run(
            decision_root, answers, evaluation_root / "clinical"
        )
        reporting_metrics = eval_reporting.run(
            reporting_root, answers, evaluation_root / "reporting"
        )
        evaluation = {
            "answer_path": str(answers.resolve()),
            "answer_sha256": sha256_file(answers),
            "separation": (
                "Post-hoc only: answers were not supplied to promotion or reporting."
            ),
            "clinical_metrics": clinical_metrics["metrics"],
            "reporting_metrics": reporting_metrics["metrics"],
        }

    parity = _parity_report(
        decision_root,
        reporting_root,
        reference_promotion_root=reference_promotion_root,
        reference_reporting_root=reference_reporting_root,
    )
    _write_parity_artifacts(output_root, parity)

    input_summary_path = input_root / "summary.json"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "answer_blind_decision_stages": True,
        "execution_boundary": (
            "Starts from frozen history-integrated clinical scorer output; "
            "upstream taxonomy/history/scorer decisions are reused, not recomputed."
        ),
        "input": {
            "root": str(input_root.resolve()),
            "summary_sha256": (
                sha256_file(input_summary_path)
                if input_summary_path.is_file() else None
            ),
        },
        "policies": {
            "promotion": {
                "path": str(promotion_policy.resolve()),
                "sha256": sha256_file(promotion_policy),
            },
            "reporting": {
                "path": str(reporting_policy.resolve()),
                "sha256": sha256_file(reporting_policy),
            },
        },
        "implementation": {
            "orchestrator": str(Path(__file__).resolve()),
            "orchestrator_sha256": sha256_file(Path(__file__)),
            "promotion_module": str(Path(promotion_v8.__file__).resolve()),
            "promotion_module_sha256": sha256_file(Path(promotion_v8.__file__)),
            "reporting_module": str(Path(reporting.__file__).resolve()),
            "reporting_module_sha256": sha256_file(Path(reporting.__file__)),
        },
        "parameters": {
            "patients": sorted(patients) if patients is not None else None,
            "rpm_threshold": rpm_threshold,
            "qc_mode": qc_mode,
        },
        "outputs": {
            "decision_root": str(decision_root.resolve()),
            "reporting_root": str(reporting_root.resolve()),
            "decision_summary_sha256": sha256_file(decision_root / "summary.json"),
            "reporting_summary_sha256": sha256_file(reporting_root / "summary.json"),
        },
        "parity_passed": parity.get("passed"),
        "evaluation_attached_post_hoc": evaluation is not None,
    }
    write_json(output_root / "run_manifest.json", manifest)

    summary = {
        "schema_version": SCHEMA_VERSION,
        "answer_blind": True,
        "execution_boundary": manifest["execution_boundary"],
        "decision_engine": decision_summary,
        "reporting": reporting_summary,
        "integrated_output": integrated_summary,
        "parity": {
            "status": parity["status"],
            "passed": parity["passed"],
            "difference_count": sum(
                len(stage.get("differences") or [])
                for stage in parity.get("stages") or []
            ),
        },
        "evaluation": evaluation,
    }
    write_json(output_root / "summary.json", summary)
    if require_parity and not parity.get("passed"):
        raise RuntimeError(
            "Integrated output failed parity; inspect parity_report.json"
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument(
        "--promotion-policy", type=Path, default=DEFAULT_PROMOTION_POLICY
    )
    parser.add_argument(
        "--reporting-policy", type=Path, default=DEFAULT_REPORTING_POLICY
    )
    parser.add_argument("--patients", type=int, nargs="*")
    parser.add_argument("--rpm-threshold", type=float)
    parser.add_argument(
        "--qc-mode", choices=("partial_or_full", "full_only")
    )
    parser.add_argument("--answers", type=Path)
    parser.add_argument("--reference-promotion-root", type=Path)
    parser.add_argument("--reference-reporting-root", type=Path)
    parser.add_argument("--require-parity", action="store_true")
    args = parser.parse_args()
    selected = set(args.patients) if args.patients else None
    print(json.dumps(
        run(
            args.input_root,
            args.output_root,
            promotion_policy=args.promotion_policy,
            reporting_policy=args.reporting_policy,
            patients=selected,
            rpm_threshold=args.rpm_threshold,
            qc_mode=args.qc_mode,
            answers=args.answers,
            reference_promotion_root=args.reference_promotion_root,
            reference_reporting_root=args.reference_reporting_root,
            require_parity=args.require_parity,
        ),
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()

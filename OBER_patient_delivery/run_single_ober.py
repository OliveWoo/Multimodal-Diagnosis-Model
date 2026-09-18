from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Iterable, Mapping

from generate_r5_final_decisions import (
    RULE_FORMULA,
    RULE_ID,
    FinalizerError,
    candidate_name,
    compose_decision,
    load_casefit_signals,
    load_json,
    normal,
    runtime_modules,
    sha256_file,
    write_outputs,
)


def _rows(value: Any) -> list[Mapping[str, Any]]:
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, Mapping)]


def _unique_names(rows: Iterable[Mapping[str, Any]]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for row in rows:
        name = candidate_name(row).strip()
        key = normal(name)
        if name and key not in seen:
            result.append(name)
            seen.add(key)
    return result


def _picked(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    deterministic = payload.get("deterministic_max") or {}
    summary = deterministic.get("best_available_summary") or {}
    rows = _rows(summary.get("picked_pathogens"))
    result: list[dict[str, Any]] = []
    for rank, row in enumerate(rows, start=1):
        name = candidate_name(row).strip()
        if not name:
            raise FinalizerError(f"picked row {rank} has no organism name")
        result.append(
            {
                "rank": rank,
                "canonical_name": name,
                "classification": row.get("classification"),
                "picked_role": row.get("picked_role"),
                "evidence_source": row.get("evidence_source"),
            }
        )
    expected = summary.get("picked_count")
    if type(expected) is int and expected != len(result):
        raise FinalizerError(
            f"picked_count mismatch: summary={expected}, rows={len(result)}"
        )
    return result


def _candidate_names(payload: Mapping[str, Any]) -> list[str]:
    deterministic = payload.get("deterministic_max") or {}
    review = payload.get("llm_missed_candidate_review") or {}
    combined: list[Mapping[str, Any]] = []
    combined.extend(_rows((deterministic.get("best_available_summary") or {}).get("picked_pathogens")))
    combined.extend(_rows(deterministic.get("pathogen_candidates")))
    for tier in (
        "review_high_priority",
        "review_context_needed",
        "review_low_specificity",
        "review_omitted_with_reason",
    ):
        combined.extend(_rows(review.get(tier)))
    return _unique_names(combined)


def run_single(
    *,
    workspace: Path,
    merge_file: Path,
    casefit_root: Path,
    clinical_config: Path,
    output: Path,
    cohort: str,
) -> dict[str, Any]:
    workspace = workspace.resolve()
    merge_file = merge_file.resolve()
    casefit_root = casefit_root.resolve()
    clinical_config = clinical_config.resolve()
    payload = load_json(merge_file)
    patient_id = str(payload.get("patient_id") or "").strip()
    if not patient_id:
        raise FinalizerError("merge file has no patient_id")

    aliases, evaluate_modules = runtime_modules(workspace)
    signals_by_patient, signal_meta = load_casefit_signals(
        casefit_root,
        clinical_config,
        aliases,
        evaluate_modules,
    )
    unexpected = sorted(set(signals_by_patient) - {patient_id})
    if unexpected:
        raise FinalizerError(f"casefit root contains other patients: {unexpected}")
    signals = signals_by_patient.get(patient_id, [])
    if not signals:
        raise FinalizerError(f"no casefit signals found for P{patient_id}")
    for signal in signals:
        if signal.get("latest_merge_sha256") != sha256_file(merge_file):
            raise FinalizerError(
                f"casefit input does not match supplied merge file for P{patient_id}"
            )

    decision = compose_decision(
        cohort=cohort,
        patient_id=patient_id,
        candidate_names=_candidate_names(payload),
        picked=_picked(payload),
        signals=signals,
        source_snapshot={
            "latest_merge_file": merge_file.name,
            "latest_merge_path": str(merge_file),
            "latest_merge_sha256": sha256_file(merge_file),
            "casefit_manifest_file": str((casefit_root / "manifest.json").resolve()),
            "casefit_manifest_sha256": sha256_file(casefit_root / "manifest.json"),
            "clinical_config_file": clinical_config.name,
            "clinical_config_sha256": sha256_file(clinical_config),
            "casefit": signal_meta,
        },
    )
    manifest = write_outputs(
        output.resolve(),
        [decision],
        {
            "mode": "single_patient_rerun",
            "rule_id": RULE_ID,
            "formula": RULE_FORMULA,
            "merge_file": str(merge_file),
            "merge_sha256": sha256_file(merge_file),
            "casefit_root": str(casefit_root),
            "casefit_manifest_sha256": sha256_file(casefit_root / "manifest.json"),
            "clinical_config": str(clinical_config),
            "clinical_config_sha256": sha256_file(clinical_config),
        },
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create an OBER R5 final decision for one merged patient JSON"
    )
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--merge-file", required=True, type=Path)
    parser.add_argument("--casefit-root", required=True, type=Path)
    parser.add_argument("--clinical-config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cohort", default="KH")
    args = parser.parse_args()
    manifest = run_single(
        workspace=args.workspace,
        merge_file=args.merge_file,
        casefit_root=args.casefit_root,
        clinical_config=args.clinical_config,
        output=args.output,
        cohort=args.cohort,
    )
    print(
        {
            "status": manifest["status"],
            "totals": manifest["totals"],
            "output": str(args.output.resolve()),
        }
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Reconstruct answer-blind run provenance from frozen merged artifacts.

This tool records what the saved JSON files prove. It intentionally does not
claim to recover the original command, environment, or historical Git commit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST_SCHEMA_VERSION = "reconstructed_pipeline_manifest_v1.0"


def load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected JSON object: {path}")
    return payload


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def patient_number(path: Path) -> int | None:
    match = re.search(r"NGS_patient_(\d+)_json", str(path))
    return int(match.group(1)) if match else None


def resolve_recorded_path(value: Any, *, patient_dir: Path) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    raw = Path(text)
    candidates = [
        raw,
        REPO_ROOT / raw,
        patient_dir / raw,
        patient_dir / "summary_outputs" / raw,
    ]
    for candidate in candidates:
        if candidate.exists() and candidate.is_file():
            return candidate.resolve()
    return None


def artifact(path: Path | None, recorded_path: Any = None) -> dict[str, Any]:
    if path is None:
        return {
            "recorded_path": str(recorded_path or "") or None,
            "resolved": False,
            "path": None,
            "sha256": None,
            "bytes": None,
        }
    return {
        "recorded_path": str(recorded_path or path),
        "resolved": True,
        "path": str(path),
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


KNOWN_RANKING_FORMULAS = (
    (0.60, 0.40, "final_score = 0.60 * code_score + 0.40 * reads_percentile"),
    (0.70, 0.30, "final_score = 0.70 * code_score + 0.30 * reads_percentile"),
)


def iter_ranking_scores(value: Any) -> Iterable[tuple[float, float, float]]:
    if isinstance(value, dict):
        if {"code_score", "reads_percentile", "final_score"}.issubset(value):
            try:
                yield (
                    float(value["code_score"]),
                    float(value["reads_percentile"]),
                    float(value["final_score"]),
                )
            except (TypeError, ValueError):
                pass
        for child in value.values():
            yield from iter_ranking_scores(child)
    elif isinstance(value, list):
        for child in value:
            yield from iter_ranking_scores(child)


def ranking_formula_observation(payload: dict[str, Any]) -> dict[str, Any]:
    metadata = payload.get("ranking_metadata") or payload.get("metadata") or {}
    if isinstance(metadata, dict):
        formula = metadata.get("score_formula")
        if formula:
            return {
                "formula": str(formula),
                "evidence_source": "metadata.score_formula",
                "informative_candidate_count": 0,
            }
        for note in metadata.get("notes") or []:
            text = str(note)
            if "final_score" in text and any(weight in text for weight in ("0.60", "0.40", "0.7", "0.3")):
                return {
                    "formula": text,
                    "evidence_source": "metadata.notes",
                    "informative_candidate_count": 0,
                }

    scores = [row for row in iter_ranking_scores(payload) if abs(row[0] - row[1]) > 1e-9]
    matches: list[tuple[str, int]] = []
    for code_weight, reads_weight, label in KNOWN_RANKING_FORMULAS:
        matching = sum(
            abs(final - (code_weight * code + reads_weight * percentile)) <= 0.0005
            for code, percentile, final in scores
        )
        if scores and matching == len(scores):
            matches.append((label, matching))
    if len(matches) == 1 and matches[0][1] >= 3:
        return {
            "formula": matches[0][0],
            "evidence_source": "inferred_from_stored_candidate_scores",
            "informative_candidate_count": matches[0][1],
        }
    return {
        "formula": None,
        "evidence_source": "not_recoverable",
        "informative_candidate_count": len(scores),
    }


def ranking_formula(payload: dict[str, Any]) -> str | None:
    return ranking_formula_observation(payload)["formula"]


def review_settings(payload: dict[str, Any]) -> dict[str, Any]:
    request = payload.get("llm_request") or {}
    if not isinstance(request, dict):
        request = {}
    return {
        "model": request.get("model"),
        "reasoning_effort": request.get("reasoning_effort"),
        "temperature_requested": request.get("temperature_requested"),
        "temperature_sent": request.get("temperature_sent"),
    }


def merge_semantics(merged: dict[str, Any], review: dict[str, Any]) -> dict[str, Any]:
    syndrome = review.get("syndrome_pattern_convergence_policy") or {}
    convergence = review.get("review_tier_convergence_policy") or {}
    nonpulmonary = review.get("nonpulmonary_hospital_evidence_guardrail") or {}
    reconciliation = merged.get("final_name_reconciliation") or {}
    return {
        "deterministic_picked_advisory_override": False,
        "name_reconciliation_enabled": bool(reconciliation.get("enabled")),
        "nonpulmonary_source_guardrail_enabled": bool(nonpulmonary.get("enabled")),
        "syndrome_convergence_version": syndrome.get("version"),
        "syndrome_promotions_enabled": syndrome.get("promotions_enabled"),
        "review_tier_convergence_version": convergence.get("version"),
        "review_tier_promotions_enabled": convergence.get("promotions_enabled"),
        "rag_visible_tiers": (review.get("review_tiering_policy") or {}).get("rag_visible_tiers"),
    }


def collect_case(merged_path: Path) -> dict[str, Any]:
    merged = load_json(merged_path)
    patient_dir = merged_path.parent.parent
    patient_id = str(merged.get("patient_id") or patient_number(merged_path) or "")
    merged_sources = merged.get("merged_source_files") or {}
    if not isinstance(merged_sources, dict):
        merged_sources = {}

    deterministic_recorded = merged_sources.get("deterministic_max")
    deterministic_path = resolve_recorded_path(deterministic_recorded, patient_dir=patient_dir)
    deterministic = load_json(deterministic_path) if deterministic_path else (merged.get("deterministic_max") or {})

    deterministic_sources = deterministic.get("source_files") or {}
    if not isinstance(deterministic_sources, dict):
        deterministic_sources = {}
    summary_recorded = deterministic_sources.get("final_summary")
    ranked_recorded = deterministic_sources.get("ranked_mngs")
    summary_path = resolve_recorded_path(summary_recorded, patient_dir=patient_dir)
    ranked_path = resolve_recorded_path(ranked_recorded, patient_dir=patient_dir)
    summary = load_json(summary_path) if summary_path else {}
    ranked = load_json(ranked_path) if ranked_path else {}

    review_recorded = merged_sources.get("llm_missed_candidate_review")
    review_path = resolve_recorded_path(review_recorded, patient_dir=patient_dir)
    review = load_json(review_path) if review_path else (merged.get("llm_missed_candidate_review") or {})
    if not isinstance(review, dict):
        review = {}

    ranking_observation = ranking_formula_observation(ranked)
    return {
        "patient_id": patient_id,
        "artifacts": {
            "merged": artifact(merged_path.resolve()),
            "deterministic": artifact(deterministic_path, deterministic_recorded),
            "summary": artifact(summary_path, summary_recorded),
            "ranked_mngs": artifact(ranked_path, ranked_recorded),
            "review": artifact(review_path, review_recorded),
        },
        "observed_configuration": {
            "scorer_rule_version": deterministic.get("rule_version"),
            "summary_rule_version": summary.get("rule_version"),
            "ranking_formula": ranking_observation["formula"],
            "ranking_formula_evidence": {
                key: value for key, value in ranking_observation.items() if key != "formula"
            },
            "review": review_settings(review),
            "merge_semantics": merge_semantics(merged, review),
        },
    }


def value_counts(cases: Iterable[dict[str, Any]], path: tuple[str, ...]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for case in cases:
        value: Any = case
        for key in path:
            value = value.get(key) if isinstance(value, dict) else None
        counts[str(value) if value is not None else "<missing>"] += 1
    return dict(sorted(counts.items()))


def build_manifest(
    patient_root: Path,
    merged_suffix: str,
    *,
    experiment_id: str,
    cohort: str,
    references: list[Path],
    patient_ids: set[str] | None = None,
) -> dict[str, Any]:
    merged_files = sorted(
        patient_root.glob(f"NGS_patient_*_json/summary_outputs/*_{merged_suffix}.json"),
        key=lambda path: (patient_number(path) is None, patient_number(path) or 0, str(path)),
    )
    if not merged_files:
        raise ValueError(f"No merged files found for suffix {merged_suffix!r} under {patient_root}")
    if patient_ids:
        merged_files = [path for path in merged_files if str(patient_number(path)) in patient_ids]
        found_ids = {str(patient_number(path)) for path in merged_files}
        missing_ids = sorted(patient_ids - found_ids, key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value))
        if missing_ids:
            raise ValueError(f"Requested patient IDs have no matching merged artifact: {missing_ids}")
    cases = [collect_case(path) for path in merged_files]
    unresolved = [
        {"patient_id": case["patient_id"], "artifact": name, "recorded_path": item["recorded_path"]}
        for case in cases
        for name, item in case["artifacts"].items()
        if name != "merged" and not item["resolved"]
    ]
    reference_records = [artifact(path.resolve()) for path in references]
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "experiment_id": experiment_id,
        "cohort": cohort,
        "manifest_status": "reconstructed_from_frozen_artifacts",
        "historical_claim_limit": (
            "This manifest proves the contents and embedded settings of currently preserved artifacts. "
            "It does not recover the original command, environment, model service state, or historical Git commit."
        ),
        "patient_root": str(patient_root.resolve()),
        "merged_suffix": merged_suffix,
        "patient_selection": sorted(
            (case["patient_id"] for case in cases),
            key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value),
        ),
        "case_count": len(cases),
        "aggregate": {
            "scorer_rule_versions": value_counts(cases, ("observed_configuration", "scorer_rule_version")),
            "summary_rule_versions": value_counts(cases, ("observed_configuration", "summary_rule_version")),
            "ranking_formulas": value_counts(cases, ("observed_configuration", "ranking_formula")),
            "review_models": value_counts(cases, ("observed_configuration", "review", "model")),
            "review_reasoning_efforts": value_counts(
                cases, ("observed_configuration", "review", "reasoning_effort")
            ),
            "unresolved_artifact_count": len(unresolved),
        },
        "metric_and_report_references": reference_records,
        "unresolved_artifacts": unresolved,
        "cases": cases,
    }


def markdown_summary(manifest: dict[str, Any]) -> str:
    aggregate = manifest["aggregate"]
    lines = [
        f"# {manifest['experiment_id']} reconstructed provenance",
        "",
        f"- Cohort: `{manifest['cohort']}`",
        f"- Cases: {manifest['case_count']}",
        f"- Merged suffix: `{manifest['merged_suffix']}`",
        f"- Status: `{manifest['manifest_status']}`",
        "",
        "## Observed configuration",
        "",
        f"- Scorer versions: `{json.dumps(aggregate['scorer_rule_versions'], ensure_ascii=False)}`",
        f"- Summary versions: `{json.dumps(aggregate['summary_rule_versions'], ensure_ascii=False)}`",
        f"- Ranking formulas: `{json.dumps(aggregate['ranking_formulas'], ensure_ascii=False)}`",
        f"- Review models: `{json.dumps(aggregate['review_models'], ensure_ascii=False)}`",
        f"- Review reasoning: `{json.dumps(aggregate['review_reasoning_efforts'], ensure_ascii=False)}`",
        f"- Unresolved recorded artifacts: {aggregate['unresolved_artifact_count']}",
        "",
        "## Interpretation limit",
        "",
        manifest["historical_claim_limit"],
        "",
        "The JSON manifest contains the per-patient artifact paths and SHA256 values.",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--merged-suffix", required=True)
    parser.add_argument("--experiment-id", required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--reference", type=Path, action="append", default=[])
    parser.add_argument(
        "--patient-id",
        action="append",
        default=[],
        help="Restrict the manifest to one patient ID; repeat for multiple IDs.",
    )
    args = parser.parse_args()

    missing_references = [path for path in args.reference if not path.exists()]
    if missing_references:
        raise SystemExit(f"Missing reference files: {missing_references}")
    manifest = build_manifest(
        args.patient_root,
        args.merged_suffix,
        experiment_id=args.experiment_id,
        cohort=args.cohort,
        references=args.reference,
        patient_ids={str(value).strip() for value in args.patient_id if str(value).strip()} or None,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(markdown_summary(manifest), encoding="utf-8")
    print(json.dumps({key: value for key, value in manifest.items() if key != "cases"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

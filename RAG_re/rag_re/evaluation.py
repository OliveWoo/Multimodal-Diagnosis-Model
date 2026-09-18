from __future__ import annotations

import math
import random
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

from .input_adapter import normalize_name
from .io_utils import load_json, sha256_file
from .versioning import EVALUATOR_FINGERPRINT, EVALUATOR_VERSION


BUILTIN_ALIASES = {
    "hsv": "hsv-1",
    "herpes simplex virus 1": "hsv-1",
    "human herpesvirus 1": "hsv-1",
    "pjp": "pneumocystis jirovecii",
    "pneumocystis jiroveci": "pneumocystis jirovecii",
    "candida albican": "candida albicans",
    "crkp": "klebsiella pneumoniae",
}


def _canonical(value: Any, aliases: Mapping[str, str]) -> str:
    normalized = normalize_name(value)
    return normalize_name(aliases.get(normalized, normalized))


def load_gold(path: Path) -> Dict[str, Any]:
    payload = load_json(path)
    if not isinstance(payload, dict):
        raise ValueError("Gold label root must be an object")

    cases: List[Dict[str, Any]] = []
    if isinstance(payload.get("cases"), list):
        for index, row in enumerate(payload["cases"]):
            if not isinstance(row, dict):
                raise ValueError(f"Gold case at index {index} must be an object")
            if row.get("patient_id") in (None, ""):
                raise ValueError(f"Gold case at index {index} is missing patient_id")
            pathogens = row.get("pathogens")
            if not isinstance(pathogens, list):
                raise ValueError(f"Gold case {row.get('patient_id')} pathogens must be a list")
            uncertain = row.get("uncertain_pathogens", [])
            if not isinstance(uncertain, list):
                raise ValueError(
                    f"Gold case {row.get('patient_id')} uncertain_pathogens must be a list"
                )
            positive_values = [str(value) for value in pathogens if str(value).strip()]
            uncertain_values = [str(value) for value in uncertain if str(value).strip()]
            overlap = {
                normalize_name(value) for value in positive_values
            } & {normalize_name(value) for value in uncertain_values}
            if overlap:
                raise ValueError(
                    f"Gold case {row.get('patient_id')} has positive/uncertain overlap: "
                    f"{sorted(overlap)}"
                )
            cases.append(
                {
                    "patient_id": str(row["patient_id"]),
                    "pathogens": positive_values,
                    "uncertain_pathogens": uncertain_values,
                    "split": row.get("split"),
                }
            )
    elif isinstance(payload.get("patients"), dict):
        for patient_id, pathogens in payload["patients"].items():
            if not isinstance(pathogens, list):
                raise ValueError(f"Gold patient {patient_id} value must be a list")
            cases.append(
                {
                    "patient_id": str(patient_id),
                    "pathogens": [str(value) for value in pathogens if str(value).strip()],
                    "uncertain_pathogens": [],
                    "split": None,
                }
            )
    else:
        raise ValueError("Gold file needs either a cases list or patients object")

    ids = [row["patient_id"] for row in cases]
    if not cases:
        raise ValueError("Gold file contains no labeled cases")
    if len(ids) != len(set(ids)):
        raise ValueError("Gold file contains duplicate patient_id values")
    aliases = dict(BUILTIN_ALIASES)
    for alias, canonical in (payload.get("aliases") or {}).items():
        aliases[normalize_name(alias)] = normalize_name(canonical)
    return {
        "schema_version": payload.get("schema_version") or "rag_re.gold.v1",
        "label_semantics": payload.get("label_semantics") or "infection_source_positive_list",
        "complete_candidate_negative_labels": bool(
            payload.get("complete_candidate_negative_labels", False)
        ),
        "development_only": bool(payload.get("development_only", False)),
        "cases": cases,
        "aliases": aliases,
        "source_sha256": sha256_file(path),
    }


def load_outputs(path: Path) -> List[Dict[str, Any]]:
    paths = [path] if path.is_file() else sorted(path.rglob("*.json"))
    outputs: List[Dict[str, Any]] = []
    for candidate in paths:
        try:
            payload = load_json(candidate)
        except (ValueError, OSError):
            continue
        if isinstance(payload, dict) and payload.get("schema_version") == "rag_re.output.v1":
            payload = dict(payload)
            payload["_source_file"] = str(candidate)
            outputs.append(payload)
    by_id: Dict[str, Dict[str, Any]] = {}
    for payload in outputs:
        patient_id = str(payload.get("patient_id") or "")
        if patient_id in by_id:
            raise ValueError(f"Duplicate RAG_re output for patient_id={patient_id}")
        by_id[patient_id] = payload
    return list(by_id.values())


def _safe_div(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def _fbeta(precision: float | None, recall: float | None, beta: float) -> float | None:
    if precision is None or recall is None:
        return None
    if precision == 0 and recall == 0:
        return 0.0
    beta_sq = beta * beta
    return (1 + beta_sq) * precision * recall / (beta_sq * precision + recall)


def _percentile(values: Sequence[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def _aggregate(case_rows: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    tp = sum(row["tp"] for row in case_rows)
    fp = sum(row["fp"] for row in case_rows)
    fn = sum(row["fn"] for row in case_rows)
    conditional_fn = sum(row["conditional_fn"] for row in case_rows if row["pool_observed"])
    conditional_tp = sum(row["tp"] for row in case_rows if row["pool_observed"])
    truth_total = sum(row["truth_count"] for row in case_rows)
    pool_truth_total = sum(row["pool_truth_count"] for row in case_rows)
    truth_in_pool = sum(
        row["truth_in_pool_count"] for row in case_rows if row["pool_observed"]
    )
    candidate_total = sum(row["candidate_count"] for row in case_rows)
    precision = _safe_div(tp, tp + fp)
    recall = _safe_div(tp, tp + fn)
    conditional_recall = _safe_div(conditional_tp, conditional_tp + conditional_fn)
    positive_cases = [row for row in case_rows if row["truth_count"] > 0]
    negative_cases = [row for row in case_rows if row["truth_count"] == 0]
    macro_recalls = [row["case_recall"] for row in positive_cases if row["case_recall"] is not None]
    return {
        "case_count": len(case_rows),
        "positive_case_count": len(positive_cases),
        "negative_case_count": len(negative_cases),
        "truth_pathogen_count": truth_total,
        "end_to_end_truth_denominator": truth_total,
        "candidate_pool_truth_denominator": pool_truth_total,
        "conditional_truth_in_pool_denominator": conditional_tp + conditional_fn,
        "output_case_coverage": _safe_div(
            sum(row["pool_observed"] for row in case_rows), len(case_rows)
        ),
        "candidate_count": candidate_total,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "agreement_precision": precision,
        "end_to_end_recall": recall,
        "conditional_classifier_recall": conditional_recall,
        "candidate_pool_recall_ceiling": _safe_div(truth_in_pool, pool_truth_total),
        "f1": _fbeta(precision, recall, 1.0),
        "f2": _fbeta(precision, recall, 2.0),
        "macro_positive_case_recall": (
            sum(macro_recalls) / len(macro_recalls) if macro_recalls else None
        ),
        "positive_case_hit_rate": _safe_div(
            sum(row["tp"] > 0 for row in positive_cases), len(positive_cases)
        ),
        "negative_case_specificity": _safe_div(
            sum(row["predicted_count"] == 0 for row in negative_cases), len(negative_cases)
        ),
        "case_exact_match": _safe_div(
            sum(row["exact_match"] for row in case_rows), len(case_rows)
        ),
        "false_positives_per_case": _safe_div(fp, len(case_rows)),
        "abstained_candidates": sum(row["abstained_count"] for row in case_rows),
        "cases_with_abstention": sum(row["abstained_count"] > 0 for row in case_rows),
        "baseline_miss_rescued": sum(row["baseline_miss_rescued"] for row in case_rows),
        "baseline_false_positives_added": sum(
            row["baseline_false_positives_added"] for row in case_rows
        ),
    }


def _bootstrap(
    case_rows: Sequence[Dict[str, Any]], iterations: int, seed: int
) -> Dict[str, Any]:
    if not case_rows or iterations <= 0:
        return {}
    randomizer = random.Random(seed)
    values: Dict[str, List[float]] = {
        "agreement_precision": [],
        "end_to_end_recall": [],
        "f2": [],
        "negative_case_specificity": [],
    }
    for _ in range(iterations):
        sample = [randomizer.choice(case_rows) for _ in range(len(case_rows))]
        metrics = _aggregate(sample)
        for key in values:
            value = metrics.get(key)
            if value is not None:
                values[key].append(float(value))
    return {
        key: {
            "lower_95": _percentile(samples, 0.025),
            "upper_95": _percentile(samples, 0.975),
        }
        for key, samples in values.items()
    }


def _bootstrap_delta(
    before_rows: Sequence[Dict[str, Any]],
    after_rows: Sequence[Dict[str, Any]],
    iterations: int,
    seed: int,
) -> Dict[str, Any]:
    if not before_rows or iterations <= 0:
        return {}
    if [row["patient_id"] for row in before_rows] != [
        row["patient_id"] for row in after_rows
    ]:
        raise ValueError("Paired bootstrap requires the same ordered patients")
    randomizer = random.Random(seed)
    keys = (
        "agreement_precision",
        "end_to_end_recall",
        "f2",
        "false_positives_per_case",
    )
    values: Dict[str, List[float]] = {key: [] for key in keys}
    count = len(before_rows)
    for _ in range(iterations):
        indices = [randomizer.randrange(count) for _ in range(count)]
        before = _aggregate([before_rows[index] for index in indices])
        after = _aggregate([after_rows[index] for index in indices])
        for key in keys:
            left, right = before.get(key), after.get(key)
            if left is not None and right is not None:
                values[key].append(float(right) - float(left))
    return {
        key: {
            "lower_95": _percentile(samples, 0.025),
            "upper_95": _percentile(samples, 0.975),
        }
        for key, samples in values.items()
    }
def _sets_for_output(
    payload: Dict[str, Any], experiment: str, aliases: Mapping[str, str]
) -> Tuple[set[str], set[str], set[str]]:
    candidate_names = {
        _canonical(row.get("organism_name"), aliases)
        for row in payload.get("candidates", [])
        if isinstance(row, dict) and _canonical(row.get("organism_name"), aliases)
    }
    experiments = payload.get("experiments") or {}
    if experiment not in experiments:
        raise ValueError(
            f"Patient {payload.get('patient_id')} output is missing experiment {experiment}"
        )
    experiment_data = experiments[experiment] or {}
    predicted = {
        _canonical(value, aliases)
        for value in experiment_data.get("predicted_pathogens", [])
        if _canonical(value, aliases)
    }
    abstained = {
        _canonical(value, aliases)
        for value in experiment_data.get("abstained_pathogens", [])
        if _canonical(value, aliases)
    }
    outside = (predicted | abstained) - candidate_names
    if outside:
        raise ValueError(
            f"Patient {payload.get('patient_id')} experiment {experiment} contains "
            f"pathogens outside its frozen candidate pool: {sorted(outside)}"
        )
    overlap = predicted & abstained
    if overlap:
        raise ValueError(
            f"Patient {payload.get('patient_id')} experiment {experiment} marks the same "
            f"pathogens positive and abstained: {sorted(overlap)}"
        )
    return candidate_names, predicted, abstained


def _case_rows(
    outputs: Mapping[str, Dict[str, Any]],
    gold_cases: Iterable[Dict[str, Any]],
    experiment: str,
    aliases: Mapping[str, str],
) -> Tuple[List[Dict[str, Any]], List[str]]:
    rows: List[Dict[str, Any]] = []
    missing: List[str] = []
    for gold in gold_cases:
        patient_id = str(gold["patient_id"])
        payload = outputs.get(patient_id)
        pool_observed = payload is not None
        if payload is None:
            missing.append(patient_id)
            payload = {
                "patient_id": patient_id,
                "candidates": [],
                "experiments": {
                    experiment: {
                        "predicted_pathogens": [],
                        "abstained_pathogens": [],
                    },
                    "E0_BASELINE": {"predicted_pathogens": []},
                },
            }
        truth = {_canonical(value, aliases) for value in gold["pathogens"]}
        uncertain = {
            _canonical(value, aliases) for value in gold.get("uncertain_pathogens", [])
        }
        canonical_overlap = truth & uncertain
        if canonical_overlap:
            raise ValueError(
                f"Gold case {patient_id} has positive/uncertain overlap after aliasing: "
                f"{sorted(canonical_overlap)}"
            )
        candidates, predicted, abstained = _sets_for_output(payload, experiment, aliases)
        uncertain_predictions = predicted & uncertain
        candidates -= uncertain
        predicted -= uncertain
        abstained -= uncertain
        baseline_data = (payload.get("experiments") or {}).get("E0_BASELINE") or {}
        baseline = {
            _canonical(value, aliases)
            for value in baseline_data.get("predicted_pathogens", [])
        }
        baseline -= uncertain
        tp_set = truth & predicted
        fp_set = predicted - truth
        fn_set = truth - predicted
        truth_in_pool = truth & candidates
        conditional_fn = truth_in_pool - predicted
        case_recall = _safe_div(len(tp_set), len(truth))
        rows.append(
            {
                "patient_id": patient_id,
                "truth": sorted(truth),
                "uncertain_pathogens_excluded": sorted(uncertain),
                "predicted_uncertain_pathogens_excluded": sorted(uncertain_predictions),
                "candidates": sorted(candidates),
                "predicted": sorted(predicted),
                "abstained": sorted(abstained),
                "tp_pathogens": sorted(tp_set),
                "fp_pathogens": sorted(fp_set),
                "fn_pathogens": sorted(fn_set),
                "truth_missing_from_candidate_pool": sorted(truth - candidates),
                "tp": len(tp_set),
                "fp": len(fp_set),
                "fn": len(fn_set),
                "conditional_fn": len(conditional_fn),
                "truth_count": len(truth),
                "truth_in_pool_count": len(truth_in_pool),
                "pool_truth_count": len(truth) if pool_observed else 0,
                "pool_observed": pool_observed,
                "candidate_count": len(candidates),
                "predicted_count": len(predicted),
                "abstained_count": len(abstained),
                "case_recall": case_recall,
                "exact_match": truth == predicted,
                "baseline_miss_rescued": len((truth - baseline) & predicted),
                "baseline_false_positives_added": len((predicted - baseline) - truth),
            }
        )
    return rows, missing


def evaluate_outputs(
    output_payloads: List[Dict[str, Any]],
    gold: Dict[str, Any],
    *,
    bootstrap_iterations: int = 1000,
    bootstrap_seed: int = 20260810,
    allow_missing_outputs: bool = False,
    allow_partial_outputs: bool = False,
    allow_mixed_provenance: bool = False,
    analysis_mode: str = "full_abc",
) -> Dict[str, Any]:
    if analysis_mode not in {"full_abc", "bc_only"}:
        raise ValueError("analysis_mode must be full_abc or bc_only")
    outputs = {str(row.get("patient_id") or ""): row for row in output_payloads}
    if len(outputs) != len(output_payloads):
        raise ValueError("RAG_re outputs contain duplicate patient_id values")
    if "" in outputs:
        raise ValueError("Every RAG_re output must contain patient_id")
    labeled_ids = {str(row["patient_id"]) for row in gold["cases"]}
    missing_ids = sorted(labeled_ids - set(outputs))
    if missing_ids and not allow_missing_outputs:
        raise ValueError(
            "Confirmatory evaluation requires every labeled output; missing patient IDs: "
            + ", ".join(missing_ids)
        )
    partial_ids = sorted(
        patient_id
        for patient_id, payload in outputs.items()
        if payload.get("run_status") != "complete"
    )
    if partial_ids and not allow_partial_outputs:
        raise ValueError(
            "Confirmatory evaluation rejects partial/non-complete artifacts: "
            + ", ".join(partial_ids)
        )
    literature_not_attempted = sorted(
        patient_id
        for patient_id, payload in outputs.items()
        if (payload.get("execution") or {}).get("literature_attempted") is not True
    )
    if analysis_mode == "full_abc" and literature_not_attempted:
        raise ValueError(
            "Full A/B/C evaluation requires Module A to have been attempted; use the "
            "explicit --bc-only analysis for patient IDs: "
            + ", ".join(literature_not_attempted)
        )
    missing_pipeline_fingerprint = sorted(
        patient_id
        for patient_id, payload in outputs.items()
        if not payload.get("pipeline_fingerprint")
    )
    if missing_pipeline_fingerprint and not allow_mixed_provenance:
        raise ValueError(
            "Confirmatory evaluation requires pipeline fingerprints; missing patient IDs: "
            + ", ".join(missing_pipeline_fingerprint)
        )
    experiment_sets = [set((payload.get("experiments") or {}).keys()) for payload in output_payloads]
    if experiment_sets and any(names != experiment_sets[0] for names in experiment_sets[1:]):
        raise ValueError("All artifacts must carry the same complete experiment set")
    provenance_fields = {
        "config_hash": {str(payload.get("config_hash")) for payload in output_payloads},
        "prompt_sha256": {str(payload.get("prompt_sha256")) for payload in output_payloads},
        "candidate_pool_provenance": {
            str((payload.get("input_meta") or {}).get("candidate_pool_provenance"))
            for payload in output_payloads
        },
        "evidence_source_config_hash": {
            str((payload.get("evidence_provenance") or {}).get("source_config_hash"))
            for payload in output_payloads
        },
        "pipeline_fingerprint": {
            str(payload.get("pipeline_fingerprint")) for payload in output_payloads
        },
        "pipeline_version": {
            str(payload.get("pipeline_version")) for payload in output_payloads
        },
    }
    mixed = {key: sorted(values) for key, values in provenance_fields.items() if len(values) > 1}
    if mixed and not allow_mixed_provenance:
        raise ValueError(f"Confirmatory evaluation rejects mixed provenance: {mixed}")
    experiment_names = sorted(
        {
            name
            for payload in output_payloads
            for name in (payload.get("experiments") or {}).keys()
        }
    )
    if analysis_mode == "bc_only":
        allowed_bc = {
            "E0_BASELINE",
            "B_MNGS",
            "C_DIRECT_SUPPORT",
            "BC_OR",
            "E0_OR_C",
            "E0_AND_B",
            "E0_AND_C",
            "E0_AND_B_OR_C",
            "E0_AND_B_AND_C",
        }
        experiment_names = [name for name in experiment_names if name in allowed_bc]
    aliases = gold["aliases"]
    results: Dict[str, Any] = {}
    all_missing = set()
    for index, experiment in enumerate(experiment_names):
        rows, missing = _case_rows(outputs, gold["cases"], experiment, aliases)
        all_missing.update(missing)
        metrics = _aggregate(rows)
        results[experiment] = {
            "semantics": next(
                (
                    (payload.get("experiments") or {}).get(experiment, {}).get("semantics")
                    for payload in output_payloads
                    if experiment in (payload.get("experiments") or {})
                ),
                None,
            ),
            "expression": next(
                (
                    (payload.get("experiments") or {}).get(experiment, {}).get("expression")
                    for payload in output_payloads
                    if experiment in (payload.get("experiments") or {})
                ),
                None,
            ),
            "metrics": metrics,
            "patient_cluster_bootstrap_95_ci": _bootstrap(
                rows, bootstrap_iterations, bootstrap_seed + index
            ),
            "cases": rows,
        }

    strict_results: Dict[str, Any] = {}
    for experiment in experiment_names:
        strict_rows, _ = _case_rows(outputs, gold["cases"], experiment, {})
        strict_results[experiment] = _aggregate(strict_rows)

    warnings = []
    if not gold["complete_candidate_negative_labels"]:
        warnings.append(
            "Gold labels are positive lists without complete candidate-level negatives; "
            "precision is reported as agreement precision, not proven clinical precision."
        )
    if all_missing:
        warnings.append(
            "Some labeled cases had no RAG_re output; they were counted as empty predictions "
            "and their gold pathogens count as false negatives."
        )
        warnings.append(
            "Candidate-pool ceiling and conditional recall exclude missing artifacts because "
            "their candidate universes are unknown; end-to-end recall still counts their FN."
        )
    if partial_ids:
        warnings.append(
            "Exploratory override included partial/non-complete artifacts: "
            + ", ".join(partial_ids)
        )
    if mixed:
        warnings.append(f"Exploratory override included mixed provenance: {mixed}")
    if missing_pipeline_fingerprint:
        warnings.append(
            "Exploratory override included artifacts without pipeline fingerprints: "
            + ", ".join(missing_pipeline_fingerprint)
        )
    if analysis_mode == "bc_only":
        warnings.append(
            "This report is explicitly restricted to E0/B/C/BC/E0+C; Module A and all "
            "A-dependent policies are not evaluated."
        )
    warnings.append(
        "Genus/family-relaxed matching is intentionally disabled; primary matching uses only "
        "exact normalized names and the frozen accepted-synonym map."
    )
    comparisons = {
        "ADD_A_TO_E0": ("E0_BASELINE", "E0_OR_A"),
        "ADD_B_AFTER_E0_A": ("E0_OR_A", "E0_OR_A_OR_B"),
        "ADD_C_AFTER_E0_AB": ("E0_OR_A_OR_B", "E0_OR_A_OR_B_OR_C"),
        "ADD_B_TO_A": ("A_LITERATURE", "AB_OR"),
        "ADD_C_TO_AB": ("AB_OR", "ABC_OR"),
        "ADD_C_TO_B": ("B_MNGS", "BC_OR"),
        "GATED_VS_E0": ("E0_BASELINE", "GATED_C_OR_A_AND_B"),
        "BASELINE_PRESERVING_GATED_VS_E0": (
            "E0_BASELINE",
            "E0_OR_C_OR_A_AND_B",
        ),
        "PRUNE_WITH_B_VS_E0": ("E0_BASELINE", "E0_AND_B"),
        "PRUNE_WITH_C_VS_E0": ("E0_BASELINE", "E0_AND_C"),
        "PRUNE_WITH_B_OR_C_VS_E0": ("E0_BASELINE", "E0_AND_B_OR_C"),
        "PRUNE_WITH_B_AND_C_VS_E0": ("E0_BASELINE", "E0_AND_B_AND_C"),
    }
    incremental: Dict[str, Any] = {}
    for label, (before_name, after_name) in comparisons.items():
        if before_name not in results or after_name not in results:
            continue
        before = results[before_name]["metrics"]
        after = results[after_name]["metrics"]

        def delta(key: str) -> float | int | None:
            left, right = before.get(key), after.get(key)
            if left is None or right is None:
                return None
            return right - left

        incremental[label] = {
            "before": before_name,
            "after": after_name,
            "delta_tp": delta("tp"),
            "delta_fp": delta("fp"),
            "delta_fn": delta("fn"),
            "delta_agreement_precision": delta("agreement_precision"),
            "delta_end_to_end_recall": delta("end_to_end_recall"),
            "delta_f2": delta("f2"),
            "delta_false_positives_per_case": delta("false_positives_per_case"),
            "paired_patient_bootstrap_delta_95_ci": _bootstrap_delta(
                results[before_name]["cases"],
                results[after_name]["cases"],
                bootstrap_iterations,
                bootstrap_seed + 10000 + len(incremental),
            ),
        }
    return {
        "schema_version": "rag_re.evaluation.v1",
        "evaluator_version": EVALUATOR_VERSION,
        "evaluator_fingerprint": EVALUATOR_FINGERPRINT,
        "analysis_mode": analysis_mode,
        "gold_sha256": gold["source_sha256"],
        "label_semantics": gold["label_semantics"],
        "name_matching_primary": "exact_species_plus_frozen_accepted_synonyms",
        "accepted_aliases": aliases,
        "labeled_case_count": len(gold["cases"]),
        "evaluated_output_count": len(labeled_ids & set(outputs)),
        "missing_output_patient_ids": sorted(all_missing),
        "unlabeled_output_patient_ids": sorted(set(outputs) - labeled_ids),
        "bootstrap": {
            "unit": "patient",
            "iterations": bootstrap_iterations,
            "seed": bootstrap_seed,
        },
        "artifact_provenance": {
            key: sorted(values) for key, values in provenance_fields.items()
        },
        "experiments": results,
        "incremental_deltas": incremental,
        "strict_exact_name_sensitivity": strict_results,
        "warnings": warnings,
    }


def render_markdown(report: Dict[str, Any]) -> str:
    lines = [
        "# RAG_re ablation evaluation",
        "",
        f"- Labeled cases: {report['labeled_case_count']}",
        f"- Evaluated outputs: {report['evaluated_output_count']}",
        f"- Analysis mode: `{report['analysis_mode']}`",
        f"- Primary name matching: `{report['name_matching_primary']}`",
        f"- Evaluator fingerprint: `{report['evaluator_fingerprint']}`",
        "",
        "| Experiment | Semantics | Precision* | End-to-end recall | Conditional recall | Pool ceiling | F2 | FP/case | Negative specificity |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]

    def fmt(value: Any, percent: bool = True) -> str:
        if value is None:
            return "NA"
        return f"{value:.1%}" if percent else f"{value:.3f}"

    for name, data in report["experiments"].items():
        metrics = data["metrics"]
        lines.append(
            "| {name} | {semantics} | {precision} | {recall} | {conditional} | "
            "{ceiling} | {f2} | {fp_case} | {specificity} |".format(
                name=name,
                semantics=data.get("semantics") or "-",
                precision=fmt(metrics.get("agreement_precision")),
                recall=fmt(metrics.get("end_to_end_recall")),
                conditional=fmt(metrics.get("conditional_classifier_recall")),
                ceiling=fmt(metrics.get("candidate_pool_recall_ceiling")),
                f2=fmt(metrics.get("f2")),
                fp_case=fmt(metrics.get("false_positives_per_case"), False),
                specificity=fmt(metrics.get("negative_case_specificity")),
            )
        )
    lines.extend(["", "\\* Precision is agreement precision unless complete negatives were adjudicated."])
    if report.get("incremental_deltas"):
        lines.extend(
            [
                "",
                "## Incremental changes",
                "",
                "| Comparison | ΔTP | ΔFP | ΔFN | Δprecision | Δrecall | ΔFP/case |",
                "|---|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for name, row in report["incremental_deltas"].items():
            lines.append(
                f"| {name} | {row['delta_tp']} | {row['delta_fp']} | {row['delta_fn']} | "
                f"{fmt(row['delta_agreement_precision'])} | "
                f"{fmt(row['delta_end_to_end_recall'])} | "
                f"{fmt(row['delta_false_positives_per_case'], False)} |"
            )
    if report.get("warnings"):
        lines.extend(["", "## Warnings", ""])
        lines.extend(f"- {warning}" for warning in report["warnings"])
    return "\n".join(lines) + "\n"

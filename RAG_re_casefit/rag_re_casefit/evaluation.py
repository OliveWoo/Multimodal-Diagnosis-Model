from __future__ import annotations

import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List


def _compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def _normal(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", _compact(value).casefold())


ALIASES = {
    "cmv": "humancytomegalovirus",
    "humancytomegalovirus": "humancytomegalovirus",
    "hsv": "herpessimplexvirustype1",
    "hsv1": "herpessimplexvirustype1",
    "herpessimplexvirus1": "herpessimplexvirustype1",
    "herpessimplexvirustype1": "herpessimplexvirustype1",
    "vzv": "varicellazostervirus",
    "varicellazostervirus": "varicellazostervirus",
    "pneumocystisjiroveci": "pneumocystisjirovecii",
    "pjp": "pneumocystisjirovecii",
    "candidaalbican": "candidaalbicans",
    "covid19": "sarscov2",
    "sarscov2": "sarscov2",
}


def canonical(value: Any) -> str:
    normalized = _normal(value)
    return ALIASES.get(normalized, normalized)


def _metrics(rows: List[Dict[str, Any]], selected: Iterable[int]) -> Dict[str, Any]:
    chosen = set(selected)
    tp = sum(index in chosen and row["gold"] for index, row in enumerate(rows))
    fp = sum(index in chosen and not row["gold"] for index, row in enumerate(rows))
    fn = sum(index not in chosen and row["gold"] for index, row in enumerate(rows))
    tn = sum(index not in chosen and not row["gold"] for index, row in enumerate(rows))
    precision = tp / (tp + fp) if tp + fp else None
    recall = tp / (tp + fn) if tp + fn else None
    specificity = tn / (tn + fp) if tn + fp else None

    def fbeta(beta: float) -> float | None:
        if precision is None or recall is None or precision + recall == 0:
            return None
        b2 = beta * beta
        return (1 + b2) * precision * recall / (b2 * precision + recall)

    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / denominator) if denominator else None
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "predicted": tp + fp,
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "f1": fbeta(1.0),
        "f0_5": fbeta(0.5),
        "mcc": mcc,
    }


def _sweep(
    rows: List[Dict[str, Any]],
    name: str,
    score: Callable[[Dict[str, Any]], int],
    eligible: Callable[[Dict[str, Any]], bool] | None = None,
) -> Dict[str, Any]:
    eligible = eligible or (lambda row: True)
    max_score = max((score(row) for row in rows), default=0)
    thresholds = []
    for k in range(0, max(10, max_score) + 1):
        chosen = [
            index for index, row in enumerate(rows)
            if eligible(row) and score(row) >= k
        ]
        thresholds.append({"k": k, **_metrics(rows, chosen)})

    def best(metric: str) -> Dict[str, Any] | None:
        valid = [row for row in thresholds if row.get(metric) is not None]
        if not valid:
            return None
        return max(
            valid,
            key=lambda row: (
                row[metric],
                row.get("precision") if row.get("precision") is not None else -1,
                row["k"],
            ),
        )

    zero_fp = [row for row in thresholds if row["fp"] == 0 and row["predicted"] > 0]
    perfect = [row for row in thresholds if row["fp"] == 0 and row["fn"] == 0]
    return {
        "name": name,
        "thresholds": thresholds,
        "best_f1": best("f1"),
        "best_f0_5": best("f0_5"),
        "best_zero_fp": max(zero_fp, key=lambda row: (row["recall"], row["k"])) if zero_fp else None,
        "perfect_separation_thresholds": [row["k"] for row in perfect],
    }


def _distribution(rows: List[Dict[str, Any]], key: str, label: bool) -> Dict[str, Any]:
    values = [int(row[key]) for row in rows if row["gold"] is label]
    return {
        "n": len(values),
        "values": sorted(values),
        "mean": statistics.mean(values) if values else None,
        "median": statistics.median(values) if values else None,
        "min": min(values) if values else None,
        "max": max(values) if values else None,
    }


def _fmt(value: Any, digits: int = 3) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def evaluate(output_dir: Path, gold_path: Path) -> Dict[str, Any]:
    gold_root = json.loads(gold_path.read_text(encoding="utf-8-sig"))
    gold_by_patient = {
        _compact(case.get("patient_id")): {canonical(value) for value in case.get("pathogens", [])}
        for case in gold_root.get("cases", [])
        if isinstance(case, dict)
    }
    rows = []
    for path in sorted(output_dir.glob("*_casefit.json")):
        artifact = json.loads(path.read_text(encoding="utf-8-sig"))
        if artifact.get("schema_version") != "rag_re_casefit.candidate_output.v1":
            raise ValueError(f"unexpected schema in {path}")
        patient_id = _compact(artifact.get("patient_id"))
        if patient_id not in gold_by_patient:
            raise ValueError(f"patient {patient_id} has no gold")
        organism = _compact(artifact.get("canonical_organism"))
        old_a = artifact.get("old_A") or {}
        a1 = artifact.get("A1_aggregate") or {}
        rows.append({
            "patient_id": patient_id,
            "organism": organism,
            "organism_key": canonical(organism),
            "gold": canonical(organism) in gold_by_patient[patient_id],
            "old_A_positive": old_a.get("positive") is True,
            "old_A_support": int(old_a.get("support_count") or 0),
            "old_A_judgeable": int(old_a.get("judgeable_count") or 0),
            "A1_strong": int(a1.get("strong_match_count") or 0),
            "A1_partial": int(a1.get("partial_match_count") or 0),
            "A1_match": int(a1.get("strong_match_count") or 0) + int(a1.get("partial_match_count") or 0),
            "A1_counter": int(a1.get("counterevidence_count") or 0),
            "A1_tier": _compact(a1.get("casefit_evidence_tier")),
            "review_tier": _compact(artifact.get("review_tier")),
            "baseline_article_source": _compact((artifact.get("provenance") or {}).get("baseline_article_source")),
            "file": path.name,
        })
    if not rows:
        raise ValueError("no candidate outputs")
    if len({(row["patient_id"], row["organism_key"]) for row in rows}) != len(rows):
        raise ValueError("duplicate patient-organism candidate")

    sweeps = [
        _sweep(rows, "old_A_support_ge_k", lambda row: row["old_A_support"]),
        _sweep(rows, "A1_strong_ge_k", lambda row: row["A1_strong"]),
        _sweep(rows, "A1_match_ge_k", lambda row: row["A1_match"]),
        _sweep(
            rows,
            "A1_counterfree_match_ge_k",
            lambda row: row["A1_match"],
            lambda row: row["A1_counter"] == 0,
        ),
        _sweep(
            rows,
            "A1_strong1_counterfree_match_ge_k",
            lambda row: row["A1_match"],
            lambda row: row["A1_strong"] >= 1 and row["A1_counter"] == 0,
        ),
    ]
    tier_policies = {}
    for name, accepted in {
        "high_only": {"high"},
        "high_or_moderate": {"high", "moderate"},
        "high_moderate_or_conflicted": {"high", "moderate", "conflicted"},
    }.items():
        tier_policies[name] = _metrics(
            rows,
            [index for index, row in enumerate(rows) if row["A1_tier"] in accepted],
        )

    by_organism: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_organism[row["organism"]].append(row)
    organism_distribution = []
    for organism, values in sorted(by_organism.items(), key=lambda item: (-len(item[1]), item[0])):
        organism_distribution.append({
            "organism": organism,
            "candidate_count": len(values),
            "gold_count": sum(row["gold"] for row in values),
            "old_A_support_values": [row["old_A_support"] for row in values],
            "A1_strong_values": [row["A1_strong"] for row in values],
            "A1_match_values": [row["A1_match"] for row in values],
            "A1_counter_values": [row["A1_counter"] for row in values],
            "A1_tiers": dict(Counter(row["A1_tier"] for row in values)),
        })

    result = {
        "schema_version": "rag_re_casefit.development_evaluation.v1",
        "scope": {
            "candidate_count": len(rows),
            "gold_candidate_count": sum(row["gold"] for row in rows),
            "non_gold_candidate_count": sum(not row["gold"] for row in rows),
            "patient_count": len({row["patient_id"] for row in rows}),
            "answer_positive_patients_only": True,
            "primary_name_matching": "exact_plus_frozen_synonyms",
            "estimand": "selected_high_context_pool_candidate_classification",
        },
        "old_A_binary": _metrics(
            rows, [index for index, row in enumerate(rows) if row["old_A_positive"]]
        ),
        "tier_policies": tier_policies,
        "sweeps": sweeps,
        "distributions": {
            key: {
                "gold": _distribution(rows, key, True),
                "non_gold": _distribution(rows, key, False),
            }
            for key in ("old_A_support", "A1_strong", "A1_match", "A1_counter")
        },
        "organism_distribution": organism_distribution,
        "candidate_rows": rows,
        "warnings": [
            "This cohort has been inspected and is development-only; selected k values are hypotheses, not test performance.",
            "The latest folder contains answer-positive patients only, so precision omits false positives from no-pathogen patients.",
            "Repeated organism names share literature evidence and are not independent biological clusters.",
            "Metrics measure separation within the preselected high/context pool, not end-to-end clinical recall.",
        ],
    }
    return result


def render_markdown(result: Dict[str, Any]) -> str:
    scope = result["scope"]
    lines = [
        "# A vs A1 development comparison",
        "",
        f"Scope: {scope['candidate_count']} non-picked high/context candidates from "
        f"{scope['patient_count']} answer-positive patients; {scope['gold_candidate_count']} gold and "
        f"{scope['non_gold_candidate_count']} non-gold candidates.",
        "",
        "## Binary/tier policies",
        "",
        "| policy | TP | FP | FN | precision | recall | F0.5 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    policies = {"old_A_positive": result["old_A_binary"], **result["tier_policies"]}
    for name, metric in policies.items():
        lines.append(
            f"| {name} | {metric['tp']} | {metric['fp']} | {metric['fn']} | "
            f"{_fmt(metric['precision'])} | {_fmt(metric['recall'])} | {_fmt(metric['f0_5'])} |"
        )
    lines += ["", "## Threshold sweep summaries", ""]
    for sweep in result["sweeps"]:
        lines.append(f"### {sweep['name']}")
        lines.append("")
        for label, key in (("best F1", "best_f1"), ("best F0.5", "best_f0_5"), ("best zero-FP", "best_zero_fp")):
            row = sweep[key]
            if row is None:
                lines.append(f"- {label}: none")
            else:
                lines.append(
                    f"- {label}: k={row['k']}, TP/FP/FN={row['tp']}/{row['fp']}/{row['fn']}, "
                    f"P={_fmt(row['precision'])}, R={_fmt(row['recall'])}, F0.5={_fmt(row['f0_5'])}"
                )
        perfect = sweep["perfect_separation_thresholds"]
        lines.append(f"- perfect separation k: {perfect if perfect else 'none'}")
        lines.append("")
    lines += [
        "## Label distributions",
        "",
        "| feature | gold values | non-gold values | gold median | non-gold median |",
        "|---|---|---|---:|---:|",
    ]
    for feature, values in result["distributions"].items():
        lines.append(
            f"| {feature} | {values['gold']['values']} | {values['non_gold']['values']} | "
            f"{values['gold']['median']} | {values['non_gold']['median']} |"
        )
    lines += [
        "",
        "## Organism distribution",
        "",
        "| organism | n | gold | old A support | A1 strong | A1 match | A1 counter | tiers |",
        "|---|---:|---:|---|---|---|---|---|",
    ]
    for row in result["organism_distribution"]:
        lines.append(
            f"| {row['organism']} | {row['candidate_count']} | {row['gold_count']} | "
            f"{row['old_A_support_values']} | {row['A1_strong_values']} | "
            f"{row['A1_match_values']} | {row['A1_counter_values']} | {row['A1_tiers']} |"
        )
    lines += ["", "## Warnings", ""]
    lines.extend(f"- {warning}" for warning in result["warnings"])
    return "\n".join(lines) + "\n"

from __future__ import annotations

import json
import math
import re
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


REPORT_SCHEMA = "rag_re_casefit.rescue_evaluation.v1"
CASEFIT_SCHEMA = "rag_re_casefit.candidate_output.v1"
CLINICAL_SCHEMA = "rag_re_clinical.output.v1"
VALID_REVIEW_TIERS = {"review_high_priority", "review_context_needed"}

FROZEN_SYNONYMS = {
    "cmv": "human cytomegalovirus",
    "hsv": "hsv-1",
    "herpes simplex virus 1": "hsv-1",
    "herpes simplex virus type 1": "hsv-1",
    "human alphaherpesvirus 1": "hsv-1",
    "human herpesvirus 1": "hsv-1",
    "human alphaherpesvirus 3": "vzv",
    "varicella-zoster virus": "vzv",
    "pjp": "pneumocystis jirovecii",
    "pneumocystis jiroveci": "pneumocystis jirovecii",
    "candida albican": "candida albicans",
    "covid-19": "sars-cov-2",
    "crkp": "klebsiella pneumoniae",
    "klebsiella pneumoniae group": "klebsiella pneumoniae",
}

GENUS_STOPWORDS = {
    "human", "herpes", "herpesvirus", "influenza", "parainfluenza",
    "respiratory", "virus", "viral",
}


class RescueEvaluationError(ValueError):
    pass


def _load(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RescueEvaluationError(f"could not read {path}: {exc}") from exc


def _normal(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).strip().lower()
    return re.sub(r"\s+", " ", text)


def _aliases(gold: Mapping[str, Any]) -> dict[str, str]:
    aliases = {_normal(k): _normal(v) for k, v in FROZEN_SYNONYMS.items()}
    raw = gold.get("aliases") or {}
    if not isinstance(raw, dict):
        raise RescueEvaluationError("gold.aliases must be an object")
    for key, value in raw.items():
        aliases[_normal(key)] = _normal(value)
    for start in list(aliases):
        current = start
        seen: set[str] = set()
        while current in aliases:
            if current in seen:
                raise RescueEvaluationError(f"alias cycle at {start!r}")
            seen.add(current)
            current = aliases[current]
        aliases[start] = current
    return aliases


def _canonical(value: Any, aliases: Mapping[str, str]) -> str:
    current = _normal(value)
    if not current:
        raise RescueEvaluationError("pathogen name must not be empty")
    seen: set[str] = set()
    while current in aliases:
        if current in seen:
            raise RescueEvaluationError(f"alias cycle for {value!r}")
        seen.add(current)
        current = aliases[current]
    return current


def _genus(value: str) -> str | None:
    match = re.match(r"^([a-z][a-z-]+)\s+([a-z][a-z0-9.-]+)", value)
    if not match or match.group(1) in GENUS_STOPWORDS:
        return None
    return match.group(1)


def _matches(left: str, right: str, *, genus_relaxed: bool) -> bool:
    if left == right:
        return True
    genus = _genus(left)
    return bool(genus_relaxed and genus and genus == _genus(right))


def _maximum_match_assignments(
    truth: Sequence[str], predictions: Sequence[str], *, genus_relaxed: bool
) -> dict[int, int]:
    edges = [
        [j for j, predicted in enumerate(predictions) if _matches(gold, predicted, genus_relaxed=genus_relaxed)]
        for gold in truth
    ]
    assigned: dict[int, int] = {}

    def augment(gold_index: int, visited: set[int]) -> bool:
        for prediction_index in edges[gold_index]:
            if prediction_index in visited:
                continue
            visited.add(prediction_index)
            if prediction_index not in assigned or augment(assigned[prediction_index], visited):
                assigned[prediction_index] = gold_index
                return True
        return False

    for index in range(len(truth)):
        augment(index, set())
    return assigned


def _maximum_match_count(
    truth: Sequence[str], predictions: Sequence[str], *, genus_relaxed: bool
) -> int:
    return len(
        _maximum_match_assignments(
            truth, predictions, genus_relaxed=genus_relaxed
        )
    )


def _metrics(
    gold_by_patient: Mapping[str, Sequence[str]],
    predictions: Mapping[str, Sequence[str]],
    *,
    genus_relaxed: bool,
) -> dict[str, Any]:
    tp = fp = fn = 0
    for patient_id, truth in gold_by_patient.items():
        predicted = tuple(predictions.get(patient_id, ()))
        matched = _maximum_match_count(truth, predicted, genus_relaxed=genus_relaxed)
        tp += matched
        fp += len(predicted) - matched
        fn += len(truth) - matched
    precision = tp / (tp + fp) if tp + fp else None
    recall = tp / (tp + fn) if tp + fn else None

    def fbeta(beta: float) -> float | None:
        if precision is None or recall is None or precision + recall == 0:
            return None
        b2 = beta * beta
        return (1 + b2) * precision * recall / (b2 * precision + recall)

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "predicted": tp + fp,
        "precision": precision,
        "recall": recall,
        "f0_5": fbeta(0.5),
        "f1": fbeta(1.0),
    }


def _strict_bool_or_none(value: Any, *, where: str) -> bool | None:
    if value is True or value is False or value is None:
        return value
    raise RescueEvaluationError(f"{where} must be boolean or null")


def _patient_id_from_merge(payload: Mapping[str, Any], path: Path) -> str:
    patient_id = str(payload.get("patient_id") or "").strip()
    if not patient_id:
        match = re.search(r"patient_(\d+)", path.name, re.IGNORECASE)
        patient_id = match.group(1) if match else ""
    if not patient_id:
        raise RescueEvaluationError(f"missing patient_id in {path}")
    return patient_id


def _clinical_rule_evaluator():
    """Load the sibling clinical rule package without copying its rule logic."""

    workspace = Path(__file__).resolve().parents[2]
    clinical_root = workspace / "RAG_re_clinical"
    if not clinical_root.is_dir():
        raise RescueEvaluationError(f"missing sibling clinical package: {clinical_root}")
    root_text = str(clinical_root)
    if root_text not in sys.path:
        sys.path.insert(0, root_text)
    try:
        from rag_re_clinical.clinical_rules import evaluate_modules
    except ImportError as exc:
        raise RescueEvaluationError(f"could not import clinical rules: {exc}") from exc
    return evaluate_modules


def _candidate_name(row: Mapping[str, Any]) -> Any:
    for key in (
        "canonical_organism_name", "organism_name", "organism", "pathogen_name", "pathogen", "name"
    ):
        if row.get(key) not in (None, ""):
            return row.get(key)
    return ""


def _load_gold_for_patients(
    gold_path: Path, patient_ids: Iterable[str]
) -> tuple[dict[str, tuple[str, ...]], dict[str, str], Mapping[str, Any]]:
    gold = _load(gold_path)
    if not isinstance(gold, dict) or not isinstance(gold.get("cases"), list):
        raise RescueEvaluationError("invalid gold schema")
    aliases = _aliases(gold)
    wanted = set(patient_ids)
    result: dict[str, tuple[str, ...]] = {}
    for case in gold["cases"]:
        if not isinstance(case, dict):
            raise RescueEvaluationError("gold case must be an object")
        patient_id = str(case.get("patient_id") or "").strip()
        if patient_id not in wanted:
            continue
        values = tuple(_canonical(value, aliases) for value in case.get("pathogens", []))
        if len(values) != len(set(values)):
            raise RescueEvaluationError(f"duplicate gold pathogen for patient {patient_id}")
        result[patient_id] = values
    missing = sorted(wanted - set(result))
    if missing:
        raise RescueEvaluationError(f"merge patients missing from gold: {missing}")
    return result, aliases, gold


def evaluate_rescue_arms(
    casefit_dir: Path,
    merge_dir: Path,
    gold_path: Path,
    clinical_config_path: Path,
) -> dict[str, Any]:
    merge_payloads: dict[str, tuple[Path, Mapping[str, Any]]] = {}
    for path in sorted(merge_dir.glob("*.json")):
        payload = _load(path)
        if not isinstance(payload, dict):
            continue
        patient_id = _patient_id_from_merge(payload, path)
        if patient_id in merge_payloads:
            raise RescueEvaluationError(f"duplicate merge patient {patient_id}")
        merge_payloads[patient_id] = (path, payload)
    if not merge_payloads:
        raise RescueEvaluationError("no merge patient JSON files")

    gold_by_patient, aliases, gold_root = _load_gold_for_patients(gold_path, merge_payloads)
    answer_positive_patients_only = all(bool(value) for value in gold_by_patient.values())

    baseline: dict[str, tuple[str, ...]] = {}
    for patient_id, (path, payload) in merge_payloads.items():
        deterministic = payload.get("deterministic_max") or {}
        summary = deterministic.get("best_available_summary") or {}
        picked = summary.get("picked_pathogens")
        if not isinstance(picked, list):
            raise RescueEvaluationError(f"missing picked_pathogens in {path}")
        names = tuple(
            _canonical(item.get("organism_name"), aliases)
            for item in picked
            if isinstance(item, dict)
        )
        if len(names) != len(picked) or len(names) != len(set(names)):
            raise RescueEvaluationError(f"invalid/duplicate picked pathogen in {path}")
        baseline[patient_id] = names

    clinical_config = _load(clinical_config_path)
    if not isinstance(clinical_config, dict) or clinical_config.get("schema_version") != "rag_re_clinical.config.v1":
        raise RescueEvaluationError(f"unexpected clinical config: {clinical_config_path}")
    evaluate_clinical_modules = _clinical_rule_evaluator()

    review_by_key: dict[tuple[str, str], tuple[str, Mapping[str, Any]]] = {}
    deterministic_by_key: dict[tuple[str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for patient_id, (_, payload) in merge_payloads.items():
        review = payload.get("llm_missed_candidate_review") or {}
        for tier in sorted(VALID_REVIEW_TIERS):
            tier_rows = review.get(tier)
            if not isinstance(tier_rows, list):
                raise RescueEvaluationError(f"patient {patient_id} missing {tier}")
            for row in tier_rows:
                if not isinstance(row, dict):
                    raise RescueEvaluationError(f"patient {patient_id} invalid {tier} row")
                organism = _canonical(_candidate_name(row), aliases)
                key = (patient_id, organism)
                if key in review_by_key:
                    raise RescueEvaluationError(f"duplicate latest review candidate {key}")
                review_by_key[key] = (tier, row)
        candidates = ((payload.get("deterministic_max") or {}).get("pathogen_candidates"))
        if not isinstance(candidates, list):
            raise RescueEvaluationError(f"patient {patient_id} missing deterministic candidates")
        for row in candidates:
            if not isinstance(row, dict):
                raise RescueEvaluationError(f"patient {patient_id} invalid deterministic candidate")
            organism = _canonical(_candidate_name(row), aliases)
            key = (patient_id, organism)
            deterministic_by_key[key].append(row)

    rows: list[dict[str, Any]] = []
    for path in sorted(casefit_dir.glob("*_casefit.json")):
        artifact = _load(path)
        if not isinstance(artifact, dict) or artifact.get("schema_version") != CASEFIT_SCHEMA:
            raise RescueEvaluationError(f"unexpected casefit schema in {path}")
        patient_id = str(artifact.get("patient_id") or "").strip()
        if patient_id not in merge_payloads:
            raise RescueEvaluationError(f"casefit patient {patient_id} is outside merge cohort")
        review_tier = str(artifact.get("review_tier") or "").strip()
        if review_tier not in VALID_REVIEW_TIERS:
            raise RescueEvaluationError(f"casefit candidate outside review scope in {path}")
        organism = _canonical(artifact.get("canonical_organism"), aliases)
        key = (patient_id, organism)
        review_match = review_by_key.get(key)
        if review_match is None:
            raise RescueEvaluationError(f"casefit candidate has no latest review match: {key}")
        if review_match[0] != review_tier:
            raise RescueEvaluationError(f"latest review tier mismatch for {key}")
        if organism in baseline[patient_id]:
            raise RescueEvaluationError(f"latest review candidate already exists in picked: {key}")
        review_row = review_match[1]
        deterministic_matches = deterministic_by_key.get(key, [])
        display_name = str(artifact.get("canonical_organism") or "").strip()
        exact_display_matches = [
            row for row in deterministic_matches
            if _normal(_candidate_name(row)) == _normal(display_name)
        ]
        if len(exact_display_matches) == 1:
            deterministic = exact_display_matches[0]
        elif len(deterministic_matches) == 1:
            deterministic = deterministic_matches[0]
        elif not deterministic_matches:
            deterministic = None
        else:
            raise RescueEvaluationError(
                f"ambiguous latest deterministic candidates for {key}: "
                f"{[_candidate_name(row) for row in deterministic_matches]}"
            )
        snapshot = review_row.get("evidence_snapshot")
        snapshot = snapshot if isinstance(snapshot, dict) else {}
        effective = dict(deterministic or snapshot)
        effective["organism_name"] = display_name
        effective["canonical_organism_name"] = display_name
        case_card = artifact.get("case_card") or {}
        source_context = {
            "site_code": case_card.get("clinical_site"),
            "target_site": case_card.get("clinical_site"),
            "specimen_site_code": case_card.get("specimen_site"),
            "specimen_class": effective.get("specimen_class"),
            "specimen_alignment": effective.get("specimen_alignment"),
            "dominant_source": case_card.get("target_syndrome"),
        }
        merged_context = {
            "source_context": source_context,
            "effective_candidate": effective,
            "deterministic_candidate": deterministic,
            "review_metadata": {"tier": review_tier, "row": review_row},
        }
        modules = evaluate_clinical_modules(merged_context, effective, clinical_config)
        b = modules.get("B_STRICT") or {}
        c = modules.get("C") or {}
        d = modules.get("D") or {}
        b_positive = _strict_bool_or_none(b.get("positive"), where=f"{key}.B_STRICT.positive")
        c_grade = str(c.get("grade") or "").strip()
        if c_grade not in {"C0", "C1", "C2", "C3", "CNEG", "CU", "CX"}:
            raise RescueEvaluationError(f"invalid C grade for {key}: {c_grade!r}")
        a1 = artifact.get("A1_aggregate") or {}
        strong = a1.get("strong_match_count")
        partial = a1.get("partial_match_count")
        if type(strong) is not int or type(partial) is not int or strong < 0 or partial < 0:
            raise RescueEvaluationError(f"invalid A1 counts in {path}")
        rows.append({
            "patient_id": patient_id,
            "organism": organism,
            "display_organism": str(artifact.get("canonical_organism") or "").strip(),
            "review_tier": review_tier,
            "A1_strong": strong,
            "A1_partial": partial,
            "A1_match": strong + partial,
            "B_strict": b_positive,
            "B_absolute_axis": _strict_bool_or_none(
                b.get("absolute_strength_axis"), where=f"{key}.B_STRICT.absolute_strength_axis"
            ),
            "B_relative_axis": _strict_bool_or_none(
                b.get("relative_strength_axis"), where=f"{key}.B_STRICT.relative_strength_axis"
            ),
            "B_site_aligned": _strict_bool_or_none(
                b.get("site_aligned"), where=f"{key}.B_STRICT.site_aligned"
            ),
            "C_grade": c_grade,
            "D_state": str(d.get("state") or "").strip(),
            "clinical_features_source": (
                "latest_merge_deterministic_candidate" if deterministic is not None
                else "latest_merge_review_evidence_snapshot"
            ),
            "gold_exact": organism in gold_by_patient[patient_id],
            "source_file": path.name,
        })
    if not rows:
        raise RescueEvaluationError("no casefit candidates")
    keys = [(row["patient_id"], row["organism"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise RescueEvaluationError("duplicate casefit patient-organism")

    selectors = {
        "R0_PICKED": lambda row: False,
        "R1_A1_STRONG2_AND_B": lambda row: row["A1_strong"] >= 2 and row["B_strict"] is True,
        "R2_R1_OR_C3": lambda row: (row["A1_strong"] >= 2 and row["B_strict"] is True) or row["C_grade"] == "C3",
        "R3_A1_MATCH7_AND_B": lambda row: row["A1_match"] >= 7 and row["B_strict"] is True,
        "R4_STRONG2_AND_B_RELATIVE": lambda row: (
            row["C_grade"] == "C3"
            or (
                row["A1_strong"] >= 2
                and row["B_site_aligned"] is True
                and row["B_relative_axis"] is True
            )
        ),
        "R5_A1_UNION_AND_B_LOOSE": lambda row: (
            row["C_grade"] == "C3"
            or (
                (row["A1_strong"] >= 2 or row["A1_match"] >= 7)
                and row["B_site_aligned"] is True
                and (
                    row["B_absolute_axis"] is True
                    or row["B_relative_axis"] is True
                )
            )
        ),
    }
    arm_predictions: dict[str, dict[str, tuple[str, ...]]] = {}
    selected_rows: dict[str, list[dict[str, Any]]] = {}
    for arm, selector in selectors.items():
        additions: dict[str, list[str]] = defaultdict(list)
        selected_rows[arm] = []
        for row in rows:
            if selector(row):
                additions[row["patient_id"]].append(row["organism"])
                selected_rows[arm].append(row)
        predictions = {}
        for patient_id, picked in baseline.items():
            combined = tuple((*picked, *additions.get(patient_id, [])))
            if len(combined) != len(set(combined)):
                raise RescueEvaluationError(f"duplicate combined prediction for patient {patient_id}")
            predictions[patient_id] = combined
        arm_predictions[arm] = predictions

    scopes = {
        "answer_positive_exact_plus_synonyms": False,
        "answer_positive_genus_relaxed": True,
    }
    results: dict[str, Any] = {}
    for scope, genus_relaxed in scopes.items():
        baseline_metrics = _metrics(gold_by_patient, arm_predictions["R0_PICKED"], genus_relaxed=genus_relaxed)
        arms = {}
        for arm, predictions in arm_predictions.items():
            metric = _metrics(gold_by_patient, predictions, genus_relaxed=genus_relaxed)
            arms[arm] = {
                **metric,
                "delta_tp_vs_R0": metric["tp"] - baseline_metrics["tp"],
                "delta_fp_vs_R0": metric["fp"] - baseline_metrics["fp"],
                "delta_fn_vs_R0": metric["fn"] - baseline_metrics["fn"],
                "delta_precision_vs_R0": (
                    metric["precision"] - baseline_metrics["precision"]
                    if metric["precision"] is not None and baseline_metrics["precision"] is not None
                    else None
                ),
                "delta_recall_vs_R0": (
                    metric["recall"] - baseline_metrics["recall"]
                    if metric["recall"] is not None and baseline_metrics["recall"] is not None
                    else None
                ),
                "added_candidate_count": len(selected_rows[arm]),
            }
        results[scope] = arms

    selected_audit = {
        arm: [
            {key: row[key] for key in (
                "patient_id", "display_organism", "review_tier", "gold_exact",
                "A1_strong", "A1_partial", "A1_match", "B_strict", "C_grade", "D_state",
            )}
            for row in selected_rows[arm]
        ]
        for arm in selectors
    }
    review_keys = {(row["patient_id"], row["organism"]) for row in rows}
    missed_gold: dict[str, list[dict[str, Any]]] = {}
    for scope, genus_relaxed in scopes.items():
        missed_rows: list[dict[str, Any]] = []
        for patient_id, truth in gold_by_patient.items():
            predicted = baseline[patient_id]
            assignments = _maximum_match_assignments(
                truth, predicted, genus_relaxed=genus_relaxed
            )
            matched_gold_indexes = set(assignments.values())
            for index, organism in enumerate(truth):
                if index not in matched_gold_indexes:
                    missed_rows.append({
                        "patient_id": patient_id,
                        "organism": organism,
                        "in_high_context_casefit_pool_exact": (patient_id, organism) in review_keys,
                    })
        missed_gold[scope] = missed_rows
    warnings = [
        "The answers in this cohort have already been inspected; do not tune a final rule and report the same cohort as untouched test performance.",
        (
            "The merge cohort contains answer-positive patients only, so precision excludes no-pathogen patient false positives."
            if answer_positive_patients_only
            else "The merge cohort includes no-pathogen patients, so their selected predictions are counted as false positives."
        ),
        "R1/R3 apply the requested formulas exactly; D guard states are reported for audit but are not hidden eligibility gates.",
        "Repeated organisms share literature evidence, so selected candidates are not independent organism clusters.",
        "B_STRICT, C3 and D are recomputed from the latest merge with the frozen clinical_v1 rules.",
    ]
    return {
        "schema_version": REPORT_SCHEMA,
        "study_status": (
            "development_only_answers_inspected"
            if gold_root.get("development_only") is True
            else "frozen_predictions_posthoc_gold_join"
        ),
        "formulas": {
            "R0_PICKED": "picked (immutable)",
            "R1_A1_STRONG2_AND_B": "R0 OR (A1_strong >= 2 AND B_STRICT is true)",
            "R2_R1_OR_C3": "R1 OR (clinical C grade == C3)",
            "R3_A1_MATCH7_AND_B": "R0 OR (A1_strong + A1_partial >= 7 AND B_STRICT is true)",
            "R4_STRONG2_AND_B_RELATIVE": "R0 OR C3 OR (A1_strong >= 2 AND site aligned AND B relative axis is true)",
            "R5_A1_UNION_AND_B_LOOSE": "R0 OR C3 OR ((A1_strong >= 2 OR A1_match >= 7) AND site aligned AND (B absolute axis OR B relative axis))",
        },
        "scope": {
            "patient_count": len(merge_payloads),
            "picked_count": sum(len(value) for value in baseline.values()),
            "casefit_candidate_count": len(rows),
            "review_high_priority_count": sum(row["review_tier"] == "review_high_priority" for row in rows),
            "review_context_needed_count": sum(row["review_tier"] == "review_context_needed" for row in rows),
            "gold_pathogen_count": sum(len(value) for value in gold_by_patient.values()),
            "answer_positive_patients_only": answer_positive_patients_only,
        },
        "metrics": results,
        "selected_candidates": selected_audit,
        "R0_missed_gold": missed_gold,
        "candidate_rows": rows,
        "warnings": warnings,
    }


def _fmt(value: Any, digits: int = 3) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def render_rescue_markdown(report: Mapping[str, Any]) -> str:
    scope = report["scope"]
    lines = [
        "# A1 + B_STRICT + C3 rescue arms",
        "",
        f"Evaluation scope: {scope['patient_count']} patients, "
        f"{scope['picked_count']} immutable picked predictions, {scope['casefit_candidate_count']} "
        "non-picked high/context candidates.",
        "",
    ]
    for scope_name, arms in report["metrics"].items():
        lines += [
            f"## {scope_name}", "",
            "| arm | added | TP | FP | FN | precision | recall | F0.5 | dTP | dFP | dP | dR |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for arm, metric in arms.items():
            lines.append(
                f"| {arm} | {metric['added_candidate_count']} | {metric['tp']} | {metric['fp']} | "
                f"{metric['fn']} | {_fmt(metric['precision'])} | {_fmt(metric['recall'])} | "
                f"{_fmt(metric['f0_5'])} | {metric['delta_tp_vs_R0']:+d} | "
                f"{metric['delta_fp_vs_R0']:+d} | {_fmt(metric['delta_precision_vs_R0'])} | "
                f"{_fmt(metric['delta_recall_vs_R0'])} |"
            )
        lines.append("")
    lines += ["## Selected candidate audit", ""]
    for arm, rows in report["selected_candidates"].items():
        lines.append(f"### {arm}")
        lines.append("")
        if not rows:
            lines.append("- none")
        for row in rows:
            lines.append(
                f"- patient {row['patient_id']}: {row['display_organism']} "
                f"(gold_exact={row['gold_exact']}, A1={row['A1_strong']} strong + "
                f"{row['A1_partial']} partial, B={row['B_strict']}, C={row['C_grade']}, "
                f"D={row['D_state']})"
            )
        lines.append("")
    lines += ["## Warnings", ""]
    lines.extend(f"- {value}" for value in report["warnings"])
    return "\n".join(lines) + "\n"

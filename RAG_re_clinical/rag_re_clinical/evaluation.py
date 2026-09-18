from __future__ import annotations

import hashlib
import json
import math
import re
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


OUTPUT_SCHEMA_VERSION = "rag_re_clinical.output.v1"
MANIFEST_SCHEMA_VERSION = "rag_re_clinical.batch_manifest.v1"
REPORT_SCHEMA_VERSION = "rag_re_clinical.evaluation.v1"
BASELINE_ARM = "CL0_BASELINE"
FIXED_ARM_NAMES = (
    BASELINE_ARM,
    "CL1_DIRECT",
    "CL2_CONVERGENT",
    "CL3_BALANCED",
    "CL4_TIERED_EXPLORATORY",
)

# These mappings are part of the evaluator, not learned from the evaluation gold.
# A gold file may add a versioned, pre-specified alias table for its own labels.
FROZEN_SYNONYMS = {
    "cmv": "human cytomegalovirus",
    "hsv": "hsv-1",
    "herpes simplex virus 1": "hsv-1",
    "herpes simplex virus type 1": "hsv-1",
    "human alphaherpesvirus 1": "hsv-1",
    "human herpesvirus 1": "hsv-1",
    "pjp": "pneumocystis jirovecii",
    "pneumocystis jiroveci": "pneumocystis jirovecii",
    "candida albican": "candida albicans",
    "crkp": "klebsiella pneumoniae",
}

GENUS_STOPWORDS = {
    "human",
    "herpes",
    "herpesvirus",
    "influenza",
    "parainfluenza",
    "respiratory",
    "virus",
    "viral",
}


class EvaluationError(ValueError):
    """Raised when clinical artifacts or gold labels violate the evaluation contract."""


@dataclass(frozen=True)
class ArmView:
    predicted: tuple[str, ...]
    manual: tuple[str, ...]
    unknown: tuple[str, ...]


@dataclass(frozen=True)
class ArtifactView:
    patient_id: str
    path: str
    candidates: tuple[str, ...]
    baseline: tuple[str, ...]
    arms: Mapping[str, ArmView]


def _normalize(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).strip().lower()
    return re.sub(r"\s+", " ", text)


def _require_name(value: Any, *, where: str) -> str:
    normalized = _normalize(value)
    if not normalized:
        raise EvaluationError(f"{where} must be a non-empty pathogen name")
    return normalized


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise EvaluationError(f"Could not read JSON from {path}: {exc}") from exc


def _build_aliases(gold_aliases: Any) -> dict[str, str]:
    aliases = {
        _require_name(key, where="frozen alias key"): _require_name(
            value, where="frozen alias value"
        )
        for key, value in FROZEN_SYNONYMS.items()
    }
    if gold_aliases is None:
        gold_aliases = {}
    if not isinstance(gold_aliases, dict):
        raise EvaluationError("gold.aliases must be an object")
    for raw_key, raw_value in gold_aliases.items():
        key = _require_name(raw_key, where="gold alias key")
        value = _require_name(raw_value, where=f"gold alias {raw_key!r}")
        existing = aliases.get(key)
        if existing is not None and existing != value:
            raise EvaluationError(
                f"Gold alias {raw_key!r} conflicts with the frozen mapping "
                f"{existing!r} != {value!r}"
            )
        aliases[key] = value

    # Resolve every chain now so cycles fail before any metrics are produced.
    for start in list(aliases):
        current = start
        visited: set[str] = set()
        while current in aliases:
            if current in visited:
                raise EvaluationError(f"Alias cycle detected at {start!r}")
            visited.add(current)
            current = aliases[current]
        aliases[start] = current
    return aliases


def _canonical(value: Any, aliases: Mapping[str, str]) -> str:
    current = _require_name(value, where="pathogen")
    visited: set[str] = set()
    while current in aliases:
        if current in visited:
            raise EvaluationError(f"Alias cycle encountered for {value!r}")
        visited.add(current)
        current = aliases[current]
    return current


def _canonical_list(
    raw: Any,
    aliases: Mapping[str, str],
    *,
    where: str,
    reject_duplicates: bool,
) -> tuple[str, ...]:
    if not isinstance(raw, list):
        raise EvaluationError(f"{where} must be a list")
    values = tuple(_canonical(value, aliases) for value in raw)
    if reject_duplicates and len(set(values)) != len(values):
        raise EvaluationError(f"{where} contains duplicate or alias-equivalent names")
    if reject_duplicates:
        return values
    return tuple(dict.fromkeys(values))


def _genus(value: str) -> str | None:
    match = re.match(r"^([a-z][a-z-]+)(?:\s+)([a-z][a-z0-9.-]+)", value)
    if not match:
        return None
    token = match.group(1)
    return None if token in GENUS_STOPWORDS else token


def _names_match(left: str, right: str, *, genus_relaxed: bool) -> bool:
    if left == right:
        return True
    left_genus = _genus(left)
    return bool(genus_relaxed and left_genus and left_genus == _genus(right))


def _maximum_match_count(
    truth: Sequence[str], predictions: Sequence[str], *, genus_relaxed: bool
) -> int:
    """Maximum one-to-one pathogen matching within one patient."""

    edges = [
        [
            index
            for index, predicted in enumerate(predictions)
            if _names_match(gold, predicted, genus_relaxed=genus_relaxed)
        ]
        for gold in truth
    ]
    assigned: dict[int, int] = {}

    def augment(gold_index: int, visited: set[int]) -> bool:
        for prediction_index in edges[gold_index]:
            if prediction_index in visited:
                continue
            visited.add(prediction_index)
            if prediction_index not in assigned or augment(
                assigned[prediction_index], visited
            ):
                assigned[prediction_index] = gold_index
                return True
        return False

    return sum(augment(index, set()) for index in range(len(truth)))


def _resolve_arm_names(arm_names: Sequence[str] | None) -> tuple[str, ...]:
    if arm_names is None:
        return FIXED_ARM_NAMES
    if isinstance(arm_names, (str, bytes)):
        raise EvaluationError("arm_names must be a sequence, not a string")
    normalized = tuple(str(value).strip() for value in arm_names)
    if not normalized or any(not value for value in normalized):
        raise EvaluationError("arm_names must contain at least one non-empty arm")
    if len(set(normalized)) != len(normalized):
        raise EvaluationError("arm_names contains duplicates")
    unsupported = sorted(set(normalized) - set(FIXED_ARM_NAMES))
    if unsupported:
        raise EvaluationError(f"Unsupported clinical arm(s): {unsupported}")
    if BASELINE_ARM not in normalized:
        normalized = (BASELINE_ARM, *normalized)
    return normalized


def _strict_bool(value: Any, *, where: str) -> bool:
    if value is True:
        return True
    if value is False:
        return False
    raise EvaluationError(f"{where} must be a boolean")


def _validate_artifact(
    payload: Any,
    *,
    path: Path,
    aliases: Mapping[str, str],
    arm_names: Sequence[str],
) -> ArtifactView:
    if not isinstance(payload, dict):
        raise EvaluationError(f"Clinical artifact {path} must be an object")
    if payload.get("schema_version") != OUTPUT_SCHEMA_VERSION:
        raise EvaluationError(
            f"Clinical artifact {path} has unsupported schema_version "
            f"{payload.get('schema_version')!r}"
        )
    patient_id = str(payload.get("patient_id") or "").strip()
    if not patient_id:
        raise EvaluationError(f"Clinical artifact {path} is missing patient_id")
    if payload.get("run_status") != "complete":
        raise EvaluationError(
            f"Clinical artifact for patient {patient_id} must have run_status='complete'"
        )

    raw_candidates = payload.get("candidates")
    if not isinstance(raw_candidates, list):
        raise EvaluationError(f"Patient {patient_id} candidates must be a list")
    candidate_names: list[str] = []
    baseline_names: list[str] = []
    candidate_arm_values: dict[str, dict[str, bool]] = {
        arm: {} for arm in arm_names
    }
    for index, candidate in enumerate(raw_candidates):
        where = f"patient {patient_id} candidate[{index}]"
        if not isinstance(candidate, dict):
            raise EvaluationError(f"{where} must be an object")
        if "organism_name" not in candidate:
            raise EvaluationError(f"{where} is missing organism_name")
        name = _canonical(
            candidate.get("canonical_organism_name") or candidate.get("organism_name"),
            aliases,
        )
        if name in candidate_names:
            raise EvaluationError(
                f"Patient {patient_id} has duplicate or alias-equivalent candidate {name!r}"
            )
        candidate_names.append(name)
        baseline = _strict_bool(
            candidate.get("baseline_selected"), where=f"{where}.baseline_selected"
        )
        if baseline:
            baseline_names.append(name)
        if "review_tier" not in candidate:
            raise EvaluationError(f"{where} is missing review_tier")
        if not isinstance(candidate.get("disposition"), str) or not str(
            candidate.get("disposition")
        ).strip():
            raise EvaluationError(f"{where}.disposition must be a non-empty string")
        if not isinstance(candidate.get("modules"), dict):
            raise EvaluationError(f"{where}.modules must be an object")
        raw_candidate_arms = candidate.get("arms")
        if not isinstance(raw_candidate_arms, dict):
            raise EvaluationError(f"{where}.arms must be an object")
        for arm in arm_names:
            if arm not in raw_candidate_arms:
                raise EvaluationError(f"{where} is missing arm {arm}")
            value = raw_candidate_arms[arm]
            if value is not True and value is not False:
                raise EvaluationError(
                    f"{where}.arms[{arm}] must be a boolean; unknown/manual candidates "
                    "belong in the root manual_review_pathogens queue"
                )
            candidate_arm_values[arm][name] = value

    candidate_set = set(candidate_names)
    baseline_set = set(baseline_names)
    raw_arms = payload.get("arms")
    if not isinstance(raw_arms, dict):
        raise EvaluationError(f"Patient {patient_id} arms must be an object")
    arms: dict[str, ArmView] = {}
    for arm in arm_names:
        if arm not in raw_arms:
            raise EvaluationError(f"Patient {patient_id} is missing root arm {arm}")
        raw_arm = raw_arms[arm]
        if not isinstance(raw_arm, dict):
            raise EvaluationError(f"Patient {patient_id} arm {arm} must be an object")
        predicted = _canonical_list(
            raw_arm.get("predicted_pathogens"),
            aliases,
            where=f"patient {patient_id} arm {arm}.predicted_pathogens",
            reject_duplicates=True,
        )
        manual = _canonical_list(
            raw_arm.get("manual_review_pathogens"),
            aliases,
            where=f"patient {patient_id} arm {arm}.manual_review_pathogens",
            reject_duplicates=True,
        )
        if type(raw_arm.get("predicted_count")) is not int or raw_arm[
            "predicted_count"
        ] != len(predicted):
            raise EvaluationError(
                f"Patient {patient_id} arm {arm} predicted_count does not match its list"
            )
        if type(raw_arm.get("manual_review_count")) is not int or raw_arm[
            "manual_review_count"
        ] != len(manual):
            raise EvaluationError(
                f"Patient {patient_id} arm {arm} manual_review_count does not match its list"
            )
        predicted_set = set(predicted)
        manual_set = set(manual)
        outside = (predicted_set | manual_set) - candidate_set
        if outside:
            raise EvaluationError(
                f"Patient {patient_id} arm {arm} contains pathogens outside candidates: "
                f"{sorted(outside)}"
            )
        overlap = predicted_set & manual_set
        if overlap:
            raise EvaluationError(
                f"Patient {patient_id} arm {arm} marks pathogens both positive and manual: "
                f"{sorted(overlap)}"
            )
        candidate_positive_set = {
            name
            for name, decision in candidate_arm_values[arm].items()
            if decision is True
        }
        if predicted_set != candidate_positive_set:
            raise EvaluationError(
                f"Patient {patient_id} arm {arm} root predictions disagree with "
                "candidate-level arm decisions"
            )
        if not baseline_set <= predicted_set:
            raise EvaluationError(
                f"Patient {patient_id} arm {arm} violates E0 immutability; missing "
                f"{sorted(baseline_set - predicted_set)}"
            )
        inferred_unknown: set[str] = set()
        if "unknown_pathogens" in raw_arm:
            explicit_unknown = set(
                _canonical_list(
                    raw_arm.get("unknown_pathogens"),
                    aliases,
                    where=f"patient {patient_id} arm {arm}.unknown_pathogens",
                    reject_duplicates=True,
                )
            )
            if not explicit_unknown <= candidate_set:
                raise EvaluationError(
                    f"Patient {patient_id} arm {arm} unknown pathogens fall outside candidates"
                )
            if explicit_unknown & (predicted_set | manual_set):
                raise EvaluationError(
                    f"Patient {patient_id} arm {arm} unknown overlaps positive/manual"
                )
            if explicit_unknown != inferred_unknown:
                raise EvaluationError(
                    f"Patient {patient_id} arm {arm} explicit unknown list disagrees "
                    "with null candidate-level decisions"
                )
        arms[arm] = ArmView(
            predicted=predicted,
            manual=manual,
            unknown=tuple(name for name in candidate_names if name in inferred_unknown),
        )

    if set(arms[BASELINE_ARM].predicted) != baseline_set:
        raise EvaluationError(
            f"Patient {patient_id} CL0_BASELINE does not exactly equal baseline_selected"
        )
    return ArtifactView(
        patient_id=patient_id,
        path=str(path.resolve()),
        candidates=tuple(candidate_names),
        baseline=tuple(baseline_names),
        arms=arms,
    )


def _load_gold(path: Path) -> tuple[dict[str, Any], dict[str, str]]:
    payload = _load_json(path)
    if not isinstance(payload, dict):
        raise EvaluationError("Gold file must contain an object")
    cases = payload.get("cases")
    if not isinstance(cases, list) or not cases:
        raise EvaluationError("gold.cases must be a non-empty list")
    aliases = _build_aliases(payload.get("aliases"))
    normalized_cases: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise EvaluationError(f"gold.cases[{index}] must be an object")
        patient_id = str(case.get("patient_id") or "").strip()
        if not patient_id:
            raise EvaluationError(f"gold.cases[{index}] is missing patient_id")
        if patient_id in seen_ids:
            raise EvaluationError(f"Gold contains duplicate patient_id {patient_id}")
        seen_ids.add(patient_id)
        truth = _canonical_list(
            case.get("pathogens"),
            aliases,
            where=f"gold patient {patient_id}.pathogens",
            reject_duplicates=False,
        )
        uncertain = _canonical_list(
            case.get("uncertain_pathogens", []),
            aliases,
            where=f"gold patient {patient_id}.uncertain_pathogens",
            reject_duplicates=False,
        )
        overlap = set(truth) & set(uncertain)
        if overlap:
            raise EvaluationError(
                f"Gold patient {patient_id} has positive/uncertain overlap: {sorted(overlap)}"
            )
        normalized_cases.append(
            {
                "patient_id": patient_id,
                "truth": truth,
                "uncertain": uncertain,
            }
        )
    return (
        {
            "schema_version": payload.get("schema_version"),
            "label_semantics": payload.get("label_semantics"),
            "complete_candidate_negative_labels": payload.get(
                "complete_candidate_negative_labels"
            ),
            "cases": normalized_cases,
            "path": str(path.resolve()),
            "sha256": _sha256(path),
        },
        aliases,
    )


def _safe_div(numerator: int | float, denominator: int | float) -> float | None:
    return numerator / denominator if denominator else None


def _fbeta(precision: float | None, recall: float | None, beta: float) -> float | None:
    if precision is None or recall is None:
        return None
    if precision == 0 and recall == 0:
        return 0.0
    beta_squared = beta * beta
    denominator = beta_squared * precision + recall
    return (
        (1 + beta_squared) * precision * recall / denominator
        if denominator
        else 0.0
    )


def _evaluate_scope(
    *,
    cases: Sequence[Mapping[str, Any]],
    artifacts: Mapping[str, ArtifactView],
    arm_names: Sequence[str],
    positive_cases_only: bool,
    genus_relaxed: bool,
) -> dict[str, Any]:
    scoped_cases = [case for case in cases if not positive_cases_only or case["truth"]]
    arm_rows: dict[str, list[dict[str, Any]]] = {arm: [] for arm in arm_names}
    for case in scoped_cases:
        patient_id = str(case["patient_id"])
        truth = tuple(case["truth"])
        uncertain = set(case["uncertain"])
        artifact = artifacts.get(patient_id)
        for arm in arm_names:
            if artifact is None:
                predicted: tuple[str, ...] = ()
                baseline: tuple[str, ...] = ()
                manual: tuple[str, ...] = ()
                unknown: tuple[str, ...] = ()
                candidate_count = 0
            else:
                predicted = tuple(
                    value for value in artifact.arms[arm].predicted if value not in uncertain
                )
                baseline = tuple(
                    value
                    for value in artifact.arms[BASELINE_ARM].predicted
                    if value not in uncertain
                )
                manual = artifact.arms[arm].manual
                unknown = artifact.arms[arm].unknown
                candidate_count = len(artifact.candidates)
            matched = _maximum_match_count(
                truth, predicted, genus_relaxed=genus_relaxed
            )
            baseline_matched = _maximum_match_count(
                truth, baseline, genus_relaxed=genus_relaxed
            )
            added_predictions = len(predicted) - len(baseline)
            added_tp = matched - baseline_matched
            added_fp = added_predictions - added_tp
            if added_predictions < 0 or added_tp < 0 or added_fp < 0:
                raise EvaluationError(
                    f"Patient {patient_id} arm {arm} produced invalid incremental counts"
                )
            arm_rows[arm].append(
                {
                    "patient_id": patient_id,
                    "artifact_missing": artifact is None,
                    "truth_count": len(truth),
                    "candidate_count": candidate_count,
                    "predicted_count": len(predicted),
                    "manual_review_count": len(manual),
                    "unknown_count": len(unknown),
                    "tp": matched,
                    "fp": len(predicted) - matched,
                    "fn": len(truth) - matched,
                    "rescue_added_tp": added_tp,
                    "rescue_added_fp": added_fp,
                    "no_path_case": not truth,
                    "no_path_correct_negative": bool(
                        not truth and artifact is not None and not predicted
                    ),
                    "no_path_predicted_positive": bool(
                        not truth and artifact is not None and predicted
                    ),
                    "predicted_pathogens": list(predicted),
                    "manual_review_pathogens": list(manual),
                    "unknown_pathogens": list(unknown),
                }
            )

    output_arms: dict[str, Any] = {}
    for arm, rows in arm_rows.items():
        tp = sum(row["tp"] for row in rows)
        fp = sum(row["fp"] for row in rows)
        fn = sum(row["fn"] for row in rows)
        precision = _safe_div(tp, tp + fp)
        recall = _safe_div(tp, tp + fn)
        added_tp = sum(row["rescue_added_tp"] for row in rows)
        added_fp = sum(row["rescue_added_fp"] for row in rows)
        no_path_total = sum(row["no_path_case"] for row in rows)
        no_path_correct = sum(row["no_path_correct_negative"] for row in rows)
        no_path_predicted_positive = sum(
            row["no_path_predicted_positive"] for row in rows
        )
        no_path_missing = sum(
            row["no_path_case"] and row["artifact_missing"] for row in rows
        )
        output_arms[arm] = {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "predicted_count": tp + fp,
            "truth_count": tp + fn,
            "precision": precision,
            "recall": recall,
            "f0_5": _fbeta(precision, recall, 0.5),
            "f1": _fbeta(precision, recall, 1.0),
            "rescue_from_baseline": {
                "added_tp": added_tp,
                "added_fp": added_fp,
                "added_predictions": added_tp + added_fp,
                "ppv": _safe_div(added_tp, added_tp + added_fp),
            },
            "no_path_specificity": _safe_div(no_path_correct, no_path_total),
            "no_path_cases": no_path_total,
            "no_path_correct_negative_cases": no_path_correct,
            "no_path_predicted_positive_cases": no_path_predicted_positive,
            "no_path_missing_artifact_cases": no_path_missing,
            "case_counts": {
                "included": len(rows),
                "artifact_present": sum(not row["artifact_missing"] for row in rows),
                "artifact_missing": sum(row["artifact_missing"] for row in rows),
                "answer_positive": sum(not row["no_path_case"] for row in rows),
                "no_path": no_path_total,
                "with_predictions": sum(row["predicted_count"] > 0 for row in rows),
                "with_manual_review": sum(
                    row["manual_review_count"] > 0 for row in rows
                ),
                "with_unknown": sum(row["unknown_count"] > 0 for row in rows),
            },
            "manual_review_count": sum(row["manual_review_count"] for row in rows),
            "manual_review_case_count": sum(
                row["manual_review_count"] > 0 for row in rows
            ),
            "unknown_count": sum(row["unknown_count"] for row in rows),
            "unknown_case_count": sum(row["unknown_count"] > 0 for row in rows),
            "cases": rows,
        }
    return {
        "case_count": len(scoped_cases),
        "positive_cases_only": positive_cases_only,
        "genus_relaxed": genus_relaxed,
        "arms": output_arms,
    }


def evaluate_directory(
    outputs_dir: str | Path,
    gold_path: str | Path,
    arm_names: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Evaluate complete clinical artifacts against frozen clinical gold.

    Missing artifacts for labeled patients are retained as empty end-to-end
    predictions, so every positive gold pathogen becomes a false negative.  A
    missing no-path artifact is reported and conservatively is not credited as a
    correct negative in no-path specificity.
    """

    outputs_path = Path(outputs_dir)
    gold_file = Path(gold_path)
    if not outputs_path.is_dir():
        raise EvaluationError(f"outputs_dir is not a directory: {outputs_path}")
    if not gold_file.is_file():
        raise EvaluationError(f"gold_path is not a file: {gold_file}")
    selected_arms = _resolve_arm_names(arm_names)
    gold, aliases = _load_gold(gold_file)

    artifacts: dict[str, ArtifactView] = {}
    artifact_hashes: dict[str, str] = {}
    manifest_paths: list[str] = []
    for path in sorted(outputs_path.glob("*.json")):
        payload = _load_json(path)
        if (
            isinstance(payload, dict)
            and payload.get("schema_version") == MANIFEST_SCHEMA_VERSION
        ):
            manifest_paths.append(str(path.resolve()))
            continue
        artifact = _validate_artifact(
            payload, path=path, aliases=aliases, arm_names=selected_arms
        )
        if artifact.patient_id in artifacts:
            raise EvaluationError(
                f"Duplicate clinical artifact for patient {artifact.patient_id}"
            )
        artifacts[artifact.patient_id] = artifact
        artifact_hashes[artifact.patient_id] = _sha256(path)

    cases = gold["cases"]
    labeled_ids = {str(case["patient_id"]) for case in cases}
    artifact_ids = set(artifacts)
    missing_ids = sorted(labeled_ids - artifact_ids)
    extra_ids = sorted(artifact_ids - labeled_ids)
    primary = _evaluate_scope(
        cases=cases,
        artifacts=artifacts,
        arm_names=selected_arms,
        positive_cases_only=False,
        genus_relaxed=False,
    )
    sensitivity = _evaluate_scope(
        cases=cases,
        artifacts=artifacts,
        arm_names=selected_arms,
        positive_cases_only=True,
        genus_relaxed=True,
    )
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "primary_scope": "all_labeled_exact_frozen_synonyms",
        "sensitivity_scope": "answer_positive_only_genus_relaxed",
        "arm_names": list(selected_arms),
        "primary_arm": "CL2_CONVERGENT" if "CL2_CONVERGENT" in selected_arms else None,
        "gold": {
            key: value for key, value in gold.items() if key != "cases"
        }
        | {"labeled_case_count": len(cases)},
        "artifacts": {
            "outputs_dir": str(outputs_path.resolve()),
            "artifact_count": len(artifacts),
            "labeled_artifact_count": len(artifact_ids & labeled_ids),
            "missing_labeled_artifact_count": len(missing_ids),
            "missing_labeled_patient_ids": missing_ids,
            "extra_unlabeled_artifact_count": len(extra_ids),
            "extra_unlabeled_patient_ids": extra_ids,
            "batch_manifest_count": len(manifest_paths),
            "batch_manifest_paths": manifest_paths,
            "sha256_by_patient": artifact_hashes,
        },
        "scopes": {
            "all_labeled_exact_frozen_synonyms": primary,
            "answer_positive_only_genus_relaxed": sensitivity,
        },
        "notes": [
            "Primary analysis includes every labeled case, exact normalized species, and frozen synonyms.",
            "Sensitivity analysis includes answer-positive cases only and permits one-to-one genus-relaxed matching.",
            "Manual-review and unknown candidates are not positive predictions but are counted separately.",
            "Missing artifacts remain in end-to-end denominators; positive labels become false negatives.",
            "A missing no-path artifact is not credited as a correct negative.",
        ],
    }


def _percent(value: Any) -> str:
    if value is None or not isinstance(value, (int, float)) or not math.isfinite(value):
        return "NA"
    return f"{value * 100:.1f}%"


def render_markdown(report: Mapping[str, Any]) -> str:
    if report.get("schema_version") != REPORT_SCHEMA_VERSION:
        raise EvaluationError("Unsupported clinical evaluation report schema")
    lines = [
        "# RAG_re clinical rescue evaluation",
        "",
        f"- Labeled cases: {report['gold']['labeled_case_count']}",
        f"- Clinical artifacts: {report['artifacts']['artifact_count']}",
        f"- Missing labeled artifacts: {report['artifacts']['missing_labeled_artifact_count']}",
    ]
    missing = report["artifacts"].get("missing_labeled_patient_ids") or []
    if missing:
        lines.append("- Missing patient IDs: " + ", ".join(map(str, missing)))
    for scope_name in (
        "all_labeled_exact_frozen_synonyms",
        "answer_positive_only_genus_relaxed",
    ):
        scope = report["scopes"][scope_name]
        lines.extend(
            [
                "",
                f"## {scope_name}",
                "",
                f"Cases: {scope['case_count']}",
                "",
                "| Arm | TP | FP | FN | Precision | Recall | F0.5 | F1 | Added TP | Added FP | Rescue PPV | No-path specificity | Manual | Unknown |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for arm in report["arm_names"]:
            metrics = scope["arms"][arm]
            rescue = metrics["rescue_from_baseline"]
            lines.append(
                f"| {arm} | {metrics['tp']} | {metrics['fp']} | {metrics['fn']} | "
                f"{_percent(metrics['precision'])} | {_percent(metrics['recall'])} | "
                f"{_percent(metrics['f0_5'])} | {_percent(metrics['f1'])} | "
                f"{rescue['added_tp']} | {rescue['added_fp']} | "
                f"{_percent(rescue['ppv'])} | {_percent(metrics['no_path_specificity'])} | "
                f"{metrics['manual_review_count']} | {metrics['unknown_count']} |"
            )
    lines.extend(["", "## Interpretation notes", ""])
    lines.extend(f"- {note}" for note in report.get("notes", []))
    return "\n".join(lines) + "\n"


def write_markdown(report: Mapping[str, Any], path: str | Path) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_markdown(report), encoding="utf-8")
    return output_path


__all__ = [
    "BASELINE_ARM",
    "FIXED_ARM_NAMES",
    "EvaluationError",
    "evaluate_directory",
    "render_markdown",
    "write_markdown",
]

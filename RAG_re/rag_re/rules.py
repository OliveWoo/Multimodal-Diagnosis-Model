from __future__ import annotations

import re
from typing import Any, Dict, Iterable, Iterator, List, Mapping


TriState = bool | None


def _norm(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", str(value or "").upper()).strip("_")


def _float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _integer(value: Any) -> int | None:
    match = re.search(r"\d+", str(value or ""))
    return int(match.group(0)) if match else None


def _site_alignment(
    candidate: Dict[str, Any],
    config: Mapping[str, Any],
    *,
    context: Mapping[str, Any] | None = None,
    infer_from_direct_evidence: bool = False,
) -> Dict[str, Any]:
    alignment = _norm(candidate.get("specimen_alignment"))
    specimen_class = _norm(candidate.get("specimen_class"))
    allowed_alignments = {
        _norm(value) for value in config.get("positive_specimen_alignments", [])
    }
    negative_alignments = {
        _norm(value) for value in config.get("negative_specimen_alignments", [])
    }
    allowed_classes = {
        _norm(value) for value in config.get("positive_specimen_classes", [])
    }
    if alignment in allowed_alignments:
        return {"aligned": True, "basis": f"specimen_alignment={alignment}"}
    if alignment in negative_alignments:
        return {"aligned": False, "basis": f"specimen_alignment={alignment}"}
    hospital_only_direct_check = (
        infer_from_direct_evidence and alignment == "HOSPITAL_ONLY"
    )
    if specimen_class in allowed_classes and not hospital_only_direct_check:
        return {"aligned": True, "basis": f"specimen_class={specimen_class}"}

    if infer_from_direct_evidence and config.get("allow_direct_evidence_site_inference"):
        site_code = _norm((context or {}).get("site_code"))
        site_values = list(_direct_evidence_site_values(candidate))
        site_patterns = {
            "LOWER_RESPIRATORY": ("LOWER_RESP", "SPUTUM", "BAL", "BALF", "BRONCH"),
            "UPPER_RESPIRATORY": ("UPPER_RESP", "NASAL", "NASOPHARY", "THROAT"),
            "BLOODSTREAM": ("BLOOD", "SERUM", "PLASMA", "SYSTEMIC"),
            "CNS": ("CSF", "CNS", "CEREBROSPINAL"),
            "URINARY": ("URINE", "URINARY"),
            "INTRA_ABDOMINAL": ("ABDOM", "ASCITES", "PERITONE"),
        }
        patterns = site_patterns.get(site_code, ())
        matched = [
            value for value in site_values if any(pattern in value for pattern in patterns)
        ]
        if matched:
            return {
                "aligned": True,
                "basis": f"direct_evidence_site={matched[0]} matches {site_code}",
            }
        if site_values:
            if not patterns:
                return {
                    "aligned": None,
                    "basis": (
                        f"direct evidence sites {site_values} were available, but target "
                        f"site {site_code or 'UNKNOWN'} was not classifiable"
                    ),
                }
            return {
                "aligned": False,
                "basis": (
                    f"direct evidence sites {site_values} did not match "
                    f"target site {site_code or 'UNKNOWN'}"
                ),
            }
    if specimen_class in allowed_classes:
        return {"aligned": True, "basis": f"specimen_class={specimen_class}"}
    return {
        "aligned": None,
        "basis": (
            f"unrecognized alignment={alignment or 'MISSING'}, "
            f"specimen_class={specimen_class or 'MISSING'}"
        ),
    }


def _direct_evidence_site_values(value: Any) -> Iterator[str]:
    if isinstance(value, dict):
        for key, nested in value.items():
            normalized_key = _norm(key)
            if normalized_key in {
                "SPECIMEN_CATEGORY",
                "SPECIMEN_TYPE",
                "SPECIMEN_SITE",
                "COLLECTION_SITE",
            } or normalized_key.endswith(("_SPECIMEN_TYPE", "_SAMPLE_TYPE")):
                normalized = _norm(nested)
                if normalized:
                    yield normalized
            yield from _direct_evidence_site_values(nested)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _direct_evidence_site_values(item)


def tri_or(*values: TriState) -> TriState:
    if any(value is True for value in values):
        return True
    if all(value is False for value in values):
        return False
    return None


def tri_and(*values: TriState) -> TriState:
    if any(value is False for value in values):
        return False
    if all(value is True for value in values):
        return True
    return None


def module_a_from_judgments(
    retrieval: Dict[str, Any],
    judgments: List[Dict[str, Any]],
    config: Mapping[str, Any],
) -> Dict[str, Any]:
    if retrieval.get("status") != "ok":
        return {
            "module": "A_LITERATURE",
            "positive": None,
            "status": retrieval.get("status") or "unavailable",
            "support_count": 0,
            "judgeable_count": 0,
            "support_ratio": None,
            "reason": "PubMed retrieval was unavailable; this is an abstention, not a negative vote.",
        }

    judgeable = [row for row in judgments if row.get("status") == "ok"]
    valid_supports: List[Dict[str, Any]] = []
    rejected_supports: List[Dict[str, Any]] = []
    for row in judgeable:
        if row.get("verdict") != "support":
            continue
        valid = True
        invalid_reasons = []
        if config.get("require_exact_species") and row.get("organism_scope") != "exact_species":
            valid = False
            invalid_reasons.append("not_exact_species")
        if config.get("require_human_clinical_evidence") and row.get(
            "human_clinical_evidence"
        ) is not True:
            valid = False
            invalid_reasons.append("not_human_clinical")
        if config.get("require_target_site_match") and row.get("target_site_match") is not True:
            valid = False
            invalid_reasons.append("site_mismatch")
        if config.get("require_evidence_span") and row.get("evidence_span_valid") is not True:
            valid = False
            invalid_reasons.append("invalid_evidence_span")
        if valid:
            valid_supports.append(row)
        else:
            rejected_supports.append(
                {"pmid": row.get("pmid"), "invalid_reasons": invalid_reasons}
            )

    support_count = len(valid_supports)
    judgeable_count = len(judgeable)
    ratio = support_count / judgeable_count if judgeable_count else 0.0
    min_judgeable = int(config.get("min_judgeable_articles", 5))
    if judgeable_count < min_judgeable:
        positive: TriState = None
        status = "insufficient_articles"
        reason = (
            f"Only {judgeable_count} judgeable articles were available; "
            f"at least {min_judgeable} are required."
        )
    else:
        min_support = int(config.get("min_support", 3))
        min_ratio = float(config.get("min_support_ratio", 0.0))
        positive = support_count >= min_support and ratio >= min_ratio
        status = "ok"
        reason = (
            f"Valid support {support_count}/{judgeable_count}; thresholds are "
            f">={min_support} articles and ratio >={min_ratio:.2f}."
        )

    verdict_counts = {
        label: sum(1 for row in judgeable if row.get("verdict") == label)
        for label in ("support", "against", "irrelevant", "unclear")
    }
    return {
        "module": "A_LITERATURE",
        "positive": positive,
        "status": status,
        "support_count": support_count,
        "judgeable_count": judgeable_count,
        "support_ratio": round(ratio, 6) if judgeable_count else None,
        "verdict_counts": verdict_counts,
        "valid_support_pmids": [str(row.get("pmid") or "") for row in valid_supports],
        "rejected_supports": rejected_supports,
        "reason": reason,
    }


def module_a_skipped() -> Dict[str, Any]:
    return {
        "module": "A_LITERATURE",
        "positive": None,
        "status": "skipped",
        "support_count": 0,
        "judgeable_count": 0,
        "support_ratio": None,
        "valid_support_pmids": [],
        "reason": "Module A was skipped by the caller.",
    }


def module_a_not_requested(reason: str) -> Dict[str, Any]:
    return {
        "module": "A_LITERATURE",
        "positive": None,
        "status": "not_requested",
        "support_count": 0,
        "judgeable_count": 0,
        "support_ratio": None,
        "valid_support_pmids": [],
        "reason": reason,
    }


def module_b_mngs(candidate: Dict[str, Any], config: Mapping[str, Any]) -> Dict[str, Any]:
    alignment = _site_alignment(candidate, config)
    aligned = alignment["aligned"]
    reasons: List[str] = []
    triggers: List[str] = []

    signal = _norm(candidate.get("mngs_signal_tier"))
    allowed_signals = {_norm(value) for value in config.get("positive_signal_tiers", [])}
    if signal in allowed_signals:
        triggers.append(f"mngs_signal_tier={candidate.get('mngs_signal_tier')}")

    reads_tier = _norm(candidate.get("reads_tier"))
    allowed_reads = {_norm(value) for value in config.get("positive_reads_tiers", [])}
    if reads_tier in allowed_reads:
        triggers.append(f"reads_tier={candidate.get('reads_tier')}")

    dominance = _integer(candidate.get("dominance_tier"))
    if dominance is not None and dominance >= int(config.get("min_dominance_tier", 2)):
        triggers.append(f"dominance_tier={candidate.get('dominance_tier')}")

    rank = _integer(candidate.get("rank_priority"))
    percentile = _float(candidate.get("reads_percentile"))
    if (
        rank is not None
        and percentile is not None
        and rank <= int(config.get("max_rank_priority", 1))
        and percentile >= float(config.get("min_reads_percentile", 0.8))
    ):
        triggers.append(f"rank_priority={rank} and reads_percentile={percentile:.4f}")

    raw_positive = bool(triggers)
    status = "ok"
    if config.get("require_site_alignment", True) and raw_positive and aligned is False:
        positive: TriState = False
        reasons.append(f"mNGS strength trigger=True, but {alignment['basis']}.")
    elif config.get("require_site_alignment", True) and raw_positive and aligned is None:
        positive = None
        status = "insufficient_input"
        reasons.append(
            f"mNGS strength trigger=True, but site alignment was unavailable: "
            f"{alignment['basis']}."
        )
    else:
        positive = raw_positive
    if triggers:
        reasons.extend(triggers)
    elif not reasons:
        reasons.append("No pre-specified mNGS strength criterion was met.")
    return {
        "module": "B_MNGS",
        "positive": positive,
        "status": status,
        "site_aligned": aligned,
        "site_alignment_basis": alignment["basis"],
        "triggers": triggers,
        "reason": "; ".join(reasons),
        "features": {
            "mngs_signal_tier": candidate.get("mngs_signal_tier"),
            "rank_priority": candidate.get("rank_priority"),
            "reads": candidate.get("reads"),
            "reads_tier": candidate.get("reads_tier"),
            "reads_percentile": candidate.get("reads_percentile"),
            "dominance_tier": candidate.get("dominance_tier"),
            "specimen_alignment": candidate.get("specimen_alignment"),
        },
    }


NEGATIVE_VALUE_MARKERS = {
    "",
    "NA",
    "N_A",
    "NONE",
    "NULL",
    "NOT_AVAILABLE",
    "UNAVAILABLE",
    "NEGATIVE",
    "NO_SUPPORT",
    "NOT_SUPPORT",
    "FALSE",
    "UNKNOWN",
    "PENDING",
    "INDETERMINATE",
    "NOT_PERFORMED",
    "NOT_DONE",
    "ORDERED",
}

SUPPORTIVE_VALUE_MARKERS = {
    "SUPPORT",
    "SUPPORTED",
    "POSITIVE",
    "DETECTED",
    "PRESENT",
    "REACTIVE",
    "CONFIRMED",
    "TRUE",
}

SUPPORT_STATUS_KEYS = {
    "STATUS",
    "RESULT",
    "VERDICT",
    "SUPPORT",
    "SUPPORTED",
    "POSITIVE",
    "DETECTED",
}


def _value_is_supportive(value: Any) -> bool:
    if value is True:
        return True
    if value in (None, False):
        return False
    if isinstance(value, (list, tuple, set)):
        return any(_value_is_supportive(item) for item in value)
    if isinstance(value, dict):
        # Dates, specimen labels and free text do not prove a positive result.  Only
        # explicit result/status keys may turn a structured object into support.
        return any(
            _value_is_supportive(item)
            for key, item in value.items()
            if _norm(key) in SUPPORT_STATUS_KEYS
        )
    normalized = _norm(value)
    if normalized in NEGATIVE_VALUE_MARKERS:
        return False
    if any(
        marker in normalized
        for marker in ("NOT_AVAILABLE", "NO_SUPPORT", "NEGATIVE", "NOT_DETECTED")
    ):
        return False
    return normalized in SUPPORTIVE_VALUE_MARKERS


def module_c_direct_support(
    candidate: Dict[str, Any],
    config: Mapping[str, Any],
    context: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    allowed = {_norm(value) for value in config.get("positive_modules", [])}
    supports: List[Dict[str, Any]] = []
    summary = candidate.get("module_support_summary")
    if isinstance(summary, dict):
        for module, value in summary.items():
            if _norm(module) in allowed and _value_is_supportive(value):
                supports.append({"module": str(module), "value": value, "source": "summary"})

    key_evidence = candidate.get("key_evidence")
    if isinstance(key_evidence, dict):
        modules = key_evidence.get("support_modules")
        if isinstance(modules, list):
            for module in modules:
                if _norm(module) in allowed:
                    supports.append(
                        {"module": str(module), "value": "listed", "source": "key_evidence"}
                    )

    # Some upstream files store module evidence as direct candidate fields.
    for field, value in candidate.items():
        normalized = _norm(field)
        if normalized in allowed and _value_is_supportive(value):
            supports.append({"module": str(field), "value": value, "source": "candidate_field"})

    deduped: List[Dict[str, Any]] = []
    seen = set()
    for support in supports:
        key = (_norm(support["module"]), str(support["source"]))
        if key not in seen:
            deduped.append(support)
            seen.add(key)
    alignment = _site_alignment(
        candidate,
        config,
        context=context,
        infer_from_direct_evidence=True,
    )
    aligned = alignment["aligned"]
    requires_alignment = bool(config.get("require_site_alignment", True))
    blocked_by_alignment = bool(deduped and requires_alignment and aligned is False)
    unknown_alignment = bool(deduped and requires_alignment and aligned is None)
    positive: TriState = (
        None
        if unknown_alignment
        else bool(deduped) and not blocked_by_alignment
    )
    status = "insufficient_input" if unknown_alignment else "ok"
    return {
        "module": "C_DIRECT_SUPPORT",
        "positive": positive,
        "status": status,
        "site_aligned": aligned,
        "site_alignment_basis": alignment["basis"],
        "supports": deduped,
        "reason": (
            "Direct non-mNGS support was present but site alignment was explicitly negative: "
            + alignment["basis"]
            if blocked_by_alignment
            else "Direct non-mNGS support was present, but site alignment was unavailable: "
            + alignment["basis"]
            if unknown_alignment
            else "Direct non-mNGS support: "
            + ", ".join(str(item["module"]) for item in deduped)
            if deduped
            else "No pre-specified same-organism non-mNGS support was found."
        ),
    }


EXPERIMENT_META: Dict[str, Dict[str, str]] = {
    "E0_BASELINE": {
        "expression": "existing best_available_summary.picked_pathogens",
        "semantics": "clinical_baseline",
    },
    "A_LITERATURE": {"expression": "A", "semantics": "screening_trigger"},
    "B_MNGS": {"expression": "B", "semantics": "screening_trigger"},
    "C_DIRECT_SUPPORT": {"expression": "C", "semantics": "screening_trigger"},
    "AB_OR": {"expression": "A OR B", "semantics": "screening_trigger"},
    "AC_OR": {"expression": "A OR C", "semantics": "screening_trigger"},
    "BC_OR": {"expression": "B OR C", "semantics": "screening_trigger"},
    "ABC_OR": {"expression": "A OR B OR C", "semantics": "screening_trigger"},
    "ABC_MAJORITY": {"expression": "at least 2 of A/B/C", "semantics": "gated_policy"},
    "GATED_C_OR_A_AND_B": {
        "expression": "C OR (A AND B)",
        "semantics": "standalone_gated_policy",
    },
    "E0_OR_A": {
        "expression": "E0 OR A",
        "semantics": "baseline_preserving_cumulative",
    },
    "E0_OR_A_OR_B": {
        "expression": "E0 OR A OR B",
        "semantics": "baseline_preserving_cumulative",
    },
    "E0_OR_A_OR_B_OR_C": {
        "expression": "E0 OR A OR B OR C",
        "semantics": "baseline_preserving_cumulative",
    },
    "E0_OR_C": {
        "expression": "E0 OR C",
        "semantics": "baseline_preserving_cumulative",
    },
    "E0_OR_C_OR_A_AND_B": {
        "expression": "E0 OR C OR (A AND B)",
        "semantics": "baseline_preserving_gated_policy",
    },
    "E0_AND_B": {
        "expression": "E0 AND B",
        "semantics": "precision_pruning_policy",
    },
    "E0_AND_C": {
        "expression": "E0 AND C",
        "semantics": "precision_pruning_policy",
    },
    "E0_AND_B_OR_C": {
        "expression": "E0 AND (B OR C)",
        "semantics": "precision_pruning_policy",
    },
    "E0_AND_B_AND_C": {
        "expression": "E0 AND B AND C",
        "semantics": "precision_pruning_policy",
    },
}


def _majority(values: Iterable[TriState]) -> TriState:
    items = list(values)
    true_count = sum(value is True for value in items)
    unknown_count = sum(value is None for value in items)
    if true_count >= 2:
        return True
    if true_count + unknown_count < 2:
        return False
    return None


def experiment_values(
    *, baseline: bool, a: TriState, b: TriState, c: TriState
) -> Dict[str, TriState]:
    return {
        "E0_BASELINE": baseline,
        "A_LITERATURE": a,
        "B_MNGS": b,
        "C_DIRECT_SUPPORT": c,
        "AB_OR": tri_or(a, b),
        "AC_OR": tri_or(a, c),
        "BC_OR": tri_or(b, c),
        "ABC_OR": tri_or(a, b, c),
        "ABC_MAJORITY": _majority([a, b, c]),
        "GATED_C_OR_A_AND_B": tri_or(c, tri_and(a, b)),
        "E0_OR_A": tri_or(baseline, a),
        "E0_OR_A_OR_B": tri_or(baseline, a, b),
        "E0_OR_A_OR_B_OR_C": tri_or(baseline, a, b, c),
        "E0_OR_C": tri_or(baseline, c),
        "E0_OR_C_OR_A_AND_B": tri_or(baseline, c, tri_and(a, b)),
        "E0_AND_B": tri_and(baseline, b),
        "E0_AND_C": tri_and(baseline, c),
        "E0_AND_B_OR_C": tri_and(baseline, tri_or(b, c)),
        "E0_AND_B_AND_C": tri_and(baseline, b, c),
    }


def build_experiment_summary(candidate_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    output: Dict[str, Any] = {}
    for experiment, meta in EXPERIMENT_META.items():
        positives = [
            row["organism_name"]
            for row in candidate_rows
            if row["experiments"].get(experiment) is True
        ]
        abstentions = [
            row["organism_name"]
            for row in candidate_rows
            if row["experiments"].get(experiment) is None
        ]
        output[experiment] = {
            **meta,
            "predicted_pathogens": positives,
            "predicted_count": len(positives),
            "abstained_pathogens": abstentions,
            "abstained_count": len(abstentions),
        }
    return output

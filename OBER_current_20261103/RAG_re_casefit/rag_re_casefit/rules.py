from __future__ import annotations

from collections import Counter
from typing import Any, Dict, Iterable, List


ALIGNMENTS = {"match", "partial", "mismatch", "unknown"}
ORGANISM_SCOPES = {"exact_species", "genus_only", "other", "unclear"}
HUMAN_CONTEXTS = {
    "causal_infection",
    "colonization_or_contamination",
    "reactivation_or_bystander",
    "detection_only",
    "nonhuman",
    "unclear",
}
EXTRA_DIMENSIONS = (
    "specimen_alignment",
    "host_alignment",
    "phenotype_alignment",
    "temporal_alignment",
)


def _compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def _validate_extraction(extraction: Dict[str, Any]) -> None:
    if extraction.get("organism_scope") not in ORGANISM_SCOPES:
        raise ValueError("invalid organism_scope")
    if extraction.get("human_clinical_context") not in HUMAN_CONTEXTS:
        raise ValueError("invalid human_clinical_context")
    for key in ("site_alignment", "syndrome_alignment", *EXTRA_DIMENSIONS):
        if extraction.get(key) not in ALIGNMENTS:
            raise ValueError(f"invalid {key}")


def derive_article_verdict(extraction: Dict[str, Any]) -> Dict[str, Any]:
    """Derive the four-state verdict without trusting an LLM Yes/No conclusion."""

    _validate_extraction(extraction)
    scope = extraction["organism_scope"]
    human = extraction["human_clinical_context"]
    site = extraction["site_alignment"]
    syndrome = extraction["syndrome_alignment"]
    reasons: List[str] = []
    counterevidence = False
    retrieval_mismatch = False

    if scope in {"genus_only", "other"}:
        reasons.append("organism_not_exact_species")
        retrieval_mismatch = True
    elif human == "nonhuman":
        reasons.append("nonhuman_evidence")
        retrieval_mismatch = True
    elif site == "mismatch" or syndrome == "mismatch":
        reasons.append("clinical_site_or_syndrome_mismatch")
        retrieval_mismatch = True
    elif human in {
        "colonization_or_contamination",
        "reactivation_or_bystander",
        "detection_only",
    }:
        reasons.append(f"noncausal_context:{human}")
        # Absence of causal attribution is neutral non-support, not evidence against
        # causality.  Colonization/contamination is explicit counterevidence at a
        # compatible site.  Reactivation/bystander only counts against the target
        # syndrome when the paper is explicitly about that same site and syndrome.
        if human == "colonization_or_contamination":
            counterevidence = scope == "exact_species" and site in {"match", "partial"}
        elif human == "reactivation_or_bystander":
            counterevidence = (
                scope == "exact_species" and site == "match" and syndrome == "match"
            )

    if reasons:
        verdict = "mismatch"
    else:
        hard_known = (
            scope == "exact_species"
            and human == "causal_infection"
            and site in {"match", "partial"}
            and syndrome in {"match", "partial"}
        )
        if not hard_known:
            verdict = "insufficient"
            reasons.append("hard_gate_information_incomplete")
        else:
            extra_values = [extraction[key] for key in EXTRA_DIMENSIONS]
            extra_match_count = sum(value == "match" for value in extra_values)
            extra_mismatch_count = sum(value == "mismatch" for value in extra_values)
            if (
                site == "match"
                and syndrome == "match"
                and extra_match_count >= 1
                and extra_mismatch_count == 0
            ):
                verdict = "strong_match"
                reasons.append("hard_gates_and_patient_specific_dimension_match")
            else:
                verdict = "partial_match"
                reasons.append("hard_gates_match_but_patient_specific_fit_is_incomplete")

    return {
        "verdict": verdict,
        "counterevidence": counterevidence,
        "retrieval_mismatch": retrieval_mismatch,
        "decision_reasons": reasons,
    }


def apply_traceability_gate(
    decision: Dict[str, Any],
    *,
    support_span_valid: bool,
    counter_span_valid: bool,
) -> Dict[str, Any]:
    """Require exact article spans before evidence enters candidate aggregation."""

    result = dict(decision)
    reasons = list(result.get("decision_reasons") or [])
    verdict = result.get("verdict")
    if verdict in {"strong_match", "partial_match"} and not support_span_valid:
        result["verdict"] = "insufficient"
        result["counterevidence"] = False
        result["retrieval_mismatch"] = False
        reasons.append("invalid_or_missing_support_span")
    if result.get("counterevidence") is True and not counter_span_valid:
        result["counterevidence"] = False
        reasons.append("invalid_or_missing_counter_span")
    result["decision_reasons"] = reasons
    result["traceability_status"] = (
        "valid"
        if not any(reason.startswith("invalid_or_missing_") for reason in reasons)
        else "invalid_required_span"
    )
    return result


def aggregate_judgments(
    judgments: Iterable[Dict[str, Any]],
    config: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    config = dict(config or {})
    min_judgeable = int(config.get("min_judgeable_articles", 5))
    min_strong_high = int(config.get("min_strong_for_high", 2))
    min_strong_moderate = int(config.get("min_strong_for_moderate", 1))
    min_partial_moderate = int(config.get("min_partial_for_moderate", 3))
    if min_judgeable < 1 or min_strong_high < 1 or min_strong_moderate < 1 or min_partial_moderate < 1:
        raise ValueError("aggregation thresholds must be positive integers")

    rows = list(judgments)
    pmids = [_compact(row.get("pmid")) for row in rows]
    nonempty_pmids = [pmid for pmid in pmids if pmid]
    if len(nonempty_pmids) != len(set(nonempty_pmids)):
        raise ValueError("duplicate PMIDs cannot be counted as independent articles")

    verdicts = Counter()
    counter_count = 0
    retrieval_mismatch_count = 0
    neutral_non_support_count = 0
    for row in rows:
        verdict = row.get("verdict")
        if verdict not in {"strong_match", "partial_match", "mismatch", "insufficient"}:
            raise ValueError("judgment has invalid verdict")
        verdicts[verdict] += 1
        counter_count += row.get("counterevidence") is True
        retrieval_mismatch_count += row.get("retrieval_mismatch") is True
        neutral_non_support_count += (
            verdict == "mismatch"
            and row.get("counterevidence") is not True
            and row.get("retrieval_mismatch") is not True
        )

    support_count = verdicts["strong_match"] + verdicts["partial_match"]
    judgeable_count = support_count + counter_count
    if judgeable_count < min_judgeable:
        tier = "insufficient"
    elif counter_count and support_count:
        tier = "conflicted"
    elif counter_count:
        tier = "low"
    elif verdicts["strong_match"] >= min_strong_high:
        tier = "high"
    elif (
        verdicts["strong_match"] >= min_strong_moderate
        or verdicts["partial_match"] >= min_partial_moderate
    ):
        tier = "moderate"
    else:
        tier = "low"

    return {
        "schema_version": "rag_re_casefit.aggregate.v1",
        "article_count": len(rows),
        "strong_match_count": verdicts["strong_match"],
        "partial_match_count": verdicts["partial_match"],
        "counterevidence_count": counter_count,
        "retrieval_mismatch_count": retrieval_mismatch_count,
        "neutral_non_support_count": neutral_non_support_count,
        "insufficient_count": verdicts["insufficient"],
        "judgeable_count": judgeable_count,
        "casefit_evidence_tier": tier,
        "auto_rescue": False,
        "interpretation": "Literature case similarity only; not patient-level causal attribution.",
        "thresholds": {
            "min_judgeable_articles": min_judgeable,
            "min_strong_for_high": min_strong_high,
            "min_strong_for_moderate": min_strong_moderate,
            "min_partial_for_moderate": min_partial_moderate,
        },
    }

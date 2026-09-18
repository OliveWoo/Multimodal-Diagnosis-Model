from __future__ import annotations

import copy
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Tuple

from .io_utils import load_json


WRAPPER_KEYS = (
    "agent_bundle",
    "infection_bundle",
    "deterministic_max",
    "bundle",
    "result",
)

MERGE_REVIEW_TIERS = (
    "review_high_priority",
    "review_context_needed",
)
MERGE_EXCLUDED_TIERS = (
    "review_low_specificity",
    "omitted_with_reason",
    "review_omitted_with_reason",
)


def normalize_name(value: Any) -> str:
    return " ".join(str(value or "").strip().lower().split())


def _find_bundles(value: Any, depth: int = 0) -> List[Dict[str, Any]]:
    if depth > 4 or not isinstance(value, dict):
        return []
    found: List[Dict[str, Any]] = []
    if isinstance(value.get("pathogen_candidates"), list):
        found.append(value)
    for key in WRAPPER_KEYS:
        nested = value.get(key)
        found.extend(_find_bundles(nested, depth + 1))
    unique: List[Dict[str, Any]] = []
    seen_ids = set()
    for bundle in found:
        if id(bundle) not in seen_ids:
            unique.append(bundle)
            seen_ids.add(id(bundle))
    return unique


def load_bundle(path: Path) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    data = load_json(path)
    if not isinstance(data, dict):
        raise ValueError(f"Input JSON root must be an object: {path}")
    bundles = _find_bundles(data)
    if not bundles:
        raise ValueError(
            f"Could not find pathogen_candidates in {path}. "
            "RAG_re accepts a direct infection bundle or an agent_bundle wrapper."
        )
    if len(bundles) > 1:
        raise ValueError(
            f"Input contains {len(bundles)} ambiguous infection bundles with pathogen_candidates: {path}"
        )
    return data, bundles[0]


def patient_id_from(path: Path, root: Dict[str, Any], bundle: Dict[str, Any]) -> str:
    def wrapper_sources(value: Any, depth: int = 0) -> Iterable[Dict[str, Any]]:
        if depth > 4 or not isinstance(value, dict):
            return
        yield value
        for key in WRAPPER_KEYS:
            yield from wrapper_sources(value.get(key), depth + 1)

    sources = [bundle]
    sources.extend(source for source in wrapper_sources(root) if source is not bundle)
    explicit_ids = {
        str(source[key]).strip()
        for source in sources
        for key in ("patient_id", "case_id")
        if source.get(key) not in (None, "") and str(source[key]).strip()
    }
    if len(explicit_ids) > 1:
        raise ValueError(f"Input contains conflicting patient/case IDs: {sorted(explicit_ids)}")
    if explicit_ids:
        return next(iter(explicit_ids))
    fallback_ids = {
        str(source["id"]).strip()
        for source in sources
        if source.get("id") not in (None, "") and str(source["id"]).strip()
    }
    if len(fallback_ids) > 1:
        raise ValueError(f"Input contains conflicting fallback IDs: {sorted(fallback_ids)}")
    if fallback_ids:
        return next(iter(fallback_ids))
    match = re.search(r"(?:NGS_)?patient[_ -]?(\d+)", path.name, re.IGNORECASE)
    return match.group(1) if match else path.stem


def candidate_name(candidate: Dict[str, Any]) -> str:
    for key in ("organism_name", "pathogen", "pathogen_name", "name", "organism"):
        value = candidate.get(key)
        if str(value or "").strip():
            return " ".join(str(value).strip().split())
    return ""


def candidate_level(candidate: Dict[str, Any]) -> int | None:
    for key in (
        "integrated_causative_level",
        "recommended_level",
        "basis_level",
        "level",
        "original_level",
    ):
        value = candidate.get(key)
        match = re.search(r"([1-5])", str(value or ""))
        if match:
            return int(match.group(1))
    return None


def _picked_rows(bundle: Dict[str, Any]) -> Iterable[Any]:
    summary = bundle.get("best_available_summary")
    if isinstance(summary, dict) and isinstance(summary.get("picked_pathogens"), list):
        yield from summary["picked_pathogens"]
    if isinstance(bundle.get("picked_pathogens"), list):
        yield from bundle["picked_pathogens"]


def baseline_picked_names(bundle: Dict[str, Any]) -> List[str]:
    names: List[str] = []
    seen = set()
    for row in _picked_rows(bundle):
        name = candidate_name(row) if isinstance(row, dict) else " ".join(str(row).split())
        key = normalize_name(name)
        if key and key not in seen:
            names.append(name)
            seen.add(key)
    return names


def source_context(bundle: Dict[str, Any], candidate: Dict[str, Any]) -> Dict[str, str]:
    specimen_class = str(candidate.get("specimen_class") or "").strip()
    specimen_alignment = str(candidate.get("specimen_alignment") or "").strip()
    dominant_source = str(bundle.get("dominant_source") or "").strip()
    explicit_clinical_site = " ".join(
        str(bundle.get(key) or "").strip()
        for key in ("clinical_target_site", "infection_site", "target_site")
        if str(bundle.get(key) or "").strip()
    )

    def infer_site(text: str) -> str | None:
        normalized = text.lower()
        if any(term in normalized for term in ("upper_resp", "nasopharyn", "nasal", "throat")):
            return "upper_respiratory"
        if any(term in normalized for term in ("blood", "systemic", "sterile", "sepsis", "bacteremia", "s1_")):
            return "bloodstream"
        if any(term in normalized for term in ("lower_resp", "respiratory", "pulmonary", "lung", "bal", "sputum", "s2_")):
            return "lower_respiratory"
        if any(term in normalized for term in ("csf", "cns", "mening", "brain")):
            return "cns"
        if any(term in normalized for term in ("urinary", "urine", "uti")):
            return "urinary"
        if any(term in normalized for term in ("abdominal", "intra_abdominal", "ascites", "peritone")):
            return "intra_abdominal"
        return None

    def evidence_site_text(value: Any) -> str:
        values: List[str] = []
        if isinstance(value, dict):
            for key, nested in value.items():
                normalized_key = normalize_name(key).replace(" ", "_")
                if normalized_key in {
                    "specimen_category",
                    "specimen_type",
                    "specimen_site",
                    "collection_site",
                } or normalized_key.endswith(("_specimen_type", "_sample_type")):
                    values.append(str(nested))
                values.append(evidence_site_text(nested))
        elif isinstance(value, list):
            values.extend(evidence_site_text(item) for item in value)
        return " ".join(value for value in values if value)

    evidence_site = infer_site(evidence_site_text(candidate.get("key_evidence")))
    specimen_site = evidence_site or infer_site(specimen_class)
    # Module A asks whether the organism is supported at the clinical disease site,
    # not merely at the mNGS collection site.  Explicit clinical site/dominant source
    # therefore wins; candidate specimen information is only a fallback.
    site = (
        infer_site(explicit_clinical_site)
        or infer_site(dominant_source)
        or specimen_site
    )
    if site == "lower_respiratory":
        target_site = "lower respiratory tract infection or pneumonia"
        query_terms = "pneumonia OR pulmonary infection OR lower respiratory tract infection"
    elif site == "upper_respiratory":
        target_site = "upper respiratory tract infection"
        query_terms = "upper respiratory tract infection"
    elif site == "bloodstream":
        target_site = "bloodstream infection, bacteremia, or sepsis"
        query_terms = "bloodstream infection OR bacteremia OR sepsis"
    elif site == "cns":
        target_site = "central nervous system infection or meningitis"
        query_terms = "central nervous system infection OR meningitis"
    elif site == "urinary":
        target_site = "urinary tract infection"
        query_terms = "urinary tract infection"
    elif site == "intra_abdominal":
        target_site = "intra-abdominal infection"
        query_terms = "intra-abdominal infection"
    else:
        target_site = "human infection at the reported specimen or clinical site"
        query_terms = "human infection"

    return {
        "target_site": target_site,
        "query_terms": query_terms,
        "specimen_class": specimen_class,
        "specimen_alignment": specimen_alignment,
        "dominant_source": dominant_source,
        "site_code": site or "unknown",
        "clinical_site_source": (
            "explicit_bundle_field"
            if infer_site(explicit_clinical_site)
            else "dominant_source"
            if infer_site(dominant_source)
            else "candidate_specimen_fallback"
            if specimen_site
            else "unknown"
        ),
        "specimen_site_code": specimen_site or "unknown",
    }


def candidate_provenance(candidate: Dict[str, Any]) -> str:
    explicit = str(candidate.get("rag_re_candidate_provenance") or "").strip()
    if explicit:
        return explicit
    key_evidence = candidate.get("key_evidence")
    if isinstance(key_evidence, dict) and key_evidence.get("manual_review_expansion") is True:
        return "prior_llm_missed_candidate_expansion"
    applied = candidate.get("applied_rules")
    applied_values = applied if isinstance(applied, list) else [applied]
    marker_values = [
        candidate.get("source_category"),
        candidate.get("rank_rule"),
        *applied_values,
    ]
    marker = " ".join(str(value or "") for value in marker_values).upper()
    if "LLM_MISSED_CANDIDATE" in marker or "MISSED_REVIEW" in marker:
        return "prior_llm_missed_candidate_expansion"
    return "original_candidate_pool"


def candidate_review_tier(candidate: Dict[str, Any]) -> str | None:
    value = str(candidate.get("rag_re_review_tier") or "").strip()
    return value or None


def candidate_literature_eligible(candidate: Dict[str, Any]) -> bool:
    value = candidate.get("rag_re_literature_eligible")
    return True if value is None else value is True


def canonical_organism_name(name: str, config: Dict[str, Any]) -> str:
    aliases = config.get("organism_normalization", {}).get("accepted_synonyms", {})
    normalized_aliases = {normalize_name(key): str(value).strip() for key, value in aliases.items()}
    return normalized_aliases.get(normalize_name(name), " ".join(str(name).split()))


def iter_candidates(bundle: Dict[str, Any]) -> Iterable[Dict[str, Any]]:
    values = bundle.get("pathogen_candidates")
    if not isinstance(values, list):
        return
    for candidate in values:
        if isinstance(candidate, dict) and candidate_name(candidate):
            yield candidate


def _merge_key(value: Any) -> str:
    if isinstance(value, dict):
        name = candidate_name(value)
        if name:
            return re.sub(r"[^a-z0-9]+", "", name.lower())
        canonical_key = str(value.get("canonical_key") or "").strip()
        if canonical_key:
            return re.sub(r"[^a-z0-9]+", "", canonical_key.lower())
    else:
        name = str(value or "")
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def _candidate_from_review(row: Dict[str, Any], tier: str) -> Dict[str, Any]:
    """Create a rule-only candidate without copying LLM rationale into Module A."""

    snapshot = row.get("evidence_snapshot")
    snapshot = snapshot if isinstance(snapshot, dict) else {}
    candidate: Dict[str, Any] = {
        "organism_name": candidate_name(row),
        "classification": row.get("classification"),
        "integrated_causative_level": row.get("integrated_causative_level")
        or row.get("observed_level"),
    }
    for key in (
        "mngs_signal_tier",
        "rank_priority",
        "reads",
        "reads_tier",
        "reads_percentile",
        "dominance_tier",
        "specimen_alignment",
        "specimen_class",
        "source_category",
        "rank_rule",
    ):
        value = snapshot.get(key, row.get(key))
        if value not in (None, ""):
            candidate[key] = value

    support_modules = snapshot.get("support_modules")
    if not isinstance(support_modules, list):
        support_modules = []
    candidate["module_support_summary"] = {
        str(module): "Support" for module in support_modules if str(module).strip()
    }
    key_evidence = {"support_modules": list(support_modules)}
    for key in (
        "non_host_support_modules",
        "related_representative_hospital_support",
        "candida_invasive_hospital_support",
        "candida_related_invasive_hospital_context",
    ):
        if key in snapshot:
            key_evidence[key] = snapshot.get(key)
    candidate["key_evidence"] = key_evidence
    candidate["rag_re_candidate_provenance"] = f"merge_{tier}"
    candidate["rag_re_review_tier"] = tier
    candidate["rag_re_literature_eligible"] = True
    return {key: value for key, value in candidate.items() if value is not None}


def select_candidate_pool(
    root: Dict[str, Any],
    bundle: Dict[str, Any],
    config: Mapping[str, Any] | None = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Select the frozen experimental candidate universe for legacy or merge inputs.

    Formal merge files contain a deterministic universe plus an LLM-created review
    tiering envelope.  Per READ/README_RAG_merge_output_zh.md, baseline picked rows
    are retained, only high/context rows are literature-eligible, and low/omitted
    rows are excluded.  Legacy infection bundles keep their previous behavior.
    """

    raw_candidates = list(iter_candidates(bundle))
    review = root.get("llm_missed_candidate_review")
    is_formal_merge = (
        root.get("deterministic_max") is bundle
        and isinstance(review, dict)
        and any(tier in review for tier in (*MERGE_REVIEW_TIERS, *MERGE_EXCLUDED_TIERS))
    )
    if not is_formal_merge:
        return raw_candidates, {
            "input_candidate_mode": "legacy_pathogen_candidates",
            "source_deterministic_candidate_count": len(raw_candidates),
            "candidate_pool_provenance": "legacy_pathogen_candidates",
            "literature_eligible_candidate_count": sum(
                candidate_literature_eligible(row) for row in raw_candidates
            ),
        }

    pool_config = dict(config or {})
    requested_tiers = pool_config.get("merge_review_tiers", list(MERGE_REVIEW_TIERS))
    if not isinstance(requested_tiers, list) or not requested_tiers:
        raise ValueError("candidate_pool.merge_review_tiers must be a non-empty list")
    unknown_tiers = sorted(set(map(str, requested_tiers)) - set(MERGE_REVIEW_TIERS))
    if unknown_tiers:
        raise ValueError(
            "candidate_pool.merge_review_tiers contains unsupported primary tiers: "
            + ", ".join(unknown_tiers)
        )

    deterministic_by_key: Dict[str, Dict[str, Any]] = {}
    for candidate in raw_candidates:
        key = _merge_key(candidate)
        if key in deterministic_by_key:
            raise ValueError(
                "Formal merge contains duplicate deterministic candidate key: "
                + candidate_name(candidate)
            )
        deterministic_by_key[key] = candidate

    selected: List[Dict[str, Any]] = []
    selected_by_key: Dict[str, Dict[str, Any]] = {}
    synthetic_names: List[str] = []

    def add_candidate(
        source: Any,
        *,
        provenance: str,
        literature_eligible: bool,
        tier: str | None = None,
    ) -> None:
        key = _merge_key(source)
        if not key:
            raise ValueError(f"Formal merge {provenance} row has no organism name")
        existing = selected_by_key.get(key)
        if existing is not None:
            if tier and not existing.get("rag_re_review_tier"):
                existing["rag_re_review_tier"] = tier
            return
        deterministic = deterministic_by_key.get(key)
        if deterministic is not None:
            candidate = copy.deepcopy(deterministic)
        elif isinstance(source, dict) and tier:
            candidate = _candidate_from_review(source, tier)
            synthetic_names.append(candidate_name(candidate))
        elif isinstance(source, dict):
            candidate = copy.deepcopy(source)
        else:
            candidate = {"organism_name": " ".join(str(source).split())}
        candidate["rag_re_candidate_provenance"] = provenance
        candidate["rag_re_literature_eligible"] = literature_eligible
        if tier:
            candidate["rag_re_review_tier"] = tier
        selected.append(candidate)
        selected_by_key[key] = candidate

    picked_rows = list(_picked_rows(bundle))
    literature_on_picked = bool(pool_config.get("literature_on_baseline_picked", False))
    for row in picked_rows:
        add_candidate(
            row,
            provenance="merge_baseline_picked",
            literature_eligible=literature_on_picked,
        )

    selected_tier_counts: Dict[str, int] = {}
    for tier in map(str, requested_tiers):
        values = review.get(tier) or []
        if not isinstance(values, list):
            raise ValueError(f"llm_missed_candidate_review.{tier} must be a list")
        selected_tier_counts[tier] = len(values)
        for row in values:
            if not isinstance(row, dict):
                raise ValueError(f"llm_missed_candidate_review.{tier} rows must be objects")
            add_candidate(
                row,
                provenance=f"merge_{tier}",
                literature_eligible=True,
                tier=tier,
            )

    excluded_tier_counts: Dict[str, int] = {}
    excluded_tier_names: Dict[str, List[str]] = {}
    for tier in MERGE_EXCLUDED_TIERS:
        values = review.get(tier) or []
        if not isinstance(values, list):
            raise ValueError(f"llm_missed_candidate_review.{tier} must be a list")
        excluded_tier_counts[tier] = len(values)
        excluded_tier_names[tier] = [
            candidate_name(row) for row in values if isinstance(row, dict) and candidate_name(row)
        ]

    selected_source_keys = set(selected_by_key)
    unselected_deterministic = [
        candidate_name(row)
        for key, row in deterministic_by_key.items()
        if key not in selected_source_keys
    ]
    tier_slug = "+".join(map(str, requested_tiers))
    return selected, {
        "input_candidate_mode": "formal_merge_review_tiered",
        "source_deterministic_candidate_count": len(raw_candidates),
        "source_baseline_picked_count": len(picked_rows),
        "selected_review_tier_counts": selected_tier_counts,
        "excluded_review_tier_counts": excluded_tier_counts,
        "excluded_review_tier_names": excluded_tier_names,
        "unselected_deterministic_candidate_count": len(unselected_deterministic),
        "unselected_deterministic_candidate_names": unselected_deterministic,
        "synthetic_review_candidate_count": len(synthetic_names),
        "synthetic_review_candidate_names": synthetic_names,
        "literature_eligible_candidate_count": sum(
            candidate_literature_eligible(row) for row in selected
        ),
        "literature_candidate_scope": (
            "all_selected_candidates"
            if literature_on_picked
            else "review_high_priority_and_context_only"
        ),
        "candidate_pool_provenance": (
            "merge_baseline_picked_plus_" + tier_slug
        ),
    }

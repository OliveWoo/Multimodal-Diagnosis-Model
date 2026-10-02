"""Central candidate evidence profiling for scorer, review, merge, and RAG.

The module is answer-blind. It converts reusable evidence axes into stable
negative-evidence flags and owns Candida evidence-strength interpretation.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Sequence

from tools import pathogen_normalization as pathogen_names


REPO_ROOT = Path(__file__).resolve().parent.parent
RULES_PATH = REPO_ROOT / "rules" / "negative_evidence_rules.json"
REVIEW_CONTEXT_RESCUE_RULES_PATH = REPO_ROOT / "rules" / "review_context_rescue_rules.json"


@lru_cache(maxsize=1)
def rules_payload() -> dict[str, Any]:
    payload = json.loads(RULES_PATH.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected object in {RULES_PATH}")
    return payload


@lru_cache(maxsize=1)
def review_context_rescue_rules_payload() -> dict[str, Any]:
    payload = json.loads(REVIEW_CONTEXT_RESCUE_RULES_PATH.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected object in {REVIEW_CONTEXT_RESCUE_RULES_PATH}")
    return payload


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def lower_text(value: Any) -> str:
    return str(value or "").strip().lower()


def _candida_rules() -> dict[str, Any]:
    return as_dict(rules_payload().get("candida_evidence"))


def _negative_rules() -> dict[str, Any]:
    return as_dict(rules_payload().get("negative_evidence"))


PROFILE_VERSION = str(rules_payload().get("profile_version") or "candidate_evidence_profile_v1")
RULE_SCHEMA_VERSION = str(rules_payload().get("schema_version") or "negative_evidence_rules_v1")
DIRECT_HOSPITAL_MODULES = tuple(
    str(value) for value in as_list(rules_payload().get("direct_hospital_modules")) if value
)


def _atypical_respiratory_serology_rule() -> dict[str, Any]:
    rules = as_dict(review_context_rescue_rules_payload().get("rules"))
    return as_dict(rules.get("atypical_respiratory_indirect_serology"))


def _recursive_lower_text(value: Any) -> str:
    if isinstance(value, dict):
        return " ".join(_recursive_lower_text(item) for item in value.values())
    if isinstance(value, list):
        return " ".join(_recursive_lower_text(item) for item in value)
    return lower_text(value)


def atypical_respiratory_serology_profile(
    organism_name: Any,
    hospital_evidence: Sequence[dict[str, Any]] = (),
) -> dict[str, Any]:
    """Profile indirect serology without treating it as direct pulmonary proof."""
    rule = _atypical_respiratory_serology_rule()
    key = pathogen_names.canonical_key(organism_name)
    allowed = {str(value) for value in as_list(rule.get("organism_keys")) if value}
    indirect_types = {lower_text(value) for value in as_list(rule.get("indirect_evidence_types")) if value}
    indirect_terms = tuple(lower_text(value) for value in as_list(rule.get("indirect_test_terms")) if value)
    direct_terms = tuple(lower_text(value) for value in as_list(rule.get("direct_test_terms")) if value)
    indirect_details: list[dict[str, Any]] = []
    direct_details: list[dict[str, Any]] = []

    for hospital_item in hospital_evidence:
        if not isinstance(hospital_item, dict):
            continue
        module_evidence = as_dict(hospital_item.get("module_evidence"))
        for module, rows in module_evidence.items():
            for row in as_list(rows):
                if not isinstance(row, dict):
                    continue
                row_text = f" {_recursive_lower_text(row)} "
                evidence_type = lower_text(row.get("evidence_type"))
                if evidence_type in indirect_types or any(term in row_text for term in indirect_terms):
                    indirect_details.append({"module": module, "evidence_type": evidence_type, "row": row})
                if any(term in row_text for term in direct_terms):
                    direct_details.append({"module": module, "evidence_type": evidence_type, "row": row})

    enabled = rule.get("enabled") is not False
    recognized = key in allowed
    indirect_found = bool(indirect_details)
    direct_found = bool(direct_details)
    return {
        "rule_id": str(rule.get("rule_id") or "R-CTX-ATYPICAL-RESP-INDIRECT-SEROLOGY"),
        "organism_key": key,
        "recognized_atypical_respiratory_pathogen": recognized,
        "indirect_serology_found": indirect_found,
        "direct_molecular_found": direct_found,
        "context_rescue_match": enabled and recognized and indirect_found and not direct_found,
        "target_tier": str(rule.get("target_tier") or "review_context_needed"),
        "maximum_tier": str(rule.get("maximum_tier") or "review_context_needed"),
        "reason": str(rule.get("reason") or "indirect serology requires case-level respiratory review"),
        "indirect_evidence": indirect_details,
        "direct_evidence": direct_details,
    }


def has_compatible_pulmonary_context(value: Any) -> bool:
    rule = _atypical_respiratory_serology_rule()
    terms = tuple(
        lower_text(term)
        for term in as_list(rule.get("compatible_pulmonary_context_terms"))
        if term
    )
    haystack = f" {_recursive_lower_text(value)} "
    return any(term in haystack for term in terms)

_STRENGTHS = as_dict(_candida_rules().get("strengths"))
CANDIDA_EVIDENCE_STRONG_INVASIVE = str(_STRENGTHS.get("strong") or "strong_invasive")
CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC = str(
    _STRENGTHS.get("intermediate") or "intermediate_systemic_or_invasive_clue"
)
CANDIDA_EVIDENCE_WEAK_NONINVASIVE = str(
    _STRENGTHS.get("weak") or "weak_noninvasive_or_colonization"
)

_STRONG_RULE = as_dict(_candida_rules().get("strong_invasive"))
_INTERMEDIATE_RULE = as_dict(_candida_rules().get("intermediate_systemic_or_invasive_clue"))
_WEAK_RULE = as_dict(_candida_rules().get("weak_noninvasive_or_colonization"))
INVASIVE_CANDIDA_SPECIMEN_CATEGORIES = {
    pathogen_names.raw_key(value) for value in as_list(_STRONG_RULE.get("specimen_category_keys"))
}
INVASIVE_CANDIDA_SPECIMEN_TERMS = tuple(lower_text(value) for value in as_list(_STRONG_RULE.get("terms")))
CANDIDA_INTERMEDIATE_EVIDENCE_TERMS = tuple(
    lower_text(value) for value in as_list(_INTERMEDIATE_RULE.get("terms"))
)


def candida_evidence_row_text(row: dict[str, Any]) -> str:
    fields = [str(value) for value in as_list(_candida_rules().get("row_text_fields"))]
    return " ".join(lower_text(row.get(field)) for field in fields)


def candida_evidence_row_strength(row: Any) -> tuple[str, str]:
    if not isinstance(row, dict):
        return CANDIDA_EVIDENCE_WEAK_NONINVASIVE, "non-dict evidence row"
    category = pathogen_names.raw_key(row.get("specimen_category"))
    haystack = candida_evidence_row_text(row)
    if category in INVASIVE_CANDIDA_SPECIMEN_CATEGORIES or any(
        term and term in haystack for term in INVASIVE_CANDIDA_SPECIMEN_TERMS
    ):
        return CANDIDA_EVIDENCE_STRONG_INVASIVE, str(
            _STRONG_RULE.get("basis") or "blood/sterile-site/tissue evidence"
        )
    if any(term and term in haystack for term in CANDIDA_INTERMEDIATE_EVIDENCE_TERMS):
        return CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC, str(
            _INTERMEDIATE_RULE.get("basis") or "systemic fungal marker"
        )
    return CANDIDA_EVIDENCE_WEAK_NONINVASIVE, str(
        _WEAK_RULE.get("basis") or "non-sterile Candida evidence"
    )


def candida_evidence_row_is_invasive(row: Any) -> bool:
    strength, _ = candida_evidence_row_strength(row)
    return strength == CANDIDA_EVIDENCE_STRONG_INVASIVE


def candida_support_details(
    evidence_item: dict[str, Any],
    *,
    strengths: set[str] | None = None,
) -> list[dict[str, Any]]:
    module_evidence = as_dict(evidence_item.get("module_evidence"))
    output: list[dict[str, Any]] = []
    for module in DIRECT_HOSPITAL_MODULES:
        rows = module_evidence.get(module)
        if not isinstance(rows, list):
            continue
        for row in rows:
            strength, basis = candida_evidence_row_strength(row)
            if strengths is not None and strength not in strengths:
                continue
            if isinstance(row, dict):
                output.append(
                    {
                        "module": module,
                        "specimen_type": row.get("specimen_type", "Unknown"),
                        "specimen_category": row.get("specimen_category", "Unknown"),
                        "quantity_tier": row.get("quantity_tier", "Unknown"),
                        "quantitation_status": row.get("quantitation_status", "Unknown"),
                        "evidence_strength": strength,
                        "evidence_strength_basis": basis,
                    }
                )
    return output


def candida_evidence_strength_profile(details: Sequence[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = {
        CANDIDA_EVIDENCE_STRONG_INVASIVE: [],
        CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC: [],
        CANDIDA_EVIDENCE_WEAK_NONINVASIVE: [],
    }
    for detail in details:
        if not isinstance(detail, dict):
            continue
        strength = str(detail.get("evidence_strength") or CANDIDA_EVIDENCE_WEAK_NONINVASIVE)
        grouped.setdefault(strength, []).append(detail)
    if grouped[CANDIDA_EVIDENCE_STRONG_INVASIVE]:
        best = CANDIDA_EVIDENCE_STRONG_INVASIVE
    elif grouped[CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC]:
        best = CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC
    elif grouped[CANDIDA_EVIDENCE_WEAK_NONINVASIVE]:
        best = CANDIDA_EVIDENCE_WEAK_NONINVASIVE
    else:
        best = "none"
    return {
        "best_strength": best,
        "strong_invasive": grouped[CANDIDA_EVIDENCE_STRONG_INVASIVE],
        "intermediate_systemic_or_invasive_clue": grouped[CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC],
        "weak_noninvasive_or_colonization": grouped[CANDIDA_EVIDENCE_WEAK_NONINVASIVE],
    }


def _nested_profile(evidence: dict[str, Any]) -> dict[str, Any]:
    profile = evidence.get("candida_evidence_strength_profile")
    if isinstance(profile, dict):
        return profile
    key_evidence = as_dict(evidence.get("key_evidence"))
    profile = key_evidence.get("candida_evidence_strength_profile")
    return profile if isinstance(profile, dict) else {}


def _module_evidence(evidence: dict[str, Any]) -> dict[str, Any]:
    for container in (evidence, as_dict(evidence.get("key_evidence"))):
        for key in ("hospital_module_evidence", "module_evidence"):
            value = container.get(key)
            if isinstance(value, dict):
                return value
    return {}


def candida_evidence_profile(evidence: dict[str, Any]) -> dict[str, Any]:
    profile = _nested_profile(evidence)
    if profile:
        return profile
    module_evidence = _module_evidence(evidence)
    if module_evidence:
        details = candida_support_details({"module_evidence": module_evidence})
        if details:
            return candida_evidence_strength_profile(details)
    invasive = evidence.get("candida_invasive_hospital_support")
    if not isinstance(invasive, list):
        invasive = as_dict(evidence.get("key_evidence")).get("candida_invasive_hospital_support")
    if isinstance(invasive, list) and invasive:
        details = []
        for item in invasive:
            detail = dict(item) if isinstance(item, dict) else {"detail": str(item)}
            detail["evidence_strength"] = CANDIDA_EVIDENCE_STRONG_INVASIVE
            details.append(detail)
        return candida_evidence_strength_profile(details)
    return {}


def candida_profile_for_candidate(
    evidence: dict[str, Any],
    exact_hospital_evidence: Sequence[dict[str, Any]] = (),
) -> dict[str, Any]:
    profile = candida_evidence_profile(evidence)
    details: list[dict[str, Any]] = []
    for item in exact_hospital_evidence:
        if isinstance(item, dict):
            details.extend(candida_support_details(item))
    if details:
        return candida_evidence_strength_profile(details)
    return profile


def candida_best_evidence_strength(evidence: dict[str, Any]) -> str:
    profile = candida_evidence_profile(evidence)
    return str(profile.get("best_strength") or "")


def candida_has_invasive_context(evidence: dict[str, Any]) -> bool:
    return candida_best_evidence_strength(evidence) == CANDIDA_EVIDENCE_STRONG_INVASIVE


def candida_has_intermediate_context(evidence: dict[str, Any]) -> bool:
    return candida_best_evidence_strength(evidence) == CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC


def _rule_enabled(code: str) -> bool:
    return as_dict(_negative_rules().get(code)).get("enabled") is not False


def _rule_detail(code: str, **values: Any) -> str:
    rule = as_dict(_negative_rules().get(code))
    if rule.get("detail"):
        return str(rule["detail"])
    template = str(rule.get("detail_template") or "")
    return template.format(**values) if template else ""


def _formal_exclusion_rules(snapshot: dict[str, Any]) -> list[str]:
    output: list[str] = []
    for item in as_list(snapshot.get("formal_pick_exclusion_rules")):
        if isinstance(item, dict):
            value = item.get("rule_id") or item.get("reason")
        else:
            value = item
        if value:
            output.append(str(value))
    single = snapshot.get("formal_pick_exclusion_rule")
    if isinstance(single, dict) and (single.get("rule_id") or single.get("reason")):
        output.append(str(single.get("rule_id") or single.get("reason")))
    return output


def build_negative_evidence(
    snapshot: dict[str, Any],
    hospital: dict[str, list[dict[str, Any]]],
    competitors: Sequence[dict[str, Any]],
    *,
    organism_name: str = "",
) -> list[dict[str, str]]:
    result: list[dict[str, str]] = []

    def add(code: str, detail: str) -> None:
        if _rule_enabled(code):
            result.append({"code": code, "detail": detail})

    exact_hospital = as_list(hospital.get("same_organism_or_alias"))
    if not exact_hospital:
        add("NO_SAME_ORGANISM_HOSPITAL_SUPPORT", _rule_detail("NO_SAME_ORGANISM_HOSPITAL_SUPPORT"))

    dominance = str(snapshot.get("dominance_tier") or "")
    dominance_rule = as_dict(_negative_rules().get("NON_DOMINANT_MNGS"))
    if dominance.startswith(tuple(str(value) for value in as_list(dominance_rule.get("tier_prefixes")))):
        add("NON_DOMINANT_MNGS", _rule_detail("NON_DOMINANT_MNGS", dominance_tier=dominance))

    reads_tier = str(snapshot.get("reads_tier") or "")
    reads_rule = as_dict(_negative_rules().get("LOW_MNGS_SIGNAL"))
    if reads_tier.startswith(tuple(str(value) for value in as_list(reads_rule.get("tier_prefixes")))):
        add("LOW_MNGS_SIGNAL", _rule_detail("LOW_MNGS_SIGNAL", reads_tier=reads_tier))

    modules = {lower_text(value) for value in as_list(snapshot.get("support_modules"))}
    non_host = as_list(snapshot.get("non_host_support_modules"))
    support_rule = as_dict(_negative_rules().get("NO_NON_HOST_CONVERGENT_SUPPORT"))
    weak_modules = {lower_text(value) for value in as_list(support_rule.get("weak_support_modules"))}
    if not non_host and modules.issubset(weak_modules):
        add("NO_NON_HOST_CONVERGENT_SUPPORT", _rule_detail("NO_NON_HOST_CONVERGENT_SUPPORT"))

    specimen_class = str(snapshot.get("specimen_class") or "")
    specimen_rule = as_dict(_negative_rules().get("NON_LOWER_RESPIRATORY_SPECIMEN"))
    accepted = {str(value) for value in as_list(specimen_rule.get("accepted_specimen_classes"))}
    if specimen_class and specimen_class not in accepted:
        add(
            "NON_LOWER_RESPIRATORY_SPECIMEN",
            _rule_detail("NON_LOWER_RESPIRATORY_SPECIMEN", specimen_class=specimen_class),
        )

    if competitors:
        prefix = str(as_dict(_negative_rules().get("COMPETING_FORMAL_PATHOGENS_PRESENT")).get("detail_prefix") or "")
        add(
            "COMPETING_FORMAL_PATHOGENS_PRESENT",
            prefix + "; ".join(str(item.get("organism_name") or "") for item in competitors),
        )

    for rule in _formal_exclusion_rules(snapshot):
        add("FORMAL_PICK_EXCLUSION_RULE", rule)

    if pathogen_names.is_candida_or_generic_yeast(organism_name):
        profile = candida_profile_for_candidate(snapshot, exact_hospital)
        best = str(profile.get("best_strength") or "none")
        if best != CANDIDA_EVIDENCE_STRONG_INVASIVE:
            add(
                "CANDIDA_NO_STRONG_INVASIVE_EVIDENCE",
                _rule_detail("CANDIDA_NO_STRONG_INVASIVE_EVIDENCE", best_strength=best),
            )
    return result


def build_candidate_evidence_profile(
    *,
    organism_name: str,
    snapshot: dict[str, Any],
    hospital: dict[str, list[dict[str, Any]]],
    competitors: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    negative = build_negative_evidence(
        snapshot,
        hospital,
        competitors,
        organism_name=organism_name,
    )
    exact = as_list(hospital.get("same_organism_or_alias"))
    related = as_list(hospital.get("same_genus_or_related"))
    candida_profile = (
        candida_profile_for_candidate(snapshot, exact)
        if pathogen_names.is_candida_or_generic_yeast(organism_name)
        else {}
    )
    return {
        "profile_version": PROFILE_VERSION,
        "rule_schema_version": RULE_SCHEMA_VERSION,
        "exact_or_alias_hospital_support": bool(exact),
        "same_genus_only_hospital_context": bool(related) and not bool(exact),
        "negative_evidence_codes": [item["code"] for item in negative],
        "candida_evidence_strength_profile": candida_profile,
    }


def clear_rule_cache() -> None:
    rules_payload.cache_clear()

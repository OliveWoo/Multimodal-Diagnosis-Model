"""Apply bounded history adjustments while preserving the analytical route.

This answer-blind shadow stage consumes test-aware deterministic patient files
and direct structured clinical-history sources.  It records both the route that
was produced without history and the route after history, moves at most one
level, and writes an evidence-consumption ledger so later promotion gates cannot
count the same history twice.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from tools.pathogen_normalization import canonical_key


DEFAULT_POLICY = Path("rules/test_aware_history_route_v4c_combined.json")
SCHEMA_VERSION = "test_aware_history_adjusted_route_shadow.v1"


def read_json(path: Path) -> Any:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if isinstance(payload, dict) and payload.get("base_policy"):
        base_path = path.parent / str(payload["base_policy"])
        base = read_json(base_path)
        if not isinstance(base, dict):
            raise ValueError(f"Base policy must be an object: {base_path}")
        actual_base_sha256 = sha256_file(base_path)
        expected_base_sha256 = str(payload.get("base_policy_sha256") or "").lower()
        if expected_base_sha256 and expected_base_sha256 != actual_base_sha256.lower():
            raise ValueError(
                "Base policy SHA-256 mismatch: "
                f"expected {expected_base_sha256}, got {actual_base_sha256} for {base_path}"
            )
        merged = copy.deepcopy(base)
        for key, value in payload.items():
            if key not in {
                "base_policy", "base_policy_sha256", "constraints_append"
            }:
                merged[key] = value
        merged["constraints"] = list(base.get("constraints") or []) + list(
            payload.get("constraints_append") or []
        )
        merged["policy_composition"] = {
            "base_policy_file": str(base_path.resolve()),
            "base_policy_sha256": actual_base_sha256,
            "overlay_policy_file": str(path.resolve()),
            "overlay_policy_sha256": sha256_file(path),
        }
        chain = list(base.get("policy_composition_chain") or [])
        if not chain and base.get("policy_composition"):
            chain.append(base["policy_composition"])
        chain.append(merged["policy_composition"])
        merged["policy_composition_chain"] = chain
        return merged
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot identify patient from {path}")
    return int(match.group(1))


def _source_record(
    patient: int, source_path: Path, json_path: str, text: Any, source_kind: str
) -> dict[str, Any] | None:
    value = " ".join(str(text or "").split())
    if not value:
        return None
    safe_path = re.sub(r"[^a-zA-Z0-9]+", "_", json_path).strip("_")
    return {
        "evidence_id": f"P{patient}:HIST:{source_kind}:{safe_path}",
        "source_kind": source_kind,
        "text": value,
        "source_ref": {
            "file": str(source_path.resolve()),
            "sha256": sha256_file(source_path),
            "json_path": json_path,
        },
    }


def structured_history_records(patient: int, patient_dir: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    admission = patient_dir / f"NGS_patient_{patient}_admission_diagnosis.json"
    if admission.is_file():
        payload = read_json(admission)
        if isinstance(payload, list):
            for index, row in enumerate(payload):
                if not isinstance(row, dict):
                    continue
                item = _source_record(
                    patient, admission, f"[{index}].diagnosis", row.get("diagnosis"),
                    "admission_diagnosis",
                )
                if item:
                    records.append(item)
                for note_index, note in enumerate(row.get("notes") or []):
                    item = _source_record(
                        patient, admission, f"[{index}].notes[{note_index}]", note,
                        "admission_diagnosis",
                    )
                    if item:
                        records.append(item)

    underlying = patient_dir / f"NGS_patient_{patient}_underlying.json"
    if underlying.is_file():
        payload = read_json(underlying)
        if isinstance(payload, list):
            for index, row in enumerate(payload):
                if not isinstance(row, dict):
                    continue
                for disease_index, disease in enumerate(row.get("underlying_diseases") or []):
                    item = _source_record(
                        patient, underlying,
                        f"[{index}].underlying_diseases[{disease_index}]", disease,
                        "underlying_condition",
                    )
                    if item:
                        records.append(item)

    unique: dict[tuple[str, str], dict[str, Any]] = {}
    for record in records:
        unique[(record["source_ref"]["file"], record["source_ref"]["json_path"])] = record
    return list(unique.values())


def _matches(text: str, terms: Iterable[str]) -> list[str]:
    lowered = text.casefold()
    return sorted({term for term in terms if term.casefold() in lowered})


def _colonization_specimen_scope(text: str, policy: dict[str, Any]) -> str:
    lowered = str(text or "").casefold()
    if any(
        str(term).casefold() in lowered
        for term in policy.get("colonization_sterile_specimen_terms") or []
    ):
        return "sterile_or_systemic"
    if any(
        str(term).casefold() in lowered
        for term in policy.get("colonization_nonsterile_respiratory_terms") or []
    ):
        return "respiratory_nonsterile"
    return "unknown"


def structured_colonization_records(
    patient: int,
    patient_dir: Path,
    history_records: list[dict[str, Any]],
    policy: dict[str, Any],
) -> list[dict[str, Any]]:
    patterns = policy.get("colonization_history_patterns") or {}
    output: list[dict[str, Any]] = []
    for record in history_records:
        categories = [
            category for category, terms in patterns.items()
            if _matches(str(record.get("text") or ""), terms)
        ]
        if not categories:
            continue
        text = str(record.get("text") or "")
        output.append({
            "evidence_id": record["evidence_id"],
            "text": text,
            "categories": sorted(categories),
            "organism_name": None,
            "organism_key": None,
            "specimen_type": None,
            "specimen_category": None,
            "specimen_scope": _colonization_specimen_scope(text, policy),
            "source_kind": record.get("source_kind"),
            "source_refs": [record.get("source_ref")],
            "direct_source_explicit": True,
        })

    for suffix in policy.get("colonization_provisional_culture_agent_suffixes") or []:
        source_path = patient_dir / "agent_outputs" / f"NGS_patient_{patient}_{suffix}"
        if not source_path.is_file():
            continue
        payload = read_json(source_path)
        observations = {
            str(row.get("observation_id")): row
            for row in payload.get("culture_observations") or []
            if isinstance(row, dict)
        }
        for index, label in enumerate(payload.get("organism_labels") or []):
            if not isinstance(label, dict):
                continue
            detail = label.get("key_evidence") or {}
            if str(detail.get("contaminant_flag") or "").casefold() != "yes":
                continue
            organism = str(label.get("organism_name") or "").strip()
            if not organism:
                continue
            observation_ids = [
                str(value) for value in label.get("evidence_observation_ids") or []
            ]
            linked = [observations[value] for value in observation_ids if value in observations]
            specimen_type = str(
                detail.get("specimen_type")
                or next((row.get("specimen_type") for row in linked if row.get("specimen_type")), "")
            )
            specimen_category = str(
                detail.get("specimen_category")
                or next((row.get("specimen_category") for row in linked if row.get("specimen_category")), "")
            )
            normalized = canonical_key(
                f"{organism}|{specimen_type}|{specimen_category}|{'|'.join(observation_ids)}"
            )
            output.append({
                "evidence_id": (
                    f"P{patient}:COLONIZATION:CULTURE_AGENT:"
                    f"{hashlib.sha256(normalized.encode('utf-8')).hexdigest()[:16]}"
                ),
                "text": (
                    f"{organism}; provisional contaminant flag; "
                    f"specimen {specimen_type or 'unknown'}"
                ),
                "categories": ["provisional_model_contaminant_label"],
                "organism_name": organism,
                "organism_key": canonical_key(organism),
                "specimen_type": specimen_type or None,
                "specimen_category": specimen_category or None,
                "specimen_scope": _colonization_specimen_scope(
                    f"{specimen_type} {specimen_category}", policy
                ),
                "source_kind": "provisional_culture_agent_label",
                "source_refs": [{
                    "file": str(source_path.resolve()),
                    "sha256": sha256_file(source_path),
                    "json_path": f"organism_labels[{index}]",
                    "observation_ids": observation_ids,
                }],
                "direct_source_explicit": False,
            })
    return output


def colonization_history_profile(
    records: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    for record in records:
        for category in record.get("categories") or []:
            counts[str(category)] += 1
    direct_categories = {
        "explicit_candidate_colonization",
        "explicit_current_not_infection",
    }
    has_direct = any(counts[category] > 0 for category in direct_categories)
    has_provisional = counts["provisional_model_contaminant_label"] > 0
    has_structural = counts["structural_airway_context_only"] > 0
    if has_direct:
        support = "direct_colonization_or_not_infection_text"
    elif has_provisional:
        support = "provisional_contaminant_label_not_deterministic"
    elif has_structural:
        support = "structural_airway_context_only"
    else:
        support = "none"
    return {
        "support_level": support,
        "has_direct_colonization_text": has_direct,
        "has_provisional_contaminant_label": has_provisional,
        "has_structural_airway_context_only": has_structural,
        "category_counts": dict(sorted(counts.items())),
        "matched_evidence": records,
    }


def candidate_colonization_history(
    candidate: dict[str, Any], history: dict[str, Any]
) -> dict[str, Any]:
    key = canonical_key(candidate.get("organism_name"))
    matched: list[dict[str, Any]] = []
    for record in history.get("matched_evidence") or []:
        record_key = str(record.get("organism_key") or "")
        text_key = canonical_key(record.get("text"))
        if key and (record_key == key or (not record_key and key in text_key)):
            matched.append(record)
    counts: Counter[str] = Counter()
    for record in matched:
        for category in record.get("categories") or []:
            counts[str(category)] += 1
    eligible_direct = [
        record for record in matched
        if record.get("direct_source_explicit")
        and record.get("specimen_scope") == "respiratory_nonsterile"
        and set(record.get("categories") or []) & {
            "explicit_candidate_colonization",
            "explicit_current_not_infection",
        }
    ]
    output = dict(history)
    output.update({
        "support_level": (
            "candidate_matched_direct_nonsterile_colonization"
            if eligible_direct
            else "candidate_matched_but_not_eligible"
            if matched else "none"
        ),
        "has_candidate_match": bool(matched),
        "has_direct_nonsterile_match": bool(eligible_direct),
        "has_explicit_colonization": any(
            "explicit_candidate_colonization" in set(row.get("categories") or [])
            for row in eligible_direct
        ),
        "has_explicit_current_not_infection": any(
            "explicit_current_not_infection" in set(row.get("categories") or [])
            for row in eligible_direct
        ),
        "has_sterile_or_systemic_label": any(
            row.get("specimen_scope") == "sterile_or_systemic" for row in matched
        ),
        "has_only_provisional_label": bool(matched) and not any(
            row.get("direct_source_explicit") for row in matched
        ),
        "consumable_categories": [
            "explicit_candidate_colonization",
            "explicit_current_not_infection",
        ],
        "category_counts": dict(sorted(counts.items())),
        "matched_evidence": matched,
        "eligible_direct_evidence": eligible_direct,
    })
    return output


def _history_claims(
    records: list[dict[str, Any]], patterns: dict[str, list[str]]
) -> tuple[list[dict[str, Any]], Counter[str]]:
    claims: dict[str, dict[str, Any]] = {}
    for record in records:
        categories: list[str] = []
        terms: list[str] = []
        for category, configured_terms in patterns.items():
            found = _matches(record["text"], configured_terms)
            if found:
                categories.append(category)
                terms.extend(found)
        if not categories:
            continue
        normalized = re.sub(r"[^a-z0-9]+", "", record["text"].casefold())
        claim = claims.setdefault(normalized, {
            "evidence_id": (
                f"{str(record['evidence_id']).split(':', 1)[0]}:HISTCLAIM:"
                f"{hashlib.sha256(normalized.encode('utf-8')).hexdigest()[:16]}"
            ),
            "text": record["text"],
            "categories": [],
            "matched_terms": [],
            "source_refs": [],
        })
        claim["categories"] = sorted(set(claim["categories"]) | set(categories))
        claim["matched_terms"] = sorted(set(claim["matched_terms"]) | set(terms))
        claim["source_refs"].append(record["source_ref"])

    matched = list(claims.values())
    category_counts: Counter[str] = Counter()
    for claim in matched:
        for category in claim["categories"]:
            category_counts[category] += 1
    return matched, category_counts


def aspiration_history_profile(
    records: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    patterns = policy["history_patterns"]
    matched, category_counts = _history_claims(records, patterns)

    explicit = category_counts["explicit_aspiration_event"] > 0
    possible = category_counts["possible_aspiration_event"] > 0
    anatomic = category_counts["durable_anatomic_risk"] > 0
    neurologic = category_counts["durable_neurologic_risk"] > 0
    if anatomic:
        support = "durable_family_specific_risk"
    elif explicit:
        support = "prior_or_current_aspiration_event_timing_unresolved"
    elif possible:
        support = "possible_aspiration_event_timing_unresolved"
    elif neurologic:
        support = "nonspecific_neurologic_risk_only"
    else:
        support = "none"
    return {
        "support_level": support,
        "has_explicit_aspiration_event": explicit,
        "has_possible_aspiration_event": possible,
        "has_durable_anatomic_risk": anatomic,
        "has_durable_neurologic_risk": neurologic,
        "qualifies_for_bounded_adjustment": anatomic or explicit or possible,
        "qualifies_for_context_to_priority": anatomic,
        "category_counts": dict(sorted(category_counts.items())),
        "matched_evidence": matched,
    }


def opportunistic_history_profile(
    records: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    patterns = policy.get("opportunistic_history_patterns") or {}
    matched, category_counts = _history_claims(records, patterns)
    strong_categories = {
        "hematologic_malignancy_or_cellular_therapy",
        "systemic_immunosuppressive_therapy",
        "cytopenia_high_risk",
    }
    context_categories = {
        "solid_tumor_or_targeted_therapy_context",
        "nonspecific_host_susceptibility_context",
    }
    has_strong = any(category_counts[category] > 0 for category in strong_categories)
    has_context = any(category_counts[category] > 0 for category in context_categories)
    if category_counts["hematologic_malignancy_or_cellular_therapy"]:
        support = "high_confidence_hematologic_or_cellular_therapy_risk"
    elif category_counts["systemic_immunosuppressive_therapy"]:
        support = "high_confidence_systemic_immunosuppressive_therapy_timing_unresolved"
    elif category_counts["cytopenia_high_risk"]:
        support = "high_confidence_cytopenia_risk_timing_unresolved"
    elif has_context:
        support = "context_only_nonspecific_host_susceptibility"
    else:
        support = "none"
    return {
        "support_level": support,
        "has_high_confidence_immunosuppression": has_strong,
        "has_context_only_host_susceptibility": has_context,
        "qualifies_for_bounded_adjustment": has_strong,
        "qualifies_for_context_to_priority": has_strong,
        "consumable_categories": sorted(strong_categories),
        "category_counts": dict(sorted(category_counts.items())),
        "matched_evidence": matched,
    }


def hap_mdr_history_profile(
    records: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    patterns = policy.get("hap_mdr_history_patterns") or {}
    matched, category_counts = _history_claims(records, patterns)
    setting_categories = {
        "hap_setting",
        "healthcare_associated_setting",
        "airway_device_exposure",
    }
    exact_mdr_categories = {
        "prior_mdr_klebsiella",
        "prior_mdr_acinetobacter",
    }
    has_setting = any(category_counts[category] > 0 for category in setting_categories)
    has_exact_mdr = any(category_counts[category] > 0 for category in exact_mdr_categories)
    has_generic_mdr = category_counts["generic_mdr_context"] > 0
    if has_exact_mdr:
        support = "organism_named_prior_mdr_history"
    elif category_counts["hap_setting"]:
        support = "hap_setting_context"
    elif category_counts["airway_device_exposure"]:
        support = "airway_device_setting_context"
    elif category_counts["healthcare_associated_setting"] or has_generic_mdr:
        support = "healthcare_exposure_context"
    else:
        support = "none"
    return {
        "support_level": support,
        "has_hap_or_airway_device_setting": has_setting,
        "has_exact_organism_mdr_history": has_exact_mdr,
        "has_generic_mdr_context": has_generic_mdr,
        "qualifies_for_bounded_adjustment": has_setting or has_exact_mdr,
        "setting_categories": sorted(setting_categories),
        "exact_mdr_categories": sorted(exact_mdr_categories),
        "consumable_categories": sorted(setting_categories | exact_mdr_categories),
        "category_counts": dict(sorted(category_counts.items())),
        "matched_evidence": matched,
    }


def reactivation_history_profile(
    records: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    patterns = policy.get("reactivation_history_patterns") or {}
    matched, category_counts = _history_claims(records, patterns)
    has_resolved = category_counts["resolved_or_asymptomatic_shedding"] > 0
    has_explicit = category_counts["explicit_current_reactivation"] > 0
    has_quantitative = category_counts["quantitative_or_tissue_evidence"] > 0
    has_pulmonary = category_counts["compatible_pulmonary_syndrome"] > 0
    has_prior = any(
        category_counts[category] > 0
        for category in {
            "prior_or_suspected_infection",
            "current_infection_timing_unresolved",
        }
    )
    if has_resolved:
        support = "resolved_or_asymptomatic_shedding"
    elif has_explicit and has_quantitative and has_pulmonary:
        support = "current_reactivation_with_independent_pulmonary_evidence"
    elif has_explicit:
        support = "explicit_reactivation_missing_independent_convergence"
    elif has_prior:
        support = "prior_or_timing_unresolved_infection_context"
    else:
        support = "none"
    return {
        "support_level": support,
        "has_resolved_or_asymptomatic_shedding": has_resolved,
        "has_explicit_current_reactivation": has_explicit,
        "has_quantitative_or_tissue_evidence": has_quantitative,
        "has_compatible_pulmonary_syndrome": has_pulmonary,
        "has_prior_or_timing_unresolved_infection": has_prior,
        "category_counts": dict(sorted(category_counts.items())),
        "matched_evidence": matched,
    }


def candidate_reactivation_history(
    candidate: dict[str, Any], history: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any]:
    key = canonical_key(candidate.get("organism_name"))
    terms = [
        str(term).casefold()
        for term in (policy.get("reactivation_candidate_terms") or {}).get(key, [])
    ]
    matched = [
        item for item in history.get("matched_evidence") or []
        if terms and any(term in str(item.get("text") or "").casefold() for term in terms)
    ]
    category_counts: Counter[str] = Counter()
    for item in matched:
        for category in item.get("categories") or []:
            category_counts[category] += 1
    has_resolved = category_counts["resolved_or_asymptomatic_shedding"] > 0
    has_explicit = category_counts["explicit_current_reactivation"] > 0
    has_quantitative = category_counts["quantitative_or_tissue_evidence"] > 0
    has_pulmonary = category_counts["compatible_pulmonary_syndrome"] > 0
    has_prior = any(
        category_counts[category] > 0
        for category in {
            "prior_or_suspected_infection",
            "current_infection_timing_unresolved",
        }
    )
    if has_resolved:
        support = "resolved_or_asymptomatic_shedding"
    elif has_explicit and has_quantitative and has_pulmonary:
        support = "current_reactivation_with_independent_pulmonary_evidence"
    elif has_explicit:
        support = "explicit_reactivation_missing_independent_convergence"
    elif has_prior:
        support = "prior_or_timing_unresolved_infection_context"
    else:
        support = "none"
    output = dict(history)
    output.update({
        "support_level": support,
        "has_resolved_or_asymptomatic_shedding": has_resolved,
        "has_explicit_current_reactivation": has_explicit,
        "has_quantitative_or_tissue_evidence": has_quantitative,
        "has_compatible_pulmonary_syndrome": has_pulmonary,
        "has_prior_or_timing_unresolved_infection": has_prior,
        "qualifies_for_context_to_audit": has_resolved,
        "qualifies_for_context_to_priority": (
            has_explicit and has_quantitative and has_pulmonary
        ),
        "consumable_categories": sorted({
            "resolved_or_asymptomatic_shedding",
            "explicit_current_reactivation",
            "quantitative_or_tissue_evidence",
            "compatible_pulmonary_syndrome",
        }),
        "category_counts": dict(sorted(category_counts.items())),
        "matched_evidence": matched,
    })
    return output


def candidate_hap_mdr_history(
    candidate: dict[str, Any], history: dict[str, Any]
) -> tuple[dict[str, Any], bool]:
    name = str(candidate.get("organism_name") or "").casefold()
    setting_categories = set(history.get("setting_categories") or [])
    exact_categories: set[str] = set()
    if "klebsiella" in name:
        exact_categories.add("prior_mdr_klebsiella")
    if "acinetobacter" in name:
        exact_categories.add("prior_mdr_acinetobacter")
    allowed = setting_categories | exact_categories
    matched = [
        item for item in history.get("matched_evidence") or []
        if set(item.get("categories") or []) & allowed
    ]
    exact_match = any(
        set(item.get("categories") or []) & exact_categories for item in matched
    )
    candidate_history = dict(history)
    candidate_history["matched_evidence"] = matched
    candidate_history["has_candidate_specific_mdr_match"] = exact_match
    candidate_history["qualifies_for_bounded_adjustment"] = bool(matched)
    candidate_history["consumable_categories"] = sorted(allowed)
    return candidate_history, exact_match


def phenotype_disagreement(
    patient: int, phenotype_root: Path | None, phenotype: str = "ASPIRATION_RISK"
) -> dict[str, Any] | None:
    if phenotype_root is None:
        return None
    path = phenotype_root / f"patient_{patient}_candidate_phenotype_shadow_v1.json"
    if not path.is_file():
        return None
    payload = read_json(path)
    rows = [
        row for row in payload.get("all_phenotype_overview") or []
        if str(row.get("phenotype") or "") == phenotype
    ]
    if not rows:
        return None
    row = rows[0]
    return {
        "phenotype": phenotype,
        "source_file": str(path.resolve()),
        "source_sha256": sha256_file(path),
        "reported_status": row.get("status"),
        "confidence": row.get("confidence"),
        "temporal_relation": row.get("temporal_relation"),
        "evidence_present": bool(row.get("evidence")),
        "interpretation": (
            "not_established_not_negative_evidence"
            if str(row.get("status") or "").upper() == "NO"
            else "phenotype_context_only"
        ),
    }


def _candidate_family(candidate: dict[str, Any]) -> str:
    return str(
        candidate.get("taxonomy_family")
        or (candidate.get("taxonomy_profile") or {}).get("primary_rule_family")
        or "unmapped_or_uncertain"
    )


def _analytic_route(candidate: dict[str, Any], policy: dict[str, Any]) -> str:
    decision = str(candidate.get("decision") or candidate.get("clinical_decision") or "hold")
    return str(policy["decision_to_analytic_route"].get(decision, "Hold"))


def _analytical_axes(candidate: dict[str, Any]) -> dict[str, Any]:
    analytical = candidate.get("analytical_profile") or {}
    rank = analytical.get("best_rank_in_retained_universe")
    tests = analytical.get("selected_positive_test_count")
    try:
        rank = int(rank)
    except (TypeError, ValueError):
        rank = None
    try:
        tests = int(tests)
    except (TypeError, ValueError):
        tests = 0
    return {
        "best_rank_in_retained_universe": rank,
        "selected_positive_test_count": tests,
        "reproducibility_axis": bool(analytical.get("reproducibility_axis")),
        "cross_molecule_selected": bool(analytical.get("cross_molecule_selected")),
        "positive_mngs_signal": tests > 0,
    }


def _candidate_specimen_context(candidate: dict[str, Any]) -> str:
    return str(
        candidate.get("specimen_context")
        or (candidate.get("source_entry") or {}).get("specimen_context")
        or "unknown"
    )


def _has_exact_hospital_support(candidate: dict[str, Any]) -> bool:
    if candidate.get("evidence_source") == "hospital_only":
        detail = candidate.get("hospital_evidence_detail") or {}
        return bool(detail.get("best_hospital_level") or detail.get("source_observations"))
    profile = candidate.get("hospital_profile") or {}
    return bool(
        profile.get("exact_or_alias_evidence")
        or profile.get("best_direct_level")
        or (profile.get("hospital_evidence_detail") or {}).get("source_observations")
    )


def strict_aspiration_anaerobe(candidate: dict[str, Any], policy: dict[str, Any]) -> bool:
    name = str(candidate.get("organism_name") or "").casefold()
    return any(
        term.casefold() in name
        for term in policy.get("strict_aspiration_anaerobe_genus_terms") or []
    )


def adjust_candidate(
    candidate: dict[str, Any], history: dict[str, Any], policy: dict[str, Any],
    disagreement: dict[str, Any] | None = None,
    opportunistic_history: dict[str, Any] | None = None,
    opportunistic_disagreement: dict[str, Any] | None = None,
    hap_mdr_history: dict[str, Any] | None = None,
    hap_mdr_disagreement: dict[str, Any] | None = None,
    reactivation_history: dict[str, Any] | None = None,
    reactivation_disagreement: dict[str, Any] | None = None,
    colonization_history: dict[str, Any] | None = None,
    colonization_disagreement: dict[str, Any] | None = None,
) -> dict[str, Any]:
    output = copy.deepcopy(candidate)
    analytic = _analytic_route(candidate, policy)
    adjusted = analytic
    family = _candidate_family(candidate)
    axes = _analytical_axes(candidate)
    rule_ids = ["HIST-S1-PRESERVE-ANALYTIC-ROUTE", "HIST-S2-MAX-ONE-LEVEL"]
    blockers: list[str] = []
    consumed: list[str] = []
    selected_history = history
    selected_disagreement = disagreement
    history_family = "none"
    opportunistic_families = set(policy.get("opportunistic_families") or [])
    hap_mdr_family = policy.get("hap_mdr_family")
    reactivation_family = policy.get("reactivation_family")
    colonization_review: dict[str, Any] | None = None

    if family == policy["aspiration_family"]:
        history_family = "aspiration"
        if not history["qualifies_for_bounded_adjustment"]:
            blockers.append("no qualifying direct-source aspiration history or durable anatomic risk")
        elif not axes["positive_mngs_signal"]:
            blockers.append("no evaluable positive mNGS test")
        elif analytic == "Hold":
            adjusted = "Audit"
            rule_ids.append("HIST-ASP1-HOLD-TO-AUDIT")
        elif analytic == "Audit":
            minimum = int(policy["analytical_thresholds"]["audit_to_context_min_positive_tests"])
            maximum_rank = int(policy["analytical_thresholds"]["audit_to_context_max_rank"])
            strong_enough = (
                (axes["best_rank_in_retained_universe"] is not None
                 and axes["best_rank_in_retained_universe"] <= maximum_rank)
                or axes["reproducibility_axis"]
                or axes["cross_molecule_selected"]
            )
            if not strict_aspiration_anaerobe(candidate, policy):
                blockers.append("not a strict aspiration-anaerobe representative")
            elif axes["selected_positive_test_count"] < minimum:
                blockers.append("insufficient evaluable positive tests for Audit-to-Context")
            elif not strong_enough:
                blockers.append("weak single-test signal is below the Audit-to-Context rank threshold")
            else:
                adjusted = "Context"
                rule_ids.append("HIST-ASP2-AUDIT-TO-CONTEXT")
        elif analytic == "Context":
            maximum_rank = int(policy["analytical_thresholds"]["context_to_priority_max_rank"])
            strong = (
                (axes["best_rank_in_retained_universe"] is not None
                 and axes["best_rank_in_retained_universe"] <= maximum_rank)
                or axes["reproducibility_axis"]
                or axes["cross_molecule_selected"]
            )
            if not history["qualifies_for_context_to_priority"]:
                blockers.append(
                    "history is a prior/uncertain aspiration event without durable anatomic risk"
                )
            elif strong:
                adjusted = "Priority"
                rule_ids.append("HIST-ASP3-CONTEXT-TO-PRIORITY")
            else:
                blockers.append("no top-three or reproducible analytical support")
        else:
            rule_ids.append("HIST-ASP4-PRIORITY-PRESERVED")
    elif family == reactivation_family and reactivation_history is not None:
        history_family = "reactivation_interpretation"
        selected_history = candidate_reactivation_history(
            candidate, reactivation_history, policy
        )
        selected_disagreement = reactivation_disagreement
        thresholds = policy["reactivation_analytical_thresholds"]
        rank = axes["best_rank_in_retained_universe"]
        reproducible = axes["reproducibility_axis"] or axes["cross_molecule_selected"]
        weak = (
            (rank is None or rank > int(thresholds["weak_signal_min_rank_exclusive"]))
            and axes["selected_positive_test_count"]
            <= int(thresholds["weak_signal_max_positive_tests"])
            and not reproducible
        )
        strong = (
            rank is not None
            and rank <= int(thresholds["context_to_priority_max_rank"])
            and axes["selected_positive_test_count"]
            >= int(thresholds["context_to_priority_min_positive_tests"])
            and reproducible
        )
        if not selected_history["matched_evidence"]:
            blockers.append("no candidate-matched direct-source reactivation history")
        elif analytic != "Context":
            blockers.append(
                "v4A only adjusts Context candidates; other routes remain unchanged"
            )
        elif selected_history["qualifies_for_context_to_audit"]:
            if weak:
                adjusted = "Audit"
                rule_ids.append("HIST-REACT1-RESOLVED-WEAK-CONTEXT-TO-AUDIT")
            else:
                blockers.append(
                    "resolved/asymptomatic history requires a weak single-test signal for Context-to-Audit"
                )
        elif selected_history["qualifies_for_context_to_priority"]:
            if strong:
                adjusted = "Priority"
                rule_ids.append("HIST-REACT2-CONVERGENT-CONTEXT-TO-PRIORITY")
            else:
                blockers.append(
                    "current reactivation requires rank at most five plus at least two reproducible positive tests"
                )
        elif selected_history["has_explicit_current_reactivation"]:
            blockers.append(
                "explicit reactivation lacks quantitative/tissue evidence or a compatible pulmonary syndrome"
            )
        else:
            blockers.append(
                "prior or timing-unresolved infection history is context only"
            )
    elif family in opportunistic_families and opportunistic_history is not None:
        history_family = "opportunistic_host_risk"
        selected_history = opportunistic_history
        selected_disagreement = opportunistic_disagreement
        thresholds = policy["opportunistic_analytical_thresholds"]
        rank = axes["best_rank_in_retained_universe"]
        reproducible = axes["reproducibility_axis"] or axes["cross_molecule_selected"]
        if not opportunistic_history["qualifies_for_bounded_adjustment"]:
            blockers.append(
                "only nonspecific host susceptibility is present; no qualifying direct-source immunosuppression"
            )
        elif not axes["positive_mngs_signal"]:
            blockers.append("no evaluable positive mNGS test")
        elif analytic == "Hold":
            adjusted = "Audit"
            rule_ids.append("HIST-OPP1-HOLD-TO-AUDIT")
        elif analytic == "Audit":
            strong = (
                rank is not None
                and rank <= int(thresholds["audit_to_context_max_rank"])
                and axes["selected_positive_test_count"]
                >= int(thresholds["audit_to_context_min_positive_tests"])
                and reproducible
            )
            if strong:
                adjusted = "Context"
                rule_ids.append("HIST-OPP2-AUDIT-TO-CONTEXT")
            else:
                blockers.append(
                    "opportunistic Audit-to-Context requires rank at most ten plus at least two reproducible positive tests"
                )
        elif analytic == "Context":
            strong = (
                rank is not None
                and rank <= int(thresholds["context_to_priority_max_rank"])
                and axes["selected_positive_test_count"]
                >= int(thresholds["context_to_priority_min_positive_tests"])
                and reproducible
            )
            if strong:
                adjusted = "Priority"
                rule_ids.append("HIST-OPP3-CONTEXT-TO-PRIORITY")
            else:
                blockers.append(
                    "opportunistic Context-to-Priority requires rank at most ten plus at least two reproducible positive tests"
                )
        else:
            rule_ids.append("HIST-OPP4-PRIORITY-PRESERVED")
    elif family == hap_mdr_family and hap_mdr_history is not None:
        history_family = "hap_mdr_exposure"
        selected_history, exact_mdr_match = candidate_hap_mdr_history(
            candidate, hap_mdr_history
        )
        selected_disagreement = hap_mdr_disagreement
        thresholds = policy["hap_mdr_analytical_thresholds"]
        rank = axes["best_rank_in_retained_universe"]
        reproducible = axes["reproducibility_axis"] or axes["cross_molecule_selected"]
        if not selected_history["qualifies_for_bounded_adjustment"]:
            blockers.append(
                "no applicable HAP, airway-device, or candidate-matched prior MDR history"
            )
        elif not axes["positive_mngs_signal"]:
            blockers.append("no evaluable positive mNGS test")
        elif analytic == "Hold":
            adjusted = "Audit"
            rule_ids.append("HIST-HAP1-HOLD-TO-AUDIT")
        elif analytic == "Audit":
            strong = (
                rank is not None
                and rank <= int(thresholds["audit_to_context_max_rank"])
                and axes["selected_positive_test_count"]
                >= int(thresholds["audit_to_context_min_positive_tests"])
                and reproducible
            )
            if strong:
                adjusted = "Context"
                rule_ids.append("HIST-HAP2-AUDIT-TO-CONTEXT")
            else:
                blockers.append(
                    "HAP/MDR Audit-to-Context requires rank at most ten plus at least two reproducible positive tests"
                )
        elif analytic == "Context":
            maximum_rank = int(
                thresholds[
                    "exact_mdr_context_to_priority_max_rank"
                    if exact_mdr_match
                    else "setting_context_to_priority_max_rank"
                ]
            )
            strong = (
                rank is not None
                and rank <= maximum_rank
                and axes["selected_positive_test_count"]
                >= int(thresholds["context_to_priority_min_positive_tests"])
                and reproducible
            )
            if strong:
                adjusted = "Priority"
                rule_ids.append(
                    "HIST-MDR3-EXACT-MATCH-CONTEXT-TO-PRIORITY"
                    if exact_mdr_match
                    else "HIST-HAP3-SETTING-CONTEXT-TO-PRIORITY"
                )
            else:
                blockers.append(
                    f"HAP/MDR Context-to-Priority requires rank at most {maximum_rank} "
                    "plus at least two reproducible positive tests"
                )
        else:
            rule_ids.append("HIST-HAP4-PRIORITY-PRESERVED")
    else:
        blockers.append("no implemented history route applies to this taxonomy family")

    if colonization_history is not None and policy.get("colonization_analytical_thresholds"):
        colonization_review = candidate_colonization_history(
            candidate, colonization_history
        )
        if colonization_review["has_candidate_match"]:
            if adjusted != analytic:
                blockers.append(
                    "candidate-matched colonization evidence conflicts with another history adjustment; existing single-level adjustment retained"
                )
            else:
                history_family = "colonization_interpretation"
                selected_history = colonization_review
                selected_disagreement = colonization_disagreement
                blockers = []
                thresholds = policy["colonization_analytical_thresholds"]
                rank = axes["best_rank_in_retained_universe"]
                nonreproducible = not (
                    axes["reproducibility_axis"] or axes["cross_molecule_selected"]
                )
                weak = bool(
                    nonreproducible
                    and (
                        rank is None
                        or rank > int(thresholds["weak_signal_min_rank_exclusive"])
                        or axes["selected_positive_test_count"]
                        <= int(thresholds["weak_signal_max_positive_tests"])
                    )
                )
                very_weak = bool(
                    nonreproducible
                    and (rank is None or rank > int(
                        thresholds["very_weak_signal_min_rank_exclusive"]
                    ))
                    and axes["selected_positive_test_count"]
                    <= int(thresholds["very_weak_signal_max_positive_tests"])
                )
                if not axes["positive_mngs_signal"]:
                    blockers.append(
                        "colonization v4B only adjusts positive mNGS candidates"
                    )
                elif _candidate_specimen_context(candidate) != "lower_respiratory":
                    blockers.append(
                        "candidate is not linked to a lower-respiratory mNGS specimen"
                    )
                elif colonization_review["has_sterile_or_systemic_label"]:
                    blockers.append(
                        "candidate-matched colonization/contaminant label is from a sterile or systemic specimen"
                    )
                elif colonization_review["has_only_provisional_label"]:
                    blockers.append(
                        "model-derived contaminant label is provisional and cannot deterministically demote"
                    )
                elif not colonization_review["has_direct_nonsterile_match"]:
                    blockers.append(
                        "no direct exact-organism colonization/not-infection evidence from the same nonsterile respiratory compartment"
                    )
                elif _has_exact_hospital_support(candidate):
                    blockers.append(
                        "exact hospital culture/PCR support blocks colonization demotion"
                    )
                elif (
                    analytic == "Priority"
                    and colonization_review["has_explicit_current_not_infection"]
                ):
                    adjusted = "Context"
                    rule_ids.append(
                        "HIST-COL1-CURRENT-NOT-INFECTION-PRIORITY-TO-CONTEXT"
                    )
                elif (
                    analytic == "Context"
                    and colonization_review["has_explicit_colonization"]
                    and weak
                ):
                    adjusted = "Audit"
                    rule_ids.append(
                        "HIST-COL2-EXPLICIT-WEAK-CONTEXT-TO-AUDIT"
                    )
                elif (
                    analytic == "Audit"
                    and colonization_review["has_explicit_colonization"]
                    and very_weak
                ):
                    adjusted = "Hold"
                    rule_ids.append(
                        "HIST-COL3-EXPLICIT-VERY-WEAK-AUDIT-TO-HOLD"
                    )
                elif analytic not in {"Priority", "Context", "Audit"}:
                    blockers.append(
                        "colonization v4B has no downward route from Hold"
                    )
                else:
                    blockers.append(
                        "candidate does not meet the route-specific colonization demotion threshold"
                    )

    route_order = list(policy["route_order"])
    if abs(route_order.index(adjusted) - route_order.index(analytic)) > 1:
        raise ValueError(f"History adjustment exceeds one level: {analytic} -> {adjusted}")
    if adjusted != analytic:
        consumable_categories = (
            {
                "explicit_aspiration_event", "possible_aspiration_event",
                "durable_anatomic_risk",
            }
            if history_family == "aspiration"
            else set(selected_history.get("consumable_categories") or [])
        )
        consumed = [
            item["evidence_id"] for item in selected_history["matched_evidence"]
            if set(item["categories"]) & consumable_categories
        ]

    candidate_history_ids = (
        [item["evidence_id"] for item in selected_history["matched_evidence"]]
        if history_family != "none" else []
    )
    original_decision = str(candidate.get("decision") or candidate.get("clinical_decision") or "hold")
    if original_decision == "picked_shadow":
        history_routed_decision = "picked_shadow"
        downstream_tier = "Picked"
    elif adjusted == "Priority":
        history_routed_decision = "review_high_priority"
        downstream_tier = "High"
    elif adjusted == "Context":
        history_routed_decision = "review_context_needed"
        downstream_tier = "Context"
    elif adjusted == "Audit":
        history_routed_decision = "review_low_specificity"
        downstream_tier = "Audit"
    else:
        history_routed_decision = "hold"
        downstream_tier = "Hold"
    downstream_blockers: list[str] = []
    if history_family == "aspiration" and downstream_tier == "High":
        downstream_blockers.extend([
            "history evidence was consumed at the route stage and cannot be counted again",
            "oral/aspiration family still requires independent same-organism or invasive evidence for Picked",
        ])
    elif history_family == "opportunistic_host_risk" and downstream_tier == "High":
        downstream_blockers.extend([
            "history evidence was consumed at the route stage and cannot be counted again",
            "opportunistic family still requires independent organism-specific evidence for Picked",
        ])
    elif history_family == "hap_mdr_exposure" and downstream_tier == "High":
        downstream_blockers.extend([
            "HAP/MDR history was consumed at the route stage and cannot be counted again",
            "healthcare setting does not identify a pathogen and cannot independently create Picked",
        ])
    elif history_family == "reactivation_interpretation" and downstream_tier == "High":
        downstream_blockers.extend([
            "reactivation evidence consumed at the route stage cannot be counted again",
            "Priority requires a later independent disease-specific Picked-gate audit",
        ])

    future_picked_review = bool(
        history_family == "reactivation_interpretation" and adjusted == "Priority"
    )

    output.update({
        "analytic_route": analytic,
        "history_adjusted_route": adjusted,
        "history_route_changed": adjusted != analytic,
        "history_rule_ids": rule_ids,
        "history_adjustment_blockers": blockers,
        "history_policy_family": history_family,
        "history_evidence_ids": candidate_history_ids,
        "consumed_evidence_ids": consumed,
        "history_evidence_consumption": [
            {
                "evidence_id": evidence_id,
                "consumed_at_stage": "history_adjusted_route",
                "downstream_reuse": "provenance_only_not_independent_support",
            }
            for evidence_id in consumed
        ],
        "history_source_disagreement": selected_disagreement,
        "colonization_interpretation": colonization_review,
        "history_routed_decision": history_routed_decision,
        "future_picked_gate_review_required": future_picked_review,
        "future_picked_gate_review_status": (
            "pending_independent_evidence_audit"
            if future_picked_review else "not_applicable"
        ),
        "future_picked_gate_review_requirements": (
            policy.get("reactivation_future_picked_gate_requirements") or []
            if future_picked_review else []
        ),
        "downstream_route": {
            "scorer_review_tier": downstream_tier,
            "projected_final_tier": downstream_tier,
            "picked_allowed_from_history_alone": False,
            "picked_promotion_blockers": downstream_blockers,
            "note": (
                "Priority enters high-priority review; Picked still requires an independent "
                "family-specific promotion route."
            ),
        },
    })
    return output


def apply_aspiration_cluster_guardrail(
    candidates: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    promoted = [
        item for item in candidates
        if item.get("analytic_route") == "Audit"
        and item.get("history_adjusted_route") == "Context"
        and _candidate_family(item) == policy["aspiration_family"]
    ]
    maximum = int(
        policy["analytical_thresholds"]["audit_to_context_max_representatives_per_patient"]
    )
    promoted.sort(key=lambda item: (
        not _analytical_axes(item)["reproducibility_axis"],
        _analytical_axes(item)["best_rank_in_retained_universe"] or 999,
        str(item.get("organism_name") or "").casefold(),
    ))
    retained = promoted[:maximum]
    suppressed = promoted[maximum:]
    for item in suppressed:
        item["history_adjusted_route"] = "Audit"
        item["history_route_changed"] = False
        item["history_routed_decision"] = "review_low_specificity"
        item["history_rule_ids"] = [
            rule for rule in item["history_rule_ids"]
            if rule != "HIST-ASP2-AUDIT-TO-CONTEXT"
        ] + ["HIST-ASP6-CLUSTER-REPRESENTATIVE-CAP"]
        item["history_adjustment_blockers"].append(
            "aspiration cluster representative cap reached; retained in Audit"
        )
        item["consumed_evidence_ids"] = []
        item["history_evidence_consumption"] = []
        item["downstream_route"] = {
            "scorer_review_tier": "Audit",
            "projected_final_tier": "Audit",
            "picked_allowed_from_history_alone": False,
            "picked_promotion_blockers": [],
            "note": "Audit candidates are retained for provenance and are not forwarded to the clinical scorer.",
        }
    return {
        "candidate_count_before_cap": len(promoted),
        "representative_cap": maximum,
        "retained_representatives": [item.get("organism_name") for item in retained],
        "suppressed_to_audit": [item.get("organism_name") for item in suppressed],
    }


def apply_opportunistic_cluster_guardrail(
    candidates: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    if not policy.get("opportunistic_analytical_thresholds"):
        return {
            "candidate_count_before_cap": 0,
            "representative_cap": 0,
            "retained_representatives": [],
            "suppressed_to_analytic_route": [],
        }
    promoted = [
        item for item in candidates
        if item.get("history_policy_family") == "opportunistic_host_risk"
        and item.get("history_route_changed")
    ]
    maximum = int(
        policy["opportunistic_analytical_thresholds"]
        ["max_adjusted_representatives_per_patient"]
    )
    route_priority = {"Priority": 0, "Context": 1, "Audit": 2, "Hold": 3}
    promoted.sort(key=lambda item: (
        route_priority.get(str(item.get("history_adjusted_route")), 9),
        not _analytical_axes(item)["reproducibility_axis"],
        _analytical_axes(item)["best_rank_in_retained_universe"] or 999,
        str(item.get("organism_name") or "").casefold(),
    ))
    retained = promoted[:maximum]
    suppressed = promoted[maximum:]
    for item in suppressed:
        analytic = str(item["analytic_route"])
        item["history_adjusted_route"] = analytic
        item["history_route_changed"] = False
        item["history_rule_ids"] = [
            rule for rule in item["history_rule_ids"]
            if rule not in {"HIST-OPP2-AUDIT-TO-CONTEXT", "HIST-OPP3-CONTEXT-TO-PRIORITY"}
        ] + ["HIST-OPP6-CLUSTER-REPRESENTATIVE-CAP"]
        item["history_adjustment_blockers"].append(
            "opportunistic cluster representative cap reached; retained at analytic route"
        )
        item["consumed_evidence_ids"] = []
        item["history_evidence_consumption"] = []
        if analytic == "Priority":
            decision, tier = "review_high_priority", "High"
        elif analytic == "Context":
            decision, tier = "review_context_needed", "Context"
        elif analytic == "Audit":
            decision, tier = "review_low_specificity", "Audit"
        else:
            decision, tier = "hold", "Hold"
        item["history_routed_decision"] = decision
        item["downstream_route"] = {
            "scorer_review_tier": tier,
            "projected_final_tier": tier,
            "picked_allowed_from_history_alone": False,
            "picked_promotion_blockers": [],
            "note": "Candidate remains at its analytic route after the representative cap.",
        }
    return {
        "candidate_count_before_cap": len(promoted),
        "representative_cap": maximum,
        "retained_representatives": [item.get("organism_name") for item in retained],
        "suppressed_to_analytic_route": [item.get("organism_name") for item in suppressed],
    }


def apply_hap_mdr_cluster_guardrail(
    candidates: list[dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    if not policy.get("hap_mdr_analytical_thresholds"):
        return {
            "candidate_count_before_cap": 0,
            "representative_cap": 0,
            "retained_representatives": [],
            "suppressed_to_analytic_route": [],
        }
    promoted = [
        item for item in candidates
        if item.get("history_policy_family") == "hap_mdr_exposure"
        and item.get("history_route_changed")
    ]
    maximum = int(
        policy["hap_mdr_analytical_thresholds"]
        ["max_adjusted_representatives_per_patient"]
    )
    route_priority = {"Priority": 0, "Context": 1, "Audit": 2, "Hold": 3}
    promoted.sort(key=lambda item: (
        route_priority.get(str(item.get("history_adjusted_route")), 9),
        not _analytical_axes(item)["reproducibility_axis"],
        _analytical_axes(item)["best_rank_in_retained_universe"] or 999,
        str(item.get("organism_name") or "").casefold(),
    ))
    retained = promoted[:maximum]
    suppressed = promoted[maximum:]
    for item in suppressed:
        analytic = str(item["analytic_route"])
        item["history_adjusted_route"] = analytic
        item["history_route_changed"] = False
        item["history_rule_ids"] = [
            rule for rule in item["history_rule_ids"]
            if rule not in {
                "HIST-HAP2-AUDIT-TO-CONTEXT",
                "HIST-HAP3-SETTING-CONTEXT-TO-PRIORITY",
                "HIST-MDR3-EXACT-MATCH-CONTEXT-TO-PRIORITY",
            }
        ] + ["HIST-HAP6-CLUSTER-REPRESENTATIVE-CAP"]
        item["history_adjustment_blockers"].append(
            "HAP/MDR cluster representative cap reached; retained at analytic route"
        )
        item["consumed_evidence_ids"] = []
        item["history_evidence_consumption"] = []
        if analytic == "Priority":
            decision, tier = "review_high_priority", "High"
        elif analytic == "Context":
            decision, tier = "review_context_needed", "Context"
        elif analytic == "Audit":
            decision, tier = "review_low_specificity", "Audit"
        else:
            decision, tier = "hold", "Hold"
        item["history_routed_decision"] = decision
        item["downstream_route"] = {
            "scorer_review_tier": tier,
            "projected_final_tier": tier,
            "picked_allowed_from_history_alone": False,
            "picked_promotion_blockers": [],
            "note": "Candidate remains at its analytic route after the representative cap.",
        }
    return {
        "candidate_count_before_cap": len(promoted),
        "representative_cap": maximum,
        "retained_representatives": [item.get("organism_name") for item in retained],
        "suppressed_to_analytic_route": [item.get("organism_name") for item in suppressed],
    }


def _all_candidates(payload: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = list(payload.get("patient_organism_decisions") or [])
    candidates.extend(payload.get("hospital_only_decisions") or [])
    unique: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        key = canonical_key(candidate.get("organism_name"))
        if key:
            unique.setdefault(key, candidate)
    return list(unique.values())


def write_csv(path: Path, patients: list[dict[str, Any]]) -> None:
    fields = [
        "patient_id", "organism_name", "taxonomy_family", "analytic_route",
        "history_adjusted_route", "history_route_changed", "scorer_review_tier",
        "best_rank", "positive_test_count", "reproducibility_axis",
        "history_policy_family", "history_rule_ids", "history_evidence_ids",
        "consumed_evidence_ids", "future_picked_gate_review_required",
        "future_picked_gate_review_status", "blockers",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for patient in patients:
            for item in patient["candidates"]:
                axes = _analytical_axes(item)
                writer.writerow({
                    "patient_id": patient["patient_id"],
                    "organism_name": item.get("organism_name"),
                    "taxonomy_family": _candidate_family(item),
                    "analytic_route": item["analytic_route"],
                    "history_adjusted_route": item["history_adjusted_route"],
                    "history_route_changed": item["history_route_changed"],
                    "scorer_review_tier": item["downstream_route"]["scorer_review_tier"],
                    "best_rank": axes["best_rank_in_retained_universe"],
                    "positive_test_count": axes["selected_positive_test_count"],
                    "reproducibility_axis": axes["reproducibility_axis"],
                    "history_policy_family": item.get("history_policy_family"),
                    "history_rule_ids": "|".join(item["history_rule_ids"]),
                    "history_evidence_ids": "|".join(item["history_evidence_ids"]),
                    "consumed_evidence_ids": "|".join(item["consumed_evidence_ids"]),
                    "future_picked_gate_review_required": item.get(
                        "future_picked_gate_review_required"
                    ),
                    "future_picked_gate_review_status": item.get(
                        "future_picked_gate_review_status"
                    ),
                    "blockers": "|".join(item["history_adjustment_blockers"]),
                })


def run(
    scorer_root: Path, patient_root: Path, output_dir: Path, *,
    policy_path: Path = DEFAULT_POLICY, phenotype_root: Path | None = None,
    patients: set[int] | None = None,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty directory: {output_dir}")
    policy = read_json(policy_path)
    paths = sorted(
        (scorer_root / "patient_outputs").glob(
            "NGS_patient_*_test_aware_deterministic_shadow.json"
        ),
        key=patient_number,
    )
    if patients is not None:
        paths = [path for path in paths if patient_number(path) in patients]
    if not paths:
        raise ValueError("No selected test-aware patient outputs found")

    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "resolved_policy_snapshot.json", policy)
    patient_output = output_dir / "patient_outputs"
    patient_output.mkdir()
    outputs: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    for source_path in paths:
        patient = patient_number(source_path)
        source = read_json(source_path)
        raw_dir = patient_root / f"NGS_patient_{patient}_json"
        records = structured_history_records(patient, raw_dir)
        history = aspiration_history_profile(records, policy)
        opportunistic_history = opportunistic_history_profile(records, policy)
        hap_mdr_history = hap_mdr_history_profile(records, policy)
        reactivation_history = reactivation_history_profile(records, policy)
        colonization_enabled = bool(policy.get("colonization_history_patterns"))
        colonization_records = (
            structured_colonization_records(patient, raw_dir, records, policy)
            if colonization_enabled else []
        )
        colonization_history = (
            colonization_history_profile(colonization_records, policy)
            if colonization_enabled else None
        )
        disagreement = phenotype_disagreement(patient, phenotype_root)
        opportunistic_disagreement = phenotype_disagreement(
            patient, phenotype_root, "IMMUNOCOMPROMISED"
        )
        hap_disagreement = phenotype_disagreement(
            patient, phenotype_root, "HAP_PHENOTYPE"
        )
        mdr_disagreement = phenotype_disagreement(
            patient, phenotype_root, "PRIOR_MDR_COLONIZATION_OR_INFECTION"
        )
        hap_mdr_disagreement = {
            "hap_phenotype": hap_disagreement,
            "prior_mdr_phenotype": mdr_disagreement,
        } if hap_disagreement or mdr_disagreement else None
        reactivation_disagreement = phenotype_disagreement(
            patient, phenotype_root, "HERPES_REACTIVATION"
        )
        colonization_airway_disagreement = phenotype_disagreement(
            patient, phenotype_root, "TRACHEOSTOMY_OR_CHRONIC_AIRWAY"
        )
        colonization_mdr_disagreement = phenotype_disagreement(
            patient, phenotype_root, "PRIOR_MDR_COLONIZATION_OR_INFECTION"
        )
        colonization_disagreement = {
            "chronic_airway_phenotype": colonization_airway_disagreement,
            "prior_mdr_colonization_or_infection_phenotype": (
                colonization_mdr_disagreement
            ),
        } if colonization_airway_disagreement or colonization_mdr_disagreement else None
        counts[f"aspiration_profile:{history['support_level']}"] += 1
        counts[f"opportunistic_profile:{opportunistic_history['support_level']}"] += 1
        counts[f"hap_mdr_profile:{hap_mdr_history['support_level']}"] += 1
        counts[f"reactivation_profile:{reactivation_history['support_level']}"] += 1
        if colonization_history is not None:
            counts[
                f"colonization_profile:{colonization_history['support_level']}"
            ] += 1
        candidates = [
            adjust_candidate(
                item,
                history,
                policy,
                disagreement,
                opportunistic_history,
                opportunistic_disagreement,
                hap_mdr_history,
                hap_mdr_disagreement,
                reactivation_history,
                reactivation_disagreement,
                colonization_history,
                colonization_disagreement,
            )
            for item in _all_candidates(source)
        ]
        cluster_guardrail = apply_aspiration_cluster_guardrail(candidates, policy)
        opportunistic_cluster_guardrail = apply_opportunistic_cluster_guardrail(
            candidates, policy
        )
        hap_mdr_cluster_guardrail = apply_hap_mdr_cluster_guardrail(candidates, policy)
        for item in candidates:
            counts[f"analytic:{item['analytic_route']}"] += 1
            counts[f"adjusted:{item['history_adjusted_route']}"] += 1
            if item["history_route_changed"]:
                counts[f"transition_family:{item['history_policy_family']}"] += 1
                counts[
                    f"transition:{item['analytic_route']}->{item['history_adjusted_route']}"
                ] += 1
            if item.get("future_picked_gate_review_required"):
                counts["future_picked_gate_review_required"] += 1
        payload = {
            "schema_version": SCHEMA_VERSION,
            "policy_id": policy["policy_id"],
            "answer_blind": True,
            "patient_id": str(patient),
            "aspiration_history_profile": history,
            "opportunistic_history_profile": opportunistic_history,
            "hap_mdr_history_profile": hap_mdr_history,
            "reactivation_history_profile": reactivation_history,
            "colonization_history_profile": colonization_history,
            "phenotype_source_disagreement": disagreement,
            "opportunistic_phenotype_source_disagreement": opportunistic_disagreement,
            "hap_mdr_phenotype_source_disagreement": hap_mdr_disagreement,
            "reactivation_phenotype_source_disagreement": reactivation_disagreement,
            "colonization_phenotype_source_disagreement": colonization_disagreement,
            "aspiration_cluster_guardrail": cluster_guardrail,
            "opportunistic_cluster_guardrail": opportunistic_cluster_guardrail,
            "hap_mdr_cluster_guardrail": hap_mdr_cluster_guardrail,
            "candidate_count": len(candidates),
            "changed_candidates": [item for item in candidates if item["history_route_changed"]],
            "clinical_scorer_forward": [
                item for item in candidates
                if item["history_adjusted_route"] in {"Priority", "Context"}
            ],
            "audit_or_hold": [
                item for item in candidates
                if item["history_adjusted_route"] in {"Audit", "Hold"}
            ],
            "candidates": candidates,
            "source": {
                "test_aware_file": str(source_path.resolve()),
                "test_aware_sha256": sha256_file(source_path),
                "patient_source_root": str(raw_dir.resolve()),
            },
            "constraints": policy["constraints"],
        }
        target = patient_output / f"NGS_patient_{patient}_history_adjusted_route_shadow.json"
        write_json(target, payload)
        outputs.append(payload)
        counts["patients"] += 1
        counts["candidates"] += len(candidates)

    write_csv(output_dir / "history_route_decisions.csv", outputs)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "policy_id": policy["policy_id"],
        "answer_blind": True,
        "patient_count": len(outputs),
        "candidate_count": sum(item["candidate_count"] for item in outputs),
        "counts": dict(sorted(counts.items())),
        "source_scorer_root": str(scorer_root.resolve()),
        "source_patient_root": str(patient_root.resolve()),
        "policy_file": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "resolved_policy_sha256": hashlib.sha256(
            json.dumps(policy, ensure_ascii=False, sort_keys=True).encode("utf-8")
        ).hexdigest(),
        "policy_composition": policy.get("policy_composition"),
        "policy_composition_chain": policy.get("policy_composition_chain") or [],
        "resolved_policy_snapshot": str(
            (output_dir / "resolved_policy_snapshot.json").resolve()
        ),
        "constraints": policy["constraints"],
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scorer_root", type=Path)
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--phenotype-root", type=Path)
    parser.add_argument("--patients", type=int, nargs="*")
    args = parser.parse_args()
    selected = set(args.patients) if args.patients else None
    print(json.dumps(run(
        args.scorer_root,
        args.patient_root,
        args.output_dir,
        policy_path=args.policy,
        phenotype_root=args.phenotype_root,
        patients=selected,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

from __future__ import annotations

import hashlib
import json
import re
from typing import Any, Dict, Iterable, List


SCHEMA_VERSION = "rag_re_casefit.case_card.v1"
SITE_CODES = {
    "lower_respiratory",
    "upper_respiratory",
    "bloodstream",
    "cns",
    "urinary",
    "intra_abdominal",
    "other",
    "unknown",
}
SYNDROME_CODES = {
    "lower_respiratory_infection",
    "upper_respiratory_infection",
    "bloodstream_infection_or_sepsis",
    "cns_infection",
    "urinary_tract_infection",
    "intra_abdominal_infection",
    "other",
    "unknown",
}
SPECIMEN_TYPES = {
    "bal_or_balf",
    "sputum_or_tracheal_aspirate",
    "blood_or_plasma",
    "csf",
    "urine",
    "tissue",
    "other",
    "unknown",
}
ONSET_PATTERNS = {"acute", "subacute", "chronic", "unknown"}
HOST_FACTORS = {
    "neutropenia",
    "hematologic_malignancy",
    "solid_tumor",
    "hematopoietic_stem_cell_transplant",
    "solid_organ_transplant",
    "immunosuppressive_therapy",
    "prolonged_corticosteroid_use",
    "hiv_or_low_cd4",
    "critical_illness_or_icu",
    "structural_lung_disease",
    "none_known",
}
CLINICAL_FEATURES = {
    "fever",
    "leukocytosis_or_leukopenia",
    "purulent_respiratory_secretions",
    "hypoxemia_or_worsening_oxygenation",
    "hemodynamic_instability",
    "neurologic_syndrome",
    "urinary_symptoms",
    "abdominal_signs_or_symptoms",
}
IMAGING_FEATURES = {
    "new_infiltrate",
    "consolidation",
    "ground_glass_opacity",
    "nodules",
    "cavitation",
    "halo_or_reverse_halo_sign",
    "abscess",
}
CASE_CARD_KEYS = {
    "schema_version",
    "patient_id",
    "target_syndrome",
    "clinical_site",
    "specimen_site",
    "specimen_type",
    "host_factors",
    "clinical_features",
    "imaging_features",
    "onset_pattern",
    "data_availability",
    "missing_fields",
    "source_sha256",
}
PROMPT_KEYS = {
    "target_syndrome",
    "clinical_site",
    "specimen_site",
    "specimen_type",
    "host_factors",
    "clinical_features",
    "imaging_features",
    "onset_pattern",
    "missing_fields",
}


def _compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def _normalized(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", _compact(value).casefold())


def _infer_site(value: Any) -> str:
    text = _compact(value).casefold()
    if any(term in text for term in ("upper_resp", "nasopharyn", "nasal", "throat")):
        return "upper_respiratory"
    if any(term in text for term in ("blood", "systemic", "sterile", "sepsis", "bacteremia", "s1_")):
        return "bloodstream"
    if any(term in text for term in ("lower_resp", "respiratory", "pulmonary", "lung", "bal", "sputum", "s2_")):
        return "lower_respiratory"
    if any(term in text for term in ("csf", "cns", "mening", "brain")):
        return "cns"
    if any(term in text for term in ("urinary", "urine", "uti")):
        return "urinary"
    if any(term in text for term in ("abdominal", "peritone", "ascites")):
        return "intra_abdominal"
    return "unknown"


def _syndrome_for_site(site: str) -> str:
    return {
        "lower_respiratory": "lower_respiratory_infection",
        "upper_respiratory": "upper_respiratory_infection",
        "bloodstream": "bloodstream_infection_or_sepsis",
        "cns": "cns_infection",
        "urinary": "urinary_tract_infection",
        "intra_abdominal": "intra_abdominal_infection",
        "other": "other",
    }.get(site, "unknown")


def _specimen_type(candidate: Dict[str, Any]) -> str:
    values: List[str] = []

    def walk(value: Any, key_hint: str = "") -> None:
        if isinstance(value, dict):
            for key, nested in value.items():
                key_text = _compact(key).casefold()
                if any(token in key_text for token in ("specimen", "sample", "collection_site")):
                    values.append(_compact(nested))
                walk(nested, key_text)
        elif isinstance(value, list):
            for nested in value:
                walk(nested, key_hint)

    walk(candidate.get("key_evidence"))
    text = " ".join(values).casefold()
    if "balf" in text or re.search(r"\bbal\b", text):
        return "bal_or_balf"
    if any(term in text for term in ("sputum", "tracheal", "endotracheal")):
        return "sputum_or_tracheal_aspirate"
    if any(term in text for term in ("blood", "plasma", "serum")):
        return "blood_or_plasma"
    if "csf" in text or "cerebrospinal" in text:
        return "csf"
    if "urine" in text:
        return "urine"
    if any(term in text for term in ("tissue", "biopsy", "histopath")):
        return "tissue"
    return "unknown"


def _candidate_name(candidate: Dict[str, Any]) -> str:
    for key in ("organism_name", "pathogen_name", "pathogen", "name", "organism"):
        value = _compact(candidate.get(key))
        if value:
            return value
    return ""


def find_candidate(root: Dict[str, Any], organism: str) -> Dict[str, Any]:
    target = _normalized(organism)
    if not target:
        raise ValueError("organism is required")
    bundle = root.get("deterministic_max") if isinstance(root.get("deterministic_max"), dict) else root
    rows: List[Dict[str, Any]] = []
    values = bundle.get("pathogen_candidates") if isinstance(bundle, dict) else None
    if isinstance(values, list):
        rows.extend(row for row in values if isinstance(row, dict))
    review = root.get("llm_missed_candidate_review")
    if isinstance(review, dict):
        for tier in ("review_high_priority", "review_context_needed", "review_low_specificity"):
            values = review.get(tier)
            if isinstance(values, list):
                rows.extend(row for row in values if isinstance(row, dict))
    matches = [row for row in rows if _normalized(_candidate_name(row)) == target]
    if not matches:
        raise ValueError(f"organism not found in candidate data: {organism}")
    deterministic = [row for row in matches if "key_evidence" in row]
    return deterministic[0] if deterministic else matches[0]


def _validate_enum_list(name: str, value: Any, allowed: Iterable[str]) -> List[str]:
    if not isinstance(value, list):
        raise ValueError(f"{name} must be a list")
    allowed_set = set(allowed)
    normalized = [_compact(item) for item in value]
    if any(not item or item not in allowed_set for item in normalized):
        raise ValueError(f"{name} contains an unsupported value")
    if len(normalized) != len(set(normalized)):
        raise ValueError(f"{name} contains duplicates")
    return normalized


def validate_case_card(card: Dict[str, Any]) -> Dict[str, Any]:
    if not isinstance(card, dict):
        raise ValueError("case card must be an object")
    extra = set(card) - CASE_CARD_KEYS
    missing = CASE_CARD_KEYS - set(card)
    if extra or missing:
        raise ValueError(f"invalid case card keys; missing={sorted(missing)}, extra={sorted(extra)}")
    if card.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("unsupported case card schema_version")
    patient_id = _compact(card.get("patient_id"))
    if not patient_id:
        raise ValueError("patient_id is required in the artifact")
    for key, allowed in (
        ("target_syndrome", SYNDROME_CODES),
        ("clinical_site", SITE_CODES),
        ("specimen_site", SITE_CODES),
        ("specimen_type", SPECIMEN_TYPES),
        ("onset_pattern", ONSET_PATTERNS),
    ):
        if card.get(key) not in allowed:
            raise ValueError(f"unsupported {key}: {card.get(key)}")
    _validate_enum_list("host_factors", card.get("host_factors"), HOST_FACTORS)
    _validate_enum_list("clinical_features", card.get("clinical_features"), CLINICAL_FEATURES)
    _validate_enum_list("imaging_features", card.get("imaging_features"), IMAGING_FEATURES)
    if not isinstance(card.get("data_availability"), dict):
        raise ValueError("data_availability must be an object")
    if not isinstance(card.get("missing_fields"), list) or any(
        not _compact(item) for item in card.get("missing_fields")
    ):
        raise ValueError("missing_fields must be a non-empty-string list")
    source_hash = _compact(card.get("source_sha256"))
    if source_hash and not re.fullmatch(r"[0-9a-f]{64}", source_hash):
        raise ValueError("source_sha256 must be empty or a lowercase SHA-256")
    return card


def prompt_case_card(card: Dict[str, Any]) -> Dict[str, Any]:
    validate_case_card(card)
    # Patient identifier, source hash, upstream tier codes, and provenance never enter
    # the model prompt. The explicit projection also drops any injected answer/gold.
    return {key: card[key] for key in sorted(PROMPT_KEYS)}


def build_case_card(
    root: Dict[str, Any],
    candidate: Dict[str, Any],
    *,
    source_bytes: bytes | None = None,
) -> Dict[str, Any]:
    if not isinstance(root, dict) or not isinstance(candidate, dict):
        raise ValueError("root and candidate must be objects")
    bundle = root.get("deterministic_max") if isinstance(root.get("deterministic_max"), dict) else root
    dominant_source = _compact(bundle.get("dominant_source"))
    explicit_site = " ".join(
        _compact(bundle.get(key))
        for key in ("clinical_target_site", "infection_site", "target_site")
        if _compact(bundle.get(key))
    )
    clinical_site = _infer_site(explicit_site) if explicit_site else _infer_site(dominant_source)
    specimen_site = _infer_site(candidate.get("specimen_class"))
    specimen_type = _specimen_type(candidate)
    host_context = bundle.get("host_context")
    host_context = host_context if isinstance(host_context, dict) else {}
    coarse_host = _compact(host_context.get("host_vulnerability_tier"))

    missing = ["clinical_features", "imaging_features", "onset_pattern"]
    if specimen_type == "unknown":
        missing.append("exact_specimen_type")
    if not coarse_host or coarse_host.casefold() == "unknown":
        missing.append("host_factors")
        host_status = "not_available"
    else:
        # V0-V3 is retained only as an availability audit. It is not interpreted or
        # shown to the LLM because its clinical definition is absent from the READ spec.
        missing.append("raw_host_factors")
        host_status = "coarse_upstream_tier_only"
    if clinical_site == "unknown":
        missing.append("clinical_site")

    if source_bytes is None:
        canonical = json.dumps(root, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        source_bytes = canonical.encode("utf-8")
    card = {
        "schema_version": SCHEMA_VERSION,
        "patient_id": _compact(root.get("patient_id") or bundle.get("patient_id")),
        "target_syndrome": _syndrome_for_site(clinical_site),
        "clinical_site": clinical_site,
        "specimen_site": specimen_site,
        "specimen_type": specimen_type,
        "host_factors": [],
        "clinical_features": [],
        "imaging_features": [],
        "onset_pattern": "unknown",
        "data_availability": {
            "clinical_site_source": "explicit_bundle_field" if explicit_site else "dominant_source" if dominant_source else "unknown",
            "host_context_status": host_status,
            "coarse_host_tier_present_but_not_prompted": bool(coarse_host and coarse_host.casefold() != "unknown"),
        },
        "missing_fields": sorted(set(missing)),
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
    }
    return validate_case_card(card)

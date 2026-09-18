"""Deterministically preserve clinically relevant evidence around LLM agents.

The helpers in this module do not decide pathogen causality or output tier.
They retain source text, assertion, timing, and missingness so downstream
rules can distinguish a real absence of evidence from an extraction loss.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping

from tools import organism_taxonomy_classifier as taxonomy


PRESERVATION_VERSION = "evidence_preservation_v1"


@dataclass(frozen=True)
class ImagingConcept:
    display_name: str
    concept_type: str
    biological_class: str
    patterns: tuple[str, ...]


IMAGING_CONCEPTS = (
    ImagingConcept(
        "Pneumocystis jirovecii",
        "organism_or_disease",
        "fungus",
        (r"\bpjp\b", r"\bpcp\s+(?:pneumonia|pneumonitis)\b", r"\bpneumocystis(?:\s+jirovecii)?\b"),
    ),
    ImagingConcept("Aspergillus spp.", "organism_or_disease", "fungus", (r"\baspergill(?:us|osis|oma)\b",)),
    ImagingConcept("Mucorales", "organism_or_disease", "fungus", (r"\bmucor(?:ales|mycosis)?\b", r"\brhizopus\b")),
    ImagingConcept("Nocardia spp.", "organism_or_disease", "bacterium", (r"\bnocardi(?:a|osis)\b",)),
    ImagingConcept(
        "Mycobacterium tuberculosis",
        "organism_or_disease",
        "bacterium",
        (r"\bmycobacterium\s+tuberculosis\b", r"\bpulmonary\s+(?:tb|tuberculosis)\b", r"\btuberculosis\b"),
    ),
    ImagingConcept("CMV", "organism_or_disease", "virus", (r"\bcmv\b", r"\bcytomegalovirus\b")),
    ImagingConcept("HSV", "organism_or_disease", "virus", (r"\bhsv(?:[-\s]?[12])?\b", r"\bherpes\s+simplex\b")),
    ImagingConcept("VZV", "organism_or_disease", "virus", (r"\bvzv\b", r"\bvaricella(?:[-\s]zoster)?\b")),
    ImagingConcept("SARS-CoV-2", "organism_or_disease", "virus", (r"\bcovid[-\s]?19\b", r"\bsars[-\s]?cov[-\s]?2\b")),
    ImagingConcept("Influenza virus", "organism_or_disease", "virus", (r"\binfluenza(?:\s+[ab])?\b",)),
    ImagingConcept("viral pneumonitis", "syndrome", "virus", (r"\bviral\s+pneumon(?:ia|itis)\b",)),
    ImagingConcept("bacterial pneumonia", "syndrome", "bacterium", (r"\bbacterial\s+pneumonia\b",)),
    ImagingConcept("fungal pneumonia", "syndrome", "fungus", (r"\bfungal\s+pneumon(?:ia|itis)\b",)),
)


MEDICATION_CLASSES: dict[str, tuple[str, tuple[str, ...]]] = {
    "prednisolone": ("systemic_corticosteroid", (r"\bprednisolone\b",)),
    "prednisone": ("systemic_corticosteroid", (r"\bprednisone\b",)),
    "methylprednisolone": ("systemic_corticosteroid", (r"\bmethylprednisolone\b",)),
    "dexamethasone": ("systemic_corticosteroid", (r"\bdexamethasone\b",)),
    "hydrocortisone": ("systemic_corticosteroid", (r"\bhydrocortisone\b",)),
    "methotrexate": ("antimetabolite_immunosuppressant", (r"\bmethotrexate\b",)),
    "azathioprine": ("antimetabolite_immunosuppressant", (r"\bazathioprine\b",)),
    "mycophenolate": ("antimetabolite_immunosuppressant", (r"\bmycophenolate\b",)),
    "cyclophosphamide": ("cytotoxic_immunosuppressant", (r"\bcyclophosphamide\b",)),
    "tacrolimus": ("calcineurin_inhibitor", (r"\btacrolimus\b",)),
    "cyclosporine": ("calcineurin_inhibitor", (r"\bcyclosporine\b", r"\bciclosporin\b")),
    "rituximab": ("biologic_immunosuppressant", (r"\brituximab\b",)),
    "infliximab": ("biologic_immunosuppressant", (r"\binfliximab\b",)),
    "adalimumab": ("biologic_immunosuppressant", (r"\badalimumab\b",)),
    "tocilizumab": ("biologic_immunosuppressant", (r"\btocilizumab\b",)),
}


NEGATION_PATTERNS = (
    r"\bno\s+(?:radiographic\s+)?evidence\s+of\b",
    r"\bnegative\s+for\b",
    r"\bnot\s+(?:compatible|consistent)\s+with\b",
    r"\bunlikely\b",
)
DIFFERENTIAL_PATTERNS = (
    r"\bdifferential\s+diagnosis\b",
    r"\bdifferential\s+includes?\b",
    r"\bcannot\s+(?:be\s+)?exclude(?:d)?\b",
    r"\bconsider\b",
    r"\bsuspect(?:ed)?\b",
    r"\bpossible\b",
    r"\br/?o\b",
)
COMPATIBLE_PATTERNS = (r"\bcompatible\s+with\b", r"\bmay\s+represent\b")
FAVORED_PATTERNS = (
    r"\bfavou?r(?:ed|s)?\b",
    r"\bconsistent\s+with\b",
    r"\bsuggestive\s+of\b",
    r"\blikely\b",
)
SUSCEPTIBILITY_PATTERNS = (
    r"\bsusceptib(?:le|ility)\b",
    r"\bresistan(?:t|ce)\b",
    r"\bminimum\s+inhibitory\s+concentration\b",
    r"\bmic\b",
    r"\bantibiogram\b",
)
PROPHYLAXIS_PATTERNS = (r"\bprophylaxis\b", r"\bprophylactic\b", r"\bprevent(?:ion|ive)\b")
DATE_PATTERN = re.compile(r"(?P<year>20\d{2})[/-](?P<month>1[0-2]|0?[1-9])(?:[/-](?P<day>3[01]|[12]\d|0?[1-9]))?")
DOSE_PATTERN = re.compile(r"\b\d+(?:\.\d+)?\s*(?:mcg|mg|g)(?:\s*/\s*(?:day|d|kg))?\b", re.IGNORECASE)

SOURCE_FIELD_CLASS_TERMS: dict[str, tuple[str, ...]] = {
    "organism_or_target": ("organism", "species", "microbe", "pathogen", "target"),
    "diagnosis_or_finding": ("diagnosis", "finding", "impression", "note"),
    "specimen_or_source": ("specimen", "sample", "source", "site"),
    "date_or_time": ("date", "time", "collected", "reported"),
    "quantity_or_result": (
        "reads",
        "sec_hit",
        "quantity",
        "value",
        "result",
        "count",
        "copies",
        "viral_load",
        "ct",
        "index",
    ),
    "unit": ("unit", "item_full"),
    "interpretation_or_status": ("interpretation", "status", "detected", "positive", "negative"),
}


def _as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _iter_text_leaves(value: Any, path: str = "$") -> Iterable[tuple[str, str]]:
    if isinstance(value, dict):
        for key, item in value.items():
            yield from _iter_text_leaves(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _iter_text_leaves(item, f"{path}[{index}]")
    elif value is not None:
        text = str(value).strip()
        if text:
            yield path, text


def _iter_leaf_fields(
    value: Any, path: str = "$", key_hint: str = ""
) -> Iterable[tuple[str, str, Any]]:
    if isinstance(value, dict):
        for key, item in value.items():
            child_path = f"{path}.{key}"
            yield from _iter_leaf_fields(item, child_path, str(key).lower())
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _iter_leaf_fields(item, f"{path}[{index}]", key_hint)
    elif value is not None and str(value).strip():
        yield path, key_hint, value


def source_field_class(key: str) -> str | None:
    for field_class, terms in SOURCE_FIELD_CLASS_TERMS.items():
        if any(term in key for term in terms):
            return field_class
    return None


def _detection_assertion(value: Any) -> str:
    text = str(value or "").strip().lower()
    if not text:
        return "unknown"
    if _matches_any(text, NEGATION_PATTERNS) or any(
        token in text for token in ("not detected", "negative", "non-reactive", "nonreactive")
    ):
        return "negative_or_excluded"
    if any(token in text for token in ("detected", "positive", "reactive", "isolated", "growth")):
        return "positive_or_detected"
    if _matches_any(text, DIFFERENTIAL_PATTERNS):
        return "uncertain_or_differential"
    return "reported_without_interpretation"


def _source_row_count(value: Any) -> int:
    if isinstance(value, list):
        return len(value)
    if isinstance(value, dict):
        records = value.get("records")
        return len(records) if isinstance(records, list) else 1
    return 0


def build_source_evidence_index(source_sections: Mapping[str, Any]) -> dict[str, Any]:
    """Retain exact critical source values without copying complete source payloads."""
    sections: list[dict[str, Any]] = []
    for section_name, payload in source_sections.items():
        fields: dict[str, list[dict[str, Any]]] = {
            field_class: [] for field_class in SOURCE_FIELD_CLASS_TERMS
        }
        repeated = Counter()
        assertions: list[dict[str, Any]] = []
        for path, key, value in _iter_leaf_fields(payload, path=str(section_name)):
            field_class = source_field_class(key)
            if field_class is None:
                continue
            fields[field_class].append({"path": path, "value": value})
            if field_class == "organism_or_target":
                normalized = re.sub(r"[^a-z0-9]+", "", str(value).lower())
                if normalized:
                    repeated[normalized] += 1
            if field_class in {"quantity_or_result", "interpretation_or_status"}:
                assertions.append(
                    {
                        "path": path,
                        "value": value,
                        "assertion": _detection_assertion(value),
                    }
                )
        sections.append(
            {
                "section": str(section_name),
                "source_available": payload is not None,
                "row_count": _source_row_count(payload),
                "field_counts": {name: len(values) for name, values in fields.items()},
                "fields": fields,
                "result_assertions": assertions,
                "repeated_organism_or_target": [
                    {"normalized_name": name, "occurrence_count": count}
                    for name, count in sorted(repeated.items())
                    if count > 1
                ],
            }
        )
    return {
        "version": PRESERVATION_VERSION,
        "lossless_for_indexed_field_classes": True,
        "sections": sections,
    }


def _matches_any(text: str, patterns: Iterable[str]) -> bool:
    return any(re.search(pattern, text, re.IGNORECASE) for pattern in patterns)


def imaging_assertion(text: str) -> str:
    if _matches_any(text, NEGATION_PATTERNS):
        return "excluded"
    if _matches_any(text, DIFFERENTIAL_PATTERNS):
        return "differential"
    if _matches_any(text, COMPATIBLE_PATTERNS):
        return "compatible"
    if _matches_any(text, FAVORED_PATTERNS):
        return "favored"
    return "mentioned"


def extract_imaging_diagnostic_mentions(raw_image: Any) -> list[dict[str, Any]]:
    rows = raw_image if isinstance(raw_image, list) else [raw_image] if isinstance(raw_image, dict) else []
    output: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    for row_index, row in enumerate(rows):
        if not isinstance(row, dict):
            continue
        exam_type = str(row.get("exam_type") or row.get("test") or "Unknown")
        collected_time = str(row.get("collected_time") or "Unknown")
        reported_time = str(row.get("reported_time") or "Unknown")
        findings = row.get("findings")
        finding_rows = findings if isinstance(findings, list) else [findings] if findings else []
        for finding_index, finding in enumerate(finding_rows):
            sentence = str(finding or "").strip()
            if not sentence:
                continue
            for concept in IMAGING_CONCEPTS:
                if not _matches_any(sentence, concept.patterns):
                    continue
                assertion = imaging_assertion(sentence)
                dedupe_key = (concept.display_name.lower(), assertion, sentence.lower())
                if dedupe_key in seen:
                    continue
                seen.add(dedupe_key)
                profile = (
                    taxonomy.classify_organism(concept.display_name, biological_class=concept.biological_class)
                    if concept.concept_type == "organism_or_disease"
                    else None
                )
                output.append(
                    {
                        "concept_name": concept.display_name,
                        "concept_type": concept.concept_type,
                        "biological_class": concept.biological_class,
                        "assertion": assertion,
                        "evidence_role": "imaging_context_only",
                        "microbiologic_confirmation": False,
                        "source_sentence": sentence,
                        "exam_type": exam_type,
                        "collected_time": collected_time,
                        "reported_time": reported_time,
                        "source_path": f"image[{row_index}].findings[{finding_index}]",
                        "taxonomy_profile": profile,
                    }
                )
    return output


def _parse_datetime(value: Any) -> datetime | None:
    text = str(value or "").strip()
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d", "%Y/%m/%d", "%Y/%m"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None


def _start_date_near(text: str, offset: int) -> tuple[str, datetime | None, bool]:
    window = text[offset : offset + 80]
    match = DATE_PATTERN.search(window)
    if not match:
        return "Unknown", None, False
    year = int(match.group("year"))
    month = int(match.group("month"))
    day = int(match.group("day") or 1)
    normalized = f"{year:04d}-{month:02d}" + (f"-{day:02d}" if match.group("day") else "")
    after = window[match.end() : match.end() + 2]
    ongoing = "-" in after
    return normalized, datetime(year, month, day), ongoing


def _medication_role(text: str, source_path: str, start: datetime | None, sample: datetime | None) -> str:
    combined = f"{source_path} {text}"
    if _matches_any(combined, SUSCEPTIBILITY_PATTERNS):
        return "susceptibility_only"
    if _matches_any(combined, PROPHYLAXIS_PATTERNS):
        return "prophylaxis"
    if start is not None and sample is not None:
        return "pre_sample_risk_medication" if start <= sample else "post_sample_treatment"
    return "medication_exposure_timing_unknown"


def extract_host_risk_medications(
    sources: Mapping[str, Any], *, sample_time: Any = None
) -> dict[str, Any]:
    sample = _parse_datetime(sample_time)
    medications: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    for source_name, payload in sources.items():
        for source_path, text in _iter_text_leaves(payload, path=source_name):
            for medication, (medication_class, patterns) in MEDICATION_CLASSES.items():
                for pattern in patterns:
                    match = re.search(pattern, text, re.IGNORECASE)
                    if not match:
                        continue
                    start_text, start, ongoing = _start_date_near(text, match.end())
                    role = _medication_role(text, source_path, start, sample)
                    dedupe_key = (medication, role, text.lower())
                    if dedupe_key in seen:
                        break
                    seen.add(dedupe_key)
                    dose_match = DOSE_PATTERN.search(text[max(0, match.start() - 30) : match.end() + 50])
                    medications.append(
                        {
                            "medication_name": medication,
                            "medication_class": medication_class,
                            "role": role,
                            "start_date": start_text,
                            "ongoing_at_source": ongoing,
                            "dose": dose_match.group(0) if dose_match else "Unknown",
                            "source_text": text,
                            "source_path": source_path,
                        }
                    )
                    break

    risk_medications = [
        item
        for item in medications
        if item["role"] in {"pre_sample_risk_medication", "medication_exposure_timing_unknown"}
    ]
    risk_classes = {item["medication_class"] for item in risk_medications}
    if len(risk_classes) >= 2:
        strength = "intermediate"
    elif risk_classes:
        strength = "limited"
    else:
        strength = "none"
    gaps: list[str] = []
    if sample is None and risk_medications:
        gaps.append("mngs_sample_time_missing_for_medication_timing")
    if any(item["medication_class"] == "systemic_corticosteroid" and item["dose"] == "Unknown" for item in risk_medications):
        gaps.append("systemic_corticosteroid_dose_missing")
    if any(item["start_date"] == "Unknown" for item in risk_medications):
        gaps.append("risk_medication_start_date_missing")
    return {
        "sample_time": str(sample_time or "Unknown"),
        "medications": medications,
        "risk_medication_count": len(risk_medications),
        "risk_classes": sorted(risk_classes),
        "host_risk_evidence_strength": strength,
        "data_gaps": gaps,
        "tier_effect": "none_evidence_only",
    }


def cbc_data_status(raw_cbc: Any) -> dict[str, Any]:
    rows = raw_cbc if isinstance(raw_cbc, list) else [raw_cbc] if isinstance(raw_cbc, dict) else []
    item_names: set[str] = set()
    dates: set[str] = set()
    for row in rows:
        if not isinstance(row, dict):
            continue
        item = str(row.get("item") or row.get("item_full") or row.get("test") or "").strip().lower()
        if item:
            item_names.add(item)
        date = str(row.get("reported_time") or row.get("collected_time") or "").strip()
        if date:
            dates.add(date)
    cbc_available = bool(rows)
    differential_terms = ("neut", "lymph", "mono", "eos", "baso", "anc", "alc")
    differential_available = any(any(term in item for term in differential_terms) for item in item_names)
    alc_available = any("alc" in item or "absolute lymph" in item for item in item_names)
    if not cbc_available:
        status = "cbc_missing"
    elif not differential_available:
        status = "cbc_available_differential_missing"
    elif not alc_available:
        status = "cbc_and_differential_available_alc_missing"
    else:
        status = "cbc_and_alc_available"
    return {
        "status": status,
        "cbc_available": cbc_available,
        "differential_available": differential_available,
        "absolute_lymphocyte_count_available": alc_available,
        "reported_times": sorted(dates),
        "observed_item_names": sorted(item_names),
    }


def _merge_dict_rows(existing: Any, additions: Iterable[dict[str, Any]], key_fields: tuple[str, ...]) -> list[dict[str, Any]]:
    output = [dict(item) for item in _as_list(existing) if isinstance(item, dict)]
    seen = {tuple(str(item.get(field) or "").lower() for field in key_fields) for item in output}
    for item in additions:
        key = tuple(str(item.get(field) or "").lower() for field in key_fields)
        if key in seen:
            continue
        seen.add(key)
        output.append(dict(item))
    return output


def enrich_image_agent_output(agent_payload: Any, raw_image: Any) -> dict[str, Any]:
    output = dict(agent_payload) if isinstance(agent_payload, dict) else {}
    mentions = extract_imaging_diagnostic_mentions(raw_image)
    output["diagnostic_mentions"] = _merge_dict_rows(
        output.get("diagnostic_mentions"), mentions, ("concept_name", "assertion", "source_sentence")
    )
    preservation = dict(_as_dict(output.get("evidence_preservation")))
    preservation.update(
        {
            "version": PRESERVATION_VERSION,
            "raw_image_available": raw_image is not None,
            "diagnostic_mention_count": len(output["diagnostic_mentions"]),
            "diagnostic_mentions_are_context_only": True,
        }
    )
    output["evidence_preservation"] = preservation
    return output


def enrich_host_agent_output(
    agent_payload: Any,
    *,
    raw_cbc: Any,
    raw_underlying: Any,
    raw_admission: Any,
    sample_time: Any = None,
) -> dict[str, Any]:
    output = dict(agent_payload) if isinstance(agent_payload, dict) else {}
    medication_profile = extract_host_risk_medications(
        {"underlying": raw_underlying, "admission_diagnosis": raw_admission},
        sample_time=sample_time,
    )
    cbc_status = cbc_data_status(raw_cbc)
    output["structured_host_evidence"] = {
        "version": PRESERVATION_VERSION,
        "medication_profile": medication_profile,
        "cbc_data_status": cbc_status,
    }
    state = dict(_as_dict(output.get("host_state")))
    state.setdefault("immunocompromise_tier", "Unknown")
    state.setdefault("acute_instability_tier", "Unknown")
    state.setdefault("host_vulnerability_tier", "Unknown")
    state.setdefault("expanded_candidate_policy", False)
    state.setdefault("opportunistic_coverage_level", "Unknown")
    flags = [str(item) for item in _as_list(state.get("key_host_flags")) if str(item).strip()]
    if medication_profile["risk_medication_count"] and "immunosuppressant_exposure" not in flags:
        flags.append("immunosuppressant_exposure")
    state["key_host_flags"] = flags
    output["host_state"] = state
    preservation = dict(_as_dict(output.get("evidence_preservation")))
    preservation.update(
        {
            "version": PRESERVATION_VERSION,
            "raw_cbc_available": raw_cbc is not None,
            "raw_underlying_available": raw_underlying is not None,
            "raw_admission_available": raw_admission is not None,
            "host_risk_evidence_strength": medication_profile["host_risk_evidence_strength"],
            "tier_changed_by_preservation": False,
        }
    )
    output["evidence_preservation"] = preservation
    return output


def mngs_sample_time_from_ranked(raw_ranked: Any) -> str:
    records = _as_dict(raw_ranked).get("records") if isinstance(raw_ranked, dict) else None
    for record in _as_list(records):
        if isinstance(record, dict) and record.get("collected_time"):
            return str(record["collected_time"])
    return "Unknown"


def patient_mngs_sample_time(patient_dir: Path) -> str:
    matches = sorted(patient_dir.glob("*mNGS_ranked_candidates.json"))
    if not matches:
        return "Unknown"
    import json

    try:
        payload = json.loads(matches[0].read_text(encoding="utf-8-sig"))
    except (OSError, ValueError):
        return "Unknown"
    return mngs_sample_time_from_ranked(payload)


def enrich_agent_payload(
    category: str,
    agent_payload: Any,
    source_sections: Mapping[str, Any],
    *,
    sample_time: Any = None,
) -> dict[str, Any]:
    normalized = str(category or "").strip().lower()
    sources = {str(key).strip().lower(): value for key, value in source_sections.items()}
    if normalized == "image":
        output = enrich_image_agent_output(agent_payload, sources.get("image"))
    elif normalized == "cbc_other_lab":
        output = enrich_host_agent_output(
            agent_payload,
            raw_cbc=sources.get("cbc"),
            raw_underlying=sources.get("underlying"),
            raw_admission=sources.get("admissiondiagnosis"),
            sample_time=sample_time,
        )
    else:
        output = dict(agent_payload) if isinstance(agent_payload, dict) else {}
    output["source_evidence_index"] = build_source_evidence_index(source_sections)
    return output

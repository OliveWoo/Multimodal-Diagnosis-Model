"""Deterministically preserve clinically relevant evidence around LLM agents.

The helpers in this module do not decide pathogen causality or output tier.
They retain source text, assertion, timing, and missingness so downstream
rules can distinguish a real absence of evidence from an extraction loss.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping

from tools import organism_taxonomy_classifier as taxonomy


PRESERVATION_VERSION = "evidence_preservation_v2_lossless_records"


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
    "decitabine": ("cytotoxic_antineoplastic", (r"\bdecitabine\b",)),
    "pemetrexed": ("cytotoxic_antineoplastic", (r"\bpemetrexed\b", r"\balimta\b")),
    "cisplatin": ("cytotoxic_antineoplastic", (r"\bcisplatin\b",)),
    "etoposide": ("cytotoxic_antineoplastic", (r"\betoposide\b",)),
    "vincristine": ("cytotoxic_antineoplastic", (r"\bvincristine\b", r"\bvincritine\b")),
    "doxorubicin": (
        "cytotoxic_antineoplastic",
        (r"\bdoxorubicin\b", r"\blipo[-\s]?doxo\b"),
    ),
    "gemcitabine": ("cytotoxic_antineoplastic", (r"\bgemcitabine\b",)),
    "cytarabine": ("cytotoxic_antineoplastic", (r"\bcytarabine\b",)),
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
DOSE_PATTERN = re.compile(
    r"\b\d+(?:\.\d+)?(?:\s*[-–~]\s*\d+(?:\.\d+)?)?\s*"
    r"(?:mcg|mg|g)(?:\s*/\s*(?:m2|day|d|kg))?\b",
    re.IGNORECASE,
)

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


HOST_RISK_PATTERNS: dict[str, tuple[tuple[str, str], ...]] = {
    "H3": (
        (
            "hematologic_malignancy",
            r"\b(?:hematologic\s+malignancy|leukemia|leukaemia|lymphoma|multiple\s+myeloma|myelodysplastic\s+syndrome|mds)\b",
        ),
        (
            "transplant_history",
            r"\b(?:solid\s+organ\s+transplant|stem\s+cell\s+transplant|hsct|kidney\s+transplant(?:ation)?|renal\s+transplant|liver\s+transplant|heart\s+transplant|lung\s+transplant)\b",
        ),
        ("chemotherapy_exposure", r"\bchemotherapy\b|\bc\s*/\s*t\s*\(\s*chemotherapy\s*\)"),
        ("biologic_exposure", r"\bbiologic(?:al)?\s+(?:therapy|treatment)\b"),
        ("severe_neutropenia_history", r"\bsevere\s+neutropenia\b|\banc\s*[<≤]\s*500\b"),
        ("profound_immune_deficiency", r"\baids\b|\bcd4\s*[<≤]\s*200\b"),
    ),
    "H2": (
        (
            "advanced_ckd",
            r"\badvanced\s+ckd\b|\bckd(?:\s*\([^)]*\))?\s+(?:stage\s*)?[45]\b",
        ),
        ("end_stage_renal_disease", r"\bend[-\s]?stage\s+renal\s+disease\b|\besrd\b"),
        ("dialysis", r"\b(?:hemo)?dialysis\b|\bregular\s+hd\b"),
        ("decompensated_cirrhosis", r"\bdecompensated\s+(?:liver\s+)?cirrhosis\b"),
        ("immunosuppressive_therapy", r"\bimmunosuppressive\s+(?:therapy|treatment)\b"),
    ),
    "H1": (
        ("diabetes", r"\bdiabetes(?:\s+mellitus)?\b|\bdm\b"),
        ("copd", r"\bcopd\b|\bchronic\s+obstructive\s+pulmonary\s+disease\b"),
        ("bronchiectasis", r"\bbronchiectasis\b"),
        ("old_tuberculosis", r"\b(?:old|previous|treated)\s+(?:tb|tuberculosis)\b"),
        ("chronic_lung_disease", r"\bchronic\s+lung\s+disease\b"),
    ),
}

HOST_UNCERTAINTY_PATTERNS = (
    r"\brule\s+out\b",
    r"\br/?o\b",
    r"\bsuspect(?:ed)?\b",
    r"\bpossible\b",
    r"\bto\s+exclude\b",
    r"\bneed\s+to\s+exclude\b",
    r"\bcause\s+to\s+be\s+determined\b",
)
HOST_NEGATION_PATTERNS = (
    r"\bno\s+(?:known\s+)?history\s+of\b",
    r"\bwithout\b",
    r"\bden(?:y|ies|ied)\b",
)


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


def _iter_source_records(value: Any, path: str) -> Iterable[tuple[str, Any]]:
    """Yield normalized records without breaking relationships between fields."""
    if isinstance(value, list):
        for index, item in enumerate(value):
            yield f"{path}[{index}]", item
        return
    if isinstance(value, dict):
        records = value.get("records")
        if isinstance(records, list):
            metadata = {key: item for key, item in value.items() if key != "records"}
            if metadata:
                yield f"{path}.metadata", metadata
            for index, item in enumerate(records):
                yield f"{path}.records[{index}]", item
            return
        yield path, value
        return
    if value is not None:
        yield path, {"value": value}


def build_source_record_index(source_sections: Mapping[str, Any]) -> dict[str, Any]:
    """Preserve complete normalized records for provenance and later re-extraction.

    This is intentionally separate from ``build_source_evidence_index``. The
    field-class index is compact and searchable, while this record index keeps
    all normalized fields from the same source row together. It does not assign
    evidence strength or change tiers.
    """
    sections: list[dict[str, Any]] = []
    for section_name, payload in source_sections.items():
        records = [
            {"record_path": record_path, "record": deepcopy(record)}
            for record_path, record in _iter_source_records(payload, str(section_name))
        ]
        sections.append(
            {
                "section": str(section_name),
                "source_available": payload is not None,
                "record_count": len(records),
                "records": records,
            }
        )
    return {
        "version": PRESERVATION_VERSION,
        "lossless_for_normalized_records": True,
        "tier_effect": "none_evidence_only",
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
    seen: set[tuple[str, str, str, str]] = set()
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


def _parsed_date(match: re.Match[str]) -> tuple[str, datetime]:
    year = int(match.group("year"))
    month = int(match.group("month"))
    day = int(match.group("day") or 1)
    normalized = f"{year:04d}-{month:02d}" + (
        f"-{day:02d}" if match.group("day") else ""
    )
    return normalized, datetime(year, month, day)


def _medication_dates(segment: str) -> tuple[str, datetime | None, str, datetime | None, bool]:
    matches = list(DATE_PATTERN.finditer(segment))
    if not matches:
        return "Unknown", None, "Unknown", None, False
    start_text, start = _parsed_date(matches[0])
    end_text = "Unknown"
    end: datetime | None = None
    if len(matches) > 1:
        end_text, end = _parsed_date(matches[-1])
    tail = segment[matches[-1].end() : matches[-1].end() + 4]
    prefix = segment[max(0, matches[0].start() - 30) : matches[0].start()]
    open_ended = bool(re.match(r"\s*[-–~](?!\s*20\d{2})", tail))
    ongoing = end is None and (open_ended or bool(re.search(r"\bsince\b", prefix, re.I)))
    return start_text, start, end_text, end, ongoing


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
            matches: list[tuple[int, int, str, str]] = []
            for medication, (medication_class, patterns) in MEDICATION_CLASSES.items():
                medication_matches = [
                    match
                    for pattern in patterns
                    for match in re.finditer(pattern, text, re.IGNORECASE)
                ]
                if medication_matches:
                    match = min(medication_matches, key=lambda item: item.start())
                    matches.append(
                        (match.start(), match.end(), medication, medication_class)
                    )
            matches.sort(key=lambda item: (item[0], item[2]))
            for index, (start_offset, end_offset, medication, medication_class) in enumerate(matches):
                next_offset = matches[index + 1][0] if index + 1 < len(matches) else len(text)
                # Keep a drug's dose/date context local so a following drug's
                # dose is not silently attributed to the preceding drug.
                segment = text[start_offset:next_offset]
                start_text, start, end_text, end, ongoing = _medication_dates(segment)
                role = _medication_role(text, source_path, start, sample)
                dedupe_key = (medication, role, source_path.lower(), text.lower())
                if dedupe_key in seen:
                    continue
                seen.add(dedupe_key)
                dose_match = DOSE_PATTERN.search(segment)
                days_from_start = (
                    (sample - start).days
                    if sample is not None and start is not None
                    else None
                )
                days_from_end = (
                    (sample - end).days
                    if sample is not None and end is not None
                    else None
                )
                if ongoing and sample is not None and start is not None and start <= sample:
                    temporal_status = "active_at_sample_from_open_ended_record"
                elif end is not None and sample is not None and end <= sample:
                    temporal_status = "ended_before_sample"
                elif start is not None and sample is not None and start <= sample:
                    temporal_status = "started_before_sample_end_unknown"
                elif start is not None and sample is not None and start > sample:
                    temporal_status = "started_after_sample"
                else:
                    temporal_status = "timing_unknown"
                medications.append(
                    {
                        "medication_name": medication,
                        "medication_class": medication_class,
                        "role": role,
                        "start_date": start_text,
                        "end_date": end_text,
                        "ongoing_at_source": ongoing,
                        "temporal_status": temporal_status,
                        "days_from_start_to_sample": days_from_start,
                        "days_from_end_to_sample": days_from_end,
                        "dose": dose_match.group(0) if dose_match else "Unknown",
                        "source_text": text,
                        "source_path": source_path,
                    }
                )

    risk_medication_observations = [
        item
        for item in medications
        if item["role"] in {"pre_sample_risk_medication", "medication_exposure_timing_unknown"}
    ]
    distinct_risk_medications = sorted(
        {item["medication_name"] for item in risk_medication_observations}
    )
    risk_classes = {item["medication_class"] for item in risk_medication_observations}
    pre_sample = [
        item
        for item in risk_medication_observations
        if item["role"] == "pre_sample_risk_medication"
    ]
    pre_sample_classes = {item["medication_class"] for item in pre_sample}
    recent_cytotoxic = any(
        item["medication_class"] == "cytotoxic_antineoplastic"
        and (
            item["temporal_status"] == "active_at_sample_from_open_ended_record"
            or (
                isinstance(
                    item["days_from_end_to_sample"]
                    if item["days_from_end_to_sample"] is not None
                    else item["days_from_start_to_sample"],
                    int,
                )
                and 0
                <= (
                    item["days_from_end_to_sample"]
                    if item["days_from_end_to_sample"] is not None
                    else item["days_from_start_to_sample"]
                )
                <= 180
            )
        )
        for item in pre_sample
    )
    active_combination = len(pre_sample_classes) >= 2 and any(
        item["temporal_status"] == "active_at_sample_from_open_ended_record"
        for item in pre_sample
    )
    if recent_cytotoxic or active_combination:
        strength = "strong"
    elif pre_sample or len(risk_classes) >= 2:
        strength = "intermediate"
    elif risk_classes:
        strength = "limited"
    else:
        strength = "none"
    gaps: list[str] = []
    if sample is None and risk_medication_observations:
        gaps.append("mngs_sample_time_missing_for_medication_timing")
    if any(
        item["medication_class"] == "systemic_corticosteroid" and item["dose"] == "Unknown"
        for item in risk_medication_observations
    ):
        gaps.append("systemic_corticosteroid_dose_missing")
    if any(item["start_date"] == "Unknown" for item in risk_medication_observations):
        gaps.append("risk_medication_start_date_missing")
    return {
        "sample_time": str(sample_time or "Unknown"),
        "medications": medications,
        "medication_observation_count": len(medications),
        "risk_medication_observation_count": len(risk_medication_observations),
        "risk_medication_count": len(distinct_risk_medications),
        "distinct_risk_medications": distinct_risk_medications,
        "risk_classes": sorted(risk_classes),
        "pre_sample_risk_classes": sorted(pre_sample_classes),
        "verified_pre_sample_observation_count": len(pre_sample),
        "host_risk_evidence_strength": strength,
        "data_gaps": gaps,
        "tier_effect": "none_evidence_only",
    }


def _host_assertion(text: str, match_start: int) -> str:
    local_prefix = text[max(0, match_start - 45) : match_start]
    if _matches_any(local_prefix, HOST_NEGATION_PATTERNS):
        return "negated"
    if _matches_any(local_prefix, HOST_UNCERTAINTY_PATTERNS):
        return "uncertain"
    return "recorded"


def host_history_observations(raw_underlying: Any, raw_admission: Any) -> list[dict[str, Any]]:
    """Preserve host-history text with stable source paths and assertions."""

    observations: list[dict[str, Any]] = []
    for source_name, payload in (
        ("Underlying", raw_underlying),
        ("AdmissionDiagnosis", raw_admission),
    ):
        for source_path, text in _iter_text_leaves(payload, path=source_name):
            path_key = source_path.casefold()
            if not any(
                token in path_key
                for token in ("underlying_diseases", ".diagnosis", ".notes", ".age")
            ):
                continue
            observations.append(
                {
                    "observation_id": f"HOST-{len(observations) + 1:04d}",
                    "source_section": source_name,
                    "source_path": source_path,
                    "source_text": text,
                    "evidence_role": "host_context_not_pathogen_identity",
                }
            )
    return observations


def _underlying_age(raw_underlying: Any) -> int | None:
    rows = raw_underlying if isinstance(raw_underlying, list) else [raw_underlying]
    for row in rows:
        if not isinstance(row, dict):
            continue
        match = re.search(r"\d+", str(row.get("age") or ""))
        if match:
            return int(match.group(0))
    return None


def deterministic_host_risk_profile(
    raw_underlying: Any,
    raw_admission: Any,
    *,
    sample_time: Any = None,
    medication_profile: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Apply the fixed H-tier dictionary without inferring pathogen identity."""

    observations = host_history_observations(raw_underlying, raw_admission)
    triggers: list[dict[str, Any]] = []
    uncertain_mentions: list[dict[str, Any]] = []
    for observation in observations:
        text = observation["source_text"]
        for tier in ("H3", "H2", "H1"):
            for rule, pattern in HOST_RISK_PATTERNS[tier]:
                for match in re.finditer(pattern, text, re.IGNORECASE):
                    assertion = _host_assertion(text, match.start())
                    item = {
                        "tier": tier,
                        "rule": rule,
                        "assertion": assertion,
                        "matched_text": match.group(0),
                        "observation_id": observation["observation_id"],
                        "source_section": observation["source_section"],
                        "source_path": observation["source_path"],
                        "source_text": text,
                    }
                    if assertion == "recorded":
                        triggers.append(item)
                    else:
                        uncertain_mentions.append(item)

    age = _underlying_age(raw_underlying)
    if age is not None and age >= 65:
        triggers.append(
            {
                "tier": "H1",
                "rule": "age_ge_65",
                "assertion": "recorded",
                "matched_text": str(age),
                "observation_id": next(
                    (
                        item["observation_id"]
                        for item in observations
                        if item["source_path"].casefold().endswith(".age")
                    ),
                    None,
                ),
                "source_section": "Underlying",
                "source_path": "Underlying.age",
                "source_text": str(age),
            }
        )

    medication = dict(medication_profile or {})
    medication_strength = str(medication.get("host_risk_evidence_strength") or "none")
    medication_trigger: dict[str, Any] | None = None
    if medication_strength == "strong":
        medication_trigger = {
            "tier": "H3",
            "rule": "verified_recent_or_combination_immunosuppressive_medication",
            "assertion": "recorded",
            "matched_text": "|".join(medication.get("distinct_risk_medications") or []),
            "observation_id": None,
            "source_section": "MedicationProfile",
            "source_path": "structured_host_evidence.medication_profile",
            "source_text": "pre-sample timing plus recent cytotoxic or active combination exposure",
        }
    elif medication_strength == "intermediate":
        pre_sample = [
            item
            for item in medication.get("medications") or []
            if item.get("role") == "pre_sample_risk_medication"
        ]
        prolonged_unknown_dose_steroid = any(
            item.get("medication_class") == "systemic_corticosteroid"
            and item.get("dose") == "Unknown"
            and isinstance(item.get("days_from_start_to_sample"), int)
            and item["days_from_start_to_sample"] >= 28
            for item in pre_sample
        )
        if prolonged_unknown_dose_steroid:
            medication_trigger = {
                "tier": "H2",
                "rule": "prolonged_steroid_duration_known_dose_unknown",
                "assertion": "provisional",
                "matched_text": "systemic corticosteroid",
                "observation_id": None,
                "source_section": "MedicationProfile",
                "source_path": "structured_host_evidence.medication_profile",
                "source_text": "duration >=28 days but dose is unavailable",
            }
    if medication_trigger:
        triggers.append(medication_trigger)

    priority = {"H3": 3, "H2": 2, "H1": 1, "H0": 0}
    tier = max(
        (str(item["tier"]) for item in triggers),
        key=lambda value: priority[value],
        default="H0" if observations else "Unknown",
    )
    tier_triggers = [item for item in triggers if item["tier"] == tier]
    provisional_only = bool(tier_triggers) and all(
        item.get("assertion") == "provisional" for item in tier_triggers
    )
    host_support_eligible = tier in {"H2", "H3"} and not provisional_only
    if host_support_eligible:
        support_status = "support"
    elif provisional_only or medication.get("risk_medication_count"):
        support_status = "insufficient"
    elif tier == "Unknown":
        support_status = "not_available"
    else:
        support_status = "insufficient"
    return {
        "profile_version": "deterministic_host_risk_v1",
        "sample_time": str(sample_time or "Unknown"),
        "immunocompromise_tier": tier,
        "trigger_count": len(triggers),
        "tier_triggers": tier_triggers,
        "all_triggers": triggers,
        "uncertain_or_negated_mentions": uncertain_mentions,
        "host_history_observations": observations,
        "host_support_eligible": host_support_eligible,
        "opportunistic_host_support_status": support_status,
        "medication_evidence_strength": medication_strength,
        "missing_is_negative": False,
        "pathogen_identity_effect": "none",
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
    differential_terms = (
        "neut",
        "lymph",
        "mono",
        "eos",
        "baso",
        "anc",
        "alc",
        "seg",
        "band",
        "metamyelo",
        "myelo",
        "promyelo",
        "blast",
    )
    differential_available = any(any(term in item for term in differential_terms) for item in item_names)
    alc_available = any("alc" in item or "absolute lymph" in item for item in item_names)
    wbc_available = any(
        item == "wbc" or "white blood cell" in item or "leukocyte" in item
        for item in item_names
    )
    lymphocyte_percentage_available = any(
        "lymph" in item and "absolute" not in item and "alc" not in item
        for item in item_names
    )
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
        "wbc_available": wbc_available,
        "lymphocyte_percentage_available": lymphocyte_percentage_available,
        "alc_potentially_derivable": (
            not alc_available and wbc_available and lymphocyte_percentage_available
        ),
        "reported_times": sorted(dates),
        "observed_item_names": sorted(item_names),
    }


def _lab_rows(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    return [value] if isinstance(value, dict) else []


def _lab_item(row: Mapping[str, Any]) -> str:
    return str(row.get("item") or row.get("test_or_item") or row.get("item_full") or "").strip()


def _lab_unit(row: Mapping[str, Any]) -> str:
    explicit = str(row.get("unit") or "").strip()
    if explicit:
        return explicit
    item_full = str(row.get("item_full") or "").strip()
    matches = re.findall(r"\(([^()]*)\)", item_full)
    return matches[-1].strip() if matches else ""


def _lab_numeric(value: Any) -> float | None:
    text = str(value or "").strip().replace(",", "")
    match = re.search(r"[-+]?\d+(?:\.\d+)?", text)
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def _lab_item_key(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").casefold())


def _lab_qc_flags(item: str, value: Any, unit: str) -> list[str]:
    """Flag impossible/ambiguous values without deleting the source row."""
    key = _lab_item_key(item)
    numeric = _lab_numeric(value)
    if numeric is None:
        return ["non_numeric_or_qualitative_value"]
    flags: list[str] = []
    if key == "ph" and not 6.5 <= numeric <= 8.0:
        flags.append("physiologically_implausible_ph")
    elif key in {"temperature", "temp"} and not 25 <= numeric <= 45:
        flags.append("physiologically_implausible_temperature")
    elif key in {"na", "sodium"} and not 80 <= numeric <= 200:
        flags.append("physiologically_implausible_sodium")
    elif key in {"k", "potassium"} and not 1 <= numeric <= 10:
        flags.append("physiologically_implausible_potassium")
    elif key in {"creatinine", "crea", "cr"} and "mg/dl" in unit.casefold() and numeric > 30:
        flags.append("physiologically_implausible_creatinine_mg_dl")
    elif key in {"biltotal", "totalbilirubin", "bilirubintotal"} and "mg/dl" in unit.casefold() and numeric > 50:
        flags.append("physiologically_implausible_bilirubin_mg_dl")
    elif key in {"hco3", "hco3serum", "bicarbonate"} and not 0 <= numeric <= 60:
        flags.append("physiologically_implausible_bicarbonate")
    elif key in {"wbc", "whitebloodcell", "leukocyte"} and not 0 <= numeric <= 500:
        flags.append("physiologically_implausible_wbc")
    elif key in {"plt", "platelet", "platelets"} and not 0 <= numeric <= 2000:
        flags.append("physiologically_implausible_platelet")
    return flags


def _lab_source_record_key(section: str, index: int, row: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        dict(row), ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    digest = hashlib.sha256(encoded).hexdigest()[:16]
    return f"{section}[{index}]#{digest}"


def lab_observations_from_sources(raw_cbc: Any, raw_other_lab: Any) -> list[dict[str, Any]]:
    """Create one deterministic, non-pathogen-specific observation per lab row."""
    observations: list[dict[str, Any]] = []
    for section, source in (("CBC", raw_cbc), ("OtherLab", raw_other_lab)):
        for source_index, row in enumerate(_lab_rows(source)):
            item = _lab_item(row)
            unit = _lab_unit(row)
            value = row.get("value")
            qc_flags = _lab_qc_flags(item, value, unit)
            observations.append(
                {
                    "observation_id": f"LAB-{len(observations) + 1:04d}",
                    "source_section": section,
                    "source_record_index": source_index,
                    "source_record_key": _lab_source_record_key(
                        section, source_index, row
                    ),
                    "test_or_item": item,
                    "item_full": str(row.get("item_full") or "").strip(),
                    "value": "" if value is None else str(value).strip(),
                    "numeric_value": _lab_numeric(value),
                    "unit": unit,
                    "collected_time": str(row.get("collected_time") or "").strip(),
                    "reported_time": str(row.get("reported_time") or "").strip(),
                    "interpretation": "recorded_non_pathogen_specific",
                    "qc_flags": qc_flags,
                    "eligible_for_acute_tiering": not any(
                        flag.startswith("physiologically_implausible_")
                        for flag in qc_flags
                    ),
                    "source_text": json.dumps(
                        row,
                        ensure_ascii=False,
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                }
            )
    return observations


def _count_scale_threshold(value: float, unit: str, *, per_ul: float, kilo_per_ul: float) -> float:
    normalized = unit.casefold().replace(" ", "")
    if any(token in normalized for token in ("x1000/ul", "10^3/ul", "k/ul")):
        return kilo_per_ul
    return per_ul


def acute_instability_profile(raw_cbc: Any, raw_other_lab: Any) -> dict[str, Any]:
    """Recompute the prompt's A-tier from valid raw values for audit and fallback.

    This is a severity/context profile only.  It never supplies pathogen identity
    and never establishes opportunistic-host support.
    """
    observations = lab_observations_from_sources(raw_cbc, raw_other_lab)
    triggers: list[dict[str, Any]] = []
    relevant_count = 0
    for observation in observations:
        if not observation["eligible_for_acute_tiering"]:
            continue
        item = _lab_item_key(observation["test_or_item"])
        value = observation["numeric_value"]
        unit = str(observation["unit"] or "")
        if value is None:
            continue
        tier = None
        rule = ""
        normalized_value = value
        if item in {"lactate", "lacticacid"}:
            relevant_count += 1
            if value >= 4:
                tier, rule = "A3", "lactate_ge_4"
            elif value >= 2:
                tier, rule = "A2", "lactate_2_to_3_9"
        elif item == "ph":
            relevant_count += 1
            if value < 7.25:
                tier, rule = "A3", "ph_lt_7_25"
            elif value <= 7.32:
                tier, rule = "A2", "ph_7_25_to_7_32"
        elif item in {"hco3", "hco3serum", "bicarbonate"}:
            relevant_count += 1
            if value < 18:
                tier, rule = "A3", "hco3_lt_18"
            elif value <= 21:
                tier, rule = "A2", "hco3_18_to_21"
        elif item in {"plt", "platelet", "platelets"}:
            relevant_count += 1
            a3 = _count_scale_threshold(value, unit, per_ul=50000, kilo_per_ul=50)
            a2 = _count_scale_threshold(value, unit, per_ul=100000, kilo_per_ul=100)
            if value < a3:
                tier, rule = "A3", "platelet_lt_50k"
            elif value < a2:
                tier, rule = "A2", "platelet_50k_to_99k"
        elif item in {"anc", "absoluteneutrophilcount"}:
            relevant_count += 1
            threshold = _count_scale_threshold(value, unit, per_ul=500, kilo_per_ul=0.5)
            if value < threshold:
                tier, rule = "A3", "anc_lt_500"
        elif item.startswith("crp") or item in {"creactiveprotein"}:
            relevant_count += 1
            if "mg/dl" in unit.casefold():
                normalized_value = value * 10
            if normalized_value >= 150:
                tier, rule = "A2", "crp_ge_150_mg_l"
            elif normalized_value >= 75:
                tier, rule = "A1", "crp_75_to_149_mg_l"
        elif item in {"wbc", "whitebloodcell", "leukocyte"}:
            relevant_count += 1
            a2 = _count_scale_threshold(value, unit, per_ul=15000, kilo_per_ul=15)
            a1 = _count_scale_threshold(value, unit, per_ul=12000, kilo_per_ul=12)
            if value >= a2:
                tier, rule = "A2", "wbc_ge_15k"
            elif value >= a1:
                tier, rule = "A1", "wbc_12k_to_14_9k"
        if tier:
            triggers.append(
                {
                    "observation_id": observation["observation_id"],
                    "tier": tier,
                    "rule": rule,
                    "test_or_item": observation["test_or_item"],
                    "value": observation["value"],
                    "normalized_value": normalized_value,
                    "unit": unit,
                    "reported_time": observation["reported_time"],
                }
            )
    rank = {"A0": 0, "A1": 1, "A2": 2, "A3": 3}
    tier = max(
        (item["tier"] for item in triggers),
        key=lambda value: rank[value],
        default="A0" if relevant_count else "Unknown",
    )
    return {
        "acute_instability_tier": tier,
        "relevant_observation_count": relevant_count,
        "trigger_count": len(triggers),
        "triggers": triggers,
        "invalid_for_tiering": [
            {
                "observation_id": item["observation_id"],
                "test_or_item": item["test_or_item"],
                "value": item["value"],
                "unit": item["unit"],
                "qc_flags": item["qc_flags"],
            }
            for item in observations
            if not item["eligible_for_acute_tiering"]
        ],
        "role": "severity_context_only_not_pathogen_specific",
        "opportunistic_host_support_effect": "none",
    }


def vulnerability_tier(immunocompromise_tier: Any, acute_instability_tier: Any) -> str:
    """Apply the fixed H/A -> V table from the host-state prompt."""
    host = str(immunocompromise_tier or "Unknown")
    acute = str(acute_instability_tier or "Unknown")
    if host == "H3" or acute == "A3" or (host == "H2" and acute in {"A2", "A3"}):
        return "V3"
    if host == "H2" or acute == "A2" or (host == "H1" and acute == "A2"):
        return "V2"
    if host == "H1" or acute == "A1":
        return "V1"
    if host == "H0" and acute == "A0":
        return "V0"
    return "Unknown"


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
    raw_other_lab: Any = None,
    sample_time: Any = None,
) -> dict[str, Any]:
    output = dict(agent_payload) if isinstance(agent_payload, dict) else {}
    medication_profile = extract_host_risk_medications(
        {"underlying": raw_underlying, "admission_diagnosis": raw_admission},
        sample_time=sample_time,
    )
    host_risk_profile = deterministic_host_risk_profile(
        raw_underlying,
        raw_admission,
        sample_time=sample_time,
        medication_profile=medication_profile,
    )
    cbc_status = cbc_data_status(raw_cbc)
    lab_observations = lab_observations_from_sources(raw_cbc, raw_other_lab)
    acute_profile = acute_instability_profile(raw_cbc, raw_other_lab)
    structured_host = dict(_as_dict(output.get("structured_host_evidence")))
    structured_host.update(
        {
            "version": PRESERVATION_VERSION,
            "medication_profile": medication_profile,
            "host_risk_profile": host_risk_profile,
            "cbc_data_status": cbc_status,
            "lab_observations": lab_observations,
            "lab_source_record_contract": {
                "cbc_source_record_count": len(_lab_rows(raw_cbc)),
                "other_lab_source_record_count": len(_lab_rows(raw_other_lab)),
                "source_record_count": len(_lab_rows(raw_cbc))
                + len(_lab_rows(raw_other_lab)),
                "observation_count": len(lab_observations),
                "one_observation_per_source_record": len(lab_observations)
                == len(_lab_rows(raw_cbc)) + len(_lab_rows(raw_other_lab)),
            },
            "acute_instability_profile": acute_profile,
            "lab_role": "non_pathogen_specific_context_only",
        }
    )
    output["structured_host_evidence"] = structured_host
    state = dict(_as_dict(output.get("host_state")))
    original_state = dict(state)
    state["immunocompromise_tier"] = host_risk_profile[
        "immunocompromise_tier"
    ]
    state["acute_instability_tier"] = acute_profile["acute_instability_tier"]
    state["host_vulnerability_tier"] = vulnerability_tier(
        state["immunocompromise_tier"], state["acute_instability_tier"]
    )
    state["expanded_candidate_policy"] = state["host_vulnerability_tier"] in {
        "V2",
        "V3",
    }
    state["opportunistic_coverage_level"] = {
        "V3": "O2",
        "V2": "O1",
        "V1": "O0",
        "V0": "O0",
    }.get(state["host_vulnerability_tier"], "Unknown")
    state["opportunistic_host_support"] = bool(
        host_risk_profile["host_support_eligible"]
    )
    state["opportunistic_host_support_status"] = host_risk_profile[
        "opportunistic_host_support_status"
    ]
    state["immunocompromise_source"] = "deterministic_source_complete_host_profile"
    state["acute_instability_source"] = "deterministic_source_complete_lab_profile"
    state["acute_severity_is_not_opportunistic_host_support"] = True
    flags = [str(item) for item in _as_list(state.get("key_host_flags")) if str(item).strip()]
    if medication_profile["risk_medication_count"] and "immunosuppressant_exposure" not in flags:
        flags.append("immunosuppressant_exposure")
    trigger_rules = {item["rule"] for item in acute_profile["triggers"]}
    if "anc_lt_500" in trigger_rules and "severe_neutropenia" not in flags:
        flags.append("severe_neutropenia")
    if "platelet_lt_50k" in trigger_rules and "thrombocytopenia" not in flags:
        flags.append("thrombocytopenia")
    if trigger_rules & {"ph_lt_7_25", "hco3_lt_18"} and "severe_acid_base_derangement" not in flags:
        flags.append("severe_acid_base_derangement")
    host_rules = {item["rule"] for item in host_risk_profile["all_triggers"]}
    if host_rules & {
        "verified_recent_or_combination_immunosuppressive_medication",
        "chemotherapy_exposure",
        "biologic_exposure",
        "prolonged_steroid_duration_known_dose_unknown",
    }:
        flags.append("immunosuppressant_exposure")
    if "hematologic_malignancy" in host_rules:
        flags.append("hematologic_malignancy")
    if "transplant_history" in host_rules:
        flags.append("transplant_history")
    if "severe_neutropenia_history" in host_rules:
        flags.append("severe_neutropenia")
    if host_rules & {
        "advanced_ckd",
        "end_stage_renal_disease",
        "dialysis",
        "decompensated_cirrhosis",
    }:
        flags.append("renal_or_hepatic_stress")
    flag_order = (
        "severe_neutropenia",
        "immunosuppressant_exposure",
        "hematologic_malignancy",
        "transplant_history",
        "severe_acid_base_derangement",
        "thrombocytopenia",
        "renal_or_hepatic_stress",
    )
    state["key_host_flags"] = [flag for flag in flag_order if flag in set(flags)]
    output["host_state"] = state
    output["infection_likelihood"] = {
        "A3": "Severe",
        "A2": "Likely",
        "A1": "Possible",
    }.get(state["acute_instability_tier"], "Unknown")
    preservation = dict(_as_dict(output.get("evidence_preservation")))
    preservation.update(
        {
            "version": PRESERVATION_VERSION,
            "raw_cbc_available": raw_cbc is not None,
            "raw_other_lab_available": raw_other_lab is not None,
            "raw_underlying_available": raw_underlying is not None,
            "raw_admission_available": raw_admission is not None,
            "host_risk_evidence_strength": medication_profile["host_risk_evidence_strength"],
            "lab_observation_count": len(lab_observations),
            "original_agent_acute_instability_tier": original_state.get(
                "acute_instability_tier", "Unknown"
            ),
            "original_agent_immunocompromise_tier": original_state.get(
                "immunocompromise_tier", "Unknown"
            ),
            "original_agent_host_vulnerability_tier": original_state.get(
                "host_vulnerability_tier", "Unknown"
            ),
            "tier_changed_by_preservation": bool(
                original_state
                and (
                    original_state.get("acute_instability_tier")
                    != state["acute_instability_tier"]
                    or original_state.get("immunocompromise_tier")
                    != state["immunocompromise_tier"]
                    or original_state.get("host_vulnerability_tier")
                    != state["host_vulnerability_tier"]
                )
            ),
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
    output["source_record_index"] = build_source_record_index(source_sections)
    return output

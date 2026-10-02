"""Audit KH-to-phenotype patient linkage without promoting links to verified.

The audit compares every KH patient with every phenotype ledger. Exact dated
laboratory observations are the primary fingerprint; dated microbiology,
diagnoses, and procedures are independent supporting axes. Patient numbers are
reported but never contribute to the similarity score.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import unicodedata
from collections import Counter
from dataclasses import dataclass, replace
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from tools import pathogen_normalization as pathogen_names


SCHEMA_VERSION = "phenotype_patient_linkage_audit_v1.1"
DATE_RE = re.compile(r"(20\d{2})[-/](\d{1,2})[-/](\d{1,2})")
PATIENT_RE = re.compile(r"(?:NGS_patient_|patient[_ ]?)(\d+)", re.I)

LAB_ALIASES = {
    "wbc": "wbc",
    "rbc": "rbc",
    "hgb": "hgb",
    "hb": "hgb",
    "hemoglobin": "hgb",
    "hct": "hct",
    "mcv": "mcv",
    "mch": "mch",
    "mchc": "mchc",
    "plt": "plt",
    "platelet": "plt",
    "platelets": "plt",
    "rdwcv": "rdwcv",
    "rdwsd": "rdwsd",
    "mpv": "mpv",
    "crp": "crp",
    "pct": "pct",
    "procalcitonin": "pct",
    "na": "na",
    "sodium": "na",
    "k": "k",
    "potassium": "k",
    "sgot": "ast",
    "ast": "ast",
    "sgpt": "alt",
    "alt": "alt",
    "urean": "bun",
    "bun": "bun",
    "creatinine": "creatinine",
    "bil_total": "bilirubin_total",
    "biltotal": "bilirubin_total",
    "bilirubintotal": "bilirubin_total",
    "albumin": "albumin",
    "ldh": "ldh",
    "lactate": "lactate",
    "temperature": "temperature",
    "ph": "ph",
    "pco2": "pco2",
    "po2": "po2",
    "hco3": "hco3",
    "totalco2": "total_co2",
    "beecf": "base_excess",
    "baseexcess": "base_excess",
    "sbc": "sbc",
    "so2c": "so2",
    "fio2": "fio2",
}


@dataclass(frozen=True)
class PatientFeatures:
    patient_number: int
    source_file: str
    labs: frozenset[str]
    microbiology: frozenset[str]
    diagnoses: frozenset[str]
    events: frozenset[str]
    age: str | None = None
    sex: str | None = None


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def patient_number(value: Any) -> int | None:
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().isdigit():
        return int(value.strip())
    match = PATIENT_RE.search(str(value or ""))
    return int(match.group(1)) if match else None


def normalize_date(value: Any) -> str | None:
    match = DATE_RE.search(str(value or ""))
    if not match:
        return None
    year, month, day = (int(part) for part in match.groups())
    try:
        return f"{year:04d}-{month:02d}-{day:02d}"
    except ValueError:
        return None


def normalize_value(value: Any) -> str | None:
    text = unicodedata.normalize("NFKC", str(value or "")).strip()
    match = re.search(r"([<>]=?)?\s*(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", text)
    if not match:
        return None
    comparator = match.group(1) or ""
    try:
        number = Decimal(match.group(2)).normalize()
    except InvalidOperation:
        return None
    rendered = format(number, "f")
    if "." in rendered:
        rendered = rendered.rstrip("0").rstrip(".")
    if rendered in {"-0", ""}:
        rendered = "0"
    return comparator + rendered


def normalize_lab_name(value: Any) -> str | None:
    text = unicodedata.normalize("NFKC", str(value or "")).strip().lower()
    text = re.sub(r"\([^)]*\)", "", text)
    compact = re.sub(r"[^a-z0-9%]+", "", text)
    if compact.startswith("crp"):
        compact = "crp"
    if compact.startswith("bilirubintotal") or compact.startswith("biltotal"):
        compact = "bilirubintotal"
    return LAB_ALIASES.get(compact)


def lab_feature(date: Any, test_name: Any, value: Any) -> str | None:
    normalized_date = normalize_date(date)
    normalized_test = normalize_lab_name(test_name)
    normalized_value = normalize_value(value)
    if not normalized_date or not normalized_test or normalized_value is None:
        return None
    return f"{normalized_date}|{normalized_test}|{normalized_value}"


def normalize_phrase(value: Any) -> str | None:
    text = unicodedata.normalize("NFKC", str(value or "")).lower()
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    text = " ".join(text.split())
    if len(text) < 4 or text in {"pneumonia", "sepsis", "shock", "infection"}:
        return None
    return text


def specimen_family(value: Any) -> str:
    text = str(value or "").lower()
    if any(token in text for token in ("balf", "lower bal", "bronch", "lavage")):
        return "lower_respiratory"
    if "sputum" in text or "endotracheal" in text:
        return "respiratory"
    if "blood" in text or "serum" in text:
        return "blood"
    if "urine" in text or "midstream" in text or "foley" in text:
        return "urine"
    if "stool" in text or "fec" in text:
        return "stool"
    return "other"


def microbiology_feature(date: Any, organism: Any) -> str | None:
    normalized_date = normalize_date(date)
    canonical = pathogen_names.canonical_key(organism)
    if not normalized_date or not canonical:
        return None
    return f"{normalized_date}|{canonical}"


def event_feature(date: Any, event: Any) -> str | None:
    normalized_date = normalize_date(date)
    text = str(event or "").lower()
    if not normalized_date:
        return None
    if any(token in text for token in ("next-generation sequencing", "next generation sequencing", "mngs")) or re.search(
        r"\bngs\b", text
    ):
        kind = "index_mngs_event"
    elif any(token in text for token in ("bronchoscopy", "bronchoscopic", "bronchial washing", "lower_respiratory")):
        kind = "lower_respiratory_sampling"
    elif "admission" in text or "admitted" in text:
        kind = "admission"
    elif "intubat" in text or "endotracheal tube insertion" in text:
        kind = "intubation"
    else:
        return None
    return f"{normalized_date}|{kind}"


def module_path(patient_dir: Path, suffix: str) -> Path:
    base = patient_dir.name.removesuffix("_json")
    return patient_dir / f"{base}_{suffix}.json"


def extract_kh_features(patient_dir: Path) -> PatientFeatures:
    number = patient_number(patient_dir.name)
    if number is None:
        raise ValueError(f"Unexpected KH directory: {patient_dir}")
    labs: set[str] = set()
    microbiology: set[str] = set()
    diagnoses: set[str] = set()
    events: set[str] = set()
    age: str | None = None
    sex: str | None = None

    for suffix in ("CBC", "other_lab"):
        for row in as_list(read_json(module_path(patient_dir, suffix))):
            item = as_dict(row)
            feature = lab_feature(
                item.get("reported_time") or item.get("collected_time"),
                item.get("item") or item.get("item_full"),
                item.get("value"),
            )
            if feature:
                labs.add(feature)

    for row in as_list(read_json(module_path(patient_dir, "underlying"))):
        item = as_dict(row)
        if age is None and item.get("age") is not None:
            age = str(item.get("age")).strip() or None
        if sex is None and item.get("gender") is not None:
            sex = str(item.get("gender")).strip().upper() or None
        for diagnosis in as_list(item.get("underlying_diseases")):
            if normalized := normalize_phrase(diagnosis):
                diagnoses.add(normalized)

    for row in as_list(read_json(module_path(patient_dir, "admission_diagnosis"))):
        if normalized := normalize_phrase(as_dict(row).get("diagnosis")):
            diagnoses.add(normalized)

    for row in as_list(read_json(module_path(patient_dir, "culture"))):
        item = as_dict(row)
        feature = microbiology_feature(
            item.get("collected_time") or item.get("reported_time"), item.get("organism")
        )
        if feature:
            microbiology.add(feature)
        if specimen_family(item.get("sample")) == "lower_respiratory":
            if event := event_feature(item.get("collected_time"), "lower_respiratory_sampling"):
                events.add(event)

    for row in as_list(read_json(module_path(patient_dir, "filmarray"))):
        item = as_dict(row)
        result = str(item.get("result") or "").lower()
        if "detected" in result and "not detected" not in result:
            feature = microbiology_feature(item.get("reported_time"), item.get("target"))
            if feature:
                microbiology.add(feature)
        if specimen_family(item.get("sample")) == "lower_respiratory":
            if event := event_feature(item.get("reported_time"), "lower_respiratory_sampling"):
                events.add(event)

    for row in as_list(read_json(module_path(patient_dir, "gm_test"))):
        item = as_dict(row)
        if specimen_family(item.get("sample")) == "lower_respiratory":
            if event := event_feature(item.get("reported_time"), "lower_respiratory_sampling"):
                events.add(event)

    return PatientFeatures(
        patient_number=number,
        source_file=str(patient_dir),
        labs=frozenset(labs),
        microbiology=frozenset(microbiology),
        diagnoses=frozenset(diagnoses),
        events=frozenset(events),
        age=age,
        sex=sex,
    )


def extract_phenotype_features(path: Path) -> PatientFeatures:
    payload = as_dict(read_json(path))
    number = patient_number(payload.get("patient_id") or path.name)
    if number is None:
        raise ValueError(f"Unexpected phenotype file: {path}")
    records = as_dict(payload.get("records"))
    labs: set[str] = set()
    microbiology: set[str] = set()
    diagnoses: set[str] = set()
    events: set[str] = set()

    for item in as_list(records.get("Lab_Results")):
        row = as_dict(as_dict(item).get("row"))
        if str(row.get("confidence") or "").upper() not in {"HIGH", "MEDIUM", ""}:
            continue
        feature = lab_feature(
            row.get("datetime") or row.get("date"), row.get("test_name"), row.get("value")
        )
        if feature:
            labs.add(feature)

    for item in as_list(records.get("Microbiology")):
        row = as_dict(as_dict(item).get("row"))
        feature = microbiology_feature(row.get("date"), row.get("organism"))
        if feature:
            microbiology.add(feature)
        if specimen_family(row.get("specimen")) == "lower_respiratory":
            if event := event_feature(row.get("date"), "lower_respiratory_sampling"):
                events.add(event)

    for item in as_list(records.get("Clinical_Facts")):
        row = as_dict(as_dict(item).get("row"))
        if str(row.get("fact_type") or "").upper() == "DIAGNOSIS":
            if normalized := normalize_phrase(row.get("description")):
                diagnoses.add(normalized)

    for item in as_list(records.get("Clinical_Timeline")):
        row = as_dict(as_dict(item).get("row"))
        if event := event_feature(
            row.get("date"), f"{row.get('event_type') or ''} {row.get('description') or ''}"
        ):
            events.add(event)

    return PatientFeatures(
        patient_number=number,
        source_file=str(path),
        labs=frozenset(labs),
        microbiology=frozenset(microbiology),
        diagnoses=frozenset(diagnoses),
        events=frozenset(events),
    )


def document_frequency(
    patients: Iterable[PatientFeatures], attribute: str
) -> Counter[str]:
    frequency: Counter[str] = Counter()
    for patient in patients:
        frequency.update(getattr(patient, attribute))
    return frequency


def weighted_overlap(
    left: frozenset[str],
    right: frozenset[str],
    frequency: Mapping[str, int],
    cohort_size: int,
    multiplier: float,
) -> tuple[list[str], float, int]:
    matches = sorted(left & right)
    score = sum(
        multiplier * (math.log((cohort_size + 1) / (frequency.get(feature, 0) + 1)) + 1)
        for feature in matches
    )
    unique = sum(1 for feature in matches if frequency.get(feature, 0) == 1)
    return matches, round(score, 3), unique


def compare_pair(
    kh: PatientFeatures,
    phenotype: PatientFeatures,
    frequencies: Mapping[str, Counter[str]],
    cohort_size: int,
) -> dict[str, Any]:
    lab_matches, lab_score, unique_labs = weighted_overlap(
        kh.labs, phenotype.labs, frequencies["labs"], cohort_size, 1.0
    )
    micro_matches, micro_score, unique_micro = weighted_overlap(
        kh.microbiology,
        phenotype.microbiology,
        frequencies["microbiology"],
        cohort_size,
        4.0,
    )
    diagnosis_matches, diagnosis_score, unique_diagnoses = weighted_overlap(
        kh.diagnoses,
        phenotype.diagnoses,
        frequencies["diagnoses"],
        cohort_size,
        2.0,
    )
    event_matches, event_score, unique_events = weighted_overlap(
        kh.events, phenotype.events, frequencies["events"], cohort_size, 1.5
    )
    index_mngs_event_matches = [
        feature for feature in event_matches if feature.endswith("|index_mngs_event")
    ]
    lower_respiratory_event_matches = [
        feature for feature in event_matches if feature.endswith("|lower_respiratory_sampling")
    ]
    # An explicit NGS date is much more discriminating than a generic admission event.
    event_score = round(event_score + 8.0 * len(index_mngs_event_matches), 3)
    return {
        "phenotype_patient_number": phenotype.patient_number,
        "score": round(lab_score + micro_score + diagnosis_score + event_score, 3),
        "evidence": {
            "exact_dated_lab_matches": lab_matches,
            "exact_dated_lab_match_count": len(lab_matches),
            "phenotype_unique_lab_match_count": unique_labs,
            "exact_dated_microbiology_matches": micro_matches,
            "exact_dated_microbiology_match_count": len(micro_matches),
            "phenotype_unique_microbiology_match_count": unique_micro,
            "diagnosis_matches": diagnosis_matches,
            "diagnosis_match_count": len(diagnosis_matches),
            "phenotype_unique_diagnosis_match_count": unique_diagnoses,
            "dated_event_matches": event_matches,
            "dated_event_match_count": len(event_matches),
            "phenotype_unique_event_match_count": unique_events,
            "explicit_index_mngs_event_matches": index_mngs_event_matches,
            "explicit_index_mngs_event_match_count": len(index_mngs_event_matches),
            "same_date_lower_respiratory_event_matches": lower_respiratory_event_matches,
            "same_date_lower_respiratory_event_match_count": len(lower_respiratory_event_matches),
        },
        "component_scores": {
            "labs": lab_score,
            "microbiology": micro_score,
            "diagnoses": diagnosis_score,
            "events": event_score,
        },
    }


def classify_link(
    *,
    same: dict[str, Any] | None,
    best: dict[str, Any],
    second: dict[str, Any] | None,
    same_identity_cluster: frozenset[int] | None = None,
) -> str:
    if same is None:
        return "phenotype_file_missing"
    evidence = as_dict(same.get("evidence"))
    lab_count = int(evidence.get("exact_dated_lab_match_count") or 0)
    unique_labs = int(evidence.get("phenotype_unique_lab_match_count") or 0)
    micro_count = int(evidence.get("exact_dated_microbiology_match_count") or 0)
    diagnosis_count = int(evidence.get("diagnosis_match_count") or 0)
    event_count = int(evidence.get("dated_event_match_count") or 0)
    margin = float(same.get("score") or 0) - float(second.get("score") or 0) if second else float(same.get("score") or 0)
    same_number = int(same.get("phenotype_patient_number"))
    best_number = int(best.get("phenotype_patient_number"))
    cluster = same_identity_cluster or frozenset({same_number})
    same_is_best = best_number == same_number
    if len(cluster) > 1 and best_number in cluster:
        return "local_same_patient_duplicate_export_cluster"
    if not same_is_best:
        best_evidence = as_dict(best.get("evidence"))
        if int(best_evidence.get("exact_dated_lab_match_count") or 0) >= 5:
            return "local_conflict_best_match_is_different_patient"
        return "local_ambiguous_best_match_is_different_patient"
    if lab_count >= 8 and unique_labs >= 4 and margin >= 8:
        return "local_strong_same_patient_support"
    if lab_count >= 3 and unique_labs >= 1 and margin >= 3:
        return "local_supportive_same_patient_match"
    if micro_count >= 1 and (diagnosis_count >= 1 or event_count >= 1) and margin >= 3:
        return "local_supportive_same_patient_match"
    return "local_insufficient_for_identity_confirmation"


def load_mngs_metadata(path: Path | None) -> dict[int, dict[str, Any]]:
    if path is None or not path.is_file():
        return {}
    payload = as_dict(read_json(path))
    result: dict[int, dict[str, Any]] = {}
    for item in as_list(payload.get("patients")):
        row = as_dict(item)
        number = patient_number(row.get("patient_id"))
        specimens = as_list(row.get("specimens"))
        if number is None:
            continue
        result[number] = {
            "source_mngs_file": row.get("source"),
            "specimens": [
                {
                    "specimen_id": specimen[0] if len(specimen) > 0 else None,
                    "collected_time": specimen[1] if len(specimen) > 1 else None,
                    "specimen_site": specimen[2] if len(specimen) > 2 else None,
                }
                for specimen in specimens
                if isinstance(specimen, list)
            ],
        }
    return result


def attach_index_mngs_events(
    patients: Sequence[PatientFeatures], metadata: Mapping[int, dict[str, Any]]
) -> list[PatientFeatures]:
    """Attach index-event anchors without changing any extracted clinical evidence."""
    enriched: list[PatientFeatures] = []
    for patient in patients:
        events = set(patient.events)
        for specimen in as_list(as_dict(metadata.get(patient.patient_number)).get("specimens")):
            row = as_dict(specimen)
            date = row.get("collected_time")
            if event := event_feature(date, "index mNGS event"):
                events.add(event)
            if specimen_family(row.get("specimen_site")) in {"lower_respiratory", "respiratory"}:
                if event := event_feature(date, "lower_respiratory_sampling"):
                    events.add(event)
        enriched.append(replace(patient, events=frozenset(events)))
    return enriched


def duplicate_export_clusters(
    patients: Sequence[PatientFeatures], *, minimum_exact_labs: int = 20
) -> tuple[list[list[int]], dict[int, frozenset[int]]]:
    """Find phenotype exports with identical, non-trivial dated lab fingerprints."""
    grouped: dict[frozenset[str], list[int]] = {}
    for patient in patients:
        if len(patient.labs) < minimum_exact_labs:
            continue
        grouped.setdefault(patient.labs, []).append(patient.patient_number)
    clusters = sorted(
        (sorted(numbers) for numbers in grouped.values() if len(numbers) > 1),
        key=lambda numbers: numbers[0],
    )
    lookup: dict[int, frozenset[int]] = {}
    for numbers in clusters:
        cluster = frozenset(numbers)
        for number in numbers:
            lookup[number] = cluster
    return clusters, lookup


def episode_anchor_strength(
    match: Mapping[str, Any], *, index_is_lower_respiratory: bool
) -> tuple[int, int]:
    evidence = as_dict(match.get("evidence"))
    return (
        int(evidence.get("explicit_index_mngs_event_match_count") or 0),
        int(evidence.get("same_date_lower_respiratory_event_match_count") or 0)
        if index_is_lower_respiratory
        else 0,
    )


def classify_episode(
    match: Mapping[str, Any] | None,
    *,
    phenotype_available: bool,
    index_is_lower_respiratory: bool,
    ambiguous: bool,
) -> str:
    if not phenotype_available:
        return "phenotype_file_missing"
    if ambiguous:
        return "multiple_exports_share_index_date_anchor"
    evidence = as_dict(as_dict(match).get("evidence"))
    if int(evidence.get("explicit_index_mngs_event_match_count") or 0) > 0:
        return "explicit_index_mngs_event_match"
    if index_is_lower_respiratory and int(
        evidence.get("same_date_lower_respiratory_event_match_count") or 0
    ) > 0:
        return "same_date_lower_respiratory_event_match"
    return "index_infection_episode_unresolved"


def compact_match(match: dict[str, Any] | None) -> dict[str, Any] | None:
    if match is None:
        return None
    evidence = as_dict(match.get("evidence"))
    return {
        "phenotype_patient_number": match.get("phenotype_patient_number"),
        "score": match.get("score"),
        "exact_dated_lab_match_count": evidence.get("exact_dated_lab_match_count"),
        "phenotype_unique_lab_match_count": evidence.get("phenotype_unique_lab_match_count"),
        "exact_dated_microbiology_match_count": evidence.get("exact_dated_microbiology_match_count"),
        "diagnosis_match_count": evidence.get("diagnosis_match_count"),
        "dated_event_match_count": evidence.get("dated_event_match_count"),
        "explicit_index_mngs_event_match_count": evidence.get(
            "explicit_index_mngs_event_match_count"
        ),
        "same_date_lower_respiratory_event_match_count": evidence.get(
            "same_date_lower_respiratory_event_match_count"
        ),
    }


def audit(
    kh_patients: Sequence[PatientFeatures],
    phenotype_patients: Sequence[PatientFeatures],
    mngs_metadata: Mapping[int, dict[str, Any]],
) -> dict[str, Any]:
    kh_patients = attach_index_mngs_events(kh_patients, mngs_metadata)
    cohort_size = len(phenotype_patients)
    frequencies = {
        attribute: document_frequency(phenotype_patients, attribute)
        for attribute in ("labs", "microbiology", "diagnoses", "events")
    }
    phenotype_by_number = {patient.patient_number: patient for patient in phenotype_patients}
    duplicate_clusters, duplicate_cluster_by_patient = duplicate_export_clusters(
        phenotype_patients
    )
    rows: list[dict[str, Any]] = []
    for kh in sorted(kh_patients, key=lambda patient: patient.patient_number):
        comparisons = sorted(
            [
                compare_pair(kh, phenotype, frequencies, cohort_size)
                for phenotype in phenotype_patients
            ],
            key=lambda item: (-float(item["score"]), int(item["phenotype_patient_number"])),
        )
        best = comparisons[0]
        same = next(
            (
                comparison
                for comparison in comparisons
                if comparison["phenotype_patient_number"] == kh.patient_number
            ),
            None,
        )
        second_for_same = next(
            (
                comparison
                for comparison in comparisons
                if comparison["phenotype_patient_number"] != kh.patient_number
            ),
            None,
        )
        identity_cluster = duplicate_cluster_by_patient.get(
            kh.patient_number, frozenset({kh.patient_number})
        )
        status = classify_link(
            same=same,
            best=best,
            second=second_for_same,
            same_identity_cluster=identity_cluster,
        )
        cluster_comparisons = [
            comparison
            for comparison in comparisons
            if int(comparison["phenotype_patient_number"]) in identity_cluster
        ]
        index_is_lower_respiratory = any(
            specimen_family(as_dict(specimen).get("specimen_site"))
            in {"lower_respiratory", "respiratory"}
            for specimen in as_list(
                as_dict(mngs_metadata.get(kh.patient_number)).get("specimens")
            )
        )
        anchor_strengths = [
            (
                episode_anchor_strength(
                    comparison,
                    index_is_lower_respiratory=index_is_lower_respiratory,
                ),
                comparison,
            )
            for comparison in cluster_comparisons
        ]
        best_anchor_strength = max(
            (strength for strength, _ in anchor_strengths), default=(0, 0)
        )
        anchored = [
            comparison
            for strength, comparison in anchor_strengths
            if strength == best_anchor_strength and strength != (0, 0)
        ]
        episode_best = anchored[0] if len(anchored) == 1 else None
        episode_status = classify_episode(
            episode_best,
            phenotype_available=same is not None,
            index_is_lower_respiratory=index_is_lower_respiratory,
            ambiguous=len(anchored) > 1,
        )
        rank = (
            next(
                index
                for index, comparison in enumerate(comparisons, start=1)
                if comparison["phenotype_patient_number"] == kh.patient_number
            )
            if kh.patient_number in phenotype_by_number
            else None
        )
        rows.append(
            {
                "pipeline_patient_number": kh.patient_number,
                "proposed_phenotype_patient_number": kh.patient_number if same else None,
                "phenotype_file_available": same is not None,
                "local_linkage_status": status,
                "patient_identity_assessment": status,
                "infection_episode_assessment": episode_status,
                "phenotype_duplicate_export_cluster": sorted(identity_cluster)
                if same is not None and len(identity_cluster) > 1
                else [],
                "same_number_match_rank": rank,
                "same_number_match": same,
                "best_match": best,
                "best_non_same_number_match": second_for_same,
                "best_episode_match_within_identity_cluster": episode_best,
                "equally_anchored_episode_candidates": [
                    compact_match(comparison) for comparison in anchored
                ],
                "same_number_score_margin_over_best_other": (
                    round(float(same["score"]) - float(second_for_same["score"]), 3)
                    if same and second_for_same
                    else None
                ),
                "top_three_matches": [compact_match(item) for item in comparisons[:3]],
                "kh_available_demographics": {"age": kh.age, "sex": kh.sex},
                "phenotype_demographics_available": False,
                "index_mngs": mngs_metadata.get(kh.patient_number, {}),
                "external_confirmation_required": True,
                "production_link_verified": False,
            }
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "scope": {
            "kh_patient_count": len(kh_patients),
            "phenotype_patient_count": len(phenotype_patients),
            "patient_number_used_in_similarity_score": False,
            "contains_direct_identifiers": False,
            "production_link_verified": False,
        },
        "method": {
            "primary_axis": "exact date + normalized lab name + normalized numeric value",
            "supporting_axes": [
                "exact dated organism",
                "diagnosis phrase",
                "dated admission/bronchoscopy/intubation event",
                "index mNGS collection date and specimen site",
            ],
            "scoring": "inverse-document-frequency weighted overlap across all phenotype patients",
            "duplicate_export_definition": (
                "Two or more phenotype ledgers with at least 20 dated laboratory features and "
                "an identical full dated-laboratory fingerprint."
            ),
            "important_boundary": (
                "Local evidence can support or challenge a proposed link but cannot establish verified identity. "
                "Verified requires an authorized mapping or source-holder confirmation of both patient and index infection episode."
            ),
        },
        "status_counts": dict(Counter(row["local_linkage_status"] for row in rows)),
        "episode_status_counts": dict(
            Counter(row["infection_episode_assessment"] for row in rows)
        ),
        "phenotype_duplicate_export_clusters": duplicate_clusters,
        "patients": rows,
    }


def write_csv_outputs(output_dir: Path, payload: dict[str, Any]) -> None:
    rows = as_list(payload.get("patients"))
    summary_fields = (
        "pipeline_patient_id",
        "patient_identity_assessment",
        "infection_episode_assessment",
        "phenotype_duplicate_export_cluster",
        "phenotype_file_available",
        "same_number_match_rank",
        "best_match_patient_id",
        "best_episode_match_patient_id",
        "equally_anchored_episode_candidate_ids",
        "same_number_score",
        "best_other_score",
        "score_margin",
        "exact_dated_lab_matches",
        "unique_lab_matches",
        "exact_dated_microbiology_matches",
        "diagnosis_matches",
        "dated_event_matches",
        "mngs_specimen_ids",
        "mngs_collected_times",
        "mngs_specimen_sites",
        "production_link_verified",
    )
    summary_rows: list[dict[str, Any]] = []
    confirmation_rows: list[dict[str, Any]] = []
    for item in rows:
        row = as_dict(item)
        same = as_dict(row.get("same_number_match"))
        same_evidence = as_dict(same.get("evidence"))
        other = as_dict(row.get("best_non_same_number_match"))
        episode_best = as_dict(row.get("best_episode_match_within_identity_cluster"))
        specimens = as_list(as_dict(row.get("index_mngs")).get("specimens"))
        summary_rows.append(
            {
                "pipeline_patient_id": f"P{row.get('pipeline_patient_number')}",
                "patient_identity_assessment": row.get("patient_identity_assessment"),
                "infection_episode_assessment": row.get("infection_episode_assessment"),
                "phenotype_duplicate_export_cluster": "; ".join(
                    f"P{number}"
                    for number in as_list(row.get("phenotype_duplicate_export_cluster"))
                ),
                "phenotype_file_available": row.get("phenotype_file_available"),
                "same_number_match_rank": row.get("same_number_match_rank"),
                "best_match_patient_id": (
                    f"P{as_dict(row.get('best_match')).get('phenotype_patient_number')}"
                    if row.get("phenotype_file_available")
                    else ""
                ),
                "best_episode_match_patient_id": (
                    f"P{episode_best.get('phenotype_patient_number')}" if episode_best else ""
                ),
                "equally_anchored_episode_candidate_ids": "; ".join(
                    f"P{as_dict(candidate).get('phenotype_patient_number')}"
                    for candidate in as_list(row.get("equally_anchored_episode_candidates"))
                ),
                "same_number_score": same.get("score"),
                "best_other_score": other.get("score"),
                "score_margin": row.get("same_number_score_margin_over_best_other"),
                "exact_dated_lab_matches": same_evidence.get("exact_dated_lab_match_count"),
                "unique_lab_matches": same_evidence.get("phenotype_unique_lab_match_count"),
                "exact_dated_microbiology_matches": same_evidence.get("exact_dated_microbiology_match_count"),
                "diagnosis_matches": same_evidence.get("diagnosis_match_count"),
                "dated_event_matches": same_evidence.get("dated_event_match_count"),
                "mngs_specimen_ids": "; ".join(str(x.get("specimen_id") or "") for x in specimens),
                "mngs_collected_times": "; ".join(str(x.get("collected_time") or "") for x in specimens),
                "mngs_specimen_sites": "; ".join(str(x.get("specimen_site") or "") for x in specimens),
                "production_link_verified": False,
            }
        )
        confirmation_rows.append(
            {
                "pipeline_patient_id": f"P{row.get('pipeline_patient_number')}",
                "proposed_phenotype_workbook": (
                    f"patient {row.get('proposed_phenotype_patient_number')}_phenotype.xlsx"
                    if row.get("proposed_phenotype_patient_number") is not None
                    else "MISSING"
                ),
                "highest_content_similarity_workbook_unverified": (
                    f"patient {as_dict(row.get('best_match')).get('phenotype_patient_number')}_phenotype.xlsx"
                    if row.get("phenotype_file_available")
                    else ""
                ),
                "mngs_specimen_ids": "; ".join(str(x.get("specimen_id") or "") for x in specimens),
                "mngs_collected_times": "; ".join(str(x.get("collected_time") or "") for x in specimens),
                "mngs_specimen_sites": "; ".join(str(x.get("specimen_site") or "") for x in specimens),
                "local_linkage_status": row.get("local_linkage_status"),
                "duplicate_export_cluster": "; ".join(
                    f"P{number}"
                    for number in as_list(row.get("phenotype_duplicate_export_cluster"))
                ),
                "best_episode_match_workbook": (
                    f"patient {episode_best.get('phenotype_patient_number')}_phenotype.xlsx"
                    if episode_best
                    else "UNRESOLVED" if row.get("phenotype_file_available") else "MISSING"
                ),
                "equally_anchored_episode_workbooks": "; ".join(
                    f"patient {as_dict(candidate).get('phenotype_patient_number')}_phenotype.xlsx"
                    for candidate in as_list(row.get("equally_anchored_episode_candidates"))
                ),
                "infection_episode_assessment": row.get("infection_episode_assessment"),
                "same_patient_confirmed_by_source_holder": "",
                "same_index_infection_episode_confirmed": "",
                "reviewer": "",
                "review_date": "",
                "notes_or_correct_phenotype_id": "",
            }
        )
    for filename, records, fields in (
        ("phenotype_patient_linkage_summary.csv", summary_rows, summary_fields),
        (
            "physician_linkage_confirmation.csv",
            confirmation_rows,
            tuple(confirmation_rows[0]) if confirmation_rows else (),
        ),
    ):
        with (output_dir / filename).open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)


def write_markdown(output_dir: Path, payload: dict[str, Any]) -> None:
    rows = as_list(payload.get("patients"))
    lines = [
        "# KH 與 phenotype 病人身分／感染事件連結 audit",
        "",
        "本 audit 不使用 P 編號參與相似度計分。它分開回答兩個問題：",
        "",
        "1. KH pipeline 與 phenotype ledger 是否可能屬於同一位病人。",
        "2. phenotype ledger 是否對應這次 mNGS 採檢所代表的感染事件。",
        "",
        "本機資料不含可直接核身的識別碼，因此所有結果都只是本地證據評估，"
        "不能標記為 production verified。正式連結仍須資料提供方確認。",
        "",
        "## 方法",
        "",
        "- 病人層：比較完整日期、標準化檢驗名稱與數值，以及日期化微生物、診斷與事件。",
        "- 事件層：比較 index mNGS 採檢日、檢體部位、明確 NGS 紀錄與同日下呼吸道採檢。",
        "- 重複匯出：至少 20 筆且完整日期化 lab fingerprint 完全相同者，標成同病人重複／分事件群組。",
        "- 邊界：短期生命徵象或治療反應不能用來證明病人身分，也不能單獨證明病原。",
        "",
        "## 整體結果",
        "",
    ]
    for status, count in sorted(as_dict(payload.get("status_counts")).items()):
        lines.append(f"- 身分 `{status}`：{count}")
    for status, count in sorted(as_dict(payload.get("episode_status_counts")).items()):
        lines.append(f"- 感染事件 `{status}`：{count}")
    clusters = as_list(payload.get("phenotype_duplicate_export_clusters"))
    if clusters:
        rendered_clusters = "、".join(
            "/".join(f"P{number}" for number in as_list(cluster))
            for cluster in clusters
        )
        lines.extend(["", f"重複／分事件 phenotype 群組：{rendered_clusters}"])
    lines.extend(
        [
            "",
            "| KH | 病人身分評估 | 感染事件評估 | 重複群組 | 同號排名 | 身分最佳 | 事件最佳 | Lab exact | Lab unique | Micro exact |",
            "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for item in rows:
        row = as_dict(item)
        same = as_dict(row.get("same_number_match"))
        evidence = as_dict(same.get("evidence"))
        best = as_dict(row.get("best_match"))
        episode_best = as_dict(row.get("best_episode_match_within_identity_cluster"))
        cluster = "/".join(
            f"P{number}"
            for number in as_list(row.get("phenotype_duplicate_export_cluster"))
        ) or "-"
        lines.append(
            "| P{patient} | {identity} | {episode} | {cluster} | {rank} | {best} | {episode_best} | {labs} | {unique} | {micro} |".format(
                patient=row.get("pipeline_patient_number"),
                identity=row.get("patient_identity_assessment"),
                episode=row.get("infection_episode_assessment"),
                cluster=cluster,
                rank=row.get("same_number_match_rank") or "-",
                best=(
                    f"P{best.get('phenotype_patient_number')}"
                    if best and row.get("phenotype_file_available")
                    else "-"
                ),
                episode_best=(
                    f"P{episode_best.get('phenotype_patient_number')}"
                    if episode_best
                    else "-"
                ),
                labs=evidence.get("exact_dated_lab_match_count") or 0,
                unique=evidence.get("phenotype_unique_lab_match_count") or 0,
                micro=evidence.get("exact_dated_microbiology_match_count") or 0,
            )
        )
    lines.extend(
        [
            "",
            "## 需要資料提供方確認",
            "",
            "特別注意：P9 的 index mNGS 為血液檢體，不能用同日支氣管鏡證明是同一次 mNGS；"
            "P25 的 P25/P26/P27 匯出均有同日呼吸道事件，無法靠日期唯一選出一份；"
            "P29/P31 與 P33/P34 各為相同 lab 指紋的群組，但部分 index 事件仍缺明確對應。",
            "",
            "請使用 `physician_linkage_confirmation.csv` 逐列確認：",
            "",
            "1. 候選 phenotype workbook 是否與 KH 病例為同一位病人。",
            "2. mNGS specimen ID、採檢時間與部位，是否對應 phenotype 所描述的同一次感染事件。",
            "",
            "收到確認前，phenotype 僅可進 shadow input，不得升級為正式 verified mapping，也不得改寫既有 Picked。",
        ]
    )
    (output_dir / "PHENOTYPE_PATIENT_LINKAGE_AUDIT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kh_root", type=Path)
    parser.add_argument("--phenotype-root", type=Path, required=True)
    parser.add_argument("--mngs-sync-report", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    kh_patients = [
        extract_kh_features(path)
        for path in args.kh_root.glob("NGS_patient_*_json")
        if path.is_dir()
    ]
    phenotype_patients = [
        extract_phenotype_features(path)
        for path in args.phenotype_root.glob("patient_*_phenotype_decision_evidence_v1.json")
        if path.is_file()
    ]
    if not kh_patients or not phenotype_patients:
        raise RuntimeError("Both KH and phenotype cohorts must contain patient files")
    payload = audit(
        kh_patients,
        phenotype_patients,
        load_mngs_metadata(args.mngs_sync_report),
    )
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / "phenotype_patient_linkage_audit.json", payload)
    write_csv_outputs(args.output, payload)
    write_markdown(args.output, payload)
    print(json.dumps({
        "schema_version": payload["schema_version"],
        "kh_patient_count": payload["scope"]["kh_patient_count"],
        "phenotype_patient_count": payload["scope"]["phenotype_patient_count"],
        "status_counts": payload["status_counts"],
        "output": str(args.output),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

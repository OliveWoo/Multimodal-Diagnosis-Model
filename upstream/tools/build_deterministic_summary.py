from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Sequence

from tools import mngs_common as mngs
from tools import evidence_preservation
from tools import normalized_agent_fallback
from tools import pathogen_normalization as pathogen_names
from utils import sanitize_filename


AGENT_OUTPUT_DIR_NAME = "agent_outputs"
SUMMARY_OUTPUT_DIR_NAME = "summary_outputs"
DEFAULT_OUTPUT_SUFFIX = "final_summary_with_filmarray_deterministic"
RULE_VERSION = "summary_deterministic_from_agents_v4_evidence_preservation"

LEVEL_RANK = {
    "Level 1": 1,
    "Level 2": 2,
    "Level 3": 3,
    "Level 4": 4,
    "Level 5": 5,
    "Unknown": 99,
    "Not_available": 99,
}

LIKELIHOOD_RANK = {
    "Severe": 1,
    "Likely": 2,
    "Possible": 3,
    "None": 90,
    "Unknown": 99,
    "Not_available": 99,
}

SOURCE_MAP = {
    "blood": "Systemic",
    "bloodstream": "Systemic",
    "systemic": "Systemic",
    "sterile": "Systemic",
    "respiratory": "Respiratory",
    "lung": "Respiratory",
    "urinary": "Urinary",
    "cns": "CNS",
    "gastrointestinal": "GI",
    "gi": "GI",
    "other": "Other",
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build final_summary_with_filmarray deterministically from small agent outputs."
    )
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--patients", nargs="*", help="Optional patient IDs, e.g. 3 10 NGS_patient_3_json.")
    parser.add_argument("--output-suffix", default=DEFAULT_OUTPUT_SUFFIX)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    return parser.parse_args(argv)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else {}


def load_json_any(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def extract_patient_identifier(patient_dir: Path) -> str:
    base_name = sanitize_filename(patient_dir.name)
    return base_name.removesuffix("_json")


def normalize_patient_id(value: str) -> str:
    text = str(value).strip()
    if text.startswith("NGS_patient_"):
        return text.removesuffix("_json")
    return f"NGS_patient_{text}"


def looks_like_patient_dir(path: Path) -> bool:
    return path.is_dir() and (path / AGENT_OUTPUT_DIR_NAME).is_dir()


def collect_patient_dirs(inputs: Sequence[Path], patients: Sequence[str] | None) -> list[Path]:
    requested = {normalize_patient_id(value) for value in patients or []}
    collected: dict[Path, None] = {}
    for entry in inputs:
        path = entry.expanduser()
        if not path.exists() or not path.is_dir():
            continue
        candidates = [path]
        candidates.extend(child for child in path.iterdir() if child.is_dir())
        for candidate in candidates:
            if looks_like_patient_dir(candidate):
                patient_id = extract_patient_identifier(candidate)
                if not requested or patient_id in requested:
                    collected[candidate.resolve()] = None
    return sorted(collected, key=lambda item: extract_patient_identifier(item))


def agent_path(patient_dir: Path, suffix: str) -> Path:
    base = extract_patient_identifier(patient_dir)
    return patient_dir / AGENT_OUTPUT_DIR_NAME / f"{base}_{suffix}.json"


def patient_source_path(patient_dir: Path, suffix: str) -> Path:
    base = extract_patient_identifier(patient_dir)
    return patient_dir / f"{base}_{suffix}.json"


def output_path(patient_dir: Path, suffix: str) -> Path:
    base = extract_patient_identifier(patient_dir)
    return patient_dir / SUMMARY_OUTPUT_DIR_NAME / f"{base}_{suffix}.json"


def level_rank(value: Any) -> int:
    return LEVEL_RANK.get(str(value or "Unknown"), 99)


def likelihood_rank(value: Any) -> int:
    return LIKELIHOOD_RANK.get(str(value or "Unknown"), 99)


def best_level(values: Sequence[Any]) -> str:
    levels = [str(value) for value in values if level_rank(value) < 99]
    if not levels:
        raw_values = {str(value or "").strip() for value in values}
        if "Unknown" in raw_values:
            return "Unknown"
        return "Not_available"
    return min(levels, key=level_rank)


def normalize_source(value: Any) -> str:
    text = str(value or "Unknown").strip()
    lowered = text.lower()
    for key, normalized in SOURCE_MAP.items():
        if key in lowered:
            return normalized
    return text if text and text not in {"NA", "N/A", "-"} else "Unknown"


def canonical_name(value: Any) -> str:
    return mngs.normalize_organism_name(value)


def classification_from_source(value: Any) -> str:
    text = str(value or "Unknown")
    if text in {"Bacterial", "Viral", "Fungal", "Parasitic", "Unknown"}:
        return text
    lowered = text.lower()
    if "virus" in lowered or "viral" in lowered:
        return "Viral"
    if "fung" in lowered or "mold" in lowered or "yeast" in lowered:
        return "Fungal"
    if "bac" in lowered or "bacteria" in lowered:
        return "Bacterial"
    if "parasite" in lowered:
        return "Parasitic"
    return "Unknown"

def lower_text(value: Any) -> str:
    return str(value or "").strip().lower()


TYPICAL_PNEUMONIA_BACTERIA_PREFIXES = (
    "streptococcuspneumoniae",
    "haemophilusinfluenzae",
    "moraxellacatarrhalis",
    "staphylococcusaureus",
    "klebsiella",
    "escherichiacoli",
    "enterobacter",
    "citrobacter",
    "serratia",
    "proteus",
    "morganella",
    "pseudomonasaeruginosa",
    "acinetobacter",
    "stenotrophomonasmaltophilia",
    "burkholderia",
    "achromobacterxylosoxidans",
    "elizabethkingia",
)

CULTURE_L2_EXCLUDED_PREFIXES = (
    "staphylococcusepidermidis",
    "staphylococcushaemolyticus",
    "staphylococcushominis",
    "staphylococcuscapitis",
    "staphylococcuswarneri",
    "staphylococcussimulans",
    "staphylococcuscohnii",
    "staphylococcuscaprae",
    "corynebacterium",
    "cutibacterium",
    "bacillus",
    "micrococcus",
)


def is_typical_pneumonia_bacterium(name: Any) -> bool:
    normalized = canonical_name(name)
    if not normalized:
        return False
    if any(normalized.startswith(prefix) for prefix in CULTURE_L2_EXCLUDED_PREFIXES):
        return False
    if mngs.is_oral_upper_airway_flora_name(name):
        return False
    return any(normalized.startswith(prefix) for prefix in TYPICAL_PNEUMONIA_BACTERIA_PREFIXES)


def adjust_culture_detected_without_quantity_level(
    *,
    module: str,
    organism_name: Any,
    classification: str,
    level: Any,
    key_evidence: Any,
) -> tuple[str, dict[str, Any] | None, bool]:
    adjusted_level = str(level or "Unknown")
    adjusted_key_evidence = dict(key_evidence) if isinstance(key_evidence, dict) else None
    if module != "culture" or not isinstance(adjusted_key_evidence, dict):
        return adjusted_level, adjusted_key_evidence, False
    if classification != "Bacterial":
        return adjusted_level, adjusted_key_evidence, False
    if not is_typical_pneumonia_bacterium(organism_name):
        return adjusted_level, adjusted_key_evidence, False

    specimen_category = str(adjusted_key_evidence.get("specimen_category") or "")
    quantity_tier = str(adjusted_key_evidence.get("quantity_tier") or "")
    growth_purity = lower_text(adjusted_key_evidence.get("growth_purity"))
    contaminant_flag = lower_text(adjusted_key_evidence.get("contaminant_flag"))
    quantitation_status = lower_text(adjusted_key_evidence.get("quantitation_status"))

    lower_resp = specimen_category == "Lower_Respiratory"
    quantity_unknown = quantity_tier in {"Unknown", "Q0", ""} or quantitation_status == "positive_detected_quantity_unknown"
    not_contaminant = contaminant_flag not in {"yes", "true"} and growth_purity not in {"contaminated", "no_growth"}
    if lower_resp and quantity_unknown and not_contaminant and level_rank(adjusted_level) > 2:
        adjusted_key_evidence["quantity_tier"] = "Unknown"
        adjusted_key_evidence["quantitation_status"] = "positive_detected_quantity_unknown"
        adjusted_key_evidence["quantity_interpretation"] = "unknown_not_low"
        adjusted_key_evidence["deterministic_adjustment"] = "respiratory_typical_pneumonia_culture_detected_quantity_unknown_level2"
        return "Level 2", adjusted_key_evidence, True
    return adjusted_level, adjusted_key_evidence, False



LOWER_RESPIRATORY_FILMARRAY_SAMPLE_TERMS = (
    "bal",
    "balf",
    "bronchoalveolar",
    "eta",
    "endotracheal",
    "ta",
    "tracheal",
    "sputum",
    "lower respiratory",
)


def iter_trace_evidence_values(item: dict[str, Any]) -> list[str]:
    values: list[str] = []
    trace = item.get("decision_trace")
    if not isinstance(trace, list):
        return values
    for step in trace:
        if not isinstance(step, dict):
            continue
        evidence = step.get("evidence")
        if not isinstance(evidence, list):
            continue
        for row in evidence:
            if not isinstance(row, dict):
                continue
            for key in ("field", "value"):
                value = row.get(key)
                if value is not None:
                    values.append(str(value))
    return values


def filmarray_detected_positive(source_item: Any) -> bool:
    if not isinstance(source_item, dict):
        return False
    key_evidence = source_item.get("key_evidence") if isinstance(source_item.get("key_evidence"), dict) else {}
    for key in ("detection_status", "result", "interpretation"):
        if lower_text(key_evidence.get(key)) in {"detected", "positive", "detected_positive"}:
            return True
    values = iter_trace_evidence_values(source_item)
    return any(lower_text(value) in {"detected", "positive", "detected_positive"} for value in values)


def filmarray_has_reported_quantity(source_item: Any) -> bool:
    if not isinstance(source_item, dict):
        return False
    key_evidence = source_item.get("key_evidence") if isinstance(source_item.get("key_evidence"), dict) else {}
    for key in ("semiquant_bin", "bin", "copy_per_ml", "copies_per_ml", "quantity", "quantity_tier"):
        value = lower_text(key_evidence.get(key))
        if value and value not in {"unknown", "not_reported", "missing", "none", "not_available", "na", "n/a"}:
            return True
    for value in iter_trace_evidence_values(source_item):
        text = lower_text(value)
        if "copy/ml" in text or "copies/ml" in text or "gc/ml" in text:
            return True
        if "10^" in text or ">=10" in text or ">10" in text:
            return True
    return False


def filmarray_lower_respiratory_sample(module_payload: Any, source_item: Any) -> bool:
    sample_values: list[str] = []
    if isinstance(module_payload, dict):
        assay = module_payload.get("assay_summary") if isinstance(module_payload.get("assay_summary"), dict) else {}
        sample_values.append(str(assay.get("filmarray_sample_type") or ""))
        sample_values.append(str(module_payload.get("probable_source") or ""))
    if isinstance(source_item, dict):
        key_evidence = source_item.get("key_evidence") if isinstance(source_item.get("key_evidence"), dict) else {}
        for key in ("sample", "sample_type", "specimen_type", "specimen_category", "specimen"):
            sample_values.append(str(key_evidence.get(key) or ""))
        values = iter_trace_evidence_values(source_item)
        for index, value in enumerate(values):
            if lower_text(value) in {"sample", "specimen", "filmarray_sample_type"} and index + 1 < len(values):
                sample_values.append(values[index + 1])
            sample_values.append(value)
    haystack = " ".join(lower_text(value) for value in sample_values)
    return any(term in haystack for term in LOWER_RESPIRATORY_FILMARRAY_SAMPLE_TERMS) or "respiratory" in haystack


def build_filmarray_key_evidence(module_payload: Any, source_item: Any) -> dict[str, Any]:
    key_evidence = dict(source_item.get("key_evidence") or {}) if isinstance(source_item, dict) else {}
    if isinstance(module_payload, dict):
        assay = module_payload.get("assay_summary") if isinstance(module_payload.get("assay_summary"), dict) else {}
        key_evidence.setdefault("filmarray_panel", assay.get("filmarray_panel", "Unknown"))
        key_evidence.setdefault("filmarray_sample_type", assay.get("filmarray_sample_type", "Unknown"))
    if isinstance(source_item, dict):
        key_evidence.setdefault("decision_trace", source_item.get("decision_trace"))
    return key_evidence


def adjust_filmarray_detected_without_quantity_level(
    *,
    module: str,
    organism_name: Any,
    classification: str,
    level: Any,
    key_evidence: Any,
    source_item: Any,
    module_payload: Any,
) -> tuple[str, dict[str, Any] | None, bool]:
    adjusted_level = str(level or "Unknown")
    if module != "filmarray_gmtest":
        return adjusted_level, dict(key_evidence) if isinstance(key_evidence, dict) else None, False
    if classification != "Bacterial" or not is_typical_pneumonia_bacterium(organism_name):
        return adjusted_level, dict(key_evidence) if isinstance(key_evidence, dict) else None, False
    if not filmarray_detected_positive(source_item):
        return adjusted_level, dict(key_evidence) if isinstance(key_evidence, dict) else None, False
    if not filmarray_lower_respiratory_sample(module_payload, source_item):
        return adjusted_level, dict(key_evidence) if isinstance(key_evidence, dict) else None, False
    if filmarray_has_reported_quantity(source_item):
        return adjusted_level, dict(key_evidence) if isinstance(key_evidence, dict) else None, False
    if level_rank(adjusted_level) <= 2:
        return adjusted_level, dict(key_evidence) if isinstance(key_evidence, dict) else None, False

    adjusted_key_evidence = build_filmarray_key_evidence(module_payload, source_item)
    adjusted_key_evidence["detection_status"] = "detected"
    adjusted_key_evidence["quantitation_status"] = "positive_detected_quantity_unknown"
    adjusted_key_evidence["quantity_interpretation"] = "unknown_not_low"
    adjusted_key_evidence["deterministic_adjustment"] = "respiratory_typical_pneumonia_filmarray_detected_quantity_unknown_level2"
    return "Level 2", adjusted_key_evidence, True


def organism_label_key(item: Any) -> str:
    if not isinstance(item, dict):
        return ""
    name = canonical_name(item.get("organism_name") or item.get("name") or item.get("target"))
    classification = classification_from_source(item.get("classification"))
    return f"{name}|{classification}" if name else ""


def labels_from_payload(payload: Any) -> list[dict[str, Any]]:
    labels = payload.get("organism_labels") if isinstance(payload, dict) else []
    return [item for item in labels if isinstance(item, dict)] if isinstance(labels, list) else []


def first_non_empty(*values: Any) -> Any:
    for value in values:
        if value not in {None, "", "Unknown", "Not_available", "NA", "N/A", "-"}:
            return value
    return None


def gm_result_is_positive(value: Any) -> bool:
    text = lower_text(value)
    if not text:
        return False
    if any(token in text for token in ("not detected", "negative", "non-reactive", "nonreactive")):
        return False
    return any(token in text for token in ("detected", "positive", "reactive"))


def gm_target_is_aspergillus(value: Any) -> bool:
    text = lower_text(value)
    return "aspergillus" in text or "galactomannan" in text


def raw_gm_source(value: Any) -> str:
    text = lower_text(value)
    if any(term in text for term in ("bal", "balf", "bronchoalveolar", "sputum", "tracheal", "respiratory")):
        return "Respiratory"
    if any(term in text for term in ("blood", "serum", "plasma")):
        return "Systemic"
    return "Unknown"


def aspergillus_gm_label_from_raw(row: dict[str, Any], index: int) -> dict[str, Any] | None:
    target = first_non_empty(row.get("target"), row.get("test"), row.get("assay"), row.get("name"))
    result = first_non_empty(row.get("result"), row.get("interpretation"), row.get("value"))
    if not gm_target_is_aspergillus(target) or not gm_result_is_positive(result):
        return None
    sample = first_non_empty(row.get("sample"), row.get("specimen"), row.get("specimen_type"), row.get("sample_type"))
    collected_time = row.get("collected_time")
    reported_time = row.get("reported_time")
    return {
        "organism_name": "Aspergillus spp.",
        "classification": "Fungal",
        "causative_level": "Level 3",
        "key_evidence": {
            "filmarray_panel": "Unknown",
            "filmarray_sample_type": "Unknown",
            "gm_test_sample_type": sample or "Unknown",
            "gm_test_target": target,
            "gm_test_result": result,
            "non_gm_test_type": "galactomannan_antigen",
            "evidence_type": "direct_detection",
            "deterministic_adjustment": "raw_gm_aspergillus_detected_level3",
            "collected_time": collected_time or "Unknown",
            "reported_time": reported_time or "Unknown",
        },
        "decision_trace": [
            {
                "step": "S2_extract_positive",
                "rule_id": "D-SUM-GM-RAW-ASPERGILLUS-L3",
                "result": "pass",
                "evidence": [
                    {"field": f"gm_test[{index}].target", "value": target},
                    {"field": f"gm_test[{index}].result", "value": result},
                    {"field": f"gm_test[{index}].sample", "value": sample or "Unknown"},
                ],
            }
        ],
    }


def raw_gm_aspergillus_labels(raw_gm: Any) -> tuple[list[dict[str, Any]], str]:
    rows = raw_gm if isinstance(raw_gm, list) else [raw_gm] if isinstance(raw_gm, dict) else []
    labels: list[dict[str, Any]] = []
    source = "Unknown"
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            continue
        label = aspergillus_gm_label_from_raw(row, index)
        if label is None:
            continue
        labels.append(label)
        source = raw_gm_source(label.get("key_evidence", {}).get("gm_test_sample_type")) or source
    return labels, source


def merge_filmarray_gmtest_payloads(
    *,
    primary: dict[str, Any] | None,
    raw_gm: Any,
) -> dict[str, Any] | None:
    if not isinstance(primary, dict) and raw_gm is None:
        return None
    merged: dict[str, Any] = dict(primary) if isinstance(primary, dict) else {}

    labels = list(labels_from_payload(merged))
    seen = {organism_label_key(item) for item in labels if organism_label_key(item)}
    auxiliary_sources: list[str] = []

    raw_labels, raw_source = raw_gm_aspergillus_labels(raw_gm)
    for item in raw_labels:
        key = organism_label_key(item)
        if key and key not in seen:
            labels.append(item)
            seen.add(key)
            auxiliary_sources.append("raw_gm_test")
        elif key:
            labels.append(item)
            auxiliary_sources.append("raw_gm_test")
    if raw_labels:
        merged.setdefault("infection_likelihood", "Possible")
        if lower_text(merged.get("infection_likelihood")) in {"", "none", "unknown", "not_available"}:
            merged["infection_likelihood"] = "Possible"
        if raw_source != "Unknown" and lower_text(merged.get("probable_source")) in {"", "unknown", "not_available"}:
            merged["probable_source"] = raw_source
        if lower_text(merged.get("probable_pathogen_type")) in {"", "unknown", "not_available"}:
            merged["probable_pathogen_type"] = "Fungal"

    if labels:
        merged["organism_labels"] = labels
    if auxiliary_sources:
        merged["merged_auxiliary_sources"] = sorted(set(auxiliary_sources))
    return merged


def module_availability(payloads: dict[str, dict[str, Any] | None]) -> dict[str, str]:
    return {key: "available" if value is not None else "missing" for key, value in payloads.items()}


def data_quality_tier(availability: dict[str, str]) -> str:
    count = sum(1 for value in availability.values() if value == "available")
    if count >= 4:
        return "Q3"
    if count >= 2:
        return "Q2"
    if count == 1:
        return "Q1"
    return "Q0"


def host_context(cbc: dict[str, Any] | None) -> dict[str, Any]:
    state = cbc.get("host_state", {}) if isinstance(cbc, dict) else {}
    if not isinstance(state, dict):
        state = {}
    structured = cbc.get("structured_host_evidence", {}) if isinstance(cbc, dict) else {}
    if not isinstance(structured, dict):
        structured = {}
    medication = structured.get("medication_profile", {})
    if not isinstance(medication, dict):
        medication = {}
    cbc_status = structured.get("cbc_data_status", {})
    if not isinstance(cbc_status, dict):
        cbc_status = {}
    tier = state.get("host_vulnerability_tier", "Unknown")
    opp = state.get("opportunistic_coverage_level", "Unknown")
    expanded = bool(state.get("expanded_candidate_policy", False))
    return {
        "host_vulnerability_tier": tier or "Unknown",
        "expanded_candidate_policy": expanded,
        "opportunistic_coverage_level": opp or "Unknown",
        "key_host_flags": list(state.get("key_host_flags") or []),
        "host_risk_evidence_strength": medication.get("host_risk_evidence_strength", "none"),
        "risk_medications": list(medication.get("medications") or []),
        "host_risk_data_gaps": list(medication.get("data_gaps") or []),
        "cbc_data_status": cbc_status,
        "evidence_preservation_changed_tier": False,
    }


def module_summaries(
    *,
    cbc: dict[str, Any] | None,
    filmarray: dict[str, Any] | None,
    image: dict[str, Any] | None,
    culture: dict[str, Any] | None,
    molecular: dict[str, Any] | None,
) -> dict[str, Any]:
    ctx = host_context(cbc)
    assay = filmarray.get("assay_summary", {}) if isinstance(filmarray, dict) else {}
    if not isinstance(assay, dict):
        assay = {}
    return {
        "cbc_other_lab": {
            "infection_likelihood": (cbc or {}).get("infection_likelihood", "Unknown"),
            "probable_source": normalize_source((cbc or {}).get("probable_source", "Unknown")),
            "probable_pathogen_type": (cbc or {}).get("probable_pathogen_type", "Unknown"),
            "support_direction": (cbc or {}).get("support_direction", "Not_available"),
            **ctx,
        },
        "filmarray_gmtest": {
            "infection_likelihood": (filmarray or {}).get("infection_likelihood", "Unknown"),
            "probable_source": normalize_source((filmarray or {}).get("probable_source", "Unknown")),
            "probable_pathogen_type": (filmarray or {}).get("probable_pathogen_type", "Unknown"),
            "panel": assay.get("filmarray_panel", "Unknown"),
        },
        "image": {
            "infection_likelihood": (image or {}).get("infection_likelihood", "Unknown"),
            "probable_source": normalize_source((image or {}).get("probable_source", "Unknown")),
            "probable_pathogen_type": (image or {}).get("probable_pathogen_type", "Unknown"),
            "support_direction": (image or {}).get("support_direction", "Not_available"),
            "diagnostic_mentions": list((image or {}).get("diagnostic_mentions") or []),
        },
        "culture": {
            "infection_likelihood": (culture or {}).get("infection_likelihood", "Unknown"),
            "probable_source": normalize_source((culture or {}).get("probable_source", "Unknown")),
            "probable_pathogen_type": (culture or {}).get("probable_pathogen_type", "Unknown"),
        },
        "molecular_microbiology": {
            "infection_likelihood": (molecular or {}).get("infection_likelihood", "Unknown"),
            "probable_source": normalize_source((molecular or {}).get("probable_source", "Unknown")),
            "probable_pathogen_type": (molecular or {}).get("probable_pathogen_type", "Unknown"),
        },
    }


def add_hospital_evidence(
    evidence_by_key: dict[str, dict[str, Any]],
    *,
    module: str,
    organism_name: Any,
    classification: Any,
    level: Any,
    source_index: int,
    key_evidence: Any = None,
    source_item: Any = None,
    module_payload: Any = None,
) -> None:
    name = str(organism_name or "").strip()
    if not name:
        return
    normalized_class = classification_from_source(classification)
    key = f"{canonical_name(name)}|{normalized_class}"
    if not key.strip("|"):
        return
    adjusted_level, adjusted_key_evidence, was_adjusted = adjust_culture_detected_without_quantity_level(
        module=module,
        organism_name=name,
        classification=normalized_class,
        level=level,
        key_evidence=key_evidence,
    )
    adjustment_rule_id = "D-SUM-CULTURE-QTY-UNKNOWN-L2" if was_adjusted else ""
    adjustment_message = (
        "typical lower-respiratory pneumonia bacterium detected by culture with unknown quantity; "
        "treated as Level 2 support, not low quantity"
        if was_adjusted
        else ""
    )
    if not was_adjusted:
        adjusted_level, adjusted_key_evidence, was_adjusted = adjust_filmarray_detected_without_quantity_level(
            module=module,
            organism_name=name,
            classification=normalized_class,
            level=level,
            key_evidence=key_evidence,
            source_item=source_item,
            module_payload=module_payload,
        )
        if was_adjusted:
            adjustment_rule_id = "D-SUM-FILMARRAY-QTY-UNKNOWN-L2"
            adjustment_message = (
                "typical lower-respiratory pneumonia bacterium detected by FilmArray without semiquant value; "
                "treated as Level 2 support, not low quantity"
            )
    evidence_item = evidence_by_key.setdefault(
        key,
        {
            "organism_name": name,
            "classification": normalized_class,
            "evidence_modules": [],
            "module_level_summary": {
                "culture": "Not_available",
                "filmarray_gmtest": "Not_available",
                "image": "Not_available",
                "molecular_microbiology": "Not_available",
                "cbc_other_lab": "Not_pathogen_specific",
            },
            "module_evidence": {},
            "evidence_trace": [],
        },
    )
    if module not in evidence_item["evidence_modules"]:
        evidence_item["evidence_modules"].append(module)
    evidence_item["module_level_summary"][module] = adjusted_level
    if isinstance(adjusted_key_evidence, dict):
        module_rows = evidence_item["module_evidence"].setdefault(module, [])
        if isinstance(module_rows, list):
            module_rows.append(adjusted_key_evidence)
    evidence_item["best_hospital_level"] = best_level(evidence_item["module_level_summary"].values())
    trace_evidence = [
        {"field": f"{module}.organism_labels[{source_index}].organism_name", "value": name},
        {"field": f"{module}.organism_labels[{source_index}].causative_level", "value": str(level or "Unknown")},
        {"field": f"{module}.organism_labels[{source_index}].classification", "value": normalized_class},
    ]
    rule_id = "D-SUM-01"
    if was_adjusted:
        rule_id = adjustment_rule_id
        trace_evidence.append(
            {
                "field": "deterministic_adjustment",
                "value": adjustment_message,
            }
        )
    evidence_item["evidence_trace"].append(
        {
            "step": "S1_extract_modules",
            "rule_id": rule_id,
            "result": "pass",
            "evidence": trace_evidence,
        }
    )


def collect_hospital_organism_evidence(
    *,
    culture: dict[str, Any] | None,
    filmarray: dict[str, Any] | None,
    image: dict[str, Any] | None,
    molecular: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    evidence_by_key: dict[str, dict[str, Any]] = {}
    for module, payload in (
        ("culture", culture),
        ("filmarray_gmtest", filmarray),
        ("molecular_microbiology", molecular),
    ):
        labels = payload.get("organism_labels", []) if isinstance(payload, dict) else []
        if not isinstance(labels, list):
            continue
        for index, item in enumerate(labels):
            if not isinstance(item, dict):
                continue
            add_hospital_evidence(
                evidence_by_key,
                module=module,
                organism_name=item.get("organism_name"),
                classification=item.get("classification"),
                level=item.get("causative_level", item.get("level", "Unknown")),
                source_index=index,
                key_evidence=item.get("key_evidence"),
                source_item=item,
                module_payload=payload,
            )
    hints = image.get("organism_hints", []) if isinstance(image, dict) else []
    if isinstance(hints, list):
        for index, item in enumerate(hints):
            if not isinstance(item, dict):
                continue
            add_hospital_evidence(
                evidence_by_key,
                module="image",
                organism_name=item.get("organism_name"),
                classification=item.get("classification", item.get("pathogen_type")),
                level=item.get("causative_level", item.get("level", "Level 3")),
                source_index=index,
            )
    output = annotate_approved_group_member_relationships(list(evidence_by_key.values()))
    output.sort(
        key=lambda item: (
            level_rank(item.get("best_hospital_level")),
            str(item.get("organism_name", "")).lower(),
        )
    )
    return output


def annotate_approved_group_member_relationships(
    evidence_items: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Link approved group/complex labels to detected member species.

    A panel group label is not an exact-species alias. Both observations remain
    visible, while the relationship is made explicit for downstream scoring and
    display. This preserves the assay's taxonomic scope and avoids duplicate
    unlinked candidates.
    """
    for group_item in evidence_items:
        group_name = group_item.get("organism_name")
        approved_members = pathogen_names.approved_group_members(group_name)
        if not approved_members:
            continue
        group_item["identity_scope"] = "approved_group_or_complex"
        group_item["group_level_support_only"] = True
        for member_item in evidence_items:
            if member_item is group_item:
                continue
            member_name = member_item.get("organism_name")
            member_key = pathogen_names.canonical_key(member_name)
            if member_key not in approved_members:
                continue
            if not pathogen_names.approved_group_member_match(group_name, member_name):
                continue

            group_links = group_item.setdefault("approved_detected_members", [])
            member_link = {
                "organism_name": str(member_name or ""),
                "canonical_key": member_key,
                "relationship": "approved_group_member_match",
            }
            if member_link not in group_links:
                group_links.append(member_link)

            related_groups = member_item.setdefault("approved_group_support", [])
            group_link = {
                "organism_name": str(group_name or ""),
                "canonical_key": pathogen_names.canonical_key(group_name),
                "relationship": "approved_group_member_match",
                "not_exact_species_support": True,
                "evidence_modules": list(group_item.get("evidence_modules") or []),
                "module_level_summary": dict(group_item.get("module_level_summary") or {}),
            }
            if group_link not in related_groups:
                related_groups.append(group_link)

            trace = member_item.setdefault("evidence_trace", [])
            trace.append(
                {
                    "step": "S1_identity_link",
                    "rule_id": "D-SUM-ID-APPROVED-GROUP-MEMBER",
                    "result": "context_only",
                    "evidence": [
                        {"field": "group_label", "value": str(group_name or "")},
                        {"field": "member_species", "value": str(member_name or "")},
                        {"field": "not_exact_species_support", "value": True},
                    ],
                }
            )
    return evidence_items


def final_infection_likelihood(evidence_items: Sequence[dict[str, Any]], modules: dict[str, Any]) -> str:
    levels = [level_rank(item.get("best_hospital_level")) for item in evidence_items]
    if any(level == 1 for level in levels):
        return "Severe"
    if any(level == 2 for level in levels):
        return "Likely"
    if any(level == 3 for level in levels):
        return "Possible"
    likelihoods = [
        module.get("infection_likelihood", "Unknown")
        for module in modules.values()
        if isinstance(module, dict)
    ]
    best = min(likelihoods, key=likelihood_rank) if likelihoods else "Unknown"
    # Keep Possible as patient-level context, but do not use it as pathogen support.
    return best if best not in {"None", "Not_available"} else "Unknown"


def dominant_pathogen_type(evidence_items: Sequence[dict[str, Any]]) -> str:
    classes = {
        str(item.get("classification", "Unknown"))
        for item in evidence_items
        if level_rank(item.get("best_hospital_level")) <= 3
        and str(item.get("classification", "Unknown")) != "Unknown"
    }
    if len(classes) == 1:
        return next(iter(classes))
    if len(classes) > 1:
        return "Mixed"
    return "Unknown"


def dominant_source(evidence_items: Sequence[dict[str, Any]], modules: dict[str, Any]) -> str:
    source_by_module = {
        "culture": modules.get("culture", {}).get("probable_source", "Unknown"),
        "filmarray_gmtest": modules.get("filmarray_gmtest", {}).get("probable_source", "Unknown"),
        "image": modules.get("image", {}).get("probable_source", "Unknown"),
        "molecular_microbiology": modules.get("molecular_microbiology", {}).get(
            "probable_source", "Unknown"
        ),
    }
    for evidence_item in evidence_items:
        if level_rank(evidence_item.get("best_hospital_level")) > 3:
            continue
        for module in evidence_item.get("evidence_modules", []):
            source = source_by_module.get(module, "Unknown")
            if source not in {"Unknown", "Not_available", "", None}:
                return normalize_source(source)
    for module in ("culture", "filmarray_gmtest", "molecular_microbiology", "image"):
        source = source_by_module.get(module, "Unknown")
        if source not in {"Unknown", "Not_available", "", None}:
            return normalize_source(source)
    return "Unknown"


def support_direction(evidence_items: Sequence[dict[str, Any]], image: dict[str, Any] | None) -> str:
    if any(
        level_rank(item.get("best_hospital_level")) <= 3
        and any(
            module in {"culture", "filmarray_gmtest", "molecular_microbiology"}
            for module in item.get("evidence_modules", [])
        )
        for item in evidence_items
    ):
        return "Support"
    if isinstance(image, dict) and image.get("support_direction") == "Support":
        return "Support"
    return "Not_available"


def amr_risk_summary(filmarray: dict[str, Any] | None) -> dict[str, Any]:
    if not isinstance(filmarray, dict):
        return {"risk_level": "Unknown", "evidence": ["not_available"]}
    findings = filmarray.get("resistance_findings")
    if not isinstance(findings, list) or not findings:
        findings = ["not_available"]
    return {
        "risk_level": filmarray.get("resistance_risk", "Unknown"),
        "evidence": findings,
    }


def merge_data_gaps(payloads: dict[str, dict[str, Any] | None]) -> list[str]:
    gaps: list[str] = []
    for module, payload in payloads.items():
        if payload is None:
            gaps.append(f"{module}: missing_agent_output")
            continue
        for item in payload.get("data_gaps", []) or []:
            text = str(item)
            if text and text not in gaps:
                gaps.append(f"{module}: {text}")
    return gaps


def build_summary(patient_dir: Path) -> dict[str, Any]:
    paths = {
        "cbc_other_lab": agent_path(patient_dir, "cbc_other_lab_agent"),
        "filmarray_gmtest": agent_path(patient_dir, "filmarray_gmtest_agent"),
        "raw_gm_test": patient_source_path(patient_dir, "gm_test"),
        "raw_culture": patient_source_path(patient_dir, "culture"),
        "raw_filmarray": patient_source_path(patient_dir, "filmarray"),
        "image": agent_path(patient_dir, "image_agent"),
        "culture": agent_path(patient_dir, "culture_agent"),
        "molecular_microbiology": agent_path(patient_dir, "molecular_microbiology_agent"),
        "raw_image": patient_source_path(patient_dir, "image"),
        "raw_cbc": patient_source_path(patient_dir, "CBC"),
        "raw_other_lab": patient_source_path(patient_dir, "other_lab"),
        "raw_underlying": patient_source_path(patient_dir, "underlying"),
        "raw_admission_diagnosis": patient_source_path(patient_dir, "admission_diagnosis"),
        "raw_mngs_ranked": patient_source_path(patient_dir, "mNGS_ranked_candidates"),
    }
    raw_path_keys = {
        "raw_gm_test",
        "raw_culture",
        "raw_filmarray",
        "raw_image",
        "raw_cbc",
        "raw_other_lab",
        "raw_underlying",
        "raw_admission_diagnosis",
        "raw_mngs_ranked",
    }
    payloads: dict[str, dict[str, Any] | None] = {
        key: load_json(path) if path.exists() else None
        for key, path in paths.items()
        if key not in raw_path_keys
    }
    agent_availability = module_availability(payloads)
    raw_sources = {
        key: load_json_any(paths[key]) if paths[key].exists() else None
        for key in raw_path_keys
    }
    deterministic_fallbacks: dict[str, str] = {}
    if payloads.get("culture") is None and raw_sources.get("raw_culture") is not None:
        payloads["culture"] = normalized_agent_fallback.culture_agent_from_source(
            raw_sources["raw_culture"]
        )
        deterministic_fallbacks["culture"] = "normalized_culture_source"
    if payloads.get("filmarray_gmtest") is None and raw_sources.get("raw_filmarray") is not None:
        payloads["filmarray_gmtest"] = normalized_agent_fallback.filmarray_agent_from_source(
            raw_sources["raw_filmarray"]
        )
        deterministic_fallbacks["filmarray_gmtest"] = "normalized_filmarray_source"
    raw_gm_payload = raw_sources.get("raw_gm_test")
    payloads["filmarray_gmtest"] = merge_filmarray_gmtest_payloads(
        primary=payloads.get("filmarray_gmtest"),
        raw_gm=raw_gm_payload,
    )
    resolved_availability = module_availability(payloads)
    sample_time = evidence_preservation.mngs_sample_time_from_ranked(raw_sources["raw_mngs_ranked"])
    payloads["image"] = evidence_preservation.enrich_image_agent_output(
        payloads.get("image"), raw_sources["raw_image"]
    )
    payloads["cbc_other_lab"] = evidence_preservation.enrich_host_agent_output(
        payloads.get("cbc_other_lab"),
        raw_cbc=raw_sources["raw_cbc"],
        raw_underlying=raw_sources["raw_underlying"],
        raw_admission=raw_sources["raw_admission_diagnosis"],
        sample_time=sample_time,
    )
    modules = module_summaries(
        cbc=payloads["cbc_other_lab"],
        filmarray=payloads["filmarray_gmtest"],
        image=payloads["image"],
        culture=payloads["culture"],
        molecular=payloads["molecular_microbiology"],
    )
    hospital_evidence = collect_hospital_organism_evidence(
        culture=payloads["culture"],
        filmarray=payloads["filmarray_gmtest"],
        image=payloads["image"],
        molecular=payloads["molecular_microbiology"],
    )
    availability = resolved_availability
    best_likelihood = final_infection_likelihood(hospital_evidence, modules)
    source = dominant_source(hospital_evidence, modules)
    pathogen_type = dominant_pathogen_type(hospital_evidence)
    image_payload = payloads["image"]
    return {
        "rule_version": RULE_VERSION,
        "final_infection_likelihood": best_likelihood,
        "dominant_source": source,
        "dominant_pathogen_type": pathogen_type,
        "support_direction": support_direction(hospital_evidence, image_payload),
        "data_quality_tier": data_quality_tier(availability),
        "module_availability": availability,
        "amr_risk_summary": amr_risk_summary(payloads["filmarray_gmtest"]),
        "host_context": host_context(payloads["cbc_other_lab"]),
        "evidence_preservation": {
            "version": evidence_preservation.PRESERVATION_VERSION,
            "agent_module_availability": agent_availability,
            "original_agent_module_availability": agent_availability,
            "resolved_module_availability": resolved_availability,
            "deterministic_source_fallbacks": deterministic_fallbacks,
            "raw_source_availability": {
                key: "available" if value is not None else "missing"
                for key, value in raw_sources.items()
            },
            "imaging_diagnostic_mentions": list(
                (payloads.get("image") or {}).get("diagnostic_mentions") or []
            ),
            "structured_host_evidence": dict(
                (payloads.get("cbc_other_lab") or {}).get("structured_host_evidence") or {}
            ),
            "source_evidence_index": evidence_preservation.build_source_evidence_index(
                {
                    key.removeprefix("raw_"): value
                    for key, value in raw_sources.items()
                }
            ),
            "tier_changed_by_preservation": False,
        },
        "time_span": (image_payload or {}).get("time_span", {"earliest": "Unknown", "latest": "Unknown"}),
        "hospital_organism_evidence": hospital_evidence,
        "pathogen_candidates": [],
        "excluded_candidates": [],
        "priority_flags": [],
        "cross_module_reasoning": [
            "D-SUM-00: deterministic summary was built directly from small agent JSON outputs.",
            "D-SUM-01: hospital_organism_evidence is pooled from culture, FilmArray/GM, molecular microbiology, and image organism_hints; this step does not decide pathogen_candidates or excluded_candidates.",
            "D-SUM-01B: FilmArray/GM evidence uses the with-FilmArray agent output and also merges raw gm_test Aspergillus galactomannan detected results when needed.",
            "D-SUM-02: negative or missing module results are recorded as Not_available/Unknown and are not used as counter-evidence.",
            f"D-SUM-03: final_infection_likelihood={best_likelihood}; dominant_source={source}; dominant_pathogen_type={pathogen_type}.",
        ],
        "module_summaries": modules,
        "data_gaps": merge_data_gaps(payloads),
        "source_files": {key: str(path) for key, path in paths.items()},
    }


def process_patient(patient_dir: Path, *, output_suffix: str, overwrite: bool, skip_existing: bool) -> Path | None:
    destination = output_path(patient_dir, output_suffix)
    if destination.exists() and skip_existing:
        return None
    if destination.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {destination}. Use --overwrite or --skip-existing.")
    summary = build_summary(patient_dir)
    write_json(destination, summary)
    return destination


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    written = []
    for patient_dir in collect_patient_dirs(args.inputs, args.patients):
        output = process_patient(
            patient_dir,
            output_suffix=args.output_suffix,
            overwrite=args.overwrite,
            skip_existing=args.skip_existing,
        )
        if output is not None:
            written.append(output)
            print(f"wrote {output}")
        else:
            print(f"skip {patient_dir.name}: output exists")
    print(f"written_count={len(written)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

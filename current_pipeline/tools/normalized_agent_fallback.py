"""Deterministic agent-shaped fallbacks for normalized microbiology sources.

These helpers are used only when a small-agent output is absent. They implement
the fixed culture, FilmArray/GM, and targeted molecular rules already documented
in the active prompts; they do not use benchmark answers or call an LLM.
"""

from __future__ import annotations

import copy
import hashlib
import json
import re
from typing import Any

from tools import organism_taxonomy_classifier as taxonomy
from tools import pathogen_normalization as pathogen_names


LEVEL_RANK = {f"Level {index}": index for index in range(1, 6)}
SPECIMEN_PRIORITY = {
    "Sterile_Site": 0,
    "Lower_Respiratory": 1,
    "NonSterile_Other": 2,
    "Unknown": 3,
}
QUANTITY_PRIORITY = {"Q4": 0, "Q3": 1, "Q2": 2, "Q1": 3, "Q0": 4, "Unknown": 5}
PNEUMONIA_AMR_TARGETS = {
    "ctxm",
    "kpc",
    "ndm",
    "oxa48like",
    "vim",
    "imp",
    "mecacandmrej",
}
BCID_AMR_TARGETS = {
    "vana",
    "vanb",
    "mcr1",
}
AMR_TARGETS = PNEUMONIA_AMR_TARGETS | BCID_AMR_TARGETS
CONTAMINANT_KEYS = {
    "coagulasenegativestaphylococcus",
    "staphylococcusepidermidis",
    "staphylococcushemolyticus",
    "staphylococcushaemolyticus",
    "cutibacteriumacnes",
    "micrococcusspp",
}
TYPICAL_RESPIRATORY_GROUP_KEYS = {
    # Curated group labels are deliberately kept distinct from exact species,
    # but their respiratory culture evidence should inherit the appropriate
    # clinical family instead of falling through to environmental/unknown.
    "acinetobactercalcoaceticusbaumanniicomplex",
}
HARD_TO_CULTURE_TERMS = (
    "legionella",
    "nocardia",
    "mycobacterium",
    "histoplasma",
    "cryptococcus",
    "aspergillus",
    "mucor",
    "rhizopus",
    "cunninghamella",
    "pneumocystis",
    "brucella",
    "bartonella",
    "francisella",
)


def _rows(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    return [value] if isinstance(value, dict) else []


def _text(value: Any) -> str:
    return str(value or "").strip()


def _positive_text(value: Any) -> bool:
    text = _text(value).lower()
    if not text or any(term in text for term in ("not detected", "negative", "no growth", "not isolated")):
        return False
    return any(term in text for term in ("detected", "positive", "reactive", "isolated", "growth"))


def _classification(name: Any) -> str:
    biological = str(taxonomy.classify_organism(name).get("biological_class") or "unknown")
    return {
        "bacterium": "Bacterial",
        "fungus": "Fungal",
        "virus": "Viral",
        "parasite": "Parasitic",
    }.get(biological, "Unknown")


def _specimen_category(value: Any) -> str:
    text = _text(value).lower()
    if any(term in text for term in ("bal", "bronchoalveolar", "sputum", "tracheal", "eta", "lower respiratory")):
        return "Lower_Respiratory"
    if any(term in text for term in ("blood", "csf", "pleural", "peritoneal", "synovial", "joint", "tissue", "bone", "sterile fluid")):
        return "Sterile_Site"
    if any(term in text for term in ("urine", "wound", "skin", "throat", "stool", "rectal", "genital")):
        return "NonSterile_Other"
    return "Unknown"


def _quantity_tier(value: Any) -> tuple[str, str]:
    text = _text(value).lower().replace("≥", ">=").replace("≤", "<=")
    if not text:
        return "Unknown", "positive_detected_quantity_unknown"
    if any(term in text for term in ("heavy", "many", "4+", "predominant")):
        return "Q4", "reported"
    if any(term in text for term in ("moderate", "3+")):
        return "Q3", "reported"
    if any(term in text for term in ("few", "2+")):
        return "Q2", "reported"
    if any(term in text for term in ("rare", "1+")):
        return "Q1", "reported"
    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*[x×]\s*10\s*\^\s*([0-9]+)", text)
    if match:
        quantity = float(match.group(1)) * (10 ** int(match.group(2)))
    else:
        match = re.search(r"10\s*\^\s*([0-9]+)", text)
        quantity = 10 ** int(match.group(1)) if match else None
    if quantity is None:
        direct = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(?:cfu|copy)", text)
        quantity = float(direct.group(1)) if direct else None
    if quantity is None:
        return "Unknown", "positive_detected_quantity_unknown"
    if quantity >= 100000:
        return "Q4", "reported"
    if quantity >= 10000:
        return "Q3", "reported"
    if quantity >= 1000:
        return "Q2", "reported"
    if quantity > 0:
        return "Q1", "reported"
    return "Q0", "reported"


def _growth_purity(row: dict[str, Any]) -> str:
    text = " ".join(_text(row.get(key)).lower() for key in ("status", "growth", "comment", "organism"))
    if any(term in text for term in ("no growth", "negative", "not isolated")):
        return "no_growth"
    if any(term in text for term in ("contamin", "mixed flora", "normal flora")):
        return "contaminated"
    if any(term in text for term in ("isolated", "pure growth", "single organism", "predominant")):
        return "isolated_pure"
    return "mixed_unspecified"


def _is_typical_bacterium(name: Any) -> bool:
    if pathogen_names.canonical_key(name) in TYPICAL_RESPIRATORY_GROUP_KEYS:
        return True
    profile = taxonomy.classify_organism(name, biological_class="Bacterial")
    return profile.get("primary_rule_family") in {
        "typical_respiratory_pathogen",
        "hospital_or_nonfermenter_gnb",
    }


def _is_exact_candida_species(name: Any) -> bool:
    """Return true for a named Candida species, not generic yeast/Candida labels."""
    profile = taxonomy.classify_organism(name, biological_class="Fungal")
    return (
        pathogen_names.genus_name(name) == "candida"
        and profile.get("taxonomic_rank") in {"species", "species_candidate"}
    )


def _is_exact_colonizer_prone_species(name: Any) -> bool:
    """Identify named colonizer-prone species that can still be true pathogens."""
    profile = taxonomy.classify_organism(name, biological_class="Bacterial")
    return (
        profile.get("primary_rule_family") == "skin_airway_colonizer_prone"
        and profile.get("taxonomic_rank") in {"species", "species_candidate"}
    )


def _is_contaminant_prone(name: Any) -> bool:
    key = pathogen_names.canonical_key(name)
    text = _text(name).lower()
    return (
        key in CONTAMINANT_KEYS
        or text in {"corynebacterium spp.", "bacillus spp.", "bacillus spp. non-anthracis"}
    )


def _culture_level(
    name: str,
    specimen: str,
    quantity: str,
    purity: str,
) -> tuple[str, list[str], str, str]:
    rules = ["R-S0-02", "R-S0-03", "R-S0-04", "R-S0-05"]
    contaminant = "Unknown"
    hard = "Yes" if any(term in name.lower() for term in HARD_TO_CULTURE_TERMS) else "No"
    sterile_strong = specimen == "Sterile_Site" and purity == "isolated_pure" and quantity != "Unknown"
    if purity in {"no_growth", "contaminated"}:
        return "Level 5", [*rules, "R-S1-01"], "Yes", hard
    if _is_contaminant_prone(name) and not sterile_strong:
        return "Level 4", [*rules, "R-S1-02"], "Yes", hard
    if pathogen_names.is_candida_or_generic_yeast(name) and specimen == "Lower_Respiratory":
        level = "Level 3" if quantity in {"Q4", "Q3"} and purity == "isolated_pure" else "Level 4"
        return level, [*rules, "R-S1-03"], "No", hard
    if specimen == "Sterile_Site" and purity == "isolated_pure" and _is_exact_candida_species(name):
        # Blood/sterile-site culture is qualitative by design: exact Candida
        # does not require a reported CFU value to be strong direct evidence.
        return "Level 1", [*rules, "R-S2-01", "R-S3-01", "R-S3-01A"], "No", hard
    if specimen == "Sterile_Site" and purity == "isolated_pure" and _is_exact_colonizer_prone_species(name):
        # An exact named species in blood/sterile-site is more informative than
        # a generic commensal label, but independence and contamination cannot
        # be resolved from one normalized record. Keep provisional support.
        return "Level 3", [*rules, "R-S1-02A", "R-S2-01"], "Unknown", hard
    if specimen == "Sterile_Site":
        level = "Level 1" if sterile_strong else "Level 2"
        rules.extend(["R-S2-01", "R-S3-01"])
    elif specimen == "Lower_Respiratory":
        if quantity == "Q4" and purity == "isolated_pure":
            level = "Level 2"
        elif quantity in {"Q3", "Q2"}:
            level = "Level 3"
        elif quantity == "Unknown" and _is_typical_bacterium(name):
            level = "Level 2"
        else:
            level = "Level 4"
        rules.append("R-S3-02")
    elif specimen == "NonSterile_Other":
        level = "Level 3" if quantity == "Q4" and purity == "isolated_pure" else "Level 4"
        rules.append("R-S3-03")
    else:
        level = "Level 3"
        rules.extend(["R-S2-02", "R-S3-04"])
    if hard == "Yes" and not pathogen_names.is_candida_or_generic_yeast(name):
        current = LEVEL_RANK[level]
        floor = 1 if specimen == "Sterile_Site" else 2 if specimen == "Lower_Respiratory" else 3
        level = f"Level {max(floor, current - 1)}"
        rules.append("R-S4-02")
    return level, rules, contaminant, hard


def _source_record_key(index: int, row: dict[str, Any]) -> str:
    encoded = json.dumps(
        row, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    digest = hashlib.sha256(encoded).hexdigest()[:16]
    return f"culture[{index}]#{digest}"


def _assay_source_record_key(
    source_name: str, index: int, row: dict[str, Any]
) -> str:
    """Stable trace key for one FilmArray/GM source row."""
    encoded = json.dumps(
        row, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    digest = hashlib.sha256(encoded).hexdigest()[:16]
    return f"{source_name}[{index}]#{digest}"


def _assay_detection_status(result: Any) -> str:
    text = _text(result).casefold()
    if any(term in text for term in ("not detected", "negative", "non-reactive", "nonreactive")):
        return "not_detected"
    if "equivocal" in text or "borderline" in text:
        return "equivocal"
    if any(term in text for term in ("invalid", "failed", "inconclusive")):
        return "invalid"
    if _positive_text(result):
        return "detected"
    return "unknown"


def _numeric_value_and_unit(value: Any) -> tuple[str, str]:
    raw = _text(value)
    if not raw:
        return "", ""
    match = re.search(r"[-+]?\d+(?:\.\d+)?", raw.replace(",", ""))
    numeric = match.group(0) if match else ""
    unit = ""
    lowered = raw.casefold()
    for candidate in ("copy/ml", "copies/ml", "gc/ml", "index", "ntu", "iu/ml"):
        if candidate in lowered:
            unit = candidate
            break
    return numeric, unit


def _float_or_none(value: Any) -> float | None:
    numeric, _ = _numeric_value_and_unit(value)
    try:
        return float(numeric) if numeric else None
    except ValueError:
        return None


def _normalized_non_gm_target(value: Any) -> str:
    text = _text(value)
    text = re.sub(
        r"\b(?:Ig[AMG]|antibody|antigen|serology|titer|Ab|Ag)\b.*$",
        "",
        text,
        flags=re.IGNORECASE,
    )
    return re.sub(r"\s+", " ", text).strip(" -_:/")


def _non_gm_test_type(value: Any) -> tuple[str, str]:
    text = _text(value).casefold()
    if "igm" in text:
        return "serology_IgM", "indirect_serology"
    if "igg" in text:
        return "serology_IgG", "indirect_serology"
    if any(term in text for term in ("antibody", "serology", " titer", " ab")):
        return "serology_antibody", "indirect_serology"
    if any(term in text for term in ("antigen", " ag")):
        return "antigen", "direct_detection"
    return "other_non_gm", "unknown"


def _quantity_parts(value: Any) -> tuple[str, str]:
    raw = _text(value)
    if not raw:
        return "", ""
    unit = "CFU/mL" if "cfu" in raw.casefold() else ""
    quantity = re.sub(r"\s*CFU\s*/\s*mL\s*$", "", raw, flags=re.IGNORECASE)
    return quantity.strip(), unit


def culture_agent_from_source(source: Any) -> dict[str, Any]:
    labels: list[dict[str, Any]] = []
    observations: list[dict[str, Any]] = []
    data_gaps: set[str] = set()
    rows = _rows(source)
    for index, row in enumerate(rows):
        name = _text(row.get("organism"))
        status = _text(row.get("status"))
        specimen_type = _text(row.get("sample") or row.get("specimen"))
        specimen = _specimen_category(specimen_type)
        raw_quantity = row.get("colony_count") or row.get("quantity") or row.get("value")
        status_lower = status.casefold()
        explicitly_negative = any(
            term in status_lower
            for term in ("not detected", "negative", "no growth", "not isolated")
        )
        # A named organism on a preliminary culture report is still a
        # positive/provisional observation.  Keep it visible and let the
        # downstream culture rules express the lower certainty.
        positive = bool(name) and not explicitly_negative
        if positive:
            quantity, quantitation = _quantity_tier(raw_quantity)
        else:
            quantity, quantitation = "Unknown", "not_applicable"
        quantity_value, quantity_unit = _quantity_parts(raw_quantity)
        purity = _growth_purity(row)
        observation_id = f"CUL-{index + 1:03d}"
        observations.append({
            "observation_id": observation_id,
            "source_record_index": index,
            "source_record_key": _source_record_key(index, row),
            "organism_name": name,
            "test": _text(row.get("test")),
            "specimen_type": specimen_type,
            "specimen_category": specimen,
            "collected_time": _text(row.get("collected_time")),
            "received_time": _text(row.get("received_time")),
            "reported_time": _text(row.get("reported_time")),
            "raw_result": status,
            "raw_quantity": _text(raw_quantity),
            "quantity_value": quantity_value,
            "quantity_unit": quantity_unit,
            "quantity_tier": quantity,
            "quantitation_status": quantitation,
            "growth_purity": purity,
            "source_text": json.dumps(
                row, ensure_ascii=False, sort_keys=True, separators=(",", ":")
            ),
        })
        if not positive:
            continue
        if specimen == "Unknown":
            data_gaps.add("missing_specimen_type")
        if quantity == "Unknown":
            data_gaps.add("missing_quantity")
        level, applied, contaminant, hard = _culture_level(name, specimen, quantity, purity)
        labels.append(
            {
                "organism_name": name,
                "classification": _classification(name),
                "causative_level": level,
                "applied_rules": applied,
                "evidence_observation_ids": [observation_id],
                "key_evidence": {
                    "specimen_type": specimen_type,
                    "specimen_category": specimen,
                    "quantity_tier": quantity,
                    "quantitation_status": quantitation,
                    "growth_purity": purity,
                    "contaminant_flag": contaminant,
                    "hard_to_culture": hard,
                    "collected_time": row.get("collected_time") or "Unknown",
                    "reported_time": row.get("reported_time") or "Unknown",
                    "source_rule": "normalized_source_deterministic_fallback",
                },
            }
        )
    labels = _dedupe_labels(labels)
    highest = min((item["causative_level"] for item in labels), key=lambda value: LEVEL_RANK[value], default="Not_available")
    likelihood = _likelihood(highest)
    top = [item for item in labels if item["causative_level"] == highest]
    source_name = _culture_source(top[0]["key_evidence"]["specimen_type"]) if top else "Unknown"
    classes = {item["classification"] for item in top if item["classification"] != "Unknown"}
    pathogen_type = next(iter(classes)) if len(classes) == 1 else "Mixed" if classes else "Unknown"
    return {
        "rule_version": "culture_normalized_source_fallback_v1",
        "infection_likelihood": likelihood,
        "probable_source": source_name,
        "probable_pathogen_type": pathogen_type,
        "culture_observations": observations,
        "organism_labels": labels,
        "module_rule_summary": {
            "highest_level": highest,
            "source_rule_id": "R-S5-03",
            "pathogen_type_rule_id": "R-S5-04",
        },
        "data_gaps": sorted(data_gaps),
        "source_record_contract": {
            "source_record_count": len(rows),
            "observation_count": len(observations),
            "one_observation_per_source_record": len(rows) == len(observations),
        },
    }


def _dedupe_labels(labels: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: dict[tuple[str, str], dict[str, Any]] = {}
    for original in labels:
        item = copy.deepcopy(original)
        key = (pathogen_names.canonical_key(item.get("organism_name")), str(item.get("classification")))
        current = selected.get(key)
        item_order = (
            LEVEL_RANK.get(str(item.get("causative_level")), 99),
            SPECIMEN_PRIORITY.get(str((item.get("key_evidence") or {}).get("specimen_category")), 99),
            QUANTITY_PRIORITY.get(str((item.get("key_evidence") or {}).get("quantity_tier")), 99),
        )
        if current is None:
            selected[key] = item
            continue
        current_order = (
            LEVEL_RANK.get(str(current.get("causative_level")), 99),
            SPECIMEN_PRIORITY.get(str((current.get("key_evidence") or {}).get("specimen_category")), 99),
            QUANTITY_PRIORITY.get(str((current.get("key_evidence") or {}).get("quantity_tier")), 99),
        )
        merged_ids = list(dict.fromkeys([
            *(current.get("evidence_observation_ids") or []),
            *(item.get("evidence_observation_ids") or []),
        ]))
        if item_order < current_order:
            if merged_ids:
                item["evidence_observation_ids"] = merged_ids
            selected[key] = item
        elif merged_ids:
            current["evidence_observation_ids"] = merged_ids
    return sorted(
        selected.values(),
        key=lambda item: (LEVEL_RANK.get(str(item.get("causative_level")), 99), str(item.get("organism_name", "")).lower()),
    )


def _likelihood(level: str) -> str:
    return {
        "Level 1": "Severe",
        "Level 2": "Likely",
        "Level 3": "Possible",
    }.get(level, "None" if level in {"Level 4", "Level 5", "Not_available"} else "Unknown")


def _culture_source(value: Any) -> str:
    category = _specimen_category(value)
    text = _text(value).lower()
    if "blood" in text:
        return "Bloodstream"
    if category == "Lower_Respiratory":
        return "Respiratory"
    if "urine" in text:
        return "Urinary"
    if any(term in text for term in ("wound", "skin")):
        return "Wound"
    if category == "Sterile_Site":
        return "Systemic"
    return "Other" if category == "NonSterile_Other" else "Unknown"


def _molecular_specimen_category(value: Any) -> str:
    """Map a targeted-molecular specimen without collapsing CNS into sterile site."""
    text = _text(value).casefold()
    if any(term in text for term in ("csf", "cerebrospinal")):
        return "CNS"
    if any(
        term in text
        for term in (
            "bal",
            "bronchoalveolar",
            "sputum",
            "endotracheal",
            "tracheal",
            "eta",
            "lower respiratory",
        )
    ):
        return "Lower_Respiratory"
    if any(
        term in text
        for term in ("nasopharyngeal", "nps", "np swab", "oropharyngeal", "throat")
    ):
        return "Upper_Respiratory"
    if any(
        term in text
        for term in (
            "blood",
            "serum",
            "plasma",
            "pleural",
            "peritoneal",
            "synovial",
            "joint",
            "tissue",
            "bone",
            "sterile fluid",
        )
    ):
        return "Sterile_Site"
    if "urine" in text:
        return "Urinary"
    if any(term in text for term in ("stool", "feces", "faeces", "rectal")):
        return "GI"
    if any(term in text for term in ("wound", "skin", "genital")):
        return "NonSterile_Other"
    return "Unknown"


def _molecular_source(category: str) -> str:
    return {
        "Lower_Respiratory": "Respiratory",
        "Upper_Respiratory": "Respiratory",
        "Sterile_Site": "Systemic",
        "Urinary": "Urinary",
        "GI": "GI",
        "CNS": "CNS",
        "NonSterile_Other": "Other",
    }.get(category, "Unknown")


def _is_molecular_resistance_target(value: Any) -> bool:
    key = pathogen_names.raw_key(value)
    if key in AMR_TARGETS:
        return True
    return bool(
        re.fullmatch(
            r"(?:bla)?(?:ctxm|kpc|ndm|vim|imp|oxa\d*[a-z]*|mec[ac]|mrej|van[ab]|mcr\d+)",
            key,
        )
    )


def _molecular_level(target: str, category: str) -> tuple[str, str]:
    """Assign the fixed molecular evidence level; molecular alone is never Level 1."""
    key = pathogen_names.raw_key(target)
    cautious_lower_respiratory = (
        "cytomegalovirus" in key
        or key == "cmv"
        or key.startswith("hsv")
        or "herpessimplex" in key
        or "candida" in key
    )
    if category in {"Sterile_Site", "CNS"}:
        return "Level 2", "R-MOL-S4-STERILE"
    if category == "Lower_Respiratory":
        if cautious_lower_respiratory:
            return "Level 3", "R-MOL-S4-LOWER-RESP-CAUTIOUS"
        return "Level 2", "R-MOL-S4-LOWER-RESP"
    return "Level 3", "R-MOL-S4-OTHER"


def molecular_microbiology_agent_from_source(source: Any) -> dict[str, Any]:
    """Build a source-complete deterministic targeted-molecular payload.

    Every source row becomes one ``molecular_observations`` record.  Only an
    explicit positive/detected result (or an unambiguously positive numeric
    load) creates an organism label.  Negative, invalid, indeterminate, and
    unknown rows remain available for audit but never become positive evidence.
    Resistance-gene targets are kept separately and never masquerade as an
    organism.
    """
    rows = _rows(source)
    observations: list[dict[str, Any]] = []
    labels: list[dict[str, Any]] = []
    resistance: list[dict[str, Any]] = []
    data_gaps: set[str] = set()

    for index, row in enumerate(rows):
        observation_id = f"MOL-{index + 1:03d}"
        test = _text(row.get("test") or row.get("assay"))
        assay_method = _text(row.get("assay_method") or row.get("method"))
        specimen = _text(
            row.get("sample")
            or row.get("specimen")
            or row.get("specimen_type")
            or row.get("sample_type")
        )
        category = _molecular_specimen_category(specimen)
        target = _text(
            row.get("target")
            or row.get("organism")
            or row.get("organism_name")
            or row.get("pathogen")
        )
        raw_result = _text(
            row.get("result") or row.get("interpretation") or row.get("status")
        )
        value = _text(row.get("value") or row.get("numeric_value"))
        viral_load = _text(
            row.get("viral_load_or_copies")
            or row.get("viral_load")
            or row.get("copies")
        )
        ct_value = _text(row.get("ct_value") or row.get("ct"))
        numeric_source = viral_load or value
        numeric_value, unit = _numeric_value_and_unit(numeric_source)
        detection_status = _assay_detection_status(raw_result)
        if detection_status == "unknown":
            detection_status = _assay_detection_status(value)
        if detection_status == "unknown" and numeric_value:
            try:
                if float(numeric_value) > 0:
                    detection_status = "detected"
            except ValueError:
                pass
        interpretation = {
            "detected": "detected_positive",
            "not_detected": "not_detected",
            "equivocal": "indeterminate",
            "invalid": "invalid",
        }.get(detection_status, "unknown")
        observations.append(
            {
                "observation_id": observation_id,
                "source_record_index": index,
                "source_record_key": _assay_source_record_key(
                    "molecular_microbiology", index, row
                ),
                "test": test,
                "assay_method": assay_method,
                "specimen_type": specimen,
                "specimen_category": category,
                "collected_time": _text(row.get("collected_time")),
                "reported_time": _text(row.get("reported_time")),
                "target": target,
                "raw_result": raw_result,
                "result_interpretation": interpretation,
                "numeric_value": numeric_value,
                "unit": unit,
                "ct_value": ct_value,
                "viral_load_or_copies": viral_load,
                "reference_or_cutoff": _text(
                    row.get("reference") or row.get("cutoff")
                ),
                "source_text": json.dumps(
                    row, ensure_ascii=False, sort_keys=True, separators=(",", ":")
                ),
            }
        )

        if not target:
            data_gaps.add("missing_target")
        if category == "Unknown":
            data_gaps.add("missing_or_unknown_specimen_type")
        if interpretation != "detected_positive" or not target:
            continue
        if _is_molecular_resistance_target(target):
            resistance.append(
                {
                    "gene": target,
                    "linked_organism": "",
                    "link_status": "Unlinked",
                    "rule_applied": "R-MOL-RES-UNLINKED",
                    "evidence_observation_ids": [observation_id],
                }
            )
            continue

        level, level_rule = _molecular_level(target, category)
        labels.append(
            {
                "organism_name": target,
                "classification": _classification(target),
                "causative_level": level,
                "result_interpretation": "detected_positive",
                "applied_rules": [
                    "R-MOL-PRESERVE-01",
                    "R-MOL-S1-01",
                    "R-MOL-S2-01",
                    "R-MOL-S3-01",
                    level_rule,
                ],
                "evidence_observation_ids": [observation_id],
                "key_evidence": {
                    "test": test,
                    "assay_method": assay_method,
                    "specimen_type": specimen,
                    "specimen_category": category,
                    "collected_time": row.get("collected_time") or "Unknown",
                    "reported_time": row.get("reported_time") or "Unknown",
                    "target": target,
                    "result": raw_result or value or viral_load,
                    "numeric_value": numeric_value,
                    "unit": unit,
                    "ct_value": ct_value,
                    "viral_load_or_copies": viral_load,
                    "source_rule": "normalized_source_deterministic_fallback",
                },
            }
        )

    labels = _dedupe_labels(labels)
    highest = min(
        (item["causative_level"] for item in labels),
        key=lambda value: LEVEL_RANK[value],
        default="Not_available",
    )
    top = [item for item in labels if item["causative_level"] == highest]
    classes = {
        item["classification"] for item in top if item["classification"] != "Unknown"
    }
    sources = {
        _molecular_source(str((item.get("key_evidence") or {}).get("specimen_category")))
        for item in top
    }
    sources.discard("Unknown")
    return {
        "rule_version": "molecular_microbiology_normalized_source_fallback_v1",
        "infection_likelihood": _likelihood(highest),
        "probable_source": next(iter(sources)) if len(sources) == 1 else "Other" if sources else "Unknown",
        "probable_pathogen_type": next(iter(classes)) if len(classes) == 1 else "Mixed" if classes else "Unknown",
        "molecular_observations": observations,
        "organism_labels": labels,
        "resistance_findings": sorted(resistance, key=lambda item: item["gene"].lower()),
        "key_findings": ["normalized_source_deterministic_fallback"] if labels else [],
        "data_gaps": sorted(data_gaps),
        "source_record_contract": {
            "source_record_count": len(rows),
            "observation_count": len(observations),
            "one_observation_per_source_record": len(rows) == len(observations),
        },
    }


def _filmarray_sample(value: Any) -> str:
    text = _text(value).lower()
    if any(term in text for term in ("bal", "bronchoalveolar")):
        return "BAL"
    if "endotracheal" in text or re.search(r"\beta\b", text):
        return "ETA"
    if "tracheal" in text or re.search(r"\bta\b", text):
        return "TA"
    if "sputum" in text:
        return "Sputum"
    if any(term in text for term in ("nps", "nasopharyngeal", "np swab")):
        return "NPS"
    if "blood" in text:
        return "Blood_culture"
    return "Unknown"


def _panel_type(rows: list[dict[str, Any]]) -> str:
    haystack = " ".join(_text(value).lower() for row in rows for value in row.values())
    target_keys = {pathogen_names.raw_key(row.get("target")) for row in rows}
    lower_respiratory = any(
        _filmarray_sample(row.get("sample")) in {"BAL", "ETA", "TA", "Sputum"}
        for row in rows
    )
    if "copy/ml" in haystack or (lower_respiratory and bool(target_keys & PNEUMONIA_AMR_TARGETS)):
        return "BioFire FilmArray Pneumonia Panel"
    if bool(target_keys & BCID_AMR_TARGETS) or any(
        term in haystack for term in ("bcid", "blood culture identification", "van a/b", "mcr-1")
    ):
        return "BioFire BCID2"
    if any(term in haystack for term in ("bordetella pertussis", "bordetella parapertussis", "respiratory panel 2.1")):
        return "BioFire Respiratory Panel 2.1"
    return "Unknown"


def _semiquant_bin(value: Any) -> str:
    text = _text(value).lower().replace("≥", ">=")
    match = re.search(r"10\s*\^\s*([4567])", text)
    if not match:
        return "not_reported"
    exponent = int(match.group(1))
    return ">=10^7" if exponent == 7 and any(term in text for term in (">", ">=", "≥")) else f"10^{exponent}"


def filmarray_agent_from_source(source: Any) -> dict[str, Any]:
    rows = _rows(source)
    panel = _panel_type(rows)
    sample = _filmarray_sample(next((_text(row.get("sample")) for row in rows if row.get("sample")), ""))
    labels: list[dict[str, Any]] = []
    resistance: list[dict[str, Any]] = []
    observations: list[dict[str, Any]] = []
    data_gaps: set[str] = set()
    for index, row in enumerate(rows):
        target = _text(row.get("target"))
        result = _text(row.get("result"))
        value = _text(row.get("value"))
        observation_id = f"ASSAY-{index + 1:03d}"
        detection_status = _assay_detection_status(result)
        if (
            detection_status == "unknown"
            and _semiquant_bin(value) != "not_reported"
        ):
            detection_status = "detected"
        numeric_value, unit = _numeric_value_and_unit(value)
        observations.append({
            "observation_id": observation_id,
            "source_record_index": index,
            "source_record_key": _assay_source_record_key(
                "filmarray", index, row
            ),
            "test_type": "filmarray",
            "panel_or_assay": panel,
            "specimen_type": _text(row.get("sample")),
            "collected_time": _text(row.get("collected_time")),
            "reported_time": _text(row.get("reported_time")),
            "target": target,
            "raw_result": result,
            "detection_status": detection_status,
            "numeric_value": numeric_value,
            "unit": unit,
            "reference_or_cutoff": _text(row.get("reference") or row.get("cutoff")),
            "semiquant_bin": _semiquant_bin(value),
            "source_text": json.dumps(
                row, ensure_ascii=False, sort_keys=True, separators=(",", ":")
            ),
        })
        if not target:
            data_gaps.add("missing_target")
            continue
        detected = detection_status == "detected"
        if not detected:
            continue
        target_key = pathogen_names.raw_key(target)
        if target_key in AMR_TARGETS:
            resistance.append(
                {
                    "gene": target,
                    "linked_organism": "",
                    "link_status": "Unlinked",
                    "rule_applied": "R-S6-RES-UNLINKED",
                    "evidence_observation_ids": [observation_id],
                }
            )
            continue
        classification = _classification(target)
        semiquant = _semiquant_bin(value)
        lower_resp = sample in {"BAL", "ETA", "TA", "Sputum"}
        if panel == "BioFire FilmArray Pneumonia Panel" and classification == "Bacterial" and semiquant != "not_reported":
            level = {">=10^7": "Level 1", "10^7": "Level 1", "10^6": "Level 2", "10^5": "Level 3", "10^4": "Level 4"}[semiquant]
            level_rule = "R-S4-01"
            quantitation = "reported"
        elif panel == "BioFire FilmArray Pneumonia Panel" and classification == "Bacterial" and lower_resp and _is_typical_bacterium(target):
            level = "Level 2"
            level_rule = "R-S4-01B"
            quantitation = "positive_detected_quantity_unknown"
            data_gaps.add("missing_semiquant_bin")
        elif panel == "BioFire FilmArray Pneumonia Panel":
            level = "Level 2" if lower_resp else "Level 3"
            level_rule = "R-S4-02"
            quantitation = "qualitative_by_design"
        elif panel == "BioFire BCID2":
            level = "Level 2"
            level_rule = "R-S4-04"
            quantitation = "qualitative_by_design"
        else:
            level = "Level 3"
            level_rule = "R-S4-03"
            quantitation = "qualitative_by_design"
        labels.append(
            {
                "organism_name": target,
                "classification": classification,
                "causative_level": level,
                "evidence_observation_ids": [observation_id],
                "key_evidence": {
                    "filmarray_panel": panel,
                    "filmarray_sample_type": sample,
                    "detection_status": "detected",
                    "semiquant_bin": semiquant if semiquant != "not_reported" else "not_applicable" if classification != "Bacterial" else "not_reported",
                    "quantitation_status": quantitation,
                    "non_gm_test_type": "none",
                    "evidence_type": "direct_detection",
                    "reported_time": row.get("reported_time") or "Unknown",
                    "source_rule": "normalized_source_deterministic_fallback",
                },
                "decision_trace": [
                    {
                        "step": "S2_extract_positive",
                        "rule_id": "R-S2-01",
                        "result": "pass",
                        "evidence": [
                            {"field": f"filmarray[{index}].target", "value": target},
                            {"field": f"filmarray[{index}].result", "value": result},
                            {"field": f"filmarray[{index}].value", "value": value},
                        ],
                    },
                    {
                        "step": "S4_filmarray_level",
                        "rule_id": level_rule,
                        "result": "pass",
                        "evidence": [{"field": "assigned_level", "value": level}],
                    },
                ],
                "causative_reasoning": [f"Deterministic normalized-source fallback applied {level_rule}."],
            }
        )
    labels = _dedupe_labels(labels)
    highest = min((item["causative_level"] for item in labels), key=lambda value: LEVEL_RANK[value], default="Not_available")
    top = [item for item in labels if item["causative_level"] == highest]
    classes = {item["classification"] for item in top if item["classification"] != "Unknown"}
    source_name = "Respiratory" if sample in {"BAL", "ETA", "TA", "Sputum", "NP", "OP", "NPS"} else "Bloodstream" if sample == "Blood_culture" else "Unknown"
    return {
        "rule_version": "filmarray_normalized_source_fallback_v2",
        "assay_summary": {
            "filmarray_panel": panel,
            "filmarray_sample_type": sample,
            "gm_test_sample_type": "Unknown",
            "gm_test_cutoff_used": {"bal_positive": 1.0, "bal_borderline": [0.5, 0.99], "serum_positive": 0.5},
        },
        "infection_likelihood": _likelihood(highest),
        "probable_source": source_name,
        "probable_pathogen_type": next(iter(classes)) if len(classes) == 1 else "Mixed" if classes else "Unknown",
        "assay_observations": observations,
        "organism_labels": labels,
        "resistance_risk": "Moderate" if resistance else "Low",
        "resistance_findings": sorted(resistance, key=lambda item: item["gene"].lower()),
        "key_findings": ["normalized_source_deterministic_fallback"],
        "data_gaps": sorted(data_gaps),
        "source_record_contract": {
            "filmarray_source_record_count": len(rows),
            "gm_source_record_count": 0,
            "source_record_count": len(rows),
            "observation_count": len(observations),
            "one_observation_per_source_record": len(rows) == len(observations),
        },
    }


def filmarray_gm_agent_from_sources(
    filmarray_source: Any,
    gm_source: Any,
) -> dict[str, Any]:
    """Build one source-complete deterministic FilmArray/GM payload.

    Every input row becomes exactly one ``assay_observations`` record.  Only
    positive FilmArray targets, positive GM values, and positive non-GM
    serology/antigen rows create organism labels.  Negative rows remain
    traceable observations and are never converted into counter-evidence.
    """
    output = filmarray_agent_from_source(filmarray_source)
    filmarray_rows = _rows(filmarray_source)
    gm_rows = _rows(gm_source)
    observations = list(output.get("assay_observations") or [])
    labels = list(output.get("organism_labels") or [])
    data_gaps = set(output.get("data_gaps") or [])
    gm_samples: set[str] = set()

    for gm_index, row in enumerate(gm_rows):
        source_index = len(filmarray_rows) + gm_index
        observation_id = f"ASSAY-{source_index + 1:03d}"
        target = _text(
            row.get("target") or row.get("test") or row.get("assay") or row.get("name")
        )
        specimen = _text(
            row.get("sample")
            or row.get("specimen")
            or row.get("specimen_type")
            or row.get("sample_type")
        )
        result = _text(row.get("result") or row.get("interpretation"))
        value = _text(row.get("value"))
        numeric_value, unit = _numeric_value_and_unit(value)
        numeric = _float_or_none(value)
        target_text = target.casefold()
        specimen_text = specimen.casefold()
        is_gm = "galactomannan" in target_text or "gm test" in target_text
        if is_gm:
            test_type = "galactomannan"
            if "bal" in specimen_text or "bronchoalveolar" in specimen_text:
                gm_sample = "BAL"
            elif any(term in specimen_text for term in ("serum", "blood", "plasma")):
                gm_sample = "Serum"
            else:
                gm_sample = "Unknown"
            gm_samples.add(gm_sample)
            positive = bool(numeric is not None and numeric >= 0.5)
            detection_status = "detected" if positive else (
                "not_detected" if numeric is not None else _assay_detection_status(result)
            )
            reference_or_cutoff = "BAL/Serum >=0.5 detection; BAL >=1.0 strong"
        else:
            non_gm_type, _ = _non_gm_test_type(target)
            test_type = (
                "serology" if non_gm_type.startswith("serology_")
                else "antigen" if non_gm_type == "antigen"
                else "other"
            )
            detection_status = _assay_detection_status(result or value)
            positive = detection_status == "detected"
            reference_or_cutoff = _text(row.get("reference") or row.get("cutoff"))

        observations.append({
            "observation_id": observation_id,
            "source_record_index": gm_index,
            "source_record_key": _assay_source_record_key(
                "gm_test", gm_index, row
            ),
            "test_type": test_type,
            "panel_or_assay": target,
            "specimen_type": specimen,
            "collected_time": _text(row.get("collected_time")),
            "reported_time": _text(row.get("reported_time")),
            "target": target,
            "raw_result": result or value,
            "detection_status": detection_status,
            "numeric_value": numeric_value,
            "unit": unit,
            "reference_or_cutoff": reference_or_cutoff,
            "semiquant_bin": "not_applicable",
            "source_text": json.dumps(
                row, ensure_ascii=False, sort_keys=True, separators=(",", ":")
            ),
        })
        if not target:
            data_gaps.add("missing_gm_or_antigen_target")
            continue
        if not positive:
            continue

        if is_gm:
            if gm_sample == "BAL" and numeric is not None and numeric >= 1.0:
                level = "Level 1"
                rule = "R-S5-02-BAL-STRONG"
            elif gm_sample == "BAL" and numeric is not None and numeric >= 0.5:
                level = "Level 2"
                rule = "R-S5-02-BAL-BORDERLINE"
            else:
                level = "Level 3"
                rule = "R-S5-02-SERUM-OR-UNKNOWN"
            organism_name = "Aspergillus spp."
            classification = "Fungal"
            non_gm_type = "none"
            evidence_type = "direct_detection"
        else:
            non_gm_type, evidence_type = _non_gm_test_type(target)
            organism_name = _normalized_non_gm_target(target)
            if not organism_name:
                data_gaps.add("non_gm_target_normalization_empty")
                continue
            classification = _classification(organism_name)
            level = "Level 4" if evidence_type == "indirect_serology" else "Level 3"
            rule = "R-S5-03-NON-GM-SEROLOGY" if level == "Level 4" else "R-S5-04-NON-GM-ANTIGEN"
            if evidence_type == "indirect_serology":
                data_gaps.update({
                    "non_gm_serology_indirect_evidence",
                    "paired_serology_or_titer_not_available",
                })

        labels.append({
            "organism_name": organism_name,
            "classification": classification,
            "causative_level": level,
            "evidence_observation_ids": [observation_id],
            "key_evidence": {
                "filmarray_panel": "Unknown",
                "filmarray_sample_type": "Unknown",
                "gm_test_sample_type": specimen or "Unknown",
                "gm_test_target": target,
                "gm_test_result": result or value,
                "detection_status": "detected",
                "semiquant_bin": "not_applicable",
                "quantitation_status": "reported" if numeric is not None else "qualitative_by_design",
                "non_gm_test_type": non_gm_type,
                "evidence_type": evidence_type,
                "reported_time": row.get("reported_time") or "Unknown",
                "source_rule": "normalized_source_deterministic_fallback",
            },
            "decision_trace": [{
                "step": "S5_gm_level",
                "rule_id": rule,
                "result": "pass",
                "evidence": [
                    {"field": f"gm_test[{gm_index}].target", "value": target},
                    {"field": f"gm_test[{gm_index}].result", "value": result or value},
                    {"field": f"gm_test[{gm_index}].sample", "value": specimen},
                ],
            }],
            "causative_reasoning": [
                f"Deterministic normalized-source fallback applied {rule}."
            ],
        })

    labels = _dedupe_labels(labels)
    highest = min(
        (item["causative_level"] for item in labels),
        key=lambda value: LEVEL_RANK[value],
        default="Not_available",
    )
    top = [item for item in labels if item["causative_level"] == highest]
    classes = {
        item["classification"] for item in top if item["classification"] != "Unknown"
    }
    output["rule_version"] = "filmarray_gm_normalized_source_fallback_v2"
    output["assay_observations"] = observations
    output["organism_labels"] = labels
    output["infection_likelihood"] = _likelihood(highest)
    output["probable_pathogen_type"] = (
        next(iter(classes)) if len(classes) == 1 else "Mixed" if classes else "Unknown"
    )
    if gm_samples == {"BAL", "Serum"}:
        gm_sample_summary = "Both"
    elif len(gm_samples) == 1:
        gm_sample_summary = next(iter(gm_samples))
    else:
        gm_sample_summary = "Unknown"
    output.setdefault("assay_summary", {})["gm_test_sample_type"] = gm_sample_summary
    output["data_gaps"] = sorted(data_gaps)
    output["source_record_contract"] = {
        "filmarray_source_record_count": len(filmarray_rows),
        "gm_source_record_count": len(gm_rows),
        "source_record_count": len(filmarray_rows) + len(gm_rows),
        "observation_count": len(observations),
        "one_observation_per_source_record": (
            len(filmarray_rows) + len(gm_rows) == len(observations)
        ),
    }
    return output

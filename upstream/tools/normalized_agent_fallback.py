"""Deterministic agent-shaped fallbacks for normalized microbiology sources.

These helpers are used only when a small-agent output is absent. They implement
the fixed culture and FilmArray rules already documented in the active prompts;
they do not use benchmark answers or call an LLM.
"""

from __future__ import annotations

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
    profile = taxonomy.classify_organism(name, biological_class="Bacterial")
    return profile.get("primary_rule_family") in {
        "typical_respiratory_pathogen",
        "hospital_or_nonfermenter_gnb",
    }


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


def culture_agent_from_source(source: Any) -> dict[str, Any]:
    labels: list[dict[str, Any]] = []
    data_gaps: set[str] = set()
    for row in _rows(source):
        name = _text(row.get("organism"))
        if not name:
            continue
        status = _text(row.get("status"))
        if not (_positive_text(status) or (not status and name)):
            continue
        specimen_type = _text(row.get("sample") or row.get("specimen"))
        specimen = _specimen_category(specimen_type)
        if specimen == "Unknown":
            data_gaps.add("missing_specimen_type")
        quantity, quantitation = _quantity_tier(
            row.get("colony_count") or row.get("quantity") or row.get("value")
        )
        if quantity == "Unknown":
            data_gaps.add("missing_quantity")
        purity = _growth_purity(row)
        level, applied, contaminant, hard = _culture_level(name, specimen, quantity, purity)
        labels.append(
            {
                "organism_name": name,
                "classification": _classification(name),
                "causative_level": level,
                "applied_rules": applied,
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
        "organism_labels": labels,
        "module_rule_summary": {
            "highest_level": highest,
            "source_rule_id": "R-S5-03",
            "pathogen_type_rule_id": "R-S5-04",
        },
        "data_gaps": sorted(data_gaps),
    }


def _dedupe_labels(labels: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: dict[tuple[str, str], dict[str, Any]] = {}
    for item in labels:
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
        if item_order < current_order:
            selected[key] = item
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
    data_gaps: set[str] = set()
    for index, row in enumerate(rows):
        target = _text(row.get("target"))
        if not target:
            continue
        result = _text(row.get("result"))
        value = _text(row.get("value"))
        detected = _positive_text(result) or (bool(_semiquant_bin(value) != "not_reported") and "not detected" not in result.lower())
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
        "rule_version": "filmarray_normalized_source_fallback_v1",
        "assay_summary": {
            "filmarray_panel": panel,
            "filmarray_sample_type": sample,
            "gm_test_sample_type": "Unknown",
            "gm_test_cutoff_used": {"bal_positive": 1.0, "bal_borderline": [0.5, 0.99], "serum_positive": 0.5},
        },
        "infection_likelihood": _likelihood(highest),
        "probable_source": source_name,
        "probable_pathogen_type": next(iter(classes)) if len(classes) == 1 else "Mixed" if classes else "Unknown",
        "organism_labels": labels,
        "resistance_risk": "Moderate" if resistance else "Low",
        "resistance_findings": sorted(resistance, key=lambda item: item["gene"].lower()),
        "key_findings": ["normalized_source_deterministic_fallback"],
        "data_gaps": sorted(data_gaps),
    }

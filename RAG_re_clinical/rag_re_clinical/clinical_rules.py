"""Pre-specified clinical precision rules for rescue-only experiments.

The functions in this module are deliberately pure: they do not read files, call
models, mutate candidates, or change the frozen E0 baseline.  They only turn an
already-merged case and candidate into auditable module results.

Important semantics
-------------------
* ``None`` means unknown/abstain.  Missing, pending, ordered and not-performed
  tests are never converted to a negative result.
* Structured, candidate-specific direct evidence is stronger than an upstream
  generic ``Aligned`` label.  An explicit site mismatch cannot be overwritten by
  that label.
* A module name or support summary is context only (C1 at most); C2/C3 require a
  structured exact-organism record with a usable specimen site.
* The engine, not this module, owns E0 immutability.
"""

from __future__ import annotations

from datetime import date, datetime
import re
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Sequence


TriState = bool | None


UNKNOWN_MARKERS = {
    "",
    "NA",
    "N_A",
    "NONE",
    "NULL",
    "MISSING",
    "UNKNOWN",
    "UNAVAILABLE",
    "NOT_AVAILABLE",
    "PENDING",
    "ORDERED",
    "NOT_DONE",
    "NOT_PERFORMED",
    "NOT_TESTED",
    "INDETERMINATE",
}

POSITIVE_MARKERS = {
    "POSITIVE",
    "DETECTED",
    "PRESENT",
    "REACTIVE",
    "CONFIRMED",
    "GROWTH",
    "ISOLATED",
    "REPORTED",
    "SUPPORT",
    "SUPPORTED",
    "TRUE",
    "PASS",
}

NEGATIVE_MARKERS = {
    "NEGATIVE",
    "NOT_DETECTED",
    "NO_GROWTH",
    "ABSENT",
    "NONREACTIVE",
    "NON_REACTIVE",
    "FALSE",
    "FAIL",
}

STATUS_KEYS = {
    "STATUS",
    "RESULT",
    "TEST_RESULT",
    "CULTURE_RESULT",
    "DETECTION_STATUS",
    "DETECTED",
    "POSITIVE",
    "VERDICT",
    "QUANTITATION_STATUS",
}

SITE_KEYS = {
    "SITE",
    "SITE_CODE",
    "SPECIMEN_SITE",
    "SPECIMEN_TYPE",
    "SPECIMEN_CATEGORY",
    "COLLECTION_SITE",
    "SAMPLE_TYPE",
    "BODY_SITE",
}

DIRECT_MODULES = {
    "culture",
    "filmarray_gmtest",
    "targeted_molecular",
    "molecular_microbiology",
    "pcr",
    "antigen",
    "gm_test",
}

HIGH_SPECIFICITY_MODULES = {
    "culture",
    "targeted_molecular",
    "molecular_microbiology",
    "pcr",
    "antigen",
}


def _record_contaminated(record: Mapping[str, Any]) -> TriState:
    """Parse only explicit contamination/poor-quality indicators."""

    for key in (
        "contaminant_flag",
        "contamination_flag",
        "is_contaminant",
        "likely_contaminant",
    ):
        if key not in record:
            continue
        value = record.get(key)
        if value is True:
            return True
        if value is False:
            return False
        token = _norm(value)
        if token in {"YES", "Y", "TRUE", "LIKELY", "CONTAMINANT", "POSITIVE"}:
            return True
        if token in {"NO", "N", "FALSE", "NOT_CONTAMINANT", "UNLIKELY"}:
            return False
    purity = _norm(record.get("growth_purity"))
    if "CONTAMIN" in purity or "POOR_QUALITY" in purity:
        return True
    return None

SKIP_EVIDENCE_SUBTREES = {
    "SUPPORT_MODULES",
    "NON_HOST_SUPPORT_MODULES",
    "MODULE_SUPPORT_SUMMARY",
    "MODULE_LEVEL_SUMMARY",
    "HOSPITAL_MODULE_LEVEL_SUMMARY",
    "RELATED_REPRESENTATIVE_HOSPITAL_SUPPORT",
    "RELATED_REPRESENTATIVE_CONTEXT",
    "HOSPITAL_EVIDENCE_TRACE",
    "LITERATURE_EVIDENCE",
}


def _norm(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", str(value or "").upper()).strip("_")


def _name_key(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]+", "", str(value or "").upper())


def _as_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_int(value: Any) -> int | None:
    match = re.search(r"-?\d+", str(value or ""))
    return int(match.group(0)) if match else None


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _first_present(mapping: Mapping[str, Any], keys: Iterable[str]) -> Any:
    normalized = {_norm(key): value for key, value in mapping.items()}
    for key in keys:
        value = normalized.get(_norm(key))
        if value not in (None, "", [], {}):
            return value
    return None


def _tri_and(*values: TriState) -> TriState:
    if any(value is False for value in values):
        return False
    if all(value is True for value in values):
        return True
    return None


def normalize_site(value: Any) -> str | None:
    """Normalize common syndrome/specimen labels to a small site vocabulary."""

    text = _norm(value)
    if text in UNKNOWN_MARKERS:
        return None

    # Specific anatomy is checked before generic words such as "sterile".
    if any(
        token in text
        for token in (
            "LOWER_RESP",
            "PULMON",
            "LUNG",
            "BRONCH",
            "SPUTUM",
            "TRACHEAL",
            "ENDOTRACHEAL",
            "BALF",
            "BAL",
            "PLEURAL",
            "PNEUMON",
            "S2_",
        )
    ):
        return "lower_respiratory"
    if any(
        token in text
        for token in ("UPPER_RESP", "NASOPHARY", "NASAL", "THROAT", "OROPHARY")
    ):
        return "upper_respiratory"
    if any(
        token in text
        for token in (
            "BLOOD",
            "SERUM",
            "PLASMA",
            "BACTEREM",
            "SEPSIS",
            "SYSTEMIC",
            "S1_",
        )
    ):
        return "bloodstream"
    if any(token in text for token in ("CSF", "CEREBROSPINAL", "MENING", "BRAIN", "CNS")):
        return "cns"
    if any(token in text for token in ("URINE", "URINARY", "FOLEY", "CATHETER_URINE", "UTI")):
        return "urinary"
    if any(token in text for token in ("INTRA_ABDOM", "ABDOM", "ASCITES", "PERITONE")):
        return "intra_abdominal"
    if any(token in text for token in ("SKIN", "WOUND", "CUTANEOUS")):
        return "skin_soft_tissue"
    if any(token in text for token in ("TISSUE", "BIOPSY", "HISTOPATH")):
        return "tissue"
    if any(token in text for token in ("STERILE_SITE", "STERILE_FLUID", "JOINT_FLUID", "SYNOVIAL")):
        return "sterile_other"
    if "RESPIRATORY" in text or text == "RESP":
        return "lower_respiratory"
    return None


def _candidate_input(frozen_candidate: Mapping[str, Any]) -> Dict[str, Any]:
    """Return a copy of the upstream rule fields without modifying the input."""

    combined: Dict[str, Any] = {}
    modules = _mapping(frozen_candidate.get("modules"))
    b_features = _mapping(_mapping(modules.get("B")).get("features"))
    combined.update(b_features)
    combined.update(_mapping(frozen_candidate.get("rule_input")))
    for key, value in frozen_candidate.items():
        if key not in {"rule_input", "modules"} and value is not None:
            combined[key] = value
    return combined


def _source_context(merged_context: Mapping[str, Any], frozen_candidate: Mapping[str, Any]) -> Dict[str, Any]:
    context: Dict[str, Any] = {}
    context.update(_mapping(merged_context.get("source_context")))
    context.update(_mapping(frozen_candidate.get("source_context")))
    # ``merged_context`` may itself be the compact source_context object.
    for key in (
        "target_site",
        "site_code",
        "specimen_site_code",
        "specimen_class",
        "specimen_alignment",
        "dominant_source",
        "clinical_site_source",
    ):
        if merged_context.get(key) not in (None, "") and key not in context:
            context[key] = merged_context.get(key)
    return context


def _candidate_name_keys(
    frozen_candidate: Mapping[str, Any], config: Mapping[str, Any] | None = None
) -> set[str]:
    candidate = _candidate_input(frozen_candidate)
    values = [
        frozen_candidate.get("organism_name"),
        frozen_candidate.get("canonical_organism_name"),
        candidate.get("organism_name"),
        candidate.get("pathogen"),
        candidate.get("pathogen_name"),
        candidate.get("canonical_key"),
    ]
    result = {_name_key(value) for value in values if _name_key(value)}
    aliases = _mapping((config or {}).get("exact_aliases"))
    changed = True
    while changed:
        changed = False
        for left, right_values in aliases.items():
            right_seq = right_values if isinstance(right_values, (list, tuple, set)) else [right_values]
            group = {_name_key(left), *(_name_key(value) for value in right_seq)}
            group.discard("")
            if result.intersection(group) and not group.issubset(result):
                result.update(group)
                changed = True
    return result


def _row_name_key(row: Mapping[str, Any]) -> str:
    for key in ("canonical_key", "canonical_organism_name", "organism_name", "pathogen_name", "pathogen", "name"):
        value = row.get(key)
        if value not in (None, ""):
            return _name_key(value)
    return ""


def _matching_candidate_sources(
    merged_context: Mapping[str, Any],
    frozen_candidate: Mapping[str, Any],
    config: Mapping[str, Any] | None = None,
) -> List[Mapping[str, Any]]:
    keys = _candidate_name_keys(frozen_candidate, config)
    sources: List[Mapping[str, Any]] = [frozen_candidate]
    rule_input = _mapping(frozen_candidate.get("rule_input"))
    if rule_input:
        sources.append(rule_input)

    def walk(value: Any) -> None:
        if isinstance(value, Mapping):
            row_key = _row_name_key(value)
            if row_key and row_key in keys:
                sources.append(value)
            for nested in value.values():
                walk(nested)
        elif isinstance(value, (list, tuple)):
            for nested in value:
                walk(nested)

    walk(merged_context)
    deduped: List[Mapping[str, Any]] = []
    seen: set[int] = set()
    for source in sources:
        marker = id(source)
        if marker not in seen:
            deduped.append(source)
            seen.add(marker)
    return deduped


def _module_from_key(key: Any) -> str | None:
    normalized = _norm(key)
    if normalized == "CULTURE" or normalized.endswith("_CULTURE"):
        return "culture"
    if "FILMARRAY" in normalized:
        return "filmarray_gmtest"
    if normalized in {"TARGETED_MOLECULAR", "TARGETED_PCR", "TARGETED_TEST"}:
        return "targeted_molecular"
    if normalized in {"MOLECULAR_MICROBIOLOGY", "MOLECULAR"}:
        return "molecular_microbiology"
    if normalized in {"PCR", "QPCR", "RT_PCR", "VIRAL_LOAD"} or normalized.endswith("_PCR"):
        return "pcr"
    if "ANTIGEN" in normalized:
        return "antigen"
    if normalized in {"GM", "GM_TEST", "GALACTOMANNAN"}:
        return "gm_test"
    return None


def _looks_like_evidence_record(record: Mapping[str, Any]) -> bool:
    keys = {_norm(key) for key in record}
    if keys.intersection(STATUS_KEYS | SITE_KEYS):
        return True
    return bool(
        keys.intersection(
            {
                "QUANTITY_TIER",
                "GROWTH_PURITY",
                "SEMIQUANT_BIN",
                "EVIDENCE_TYPE",
                "FILMARRAY_PANEL",
                "TEST_TYPE",
                "ASSAY",
                "CT_VALUE",
                "VIRAL_LOAD",
            }
        )
    )


def _status_from_record(record: Mapping[str, Any], path: Sequence[str]) -> str:
    statuses: List[str] = []
    for key, value in record.items():
        normalized_key = _norm(key)
        if normalized_key not in STATUS_KEYS:
            continue
        if value is True:
            statuses.append("POSITIVE")
        elif value is False:
            statuses.append("NEGATIVE")
        else:
            statuses.append(_norm(value))

    def classified(value: str) -> str | None:
        if value in UNKNOWN_MARKERS or any(
            token in value
            for token in ("NOT_AVAILABLE", "NOT_PERFORMED", "NOT_DONE", "PENDING", "ORDERED", "MISSING")
        ):
            return "unknown"
        # Negative phrases must be tested before the substring DETECTED.
        if value in NEGATIVE_MARKERS or any(
            token in value for token in ("NOT_DETECTED", "NO_GROWTH", "NEGATIVE", "NONREACTIVE")
        ):
            return "negative"
        if value in POSITIVE_MARKERS or any(
            token in value for token in ("POSITIVE", "DETECTED", "ISOLATED", "CONFIRMED")
        ):
            return "positive"
        return None

    classes = [classified(value) for value in statuses]
    classes = [value for value in classes if value]
    if "positive" in classes and "negative" in classes:
        return "discordant"
    if "positive" in classes:
        return "positive"
    if "negative" in classes:
        return "negative"
    if "unknown" in classes:
        return "unknown"

    # In current merge files, a row under candidate-specific hospital evidence is
    # already a matched positive row; culture rows often have no explicit result
    # key.  This structural inference is intentionally unavailable to summaries.
    structural_parent = any(
        part in {"HOSPITAL_MODULE_EVIDENCE", "HOSPITAL_EVIDENCE_DETAIL"}
        for part in path
    )
    if structural_parent and record:
        return "positive"
    return "unknown"


def _record_sites(record: Mapping[str, Any]) -> tuple[List[str], List[str]]:
    sites: List[str] = []
    raw_values: List[str] = []
    for key, value in record.items():
        normalized_key = _norm(key)
        if normalized_key in SITE_KEYS or normalized_key.endswith(("_SPECIMEN_TYPE", "_SAMPLE_TYPE", "_SPECIMEN_SITE")):
            raw = str(value or "").strip()
            site = normalize_site(raw)
            if raw:
                raw_values.append(raw)
            if site and site not in sites:
                sites.append(site)
    return sites, raw_values


def _record_is_sterile(record: Mapping[str, Any], sites: Sequence[str]) -> bool:
    text = " ".join(_norm(value) for value in record.values() if not isinstance(value, (dict, list)))
    if any(token in text for token in ("NONSTERILE", "NON_STERILE", "LOWER_RESPIRATORY", "SPUTUM", "BAL")):
        return False
    if any(site in {"bloodstream", "cns", "tissue", "sterile_other"} for site in sites):
        return True
    return any(
        token in text
        for token in ("STERILE_SITE", "TISSUE", "BIOPSY", "HISTOPATH", "PLEURAL_FLUID", "ABSCESS")
    )


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return tuple(sorted((_norm(key), _freeze(nested)) for key, nested in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return str(value)


def _structured_direct_records(
    merged_context: Mapping[str, Any],
    frozen_candidate: Mapping[str, Any],
    config: Mapping[str, Any] | None = None,
) -> List[Dict[str, Any]]:
    name_keys = _candidate_name_keys(frozen_candidate, config)
    records: List[Dict[str, Any]] = []

    def walk(value: Any, path: tuple[str, ...], active_module: str | None) -> None:
        if isinstance(value, Mapping):
            support_type = _norm(value.get("support_type"))
            related = "RELATED_REPRESENTATIVE" in support_type or any(
                "RELATED_REPRESENTATIVE" in part for part in path
            )
            if related:
                return
            if active_module and _looks_like_evidence_record(value):
                row_name = _row_name_key(value)
                exact = not row_name or row_name in name_keys
                sites, raw_sites = _record_sites(value)
                records.append(
                    {
                        "module": active_module,
                        "status": _status_from_record(value, path),
                        "exact_organism": exact,
                        "organism_explicit": bool(row_name) and exact,
                        "sites": sites,
                        "raw_sites": raw_sites,
                        "sterile": _record_is_sterile(value, sites),
                        "contaminated": _record_contaminated(value),
                        "path": ".".join(path),
                        "record": dict(value),
                    }
                )
                return
            for key, nested in value.items():
                normalized_key = _norm(key)
                if normalized_key in SKIP_EVIDENCE_SUBTREES or "RELATED_REPRESENTATIVE" in normalized_key:
                    continue
                module = _module_from_key(key) or active_module
                walk(nested, (*path, normalized_key), module)
        elif isinstance(value, (list, tuple)):
            for index, nested in enumerate(value):
                walk(nested, (*path, str(index)), active_module)

    for source in _matching_candidate_sources(merged_context, frozen_candidate, config):
        walk(source, (), None)

    deduped: List[Dict[str, Any]] = []
    seen = set()
    for row in records:
        marker = (
            row["module"],
            row["status"],
            row["exact_organism"],
            row["organism_explicit"],
            tuple(row["sites"]),
            row["sterile"],
            row["contaminated"],
            _freeze(row["record"]),
        )
        if marker not in seen:
            deduped.append(row)
            seen.add(marker)
    return deduped


def _target_site(
    merged_context: Mapping[str, Any],
    frozen_candidate: Mapping[str, Any],
    config: Mapping[str, Any] | None = None,
) -> str | None:
    override = _first_present(config or {}, ("target_site", "site_code"))
    if normalize_site(override):
        return normalize_site(override)
    context = _source_context(merged_context, frozen_candidate)
    for key in ("site_code", "target_site", "dominant_source"):
        site = normalize_site(context.get(key))
        if site:
            return site
    candidate = _candidate_input(frozen_candidate)
    for container in (merged_context, _mapping(merged_context.get("deterministic_max")), candidate):
        for key in ("clinical_target_site", "infection_site", "target_site", "dominant_source"):
            site = normalize_site(container.get(key))
            if site:
                return site
    return None


def _mngs_site_alignment(
    merged_context: Mapping[str, Any], frozen_candidate: Mapping[str, Any]
) -> Dict[str, Any]:
    candidate = _candidate_input(frozen_candidate)
    context = _source_context(merged_context, frozen_candidate)
    target = _target_site(merged_context, frozen_candidate)
    specimen = normalize_site(context.get("specimen_site_code")) or normalize_site(
        candidate.get("specimen_class")
    )
    alignment = _norm(candidate.get("specimen_alignment") or context.get("specimen_alignment"))
    if target and specimen:
        aligned = target == specimen
        return {
            "aligned": aligned,
            "target_site": target,
            "specimen_site": specimen,
            "basis": f"explicit mNGS specimen site {specimen} {'matches' if aligned else 'does not match'} {target}",
        }
    if alignment in {"NOT_ALIGNED", "MISMATCH", "MISALIGNED", "NEGATIVE"}:
        return {
            "aligned": False,
            "target_site": target,
            "specimen_site": specimen,
            "basis": f"specimen_alignment={alignment}",
        }
    if alignment in {"ALIGNED", "MATCHED", "STERILE_OR_SYSTEMIC"}:
        return {
            "aligned": True,
            "target_site": target,
            "specimen_site": specimen,
            "basis": f"specimen_alignment={alignment}",
        }
    return {
        "aligned": None,
        "target_site": target,
        "specimen_site": specimen,
        "basis": "mNGS target/specimen site or alignment is missing/unrecognized",
    }


def _guardrail_values(candidate: Mapping[str, Any]) -> List[str]:
    values: List[str] = []

    def is_guard_rule(value: Any) -> bool:
        text = _norm(value)
        return any(
            token in text
            for token in (
                "GUARD",
                "EXCLUSION",
                "CONTAMIN",
                "COLONIZER",
                "BACKGROUND",
                "CANDIDA_RESP",
                "HERPES_REACTIVATION",
                "SKIN_FLORA",
                "ORAL_ASPIRATION",
                "WATER_ENV",
                "IMPOSSIBLE",
                "NOT_PICKED",
                "OMITTED",
            )
        )

    def walk(value: Any, key_hint: str = "") -> None:
        if isinstance(value, Mapping):
            for key, nested in value.items():
                normalized = _norm(key)
                exclusion_key = "EXCLUSION_RULE" in normalized
                state_flag_key = any(
                    token in normalized
                    for token in ("CONTAMIN", "COLONIZER", "BACKGROUND")
                )
                guard_key = "GUARD" in normalized
                if exclusion_key or state_flag_key or guard_key:
                    if nested is True:
                        values.append(str(nested))
                    elif (
                        not isinstance(nested, bool)
                        and isinstance(nested, (str, int, float))
                        and str(nested).strip()
                    ):
                        token = _norm(nested)
                        negative_or_unknown = token in {
                            "NO",
                            "N",
                            "FALSE",
                            "NEGATIVE",
                            "NONE",
                            "UNKNOWN",
                            "NOT_AVAILABLE",
                            "NOT_APPLICABLE",
                        }
                        if not negative_or_unknown and (
                            exclusion_key
                            or is_guard_rule(nested)
                            or (
                                state_flag_key
                                and token
                                in {"YES", "Y", "TRUE", "LIKELY", "POSITIVE"}
                            )
                        ):
                            values.append(str(nested))
                walk(nested, normalized)
        elif isinstance(value, (list, tuple)) and any(
            token in key_hint for token in ("RULE", "GUARD", "APPLIED")
        ):
            for nested in value:
                if str(nested).strip() and is_guard_rule(nested):
                    values.append(str(nested))

    walk(candidate)
    return values


def specimen_colonization_gate(
    merged_context: Mapping[str, Any],
    frozen_candidate: Mapping[str, Any],
    config: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    """Module D: require site coherence and surface colonization guardrails."""

    cfg = config or {}
    candidate = _candidate_input(frozen_candidate)
    target = _target_site(merged_context, frozen_candidate, cfg)
    generic = _mngs_site_alignment(merged_context, frozen_candidate)
    records = _structured_direct_records(merged_context, frozen_candidate, cfg)
    positive_direct_sites = sorted(
        {
            site
            for row in records
            if row["status"] == "positive" and row["exact_organism"]
            for site in row["sites"]
        }
    )

    # The frozen RAG_re artifact intentionally carries a narrow rule-input
    # allowlist.  Colonizer/background flags and formal exclusion rules still
    # live in the exact-matched deterministic candidate, so D must inspect the
    # adapter's candidate-specific layers as well.  This remains patient/name
    # scoped and never uses genus/fuzzy transfer.
    guard_sources = _matching_candidate_sources(
        merged_context, frozen_candidate, cfg
    )
    guard_values = list(
        dict.fromkeys(
            value
            for source in guard_sources
            for value in _guardrail_values(source)
        )
    )
    guard_text = " ".join(_norm(value) for value in guard_values)
    organism = _norm(
        frozen_candidate.get("canonical_organism_name")
        or frozen_candidate.get("organism_name")
        or candidate.get("organism_name")
    )
    guard_reasons: List[str] = []
    hard_block = False
    raw_colonizer_flag = any(
        source.get("is_likely_colonizer_or_background") is True
        for source in guard_sources
    )
    if raw_colonizer_flag:
        guard_reasons.append("upstream colonizer/background flag")
    if any(token in guard_text for token in ("IMPOSSIBLE", "CONTAMINANT", "SKIN_FLORA_HARD", "CONS_SKIN_FLORA_HARD")):
        hard_block = True
        guard_reasons.append("hard contamination/impossible-source rule")
    elif guard_text:
        guard_reasons.append("upstream guardrail rule")
    if target == "lower_respiratory" and ("CANDIDA" in organism or organism in {"YEAST", "GENERICYEAST"}):
        guard_reasons.append("respiratory Candida/yeast colonization prior")
    if target == "lower_respiratory" and any(
        token in organism for token in ("HERPESSIMPLEX", "HSV1", "HSV2", "CYTOMEGALOVIRUS", "CMV", "EPSTEINBARR", "EBV")
    ):
        guard_reasons.append("respiratory herpesvirus reactivation prior")

    matching = [site for site in positive_direct_sites if target and site == target]
    mismatching = [site for site in positive_direct_sites if target and site != target]
    if hard_block:
        state = "D_BLOCK"
        aligned: TriState = False
        basis = "; ".join(guard_reasons)
    elif target and positive_direct_sites and not matching:
        # Direct structured sites outrank a generic upstream Aligned label.
        state = "D_BLOCK"
        aligned = False
        basis = f"structured direct evidence sites {positive_direct_sites} do not match target {target}"
    elif target and matching and mismatching:
        state = "D_GUARDED"
        aligned = None
        basis = f"structured direct evidence is site-discordant: match={matching}, mismatch={mismatching}"
        guard_reasons.append("mixed direct-evidence sites")
    elif target and matching:
        aligned = True
        basis = f"structured direct evidence site {matching[0]} matches target {target}"
        state = "D_GUARDED" if guard_reasons else "D_PASS"
    elif generic["aligned"] is False:
        state = "D_BLOCK"
        aligned = False
        basis = generic["basis"]
    elif generic["aligned"] is True:
        aligned = True
        basis = generic["basis"]
        state = "D_GUARDED" if guard_reasons else "D_PASS"
    else:
        state = "D_UNKNOWN"
        aligned = None
        basis = generic["basis"]

    return {
        "module": "D_SPECIMEN_COLONIZATION",
        "state": state,
        "status": (
            "unknown"
            if state == "D_UNKNOWN"
            else "guarded"
            if state == "D_GUARDED"
            else "ok"
        ),
        "positive": True if state == "D_PASS" else False if state == "D_BLOCK" else None,
        "site_aligned": aligned,
        "target_site": target,
        "mngs_specimen_site": generic["specimen_site"],
        "direct_evidence_sites": positive_direct_sites,
        "direct_site_priority_applied": bool(positive_direct_sites),
        "guard_reasons": guard_reasons,
        "reason": basis,
    }


def evaluate_b_strict(
    merged_context: Mapping[str, Any],
    frozen_candidate: Mapping[str, Any],
    config: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    """Module B_strict using only mNGS strength and mNGS/site fields."""

    cfg = config or {}
    candidate = _candidate_input(frozen_candidate)
    site = _mngs_site_alignment(merged_context, frozen_candidate)

    signal = _norm(candidate.get("mngs_signal_tier"))
    reads_tier = _norm(candidate.get("reads_tier"))
    dominance = _as_int(candidate.get("dominance_tier"))
    rank = _as_int(candidate.get("rank_priority"))
    percentile = _as_float(candidate.get("reads_percentile"))

    allowed_signals = {
        _norm(value)
        for value in cfg.get("positive_signal_tiers", ("M1_STRONG", "M2_MODERATE", "M1", "M2"))
    }
    allowed_reads = {
        _norm(value)
        for value in cfg.get("positive_reads_tiers", ("R3_HIGH", "R4_VERY_HIGH", "R3", "R4"))
    }
    absolute_inputs = [signal, reads_tier]
    absolute_present = any(value and value not in UNKNOWN_MARKERS for value in absolute_inputs)
    absolute: TriState = (
        signal in allowed_signals or reads_tier in allowed_reads if absolute_present else None
    )

    relative_present = any(value is not None for value in (dominance, rank, percentile))
    dominance_ok = dominance is not None and dominance >= int(cfg.get("min_dominance_tier", 2))
    rank_percentile_ok = (
        rank is not None
        and percentile is not None
        and rank <= int(cfg.get("max_rank_priority", 1))
        and percentile >= float(cfg.get("min_reads_percentile", 0.8))
    )
    relative: TriState = dominance_ok or rank_percentile_ok if relative_present else None

    ntc_ratio = _as_float(
        _first_present(
            candidate,
            ("sample_to_ntc_ratio", "sample_ntc_ratio", "ntc_ratio", "rk_ntc_ratio"),
        )
    )
    rpm = _as_float(
        _first_present(candidate, ("rpm", "reads_per_million", "organism_rpm", "sample_rpm"))
    )
    quality_gate: TriState = True
    quality_reasons: List[str] = []
    if cfg.get("require_ntc_ratio"):
        quality_gate = _tri_and(
            quality_gate,
            None if ntc_ratio is None else ntc_ratio >= float(cfg.get("min_ntc_ratio", 1.0)),
        )
        quality_reasons.append("NTC ratio required")
    if cfg.get("require_rpm"):
        quality_gate = _tri_and(
            quality_gate,
            None if rpm is None else rpm >= float(cfg.get("min_rpm", 0.0)),
        )
        quality_reasons.append("RPM required")

    positive = _tri_and(absolute, relative, site["aligned"], quality_gate)
    missing_features: List[str] = []
    if ntc_ratio is None:
        missing_features.append("ntc_control_ratio")
    if rpm is None:
        missing_features.append("rpm")
    core_missing: List[str] = []
    if absolute is None:
        core_missing.append("absolute_mngs_strength")
    if relative is None:
        core_missing.append("relative_mngs_strength")
    if site["aligned"] is None:
        core_missing.append("site_alignment")
    if core_missing:
        completeness = "missing"
        status = "missing"
    elif missing_features:
        completeness = "partial"
        status = "partial"
    else:
        completeness = "complete"
        status = "ok"

    return {
        "module": "B_STRICT_MNGS",
        "positive": positive,
        "status": status,
        "evidence_completeness": completeness,
        "missing_core_features": core_missing,
        "missing_quality_features": missing_features,
        "absolute_strength_axis": absolute,
        "relative_strength_axis": relative,
        "site_aligned": site["aligned"],
        "site_alignment_basis": site["basis"],
        "quality_gate": quality_gate,
        "reason": (
            f"absolute={absolute}, relative={relative}, site={site['aligned']}, quality={quality_gate}; "
            f"completeness={completeness}"
        ),
        "features": {
            "mngs_signal_tier": candidate.get("mngs_signal_tier"),
            "reads_tier": candidate.get("reads_tier"),
            "dominance_tier": candidate.get("dominance_tier"),
            "rank_priority": candidate.get("rank_priority"),
            "reads_percentile": candidate.get("reads_percentile"),
            "reads": candidate.get("reads"),
            "ntc_control_ratio": ntc_ratio,
            "rpm": rpm,
            "specimen_alignment": candidate.get("specimen_alignment"),
            "specimen_class": candidate.get("specimen_class"),
        },
        "quality_reasons": quality_reasons,
    }


def _unstructured_support_modules(candidate: Mapping[str, Any]) -> List[str]:
    found: List[str] = []
    key_evidence = _mapping(candidate.get("key_evidence"))
    for container in (candidate, key_evidence):
        modules = container.get("support_modules")
        if isinstance(modules, (list, tuple, set)):
            for module in modules:
                canonical = _module_from_key(module)
                if canonical and canonical in DIRECT_MODULES:
                    found.append(canonical)
    summary = _mapping(candidate.get("module_support_summary"))
    for module, value in summary.items():
        canonical = _module_from_key(module)
        normalized_value = _norm(value)
        if canonical and canonical in DIRECT_MODULES and normalized_value in POSITIVE_MARKERS:
            found.append(canonical)
    return sorted(set(found))


def grade_direct_evidence(
    merged_context: Mapping[str, Any],
    frozen_candidate: Mapping[str, Any],
    config: Mapping[str, Any] | None = None,
    d_result: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    """Module C: grade same-organism non-mNGS evidence from C0 to C3/CNEG."""

    cfg = config or {}
    candidate = _candidate_input(frozen_candidate)
    target = _target_site(merged_context, frozen_candidate, cfg)
    records = _structured_direct_records(merged_context, frozen_candidate, cfg)
    exact = [row for row in records if row["exact_organism"]]
    positive = [row for row in exact if row["status"] == "positive"]
    negative = [row for row in exact if row["status"] == "negative"]
    unknown = [row for row in exact if row["status"] in {"unknown", "discordant"}]
    aligned_positive = [row for row in positive if target and target in row["sites"]]
    mismatch_positive = [row for row in positive if row["sites"] and (not target or target not in row["sites"])]
    aligned_negative = [row for row in negative if target and target in row["sites"]]
    unstructured = _unstructured_support_modules(candidate)

    organism = _norm(
        frozen_candidate.get("canonical_organism_name")
        or frozen_candidate.get("organism_name")
        or candidate.get("organism_name")
    )
    respiratory_candida = target == "lower_respiratory" and (
        "CANDIDA" in organism or organism in {"YEAST", "GENERICYEAST"}
    )
    positive_modules = {row["module"] for row in positive}
    candida_culture_only = respiratory_candida and bool(
        positive_modules or unstructured
    ) and (positive_modules | set(unstructured)).issubset({"culture"})

    discordant = bool(aligned_positive and aligned_negative) or any(
        row["status"] == "discordant" for row in exact
    )
    manual_review = False
    positive_value: TriState
    evidence_state: str
    if discordant:
        grade = "C1"
        positive_value = None
        evidence_state = "discordant"
        manual_review = True
        reason = "Aligned structured positive and negative evidence conflict."
    elif candida_culture_only:
        grade = "C1"
        positive_value = None
        evidence_state = "context_only"
        manual_review = True
        reason = "Respiratory Candida/yeast supported only by culture is C1 and requires manual review."
    else:
        qualifying = [
            row
            for row in aligned_positive
            if row.get("contaminated") is not True
            and (
                row["module"] in HIGH_SPECIFICITY_MODULES
                or (
                    row["module"] == "filmarray_gmtest"
                    and row.get("organism_explicit") is True
                )
            )
        ]
        sterile_qualifying = [row for row in qualifying if row["sterile"]]
        if sterile_qualifying:
            grade = "C3"
            positive_value = True
            evidence_state = "positive"
            reason = "Exact-organism high-specificity structured evidence is site-aligned and sterile/tissue based."
        elif qualifying:
            grade = "C2"
            positive_value = True
            evidence_state = "positive"
            reason = "Exact-organism high-specificity structured evidence has an explicit matching site."
        elif positive or unstructured:
            grade = "C1"
            positive_value = None
            evidence_state = "context_only"
            manual_review = True
            if mismatch_positive:
                reason = "Structured direct evidence is at a non-target site; generic Aligned cannot override it."
            elif any(row.get("contaminated") is True for row in positive):
                reason = "Structured positive evidence is explicitly flagged as contaminant or poor quality."
            elif positive:
                reason = "Structured positive evidence lacks a usable matching site or exact high-specificity assay."
            else:
                reason = "Only module-level/support-summary context is available; it cannot establish C2/C3."
        elif aligned_negative:
            grade = "CNEG"
            positive_value = False
            evidence_state = "explicit_negative"
            reason = "An exact-organism structured test at the target site is explicitly negative."
        else:
            grade = "C0"
            positive_value = None
            evidence_state = "unknown" if unknown else "no_direct_evidence"
            reason = (
                "Direct-test state is pending/not done/missing/unknown; this is not negative."
                if unknown
                else "No exact-organism structured direct evidence was found."
            )

    return {
        "module": "C_DIRECT_EVIDENCE_GRADE",
        "grade": grade,
        "positive": positive_value,
        "status": evidence_state,
        "manual_review": manual_review,
        "target_site": target,
        "structured_record_count": len(records),
        "exact_record_count": len(exact),
        "positive_record_count": len(positive),
        "negative_record_count": len(negative),
        "unknown_record_count": len(unknown),
        "aligned_positive_count": len(aligned_positive),
        "mismatched_positive_count": len(mismatch_positive),
        "support_modules_context_only": unstructured,
        "related_representative_is_exact": False,
        "respiratory_candida_culture_only": candida_culture_only,
        "reason": reason,
        "records": records,
    }


DATE_KEYS_INDEX = {
    "INDEX_DATE",
    "SYNDROME_ONSET",
    "INFECTION_ONSET",
    "CLINICAL_ONSET",
    "SYMPTOM_ONSET",
    "ADMISSION_DATE",
}
DATE_KEYS_SPECIMEN = {
    "SPECIMEN_DATE",
    "COLLECTION_DATE",
    "COLLECTED_AT",
    "SAMPLED_AT",
    "TEST_DATE",
    "SAMPLING_DATE",
}
DATE_KEYS_TREATMENT = {
    "ANTIMICROBIAL_START",
    "ANTIBIOTIC_START",
    "ANTIVIRAL_START",
    "ANTIFUNGAL_START",
    "TREATMENT_START",
}


def _find_values_by_keys(value: Any, keys: set[str]) -> List[Any]:
    found: List[Any] = []
    if isinstance(value, Mapping):
        for key, nested in value.items():
            if _norm(key) in keys and nested not in (None, ""):
                found.append(nested)
            if isinstance(nested, (Mapping, list, tuple)):
                found.extend(_find_values_by_keys(nested, keys))
    elif isinstance(value, (list, tuple)):
        for nested in value:
            found.extend(_find_values_by_keys(nested, keys))
    return found


def _parse_date(value: Any) -> date | None:
    if isinstance(value, datetime):
        return value.date()
    if isinstance(value, date):
        return value
    text = str(value or "").strip()
    if not text or _norm(text) in UNKNOWN_MARKERS:
        return None
    text = text.replace("Z", "+00:00")
    try:
        return datetime.fromisoformat(text).date()
    except ValueError:
        pass
    for pattern in ("%Y/%m/%d", "%Y-%m-%d", "%d/%m/%Y", "%Y%m%d"):
        try:
            return datetime.strptime(text[:10], pattern).date()
        except ValueError:
            continue
    return None


def _first_date(values: Iterable[Any]) -> date | None:
    dates = [_parse_date(value) for value in values]
    usable = [value for value in dates if value is not None]
    return min(usable) if usable else None


def evaluate_temporal_coherence(
    merged_context: Mapping[str, Any],
    frozen_candidate: Mapping[str, Any],
    config: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    """Module F: compare candidate specimen timing with the clinical index date."""

    cfg = config or {}
    sources = _matching_candidate_sources(merged_context, frozen_candidate, cfg)
    index_date = _first_date(_find_values_by_keys(merged_context, DATE_KEYS_INDEX))
    specimen_date = _first_date(
        value for source in sources for value in _find_values_by_keys(source, DATE_KEYS_SPECIMEN)
    ) or _first_date(_find_values_by_keys(merged_context, DATE_KEYS_SPECIMEN))
    treatment_date = _first_date(_find_values_by_keys(merged_context, DATE_KEYS_TREATMENT))

    before = int(cfg.get("window_before_days", 3))
    after = int(cfg.get("window_after_days", 3))
    delta: int | None = None
    if index_date and specimen_date:
        delta = (specimen_date - index_date).days
        coherent = -before <= delta <= after
        state = "F_COHERENT" if coherent else "F_MISMATCH"
        reason = (
            f"specimen is {delta:+d} day(s) from the clinical index; allowed window is "
            f"-{before} to +{after} days"
        )
    else:
        coherent = None
        state = "F_UNKNOWN"
        missing = []
        if not index_date:
            missing.append("clinical_index_date")
        if not specimen_date:
            missing.append("candidate_specimen_date")
        reason = f"Missing temporal fields: {', '.join(missing)}; unknown is not negative."

    pre_antimicrobial: TriState = (
        specimen_date <= treatment_date if specimen_date and treatment_date else None
    )
    return {
        "module": "F_TEMPORAL_COHERENCE",
        "state": state,
        "positive": coherent,
        "status": "ok" if coherent is not None else "unknown",
        "index_date": index_date.isoformat() if index_date else None,
        "specimen_date": specimen_date.isoformat() if specimen_date else None,
        "treatment_start_date": treatment_date.isoformat() if treatment_date else None,
        "delta_days": delta,
        "window_before_days": before,
        "window_after_days": after,
        "pre_antimicrobial": pre_antimicrobial,
        "post_treatment_specimen": pre_antimicrobial is False,
        "reason": reason,
    }


def evaluate_convergence(
    d_result: Mapping[str, Any],
    b_result: Mapping[str, Any],
    c_result: Mapping[str, Any],
    a_positive: TriState = None,
    f_result: Mapping[str, Any] | None = None,
    config: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    """Module E: conservative convergence policy; never mutates E0."""

    cfg = config or {}
    d_state = d_result.get("state")
    b = b_result.get("positive")
    c_grade = c_result.get("grade")
    f_state = (f_result or {}).get("state")

    auto = False
    review = False
    if d_state == "D_BLOCK":
        tier = "E_BLOCKED"
        disposition = "do_not_rescue"
        reason = "Specimen/colonization gate blocks automatic rescue."
    elif c_grade == "CNEG":
        tier = "E_DISCORDANT_NEGATIVE"
        disposition = "manual_review" if b is True or a_positive is True else "do_not_rescue"
        review = disposition == "manual_review"
        reason = "Explicit target-site direct negative evidence prevents automatic rescue."
    elif d_state == "D_UNKNOWN":
        tier = "E_UNKNOWN_SITE"
        disposition = "manual_review" if b is True or c_grade in {"C1", "C2", "C3"} else "do_not_rescue"
        review = disposition == "manual_review"
        reason = "Required site semantics are unknown."
    elif d_state == "D_GUARDED" and c_grade != "C3":
        tier = "E_GUARDED_REVIEW"
        disposition = "manual_review" if b is True or c_grade in {"C1", "C2"} or a_positive is True else "do_not_rescue"
        review = disposition == "manual_review"
        reason = "A colonization/reactivation guardrail requires manual review unless C3 is present."
    elif c_grade == "C3" and d_state in {"D_PASS", "D_GUARDED"}:
        tier = "E3_DIRECT_HIGH"
        disposition = "auto_rescue"
        auto = True
        reason = "Site-coherent C3 direct evidence satisfies the highest convergence tier."
    elif d_state == "D_PASS" and c_grade == "C2" and b is True:
        tier = "E2_CONVERGENT_HIGH"
        disposition = "auto_rescue"
        auto = True
        reason = "Strict mNGS and exact site-aligned direct evidence converge."
    elif d_state == "D_PASS" and b is True and a_positive is True:
        tier = "E1_BIOLOGIC_PLUS_SIGNAL"
        auto = bool(cfg.get("allow_literature_mngs_auto", False))
        disposition = "auto_rescue" if auto else "manual_review"
        review = not auto
        reason = "Literature plus strict mNGS is supportive, but literature is not patient-specific."
    elif any((b is True, c_grade in {"C1", "C2"}, a_positive is True)):
        tier = "E_SINGLE_AXIS_REVIEW"
        disposition = "manual_review"
        review = True
        reason = "Evidence is present without a pre-specified automatic convergence pattern."
    else:
        tier = "E_NO_CONVERGENCE"
        disposition = "do_not_rescue"
        reason = "No pre-specified clinical convergence pattern was met."

    if auto and f_state == "F_MISMATCH":
        auto = False
        review = True
        tier = "E_TEMPORAL_MISMATCH_REVIEW"
        disposition = "manual_review"
        reason = "The evidence converges, but its timing is outside the pre-specified clinical window."
    elif auto and f_state == "F_UNKNOWN" and cfg.get("require_temporal_coherence_for_auto", False):
        auto = False
        review = True
        tier = "E_TEMPORAL_UNKNOWN_REVIEW"
        disposition = "manual_review"
        reason = "Temporal coherence is required for this arm but is unknown."

    return {
        "module": "E_CONVERGENCE",
        "state": tier,
        "tier": tier,
        "disposition": disposition,
        "auto_rescue": auto,
        "review_required": review,
        "e0_immutable": True,
        "inputs": {
            "D": d_state,
            "B_strict": b,
            "C": c_grade,
            "A": a_positive,
            "F": f_state,
        },
        "reason": reason,
    }


def _config_section(config: Mapping[str, Any], *names: str) -> Mapping[str, Any]:
    for name in names:
        value = config.get(name)
        if isinstance(value, Mapping):
            return value
    return {}


def evaluate_modules(
    merged_context: Mapping[str, Any],
    frozen_candidate: Mapping[str, Any],
    config: Mapping[str, Any] | None = None,
) -> Dict[str, Any]:
    """Evaluate D, B_strict, C, F and E without changing the frozen candidate."""

    raw_config = config or {}
    cfg = _mapping(raw_config.get("clinical_rules")) or raw_config
    d_cfg = _config_section(cfg, "d_specimen_colonization", "module_d")
    b_cfg = _config_section(cfg, "b_strict", "module_b_strict")
    c_cfg = _config_section(cfg, "c_evidence", "module_c")
    f_cfg = _config_section(cfg, "f_temporal", "module_f")
    e_cfg = _config_section(cfg, "e_convergence", "decision")

    d_result = specimen_colonization_gate(merged_context, frozen_candidate, d_cfg)
    b_result = evaluate_b_strict(merged_context, frozen_candidate, b_cfg)
    c_result = grade_direct_evidence(
        merged_context, frozen_candidate, c_cfg, d_result=d_result
    )
    f_result = evaluate_temporal_coherence(merged_context, frozen_candidate, f_cfg)
    modules = _mapping(frozen_candidate.get("modules"))
    a_module = _mapping(modules.get("A"))
    a_positive = a_module.get("positive") if a_module.get("positive") in (True, False) else None
    e_result = evaluate_convergence(
        d_result,
        b_result,
        c_result,
        a_positive=a_positive,
        f_result=f_result,
        config=e_cfg,
    )
    return {
        "D": d_result,
        "B_STRICT": b_result,
        "C": c_result,
        "E": e_result,
        "F": f_result,
        "contract": {
            "pure": True,
            "e0_immutable": True,
            "unknown_is_negative": False,
        },
    }

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from tools.llm_requester import load_prompt_template, send_to_llm

DEFAULT_SOURCE = Path("outputs") / "mngs_candidate_microbes.json"
DEFAULT_METADATA = Path("outputs") / "specimen_lookup_summary_mNGS_grouped_with_type_filtered.json"
DEFAULT_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs.json"
PROMPT_TASK = "mngs_to_specimen_agent"
NAME_SANITIZER = re.compile(r"[^a-z0-9]+")

MATCH_ALIASES = {
    "cmv": "humanbetaherpesvirus5",
    "humanbetaherpesvirus5": "humanbetaherpesvirus5",
    "ebv": "humangammaherpesvirus4",
    "humangammaherpesvirus4": "humangammaherpesvirus4",
    "hhv6b": "humanbetaherpesvirus6b",
    "humanbetaherpesvirus6b": "humanbetaherpesvirus6b",
    "hhv7": "humanbetaherpesvirus7",
    "humanbetaherpesvirus7": "humanbetaherpesvirus7",
    "hsv1": "humanalphaherpesvirus1",
    "humanalphaherpesvirus1": "humanalphaherpesvirus1",
    "hsv2": "humanalphaherpesvirus2",
    "humanalphaherpesvirus2": "humanalphaherpesvirus2",
    "metamycoplasmaorale": "mycoplasmaorale",
    "mycoplasmaorale": "mycoplasmaorale",
    "metamycoplasmasalivarium": "mycoplasmasalivarium",
    "mycoplasmasalivarium": "mycoplasmasalivarium",
    "metamycoplasmahominis": "mycoplasmahominis",
    "mycoplasmahominis": "mycoplasmahominis",
    "nakaseomycesglabratus": "candidaglabrata",
    "candidaglabrata": "candidaglabrata",
    "candidaduobushaemulonis": "candidahaemuloniicomplex",
    "candidahaemuloniicomplex": "candidahaemuloniicomplex",
    "legionellapneumophilasubsppneumophila": "legionellapneumophila",
    "legionellapneumophila": "legionellapneumophila",
}

DISPLAY_NAME_MAP = {
    "Metamycoplasma orale": "Mycoplasma orale",
    "Metamycoplasma salivarium": "Mycoplasma salivarium",
    "Metamycoplasma hominis": "Mycoplasma hominis",
    "Nakaseomyces glabratus": "Candida glabrata",
    "Human betaherpesvirus 5": "CMV",
    "Human gammaherpesvirus 4": "EBV",
    "Human betaherpesvirus 7": "HHV-7",
    "Legionella pneumophila subsp. pneumophila": "Legionella pneumophila",
    "[Candida] duobushaemulonis": "Candida haemulonii complex",
}

FORCE_KEEP_BACTERIA = {
    "enterococcusfaecium",
    "acinetobacterbaumannii",
    "corynebacteriumstriatum",
    "staphylococcushaemolyticus",
    "pseudomonasaeruginosa",
    "klebsiellapneumoniae",
    "mycoplasmasalivarium",
    "legionellapneumophila",
}
FORCE_KEEP_VIRUSES = {
    "humanbetaherpesvirus5",
    "humangammaherpesvirus4",
    "humanbetaherpesvirus6b",
    "humanbetaherpesvirus7",
    "humanpapillomavirus11",
}
PNEUMOCYSTIS = "pneumocystisjirovecii"
ORAL_ANAEROBE_PRIORITY = {
    "parvimonasmicra": 0,
    "prevotellanigrescens": 1,
    "prevotellahisticola": 1,
    "prevotellamelaninogenica": 2,
    "porphyromonasendodontalis": 3,
    "fusobacteriumnucleatum": 4,
}
ORAL_COLONIZER_GENERA = ("parvimonas", "prevotella", "porphyromonas", "fusobacterium")
CATEGORY_TO_GROUP = {"1.Bac": "bacterial", "2.Fungi": "fungal", "3.Virus": "viral"}
GROUP_KEYS = ("bacterial", "viral", "fungal", "others")
LIKELIHOOD_LEVEL_KEYS = ("high", "medium", "low_colonizer")
HIGH_READ_THRESHOLDS = {"bacterial": 120, "viral": 80, "fungal": 80, "others": 120}
MEDIUM_READ_THRESHOLDS = {"bacterial": 20, "viral": 20, "fungal": 20, "others": 30}
BACKFILL_MEDIUM_LIMIT = {"bacterial": 4, "viral": 3, "fungal": 3, "others": 2}
BACKFILL_LOW_LIMIT = {"bacterial": 3, "viral": 2, "fungal": 2, "others": 2}
BALF_ORAL_MAX_KEEP = 3
AGGRESSIVE_ALL_LEVELS_RECALL = True
LOW_COLONIZER_MIN_READS = 5
LOW_COLONIZER_BACKGROUND_PREFIXES = (
    "anelloviridae",
    "torquetenovirus",
    "betatorquevirus",
    "ttv",
)
LOW_COLONIZER_LOW_READ_ALLOW_PREFIXES = (
    "humanbetaherpesvirus5",
    "humangammaherpesvirus4",
    "humanbetaherpesvirus6b",
    "humanbetaherpesvirus7",
    "humanalphaherpesvirus",
    "humanpapillomavirus",
    "pneumocystisjirovecii",
    "candida",
    "enterobacterales",
    "staphylococcus",
    "pseudomonas",
    "xanthomonadaceae",
    "chryseobacterium",
    "enterococcusfaecium",
    "escherichiacoli",
)
LOW_COLONIZER_EXACT_EXCLUDE = {
    # High-frequency output-only species in this project that never matched Excel labels.
    "acinetobacterbereziniae",
    "acinetobacterjohnsonii",
    "anelloviridaesp",
    "prevotellajejuni",
    "rothiamucilaginosa",
    "scardoviawiggsiae",
    "streptococcusparasanguinis",
    "talaromycespinophilus",
    "undibacteriumoligocarboniphilum",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the mNGS-to-SPECIMEN prompt on one patient or all patients."
    )
    parser.add_argument(
        "--source-json",
        type=Path,
        default=DEFAULT_SOURCE,
        help=f"Path to mngs candidate JSON (default: {DEFAULT_SOURCE})",
    )
    parser.add_argument(
        "--metadata-json",
        type=Path,
        default=DEFAULT_METADATA,
        help=(
            "Optional specimen-style metadata JSON used to fill specimen_code/collected_time/"
            f"specimen_site (default: {DEFAULT_METADATA})"
        ),
    )
    parser.add_argument(
        "--patient-id",
        help="Run only one patient id. If omitted, run all patients.",
    )
    parser.add_argument(
        "--model",
        default="gpt-5",
        help="OpenAI model name (default: gpt-5)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output path for aggregated JSON results (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the rendered prompt payload for the target patient(s) without calling the API.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict payload from {path}, got {type(payload)}")
    return payload


def get_metadata_entry(metadata_payload: dict[str, Any], patient_id: str) -> dict[str, Any]:
    entries = metadata_payload.get(patient_id)
    if isinstance(entries, list) and entries:
        entry = entries[0]
        if isinstance(entry, dict):
            return entry
    return {}


def get_metadata_entries(metadata_payload: dict[str, Any], patient_id: str) -> list[dict[str, Any]]:
    entries = metadata_payload.get(patient_id)
    if isinstance(entries, list):
        return [entry for entry in entries if isinstance(entry, dict)]
    return []


def normalize_organism_name(name: Any) -> str:
    text = str(name or "").strip().lower()
    if not text:
        return ""
    normalized = NAME_SANITIZER.sub("", text)
    return MATCH_ALIASES.get(normalized, normalized)


def canonical_display_name(name: str) -> str:
    return DISPLAY_NAME_MAP.get(str(name or "").strip(), str(name or "").strip())


def should_force_keep_candida(normalized: str, sec_hit: Any) -> bool:
    numeric = sec_hit or 0
    if normalized == "candidaalbicans":
        return True
    if normalized == "candidahaemuloniicomplex":
        return numeric >= 10
    if normalized in {"candidatropicalis", "candidaglabrata"}:
        return numeric >= 50
    return False


def should_force_keep(candidate: dict[str, Any]) -> bool:
    normalized = normalize_organism_name(candidate.get("organism_name"))
    category = str(candidate.get("source_category", ""))
    sec_hit = candidate.get("sec_hit", 0) or 0
    if not normalized:
        return False
    if category == "2.Fungi":
        return normalized == PNEUMOCYSTIS or should_force_keep_candida(normalized, sec_hit)
    if category == "3.Virus":
        return normalized in FORCE_KEEP_VIRUSES
    if category == "1.Bac":
        if normalized == "stenotrophomonasmaltophilia":
            return sec_hit >= 50
        return normalized in FORCE_KEEP_BACTERIA
    return False


def _to_float(value: Any) -> float:
    if value is None:
        return 0.0
    if isinstance(value, bool):
        return float(int(value))
    try:
        return float(value)
    except (TypeError, ValueError):
        text = str(value).strip()
        if not text:
            return 0.0
        try:
            return float(text)
        except ValueError:
            return 0.0


def _to_output_reads(value: Any) -> int | float:
    numeric = _to_float(value)
    return int(numeric) if float(numeric).is_integer() else numeric


def _is_oral_colonizer_name(normalized: str) -> bool:
    if not normalized:
        return False
    if normalized in ORAL_ANAEROBE_PRIORITY:
        return True
    return any(normalized.startswith(genus) for genus in ORAL_COLONIZER_GENERA)


def _is_force_keep_name(group_name: str, normalized: str, reads: Any) -> bool:
    numeric_reads = _to_float(reads)
    if group_name == "fungal":
        return normalized == PNEUMOCYSTIS or should_force_keep_candida(normalized, numeric_reads)
    if group_name == "viral":
        return normalized in FORCE_KEEP_VIRUSES
    if group_name == "bacterial":
        if normalized == "stenotrophomonasmaltophilia":
            return numeric_reads >= 50
        return normalized in FORCE_KEEP_BACTERIA
    return False


def _is_high_value_non_oral_bacteria(normalized: str) -> bool:
    return normalized in {
        "pseudomonasaeruginosa",
        "klebsiellapneumoniae",
        "enterococcusfaecium",
        "acinetobacterbaumannii",
    }


def _build_related_support_set(related_specimens: list[dict[str, Any]]) -> set[str]:
    supported: set[str] = set()
    for specimen in related_specimens:
        if not isinstance(specimen, dict):
            continue
        for candidate in specimen.get("related_candidates", []):
            if not isinstance(candidate, dict):
                continue
            normalized = normalize_organism_name(candidate.get("organism_name"))
            if normalized:
                supported.add(normalized)
    return supported


def _normalize_group_pathogens(raw_group: Any) -> dict[str, list[dict[str, Any]]]:
    normalized = {key: [] for key in GROUP_KEYS}
    if not isinstance(raw_group, dict):
        return normalized
    for group_name in GROUP_KEYS:
        values = raw_group.get(group_name)
        if isinstance(values, list):
            normalized[group_name] = [item for item in values if isinstance(item, dict)]
    return normalized


def _add_unique_item(
    target: dict[str, list[dict[str, Any]]],
    seen: set[tuple[str, str]],
    group_name: str,
    item: dict[str, Any],
) -> bool:
    normalized = normalize_organism_name(item.get("name"))
    if not normalized:
        return False
    key = (group_name, normalized)
    if key in seen:
        return False
    target[group_name].append(
        {
            "name": str(item.get("name", "")).strip(),
            "reads": _to_output_reads(item.get("reads", 0)),
        }
    )
    seen.add(key)
    return True


def _sort_pathogen_lists(pathogens: dict[str, list[dict[str, Any]]]) -> None:
    for group_name in GROUP_KEYS:
        pathogens[group_name] = sorted(
            [item for item in pathogens[group_name] if isinstance(item, dict)],
            key=lambda item: (-_to_float(item.get("reads", 0)), str(item.get("name", ""))),
        )


def _is_low_colonizer_low_read_allowed(normalized: str) -> bool:
    return any(normalized.startswith(prefix) for prefix in LOW_COLONIZER_LOW_READ_ALLOW_PREFIXES)


def _prune_low_colonizer_noise(tiers: dict[str, dict[str, list[dict[str, Any]]]]) -> None:
    low_group = tiers.get("low_colonizer")
    if not isinstance(low_group, dict):
        return
    for group_name in GROUP_KEYS:
        values = low_group.get(group_name)
        if not isinstance(values, list):
            continue
        kept: list[dict[str, Any]] = []
        for item in values:
            if not isinstance(item, dict):
                continue
            normalized = normalize_organism_name(item.get("name"))
            if not normalized:
                continue
            if normalized in LOW_COLONIZER_EXACT_EXCLUDE:
                continue
            if any(normalized.startswith(prefix) for prefix in LOW_COLONIZER_BACKGROUND_PREFIXES):
                continue
            reads = _to_float(item.get("reads", 0))
            if reads < LOW_COLONIZER_MIN_READS and not _is_low_colonizer_low_read_allowed(normalized):
                continue
            kept.append(item)
        low_group[group_name] = kept


def _choose_likelihood_tier(
    *,
    group_name: str,
    normalized: str,
    reads: Any,
    specimen_site: str,
    related_supported: bool,
    has_strong_non_oral_bacteria: bool,
) -> str:
    read_value = _to_float(reads)
    site = str(specimen_site or "").strip().lower()

    if _is_force_keep_name(group_name, normalized, read_value):
        return "high"
    if related_supported and read_value >= 5:
        return "high"
    if read_value >= HIGH_READ_THRESHOLDS.get(group_name, 120):
        return "high"

    if group_name == "bacterial" and site == "balf" and _is_oral_colonizer_name(normalized):
        if has_strong_non_oral_bacteria and read_value < 120:
            return "low_colonizer"
        if read_value >= 120:
            return "medium"
        return "low_colonizer"

    if related_supported or read_value >= MEDIUM_READ_THRESHOLDS.get(group_name, 20):
        return "medium"

    if group_name == "bacterial" and site == "balf":
        return "low_colonizer"
    return "medium"


def _should_backfill_candidate(
    *,
    candidate: dict[str, Any],
    group_name: str,
    specimen_site: str,
    related_supported: bool,
    has_strong_non_oral_bacteria: bool,
) -> str | None:
    normalized = normalize_organism_name(candidate.get("organism_name"))
    if not normalized:
        return None
    reads = _to_float(candidate.get("sec_hit", 0))
    site = str(specimen_site or "").strip().lower()

    if _is_force_keep_name(group_name, normalized, reads):
        return "high"
    if group_name == "bacterial" and site == "balf" and _is_oral_colonizer_name(normalized):
        if has_strong_non_oral_bacteria and reads < 120:
            return "low_colonizer"
        return "medium" if reads >= MEDIUM_READ_THRESHOLDS["bacterial"] else "low_colonizer"
    if related_supported and reads >= 5:
        return "medium"
    if reads >= MEDIUM_READ_THRESHOLDS.get(group_name, 20):
        return "medium"
    return None


def prune_balf_oral_anaerobes(
    curated_output: dict[str, Any],
    specimen_site: str,
    current_candidates: list[dict[str, Any]],
) -> dict[str, Any]:
    if str(specimen_site or "").strip().lower() != "balf":
        return curated_output

    pathogens = curated_output.get("pathogens", {})
    bacteria = pathogens.get("bacterial")
    if not isinstance(bacteria, list):
        return curated_output

    candidate_norms = {
        normalize_organism_name(candidate.get("organism_name"))
        for candidate in current_candidates
        if isinstance(candidate, dict)
    }
    has_strong_non_anaerobe = any(
        name in candidate_norms
        for name in {
            "pseudomonasaeruginosa",
            "klebsiellapneumoniae",
            "enterococcusfaecium",
            "acinetobacterbaumannii",
        }
    )
    if not has_strong_non_anaerobe:
        return curated_output

    oral_items: list[dict[str, Any]] = []
    kept_items: list[dict[str, Any]] = []
    for item in bacteria:
        normalized = normalize_organism_name(item.get("name"))
        if normalized in ORAL_ANAEROBE_PRIORITY:
            oral_items.append(item)
        else:
            kept_items.append(item)

    if len(oral_items) <= BALF_ORAL_MAX_KEEP:
        return curated_output

    oral_items = sorted(
        oral_items,
        key=lambda item: (
            ORAL_ANAEROBE_PRIORITY.get(normalize_organism_name(item.get("name")), 99),
            -_to_float(item.get("reads", 0)),
            str(item.get("name", "")),
        ),
    )
    pathogens["bacterial"] = kept_items + oral_items[:BALF_ORAL_MAX_KEEP]
    return curated_output


def augment_curated_output(
    curated_output: dict[str, Any],
    current_candidates: list[dict[str, Any]],
    related_specimens: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    pathogens = _normalize_group_pathogens(curated_output.get("pathogens"))
    curated_output["pathogens"] = pathogens

    existing: set[tuple[str, str]] = set()
    for group_name in GROUP_KEYS:
        for item in pathogens[group_name]:
            normalized = normalize_organism_name(item.get("name"))
            if normalized:
                existing.add((group_name, normalized))

    for candidate in current_candidates:
        if not isinstance(candidate, dict) or not should_force_keep(candidate):
            continue
        group_name = CATEGORY_TO_GROUP.get(str(candidate.get("source_category", "")))
        normalized = normalize_organism_name(candidate.get("organism_name"))
        if not group_name or not normalized or (group_name, normalized) in existing:
            continue
        pathogens[group_name].append(
            {
                "name": canonical_display_name(str(candidate.get("organism_name", ""))),
                "reads": _to_output_reads(candidate.get("sec_hit", 0)),
            }
        )
        existing.add((group_name, normalized))

    _sort_pathogen_lists(pathogens)
    prune_balf_oral_anaerobes(
        curated_output,
        curated_output.get("specimen_site", ""),
        current_candidates,
    )
    pathogens = _normalize_group_pathogens(curated_output.get("pathogens"))
    curated_output["pathogens"] = pathogens

    related_support_set = _build_related_support_set(related_specimens or [])
    bacterial_candidate_names = {
        normalize_organism_name(candidate.get("organism_name"))
        for candidate in current_candidates
        if isinstance(candidate, dict) and CATEGORY_TO_GROUP.get(str(candidate.get("source_category", ""))) == "bacterial"
    }
    has_strong_non_oral_bacteria = any(
        _is_high_value_non_oral_bacteria(name) for name in bacterial_candidate_names if name
    )

    tiers: dict[str, dict[str, list[dict[str, Any]]]] = {
        tier_name: {group_name: [] for group_name in GROUP_KEYS}
        for tier_name in LIKELIHOOD_LEVEL_KEYS
    }
    tier_seen: set[tuple[str, str]] = set()
    specimen_site = str(curated_output.get("specimen_site", ""))

    for group_name in GROUP_KEYS:
        for item in pathogens[group_name]:
            normalized = normalize_organism_name(item.get("name"))
            if not normalized:
                continue
            tier_name = _choose_likelihood_tier(
                group_name=group_name,
                normalized=normalized,
                reads=item.get("reads", 0),
                specimen_site=specimen_site,
                related_supported=normalized in related_support_set,
                has_strong_non_oral_bacteria=has_strong_non_oral_bacteria,
            )
            _add_unique_item(tiers[tier_name], tier_seen, group_name, item)

    added_counts: dict[str, dict[str, int]] = {
        tier_name: {group_name: 0 for group_name in GROUP_KEYS}
        for tier_name in LIKELIHOOD_LEVEL_KEYS
    }
    sorted_candidates = sorted(
        [candidate for candidate in current_candidates if isinstance(candidate, dict)],
        key=lambda item: (
            str(item.get("source_category", "")),
            -_to_float(item.get("sec_hit", 0)),
            str(item.get("organism_name", "")),
        ),
    )
    for candidate in sorted_candidates:
        group_name = CATEGORY_TO_GROUP.get(str(candidate.get("source_category", "")))
        normalized = normalize_organism_name(candidate.get("organism_name"))
        if not group_name or not normalized or (group_name, normalized) in tier_seen:
            continue

        tier_name = _should_backfill_candidate(
            candidate=candidate,
            group_name=group_name,
            specimen_site=specimen_site,
            related_supported=normalized in related_support_set,
            has_strong_non_oral_bacteria=has_strong_non_oral_bacteria,
        )
        if tier_name is None:
            continue
        if tier_name == "medium" and added_counts[tier_name][group_name] >= BACKFILL_MEDIUM_LIMIT[group_name]:
            continue
        if tier_name == "low_colonizer" and added_counts[tier_name][group_name] >= BACKFILL_LOW_LIMIT[group_name]:
            continue

        item = {
            "name": canonical_display_name(str(candidate.get("organism_name", ""))),
            "reads": _to_output_reads(candidate.get("sec_hit", 0)),
        }
        if _add_unique_item(tiers[tier_name], tier_seen, group_name, item):
            added_counts[tier_name][group_name] += 1

    if AGGRESSIVE_ALL_LEVELS_RECALL:
        for candidate in sorted_candidates:
            group_name = CATEGORY_TO_GROUP.get(str(candidate.get("source_category", "")))
            normalized = normalize_organism_name(candidate.get("organism_name"))
            if not group_name or not normalized or (group_name, normalized) in tier_seen:
                continue
            item = {
                "name": canonical_display_name(str(candidate.get("organism_name", ""))),
                "reads": _to_output_reads(candidate.get("sec_hit", 0)),
            }
            _add_unique_item(tiers["low_colonizer"], tier_seen, group_name, item)

    _prune_low_colonizer_noise(tiers)
    for tier_name in LIKELIHOOD_LEVEL_KEYS:
        _sort_pathogen_lists(tiers[tier_name])

    curated_output["pathogens_by_likelihood"] = tiers

    primary_pathogens = {group_name: [] for group_name in GROUP_KEYS}
    primary_seen: set[tuple[str, str]] = set()
    for tier_name in ("high", "medium"):
        for group_name in GROUP_KEYS:
            for item in tiers[tier_name][group_name]:
                _add_unique_item(primary_pathogens, primary_seen, group_name, item)
    _sort_pathogen_lists(primary_pathogens)
    curated_output["pathogens"] = primary_pathogens
    return curated_output


def build_candidates(
    patient_payload: dict[str, Any],
    specimen_code: str | None = None,
) -> list[dict[str, Any]]:
    deduped: dict[tuple[str, str], dict[str, Any]] = {}
    for record in patient_payload.get("records", []):
        if not isinstance(record, dict):
            continue
        record_specimen_code = str(record.get("specimen_code", ""))
        if specimen_code and record_specimen_code != specimen_code:
            continue
        for candidate in record.get("candidates", []):
            if not isinstance(candidate, dict):
                continue
            organism_name = str(candidate.get("organism_name", ""))
            source_category = str(candidate.get("source_category", ""))
            sec_hit = candidate.get("sec_hit", 0)
            key = (organism_name, source_category)
            current = deduped.get(key)
            if current is None or (sec_hit or 0) > (current.get("sec_hit", 0) or 0):
                deduped[key] = {
                    "organism_name": organism_name,
                    "sec_hit": sec_hit,
                    "source_category": source_category,
                }
    return sorted(
        deduped.values(),
        key=lambda item: (str(item.get("source_category", "")), -(item.get("sec_hit", 0) or 0), str(item.get("organism_name", ""))),
    )


def derive_specimen_entries_from_source(patient_payload: dict[str, Any]) -> list[dict[str, Any]]:
    seen: set[str] = set()
    entries: list[dict[str, Any]] = []
    for record in patient_payload.get("records", []):
        if not isinstance(record, dict):
            continue
        specimen_code = str(record.get("specimen_code", "")).strip()
        if not specimen_code or specimen_code in seen:
            continue
        seen.add(specimen_code)
        entries.append(
            {
                "specimen_code": specimen_code,
                "collected_time": "",
                "specimen_site": "",
            }
        )
    return entries


def extract_related_specimens(
    patient_payload: dict[str, Any],
    current_specimen_code: str,
    metadata_entries: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    related: list[dict[str, Any]] = []
    for entry in metadata_entries or derive_specimen_entries_from_source(patient_payload):
        specimen_code = str(entry.get("specimen_code", ""))
        if specimen_code == current_specimen_code:
            continue
        if not specimen_code:
            continue
        related.append(
            {
                "specimen_code": specimen_code,
                "collected_time": entry.get("collected_time", ""),
                "specimen_site": entry.get("specimen_site", ""),
                "related_candidates": build_candidates(patient_payload, specimen_code=specimen_code),
            }
        )
    return related


def build_specimen_payload(
    patient_id: str,
    patient_payload: dict[str, Any],
    metadata_entry: dict[str, Any] | None,
    metadata_payload: dict[str, Any] | None,
) -> dict[str, Any]:
    metadata_entry = metadata_entry or {}
    specimen_code = str(metadata_entry.get("specimen_code", ""))
    metadata_entries = get_metadata_entries(metadata_payload or {}, patient_id)
    candidates = build_candidates(patient_payload, specimen_code=specimen_code or None)
    return {
        "patient_id": patient_id,
        "specimen_code": specimen_code,
        "collected_time": metadata_entry.get("collected_time", ""),
        "specimen_site": metadata_entry.get("specimen_site", ""),
        "candidates": candidates,
        "related_specimens": extract_related_specimens(
            patient_payload,
            specimen_code,
            metadata_entries=metadata_entries or None,
        ),
    }


def render_prompt(patient_payload: dict[str, Any]) -> str:
    template = load_prompt_template(PROMPT_TASK)
    return (
        f"{template.rstrip()}\n\n"
        "### Real input payload ###\n"
        f"{json.dumps(patient_payload, ensure_ascii=False, indent=2)}\n"
    )


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def _sortable_patient_key(value: str) -> tuple[int, str]:
    text = str(value).strip()
    if text.isdigit():
        return (0, f"{int(text):09d}")
    return (1, text)


def main() -> None:
    args = parse_args()
    source_payload = load_json(args.source_json)
    patients = source_payload.get("patients", {})
    if not isinstance(patients, dict):
        raise TypeError("Expected source payload to contain a dict at 'patients'")

    metadata_payload: dict[str, Any] | None = None
    if args.metadata_json.exists():
        metadata_payload = load_json(args.metadata_json)

    if args.patient_id:
        patient_ids = [str(args.patient_id)]
    else:
        patient_ids = [
            str(patient_id)
            for patient_id in sorted(patients.keys(), key=_sortable_patient_key)
            if str(patient_id).strip().isdigit()
        ]
    results: dict[str, Any] = {}

    for patient_id in patient_ids:
        patient_payload = patients.get(patient_id)
        if not isinstance(patient_payload, dict):
            raise KeyError(f"Patient id not found or not a dict: {patient_id}")

        metadata_entries = get_metadata_entries(metadata_payload or {}, patient_id)
        if not metadata_entries:
            metadata_entries = derive_specimen_entries_from_source(patient_payload)
        if not metadata_entries:
            metadata_entries = [get_metadata_entry(metadata_payload or {}, patient_id)]

        patient_results: list[dict[str, Any]] = []
        for metadata_entry in metadata_entries:
            llm_input = build_specimen_payload(
                patient_id,
                patient_payload,
                metadata_entry,
                metadata_payload,
            )
            prompt = render_prompt(llm_input)

            if args.dry_run:
                print(prompt)
                continue

            output_text, usage = send_to_llm(prompt, model=args.model)
            try:
                patient_results.append(
                    augment_curated_output(
                        json.loads(output_text),
                        llm_input.get("candidates", []),
                        llm_input.get("related_specimens", []),
                    )
                )
            except json.JSONDecodeError as exc:
                specimen_code = metadata_entry.get("specimen_code", "")
                raise RuntimeError(
                    f"LLM output for patient {patient_id} specimen {specimen_code} was not valid JSON"
                ) from exc

            specimen_code = metadata_entry.get("specimen_code", "")
            if usage:
                print(
                    f"patient={patient_id} specimen={specimen_code} "
                    f"prompt_tokens={usage.get('prompt_tokens')} "
                    f"completion_tokens={usage.get('completion_tokens')} "
                    f"total_tokens={usage.get('total_tokens')}"
                )
            else:
                print(f"patient={patient_id} specimen={specimen_code} completed")

        if not args.dry_run:
            results[patient_id] = patient_results

    if not args.dry_run:
        write_json(args.output, results)
        print(f"Wrote results to {args.output}")


if __name__ == "__main__":
    main()

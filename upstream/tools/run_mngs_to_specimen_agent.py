from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from tools import pathogen_normalization as pathogen_names
from tools.llm_requester import load_prompt_template, send_to_llm

DEFAULT_SOURCE = Path("outputs") / "mngs_candidate_microbes.json"
DEFAULT_METADATA = Path("outputs") / "specimen_lookup_summary_mNGS_grouped_with_type_filtered.json"
DEFAULT_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs.json"
PROMPT_TASK = "mngs_to_specimen_agent"
NAME_SANITIZER = pathogen_names.NAME_SANITIZER
MATCH_ALIASES = pathogen_names.alias_mapping()
DISPLAY_NAME_MAP = pathogen_names.display_mapping()

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
ORAL_UPPER_AIRWAY_GENERA = ORAL_COLONIZER_GENERA + (
    "actinomyces",
    "alloprevotella",
    "capnocytophaga",
    "dialister",
    "granulicatella",
    "lacticaseibacillus",
    "lactobacillus",
    "lancefieldella",
    "ligilactobacillus",
    "limosilactobacillus",
    "megasphaera",
    "moraxella",
    "neisseria",
    "rothia",
    "schaalia",
    "segatella",
    "streptococcus",
    "tannerella",
    "treponema",
    "veillonella",
)
CATEGORY_TO_GROUP = {"1.Bac": "bacterial", "2.Fungi": "fungal", "3.Virus": "viral"}
GROUP_KEYS = ("bacterial", "viral", "fungal", "others")
LIKELIHOOD_LEVEL_KEYS = ("high", "medium", "low_colonizer")
HIGH_READ_THRESHOLDS = {"bacterial": 120, "viral": 80, "fungal": 80, "others": 120}
MEDIUM_READ_THRESHOLDS = {"bacterial": 20, "viral": 20, "fungal": 20, "others": 30}
BACKFILL_MEDIUM_LIMIT = {"bacterial": 4, "viral": 3, "fungal": 3, "others": 2}
BACKFILL_LOW_LIMIT = {"bacterial": 3, "viral": 2, "fungal": 2, "others": 2}
BALF_ORAL_MAX_KEEP = 3
LOW_COLONIZER_ORAL_MAX_KEEP = 2
AGGRESSIVE_ALL_LEVELS_RECALL = True
LOW_COLONIZER_MIN_READS = 5
LOW_COLONIZER_BACKGROUND_PREFIXES = (
    "anelloviridae",
    "alphatorquevirus",
    "torquetenovirus",
    "torqueteno",
    "betatorquevirus",
    "gammatorquevirus",
    "ttv",
)
BACKGROUND_VIRUS_PREFIXES = LOW_COLONIZER_BACKGROUND_PREFIXES + (
    "cressdnaviricota",
    "circoviridae",
    "genomoviridae",
)
ENVIRONMENTAL_LOW_PRIORITY_PREFIXES = (
    "agrobacterium",
    "anthropogastromicrobium",
    "arcicella",
    "brachybacterium",
    "cloacibacterium",
    "delftia",
    "flavobacterium",
    "flectobacillus",
    "fluviibacter",
    "janibacter",
    "kocuria",
    "limnohabitans",
    "methylobacterium",
    "microbacterium",
    "micrococcus",
    "novosphingobium",
    "ornithinimicrobium",
    "paracoccus",
    "phenylobacterium",
    "pseudozyma",
    "ralstonia",
    "rhodococcus",
    "sphingobium",
    "sphingomonas",
    "talaromyces",
    "thermomonas",
    "undibacterium",
)
LOW_PRIORITY_FUNGAL_PREFIXES = (
    "alternaria",
    "cladosporium",
    "irpex",
    "malassezia",
    "starmerella",
)
LOW_COLONIZER_SKIN_PROTECTED_EXACT = {
    "corynebacteriumstriatum",
    "staphylococcusaureus",
    "staphylococcushaemolyticus",
}
LOW_COLONIZER_ORAL_CAP_EXEMPT_EXACT = {
    "mycoplasmasalivarium",
    "streptococcuspneumoniae",
}
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
        "--include-non-digit-patient-ids",
        action="store_true",
        help=(
            "When running all patients, include IDs such as K041/T001. "
            "Default keeps the legacy numeric-only behavior."
        ),
    )
    parser.add_argument(
        "--model",
        default="gpt-5.4",
        help="OpenAI model name (default: gpt-5.4)",
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
    return pathogen_names.canonical_key(name)


def canonical_display_name(name: str) -> str:
    return pathogen_names.display_name(name)


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


def _is_oral_upper_airway_name(normalized: str) -> bool:
    if not normalized:
        return False
    return any(normalized.startswith(genus) for genus in ORAL_UPPER_AIRWAY_GENERA)


def _is_background_virus_name(normalized: str) -> bool:
    if not normalized:
        return False
    return any(normalized.startswith(prefix) for prefix in BACKGROUND_VIRUS_PREFIXES)


def _is_environmental_low_priority_name(normalized: str) -> bool:
    if not normalized:
        return False
    return any(normalized.startswith(prefix) for prefix in ENVIRONMENTAL_LOW_PRIORITY_PREFIXES)


def _is_low_read_manual_style_name(group_name: str, normalized: str, read_value: float) -> bool:
    if group_name == "viral" and (
        normalized in FORCE_KEEP_VIRUSES or normalized.startswith("humanalphaherpesvirus")
    ):
        return read_value < 20
    if group_name == "fungal" and normalized.startswith("candida"):
        return read_value < 10
    if normalized == "mycobacteriumtuberculosis":
        return read_value < 10
    return False


def _is_non_priority_staph_coryne_name(normalized: str) -> bool:
    if normalized in {"corynebacteriumstriatum", "staphylococcushaemolyticus"}:
        return False
    return normalized.startswith("corynebacterium") or normalized.startswith("staphylococcus")


def _is_non_aeruginosa_pseudomonas_name(normalized: str) -> bool:
    return normalized.startswith("pseudomonas") and normalized != "pseudomonasaeruginosa"


def _is_non_baumannii_acinetobacter_name(normalized: str) -> bool:
    return normalized.startswith("acinetobacter") and normalized != "acinetobacterbaumannii"


def _is_skin_colonizer_name(normalized: str) -> bool:
    return (
        normalized.startswith("corynebacterium")
        or normalized.startswith("staphylococcus")
        or normalized.startswith("cutibacterium")
    )


def _is_low_priority_fungal_name(normalized: str) -> bool:
    if not normalized:
        return False
    return any(normalized.startswith(prefix) for prefix in LOW_PRIORITY_FUNGAL_PREFIXES)


def _is_low_colonizer_protected_name(group_name: str, normalized: str) -> bool:
    if normalized in {PNEUMOCYSTIS, "mycobacteriumtuberculosis"}:
        return True
    if group_name == "viral" and (
        normalized in FORCE_KEEP_VIRUSES or normalized.startswith("humanalphaherpesvirus")
    ):
        return True
    if group_name == "bacterial" and normalized in LOW_COLONIZER_SKIN_PROTECTED_EXACT:
        return True
    return False


def _is_oral_cap_candidate_name(normalized: str) -> bool:
    return _is_oral_upper_airway_name(normalized) and normalized not in LOW_COLONIZER_ORAL_CAP_EXEMPT_EXACT


def _oral_low_colonizer_priority(normalized: str) -> int:
    if normalized.startswith(("parvimonas", "prevotella", "fusobacterium", "veillonella", "segatella")):
        return 4
    if normalized.startswith(("actinomyces", "capnocytophaga", "tannerella", "treponema")):
        return 3
    if normalized.startswith(("lactobacillus", "ligilactobacillus", "limosilactobacillus", "lacticaseibacillus")):
        return 1
    if normalized.startswith(("neisseria", "moraxella", "rothia", "streptococcus")):
        return 1
    return 2


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


def _rank_priority_from_info(rank_info: Any) -> int | None:
    if not isinstance(rank_info, dict):
        return None
    try:
        return int(rank_info.get("rank_priority"))
    except (TypeError, ValueError):
        return None


def _possibility_level_from_info(rank_info: Any) -> str:
    if not isinstance(rank_info, dict):
        return ""
    return str(rank_info.get("possibility_level") or "").strip().lower()


def _has_strong_rank(rank_priority: int | None) -> bool:
    return rank_priority is not None and rank_priority <= 2


def _has_acceptable_rank(rank_priority: int | None) -> bool:
    return rank_priority is not None and rank_priority <= 3


def _has_weak_rank(rank_priority: int | None) -> bool:
    return rank_priority is not None and rank_priority >= 4


def _very_high_read_threshold(group_name: str) -> float:
    return {
        "bacterial": 1000.0,
        "viral": 1000.0,
        "fungal": 500.0,
        "others": 1000.0,
    }.get(group_name, 1000.0)


def _build_group_top_reads(current_candidates: list[dict[str, Any]]) -> dict[str, float]:
    top_reads = {group_name: 0.0 for group_name in GROUP_KEYS}
    for candidate in current_candidates:
        if not isinstance(candidate, dict):
            continue
        group_name = CATEGORY_TO_GROUP.get(str(candidate.get("source_category", "")))
        if not group_name:
            continue
        top_reads[group_name] = max(top_reads[group_name], _to_float(candidate.get("sec_hit", 0)))
    return top_reads


def _has_high_read_dominance(group_name: str, reads: Any, group_top_reads: dict[str, float]) -> bool:
    read_value = _to_float(reads)
    top_read = group_top_reads.get(group_name, 0.0)
    if top_read <= 0:
        return False
    if read_value >= top_read:
        return True
    return read_value >= max(HIGH_READ_THRESHOLDS.get(group_name, 120), top_read * 0.25)


def _build_candidate_evidence_map(
    current_candidates: list[dict[str, Any]],
) -> dict[tuple[str, str], dict[str, Any]]:
    evidence: dict[tuple[str, str], dict[str, Any]] = {}

    def evidence_key(candidate: dict[str, Any]) -> tuple[float, float, float]:
        rank_info = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
        rank_priority = _rank_priority_from_info(rank_info)
        rank_score = -(rank_priority if rank_priority is not None else 999.0)
        possibility_score = {"high": 3.0, "medium": 2.0, "low": 1.0}.get(
            _possibility_level_from_info(rank_info),
            0.0,
        )
        reads = _to_float(candidate.get("sec_hit", candidate.get("reads", 0)))
        return rank_score, possibility_score, reads

    for candidate in current_candidates:
        if not isinstance(candidate, dict):
            continue
        group_name = CATEGORY_TO_GROUP.get(str(candidate.get("source_category", "")))
        normalized = normalize_organism_name(candidate.get("organism_name"))
        if not group_name or not normalized:
            continue
        key = (group_name, normalized)
        current = evidence.get(key)
        if current is None or evidence_key(candidate) > evidence_key(current):
            evidence[key] = candidate
    return evidence


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


def _genus_key_from_name(name: Any) -> str:
    text = str(name or "").strip()
    if not text:
        return ""
    first_token = text.replace("[", " ").replace("]", " ").split()[0] if text.split() else ""
    return NAME_SANITIZER.sub("", first_token.lower())


def _low_colonizer_merge_key(group_name: str, item: dict[str, Any], normalized: str) -> str:
    if group_name not in {"bacterial", "fungal"}:
        return ""
    if _is_low_colonizer_protected_name(group_name, normalized):
        return ""
    return _genus_key_from_name(item.get("name"))


def _is_low_colonizer_merge_candidate(group_name: str, normalized: str) -> bool:
    if group_name == "bacterial":
        return (
            _is_oral_upper_airway_name(normalized)
            or _is_environmental_low_priority_name(normalized)
            or _is_skin_colonizer_name(normalized)
            or _is_non_aeruginosa_pseudomonas_name(normalized)
        )
    if group_name == "fungal":
        return _is_low_priority_fungal_name(normalized) or normalized.startswith("candida")
    return False


def _prune_low_colonizer_noise(
    tiers: dict[str, dict[str, list[dict[str, Any]]]],
    candidate_evidence: dict[tuple[str, str], dict[str, Any]] | None = None,
) -> None:
    low_group = tiers.get("low_colonizer")
    if not isinstance(low_group, dict):
        return

    evidence_map = candidate_evidence or {}
    primary_genus_by_group = {group_name: set() for group_name in GROUP_KEYS}
    for tier_name in ("high", "medium"):
        tier_group = tiers.get(tier_name)
        if not isinstance(tier_group, dict):
            continue
        for group_name in GROUP_KEYS:
            for item in tier_group.get(group_name, []):
                if not isinstance(item, dict):
                    continue
                normalized = normalize_organism_name(item.get("name"))
                if not normalized:
                    continue
                genus_key = _low_colonizer_merge_key(group_name, item, normalized)
                if genus_key:
                    primary_genus_by_group[group_name].add(genus_key)

    def _rank_priority_for(group_name: str, normalized: str) -> int | None:
        evidence = evidence_map.get((group_name, normalized), {})
        ranking = evidence.get("ranking") if isinstance(evidence, dict) else None
        return _rank_priority_from_info(ranking)

    def _item_keep_score(group_name: str, item: dict[str, Any]) -> tuple[int, int, float]:
        normalized = normalize_organism_name(item.get("name"))
        rank_priority = _rank_priority_for(group_name, normalized)
        has_rank = 1 if rank_priority is not None else 0
        rank_score = -(rank_priority if rank_priority is not None else 999)
        return has_rank, rank_score, _to_float(item.get("reads", 0))

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
            reads = _to_float(item.get("reads", 0))
            rank_priority = _rank_priority_for(group_name, normalized)
            background_virus_allowed = False

            if group_name == "viral" and _is_background_virus_name(normalized):
                if any(normalized.startswith(prefix) for prefix in LOW_COLONIZER_BACKGROUND_PREFIXES):
                    if reads < 100:
                        continue
                    background_virus_allowed = _has_strong_rank(rank_priority)
                elif normalized.startswith(("cressdnaviricota", "circoviridae", "genomoviridae")):
                    background_virus_allowed = _has_strong_rank(rank_priority) and reads >= 100
                else:
                    background_virus_allowed = _has_acceptable_rank(rank_priority) and reads >= 100
                if not background_virus_allowed:
                    continue

            if _is_environmental_low_priority_name(normalized) and rank_priority is not None:
                if rank_priority >= 3 and reads < 100:
                    continue

            if group_name == "bacterial" and _is_non_baumannii_acinetobacter_name(normalized):
                if rank_priority is not None and rank_priority >= 3 and reads < 100:
                    continue

            if group_name == "bacterial" and _is_non_aeruginosa_pseudomonas_name(normalized):
                if rank_priority is not None and rank_priority >= 3 and reads < 200:
                    continue

            if group_name == "bacterial" and _is_skin_colonizer_name(normalized):
                if normalized in LOW_COLONIZER_SKIN_PROTECTED_EXACT:
                    if reads < LOW_COLONIZER_MIN_READS and not _has_acceptable_rank(rank_priority):
                        continue
                elif _has_strong_rank(rank_priority) and reads >= 50:
                    pass
                elif rank_priority == 3 and reads >= 200:
                    pass
                elif rank_priority is None and reads >= 500:
                    pass
                else:
                    continue

            if group_name == "fungal":
                if normalized == PNEUMOCYSTIS:
                    pass
                elif normalized.startswith("candida"):
                    if reads < 3:
                        continue
                    if rank_priority is not None and rank_priority > 3 and reads < 10:
                        continue
                elif _is_low_priority_fungal_name(normalized):
                    if not (_has_strong_rank(rank_priority) and reads >= 100):
                        continue

            merge_key = _low_colonizer_merge_key(group_name, item, normalized)
            if (
                merge_key
                and merge_key in primary_genus_by_group[group_name]
                and not _is_low_colonizer_protected_name(group_name, normalized)
            ):
                continue

            if (
                reads < LOW_COLONIZER_MIN_READS
                and not background_virus_allowed
                and not _is_low_colonizer_low_read_allowed(normalized)
            ):
                continue
            kept.append(item)

        passthrough: list[dict[str, Any]] = []
        best_by_merge_key: dict[str, dict[str, Any]] = {}
        for item in kept:
            normalized = normalize_organism_name(item.get("name"))
            merge_key = _low_colonizer_merge_key(group_name, item, normalized)
            if merge_key and _is_low_colonizer_merge_candidate(group_name, normalized):
                current = best_by_merge_key.get(merge_key)
                if current is None or _item_keep_score(group_name, item) > _item_keep_score(group_name, current):
                    best_by_merge_key[merge_key] = item
            else:
                passthrough.append(item)
        final_items = passthrough + list(best_by_merge_key.values())
        if group_name == "bacterial":
            oral_items = [
                item
                for item in final_items
                if _is_oral_cap_candidate_name(normalize_organism_name(item.get("name")))
            ]
            non_oral_items = [
                item
                for item in final_items
                if not _is_oral_cap_candidate_name(normalize_organism_name(item.get("name")))
            ]
            oral_items = sorted(
                oral_items,
                key=lambda item: (
                    _item_keep_score(group_name, item)[:2],
                    _oral_low_colonizer_priority(normalize_organism_name(item.get("name"))),
                    _to_float(item.get("reads", 0)),
                ),
                reverse=True,
            )[:LOW_COLONIZER_ORAL_MAX_KEEP]
            final_items = non_oral_items + oral_items
        low_group[group_name] = final_items


def _choose_likelihood_tier(
    *,
    group_name: str,
    normalized: str,
    reads: Any,
    specimen_site: str,
    related_supported: bool,
    has_strong_non_oral_bacteria: bool,
    dominance_supported: bool = True,
    rank_info: dict[str, Any] | None = None,
) -> str:
    read_value = _to_float(reads)
    rank_priority = _rank_priority_from_info(rank_info)
    high_read_threshold = HIGH_READ_THRESHOLDS.get(group_name, 120)
    medium_read_threshold = MEDIUM_READ_THRESHOLDS.get(group_name, 20)
    high_read = read_value >= high_read_threshold
    medium_read = read_value >= medium_read_threshold
    very_high_read = read_value >= _very_high_read_threshold(group_name)
    manual_style_protected = _is_force_keep_name(group_name, normalized, read_value)
    is_background_virus = group_name == "viral" and _is_background_virus_name(normalized)
    is_oral_airway = group_name == "bacterial" and _is_oral_upper_airway_name(normalized)
    is_environmental = _is_environmental_low_priority_name(normalized)
    is_skin_or_generic_bacterial = group_name == "bacterial" and (
        _is_non_priority_staph_coryne_name(normalized) or _is_non_aeruginosa_pseudomonas_name(normalized)
    )
    is_low_priority_fungus = group_name == "fungal" and normalized.startswith("malassezia")
    high_disallowed = (
        is_background_virus
        or is_oral_airway
        or is_environmental
        or is_skin_or_generic_bacterial
        or is_low_priority_fungus
        or _is_low_read_manual_style_name(group_name, normalized, read_value)
    )

    if is_background_virus:
        return "low_colonizer"

    if is_oral_airway:
        if has_strong_non_oral_bacteria and read_value < max(120.0, high_read_threshold):
            return "low_colonizer"
        if _has_strong_rank(rank_priority) and read_value >= max(120.0, medium_read_threshold):
            return "medium"
        if rank_priority == 3 and very_high_read:
            return "medium"
        if related_supported and read_value >= high_read_threshold:
            return "medium"
        return "low_colonizer"

    if is_environmental:
        if _has_strong_rank(rank_priority) and very_high_read:
            return "medium"
        return "low_colonizer"

    if is_low_priority_fungus:
        if _has_strong_rank(rank_priority) and read_value >= _very_high_read_threshold(group_name):
            return "medium"
        return "low_colonizer"

    if is_skin_or_generic_bacterial:
        if _has_strong_rank(rank_priority) and medium_read:
            return "medium"
        if rank_priority == 3 and read_value >= _very_high_read_threshold(group_name):
            return "medium"
        return "low_colonizer"

    if _is_low_read_manual_style_name(group_name, normalized, read_value):
        return "low_colonizer"

    # High is evidence-first: rank 1-2 plus meaningful reads, with no background/colonizer flags.
    if not high_disallowed and _has_strong_rank(rank_priority) and high_read and dominance_supported:
        return "high"

    # Manual-style organisms are recall guardrails, not direct high-confidence evidence.
    if manual_style_protected:
        if _has_acceptable_rank(rank_priority) and medium_read:
            return "medium"
        if related_supported and medium_read:
            return "medium"
        return "low_colonizer"

    if related_supported and medium_read and not high_disallowed:
        return "medium"
    if _has_acceptable_rank(rank_priority) and medium_read:
        return "medium"
    if very_high_read and not _has_weak_rank(rank_priority):
        return "medium"

    return "low_colonizer"


def _should_backfill_candidate(
    *,
    candidate: dict[str, Any],
    group_name: str,
    specimen_site: str,
    related_supported: bool,
    has_strong_non_oral_bacteria: bool,
    dominance_supported: bool = True,
) -> str | None:
    normalized = normalize_organism_name(candidate.get("organism_name"))
    if not normalized:
        return None
    reads = _to_float(candidate.get("sec_hit", 0))
    tier_name = _choose_likelihood_tier(
        group_name=group_name,
        normalized=normalized,
        reads=reads,
        specimen_site=specimen_site,
        related_supported=related_supported,
        has_strong_non_oral_bacteria=has_strong_non_oral_bacteria,
        dominance_supported=dominance_supported,
        rank_info=candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else None,
    )
    if tier_name in {"high", "medium"}:
        return tier_name
    if _is_force_keep_name(group_name, normalized, reads) or tier_name == "low_colonizer":
        return "low_colonizer"
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
    candidate_evidence = _build_candidate_evidence_map(current_candidates)
    group_top_reads = _build_group_top_reads(current_candidates)
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
                dominance_supported=_has_high_read_dominance(group_name, item.get("reads", 0), group_top_reads),
                rank_info=(
                    candidate_evidence.get((group_name, normalized), {}).get("ranking")
                    if isinstance(candidate_evidence.get((group_name, normalized), {}).get("ranking"), dict)
                    else None
                ),
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
            dominance_supported=_has_high_read_dominance(group_name, candidate.get("sec_hit", 0), group_top_reads),
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

    _prune_low_colonizer_noise(tiers, candidate_evidence)
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
    def _merged_selection_reasons(*values: Any) -> list[str]:
        merged: list[str] = []
        seen: set[str] = set()
        for value in values:
            if isinstance(value, list):
                items = value
            elif value:
                items = [value]
            else:
                items = []
            for item in items:
                text = str(item).strip()
                if text and text not in seen:
                    merged.append(text)
                    seen.add(text)
        return merged

    def _candidate_keep_key(candidate: dict[str, Any]) -> tuple[float, float, float, float, float]:
        ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
        rank_priority_raw = ranking.get("rank_priority")
        rank_priority = _to_float(rank_priority_raw) if rank_priority_raw not in (None, "") else 999.0
        has_ranking = 1.0 if ranking else 0.0
        possibility_order = {"high": 3.0, "medium": 2.0, "low": 1.0}
        possibility_score = possibility_order.get(str(ranking.get("possibility_level", "")).lower(), 0.0)
        sec_hit_value = _to_float(candidate.get("sec_hit", candidate.get("reads", 0)))
        has_reasons = 1.0 if isinstance(candidate.get("selection_reasons"), list) and candidate.get("selection_reasons") else 0.0
        richness = float(
            sum(
                1
                for key in ("reads", "selection_reasons", "ranking")
                if key in candidate and candidate.get(key) not in (None, "", [], {})
            )
        )
        # Higher is better. Lower rank_priority means stronger ranked evidence, so use negative value.
        return (has_ranking, -rank_priority, possibility_score, sec_hit_value, has_reasons + richness)

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
            organism_name = str(candidate.get("organism_name", candidate.get("name", "")))
            source_category = str(candidate.get("source_category", candidate.get("type", "")))
            sec_hit = candidate.get("sec_hit", candidate.get("reads", 0))
            key = (organism_name, source_category)
            current = deduped.get(key)
            selection_reasons = _merged_selection_reasons(candidate.get("selection_reasons"))
            if current is None or _candidate_keep_key(candidate) > _candidate_keep_key(current):
                if current is not None:
                    selection_reasons = _merged_selection_reasons(
                        current.get("selection_reasons"),
                        selection_reasons,
                    )
                normalized_candidate: dict[str, Any] = {
                    "organism_name": organism_name,
                    "sec_hit": sec_hit,
                    "source_category": source_category,
                    "selection_reasons": selection_reasons,
                }
                if "reads" in candidate and candidate.get("reads") is not None:
                    normalized_candidate["reads"] = candidate.get("reads")
                if isinstance(candidate.get("ranking"), dict):
                    normalized_candidate["ranking"] = dict(candidate.get("ranking"))
                deduped[key] = normalized_candidate
            elif current is not None:
                current["selection_reasons"] = _merged_selection_reasons(
                    current.get("selection_reasons"),
                    selection_reasons,
                )
    return sorted(
        deduped.values(),
        key=lambda item: (
            str(item.get("source_category", "")),
            _to_float((item.get("ranking") or {}).get("rank_priority", 999)),
            -{"high": 3.0, "medium": 2.0, "low": 1.0}.get(
                str((item.get("ranking") or {}).get("possibility_level", "")).lower(),
                0.0,
            ),
            -_to_float(item.get("sec_hit", item.get("reads", 0))),
            str(item.get("organism_name", "")),
        ),
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
                "collected_time": record.get("collected_time", ""),
                "specimen_site": record.get("specimen_site", ""),
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
        patient_ids = []
        for patient_id in sorted(patients.keys(), key=_sortable_patient_key):
            normalized_patient_id = str(patient_id).strip()
            if not normalized_patient_id:
                continue
            if normalized_patient_id.isdigit() or args.include_non_digit_patient_ids:
                patient_ids.append(normalized_patient_id)
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

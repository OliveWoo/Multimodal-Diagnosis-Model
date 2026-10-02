from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Sequence

from tools.run_mngs_to_specimen_agent import (
    FORCE_KEEP_BACTERIA,
    FORCE_KEEP_VIRUSES,
    ORAL_ANAEROBE_PRIORITY,
    PNEUMOCYSTIS,
    normalize_organism_name,
    should_force_keep_candida,
)

DEFAULT_INPUT = Path("outputs") / "mngs_to_specimen_outputs_main_baseline.json"
DEFAULT_RANKED = Path("mngs_candidate_microbes_ranked.json")
DEFAULT_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs_main_rank_opt.json"
DEFAULT_SUMMARY = Path("outputs") / "mngs_rank_read_precision_opt_main.json"

GROUP_KEYS = ("bacterial", "viral", "fungal", "others")
TIER_KEYS = ("high", "medium", "low_colonizer")

GENUS_SANITIZER = re.compile(r"[^a-z]+")
GENERIC_TAXON_SUFFIXES = ("aceae", "ales")
BACKGROUND_VIRUS_PREFIXES = (
    "anelloviridae",
    "torqueteno",
    "torquetenovirus",
    "betatorquevirus",
    "alphatorquevirus",
    "gammatorquevirus",
    "ttvirus",
    "genomoviridae",
    "cressdnaviricota",
    "circoviridae",
    "gemykibivirus",
)
ORAL_UPPER_AIRWAY_PREFIXES = (
    "parvimonas",
    "prevotella",
    "porphyromonas",
    "fusobacterium",
    "veillonella",
    "rothia",
    "oribacterium",
    "slackia",
    "leptotrichia",
    "capnocytophaga",
    "actinomyces",
    "aggregatibacter",
    "treponema",
)
ENVIRONMENTAL_OR_LOW_PRIORITY_PREFIXES = (
    "acidovorax",
    "brachybacterium",
    "brevundimonas",
    "chryseobacterium",
    "cutibacterium",
    "janibacter",
    "kocuria",
    "microbacterium",
    "novosphingobium",
    "paracoccus",
    "pseudoxanthomonas",
    "sphingobium",
    "sphingomonas",
    "stutzerimonas",
    "talaromyces",
    "undibacterium",
)
ENVIRONMENTAL_EXACT = {
    "acinetobacterjohnsonii",
    "pseudomonassoli",
    "pseudomonasalcaligenes",
    "pseudomonasjaponica",
    "pseudomonaskhazarica",
    "pseudomonasceruminis",
    "xanthomonadaceae",
}
COMMON_MANUAL_LIKE = {
    "achromobacterxylosoxidans",
    "acinetobacterbaumannii",
    "burkholderiacenocepacia",
    "candidaalbicans",
    "candidaglabrata",
    "candidahaemuloniicomplex",
    "candidatropicalis",
    "corynebacteriumresistens",
    "corynebacteriumstriatum",
    "enterococcusfaecium",
    "humanalphaherpesvirus1",
    "humanalphaherpesvirus2",
    "humanbetaherpesvirus5",
    "humanbetaherpesvirus6b",
    "humanbetaherpesvirus7",
    "humangammaherpesvirus4",
    "klebsiellapneumoniae",
    "legionellapneumophila",
    "mycoplasmasalivarium",
    "pneumocystisjirovecii",
    "pseudomonasaeruginosa",
    "staphylococcusaureus",
    "staphylococcushaemolyticus",
    "stenotrophomonasmaltophilia",
}

# Pulmonary coarse filter design:
# - keep aspiration-relevant oral/gut anaerobes for later max/LLM review;
# - remove taxa that are usually skin/environment/background and have weak direct
#   pulmonary-pathogen plausibility. This is intentionally a coarse pre-filter,
#   not a final clinical decision.
PULMONARY_CONDITIONAL_ANAEROBE_PREFIXES = (
    "acidaminococcus",
    "alistipes",
    "anaerococcus",
    "bacteroides",
    "fusobacterium",
    "parvimonas",
    "peptostreptococcus",
    "phocaeicola",
    "porphyromonas",
    "prevotella",
    "ruthenibacterium",
    "veillonella",
)
PULMONARY_BACKGROUND_PREFIXES = (
    "acidovorax",
    "akkermansia",
    "brachybacterium",
    "brevundimonas",
    "cutibacterium",
    "dolosigranulum",
    "janibacter",
    "kocuria",
    "malassezia",
    "microbacterium",
    "micrococcus",
    "novosphingobium",
    "paracoccus",
    "pseudoxanthomonas",
    "sphingobium",
    "sphingomonas",
    "stutzerimonas",
    "undibacterium",
)
PULMONARY_LOW_PRIORITY_FUNGAL_PREFIXES = (
    "exophiala",
    "penicillium",
    "talaromyces",
)
PULMONARY_BACKGROUND_VIRUS_EXACT = {
    "chickenassociatedgemycircularvirus1",
    "cucumbergreenmottlemosaicvirus",
    "humanpolyomavirus6",
    "humangammaherpesvirus4",
    "humanbetaherpesvirus6b",
    "humanbetaherpesvirus7",
    "tomatobrownrugosefruitvirus",
}
PULMONARY_BACKGROUND_VIRUS_PREFIXES = (
    "gemycircularvirus",
    "gemykibivirus",
)
PULMONARY_BACKGROUND_EXEMPT = {
    "corynebacteriumresistens",
    "corynebacteriumstriatum",
    "humanalphaherpesvirus1",
    "humanalphaherpesvirus2",
    "humanbetaherpesvirus5",
    "pneumocystisjirovecii",
    "staphylococcusaureus",
    "staphylococcushaemolyticus",
    "staphylococcuslugdunensis",
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "No-answer mNGS-to-specimen optimizer using rank_priority, reads, tier, "
            "and manual-style pathogen protection."
        )
    )
    parser.add_argument("--input-json", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--ranked-json", type=Path, default=DEFAULT_RANKED)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-json", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--remove-medium", action="store_true", default=True)
    parser.add_argument("--disable-remove-medium", action="store_true")
    parser.add_argument("--remove-low", action="store_true", default=True)
    parser.add_argument("--disable-remove-low", action="store_true")
    parser.add_argument("--demote-medium", action="store_true", default=True)
    parser.add_argument("--disable-demote-medium", action="store_true")
    parser.add_argument("--background-virus-max-reads", type=float, default=100.0)
    parser.add_argument("--medium-env-max-reads", type=float, default=250.0)
    parser.add_argument("--medium-low-priority-max-reads", type=float, default=50.0)
    parser.add_argument("--low-env-max-reads", type=float, default=1000.0)
    parser.add_argument("--low-rank-remove-cutoff", type=int, default=4)
    parser.add_argument("--low-rank-remove-max-reads", type=float, default=120.0)
    parser.add_argument("--low-possibility-remove-max-reads", type=float, default=3.0)
    parser.add_argument("--balf-medium-oral-max-keep", type=int, default=-1)
    parser.add_argument("--balf-low-oral-max-keep", type=int, default=-1)
    parser.add_argument("--balf-non-oral-strong-read", type=float, default=30.0)
    parser.add_argument("--protect-rank-priority", type=int, default=2)
    parser.add_argument("--protect-read-bacterial", type=float, default=1000.0)
    parser.add_argument("--protect-read-viral", type=float, default=1000.0)
    parser.add_argument("--protect-read-fungal", type=float, default=500.0)
    parser.add_argument("--protect-read-others", type=float, default=1000.0)
    parser.add_argument(
        "--pulmonary-coarse-filter",
        action="store_true",
        help=(
            "Apply an additional pulmonary relevance coarse filter before the usual "
            "rank/read optimizer. This removes obvious non-pulmonary background taxa "
            "but keeps aspiration-relevant oral/gut anaerobes for downstream review."
        ),
    )
    return parser.parse_args(argv)


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict JSON from {path}, got {type(payload)}")
    return payload


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def _to_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _rank_priority(rank_info: dict[str, Any]) -> int | None:
    try:
        return int(rank_info.get("rank_priority"))
    except (TypeError, ValueError):
        return None


def _possibility_level(rank_info: dict[str, Any]) -> str:
    return str(rank_info.get("possibility_level") or "").strip().lower()


def _rank_sort_key(candidate: dict[str, Any]) -> tuple[int, float]:
    ranking = candidate.get("ranking") if isinstance(candidate.get("ranking"), dict) else {}
    priority = _rank_priority(ranking)
    if priority is None:
        priority = 999
    return priority, -_to_float(candidate.get("sec_hit", candidate.get("reads", 0)))


def _build_rank_lookup(payload: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    lookup: dict[tuple[str, str], dict[str, Any]] = {}
    patients = payload.get("patients")
    if not isinstance(patients, dict):
        return lookup

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for patient_payload in patients.values():
        if not isinstance(patient_payload, dict):
            continue
        for record in patient_payload.get("records") or []:
            if not isinstance(record, dict):
                continue
            specimen_code = str(record.get("specimen_code") or "").strip()
            if not specimen_code:
                continue
            for candidate in record.get("candidates") or []:
                if not isinstance(candidate, dict):
                    continue
                normalized = normalize_organism_name(candidate.get("organism_name") or candidate.get("name"))
                if normalized:
                    grouped[(specimen_code, normalized)].append(candidate)

    for key, values in grouped.items():
        lookup[key] = sorted(values, key=_rank_sort_key)[0]
    return lookup


def _normalize_tiers(raw: Any) -> dict[str, dict[str, list[dict[str, Any]]]]:
    tiers = {tier: {group: [] for group in GROUP_KEYS} for tier in TIER_KEYS}
    if not isinstance(raw, dict):
        return tiers
    for tier in TIER_KEYS:
        level = raw.get(tier)
        if tier == "low_colonizer" and not isinstance(level, dict):
            level = raw.get("low")
        if not isinstance(level, dict):
            continue
        for group in GROUP_KEYS:
            values = level.get(group)
            if isinstance(values, list):
                tiers[tier][group] = [item for item in values if isinstance(item, dict)]
    return tiers


def _item_sort_key(item: dict[str, Any]) -> tuple[float, str]:
    return -_to_float(item.get("reads", 0)), str(item.get("name", ""))


def _dedupe(values: list[dict[str, Any]]) -> list[dict[str, Any]]:
    best: dict[str, dict[str, Any]] = {}
    for item in values:
        normalized = normalize_organism_name(item.get("name"))
        if not normalized:
            continue
        current = best.get(normalized)
        if current is None or _to_float(item.get("reads", 0)) > _to_float(current.get("reads", 0)):
            best[normalized] = item
    return sorted(best.values(), key=_item_sort_key)


def _rebuild_primary_pathogens(tiers: dict[str, dict[str, list[dict[str, Any]]]]) -> dict[str, list[dict[str, Any]]]:
    primary = {group: [] for group in GROUP_KEYS}
    seen: set[tuple[str, str]] = set()
    for tier in ("high", "medium"):
        for group in GROUP_KEYS:
            for item in tiers[tier][group]:
                normalized = normalize_organism_name(item.get("name"))
                key = (group, normalized)
                if normalized and key not in seen:
                    seen.add(key)
                    primary[group].append(item)
    for group in GROUP_KEYS:
        primary[group] = _dedupe(primary[group])
    return primary


def _extract_genus(display_name: str, normalized: str) -> str:
    text = str(display_name or "").strip().lower()
    if text:
        token = text.split()[0]
        genus = GENUS_SANITIZER.sub("", token)
        if genus:
            return genus
    return GENUS_SANITIZER.sub("", normalized[:25])


def _is_generic_taxon_label(display_name: str, normalized: str) -> bool:
    text = " ".join(str(display_name or "").strip().lower().split())
    if not text:
        return False
    tokens = text.split(" ")
    first = GENUS_SANITIZER.sub("", tokens[0]) if tokens else ""
    if re.search(r"\bspp?\.?$", text):
        return True
    if len(tokens) == 1 and first:
        return True
    if first.endswith(GENERIC_TAXON_SUFFIXES):
        return True
    return normalized.endswith("sp") or normalized.endswith("spp")


def _is_background_virus(normalized: str) -> bool:
    return any(normalized.startswith(prefix) for prefix in BACKGROUND_VIRUS_PREFIXES)


def _is_oral_upper_airway(normalized: str) -> bool:
    if normalized in ORAL_ANAEROBE_PRIORITY:
        return True
    return any(normalized.startswith(prefix) for prefix in ORAL_UPPER_AIRWAY_PREFIXES)


def _is_environmental_or_low_priority(normalized: str) -> bool:
    if normalized in ENVIRONMENTAL_EXACT:
        return True
    return any(normalized.startswith(prefix) for prefix in ENVIRONMENTAL_OR_LOW_PRIORITY_PREFIXES)


def _is_pulmonary_conditional_anaerobe(normalized: str) -> bool:
    return any(normalized.startswith(prefix) for prefix in PULMONARY_CONDITIONAL_ANAEROBE_PREFIXES)


def _pulmonary_coarse_reasons(*, group: str, tier: str, normalized: str, reads: float) -> list[str]:
    if not normalized or normalized in PULMONARY_BACKGROUND_EXEMPT:
        return []
    if _is_pulmonary_conditional_anaerobe(normalized):
        return []

    reasons: list[str] = []
    if group == "viral" and (
        _is_background_virus(normalized)
        or normalized in PULMONARY_BACKGROUND_VIRUS_EXACT
        or any(normalized.startswith(prefix) for prefix in PULMONARY_BACKGROUND_VIRUS_PREFIXES)
    ):
        reasons.append("pulmonary_background_virus")

    if any(normalized.startswith(prefix) for prefix in PULMONARY_BACKGROUND_PREFIXES):
        reasons.append("pulmonary_background_or_unrelated_taxon")

    if (
        group == "fungal"
        and tier != "high"
        and any(normalized.startswith(prefix) for prefix in PULMONARY_LOW_PRIORITY_FUNGAL_PREFIXES)
    ):
        reasons.append("pulmonary_low_priority_fungal_taxon")

    if normalized.startswith("corynebacterium") and tier != "high":
        reasons.append("pulmonary_skin_or_environmental_taxon")
    if normalized.startswith("staphylococcus") and tier != "high":
        reasons.append("pulmonary_skin_or_environmental_taxon")

    return list(dict.fromkeys(reasons))


def _is_common_manual_like(group: str, normalized: str, reads: float) -> bool:
    if normalized.startswith("candida"):
        return True
    if normalized.startswith(("klebsiella", "enterobacter", "serratia", "escherichia", "citrobacter")):
        return True
    if normalized in COMMON_MANUAL_LIKE:
        if normalized == "stenotrophomonasmaltophilia":
            return reads >= 50
        if normalized == "corynebacteriumresistens":
            return reads >= 100
        return True
    if group == "bacterial":
        if normalized == "stenotrophomonasmaltophilia":
            return reads >= 50
        return normalized in FORCE_KEEP_BACTERIA
    if group == "viral":
        return normalized in FORCE_KEEP_VIRUSES or normalized.startswith("humanalphaherpesvirus")
    if group == "fungal":
        return normalized == PNEUMOCYSTIS or should_force_keep_candida(normalized, reads)
    return False


def _read_protect_threshold(group: str, args: argparse.Namespace) -> float:
    return {
        "bacterial": args.protect_read_bacterial,
        "viral": args.protect_read_viral,
        "fungal": args.protect_read_fungal,
        "others": args.protect_read_others,
    }[group]


def _protection_signals(
    *,
    group: str,
    normalized: str,
    reads: float,
    tier: str,
    rank_info: dict[str, Any],
    args: argparse.Namespace,
) -> list[str]:
    signals: list[str] = []
    priority = _rank_priority(rank_info)
    if tier == "high":
        signals.append("tier_high")
    if _is_common_manual_like(group, normalized, reads):
        signals.append("common_manual_like")
    if (
        priority is not None
        and priority <= args.protect_rank_priority
        and reads >= 100
        and not _is_environmental_or_low_priority(normalized)
        and not _is_background_virus(normalized)
    ):
        signals.append("top_rank_with_reads")
    if reads >= _read_protect_threshold(group, args):
        signals.append("very_high_reads")
    return signals


def _strong_non_oral_balf(tiers: dict[str, dict[str, list[dict[str, Any]]]], args: argparse.Namespace) -> bool:
    for tier in ("high", "medium"):
        for item in tiers[tier]["bacterial"]:
            normalized = normalize_organism_name(item.get("name"))
            reads = _to_float(item.get("reads", 0))
            if not normalized or _is_oral_upper_airway(normalized):
                continue
            if _is_common_manual_like("bacterial", normalized, reads) or reads >= args.balf_non_oral_strong_read:
                return True
    return False


def _specific_genus_set(tiers: dict[str, dict[str, list[dict[str, Any]]]], group: str) -> set[str]:
    genus_set: set[str] = set()
    for tier in TIER_KEYS:
        for item in tiers[tier][group]:
            name = str(item.get("name") or "")
            normalized = normalize_organism_name(name)
            if not normalized or _is_generic_taxon_label(name, normalized):
                continue
            genus = _extract_genus(name, normalized)
            if genus:
                genus_set.add(genus)
    return genus_set


def _noise_reasons(
    *,
    tier: str,
    group: str,
    name: str,
    normalized: str,
    reads: float,
    rank_info: dict[str, Any],
    specific_genus: set[str],
    strong_non_oral_balf: bool,
    specimen_site: str,
    args: argparse.Namespace,
) -> list[str]:
    reasons: list[str] = []
    priority = _rank_priority(rank_info)
    possibility = _possibility_level(rank_info)
    is_generic = _is_generic_taxon_label(name, normalized)
    genus = _extract_genus(name, normalized)

    if group == "viral" and _is_background_virus(normalized) and reads <= args.background_virus_max_reads:
        reasons.append("background_virus")
    if is_generic and genus in specific_genus:
        reasons.append("generic_shadowed_by_specific_species")
    elif is_generic and tier == "low_colonizer":
        reasons.append("low_generic_label")
    if _is_environmental_or_low_priority(normalized):
        reasons.append("environmental_or_low_priority_taxon")
    if group == "bacterial" and _is_oral_upper_airway(normalized):
        reasons.append("oral_upper_airway_flora")
    if priority is not None and priority >= 4:
        reasons.append("rank_priority_4_or_5")
    if priority == 3 and reads < 20:
        reasons.append("rank_priority_3_low_reads")
    if possibility == "low" and reads <= args.low_possibility_remove_max_reads:
        reasons.append("low_possibility_low_reads")
    if tier == "medium" and reads <= args.medium_low_priority_max_reads:
        reasons.append("medium_low_reads")
    if tier == "low_colonizer" and reads <= args.low_rank_remove_max_reads:
        reasons.append("low_colonizer_low_reads")
    if (
        group == "bacterial"
        and "balf" in specimen_site.lower()
        and strong_non_oral_balf
        and _is_oral_upper_airway(normalized)
    ):
        reasons.append("balf_oral_with_strong_non_oral")
    if normalized.startswith("corynebacterium") and normalized not in {
        "corynebacteriumstriatum",
        "corynebacteriumresistens",
    }:
        reasons.append("non_priority_corynebacterium")
    if normalized.startswith("staphylococcus") and normalized not in {
        "staphylococcusaureus",
        "staphylococcushaemolyticus",
    }:
        reasons.append("non_priority_staphylococcus")
    return reasons


def _should_remove(
    *,
    tier: str,
    group: str,
    normalized: str,
    reads: float,
    reasons: list[str],
    protected: bool,
    args: argparse.Namespace,
) -> bool:
    if protected:
        return False
    reason_set = set(reasons)

    if tier == "high":
        return False
    if tier == "medium" and not args.remove_medium:
        return False
    if tier == "low_colonizer" and not args.remove_low:
        return False

    if "background_virus" in reason_set:
        return True
    if tier == "medium":
        if "environmental_or_low_priority_taxon" in reason_set and reads <= args.medium_env_max_reads:
            return True
        if "generic_shadowed_by_specific_species" in reason_set:
            return True
        if "rank_priority_4_or_5" in reason_set and "medium_low_reads" in reason_set and reads <= 20:
            return True
        return False

    if tier == "low_colonizer":
        if "environmental_or_low_priority_taxon" in reason_set and reads <= args.low_env_max_reads:
            return True
        if "low_generic_label" in reason_set:
            return True
        if "generic_shadowed_by_specific_species" in reason_set:
            return True
        if (
            "rank_priority_4_or_5" in reason_set
            and "low_possibility_low_reads" in reason_set
            and reads <= args.low_possibility_remove_max_reads
        ):
            return True
        if "non_priority_corynebacterium" in reason_set and reads < 10:
            return True
        if "non_priority_staphylococcus" in reason_set and reads < 10:
            return True
        return False
    return False


def _should_demote_medium(
    *,
    group: str,
    reads: float,
    reasons: list[str],
    protected: bool,
    args: argparse.Namespace,
) -> bool:
    if protected or not args.demote_medium:
        return False
    reason_set = set(reasons)
    if "environmental_or_low_priority_taxon" in reason_set and reads <= 1000:
        return True
    if "rank_priority_4_or_5" in reason_set:
        return True
    if group == "bacterial" and "balf_oral_with_strong_non_oral" in reason_set:
        return True
    return False


def _record(
    entries: list[dict[str, Any]],
    counters: dict[str, Counter[str]],
    *,
    action: str,
    specimen_code: str,
    specimen_site: str,
    tier: str,
    group: str,
    name: str,
    normalized: str,
    reads: float,
    rank_info: dict[str, Any],
    reasons: list[str],
    protection: list[str],
) -> None:
    for reason in reasons or ["unspecified"]:
        counters["by_reason"][reason] += 1
    counters["by_action"][action] += 1
    counters["by_tier"][tier] += 1
    counters["by_group"][group] += 1
    counters["by_name"][name or normalized] += 1
    counters["by_specimen"][specimen_code or "__missing_specimen_code__"] += 1
    entries.append(
        {
            "action": action,
            "specimen_code": specimen_code,
            "specimen_site": specimen_site,
            "tier": tier,
            "group": group,
            "name": name,
            "normalized": normalized,
            "reads": reads,
            "rank_priority": _rank_priority(rank_info),
            "rank_rule": rank_info.get("rank_rule"),
            "possibility_level": rank_info.get("possibility_level"),
            "reasons": reasons,
            "protection_signals": protection,
        }
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.disable_remove_medium:
        args.remove_medium = False
    if args.disable_remove_low:
        args.remove_low = False
    if args.disable_demote_medium:
        args.demote_medium = False

    payload = _load_json(args.input_json)
    ranked_lookup = _build_rank_lookup(_load_json(args.ranked_json))

    counters = {
        "by_action": Counter(),
        "by_reason": Counter(),
        "by_tier": Counter(),
        "by_group": Counter(),
        "by_name": Counter(),
        "by_specimen": Counter(),
    }
    changed_entries = 0
    total_entries = 0
    skipped_entries_without_tiers = 0
    change_entries: list[dict[str, Any]] = []
    updated: dict[str, Any] = {}

    for patient_id, entries in payload.items():
        if not isinstance(entries, list):
            updated[str(patient_id)] = entries
            continue
        out_entries: list[Any] = []
        for entry in entries:
            if not isinstance(entry, dict):
                out_entries.append(entry)
                continue
            total_entries += 1
            raw_tiers = entry.get("pathogens_by_likelihood")
            if not isinstance(raw_tiers, dict):
                skipped_entries_without_tiers += 1
                out_entries.append(entry)
                continue

            specimen_code = str(entry.get("specimen_code") or "").strip()
            specimen_site = str(entry.get("specimen_site") or "").strip()
            tiers = _normalize_tiers(raw_tiers)
            before = json.dumps(tiers, ensure_ascii=False, sort_keys=True)
            strong_non_oral = _strong_non_oral_balf(tiers, args)

            for group in GROUP_KEYS:
                specific_genus = _specific_genus_set(tiers, group)

                for tier in ("high", "medium", "low_colonizer"):
                    keep: list[dict[str, Any]] = []
                    demote_to_low: list[dict[str, Any]] = []

                    for item in tiers[tier][group]:
                        name = str(item.get("name") or "").strip()
                        normalized = normalize_organism_name(name)
                        reads = _to_float(item.get("reads", 0))
                        if not normalized:
                            _record(
                                change_entries,
                                counters,
                                action="remove",
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier=tier,
                                group=group,
                                name=name,
                                normalized=normalized,
                                reads=reads,
                                rank_info={},
                                reasons=["invalid_name"],
                                protection=[],
                            )
                            continue

                        candidate = ranked_lookup.get((specimen_code, normalized), {})
                        rank_info = (
                            candidate.get("ranking")
                            if isinstance(candidate.get("ranking"), dict)
                            else {}
                        )
                        protection = _protection_signals(
                            group=group,
                            normalized=normalized,
                            reads=reads,
                            tier=tier,
                            rank_info=rank_info,
                            args=args,
                        )
                        protected = bool(protection)
                        reasons = _noise_reasons(
                            tier=tier,
                            group=group,
                            name=name,
                            normalized=normalized,
                            reads=reads,
                            rank_info=rank_info,
                            specific_genus=specific_genus,
                            strong_non_oral_balf=strong_non_oral,
                            specimen_site=specimen_site,
                            args=args,
                        )

                        if args.pulmonary_coarse_filter:
                            pulmonary_reasons = _pulmonary_coarse_reasons(
                                group=group,
                                tier=tier,
                                normalized=normalized,
                                reads=reads,
                            )
                            if pulmonary_reasons:
                                reasons = list(dict.fromkeys(reasons + pulmonary_reasons))
                                _record(
                                    change_entries,
                                    counters,
                                    action="remove",
                                    specimen_code=specimen_code,
                                    specimen_site=specimen_site,
                                    tier=tier,
                                    group=group,
                                    name=name,
                                    normalized=normalized,
                                    reads=reads,
                                    rank_info=rank_info,
                                    reasons=reasons,
                                    protection=protection,
                                )
                                continue

                        if _should_remove(
                            tier=tier,
                            group=group,
                            normalized=normalized,
                            reads=reads,
                            reasons=reasons,
                            protected=protected,
                            args=args,
                        ):
                            _record(
                                change_entries,
                                counters,
                                action="remove",
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier=tier,
                                group=group,
                                name=name,
                                normalized=normalized,
                                reads=reads,
                                rank_info=rank_info,
                                reasons=reasons,
                                protection=protection,
                            )
                            continue

                        if tier == "medium" and _should_demote_medium(
                            group=group,
                            reads=reads,
                            reasons=reasons,
                            protected=protected,
                            args=args,
                        ):
                            demote_to_low.append(item)
                            _record(
                                change_entries,
                                counters,
                                action="demote_medium_to_low",
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier=tier,
                                group=group,
                                name=name,
                                normalized=normalized,
                                reads=reads,
                                rank_info=rank_info,
                                reasons=reasons,
                                protection=protection,
                            )
                            continue

                        keep.append(item)

                    tiers[tier][group] = _dedupe(keep)
                    if tier == "medium" and demote_to_low:
                        tiers["low_colonizer"][group] = _dedupe(
                            tiers["low_colonizer"][group] + demote_to_low
                        )

            if "balf" in specimen_site.lower() and strong_non_oral:
                for tier, max_keep in (
                    ("medium", args.balf_medium_oral_max_keep),
                    ("low_colonizer", args.balf_low_oral_max_keep),
                ):
                    if max_keep < 0:
                        continue
                    oral_items: list[dict[str, Any]] = []
                    non_oral_items: list[dict[str, Any]] = []
                    for item in tiers[tier]["bacterial"]:
                        normalized = normalize_organism_name(item.get("name"))
                        if normalized and _is_oral_upper_airway(normalized):
                            oral_items.append(item)
                        else:
                            non_oral_items.append(item)
                    if len(oral_items) > max_keep:
                        oral_items.sort(
                            key=lambda item: (
                                ORAL_ANAEROBE_PRIORITY.get(
                                    normalize_organism_name(item.get("name")),
                                    99,
                                ),
                                -_to_float(item.get("reads", 0)),
                                str(item.get("name") or ""),
                            )
                        )
                        keep_oral = oral_items[:max_keep]
                        drop_oral = oral_items[max_keep:]
                        for item in drop_oral:
                            name = str(item.get("name") or "").strip()
                            normalized = normalize_organism_name(name)
                            reads = _to_float(item.get("reads", 0))
                            candidate = ranked_lookup.get((specimen_code, normalized), {})
                            rank_info = (
                                candidate.get("ranking")
                                if isinstance(candidate.get("ranking"), dict)
                                else {}
                            )
                            protection = _protection_signals(
                                group="bacterial",
                                normalized=normalized,
                                reads=reads,
                                tier=tier,
                                rank_info=rank_info,
                                args=args,
                            )
                            if protection:
                                keep_oral.append(item)
                                continue
                            _record(
                                change_entries,
                                counters,
                                action="remove",
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier=tier,
                                group="bacterial",
                                name=name,
                                normalized=normalized,
                                reads=reads,
                                rank_info=rank_info,
                                reasons=["balf_oral_max_keep"],
                                protection=protection,
                            )
                        tiers[tier]["bacterial"] = _dedupe(non_oral_items + keep_oral)

            for tier in TIER_KEYS:
                for group in GROUP_KEYS:
                    tiers[tier][group] = _dedupe(tiers[tier][group])

            after = json.dumps(tiers, ensure_ascii=False, sort_keys=True)
            if after != before:
                changed_entries += 1

            new_entry = dict(entry)
            new_entry["pathogens_by_likelihood"] = tiers
            new_entry["pathogens"] = _rebuild_primary_pathogens(tiers)
            out_entries.append(new_entry)

        updated[str(patient_id)] = out_entries

    _write_json(args.output_json, updated)
    summary = {
        "mode": "rank-read-pathogen-aware-no-answer-opt",
        "source_input_json": str(args.input_json),
        "ranked_json": str(args.ranked_json),
        "optimized_output_json": str(args.output_json),
        "summary_json": str(args.summary_json),
        "rules": {
            "remove_medium": args.remove_medium,
            "remove_low": args.remove_low,
            "demote_medium": args.demote_medium,
            "background_virus_max_reads": args.background_virus_max_reads,
            "medium_env_max_reads": args.medium_env_max_reads,
            "medium_low_priority_max_reads": args.medium_low_priority_max_reads,
            "low_env_max_reads": args.low_env_max_reads,
            "low_rank_remove_cutoff": args.low_rank_remove_cutoff,
            "low_rank_remove_max_reads": args.low_rank_remove_max_reads,
            "low_possibility_remove_max_reads": args.low_possibility_remove_max_reads,
            "balf_medium_oral_max_keep": args.balf_medium_oral_max_keep,
            "balf_low_oral_max_keep": args.balf_low_oral_max_keep,
            "protect_rank_priority": args.protect_rank_priority,
            "protect_read_bacterial": args.protect_read_bacterial,
            "protect_read_viral": args.protect_read_viral,
            "protect_read_fungal": args.protect_read_fungal,
            "protect_read_others": args.protect_read_others,
            "pulmonary_coarse_filter": args.pulmonary_coarse_filter,
        },
        "stats": {
            "total_entries": total_entries,
            "skipped_entries_without_tiers": skipped_entries_without_tiers,
            "changed_entries": changed_entries,
            "changed_item_count": len(change_entries),
            "by_action": dict(counters["by_action"]),
            "by_reason": dict(counters["by_reason"]),
            "by_tier": dict(counters["by_tier"]),
            "by_group": dict(counters["by_group"]),
            "top_changed_species": [
                {"name": name, "count": count}
                for name, count in counters["by_name"].most_common(80)
            ],
            "top_changed_specimens": [
                {"specimen_code": code, "count": count}
                for code, count in counters["by_specimen"].most_common(100)
            ],
        },
        "changed_entries": change_entries,
    }
    _write_json(args.summary_json, summary)

    print(
        "Rank/read/pathogen-aware OPT complete. "
        f"changed_entries={changed_entries} changed_items={len(change_entries)}"
    )
    print(f"Wrote optimized output to {args.output_json}")
    print(f"Wrote optimization summary to {args.summary_json}")


if __name__ == "__main__":
    main()

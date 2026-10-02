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
    LOW_COLONIZER_BACKGROUND_PREFIXES,
    LOW_COLONIZER_EXACT_EXCLUDE,
    LOW_COLONIZER_LOW_READ_ALLOW_PREFIXES,
    ORAL_ANAEROBE_PRIORITY,
    PNEUMOCYSTIS,
    normalize_organism_name,
    should_force_keep_candida,
)

DEFAULT_INPUT = Path("outputs") / "mngs_to_specimen_outputs.json"
DEFAULT_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs_optimized_lite.json"
DEFAULT_SUMMARY = Path("outputs") / "mngs_low_colonizer_optimization_lite.json"

GROUP_KEYS = ("bacterial", "viral", "fungal", "others")
ORAL_COLONIZER_GENERA = ("parvimonas", "prevotella", "porphyromonas", "fusobacterium")
GENUS_SANITIZER = re.compile(r"[^a-z]+")
BACKGROUND_VIRUS_PREFIXES = (
    "anelloviridae",
    "torqueteno",
    "genomoviridae",
    "cressdnaviricota",
    "circoviridae",
    "gemykibivirus",
)
GENERIC_TAXON_SUFFIXES = ("aceae", "ales")
STERILE_SITE_HINTS = (
    "blood",
    "csf",
    "cerebrospinal",
    "pleural",
    "ascites",
    "peritoneal",
    "synovial",
    "bone",
    "tissue",
    "abscess",
)
PROTECT_PREFIXES = (
    "pseudomonas",
    "klebsiella",
    "acinetobacter",
    "enterococcus",
    "staphylococcus",
    "streptococcus",
    "escherichia",
    "serratia",
    "enterobacter",
    "legionella",
    "mycoplasma",
    "candida",
    "aspergillus",
    "pneumocystisjirovecii",
    "humanbetaherpesvirus",
    "humangammaherpesvirus",
    "humanalphaherpesvirus",
)
PROTECT_READ_THRESHOLDS = {"bacterial": 20.0, "viral": 20.0, "fungal": 12.0, "others": 25.0}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rule-based OPT-lite (no Excel). "
            "Includes demote-before-delete, expanded noise removal, and protection score."
        )
    )
    parser.add_argument("--input-json", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-json", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--include-medium", action="store_true")
    parser.add_argument("--min-low-reads", type=float, default=5.0)
    parser.add_argument("--min-medium-reads", type=float, default=8.0)
    parser.add_argument("--background-max-reads", type=float, default=20.0)
    parser.add_argument("--background-virus-hard-filter", action="store_true")
    parser.add_argument("--background-virus-max-reads", type=float, default=10.0)
    parser.add_argument("--min-medium-demote-conditions", type=int, default=1)
    parser.add_argument("--min-remove-conditions", type=int, default=1)
    parser.add_argument("--protection-score-keep", type=int, default=2)
    parser.add_argument("--oral-low-reads", type=float, default=10.0)
    parser.add_argument("--viral-trace-reads", type=float, default=3.0)
    parser.add_argument("--fungal-trace-reads", type=float, default=4.0)
    parser.add_argument("--others-trace-reads", type=float, default=6.0)
    parser.add_argument("--low-rank-cutoff", type=int, default=4)
    parser.add_argument("--low-rank-max-reads", type=float, default=8.0)
    parser.add_argument("--balf-oral-max-keep", type=int, default=3)
    parser.add_argument("--balf-non-oral-strong-read", type=float, default=30.0)
    parser.add_argument("--max-low-bacterial", type=int, default=8)
    parser.add_argument("--max-low-viral", type=int, default=4)
    parser.add_argument("--max-low-fungal", type=int, default=4)
    parser.add_argument("--max-low-others", type=int, default=3)
    parser.add_argument("--generic-label-prune", action="store_true")
    parser.add_argument("--generic-label-max-reads", type=float, default=25.0)
    parser.add_argument("--balf-medium-oral-max-keep", type=int, default=2)
    parser.add_argument("--genus-shadow-prune", action="store_true")
    parser.add_argument("--shadow-read-min", type=float, default=20.0)
    parser.add_argument("--shadow-ratio", type=float, default=2.5)
    parser.add_argument("--use-project-exclude-list", action="store_true")
    return parser.parse_args(argv)


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict payload from {path}, got {type(payload)}")
    return payload


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


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


def _normalize_tiers(raw_tiers: Any) -> dict[str, dict[str, list[dict[str, Any]]]]:
    tiers = {
        "high": {group: [] for group in GROUP_KEYS},
        "medium": {group: [] for group in GROUP_KEYS},
        "low_colonizer": {group: [] for group in GROUP_KEYS},
    }
    if not isinstance(raw_tiers, dict):
        return tiers
    for tier_name in ("high", "medium", "low_colonizer"):
        level = raw_tiers.get(tier_name)
        if tier_name == "low_colonizer" and not isinstance(level, dict):
            level = raw_tiers.get("low")
        if not isinstance(level, dict):
            continue
        for group_name in GROUP_KEYS:
            values = level.get(group_name)
            if isinstance(values, list):
                tiers[tier_name][group_name] = [item for item in values if isinstance(item, dict)]
    return tiers


def _item_sort_key(item: dict[str, Any]) -> tuple[float, str]:
    return (-_to_float(item.get("reads", 0)), str(item.get("name", "")))


def _dedupe_items_keep_max(values: list[dict[str, Any]]) -> list[dict[str, Any]]:
    deduped: dict[str, dict[str, Any]] = {}
    for item in values:
        normalized = normalize_organism_name(item.get("name"))
        if not normalized:
            continue
        current = deduped.get(normalized)
        if current is None or _to_float(item.get("reads", 0)) > _to_float(current.get("reads", 0)):
            deduped[normalized] = item
    return sorted(deduped.values(), key=_item_sort_key)


def _is_force_keep_name(group_name: str, normalized: str, reads: float) -> bool:
    if group_name == "fungal":
        return normalized == PNEUMOCYSTIS or should_force_keep_candida(normalized, reads)
    if group_name == "viral":
        return normalized in FORCE_KEEP_VIRUSES
    if group_name == "bacterial":
        if normalized == "stenotrophomonasmaltophilia":
            return reads >= 50
        return normalized in FORCE_KEEP_BACTERIA
    return False


def _is_low_read_allowed(normalized: str) -> bool:
    return any(normalized.startswith(prefix) for prefix in LOW_COLONIZER_LOW_READ_ALLOW_PREFIXES)


def _is_background_prefix(normalized: str) -> bool:
    return any(normalized.startswith(prefix) for prefix in LOW_COLONIZER_BACKGROUND_PREFIXES)


def _is_background_virus_name(normalized: str) -> bool:
    if not normalized:
        return False
    return any(normalized.startswith(prefix) for prefix in BACKGROUND_VIRUS_PREFIXES)


def _is_oral_colonizer_name(normalized: str) -> bool:
    if not normalized:
        return False
    if normalized in ORAL_ANAEROBE_PRIORITY:
        return True
    return any(normalized.startswith(prefix) for prefix in ORAL_COLONIZER_GENERA)


def _extract_genus(display_name: str, normalized_name: str) -> str:
    text = str(display_name or "").strip().lower()
    if text:
        token = text.split()[0]
        genus = GENUS_SANITIZER.sub("", token)
        if genus:
            return genus
    if not normalized_name:
        return ""
    return GENUS_SANITIZER.sub("", normalized_name[:25])


def _is_generic_taxon_label(display_name: str, normalized_name: str) -> tuple[bool, str]:
    text = str(display_name or "").strip().lower()
    if not text:
        return False, ""
    text = " ".join(text.split())
    tokens = text.split(" ")
    first = GENUS_SANITIZER.sub("", tokens[0]) if tokens else ""

    if re.search(r"\bspp?\.?$", text):
        return True, first
    if first.endswith(GENERIC_TAXON_SUFFIXES):
        return True, first
    if len(tokens) == 1 and first:
        return True, first
    if normalized_name.endswith("sp") or normalized_name.endswith("spp"):
        return True, first
    return False, ""


def _build_rank_index(values: list[dict[str, Any]]) -> dict[str, int]:
    best_reads: dict[str, float] = {}
    for item in values:
        normalized = normalize_organism_name(item.get("name"))
        if not normalized:
            continue
        best_reads[normalized] = max(best_reads.get(normalized, 0.0), _to_float(item.get("reads", 0)))
    ordered = sorted(best_reads.items(), key=lambda kv: (-kv[1], kv[0]))
    return {normalized: idx + 1 for idx, (normalized, _) in enumerate(ordered)}


def _build_high_medium_genus_max(tiers: dict[str, dict[str, list[dict[str, Any]]]]) -> dict[str, float]:
    genus_max: dict[str, float] = defaultdict(float)
    for tier_name in ("high", "medium"):
        for group_name in GROUP_KEYS:
            for item in tiers[tier_name][group_name]:
                normalized = normalize_organism_name(item.get("name"))
                genus = _extract_genus(item.get("name", ""), normalized)
                if genus:
                    genus_max[genus] = max(genus_max.get(genus, 0.0), _to_float(item.get("reads", 0)))
    return genus_max


def _build_specific_genus_sets(
    tiers: dict[str, dict[str, list[dict[str, Any]]]]
) -> dict[str, set[str]]:
    by_group: dict[str, set[str]] = {group: set() for group in GROUP_KEYS}
    for tier_name in ("high", "medium", "low_colonizer"):
        for group_name in GROUP_KEYS:
            for item in tiers[tier_name][group_name]:
                normalized = normalize_organism_name(item.get("name"))
                name = str(item.get("name", "")).strip()
                if not normalized:
                    continue
                is_generic, _ = _is_generic_taxon_label(name, normalized)
                if is_generic:
                    continue
                genus = _extract_genus(name, normalized)
                if genus:
                    by_group[group_name].add(genus)
    return by_group


def _has_sterile_site(specimen_site: str) -> bool:
    site = str(specimen_site or "").strip().lower()
    return bool(site) and any(token in site for token in STERILE_SITE_HINTS)


def _protection_score(
    *,
    normalized: str,
    reads: float,
    group_name: str,
    tier_name: str,
    specimen_site: str,
    rank_index: dict[str, int],
) -> tuple[int, list[str]]:
    score = 0
    signals: list[str] = []
    if tier_name == "medium":
        score += 1
        signals.append("medium_default")
    if reads >= PROTECT_READ_THRESHOLDS.get(group_name, 20.0):
        score += 1
        signals.append("read_threshold")
    rank = rank_index.get(normalized)
    if rank is not None and rank <= 2:
        score += 1
        signals.append("top_rank")
    if _has_sterile_site(specimen_site) and group_name in {"bacterial", "fungal"}:
        score += 1
        signals.append("sterile_site")
    if any(normalized.startswith(prefix) for prefix in PROTECT_PREFIXES):
        score += 1
        signals.append("priority_prefix")
    return score, signals


def _noise_conditions(
    *,
    normalized: str,
    reads: float,
    tier_name: str,
    group_name: str,
    display_name: str,
    rank: int | None,
    high_medium_genus_max: dict[str, float],
    specific_genus_set: set[str] | None,
    args: argparse.Namespace,
) -> list[str]:
    conditions: list[str] = []
    if args.use_project_exclude_list and normalized in LOW_COLONIZER_EXACT_EXCLUDE:
        conditions.append("project_exclude_list")
    if _is_background_prefix(normalized) and reads <= args.background_max_reads:
        conditions.append("background_prefix_low_reads")
    if tier_name == "medium" and reads < args.min_medium_reads and not _is_low_read_allowed(normalized):
        conditions.append("medium_below_min_reads")
    if tier_name == "low_colonizer" and reads < args.min_low_reads and not _is_low_read_allowed(normalized):
        conditions.append("low_below_min_reads")
    if group_name == "bacterial" and _is_oral_colonizer_name(normalized) and reads <= args.oral_low_reads:
        conditions.append("oral_colonizer_low_reads")
    if group_name == "viral" and reads <= args.viral_trace_reads and not _is_low_read_allowed(normalized):
        conditions.append("viral_trace_reads")
    if (
        group_name == "fungal"
        and reads <= args.fungal_trace_reads
        and normalized != PNEUMOCYSTIS
        and not normalized.startswith("candida")
    ):
        conditions.append("fungal_trace_reads")
    if group_name == "others" and reads <= args.others_trace_reads:
        conditions.append("others_trace_reads")
    if rank is not None and rank > args.low_rank_cutoff and reads <= args.low_rank_max_reads:
        conditions.append("low_rank_low_reads")
    if args.genus_shadow_prune:
        genus = _extract_genus(display_name, normalized)
        higher = high_medium_genus_max.get(genus, 0.0) if genus else 0.0
        if genus and higher >= args.shadow_read_min and higher >= max(reads, 0.1) * args.shadow_ratio:
            conditions.append("genus_shadowed_by_higher_tier")
    if args.generic_label_prune:
        is_generic, generic_root = _is_generic_taxon_label(display_name, normalized)
        if is_generic:
            if generic_root and generic_root in (specific_genus_set or set()):
                conditions.append("generic_label_shadowed_by_species")
            if reads <= args.generic_label_max_reads and not _is_low_read_allowed(normalized):
                conditions.append("generic_label_low_reads")
    return conditions


def _record_change(
    *,
    reason_counter: Counter[str],
    group_counter: Counter[str],
    tier_counter: Counter[str],
    name_counter: Counter[str],
    specimen_counter: Counter[str],
    entries: list[dict[str, Any]],
    specimen_code: str,
    specimen_site: str,
    tier_name: str,
    group_name: str,
    name: str,
    normalized: str,
    reads: float,
    reasons: list[str],
    protection_score: int,
    protection_signals: list[str],
) -> None:
    reasons = reasons or ["unspecified"]
    for reason in reasons:
        reason_counter[reason] += 1
    group_counter[group_name] += 1
    tier_counter[tier_name] += 1
    name_counter[name or normalized] += 1
    specimen_counter[specimen_code or "__missing_specimen_code__"] += 1
    entries.append(
        {
            "specimen_code": specimen_code,
            "specimen_site": specimen_site,
            "tier": tier_name,
            "group": group_name,
            "name": name,
            "reads": reads,
            "reasons": reasons,
            "protection_score": protection_score,
            "protection_signals": protection_signals,
        }
    )


def _rebuild_primary_pathogens(tiers: dict[str, dict[str, list[dict[str, Any]]]]) -> dict[str, list[dict[str, Any]]]:
    primary = {group: [] for group in GROUP_KEYS}
    seen: dict[str, set[str]] = {group: set() for group in GROUP_KEYS}
    for tier_name in ("high", "medium"):
        for group_name in GROUP_KEYS:
            for item in tiers[tier_name][group_name]:
                normalized = normalize_organism_name(item.get("name"))
                if not normalized or normalized in seen[group_name]:
                    continue
                seen[group_name].add(normalized)
                primary[group_name].append(item)
    for group_name in GROUP_KEYS:
        primary[group_name].sort(key=_item_sort_key)
    return primary


def _has_strong_non_oral_bacteria(
    tiers: dict[str, dict[str, list[dict[str, Any]]]], args: argparse.Namespace
) -> bool:
    for tier_name in ("high", "medium"):
        for item in tiers[tier_name]["bacterial"]:
            normalized = normalize_organism_name(item.get("name"))
            if not normalized or _is_oral_colonizer_name(normalized):
                continue
            reads = _to_float(item.get("reads", 0))
            if normalized in FORCE_KEEP_BACTERIA or reads >= args.balf_non_oral_strong_read:
                return True
    return False


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    payload = _load_json(args.input_json)

    removed_by_reason: Counter[str] = Counter()
    removed_by_group: Counter[str] = Counter()
    removed_by_tier: Counter[str] = Counter()
    removed_by_name: Counter[str] = Counter()
    removed_entries: list[dict[str, Any]] = []
    per_specimen_removed_count: Counter[str] = Counter()

    demoted_by_reason: Counter[str] = Counter()
    demoted_by_group: Counter[str] = Counter()
    demoted_by_tier: Counter[str] = Counter()
    demoted_by_name: Counter[str] = Counter()
    demoted_entries: list[dict[str, Any]] = []
    per_specimen_demoted_count: Counter[str] = Counter()

    total_entries = 0
    changed_entries = 0
    skipped_entries_without_tiers = 0
    updated_payload: dict[str, Any] = {}

    for patient_id, entries in payload.items():
        if not isinstance(entries, list):
            updated_payload[patient_id] = entries
            continue
        out_entries: list[Any] = []
        for entry in entries:
            if not isinstance(entry, dict):
                out_entries.append(entry)
                continue
            total_entries += 1
            specimen_code = str(entry.get("specimen_code", "")).strip()
            specimen_site = str(entry.get("specimen_site", "")).strip()
            raw_tiers = entry.get("pathogens_by_likelihood")
            if not isinstance(raw_tiers, dict):
                skipped_entries_without_tiers += 1
                out_entries.append(entry)
                continue

            tiers = _normalize_tiers(raw_tiers)
            before = json.dumps(tiers, ensure_ascii=False, sort_keys=True)
            specific_genus_by_group = _build_specific_genus_sets(tiers)

            if args.include_medium:
                genus_max = _build_high_medium_genus_max(tiers)
                for group_name in GROUP_KEYS:
                    rank = _build_rank_index(tiers["medium"][group_name])
                    keep_medium: list[dict[str, Any]] = []
                    demote_low: list[dict[str, Any]] = []
                    for item in tiers["medium"][group_name]:
                        name = str(item.get("name", "")).strip()
                        normalized = normalize_organism_name(name)
                        reads = _to_float(item.get("reads", 0))
                        if not normalized:
                            _record_change(
                                reason_counter=removed_by_reason,
                                group_counter=removed_by_group,
                                tier_counter=removed_by_tier,
                                name_counter=removed_by_name,
                                specimen_counter=per_specimen_removed_count,
                                entries=removed_entries,
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name="medium",
                                group_name=group_name,
                                name=name,
                                normalized=normalized,
                                reads=reads,
                                reasons=["invalid_name"],
                                protection_score=0,
                                protection_signals=[],
                            )
                            continue
                        if _is_force_keep_name(group_name, normalized, reads):
                            keep_medium.append(item)
                            continue
                        if (
                            args.background_virus_hard_filter
                            and group_name == "viral"
                            and _is_background_virus_name(normalized)
                            and reads <= args.background_virus_max_reads
                        ):
                            _record_change(
                                reason_counter=removed_by_reason,
                                group_counter=removed_by_group,
                                tier_counter=removed_by_tier,
                                name_counter=removed_by_name,
                                specimen_counter=per_specimen_removed_count,
                                entries=removed_entries,
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name="medium",
                                group_name=group_name,
                                name=name,
                                normalized=normalized,
                                reads=reads,
                                reasons=["background_virus_hard_filter"],
                                protection_score=0,
                                protection_signals=[],
                            )
                            continue
                        conds = _noise_conditions(
                            normalized=normalized,
                            reads=reads,
                            tier_name="medium",
                            group_name=group_name,
                            display_name=name,
                            rank=rank.get(normalized),
                            high_medium_genus_max=genus_max,
                            specific_genus_set=specific_genus_by_group.get(group_name, set()),
                            args=args,
                        )
                        score, signals = _protection_score(
                            normalized=normalized,
                            reads=reads,
                            group_name=group_name,
                            tier_name="medium",
                            specimen_site=specimen_site,
                            rank_index=rank,
                        )
                        if (
                            args.generic_label_prune
                            and "generic_label_shadowed_by_species" in conds
                            and score < args.protection_score_keep
                        ):
                            _record_change(
                                reason_counter=removed_by_reason,
                                group_counter=removed_by_group,
                                tier_counter=removed_by_tier,
                                name_counter=removed_by_name,
                                specimen_counter=per_specimen_removed_count,
                                entries=removed_entries,
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name="medium",
                                group_name=group_name,
                                name=name,
                                normalized=normalized,
                                reads=reads,
                                reasons=["generic_label_shadowed_by_species"],
                                protection_score=score,
                                protection_signals=signals,
                            )
                            continue
                        # 1) demote first, do not directly delete medium.
                        if len(conds) >= args.min_medium_demote_conditions and score < args.protection_score_keep:
                            demote_low.append(item)
                            _record_change(
                                reason_counter=demoted_by_reason,
                                group_counter=demoted_by_group,
                                tier_counter=demoted_by_tier,
                                name_counter=demoted_by_name,
                                specimen_counter=per_specimen_demoted_count,
                                entries=demoted_entries,
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name="medium",
                                group_name=group_name,
                                name=name,
                                normalized=normalized,
                                reads=reads,
                                reasons=conds,
                                protection_score=score,
                                protection_signals=signals,
                            )
                        else:
                            keep_medium.append(item)
                    tiers["medium"][group_name] = _dedupe_items_keep_max(keep_medium)
                    tiers["low_colonizer"][group_name] = _dedupe_items_keep_max(
                        tiers["low_colonizer"][group_name] + demote_low
                    )

            if args.background_virus_hard_filter or args.generic_label_prune:
                specific_genus_by_group = _build_specific_genus_sets(tiers)
                for group_name in GROUP_KEYS:
                    rank_high = _build_rank_index(tiers["high"][group_name])
                    keep_high: list[dict[str, Any]] = []
                    for item in tiers["high"][group_name]:
                        name = str(item.get("name", "")).strip()
                        normalized = normalize_organism_name(name)
                        reads = _to_float(item.get("reads", 0))
                        if not normalized:
                            _record_change(
                                reason_counter=removed_by_reason,
                                group_counter=removed_by_group,
                                tier_counter=removed_by_tier,
                                name_counter=removed_by_name,
                                specimen_counter=per_specimen_removed_count,
                                entries=removed_entries,
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name="high",
                                group_name=group_name,
                                name=name,
                                normalized=normalized,
                                reads=reads,
                                reasons=["invalid_name"],
                                protection_score=0,
                                protection_signals=[],
                            )
                            continue
                        if _is_force_keep_name(group_name, normalized, reads):
                            keep_high.append(item)
                            continue
                        remove_reasons: list[str] = []
                        if (
                            args.background_virus_hard_filter
                            and group_name == "viral"
                            and _is_background_virus_name(normalized)
                            and reads <= args.background_virus_max_reads
                        ):
                            remove_reasons.append("background_virus_hard_filter")
                        if args.generic_label_prune:
                            is_generic, generic_root = _is_generic_taxon_label(name, normalized)
                            if generic_root and is_generic and generic_root in specific_genus_by_group.get(
                                group_name, set()
                            ):
                                remove_reasons.append("generic_label_shadowed_by_species")
                        if remove_reasons:
                            score, signals = _protection_score(
                                normalized=normalized,
                                reads=reads,
                                group_name=group_name,
                                tier_name="high",
                                specimen_site=specimen_site,
                                rank_index=rank_high,
                            )
                            if score < args.protection_score_keep:
                                _record_change(
                                    reason_counter=removed_by_reason,
                                    group_counter=removed_by_group,
                                    tier_counter=removed_by_tier,
                                    name_counter=removed_by_name,
                                    specimen_counter=per_specimen_removed_count,
                                    entries=removed_entries,
                                    specimen_code=specimen_code,
                                    specimen_site=specimen_site,
                                    tier_name="high",
                                    group_name=group_name,
                                    name=name,
                                    normalized=normalized,
                                    reads=reads,
                                    reasons=remove_reasons,
                                    protection_score=score,
                                    protection_signals=signals,
                                )
                                continue
                        keep_high.append(item)
                    tiers["high"][group_name] = _dedupe_items_keep_max(keep_high)

            if args.include_medium and "balf" in specimen_site.lower():
                strong_non_oral_medium = _has_strong_non_oral_bacteria(tiers, args)
                if strong_non_oral_medium:
                    rank_medium_bac = _build_rank_index(tiers["medium"]["bacterial"])
                    oral_medium: list[dict[str, Any]] = []
                    non_oral_medium: list[dict[str, Any]] = []
                    for item in tiers["medium"]["bacterial"]:
                        name = str(item.get("name", "")).strip()
                        normalized = normalize_organism_name(name)
                        reads = _to_float(item.get("reads", 0))
                        if normalized and _is_oral_colonizer_name(normalized):
                            score, signals = _protection_score(
                                normalized=normalized,
                                reads=reads,
                                group_name="bacterial",
                                tier_name="medium",
                                specimen_site=specimen_site,
                                rank_index=rank_medium_bac,
                            )
                            oral_medium.append(
                                {
                                    "item": item,
                                    "name": name,
                                    "normalized": normalized,
                                    "reads": reads,
                                    "score": score,
                                    "signals": signals,
                                }
                            )
                        else:
                            non_oral_medium.append(item)
                    if len(oral_medium) > args.balf_medium_oral_max_keep:
                        oral_medium.sort(
                            key=lambda c: (
                                ORAL_ANAEROBE_PRIORITY.get(c["normalized"], 99),
                                -c["reads"],
                                c["name"],
                            )
                        )
                        overflow = len(oral_medium) - args.balf_medium_oral_max_keep
                        removable = [c for c in oral_medium if c["score"] < args.protection_score_keep]
                        drop_ids = {id(c) for c in removable[:overflow]}
                        if drop_ids:
                            for dropped in [c for c in oral_medium if id(c) in drop_ids]:
                                _record_change(
                                    reason_counter=removed_by_reason,
                                    group_counter=removed_by_group,
                                    tier_counter=removed_by_tier,
                                    name_counter=removed_by_name,
                                    specimen_counter=per_specimen_removed_count,
                                    entries=removed_entries,
                                    specimen_code=specimen_code,
                                    specimen_site=specimen_site,
                                    tier_name="medium",
                                    group_name="bacterial",
                                    name=dropped["name"],
                                    normalized=dropped["normalized"],
                                    reads=dropped["reads"],
                                    reasons=["balf_medium_oral_converge"],
                                    protection_score=dropped["score"],
                                    protection_signals=dropped["signals"],
                                )
                            oral_medium = [c for c in oral_medium if id(c) not in drop_ids]
                        tiers["medium"]["bacterial"] = _dedupe_items_keep_max(
                            non_oral_medium + [c["item"] for c in oral_medium]
                        )

            specific_genus_by_group = _build_specific_genus_sets(tiers)
            genus_max = _build_high_medium_genus_max(tiers)
            strong_non_oral = _has_strong_non_oral_bacteria(tiers, args)
            max_keep_by_group = {
                "bacterial": args.max_low_bacterial,
                "viral": args.max_low_viral,
                "fungal": args.max_low_fungal,
                "others": args.max_low_others,
            }

            for group_name in GROUP_KEYS:
                rank = _build_rank_index(tiers["low_colonizer"][group_name])
                keep_candidates: list[dict[str, Any]] = []
                for item in tiers["low_colonizer"][group_name]:
                    name = str(item.get("name", "")).strip()
                    normalized = normalize_organism_name(name)
                    reads = _to_float(item.get("reads", 0))
                    if not normalized:
                        _record_change(
                            reason_counter=removed_by_reason,
                            group_counter=removed_by_group,
                            tier_counter=removed_by_tier,
                            name_counter=removed_by_name,
                            specimen_counter=per_specimen_removed_count,
                            entries=removed_entries,
                            specimen_code=specimen_code,
                            specimen_site=specimen_site,
                            tier_name="low_colonizer",
                            group_name=group_name,
                            name=name,
                            normalized=normalized,
                            reads=reads,
                            reasons=["invalid_name"],
                            protection_score=0,
                            protection_signals=[],
                        )
                        continue
                    if _is_force_keep_name(group_name, normalized, reads):
                        keep_candidates.append(
                            {
                                "item": item,
                                "name": name,
                                "normalized": normalized,
                                "reads": reads,
                                "score": 99,
                                "signals": ["force_keep"],
                                "conds": [],
                            }
                        )
                        continue
                    conds = _noise_conditions(
                        normalized=normalized,
                        reads=reads,
                        tier_name="low_colonizer",
                        group_name=group_name,
                        display_name=name,
                        rank=rank.get(normalized),
                        high_medium_genus_max=genus_max,
                        specific_genus_set=specific_genus_by_group.get(group_name, set()),
                        args=args,
                    )
                    score, signals = _protection_score(
                        normalized=normalized,
                        reads=reads,
                        group_name=group_name,
                        tier_name="low_colonizer",
                        specimen_site=specimen_site,
                        rank_index=rank,
                    )
                    # 2) delete when noise condition threshold is met.
                    if len(conds) >= args.min_remove_conditions and score < args.protection_score_keep:
                        _record_change(
                            reason_counter=removed_by_reason,
                            group_counter=removed_by_group,
                            tier_counter=removed_by_tier,
                            name_counter=removed_by_name,
                            specimen_counter=per_specimen_removed_count,
                            entries=removed_entries,
                            specimen_code=specimen_code,
                            specimen_site=specimen_site,
                            tier_name="low_colonizer",
                            group_name=group_name,
                            name=name,
                            normalized=normalized,
                            reads=reads,
                            reasons=conds,
                            protection_score=score,
                            protection_signals=signals,
                        )
                        continue
                    keep_candidates.append(
                        {
                            "item": item,
                            "name": name,
                            "normalized": normalized,
                            "reads": reads,
                            "score": score,
                            "signals": signals,
                            "conds": conds,
                        }
                    )

                # 5) BALF oral only convergence (not blanket delete).
                if group_name == "bacterial" and "balf" in specimen_site.lower() and strong_non_oral:
                    oral = [c for c in keep_candidates if _is_oral_colonizer_name(c["normalized"])]
                    non_oral = [c for c in keep_candidates if not _is_oral_colonizer_name(c["normalized"])]
                    if len(oral) > args.balf_oral_max_keep:
                        oral.sort(
                            key=lambda c: (
                                ORAL_ANAEROBE_PRIORITY.get(c["normalized"], 99),
                                -c["reads"],
                                c["name"],
                            )
                        )
                        dropped_oral = oral[args.balf_oral_max_keep :]
                        oral = oral[: args.balf_oral_max_keep]
                        for dropped in dropped_oral:
                            _record_change(
                                reason_counter=removed_by_reason,
                                group_counter=removed_by_group,
                                tier_counter=removed_by_tier,
                                name_counter=removed_by_name,
                                specimen_counter=per_specimen_removed_count,
                                entries=removed_entries,
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name="low_colonizer",
                                group_name=group_name,
                                name=dropped["name"],
                                normalized=dropped["normalized"],
                                reads=dropped["reads"],
                                reasons=(dropped.get("conds", []) + ["balf_oral_converge"]),
                                protection_score=dropped["score"],
                                protection_signals=dropped["signals"],
                            )
                        keep_candidates = non_oral + oral

                # Soft cap: remove weakest unprotected first; do not force-remove protected.
                max_keep = max_keep_by_group[group_name]
                if len(keep_candidates) > max_keep:
                    overflow = len(keep_candidates) - max_keep
                    removable = sorted(
                        (
                            c
                            for c in keep_candidates
                            if c["score"] < args.protection_score_keep
                            and len(c.get("conds", [])) >= args.min_remove_conditions
                        ),
                        key=lambda c: (c["score"], c["reads"], c["name"]),
                    )
                    drop_ids = {id(c) for c in removable[:overflow]}
                    if drop_ids:
                        for dropped in [c for c in keep_candidates if id(c) in drop_ids]:
                            _record_change(
                                reason_counter=removed_by_reason,
                                group_counter=removed_by_group,
                                tier_counter=removed_by_tier,
                                name_counter=removed_by_name,
                                specimen_counter=per_specimen_removed_count,
                                entries=removed_entries,
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name="low_colonizer",
                                group_name=group_name,
                                name=dropped["name"],
                                normalized=dropped["normalized"],
                                reads=dropped["reads"],
                                reasons=(dropped.get("conds", []) + ["low_group_soft_cap"]),
                                protection_score=dropped["score"],
                                protection_signals=dropped["signals"],
                            )
                        keep_candidates = [c for c in keep_candidates if id(c) not in drop_ids]

                tiers["low_colonizer"][group_name] = _dedupe_items_keep_max(
                    [c["item"] for c in keep_candidates]
                )

            for tier_name in ("high", "medium", "low_colonizer"):
                for group_name in GROUP_KEYS:
                    tiers[tier_name][group_name] = _dedupe_items_keep_max(tiers[tier_name][group_name])

            after = json.dumps(tiers, ensure_ascii=False, sort_keys=True)
            if after != before:
                changed_entries += 1

            new_entry = dict(entry)
            new_entry["pathogens_by_likelihood"] = tiers
            new_entry["pathogens"] = _rebuild_primary_pathogens(tiers)
            out_entries.append(new_entry)
        updated_payload[str(patient_id)] = out_entries

    _write_json(args.output_json, updated_payload)

    summary = {
        "source_input_json": str(args.input_json),
        "optimized_output_json": str(args.output_json),
        "summary_json": str(args.summary_json),
        "mode": "opt-lite-no-excel",
        "include_medium": args.include_medium,
        "rules": {
            "min_low_reads": args.min_low_reads,
            "min_medium_reads": args.min_medium_reads,
            "background_max_reads": args.background_max_reads,
            "background_virus_hard_filter": args.background_virus_hard_filter,
            "background_virus_max_reads": args.background_virus_max_reads,
            "min_medium_demote_conditions": args.min_medium_demote_conditions,
            "min_remove_conditions": args.min_remove_conditions,
            "protection_score_keep": args.protection_score_keep,
            "oral_low_reads": args.oral_low_reads,
            "viral_trace_reads": args.viral_trace_reads,
            "fungal_trace_reads": args.fungal_trace_reads,
            "others_trace_reads": args.others_trace_reads,
            "low_rank_cutoff": args.low_rank_cutoff,
            "low_rank_max_reads": args.low_rank_max_reads,
            "balf_oral_max_keep": args.balf_oral_max_keep,
            "balf_non_oral_strong_read": args.balf_non_oral_strong_read,
            "max_low_bacterial": args.max_low_bacterial,
            "max_low_viral": args.max_low_viral,
            "max_low_fungal": args.max_low_fungal,
            "max_low_others": args.max_low_others,
            "generic_label_prune": args.generic_label_prune,
            "generic_label_max_reads": args.generic_label_max_reads,
            "balf_medium_oral_max_keep": args.balf_medium_oral_max_keep,
            "genus_shadow_prune": args.genus_shadow_prune,
            "shadow_read_min": args.shadow_read_min,
            "shadow_ratio": args.shadow_ratio,
            "use_project_exclude_list": args.use_project_exclude_list,
        },
        "stats": {
            "total_entries": total_entries,
            "skipped_entries_without_tiers": skipped_entries_without_tiers,
            "changed_entries": changed_entries,
            "demoted_item_count": len(demoted_entries),
            "removed_item_count": len(removed_entries),
            "demoted_by_tier": dict(demoted_by_tier),
            "demoted_by_group": dict(demoted_by_group),
            "demoted_by_reason": dict(demoted_by_reason),
            "removed_by_tier": dict(removed_by_tier),
            "removed_by_group": dict(removed_by_group),
            "removed_by_reason": dict(removed_by_reason),
            "top_demoted_species": [{"name": n, "count": c} for n, c in demoted_by_name.most_common(50)],
            "top_removed_species": [{"name": n, "count": c} for n, c in removed_by_name.most_common(50)],
            "top_specimen_demoted_counts": [
                {"specimen_code": k, "count": v} for k, v in per_specimen_demoted_count.most_common(100)
            ],
            "top_specimen_removed_counts": [
                {"specimen_code": k, "count": v} for k, v in per_specimen_removed_count.most_common(100)
            ],
        },
        "demoted_entries": demoted_entries,
        "removed_entries": removed_entries,
    }
    _write_json(args.summary_json, summary)

    print(
        "OPT-lite complete. "
        f"changed_entries={changed_entries} demoted_items={len(demoted_entries)} "
        f"removed_items={len(removed_entries)} include_medium={args.include_medium}"
    )
    print(f"Wrote optimized output to {args.output_json}")
    print(f"Wrote optimization summary to {args.summary_json}")


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Sequence

import tools.compare_mngs_to_specimen_with_excel as cmp
from tools.run_mngs_to_specimen_agent import (
    FORCE_KEEP_BACTERIA,
    FORCE_KEEP_VIRUSES,
    ORAL_ANAEROBE_PRIORITY,
    PNEUMOCYSTIS,
    normalize_organism_name,
    should_force_keep_candida,
)

DEFAULT_INPUT = Path("outputs") / "mngs_to_specimen_outputs_main_opt_lite_med_vnext2.json"
DEFAULT_OUTPUT = Path("outputs") / "mngs_to_specimen_outputs_main_rerank_guardrail.json"
DEFAULT_SUMMARY = Path("outputs") / "mngs_two_stage_rerank_summary.json"
DEFAULT_GUARDRAIL_REPORT = Path("outputs") / "mngs_two_stage_guardrail_report.json"

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


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Two-stage re-rank for mNGS tiered outputs with optional recall guardrail. "
            "Stage1: hard filters (background viral / generic labels / BALF medium oral). "
            "Stage2: high+medium budget rerank."
        )
    )
    parser.add_argument("--input-json", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-json", type=Path, default=DEFAULT_SUMMARY)

    # Stage 1: hard filters
    parser.add_argument("--background-virus-hard-filter", action="store_true")
    parser.add_argument("--background-virus-max-reads", type=float, default=10.0)
    parser.add_argument("--generic-label-prune", action="store_true")
    parser.add_argument("--generic-label-max-reads", type=float, default=25.0)
    parser.add_argument("--generic-label-medium-block", action="store_true", default=True)
    parser.add_argument("--disable-generic-label-medium-block", action="store_true")
    parser.add_argument("--balf-medium-oral-max-keep", type=int, default=1)
    parser.add_argument("--balf-medium-oral-min-reads", type=float, default=8.0)
    parser.add_argument("--balf-non-oral-strong-read", type=float, default=30.0)

    # Stage 2: rerank budget
    parser.add_argument("--enable-stage2-rerank", action="store_true", default=True)
    parser.add_argument("--disable-stage2-rerank", action="store_true")
    parser.add_argument("--stage2-medium-first", action="store_true", default=True)
    parser.add_argument("--disable-stage2-medium-first", action="store_true")
    parser.add_argument("--hm-budget-bacterial", type=int, default=8)
    parser.add_argument("--hm-budget-viral", type=int, default=4)
    parser.add_argument("--hm-budget-fungal", type=int, default=4)
    parser.add_argument("--hm-budget-others", type=int, default=2)
    parser.add_argument("--hm-strict-budget-bacterial", type=int, default=6)
    parser.add_argument("--hm-strict-budget-viral", type=int, default=3)
    parser.add_argument("--hm-strict-budget-fungal", type=int, default=3)
    parser.add_argument("--hm-strict-budget-others", type=int, default=1)
    parser.add_argument("--strict-when-anchor", action="store_true")
    parser.add_argument("--anchor-bacterial-reads", type=float, default=40.0)
    parser.add_argument("--anchor-viral-reads", type=float, default=60.0)
    parser.add_argument("--anchor-fungal-reads", type=float, default=30.0)
    parser.add_argument("--anchor-others-reads", type=float, default=80.0)

    # Stage 2 scoring
    parser.add_argument("--score-high-bonus", type=float, default=2.0)
    parser.add_argument("--score-medium-bonus", type=float, default=0.5)
    parser.add_argument("--score-reads-weight", type=float, default=1.0)
    parser.add_argument("--score-priority-prefix-bonus", type=float, default=0.8)
    parser.add_argument("--score-sterile-site-bonus", type=float, default=0.5)
    parser.add_argument("--score-generic-penalty", type=float, default=1.0)
    parser.add_argument("--score-background-virus-penalty", type=float, default=1.2)
    parser.add_argument("--score-balf-oral-penalty", type=float, default=0.8)
    parser.add_argument("--force-keep-score", type=float, default=99.0)

    # Guardrail
    parser.add_argument("--specimen-excel", type=Path, default=None)
    parser.add_argument(
        "--guardrail-scope",
        choices=("by_specimen", "by_patient"),
        default="by_specimen",
    )
    parser.add_argument(
        "--guardrail-view",
        choices=("high", "high_medium", "all_levels"),
        default="all_levels",
    )
    parser.add_argument("--guardrail-recall-floor", type=float, default=0.85)
    parser.add_argument("--guardrail-max-drop", type=float, default=0.01)
    parser.add_argument("--guardrail-baseline-json", type=Path, default=None)
    parser.add_argument("--guardrail-rollback-on-fail", action="store_true", default=True)
    parser.add_argument("--guardrail-report-json", type=Path, default=DEFAULT_GUARDRAIL_REPORT)

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


def _has_sterile_site(specimen_site: str) -> bool:
    site = str(specimen_site or "").strip().lower()
    return bool(site) and any(token in site for token in STERILE_SITE_HINTS)


def _build_specific_genus_sets(
    tiers: dict[str, dict[str, list[dict[str, Any]]]]
) -> dict[str, set[str]]:
    by_group: dict[str, set[str]] = {group: set() for group in GROUP_KEYS}
    for tier_name in ("high", "medium", "low_colonizer"):
        for group_name in GROUP_KEYS:
            for item in tiers[tier_name][group_name]:
                normalized = normalize_organism_name(item.get("name"))
                if not normalized:
                    continue
                name = str(item.get("name", "")).strip()
                is_generic, _ = _is_generic_taxon_label(name, normalized)
                if is_generic:
                    continue
                genus = _extract_genus(name, normalized)
                if genus:
                    by_group[group_name].add(genus)
    return by_group


def _has_strong_non_oral_bacteria_hm(
    tiers: dict[str, dict[str, list[dict[str, Any]]]],
    args: argparse.Namespace,
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


def _score_candidate(
    *,
    normalized: str,
    name: str,
    reads: float,
    tier_name: str,
    group_name: str,
    force_keep: bool,
    specimen_site: str,
    strong_non_oral_balf: bool,
    args: argparse.Namespace,
) -> float:
    if force_keep:
        return args.force_keep_score
    score = math.log10(max(reads, 0.0) + 1.0) * args.score_reads_weight
    if tier_name == "high":
        score += args.score_high_bonus
    elif tier_name == "medium":
        score += args.score_medium_bonus
    if any(normalized.startswith(prefix) for prefix in PROTECT_PREFIXES):
        score += args.score_priority_prefix_bonus
    if _has_sterile_site(specimen_site) and group_name in {"bacterial", "fungal"}:
        score += args.score_sterile_site_bonus
    is_generic, _ = _is_generic_taxon_label(name, normalized)
    if is_generic:
        score -= args.score_generic_penalty
    if group_name == "viral" and _is_background_virus_name(normalized):
        score -= args.score_background_virus_penalty
    if (
        group_name == "bacterial"
        and tier_name == "medium"
        and "balf" in specimen_site.lower()
        and strong_non_oral_balf
        and _is_oral_colonizer_name(normalized)
    ):
        score -= args.score_balf_oral_penalty
    return score


def _build_output_lookup_by_level_from_payload(
    payload: dict[str, Any],
) -> tuple[dict[str, dict[str, set[str]]], dict[str, dict[str, set[str]]]]:
    by_specimen: dict[str, dict[str, set[str]]] = {
        key: defaultdict(set) for key in cmp.LIKELIHOOD_KEYS
    }
    by_patient: dict[str, dict[str, set[str]]] = {
        key: defaultdict(set) for key in cmp.LIKELIHOOD_KEYS
    }

    def extract_names(pathogens_group: Any) -> set[str]:
        names: set[str] = set()
        if not isinstance(pathogens_group, dict):
            return names
        for values in pathogens_group.values():
            if not isinstance(values, list):
                continue
            for item in values:
                if not isinstance(item, dict):
                    continue
                normalized = cmp._normalize_name(item.get("name"))
                if normalized:
                    names.add(normalized)
        return names

    for patient_id, entries in payload.items():
        if not isinstance(entries, list):
            continue
        patient_key = cmp._normalize_identifier(patient_id)
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            specimen_code = cmp._normalize_identifier(entry.get("specimen_code", ""))
            tiered = entry.get("pathogens_by_likelihood")
            if isinstance(tiered, dict):
                high_names = extract_names(tiered.get("high", {}))
                medium_names = extract_names(tiered.get("medium", {}))
                low_names = extract_names(tiered.get("low_colonizer", tiered.get("low", {})))
            else:
                legacy = extract_names(entry.get("pathogens", {}))
                high_names = legacy
                medium_names = set()
                low_names = set()

            names_by_level = {
                "high": high_names,
                "high_medium": high_names | medium_names,
                "all_levels": high_names | medium_names | low_names,
            }
            for level_key, names in names_by_level.items():
                if specimen_code:
                    by_specimen[level_key][specimen_code].update(names)
                if patient_key:
                    by_patient[level_key][patient_key].update(names)
    return (
        {level: dict(values) for level, values in by_specimen.items()},
        {level: dict(values) for level, values in by_patient.items()},
    )


def _derive_excel_by_patient_from_payload(
    payload: dict[str, Any],
    excel_by_specimen: dict[str, set[str]],
) -> dict[str, set[str]]:
    by_patient: dict[str, set[str]] = defaultdict(set)
    for patient_id, entries in payload.items():
        if not isinstance(entries, list):
            continue
        normalized_patient = cmp._normalize_identifier(patient_id)
        if not normalized_patient:
            continue
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            specimen_code = cmp._normalize_identifier(entry.get("specimen_code", ""))
            if not specimen_code:
                continue
            by_patient[normalized_patient].update(excel_by_specimen.get(specimen_code, set()))
    return dict(by_patient)


def _compute_recall_summary(
    payload: dict[str, Any],
    specimen_excel: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], str]:
    excel_by_specimen, excel_by_patient = cmp._load_excel_lookup(specimen_excel)
    output_by_specimen, output_by_patient = _build_output_lookup_by_level_from_payload(payload)
    by_specimen_report = cmp._build_comparison_by_level(output_by_specimen, excel_by_specimen)

    patient_alignment_mode = "excel_patient_id_direct"
    output_patient_keys = set(output_by_patient.get("all_levels", {}))
    excel_patient_keys = set(excel_by_patient)
    shared_patient_keys = len(output_patient_keys & excel_patient_keys)
    shared_specimen_keys = int(by_specimen_report.get("all_levels", {}).get("shared_keys", 0))
    if shared_patient_keys == 0 and shared_specimen_keys > 0:
        excel_by_patient = _derive_excel_by_patient_from_payload(payload, excel_by_specimen)
        patient_alignment_mode = "derived_from_output_specimen_codes"

    by_patient_report = cmp._build_comparison_by_level(output_by_patient, excel_by_patient)
    summary = cmp._build_summary(by_specimen_report, by_patient_report)
    return summary, by_specimen_report, by_patient_report, patient_alignment_mode


def _get_budget(group_name: str, strict: bool, args: argparse.Namespace) -> int:
    if strict:
        return {
            "bacterial": args.hm_strict_budget_bacterial,
            "viral": args.hm_strict_budget_viral,
            "fungal": args.hm_strict_budget_fungal,
            "others": args.hm_strict_budget_others,
        }[group_name]
    return {
        "bacterial": args.hm_budget_bacterial,
        "viral": args.hm_budget_viral,
        "fungal": args.hm_budget_fungal,
        "others": args.hm_budget_others,
    }[group_name]


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.disable_stage2_rerank:
        args.enable_stage2_rerank = False
    if args.disable_stage2_medium_first:
        args.stage2_medium_first = False
    if args.disable_generic_label_medium_block:
        args.generic_label_medium_block = False

    input_payload = _load_json(args.input_json)
    output_payload: dict[str, Any] = {}

    removed_by_reason: Counter[str] = Counter()
    removed_by_stage: Counter[str] = Counter()
    removed_by_tier: Counter[str] = Counter()
    removed_by_group: Counter[str] = Counter()
    removed_by_name: Counter[str] = Counter()
    per_specimen_removed: Counter[str] = Counter()
    removed_entries: list[dict[str, Any]] = []
    changed_entries = 0
    total_entries = 0
    skipped_entries_without_tiers = 0

    def record_remove(
        *,
        stage: str,
        reason: str,
        specimen_code: str,
        specimen_site: str,
        tier_name: str,
        group_name: str,
        name: str,
        normalized: str,
        reads: float,
    ) -> None:
        removed_by_reason[reason] += 1
        removed_by_stage[stage] += 1
        removed_by_tier[tier_name] += 1
        removed_by_group[group_name] += 1
        removed_by_name[name or normalized] += 1
        per_specimen_removed[specimen_code or "__missing_specimen_code__"] += 1
        removed_entries.append(
            {
                "stage": stage,
                "reason": reason,
                "specimen_code": specimen_code,
                "specimen_site": specimen_site,
                "tier": tier_name,
                "group": group_name,
                "name": name,
                "reads": reads,
            }
        )

    for patient_id, entries in input_payload.items():
        if not isinstance(entries, list):
            output_payload[str(patient_id)] = entries
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

            # Stage 1A: high/medium hard filters.
            for tier_name in ("high", "medium"):
                for group_name in GROUP_KEYS:
                    keep: list[dict[str, Any]] = []
                    for item in tiers[tier_name][group_name]:
                        name = str(item.get("name", "")).strip()
                        normalized = normalize_organism_name(name)
                        reads = _to_float(item.get("reads", 0))
                        if not normalized:
                            record_remove(
                                stage="stage1_hard_filter",
                                reason="invalid_name",
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name=tier_name,
                                group_name=group_name,
                                name=name,
                                normalized=normalized,
                                reads=reads,
                            )
                            continue
                        if _is_force_keep_name(group_name, normalized, reads):
                            keep.append(item)
                            continue

                        reasons: list[str] = []
                        if (
                            args.background_virus_hard_filter
                            and group_name == "viral"
                            and _is_background_virus_name(normalized)
                            and reads <= args.background_virus_max_reads
                        ):
                            reasons.append("background_virus_hard_filter")

                        if args.generic_label_prune:
                            is_generic, generic_root = _is_generic_taxon_label(name, normalized)
                            if is_generic:
                                if generic_root and generic_root in specific_genus_by_group.get(group_name, set()):
                                    reasons.append("generic_label_shadowed_by_species")
                                elif args.generic_label_medium_block and tier_name == "medium":
                                    reasons.append("generic_label_medium_block")
                                elif reads <= args.generic_label_max_reads:
                                    reasons.append("generic_label_low_reads")

                        if reasons:
                            for reason in reasons:
                                record_remove(
                                    stage="stage1_hard_filter",
                                    reason=reason,
                                    specimen_code=specimen_code,
                                    specimen_site=specimen_site,
                                    tier_name=tier_name,
                                    group_name=group_name,
                                    name=name,
                                    normalized=normalized,
                                    reads=reads,
                                )
                            continue
                        keep.append(item)
                    tiers[tier_name][group_name] = _dedupe_items_keep_max(keep)

            # Stage 1B: BALF medium oral convergence.
            if args.balf_medium_oral_max_keep >= 0 and "balf" in specimen_site.lower():
                strong_non_oral = _has_strong_non_oral_bacteria_hm(tiers, args)
                if strong_non_oral:
                    oral_medium: list[dict[str, Any]] = []
                    non_oral_medium: list[dict[str, Any]] = []
                    for item in tiers["medium"]["bacterial"]:
                        normalized = normalize_organism_name(item.get("name"))
                        if normalized and _is_oral_colonizer_name(normalized):
                            oral_medium.append(item)
                        else:
                            non_oral_medium.append(item)
                    if args.balf_medium_oral_min_reads > 0:
                        kept_oral: list[dict[str, Any]] = []
                        for item in oral_medium:
                            name = str(item.get("name", "")).strip()
                            normalized = normalize_organism_name(name)
                            reads = _to_float(item.get("reads", 0))
                            if _is_force_keep_name("bacterial", normalized, reads):
                                kept_oral.append(item)
                                continue
                            if reads < args.balf_medium_oral_min_reads:
                                record_remove(
                                    stage="stage1_balf_medium_oral",
                                    reason="balf_medium_oral_low_reads",
                                    specimen_code=specimen_code,
                                    specimen_site=specimen_site,
                                    tier_name="medium",
                                    group_name="bacterial",
                                    name=name,
                                    normalized=normalized,
                                    reads=reads,
                                )
                                continue
                            kept_oral.append(item)
                        oral_medium = kept_oral
                    if len(oral_medium) > args.balf_medium_oral_max_keep:
                        oral_medium.sort(
                            key=lambda x: (
                                ORAL_ANAEROBE_PRIORITY.get(normalize_organism_name(x.get("name")), 99),
                                -_to_float(x.get("reads", 0)),
                                str(x.get("name", "")),
                            )
                        )
                        keep_oral = oral_medium[: args.balf_medium_oral_max_keep]
                        drop_oral = oral_medium[args.balf_medium_oral_max_keep :]
                        for item in drop_oral:
                            name = str(item.get("name", "")).strip()
                            normalized = normalize_organism_name(name)
                            reads = _to_float(item.get("reads", 0))
                            if _is_force_keep_name("bacterial", normalized, reads):
                                keep_oral.append(item)
                                continue
                            record_remove(
                                stage="stage1_balf_medium_oral",
                                reason="balf_medium_oral_converge",
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name="medium",
                                group_name="bacterial",
                                name=name,
                                normalized=normalized,
                                reads=reads,
                            )
                        tiers["medium"]["bacterial"] = _dedupe_items_keep_max(non_oral_medium + keep_oral)

            # Stage 2: high+medium rerank budget by group.
            if args.enable_stage2_rerank:
                for group_name in GROUP_KEYS:
                    high_candidates: list[dict[str, Any]] = []
                    medium_candidates: list[dict[str, Any]] = []
                    for item in tiers["high"][group_name]:
                        name = str(item.get("name", "")).strip()
                        normalized = normalize_organism_name(name)
                        if not normalized:
                            continue
                        reads = _to_float(item.get("reads", 0))
                        high_candidates.append(
                            {
                                "tier": "high",
                                "group": group_name,
                                "item": item,
                                "name": name,
                                "normalized": normalized,
                                "reads": reads,
                                "force_keep": _is_force_keep_name(group_name, normalized, reads),
                            }
                        )
                    for item in tiers["medium"][group_name]:
                        name = str(item.get("name", "")).strip()
                        normalized = normalize_organism_name(name)
                        if not normalized:
                            continue
                        reads = _to_float(item.get("reads", 0))
                        medium_candidates.append(
                            {
                                "tier": "medium",
                                "group": group_name,
                                "item": item,
                                "name": name,
                                "normalized": normalized,
                                "reads": reads,
                                "force_keep": _is_force_keep_name(group_name, normalized, reads),
                            }
                        )

                    all_candidates = high_candidates + medium_candidates
                    if not all_candidates:
                        continue

                    strict = False
                    if args.strict_when_anchor:
                        threshold = {
                            "bacterial": args.anchor_bacterial_reads,
                            "viral": args.anchor_viral_reads,
                            "fungal": args.anchor_fungal_reads,
                            "others": args.anchor_others_reads,
                        }[group_name]
                        strict = any(c["force_keep"] or c["reads"] >= threshold for c in all_candidates)
                    budget = _get_budget(group_name, strict, args)
                    if budget < 0 or len(all_candidates) <= budget:
                        continue

                    strong_non_oral_balf = (
                        group_name == "bacterial"
                        and "balf" in specimen_site.lower()
                        and _has_strong_non_oral_bacteria_hm(tiers, args)
                    )
                    for c in all_candidates:
                        c["score"] = _score_candidate(
                            normalized=c["normalized"],
                            name=c["name"],
                            reads=c["reads"],
                            tier_name=c["tier"],
                            group_name=group_name,
                            force_keep=c["force_keep"],
                            specimen_site=specimen_site,
                            strong_non_oral_balf=strong_non_oral_balf,
                            args=args,
                        )

                    keep_set: set[int] = set()
                    if args.stage2_medium_first:
                        # Medium-first rerank: keep high first, trim medium first.
                        high_force_keep = [c for c in high_candidates if c["force_keep"]]
                        high_non_force = [c for c in high_candidates if not c["force_keep"]]
                        medium_force_keep = [c for c in medium_candidates if c["force_keep"]]
                        medium_non_force = [c for c in medium_candidates if not c["force_keep"]]

                        keep_candidates_ordered: list[dict[str, Any]] = []
                        keep_candidates_ordered.extend(high_force_keep)

                        high_slots_left = max(0, budget - len(keep_candidates_ordered))
                        if high_slots_left > 0:
                            high_non_force_sorted = sorted(
                                high_non_force,
                                key=lambda x: (-x["score"], -x["reads"], x["name"]),
                            )
                            keep_candidates_ordered.extend(high_non_force_sorted[:high_slots_left])

                        slots_after_high = max(0, budget - len(keep_candidates_ordered))
                        if slots_after_high > 0:
                            keep_candidates_ordered.extend(medium_force_keep[:slots_after_high])
                            slots_after_high = max(0, budget - len(keep_candidates_ordered))
                        if slots_after_high > 0:
                            medium_non_force_sorted = sorted(
                                medium_non_force,
                                key=lambda x: (-x["score"], -x["reads"], x["name"]),
                            )
                        keep_candidates_ordered.extend(medium_non_force_sorted[:slots_after_high])

                        keep_set = {id(c) for c in keep_candidates_ordered}
                    else:
                        force_keep_items = [c for c in all_candidates if c["force_keep"]]
                        if len(force_keep_items) >= budget:
                            keep_set = {id(c) for c in force_keep_items}
                        else:
                            slots = budget - len(force_keep_items)
                            others = sorted(
                                (c for c in all_candidates if not c["force_keep"]),
                                key=lambda x: (-x["score"], -x["reads"], x["name"]),
                            )
                            keep_set = {id(c) for c in force_keep_items + others[:slots]}

                    kept_high: list[dict[str, Any]] = []
                    kept_medium: list[dict[str, Any]] = []
                    for c in all_candidates:
                        if id(c) in keep_set:
                            if c["tier"] == "high":
                                kept_high.append(c["item"])
                            else:
                                kept_medium.append(c["item"])
                        else:
                            remove_reason = (
                                "hm_budget_rerank_medium_first"
                                if args.stage2_medium_first
                                else "hm_budget_rerank"
                            )
                            record_remove(
                                stage="stage2_budget_rerank",
                                reason=remove_reason,
                                specimen_code=specimen_code,
                                specimen_site=specimen_site,
                                tier_name=c["tier"],
                                group_name=group_name,
                                name=c["name"],
                                normalized=c["normalized"],
                                reads=c["reads"],
                            )
                    tiers["high"][group_name] = _dedupe_items_keep_max(kept_high)
                    tiers["medium"][group_name] = _dedupe_items_keep_max(kept_medium)

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

        output_payload[str(patient_id)] = out_entries

    final_payload = output_payload
    guardrail_report: dict[str, Any] = {"enabled": False}

    if args.specimen_excel:
        baseline_payload: dict[str, Any]
        baseline_source = str(args.input_json)
        if args.guardrail_baseline_json and args.guardrail_baseline_json.exists():
            baseline_payload = _load_json(args.guardrail_baseline_json)
            baseline_source = str(args.guardrail_baseline_json)
        else:
            baseline_payload = input_payload

        baseline_summary, _, _, baseline_align = _compute_recall_summary(
            baseline_payload, args.specimen_excel
        )
        candidate_summary, _, _, candidate_align = _compute_recall_summary(
            output_payload, args.specimen_excel
        )
        baseline_recall = float(
            baseline_summary[args.guardrail_scope][args.guardrail_view]["recall"]  # type: ignore[index]
        )
        candidate_recall = float(
            candidate_summary[args.guardrail_scope][args.guardrail_view]["recall"]  # type: ignore[index]
        )
        min_allowed = args.guardrail_recall_floor
        passed = candidate_recall >= min_allowed
        rolled_back = False
        if not passed and args.guardrail_rollback_on_fail:
            final_payload = baseline_payload
            rolled_back = True

        guardrail_report = {
            "enabled": True,
            "scope": args.guardrail_scope,
            "view": args.guardrail_view,
            "baseline_source": baseline_source,
            "baseline_alignment_mode": baseline_align,
            "candidate_alignment_mode": candidate_align,
            "baseline_recall": baseline_recall,
            "candidate_recall": candidate_recall,
            "recall_floor": args.guardrail_recall_floor,
            "max_drop": args.guardrail_max_drop,
            "min_allowed_recall": min_allowed,
            "passed": passed,
            "rolled_back": rolled_back,
        }
        if args.guardrail_report_json:
            _write_json(args.guardrail_report_json, guardrail_report)

    _write_json(args.output_json, final_payload)

    summary = {
        "source_input_json": str(args.input_json),
        "output_json": str(args.output_json),
        "summary_json": str(args.summary_json),
        "rules": {
            "background_virus_hard_filter": args.background_virus_hard_filter,
            "background_virus_max_reads": args.background_virus_max_reads,
            "generic_label_prune": args.generic_label_prune,
            "generic_label_max_reads": args.generic_label_max_reads,
            "generic_label_medium_block": args.generic_label_medium_block,
            "balf_medium_oral_max_keep": args.balf_medium_oral_max_keep,
            "balf_medium_oral_min_reads": args.balf_medium_oral_min_reads,
            "balf_non_oral_strong_read": args.balf_non_oral_strong_read,
            "enable_stage2_rerank": args.enable_stage2_rerank,
            "stage2_medium_first": args.stage2_medium_first,
            "strict_when_anchor": args.strict_when_anchor,
            "hm_budget_bacterial": args.hm_budget_bacterial,
            "hm_budget_viral": args.hm_budget_viral,
            "hm_budget_fungal": args.hm_budget_fungal,
            "hm_budget_others": args.hm_budget_others,
            "hm_strict_budget_bacterial": args.hm_strict_budget_bacterial,
            "hm_strict_budget_viral": args.hm_strict_budget_viral,
            "hm_strict_budget_fungal": args.hm_strict_budget_fungal,
            "hm_strict_budget_others": args.hm_strict_budget_others,
            "anchor_bacterial_reads": args.anchor_bacterial_reads,
            "anchor_viral_reads": args.anchor_viral_reads,
            "anchor_fungal_reads": args.anchor_fungal_reads,
            "anchor_others_reads": args.anchor_others_reads,
        },
        "stats": {
            "total_entries": total_entries,
            "skipped_entries_without_tiers": skipped_entries_without_tiers,
            "changed_entries": changed_entries,
            "removed_item_count": len(removed_entries),
            "removed_by_reason": dict(removed_by_reason),
            "removed_by_stage": dict(removed_by_stage),
            "removed_by_tier": dict(removed_by_tier),
            "removed_by_group": dict(removed_by_group),
            "top_removed_species": [{"name": n, "count": c} for n, c in removed_by_name.most_common(60)],
            "top_specimen_removed_counts": [
                {"specimen_code": k, "count": v} for k, v in per_specimen_removed.most_common(100)
            ],
        },
        "guardrail": guardrail_report,
        "removed_entries": removed_entries,
    }
    _write_json(args.summary_json, summary)

    print(
        "Two-stage rerank complete. "
        f"changed_entries={changed_entries} removed_items={len(removed_entries)} "
        f"guardrail_enabled={bool(args.specimen_excel)}"
    )
    if guardrail_report.get("enabled"):
        print(
            "Guardrail: "
            f"baseline_recall={guardrail_report['baseline_recall']:.4f} "
            f"candidate_recall={guardrail_report['candidate_recall']:.4f} "
            f"min_allowed={guardrail_report['min_allowed_recall']:.4f} "
            f"passed={guardrail_report['passed']} rolled_back={guardrail_report['rolled_back']}"
        )
    print(f"Wrote reranked output to {args.output_json}")
    print(f"Wrote rerank summary to {args.summary_json}")
    if args.guardrail_report_json and guardrail_report.get("enabled"):
        print(f"Wrote guardrail report to {args.guardrail_report_json}")


if __name__ == "__main__":
    main()

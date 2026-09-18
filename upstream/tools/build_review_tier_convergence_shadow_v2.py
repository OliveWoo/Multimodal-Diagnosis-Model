"""Build a species-aware review-tier convergence shadow report.

This tool is intentionally read-only with respect to patient outputs. It runs
after merge, proposes conservative high/context-to-low tier changes, and
reports both the historical genus-relaxed benchmark and a strict species-level
sensitivity analysis.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import build_context_convergence_shadow as context_v1  # noqa: E402
from tools import build_syndrome_pattern_shadow_audit as syndrome  # noqa: E402
from tools import compare_aliasclean_baseline as compare  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402


REVIEW_TIERS = ("review_high_priority", "review_context_needed")
RAG_VISIBLE_TIERS = {"picked", *REVIEW_TIERS}
DIRECT_HOSPITAL_MODULES = {"culture", "filmarray_gmtest", "molecular_microbiology"}

FIELDS = [
    *context_v1.FIELDS,
    "same_genus_distinct_species_answer",
    "answer_relationship",
    "same_genus_picked_representative",
    "same_genus_picked_direct_hospital_support",
    "shadow_rule_id",
]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build review-tier convergence shadow v2.")
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--merged-suffix", required=True)
    parser.add_argument("--answer-csv", type=Path, required=True)
    parser.add_argument("--answer-column", default="auto")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def strict_and_related_match(
    name: str,
    answers: list[compare.NameEntry],
) -> tuple[str, str, str, str, str]:
    entry = compare.name_entry(name)
    if not entry:
        return "no", "", "no", "", "none"
    for answer in answers:
        match_type = pathogen_names.name_match_type(entry.raw, answer.raw, allow_genus_relaxed=False)
        if match_type:
            return "yes", answer.raw, "no", "", match_type
    for answer in answers:
        if not entry.genus or entry.genus != answer.genus:
            continue
        if pathogen_names.is_generic_representative_label(name) or pathogen_names.is_generic_representative_label(
            answer.raw
        ):
            relationship = "generic_to_species_same_genus"
        else:
            relationship = "distinct_species_same_genus"
        return "no", "", "yes", answer.raw, relationship
    return "no", "", "no", "", "none"


def picked_items(payload: dict[str, Any]) -> list[dict[str, Any]]:
    deterministic = as_dict(payload.get("deterministic_max"))
    best = as_dict(deterministic.get("best_available_summary"))
    return [item for item in as_list(best.get("picked_pathogens")) if isinstance(item, dict)]


def item_direct_hospital_support(item: dict[str, Any]) -> bool:
    modules = {str(value) for value in as_list(item.get("support_modules"))}
    detail = as_dict(item.get("hospital_evidence_detail"))
    modules.update(str(value) for value in as_list(detail.get("support_modules")))
    return bool(modules & DIRECT_HOSPITAL_MODULES)


def stronger_same_genus_picked(name: str, picked: list[dict[str, Any]]) -> tuple[str, str]:
    key = pathogen_names.canonical_key(name)
    genus = pathogen_names.genus_name(name)
    if not key or not genus:
        return "", "no"
    for item in picked:
        picked_name = syndrome.clean_name(item)
        if (
            pathogen_names.genus_name(picked_name) == genus
            and pathogen_names.canonical_key(picked_name) != key
            and not pathogen_names.is_generic_representative_label(picked_name)
            and item_direct_hospital_support(item)
        ):
            return picked_name, "yes"
    return "", "no"


def recommend_tier(row: dict[str, Any], item: dict[str, Any]) -> tuple[str, str, str]:
    original = row["original_tier"]
    if original != "review_context_needed":
        return original, "RTCV2-HOLD-HIGH", "high-priority items are not demoted by this conservative shadow"

    v1_tier, v1_reason = context_v1.suggested_tier(row, item)
    if v1_tier == "review_low_specificity":
        group = row["clinical_ecology_group"]
        if group == "candida_or_generic_yeast":
            rule_id = "RTCV2-WEAK-CANDIDA-NONINVASIVE"
        else:
            rule_id = "RTCV2-NONSTRIATUM-CORYNE-WEAK"
        return v1_tier, rule_id, v1_reason

    dominance = str(row.get("dominance_tier") or "")
    no_exact_support = row.get("same_organism_hospital_support") != "yes"
    no_non_host_support = not str(row.get("non_host_support_modules") or "").strip()
    has_stronger_picked = row.get("same_genus_picked_direct_hospital_support") == "yes"
    if dominance.startswith(("D0", "D1")) and no_exact_support and no_non_host_support and has_stronger_picked:
        return (
            "review_low_specificity",
            "RTCV2-DISTINCT-SPECIES-STRONGER-PICKED",
            "distinct species lacks exact support and dominance while another species in the same genus is already picked with direct hospital evidence",
        )
    return original, "RTCV2-HOLD-CONTEXT", "no agreed low-risk demotion rule was met"


def output_entries(
    payload: dict[str, Any],
    demoted_keys: set[tuple[str, str]],
    included_review_tiers: set[str],
) -> list[compare.NameEntry]:
    entries: list[compare.NameEntry] = []
    for item in picked_items(payload):
        entry = compare.name_entry(syndrome.clean_name(item))
        if entry:
            entries.append(entry)
    review = as_dict(payload.get("llm_missed_candidate_review"))
    for tier in REVIEW_TIERS:
        if tier not in included_review_tiers:
            continue
        for item in as_list(review.get(tier)):
            if not isinstance(item, dict):
                continue
            name = syndrome.clean_name(item)
            if (tier, pathogen_names.canonical_key(name)) in demoted_keys:
                continue
            entry = compare.name_entry(name)
            if entry:
                entries.append(entry)
    return entries


def strict_pairwise_compare(
    output_entries: list[compare.NameEntry], answer_entries: list[compare.NameEntry]
) -> dict[str, Any]:
    unmatched_answers = set(range(len(answer_entries)))
    matched: list[dict[str, str]] = []
    output_only: list[str] = []
    for output in output_entries:
        matched_index = next(
            (
                index
                for index in sorted(unmatched_answers)
                if pathogen_names.name_match_type(
                    output.raw,
                    answer_entries[index].raw,
                    allow_genus_relaxed=False,
                )
            ),
            None,
        )
        if matched_index is None:
            output_only.append(output.raw)
            continue
        unmatched_answers.remove(matched_index)
        matched.append(
            {
                "output": output.raw,
                "answer": answer_entries[matched_index].raw,
                "match_type": pathogen_names.name_match_type(
                    output.raw,
                    answer_entries[matched_index].raw,
                    allow_genus_relaxed=False,
                ),
            }
        )
    return {
        "matched": matched,
        "output_only": output_only,
        "answer_only": [answer_entries[index].raw for index in sorted(unmatched_answers)],
    }


def metric_summary(
    patient_payloads: dict[str, dict[str, Any]],
    answers: dict[str, list[compare.NameEntry]],
    demotions: dict[str, set[tuple[str, str]]],
    *,
    strict: bool,
    included_review_tiers: set[str],
) -> dict[str, Any]:
    matched = 0
    output_count = 0
    answer_total = sum(len(items) for items in answers.values())
    answer_only: list[dict[str, str]] = []
    match_types: Counter[str] = Counter()
    for pid, answer_entries in answers.items():
        payload = patient_payloads.get(pid)
        if not payload:
            continue
        entries = output_entries(payload, demotions.get(pid, set()), included_review_tiers)
        result = strict_pairwise_compare(entries, answer_entries) if strict else compare.pairwise_compare(entries, answer_entries)
        output_count += len(entries)
        matched += len(result["matched"])
        match_types.update(item["match_type"] for item in result["matched"])
        answer_only.extend({"patient_id": pid, "organism_name": name} for name in result["answer_only"])
    return {
        "matched": matched,
        "output": output_count,
        "answer_total": answer_total,
        "precision": matched / output_count if output_count else 0.0,
        "recall": matched / answer_total if answer_total else 0.0,
        "match_type_counts": dict(match_types),
        "answer_only": answer_only,
    }


def metric_text(metric: dict[str, Any]) -> str:
    return (
        f"{metric['matched']}/{metric['output']} = {metric['precision']:.3f} precision; "
        f"{metric['matched']}/{metric['answer_total']} = {metric['recall']:.3f} recall"
    )


def write_results(path: Path, summary: dict[str, Any]) -> None:
    changes = summary["suggested_changes"]
    scope_metrics = summary["scope_metrics"]
    lines = [
        "# KH 0728 2Days Review Tier Convergence Shadow v2",
        "",
        "- No GPT calls.",
        "- Patient JSON outputs were not modified by this shadow step.",
        "- Distinct species in the same genus are not strict answer matches.",
        "",
        "## Scope metrics",
        "",
        "| Scope | Historical current | Historical shadow | Strict current | Strict shadow |",
        "|---|---|---|---|---|",
    ]
    for scope in ("picked", "picked_plus_high", "combined"):
        metrics = scope_metrics[scope]
        lines.append(
            f"| {scope} | {metric_text(metrics['historical_current'])} | "
            f"{metric_text(metrics['historical_shadow'])} | "
            f"{metric_text(metrics['strict_current'])} | "
            f"{metric_text(metrics['strict_shadow'])} |"
        )
    lines.extend(
        [
        "",
        "## Historical benchmark (genus-relaxed, retained for continuity)",
        "",
        f"- Current: {metric_text(summary['metrics']['historical_current'])}.",
        f"- Shadow: {metric_text(summary['metrics']['historical_shadow'])}.",
        "",
        "## Strict species/alias sensitivity analysis",
        "",
        f"- Current: {metric_text(summary['metrics']['strict_current'])}.",
        f"- Shadow: {metric_text(summary['metrics']['strict_shadow'])}.",
        "",
        f"## Suggested demotions ({len(changes)})",
        "",
        ]
    )
    lines.extend(
        f"- P{row['patient_id']} {row['organism_name']}: {row['shadow_rule_id']}."
        for row in changes
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    answers, _raw_answers, answer_column, answer_row_count = compare.load_answers(
        args.answer_csv, args.answer_column
    )
    rows: list[dict[str, Any]] = []
    patient_payloads: dict[str, dict[str, Any]] = {}
    demotions: dict[str, set[tuple[str, str]]] = {}
    missing_patients: list[str] = []

    for patient_dir in sorted(args.patient_root.glob("NGS_patient_*_json"), key=compare.patient_sort_key):
        pid = compare.patient_id(patient_dir)
        path = syndrome.patient_output_path(patient_dir, args.merged_suffix)
        if not path.exists():
            missing_patients.append(pid)
            continue
        payload = read_json(path)
        patient_payloads[pid] = payload
        candidates = syndrome.candidate_index(payload)
        picked = picked_items(payload)
        review = as_dict(payload.get("llm_missed_candidate_review"))
        for tier in REVIEW_TIERS:
            for item in as_list(review.get(tier)):
                if not isinstance(item, dict):
                    continue
                name = syndrome.clean_name(item)
                candidate = candidates.get(pathogen_names.canonical_key(name), {})
                base = syndrome.build_row(
                    label="review_tier_convergence_shadow_v2",
                    pid=pid,
                    tier=tier,
                    item=item,
                    candidate=candidate,
                    answers_text="",
                    answer_count=len(answers.get(pid, [])),
                )
                strict, strict_answer, related, related_answer, relationship = strict_and_related_match(
                    name, answers.get(pid, [])
                )
                representative, representative_support = stronger_same_genus_picked(name, picked)
                row = {
                    **base,
                    "original_tier": tier,
                    "strict_answer_hit": strict,
                    "strict_matched_answer": strict_answer,
                    "related_genus_match": related,
                    "related_genus_answer": related_answer,
                    "same_genus_distinct_species_answer": related_answer
                    if relationship == "distinct_species_same_genus"
                    else "",
                    "answer_relationship": relationship,
                    "same_genus_picked_representative": representative,
                    "same_genus_picked_direct_hospital_support": representative_support,
                    "nonpulmonary_hospital_evidence_only": "yes"
                    if isinstance(item.get("nonpulmonary_hospital_evidence_guardrail"), dict)
                    and item["nonpulmonary_hospital_evidence_guardrail"].get("nonpulmonary_only")
                    else "no",
                }
                negative, rescue = context_v1.evidence_flags(row, item)
                row["negative_evidence_flags"] = "; ".join(negative)
                row["rescue_evidence_flags"] = "; ".join(rescue)
                row["suggested_tier"], row["shadow_rule_id"], row["suggested_reason"] = recommend_tier(
                    row, item
                )
                if row["suggested_tier"] not in RAG_VISIBLE_TIERS:
                    demotions.setdefault(pid, set()).add((tier, pathogen_names.canonical_key(name)))
                rows.append(row)

    changes = [row for row in rows if row["suggested_tier"] != row["original_tier"]]
    no_demotions: dict[str, set[tuple[str, str]]] = {}
    scope_definitions = {
        "picked": set(),
        "picked_plus_high": {"review_high_priority"},
        "combined": {"review_high_priority", "review_context_needed"},
    }
    scope_metrics: dict[str, dict[str, Any]] = {}
    for scope, included_tiers in scope_definitions.items():
        scope_metrics[scope] = {
            "historical_current": metric_summary(
                patient_payloads,
                answers,
                no_demotions,
                strict=False,
                included_review_tiers=included_tiers,
            ),
            "historical_shadow": metric_summary(
                patient_payloads,
                answers,
                demotions,
                strict=False,
                included_review_tiers=included_tiers,
            ),
            "strict_current": metric_summary(
                patient_payloads,
                answers,
                no_demotions,
                strict=True,
                included_review_tiers=included_tiers,
            ),
            "strict_shadow": metric_summary(
                patient_payloads,
                answers,
                demotions,
                strict=True,
                included_review_tiers=included_tiers,
            ),
        }
    historical_current = scope_metrics["combined"]["historical_current"]
    historical_shadow = scope_metrics["combined"]["historical_shadow"]
    strict_current = scope_metrics["combined"]["strict_current"]
    strict_shadow = scope_metrics["combined"]["strict_shadow"]
    summary = {
        "schema_version": "review_tier_convergence_shadow_v2",
        "patient_root": str(args.patient_root),
        "merged_suffix": args.merged_suffix,
        "answer_csv": str(args.answer_csv),
        "answer_column": answer_column,
        "answer_csv_row_count": answer_row_count,
        "review_candidate_count": len(rows),
        "suggested_change_count": len(changes),
        "suggested_change_strict_answer_hit_count": sum(
            row["strict_answer_hit"] == "yes" for row in changes
        ),
        "suggested_change_same_genus_only_count": sum(
            row["related_genus_match"] == "yes" for row in changes
        ),
        "suggested_tier_counts": dict(Counter(row["suggested_tier"] for row in rows)),
        "shadow_rule_counts": dict(Counter(row["shadow_rule_id"] for row in changes)),
        "suggested_changes": [
            {
                "patient_id": row["patient_id"],
                "organism_name": row["organism_name"],
                "original_tier": row["original_tier"],
                "suggested_tier": row["suggested_tier"],
                "shadow_rule_id": row["shadow_rule_id"],
                "strict_answer_hit": row["strict_answer_hit"],
                "related_genus_answer": row["related_genus_answer"],
                "same_genus_picked_representative": row["same_genus_picked_representative"],
                "suggested_reason": row["suggested_reason"],
            }
            for row in changes
        ],
        "metrics": {
            "historical_current": historical_current,
            "historical_shadow": historical_shadow,
            "strict_current": strict_current,
            "strict_shadow": strict_shadow,
        },
        "scope_metrics": scope_metrics,
        "historical_incremental_hit_loss": historical_current["matched"] - historical_shadow["matched"],
        "strict_incremental_hit_loss": strict_current["matched"] - strict_shadow["matched"],
        "missing_patients": missing_patients,
        "formal_outputs_changed": False,
    }

    write_csv(args.output_dir / "review_tier_convergence_shadow_v2.csv", rows)
    write_csv(args.output_dir / "review_tier_convergence_shadow_v2_suggested_changes.csv", changes)
    write_json(args.output_dir / "review_tier_convergence_shadow_v2_summary.json", summary)
    write_results(args.output_dir / "RESULTS.md", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

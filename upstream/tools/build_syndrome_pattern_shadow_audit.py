"""Build a shadow audit for syndrome-pattern interpretation of merged outputs.

This report is read-only with respect to patient outputs. It inspects the
current picked/high/context/low tiers and adds patient-level pattern labels so
we can discuss rules without changing the production merge files.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import compare_aliasclean_baseline as compare  # noqa: E402
from tools import deterministic_mngs_max_scorer as scorer  # noqa: E402
from tools import mngs_common as mngs  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402
from tools import review_missed_mngs_candidates as review_rules  # noqa: E402


TIERS = ("picked", "review_high_priority", "review_context_needed", "review_low_specificity")
RAG_VISIBLE = {"picked", "review_high_priority", "review_context_needed"}
DIRECT_HOSPITAL_MODULES = {"culture", "filmarray_gmtest", "molecular_microbiology"}

DETAIL_FIELDS = [
    "label",
    "patient_id",
    "answers",
    "answer_count",
    "tier",
    "organism_name",
    "canonical_key",
    "genus",
    "classification",
    "clinical_ecology_group",
    "syndrome_pattern",
    "pattern_strength",
    "pattern_role",
    "suggested_shadow_tier",
    "suggested_shadow_reason",
    "answer_hit",
    "matched_answer",
    "match_type",
    "rank_priority",
    "reads",
    "reads_tier",
    "reads_percentile",
    "dominance_tier",
    "specimen_class",
    "same_organism_hospital_support",
    "same_genus_or_related_hospital_support",
    "support_modules",
    "non_host_support_modules",
    "review_tier_reason",
    "rationale_zh",
]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build syndrome-pattern shadow audit from merged outputs.")
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--merged-suffix", required=True)
    parser.add_argument("--answer-csv", type=Path, required=True)
    parser.add_argument("--answer-column", default="auto")
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def clean_name(item: dict[str, Any]) -> str:
    return pathogen_names.clean_display_text(item.get("organism_name") or item.get("name"))


def patient_output_path(patient_dir: Path, suffix: str) -> Path:
    pid = compare.patient_id(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{pid}_{suffix}.json"


def evidence_snapshot(item: dict[str, Any]) -> dict[str, Any]:
    return as_dict(item.get("evidence_snapshot"))


def candidate_index(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    det = as_dict(payload.get("deterministic_max"))
    output: dict[str, dict[str, Any]] = {}
    for candidate in as_list(det.get("pathogen_candidates")):
        if not isinstance(candidate, dict):
            continue
        key = pathogen_names.canonical_key(candidate.get("organism_name") or candidate.get("name"))
        if key and key not in output:
            output[key] = candidate
    return output


def best_summary(payload: dict[str, Any]) -> dict[str, Any]:
    det = as_dict(payload.get("deterministic_max"))
    return as_dict(det.get("best_available_summary"))


def merged_items(payload: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    review = as_dict(payload.get("llm_missed_candidate_review"))
    rows: list[tuple[str, dict[str, Any]]] = []
    for item in as_list(best_summary(payload).get("picked_pathogens")):
        if isinstance(item, dict):
            rows.append(("picked", item))
    for tier in ("review_high_priority", "review_context_needed", "review_low_specificity"):
        for item in as_list(review.get(tier)):
            if isinstance(item, dict):
                rows.append((tier, item))
    return rows


def item_entry(item: dict[str, Any]) -> compare.NameEntry | None:
    return compare.name_entry(clean_name(item))


def support_modules(item: dict[str, Any], candidate: dict[str, Any]) -> list[str]:
    modules: list[str] = []
    modules.extend(str(value) for value in as_list(item.get("support_modules")))
    modules.extend(str(value) for value in as_list(evidence_snapshot(item).get("support_modules")))
    modules.extend(str(value) for value in as_list(candidate.get("support_modules")))
    summary = as_dict(candidate.get("module_support_summary"))
    modules.extend(str(key) for key, value in summary.items() if value == "Support")
    output: list[str] = []
    seen: set[str] = set()
    for module in modules:
        if not module or module in seen:
            continue
        seen.add(module)
        output.append(module)
    return output


def non_host_support_modules(item: dict[str, Any], candidate: dict[str, Any]) -> list[str]:
    modules: list[str] = []
    modules.extend(str(value) for value in as_list(evidence_snapshot(item).get("non_host_support_modules")))
    modules.extend(module for module in support_modules(item, candidate) if module in DIRECT_HOSPITAL_MODULES)
    output: list[str] = []
    seen: set[str] = set()
    for module in modules:
        if not module or module == "host" or module in seen:
            continue
        seen.add(module)
        output.append(module)
    return output


def first_present(*values: Any) -> Any:
    for value in values:
        if value not in (None, "", [], {}):
            return value
    return ""


def clinical_ecology_group(name: str, classification: str) -> str:
    normalized = mngs.normalize_organism_name(name)
    if scorer.is_candida_or_generic_yeast_name(name):
        return "candida_or_generic_yeast"
    if scorer.is_rare_opportunistic_yeast_name(name):
        return "rare_opportunistic_yeast"
    if normalized in getattr(scorer, "REACTIVATION_VIRUS_NAMES", set()):
        return "herpes_reactivation_virus"
    if scorer.is_core_answer_sensitive_respiratory_virus_name(name):
        return "core_respiratory_virus"
    if review_rules.is_oral_aspiration_flora_name(name):
        return review_rules.oral_aspiration_flora_category(name)
    if scorer.is_water_environmental_gnb_name(name) or scorer.is_low_pulmonary_specificity_environmental_name(name):
        return "environmental_water_soil_low_specificity_gnb"
    if scorer.is_respiratory_commensal_flora_name(name) or scorer.is_coagulase_negative_staph_background_name(name):
        return "skin_or_airway_background"
    if normalized.startswith("corynebacterium"):
        return "potentially_pathogenic_nondiphtherial_corynebacterium"
    if any(token in normalized for token in ("klebsiella", "escherichia", "enterobacter", "serratia", "proteus", "cronobacter")):
        return "enterobacterales_or_hospital_gnb"
    if any(token in normalized for token in ("pseudomonas", "acinetobacter", "stenotrophomonas")):
        return "nonfermenter_hospital_gnb"
    cls = classification.lower()
    if cls == "viral":
        return "other_virus"
    if cls == "fungal":
        return "other_fungus"
    if cls == "bacterial":
        return "other_bacterial_or_unclassified"
    return "other_or_unclassified"


def patient_pattern_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    groups = Counter(row["clinical_ecology_group"] for row in rows)
    oral_count = sum(
        groups.get(group, 0)
        for group in (
            "strict_aspiration_anaerobe",
            "broad_oral_upper_airway_commensal",
            "atypical_oral_associated",
        )
    )
    hospital_gnb_count = groups.get("enterobacterales_or_hospital_gnb", 0) + groups.get("nonfermenter_hospital_gnb", 0)
    skin_count = groups.get("skin_or_airway_background", 0) + groups.get(
        "potentially_pathogenic_nondiphtherial_corynebacterium", 0
    )
    patterns: list[str] = []
    if oral_count >= 2 or groups.get("strict_aspiration_anaerobe", 0):
        patterns.append("aspiration_or_anaerobe_cluster")
    if hospital_gnb_count >= 2:
        patterns.append("hospital_gnb_cluster")
    if groups.get("herpes_reactivation_virus", 0):
        patterns.append("herpes_reactivation_context")
    if groups.get("candida_or_generic_yeast", 0) or groups.get("rare_opportunistic_yeast", 0):
        patterns.append("yeast_context")
    if skin_count:
        patterns.append("skin_airway_background_context")
    return {
        "patterns": patterns,
        "group_counts": dict(groups),
        "oral_count": oral_count,
        "hospital_gnb_count": hospital_gnb_count,
        "skin_airway_count": skin_count,
    }


def pattern_strength(row: dict[str, Any], patient_patterns: dict[str, Any]) -> tuple[str, str]:
    group = row["clinical_ecology_group"]
    tier = row["tier"]
    d = row["dominance_tier"]
    reads_tier = row["reads_tier"]
    has_same = row["same_organism_hospital_support"] == "yes"
    has_related = row["same_genus_or_related_hospital_support"] == "yes"
    lower_resp = row["specimen_class"] == "S2_lower_respiratory"
    high_reads = reads_tier in {"R3_high", "R4_very_high"}
    dominant = d in {"D2_moderate", "D3_dominant"}
    patterns = set(patient_patterns.get("patterns") or [])

    if tier == "picked":
        return "formal", "lead_or_selected_pathogen"
    if group == "strict_aspiration_anaerobe":
        if has_same or dominant or ("aspiration_or_anaerobe_cluster" in patterns and lower_resp and high_reads):
            return "moderate", "syndrome_representative_context"
        return "weak", "single_anaerobe_context"
    if group in {"broad_oral_upper_airway_commensal", "atypical_oral_associated"}:
        if "aspiration_or_anaerobe_cluster" in patterns and (has_same or dominant):
            return "moderate", "supporting_syndrome_member"
        return "weak", "background_oral_member"
    if group == "potentially_pathogenic_nondiphtherial_corynebacterium":
        if has_same and (high_reads or dominant):
            return "moderate", "possible_device_or_airway_pathogen"
        if high_reads or dominant:
            return "weak_to_moderate", "high-burden nondiphtherial Corynebacterium needs case-level review"
        return "weak", "background_corynebacterium_signal"
    if group == "skin_or_airway_background":
        if has_same and (high_reads or dominant):
            return "moderate", "possible_device_or_airway_pathogen"
        if has_related and high_reads:
            return "weak_to_moderate", "related_background_signal"
        return "weak", "background_skin_airway_signal"
    if group == "environmental_water_soil_low_specificity_gnb":
        if has_same or dominant:
            return "moderate", "environmental_gnb_with_convergence"
        return "weak", "environmental_low_specificity_signal"
    if group == "candida_or_generic_yeast":
        if has_same and lower_resp and high_reads:
            return "moderate", "candida_respiratory_context"
        return "weak_to_moderate", "candida_colonization_vs_infection_context"
    if group == "herpes_reactivation_virus":
        if dominant or has_same:
            return "moderate", "viral_reactivation_with_direct_support"
        return "weak_to_moderate", "viral_reactivation_context"
    if group in {"enterobacterales_or_hospital_gnb", "nonfermenter_hospital_gnb", "core_respiratory_virus"}:
        if has_same or dominant:
            return "strong", "actionable_respiratory_pathogen_context"
        return "moderate", "actionable_but_incomplete_context"
    return "unknown", "needs_rule_review"


def suggested_tier(row: dict[str, Any]) -> tuple[str, str]:
    group = row["clinical_ecology_group"]
    strength = row["pattern_strength"]
    current = row["tier"]
    if current == "picked":
        return "picked", "formal picked is not changed by this review-tier shadow audit"
    if current == "review_low_specificity":
        return "review_low_specificity", "pure convergence shadow: audit-only organisms are not promoted into RAG-visible tiers"
    if strength == "strong":
        return "review_high_priority", "strong local convergence or actionable respiratory pathogen context"
    if strength == "moderate":
        return "review_context_needed", "pattern-level evidence is plausible but incomplete"
    if strength == "weak_to_moderate":
        return "review_context_needed", "clinically important context, but should not be upgraded without direct support"
    if group in {
        "broad_oral_upper_airway_commensal",
        "atypical_oral_associated",
        "environmental_water_soil_low_specificity_gnb",
        "skin_or_airway_background",
    }:
        return "review_low_specificity", "low-specificity single-organism signal without convergent syndrome evidence"
    return "review_low_specificity", "weak pattern signal; keep in audit rather than RAG-visible by default"


def match_rows_to_answers(rows: list[dict[str, Any]], answers: list[compare.NameEntry]) -> None:
    unmatched = set(range(len(answers)))
    for row in rows:
        row["answer_hit"] = "no"
        row["matched_answer"] = ""
        row["match_type"] = ""
        entry = compare.name_entry(row["organism_name"])
        if not entry:
            continue
        best_idx = None
        best_type = ""
        for idx in list(unmatched):
            match_type = compare.names_match(entry, answers[idx])
            if not match_type:
                continue
            if match_type == "exact_or_alias":
                best_idx = idx
                best_type = match_type
                break
            if best_idx is None or (
                match_type == "approved_group_member_match" and best_type == "genus_relaxed"
            ):
                best_idx = idx
                best_type = match_type
        if best_idx is None:
            continue
        unmatched.remove(best_idx)
        row["answer_hit"] = "yes"
        row["matched_answer"] = answers[best_idx].raw
        row["match_type"] = best_type


def counter_to_rows(counter: Counter[tuple[str, ...]], fields: Sequence[str]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for key, count in counter.most_common():
        row = {field: value for field, value in zip(fields, key)}
        row["count"] = count
        output.append(row)
    return output


def safe_ratio(num: int, den: int) -> float:
    return num / den if den else 0.0


def build_row(
    *,
    label: str,
    pid: str,
    tier: str,
    item: dict[str, Any],
    candidate: dict[str, Any],
    answers_text: str,
    answer_count: int,
) -> dict[str, Any]:
    name = clean_name(item)
    snapshot = evidence_snapshot(item)
    classification = str(first_present(item.get("classification"), candidate.get("classification"), snapshot.get("classification")))
    modules = support_modules(item, candidate)
    non_host = non_host_support_modules(item, candidate)
    same_genus_related = first_present(
        item.get("same_genus_or_related_hospital_evidence"),
        item.get("related_representative_hospital_support"),
        candidate.get("related_representative_hospital_support"),
        as_dict(candidate.get("key_evidence")).get("related_representative_hospital_support"),
    )
    row = {
        "label": label,
        "patient_id": pid,
        "answers": answers_text,
        "answer_count": answer_count,
        "tier": tier,
        "organism_name": name,
        "canonical_key": pathogen_names.canonical_key(name),
        "genus": pathogen_names.genus_name(name),
        "classification": classification,
        "clinical_ecology_group": clinical_ecology_group(name, classification),
        "rank_priority": str(first_present(item.get("rank_priority"), snapshot.get("rank_priority"), candidate.get("rank_priority"))),
        "reads": str(first_present(item.get("reads"), snapshot.get("reads"), candidate.get("reads"))),
        "reads_tier": str(first_present(item.get("reads_tier"), snapshot.get("reads_tier"), candidate.get("reads_tier"))),
        "reads_percentile": str(
            first_present(item.get("reads_percentile"), snapshot.get("reads_percentile"), candidate.get("reads_percentile"))
        ),
        "dominance_tier": str(first_present(item.get("dominance_tier"), snapshot.get("dominance_tier"), candidate.get("dominance_tier"))),
        "specimen_class": str(first_present(item.get("specimen_class"), snapshot.get("specimen_class"), candidate.get("specimen_class"))),
        "same_organism_hospital_support": "yes" if set(non_host) & DIRECT_HOSPITAL_MODULES else "no",
        "same_genus_or_related_hospital_support": "yes" if same_genus_related else "no",
        "support_modules": "; ".join(modules),
        "non_host_support_modules": "; ".join(non_host),
        "review_tier_reason": str(item.get("review_tier_reason") or ""),
        "rationale_zh": str(item.get("rationale_zh") or ""),
    }
    return row


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    answers, raw_answers, answer_column, answer_row_count = compare.load_answers(args.answer_csv, args.answer_column)
    all_rows: list[dict[str, Any]] = []
    pattern_rows: list[dict[str, Any]] = []
    missing_patients: list[str] = []

    for patient_dir in sorted(args.patient_root.glob("NGS_patient_*_json"), key=compare.patient_sort_key):
        pid = compare.patient_id(patient_dir)
        path = patient_output_path(patient_dir, args.merged_suffix)
        if not path.exists():
            missing_patients.append(pid)
            continue
        payload = read_json(path)
        candidates = candidate_index(payload)
        rows: list[dict[str, Any]] = []
        for tier, item in merged_items(payload):
            name = clean_name(item)
            candidate = candidates.get(pathogen_names.canonical_key(name), {})
            rows.append(
                build_row(
                    label=args.label,
                    pid=pid,
                    tier=tier,
                    item=item,
                    candidate=candidate,
                    answers_text=raw_answers.get(pid, ""),
                    answer_count=len(answers.get(pid, [])),
                )
            )
        pattern = patient_pattern_summary(rows)
        match_rows_to_answers(rows, answers.get(pid, []))
        for row in rows:
            row["syndrome_pattern"] = "; ".join(pattern["patterns"])
            strength, role = pattern_strength(row, pattern)
            row["pattern_strength"] = strength
            row["pattern_role"] = role
            suggested, reason = suggested_tier(row)
            row["suggested_shadow_tier"] = suggested
            row["suggested_shadow_reason"] = reason
        all_rows.extend(rows)
        pattern_rows.append(
            {
                "patient_id": pid,
                "answers": raw_answers.get(pid, ""),
                "answer_count": len(answers.get(pid, [])),
                "patterns": "; ".join(pattern["patterns"]),
                "oral_count": pattern["oral_count"],
                "hospital_gnb_count": pattern["hospital_gnb_count"],
                "skin_airway_count": pattern["skin_airway_count"],
                "group_counts": json.dumps(pattern["group_counts"], ensure_ascii=False, sort_keys=True),
            }
        )

    answer_rows = [row for row in all_rows if row["answer_count"]]
    visible_rows = [row for row in answer_rows if row["tier"] in RAG_VISIBLE]
    suggested_visible_rows = [row for row in answer_rows if row["suggested_shadow_tier"] in RAG_VISIBLE]
    baseline_hits = sum(1 for row in visible_rows if row["answer_hit"] == "yes")
    suggested_hits = sum(1 for row in suggested_visible_rows if row["answer_hit"] == "yes")
    answer_total = sum(len(value) for value in answers.values())

    summary = {
        "label": args.label,
        "patient_root": str(args.patient_root),
        "merged_suffix": args.merged_suffix,
        "answer_csv": str(args.answer_csv),
        "answer_column": answer_column,
        "answer_csv_row_count": answer_row_count,
        "answer_patient_count": len(answers),
        "answer_organism_count": answer_total,
        "missing_patients": missing_patients,
        "all_rows": len(all_rows),
        "answer_patient_rows": len(answer_rows),
        "current_rag_visible": {
            "output": len(visible_rows),
            "matched": baseline_hits,
            "precision": safe_ratio(baseline_hits, len(visible_rows)),
            "recall": safe_ratio(baseline_hits, answer_total),
        },
        "shadow_rag_visible": {
            "output": len(suggested_visible_rows),
            "matched": suggested_hits,
            "precision": safe_ratio(suggested_hits, len(suggested_visible_rows)),
            "recall": safe_ratio(suggested_hits, answer_total),
        },
        "tier_counts": dict(Counter(row["tier"] for row in answer_rows)),
        "suggested_tier_counts": dict(Counter(row["suggested_shadow_tier"] for row in answer_rows)),
        "clinical_ecology_counts": dict(Counter(row["clinical_ecology_group"] for row in answer_rows)),
        "clinical_ecology_answer_hits": {
            group: {
                "hit": sum(1 for row in rows if row["answer_hit"] == "yes"),
                "total": len(rows),
            }
            for group, rows in sorted(
                ((group, [row for row in answer_rows if row["clinical_ecology_group"] == group]) for group in {row["clinical_ecology_group"] for row in answer_rows}),
                key=lambda item: item[0],
            )
        },
        "demotions_to_low": [
            {
                "patient_id": row["patient_id"],
                "organism_name": row["organism_name"],
                "tier": row["tier"],
                "suggested_shadow_tier": row["suggested_shadow_tier"],
                "clinical_ecology_group": row["clinical_ecology_group"],
                "answer_hit": row["answer_hit"],
                "matched_answer": row["matched_answer"],
                "suggested_shadow_reason": row["suggested_shadow_reason"],
            }
            for row in answer_rows
            if row["tier"] in RAG_VISIBLE and row["suggested_shadow_tier"] == "review_low_specificity"
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "syndrome_pattern_shadow_detail.csv", all_rows, DETAIL_FIELDS)
    write_csv(
        args.output_dir / "syndrome_pattern_patient_summary.csv",
        pattern_rows,
        ("patient_id", "answers", "answer_count", "patterns", "oral_count", "hospital_gnb_count", "skin_airway_count", "group_counts"),
    )
    write_csv(
        args.output_dir / "syndrome_pattern_shadow_counts.csv",
        counter_to_rows(Counter((row["clinical_ecology_group"], row["tier"], row["suggested_shadow_tier"], row["answer_hit"]) for row in answer_rows), ("clinical_ecology_group", "tier", "suggested_shadow_tier", "answer_hit")),
        ("clinical_ecology_group", "tier", "suggested_shadow_tier", "answer_hit", "count"),
    )
    write_json(args.output_dir / "syndrome_pattern_shadow_summary.json", summary)
    print(f"rows={len(all_rows)}")
    print(f"missing_patients={len(missing_patients)}")
    print(
        "current_rag_visible="
        f"{baseline_hits}/{len(visible_rows)} precision={safe_ratio(baseline_hits, len(visible_rows)):.3f} "
        f"recall={baseline_hits}/{answer_total}={safe_ratio(baseline_hits, answer_total):.3f}"
    )
    print(
        "shadow_rag_visible="
        f"{suggested_hits}/{len(suggested_visible_rows)} precision={safe_ratio(suggested_hits, len(suggested_visible_rows)):.3f} "
        f"recall={suggested_hits}/{answer_total}={safe_ratio(suggested_hits, answer_total):.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

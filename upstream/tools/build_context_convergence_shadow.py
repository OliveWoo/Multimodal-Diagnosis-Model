"""Build a KH-only shadow audit for review_context_needed convergence.

The tool never changes patient outputs. It reports evidence axes, conservative
negative/rescue flags, and a suggested tier for later rule discussion.
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

from tools import build_syndrome_pattern_shadow_audit as syndrome  # noqa: E402
from tools import compare_aliasclean_baseline as compare  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402


FIELDS = [
    "patient_id",
    "organism_name",
    "canonical_key",
    "clinical_ecology_group",
    "original_tier",
    "suggested_tier",
    "suggested_reason",
    "strict_answer_hit",
    "strict_matched_answer",
    "related_genus_match",
    "related_genus_answer",
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
    "negative_evidence_flags",
    "rescue_evidence_flags",
    "nonpulmonary_hospital_evidence_only",
    "review_tier_reason",
]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a context convergence shadow report.")
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


def evidence(item: dict[str, Any]) -> dict[str, Any]:
    value = item.get("evidence_snapshot")
    return value if isinstance(value, dict) else item


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def numeric(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def strict_and_related_match(
    name: str,
    answers: list[compare.NameEntry],
) -> tuple[str, str, str, str]:
    entry = compare.name_entry(name)
    if not entry:
        return "no", "", "no", ""
    for answer in answers:
        if pathogen_names.name_match_type(entry.raw, answer.raw, allow_genus_relaxed=False):
            return "yes", answer.raw, "no", ""
    for answer in answers:
        if entry.genus and entry.genus == answer.genus:
            return "no", "", "yes", answer.raw
    return "no", "", "no", ""


def evidence_flags(row: dict[str, Any], item: dict[str, Any]) -> tuple[list[str], list[str]]:
    negative: list[str] = []
    rescue: list[str] = []
    dominance = row["dominance_tier"]
    reads_tier = row["reads_tier"]
    same = row["same_organism_hospital_support"] == "yes"
    related = row["same_genus_or_related_hospital_support"] == "yes"
    non_host = bool(row["non_host_support_modules"])
    guardrail = item.get("nonpulmonary_hospital_evidence_guardrail")

    if not same:
        negative.append("NO_EXACT_SPECIES_PULMONARY_HOSPITAL_SUPPORT")
    if dominance.startswith(("D0", "D1")):
        negative.append("NON_DOMINANT_D0_D1")
    if reads_tier.startswith(("R0", "R1")):
        negative.append("LOW_SIGNAL_R0_R1")
    if not non_host:
        negative.append("NO_NON_HOST_CONVERGENT_SUPPORT")
    if isinstance(guardrail, dict) and guardrail.get("nonpulmonary_only"):
        negative.append("NONPULMONARY_HOSPITAL_EVIDENCE_ONLY")

    if same:
        rescue.append("EXACT_SPECIES_HOSPITAL_SUPPORT")
    if related:
        rescue.append("RELATED_GENUS_HOSPITAL_CONTEXT")
    if dominance.startswith(("D2", "D3")):
        rescue.append("D2_D3_DOMINANCE")
    if reads_tier.startswith("R4"):
        rescue.append("R4_VERY_HIGH_SIGNAL")
    elif reads_tier.startswith("R3"):
        rescue.append("R3_HIGH_SIGNAL")
    profile = evidence(item).get("candida_evidence_strength_profile")
    if isinstance(profile, dict) and profile.get("best_strength") == "strong_invasive":
        rescue.append("STRONG_INVASIVE_CANDIDA_EVIDENCE")
    return negative, rescue


def suggested_tier(row: dict[str, Any], item: dict[str, Any]) -> tuple[str, str]:
    group = row["clinical_ecology_group"]
    key = row["canonical_key"]
    dominance = row["dominance_tier"]
    reads_tier = row["reads_tier"]
    same = row["same_organism_hospital_support"] == "yes"
    related = row["same_genus_or_related_hospital_support"] == "yes"
    non_host = bool(row["non_host_support_modules"])
    negatives = set(row["negative_evidence_flags"].split("; "))
    rescues = set(row["rescue_evidence_flags"].split("; "))

    if "NONPULMONARY_HOSPITAL_EVIDENCE_ONLY" in negatives and numeric(row["reads"]) <= 0:
        return "review_omitted_with_reason", "non-pulmonary hospital-only evidence cannot support pulmonary causality"
    if "STRONG_INVASIVE_CANDIDA_EVIDENCE" in rescues:
        return "review_high_priority", "strong invasive Candida evidence is a high-priority rescue"
    if group == "candida_or_generic_yeast":
        weak_signal = reads_tier.startswith(("R0", "R1", "R2"))
        if dominance.startswith(("D0", "D1")) and weak_signal and not same and not related and not non_host:
            return "review_low_specificity", "weak non-dominant respiratory Candida without convergent hospital or invasive evidence"
    if group == "environmental_water_soil_low_specificity_gnb":
        if dominance.startswith(("D0", "D1")) and not same and not non_host:
            return "review_low_specificity", "environmental/water GNB is non-dominant and lacks same-organism hospital support"
    if group == "potentially_pathogenic_nondiphtherial_corynebacterium" and key != "corynebacteriumstriatum":
        if dominance.startswith(("D0", "D1")) and not same and not non_host:
            return "review_low_specificity", "non-striatum Corynebacterium lacks dominance and exact-species hospital support"
    return "review_context_needed", "hold current tier pending category-specific convergence discussion"


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    answers, _raw, answer_column, answer_row_count = compare.load_answers(args.answer_csv, args.answer_column)
    rows: list[dict[str, Any]] = []
    missing_patients: list[str] = []

    for patient_dir in sorted(args.patient_root.glob("NGS_patient_*_json"), key=compare.patient_sort_key):
        pid = compare.patient_id(patient_dir)
        path = syndrome.patient_output_path(patient_dir, args.merged_suffix)
        if not path.exists():
            missing_patients.append(pid)
            continue
        payload = read_json(path)
        candidates = syndrome.candidate_index(payload)
        review = payload.get("llm_missed_candidate_review") or {}
        for item in as_list(review.get("review_context_needed")):
            if not isinstance(item, dict):
                continue
            name = syndrome.clean_name(item)
            candidate = candidates.get(pathogen_names.canonical_key(name), {})
            base = syndrome.build_row(
                label="context_convergence_shadow",
                pid=pid,
                tier="review_context_needed",
                item=item,
                candidate=candidate,
                answers_text="",
                answer_count=len(answers.get(pid, [])),
            )
            strict, strict_answer, related, related_answer = strict_and_related_match(name, answers.get(pid, []))
            row = {
                **base,
                "original_tier": "review_context_needed",
                "strict_answer_hit": strict,
                "strict_matched_answer": strict_answer,
                "related_genus_match": related,
                "related_genus_answer": related_answer,
                "nonpulmonary_hospital_evidence_only": "yes"
                if isinstance(item.get("nonpulmonary_hospital_evidence_guardrail"), dict)
                and item["nonpulmonary_hospital_evidence_guardrail"].get("nonpulmonary_only")
                else "no",
            }
            negative, rescue = evidence_flags(row, item)
            row["negative_evidence_flags"] = "; ".join(negative)
            row["rescue_evidence_flags"] = "; ".join(rescue)
            row["suggested_tier"], row["suggested_reason"] = suggested_tier(row, item)
            rows.append(row)

    write_csv(args.output_dir / "context_convergence_shadow.csv", rows)
    recommendation_counts = Counter(row["suggested_tier"] for row in rows)
    category_counts = Counter(row["clinical_ecology_group"] for row in rows)
    demoted = [row for row in rows if row["suggested_tier"] != row["original_tier"]]
    summary = {
        "schema_version": "context_convergence_shadow_v1",
        "patient_root": str(args.patient_root),
        "merged_suffix": args.merged_suffix,
        "answer_csv": str(args.answer_csv),
        "answer_column": answer_column,
        "answer_csv_row_count": answer_row_count,
        "context_candidate_count": len(rows),
        "strict_answer_hit_count": sum(row["strict_answer_hit"] == "yes" for row in rows),
        "related_genus_match_count": sum(row["related_genus_match"] == "yes" for row in rows),
        "suggested_tier_counts": dict(recommendation_counts),
        "category_counts": dict(category_counts),
        "suggested_change_count": len(demoted),
        "suggested_change_strict_answer_hit_count": sum(row["strict_answer_hit"] == "yes" for row in demoted),
        "missing_patients": missing_patients,
        "formal_outputs_changed": False,
    }
    write_json(args.output_dir / "context_convergence_shadow_summary.json", summary)
    write_csv(args.output_dir / "context_convergence_shadow_suggested_changes.csv", demoted)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Build an output-only false-positive audit for merged mNGS outputs.

The report focuses on answer patients. Patients without benchmark answers are
counted in the JSON summary as unlabeled visible outputs, but their organisms
are not labeled false positives.
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

from tools import compare_aliasclean_baseline as compare  # noqa: E402
from tools import deterministic_mngs_max_scorer as scorer  # noqa: E402
from tools import mngs_common as mngs  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402


DIRECT_HOSPITAL_MODULES = {"culture", "filmarray_gmtest", "molecular_microbiology"}

DETAIL_FIELDS = [
    "label",
    "patient_id",
    "data_richness",
    "answer_count",
    "answers",
    "output_tier",
    "organism_name",
    "canonical_key",
    "genus",
    "classification",
    "fp_category",
    "primary_rule_id",
    "guardrail_rule",
    "level_rule",
    "formal_pick_exclusion_rule",
    "applied_rules",
    "review_reasons",
    "review_tier_reason",
    "basis_level",
    "mngs_signal_tier",
    "rank_priority",
    "reads",
    "reads_tier",
    "reads_percentile",
    "dominance_tier",
    "specimen_class",
    "is_lower_respiratory",
    "source_category",
    "evidence_source",
    "same_organism_hospital_support",
    "hospital_support_modules",
    "support_modules",
    "non_host_support_modules",
    "related_representative_hospital_support",
    "is_protected_pathogen",
    "is_likely_colonizer_or_background",
    "original_review_section",
    "matched_status",
]


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build merged-output output-only false-positive audit.")
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--max-suffix", required=True)
    parser.add_argument("--answer-csv", type=Path, required=True)
    parser.add_argument("--answer-column", default="auto")
    parser.add_argument("--patient-database-csv", type=Path)
    parser.add_argument("--dataset", default="")
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--summary-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    return parser.parse_args(argv)


def load_richness(path: Path | None, dataset: str) -> dict[str, str]:
    if path is None or not path.exists():
        return {}
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if dataset and str(row.get("dataset") or "") != dataset:
                continue
            pid = str(row.get("patient_id") or "").strip()
            if pid:
                out[pid] = str(row.get("data_richness") or "").strip()
    return out


def item_name(item: dict[str, Any]) -> str:
    return str(item.get("organism_name") or item.get("name") or "").strip()


def item_entry(item: dict[str, Any]) -> compare.NameEntry | None:
    return compare.name_entry(item_name(item))


def merged_output_path(patient_dir: Path, suffix: str) -> Path:
    pid = compare.patient_id(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{pid}_{suffix}.json"


def candidate_key(name: Any) -> str:
    return pathogen_names.canonical_key(name)


def deterministic_candidate_index(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    det = payload.get("deterministic_max") if isinstance(payload.get("deterministic_max"), dict) else payload
    out: dict[str, dict[str, Any]] = {}
    for candidate in det.get("pathogen_candidates") or []:
        if not isinstance(candidate, dict):
            continue
        key = candidate_key(candidate.get("organism_name") or candidate.get("name"))
        if key and key not in out:
            out[key] = candidate
    return out


def visible_outputs(payload: dict[str, Any]) -> list[dict[str, Any]]:
    det = payload.get("deterministic_max") if isinstance(payload.get("deterministic_max"), dict) else payload
    review = payload.get("llm_missed_candidate_review") if isinstance(payload.get("llm_missed_candidate_review"), dict) else {}
    best = det.get("best_available_summary") if isinstance(det, dict) else {}
    outputs: list[dict[str, Any]] = []
    for item in best.get("picked_pathogens") or []:
        if isinstance(item, dict):
            outputs.append({"output_tier": "picked", "item": item})
    for tier in ("review_high_priority", "review_context_needed"):
        for item in review.get(tier) or []:
            if isinstance(item, dict):
                outputs.append({"output_tier": tier, "item": item})
    return outputs


def match_output_indices(
    outputs: list[dict[str, Any]],
    answer_entries: list[compare.NameEntry],
) -> tuple[set[int], list[dict[str, str]], list[str]]:
    unmatched_answers = set(range(len(answer_entries)))
    matched_indices: set[int] = set()
    matched: list[dict[str, str]] = []
    for idx, output in enumerate(outputs):
        entry = item_entry(output["item"])
        if entry is None:
            continue
        best_idx = None
        best_type = ""
        for answer_idx in list(unmatched_answers):
            match_type = compare.names_match(entry, answer_entries[answer_idx])
            if not match_type:
                continue
            if match_type == "exact_or_alias":
                best_idx = answer_idx
                best_type = match_type
                break
            if best_idx is None or (
                match_type == "approved_group_member_match" and best_type == "genus_relaxed"
            ):
                best_idx = answer_idx
                best_type = match_type
        if best_idx is None:
            continue
        unmatched_answers.remove(best_idx)
        matched_indices.add(idx)
        matched.append({"output": entry.raw, "answer": answer_entries[best_idx].raw, "match_type": best_type})
    return matched_indices, matched, [answer_entries[idx].raw for idx in sorted(unmatched_answers)]


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def joined(values: Sequence[Any]) -> str:
    return "; ".join(str(value) for value in values if str(value) != "")


def evidence_snapshot(item: dict[str, Any]) -> dict[str, Any]:
    snapshot = item.get("evidence_snapshot")
    return snapshot if isinstance(snapshot, dict) else {}


def key_evidence(candidate: dict[str, Any]) -> dict[str, Any]:
    value = candidate.get("key_evidence")
    return value if isinstance(value, dict) else {}


def first_present(*values: Any) -> Any:
    for value in values:
        if value not in (None, "", [], {}):
            return value
    return ""


def support_modules(item: dict[str, Any], candidate: dict[str, Any]) -> list[str]:
    snapshot = evidence_snapshot(item)
    modules: list[str] = []
    modules.extend(str(value) for value in as_list(item.get("support_modules")))
    modules.extend(str(value) for value in as_list(snapshot.get("support_modules")))
    summary = candidate.get("module_support_summary") if isinstance(candidate, dict) else {}
    if isinstance(summary, dict):
        modules.extend(str(key) for key, value in summary.items() if value == "Support")
    out: list[str] = []
    seen: set[str] = set()
    for module in modules:
        if not module or module in seen:
            continue
        seen.add(module)
        out.append(module)
    return out


def non_host_support_modules(item: dict[str, Any], candidate: dict[str, Any]) -> list[str]:
    snapshot = evidence_snapshot(item)
    modules = [str(value) for value in as_list(snapshot.get("non_host_support_modules"))]
    modules.extend(module for module in support_modules(item, candidate) if module in DIRECT_HOSPITAL_MODULES)
    out: list[str] = []
    seen: set[str] = set()
    for module in modules:
        if not module or module in seen:
            continue
        seen.add(module)
        out.append(module)
    return out


def rule_values(item: dict[str, Any], candidate: dict[str, Any]) -> dict[str, str]:
    evidence = key_evidence(candidate)
    snapshot = evidence_snapshot(item)
    formal = candidate.get("formal_pick_exclusion_rule")
    formal_rule = ""
    if isinstance(formal, dict):
        formal_rule = str(formal.get("rule_id") or formal.get("rule") or "")
    applied = [str(value) for value in as_list(candidate.get("applied_rules"))]
    review_reasons = [str(value) for value in as_list(snapshot.get("review_reasons"))]
    guardrail = str(first_present(evidence.get("guardrail_rule"), snapshot.get("guardrail_rule")))
    level_rule = str(first_present(evidence.get("level_rule"), snapshot.get("level_rule")))
    primary = str(
        first_present(
            guardrail,
            formal_rule,
            review_reasons[0] if review_reasons else "",
            level_rule,
            applied[0] if applied else "",
            item.get("review_tier_reason"),
        )
    )
    return {
        "primary_rule_id": primary,
        "guardrail_rule": guardrail,
        "level_rule": level_rule,
        "formal_pick_exclusion_rule": formal_rule,
        "applied_rules": joined(applied),
        "review_reasons": joined(review_reasons),
    }


def fp_category(name: str, item: dict[str, Any], candidate: dict[str, Any]) -> str:
    classification = str(first_present(item.get("classification"), candidate.get("classification"), "")).lower()
    normalized = mngs.normalize_organism_name(name)
    if scorer.is_candida_or_generic_yeast_name(name):
        return "Candida/generic yeast"
    if scorer.is_rare_opportunistic_yeast_name(name):
        return "rare/opportunistic yeast"
    if normalized in getattr(scorer, "REACTIVATION_VIRUS_NAMES", set()):
        return "herpes/reactivation virus"
    if scorer.is_water_environmental_gnb_name(name) or scorer.is_low_pulmonary_specificity_environmental_name(name):
        return "water/environmental GNB"
    if mngs.is_oral_upper_airway_flora_name(name):
        return "oral/upper-airway flora"
    if scorer.is_respiratory_commensal_flora_name(name) or scorer.is_coagulase_negative_staph_background_name(name):
        return "skin/respiratory commensal"
    if classification == "fungal":
        return "other fungus/mold"
    if classification == "viral":
        return "other virus"
    if classification == "bacterial":
        return "other bacterial"
    return "other/unknown"


def audit_row(
    *,
    label: str,
    pid: str,
    data_richness: str,
    answer_text: str,
    answer_count: int,
    output_tier: str,
    item: dict[str, Any],
    candidate: dict[str, Any],
) -> dict[str, Any]:
    name = item_name(item)
    snapshot = evidence_snapshot(item)
    rules = rule_values(item, candidate)
    modules = support_modules(item, candidate)
    non_host = non_host_support_modules(item, candidate)
    specimen_class = str(first_present(item.get("specimen_class"), snapshot.get("specimen_class"), candidate.get("specimen_class")))
    related = first_present(
        item.get("related_representative_hospital_support"),
        candidate.get("related_representative_hospital_support"),
        key_evidence(candidate).get("related_representative_hospital_support"),
    )
    if isinstance(related, (list, dict)):
        related_text = json.dumps(related, ensure_ascii=False)
    else:
        related_text = str(related or "")
    classification = str(first_present(item.get("classification"), candidate.get("classification")))
    return {
        "label": label,
        "patient_id": pid,
        "data_richness": data_richness,
        "answer_count": answer_count,
        "answers": answer_text,
        "output_tier": output_tier,
        "organism_name": name,
        "canonical_key": candidate_key(name),
        "genus": pathogen_names.genus_name(name),
        "classification": classification,
        "fp_category": fp_category(name, item, candidate),
        **rules,
        "review_tier_reason": str(item.get("review_tier_reason") or ""),
        "basis_level": str(first_present(item.get("basis_level"), candidate.get("integrated_causative_level"))),
        "mngs_signal_tier": str(first_present(item.get("mngs_signal_tier"), candidate.get("mngs_signal_tier"))),
        "rank_priority": str(first_present(item.get("rank_priority"), snapshot.get("rank_priority"), candidate.get("rank_priority"))),
        "reads": str(first_present(item.get("reads"), snapshot.get("reads"), candidate.get("reads"))),
        "reads_tier": str(first_present(item.get("reads_tier"), snapshot.get("reads_tier"), candidate.get("reads_tier"))),
        "reads_percentile": str(
            first_present(item.get("reads_percentile"), snapshot.get("reads_percentile"), candidate.get("reads_percentile"))
        ),
        "dominance_tier": str(first_present(item.get("dominance_tier"), snapshot.get("dominance_tier"), candidate.get("dominance_tier"))),
        "specimen_class": specimen_class,
        "is_lower_respiratory": str(specimen_class == "S2_lower_respiratory"),
        "source_category": str(first_present(item.get("source_category"), snapshot.get("source_category"), candidate.get("source_category"))),
        "evidence_source": str(first_present(item.get("evidence_source"), candidate.get("evidence_source"))),
        "same_organism_hospital_support": str(bool(set(non_host) & DIRECT_HOSPITAL_MODULES)),
        "hospital_support_modules": joined(module for module in non_host if module in DIRECT_HOSPITAL_MODULES),
        "support_modules": joined(modules),
        "non_host_support_modules": joined(non_host),
        "related_representative_hospital_support": related_text,
        "is_protected_pathogen": str(bool(first_present(item.get("is_protected_pathogen"), candidate.get("is_protected_pathogen")))),
        "is_likely_colonizer_or_background": str(bool(candidate.get("is_likely_colonizer_or_background"))),
        "original_review_section": str(item.get("original_review_section") or ""),
        "matched_status": "output_only",
    }


def counter_rows(rows: list[dict[str, Any]], key: str) -> list[dict[str, Any]]:
    counter = Counter(str(row.get(key) or "Unknown") for row in rows)
    return [{"grouping": key, "value": value, "count": count} for value, count in counter.most_common()]


def write_detail_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DETAIL_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["grouping", "value", "count"])
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    answers, raw_answers, answer_column, answer_row_count = compare.load_answers(args.answer_csv, args.answer_column)
    richness = load_richness(args.patient_database_csv, args.dataset)
    rows: list[dict[str, Any]] = []
    missing_outputs: list[str] = []
    unlabeled_visible_output_count = 0
    matched_total = 0
    output_total_answer_patients = 0
    answer_total = 0

    for patient_dir in sorted(args.patient_root.glob("NGS_patient_*_json"), key=compare.patient_sort_key):
        pid = compare.patient_id(patient_dir)
        path = merged_output_path(patient_dir, args.max_suffix)
        if not path.exists():
            missing_outputs.append(pid)
            continue
        payload = read_json(path)
        candidate_index = deterministic_candidate_index(payload)
        outputs = visible_outputs(payload)
        answer_entries = answers.get(pid, [])
        if not answer_entries:
            unlabeled_visible_output_count += len(outputs)
            continue
        matched_indices, matched, _missing_answers = match_output_indices(outputs, answer_entries)
        matched_total += len(matched)
        output_total_answer_patients += len(outputs)
        answer_total += len(answer_entries)
        for idx, output in enumerate(outputs):
            if idx in matched_indices:
                continue
            item = output["item"]
            name = item_name(item)
            candidate = candidate_index.get(candidate_key(name), {})
            rows.append(
                audit_row(
                    label=args.label,
                    pid=pid,
                    data_richness=richness.get(pid, ""),
                    answer_text=raw_answers.get(pid, ""),
                    answer_count=len(answer_entries),
                    output_tier=output["output_tier"],
                    item=item,
                    candidate=candidate,
                )
            )

    summary_rows: list[dict[str, Any]] = []
    for key in [
        "output_tier",
        "fp_category",
        "data_richness",
        "primary_rule_id",
        "classification",
        "same_organism_hospital_support",
        "is_lower_respiratory",
        "reads_tier",
        "dominance_tier",
        "specimen_class",
        "source_category",
        "hospital_support_modules",
    ]:
        summary_rows.extend(counter_rows(rows, key))

    write_detail_csv(args.output_csv, rows)
    write_summary_csv(args.summary_csv, summary_rows)
    summary = {
        "label": args.label,
        "patient_root": str(args.patient_root),
        "max_suffix": args.max_suffix,
        "answer_csv": str(args.answer_csv),
        "answer_column": answer_column,
        "answer_csv_row_count": answer_row_count,
        "answer_patient_count": len(answers),
        "answer_organism_count": answer_total,
        "visible_output_count_answer_patients": output_total_answer_patients,
        "matched_output_count_answer_patients": matched_total,
        "output_only_count_answer_patients": len(rows),
        "precision": matched_total / output_total_answer_patients if output_total_answer_patients else 0.0,
        "recall": matched_total / answer_total if answer_total else 0.0,
        "unlabeled_visible_output_count_no_answer_patients": unlabeled_visible_output_count,
        "missing_output_patient_ids": missing_outputs,
        "detail_csv": str(args.output_csv),
        "summary_csv": str(args.summary_csv),
        "distributions": {
            key: dict(Counter(str(row.get(key) or "Unknown") for row in rows).most_common())
            for key in [
                "output_tier",
                "fp_category",
                "data_richness",
                "primary_rule_id",
                "classification",
                "same_organism_hospital_support",
                "is_lower_respiratory",
                "reads_tier",
                "dominance_tier",
            ]
        },
    }
    write_json(args.output_json, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

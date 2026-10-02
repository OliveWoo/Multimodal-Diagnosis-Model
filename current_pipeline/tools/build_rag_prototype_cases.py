"""Build prototype RAG inputs for benchmark answers missed by formal picked output.

This script does not perform literature search. It packages local evidence and
a fixed evidence-card schema so downstream RAG can judge cases consistently.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import compare_aliasclean_baseline as compare  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402

RAG_VISIBLE_TIERS = {"review_high_priority", "review_context_needed"}
REVIEW_TIERS = (
    "review_high_priority",
    "review_context_needed",
    "review_low_specificity",
    "omitted_with_reason",
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build prototype RAG evidence-card inputs.")
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


def write_jsonl(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def patient_output_path(patient_dir: Path, suffix: str) -> Path:
    pid = compare.patient_id(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{pid}_{suffix}.json"


def final_summary_path(patient_dir: Path) -> Path:
    pid = compare.patient_id(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{pid}_final_summary_with_filmarray_deterministic.json"


def item_name(item: dict[str, Any]) -> str:
    return str(item.get("organism_name") or item.get("name") or "").strip()


def item_entry(item: dict[str, Any]) -> compare.NameEntry | None:
    return compare.name_entry(item_name(item))


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def joined(values: Sequence[Any]) -> str:
    return "; ".join(str(value) for value in values if str(value).strip())


def deterministic_payload(merged: dict[str, Any]) -> dict[str, Any]:
    value = merged.get("deterministic_max")
    return value if isinstance(value, dict) else merged


def review_payload(merged: dict[str, Any]) -> dict[str, Any]:
    value = merged.get("llm_missed_candidate_review")
    return value if isinstance(value, dict) else {}


def picked_items(det: dict[str, Any]) -> list[dict[str, Any]]:
    best = det.get("best_available_summary") if isinstance(det, dict) else {}
    return [item for item in (best.get("picked_pathogens") or []) if isinstance(item, dict)]


def review_items(review: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    out: list[tuple[str, dict[str, Any]]] = []
    for tier in REVIEW_TIERS:
        for item in review.get(tier) or []:
            if isinstance(item, dict):
                out.append((tier, item))
    return out


def match_in_items(answer: compare.NameEntry, items: Sequence[dict[str, Any]]) -> tuple[dict[str, Any], str] | None:
    for item in items:
        entry = item_entry(item)
        if entry is None:
            continue
        match_type = compare.names_match(entry, answer)
        if match_type:
            return item, match_type
    return None


def match_in_review(answer: compare.NameEntry, items: Sequence[tuple[str, dict[str, Any]]]) -> tuple[str, dict[str, Any], str] | None:
    for tier, item in items:
        entry = item_entry(item)
        if entry is None:
            continue
        match_type = compare.names_match(entry, answer)
        if match_type:
            return tier, item, match_type
    return None


def candidate_index(det: dict[str, Any]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for item in det.get("pathogen_candidates") or []:
        if not isinstance(item, dict):
            continue
        key = pathogen_names.canonical_key(item.get("organism_name") or item.get("name"))
        if key and key not in out:
            out[key] = item
    return out


def find_candidate(answer: compare.NameEntry, det: dict[str, Any], matched_review_item: dict[str, Any] | None) -> dict[str, Any]:
    index = candidate_index(det)
    if matched_review_item:
        key = pathogen_names.canonical_key(item_name(matched_review_item))
        if key in index:
            return index[key]
    if answer.canonical in index:
        return index[answer.canonical]
    for item in index.values():
        entry = compare.name_entry(item.get("organism_name"))
        if entry and compare.names_match(entry, answer):
            return item
    return {}


def matching_hospital_evidence(summary: dict[str, Any], answer: compare.NameEntry) -> dict[str, list[dict[str, Any]]]:
    exact: list[dict[str, Any]] = []
    same_genus: list[dict[str, Any]] = []
    for item in summary.get("hospital_organism_evidence") or []:
        if not isinstance(item, dict):
            continue
        entry = compare.name_entry(item.get("organism_name"))
        if entry is None:
            continue
        if compare.names_match(entry, answer):
            exact.append(item)
        elif entry.genus and answer.genus and entry.genus == answer.genus:
            same_genus.append(item)
    return {"same_organism_or_alias": exact, "same_genus_or_related": same_genus}


def support_modules(candidate: dict[str, Any], review_item: dict[str, Any] | None) -> list[str]:
    modules: list[str] = []
    key_evidence = candidate.get("key_evidence") if isinstance(candidate.get("key_evidence"), dict) else {}
    modules.extend(str(value) for value in as_list(key_evidence.get("support_modules")))
    summary = candidate.get("module_support_summary") if isinstance(candidate.get("module_support_summary"), dict) else {}
    modules.extend(str(key) for key, value in summary.items() if value == "Support")
    if review_item:
        snapshot = review_item.get("evidence_snapshot") if isinstance(review_item.get("evidence_snapshot"), dict) else {}
        modules.extend(str(value) for value in as_list(snapshot.get("support_modules")))
    out: list[str] = []
    seen: set[str] = set()
    for module in modules:
        if module and module not in seen:
            out.append(module)
            seen.add(module)
    return out


def candidate_snapshot(candidate: dict[str, Any], review_item: dict[str, Any] | None) -> dict[str, Any]:
    snapshot = review_item.get("evidence_snapshot") if review_item and isinstance(review_item.get("evidence_snapshot"), dict) else {}
    return {
        "organism_name": candidate.get("organism_name") or (item_name(review_item) if review_item else ""),
        "classification": candidate.get("classification") or (review_item.get("classification") if review_item else ""),
        "integrated_causative_level": candidate.get("integrated_causative_level"),
        "mngs_signal_tier": candidate.get("mngs_signal_tier"),
        "rank_priority": candidate.get("rank_priority") or snapshot.get("rank_priority"),
        "reads": candidate.get("reads") if candidate.get("reads") is not None else snapshot.get("reads"),
        "reads_tier": candidate.get("reads_tier") or snapshot.get("reads_tier"),
        "reads_percentile": candidate.get("reads_percentile") if candidate.get("reads_percentile") is not None else snapshot.get("reads_percentile"),
        "dominance_tier": candidate.get("dominance_tier") or snapshot.get("dominance_tier"),
        "specimen_class": candidate.get("specimen_class") or snapshot.get("specimen_class"),
        "specimen_alignment": candidate.get("specimen_alignment"),
        "rank_rule": candidate.get("rank_rule") or snapshot.get("rank_rule"),
        "support_modules": support_modules(candidate, review_item),
        "non_host_support_modules": as_list(snapshot.get("non_host_support_modules")),
        "applied_rules": as_list(candidate.get("applied_rules")),
        "key_evidence": candidate.get("key_evidence") if isinstance(candidate.get("key_evidence"), dict) else {},
        "integrated_reasoning": as_list(candidate.get("integrated_reasoning")),
    }


def pathogen_category(name: str, classification: str = "") -> str:
    key = pathogen_names.canonical_key(name)
    lower = name.lower()
    if "candida" in key or key == "yeast":
        return "Candida/yeast"
    if any(token in key for token in ("hsv", "herpessimplex", "cytomegalovirus", "cmv", "vzv", "varicellazoster", "ebv", "epsteinbarr")):
        return "Herpesvirus/reactivation"
    if any(token in lower for token in ("covid", "sars-cov-2", "influenza", "adenovirus", "rsv", "parainfluenza", "metapneumovirus")):
        return "Respiratory virus"
    if any(token in key for token in ("bacteroides", "prevotella", "fusobacterium", "porphyromonas", "parvimonas", "peptostreptococcus", "finegoldia")):
        return "Strict anaerobe/aspiration flora"
    if any(token in key for token in ("corynebacterium", "cutibacterium", "staphylococcushaemolyticus", "staphylococcusepidermidis")):
        return "Skin/airway colonizer Gram-positive"
    if any(token in key for token in ("chryseobacterium", "sphingomonas", "ralstonia", "achromobacter")):
        return "Environmental/water-associated GNB"
    if any(token in key for token in ("enterococcus", "clostridioidesdifficile", "clostridiumdifficile")):
        return "GI/urinary/non-pulmonary-prone organism"
    if any(token in key for token in ("escherichia", "klebsiella", "enterobacter", "cronobacter", "serratia", "proteus")):
        return "Enterobacterales/hospital GNB"
    return classification or "Other/unclassified"


def search_queries(name: str, category: str) -> list[str]:
    queries = [
        f'"{name}" pneumonia adult',
        f'"{name}" lower respiratory tract infection',
        f'"{name}" respiratory colonization contaminant',
    ]
    if category == "Candida/yeast":
        queries.extend([
            f'"{name}" Candida pneumonia BAL colonization',
            "Candida respiratory secretions colonization pneumonia guideline",
            "IDSA candidiasis guideline respiratory secretions Candida pneumonia",
        ])
    elif category == "Herpesvirus/reactivation":
        queries.extend([
            f'"{name}" BAL pneumonia immunocompromised viral load',
            f'"{name}" reactivation shedding lower respiratory tract',
        ])
    elif category == "Strict anaerobe/aspiration flora":
        queries.extend([
            f'"{name}" aspiration pneumonia lung abscess empyema',
            f'"{name}" anaerobic pneumonia metagenomic sequencing',
        ])
    elif category == "Skin/airway colonizer Gram-positive":
        queries.extend([
            f'"{name}" ventilator associated pneumonia case series',
            f'"{name}" respiratory tract colonization infection criteria',
        ])
    elif category == "Environmental/water-associated GNB":
        queries.extend([
            f'"{name}" hospital acquired pneumonia immunocompromised',
            f'"{name}" water environmental organism respiratory infection',
        ])
    return queries


def rag_questions() -> list[str]:
    return [
        "Does reliable literature support this organism as a cause of adult pneumonia, VAP/HAP, aspiration pneumonia, lung abscess, empyema, or pulmonary infection in immunocompromised hosts?",
        "In respiratory specimens, is this organism more often interpreted as a pathogen, colonizer, contaminant, environmental background, or bystander?",
        "What evidence is usually required to call this a true pulmonary pathogen: repeated lower-respiratory culture, BAL quantitative burden, blood/sterile-site evidence, PCR/viral load, histopathology, compatible imaging, or treatment response?",
        "Does this patient's local evidence meet those criteria? Which details support infection, and which argue against it?",
        "What RAG tier is recommended: review_high_priority/keep, review_context_needed, review_low_specificity, or audit_reject?",
        "If literature is sparse or conflicting, state uncertainty and avoid upgrading based only on a single weak case report.",
    ]


def evidence_card_template() -> dict[str, Any]:
    return {
        "literature_search": {
            "databases": ["PubMed", "guidelines when applicable", "PMC full text when available"],
            "target_paper_count": 10,
            "prefer_guideline_or_review_first": True,
        },
        "papers": [
            {
                "rank": None,
                "title": "",
                "year": None,
                "pmid_or_doi": "",
                "study_type": "guideline/review/case series/case report/diagnostic study",
                "population_similarity": "high/medium/low",
                "specimen_similarity": "BAL/sputum/ET aspirate/blood/sterile site/other",
                "support_direction": "support_strong/support_weak/neutral/against/irrelevant",
                "quality": "high/medium/low",
                "key_point": "",
            }
        ],
        "paper_vote_summary": {
            "support_strong": 0,
            "support_weak": 0,
            "neutral": 0,
            "against": 0,
            "irrelevant": 0,
        },
        "local_evidence_assessment": {
            "supports_true_pathogen": [],
            "argues_against_true_pathogen": [],
            "missing_key_evidence": [],
        },
        "rag_decision": {
            "recommended_tier": "review_high_priority/review_context_needed/review_low_specificity/audit_reject",
            "should_send_to_doctor": None,
            "confidence": "high/medium/low",
            "one_sentence_reason": "",
        },
    }


def safe_case_id(patient_id: str, organism: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", organism).strip("_")
    return f"P{patient_id}_{slug}"


def build_case(
    label: str,
    patient_id: str,
    answer: compare.NameEntry,
    det: dict[str, Any],
    review: dict[str, Any],
    final_summary: dict[str, Any],
) -> dict[str, Any] | None:
    picked = picked_items(det)
    if match_in_items(answer, picked):
        return None

    review_match = match_in_review(answer, review_items(review))
    current_tier = "not_found_in_current_review"
    review_item: dict[str, Any] | None = None
    match_type = ""
    matched_output = ""
    if review_match:
        current_tier, review_item, match_type = review_match
        matched_output = item_name(review_item)

    candidate = find_candidate(answer, det, review_item)
    snap = candidate_snapshot(candidate, review_item)
    display_name = matched_output or answer.raw
    category = pathogen_category(display_name, str(snap.get("classification") or ""))
    hospital = matching_hospital_evidence(final_summary, answer)
    is_rag_visible = current_tier in RAG_VISIBLE_TIERS
    reason = ""
    if review_item:
        reason = str(review_item.get("rationale_zh") or review_item.get("review_tier_reason") or "")

    if is_rag_visible:
        suggested_use = "primary_rag_input"
    elif current_tier in {"review_low_specificity", "omitted_with_reason"}:
        suggested_use = "do_not_send_by_default_but_useful_for_validation"
    else:
        suggested_use = "outside_current_rag_queue_or_missing_source_evidence"

    return {
        "case_id": safe_case_id(patient_id, answer.raw),
        "dataset_label": label,
        "patient_id": patient_id,
        "benchmark_answer": answer.raw,
        "matched_output_organism": matched_output,
        "match_type": match_type,
        "picked_status": "answer_missed_by_formal_picked",
        "current_review_tier": current_tier,
        "current_rag_visible": is_rag_visible,
        "suggested_rag_use": suggested_use,
        "pathogen_category": category,
        "current_formal_picked_pathogens": [item_name(item) for item in picked if item_name(item)],
        "local_evidence": {
            "final_infection_likelihood": final_summary.get("final_infection_likelihood"),
            "dominant_source": final_summary.get("dominant_source"),
            "dominant_pathogen_type": final_summary.get("dominant_pathogen_type"),
            "data_quality_tier": final_summary.get("data_quality_tier"),
            "host_context": final_summary.get("host_context") or det.get("host_context"),
            "candidate_snapshot": snap,
            "hospital_evidence_matching_answer": hospital,
            "current_review_reason": reason,
            "risk_if_ignored": review_item.get("risk_if_ignored") if review_item else "",
            "why_not_auto_upgrade": review_item.get("why_not_auto_upgrade") if review_item else "",
            "suggested_next_check_from_review": review_item.get("suggested_next_check") if review_item else "",
        },
        "rag_task": {
            "fixed_questions": rag_questions(),
            "suggested_search_queries": search_queries(display_name, category),
            "evidence_card_template": evidence_card_template(),
        },
    }


CSV_FIELDS = [
    "case_id",
    "patient_id",
    "benchmark_answer",
    "matched_output_organism",
    "current_review_tier",
    "current_rag_visible",
    "suggested_rag_use",
    "pathogen_category",
    "picked_pathogens",
    "rank_priority",
    "reads",
    "reads_tier",
    "reads_percentile",
    "dominance_tier",
    "specimen_class",
    "support_modules",
    "same_organism_hospital_evidence_count",
    "same_genus_hospital_evidence_count",
    "current_review_reason",
]


def csv_row(case: dict[str, Any]) -> dict[str, Any]:
    snap = case["local_evidence"]["candidate_snapshot"]
    hospital = case["local_evidence"]["hospital_evidence_matching_answer"]
    return {
        "case_id": case["case_id"],
        "patient_id": case["patient_id"],
        "benchmark_answer": case["benchmark_answer"],
        "matched_output_organism": case["matched_output_organism"],
        "current_review_tier": case["current_review_tier"],
        "current_rag_visible": case["current_rag_visible"],
        "suggested_rag_use": case["suggested_rag_use"],
        "pathogen_category": case["pathogen_category"],
        "picked_pathogens": joined(case["current_formal_picked_pathogens"]),
        "rank_priority": snap.get("rank_priority"),
        "reads": snap.get("reads"),
        "reads_tier": snap.get("reads_tier"),
        "reads_percentile": snap.get("reads_percentile"),
        "dominance_tier": snap.get("dominance_tier"),
        "specimen_class": snap.get("specimen_class"),
        "support_modules": joined(snap.get("support_modules") or []),
        "same_organism_hospital_evidence_count": len(hospital.get("same_organism_or_alias") or []),
        "same_genus_hospital_evidence_count": len(hospital.get("same_genus_or_related") or []),
        "current_review_reason": case["local_evidence"]["current_review_reason"],
    }


def markdown_case(case: dict[str, Any]) -> str:
    snap = case["local_evidence"]["candidate_snapshot"]
    hospital = case["local_evidence"]["hospital_evidence_matching_answer"]
    questions = "\n".join(f"{idx}. {question}" for idx, question in enumerate(case["rag_task"]["fixed_questions"], 1))
    queries = "\n".join(f"- {query}" for query in case["rag_task"]["suggested_search_queries"])
    picked = joined(case["current_formal_picked_pathogens"]) or "-"
    return f"""# {case['case_id']}

## Case
- patient_id: {case['patient_id']}
- benchmark_answer: {case['benchmark_answer']}
- matched_output_organism: {case['matched_output_organism'] or '-'}
- current_review_tier: {case['current_review_tier']}
- current_rag_visible: {case['current_rag_visible']}
- suggested_rag_use: {case['suggested_rag_use']}
- pathogen_category: {case['pathogen_category']}
- formal_picked_pathogens: {picked}

## Local Evidence
- final_infection_likelihood: {case['local_evidence']['final_infection_likelihood']}
- dominant_source: {case['local_evidence']['dominant_source']}
- dominant_pathogen_type: {case['local_evidence']['dominant_pathogen_type']}
- host_context: {json.dumps(case['local_evidence']['host_context'], ensure_ascii=False)}
- mNGS: rank={snap.get('rank_priority')}, reads={snap.get('reads')}, reads_tier={snap.get('reads_tier')}, percentile={snap.get('reads_percentile')}, dominance={snap.get('dominance_tier')}, specimen={snap.get('specimen_class')}
- support_modules: {joined(snap.get('support_modules') or []) or '-'}
- same_organism_hospital_evidence_count: {len(hospital.get('same_organism_or_alias') or [])}
- same_genus_hospital_evidence_count: {len(hospital.get('same_genus_or_related') or [])}
- current_review_reason: {case['local_evidence']['current_review_reason'] or '-'}

## RAG Questions
{questions}

## Suggested Search Queries
{queries}

## Evidence Card To Fill
```json
{json.dumps(case['rag_task']['evidence_card_template'], ensure_ascii=False, indent=2)}
```
"""


def write_case_cards(output_dir: Path, cases: Sequence[dict[str, Any]]) -> None:
    card_dir = output_dir / "case_cards_md"
    card_dir.mkdir(parents=True, exist_ok=True)
    for case in cases:
        (card_dir / f"{case['case_id']}.md").write_text(markdown_case(case), encoding="utf-8")


def write_readme(output_dir: Path, label: str, merged_suffix: str, total_cases: int, rag_visible_cases: int) -> None:
    text = f"""# RAG prototype input

Source label: `{label}`
Merged suffix: `{merged_suffix}`

This folder packages benchmark-answer organisms that were missed by formal picked output.
It is meant for downstream RAG testing, not as a final medical decision.

Files:
- `rag_cases_rag_visible_only.jsonl`: first batch to send to RAG. These are currently in `review_high_priority` or `review_context_needed`.
- `rag_cases_all_picked_missed_answers.json`: all picked-missed answer cases, including low/audit/not-found cases for validation.
- `rag_cases_summary.csv`: compact table for human inspection.
- `case_cards_md/`: one Markdown evidence-card draft per case.

RAG should judge only the review/RAG organisms. It should return one of:
- `review_high_priority` / keep
- `review_context_needed`
- `review_low_specificity`
- `audit_reject`

Counts:
- picked-missed answer cases: {total_cases}
- current RAG-visible cases: {rag_visible_cases}
"""
    (output_dir / "README_for_RAG.md").write_text(text, encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    answers, _, answer_column, answer_row_count = compare.load_answers(args.answer_csv, args.answer_column)
    cases: list[dict[str, Any]] = []
    missing_outputs: list[str] = []

    for patient_dir in sorted(args.patient_root.glob("NGS_patient_*_json"), key=compare.patient_sort_key):
        pid = compare.patient_id(patient_dir)
        if pid not in answers:
            continue
        merged_path = patient_output_path(patient_dir, args.merged_suffix)
        if not merged_path.exists():
            missing_outputs.append(str(merged_path))
            continue
        merged = read_json(merged_path)
        summary_path = final_summary_path(patient_dir)
        final_summary = read_json(summary_path) if summary_path.exists() else {}
        det = deterministic_payload(merged)
        review = review_payload(merged)
        for answer in answers.get(pid, []):
            case = build_case(args.label, pid, answer, det, review, final_summary)
            if case is not None:
                cases.append(case)

    rag_visible = [case for case in cases if case["current_rag_visible"]]
    category_counts = Counter(case["pathogen_category"] for case in cases)
    tier_counts = Counter(case["current_review_tier"] for case in cases)
    summary = {
        "label": args.label,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "patient_root": str(args.patient_root),
        "merged_suffix": args.merged_suffix,
        "answer_csv": str(args.answer_csv),
        "answer_column": answer_column,
        "answer_csv_row_count": answer_row_count,
        "picked_missed_answer_case_count": len(cases),
        "rag_visible_case_count": len(rag_visible),
        "non_rag_visible_case_count": len(cases) - len(rag_visible),
        "category_counts": dict(category_counts),
        "tier_counts": dict(tier_counts),
        "missing_outputs": missing_outputs,
        "cases": cases,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json(args.output_dir / "rag_cases_all_picked_missed_answers.json", summary)
    write_jsonl(args.output_dir / "rag_cases_rag_visible_only.jsonl", rag_visible)
    write_csv(args.output_dir / "rag_cases_summary.csv", [csv_row(case) for case in cases], CSV_FIELDS)
    write_case_cards(args.output_dir, cases)
    write_readme(args.output_dir, args.label, args.merged_suffix, len(cases), len(rag_visible))

    print(f"output_dir={args.output_dir}")
    print(f"picked_missed_answer_cases={len(cases)}")
    print(f"rag_visible_cases={len(rag_visible)}")
    print(f"non_rag_visible_cases={len(cases) - len(rag_visible)}")
    print("tier_counts=" + json.dumps(dict(tier_counts), ensure_ascii=False, sort_keys=True))
    print("category_counts=" + json.dumps(dict(category_counts), ensure_ascii=False, sort_keys=True))
    print(f"missing_outputs={len(missing_outputs)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

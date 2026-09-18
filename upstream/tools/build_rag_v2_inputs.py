"""Build answer-blind RAG v2 production inputs and private acceptance reports."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import build_rag_prototype_cases as prototype  # noqa: E402
from tools import candidate_evidence_profile as evidence_profiles  # noqa: E402
from tools import compare_aliasclean_baseline as compare  # noqa: E402
from tools import organism_taxonomy_classifier as organism_taxonomy  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402

PRODUCTION_SCHEMA_VERSION = "rag_production_case_v2.0"
ADJUDICATION_SCHEMA_VERSION = "rag_adjudication_v2.0"
PRODUCTION_TIERS = ("review_high_priority", "review_context_needed")
TIER_PRIORITY = {"review_high_priority": 2, "review_context_needed": 1}
FORBIDDEN_PRODUCTION_KEYS = {
    "answer",
    "answers",
    "answer_hit",
    "answer_organism",
    "benchmark_answer",
    "expected_in_production_queue",
    "ground_truth",
    "is_answer",
    "match_type",
    "matched_output_organism",
    "picked_status",
}


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def clean_name(item: dict[str, Any]) -> str:
    return pathogen_names.clean_display_text(item.get("organism_name") or item.get("name"))


def safe_case_id(patient_id: str, organism_name: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", organism_name).strip("_")
    return f"P{patient_id}_{slug}"


def merged_path(patient_dir: Path, suffix: str) -> Path:
    patient_id = compare.patient_id(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{patient_id}_{suffix}.json"


def final_summary_path(patient_dir: Path) -> Path:
    patient_id = compare.patient_id(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{patient_id}_final_summary_with_filmarray_deterministic.json"


def deterministic_payload(merged: dict[str, Any]) -> dict[str, Any]:
    return as_dict(merged.get("deterministic_max")) or merged


def review_payload(merged: dict[str, Any]) -> dict[str, Any]:
    return as_dict(merged.get("llm_missed_candidate_review"))


def picked_items(deterministic: dict[str, Any]) -> list[dict[str, Any]]:
    best = as_dict(deterministic.get("best_available_summary"))
    return [item for item in as_list(best.get("picked_pathogens")) if isinstance(item, dict)]


def candidate_index(deterministic: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for item in as_list(deterministic.get("pathogen_candidates")):
        if not isinstance(item, dict):
            continue
        key = pathogen_names.canonical_key(clean_name(item))
        if key and key not in result:
            result[key] = item
    return result


def candidate_for_review(
    review_item: dict[str, Any],
    candidates: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    target = compare.name_entry(clean_name(review_item))
    if target is None:
        return {}
    if target.canonical in candidates:
        return candidates[target.canonical]
    for candidate in candidates.values():
        candidate_entry = compare.name_entry(clean_name(candidate))
        if candidate_entry and candidate_entry.canonical == target.canonical:
            return candidate
    return {}


def hospital_evidence_for_name(
    summary: dict[str, Any],
    organism_name: str,
) -> dict[str, list[dict[str, Any]]]:
    """Keep exact aliases and same-genus context separate."""
    target = compare.name_entry(organism_name)
    exact: list[dict[str, Any]] = []
    same_genus: list[dict[str, Any]] = []
    if target is None:
        return {"same_organism_or_alias": exact, "same_genus_or_related": same_genus}
    for item in as_list(summary.get("hospital_organism_evidence")):
        if not isinstance(item, dict):
            continue
        entry = compare.name_entry(item.get("organism_name"))
        if entry is None:
            continue
        if entry.canonical == target.canonical:
            exact.append(item)
        elif entry.genus and target.genus and entry.genus == target.genus:
            same_genus.append(item)
    return {"same_organism_or_alias": exact, "same_genus_or_related": same_genus}


def support_modules(candidate: dict[str, Any], review_item: dict[str, Any]) -> list[str]:
    modules: list[str] = []
    key_evidence = as_dict(candidate.get("key_evidence"))
    snapshot = as_dict(review_item.get("evidence_snapshot"))
    modules.extend(str(value) for value in as_list(key_evidence.get("support_modules")))
    modules.extend(str(value) for value in as_list(snapshot.get("support_modules")))
    module_summary = as_dict(candidate.get("module_support_summary"))
    modules.extend(str(key) for key, value in module_summary.items() if value == "Support")
    return list(dict.fromkeys(module for module in modules if module))


def candidate_snapshot(candidate: dict[str, Any], review_item: dict[str, Any]) -> dict[str, Any]:
    snapshot = as_dict(review_item.get("evidence_snapshot"))
    key_evidence = as_dict(candidate.get("key_evidence"))
    return {
        "organism_name": clean_name(candidate) or clean_name(review_item),
        "classification": candidate.get("classification") or review_item.get("classification"),
        "integrated_causative_level": candidate.get("integrated_causative_level"),
        "mngs_signal_tier": candidate.get("mngs_signal_tier"),
        "rank_priority": candidate.get("rank_priority") or snapshot.get("rank_priority"),
        "reads": candidate.get("reads") if candidate.get("reads") is not None else snapshot.get("reads"),
        "reads_tier": candidate.get("reads_tier") or snapshot.get("reads_tier"),
        "reads_percentile": (
            candidate.get("reads_percentile")
            if candidate.get("reads_percentile") is not None
            else snapshot.get("reads_percentile")
        ),
        "dominance_tier": candidate.get("dominance_tier") or snapshot.get("dominance_tier"),
        "specimen_class": candidate.get("specimen_class") or snapshot.get("specimen_class"),
        "specimen_alignment": candidate.get("specimen_alignment"),
        "rank_rule": candidate.get("rank_rule") or snapshot.get("rank_rule"),
        "support_modules": support_modules(candidate, review_item),
        "non_host_support_modules": as_list(snapshot.get("non_host_support_modules")),
        "applied_rules": as_list(candidate.get("applied_rules")),
        "formal_pick_exclusion_rules": as_list(
            key_evidence.get("formal_pick_exclusion_rules")
        ),
        "formal_pick_exclusion_rule": key_evidence.get("formal_pick_exclusion_rule"),
        "candida_evidence_strength_profile": (
            key_evidence.get("candida_evidence_strength_profile")
            or snapshot.get("candida_evidence_strength_profile")
            or {}
        ),
        "candida_invasive_hospital_support": (
            key_evidence.get("candida_invasive_hospital_support")
            or snapshot.get("candida_invasive_hospital_support")
            or []
        ),
    }


def competing_pathogens(
    picked: Sequence[dict[str, Any]],
    target_key: str,
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for item in picked:
        name = clean_name(item)
        if not name or pathogen_names.canonical_key(name) == target_key:
            continue
        result.append(
            {
                "organism_name": name,
                "classification": item.get("classification"),
                "integrated_causative_level": item.get("integrated_causative_level"),
                "mngs_signal_tier": item.get("mngs_signal_tier"),
                "rank_priority": item.get("rank_priority"),
                "reads": item.get("reads"),
                "reads_tier": item.get("reads_tier"),
                "reads_percentile": item.get("reads_percentile"),
                "dominance_tier": item.get("dominance_tier"),
                "support_modules": as_list(as_dict(item.get("key_evidence")).get("support_modules")),
            }
        )
    return result


def negative_evidence(
    snapshot: dict[str, Any],
    hospital: dict[str, list[dict[str, Any]]],
    competitors: Sequence[dict[str, Any]],
    *,
    organism_name: str = "",
) -> list[dict[str, str]]:
    return evidence_profiles.build_negative_evidence(
        snapshot,
        hospital,
        competitors,
        organism_name=organism_name,
    )


def review_policy_version(review: dict[str, Any]) -> str | None:
    policy = as_dict(review.get("review_tiering_policy"))
    value = policy.get("version")
    return str(value) if value is not None else None


def build_case(
    *,
    label: str,
    merged_suffix: str,
    patient_id: str,
    tier: str,
    review_item: dict[str, Any],
    deterministic: dict[str, Any],
    review: dict[str, Any],
    summary: dict[str, Any],
) -> dict[str, Any] | None:
    organism_name = clean_name(review_item)
    canonical_key = pathogen_names.canonical_key(organism_name)
    if not organism_name or not canonical_key:
        return None
    candidates = candidate_index(deterministic)
    candidate = candidate_for_review(review_item, candidates)
    snapshot = candidate_snapshot(candidate, review_item)
    classification = str(snapshot.get("classification") or "")
    taxonomy_profile = review_item.get("taxonomy_profile")
    if not isinstance(taxonomy_profile, dict):
        taxonomy_profile = candidate.get("taxonomy_profile")
    if not isinstance(taxonomy_profile, dict):
        taxonomy_profile = organism_taxonomy.classify_organism(
            organism_name,
            biological_class=classification,
        )
    category = prototype.pathogen_category(organism_name, classification)
    hospital = hospital_evidence_for_name(summary, organism_name)
    picked = picked_items(deterministic)
    competitors = competing_pathogens(picked, canonical_key)
    evidence_profile = evidence_profiles.build_candidate_evidence_profile(
        organism_name=organism_name,
        snapshot=snapshot,
        hospital=hospital,
        competitors=competitors,
    )
    negatives = negative_evidence(
        snapshot,
        hospital,
        competitors,
        organism_name=organism_name,
    )
    policy = as_dict(review.get("review_tiering_policy"))
    return {
        "schema_version": PRODUCTION_SCHEMA_VERSION,
        "case_id": safe_case_id(patient_id, organism_name),
        "dataset_label": label,
        "patient_id": patient_id,
        "organism": {
            "display_name": organism_name,
            "canonical_key": canonical_key,
            "classification": classification or None,
            "pathogen_category": category,
            "taxonomy_unknown": taxonomy_profile.get("mapping_status") == "unmapped",
            "taxonomy_profile": taxonomy_profile,
        },
        "current_model_state": {
            "source_tier": tier,
            "formal_picked": False,
            "current_formal_picked_pathogens": [clean_name(item) for item in picked if clean_name(item)],
            "competing_pathogens": competitors,
            "current_review_reason": str(
                review_item.get("rationale_zh") or review_item.get("review_tier_reason") or ""
            ),
            "risk_if_ignored": str(review_item.get("risk_if_ignored") or ""),
            "why_not_auto_upgrade": str(review_item.get("why_not_auto_upgrade") or ""),
            "suggested_next_check": str(review_item.get("suggested_next_check") or ""),
        },
        "local_evidence": {
            "case_context": {
                "final_infection_likelihood": summary.get("final_infection_likelihood")
                or deterministic.get("final_infection_likelihood"),
                "dominant_source": summary.get("dominant_source") or deterministic.get("dominant_source"),
                "dominant_pathogen_type": summary.get("dominant_pathogen_type")
                or deterministic.get("dominant_pathogen_type"),
                "data_quality_tier": summary.get("data_quality_tier"),
                "host_context": summary.get("host_context") or deterministic.get("host_context") or {},
                "time_span": summary.get("time_span") or {},
                "priority_flags": as_list(summary.get("priority_flags")),
            },
            "mngs_evidence": snapshot,
            "same_organism_hospital_evidence": hospital["same_organism_or_alias"],
            "same_genus_or_related_hospital_evidence": hospital["same_genus_or_related"],
            "negative_evidence": negatives,
            "candidate_evidence_profile": evidence_profile,
            "data_gaps": as_list(summary.get("data_gaps")) or as_list(deterministic.get("data_gaps")),
        },
        "rag_task": {
            "knowledge_card_key": f"adult_pulmonary:{canonical_key}",
            "retrieval_priority": "high" if tier == "review_high_priority" else "standard",
            "fixed_questions": prototype.rag_questions(),
            "suggested_search_queries": prototype.search_queries(organism_name, category),
            "required_output_schema": ADJUDICATION_SCHEMA_VERSION,
        },
        "provenance": {
            "merged_suffix": merged_suffix,
            "deterministic_rule_version": deterministic.get("rule_version"),
            "review_policy_version": review_policy_version(review) or policy.get("policy_version"),
            "candidate_evidence_profile_version": evidence_profiles.PROFILE_VERSION,
        },
    }


def find_forbidden_keys(value: Any, path: str = "$") -> list[str]:
    findings: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            child_path = f"{path}.{key}"
            if str(key).lower() in FORBIDDEN_PRODUCTION_KEYS:
                findings.append(child_path)
            findings.extend(find_forbidden_keys(child, child_path))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            findings.extend(find_forbidden_keys(child, f"{path}[{index}]"))
    return findings


def validate_production_case(case: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    required = {
        "schema_version",
        "case_id",
        "dataset_label",
        "patient_id",
        "organism",
        "current_model_state",
        "local_evidence",
        "rag_task",
        "provenance",
    }
    missing = sorted(required - set(case))
    if missing:
        errors.append("missing required keys: " + ", ".join(missing))
    if case.get("schema_version") != PRODUCTION_SCHEMA_VERSION:
        errors.append("invalid schema_version")
    tier = as_dict(case.get("current_model_state")).get("source_tier")
    if tier not in PRODUCTION_TIERS:
        errors.append(f"invalid source_tier: {tier}")
    if as_dict(case.get("current_model_state")).get("formal_picked") is not False:
        errors.append("formal_picked must be false for the review-only production queue")
    errors.extend("forbidden benchmark key: " + path for path in find_forbidden_keys(case))
    return errors


def build_knowledge_requests(cases: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    requests: dict[str, dict[str, Any]] = {}
    for case in cases:
        task = case["rag_task"]
        key = task["knowledge_card_key"]
        if key in requests:
            continue
        requests[key] = {
            "schema_version": "rag_knowledge_request_v2.0",
            "knowledge_card_key": key,
            "organism": case["organism"],
            "clinical_domain": "adult pulmonary infection",
            "suggested_search_queries": task["suggested_search_queries"],
            "retrieval_order": [
                "guideline",
                "systematic_review_or_review",
                "cohort_or_case_series",
                "diagnostic_study",
                "case_report_for_rare_organisms_only",
                "colonization_contamination_or_microbiome_counterevidence",
            ],
        }
    return list(requests.values())


def production_readme(manifest: dict[str, Any]) -> str:
    return f"""# RAG v2 production input

這個資料夾是正式 RAG 輸入，不含 benchmark answer 或 answer-hit 標記。

來源 merge suffix：`{manifest['merged_suffix']}`

## 要給 RAG 的檔案

- `production_queue.jsonl`：逐一判斷病人與候選菌的病例卡，只含 review_high_priority 與 review_context_needed。
- `knowledge_card_requests.jsonl`：依菌種去重的文獻檢索工作，相同菌種跨病人共用 knowledge card。
- `schemas/rag_v2_adjudication_output.schema.json`：RAG 回傳格式。
- `production_queue_summary.csv`：人工快速檢查用，不應取代完整 JSONL。

## 不可混入 production 的資料

`evaluation_private` 內含 benchmark labels，只能在 RAG 完成後做驗證，不能提供給檢索或判斷模型。

## 數量

- production cases：{manifest['production_case_count']}
- knowledge card requests：{manifest['knowledge_request_count']}
- review_high_priority：{manifest['tier_counts'].get('review_high_priority', 0)}
- review_context_needed：{manifest['tier_counts'].get('review_context_needed', 0)}

RAG 不應直接覆寫 deterministic picked。請把結果寫成獨立 `rag_adjudication`，依 schema 回傳 primary_pathogen、secondary_pathogen、context_only、reject 或 insufficient_evidence。
"""


def build_production(args: argparse.Namespace) -> int:
    patient_dirs = sorted(args.patient_root.glob("NGS_patient_*_json"), key=compare.patient_sort_key)
    cases_by_key: dict[tuple[str, str], dict[str, Any]] = {}
    missing_merge: list[str] = []
    missing_summary: list[str] = []

    for patient_dir in patient_dirs:
        patient_id = compare.patient_id(patient_dir)
        merge_file = merged_path(patient_dir, args.merged_suffix)
        if not merge_file.exists():
            missing_merge.append(patient_id)
            continue
        merged = read_json(merge_file)
        deterministic = deterministic_payload(merged)
        review = review_payload(merged)
        summary_file = final_summary_path(patient_dir)
        if summary_file.exists():
            summary = read_json(summary_file)
        else:
            summary = {}
            missing_summary.append(patient_id)
        for tier in PRODUCTION_TIERS:
            for item in as_list(review.get(tier)):
                if not isinstance(item, dict):
                    continue
                case = build_case(
                    label=args.label,
                    merged_suffix=args.merged_suffix,
                    patient_id=patient_id,
                    tier=tier,
                    review_item=item,
                    deterministic=deterministic,
                    review=review,
                    summary=summary,
                )
                if case is None:
                    continue
                key = (patient_id, case["organism"]["canonical_key"])
                previous = cases_by_key.get(key)
                if previous is None or TIER_PRIORITY[tier] > TIER_PRIORITY[previous["current_model_state"]["source_tier"]]:
                    cases_by_key[key] = case

    cases = sorted(
        cases_by_key.values(),
        key=lambda case: (
            int(case["patient_id"]) if str(case["patient_id"]).isdigit() else 10**9,
            -TIER_PRIORITY[case["current_model_state"]["source_tier"]],
            case["organism"]["canonical_key"],
        ),
    )
    validation_errors: list[dict[str, Any]] = []
    for case in cases:
        errors = validate_production_case(case)
        if errors:
            validation_errors.append({"case_id": case.get("case_id"), "errors": errors})
    leakage_findings = [
        {"case_id": case["case_id"], "paths": findings}
        for case in cases
        if (findings := find_forbidden_keys(case))
    ]
    if validation_errors or leakage_findings:
        raise RuntimeError(
            f"Production validation failed: errors={len(validation_errors)}, leakage={len(leakage_findings)}"
        )

    knowledge_requests = build_knowledge_requests(cases)
    tier_counts = Counter(case["current_model_state"]["source_tier"] for case in cases)
    category_counts = Counter(case["organism"]["pathogen_category"] for case in cases)
    manifest = {
        "schema_version": "rag_production_manifest_v2.0",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "dataset_label": args.label,
        "patient_root": str(args.patient_root),
        "merged_suffix": args.merged_suffix,
        "patient_directory_count": len(patient_dirs),
        "production_case_count": len(cases),
        "knowledge_request_count": len(knowledge_requests),
        "tier_counts": dict(tier_counts),
        "category_counts": dict(category_counts),
        "missing_merge_patient_ids": missing_merge,
        "missing_summary_patient_ids": missing_summary,
        "production_validation": {
            "valid_case_count": len(cases),
            "invalid_case_count": 0,
            "forbidden_benchmark_key_count": 0,
            "answer_source_read": False,
        },
        "schema_files": {
            "production_input": "schemas/rag_v2_production_input.schema.json",
            "adjudication_output": "schemas/rag_v2_adjudication_output.schema.json",
        },
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_jsonl(args.output_dir / "production_queue.jsonl", cases)
    write_jsonl(args.output_dir / "knowledge_card_requests.jsonl", knowledge_requests)
    write_json(args.output_dir / "production_manifest.json", manifest)
    write_json(
        args.output_dir / "schemas" / "rag_v2_production_input.schema.json",
        read_json(REPO_ROOT / "rules" / "rag_v2_production_input.schema.json"),
    )
    write_json(
        args.output_dir / "schemas" / "rag_v2_adjudication_output.schema.json",
        read_json(REPO_ROOT / "rules" / "rag_v2_adjudication_output.schema.json"),
    )
    (args.output_dir / "README_繁中.md").write_text(production_readme(manifest), encoding="utf-8")
    write_csv(
        args.output_dir / "production_queue_summary.csv",
        [
            {
                "case_id": case["case_id"],
                "patient_id": case["patient_id"],
                "organism_name": case["organism"]["display_name"],
                "source_tier": case["current_model_state"]["source_tier"],
                "pathogen_category": case["organism"]["pathogen_category"],
                "reads": case["local_evidence"]["mngs_evidence"].get("reads"),
                "reads_tier": case["local_evidence"]["mngs_evidence"].get("reads_tier"),
                "dominance_tier": case["local_evidence"]["mngs_evidence"].get("dominance_tier"),
                "same_organism_hospital_evidence_count": len(
                    case["local_evidence"]["same_organism_hospital_evidence"]
                ),
                "negative_evidence_count": len(case["local_evidence"]["negative_evidence"]),
                "knowledge_card_key": case["rag_task"]["knowledge_card_key"],
            }
            for case in cases
        ],
        (
            "case_id",
            "patient_id",
            "organism_name",
            "source_tier",
            "pathogen_category",
            "reads",
            "reads_tier",
            "dominance_tier",
            "same_organism_hospital_evidence_count",
            "negative_evidence_count",
            "knowledge_card_key",
        ),
    )
    print(f"production_cases={len(cases)}")
    print(f"knowledge_requests={len(knowledge_requests)}")
    print("tier_counts=" + json.dumps(dict(tier_counts), ensure_ascii=False, sort_keys=True))
    print(f"missing_merge={len(missing_merge)}")
    print(f"forbidden_benchmark_keys=0")
    return 0


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def evaluation_match(case: dict[str, Any], patient_id: str, expected_name: str) -> bool:
    if str(case.get("patient_id")) != patient_id:
        return False
    actual = compare.name_entry(as_dict(case.get("organism")).get("display_name"))
    expected = compare.name_entry(expected_name)
    return bool(actual and expected and compare.names_match(actual, expected))


def run_evaluation(args: argparse.Namespace) -> int:
    production_cases = load_jsonl(args.production_queue)
    labels = read_json(args.acceptance_cases)
    rows: list[dict[str, Any]] = []
    for expected in as_list(labels.get("cases")):
        if not isinstance(expected, dict):
            continue
        patient_id = str(expected.get("patient_id") or "")
        organism = str(expected.get("answer_organism") or "")
        expected_present = bool(expected.get("expected_in_production_queue"))
        matches = [
            case for case in production_cases if evaluation_match(case, patient_id, organism)
        ]
        actual_present = bool(matches)
        rows.append(
            {
                "patient_id": patient_id,
                "answer_organism": organism,
                "expected_in_production_queue": expected_present,
                "actual_in_production_queue": actual_present,
                "acceptance_pass": actual_present == expected_present,
                "matched_case_ids": "; ".join(case["case_id"] for case in matches),
                "actual_tiers": "; ".join(
                    case["current_model_state"]["source_tier"] for case in matches
                ),
            }
        )
    passed = sum(bool(row["acceptance_pass"]) for row in rows)
    report = {
        "schema_version": "rag_acceptance_report_v2.0",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "production_queue": str(args.production_queue),
        "acceptance_cases": str(args.acceptance_cases),
        "production_case_count": len(production_cases),
        "acceptance_case_count": len(rows),
        "passed_count": passed,
        "failed_count": len(rows) - passed,
        "all_passed": passed == len(rows),
        "production_forbidden_benchmark_key_count": sum(
            len(find_forbidden_keys(case)) for case in production_cases
        ),
        "cases": rows,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json(args.output_dir / "acceptance_report.private.json", report)
    write_csv(
        args.output_dir / "acceptance_report.private.csv",
        rows,
        (
            "patient_id",
            "answer_organism",
            "expected_in_production_queue",
            "actual_in_production_queue",
            "acceptance_pass",
            "matched_case_ids",
            "actual_tiers",
        ),
    )
    print(f"acceptance_passed={passed}/{len(rows)}")
    print(f"production_forbidden_benchmark_keys={report['production_forbidden_benchmark_key_count']}")
    return 0 if report["all_passed"] and report["production_forbidden_benchmark_key_count"] == 0 else 1


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    production = subparsers.add_parser("production", help="Build answer-blind high/context RAG inputs.")
    production.add_argument("patient_root", type=Path)
    production.add_argument("--merged-suffix", required=True)
    production.add_argument("--label", required=True)
    production.add_argument("--output-dir", type=Path, required=True)

    evaluate = subparsers.add_parser("evaluate", help="Evaluate a production queue with private labels.")
    evaluate.add_argument("--production-queue", type=Path, required=True)
    evaluate.add_argument("--acceptance-cases", type=Path, required=True)
    evaluate.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.command == "production":
        return build_production(args)
    return run_evaluation(args)


if __name__ == "__main__":
    raise SystemExit(main())

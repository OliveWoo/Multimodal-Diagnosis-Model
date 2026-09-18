"""Shadow-mode review/RAG triage from a central rule table.

This tool reads existing missed-candidate review queues and deterministic max
outputs, applies rules/review_triage_rules.json, and writes a comparison report.
It does not modify patient JSONs or final merged outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import mngs_common as mngs  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402


DEFAULT_RULES_PATH = REPO_ROOT / "rules" / "review_triage_rules.json"
DEFAULT_QUEUE_SUFFIX = "mNGS_missed_candidate_review_queue"
DEFAULT_MAX_SUFFIX = "mNGS_max_deterministic_resp_commensal_dominance_guardrail_opt_chosen_full"

DIRECT_REVIEW_TIERS = (
    "review_high_priority",
    "review_context_needed",
    "review_low_specificity",
    "review_omitted_with_reason",
    "omitted_with_reason",
)
TIER_PRIORITY = {
    "review_high_priority": 400,
    "review_context_needed": 300,
    "review_low_specificity": 200,
    "review_omitted_with_reason": 100,
    "omitted_with_reason": 100,
}
CSV_FIELDS = (
    "patient_id",
    "organism_name",
    "canonical_key",
    "genus",
    "shadow_tier",
    "shadow_status",
    "shadow_rag_visible",
    "matched_rule_ids",
    "shadow_reason",
    "existing_review_tier",
    "existing_rag_visible",
    "tier_changed",
    "rag_visibility_changed",
    "rank_priority",
    "reads",
    "reads_tier",
    "reads_percentile",
    "dominance_tier",
    "specimen_class",
    "source_category",
    "support_modules",
    "review_reasons",
)


def resolve_path(value: str | Path | None, *, default: Path | None = None) -> Path | None:
    if value is None:
        return default
    path = Path(value)
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})


def parse_patient_ids(values: Sequence[str] | None) -> set[str]:
    output: set[str] = set()
    for value in values or []:
        text = str(value).strip()
        if not text:
            continue
        output.add(text if text.startswith("NGS_patient_") else f"NGS_patient_{text}")
    return output


def patient_sort_key(patient_dir: Path) -> tuple[int, str]:
    match = re.search(r"NGS_patient_(\d+)", patient_dir.name)
    if not match:
        return (10**9, patient_dir.name)
    return (int(match.group(1)), patient_dir.name)


def patient_dirs(patient_root: Path, patient_ids: set[str]) -> list[Path]:
    dirs = [path for path in patient_root.iterdir() if path.is_dir() and path.name.startswith("NGS_patient_")]
    if patient_ids:
        dirs = [path for path in dirs if path.name.removesuffix("_json") in patient_ids or path.name in patient_ids]
    return sorted(dirs, key=patient_sort_key)


def patient_id_from_dir(patient_dir: Path) -> str:
    return mngs.extract_patient_identifier(patient_dir).replace("NGS_patient_", "")


def path_for_suffix(patient_dir: Path, suffix: str) -> Path:
    base = mngs.extract_patient_identifier(patient_dir)
    return patient_dir / mngs.SUMMARY_OUTPUT_DIR_NAME / f"{base}_{suffix}.json"


def numeric(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def evidence_from_item(item: dict[str, Any]) -> dict[str, Any]:
    evidence = item.get("evidence_snapshot")
    return evidence if isinstance(evidence, dict) else item


def canonical_key(value: Any) -> str:
    return pathogen_names.canonical_key(value)


def display_name(value: Any) -> str:
    return pathogen_names.display_name(value)


def genus_of(value: Any) -> str:
    return pathogen_names.genus_name(value)


def picked_name_sets(deterministic_max: dict[str, Any]) -> tuple[set[str], set[str]]:
    best = deterministic_max.get("best_available_summary") if isinstance(deterministic_max, dict) else {}
    picked = best.get("picked_pathogens") if isinstance(best, dict) else []
    keys: set[str] = set()
    genera: set[str] = set()
    for item in picked if isinstance(picked, list) else []:
        if not isinstance(item, dict):
            continue
        name = item.get("organism_name") or item.get("name")
        key = canonical_key(name)
        genus = genus_of(name)
        if key:
            keys.add(key)
        if genus:
            genera.add(genus)
    return keys, genera


def item_text_blob(item: dict[str, Any]) -> str:
    return json.dumps(item, ensure_ascii=False, sort_keys=True).lower()


def rule_sets(rule_table: dict[str, Any]) -> dict[str, set[str]]:
    sets = rule_table.get("organism_sets")
    if not isinstance(sets, dict):
        return {}
    return {str(name): {canonical_key(value) for value in values or []} for name, values in sets.items()}


def list_values_match_any(actual: Any, expected: Sequence[Any]) -> bool:
    actual_values = {str(value).lower() for value in as_list(actual)}
    expected_values = {str(value).lower() for value in expected}
    return bool(actual_values & expected_values)


def support_modules_only(item: dict[str, Any], expected: Sequence[Any]) -> bool:
    modules = item.get("support_modules") or item.get("non_host_support_modules") or []
    actual = {str(value).lower() for value in as_list(modules)}
    wanted = {str(value).lower() for value in expected}
    return actual == wanted


def support_modules_match_any(item: dict[str, Any], evidence: dict[str, Any], expected: Sequence[Any]) -> bool:
    modules = []
    modules.extend(as_list(item.get("support_modules")))
    modules.extend(as_list(item.get("non_host_support_modules")))
    modules.extend(as_list(evidence.get("support_modules")))
    modules.extend(as_list(evidence.get("non_host_support_modules")))
    return list_values_match_any(modules, expected)


def condition_matches(
    conditions: dict[str, Any],
    *,
    item: dict[str, Any],
    evidence: dict[str, Any],
    key: str,
    genus: str,
    picked_keys: set[str],
    picked_genera: set[str],
    organism_sets: dict[str, set[str]],
) -> bool:
    for field, expected in conditions.items():
        if field == "canonical_in_picked":
            if (key in picked_keys) != bool(expected):
                return False
        elif field == "genus_in_picked":
            if (bool(genus) and genus in picked_genera) != bool(expected):
                return False
        elif field == "canonical_key_equals":
            if key != canonical_key(expected):
                return False
        elif field == "canonical_key_in":
            if key not in {canonical_key(value) for value in as_list(expected)}:
                return False
        elif field == "canonical_key_in_set":
            if key not in organism_sets.get(str(expected), set()):
                return False
        elif field == "genus_equals":
            if genus != str(expected).strip().lower():
                return False
        elif field == "source_category_equals":
            if str(evidence.get("source_category") or item.get("source_category") or "") != str(expected):
                return False
        elif field == "specimen_class_equals":
            if str(evidence.get("specimen_class") or item.get("specimen_class") or "") != str(expected):
                return False
        elif field == "dominance_tier_equals":
            if str(evidence.get("dominance_tier") or item.get("dominance_tier") or "") != str(expected):
                return False
        elif field == "dominance_tier_in":
            actual = str(evidence.get("dominance_tier") or item.get("dominance_tier") or "")
            if actual not in {str(value) for value in as_list(expected)}:
                return False
        elif field == "reads_tier_in":
            actual = str(evidence.get("reads_tier") or item.get("reads_tier") or "")
            if actual not in {str(value) for value in as_list(expected)}:
                return False
        elif field == "confidence_equals":
            if str(item.get("confidence") or "").lower() != str(expected).lower():
                return False
        elif field == "reads_gte":
            if numeric(evidence.get("reads") or item.get("reads")) < numeric(expected):
                return False
        elif field == "reads_lt":
            if numeric(evidence.get("reads") or item.get("reads")) >= numeric(expected):
                return False
        elif field == "reads_percentile_gte":
            if numeric(evidence.get("reads_percentile") or item.get("reads_percentile")) < numeric(expected):
                return False
        elif field == "review_reasons_any":
            if not list_values_match_any(item.get("review_reasons"), as_list(expected)):
                return False
        elif field == "applied_rules_any":
            if not list_values_match_any(item.get("applied_rules"), as_list(expected)):
                return False
        elif field == "guardrail_rule_in":
            if str(item.get("guardrail_rule") or evidence.get("guardrail_rule") or "") not in {str(value) for value in as_list(expected)}:
                return False
        elif field == "is_likely_colonizer_or_background":
            if bool(item.get("is_likely_colonizer_or_background")) != bool(expected):
                return False
        elif field == "support_modules_only":
            if not support_modules_only(item, as_list(expected)):
                return False
        elif field == "support_modules_any":
            if not support_modules_match_any(item, evidence, as_list(expected)):
                return False
        elif field == "text_contains_any":
            blob = item_text_blob(item)
            if not any(str(value).lower() in blob for value in as_list(expected)):
                return False
        else:
            raise ValueError(f"Unsupported rule condition: {field}")
    return True


def classify_item(
    item: dict[str, Any],
    *,
    rule_table: dict[str, Any],
    picked_keys: set[str],
    picked_genera: set[str],
    organism_sets: dict[str, set[str]],
) -> dict[str, Any]:
    name = item.get("organism_name") or item.get("name")
    key = canonical_key(name)
    genus = genus_of(name)
    evidence = evidence_from_item(item)
    rules = sorted(rule_table.get("rules") or [], key=lambda rule: numeric(rule.get("priority")), reverse=True)
    for rule in rules:
        if not isinstance(rule, dict):
            continue
        conditions = rule.get("when") or {}
        if not isinstance(conditions, dict):
            continue
        if condition_matches(
            conditions,
            item=item,
            evidence=evidence,
            key=key,
            genus=genus,
            picked_keys=picked_keys,
            picked_genera=picked_genera,
            organism_sets=organism_sets,
        ):
            tier = str(rule.get("tier"))
            return {
                "tier": tier,
                "status": str(rule.get("status") or "matched"),
                "matched_rule_ids": [str(rule.get("id") or "")],
                "reason": str(rule.get("reason") or ""),
                "canonical_key": key,
                "display_name": display_name(name),
                "genus": genus,
                "rag_visible": tier in set(rule_table.get("rag_visible_tiers") or []),
            }

    defaults = rule_table.get("defaults") if isinstance(rule_table.get("defaults"), dict) else {}
    tier = str(defaults.get("tier") or "review_low_specificity")
    return {
        "tier": tier,
        "status": str(defaults.get("status") or "needs_rule_review"),
        "matched_rule_ids": [str(defaults.get("rule_id") or "default_needs_rule_review")],
        "reason": str(defaults.get("reason") or "No explicit rule matched."),
        "canonical_key": key,
        "display_name": display_name(name),
        "genus": genus,
        "rag_visible": tier in set(rule_table.get("rag_visible_tiers") or []),
    }


def review_tier_index(review: dict[str, Any]) -> dict[str, dict[str, Any]]:
    index: dict[str, dict[str, Any]] = {}
    for tier in DIRECT_REVIEW_TIERS:
        canonical_tier = "review_omitted_with_reason" if tier == "omitted_with_reason" else tier
        items = review.get(tier)
        if not isinstance(items, list):
            continue
        for item in items:
            if not isinstance(item, dict):
                continue
            key = canonical_key(item.get("organism_name") or item.get("name"))
            if not key:
                continue
            previous = index.get(key)
            current_priority = TIER_PRIORITY.get(canonical_tier, 0)
            previous_priority = TIER_PRIORITY.get(previous.get("tier"), 0) if previous else -1
            if previous is None or current_priority > previous_priority:
                index[key] = {
                    "tier": canonical_tier,
                    "confidence": item.get("confidence"),
                    "reason": item.get("review_tier_reason") or item.get("rationale_zh"),
                }
    return index


def joined(value: Any) -> str:
    if isinstance(value, list):
        return "; ".join(str(item) for item in value)
    return "" if value is None else str(value)


def row_for_item(
    *,
    patient_id: str,
    item: dict[str, Any],
    shadow: dict[str, Any],
    existing: dict[str, Any] | None,
    rag_visible_tiers: set[str],
) -> dict[str, Any]:
    evidence = evidence_from_item(item)
    existing_tier = str(existing.get("tier")) if existing else ""
    shadow_tier = str(shadow.get("tier"))
    return {
        "patient_id": patient_id,
        "organism_name": shadow.get("display_name") or item.get("organism_name") or item.get("name"),
        "canonical_key": shadow.get("canonical_key"),
        "genus": shadow.get("genus"),
        "shadow_tier": shadow_tier,
        "shadow_status": shadow.get("status"),
        "shadow_rag_visible": bool(shadow.get("rag_visible")),
        "matched_rule_ids": joined(shadow.get("matched_rule_ids")),
        "shadow_reason": shadow.get("reason"),
        "existing_review_tier": existing_tier,
        "existing_rag_visible": existing_tier in rag_visible_tiers if existing_tier else "",
        "tier_changed": (existing_tier != shadow_tier) if existing_tier else "",
        "rag_visibility_changed": ((existing_tier in rag_visible_tiers) != bool(shadow.get("rag_visible"))) if existing_tier else "",
        "rank_priority": evidence.get("rank_priority"),
        "reads": evidence.get("reads"),
        "reads_tier": evidence.get("reads_tier"),
        "reads_percentile": evidence.get("reads_percentile"),
        "dominance_tier": evidence.get("dominance_tier"),
        "specimen_class": evidence.get("specimen_class"),
        "source_category": evidence.get("source_category"),
        "support_modules": joined(item.get("support_modules") or evidence.get("support_modules")),
        "review_reasons": joined(item.get("review_reasons")),
    }


def process_patient(
    patient_dir: Path,
    *,
    rule_table: dict[str, Any],
    organism_sets: dict[str, set[str]],
    queue_suffix: str,
    max_suffix: str,
    review_suffix: str | None,
) -> tuple[list[dict[str, Any]], dict[str, Any] | None]:
    pid = patient_id_from_dir(patient_dir)
    queue_path = path_for_suffix(patient_dir, queue_suffix)
    max_path = path_for_suffix(patient_dir, max_suffix)
    if not queue_path.exists():
        return [], {"patient_id": pid, "reason": "queue_missing", "path": str(queue_path)}
    if not max_path.exists():
        return [], {"patient_id": pid, "reason": "deterministic_max_missing", "path": str(max_path)}

    queue_payload = load_json(queue_path)
    deterministic_max = load_json(max_path)
    picked_keys, picked_genera = picked_name_sets(deterministic_max)
    existing_index: dict[str, dict[str, Any]] = {}
    if review_suffix:
        review_path = path_for_suffix(patient_dir, review_suffix)
        if review_path.exists():
            existing_index = review_tier_index(load_json(review_path))

    rows: list[dict[str, Any]] = []
    rag_visible_tiers = set(rule_table.get("rag_visible_tiers") or [])
    for item in queue_payload.get("review_queue") or []:
        if not isinstance(item, dict):
            continue
        shadow = classify_item(
            item,
            rule_table=rule_table,
            picked_keys=picked_keys,
            picked_genera=picked_genera,
            organism_sets=organism_sets,
        )
        existing = existing_index.get(str(shadow.get("canonical_key") or ""))
        rows.append(
            row_for_item(
                patient_id=pid,
                item=item,
                shadow=shadow,
                existing=existing,
                rag_visible_tiers=rag_visible_tiers,
            )
        )
    return rows, None


def counter_dict(rows: list[dict[str, Any]], field: str) -> dict[str, int]:
    return dict(sorted(Counter(str(row.get(field) or "") for row in rows).items()))


def build_report(
    rows: list[dict[str, Any]],
    skipped: list[dict[str, Any]],
    *,
    processed_patients: int,
    rule_table: dict[str, Any],
    args: argparse.Namespace,
) -> dict[str, Any]:
    changed = [row for row in rows if row.get("tier_changed") is True]
    rag_visibility_changed = [row for row in rows if row.get("rag_visibility_changed") is True]
    needs_rule_review = [row for row in rows if row.get("shadow_status") == "needs_rule_review"]
    return {
        "report_version": "review_triage_rule_table_shadow_report_v1",
        "mode": "shadow_only",
        "rule_table_schema_version": rule_table.get("schema_version"),
        "inputs": {
            "patient_root": str(resolve_path(args.patient_root)),
            "queue_suffix": args.queue_suffix,
            "deterministic_max_suffix": args.deterministic_max_suffix,
            "existing_review_suffix": args.review_suffix,
            "rules_path": str(resolve_path(args.rules_path, default=DEFAULT_RULES_PATH)),
        },
        "summary": {
            "processed_patients": processed_patients,
            "patients_with_rows": len({row.get("patient_id") for row in rows}),
            "skipped_patients": len(skipped),
            "candidate_rows": len(rows),
            "shadow_tier_counts": counter_dict(rows, "shadow_tier"),
            "shadow_status_counts": counter_dict(rows, "shadow_status"),
            "existing_review_tier_counts": counter_dict(rows, "existing_review_tier"),
            "shadow_rag_visible_count": sum(1 for row in rows if row.get("shadow_rag_visible") is True),
            "existing_rag_visible_count": sum(1 for row in rows if row.get("existing_rag_visible") is True),
            "tier_changed_count": len(changed),
            "rag_visibility_changed_count": len(rag_visibility_changed),
            "needs_rule_review_count": len(needs_rule_review),
        },
        "skipped": skipped,
        "tier_changed_rows": changed,
        "rag_visibility_changed_rows": rag_visibility_changed,
        "needs_rule_review_rows": needs_rule_review,
        "rows": rows,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run shadow-mode LLM-missing/RAG triage from a rule table.")
    parser.add_argument("--patient-root", required=True)
    parser.add_argument("--patient-ids", nargs="*")
    parser.add_argument("--rules-path", default=str(DEFAULT_RULES_PATH))
    parser.add_argument("--queue-suffix", default=DEFAULT_QUEUE_SUFFIX)
    parser.add_argument("--deterministic-max-suffix", default=DEFAULT_MAX_SUFFIX)
    parser.add_argument("--review-suffix", default=None, help="Optional existing LLM review suffix to compare against.")
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-csv", default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    patient_root = resolve_path(args.patient_root)
    rules_path = resolve_path(args.rules_path, default=DEFAULT_RULES_PATH)
    if patient_root is None or rules_path is None:
        raise RuntimeError("patient root and rules path are required")
    rule_table = load_json(rules_path)
    sets = rule_sets(rule_table)
    selected_ids = parse_patient_ids(args.patient_ids)
    rows: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    dirs = patient_dirs(patient_root, selected_ids)
    for patient_dir in dirs:
        patient_rows, skip = process_patient(
            patient_dir,
            rule_table=rule_table,
            organism_sets=sets,
            queue_suffix=args.queue_suffix,
            max_suffix=args.deterministic_max_suffix,
            review_suffix=args.review_suffix,
        )
        rows.extend(patient_rows)
        if skip:
            skipped.append(skip)

    report = build_report(rows, skipped, processed_patients=len(dirs), rule_table=rule_table, args=args)
    output_json = resolve_path(args.output_json)
    if output_json is None:
        raise RuntimeError("output JSON path is required")
    write_json(output_json, report)
    output_csv = resolve_path(args.output_csv) if args.output_csv else output_json.with_suffix(".csv")
    if output_csv is not None:
        write_csv(output_csv, rows)
    print(
        "wrote shadow report "
        f"rows={len(rows)} shadow_rag_visible={report['summary']['shadow_rag_visible_count']} "
        f"needs_rule_review={report['summary']['needs_rule_review_count']} json={output_json} csv={output_csv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

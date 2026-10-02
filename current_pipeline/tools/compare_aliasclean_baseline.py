"""Compare pre-alias-clean and alias-clean merged mNGS max outputs.

This is a baseline/freeze check. It verifies that alias cleanup reduces visible
duplicates without dropping answer hits, while preserving collapsed items in
related_evidence for later RAG review.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from tools import pathogen_normalization as pathogen_names

DEFAULT_ANSWER_CSV = Path("outputs") / "filtered_guardrail_20260804_vs_previous_two_comparison.csv"


@dataclass(frozen=True)
class NameEntry:
    raw: str
    canonical: str
    genus: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build an old-vs-aliasclean comparison for merged mNGS max outputs."
    )
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--before-suffix", required=True)
    parser.add_argument("--after-suffix", required=True)
    parser.add_argument("--answer-csv", type=Path, default=DEFAULT_ANSWER_CSV)
    parser.add_argument(
        "--answer-column",
        default="auto",
        help="Answer column in --answer-csv. Use 'auto' to prefer answer, then answers.",
    )
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    return parser.parse_args(argv)


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def patient_id(patient_dir: Path) -> str:
    return patient_dir.name.removeprefix("NGS_patient_").removesuffix("_json")


def patient_sort_key(patient_dir: Path) -> tuple[int, str]:
    pid = patient_id(patient_dir)
    return (int(pid), pid) if pid.isdigit() else (10**9, pid)


def output_path(patient_dir: Path, suffix: str) -> Path:
    pid = patient_id(patient_dir)
    return patient_dir / "summary_outputs" / f"NGS_patient_{pid}_{suffix}.json"


def load_sections(patient_dir: Path, suffix: str) -> dict[str, list[dict[str, Any]]] | None:
    path = output_path(patient_dir, suffix)
    if not path.is_file():
        return None
    payload = load_json(path)
    deterministic = payload.get("deterministic_max") if isinstance(payload.get("deterministic_max"), dict) else payload
    review = (
        payload.get("llm_missed_candidate_review")
        if isinstance(payload.get("llm_missed_candidate_review"), dict)
        else payload
    )
    best = deterministic.get("best_available_summary") if isinstance(deterministic, dict) else {}
    high_priority = [
        item for item in ((review or {}).get("review_high_priority") or []) if isinstance(item, dict)
    ]
    context_needed = [
        item for item in ((review or {}).get("review_context_needed") or []) if isinstance(item, dict)
    ]
    if high_priority or context_needed:
        possible_items = high_priority
        watch_items = context_needed
    else:
        possible_items = [
            item for item in ((review or {}).get("possible_missed_pathogens") or []) if isinstance(item, dict)
        ]
        watch_items = [item for item in ((review or {}).get("watchlist_candidates") or []) if isinstance(item, dict)]
    return {
        "picked": [item for item in (best.get("picked_pathogens") or []) if isinstance(item, dict)],
        "possible": possible_items,
        "watch": watch_items,
    }


def display_item_name(item: dict[str, Any]) -> str:
    return str(item.get("organism_name") or item.get("name") or "").strip()


def name_entry(name: Any) -> NameEntry | None:
    raw = pathogen_names.clean_display_text(name)
    if not raw or raw == "-":
        return None
    canonical = pathogen_names.canonical_key(raw)
    if not canonical or canonical.isdigit():
        return None
    return NameEntry(raw=raw, canonical=canonical, genus=pathogen_names.genus_name(raw))


def extract_name_entries(value: Any) -> list[NameEntry]:
    if value is None:
        return []
    text = str(value).strip()
    if not text or text in {"-", "NA"}:
        return []
    entries: list[NameEntry] = []
    for part in re.split(r"[;\n\r\uFF1B\u3001]+", text):
        part = part.strip()
        if not part or part in {"-", "NA"}:
            continue
        entry = name_entry(part)
        if entry:
            entries.append(entry)
    return entries


def section_entries(items: list[dict[str, Any]]) -> list[NameEntry]:
    output: list[NameEntry] = []
    for item in items:
        entry = name_entry(display_item_name(item))
        if entry:
            output.append(entry)
    return output


def names_match(output: NameEntry, answer: NameEntry) -> str:
    return pathogen_names.name_match_type(output.raw, answer.raw)


def pairwise_compare(output_entries: list[NameEntry], answer_entries: list[NameEntry]) -> dict[str, Any]:
    unmatched_answers = set(range(len(answer_entries)))
    matched: list[dict[str, str]] = []
    output_only: list[str] = []
    for output in output_entries:
        best_idx = None
        best_type = ""
        for idx in list(unmatched_answers):
            match_type = names_match(output, answer_entries[idx])
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
            output_only.append(output.raw)
            continue
        unmatched_answers.remove(best_idx)
        matched.append(
            {
                "output": output.raw,
                "answer": answer_entries[best_idx].raw,
                "match_type": best_type,
            }
        )
    return {
        "matched": matched,
        "output_only": output_only,
        "answer_only": [answer_entries[idx].raw for idx in sorted(unmatched_answers)],
    }


def ratio_text(matched: int, total: int) -> str:
    return f"{matched}/{total}" if total else ""


def safe_ratio(matched: int, total: int) -> float:
    return matched / total if total else 0.0


def related_evidence_pairs(items: list[dict[str, Any]]) -> list[str]:
    pairs: list[str] = []
    for item in items:
        retained = display_item_name(item)
        related = item.get("related_evidence") or []
        if not isinstance(related, list):
            continue
        for rel in related:
            if not isinstance(rel, dict):
                continue
            original = str(rel.get("original_name") or "").strip()
            relationship = str(rel.get("relationship") or "").strip()
            source = str(rel.get("source_section") or "").strip()
            if original:
                suffix = f" [{relationship}" + (f", {source}]" if source else "]")
                pairs.append(f"{retained} <= {original}{suffix}")
    return pairs


def duplicate_group_details(items: list[dict[str, Any]]) -> list[str]:
    groups: dict[str, list[str]] = defaultdict(list)
    for item in items:
        name = display_item_name(item)
        if not name:
            continue
        key = pathogen_names.representative_group_key(name)
        groups[key].append(name)
    details: list[str] = []
    for key, names in sorted(groups.items()):
        if len(names) <= 1:
            continue
        counts = Counter(names)
        rendered = []
        for name, count in sorted(counts.items()):
            rendered.append(f"{name} x{count}" if count > 1 else name)
        details.append(f"{key}: " + " | ".join(rendered))
    return details


def multiset_removed_added(before: list[dict[str, Any]], after: list[dict[str, Any]]) -> tuple[list[str], list[str]]:
    before_names = Counter(display_item_name(item) for item in before if display_item_name(item))
    after_names = Counter(display_item_name(item) for item in after if display_item_name(item))
    removed: list[str] = []
    added: list[str] = []
    for name, count in sorted((before_names - after_names).items()):
        removed.extend([name] * count)
    for name, count in sorted((after_names - before_names).items()):
        added.extend([name] * count)
    return removed, added


def resolve_answer_column(fieldnames: Sequence[str], requested: str) -> str:
    if requested != "auto":
        if requested not in fieldnames:
            raise ValueError(f"Answer column {requested!r} not found in CSV fields: {fieldnames}")
        return requested
    for candidate in ("answer", "answers"):
        if candidate in fieldnames:
            return candidate
    raise ValueError(f"Cannot auto-detect answer column from CSV fields: {fieldnames}")


def load_answers(path: Path, answer_column: str) -> tuple[dict[str, list[NameEntry]], dict[str, str], str, int]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return {}, {}, answer_column, 0
        column = resolve_answer_column(reader.fieldnames, answer_column)
        answers: dict[str, list[NameEntry]] = {}
        raw_answers: dict[str, str] = {}
        row_count = 0
        for row in reader:
            row_count += 1
            pid = str(row.get("patient_id") or "").strip()
            if not pid:
                continue
            raw = str(row.get(column) or "").strip()
            raw_answers[pid] = raw
            entries = extract_name_entries(raw)
            if entries:
                answers[pid] = entries
        return answers, raw_answers, column, row_count


def totals_template() -> dict[str, int]:
    return {
        "picked_output": 0,
        "picked_matched": 0,
        "picked_exact_or_alias": 0,
        "picked_approved_group_member_match": 0,
        "picked_genus_relaxed": 0,
        "combined_output": 0,
        "combined_matched": 0,
        "combined_exact_or_alias": 0,
        "combined_approved_group_member_match": 0,
        "combined_genus_relaxed": 0,
        "answer_organisms": 0,
    }


def add_match_totals(
    totals: dict[str, int],
    prefix: str,
    output_entries: list[NameEntry],
    answers: list[NameEntry],
    comparison: dict[str, Any],
) -> None:
    totals[f"{prefix}_output"] += len(output_entries)
    totals[f"{prefix}_matched"] += len(comparison["matched"])
    totals[f"{prefix}_exact_or_alias"] += sum(1 for item in comparison["matched"] if item["match_type"] == "exact_or_alias")
    totals[f"{prefix}_approved_group_member_match"] += sum(
        1 for item in comparison["matched"] if item["match_type"] == "approved_group_member_match"
    )
    totals[f"{prefix}_genus_relaxed"] += sum(1 for item in comparison["matched"] if item["match_type"] == "genus_relaxed")
    if prefix == "picked":
        totals["answer_organisms"] += len(answers)


def summarize_totals(totals: dict[str, int]) -> dict[str, Any]:
    answer_total = totals["answer_organisms"]
    return {
        **totals,
        "picked_precision": safe_ratio(totals["picked_matched"], totals["picked_output"]),
        "picked_recall": safe_ratio(totals["picked_matched"], answer_total),
        "combined_precision": safe_ratio(totals["combined_matched"], totals["combined_output"]),
        "combined_recall": safe_ratio(totals["combined_matched"], answer_total),
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    answers, raw_answers, answer_column, answer_row_count = load_answers(args.answer_csv, args.answer_column)
    rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {
        "label": args.label,
        "patient_root": str(args.patient_root),
        "before_suffix": args.before_suffix,
        "after_suffix": args.after_suffix,
        "answer_csv": str(args.answer_csv),
        "answer_column": answer_column,
        "answer_csv_row_count": answer_row_count,
        "answer_patient_count": len(answers),
        "answer_organism_count": sum(len(value) for value in answers.values()),
        "missing_before": [],
        "missing_after": [],
        "all_patients": {
            "before_picked": 0,
            "after_picked": 0,
            "before_possible": 0,
            "after_possible": 0,
            "before_watch": 0,
            "after_watch": 0,
            "before_combined": 0,
            "after_combined": 0,
            "before_visible_duplicate_groups": 0,
            "after_visible_duplicate_groups": 0,
            "before_visible_duplicate_extra_items": 0,
            "after_visible_duplicate_extra_items": 0,
            "after_related_evidence_count": 0,
        },
        "answer_patients": {
            "before": totals_template(),
            "after": totals_template(),
        },
    }

    for patient_dir in sorted(args.patient_root.glob("NGS_patient_*_json"), key=patient_sort_key):
        pid = patient_id(patient_dir)
        before = load_sections(patient_dir, args.before_suffix)
        after = load_sections(patient_dir, args.after_suffix)
        if before is None:
            summary["missing_before"].append(pid)
            before = {"picked": [], "possible": [], "watch": []}
        if after is None:
            summary["missing_after"].append(pid)
            after = {"picked": [], "possible": [], "watch": []}

        before_combined = before["picked"] + before["possible"] + before["watch"]
        after_combined = after["picked"] + after["possible"] + after["watch"]
        before_duplicates = duplicate_group_details(before_combined)
        after_duplicates = duplicate_group_details(after_combined)
        after_related = related_evidence_pairs(after_combined)
        removed_visible, added_visible = multiset_removed_added(before_combined, after_combined)

        all_summary = summary["all_patients"]
        all_summary["before_picked"] += len(before["picked"])
        all_summary["after_picked"] += len(after["picked"])
        all_summary["before_possible"] += len(before["possible"])
        all_summary["after_possible"] += len(after["possible"])
        all_summary["before_watch"] += len(before["watch"])
        all_summary["after_watch"] += len(after["watch"])
        all_summary["before_combined"] += len(before_combined)
        all_summary["after_combined"] += len(after_combined)
        all_summary["before_visible_duplicate_groups"] += len(before_duplicates)
        all_summary["after_visible_duplicate_groups"] += len(after_duplicates)
        all_summary["before_visible_duplicate_extra_items"] += sum(
            max(0, detail.count("|")) for detail in before_duplicates
        )
        all_summary["after_visible_duplicate_extra_items"] += sum(max(0, detail.count("|")) for detail in after_duplicates)
        all_summary["after_related_evidence_count"] += len(after_related)

        answer_entries = answers.get(pid, [])
        before_picked_entries = section_entries(before["picked"])
        after_picked_entries = section_entries(after["picked"])
        before_combined_entries = section_entries(before_combined)
        after_combined_entries = section_entries(after_combined)

        before_picked_cmp = pairwise_compare(before_picked_entries, answer_entries)
        after_picked_cmp = pairwise_compare(after_picked_entries, answer_entries)
        before_combined_cmp = pairwise_compare(before_combined_entries, answer_entries)
        after_combined_cmp = pairwise_compare(after_combined_entries, answer_entries)

        if answer_entries:
            add_match_totals(summary["answer_patients"]["before"], "picked", before_picked_entries, answer_entries, before_picked_cmp)
            add_match_totals(summary["answer_patients"]["after"], "picked", after_picked_entries, answer_entries, after_picked_cmp)
            add_match_totals(summary["answer_patients"]["before"], "combined", before_combined_entries, answer_entries, before_combined_cmp)
            add_match_totals(summary["answer_patients"]["after"], "combined", after_combined_entries, answer_entries, after_combined_cmp)

        rows.append(
            {
                "label": args.label,
                "patient_id": pid,
                "answer": raw_answers.get(pid, ""),
                "answer_count": len(answer_entries),
                "before_picked_count": len(before["picked"]),
                "after_picked_count": len(after["picked"]),
                "picked_count_delta": len(after["picked"]) - len(before["picked"]),
                "before_possible_count": len(before["possible"]),
                "after_possible_count": len(after["possible"]),
                "possible_count_delta": len(after["possible"]) - len(before["possible"]),
                "before_watch_count": len(before["watch"]),
                "after_watch_count": len(after["watch"]),
                "watch_count_delta": len(after["watch"]) - len(before["watch"]),
                "before_combined_count": len(before_combined),
                "after_combined_count": len(after_combined),
                "combined_count_delta": len(after_combined) - len(before_combined),
                "before_picked_ratio": ratio_text(len(before_picked_cmp["matched"]), len(answer_entries)),
                "after_picked_ratio": ratio_text(len(after_picked_cmp["matched"]), len(answer_entries)),
                "picked_hit_delta": len(after_picked_cmp["matched"]) - len(before_picked_cmp["matched"]),
                "before_combined_ratio": ratio_text(len(before_combined_cmp["matched"]), len(answer_entries)),
                "after_combined_ratio": ratio_text(len(after_combined_cmp["matched"]), len(answer_entries)),
                "combined_hit_delta": len(after_combined_cmp["matched"]) - len(before_combined_cmp["matched"]),
                "before_visible_duplicate_groups": " ; ".join(before_duplicates),
                "after_visible_duplicate_groups": " ; ".join(after_duplicates),
                "duplicate_group_delta": len(after_duplicates) - len(before_duplicates),
                "after_related_evidence_count": len(after_related),
                "after_related_evidence_pairs": " ; ".join(after_related),
                "removed_visible_outputs": " ; ".join(removed_visible),
                "added_visible_outputs": " ; ".join(added_visible),
                "before_picked_missing_answers": " ; ".join(before_picked_cmp["answer_only"]),
                "after_picked_missing_answers": " ; ".join(after_picked_cmp["answer_only"]),
                "before_combined_missing_answers": " ; ".join(before_combined_cmp["answer_only"]),
                "after_combined_missing_answers": " ; ".join(after_combined_cmp["answer_only"]),
            }
        )

    summary["all_patients"]["picked_delta"] = (
        summary["all_patients"]["after_picked"] - summary["all_patients"]["before_picked"]
    )
    summary["all_patients"]["possible_delta"] = (
        summary["all_patients"]["after_possible"] - summary["all_patients"]["before_possible"]
    )
    summary["all_patients"]["watch_delta"] = summary["all_patients"]["after_watch"] - summary["all_patients"]["before_watch"]
    summary["all_patients"]["combined_delta"] = (
        summary["all_patients"]["after_combined"] - summary["all_patients"]["before_combined"]
    )
    summary["answer_patients"]["before"] = summarize_totals(summary["answer_patients"]["before"])
    summary["answer_patients"]["after"] = summarize_totals(summary["answer_patients"]["after"])
    summary["answer_patients"]["picked_hit_delta"] = (
        summary["answer_patients"]["after"]["picked_matched"]
        - summary["answer_patients"]["before"]["picked_matched"]
    )
    summary["answer_patients"]["combined_hit_delta"] = (
        summary["answer_patients"]["after"]["combined_matched"]
        - summary["answer_patients"]["before"]["combined_matched"]
    )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else [])
        writer.writeheader()
        writer.writerows(rows)
    write_json(args.output_json, summary)
    print(f"wrote {args.output_csv}")
    print(f"wrote {args.output_json}")
    print(json.dumps(summary["all_patients"], ensure_ascii=False, indent=2))
    print(json.dumps(summary["answer_patients"], ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

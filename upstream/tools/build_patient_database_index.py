from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import pathogen_normalization as pathogen_names  # noqa: E402

PATIENT_DIR_RE = re.compile(r"^NGS_patient_(\d+)_json$")
NAME_SANITIZER = re.compile(r"[^a-z0-9]+")
SPLIT_RE = re.compile(r"[;；,，\n\r]+")

DEFAULT_MERGE_SUFFIX = "mNGS_max_deterministic_pulmonary_opt_with_missed_candidate_review"
DEFAULT_FILTERED_MERGE_SUFFIX = "mNGS_max_deterministic_pulmonary_opt_guardrail_20260804_aliasclean_20260805_with_missed_candidate_review"
DEFAULT_KH_MERGE_SUFFIX = "mNGS_max_deterministic_pulmonary_opt_guardrail_20260804_aliasclean_20260805_direct_triage"
DEFAULT_TWO_HOSPITALS_MERGE_SUFFIX = "mNGS_max_deterministic_pulmonary_opt_with_missed_candidate_review"
DEFAULT_FINAL_SUMMARY_SUFFIXES = [
    "final_summary_with_filmarray_deterministic",
    "final_summary_with_underlying_with_filmarray",
    "final_summary",
]

NON_ANSWER = {"", "-", "none", "na", "no", "nopathogen", "nopatogen", "notinfectioncase", "negative"}

COLUMNS = [
    "patient_id",
    "specimen_id",
    "dataset",
    "specimen_type",
    "data_richness",
    "available_modules",
    "missing_modules",
    "hospital_detected_organisms",
    "mngs_raw_count",
    "picked_count",
    "possible_missed_count",
    "watchlist_count",
    "picked_pathogens",
    "possible_missed_pathogens",
    "watchlist_pathogens",
    "hospital_answer",
    "picked_answer_hit_ratio",
    "combined_answer_hit_ratio",
    "needs_review",
    "review_reason",
    "final_json_path",
    "patient_folder_path",
]


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def list_patient_dirs(root: Path) -> list[tuple[int, Path]]:
    result = []
    for path in root.iterdir() if root.exists() else []:
        if not path.is_dir():
            continue
        match = PATIENT_DIR_RE.match(path.name)
        if match:
            result.append((int(match.group(1)), path))
    return sorted(result)


def first_existing_summary(patient_dir: Path, patient_id: int) -> Path | None:
    summary_dir = patient_dir / "summary_outputs"
    for suffix in DEFAULT_FINAL_SUMMARY_SUFFIXES:
        path = summary_dir / f"NGS_patient_{patient_id}_{suffix}.json"
        if path.exists():
            return path
    return None


def merge_path(patient_dir: Path, patient_id: int, suffix: str) -> Path:
    return patient_dir / "summary_outputs" / f"NGS_patient_{patient_id}_{suffix}.json"


def canonical_name(value: Any) -> str:
    text = str(value or "").strip().lower()
    if not text or text == "-":
        return ""
    norm = pathogen_names.raw_key(text)
    if norm in NON_ANSWER:
        return ""
    return pathogen_names.canonical_key(text)


def genus_name(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text == "-":
        return ""
    return pathogen_names.genus_name(text)


def split_names(value: Any) -> list[str]:
    if value is None:
        return []
    raw = str(value).strip()
    if not raw or raw == "-":
        return []
    output = []
    for part in SPLIT_RE.split(raw):
        name = part.strip()
        if not name or canonical_name(name) in NON_ANSWER:
            continue
        if canonical_name(name):
            output.append(name)
    return output


def names_match(answer: str, candidate: str) -> bool:
    return bool(pathogen_names.name_match_type(answer, candidate))


def ratio_text(hit_count: int, answer_count: int) -> str:
    if answer_count <= 0:
        return "NA"
    return f"{hit_count}/{answer_count}"


def load_answer_lookup(paths: list[tuple[str, Path]]) -> dict[tuple[str, int], list[str]]:
    lookup: dict[tuple[str, int], list[str]] = {}
    for dataset, path in paths:
        if not path.exists():
            print(f"warning: answer CSV not found for {dataset}: {path}", file=sys.stderr)
            continue
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                try:
                    pid = int(str(row.get("patient_id") or "").strip())
                except ValueError:
                    continue
                answers = split_names(row.get("answer") or row.get("answers") or row.get("hospital_answer") or "")
                lookup[(dataset, pid)] = answers
    return lookup


def dedupe_texts(values: list[str]) -> list[str]:
    seen = set()
    output = []
    for value in values:
        text = str(value or "").strip()
        if not text or text == "-":
            continue
        key = text.lower()
        if key in seen:
            continue
        seen.add(key)
        output.append(text)
    return output


def unique_join(values: list[str]) -> str:
    output = dedupe_texts(values)
    return "; ".join(output) if output else ""


def multiline_join(values: list[str]) -> str:
    output = dedupe_texts(values)
    return "\n".join(output) if output else ""


def iter_raw_mngs_items(patient_dir: Path, patient_id: int) -> tuple[bool, list[dict[str, Any]], list[str], list[str]]:
    path = patient_dir / f"NGS_patient_{patient_id}_all_RK_NTC_microbes.json"
    if not path.exists():
        return False, [], [], []
    try:
        payload = read_json(path)
    except Exception:
        return True, [], [], []
    records = payload if isinstance(payload, list) else payload.get("records", []) if isinstance(payload, dict) else []
    specimen_ids = []
    specimen_types = []
    items: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        specimen_ids.append(str(record.get("specimen_code") or record.get("seq_id") or "").strip())
        specimen_types.append(str(record.get("specimen_site") or record.get("specimen_type") or "").strip())
        pathogens = record.get("pathogens") or {}
        if isinstance(pathogens, dict):
            groups = pathogens.values()
        elif isinstance(pathogens, list):
            groups = [pathogens]
        else:
            groups = []
        for group in groups:
            if not isinstance(group, list):
                continue
            for item in group:
                if isinstance(item, dict):
                    items.append(item)
    return True, items, specimen_ids, specimen_types


def raw_mngs_info(patient_dir: Path, patient_id: int) -> tuple[bool, int, str, str]:
    exists, items, specimen_ids, specimen_types = iter_raw_mngs_items(patient_dir, patient_id)
    names = [str(item.get("name") or item.get("organism_name") or "") for item in items if isinstance(item, dict)]
    return exists, len({canonical_name(n) or n.lower() for n in names if n}), unique_join(specimen_ids), unique_join(specimen_types)


def parse_reads(value: Any) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(float(str(value).replace(",", "").strip()))
    except (TypeError, ValueError):
        return None


def item_reads(item: dict[str, Any]) -> int | None:
    for key in ["reads", "mngs_reads", "read_count", "mNGS Read", "mNGS_reads"]:
        reads = parse_reads(item.get(key))
        if reads is not None:
            return reads
    snapshot = item.get("evidence_snapshot")
    if isinstance(snapshot, dict):
        for key in ["reads", "mngs_reads", "read_count"]:
            reads = parse_reads(snapshot.get(key))
            if reads is not None:
                return reads
    return None


def raw_reads_by_name(patient_dir: Path, patient_id: int) -> dict[str, int]:
    _, items, _, _ = iter_raw_mngs_items(patient_dir, patient_id)
    output: dict[str, int] = {}
    for item in items:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name") or item.get("organism_name") or "").strip()
        key = canonical_name(name)
        reads = item_reads(item)
        if not key or reads is None:
            continue
        output[key] = max(output.get(key, 0), reads)
    return output


def format_name_with_reads(name: str, reads: int | None) -> str:
    name = str(name or "").strip()
    if not name:
        return ""
    if reads is None:
        return name
    return f"{name} ({reads})"


def format_candidate_items(items: list[dict[str, Any]]) -> list[str]:
    output = []
    for item in items:
        name = str(item.get("organism_name") or item.get("name") or "").strip()
        if name:
            output.append(format_name_with_reads(name, item_reads(item)))
    return output


def format_answers_with_reads(answers: list[str], candidate_items: list[dict[str, Any]], raw_reads: dict[str, int]) -> list[str]:
    output = []
    for answer in answers:
        reads_values: list[int] = []
        for item in candidate_items:
            name = str(item.get("organism_name") or item.get("name") or "").strip()
            reads = item_reads(item)
            if reads is not None and names_match(answer, name):
                reads_values.append(reads)
        key = canonical_name(answer)
        if key in raw_reads:
            reads_values.append(raw_reads[key])
        reads = max(reads_values) if reads_values else 0
        output.append(format_name_with_reads(answer, reads))
    return output


def module_available(summary: dict[str, Any], module: str) -> bool:
    availability = summary.get("module_availability") or {}
    status = str(availability.get(module) or "").strip().lower()
    if status == "available":
        return True
    return False


def file_exists(patient_dir: Path, patient_id: int, suffix: str) -> bool:
    return (patient_dir / f"NGS_patient_{patient_id}_{suffix}.json").exists() or (
        patient_dir / "agent_outputs" / f"NGS_patient_{patient_id}_{suffix}_agent.json"
    ).exists()


def hospital_evidence(summary: dict[str, Any]) -> tuple[list[str], set[str]]:
    evidence = summary.get("hospital_organism_evidence") or []
    organisms = []
    modules_with_organisms: set[str] = set()
    if not isinstance(evidence, list):
        return organisms, modules_with_organisms
    for item in evidence:
        if not isinstance(item, dict):
            continue
        name = str(item.get("organism_name") or "").strip()
        if not name:
            continue
        modules = [str(m) for m in item.get("evidence_modules") or [] if m]
        modules_with_organisms.update(modules)
        level = str(item.get("best_hospital_level") or "").strip()
        detail = "+".join(modules) if modules else "hospital"
        if level:
            detail = f"{detail}; {level}"
        organisms.append(f"{name} [{detail}]")
    return organisms, modules_with_organisms


def available_and_missing(dataset: str, patient_dir: Path, patient_id: int, summary: dict[str, Any], mngs_exists: bool, modules_with_organisms: set[str]) -> tuple[list[str], list[str], int, int]:
    available = []
    if mngs_exists:
        available.append("mNGS")

    common_expected = ["culture", "cbc_other_lab", "image", "underlying", "admission_diagnosis"]
    if dataset in {"filtered", "KH"}:
        expected = ["mNGS", "culture", "filmarray_gmtest", "cbc_other_lab", "image", "underlying", "admission_diagnosis"]
    else:
        expected = ["mNGS", "culture", "molecular_microbiology", "cbc_other_lab", "image", "underlying", "admission_diagnosis"]

    for module in ["culture", "filmarray_gmtest", "molecular_microbiology", "cbc_other_lab", "image"]:
        if module_available(summary, module) or module in modules_with_organisms:
            available.append(module)
    if file_exists(patient_dir, patient_id, "underlying"):
        available.append("underlying")
    if file_exists(patient_dir, patient_id, "admission_diagnosis"):
        available.append("admission_diagnosis")

    # For two_hospitals, FilmArray/GM is optional but should appear if it has real organism evidence.
    if dataset == "two_hospitals" and "filmarray_gmtest" in modules_with_organisms and "filmarray_gmtest" not in available:
        available.append("filmarray_gmtest")

    available = list(dict.fromkeys(available))
    missing = [m for m in expected if m not in available]
    hospital_count = len([m for m in available if m != "mNGS"])
    # data_richness uses effective strong evidence modules, not mere file/module presence.
    strong_count = len([m for m in ["culture", "filmarray_gmtest", "molecular_microbiology"] if m in modules_with_organisms])
    return available, missing, strong_count, hospital_count


def data_richness(mngs_exists: bool, strong_count: int, hospital_count: int) -> str:
    if not mngs_exists:
        return "mNGS-missing"
    if hospital_count == 0:
        return "mNGS-only"
    if strong_count >= 2:
        return "rich"
    if strong_count >= 1 and hospital_count >= 3:
        return "moderate"
    return "sparse"


def extract_merge_items(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    if not path.exists():
        return [], [], []
    try:
        payload = read_json(path)
    except Exception:
        return [], [], []
    deterministic = payload.get("deterministic_max") or payload
    best = deterministic.get("best_available_summary") or {}
    picked = [item for item in best.get("picked_pathogens") or [] if isinstance(item, dict) and str(item.get("organism_name") or "").strip()]
    review = payload.get("llm_missed_candidate_review") or {}
    high_priority = [item for item in review.get("review_high_priority") or [] if isinstance(item, dict) and str(item.get("organism_name") or "").strip()]
    context_needed = [item for item in review.get("review_context_needed") or [] if isinstance(item, dict) and str(item.get("organism_name") or "").strip()]
    if high_priority or context_needed:
        possible = high_priority
        watch = context_needed
    else:
        possible = [item for item in review.get("possible_missed_pathogens") or [] if isinstance(item, dict) and str(item.get("organism_name") or "").strip()]
        watch = [item for item in review.get("watchlist_candidates") or [] if isinstance(item, dict) and str(item.get("organism_name") or "").strip()]
    return picked, possible, watch


def item_names(items: list[dict[str, Any]]) -> list[str]:
    return [str(item.get("organism_name") or item.get("name") or "").strip() for item in items if str(item.get("organism_name") or item.get("name") or "").strip()]


def count_answer_hits(answers: list[str], candidates: list[str]) -> int:
    hits = 0
    used = set()
    for answer in answers:
        for idx, candidate in enumerate(candidates):
            if idx in used:
                continue
            if names_match(answer, candidate):
                hits += 1
                used.add(idx)
                break
    return hits


def build_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    answer_lookup = load_answer_lookup([
        ("filtered", args.filtered_answer_csv),
        ("KH", args.kh_answer_csv),
        ("two_hospitals", args.two_hospitals_answer_csv),
    ])
    rows = []
    roots = [
        ("filtered", args.filtered_root, args.filtered_merge_suffix or args.merge_suffix),
        ("KH", args.kh_root, args.kh_merge_suffix or args.merge_suffix),
        ("two_hospitals", args.two_hospitals_root, args.two_hospitals_merge_suffix or args.merge_suffix),
    ]
    for dataset, root, merge_suffix in roots:
        for patient_id, patient_dir in list_patient_dirs(root):
            final_summary_path = first_existing_summary(patient_dir, patient_id)
            summary = read_json(final_summary_path) if final_summary_path and final_summary_path.exists() else {}
            mngs_exists, raw_count, specimen_id, specimen_type = raw_mngs_info(patient_dir, patient_id)
            hospital_orgs, modules_with_organisms = hospital_evidence(summary)
            available, missing, strong_count, hospital_count = available_and_missing(
                dataset, patient_dir, patient_id, summary, mngs_exists, modules_with_organisms
            )
            final_path = merge_path(patient_dir, patient_id, merge_suffix)
            picked_items, possible_items, watch_items = extract_merge_items(final_path)
            picked = item_names(picked_items)
            possible = item_names(possible_items)
            watch = item_names(watch_items)
            answers = answer_lookup.get((dataset, patient_id), [])
            picked_hits = count_answer_hits(answers, picked)
            combined_hits = count_answer_hits(answers, picked + possible + watch)
            raw_reads = raw_reads_by_name(patient_dir, patient_id)
            row = {
                "patient_id": patient_id,
                "specimen_id": specimen_id,
                "dataset": dataset,
                "specimen_type": specimen_type,
                "data_richness": data_richness(mngs_exists, strong_count, hospital_count),
                "available_modules": unique_join(available),
                "missing_modules": unique_join(missing),
                "hospital_detected_organisms": multiline_join(hospital_orgs),
                "mngs_raw_count": raw_count,
                "picked_count": len(picked),
                "possible_missed_count": len(possible),
                "watchlist_count": len(watch),
                "picked_pathogens": multiline_join(format_candidate_items(picked_items)),
                "possible_missed_pathogens": multiline_join(format_candidate_items(possible_items)),
                "watchlist_pathogens": multiline_join(format_candidate_items(watch_items)),
                "hospital_answer": multiline_join(format_answers_with_reads(answers, picked_items + possible_items + watch_items, raw_reads)),
                "picked_answer_hit_ratio": ratio_text(picked_hits, len(answers)),
                "combined_answer_hit_ratio": ratio_text(combined_hits, len(answers)),
                "needs_review": "",
                "review_reason": "",
                "final_json_path": str(final_path.resolve()) if final_path.exists() else "",
                "patient_folder_path": str(patient_dir.resolve()),
            }
            rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a patient-level CSV/JSON database index for mNGS review outputs.")
    parser.add_argument("--filtered-root", type=Path, default=Path(r"D:\CSIE_PROJECT\outputs\patient_info_filtered"))
    parser.add_argument("--kh-root", type=Path, default=Path(r"D:\CSIE_PROJECT\outputs\patient_info_KH_0726_2Days"))
    parser.add_argument("--two-hospitals-root", type=Path, default=Path(r"D:\CSIE_PROJECT\outputs\patient_info_two_hospitals_0611"))
    parser.add_argument("--filtered-answer-csv", type=Path, default=Path(r"D:\CSIE_PROJECT\outputs\kmuh_filtered_latest_pulmonary_opt_answer_comparison_after_queue_review_update.csv"))
    parser.add_argument("--kh-answer-csv", type=Path, default=Path(r"D:\CSIE_PROJECT\outputs\runs\legacy_reports\filtered\20260804_guardrail\filtered_guardrail_20260804_vs_previous_two_comparison.csv"))
    parser.add_argument(
        "--two-hospitals-answer-csv",
        type=Path,
        default=Path(
            r"D:\CSIE_PROJECT\outputs\runs\legacy_reports\two_hospitals_0611\answer_match\two_hospitals_0611_all_patients_answer_match_clean_latest.csv"
        ),
    )
    parser.add_argument("--merge-suffix", default=DEFAULT_MERGE_SUFFIX)
    parser.add_argument("--filtered-merge-suffix", default=DEFAULT_FILTERED_MERGE_SUFFIX)
    parser.add_argument("--kh-merge-suffix", default=DEFAULT_KH_MERGE_SUFFIX)
    parser.add_argument("--two-hospitals-merge-suffix", default=DEFAULT_TWO_HOSPITALS_MERGE_SUFFIX)
    parser.add_argument("--output-dir", type=Path, default=Path(r"D:\CSIE_PROJECT\outputs\patient_database_latest"))
    args = parser.parse_args()

    rows = build_rows(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "patient_database_index.csv"
    json_path = args.output_dir / "patient_database_index.json"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    write_json(json_path, {"rows": rows, "row_count": len(rows), "columns": COLUMNS})
    print(f"wrote {len(rows)} rows")
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()


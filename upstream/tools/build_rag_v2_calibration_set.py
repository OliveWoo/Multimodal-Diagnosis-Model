"""Build a private positive/negative calibration subset for RAG v2."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import build_rag_v2_inputs as rag_inputs  # noqa: E402
from tools import compare_aliasclean_baseline as compare  # noqa: E402

DEFAULT_TARGET_NEGATIVE_COUNTS = {
    "Candida/yeast": 3,
    "Viral": 3,
    "Skin/airway colonizer Gram-positive": 3,
    "Strict anaerobe/aspiration flora": 2,
    "Bacterial": 2,
    "Enterobacterales/hospital GNB": 1,
    "Respiratory virus": 1,
    "Other/unclassified": 1,
}


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


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
        writer.writerows(rows)


def dataset_key(case: dict[str, Any]) -> str:
    label = str(case.get("dataset_label") or "").lower()
    if "two_hospitals" in label or "two-hospital" in label:
        return "two_hospitals"
    if "kh" in label:
        return "KH"
    if "filtered" in label:
        return "filtered"
    return "unknown"


def answer_lookup(args: argparse.Namespace) -> dict[tuple[str, str], list[compare.NameEntry]]:
    lookup: dict[tuple[str, str], list[compare.NameEntry]] = {}
    for dataset, path in (
        ("KH", args.kh_answer_csv),
        ("two_hospitals", args.two_hospitals_answer_csv),
        ("filtered", args.filtered_answer_csv),
    ):
        if path is None or not path.exists():
            continue
        answers, _, _, _ = compare.load_answers(path, "answers")
        for patient_id, entries in answers.items():
            lookup[(dataset, str(patient_id))] = entries
    return lookup


def case_entry(case: dict[str, Any]) -> compare.NameEntry | None:
    organism = case.get("organism")
    name = organism.get("display_name") if isinstance(organism, dict) else ""
    return compare.name_entry(name)


def case_matches_answers(case: dict[str, Any], answers: list[compare.NameEntry]) -> tuple[bool, str, str]:
    entry = case_entry(case)
    if entry is None:
        return False, "", ""
    for answer in answers:
        match_type = compare.names_match(entry, answer)
        if match_type:
            return True, answer.raw, match_type
    return False, "", ""


def find_positive_case(
    cases: Sequence[dict[str, Any]],
    expected: dict[str, Any],
    dataset_default: str,
) -> dict[str, Any] | None:
    patient_id = str(expected.get("patient_id") or "")
    expected_entry = compare.name_entry(expected.get("answer_organism"))
    dataset = str(expected.get("dataset") or dataset_default)
    if expected_entry is None:
        return None
    matches: list[dict[str, Any]] = []
    for case in cases:
        if str(case.get("patient_id") or "") != patient_id:
            continue
        if dataset_key(case) != dataset:
            continue
        entry = case_entry(case)
        if entry and compare.names_match(entry, expected_entry):
            matches.append(case)
    if not matches:
        return None
    tier_priority = {"review_high_priority": 2, "review_context_needed": 1}
    return max(
        matches,
        key=lambda case: tier_priority.get(str((case.get("current_model_state") or {}).get("source_tier")), 0),
    )


def negative_sort_key(case: dict[str, Any]) -> tuple[int, int, int, str]:
    current = case.get("current_model_state") if isinstance(case.get("current_model_state"), dict) else {}
    local = case.get("local_evidence") if isinstance(case.get("local_evidence"), dict) else {}
    mngs = local.get("mngs_evidence") if isinstance(local.get("mngs_evidence"), dict) else {}
    hospital = local.get("same_organism_hospital_evidence")
    negatives = local.get("negative_evidence")
    tier = str(current.get("source_tier") or "")
    tier_score = 2 if tier == "review_high_priority" else 1
    hospital_score = 1 if isinstance(hospital, list) and hospital else 0
    negative_count = len(negatives) if isinstance(negatives, list) else 0
    reads = int(float(mngs.get("reads") or 0))
    return (-tier_score, -hospital_score, -negative_count, f"{reads:012d}_{case.get('case_id')}")


def label_row(
    case: dict[str, Any],
    *,
    label_type: str,
    benchmark_answer: str,
    match_type: str,
    selection_reason: str,
    patient_answers: str = "",
) -> dict[str, Any]:
    organism = case.get("organism") if isinstance(case.get("organism"), dict) else {}
    state = case.get("current_model_state") if isinstance(case.get("current_model_state"), dict) else {}
    local = case.get("local_evidence") if isinstance(case.get("local_evidence"), dict) else {}
    mngs = local.get("mngs_evidence") if isinstance(local.get("mngs_evidence"), dict) else {}
    hospital = local.get("same_organism_hospital_evidence")
    negatives = local.get("negative_evidence")
    return {
        "case_id": case.get("case_id"),
        "dataset": dataset_key(case),
        "dataset_label": case.get("dataset_label"),
        "patient_id": case.get("patient_id"),
        "organism_name": organism.get("display_name"),
        "canonical_key": organism.get("canonical_key"),
        "pathogen_category": organism.get("pathogen_category"),
        "source_tier": state.get("source_tier"),
        "label_type": label_type,
        "expected_rag_behavior": "retain_or_show_context" if label_type == "positive_control" else "reject_or_do_not_show",
        "benchmark_answer_private": benchmark_answer,
        "patient_answers_private": patient_answers or benchmark_answer,
        "match_type_private": match_type,
        "selection_reason_private": selection_reason,
        "reads": mngs.get("reads"),
        "reads_tier": mngs.get("reads_tier"),
        "dominance_tier": mngs.get("dominance_tier"),
        "same_organism_hospital_evidence_count": len(hospital) if isinstance(hospital, list) else 0,
        "negative_evidence_count": len(negatives) if isinstance(negatives, list) else 0,
    }


def parse_target_counts(values: Sequence[str] | None) -> dict[str, int]:
    counts = dict(DEFAULT_TARGET_NEGATIVE_COUNTS)
    for value in values or []:
        if "=" not in value:
            raise ValueError(f"Invalid --negative-count value: {value!r}; expected Category=Count")
        key, raw_count = value.split("=", 1)
        counts[key.strip()] = int(raw_count)
    return counts


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("production_queue", type=Path)
    parser.add_argument("--positive-cases", type=Path, required=True)
    parser.add_argument("--kh-answer-csv", type=Path, required=True)
    parser.add_argument("--two-hospitals-answer-csv", type=Path, required=True)
    parser.add_argument("--filtered-answer-csv", type=Path)
    parser.add_argument("--schemas-dir", type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--negative-count", action="append")
    parser.add_argument("--positive-dataset-default", default="KH")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    cases = read_jsonl(args.production_queue)
    answers = answer_lookup(args)
    selected: dict[str, dict[str, Any]] = {}
    labels: list[dict[str, Any]] = []
    missing_positive: list[dict[str, Any]] = []

    positive_payload = read_json(args.positive_cases)
    for expected in positive_payload.get("cases") or []:
        if not isinstance(expected, dict):
            continue
        case = find_positive_case(cases, expected, args.positive_dataset_default)
        if case is None:
            missing_positive.append(expected)
            continue
        selected[str(case["case_id"])] = case
        labels.append(
            label_row(
                case,
                label_type="positive_control",
                benchmark_answer=str(expected.get("answer_organism") or ""),
                match_type="expected_positive",
                selection_reason="previous picked-missed answer-hit case expected to remain in RAG production queue",
            )
        )

    target_counts = parse_target_counts(args.negative_count)
    candidates_by_category: dict[str, list[tuple[dict[str, Any], str, str]]] = defaultdict(list)
    skipped_without_answer = 0
    for case in cases:
        case_id = str(case.get("case_id") or "")
        if case_id in selected:
            continue
        dataset = dataset_key(case)
        patient_id = str(case.get("patient_id") or "")
        patient_answers = answers.get((dataset, patient_id), [])
        if not patient_answers:
            skipped_without_answer += 1
            continue
        matched, _answer, _match_type = case_matches_answers(case, patient_answers)
        if matched:
            continue
        patient_answers_text = "; ".join(answer.raw for answer in patient_answers)
        organism = case.get("organism") if isinstance(case.get("organism"), dict) else {}
        category = str(organism.get("pathogen_category") or "Other/unclassified")
        candidates_by_category[category].append((case, patient_answers_text, "no_match"))

    negative_labels: list[dict[str, Any]] = []
    for category, target_count in target_counts.items():
        rows = sorted(candidates_by_category.get(category, []), key=lambda row: negative_sort_key(row[0]))
        for case, answer, match_type in rows[: max(target_count, 0)]:
            selected[str(case["case_id"])] = case
            negative_labels.append(
                label_row(
                    case,
                    label_type="negative_control",
                    benchmark_answer=answer,
                    match_type=match_type,
                    patient_answers=answer,
                    selection_reason=f"non-answer high/context candidate selected from category: {category}",
                )
            )
    labels.extend(negative_labels)

    selected_cases = sorted(
        selected.values(),
        key=lambda case: (
            0 if any(row["case_id"] == case.get("case_id") and row["label_type"] == "positive_control" for row in labels) else 1,
            int(case.get("patient_id")) if str(case.get("patient_id")).isdigit() else 10**9,
            str(case.get("case_id") or ""),
        ),
    )
    selected_ids = {str(case.get("case_id") or "") for case in selected_cases}
    selected_labels = [row for row in labels if str(row.get("case_id") or "") in selected_ids]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_jsonl(args.output_dir / "production_queue.jsonl", selected_cases)
    write_jsonl(args.output_dir / "knowledge_card_requests.jsonl", rag_inputs.build_knowledge_requests(selected_cases))
    if args.schemas_dir and args.schemas_dir.exists():
        shutil.copytree(args.schemas_dir, args.output_dir / "schemas", dirs_exist_ok=True)
    write_json(
        args.output_dir / "calibration_labels.private.json",
        {
            "schema_version": "rag_calibration_labels_v1.0",
            "description": "Private labels. Do not provide this file to retrieval or adjudication models.",
            "cases": selected_labels,
        },
    )
    write_csv(
        args.output_dir / "calibration_labels.private.csv",
        selected_labels,
        (
            "case_id",
            "dataset",
            "patient_id",
            "organism_name",
            "pathogen_category",
            "source_tier",
            "label_type",
            "expected_rag_behavior",
            "benchmark_answer_private",
            "patient_answers_private",
            "match_type_private",
            "selection_reason_private",
            "reads",
            "reads_tier",
            "dominance_tier",
            "same_organism_hospital_evidence_count",
            "negative_evidence_count",
        ),
    )
    write_csv(
        args.output_dir / "production_queue_summary.csv",
        [
            {
                "case_id": case.get("case_id"),
                "dataset_label": case.get("dataset_label"),
                "patient_id": case.get("patient_id"),
                "organism_name": (case.get("organism") or {}).get("display_name"),
                "source_tier": (case.get("current_model_state") or {}).get("source_tier"),
                "pathogen_category": (case.get("organism") or {}).get("pathogen_category"),
                "knowledge_card_key": (case.get("rag_task") or {}).get("knowledge_card_key"),
            }
            for case in selected_cases
        ],
        (
            "case_id",
            "dataset_label",
            "patient_id",
            "organism_name",
            "source_tier",
            "pathogen_category",
            "knowledge_card_key",
        ),
    )
    category_counts = Counter(row["pathogen_category"] for row in selected_labels)
    label_counts = Counter(row["label_type"] for row in selected_labels)
    manifest = {
        "schema_version": "rag_calibration_manifest_v1.0",
        "label": args.label,
        "source_production_queue": str(args.production_queue),
        "production_case_count": len(selected_cases),
        "knowledge_request_count": len(rag_inputs.build_knowledge_requests(selected_cases)),
        "label_counts": dict(label_counts),
        "category_counts": dict(category_counts),
        "negative_candidate_counts_by_category": {
            category: len(rows) for category, rows in sorted(candidates_by_category.items())
        },
        "target_negative_counts": target_counts,
        "missing_positive_cases": missing_positive,
        "skipped_cases_without_private_answers": skipped_without_answer,
        "production_queue_answer_blind": True,
    }
    write_json(args.output_dir / "calibration_manifest.json", manifest)
    (args.output_dir / "README_繁中.md").write_text(
        (
            "# RAG v2 calibration set\n\n"
            "這是 RAG adjudication 校準用的小型 pilot set。\n\n"
            "- `production_queue.jsonl`：給 RAG 使用，answer-blind。\n"
            "- `knowledge_card_requests.jsonl`：PubMed/retrieval 需要查的去重菌種卡。\n"
            "- `calibration_labels.private.*`：私有評估標籤，不可給 RAG 或 adjudication model。\n\n"
            "目標是同時測 positive retention 與 negative rejection，而不是只測能不能保留答案 hit。\n"
        ),
        encoding="utf-8",
    )
    clean_readme_text = (
        "# RAG v2 calibration set\n\n"
        "這是 RAG adjudication 的正負案例校準集，用來測試 RAG 是否能保留重要候選並排除低價值候選。\n\n"
        "- `production_queue.jsonl`：給 RAG / adjudication model 使用，answer-blind，不含標準答案。\n"
        "- `knowledge_card_requests.jsonl`：給 PubMed / guideline retrieval 使用，用來抓文獻知識卡。\n"
        "- `calibration_labels.private.*`：私有答案與正負案例標籤，只能用於最後評估，不可以交給 RAG 或 adjudication model。\n\n"
        "評估重點是 positive controls 是否被保留，以及 negative controls 是否被拒絕或降到不給醫師看的層級。\n"
    )
    (args.output_dir / "README_繁中.md").write_text(clean_readme_text, encoding="utf-8")
    final_readme_text = (
        "# RAG v2 calibration set\n\n"
        "\u9019\u662f RAG adjudication \u7684\u6b63\u8ca0\u6848\u4f8b\u6821\u6e96\u96c6\uff0c"
        "\u7528\u4f86\u6e2c\u8a66 RAG \u662f\u5426\u80fd\u4fdd\u7559\u91cd\u8981\u5019\u9078\u4e26\u6392\u9664\u4f4e\u50f9\u503c\u5019\u9078\u3002\n\n"
        "- `production_queue.jsonl`\uff1a\u7d66 RAG / adjudication model \u4f7f\u7528\uff0canswer-blind\uff0c\u4e0d\u542b\u6a19\u6e96\u7b54\u6848\u3002\n"
        "- `knowledge_card_requests.jsonl`\uff1a\u7d66 PubMed / guideline retrieval \u4f7f\u7528\uff0c\u7528\u4f86\u6293\u6587\u737b\u77e5\u8b58\u5361\u3002\n"
        "- `calibration_labels.private.*`\uff1a\u79c1\u6709\u7b54\u6848\u8207\u6b63\u8ca0\u6848\u4f8b\u6a19\u7c64\uff0c"
        "\u53ea\u80fd\u7528\u65bc\u6700\u5f8c\u8a55\u4f30\uff0c\u4e0d\u53ef\u4ee5\u4ea4\u7d66 RAG \u6216 adjudication model\u3002\n\n"
        "\u8a55\u4f30\u91cd\u9ede\u662f positive controls \u662f\u5426\u88ab\u4fdd\u7559\uff0c"
        "\u4ee5\u53ca negative controls \u662f\u5426\u88ab\u62d2\u7d55\u6216\u964d\u5230\u4e0d\u7d66\u91ab\u5e2b\u770b\u7684\u5c64\u7d1a\u3002\n"
    )
    (args.output_dir / "README_\u7e41\u4e2d.md").write_text(final_readme_text, encoding="utf-8")
    print(f"selected_cases={len(selected_cases)}")
    print("label_counts=" + json.dumps(dict(label_counts), ensure_ascii=False, sort_keys=True))
    print("category_counts=" + json.dumps(dict(category_counts), ensure_ascii=False, sort_keys=True))
    print(f"missing_positive={len(missing_positive)}")
    print(f"skipped_without_answers={skipped_without_answer}")
    return 1 if missing_positive else 0


if __name__ == "__main__":
    raise SystemExit(main())

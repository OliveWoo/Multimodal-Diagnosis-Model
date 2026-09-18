from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Sequence


DEFAULT_SOURCE_ROOT = Path("0511跑") / "其他兩院資料_spilt" / "rich"
DEFAULT_COMBINED_OUTPUT = Path("outputs") / "mngs_candidate_microbes_rich_ranked.json"
DEFAULT_REPORT = Path("outputs") / "mngs_candidate_microbes_rich_ranked_report.json"

DUPLICATE_SUFFIX_RE = re.compile(r"\s+\(\d+\)(?=\.json$|$)")
PATIENT_ID_RE = re.compile(r"NGS_patient_(?P<id>\d+)(?:_json|_all_RK_NTC_microbes|_mNGS_grouped|$)")

GROUP_TO_SOURCE_CATEGORY = {
    "bacterial": "1.Bac",
    "fungal": "2.Fungi",
    "viral": "3.Virus",
    "others": "4.Others",
}

CODE_SCORE = {
    1: 1.0,
    2: 0.85,
    3: 0.65,
    4: 0.35,
    5: 0.2,
}


def clean_duplicate_suffix(name: str) -> str:
    previous = None
    current = name
    while current != previous:
        previous = current
        current = DUPLICATE_SUFFIX_RE.sub("", current)
    return current


def extract_patient_id(path: Path) -> str | None:
    match = PATIENT_ID_RE.search(clean_duplicate_suffix(path.name))
    if not match:
        return None
    return match.group("id")


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def normalize_code(raw: Any) -> str:
    text = str(raw or "").strip()
    if not text:
        return ""
    text = text.replace("Code_RK_NTC=RK08", "Code_RK_NTC=RK8")
    text = text.replace("Code_RK_NTC=RK00", "Code_RK_NTC=RK0")
    return text


def rank_from_codes(raw_codes: Any) -> tuple[int | None, str]:
    codes = {normalize_code(code) for code in (raw_codes or [])}
    codes.discard("")
    has_ntc00 = "Code_NTC=00" in codes
    has_rk0 = "Code_RK_NTC=RK0" in codes
    has_rk8 = "Code_RK_NTC=RK8" in codes

    if has_ntc00 and has_rk0:
        return 1, "Code_NTC=00 + Code_RK_NTC=RK00"
    if has_ntc00 and has_rk8:
        return 2, "Code_RK_NTC=RK8 + Code_NTC=00"
    if has_rk0:
        return 3, "Code_RK_NTC=RK00"
    if has_rk8:
        return 4, "Code_RK_NTC=RK8"
    if has_ntc00:
        return 5, "Code_NTC=00"
    return None, "Unknown"


def to_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def numeric_for_json(value: float) -> int | float:
    if float(value).is_integer():
        return int(value)
    return value


def possibility_level(rank_priority: int | None, reads: float) -> str:
    if reads <= 0 or rank_priority is None:
        return "low"
    if rank_priority in {1, 2}:
        return "high"
    if rank_priority == 3:
        return "medium"
    return "low"


def build_ranked_payload(grouped_path: Path, patient_id: str) -> dict[str, Any]:
    grouped = load_json(grouped_path)
    if not isinstance(grouped, list):
        raise ValueError(f"Grouped mNGS JSON must be a list: {grouped_path}")

    records: list[dict[str, Any]] = []
    for record in grouped:
        if not isinstance(record, dict):
            continue
        pathogens = record.get("pathogens")
        if not isinstance(pathogens, dict):
            continue

        candidates: list[dict[str, Any]] = []
        for group_name, source_category in GROUP_TO_SOURCE_CATEGORY.items():
            values = pathogens.get(group_name)
            if not isinstance(values, list):
                continue
            for item in values:
                if not isinstance(item, dict):
                    continue
                name = str(item.get("name") or "").strip()
                if not name:
                    continue
                reads = to_float(item.get("reads"))
                rank_priority, rank_rule = rank_from_codes(item.get("rk_ntc"))
                selection_reasons = [
                    normalized
                    for normalized in (normalize_code(code) for code in (item.get("rk_ntc") or []))
                    if normalized
                ]
                candidates.append(
                    {
                        "organism_name": name,
                        "source_category": source_category,
                        "reads": reads,
                        "sec_hit": reads,
                        "selection_reasons": selection_reasons,
                        "ranking": {
                            "rank_priority": rank_priority,
                            "rank_rule": rank_rule,
                            "code_score": CODE_SCORE.get(rank_priority or 0, 0.0),
                            "reads_percentile": 0.0,
                            "final_score": 0.0,
                            "possibility_level": possibility_level(rank_priority, reads),
                        },
                    }
                )

        if not candidates:
            continue

        sorted_reads = sorted(to_float(candidate.get("reads")) for candidate in candidates)
        total = len(sorted_reads)
        for candidate in candidates:
            reads = to_float(candidate.get("reads"))
            percentile = sum(1 for value in sorted_reads if value <= reads) / total
            ranking = candidate["ranking"]
            ranking["reads_percentile"] = round(percentile, 4)
            ranking["final_score"] = round(
                (to_float(ranking.get("code_score")) * 0.7) + (ranking["reads_percentile"] * 0.3),
                4,
            )
            candidate["reads"] = numeric_for_json(reads)
            candidate["sec_hit"] = numeric_for_json(reads)

        records.append(
            {
                "specimen_code": record.get("specimen_code"),
                "collected_time": record.get("collected_time"),
                "specimen_site": record.get("specimen_site"),
                "seq_id": record.get("seq_id") or record.get("specimen_code"),
                "status": "generated_from_mNGS_grouped",
                "source_grouped_path": str(grouped_path),
                "candidates": sorted(
                    candidates,
                    key=lambda item: (
                        item["ranking"]["rank_priority"] or 99,
                        -to_float(item["ranking"].get("final_score")),
                        -to_float(item.get("reads")),
                        str(item.get("organism_name")),
                    ),
                ),
            }
        )

    return {
        "patient_id": patient_id,
        "source_grouped_path": str(grouped_path),
        "records": records,
        "ranking_metadata": {
            "source": "generated_from_mNGS_grouped_folder",
            "rank_priority_order": [
                "1: Code_NTC=00 + Code_RK_NTC=RK00",
                "2: Code_RK_NTC=RK8 + Code_NTC=00",
                "3: Code_RK_NTC=RK00",
                "4: Code_RK_NTC=RK8",
                "5: Code_NTC=00",
            ],
            "score_formula": "final_score = code_score * 0.7 + reads_percentile * 0.3",
            "notes": [
                "Code_RK_NTC=RK0 is normalized to Code_RK_NTC=RK00.",
                "reads_percentile is computed within each specimen record.",
                "Candidates are sorted by rank_priority, final_score, reads, organism_name.",
            ],
        },
    }


def find_source_file(patient_dir: Path, patient_id: str, source_kind: str) -> Path | None:
    if source_kind == "all_rk_ntc":
        expected_clean = f"NGS_patient_{patient_id}_all_RK_NTC_microbes.json"
        pattern = f"NGS_patient_{patient_id}_all_RK_NTC_microbes*.json"
    elif source_kind == "mngs_grouped":
        expected_clean = f"NGS_patient_{patient_id}_mNGS_grouped.json"
        pattern = f"NGS_patient_{patient_id}_mNGS_grouped*.json"
    else:
        raise ValueError(f"Unsupported source_kind: {source_kind}")

    candidates = [
        path
        for path in patient_dir.glob(pattern)
        if clean_duplicate_suffix(path.name) == expected_clean
    ]
    if not candidates:
        return None
    return sorted(candidates, key=lambda path: (len(path.name), path.name))[0]


def source_file_pattern(source_kind: str) -> str:
    if source_kind == "all_rk_ntc":
        return "NGS_patient_*_all_RK_NTC_microbes*.json"
    if source_kind == "mngs_grouped":
        return "NGS_patient_*_mNGS_grouped*.json"
    raise ValueError(f"Unsupported source_kind: {source_kind}")


def output_patient_dir(source_root: Path, patient_id: str, patient_root: Path | None) -> Path:
    if patient_root is not None:
        return patient_root / f"NGS_patient_{patient_id}_json"
    return source_root / f"NGS_patient_{patient_id}_json"


def build_folder(
    source_root: Path,
    combined_output: Path,
    report_output: Path,
    *,
    source_kind: str = "mngs_grouped",
    patient_root: Path | None = None,
) -> dict[str, Any]:
    if not source_root.exists():
        raise FileNotFoundError(f"Source root does not exist: {source_root}")

    patients: dict[str, dict[str, Any]] = {}
    report_rows: list[dict[str, Any]] = []
    skipped: list[str] = []

    source_rows: list[tuple[str, Path, Path]] = []

    for patient_dir in sorted(path for path in source_root.iterdir() if path.is_dir()):
        patient_id = extract_patient_id(patient_dir)
        if patient_id is None:
            skipped.append(str(patient_dir))
            continue
        grouped_path = find_source_file(patient_dir, patient_id, source_kind)
        if grouped_path is None:
            skipped.append(f"{patient_dir}: {source_kind} file not found")
            continue
        source_rows.append((patient_id, grouped_path, patient_dir))

    for source_file in sorted(source_root.glob(source_file_pattern(source_kind))):
        patient_id = extract_patient_id(source_file)
        if patient_id is None:
            skipped.append(str(source_file))
            continue
        destination_dir = output_patient_dir(source_root, patient_id, patient_root)
        source_rows.append((patient_id, source_file, destination_dir))

    seen_outputs: set[Path] = set()
    for patient_id, grouped_path, patient_dir in source_rows:
        payload = build_ranked_payload(grouped_path, patient_id)
        output_path = patient_dir / f"NGS_patient_{patient_id}_mNGS_ranked_candidates.json"
        if output_path in seen_outputs:
            skipped.append(f"{grouped_path}: duplicate output {output_path}")
            continue
        seen_outputs.add(output_path)
        write_json(output_path, payload)
        patients[patient_id] = payload
        candidate_count = sum(len(record.get("candidates", [])) for record in payload.get("records", []))
        report_rows.append(
            {
                "patient_id": patient_id,
                "source_grouped_path": str(grouped_path),
                "ranked_candidates_path": str(output_path),
                "record_count": len(payload.get("records", [])),
                "candidate_count": candidate_count,
            }
        )

    combined = {
        "patients": patients,
        "ranking_metadata": {
            "source_root": str(source_root),
        "source": f"generated_from_per_patient_{source_kind}_json",
        "source_kind": source_kind,
            "patient_count": len(patients),
        },
    }
    write_json(combined_output, combined)

    report = {
        "source_root": str(source_root),
        "source_kind": source_kind,
        "combined_output": str(combined_output),
        "patient_root": str(patient_root) if patient_root else None,
        "patient_count": len(patients),
        "skipped_count": len(skipped),
        "patients": report_rows,
        "skipped": skipped,
    }
    write_json(report_output, report)
    return report


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build per-patient ranked mNGS candidates from *_mNGS_grouped*.json files."
    )
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument(
        "--source-kind",
        choices=["mngs_grouped", "all_rk_ntc"],
        default="mngs_grouped",
        help="Which per-patient source JSON to convert into ranked candidates.",
    )
    parser.add_argument("--combined-output", type=Path, default=DEFAULT_COMBINED_OUTPUT)
    parser.add_argument("--report-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--patient-root",
        type=Path,
        default=None,
        help=(
            "Optional output root containing NGS_patient_<id>_json folders. "
            "Use this when --source-root is a flat folder of *_all_RK_NTC_microbes.json files."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report = build_folder(
        args.source_root,
        args.combined_output,
        args.report_output,
        source_kind=args.source_kind,
        patient_root=args.patient_root,
    )
    print(f"patient_count={report['patient_count']}")
    print(f"skipped_count={report['skipped_count']}")
    print(f"combined_output={report['combined_output']}")
    print(f"report_output={args.report_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

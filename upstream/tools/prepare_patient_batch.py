from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path
from typing import Any, Sequence


DEFAULT_SOURCE_ROOT = Path("0511跑") / "其他兩院資料_spilt" / "rich"
DEFAULT_OUTPUT_ROOT = Path("outputs") / "patient_info_rich_normalized"
DUPLICATE_SUFFIX_RE = re.compile(r"\s+\(\d+\)(?=\.json$|$)")
PATIENT_DIR_RE = re.compile(r"NGS_patient_(?P<id>\d+)_json")

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


def clean_name(name: str) -> str:
    previous = None
    current = name
    while current != previous:
        previous = current
        current = DUPLICATE_SUFFIX_RE.sub("", current)
    return current


def extract_patient_id(path: Path) -> str | None:
    match = PATIENT_DIR_RE.search(clean_name(path.name))
    if not match:
        return None
    return match.group("id")


def copy_tree_normalized(source: Path, destination: Path) -> dict[str, Any]:
    copied_files = 0
    collisions: list[str] = []

    for item in sorted(source.rglob("*")):
        rel = item.relative_to(source)
        cleaned_parts = [clean_name(part) for part in rel.parts]
        target = destination.joinpath(*cleaned_parts)

        if item.is_dir():
            target.mkdir(parents=True, exist_ok=True)
            continue

        if target.exists():
            collisions.append(str(target))
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(item, target)
        copied_files += 1

    return {
        "source": str(source),
        "destination": str(destination),
        "copied_files": copied_files,
        "collisions": collisions,
    }


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


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


def possibility_level(rank_priority: int | None, reads: float) -> str:
    if reads <= 0 or rank_priority is None:
        return "low"
    if rank_priority in {1, 2}:
        return "high"
    if rank_priority == 3:
        return "medium"
    return "low"


def build_ranked_candidates(patient_dir: Path, patient_id: str) -> dict[str, Any] | None:
    base_name = f"NGS_patient_{patient_id}"
    grouped_path = patient_dir / f"{base_name}_mNGS_grouped.json"
    if not grouped_path.exists():
        return None

    grouped = load_json(grouped_path)
    if not isinstance(grouped, list):
        return None

    records: list[dict[str, Any]] = []
    for record in grouped:
        if not isinstance(record, dict):
            continue
        flat: list[dict[str, Any]] = []
        pathogens = record.get("pathogens")
        if not isinstance(pathogens, dict):
            continue

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
                flat.append(
                    {
                        "organism_name": name,
                        "source_category": source_category,
                        "reads": reads,
                        "sec_hit": reads,
                        "selection_reasons": [
                            normalize_code(code)
                            for code in (item.get("rk_ntc") or [])
                            if normalize_code(code)
                        ],
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

        if not flat:
            continue

        sorted_reads = sorted(candidate["reads"] for candidate in flat)
        total = len(sorted_reads)
        for candidate in flat:
            reads = candidate["reads"]
            percentile = sum(1 for value in sorted_reads if value <= reads) / total
            ranking = candidate["ranking"]
            ranking["reads_percentile"] = round(percentile, 4)
            ranking["final_score"] = round(
                (ranking["code_score"] * 0.7) + (ranking["reads_percentile"] * 0.3),
                4,
            )
            if reads.is_integer():
                candidate["reads"] = int(reads)
                candidate["sec_hit"] = int(reads)

        records.append(
            {
                "specimen_code": record.get("specimen_code"),
                "collected_time": record.get("collected_time"),
                "specimen_site": record.get("specimen_site"),
                "seq_id": record.get("seq_id") or record.get("specimen_code"),
                "status": "generated_from_mNGS_grouped",
                "candidates": sorted(
                    flat,
                    key=lambda item: (
                        item["ranking"]["rank_priority"] or 99,
                        -to_float(item.get("reads")),
                        str(item.get("organism_name")),
                    ),
                ),
            }
        )

    return {
        "patient_id": patient_id,
        "records": records,
        "ranking_metadata": {
            "source": "generated_from_patient_mNGS_grouped",
            "rank_priority_order": [
                "1: Code_NTC=00 + Code_RK_NTC=RK00",
                "2: Code_RK_NTC=RK8 + Code_NTC=00",
                "3: Code_RK_NTC=RK00",
                "4: Code_RK_NTC=RK8",
                "5: Code_NTC=00",
            ],
            "notes": [
                "Code_RK_NTC=RK0 in rich source data is normalized to Code_RK_NTC=RK00.",
                "reads_percentile is computed within each specimen record from mNGS_grouped reads.",
            ],
        },
    }


def prepare_batch(source_root: Path, output_root: Path, *, overwrite: bool) -> dict[str, Any]:
    if not source_root.exists():
        raise FileNotFoundError(f"Source root does not exist: {source_root}")
    if output_root.exists() and overwrite:
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    patient_reports: list[dict[str, Any]] = []
    skipped: list[str] = []
    for source_patient in sorted(path for path in source_root.iterdir() if path.is_dir()):
        patient_id = extract_patient_id(source_patient)
        if patient_id is None:
            skipped.append(str(source_patient))
            continue
        destination = output_root / f"NGS_patient_{patient_id}_json"
        if destination.exists() and not overwrite:
            skipped.append(f"{source_patient} -> {destination} exists")
            continue
        if destination.exists():
            shutil.rmtree(destination)
        copy_report = copy_tree_normalized(source_patient, destination)
        ranked_payload = build_ranked_candidates(destination, patient_id)
        ranked_path = None
        ranked_count = 0
        if ranked_payload is not None:
            ranked_path = destination / f"NGS_patient_{patient_id}_mNGS_ranked_candidates.json"
            write_json(ranked_path, ranked_payload)
            ranked_count = sum(len(record.get("candidates", [])) for record in ranked_payload["records"])
        patient_reports.append(
            {
                "patient_id": patient_id,
                "source": str(source_patient),
                "destination": str(destination),
                "copied_files": copy_report["copied_files"],
                "collisions": copy_report["collisions"],
                "ranked_candidates_path": str(ranked_path) if ranked_path else None,
                "ranked_candidate_count": ranked_count,
            }
        )

    report = {
        "source_root": str(source_root),
        "output_root": str(output_root),
        "patient_count": len(patient_reports),
        "skipped_count": len(skipped),
        "patients": patient_reports,
        "skipped": skipped,
    }
    report_path = output_root.parent / f"{output_root.name}_prepare_report.json"
    write_json(report_path, report)
    report["report_path"] = str(report_path)
    return report


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize external patient folders and generate local ranked mNGS candidates."
    )
    parser.add_argument("source_root", nargs="?", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report = prepare_batch(args.source_root, args.output_root, overwrite=args.overwrite)
    print(f"patient_count={report['patient_count']}")
    print(f"skipped_count={report['skipped_count']}")
    print(f"output_root={report['output_root']}")
    print(f"report={report['report_path']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

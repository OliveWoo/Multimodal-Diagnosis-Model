"""Build an event-aligned culture overlay from the three-hospital workbook.

The source workbook stores up to six culture events in repeated column blocks.
This adapter keeps every source event in an audit file, but only places culture
observations collected within the fixed index-event window in the decision
input.  Resistance labels such as CRKP are normalized to organism identity
(``Klebsiella pneumoniae``) plus a separate resistance phenotype.

Benchmark/clinical-diagnosis columns are intentionally never read.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from datetime import date, datetime
from pathlib import Path
from typing import Any, Iterable

import openpyxl
from openpyxl.utils import get_column_letter

from tools import build_deterministic_summary
from tools import pathogen_normalization
from tools.clinical_event_window import (
    DEFAULT_CALENDAR_DAYS_EACH_SIDE,
    PRIMARY_WINDOW_MODE,
    signed_calendar_day_delta,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_ROOT = (
    ROOT / "outputs/runs/2026-10-07_two_hospital_new_mngs_11_v2_6_input_v2"
)
DEFAULT_WORKBOOK = ROOT / "三院-final-mNGS merged file.xlsx"
DEFAULT_OUTPUT_ROOT = (
    ROOT
    / "outputs/runs/2026-10-08_two_hospital_A11_temporal_culture_overlay_v1"
)
DEFAULT_SUMMARY_SUFFIX = "release_input_new_mngs_11_temporal_culture_v1"
WORKSHEET_NAME = "微生物診斷資料"
EVENT_WINDOW_HOURS = 48.0
CULTURE_SLOT_START_COLUMNS = (4, 16, 28, 40, 52, 64)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_datetime(value: Any) -> datetime | None:
    if isinstance(value, datetime):
        return value.replace(tzinfo=None)
    if isinstance(value, date):
        return datetime.combine(value, datetime.min.time())
    text = str(value or "").strip()
    if not text:
        return None
    normalized = text.replace("/", "-")
    for pattern in (
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d",
    ):
        try:
            return datetime.strptime(normalized, pattern)
        except ValueError:
            continue
    raise ValueError(f"Unsupported date/time value: {value!r}")


def iso_date(value: Any) -> str:
    parsed = parse_datetime(value)
    return parsed.strftime("%Y-%m-%d") if parsed else ""


def split_organisms(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    return [
        item.strip()
        for item in re.split(r"[,;\n]+", text)
        if item.strip()
    ]


def alignment(index_time: datetime, observed_time: datetime) -> dict[str, Any]:
    delta_hours = (observed_time - index_time).total_seconds() / 3600.0
    calendar_delta = signed_calendar_day_delta(observed_time, index_time)
    aligned = abs(calendar_delta) <= DEFAULT_CALENDAR_DAYS_EACH_SIDE
    if aligned:
        status = "within_plus_minus_two_calendar_dates"
    elif calendar_delta < 0:
        status = "pre_index_outside_window"
    else:
        status = "post_index_outside_window"
    return {
        "event_aligned": aligned,
        "event_alignment_status": status,
        "signed_delta_hours": delta_hours,
        "signed_calendar_day_delta": calendar_delta,
        "event_window_mode": PRIMARY_WINDOW_MODE,
        "event_window_calendar_days_each_side": DEFAULT_CALENDAR_DAYS_EACH_SIDE,
        "event_window_hours": EVENT_WINDOW_HOURS,
        "event_window_hours_role": "audit_sensitivity_only",
    }


def find_source_row(
    sheet: Any, *, hospital_case_id: str, specimen_id: str
) -> tuple[int, tuple[Any, ...]]:
    matches: list[tuple[int, tuple[Any, ...]]] = []
    for row_number, row in enumerate(sheet.iter_rows(values_only=True), start=1):
        if row_number == 1:
            continue
        case_value = str(row[0] or "").strip() if len(row) >= 1 else ""
        specimen_value = str(row[1] or "").strip() if len(row) >= 2 else ""
        if case_value == hospital_case_id and specimen_value == specimen_id:
            matches.append((row_number, tuple(row)))
    if len(matches) != 1:
        raise ValueError(
            f"Expected one workbook row for {hospital_case_id}/{specimen_id}; "
            f"found {len(matches)}"
        )
    return matches[0]


def extract_culture_records(
    row_number: int,
    row: tuple[Any, ...],
    *,
    index_time: datetime,
    workbook_path: Path,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for slot_number, start_column in enumerate(CULTURE_SLOT_START_COLUMNS, start=1):
        start_index = start_column - 1
        if start_index + 10 >= len(row):
            continue
        specimen = str(row[start_index] or "").strip()
        organisms_cell = row[start_index + 2]
        collected_value = row[start_index + 9]
        reported_value = row[start_index + 10]
        resistance_flag = str(row[start_index + 11] or "").strip()
        organisms = split_organisms(organisms_cell)
        if not organisms:
            continue
        collected_time = parse_datetime(collected_value)
        if collected_time is None:
            raise ValueError(
                f"Culture slot {slot_number} row {row_number} has organisms but no collection date"
            )
        timing = alignment(index_time, collected_time)
        for source_label in organisms:
            canonical_key = pathogen_normalization.canonical_key(source_label)
            canonical_display = (
                pathogen_normalization.display_name(source_label) or source_label
            )
            resistance_note = pathogen_normalization.resistance_note(source_label)
            records.append(
                {
                    "test": "Microbiology culture",
                    "sample": specimen,
                    "collected_time": iso_date(collected_value),
                    "received_time": "",
                    "reported_time": iso_date(reported_value),
                    "organism": canonical_display,
                    "status": "Detected",
                    "antimicrobial_susceptibility_test": [],
                    "source_organism_label": source_label,
                    "canonical_organism_key": canonical_key,
                    "resistance_phenotype": (
                        {
                            "label": resistance_note,
                            "source_label": source_label,
                            "relationship": "resistance_phenotype",
                            "organism_identity": canonical_display,
                            "status": "reported",
                        }
                        if resistance_note
                        else None
                    ),
                    "source_resistance_flag": resistance_flag,
                    "source_provenance": {
                        "workbook": str(workbook_path.resolve()),
                        "worksheet": WORKSHEET_NAME,
                        "row": row_number,
                        "culture_slot": slot_number,
                        "organism_cell": f"{get_column_letter(start_column + 2)}{row_number}",
                        "collected_time_cell": f"{get_column_letter(start_column + 9)}{row_number}",
                        "reported_time_cell": f"{get_column_letter(start_column + 10)}{row_number}",
                    },
                    **timing,
                }
            )
    return records


def decision_culture_row(record: dict[str, Any]) -> dict[str, Any]:
    return {
        key: record[key]
        for key in (
            "test",
            "sample",
            "collected_time",
            "received_time",
            "reported_time",
            "organism",
            "status",
            "antimicrobial_susceptibility_test",
            "source_organism_label",
            "canonical_organism_key",
            "resistance_phenotype",
            "source_resistance_flag",
            "source_provenance",
            "event_aligned",
            "event_alignment_status",
            "signed_delta_hours",
            "event_window_hours",
        )
    }


def patient_documents(
    root: Path, patient_id: int
) -> tuple[Path, Path, dict[str, Any]]:
    patient_dir = root / "patients" / f"NGS_patient_{patient_id}_json"
    multi_path = (
        root
        / "multi_assay_source"
        / f"NGS_patient_{patient_id}_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
    )
    payload = read_json(multi_path)
    if int(payload.get("patient_id") or -1) != patient_id:
        raise ValueError(f"mNGS patient mismatch in {multi_path}")
    return patient_dir, multi_path, payload


def case_identity(payload: dict[str, Any]) -> tuple[str, str, datetime]:
    cases = payload.get("cases") or []
    if len(cases) != 1:
        raise ValueError(f"Expected exactly one mNGS case; found {len(cases)}")
    case = cases[0]
    hospital_case_id = str(payload.get("hospital_case_id") or "").strip()
    specimen_id = str(case.get("specimen_code") or case.get("case_review_id") or "").strip()
    index_time = parse_datetime(case.get("collected_time"))
    if not hospital_case_id or not specimen_id or index_time is None:
        raise ValueError("Missing hospital case, specimen, or index-event date")
    return hospital_case_id, specimen_id, index_time


def link_patient(
    workbook: Any,
    output_root: Path,
    workbook_path: Path,
    patient_id: int,
) -> dict[str, Any]:
    patient_dir, multi_path, multi_payload = patient_documents(output_root, patient_id)
    hospital_case_id, specimen_id, index_time = case_identity(multi_payload)
    sheet = workbook[WORKSHEET_NAME]
    row_number, row = find_source_row(
        sheet,
        hospital_case_id=hospital_case_id,
        specimen_id=specimen_id,
    )
    records = extract_culture_records(
        row_number,
        row,
        index_time=index_time,
        workbook_path=workbook_path,
    )
    aligned = [item for item in records if item["event_aligned"]]
    outside = [item for item in records if not item["event_aligned"]]
    culture_path = patient_dir / f"NGS_patient_{patient_id}_culture.json"
    write_json(culture_path, [decision_culture_row(item) for item in aligned])
    all_path = patient_dir / f"NGS_patient_{patient_id}_culture_all_events.json"
    outside_path = (
        patient_dir / f"NGS_patient_{patient_id}_culture_outside_event_window.json"
    )
    write_json(all_path, records)
    write_json(outside_path, outside)

    phenotype_links = [
        {
            "source_organism_label": item["source_organism_label"],
            "organism_identity": item["organism"],
            "canonical_organism_key": item["canonical_organism_key"],
            "resistance_phenotype": item["resistance_phenotype"],
            "collected_time": item["collected_time"],
            "event_alignment_status": item["event_alignment_status"],
            "signed_delta_hours": item["signed_delta_hours"],
            "decision_effect": (
                "eligible_as_event_aligned_culture_evidence"
                if item["event_aligned"]
                else "audit_only_not_direct_hospital_evidence"
            ),
            "source_provenance": item["source_provenance"],
        }
        for item in records
        if item.get("resistance_phenotype")
    ]
    linkage_path = (
        patient_dir / f"NGS_patient_{patient_id}_resistance_phenotype_linkage.json"
    )
    linkage = {
        "schema_version": "temporal_resistance_phenotype_linkage.v1",
        "patient_id": patient_id,
        "answer_blind": True,
        "index_event": {
            "hospital_case_id": hospital_case_id,
            "specimen_id": specimen_id,
            "collected_time": index_time.strftime("%Y-%m-%d"),
            "window_hours": EVENT_WINDOW_HOURS,
        },
        "links": phenotype_links,
        "constraint": (
            "Organism identity and resistance phenotype are linked, but only "
            "event-aligned microbiology may act as direct hospital evidence."
        ),
    }
    write_json(linkage_path, linkage)
    return {
        "patient_id": patient_id,
        "hospital_case_id": hospital_case_id,
        "specimen_id": specimen_id,
        "mngs_source": str(multi_path.resolve()),
        "index_time": index_time.strftime("%Y-%m-%d"),
        "workbook_row": row_number,
        "all_culture_observation_count": len(records),
        "event_aligned_observation_count": len(aligned),
        "outside_window_observation_count": len(outside),
        "resistance_phenotype_link_count": len(phenotype_links),
        "event_aligned_culture": str(culture_path.resolve()),
        "all_event_audit": str(all_path.resolve()),
        "outside_window_audit": str(outside_path.resolve()),
        "phenotype_linkage": str(linkage_path.resolve()),
    }


def attach_audit_to_summary(
    patient_dir: Path,
    patient_id: int,
    summary_suffix: str,
    patient_manifest: dict[str, Any],
) -> None:
    summary_path = (
        patient_dir
        / "summary_outputs"
        / f"NGS_patient_{patient_id}_{summary_suffix}.json"
    )
    summary = read_json(summary_path)
    outside = read_json(Path(patient_manifest["outside_window_audit"]))
    linkage = read_json(Path(patient_manifest["phenotype_linkage"]))
    summary["temporal_hospital_evidence_audit"] = {
        "event_window_hours": EVENT_WINDOW_HOURS,
        "outside_window_evidence_count": len(outside),
        "outside_window_evidence": outside,
        "resistance_phenotype_linkage": linkage,
        "decision_effect": "outside-window records are audit-visible and cannot provide direct hospital evidence",
    }
    write_json(summary_path, summary)


def run(
    source_root: Path,
    workbook_path: Path,
    output_root: Path,
    patient_ids: Iterable[int],
    summary_suffix: str,
) -> dict[str, Any]:
    if output_root.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {output_root}")
    if not workbook_path.is_file():
        raise FileNotFoundError(workbook_path)
    shutil.copytree(source_root, output_root)
    workbook = openpyxl.load_workbook(workbook_path, read_only=True, data_only=True)
    if WORKSHEET_NAME not in workbook.sheetnames:
        raise ValueError(f"Missing worksheet {WORKSHEET_NAME!r}")
    patient_manifests = [
        link_patient(workbook, output_root, workbook_path, int(patient_id))
        for patient_id in sorted(set(patient_ids))
    ]
    workbook.close()

    all_patient_dirs = sorted(
        (output_root / "patients").glob("NGS_patient_*_json"),
        key=lambda path: int(re.search(r"(\d+)", path.name).group(1)),
    )
    for patient_dir in all_patient_dirs:
        build_deterministic_summary.process_patient(
            patient_dir,
            output_suffix=summary_suffix,
            agent_suffix_tag=None,
            overwrite=True,
            skip_existing=False,
        )
    manifest_by_patient = {
        int(item["patient_id"]): item for item in patient_manifests
    }
    for patient_id, patient_manifest in manifest_by_patient.items():
        patient_dir = output_root / "patients" / f"NGS_patient_{patient_id}_json"
        attach_audit_to_summary(
            patient_dir,
            patient_id,
            summary_suffix,
            patient_manifest,
        )

    manifest = {
        "schema_version": "two_hospital_temporal_culture_overlay.v1",
        "answer_blind": True,
        "source_root": str(source_root.resolve()),
        "output_root": str(output_root.resolve()),
        "source_workbook": str(workbook_path.resolve()),
        "source_workbook_sha256": sha256_file(workbook_path),
        "worksheet": WORKSHEET_NAME,
        "clinical_diagnosis_columns_read": False,
        "event_window_hours": EVENT_WINDOW_HOURS,
        "summary_suffix": summary_suffix,
        "patient_count": len(patient_manifests),
        "patients": patient_manifests,
        "rules": [
            "Resistance phenotype labels are separated from canonical organism identity.",
            "Only collection times within plus/minus 48 hours enter decision culture input.",
            "All other culture events remain available in explicit audit files.",
            "Clinical diagnosis/answer columns are never read by this adapter.",
        ],
    }
    write_json(output_root / "temporal_culture_overlay_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--patients", nargs="+", type=int, default=[62])
    parser.add_argument("--summary-suffix", default=DEFAULT_SUMMARY_SUFFIX)
    args = parser.parse_args()
    print(
        json.dumps(
            run(
                args.source_root,
                args.workbook,
                args.output_root,
                args.patients,
                args.summary_suffix,
            ),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

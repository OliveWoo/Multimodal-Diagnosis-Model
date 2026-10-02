from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timedelta
from functools import lru_cache
from pathlib import Path
from typing import Any

from tools.evidence_preservation import cbc_data_status, extract_host_risk_medications, extract_imaging_diagnostic_mentions


DATE_FORMATS = (
    ("%Y-%m-%d %H:%M:%S", "second"),
    ("%Y-%m-%d %H:%M", "minute"),
    ("%Y/%m/%d %H:%M:%S", "second"),
    ("%Y/%m/%d %H:%M", "minute"),
    ("%Y-%m-%d", "day"),
    ("%Y/%m/%d", "day"),
    ("%Y-%m", "month"),
    ("%Y/%m", "month"),
)
EPISODE_WINDOW = timedelta(days=2)


def time_interval(value: Any) -> tuple[datetime, datetime, str] | None:
    text = str(value or "").strip()
    if text.lower() in {"", "unknown", "none", "-"}:
        return None
    for fmt, precision in DATE_FORMATS:
        try:
            start = datetime.strptime(text, fmt)
        except ValueError:
            continue
        if precision == "month":
            end = datetime(start.year + (start.month == 12), start.month % 12 + 1, 1)
        else:
            end = start + {"day": timedelta(days=1), "minute": timedelta(minutes=1), "second": timedelta(seconds=1)}[precision]
        return start, end, precision
    return None


def temporal_relation(value: Any, index_value: Any) -> str:
    event = time_interval(value)
    index = time_interval(index_value)
    if event is None or index is None:
        return "unknown"
    if event[1] <= index[0]:
        return "before_index"
    if event[0] >= index[1]:
        return "after_index"
    return "overlaps_index_precision_uncertain"


def episode_window_relation(value: Any, index_value: Any) -> str:
    observation = time_interval(value)
    index = time_interval(index_value)
    if observation is None or index is None:
        return "unknown"
    earliest_start = index[0] - EPISODE_WINDOW
    latest_end = index[1] + EPISODE_WINDOW
    if observation[1] <= earliest_start or observation[0] >= latest_end:
        return "outside_window"
    if observation[2] in {"minute", "second"} and index[2] in {"minute", "second"}:
        return "within_window" if earliest_start <= observation[0] <= index[0] + EPISODE_WINDOW else "outside_window"
    if observation[0] >= index[1] - EPISODE_WINDOW and observation[1] <= index[0] + EPISODE_WINDOW:
        return "within_window"
    return "boundary_uncertain"


def episode_evidence_status(role: str, observation_window: str, report_window: str) -> str:
    if role in {"host_history", "host_medication"}:
        if observation_window == "outside_window" and role == "host_medication":
            return "host_exposure_timing_review"
        return "host_context_timing_review"
    if observation_window == "within_window":
        return "episode_observation"
    if observation_window == "outside_window":
        return "outside_episode_window"
    if observation_window == "boundary_uncertain":
        return "episode_boundary_review"
    if report_window == "within_window":
        return "report_only_observation_time_unverified"
    return "observation_time_unknown_review"


def availability_at_collection(role: str, observed_at: Any, reported_at: Any, index_time: Any) -> str:
    observed = temporal_relation(observed_at, index_time)
    reported = temporal_relation(reported_at, index_time)
    if observed == "after_index" and reported == "before_index":
        return "time_conflict_review"
    if observed == "after_index" or reported == "after_index":
        return "post_index_not_available"
    if observed == "overlaps_index_precision_uncertain" or reported == "overlaps_index_precision_uncertain":
        return "timing_uncertain"
    if role == "host_history":
        return "preexisting_context_documentation_unverified"
    if role == "host_medication":
        return "preexisting_context_documentation_unverified" if observed == "before_index" else "timing_unknown"
    if reported == "before_index":
        return "confirmed_available_before_index"
    if observed == "before_index":
        return "report_time_missing_availability_unverified"
    return "timing_unknown"


@lru_cache(maxsize=128)
def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_source(patient_dir: Path, patient_id: str, suffix: str) -> tuple[Path | None, Any]:
    path = patient_dir / f"NGS_patient_{patient_id}_{suffix}.json"
    if not path.is_file():
        return None, None
    return path, json.loads(path.read_text(encoding="utf-8-sig"))


def source_pointer(path: Path, json_path: str, *, claim_path: str | None = None) -> dict[str, str]:
    source_sha = sha256_file(path)
    return {
        "file": str(path.resolve()),
        "sha256": source_sha,
        "json_path": json_path,
        "source_claim_group": f"{source_sha}:{claim_path or json_path}",
    }


def build_patient(patient_dir: Path, multi_assay_path: Path | None = None) -> dict[str, Any]:
    match = re.fullmatch(r"NGS_patient_(\d+)_json", patient_dir.name)
    if match is None:
        raise ValueError(f"Not a patient directory: {patient_dir}")
    patient_id = match.group(1)
    events = []
    seen: set[tuple[str, str, str]] = set()
    if multi_assay_path is not None:
        multi_assay = json.loads(multi_assay_path.read_text(encoding="utf-8-sig"))
        if str(int(multi_assay.get("patient_id"))) != str(int(patient_id)):
            raise ValueError(f"Multi-assay patient mismatch: {multi_assay_path}")
        cases = multi_assay.get("cases") or []
        if not isinstance(cases, list) or not cases:
            raise ValueError(f"No multi-assay cases for P{patient_id}")
        for i, case in enumerate(cases):
            key = (
                str(case.get("specimen_code") or case.get("case_review_id") or ""),
                str(case.get("collected_time") or ""),
                str(case.get("specimen_site") or ""),
            )
            if key in seen:
                continue
            seen.add(key)
            tests = case.get("tests") or []
            seq_ids = sorted({str(item.get("seq_id")) for item in tests if item.get("seq_id")})
            events.append({
                "event_id": f"P{patient_id}:mngs:{len(events) + 1}",
                "case_review_id": case.get("case_review_id"),
                "specimen_code": case.get("specimen_code") or case.get("case_review_id"),
                "specimen_site": case.get("specimen_site"),
                "collected_time": case.get("collected_time"),
                "seq_ids": seq_ids,
                "source_ref": source_pointer(multi_assay_path, f"cases[{i}]"),
            })
        schema_version = "index_event_clinical_timeline_shadow.v4"
        index_event_source = "multi_assay_case_inventory"
        multi_assay_status = "joined_by_patient_and_case_source_contract"
    else:
        ranked_path, ranked = read_source(patient_dir, patient_id, "mNGS_ranked_candidates")
        if ranked_path is None or not isinstance(ranked, dict):
            raise FileNotFoundError(f"Index mNGS ranked records missing for P{patient_id}")
        records = ranked.get("records") or []
        if not isinstance(records, list) or not records:
            raise ValueError(f"No index mNGS events for P{patient_id}")
        for i, record in enumerate(records):
            key = (str(record.get("specimen_code") or ""), str(record.get("collected_time") or ""), str(record.get("specimen_site") or ""))
            if key in seen:
                continue
            seen.add(key)
            events.append({
                "event_id": f"P{patient_id}:mngs:{len(events) + 1}",
                "case_review_id": record.get("specimen_code"),
                "specimen_code": record.get("specimen_code"),
                "specimen_site": record.get("specimen_site"),
                "collected_time": record.get("collected_time"),
                "seq_id": record.get("seq_id"),
                "source_ref": source_pointer(ranked_path, f"records[{i}]"),
            })
        schema_version = "index_event_clinical_timeline_shadow.v3"
        index_event_source = "current_ranked_mngs_records_only"
        multi_assay_status = "not_joined_pending_lab_field_definitions"
    sources = {name: read_source(patient_dir, patient_id, suffix) for name, suffix in (
        ("underlying", "underlying"),
        ("admission_diagnosis", "admission_diagnosis"),
        ("cbc", "CBC"),
        ("other_lab", "other_lab"),
        ("image", "image"),
    )}
    evidence: list[dict[str, Any]] = []

    def add(
        role: str, pointer: dict[str, str], *, label: str, assertion: str,
        observed_at: Any = None, reported_at: Any = None, time_basis: str = "unknown",
        source_text: Any = None, value: Any = None, unit: Any = None,
        diagnostic_mentions: list[dict[str, Any]] | None = None,
        details: dict[str, Any] | None = None, qc_flags: list[str] | None = None,
    ) -> None:
        observation_relation = {event["event_id"]: temporal_relation(observed_at, event["collected_time"]) for event in events}
        report_relation = {event["event_id"]: temporal_relation(reported_at, event["collected_time"]) for event in events}
        availability = {event["event_id"]: availability_at_collection(
            role, observed_at, reported_at, event["collected_time"]
        ) for event in events}
        observation_window = {event["event_id"]: episode_window_relation(
            observed_at, event["collected_time"]
        ) for event in events}
        report_window = {event["event_id"]: episode_window_relation(
            reported_at, event["collected_time"]
        ) for event in events}
        episode_status = {event_id: episode_evidence_status(
            role, relation, report_window[event_id]
        ) for event_id, relation in observation_window.items()}
        evidence.append({
            "evidence_id": f"P{patient_id}:E{len(evidence) + 1:05d}",
            "evidence_role": role,
            "label": label,
            "assertion": assertion,
            "microbiologic_confirmation": False,
            "source_ref": pointer,
            "source_text": source_text,
            "value": value,
            "unit": unit,
            "observed_at": observed_at,
            "reported_at": reported_at,
            "time_basis": time_basis,
            "observation_relation_to_index": observation_relation,
            "report_relation_to_index": report_relation,
            "available_at_index": availability,
            "observation_window_relation_to_index": observation_window,
            "report_window_relation_to_index": report_window,
            "episode_evidence_status": episode_status,
            "event_link_status": "patient_level_temporal_only_not_proven_same_infection_event",
            "relative_to_treatment": "unresolved_no_verified_treatment_timeline",
            "diagnostic_mentions": diagnostic_mentions or [],
            "details": details or {},
            "qc_flags": qc_flags or [],
            "derived_from": [],
        })

    underlying_path, underlying = sources["underlying"]
    if underlying_path is not None and isinstance(underlying, list):
        for i, row in enumerate(underlying):
            if not isinstance(row, dict):
                continue
            for j, disease in enumerate(row.get("underlying_diseases") or []):
                if not str(disease).strip():
                    continue
                add("host_history", source_pointer(underlying_path, f"[{i}].underlying_diseases[{j}]") ,
                    label=str(disease), assertion="recorded_history", source_text=str(disease),
                    qc_flags=["history_onset_time_unknown"])

    medication_inputs = {
        name: payload for name in ("underlying", "admission_diagnosis")
        if (payload := sources[name][1]) is not None
    }
    medication_profile = extract_host_risk_medications(
        medication_inputs,
        sample_time=events[0].get("collected_time") if events else None,
    )
    for medication in medication_profile["medications"]:
        source_path = str(medication["source_path"])
        source_name = source_path.split("[", 1)[0]
        file_path = sources[source_name][0]
        if file_path is None:
            continue
        # Multiple drugs in one note share a claim group, not independent support.
        add("host_medication", source_pointer(file_path, source_path),
            label=medication["medication_name"],
            assertion=(
                "recorded_pre_sample_exposure_not_administration_verified"
                if medication["role"] == "pre_sample_risk_medication"
                else "recorded_exposure_not_administration_verified"
            ),
            observed_at=None if medication["start_date"] == "Unknown" else medication["start_date"],
            time_basis="recorded_start_date_not_administration_time", source_text=medication["source_text"],
            details={
                "medication_class": medication["medication_class"],
                "record_role": medication["role"],
                "dose": medication["dose"],
                "end_date": medication["end_date"],
                "ongoing_at_source": medication["ongoing_at_source"],
                "temporal_status": medication["temporal_status"],
                "days_from_start_to_sample": medication[
                    "days_from_start_to_sample"
                ],
                "days_from_end_to_sample": medication[
                    "days_from_end_to_sample"
                ],
            }, qc_flags=(
                (["dose_unknown"] if medication["dose"] == "Unknown" else [])
                + (["medication_timing_unknown"] if medication["temporal_status"] == "timing_unknown" else [])
            ))

    for source_name in ("cbc", "other_lab"):
        path, rows = sources[source_name]
        if path is None or not isinstance(rows, list):
            continue
        for i, row in enumerate(rows):
            if not isinstance(row, dict):
                continue
            label = str(row.get("item_full") or row.get("item") or "Unknown")
            reported_at = row.get("reported_time")
            collected_at = row.get("collected_time")
            add("illness_lab", source_pointer(path, f"[{i}]"), label=label, assertion="recorded_result",
                observed_at=collected_at, reported_at=reported_at,
                time_basis="collected_time" if collected_at else "reported_time_only_observation_time_unknown",
                value=row.get("value"), unit=row.get("unit"),
                details={"module": source_name, "test": row.get("test")},
                qc_flags=[] if collected_at else ["collection_time_missing"])

    image_path, image_rows = sources["image"]
    if image_path is not None and isinstance(image_rows, list):
        mentions_by_path: dict[str, list[dict[str, Any]]] = {}
        for mention in extract_imaging_diagnostic_mentions(image_rows):
            mentions_by_path.setdefault(str(mention["source_path"]), []).append({
                "concept_name": mention["concept_name"],
                "concept_type": mention["concept_type"],
                "assertion": mention["assertion"],
            })
        for i, row in enumerate(image_rows):
            if not isinstance(row, dict):
                continue
            for j, finding in enumerate(row.get("findings") or []):
                if not isinstance(finding, str) or not finding.strip():
                    continue
                json_path = f"[{i}].findings[{j}]"
                add("imaging_context", source_pointer(image_path, json_path),
                    label=str(row.get("exam_type") or row.get("test") or "image"),
                    assertion="described_finding", observed_at=row.get("collected_time"),
                    reported_at=row.get("reported_time"), time_basis="examination_collected_time",
                    source_text=finding, diagnostic_mentions=mentions_by_path.get(f"image[{i}].findings[{j}]", []),
                    details={"exam_type": row.get("exam_type")})

    return {
        "schema_version": schema_version,
        "patient_id": patient_id,
        "answer_blind": True,
        "decision_effect": "none_shadow_only",
        "primary_analysis": "retrospective_infection_episode_window",
        "episode_window_anchor": "index_mngs_specimen_collection_time",
        "episode_window_hours_each_side": 48,
        "primary_view_is_at_collection_prediction": False,
        "index_event_source": index_event_source,
        "multi_assay_event_link_status": multi_assay_status,
        "index_events": events,
        "source_availability": {name: path is not None for name, (path, _) in sources.items()},
        "clinical_data_status": {"cbc": cbc_data_status(sources["cbc"][1])},
        "phenotype_link_status": "not_attached_patient_and_event_link_unverified",
        "clinical_evidence": evidence,
    }


def build_shadow(
    patient_root: Path,
    output_dir: Path,
    patient_ids: list[str] | None = None,
    multi_assay_root: Path | None = None,
) -> dict[str, Any]:
    directories = sorted(patient_root.glob("NGS_patient_*_json"))
    if patient_ids is not None:
        wanted = {str(int(item)) for item in patient_ids}
        directories = [path for path in directories if re.fullmatch(r"NGS_patient_(\d+)_json", path.name)
                       and path.name.split("_")[2] in wanted]
        if {path.name.split("_")[2] for path in directories} != wanted:
            raise FileNotFoundError("One or more requested patient directories are missing")
    if not directories:
        raise FileNotFoundError("No patient directories found")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Shadow output already exists: {output_dir}")
    packets = []
    for path in directories:
        patient_id = path.name.split("_")[2]
        multi_assay_path = None
        if multi_assay_root is not None:
            multi_assay_path = multi_assay_root / f"NGS_patient_{patient_id}_multi_assay_mngs_evidence.json"
            if not multi_assay_path.is_file():
                raise FileNotFoundError(multi_assay_path)
        packets.append(build_patient(path, multi_assay_path))
    output_dir.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    for packet in packets:
        counts["patients"] += 1
        counts["index_events"] += len(packet["index_events"])
        counts.update(row["evidence_role"] for row in packet["clinical_evidence"])
        counts.update(
            f"availability:{status}"
            for row in packet["clinical_evidence"]
            for status in row["available_at_index"].values()
        )
        counts.update(
            f"episode:{status}"
            for row in packet["clinical_evidence"]
            for status in row["episode_evidence_status"].values()
        )
        (output_dir / f"NGS_patient_{packet['patient_id']}_clinical_timeline_shadow.json").write_text(
            json.dumps(packet, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
    summary = {"schema_version": (
                   "index_event_clinical_timeline_shadow_summary.v4"
                   if multi_assay_root is not None else
                   "index_event_clinical_timeline_shadow_summary.v3"
               ),
               "source_root": str(patient_root.resolve()), "counts": dict(counts),
               "multi_assay_event_root": str(multi_assay_root.resolve()) if multi_assay_root else None,
               "primary_analysis": "retrospective_infection_episode_window",
               "episode_window_anchor": "index_mngs_specimen_collection_time",
               "episode_window_hours_each_side": 48,
               "primary_view_is_at_collection_prediction": False,
               "decision_effect": "none_shadow_only", "phenotype_attached": False,
               "treatment_response_inferred": False}
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a source-traceable clinical timeline without scoring")
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--patients", nargs="+", help="Optional patient numbers, e.g. 9 14")
    parser.add_argument("--multi-assay-root", type=Path,
                        help="Optional patient-level multi-assay inventory directory used as the event source.")
    args = parser.parse_args()
    print(json.dumps(build_shadow(
        args.patient_root,
        args.output_dir,
        args.patients,
        multi_assay_root=args.multi_assay_root,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build per-patient candidate-window microbe evidence bundles.

The candidate list is intentionally separated from the evidence scope:

* candidate_window_days controls which positive microbiology records can create
  candidate microbes.
* evidence_scope=all_time keeps all available records from the input patient
  directory as evidence for those candidates.

The script is written for roots containing structured NGS_patient_*_json
folders, such as data/patients.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable


PATIENT_DIR_RE = re.compile(r"NGS_patient_(\d+)_json$")
PATIENT_FILE_RE = re.compile(r"NGS_patient_(\d+)_(.+)\.json$", re.IGNORECASE)
DATE_RE = re.compile(
    r"\d{4}[-/]\d{1,2}[-/]\d{1,2}(?:[ T]\d{1,2}:\d{2}(?::\d{2})?)?"
)

DATE_FIELDS = ("collected_time", "reported_time", "received_time")
MICROBE_RECORD_TYPES = {
    "culture",
    "filmarray",
    "molecular_microbiology",
    "mngs_grouped",
    "mngs_name_reads",
    "all_RK_NTC_microbes",
}
CONTEXT_RECORD_TYPES = {
    "admission_diagnosis",
    "underlying",
    "image",
    "CBC",
    "other_lab",
    "gm_test",
}
RESISTANCE_MARKERS = {
    "ctx-m",
    "ctxm",
    "imp",
    "kpc",
    "ndm",
    "oxa-48-like",
    "oxa 48 like",
    "vim",
    "meca/c and mreJ".lower(),
    "meca c and mrej",
    "mecA/C and MREJ".lower(),
}
NON_SPECIFIC_MICROBE_PATTERNS = [
    "identification to follow",
    "preliminary report",
    "final report",
    "organism isolated",
    "no growth",
    "normal flora",
    "mixed flora",
    "commensal flora",
    "not detected",
    "negative",
    "g(+) bacteria",
    "g(+) bacilli",
    "g(-) bacteria",
    "g(-) bacilli",
]

ALIASES = {
    "human respirovirus 3": {"parainfluenza virus", "human parainfluenza virus 3"},
    "parainfluenza virus": {"human respirovirus 3", "human parainfluenza virus 3"},
    "human betaherpesvirus 5": {"cmv", "cytomegalovirus"},
    "cmv": {"human betaherpesvirus 5", "cytomegalovirus"},
    "human gammaherpesvirus 4": {"ebv", "epstein barr virus"},
    "ebv": {"human gammaherpesvirus 4", "epstein barr virus"},
    "human alphaherpesvirus 1": {"hsv-1", "herpes simplex virus 1"},
    "hsv-1": {"human alphaherpesvirus 1", "herpes simplex virus 1"},
    "human alphaherpesvirus 2": {"hsv-2", "herpes simplex virus 2"},
    "hsv-2": {"human alphaherpesvirus 2", "herpes simplex virus 2"},
    "severe acute respiratory syndrome related coronavirus": {
        "sars-cov-2",
        "coronavirus",
    },
    "influenza a virus": {"influenza a"},
    "rhinovirus a89": {"human rhinovirus/enterovirus", "human rhinovirus enterovirus"},
    "acinetobacter baumannii": {
        "acineto calc baumannii complex",
        "acinetobacter baumannii complex",
    },
    "acineto calc baumannii complex": {
        "acinetobacter baumannii",
        "acinetobacter baumannii complex",
    },
    "klebsiella pneumoniae group": {"klebsiella pneumoniae"},
}

CANONICAL_NAME_BY_NORM = {
    "acineto calc baumannii complex": "acinetobacter baumannii",
    "acinetobacter baumannii complex": "acinetobacter baumannii",
    "klebsiella pneumoniae group": "klebsiella pneumoniae",
    "influenza a virus": "influenza a",
    "human rhinovirus enterovirus": "human rhinovirus/enterovirus",
    "rhinovirus a89": "human rhinovirus/enterovirus",
    "human respirovirus 3": "parainfluenza virus",
    "human parainfluenza virus 3": "parainfluenza virus",
    "human betaherpesvirus 5": "cmv",
    "cytomegalovirus": "cmv",
    "human gammaherpesvirus 4": "ebv",
    "epstein barr virus": "ebv",
    "human alphaherpesvirus 1": "hsv-1",
    "herpes simplex virus 1": "hsv-1",
    "human alphaherpesvirus 2": "hsv-2",
    "herpes simplex virus 2": "hsv-2",
}


@dataclass
class TimedRecord:
    record_type: str
    source_file: str
    index: int
    record: dict[str, Any]
    record_time: datetime | None
    time_field: str | None
    microbe_name: str | None
    normalized_name: str | None
    result: str | None
    is_positive: bool
    is_candidate_eligible: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build candidate-window/all-time evidence bundles from patient JSON dirs."
    )
    parser.add_argument(
        "--input-root",
        required=True,
        type=Path,
        help="Root containing NGS_patient_*_json directories.",
    )
    parser.add_argument(
        "--output-root",
        required=True,
        type=Path,
        help="Directory where per-patient bundles and manifest are written.",
    )
    parser.add_argument(
        "--candidate-window-days",
        type=float,
        default=2,
        help="Days before/after anchor time used to select candidate microbes.",
    )
    parser.add_argument(
        "--window-mode",
        choices=("elapsed", "date"),
        default="elapsed",
        help=(
            "elapsed uses exact +/- N*24 hours. date uses the anchor calendar "
            "date +/- N whole dates, inclusive."
        ),
    )
    parser.add_argument(
        "--evidence-scope",
        choices=("all_time",),
        default="all_time",
        help="Evidence records to include for selected candidates.",
    )
    parser.add_argument(
        "--mngs-root",
        type=Path,
        default=None,
        help=(
            "Optional root containing NGS_patient_*_all_RK_NTC_microbes.json "
            "or mNGS_grouped files. Used only for anchor/candidate supplementation."
        ),
    )
    parser.add_argument(
        "--anchor-source",
        choices=("auto", "input_dates", "mngs"),
        default="auto",
        help="How to choose the candidate window anchor.",
    )
    parser.add_argument(
        "--include-mngs-pathogens",
        action="store_true",
        help=(
            "Also add every pathogen from the selected mNGS record as a candidate. "
            "By default, mNGS is used only to choose the anchor time."
        ),
    )
    parser.add_argument(
        "--selected-root",
        type=Path,
        default=None,
        help=(
            "Optional directory containing NGS_patient_<ID>_selected_pathogens.json "
            "files produced by extract_selected_pathogens.py."
        ),
    )
    parser.add_argument(
        "--candidate-source",
        choices=("window", "selected", "upstream_all", "union"),
        default="window",
        help=(
            "window uses positive tests in the anchor window; selected uses only "
            "upstream picked pathogens; upstream_all uses every frozen upstream "
            "candidate; union combines window and upstream_all."
        ),
    )
    parser.add_argument(
        "--selected-missing-policy",
        choices=("error", "skip"),
        default="error",
        help="How to handle an input patient without a selected-pathogens file.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output.",
    )
    return parser.parse_args()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as f:
        return json.load(f)


def selected_file(selected_root: Path, patient_id: int) -> Path:
    return selected_root / f"NGS_patient_{patient_id}_selected_pathogens.json"


def load_selected_bundle(selected_root: Path, patient_id: int) -> dict[str, Any]:
    path = selected_file(selected_root, patient_id)
    if not path.is_file():
        raise FileNotFoundError(f"selected-pathogens input is missing: {path}")
    data = load_json(path)
    if not isinstance(data, dict):
        raise ValueError(f"selected-pathogens input must be an object: {path}")
    if str(data.get("patient_id") or "") != str(patient_id):
        raise ValueError(f"selected-pathogens patient_id does not match folder: {path}")
    if not isinstance(data.get("selected_pathogens"), list):
        raise ValueError(f"selected_pathogens must be an array: {path}")
    if not isinstance(data.get("candidate_pathogens"), list):
        raise ValueError(f"candidate_pathogens must be an array: {path}")
    data["_source_file"] = str(path.resolve())
    return data


def write_json(path: Path, data: Any, pretty: bool) -> None:
    with path.open("w", encoding="utf-8") as f:
        if pretty:
            json.dump(data, f, ensure_ascii=False, indent=2)
            f.write("\n")
        else:
            json.dump(data, f, ensure_ascii=False, separators=(",", ":"))
            f.write("\n")


def parse_datetime(value: Any) -> datetime | None:
    if not isinstance(value, str):
        return None
    match = DATE_RE.search(value)
    if not match:
        return None
    text = match.group(0).replace("/", "-").replace("T", " ")
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None


def normalize_name(name: Any) -> str | None:
    if not isinstance(name, str):
        return None
    value = name.strip()
    if not value:
        return None
    value = value.replace("[", "").replace("]", "")
    value = re.sub(r"\([^)]*\)", " ", value)
    value = re.sub(r"[^A-Za-z0-9/+.-]+", " ", value)
    value = re.sub(r"\s+", " ", value).strip().lower()
    return value or None


def canonical_key(name: str) -> str:
    norm = normalize_name(name) or name.lower()
    return CANONICAL_NAME_BY_NORM.get(norm, norm)


def is_resistance_marker(name: str | None) -> bool:
    if not name:
        return False
    norm = normalize_name(name) or ""
    return norm in RESISTANCE_MARKERS


def is_specific_microbe(name: str | None) -> bool:
    if not name:
        return False
    norm = normalize_name(name) or ""
    if not norm or is_resistance_marker(name):
        return False
    return not any(pattern in norm for pattern in NON_SPECIFIC_MICROBE_PATTERNS)


def is_detected_result(result: Any) -> bool:
    if result is None:
        return False
    text = str(result).strip().lower()
    if not text:
        return False
    negative_terms = ("not detected", "negative", "not isolated", "no growth")
    if any(term in text for term in negative_terms):
        return False
    positive_terms = ("detected", "isolated", "positive", "present")
    return any(term in text for term in positive_terms)


def choose_record_time(record: dict[str, Any]) -> tuple[datetime | None, str | None]:
    for field in DATE_FIELDS:
        parsed = parse_datetime(record.get(field))
        if parsed:
            return parsed, field
    return None, None


def iter_patient_dirs(input_root: Path) -> Iterable[tuple[int, Path]]:
    dirs = []
    for path in input_root.iterdir():
        if not path.is_dir():
            continue
        match = PATIENT_DIR_RE.match(path.name)
        if match:
            dirs.append((int(match.group(1)), path))
    yield from sorted(dirs, key=lambda item: item[0])


def record_type_from_path(path: Path, patient_id: int) -> str:
    match = PATIENT_FILE_RE.match(path.name)
    if not match:
        return path.stem
    prefix = f"NGS_patient_{patient_id}_"
    stem = path.stem
    if stem.startswith(prefix):
        return stem[len(prefix) :]
    return match.group(2)


def flatten_records(data: Any) -> list[dict[str, Any]]:
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    if isinstance(data, dict):
        return [data]
    return []


def extract_mngs_pathogens(record: dict[str, Any]) -> list[dict[str, Any]]:
    pathogens = record.get("pathogens")
    if not isinstance(pathogens, dict):
        return []
    out = []
    for category, items in pathogens.items():
        if not isinstance(items, list):
            continue
        for item in items:
            if not isinstance(item, dict):
                continue
            name = item.get("name")
            if not is_specific_microbe(name):
                continue
            out.append(
                {
                    "name": name,
                    "category": category,
                    "reads": item.get("reads"),
                    "rk_ntc": item.get("rk_ntc"),
                }
            )
    return out


def extract_timed_records(patient_dir: Path, patient_id: int) -> tuple[list[TimedRecord], dict[str, Any]]:
    timed_records: list[TimedRecord] = []
    raw_by_type: dict[str, Any] = {}
    for path in sorted(patient_dir.glob("*.json")):
        record_type = record_type_from_path(path, patient_id)
        data = load_json(path)
        raw_by_type[record_type] = data
        records = flatten_records(data)
        for idx, record in enumerate(records):
            record_time, time_field = choose_record_time(record)
            microbe_name: str | None = None
            result: str | None = None
            is_positive = False
            is_candidate_eligible = False

            if record_type == "culture":
                microbe_name = record.get("organism")
                result = record.get("status")
                is_positive = is_specific_microbe(microbe_name) and (
                    str(record.get("status", "")).strip().lower() == "isolated"
                    or bool(record.get("colony_count"))
                )
                is_candidate_eligible = is_positive
            elif record_type in ("filmarray", "molecular_microbiology"):
                microbe_name = record.get("target")
                result = record.get("result")
                is_positive = is_specific_microbe(microbe_name) and is_detected_result(result)
                is_candidate_eligible = is_positive
            elif record_type in ("mngs_name_reads",):
                microbe_name = record.get("name")
                result = "mNGS"
                is_positive = is_specific_microbe(microbe_name)
                is_candidate_eligible = is_positive
            elif record_type in ("mNGS_grouped", "mngs_grouped", "all_RK_NTC_microbes"):
                # Keep the grouped record itself as a timed record. Individual
                # pathogens are expanded elsewhere for candidate extraction.
                result = "mNGS"
                is_positive = True
                is_candidate_eligible = False

            timed_records.append(
                TimedRecord(
                    record_type=record_type,
                    source_file=path.name,
                    index=idx,
                    record=record,
                    record_time=record_time,
                    time_field=time_field,
                    microbe_name=microbe_name,
                    normalized_name=normalize_name(microbe_name),
                    result=result,
                    is_positive=is_positive,
                    is_candidate_eligible=is_candidate_eligible,
                )
            )
    return timed_records, raw_by_type


def load_mngs_records(mngs_root: Path | None) -> list[dict[str, Any]]:
    if not mngs_root or not mngs_root.exists():
        return []
    records: list[dict[str, Any]] = []
    for path in sorted(mngs_root.glob("*.json")):
        match = re.search(r"NGS_patient_(\d+)_", path.name)
        if not match:
            continue
        source_patient_id = int(match.group(1))
        try:
            data = load_json(path)
        except json.JSONDecodeError:
            continue
        for idx, record in enumerate(flatten_records(data)):
            record_time, time_field = choose_record_time(record)
            records.append(
                {
                    "source_patient_id": source_patient_id,
                    "source_file": str(path),
                    "source_file_name": path.name,
                    "index": idx,
                    "record": record,
                    "record_time": record_time,
                    "time_field": time_field,
                    "pathogens": extract_mngs_pathogens(record),
                }
            )
    return records


def collect_positive_microbiology_dates(records: list[TimedRecord]) -> list[datetime]:
    return sorted(
        record.record_time
        for record in records
        if record.record_time
        and record.is_positive
        and record.record_type in MICROBE_RECORD_TYPES
    )


def midpoint_time(times: list[datetime]) -> datetime | None:
    if not times:
        return None
    sorted_times = sorted(times)
    if len(sorted_times) % 2 == 1:
        return sorted_times[len(sorted_times) // 2]
    left = sorted_times[len(sorted_times) // 2 - 1]
    right = sorted_times[len(sorted_times) // 2]
    return left + (right - left) / 2


def find_mngs_anchor(
    patient_records: list[TimedRecord],
    all_mngs_records: list[dict[str, Any]],
    patient_id: int,
) -> tuple[datetime | None, dict[str, Any] | None, str]:
    patient_times = [r.record_time for r in patient_records if r.record_time]
    if not patient_times or not all_mngs_records:
        return None, None, "none"
    min_time = min(patient_times) - timedelta(days=1)
    max_time = max(patient_times) + timedelta(days=1)

    candidates: list[tuple[int, float, dict[str, Any]]] = []
    median = midpoint_time([t for t in patient_times if t])
    for mngs in all_mngs_records:
        if mngs.get("source_patient_id") != patient_id:
            continue
        record_time = mngs.get("record_time")
        if not isinstance(record_time, datetime):
            continue
        if min_time <= record_time <= max_time:
            distance = abs((record_time - median).total_seconds()) if median else 0
            candidates.append((0, distance, mngs))
    if not candidates:
        return None, None, "none"
    candidates.sort(key=lambda item: (item[0], item[1]))
    selected = candidates[0][2]
    return selected.get("record_time"), selected, "mngs_date_overlap"


def choose_anchor(
    records: list[TimedRecord],
    all_mngs_records: list[dict[str, Any]],
    patient_id: int,
    anchor_source: str,
) -> tuple[datetime | None, str, dict[str, Any] | None]:
    if anchor_source in ("auto", "mngs"):
        anchor, mngs_record, source = find_mngs_anchor(records, all_mngs_records, patient_id)
        if anchor or anchor_source == "mngs":
            return anchor, source, mngs_record

    positive_micro_dates = collect_positive_microbiology_dates(records)
    anchor = midpoint_time(positive_micro_dates)
    if anchor:
        return anchor, "positive_microbiology_median", None

    all_dates = sorted(r.record_time for r in records if r.record_time)
    anchor = midpoint_time(all_dates)
    if anchor:
        return anchor, "all_record_dates_median", None

    return None, "none", None


def date_window_day_count(days: float) -> int:
    whole_days = int(days)
    if days != whole_days:
        raise ValueError("--window-mode date requires an integer --candidate-window-days")
    return whole_days


def window_bounds(
    anchor: datetime | None, days: float, window_mode: str
) -> tuple[datetime | None, datetime | None]:
    if anchor is None:
        return None, None
    if window_mode == "date":
        whole_days = date_window_day_count(days)
        start_date = anchor.date() - timedelta(days=whole_days)
        end_date = anchor.date() + timedelta(days=whole_days)
        return (
            datetime.combine(start_date, datetime.min.time()),
            datetime.combine(end_date, datetime.max.time()),
        )
    return anchor - timedelta(days=days), anchor + timedelta(days=days)


def in_window(
    record_time: datetime | None,
    anchor: datetime | None,
    days: float,
    window_mode: str,
) -> bool:
    if record_time is None or anchor is None:
        return False
    start, end = window_bounds(anchor, days, window_mode)
    return start is not None and end is not None and start <= record_time <= end


def make_evidence_ref(record: TimedRecord) -> dict[str, Any]:
    out: dict[str, Any] = {
        "record_type": record.record_type,
        "source_file": record.source_file,
        "record_index": record.index,
        "time": record.record_time.isoformat(sep=" ") if record.record_time else None,
        "time_field": record.time_field,
        "microbe_name": record.microbe_name,
        "result": record.result,
        "record": record.record,
    }
    return out


def make_mngs_evidence_ref(mngs_record: dict[str, Any], pathogen: dict[str, Any] | None = None) -> dict[str, Any]:
    record_time = mngs_record.get("record_time")
    raw_record = mngs_record.get("record")
    evidence_record = dict(pathogen) if pathogen else raw_record
    if pathogen and isinstance(raw_record, dict):
        for field in (
            "specimen_code",
            "specimen_site",
            "test_id",
            "accession_number",
            "accession_no",
            "sample_id",
            "order_id",
            "lab_no",
        ):
            if field in raw_record and field not in evidence_record:
                evidence_record[field] = raw_record[field]
    return {
        "record_type": "mngs",
        "source_file": mngs_record.get("source_file_name"),
        "source_patient_id": mngs_record.get("source_patient_id"),
        "record_index": mngs_record.get("index"),
        "time": record_time.isoformat(sep=" ") if isinstance(record_time, datetime) else None,
        "time_field": mngs_record.get("time_field"),
        "microbe_name": pathogen.get("name") if pathogen else None,
        "result": "mNGS",
        "record": evidence_record,
    }


def names_match(a: str | None, b: str | None) -> bool:
    na = normalize_name(a)
    nb = normalize_name(b)
    if not na or not nb:
        return False
    if na == nb:
        return True
    if na in ALIASES and nb in ALIASES[na]:
        return True
    if nb in ALIASES and na in ALIASES[nb]:
        return True
    if len(na) >= 8 and len(nb) >= 8 and (na in nb or nb in na):
        return True
    return False


def gm_related_to_candidate(candidate_name: str) -> bool:
    norm = normalize_name(candidate_name) or ""
    return "aspergillus" in norm or "aspergill" in norm


def build_candidates(
    records: list[TimedRecord],
    anchor: datetime | None,
    window_days: float,
    window_mode: str,
    selected_mngs_record: dict[str, Any] | None,
    selected_bundle: dict[str, Any] | None,
    candidate_source: str,
) -> dict[str, dict[str, Any]]:
    candidates: dict[str, dict[str, Any]] = {}

    if candidate_source in {"window", "union"}:
        for record in records:
            if not record.is_candidate_eligible or not in_window(
                record.record_time, anchor, window_days, window_mode
            ):
                continue
            if not record.microbe_name or not record.normalized_name:
                continue
            key = canonical_key(record.microbe_name)
            candidate = candidates.setdefault(
                key,
                {
                    "name": record.microbe_name,
                    "normalized_name": key,
                    "observed_names": [],
                    "candidate_sources": [],
                    "all_time_evidence": {},
                    "selection_status": "not_linked_to_upstream",
                    "selection_provenance": None,
                },
            )
            if record.microbe_name not in candidate["observed_names"]:
                candidate["observed_names"].append(record.microbe_name)
            candidate["candidate_sources"].append(make_evidence_ref(record))

    if candidate_source in {"window", "union"} and selected_mngs_record and anchor:
        for pathogen in selected_mngs_record.get("pathogens", []):
            name = pathogen.get("name")
            if not is_specific_microbe(name):
                continue
            key = canonical_key(name)
            candidate = candidates.setdefault(
                key,
                {
                    "name": name,
                    "normalized_name": key,
                    "observed_names": [],
                    "candidate_sources": [],
                    "all_time_evidence": {},
                    "selection_status": "not_linked_to_upstream",
                    "selection_provenance": None,
                },
            )
            if name not in candidate["observed_names"]:
                candidate["observed_names"].append(name)
            candidate["candidate_sources"].append(
                make_mngs_evidence_ref(selected_mngs_record, pathogen)
            )

    if candidate_source in {"selected", "upstream_all", "union"}:
        if selected_bundle is None:
            raise ValueError(f"candidate_source={candidate_source} requires --selected-root")
        source_key = (
            "selected_pathogens" if candidate_source == "selected" else "candidate_pathogens"
        )
        for upstream in selected_bundle.get(source_key, []):
            if not isinstance(upstream, dict):
                continue
            name = upstream.get("canonical_name") or upstream.get("original_name")
            if not isinstance(name, str) or not is_specific_microbe(name):
                continue
            key = canonical_key(name)
            candidate = candidates.setdefault(
                key,
                {
                    "name": name,
                    "normalized_name": key,
                    "observed_names": [],
                    "candidate_sources": [],
                    "all_time_evidence": {},
                    "selection_status": "not_linked_to_upstream",
                    "selection_provenance": None,
                },
            )
            original_name = upstream.get("original_name")
            if isinstance(original_name, str) and original_name not in candidate["observed_names"]:
                candidate["observed_names"].append(original_name)
            candidate["selection_status"] = upstream.get("selection_status") or "unknown"
            candidate["selection_provenance"] = upstream

    return candidates


def group_matching_evidence(
    candidate_name: str,
    records: list[TimedRecord],
    selected_mngs_record: dict[str, Any] | None,
) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        if record.microbe_name and names_match(candidate_name, record.microbe_name):
            grouped[record.record_type].append(make_evidence_ref(record))
        elif record.record_type in {"mngs_grouped", "mNGS_grouped", "all_RK_NTC_microbes"}:
            for pathogen in extract_mngs_pathogens(record.record):
                if not names_match(candidate_name, pathogen.get("name")):
                    continue
                evidence = make_evidence_ref(record)
                evidence["microbe_name"] = pathogen.get("name")
                evidence_record = dict(pathogen)
                for field in (
                    "specimen_code",
                    "specimen_site",
                    "test_id",
                    "accession_number",
                    "accession_no",
                    "sample_id",
                    "order_id",
                    "lab_no",
                ):
                    if field in record.record and field not in evidence_record:
                        evidence_record[field] = record.record[field]
                evidence["record"] = evidence_record
                grouped[record.record_type].append(evidence)
        elif record.record_type == "gm_test" and gm_related_to_candidate(candidate_name):
            grouped[record.record_type].append(make_evidence_ref(record))

    if selected_mngs_record:
        for pathogen in selected_mngs_record.get("pathogens", []):
            if names_match(candidate_name, pathogen.get("name")):
                grouped["mngs"].append(make_mngs_evidence_ref(selected_mngs_record, pathogen))

    return dict(sorted(grouped.items(), key=lambda item: item[0]))


def summarize_context(raw_by_type: dict[str, Any]) -> dict[str, Any]:
    context: dict[str, Any] = {}
    for record_type in sorted(CONTEXT_RECORD_TYPES):
        if record_type not in raw_by_type:
            continue
        data = raw_by_type[record_type]
        records = flatten_records(data)
        context[record_type] = {
            "available": True,
            "record_count": len(records),
            "records": data,
        }
    return context


def summarize_positive_microbiology(records: list[TimedRecord]) -> list[dict[str, Any]]:
    out = []
    for record in records:
        if record.is_positive and record.microbe_name:
            out.append(make_evidence_ref(record))
    return out


def build_patient_bundle(
    patient_id: int,
    patient_dir: Path,
    all_mngs_records: list[dict[str, Any]],
    selected_bundle: dict[str, Any] | None,
    args: argparse.Namespace,
) -> tuple[dict[str, Any], dict[str, Any]]:
    records, raw_by_type = extract_timed_records(patient_dir, patient_id)
    anchor, anchor_source, selected_mngs_record = choose_anchor(
        records, all_mngs_records, patient_id, args.anchor_source
    )
    candidates = build_candidates(
        records,
        anchor,
        args.candidate_window_days,
        args.window_mode,
        selected_mngs_record if args.include_mngs_pathogens else None,
        selected_bundle,
        args.candidate_source,
    )
    window_start, window_end = window_bounds(anchor, args.candidate_window_days, args.window_mode)

    for candidate in candidates.values():
        candidate["all_time_evidence"] = group_matching_evidence(
            candidate["name"], records, selected_mngs_record
        )
        candidate["all_time_evidence_counts"] = {
            record_type: len(items)
            for record_type, items in candidate["all_time_evidence"].items()
        }

    sorted_candidates = sorted(
        candidates.values(),
        key=lambda item: (
            -len(item["candidate_sources"]),
            item["name"].lower(),
        ),
    )
    record_type_counts = Counter(record.record_type for record in records)
    positive_micro = summarize_positive_microbiology(records)

    bundle: dict[str, Any] = {
        "patient_id": patient_id,
        "input_patient_dir": str(patient_dir),
        "parameters": {
            "candidate_window_days": args.candidate_window_days,
            "evidence_scope": args.evidence_scope,
            "anchor_source_requested": args.anchor_source,
            "include_mngs_pathogens": args.include_mngs_pathogens,
            "window_mode": args.window_mode,
            "candidate_source": args.candidate_source,
        },
        "selected_input": (
            {
                "source_file": selected_bundle.get("_source_file"),
                "hospital": selected_bundle.get("hospital"),
                "dataset": selected_bundle.get("dataset"),
                "selection_source": selected_bundle.get("selection_source"),
                "pipeline_version": selected_bundle.get("pipeline_version"),
            }
            if selected_bundle
            else None
        ),
        "anchor": {
            "time": anchor.isoformat(sep=" ") if anchor else None,
            "source": anchor_source,
            "mngs_source_file": (
                selected_mngs_record.get("source_file")
                if selected_mngs_record
                else None
            ),
            "mngs_source_patient_id": (
                selected_mngs_record.get("source_patient_id")
                if selected_mngs_record
                else None
            ),
        },
        "candidate_window": {
            "start": window_start.isoformat(sep=" ") if window_start else None,
            "end": window_end.isoformat(sep=" ") if window_end else None,
        },
        "record_type_counts": dict(sorted(record_type_counts.items())),
        "positive_microbiology_records": positive_micro,
        "all_time_context": summarize_context(raw_by_type),
        "candidates": sorted_candidates,
    }

    patient_summary = {
        "patient_id": patient_id,
        "status": "success",
        "anchor_time": bundle["anchor"]["time"],
        "anchor_source": anchor_source,
        "mngs_source_patient_id": bundle["anchor"]["mngs_source_patient_id"],
        "candidate_count": len(sorted_candidates),
        "candidate_names": [item["name"] for item in sorted_candidates],
        "picked_count": sum(
            item.get("selection_status") == "picked" for item in sorted_candidates
        ),
        "positive_microbiology_record_count": len(positive_micro),
        "record_type_counts": bundle["record_type_counts"],
    }
    return bundle, patient_summary


def main() -> int:
    args = parse_args()
    input_root = args.input_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    mngs_records = load_mngs_records(args.mngs_root.resolve() if args.mngs_root else None)
    selected_root = args.selected_root.resolve() if args.selected_root else None
    if args.candidate_source != "window" and selected_root is None:
        raise SystemExit(f"candidate-source {args.candidate_source} requires --selected-root")
    if selected_root is not None and not selected_root.is_dir():
        raise SystemExit(f"selected root is not a directory: {selected_root}")
    patient_summaries: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []

    for patient_id, patient_dir in iter_patient_dirs(input_root):
        try:
            selected_bundle = None
            if args.candidate_source != "window":
                assert selected_root is not None
                path = selected_file(selected_root, patient_id)
                if not path.is_file() and args.selected_missing_policy == "skip":
                    skipped.append(
                        {
                            "patient_id": patient_id,
                            "patient_dir": str(patient_dir),
                            "reason": "selected_pathogens_file_missing",
                        }
                    )
                    continue
                selected_bundle = load_selected_bundle(selected_root, patient_id)
            bundle, summary = build_patient_bundle(
                patient_id, patient_dir, mngs_records, selected_bundle, args
            )
            output_file = output_root / (
                f"NGS_patient_{patient_id}_candidate_window_"
                f"{args.candidate_window_days:g}d_{args.evidence_scope}_evidence.json"
            )
            write_json(output_file, bundle, args.pretty)
            summary["output_file"] = str(output_file)
            patient_summaries.append(summary)
        except Exception as exc:  # pragma: no cover - manifest captures data issues.
            failures.append(
                {
                    "patient_id": patient_id,
                    "patient_dir": str(patient_dir),
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input_root": str(input_root),
        "output_root": str(output_root),
        "mngs_root": str(args.mngs_root.resolve()) if args.mngs_root else None,
        "selected_root": str(selected_root) if selected_root else None,
        "parameters": {
            "candidate_window_days": args.candidate_window_days,
            "evidence_scope": args.evidence_scope,
            "anchor_source": args.anchor_source,
            "include_mngs_pathogens": args.include_mngs_pathogens,
            "window_mode": args.window_mode,
            "candidate_source": args.candidate_source,
            "selected_missing_policy": args.selected_missing_policy,
        },
        "summary": {
            "patients_processed": len(patient_summaries),
            "patients_failed": len(failures),
            "patients_skipped": len(skipped),
            "total_candidates": sum(item["candidate_count"] for item in patient_summaries),
            "total_picked": sum(item["picked_count"] for item in patient_summaries),
        },
        "patients": patient_summaries,
        "skipped": skipped,
        "failures": failures,
    }
    write_json(output_root / "candidate_window_evidence_manifest.json", manifest, True)

    csv_path = output_root / "candidate_window_evidence_summary.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "patient_id",
                "anchor_time",
                "anchor_source",
                "mngs_source_patient_id",
                "candidate_count",
                "candidate_names",
                "picked_count",
                "positive_microbiology_record_count",
                "output_file",
            ],
        )
        writer.writeheader()
        for item in patient_summaries:
            row = dict(item)
            row["candidate_names"] = "; ".join(item["candidate_names"])
            row.pop("status", None)
            row.pop("record_type_counts", None)
            writer.writerow(row)

    print(
        "processed={processed} failed={failed} total_candidates={candidates} output={output}".format(
            processed=len(patient_summaries),
            failed=len(failures),
            candidates=manifest["summary"]["total_candidates"],
            output=output_root,
        )
    )
    if failures:
        print(f"failures={len(failures)}; see candidate_window_evidence_manifest.json")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

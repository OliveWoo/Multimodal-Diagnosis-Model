"""Build an answer-blind small-agent workspace from a multi-cohort bundle.

Layer 01 contains the 2026-09-25 supplement. Some current handoff bundles have
already reconciled those records into layer 02, while older copies have not.
This builder keeps the normalized layer-02 records and adds only supplemental
v17 ``clinical_result`` records that are not already represented by the same
parser file and JSON pointer. It then applies the inclusive +/-2 calendar-date
event window and writes an isolated flat patient workspace expected by
``tools.analyze_llm_agent``.

Ground truth and prior derived outputs are never read or copied.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import shutil
import unicodedata
from collections import Counter
from datetime import date, datetime
from pathlib import Path
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "2026-10-08_44patients_pipeline_input_with_ground_truth"
DEFAULT_OUTPUT = ROOT / "outputs" / "runs" / "2026-10-09_44patients_small_agents_stage1_v5"

PATIENT_RE = re.compile(r"^NGS_patient_(\d+)$")
RAW_RE = re.compile(r"^NGS_patient_\d+_(.+)_Raw\.json$")
SUPPLEMENT_HEADER = ["925_source_sheet", "925_source_row", "v17_class", "v17_records_json"]
SUPPLEMENT_END = "[END 925 SUPPLEMENT]"

MODULES = (
    "admission_diagnosis",
    "CBC",
    "culture",
    "filmarray",
    "gm_test",
    "image",
    "molecular_microbiology",
    "other_lab",
    "underlying",
)
TEMPORAL_MODULES = {
    "CBC",
    "culture",
    "filmarray",
    "gm_test",
    "image",
    "molecular_microbiology",
    "other_lab",
}
TIME_FIELDS = ("collected_time", "reported_time", "source_time")
SEMANTIC_IGNORED_FIELDS = {
    "source",
    "window_match",
    "record_origin",
    "supplement_provenance",
    "temporal_quality",
    "validation_status",
    "qc",
    "qc_NOTE",
}
PROMPT_OMITTED_FIELDS = {
    "source",
    "supplement_provenance",
    "original_row",
    "qc_assessment",
    "record_class",
    "parser_reference",
    "record_origin",
    "temporal_quality",
}
CBC_HOST_STATE_PATTERN = re.compile(
    r"^(?:wbc|hgb|hemoglobin|hct|hematocrit|plt|platelet|anc|absolute neut|neut|seg|band|lym|lymph|alc|absolute lymph|blast)",
    re.IGNORECASE,
)
HOST_STATE_OTHER_LAB_PATTERNS = tuple(
    re.compile(pattern, re.IGNORECASE)
    for pattern in (
        r"\bcrp\b|c[- ]?reactive|procalcitonin|\bpct\b|lactat|lactic",
        r"creatinine|\burea\b|urea n|\bbun\b|\begfr\b|cystatin",
        r"\bsgot\b|\bsgpt\b|\bast\b|\balt\b|bilirubin|\bbil\b|albumin|\balp\b|\bggt\b|ammonia",
        r"hco3|bicarbonate|pco2|po2|\bph\b|base excess|beecf|\bsbc\b|total co2|fio2|so2|oxygen saturation",
        r"sodium|\bna\b|potassium|\bk\b|chloride|\bcl\b|calcium|ca\+\+|magnesium|\bmg\b|phosph|\bp\b",
        r"pt\s*\(?inr|\binr\b|\bpt\b|\bptt\b|\baptt\b|d[- ]?dimer|fibrinogen",
        r"\bldh\b|ferritin|troponin|tropo|ntpro[- ]?bnp|\bbnp\b|ck[- ]?mb|\bcpk\b",
        r"temperature|fever|glucose|sugar|\bglu\b",
    )
)


def load_json(path: Path) -> Any:
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


def parse_date(value: Any) -> date | None:
    if value in (None, ""):
        return None
    text = str(value).strip().replace("/", "-")
    if not text:
        return None
    try:
        return datetime.fromisoformat(text.replace("Z", "+00:00")).date()
    except ValueError:
        pass
    for fmt in ("%Y-%m-%d", "%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(text, fmt).date()
        except ValueError:
            continue
    return None


def parse_datetime(value: Any) -> datetime | None:
    if value in (None, ""):
        return None
    text = str(value).strip().replace("/", "-")
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
        return parsed.replace(tzinfo=None) if parsed.tzinfo is not None else parsed
    except ValueError:
        pass
    for fmt in ("%Y-%m-%d", "%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None


def patient_dirs(source_root: Path) -> list[tuple[str, int, Path]]:
    rows: list[tuple[str, int, Path]] = []
    patients_root = source_root / "patients"
    for cohort_dir in sorted(path for path in patients_root.iterdir() if path.is_dir()):
        for patient_dir in sorted(path for path in cohort_dir.iterdir() if path.is_dir()):
            match = PATIENT_RE.fullmatch(patient_dir.name)
            if match:
                rows.append((cohort_dir.name, int(match.group(1)), patient_dir))
    return sorted(rows, key=lambda row: (row[1], row[0]))


def event_anchors(mngs_path: Path) -> list[dict[str, Any]]:
    payload = load_json(mngs_path)
    anchors: list[dict[str, Any]] = []
    for index, case in enumerate(payload.get("cases") or [], start=1):
        if not isinstance(case, dict):
            continue
        collected_time = case.get("collected_time")
        collected_date = parse_date(collected_time)
        if collected_date is None:
            raise ValueError(f"Missing/unparseable mNGS collected_time in {mngs_path}: {collected_time!r}")
        anchors.append({
            "anchor_id": f"mngs_{index}",
            "case_review_id": case.get("case_review_id"),
            "collected_time": collected_time,
            "date": collected_date,
            "datetime": parse_datetime(collected_time),
        })
    if not anchors:
        raise ValueError(f"No mNGS cases found in {mngs_path}")
    return anchors


def record_event_date(record: dict[str, Any]) -> tuple[str | None, date | None]:
    for field in TIME_FIELDS:
        parsed = parse_date(record.get(field))
        if parsed is not None:
            return field, parsed
    return None, None


def match_window(record_date: date, anchors: Sequence[dict[str, Any]]) -> dict[str, Any] | None:
    matches: list[tuple[dict[str, Any], int]] = []
    for anchor in anchors:
        offset = (record_date - anchor["date"]).days
        if -2 <= offset <= 2:
            matches.append((anchor, offset))
    if not matches:
        return None
    return {
        "policy": "inclusive_plus_minus_2_calendar_dates",
        "event_date": record_date.isoformat(),
        "anchor_ids": [anchor["anchor_id"] for anchor, _ in matches],
        "anchor_case_review_ids": [anchor.get("case_review_id") for anchor, _ in matches],
        "calendar_day_offsets": [offset for _, offset in matches],
    }


def semantic_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: semantic_value(item)
            for key, item in sorted(value.items())
            if key not in SEMANTIC_IGNORED_FIELDS
        }
    if isinstance(value, list):
        return [semantic_value(item) for item in value]
    if isinstance(value, str):
        return " ".join(value.split())
    return value


def semantic_key(record: dict[str, Any]) -> str:
    return json.dumps(semantic_value(record), ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def parser_source_key(record: dict[str, Any]) -> tuple[str, str] | None:
    """Return the stable parser identity shared by layer-01 and layer-02 records."""
    for field in ("parser_reference", "supplement_provenance"):
        reference = record.get(field)
        if not isinstance(reference, dict):
            continue
        parser_file = str(reference.get("parser_file") or "").strip()
        json_pointer = str(reference.get("json_pointer") or "").strip()
        if parser_file and json_pointer:
            return parser_file, json_pointer
    return None


def normalized_analyte(record: dict[str, Any]) -> str:
    value = (
        record.get("item")
        or record.get("item_full")
        or record.get("target")
        or record.get("test")
        or ""
    )
    text = unicodedata.normalize("NFKC", str(value)).lower()
    return " ".join(text.split())


def host_state_relevant_other_lab(record: dict[str, Any]) -> bool:
    text = " ".join(
        normalized_analyte(record).split()
        + [str(record.get("item_full") or "").lower(), str(record.get("unit") or "").lower()]
    )
    return any(pattern.search(text) for pattern in HOST_STATE_OTHER_LAB_PATTERNS)


def numeric_value(record: dict[str, Any]) -> float | None:
    value = record.get("value")
    if value in (None, ""):
        return None
    match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)", str(value).replace(",", ""))
    if match is None:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def record_datetime(record: dict[str, Any]) -> datetime | None:
    for field in TIME_FIELDS:
        parsed = parse_datetime(record.get(field))
        if parsed is not None:
            return parsed
    return None


def trajectory_indices(
    records: Sequence[dict[str, Any]], anchors: Sequence[dict[str, Any]]
) -> set[int]:
    if len(records) <= 4:
        return set(range(len(records)))
    selected: set[int] = set()
    numeric = [(index, numeric_value(record)) for index, record in enumerate(records)]
    numeric = [(index, value) for index, value in numeric if value is not None]
    if numeric:
        selected.add(min(numeric, key=lambda row: row[1])[0])
        selected.add(max(numeric, key=lambda row: row[1])[0])
    timed = [(index, record_datetime(record)) for index, record in enumerate(records)]
    timed = [(index, value) for index, value in timed if value is not None]
    for anchor in anchors:
        anchor_dt = anchor.get("datetime")
        if anchor_dt is None or not timed:
            continue
        before = [row for row in timed if row[1] <= anchor_dt]
        after = [row for row in timed if row[1] >= anchor_dt]
        if before:
            selected.add(max(before, key=lambda row: row[1])[0])
        if after:
            selected.add(min(after, key=lambda row: row[1])[0])
    if not selected:
        selected.update({0, len(records) - 1})
    return selected


def select_lab_records_for_prompt(
    records: Sequence[Any], module: str, anchors: Sequence[dict[str, Any]]
) -> tuple[list[Any], list[dict[str, Any]]]:
    if module not in {"CBC", "other_lab"}:
        return list(records), []
    eligible: list[tuple[int, dict[str, Any]]] = []
    omitted: list[dict[str, Any]] = []
    for index, item in enumerate(records):
        if not isinstance(item, dict):
            omitted.append({"record_index": index, "reason": "non_object_record", "record": item})
            continue
        if module == "CBC" and not CBC_HOST_STATE_PATTERN.search(normalized_analyte(item)):
            omitted.append({"record_index": index, "reason": "outside_host_state_cbc_contract", "record": item})
            continue
        if module == "other_lab" and not host_state_relevant_other_lab(item):
            omitted.append({"record_index": index, "reason": "outside_host_state_lab_contract", "record": item})
            continue
        eligible.append((index, item))

    groups: dict[str, list[tuple[int, dict[str, Any]]]] = {}
    for original_index, record in eligible:
        key = f"{normalized_analyte(record)}|{str(record.get('unit') or '').lower()}"
        groups.setdefault(key, []).append((original_index, record))
    selected_indices: set[int] = set()
    for rows in groups.values():
        rows.sort(key=lambda row: (record_datetime(row[1]) or datetime.min, row[0]))
        local_records = [record for _, record in rows]
        for local_index in trajectory_indices(local_records, anchors):
            selected_indices.add(rows[local_index][0])
    selected: list[Any] = []
    for index, item in enumerate(records):
        if index in selected_indices:
            selected.append(item)
        elif any(index == row[0] for row in eligible):
            omitted.append({"record_index": index, "reason": "trajectory_compaction", "record": item})
    return selected, omitted


def compact_value_for_prompt(value: Any) -> Any:
    if isinstance(value, list):
        return [compact_value_for_prompt(item) for item in value]
    if not isinstance(value, dict):
        return value
    compact: dict[str, Any] = {}
    for key, item in value.items():
        if key in PROMPT_OMITTED_FIELDS:
            continue
        if key == "window_match" and isinstance(item, dict):
            offsets = item.get("calendar_day_offsets")
            if offsets not in (None, [], ""):
                compact["relative_days_to_mngs"] = offsets
            continue
        compact[key] = compact_value_for_prompt(item)
    return compact


def compact_records_for_prompt(
    records: Sequence[Any], patient_name: str, module: str
) -> tuple[list[Any], list[dict[str, Any]]]:
    compact_records: list[Any] = []
    provenance: list[dict[str, Any]] = []
    for index, item in enumerate(records, start=1):
        if not isinstance(item, dict):
            compact_records.append(item)
            continue
        source_ref = f"P{patient_name.rsplit('_', 1)[-1]}:{module[:3]}:{index:04d}"
        compact = compact_value_for_prompt(item)
        compact["source_ref"] = source_ref
        compact_records.append(compact)
        provenance.append({
            "source_ref": source_ref,
            "record_origin": item.get("record_origin", "existing_pipeline_direct_input"),
            "source": item.get("source"),
            "parser_reference": item.get("parser_reference"),
            "supplement_provenance": item.get("supplement_provenance"),
            "window_match": item.get("window_match"),
            "merge_conflicts": item.get("merge_conflicts"),
        })
    return compact_records, provenance


def iter_supplement_rows(raw_path: Path) -> Iterable[tuple[int, list[Any]]]:
    payload = load_json(raw_path)
    rows = payload.get("raw_data") if isinstance(payload, dict) else None
    if not isinstance(rows, list):
        return
    inside = False
    for index, row in enumerate(rows):
        if not isinstance(row, list) or not row:
            continue
        if row[:4] == SUPPLEMENT_HEADER:
            inside = True
            continue
        if row[0] == SUPPLEMENT_END:
            break
        if inside and len(row) >= 4:
            yield index, row


def supplement_records(
    raw_path: Path,
    module: str,
    anchors: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]], Counter[str]]:
    accepted: list[dict[str, Any]] = []
    audit: dict[str, list[dict[str, Any]]] = {
        "undated": [],
        "outside_window": [],
        "excluded_class": [],
        "parse_error": [],
    }
    counts: Counter[str] = Counter()
    for raw_index, row in iter_supplement_rows(raw_path):
        row_class = str(row[2])
        try:
            parsed_entries = json.loads(row[3])
        except (TypeError, json.JSONDecodeError) as exc:
            audit["parse_error"].append({"raw_index": raw_index, "error": str(exc)})
            counts["parse_error"] += 1
            continue
        if not isinstance(parsed_entries, list):
            parsed_entries = [parsed_entries]
        for entry in parsed_entries:
            if not isinstance(entry, dict):
                counts["non_object_entry"] += 1
                continue
            entry_class = str(entry.get("class") or row_class)
            if entry_class != "clinical_result":
                audit["excluded_class"].append({
                    "raw_index": raw_index,
                    "row_class": row_class,
                    "entry_class": entry_class,
                    "json_pointer": entry.get("json_pointer"),
                })
                counts[f"excluded_{entry_class}"] += 1
                continue
            source_record = entry.get("record")
            if not isinstance(source_record, dict):
                counts["missing_record"] += 1
                continue
            record = copy.deepcopy(source_record)
            time_field, event_date = record_event_date(record)
            if module in TEMPORAL_MODULES:
                if event_date is None:
                    audit["undated"].append({
                        "raw_index": raw_index,
                        "row_class": row_class,
                        "json_pointer": entry.get("json_pointer"),
                        "record": record,
                    })
                    counts["undated"] += 1
                    continue
                window = match_window(event_date, anchors)
                if window is None:
                    audit["outside_window"].append({
                        "raw_index": raw_index,
                        "row_class": row_class,
                        "json_pointer": entry.get("json_pointer"),
                        "time_field": time_field,
                        "event_date": event_date.isoformat(),
                    })
                    counts["outside_window"] += 1
                    continue
                window["time_field"] = time_field
                record["window_match"] = window
                if time_field == "source_time":
                    record["temporal_quality"] = "explicit_date_role_unspecified_event_context_only"
            record["record_origin"] = "925_v17_supplement"
            record["supplement_provenance"] = {
                "raw_file": raw_path.name,
                "raw_data_index": raw_index,
                "source_sheet": row[0],
                "source_row": row[1],
                "parser_file": entry.get("parser_file"),
                "json_pointer": entry.get("json_pointer"),
            }
            accepted.append(record)
            counts["accepted_before_dedup"] += 1
    return accepted, audit, counts


def merge_records(existing: list[Any], additions: Iterable[dict[str, Any]]) -> tuple[list[Any], int]:
    merged = list(existing)
    seen_semantic = {semantic_key(item) for item in existing if isinstance(item, dict)}
    seen_parser_sources = {
        source_key
        for item in existing
        if isinstance(item, dict)
        for source_key in [parser_source_key(item)]
        if source_key is not None
    }
    duplicate_count = 0
    for record in additions:
        source_key = parser_source_key(record)
        key = semantic_key(record)
        if (source_key is not None and source_key in seen_parser_sources) or key in seen_semantic:
            duplicate_count += 1
            continue
        seen_semantic.add(key)
        if source_key is not None:
            seen_parser_sources.add(source_key)
        merged.append(record)
    return merged, duplicate_count


def _image_text_tokens(record: dict[str, Any]) -> set[str]:
    findings = record.get("findings")
    values = list(findings) if isinstance(findings, list) else [findings] if findings else []
    impression = record.get("impression")
    if impression:
        values.append(impression)
    text = unicodedata.normalize("NFKC", " ".join(str(value) for value in values)).casefold()
    return set(re.findall(r"\w+", text, flags=re.UNICODE))


def source_aware_image_prompt_records(
    records: Sequence[Any],
) -> tuple[list[Any], list[dict[str, Any]]]:
    """Drop source-time-only copies when a richer record describes the same study."""

    object_rows = [item for item in records if isinstance(item, dict)]
    omitted_ids: set[int] = set()
    omitted: list[dict[str, Any]] = []
    for low in object_rows:
        if low.get("collected_time") or low.get("reported_time"):
            continue
        low_date = parse_date(low.get("source_time"))
        low_tokens = _image_text_tokens(low)
        if low_date is None or not low_tokens:
            continue
        low_exam = " ".join(
            unicodedata.normalize("NFKC", str(low.get("exam_type") or ""))
            .casefold()
            .split()
        )
        for rich in object_rows:
            if rich is low or not (rich.get("collected_time") or rich.get("reported_time")):
                continue
            _, rich_date = record_event_date(rich)
            if rich_date != low_date:
                continue
            rich_exam = " ".join(
                unicodedata.normalize("NFKC", str(rich.get("exam_type") or ""))
                .casefold()
                .split()
            )
            if low_exam and rich_exam and low_exam != rich_exam:
                continue
            rich_tokens = _image_text_tokens(rich)
            if not rich_tokens:
                continue
            containment = len(low_tokens & rich_tokens) / min(len(low_tokens), len(rich_tokens))
            if containment < 0.8:
                continue
            omitted_ids.add(id(low))
            omitted.append(
                {
                    "reason": "source_time_only_duplicate_of_richer_image_record",
                    "source_time": low.get("source_time"),
                    "exam_type": low.get("exam_type"),
                    "matched_collected_time": rich.get("collected_time"),
                    "matched_reported_time": rich.get("reported_time"),
                    "content_token_containment": round(containment, 4),
                    "record": low,
                }
            )
            break
    return [item for item in records if id(item) not in omitted_ids], omitted


def build_patient(
    cohort: str,
    patient_id: int,
    source_patient: Path,
    output_root: Path,
) -> dict[str, Any]:
    patient_name = f"NGS_patient_{patient_id}"
    layer_01 = source_patient / "01_upstream_raw_clinical_sheet_exports"
    layer_02 = source_patient / "02_pipeline_direct_clinical_inputs"
    layer_03 = source_patient / "03_raw_mngs_per_test"
    mngs_path = layer_03 / f"{patient_name}_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
    if not mngs_path.is_file():
        raise FileNotFoundError(mngs_path)
    anchors = event_anchors(mngs_path)

    output_patient = output_root / "patients" / f"{patient_name}_json"
    output_patient.mkdir(parents=True, exist_ok=False)
    (output_patient / "agent_outputs").mkdir()
    (output_patient / "summary_outputs").mkdir()
    shutil.copy2(mngs_path, output_patient / mngs_path.name)

    module_rows: dict[str, Any] = {}
    patient_audit: dict[str, Any] = {}
    patient_provenance: dict[str, list[dict[str, Any]]] = {}
    patient_counts: Counter[str] = Counter()
    for module in MODULES:
        direct_path = layer_02 / f"{patient_name}_{module}.json"
        existing = load_json(direct_path) if direct_path.is_file() else []
        if not isinstance(existing, list):
            raise TypeError(f"Expected list in {direct_path}")

        additions: list[dict[str, Any]] = []
        audit = {"undated": [], "outside_window": [], "excluded_class": [], "parse_error": []}
        raw_path = layer_01 / f"{patient_name}_{module}_Raw.json"
        counts: Counter[str] = Counter()
        if raw_path.is_file():
            additions, audit, counts = supplement_records(raw_path, module, anchors)
        merged, duplicate_count = merge_records(existing, additions)
        pre_prompt_records: list[Any] = merged
        source_duplicate_omitted: list[dict[str, Any]] = []
        if module == "image":
            pre_prompt_records, source_duplicate_omitted = source_aware_image_prompt_records(
                merged
            )
        prompt_records, prompt_omitted = select_lab_records_for_prompt(
            pre_prompt_records, module, anchors
        )
        prompt_omitted = [*source_duplicate_omitted, *prompt_omitted]
        compact_merged, provenance = compact_records_for_prompt(prompt_records, patient_name, module)
        output_path = output_patient / f"{patient_name}_{module}.json"
        write_json(output_path, compact_merged)
        patient_provenance[module] = provenance
        module_rows[module] = {
            "existing_record_count": len(existing),
            "supplement_accepted_before_dedup": len(additions),
            "supplement_duplicate_count": duplicate_count,
            "output_record_count": len(merged),
            "prompt_record_count": len(compact_merged),
            "prompt_omitted_record_count": len(prompt_omitted),
            "source_aware_duplicate_count": len(source_duplicate_omitted),
            "output_sha256": sha256_file(output_path),
        }
        module_rows[module].update(dict(sorted(counts.items())))
        if any(audit.values()):
            patient_audit[module] = audit
        if prompt_omitted:
            patient_audit.setdefault(module, {})["prompt_omitted"] = prompt_omitted
        patient_counts["existing_records"] += len(existing)
        patient_counts["supplement_accepted_before_dedup"] += len(additions)
        patient_counts["supplement_duplicates"] += duplicate_count
        patient_counts["output_records"] += len(merged)
        patient_counts["prompt_records"] += len(compact_merged)
        patient_counts["prompt_omitted_records"] += len(prompt_omitted)
        patient_counts.update(counts)

    audit_path = output_root / "audit" / f"{patient_name}_supplement_filter_audit.json"
    write_json(audit_path, {
        "schema_version": "small_agent_supplement_filter_audit.v1",
        "patient_id": patient_id,
        "cohort": cohort,
        "window_policy": "inclusive_plus_minus_2_calendar_dates",
        "modules": patient_audit,
    })
    provenance_path = output_root / "audit" / f"{patient_name}_input_provenance.json"
    write_json(provenance_path, {
        "schema_version": "small_agent_input_provenance.v1",
        "patient_id": patient_id,
        "cohort": cohort,
        "interpretation": "source_ref resolves prompt records to omitted source metadata; this file is not sent to the LLM",
        "modules": patient_provenance,
    })
    manifest = {
        "schema_version": "small_agent_patient_input.v1",
        "patient_id": patient_id,
        "cohort": cohort,
        "source_patient": str(source_patient.resolve()),
        "window_policy": "inclusive_plus_minus_2_calendar_dates",
        "mngs_anchors": [
            {key: value for key, value in anchor.items() if key not in {"date", "datetime"}}
            for anchor in anchors
        ],
        "ground_truth_accessed": False,
        "derived_output_reference_accessed": False,
        "modules": module_rows,
        "audit_path": str(audit_path.resolve()),
        "provenance_path": str(provenance_path.resolve()),
    }
    write_json(output_patient / f"{patient_name}_small_agent_input_manifest.json", manifest)
    return {
        "patient_id": patient_id,
        "cohort": cohort,
        "patient_dir": str(output_patient.resolve()),
        "mngs_case_count": len(anchors),
        **dict(patient_counts),
    }


def ground_truth_firewall(output_root: Path) -> list[str]:
    violations: list[str] = []
    for path in output_root.rglob("*"):
        relative = str(path.relative_to(output_root)).lower()
        if "ground_truth" in relative or "clinical_pathogen" in relative:
            violations.append(str(path.resolve()))
    return violations


def run(source_root: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output root: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    rows = [
        build_patient(cohort, patient_id, source_patient, output_root)
        for cohort, patient_id, source_patient in patient_dirs(source_root)
    ]
    violations = ground_truth_firewall(output_root)
    if violations:
        raise RuntimeError(f"Ground-truth firewall violation: {violations}")
    totals: Counter[str] = Counter()
    for row in rows:
        for key, value in row.items():
            if isinstance(value, int) and key != "patient_id":
                totals[key] += value
    summary = {
        "schema_version": "small_agent_workspace_summary.v1",
        "source_root": str(source_root.resolve()),
        "output_root": str(output_root.resolve()),
        "patient_count": len(rows),
        "cohort_counts": dict(Counter(row["cohort"] for row in rows)),
        "patient_ids": [row["patient_id"] for row in rows],
        "window_policy": "inclusive_plus_minus_2_calendar_dates",
        "ground_truth_firewall": "passed",
        "ground_truth_files_copied": 0,
        "derived_reference_files_copied": 0,
        "totals": dict(sorted(totals.items())),
        "patients": rows,
        "next_stage": "small_agent_prompt_preflight",
    }
    write_json(output_root / "workspace_summary.json", summary)
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    print(json.dumps(run(args.source_root, args.output_root), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

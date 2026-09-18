from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import re
import sys
import time
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from input_identity import load_mapping, normalized_code, read_frozen_json, sha256_file


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_MANIFEST = SCRIPT_DIR / "cohort_manifest.json"
DEFAULT_PROMPT = SCRIPT_DIR / "prompt_v1.txt"
DEFAULT_AUDIT = SCRIPT_DIR / "audit" / "direct_raw_input_audit.json"
DEFAULT_OUTPUT_ROOT = SCRIPT_DIR / "outputs"

CLINICAL_SECTIONS = (
    "admission_diagnosis",
    "underlying",
    "CBC",
    "other_lab",
    "culture",
    "filmarray",
    "gm_test",
    "molecular_microbiology",
    "image",
)

DROP_KEYS = {
    "birth_date",
    "date_of_birth",
    "dob",
    "patient_name",
    "name_of_patient",
    "medical_record_number",
    "mrn",
    "chart_number",
    "national_id",
    "address",
    "phone",
    "telephone",
    "email",
    "gold",
    "gold_answer",
    "ground_truth",
    "answer",
    "final_answer",
    "picked",
    "ober",
    "llm_missed_candidate_review",
    "review_high_priority",
    "review_context_needed",
    "a1",
    "c3",
    "site_alignment",
    "site_align",
}

DATE_KEYS = {
    "admission_date",
    "collected_time",
    "received_time",
    "reported_time",
    "collection_date",
    "report_date",
}

DATE_PATTERNS = (
    re.compile(r"\b(?:19|20)\d{2}[-/](?:0?[1-9]|1[0-2])[-/](?:0?[1-9]|[12]\d|3[01])\b"),
    re.compile(r"\b(?:19|20)\d{6}\b"),
)

DATE_RANGE_PATTERN = re.compile(
    r"\b((?:19|20)\d{2})[-/](0?[1-9]|1[0-2])[-/](0?[1-9]|[12]\d|3[01])"
    r"\s*[-~至]\s*(0?[1-9]|1[0-2])[-/](0?[1-9]|[12]\d|3[01])\b"
)

OUTPUT_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "patient_id": {"type": "string"},
        "infection_assessment": {
            "type": "string",
            "enum": ["supported", "not_supported", "uncertain"],
        },
        "predicted_pathogens": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "organism": {"type": "string"},
                    "confidence": {
                        "type": "string",
                        "enum": ["high", "moderate", "low"],
                    },
                    "supporting_evidence": {
                        "type": "array",
                        "items": {"type": "string"},
                    },
                    "counterevidence": {
                        "type": "array",
                        "items": {"type": "string"},
                    },
                },
                "required": [
                    "organism",
                    "confidence",
                    "supporting_evidence",
                    "counterevidence",
                ],
            },
        },
        "excluded_detected_organisms": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "organism": {"type": "string"},
                    "reason": {"type": "string"},
                },
                "required": ["organism", "reason"],
            },
        },
        "no_supported_pathogen": {"type": "boolean"},
        "limitations": {
            "type": "array",
            "items": {"type": "string"},
        },
    },
    "required": [
        "patient_id",
        "infection_assessment",
        "predicted_pathogens",
        "excluded_detected_organisms",
        "no_supported_pathogen",
        "limitations",
    ],
}


@dataclass(frozen=True)
class PatientSpec:
    patient_id: str
    cohort: str
    clinical_root: Path
    raw_mngs_root: Path
    raw_mngs_patient_id: str
    expected_case_code: str
    expected_specimen_code: str
    raw_mngs_sha256: str
    clinical_hashes: tuple[tuple[str, str | None], ...]
    mapping_sha256: str

    @property
    def numeric_id(self) -> str:
        return self.patient_id.removeprefix("P")

    @property
    def clinical_dir(self) -> Path:
        return self.clinical_root / f"NGS_patient_{self.numeric_id}_json"

    @property
    def raw_mngs_file(self) -> Path:
        raw_id = self.raw_mngs_patient_id.removeprefix("P")
        return self.raw_mngs_root / f"NGS_patient_{raw_id}_all_RK_NTC_microbes.json"


@dataclass
class SanitizeStats:
    dropped_keys: dict[str, int] = field(default_factory=dict)
    pseudonymized_specimen_codes: int = 0
    converted_dates: int = 0
    unparseable_dates_removed: int = 0
    free_text_dates_converted: int = 0

    def add_drop(self, key: str) -> None:
        self.dropped_keys[key] = self.dropped_keys.get(key, 0) + 1


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def stable_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def normalize_key(key: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", key.casefold()).strip("_")


def parse_date_value(value: Any) -> date | None:
    if not isinstance(value, str):
        return None
    text = value.strip()
    if not text:
        return None
    candidates = [text, text[:10].strip()]
    formats = (
        "%Y-%m-%d",
        "%Y/%m/%d",
        "%Y%m%d",
        "%Y/%m/%d %H:%M",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d %H:%M:%S",
        "%Y/%m/%d %H:%M:%S",
    )
    for candidate in candidates:
        for fmt in formats:
            try:
                return datetime.strptime(candidate, fmt).date()
            except ValueError:
                continue
    return None


def find_index_date(raw_mngs: Any) -> date | None:
    found: list[date] = []

    def walk(value: Any) -> None:
        if isinstance(value, dict):
            for key, item in value.items():
                if normalize_key(str(key)) == "collected_time":
                    parsed = parse_date_value(item)
                    if parsed is not None:
                        found.append(parsed)
                walk(item)
        elif isinstance(value, list):
            for item in value:
                walk(item)

    walk(raw_mngs)
    return min(found) if found else None


def relative_day_label(value_date: date, index_date: date) -> str:
    delta = (value_date - index_date).days
    return f"D{delta:+d}"


def replace_free_text_dates(text: str, index_date: date | None, stats: SanitizeStats) -> str:
    def range_replacement(match: re.Match[str]) -> str:
        year, start_month, start_day, end_month, end_day = match.groups()
        start = parse_date_value(f"{year}/{start_month}/{start_day}")
        end = parse_date_value(f"{year}/{end_month}/{end_day}")
        stats.free_text_dates_converted += 2
        if start is None or end is None or index_date is None:
            return "[date range removed]"
        return f"{relative_day_label(start, index_date)} to {relative_day_label(end, index_date)}"

    result = DATE_RANGE_PATTERN.sub(range_replacement, text)
    for pattern in DATE_PATTERNS:
        def replacement(match: re.Match[str]) -> str:
            parsed = parse_date_value(match.group(0))
            stats.free_text_dates_converted += 1
            if parsed is None or index_date is None:
                return "[date removed]"
            return relative_day_label(parsed, index_date)

        previous = None
        while previous != result:
            previous = result
            result = pattern.sub(replacement, result)
    return result


def sanitize_value(
    value: Any,
    *,
    index_date: date | None,
    specimen_map: dict[str, str],
    stats: SanitizeStats,
) -> Any:
    if isinstance(value, dict):
        cleaned: dict[str, Any] = {}
        for raw_key, raw_item in value.items():
            key = str(raw_key)
            normalized = normalize_key(key)
            if normalized in DROP_KEYS:
                stats.add_drop(normalized)
                continue
            if normalized == "specimen_code":
                raw_code = str(raw_item)
                if raw_code not in specimen_map:
                    specimen_map[raw_code] = f"specimen_{len(specimen_map) + 1}"
                cleaned[key] = specimen_map[raw_code]
                stats.pseudonymized_specimen_codes += 1
                continue
            if normalized in DATE_KEYS:
                parsed = parse_date_value(raw_item)
                if parsed is not None and index_date is not None:
                    cleaned[key] = relative_day_label(parsed, index_date)
                    stats.converted_dates += 1
                else:
                    cleaned[key] = "[date unavailable]"
                    stats.unparseable_dates_removed += 1
                continue
            cleaned[key] = sanitize_value(
                raw_item,
                index_date=index_date,
                specimen_map=specimen_map,
                stats=stats,
            )
        return cleaned
    if isinstance(value, list):
        return [
            sanitize_value(
                item,
                index_date=index_date,
                specimen_map=specimen_map,
                stats=stats,
            )
            for item in value
        ]
    if isinstance(value, str):
        return replace_free_text_dates(value, index_date, stats)
    return copy.deepcopy(value)


def exact_date_hits(value: Any) -> list[str]:
    serialized = stable_json(value)
    hits: list[str] = []
    for pattern in DATE_PATTERNS:
        hits.extend(match.group(0) for match in pattern.finditer(serialized))
    return sorted(set(hits))


def resolve_relative(base: Path, configured: str) -> Path:
    path = Path(configured)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def load_manifest(path: Path) -> tuple[dict[str, Any], list[PatientSpec]]:
    manifest = load_json(path)
    base = path.resolve().parent
    raw_mngs_root = resolve_relative(base, manifest["raw_mngs_source"])
    mappings, mapping_hash = load_mapping(manifest, base, CLINICAL_SECTIONS)
    patients: list[PatientSpec] = []
    seen: set[str] = set()
    for cohort in manifest["cohorts"]:
        clinical_root = resolve_relative(base, cohort["clinical_root"])
        for patient_id in cohort["patient_ids"]:
            normalized = patient_id.upper()
            if normalized in seen:
                raise ValueError(f"Duplicate patient in manifest: {normalized}")
            if not re.fullmatch(r"P\d+", normalized):
                raise ValueError(f"Invalid patient ID: {patient_id}")
            seen.add(normalized)
            identity = mappings[(cohort["name"], normalized)]
            patients.append(
                PatientSpec(
                    patient_id=normalized,
                    cohort=cohort["name"],
                    clinical_root=clinical_root,
                    raw_mngs_root=raw_mngs_root,
                    raw_mngs_patient_id=identity["raw_mngs_patient_id"],
                    expected_case_code=identity["case_code"],
                    expected_specimen_code=identity["specimen_code"],
                    raw_mngs_sha256=identity["raw_mngs_sha256"],
                    clinical_hashes=tuple(identity["clinical_files_sha256"].items()),
                    mapping_sha256=mapping_hash,
                )
            )
    return manifest, patients


def clinical_file(spec: PatientSpec, section: str) -> Path:
    return spec.clinical_dir / f"NGS_patient_{spec.numeric_id}_{section}.json"


def build_patient_payload(spec: PatientSpec) -> tuple[dict[str, Any], dict[str, Any]]:
    if not spec.clinical_dir.is_dir():
        raise FileNotFoundError(f"Clinical folder not found: {spec.clinical_dir}")
    if not spec.raw_mngs_file.is_file():
        raise FileNotFoundError(f"Raw mNGS file not found: {spec.raw_mngs_file}")

    raw_mngs = read_frozen_json(spec.raw_mngs_file, spec.raw_mngs_sha256)
    if not isinstance(raw_mngs, list) or len(raw_mngs) != 1 or not isinstance(raw_mngs[0], dict):
        raise ValueError("An explicit single specimen is required; multiple/empty specimens cannot be silently merged")
    if normalized_code(raw_mngs[0].get("specimen_code")) != normalized_code(spec.expected_specimen_code):
        raise ValueError(f"Raw specimen identity mismatch for {spec.patient_id}")
    index_date = find_index_date(raw_mngs)
    raw_sections: dict[str, Any] = {"all_RK_NTC_microbes": raw_mngs}
    missing_sections: list[str] = []
    empty_sections: list[str] = []

    if raw_mngs in ([], {}, None):
        empty_sections.append("all_RK_NTC_microbes")

    clinical_hashes = dict(spec.clinical_hashes)
    for section in CLINICAL_SECTIONS:
        path = clinical_file(spec, section)
        expected_hash = clinical_hashes[section]
        if path.exists() != (expected_hash is not None):
            raise ValueError(f"Clinical file availability changed since mapping freeze: {path.name}")
        if not path.is_file():
            missing_sections.append(section)
            raw_sections[section] = []
            continue
        content = read_frozen_json(path, expected_hash)
        raw_sections[section] = content
        if content in ([], {}, None):
            empty_sections.append(section)

    stats = SanitizeStats()
    specimen_map: dict[str, str] = {}
    sanitized_sections = sanitize_value(
        raw_sections,
        index_date=index_date,
        specimen_map=specimen_map,
        stats=stats,
    )
    payload = {
        "patient_id": spec.patient_id,
        "data_notes": {
            "task": "Identify the clinical pathogen(s) responsible for the current pulmonary/lower-respiratory infection episode.",
            "empty_array_semantics": "No available record; not a negative test.",
            "date_semantics": "D+0 is the earliest raw mNGS collection date; all available dates are relative days.",
            "mngs_semantics": "Unfiltered all_RK_NTC_microbes; no upstream mNGS_grouped filtering is included.",
        },
        "clinical_record": sanitized_sections,
    }
    date_hits = exact_date_hits(payload)
    audit = {
        "patient_id": spec.patient_id,
        "cohort": spec.cohort,
        "clinical_folder_exists": spec.clinical_dir.is_dir(),
        "raw_mngs_file_exists": spec.raw_mngs_file.is_file(),
        "raw_mngs_nonempty": raw_mngs not in ([], {}, None),
        "index_date_available": index_date is not None,
        "missing_sections": missing_sections,
        "empty_sections": sorted(empty_sections),
        "dropped_keys": dict(sorted(stats.dropped_keys.items())),
        "pseudonymized_specimen_codes": stats.pseudonymized_specimen_codes,
        "converted_dates": stats.converted_dates,
        "unparseable_dates_removed": stats.unparseable_dates_removed,
        "free_text_dates_converted": stats.free_text_dates_converted,
        "remaining_exact_date_hits": date_hits,
        "payload_bytes": len(stable_json(payload).encode("utf-8")),
        "payload_sha256": sha256_text(stable_json(payload)),
        # Local audit only: identity metadata is NEVER inserted into the API payload.
        "input_identity": {
            "clinical_patient_id": spec.patient_id,
            "case_code": spec.expected_case_code,
            "raw_mngs_patient_id": spec.raw_mngs_patient_id,
            "raw_mngs_sha256": spec.raw_mngs_sha256,
            "clinical_files_sha256": clinical_hashes,
            "mapping_sha256": spec.mapping_sha256,
            "specimen_verified": True,
        },
    }
    return payload, audit


def select_patients(all_patients: list[PatientSpec], requested: Iterable[str] | None) -> list[PatientSpec]:
    if not requested:
        return all_patients
    normalized = {value.upper() if value.upper().startswith("P") else f"P{value}" for value in requested}
    available = {spec.patient_id for spec in all_patients}
    unknown = sorted(normalized - available)
    if unknown:
        raise ValueError(f"Patients not in frozen cohort: {', '.join(unknown)}")
    return [spec for spec in all_patients if spec.patient_id in normalized]


def build_audit(manifest_path: Path, requested: Iterable[str] | None = None) -> dict[str, Any]:
    manifest, all_patients = load_manifest(manifest_path)
    patients = select_patients(all_patients, requested)
    records: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    for spec in patients:
        try:
            _, audit = build_patient_payload(spec)
            records.append(audit)
        except Exception as exc:  # audit must report every patient rather than stop early
            failures.append({"patient_id": spec.patient_id, "error": str(exc)})

    blocked = [
        item["patient_id"]
        for item in records
        if not item["raw_mngs_nonempty"]
        or not item["index_date_available"]
        or item["remaining_exact_date_hits"]
    ]
    return {
        "schema_version": "direct_raw_benchmark.input_audit.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "manifest_sha256": sha256_text(stable_json(manifest)),
        "patient_count": len(patients),
        "success_count": len(records),
        "failure_count": len(failures),
        "blocked_patient_ids": blocked,
        "records": records,
        "failures": failures,
    }


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def load_env_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Environment file not found: {path}")
    for raw_line in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if key.startswith("$env:"):
            key = key[5:]
        value = value.strip().strip('"').strip("'")
        if key and key not in os.environ:
            os.environ[key] = value


def usage_to_dict(usage: Any) -> dict[str, Any] | None:
    if usage is None:
        return None
    if hasattr(usage, "model_dump"):
        return usage.model_dump()
    if isinstance(usage, dict):
        return usage
    return {"value": str(usage)}


def validate_prediction(patient_id: str, prediction: dict[str, Any]) -> None:
    if prediction.get("patient_id") != patient_id:
        raise ValueError(
            f"Response patient_id mismatch: expected {patient_id}, got {prediction.get('patient_id')}"
        )
    predicted = prediction.get("predicted_pathogens")
    if not isinstance(predicted, list):
        raise ValueError("predicted_pathogens is not a list")
    no_supported = prediction.get("no_supported_pathogen")
    if bool(predicted) == bool(no_supported):
        raise ValueError(
            "no_supported_pathogen must be true exactly when predicted_pathogens is empty"
        )


def validate_cached_output(record: dict[str, Any], spec: PatientSpec, audit: dict[str, Any],
                           prompt_hash: str, model: str, effort: str, max_tokens: int = 4000) -> None:
    expected = {"patient_id": spec.patient_id, "cohort": spec.cohort,
                "payload_sha256": audit["payload_sha256"], "prompt_sha256": prompt_hash,
                "model_requested": model, "reasoning_effort": effort}
    for field, value in expected.items():
        if record.get(field) != value:
            raise ValueError(f"Cannot reuse stale/mismatched output for {spec.patient_id}: {field}")
    if record.get("max_output_tokens", 4000) != max_tokens:
        raise ValueError("Cached output used a different output-token cap")
    validate_prediction(spec.patient_id, record["prediction"])


def call_openai(
    *,
    client: Any,
    patient_id: str,
    payload: dict[str, Any],
    prompt: str,
    model: str,
    reasoning_effort: str,
    max_output_tokens: int,
    retries: int,
) -> Any:
    request_input = (
        "以下是單一病人的 Direct-Raw 病例資料。請勿使用資料外的病例資訊。\n\n"
        + json.dumps(payload, ensure_ascii=False, indent=2)
    )
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            return client.responses.create(
                model=model,
                instructions=prompt,
                input=request_input,
                reasoning={"effort": reasoning_effort},
                text={
                    "format": {
                        "type": "json_schema",
                        "name": "direct_raw_clinical_pathogen_prediction",
                        "strict": True,
                        "schema": OUTPUT_SCHEMA,
                    }
                },
                max_output_tokens=max_output_tokens,
                store=False,
            )
        except Exception as exc:
            last_error = exc
            if attempt >= retries:
                break
            time.sleep(min(2**attempt, 8))
    raise RuntimeError(f"OpenAI request failed for {patient_id}: {last_error}") from last_error


def run_patients(args: argparse.Namespace) -> int:
    if args.env_file:
        load_env_file(Path(args.env_file).resolve())
    if not os.getenv("OPENAI_API_KEY"):
        raise RuntimeError(
            "OPENAI_API_KEY is not set. Pass --env-file or set the environment variable."
        )

    try:
        from openai import OpenAI
    except ImportError as exc:
        raise RuntimeError("Install dependencies with: pip install -r requirements.txt") from exc

    manifest, all_patients = load_manifest(Path(args.manifest).resolve())
    patients = select_patients(all_patients, args.patients)
    prompt_path = Path(args.prompt).resolve()
    prompt = prompt_path.read_text(encoding="utf-8-sig").strip()
    prompt_sha256 = sha256_text(prompt)
    model = args.model or manifest["model"]
    reasoning_effort = args.reasoning_effort or manifest["reasoning_effort"]
    run_name = args.run_name or datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
    output_dir = Path(args.output_root).resolve() / run_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate EVERY input and existing output before the first paid request.
    prepared = []
    for spec in patients:
        payload, audit = build_patient_payload(spec)
        if audit["remaining_exact_date_hits"] or not audit["index_date_available"]:
            raise ValueError(f"Input sanitization/index date failed for {spec.patient_id}")
        output_path = output_dir / f"{spec.patient_id}_direct_raw.json"
        if output_path.exists() and not args.overwrite:
            validate_cached_output(load_json(output_path), spec, audit, prompt_sha256, model, reasoning_effort, args.max_output_tokens)
        prepared.append((spec, payload, audit))
    write_json(output_dir / "input_audit_v2.json", {
        "schema_version": "direct_raw_benchmark.run_inputs.v2",
        "manifest_sha256": sha256_file(Path(args.manifest).resolve()),
        "mapping_sha256": patients[0].mapping_sha256 if patients else None,
        "records": [audit for _, _, audit in prepared],
    })
    write_json(output_dir / "run_configuration.json", {
        "model": model, "reasoning_effort": reasoning_effort,
        "max_output_tokens": args.max_output_tokens, "store": False, "tools": [],
        "prompt_sha256": prompt_sha256, "runner_sha256": sha256_file(Path(__file__)),
        "output_schema_sha256": sha256_text(stable_json(OUTPUT_SCHEMA)),
        "mapping_sha256": patients[0].mapping_sha256 if patients else None,
    })
    client = OpenAI(max_retries=0, timeout=180.0)
    successes: list[str] = []
    reused: list[str] = []
    failures: list[dict[str, str]] = []
    for spec, payload, audit in prepared:
        output_path = output_dir / f"{spec.patient_id}_direct_raw.json"
        if output_path.exists() and not args.overwrite:
            successes.append(spec.patient_id)
            reused.append(spec.patient_id)
            print(f"REUSE {spec.patient_id}: matching payload/prompt/model verified", flush=True)
            continue
        try:
            if audit["remaining_exact_date_hits"]:
                raise ValueError(f"Exact dates remain after sanitization: {audit['remaining_exact_date_hits']}")
            if not audit["raw_mngs_nonempty"] or not audit["index_date_available"]:
                raise ValueError("Raw mNGS input is empty or lacks an index collection date")
            response = call_openai(
                client=client,
                patient_id=spec.patient_id,
                payload=payload,
                prompt=prompt,
                model=model,
                reasoning_effort=reasoning_effort,
                max_output_tokens=args.max_output_tokens,
                retries=args.retries,
            )
            if getattr(response, "status", "completed") != "completed":
                write_json(output_dir / f"{spec.patient_id}_incomplete_response.json", {
                    "response_id": getattr(response, "id", None), "status": getattr(response, "status", None),
                    "usage": usage_to_dict(getattr(response, "usage", None)),
                    "output_text": getattr(response, "output_text", None),
                })
                raise ValueError(f"Response is not completed for {spec.patient_id}; saved for review")
            prediction = json.loads(response.output_text)
            validate_prediction(spec.patient_id, prediction)
            record = {
                "schema_version": "direct_raw_benchmark.output.v2",
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "patient_id": spec.patient_id,
                "cohort": spec.cohort,
                "model_requested": model,
                "model_returned": getattr(response, "model", None),
                "reasoning_effort": reasoning_effort,
                "prompt_sha256": prompt_sha256,
                "payload_sha256": audit["payload_sha256"],
                "input_identity": audit["input_identity"],
                "max_output_tokens": args.max_output_tokens,
                "response_id": getattr(response, "id", None),
                "usage": usage_to_dict(getattr(response, "usage", None)),
                "prediction": prediction,
            }
            write_json(output_path, record)
            successes.append(spec.patient_id)
            print(f"OK {spec.patient_id}", flush=True)
        except Exception as exc:
            failures.append({"patient_id": spec.patient_id, "error": str(exc)})
            print(f"ERROR {spec.patient_id}: {exc}", file=sys.stderr, flush=True)
            if args.fail_fast:
                break

    summary = {
        "schema_version": "direct_raw_benchmark.run_summary.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_name": run_name,
        "model": model,
        "reasoning_effort": reasoning_effort,
        "prompt_sha256": prompt_sha256,
        "requested_patient_ids": [spec.patient_id for spec in patients],
        "success_patient_ids": successes,
        "reused_patient_ids": reused,
        "new_success_patient_ids": [p for p in successes if p not in reused],
        "input_mapping_sha256": patients[0].mapping_sha256 if patients else None,
        "failures": failures,
    }
    write_json(output_dir / "run_summary.json", summary)
    return 1 if failures else 0


def audit_command(args: argparse.Namespace) -> int:
    audit = build_audit(Path(args.manifest).resolve(), args.patients)
    output = Path(args.output).resolve()
    write_json(output, audit)
    print(f"Audit written: {output}")
    print(
        f"patients={audit['patient_count']} success={audit['success_count']} "
        f"failures={audit['failure_count']} blocked={len(audit['blocked_patient_ids'])}"
    )
    if audit["blocked_patient_ids"]:
        print("blocked_patient_ids=" + ",".join(audit["blocked_patient_ids"]))
    return 1 if audit["failure_count"] or audit["blocked_patient_ids"] else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the frozen 41-patient Direct-Raw GPT benchmark one patient folder at a time."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    audit_parser = subparsers.add_parser("audit", help="Build sanitized payloads in memory and audit them without API calls.")
    audit_parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    audit_parser.add_argument("--patients", nargs="*")
    audit_parser.add_argument("--output", default=str(DEFAULT_AUDIT))
    audit_parser.set_defaults(func=audit_command)

    run_parser = subparsers.add_parser("run", help="Run independent Responses API calls for selected patients.")
    run_parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    run_parser.add_argument("--prompt", default=str(DEFAULT_PROMPT))
    run_parser.add_argument("--patients", nargs="*")
    run_parser.add_argument("--model")
    run_parser.add_argument("--reasoning-effort", choices=["none", "low", "medium", "high", "xhigh", "max"])
    run_parser.add_argument("--max-output-tokens", type=int, default=4000)
    run_parser.add_argument("--retries", type=int, default=2)
    run_parser.add_argument("--env-file")
    run_parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    run_parser.add_argument("--run-name")
    run_parser.add_argument("--overwrite", action="store_true")
    run_parser.add_argument("--fail-fast", action="store_true")
    run_parser.set_defaults(func=run_patients)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        return args.func(args)
    except Exception as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from tools.time_window import cutoff_entries_from_manifest, filter_entries_by_mngs_window


DEFAULT_PATIENT_EXPORT_ROOT = Path("outputs") / "excel_patient_exports_demo_v2"
DEFAULT_OUTPUT_NAME = "ppt_link_summary.json"

DATE_TOKEN_PATTERN = re.compile(r"(20\d{2})(\d{2})(\d{2})")
CASE_DIGITS_PATTERN = re.compile(r"(\d+)$")
RADIOLOGY_CHUNK_PATTERN = re.compile(
    r"(?P<date>20\d{6})\s*Imaging findings\s*:?\s*(?P<body>.*?)(?=(?:20\d{6}\s*Imaging findings)|\Z)",
    re.IGNORECASE | re.DOTALL,
)
RADIOLOGY_DATE_PATTERN = re.compile(r"(20\d{2})[/-](\d{1,2})[/-](\d{1,2})")
RADIOLOGY_EXAM_TITLE_PATTERN = re.compile(
    r"\b(CHEST|CXR|X-?RAY|CT|MRI|SONO|ULTRASOUND|KUB|PA VIEW|AP VIEW|PORTABLE)\b",
    re.IGNORECASE,
)


def _read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _read_list_json(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    payload = _read_json(path)
    if not isinstance(payload, list):
        return []
    return [entry for entry in payload if isinstance(entry, dict)]


def _load_mngs_window_entries(patient_dir: Path, patient_key: str) -> list[dict[str, Any]]:
    mngs_entries = _read_list_json(patient_dir / f"{patient_key}_mNGS_grouped.json")
    if mngs_entries:
        return mngs_entries
    manifest_path = patient_dir / f"{patient_key}_patient_manifest.json"
    if not manifest_path.exists():
        return []
    manifest = _read_json(manifest_path)
    return cutoff_entries_from_manifest(manifest)


def _text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def _normalize_hospital_id(value: Any) -> str:
    text = _text(value).replace("\u3000", " ")
    return re.sub(r"\s+", "", text)


def _normalize_case_code(value: Any) -> str:
    text = _text(value).upper().replace("-", "")
    match = CASE_DIGITS_PATTERN.search(text)
    if not match:
        return text
    digits = match.group(1).zfill(3)
    if text.startswith("TS"):
        return f"T{digits}"
    prefix = re.sub(r"\d+$", "", text)
    return f"{prefix}{digits}"


def _fallback_hospital_id_from_source_file(value: Any) -> str:
    source_file = Path(_text(value)).stem
    return _normalize_hospital_id(source_file)


def _normalize_compact_date(value: str) -> str:
    match = DATE_TOKEN_PATTERN.fullmatch(value)
    if not match:
        return value
    year, month, day = match.groups()
    return f"{year}-{month}-{day}"


def _normalize_slashed_date(value: str) -> str:
    match = RADIOLOGY_DATE_PATTERN.search(value)
    if not match:
        return ""
    year, month, day = match.groups()
    return f"{year}-{int(month):02d}-{int(day):02d}"


def _load_patient_index(patient_export_root: Path) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    index_path = patient_export_root / "patient_index.json"
    payload = _read_json(index_path)
    by_hospital_id: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in payload:
        hospital_id = _normalize_hospital_id(row.get("hospital_id"))
        if hospital_id:
            by_hospital_id[hospital_id].append(row)
    return payload, by_hospital_id


def _resolve_patients_for_slide(
    slide_payload: Mapping[str, Any],
    patient_index: Sequence[Mapping[str, Any]],
    by_hospital_id: Mapping[str, list[Mapping[str, Any]]],
) -> tuple[list[str], str | None]:
    patient_id = _normalize_hospital_id(slide_payload.get("patient_id"))
    case_id = _normalize_case_code(slide_payload.get("case_id"))
    fallback_patient_id = _fallback_hospital_id_from_source_file(slide_payload.get("source_file"))

    hospital_id_candidates: list[str] = []
    for candidate in (patient_id, fallback_patient_id):
        if candidate and candidate not in hospital_id_candidates:
            hospital_id_candidates.append(candidate)

    for hospital_id in hospital_id_candidates:
        matches = list(by_hospital_id.get(hospital_id, []))
        if len(matches) == 1:
            return [_text(matches[0].get("patient_key"))], None
        if len(matches) > 1:
            for row in matches:
                if _normalize_case_code(row.get("case_code")) == case_id:
                    return [_text(row.get("patient_key"))], None
            # When one hospital_id maps to multiple specimen-level patient keys,
            # fan out the PPT to every matching key instead of dropping it.
            fanout_keys = [_text(row.get("patient_key")) for row in matches if _text(row.get("patient_key"))]
            if fanout_keys:
                return fanout_keys, f"fanout_hospital_id:{hospital_id}"
            return [], f"ambiguous_hospital_id:{hospital_id}"

    if case_id:
        exact_case_matches = [row for row in patient_index if _normalize_case_code(row.get("case_code")) == case_id]
        if len(exact_case_matches) == 1:
            return [_text(exact_case_matches[0].get("patient_key"))], None
        if len(exact_case_matches) > 1:
            return [], f"ambiguous_case_code:{case_id}"

    return [], "no_matching_patient"


def _clean_finding_line(line: str) -> str:
    cleaned = line.strip()
    while True:
        updated = re.sub(r"^(?:[>\.\-\u2022]+\s*)+", "", cleaned)
        if updated == cleaned:
            break
        cleaned = updated.strip()
    cleaned = re.sub(r"\s+", " ", cleaned)
    return cleaned.strip()


def _looks_like_radiology_admin_line(line: str) -> bool:
    lowered = line.lower().strip()
    if not lowered:
        return True
    prefixes = (
        "date of examination",
        "date of rcp",
        "reported by",
        "on:",
        "license#:",
    )
    if lowered.startswith(prefixes):
        return True
    if "at:report" in lowered:
        return True
    if lowered.startswith("date:"):
        return True
    return False


def _looks_like_radiology_exam_title(line: str) -> bool:
    return bool(RADIOLOGY_EXAM_TITLE_PATTERN.search(line.strip()))


def _extract_radiology_report_fallback(slide_payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    line_candidates = (
        slide_payload.get("extracted_findings_text")
        or slide_payload.get("analysis_text_blocks")
        or slide_payload.get("image_derived_text_blocks")
        or []
    )
    flattened_lines: list[str] = []
    for block in line_candidates:
        text = _text(block)
        if not text:
            continue
        flattened_lines.extend(text.splitlines())

    if not flattened_lines:
        return []

    exam_title = ""
    findings: list[str] = []
    saw_admin_hint = False
    saw_bullet_finding = False
    for raw_line in flattened_lines:
        stripped = raw_line.strip()
        if not stripped:
            continue
        if _looks_like_radiology_admin_line(stripped):
            saw_admin_hint = True
            continue
        if not exam_title and not stripped.startswith(">") and not re.match(r"^[\.\-\u2022]\s*", stripped):
            exam_title = _clean_finding_line(stripped)
            continue
        cleaned = _clean_finding_line(stripped)
        if not cleaned:
            continue
        is_bullet = stripped.startswith(">") or bool(re.match(r"^[\.\-\u2022]\s*", stripped))
        if is_bullet:
            saw_bullet_finding = True
        if is_bullet or not findings:
            findings.append(cleaned)
        else:
            findings[-1] = f"{findings[-1]} {cleaned}".strip()

    if not findings:
        return []
    if not saw_bullet_finding:
        return []
    if not (saw_admin_hint or _looks_like_radiology_exam_title(exam_title)):
        return []

    date_candidates = [
        _text(slide_payload.get("date")),
        *[_text(value) for value in slide_payload.get("detected_dates") or []],
    ]
    normalized_date = ""
    for candidate in date_candidates:
        if not candidate:
            continue
        normalized_date = _normalize_slashed_date(candidate) or candidate
        if normalized_date:
            break

    entry = {
        "test": "Radiology / Special Examination",
        "exam_type": exam_title or "Imaging findings from PPT",
        "collected_time": normalized_date,
        "reported_time": normalized_date,
        "findings": findings,
        "source": {
            "source_file": _text(slide_payload.get("source_file")),
            "slide_index": slide_payload.get("slide_index"),
        },
    }
    return [entry]


def extract_radiology_entries(slide_payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    modalities = set(slide_payload.get("candidate_modalities") or [])

    entries: list[dict[str, Any]] = []
    source_blocks = slide_payload.get("extracted_findings_text") or slide_payload.get("detected_text_blocks") or []

    for block in source_blocks:
        text = _text(block)
        if not text or "Imaging findings" not in text:
            continue
        for match in RADIOLOGY_CHUNK_PATTERN.finditer(text):
            date = _normalize_compact_date(match.group("date"))
            body = match.group("body")
            findings: list[str] = []
            skipping_note = False
            for raw_line in body.splitlines():
                stripped = raw_line.strip()
                if not stripped:
                    continue
                if stripped.startswith("("):
                    skipping_note = True
                    continue
                if skipping_note:
                    continue
                line = _clean_finding_line(raw_line)
                if not line or line.lower().startswith("imaging findings"):
                    continue
                is_bullet = bool(re.match(r"^[\.\-\u2022]\s*", stripped))
                if is_bullet or not findings:
                    findings.append(line)
                else:
                    findings[-1] = f"{findings[-1]} {line}".strip()
            if not findings:
                continue
            entry = {
                "test": "Radiology / Special Examination",
                "exam_type": "Imaging findings from PPT",
                "collected_time": date,
                "reported_time": date,
                "findings": findings,
                "source": {
                    "source_file": _text(slide_payload.get("source_file")),
                    "slide_index": slide_payload.get("slide_index"),
                },
            }
            entries.append(entry)

    deduped: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for entry in entries:
        key = (entry["collected_time"], " | ".join(entry["findings"]))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(entry)
    if deduped:
        return deduped
    return _extract_radiology_report_fallback(slide_payload)


def _merge_image_entries(existing_entries: list[dict[str, Any]], new_entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    merged = list(existing_entries)
    seen = {
        (
            _text(entry.get("collected_time")),
            " | ".join(entry.get("findings", [])),
        )
        for entry in merged
        if isinstance(entry, dict)
    }
    for entry in new_entries:
        key = (_text(entry.get("collected_time")), " | ".join(entry.get("findings", [])))
        if key in seen:
            continue
        seen.add(key)
        merged.append(entry)
    return merged


def _drop_entries_from_same_slide(
    existing_entries: list[dict[str, Any]],
    *,
    source_file: str,
    slide_index: Any,
) -> list[dict[str, Any]]:
    filtered: list[dict[str, Any]] = []
    for entry in existing_entries:
        if not isinstance(entry, dict):
            continue
        source = entry.get("source")
        if isinstance(source, dict):
            same_source = _text(source.get("source_file")) == source_file
            same_slide = source.get("slide_index") == slide_index
            if same_source and same_slide:
                continue
        filtered.append(entry)
    return filtered


def _scan_slide_jsons(ppt_root: Path) -> list[Path]:
    return sorted(ppt_root.glob("**/slides/slide_*/parsed.json"))


def link_ppt_to_patient_exports(
    ppt_root: Path,
    patient_export_root: Path,
) -> dict[str, Any]:
    patient_index, by_hospital_id = _load_patient_index(patient_export_root)
    slide_paths = _scan_slide_jsons(ppt_root)

    linked_by_patient: dict[str, list[dict[str, Any]]] = defaultdict(list)
    unmatched: list[dict[str, Any]] = []

    for slide_path in slide_paths:
        slide_payload = _read_json(slide_path)
        patient_keys, reason = _resolve_patients_for_slide(slide_payload, patient_index, by_hospital_id)
        bundle_entry = {
            "source_file": _text(slide_payload.get("source_file")),
            "slide_index": slide_payload.get("slide_index"),
            "case_id": _text(slide_payload.get("case_id")),
            "patient_id": _normalize_hospital_id(slide_payload.get("patient_id")),
            "date": _text(slide_payload.get("date")),
            "candidate_modalities": list(slide_payload.get("candidate_modalities") or []),
            "review_reasons": list(slide_payload.get("review_reasons") or []),
            "parsed_json_path": str(slide_path.resolve()),
        }
        if not patient_keys:
            unmatched.append({**bundle_entry, "reason": reason})
            continue
        for patient_key in patient_keys:
            linked_by_patient[patient_key].append(bundle_entry)

            patient_dir = patient_export_root / patient_key
            bundle_path = patient_dir / f"{patient_key}_ppt_slide_bundle.json"
            existing_bundle = _read_json(bundle_path) if bundle_path.exists() else []
            existing_bundle = [entry for entry in existing_bundle if isinstance(entry, dict)]
            deduped_bundle = []
            seen_bundle_keys: set[tuple[str, Any]] = set()
            for entry in existing_bundle + [bundle_entry]:
                bundle_key = (_text(entry.get("source_file")), entry.get("slide_index"))
                if bundle_key in seen_bundle_keys:
                    continue
                seen_bundle_keys.add(bundle_key)
                deduped_bundle.append(entry)
            _write_json(bundle_path, deduped_bundle)

            radiology_entries = extract_radiology_entries(slide_payload)
            if radiology_entries:
                image_path = patient_dir / f"{patient_key}_image.json"
                mngs_entries = _load_mngs_window_entries(patient_dir, patient_key)
                radiology_entries = filter_entries_by_mngs_window(radiology_entries, mngs_entries)
                existing_image = _read_list_json(image_path)
                existing_image = _drop_entries_from_same_slide(
                    existing_image,
                    source_file=_text(slide_payload.get("source_file")),
                    slide_index=slide_payload.get("slide_index"),
                )
                merged_image = _merge_image_entries(existing_image, radiology_entries)
                _write_json(image_path, merged_image)

    summary = {
        "ppt_root": str(ppt_root.resolve()),
        "patient_export_root": str(patient_export_root.resolve()),
        "slides_scanned": len(slide_paths),
        "slides_linked": sum(len(entries) for entries in linked_by_patient.values()),
        "patients_with_linked_slides": len(linked_by_patient),
        "unmatched_slides": unmatched,
    }
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Link PPT slide-level parsed outputs to the exported Excel patient folders.",
    )
    parser.add_argument("ppt_root", type=Path, help="Path to the PPT parser output root.")
    parser.add_argument(
        "--patient-export-root",
        type=Path,
        default=DEFAULT_PATIENT_EXPORT_ROOT,
        help=f"Path to the per-patient Excel export root (default: {DEFAULT_PATIENT_EXPORT_ROOT}).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    ppt_root = args.ppt_root.expanduser().resolve()
    patient_export_root = args.patient_export_root.expanduser().resolve()

    summary = link_ppt_to_patient_exports(ppt_root, patient_export_root)
    summary_path = patient_export_root / DEFAULT_OUTPUT_NAME
    _write_json(summary_path, summary)
    print(f"Linked {summary['slides_linked']} slides across {summary['patients_with_linked_slides']} patients.")
    print(f"Wrote summary: {summary_path}")


if __name__ == "__main__":
    main()

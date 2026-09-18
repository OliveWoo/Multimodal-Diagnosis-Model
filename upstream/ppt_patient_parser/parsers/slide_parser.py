from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from ..modality_classifier import classify_modalities

CASE_ID_PATTERN = re.compile(r"\b[A-Z]{1,4}[-_ ]?\d{1,4}\b")
PATIENT_ID_PATTERN = re.compile(r"\b\d{6,10}\b")
DATE_PATTERN = re.compile(r"\b((?:20\d{2}[-/]\d{1,2}[-/]\d{1,2})|(?:20\d{6}))\b")


def _dedupe_preserve_order(values: list[str]) -> list[str]:
    deduped: list[str] = []
    seen: set[str] = set()
    for value in values:
        normalized = value.strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(normalized)
    return deduped


def _normalize_date(raw_value: str) -> str | None:
    digits_only = re.sub(r"\D", "", raw_value)
    if len(digits_only) == 8 and digits_only.startswith("20"):
        return f"{digits_only[:4]}-{digits_only[4:6]}-{digits_only[6:8]}"

    match = re.fullmatch(r"(20\d{2})[-/](\d{1,2})[-/](\d{1,2})", raw_value)
    if not match:
        return None
    year, month, day = match.groups()
    return f"{year}-{int(month):02d}-{int(day):02d}"


def _guess_case_id(text_blocks: list[str]) -> str | None:
    for block in text_blocks:
        match = CASE_ID_PATTERN.search(block)
        if match:
            return match.group(0).replace(" ", "").replace("_", "-")
    return None


def _guess_patient_id(text_blocks: list[str]) -> str | None:
    for block in text_blocks:
        match = PATIENT_ID_PATTERN.search(block)
        if match:
            return match.group(0)
    return None


def _guess_date(text_blocks: list[str]) -> str | None:
    for block in text_blocks:
        match = DATE_PATTERN.search(block)
        if match:
            return _normalize_date(match.group(1))
    return None


def _extract_dates(text_blocks: list[str]) -> list[str]:
    collected: list[str] = []
    for block in text_blocks:
        for match in DATE_PATTERN.findall(block):
            normalized = _normalize_date(match)
            if normalized:
                collected.append(normalized)
    return _dedupe_preserve_order(collected)


def _geometry(payload: dict[str, Any]) -> dict[str, int]:
    geometry = payload.get("geometry") or payload.get("display_geometry") or {}
    left = int(geometry.get("left", 0) or 0)
    top = int(geometry.get("top", 0) or 0)
    width = int(geometry.get("width", 0) or 0)
    height = int(geometry.get("height", 0) or 0)
    return {
        "left": left,
        "top": top,
        "width": width,
        "height": height,
        "right": int(geometry.get("right", left + width) or (left + width)),
        "bottom": int(geometry.get("bottom", top + height) or (top + height)),
        "center_x": int(geometry.get("center_x", left + (width // 2)) or (left + (width // 2))),
        "center_y": int(geometry.get("center_y", top + (height // 2)) or (top + (height // 2))),
        "shape_order": int(geometry.get("shape_order", 0) or 0),
    }


def _vertical_gap(text_geometry: dict[str, int], image_geometry: dict[str, int]) -> int:
    if text_geometry["bottom"] < image_geometry["top"]:
        return image_geometry["top"] - text_geometry["bottom"]
    if image_geometry["bottom"] < text_geometry["top"]:
        return text_geometry["top"] - image_geometry["bottom"]
    return 0


def _horizontal_overlap_ratio(text_geometry: dict[str, int], image_geometry: dict[str, int]) -> float:
    overlap = max(
        0,
        min(text_geometry["right"], image_geometry["right"]) - max(text_geometry["left"], image_geometry["left"]),
    )
    baseline = min(
        max(1, text_geometry["width"]),
        max(1, image_geometry["width"]),
    )
    return overlap / baseline


def _nearby_text_records(image: dict[str, Any], text_records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    image_geometry = _geometry(image)
    if image_geometry["width"] <= 0 and image_geometry["height"] <= 0:
        return []

    image_area = max(1, image_geometry["width"] * image_geometry["height"])
    caption_like_candidates: list[tuple[tuple[int, int, int], dict[str, Any]]] = []
    general_candidates: list[tuple[tuple[int, int, int], dict[str, Any]]] = []
    for record in text_records:
        text_geometry = _geometry(record)
        if text_geometry["width"] <= 0 and text_geometry["height"] <= 0:
            continue

        center_gap = abs(text_geometry["center_x"] - image_geometry["center_x"])
        vertical_gap = _vertical_gap(text_geometry, image_geometry)
        overlap_ratio = _horizontal_overlap_ratio(text_geometry, image_geometry)
        aligned_horizontally = overlap_ratio >= 0.1 or center_gap <= max(
            image_geometry["width"],
            text_geometry["width"],
        )
        max_vertical_gap = max(
            image_geometry["height"] // 2,
            text_geometry["height"] * 3,
            1,
        )
        if aligned_horizontally and vertical_gap <= max_vertical_gap:
            order_gap = abs(text_geometry["shape_order"] - image_geometry["shape_order"])
            text_area = max(1, text_geometry["width"] * text_geometry["height"])
            sort_key = (vertical_gap, center_gap, order_gap)
            general_candidates.append((sort_key, record))
            if (
                text_area <= image_area // 2
                or text_geometry["height"] <= max(1, image_geometry["height"] // 3)
            ):
                caption_like_candidates.append((sort_key, record))

    selected_candidates = caption_like_candidates or general_candidates
    selected_candidates.sort(key=lambda item: item[0])
    return [record for _, record in selected_candidates[:3]]


def _enrich_image_records(
    image_records: list[dict[str, Any]],
    text_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    enriched_records: list[dict[str, Any]] = []
    for image in image_records:
        enriched = dict(image)
        nearby_records = _nearby_text_records(enriched, text_records)
        enriched["nearby_text_blocks"] = [record["text"] for record in nearby_records]
        enriched["nearby_text_records"] = nearby_records
        enriched_records.append(enriched)
    return enriched_records


def _image_vision_status(image: dict[str, Any]) -> str:
    vision = image.get("vision") or {}
    status = vision.get("status")
    if isinstance(status, str) and status:
        return status
    if vision.get("raw_text") or vision.get("findings"):
        return "success"
    return "not_requested"


def _extract_image_text(image_records: list[dict[str, Any]]) -> tuple[list[str], list[str], dict[str, int]]:
    image_text_blocks: list[str] = []
    image_findings: list[str] = []
    status_summary = {
        "success": 0,
        "empty": 0,
        "not_requested": 0,
        "error": 0,
    }
    for image in image_records:
        vision = image.get("vision") or {}
        status = _image_vision_status(image)
        if status not in status_summary:
            status_summary[status] = 0
        status_summary[status] += 1
        raw_text = vision.get("raw_text") or []
        findings = vision.get("findings") or []
        if isinstance(raw_text, list):
            image_text_blocks.extend(str(item) for item in raw_text if str(item).strip())
        if isinstance(findings, list):
            image_findings.extend(str(item) for item in findings if str(item).strip())
    return (
        _dedupe_preserve_order(image_text_blocks),
        _dedupe_preserve_order(image_findings),
        status_summary,
    )


def _image_has_extracted_content(image: dict[str, Any]) -> bool:
    vision = image.get("vision") or {}
    raw_text = vision.get("raw_text") or []
    findings = vision.get("findings") or []
    nearby_text = image.get("nearby_text_blocks") or []
    return bool(raw_text or findings or nearby_text)


def parse_slide(
    *,
    source_file: Path,
    slide_index: int,
    text_blocks: list[str],
    text_records: list[dict[str, Any]] | None = None,
    image_records: list[dict],
) -> dict:
    text_records = text_records or []
    enriched_images = _enrich_image_records(image_records, text_records)
    image_derived_text_blocks, image_findings, image_extraction_status_summary = _extract_image_text(enriched_images)
    analysis_text_blocks = _dedupe_preserve_order(text_blocks + image_derived_text_blocks + image_findings)

    candidate_modalities, content_type_guess = classify_modalities(analysis_text_blocks)
    case_id = _guess_case_id(analysis_text_blocks)
    patient_id = _guess_patient_id(analysis_text_blocks)
    date = _guess_date(analysis_text_blocks)
    detected_dates = _extract_dates(analysis_text_blocks)

    review_reasons: list[str] = []
    if not text_blocks:
        review_reasons.append("no_slide_text_detected")
    if not enriched_images:
        review_reasons.append("no_embedded_images_detected")
    if any(image["preprocessing"].get("suspicious_distortion") for image in enriched_images):
        review_reasons.append("contains_suspiciously_distorted_images")
    if not candidate_modalities:
        review_reasons.append("no_modality_keywords_detected")
    if case_id is None:
        review_reasons.append("case_id_not_detected")
    if enriched_images and any(_image_vision_status(image) == "error" for image in enriched_images):
        review_reasons.append("image_text_extraction_failed")
    if enriched_images and not any(_image_has_extracted_content(image) for image in enriched_images):
        review_reasons.append("embedded_images_without_extracted_content")

    extracted_findings_text: list[str] = []
    if content_type_guess["contains_radiology_text"]:
        extracted_findings_text.extend(block for block in analysis_text_blocks if len(block) > 20)

    lung_related_hints: list[str] = []
    lowered = "\n".join(analysis_text_blocks).lower()
    for keyword in ("opacity", "ground glass", "ground-glass", "pleural effusion", "consolidation", "pneumonia"):
        if keyword in lowered:
            lung_related_hints.append(keyword)

    severity_hints: list[str] = []
    for keyword in ("shock", "intubation", "icu", "vasopressor", "severe"):
        if keyword in lowered:
            severity_hints.append(keyword)

    return {
        "source_file": str(source_file),
        "slide_index": slide_index,
        "case_id": case_id,
        "patient_id": patient_id,
        "date": date,
        "detected_dates": detected_dates,
        "detected_text_blocks": text_blocks,
        "image_derived_text_blocks": image_derived_text_blocks,
        "image_findings": image_findings,
        "analysis_text_blocks": analysis_text_blocks,
        "images": enriched_images,
        "image_extraction_status_summary": image_extraction_status_summary,
        "candidate_modalities": candidate_modalities,
        "content_type_guess": content_type_guess,
        "extracted_findings_text": _dedupe_preserve_order(extracted_findings_text + image_findings),
        "severity_hints": sorted(set(severity_hints)),
        "lung_related_hints": sorted(set(lung_related_hints)),
        "needs_review": bool(review_reasons),
        "review_reasons": review_reasons,
    }

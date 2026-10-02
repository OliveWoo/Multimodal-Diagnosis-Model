"""Source provenance and anchor extraction for Raw-to-structured conversion."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable


PROVENANCE_SCHEMA_VERSION = "source_provenance_v1"
DATE_PATTERN = re.compile(
    r"(?<!\d)(20\d{2})[-/](\d{1,2})[-/](\d{1,2})(?!\d)|(?<!\d)(20\d{6})(?!\d)"
)
QUANTITY_PATTERN = re.compile(
    r"(?i)(?:[<>]=?|≥|≤)?\s*\d+(?:\.\d+)?"
    r"(?:\s*[x×]\s*10\s*\^?\s*\d+)?\s*"
    r"(?:cfu(?:/ml)?|copies(?:/ml)?|iu(?:/ml)?|pg/ml|ng/ml|mg/dl|g/dl|"
    r"cells?/u?l|index|s/co|%|mmol/l|meq/l)"
)
SPECIMEN_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("balf", re.compile(r"(?i)\b(?:balf|bal|bronchoalveolar\s+lavage)\b")),
    ("sputum", re.compile(r"(?i)\bsputum\b|痰")),
    ("endotracheal_aspirate", re.compile(r"(?i)\b(?:et|endotracheal)\s*(?:aspirate|suction)?\b")),
    ("blood", re.compile(r"(?i)\b(?:blood|serum|plasma)\b|血液|血清|血漿")),
    ("urine", re.compile(r"(?i)\b(?:urine|foley)\b|尿液|導尿")),
    ("stool", re.compile(r"(?i)\b(?:stool|feces|faeces)\b|糞便")),
    ("pleural_fluid", re.compile(r"(?i)\bpleural\s+fluid\b|胸水")),
    ("csf", re.compile(r"(?i)\b(?:csf|cerebrospinal\s+fluid)\b|腦脊髓液")),
    ("tissue", re.compile(r"(?i)\b(?:tissue|biopsy)\b|組織|切片")),
)


def stable_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def payload_sha256(value: Any) -> str:
    return hashlib.sha256(stable_json_bytes(value)).hexdigest()


def iter_nonempty_raw_cells(payload: Any) -> Iterable[dict[str, Any]]:
    if not isinstance(payload, dict):
        return
    rows = payload.get("raw_data")
    if not isinstance(rows, list):
        return
    for row_index, row in enumerate(rows, start=1):
        values = row if isinstance(row, list) else [row]
        for column_index, value in enumerate(values, start=1):
            text = str(value).strip() if value is not None else ""
            if text:
                yield {"row": row_index, "column": column_index, "value": text}


def normalize_date_token(value: str) -> str | None:
    match = DATE_PATTERN.search(value)
    if not match:
        return None
    if match.group(4):
        digits = match.group(4)
        return f"{digits[:4]}-{digits[4:6]}-{digits[6:8]}"
    return f"{match.group(1)}-{int(match.group(2)):02d}-{int(match.group(3)):02d}"


def extract_source_anchors(payload: Any) -> dict[str, list[dict[str, Any]]]:
    anchors: dict[str, list[dict[str, Any]]] = {
        "dates": [], "quantities": [], "specimens": []
    }
    seen: dict[str, set[tuple[Any, ...]]] = {key: set() for key in anchors}
    for cell in iter_nonempty_raw_cells(payload):
        text = cell["value"]
        for match in DATE_PATTERN.finditer(text):
            normalized = normalize_date_token(match.group(0))
            key = (normalized, cell["row"], cell["column"])
            if normalized and key not in seen["dates"]:
                seen["dates"].add(key)
                anchors["dates"].append(
                    {**cell, "raw": match.group(0), "normalized": normalized}
                )
        for match in QUANTITY_PATTERN.finditer(text):
            raw = re.sub(r"\s+", " ", match.group(0)).strip()
            key = (raw.lower(), cell["row"], cell["column"])
            if key not in seen["quantities"]:
                seen["quantities"].add(key)
                anchors["quantities"].append({**cell, "raw": raw})
        for specimen, pattern in SPECIMEN_PATTERNS:
            if pattern.search(text):
                key = (specimen, cell["row"], cell["column"])
                if key not in seen["specimens"]:
                    seen["specimens"].add(key)
                    anchors["specimens"].append({**cell, "canonical": specimen})
    return anchors


def build_source_provenance(payload: Any, source_path: Path | None = None) -> dict[str, Any]:
    cells = list(iter_nonempty_raw_cells(payload))
    rows = payload.get("raw_data") if isinstance(payload, dict) else None
    return {
        "schema_version": PROVENANCE_SCHEMA_VERSION,
        "source_raw_json": str(source_path.resolve()) if source_path else None,
        "source_file": payload.get("source_file") if isinstance(payload, dict) else None,
        "sheet": payload.get("sheet") if isinstance(payload, dict) else None,
        "raw_payload_sha256": payload_sha256(payload),
        "raw_row_count": len(rows) if isinstance(rows, list) else 0,
        "raw_nonempty_cell_count": len(cells),
        "critical_anchor_index": extract_source_anchors(payload),
        "raw_payload_retained_externally": source_path is not None,
        "tier_effect": "none_provenance_only",
    }


def attach_source_provenance(
    parsed: dict[str, Any], raw_payload: Any, source_path: Path | None = None
) -> dict[str, Any]:
    enriched = dict(parsed)
    enriched["_source_provenance"] = build_source_provenance(raw_payload, source_path)
    return enriched

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Sequence

from tools.time_window import cutoff_entries_from_manifest, filter_entries_by_mngs_window


DEFAULT_PATIENT_EXPORT_ROOT = Path("outputs") / "excel_patient_exports"
SUMMARY_NAME = "ppt_lab_import_summary.json"

DATETIME_PATTERN = re.compile(r"^(?:20)?\d{2}[/-]\d{2}[/-]\d{2}(?:\s+\d{2}:\d{2}(?::\d{2})?)?$")
VALUE_PATTERN = re.compile(r"^[<>]?\d+(?:\.\d+)?$")
RANGE_PATTERN = re.compile(r"^(?:[<>]?\d+(?:\.\d+)?)?\s*-\s*[<>]?\d+(?:\.\d+)?$")
UNIT_PATTERN = re.compile(r"^[A-Za-z%/^0-9_.-]+(?:/[A-Za-z%/^0-9_.-]+)?$")
QUOTED_TOKEN_PATTERN = re.compile(r'"((?:\\.|[^"])*)"')
NORMALIZED_SPACE_PATTERN = re.compile(r"\s+")
ROW_TOKEN_PATTERN = re.compile(
    r"^(?P<date>(?:20)?\d{2}-\d{2}-\d{2})\s+(?P<time>\d{2}:\d{2}(?::\d{2})?)\s+(?P<body>.+)$"
)

CBC_ITEM_PREFIXES = {
    "wbc",
    "rbc",
    "hb",
    "hgb",
    "hct",
    "plt",
    "mcv",
    "mch",
    "mchc",
    "rdw-sd",
    "rdw-cv",
    "neutrophil",
    "lymphocyte",
    "monocyte",
    "eosinophil",
    "basophil",
    "band",
}
KNOWN_ALPHA_UNITS = {"mmhg", "fl", "pg", "iu"}
KNOWN_SPECIMEN_PREFIXES = ("blood", "urine", "sputum", "balf", "serum", "plasma", "csf")
ROW_MAJOR_CBC_HEADER_MAP = {
    "WBC": ("WBC", "CBC"),
    "RBC": ("RBC", "CBC"),
    "HGB": ("Hgb", "CBC"),
    "HCT": ("Hct", "CBC"),
    "MCV": ("MCV", "CBC"),
    "MCH": ("MCH", "CBC"),
    "MCHC": ("MCHC", "CBC"),
    "RDW": ("RDW", "CBC"),
    "PLT": ("PLT", "CBC"),
    "BAND": ("Band", "CBC"),
    "SEG": ("Seg", "CBC"),
    "LYM": ("Lymphocyte", "CBC"),
    "MONO": ("Monocyte", "CBC"),
    "EOS": ("Eosinophil", "CBC"),
    "BASO": ("Basophil", "CBC"),
    "INR(PT)": ("INR(PT)", "Other Lab Tests"),
    "PT": ("PT", "Other Lab Tests"),
    "APTT": ("APTT", "Other Lab Tests"),
    "D-DIMER": ("D-dimer", "Other Lab Tests"),
    "FDP": ("FDP", "Other Lab Tests"),
    "FIBRINOGEN": ("Fibrinogen", "Other Lab Tests"),
}
ROW_MAJOR_ABG_WITH_FIO2_HEADERS = (
    "Blood gas FIO2",
    "Blood gas PH",
    "Blood gas PO2",
    "Blood gas PCO2",
    "Blood gas HCO3",
    "Blood gas BE",
    "Blood gas SO2",
    "Blood gas TCO2",
)
ROW_MAJOR_ABG_NO_FIO2_HEADERS = (
    "Blood gas PH",
    "Blood gas PO2",
    "Blood gas PCO2",
    "Blood gas HCO3",
    "Blood gas SO2",
    "Blood gas TCO2",
)
ROW_MAJOR_COAG_HEADERS = (
    ("INR(PT)", "Other Lab Tests"),
    ("PT", "Other Lab Tests"),
    ("APTT", "Other Lab Tests"),
    ("Fibrinogen", "Other Lab Tests"),
)
CANONICAL_ITEM_PREFIXES = (
    ("alkaline phosphatase", "Alkaline Phosphatase"),
    ("direct bilirubin", "Direct Bilirubin"),
    ("total bilirubin", "Total Bilirubin"),
    ("free calcium", "Free Calcium"),
    ("immunoglobulin d", "Immunoglobulin D"),
    ("troponin i", "Troponin I"),
    ("glucose (pc/dextro)", "Glucose (PC/DEXTRO)"),
    ("glucose", "Glucose"),
    ("blood gas be(b)", "Blood gas BE(B)"),
    ("blood gas fcohb", "Blood gas FCOHb"),
    ("blood gas glu", "Blood gas GLU"),
    ("blood gas hco3-act", "Blood gas HCO3-act"),
    ("blood gas hct", "Blood gas HCT"),
    ("blood gas ica", "Blood gas iCa"),
    ("blood gas k", "Blood gas K"),
    ("blood gas lac", "Blood gas Lac"),
    ("blood gas na", "Blood gas Na"),
    ("blood gas pco2", "Blood gas PCO2"),
    ("blood gas ph", "Blood gas PH"),
    ("blood gas po2", "Blood gas PO2"),
    ("blood gas so2", "Blood gas so2"),
    ("ferritin", "Ferritin"),
    ("esr", "ESR"),
    ("tsh", "TSH"),
    ("free t4", "Free T4"),
    ("ldh", "LDH"),
    ("r-gt", "r-GT"),
    ("amylase", "Amylase"),
    ("lipase", "Lipase"),
    ("bun", "BUN"),
    ("creatinine", "Creatinine"),
    ("ast", "AST"),
    ("alt", "ALT"),
    ("na", "Na"),
    ("k", "K"),
    ("cl", "Cl"),
    ("magnesium", "Magnesium"),
    ("crp", "CRP"),
    ("ck", "CK"),
    ("wbc", "WBC"),
    ("rbc", "RBC"),
    ("hb", "Hb"),
    ("hgb", "Hgb"),
    ("hct", "HCT"),
    ("plt", "PLT"),
    ("mcv", "MCV"),
    ("mchc", "MCHC"),
    ("mch", "MCH"),
    ("rdw-sd", "RDW-SD"),
    ("rdw-cv", "RDW-CV"),
    ("neutrophil", "Neutrophil"),
    ("lymphocyte", "Lymphocyte"),
    ("monocyte", "Monocyte"),
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


def _load_mngs_window_entries(patient_dir: Path) -> list[dict[str, Any]]:
    mngs_candidates = list(patient_dir.glob("*_mNGS_grouped.json"))
    mngs_entries = _read_list_json(mngs_candidates[0]) if mngs_candidates else []
    if mngs_entries:
        return mngs_entries

    manifest_candidates = list(patient_dir.glob("*_patient_manifest.json"))
    if not manifest_candidates:
        return []
    manifest = _read_json(manifest_candidates[0])
    return cutoff_entries_from_manifest(manifest)


def _text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def _normalize_datetime_token(token: str) -> str:
    token = _text(token).replace("/", "-")
    short_match = re.fullmatch(
        r"(?P<yy>\d{2})-(?P<mm>\d{2})-(?P<dd>\d{2})(?P<rest>(?:\s+\d{2}:\d{2}(?::\d{2})?)?)",
        token,
    )
    if short_match:
        year = short_match.group("yy")
        month = short_match.group("mm")
        day = short_match.group("dd")
        rest = short_match.group("rest")
        return f"20{year}-{month}-{day}{rest}"
    if len(token) == 10:
        return token
    if len(token) == 16:
        return token
    return token


def _is_datetime(token: str) -> bool:
    return bool(DATETIME_PATTERN.fullmatch(_normalize_datetime_token(token)))


def _looks_like_value(token: str) -> bool:
    token = _text(token).replace(",", "")
    return bool(VALUE_PATTERN.fullmatch(token))


def _looks_like_range(token: str) -> bool:
    token = _text(token).replace(",", "")
    if RANGE_PATTERN.fullmatch(token):
        return True
    return bool(re.fullmatch(r"^-\s*\d+(?:\.\d+)?$", token))


def _looks_like_unit(token: str) -> bool:
    token = _text(token)
    if not token or " " in token:
        return False
    if _looks_like_value(token):
        return False
    if token.lower() in KNOWN_ALPHA_UNITS:
        return True
    if token.isalpha():
        return False
    return bool(UNIT_PATTERN.fullmatch(token))


def _looks_like_item(token: str) -> bool:
    token = _text(token)
    if not token:
        return False
    if _is_datetime(token) or _looks_like_value(token) or _looks_like_range(token) or _looks_like_unit(token):
        return False
    return any(ch.isalpha() for ch in token)


def _short_item_name(item_full: str) -> str:
    cleaned = _canonical_item_label(item_full)
    if not cleaned:
        return ""
    return cleaned


def _canonical_item_label(item_full: str) -> str:
    cleaned = NORMALIZED_SPACE_PATTERN.sub(" ", _text(item_full)).strip()
    if not cleaned:
        return ""
    lower = cleaned.lower()
    if "blood gas" in lower:
        cleaned = cleaned[lower.index("blood gas") :]
        lower = cleaned.lower()
    normalized = re.sub(r"[^a-z0-9()/%.+-]+", " ", lower).strip()
    for prefix, canonical in CANONICAL_ITEM_PREFIXES:
        if normalized.startswith(prefix):
            return canonical
    return cleaned


def _is_cbc_item(item_full: str) -> bool:
    short = _short_item_name(item_full).lower()
    return short in CBC_ITEM_PREFIXES


def _flatten_raw_text(value: Any) -> list[str]:
    flattened: list[str] = []
    if isinstance(value, list):
        for item in value:
            flattened.extend(_flatten_raw_text(item))
        return flattened
    if isinstance(value, dict):
        flattened.extend(_flatten_raw_text(value.get("raw_text", [])))
        return flattened

    text = _text(value)
    if not text:
        return flattened

    if '"raw_text"' in text:
        try:
            payload = json.loads(text)
        except json.JSONDecodeError:
            quoted = [bytes(match.group(1), "utf-8").decode("unicode_escape") for match in QUOTED_TOKEN_PATTERN.finditer(text)]
            if quoted:
                if quoted and quoted[0] == "raw_text":
                    quoted = quoted[1:]
                if "findings" in quoted:
                    quoted = quoted[: quoted.index("findings")]
                return [token for token in quoted if token not in {"raw_text", "findings", "confidence"}]
        else:
            return _flatten_raw_text(payload)

    flattened.append(re.sub(r"\s+", " ", text).strip())
    return flattened


def _table_header_key(token: str) -> str:
    normalized = re.sub(r"\s+", "", _text(token).upper())
    normalized = normalized.replace("O2SAT", "SO2")
    return normalized


def _safe_float(token: str) -> float | None:
    text = _text(token).replace(",", "")
    try:
        return float(text)
    except ValueError:
        return None


def _looks_like_specimen_token(token: str) -> bool:
    normalized = re.sub(r"[^a-z]", "", _text(token).lower())
    if not normalized:
        return False
    return any(normalized.startswith(prefix) for prefix in KNOWN_SPECIMEN_PREFIXES)


def _extract_row_major_cbc_header_specs(tokens: Sequence[str], first_datetime_index: int) -> list[tuple[str, str]]:
    recognized_headers: list[tuple[int, tuple[str, str]]] = []
    for index, token in enumerate(tokens[:first_datetime_index]):
        spec = ROW_MAJOR_CBC_HEADER_MAP.get(_table_header_key(token))
        if spec:
            recognized_headers.append((index, spec))

    if not recognized_headers:
        return []

    header_blocks: list[list[tuple[str, str]]] = []
    current_block: list[tuple[str, str]] = []
    last_index: int | None = None
    for index, spec in recognized_headers:
        if last_index is None or index == last_index + 1:
            current_block.append(spec)
        else:
            header_blocks.append(current_block)
            current_block = [spec]
        last_index = index

    if current_block:
        header_blocks.append(current_block)

    viable_blocks = [block for block in header_blocks if len(block) >= 8]
    if viable_blocks:
        return viable_blocks[-1]
    return header_blocks[-1]


def _is_row_token(token: str) -> bool:
    return bool(ROW_TOKEN_PATTERN.match(_normalize_datetime_token(token)))


def _first_row_major_data_index(tokens: Sequence[str]) -> int | None:
    for index, token in enumerate(tokens):
        if _is_datetime(token) or _is_row_token(token):
            return index
    return None


def _extract_row_major_cbc_entries(tokens: Sequence[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    first_datetime_index = _first_row_major_data_index(tokens)
    if first_datetime_index is None:
        return [], []

    header_specs = _extract_row_major_cbc_header_specs(tokens, first_datetime_index)
    if len(header_specs) < 8:
        return [], []

    cbc_entries: list[dict[str, Any]] = []
    other_entries: list[dict[str, Any]] = []
    index = first_datetime_index
    while index < len(tokens):
        row_token = _text(tokens[index])
        if _is_row_token(row_token):
            normalized_row_token = _normalize_datetime_token(row_token)
            match = ROW_TOKEN_PATTERN.match(normalized_row_token)
            if not match:
                index += 1
                continue

            reported_time = _normalize_datetime_token(f"{match.group('date')} {match.group('time')}")
            row_values = match.group("body").split()
            if row_values and _looks_like_specimen_token(row_values[0]):
                row_values = row_values[1:]
            row_header_specs = header_specs
            if len(row_values) > len(header_specs):
                index += 1
                continue
            if len(row_values) < len(header_specs):
                if len(row_values) >= 8:
                    row_header_specs = header_specs[: len(row_values)]
                elif 1 <= len(row_values) <= len(ROW_MAJOR_COAG_HEADERS):
                    row_header_specs = list(ROW_MAJOR_COAG_HEADERS[: len(row_values)])
                else:
                    index += 1
                    continue

            for (item_name, test_name), value in zip(row_header_specs, row_values):
                normalized_value = _text(value)
                if normalized_value in {"", "-"}:
                    continue
                entry = {
                    "reported_time": reported_time,
                    "item": item_name,
                    "item_full": item_name,
                    "value": normalized_value,
                    "test": test_name,
                }
                if test_name == "CBC":
                    cbc_entries.append(entry)
                else:
                    other_entries.append(entry)

            index += 1
            continue

        if not _is_datetime(tokens[index]):
            index += 1
            continue

        reported_time = _normalize_datetime_token(tokens[index])
        row_start = index + 1
        if row_start < len(tokens) and _looks_like_specimen_token(tokens[row_start]):
            row_start += 1

        next_index = next((probe for probe in range(row_start, len(tokens)) if _is_datetime(tokens[probe])), len(tokens))
        row_values = list(tokens[row_start:next_index])
        row_header_specs = header_specs
        if len(row_values) > len(header_specs):
            index = next_index
            continue
        if len(row_values) < len(header_specs):
            if len(row_values) >= 8:
                row_header_specs = header_specs[: len(row_values)]
            elif 1 <= len(row_values) <= len(ROW_MAJOR_COAG_HEADERS):
                row_header_specs = list(ROW_MAJOR_COAG_HEADERS[: len(row_values)])
            else:
                index = next_index
                continue

        for (item_name, test_name), value in zip(row_header_specs, row_values):
            normalized_value = _text(value)
            if normalized_value in {"", "-"}:
                continue
            entry = {
                "reported_time": reported_time,
                "item": item_name,
                "item_full": item_name,
                "value": normalized_value,
                "test": test_name,
            }
            if test_name == "CBC":
                cbc_entries.append(entry)
            else:
                other_entries.append(entry)

        index = next_index

    return cbc_entries, other_entries


def _extract_row_major_smac_entries(tokens: Sequence[str]) -> list[dict[str, Any]]:
    first_datetime_index = _first_row_major_data_index(tokens)
    if first_datetime_index is None:
        return []

    header_names = [
        _canonical_item_label(token)
        for token in tokens[:first_datetime_index]
        if re.search(r"[A-Za-z]", _text(token))
        and not _is_datetime(token)
        and not _looks_like_value(token)
        and not _looks_like_range(token)
    ]
    header_names = [name for name in header_names if name]
    if len(header_names) < 8:
        return []

    other_entries: list[dict[str, Any]] = []
    index = first_datetime_index
    while index < len(tokens):
        token = _text(tokens[index])
        if _is_row_token(token):
            normalized_row_token = _normalize_datetime_token(token)
            match = ROW_TOKEN_PATTERN.match(normalized_row_token)
            if not match:
                index += 1
                continue
            reported_time = _normalize_datetime_token(f"{match.group('date')} {match.group('time')}")
            row_values = match.group("body").split()
            if row_values and _looks_like_specimen_token(row_values[0]):
                row_values = row_values[1:]
            next_index = index + 1
        elif _is_datetime(token):
            reported_time = _normalize_datetime_token(token)
            row_start = index + 1
            if row_start < len(tokens) and _looks_like_specimen_token(tokens[row_start]):
                row_start += 1
            next_index = next(
                (
                    probe
                    for probe in range(row_start, len(tokens))
                    if _is_datetime(tokens[probe]) or _is_row_token(tokens[probe])
                ),
                len(tokens),
            )
            row_values = list(tokens[row_start:next_index])
        else:
            index += 1
            continue

        if not row_values:
            index = next_index
            continue

        row_header_names = header_names
        if len(row_values) > len(header_names):
            index = next_index
            continue
        if len(row_values) < len(header_names):
            if len(row_values) < 8:
                index = next_index
                continue
            row_header_names = header_names[: len(row_values)]

        for item_name, value in zip(row_header_names, row_values):
            normalized_value = _text(value)
            if normalized_value in {"", "-"}:
                continue
            other_entries.append(
                {
                    "reported_time": reported_time,
                    "item": _short_item_name(item_name),
                    "item_full": item_name,
                    "value": normalized_value,
                    "test": "Other Lab Tests",
                }
            )

        index = next_index

    return other_entries


def _extract_row_major_abg_entries(tokens: Sequence[str]) -> list[dict[str, Any]]:
    row_tokens = [token for token in tokens if ROW_TOKEN_PATTERN.match(_normalize_datetime_token(token))]
    if not row_tokens:
        return []

    other_entries: list[dict[str, Any]] = []
    for token in row_tokens:
        normalized_token = _normalize_datetime_token(token)
        match = ROW_TOKEN_PATTERN.match(normalized_token)
        if not match:
            continue
        reported_time = _normalize_datetime_token(f"{match.group('date')} {match.group('time')}")
        values = match.group("body").split()
        if not values:
            continue

        first_number = _safe_float(values[0])
        if first_number is not None and first_number > 20:
            header_names = ROW_MAJOR_ABG_WITH_FIO2_HEADERS[: len(values)]
        else:
            header_names = ROW_MAJOR_ABG_NO_FIO2_HEADERS[: len(values)]

        if not header_names:
            continue

        for item_name, value in zip(header_names, values):
            normalized_value = _text(value)
            if normalized_value in {"", "-"}:
                continue
            other_entries.append(
                {
                    "reported_time": reported_time,
                    "item": item_name,
                    "item_full": item_name,
                    "value": normalized_value,
                    "test": "Other Lab Tests",
                }
            )

    return other_entries


def _collect_slide_tokens(slide_payload: dict[str, Any]) -> list[str]:
    tokens: list[str] = []
    for image in slide_payload.get("images", []):
        if not isinstance(image, dict):
            continue
        vision = image.get("vision")
        if not isinstance(vision, dict):
            continue
        tokens.extend(_flatten_raw_text(vision.get("raw_text", [])))
    return [token for token in tokens if token]


def parse_lab_entries_from_tokens(tokens: Sequence[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    cbc_entries: list[dict[str, Any]] = []
    other_lab_entries: list[dict[str, Any]] = []
    current_time = ""
    index = 0

    while index < len(tokens):
        token = _text(tokens[index])
        if not token:
            index += 1
            continue

        if _is_datetime(token):
            current_time = _normalize_datetime_token(token)
            index += 1
            continue

        if not _looks_like_item(token):
            index += 1
            continue

        item_parts = [token]
        value_index = None
        probe = index + 1
        while probe < min(len(tokens), index + 4):
            probe_token = _text(tokens[probe])
            if _is_datetime(probe_token):
                break
            if _looks_like_value(probe_token):
                value_index = probe
                break
            if _looks_like_item(probe_token) and not _looks_like_unit(probe_token):
                item_parts.append(probe_token)
            probe += 1

        if value_index is None:
            index += 1
            continue

        item_full = re.sub(r"\s+", " ", " ".join(item_parts)).strip()
        value = _text(tokens[value_index])
        unit = ""
        trailing_time = ""

        next_index = value_index + 1
        if next_index < len(tokens) and _looks_like_range(tokens[next_index]):
            next_index += 1
        if next_index < len(tokens) and _looks_like_unit(tokens[next_index]):
            unit = _text(tokens[next_index])
            next_index += 1
        if next_index < len(tokens) and _is_datetime(tokens[next_index]):
            trailing_time = _normalize_datetime_token(tokens[next_index])
            current_time = trailing_time
            next_index += 1

        reported_time = trailing_time or current_time
        if not reported_time:
            index = value_index + 1
            continue

        entry = {
            "reported_time": reported_time,
            "item": _short_item_name(item_full),
            "item_full": _canonical_item_label(item_full) if not unit else f"{_canonical_item_label(item_full)} ({unit})",
            "value": value,
        }

        if _is_cbc_item(item_full):
            entry["test"] = "CBC"
            cbc_entries.append(entry)
        else:
            entry["test"] = "Other Lab Tests"
            other_lab_entries.append(entry)

        index = max(index + 1, next_index)

    return cbc_entries, other_lab_entries


def extract_lab_entries_from_slide(slide_payload: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    modalities = set(slide_payload.get("candidate_modalities") or [])
    if not {"cbc", "smac", "abg"} & modalities:
        return [], []

    tokens = _collect_slide_tokens(slide_payload)
    row_major_cbc_entries: list[dict[str, Any]] = []
    row_major_other_entries: list[dict[str, Any]] = []
    row_major_smac_entries: list[dict[str, Any]] = []
    if "cbc" in modalities:
        row_major_cbc_entries, row_major_other_entries = _extract_row_major_cbc_entries(tokens)
    if "smac" in modalities:
        row_major_smac_entries = _extract_row_major_smac_entries(tokens)

    if row_major_cbc_entries:
        cbc_entries = row_major_cbc_entries
        other_lab_entries = row_major_other_entries
    elif row_major_smac_entries:
        cbc_entries = []
        other_lab_entries = row_major_smac_entries
    else:
        cbc_entries, other_lab_entries = parse_lab_entries_from_tokens(tokens)

    if "abg" in modalities:
        row_major_abg_entries = _extract_row_major_abg_entries(tokens)
        if row_major_abg_entries:
            other_lab_entries.extend(row_major_abg_entries)

    source = {
        "source_file": _text(slide_payload.get("source_file")),
        "slide_index": slide_payload.get("slide_index"),
    }
    for entry in cbc_entries:
        entry["source"] = source
    for entry in other_lab_entries:
        entry["source"] = source

    return cbc_entries, other_lab_entries


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


def _merge_entries(existing_entries: list[dict[str, Any]], new_entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    merged = list(existing_entries)
    seen = {
        (
            _text(entry.get("reported_time")),
            _text(entry.get("item")),
            _text(entry.get("value")),
            _text((entry.get("source") or {}).get("source_file")),
            (entry.get("source") or {}).get("slide_index"),
        )
        for entry in merged
        if isinstance(entry, dict)
    }
    for entry in new_entries:
        key = (
            _text(entry.get("reported_time")),
            _text(entry.get("item")),
            _text(entry.get("value")),
            _text((entry.get("source") or {}).get("source_file")),
            (entry.get("source") or {}).get("slide_index"),
        )
        if key in seen:
            continue
        seen.add(key)
        merged.append(entry)
    return merged


def import_ppt_labs(patient_export_root: Path) -> dict[str, Any]:
    patient_dirs = [path for path in patient_export_root.iterdir() if path.is_dir()]
    slides_examined = 0
    slides_with_lab = 0
    cbc_written = 0
    other_lab_written = 0
    cbc_filtered_out = 0
    other_lab_filtered_out = 0

    for patient_dir in patient_dirs:
        bundle_candidates = list(patient_dir.glob("*_ppt_slide_bundle.json"))
        if not bundle_candidates:
            continue
        bundle = _read_json(bundle_candidates[0])
        if not isinstance(bundle, list):
            continue

        cbc_path = next(patient_dir.glob("*_CBC.json"))
        other_lab_path = next(patient_dir.glob("*_other_lab.json"))
        mngs_entries = _load_mngs_window_entries(patient_dir)
        existing_cbc = _read_list_json(cbc_path)
        existing_other = _read_list_json(other_lab_path)

        for slide_entry in bundle:
            parsed_json_path = Path(_text(slide_entry.get("parsed_json_path")))
            if not parsed_json_path.exists():
                continue
            slides_examined += 1
            slide_payload = _read_json(parsed_json_path)
            cbc_entries, other_entries = extract_lab_entries_from_slide(slide_payload)
            if not cbc_entries and not other_entries:
                continue

            slides_with_lab += 1
            filtered_cbc_entries = filter_entries_by_mngs_window(cbc_entries, mngs_entries)
            filtered_other_entries = filter_entries_by_mngs_window(other_entries, mngs_entries)
            cbc_filtered_out += len(cbc_entries) - len(filtered_cbc_entries)
            other_lab_filtered_out += len(other_entries) - len(filtered_other_entries)

            source_file = _text(slide_payload.get("source_file"))
            slide_index = slide_payload.get("slide_index")
            existing_cbc = _drop_entries_from_same_slide(existing_cbc, source_file=source_file, slide_index=slide_index)
            existing_other = _drop_entries_from_same_slide(existing_other, source_file=source_file, slide_index=slide_index)

            existing_cbc = _merge_entries(existing_cbc, filtered_cbc_entries)
            existing_other = _merge_entries(existing_other, filtered_other_entries)
            cbc_written += len(filtered_cbc_entries)
            other_lab_written += len(filtered_other_entries)

        _write_json(cbc_path, existing_cbc)
        _write_json(other_lab_path, existing_other)

    return {
        "patient_export_root": str(patient_export_root.resolve()),
        "slides_examined": slides_examined,
        "slides_with_lab_content": slides_with_lab,
        "cbc_entries_written": cbc_written,
        "other_lab_entries_written": other_lab_written,
        "cbc_entries_filtered_out_by_mngs_window": cbc_filtered_out,
        "other_lab_entries_filtered_out_by_mngs_window": other_lab_filtered_out,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Import CBC and other lab entries from PPT slide OCR into the patient export folders.",
    )
    parser.add_argument(
        "--patient-export-root",
        type=Path,
        default=DEFAULT_PATIENT_EXPORT_ROOT,
        help=f"Path to the per-patient export root (default: {DEFAULT_PATIENT_EXPORT_ROOT}).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    patient_export_root = args.patient_export_root.expanduser().resolve()
    summary = import_ppt_labs(patient_export_root)
    summary_path = patient_export_root / SUMMARY_NAME
    _write_json(summary_path, summary)
    print(
        f"Imported {summary['cbc_entries_written']} CBC entries and "
        f"{summary['other_lab_entries_written']} other lab entries."
    )
    print(f"Wrote summary: {summary_path}")


if __name__ == "__main__":
    main()

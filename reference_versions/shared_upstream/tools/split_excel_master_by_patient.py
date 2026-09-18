from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd

from tools.time_window import cutoff_entries_from_dates, filter_entries_by_mngs_window, infer_date_from_specimen_id
from utils import sanitize_filename

BASIC_SHEET = "基本資料"
MICRO_SHEET = "微生物診斷資料"
MNGS_SHEET = "mNGS"
READ_COUNT_SHEET = "Read count"

DEFAULT_EXCEL_DIR = Path("excel")
DEFAULT_OUTPUT_DIR = Path("outputs") / "excel_patient_exports"

STANDARD_SUFFIXES = (
    "underlying",
    "admission_diagnosis",
    "culture",
    "filmarray",
    "gm_test",
    "molecular_microbiology",
    "image",
    "CBC",
    "other_lab",
    "mNGS_grouped",
)

CULTURE_SLOT_MAP = {
    1: ("檢體1", "檢體1說明", "病原菌1", "採檢日1", "報告日1", "抗藥"),
    2: ("檢體2", "檢體2說明", "病原菌2", "採檢日2", "報告日2", "抗藥.1"),
    3: ("檢體3", "檢體3說明", "病原菌3", "採檢日3", "報告日3", "無"),
    4: ("檢體4", "檢體4說明", "病原菌4", "採檢日4", "報告日4", "抗藥.2"),
    5: ("檢體5", "檢體5說明", "病原菌5", "採檢日5", "報告日5", "抗藥.3"),
    6: ("檢體6", "檢體6說明", "病原菌6", "採檢日6", "報告日6", "抗藥.4"),
}

NON_CULTURE_SLOT_MAP = {
    1: ("檢體1.1", "檢體1說明.1", "非培養病原菌1", "採檢日1.1", "報告日1.1"),
    2: ("檢體2.1", "檢體2說明.1", "非培養病原菌2", "採檢日2.1", "報告日2.1"),
    3: ("檢體3.1", "檢體3說明.1", "非培養病原菌3", "採檢日3.1", "報告日3.1"),
    4: ("檢體4.1", "檢體4說明.1", "非培養病原菌4", "採檢日4.1", "報告日4.1"),
    5: ("檢體5.1", "檢體5說明.1", "非培養病原菌5", "採檢日5.1", "報告日5.1"),
    6: ("檢體6.1", "檢體6說明.1", "非培養病原菌6", "採檢日6.1", "報告日6.1"),
    7: ("檢體7", "檢體7說明", "非培養病原菌7", "採檢日7", "報告日7"),
    8: ("檢體8", "檢體8說明˙", "非培養病原菌8", "採檢日8", "報告日8"),
}

MNGS_PATHOGEN_COLUMNS = [f"pathogen {index}" for index in range(1, 11)]
READ_COUNT_COLUMNS = [(f"report_pathogens_{index}", f"Read_{index}") for index in range(1, 12)]
CLINICAL_DIAGNOSIS_COLUMNS = [f"臨床病原診斷{index}" for index in range(1, 9)]

SKIP_VALUES = {"", "-", "nan", "none", "病原總數", "blood", "sputum", "balf", "urine"}
VIRUS_HINTS = (
    "virus",
    "viral",
    "cmv",
    "ebv",
    "hbv",
    "hcv",
    "hsv",
    "hhv",
    "influenza",
    "coronavirus",
    "rsv",
    "adenovirus",
    "rhinovirus",
    "enterovirus",
    "parainfluenza",
    "metapneumovirus",
    "bocavirus",
    "norovirus",
    "rotavirus",
    "astrovirus",
)
FUNGAL_HINTS = (
    "candida",
    "asperg",
    "pneumocystis",
    "cryptococ",
    "mucor",
    "rhizopus",
    "cunninghamella",
    "fusarium",
    "trichosporon",
    "yeast",
    "fung",
)
GM_TEST_HINTS = (
    "galactomannan",
    " antigen",
    " ag",
    " ab",
    "igm",
    "igg",
    "legionella",
    "chlamydophila",
    "mycoplasma",
)
FILMARRAY_HINTS = (
    "-pn",
    " panel",
    "ctx-m",
    "kpc",
    "ndm",
    "oxa-48",
    "vim",
    "imp",
)
PATIENT_CODE_PATTERN = re.compile(r"^[A-Za-z]{1,3}\d{3}$")


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    text = str(value).strip()
    if text.endswith(".0") and re.fullmatch(r"-?\d+\.0", text):
        return text[:-2]
    return text


def _format_datetime(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, pd.Timestamp):
        if pd.isna(value):
            return ""
        value = value.to_pydatetime()
    if hasattr(value, "strftime"):
        if getattr(value, "hour", 0) == 0 and getattr(value, "minute", 0) == 0 and getattr(value, "second", 0) == 0:
            return value.strftime("%Y-%m-%d")
        if getattr(value, "second", 0) == 0 and getattr(value, "microsecond", 0) == 0:
            return value.strftime("%Y-%m-%d %H:%M")
        return value.strftime("%Y-%m-%d %H:%M:%S")
    return _text(value)


def _is_yes(value: Any) -> bool:
    return _text(value).lower() in {"yes", "y", "true", "1", "是", "有"}


def _normalize_gender(value: Any) -> str:
    text = _text(value).lower()
    if text in {"男", "m", "male"}:
        return "M"
    if text in {"女", "f", "female"}:
        return "F"
    return _text(value)


def _compose_sample(sample: Any, description: Any) -> str:
    sample_text = _text(sample)
    description_text = _text(description)
    if description_text and description_text.lower() != sample_text.lower():
        if sample_text and sample_text != "其他":
            return f"{sample_text} ({description_text})"
        return description_text
    return sample_text


def _split_organisms(raw_value: Any) -> list[str]:
    text = _text(raw_value)
    if not text:
        return []
    parts = re.split(r"[;,、]", text)
    return [part.strip() for part in parts if part.strip()]


def _split_detected_targets(raw_value: Any) -> list[tuple[str, str | None]]:
    targets: list[tuple[str, str | None]] = []
    for part in _split_organisms(raw_value):
        name, value = _parse_name_with_optional_value(part)
        if name:
            targets.append((name, value))
    return targets


def _parse_name_with_optional_value(raw_value: Any) -> tuple[str, str | None]:
    part = _text(raw_value)
    lowered = part.lower()
    if not part or lowered in SKIP_VALUES:
        return "", None

    match = re.match(r"^(?P<name>.+?)\((?P<value>[^()]+)\)$", part)
    if match:
        name = match.group("name").strip()
        value = match.group("value").strip()
    else:
        name = part.strip()
        value = None

    if name.lower().endswith("-pn"):
        name = name[:-3].rstrip(" -")
    return name, value


def _parse_optional_int(raw_value: Any) -> int | None:
    text = _text(raw_value).replace(",", "")
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _guess_mngs_category(pathogen_name: str) -> str:
    lowered = pathogen_name.lower()
    if "haemophilus influenzae" in lowered:
        return "bacterial"
    if any(hint in lowered for hint in VIRUS_HINTS):
        return "viral"
    if any(hint in lowered for hint in FUNGAL_HINTS):
        return "fungal"
    if "mycobacter" in lowered:
        return "others"
    return "bacterial"


def _is_gm_test_target(target: str) -> bool:
    lowered = f" {target.lower()} "
    return any(hint in lowered for hint in GM_TEST_HINTS)


def _is_filmarray_target(target: str, description: str) -> bool:
    lowered = f"{target} {description}".lower()
    return any(hint in lowered for hint in FILMARRAY_HINTS)


def _normalize_patient_identifier(value: Any) -> str:
    return _text(value)


def _patient_slug(case_code: str, specimen_id: str, hospital_id: str) -> str:
    slug = "__".join(part for part in (case_code, specimen_id, hospital_id) if part)
    return sanitize_filename(slug) or "unknown_patient"


def _is_patient_code(value: Any) -> bool:
    return bool(PATIENT_CODE_PATTERN.fullmatch(_text(value)))


def _read_sheet(excel_path: Path, sheet_name: str) -> pd.DataFrame:
    frame = pd.read_excel(excel_path, sheet_name=sheet_name, dtype=object)
    return frame.fillna("")


def load_master_workbook(excel_path: Path) -> dict[str, pd.DataFrame]:
    return {
        BASIC_SHEET: _read_sheet(excel_path, BASIC_SHEET),
        MICRO_SHEET: _read_sheet(excel_path, MICRO_SHEET),
        MNGS_SHEET: _read_sheet(excel_path, MNGS_SHEET),
        READ_COUNT_SHEET: _read_sheet(excel_path, READ_COUNT_SHEET),
    }


def build_underlying_json(basic_row: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not basic_row:
        return []

    diseases: list[str] = []
    if _is_yes(basic_row.get("糖尿病")):
        diseases.append("Diabetes mellitus")
    if _is_yes(basic_row.get("慢性腎衰竭")):
        diseases.append("Chronic kidney disease")
    if _is_yes(basic_row.get("ESRD洗腎")):
        diseases.append("ESRD on dialysis")
    if _is_yes(basic_row.get("肝硬化")):
        diseases.append("Liver cirrhosis")
    if _is_yes(basic_row.get("自體免疫疾病")):
        diseases.append("Autoimmune disease")
    if _is_yes(basic_row.get("活動性腫瘤")):
        tumor_site = _text(basic_row.get("腫瘤部位"))
        if tumor_site:
            diseases.append(f"Active malignancy ({tumor_site})")
        else:
            diseases.append("Active malignancy")
    if _is_yes(basic_row.get("類風溼性關節炎")):
        diseases.append("Rheumatoid arthritis")
    if _is_yes(basic_row.get("紅斑性狼瘡")):
        diseases.append("Systemic lupus erythematosus")
    other_autoimmune = _text(basic_row.get("其他自體免疫疾病"))
    if other_autoimmune and other_autoimmune not in {"否", "無"}:
        diseases.append(other_autoimmune)
    if _is_yes(basic_row.get("全身性類固醇")):
        diseases.append("Systemic steroid exposure")
    if _is_yes(basic_row.get("抽菸史")):
        diseases.append("Smoking history")

    medical_history = {
        "diabetes_mellitus": _text(basic_row.get("糖尿病")),
        "chronic_kidney_disease": _text(basic_row.get("慢性腎衰竭")),
        "esrd_on_dialysis": _text(basic_row.get("ESRD洗腎")),
        "liver_cirrhosis": _text(basic_row.get("肝硬化")),
        "autoimmune_disease": _text(basic_row.get("自體免疫疾病")),
        "active_malignancy": _text(basic_row.get("活動性腫瘤")),
        "tumor_site": _text(basic_row.get("腫瘤部位")),
        "rheumatoid_arthritis": _text(basic_row.get("類風溼性關節炎")),
        "systemic_lupus_erythematosus": _text(basic_row.get("紅斑性狼瘡")),
        "other_autoimmune_disease": _text(basic_row.get("其他自體免疫疾病")),
        "systemic_steroid_exposure": _text(basic_row.get("全身性類固醇")),
        "smoking_history": _text(basic_row.get("抽菸史")),
    }
    medical_history = {key: value for key, value in medical_history.items() if value}

    clinical_context = {
        "ward_unit": _text(basic_row.get("病房單位")),
        "admission_date": _format_datetime(basic_row.get("住院日")),
        "ards_on_specimen_day": _text(basic_row.get("送檢日ARDS")),
        "pneumonia_type": _text(basic_row.get("肺炎種類")),
        "pf_ratio_on_specimen_day": _text(basic_row.get("送檢日PF ratio")),
    }
    clinical_context = {key: value for key, value in clinical_context.items() if value}

    payload = {
        "test": "Underlying Conditions",
        "age": _text(basic_row.get("年齡")),
        "gender": _normalize_gender(basic_row.get("性別")),
        "birth_date": _format_datetime(basic_row.get("生日")),
        "height": _text(basic_row.get("身高")),
        "weight": _text(basic_row.get("體重")),
        "underlying_diseases": diseases,
    }
    if medical_history:
        payload["medical_history"] = medical_history
    if clinical_context:
        payload["clinical_context"] = clinical_context

    return [
        {key: value for key, value in payload.items() if value is not None and value != ""}
    ]


def build_culture_json(micro_row: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not micro_row:
        return []

    entries: list[dict[str, Any]] = []
    for sample_col, desc_col, organism_col, collected_col, reported_col, resistance_col in CULTURE_SLOT_MAP.values():
        sample = _compose_sample(micro_row.get(sample_col), micro_row.get(desc_col))
        organisms = _split_organisms(micro_row.get(organism_col))
        if not sample or not organisms:
            continue
        resistance = _text(micro_row.get(resistance_col))
        for organism in organisms:
            entry = {
                "test": "Microbiology culture",
                "sample": sample,
                "collected_time": _format_datetime(micro_row.get(collected_col)),
                "received_time": "",
                "reported_time": _format_datetime(micro_row.get(reported_col)),
                "organism": organism,
                "status": "Detected",
                "antimicrobial_susceptibility_test": [],
            }
            if resistance and resistance not in {"無", "否"}:
                entry["notes"] = [f"Resistance marker noted in source sheet: {resistance}"]
            entries.append(entry)
    return entries


def build_non_culture_json(
    micro_row: Mapping[str, Any] | None,
    *,
    case_code: str = "",
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    if not micro_row:
        return [], [], []

    filmarray_entries: list[dict[str, Any]] = []
    gm_entries: list[dict[str, Any]] = []
    molecular_entries: list[dict[str, Any]] = []
    case_prefix = _text(case_code).upper()[:1]
    strict_filmarray_classification = case_prefix in {"T", "V"}

    for sample_col, desc_col, target_col, collected_col, reported_col in NON_CULTURE_SLOT_MAP.values():
        sample = _compose_sample(micro_row.get(sample_col), micro_row.get(desc_col))
        description = _text(micro_row.get(desc_col))
        reported_time = _format_datetime(micro_row.get(reported_col))
        collected_time = _format_datetime(micro_row.get(collected_col))
        for target, value in _split_detected_targets(micro_row.get(target_col)):
            is_gm_target = _is_gm_test_target(target)
            is_filmarray_target = _is_filmarray_target(target, description)
            if is_gm_target:
                test_name = "GM / Antigen / Serology"
            elif is_filmarray_target:
                test_name = "FilmArray / Molecular Panel"
            else:
                test_name = "Non-panel Molecular Test"

            entry = {
                "test": test_name,
                "sample": sample,
                "collected_time": collected_time,
                "reported_time": reported_time or collected_time,
                "target": target,
                "result": "Detected",
            }
            if value:
                entry["value"] = value
            if is_gm_target:
                gm_entries.append(entry)
            elif strict_filmarray_classification and not is_filmarray_target:
                molecular_entries.append(entry)
            else:
                filmarray_entries.append(entry)

    return filmarray_entries, gm_entries, molecular_entries


def build_admission_diagnosis_json(basic_row: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not basic_row:
        return []

    entries: list[dict[str, Any]] = []
    pneumonia_type = _text(basic_row.get("肺炎種類"))
    if pneumonia_type:
        entries.append({"test": "Admission diagnosis", "diagnosis": f"Pneumonia type: {pneumonia_type}"})

    ards = _text(basic_row.get("送檢日ARDS"))
    if ards and ards not in {"無", "否"}:
        entries.append({"test": "Admission diagnosis", "diagnosis": f"ARDS on specimen day: {ards}"})

    ward = _text(basic_row.get("病房單位"))
    if ward:
        entries.append({"test": "Admission diagnosis", "diagnosis": f"Ward / unit: {ward}"})

    return entries


def _build_read_count_lookup(read_count_frame: pd.DataFrame) -> tuple[dict[tuple[str, str], dict[str, int]], dict[str, dict[str, int]]]:
    exact_lookup: dict[tuple[str, str], dict[str, int]] = {}
    specimen_lookup: dict[str, dict[str, int]] = {}

    for _, row in read_count_frame.iterrows():
        specimen_id = _normalize_patient_identifier(row.get("specimen_id"))
        hospital_id = _normalize_patient_identifier(row.get("hospital_patient_id"))
        if not specimen_id:
            continue

        exact_bucket = exact_lookup.setdefault((specimen_id, hospital_id), defaultdict(int))
        specimen_bucket = specimen_lookup.setdefault(specimen_id, defaultdict(int))
        for pathogen_col, read_col in READ_COUNT_COLUMNS:
            pathogen, _ = _parse_name_with_optional_value(row.get(pathogen_col))
            reads = _parse_optional_int(row.get(read_col))
            if not pathogen or reads is None:
                continue
            # Read count sheets can repeat the same pathogen across duplicate rows for the
            # same specimen/hospital pair. Keeping the maximum avoids inflating reads by
            # summing duplicated report rows.
            exact_bucket[pathogen] = max(exact_bucket[pathogen], reads)
            specimen_bucket[pathogen] = max(specimen_bucket[pathogen], reads)

    return exact_lookup, specimen_lookup


def build_mngs_grouped_json(
    mngs_rows: Sequence[Mapping[str, Any]],
    read_count_frame: pd.DataFrame,
) -> list[dict[str, Any]]:
    exact_lookup, specimen_lookup = _build_read_count_lookup(read_count_frame)
    entries: list[dict[str, Any]] = []
    seen_keys: set[tuple[str, str, str]] = set()

    for row in mngs_rows:
        specimen_id = _normalize_patient_identifier(row.get("specimen_id"))
        hospital_id = _normalize_patient_identifier(row.get("hospital_patient_id"))
        specimen_type = _text(row.get("specimen_type"))
        collected_time = _format_datetime(row.get("arrival_date"))
        key = (specimen_id, hospital_id, collected_time)
        if not specimen_id or key in seen_keys:
            continue
        seen_keys.add(key)

        exact_reads = exact_lookup.get((specimen_id, hospital_id), {})
        fallback_reads = specimen_lookup.get(specimen_id, {})
        grouped = {"bacterial": [], "viral": [], "fungal": [], "others": []}

        for pathogen_col in MNGS_PATHOGEN_COLUMNS:
            pathogen, inline_value = _parse_name_with_optional_value(row.get(pathogen_col))
            if not pathogen:
                continue
            inline_reads = _parse_optional_int(inline_value)
            reads = exact_reads.get(pathogen, fallback_reads.get(pathogen, inline_reads or 0))
            category = _guess_mngs_category(pathogen)
            grouped[category].append({"name": pathogen, "reads": reads})

        if not any(grouped.values()):
            continue

        entries.append(
            {
                "specimen_code": specimen_id,
                "collected_time": collected_time,
                "specimen_site": specimen_type,
                "pathogens": grouped,
            }
        )

    return entries


def build_patient_manifest(
    workbook_name: str,
    patient_key: str,
    case_code: str,
    specimen_id: str,
    hospital_id: str,
    basic_row: Mapping[str, Any] | None,
    micro_row: Mapping[str, Any] | None,
    mngs_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    clinical_pathogen_diagnoses: list[str] = []
    if micro_row:
        for column in CLINICAL_DIAGNOSIS_COLUMNS:
            diagnosis = _text(micro_row.get(column))
            if diagnosis and diagnosis not in {"無", "不明"}:
                clinical_pathogen_diagnoses.append(diagnosis)
    mngs_collected_dates = sorted(
        {
            formatted
            for row in mngs_rows
            if (formatted := _format_datetime(row.get("arrival_date")))
        }
    )

    return {
        "patient_key": patient_key,
        "source_workbook": workbook_name,
        "mngs_collected_dates": mngs_collected_dates,
        "identifiers": {
            "case_code": case_code,
            "specimen_id": specimen_id,
            "hospital_id": hospital_id,
            "hospital_code": _text((basic_row or {}).get("醫院")),
            "hospital_name": _text(mngs_rows[0].get("hospital")) if mngs_rows else "",
        },
        "clinical_context": {
            "ward_unit": _text((basic_row or {}).get("病房單位")),
            "admission_date": _format_datetime((basic_row or {}).get("住院日")),
            "pneumonia_type": _text((basic_row or {}).get("肺炎種類")),
            "ards_on_specimen_day": _text((basic_row or {}).get("送檢日ARDS")),
            "pf_ratio_on_specimen_day": _text((basic_row or {}).get("送檢日PF ratio")),
            "treatment_changed": _text((micro_row or {}).get("微生物治療改變")),
            "treatment_change_mode": _text((micro_row or {}).get("治療變化方式")),
            "mngs_impacted_management": _text((micro_row or {}).get("mNGS是否影響")),
            "clinical_pathogen_diagnoses": clinical_pathogen_diagnoses,
            "ventilator_within_3_days": _text((micro_row or {}).get("送檢日三日內呼吸器使用")),
            "ventilator_end_date": _format_datetime((micro_row or {}).get("呼吸器使用結束日")),
            "discharge_date": _format_datetime((micro_row or {}).get("出院日")),
            "discharge_status": _text((micro_row or {}).get("出院狀態")),
        },
        "source_row_counts": {
            "basic_rows": 1 if basic_row else 0,
            "micro_rows": 1 if micro_row else 0,
            "mngs_rows": len(mngs_rows),
        },
    }


def _group_frame_rows(frame: pd.DataFrame, key_columns: Sequence[str]) -> dict[tuple[str, str, str], list[dict[str, Any]]]:
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for _, row in frame.iterrows():
        payload = {str(column): row[column] for column in frame.columns}
        if not _is_patient_code(payload.get(key_columns[0])):
            continue
        key = tuple(_normalize_patient_identifier(payload.get(column)) for column in key_columns)
        if not any(key):
            continue
        grouped[key].append(payload)
    return grouped


def export_patient_workbook(
    excel_path: Path,
    output_dir: Path,
) -> tuple[list[Path], Path]:
    sheets = load_master_workbook(excel_path)
    basic_groups = _group_frame_rows(sheets[BASIC_SHEET], ("編碼", "Specimen_ID", "Hospital_ID"))
    micro_groups = _group_frame_rows(sheets[MICRO_SHEET], ("編碼", "Specimen_ID", "Hospital_ID"))
    mngs_groups = _group_frame_rows(sheets[MNGS_SHEET], ("編碼", "specimen_id", "hospital_patient_id"))

    all_keys = sorted(set(basic_groups) | set(micro_groups) | set(mngs_groups))
    output_dir.mkdir(parents=True, exist_ok=True)

    index_payload: list[dict[str, Any]] = []
    written_dirs: list[Path] = []

    for case_code, specimen_id, hospital_id in all_keys:
        patient_key = _patient_slug(case_code, specimen_id, hospital_id)
        patient_dir = output_dir / patient_key
        patient_dir.mkdir(parents=True, exist_ok=True)

        basic_row = basic_groups.get((case_code, specimen_id, hospital_id), [{}])[0]
        micro_row = micro_groups.get((case_code, specimen_id, hospital_id), [{}])[0]
        mngs_rows = mngs_groups.get((case_code, specimen_id, hospital_id), [])

        payloads: dict[str, Any] = {
            "underlying": build_underlying_json(basic_row if basic_row else None),
            "admission_diagnosis": build_admission_diagnosis_json(basic_row if basic_row else None),
            "culture": build_culture_json(micro_row if micro_row else None),
            "image": [],
            "CBC": [],
            "other_lab": [],
            "mNGS_grouped": build_mngs_grouped_json(mngs_rows, sheets[READ_COUNT_SHEET]),
        }
        filmarray_entries, gm_entries, molecular_entries = build_non_culture_json(
            micro_row if micro_row else None,
            case_code=case_code,
        )
        payloads["filmarray"] = filmarray_entries
        payloads["gm_test"] = gm_entries
        payloads["molecular_microbiology"] = molecular_entries

        mngs_window_entries = cutoff_entries_from_dates(_format_datetime(row.get("arrival_date")) for row in mngs_rows)
        if not mngs_window_entries:
            inferred_specimen_date = infer_date_from_specimen_id(specimen_id)
            if inferred_specimen_date is not None:
                mngs_window_entries = [{"collected_time": inferred_specimen_date.isoformat()}]
        for suffix in ("culture", "filmarray", "gm_test", "molecular_microbiology", "image", "CBC", "other_lab"):
            payloads[suffix] = filter_entries_by_mngs_window(payloads[suffix], mngs_window_entries)

        manifest = build_patient_manifest(
            workbook_name=excel_path.name,
            patient_key=patient_key,
            case_code=case_code,
            specimen_id=specimen_id,
            hospital_id=hospital_id,
            basic_row=basic_row if basic_row else None,
            micro_row=micro_row if micro_row else None,
            mngs_rows=mngs_rows,
        )

        for suffix in STANDARD_SUFFIXES:
            output_path = patient_dir / f"{patient_key}_{suffix}.json"
            output_path.write_text(json.dumps(payloads[suffix], ensure_ascii=False, indent=2), encoding="utf-8")

        manifest_path = patient_dir / f"{patient_key}_patient_manifest.json"
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

        index_payload.append(
            {
                "patient_key": patient_key,
                "case_code": case_code,
                "specimen_id": specimen_id,
                "hospital_id": hospital_id,
                "output_dir": str(patient_dir),
            }
        )
        written_dirs.append(patient_dir)

    index_path = output_dir / "patient_index.json"
    index_path.write_text(json.dumps(index_payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return written_dirs, index_path


def _resolve_default_excel_path(excel_path: Path | None) -> Path:
    if excel_path is not None:
        return excel_path.expanduser().resolve()
    candidates = sorted(DEFAULT_EXCEL_DIR.glob("*.xlsx"))
    if not candidates:
        raise FileNotFoundError(f"No .xlsx files found under {DEFAULT_EXCEL_DIR}")
    return candidates[0].resolve()


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split the Excel master workbook into per-patient JSON files aligned with the existing project layout.",
    )
    parser.add_argument(
        "--excel",
        type=Path,
        default=None,
        help="Path to the master workbook. Defaults to the first .xlsx file under ./excel.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Where to write the per-patient folders (default: {DEFAULT_OUTPUT_DIR}).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    excel_path = _resolve_default_excel_path(args.excel)
    output_dir = args.output.expanduser().resolve()

    if not excel_path.exists():
        raise FileNotFoundError(f"Excel workbook not found: {excel_path}")

    written_dirs, index_path = export_patient_workbook(excel_path, output_dir)
    print(f"Exported {len(written_dirs)} patient folders to {output_dir}")
    print(f"Wrote patient index: {index_path}")


if __name__ == "__main__":
    main()

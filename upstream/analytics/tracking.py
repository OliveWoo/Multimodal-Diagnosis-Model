from __future__ import annotations

import csv
import logging
import re
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence

LOG = logging.getLogger(__name__)

ANALYTICS_DIR = Path(__file__).resolve().parent
DATA_DIR = ANALYTICS_DIR / "data"
FIGURES_DIR = ANALYTICS_DIR / "figures"

ANALYZE_DATASET = DATA_DIR / "analyze_agent_results.csv"
SUMMARY_DATASET = DATA_DIR / "summarize_agent_results.csv"

# Known agent suffixes used by tools/analyze_llm_agent.py
AGENT_SUFFIXES = {
    "underlying": "_underlying_agent.json",
    "cbc_other_lab": "_cbc_other_lab_agent.json",
    "filmarray_gmtest": "_filmarray_gmtest_agent.json",
    "image": "_image_agent.json",
    "culture": "_culture_agent.json",
}

CATEGORY_FIELDNAMES = [
    "timestamp",
    "patient_id",
    "category",
    "rule_version",
    "infection_likelihood",
    "probable_source",
    "probable_pathogen_type",
    "support_direction",
    "host_vulnerability_tier",
    "opportunistic_coverage_level",
    "candidate_count",
    "top_candidate_level",
]

SUMMARY_FIELDNAMES = [
    "timestamp",
    "patient_id",
    "rule_version",
    "final_infection_likelihood",
    "dominant_source",
    "dominant_pathogen_type",
    "support_direction",
    "data_quality_tier",
    "amr_risk_level",
    "host_vulnerability_tier",
    "opportunistic_coverage_level",
    "candidate_count",
    "top_candidate_level",
]


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _append_row(dataset: Path, fieldnames: Sequence[str], row: dict) -> Path:
    _ensure_parent(dataset)
    file_exists = dataset.exists()
    with dataset.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)
    return dataset


def _parse_level_number(value: object) -> int | None:
    if value is None:
        return None
    if isinstance(value, int):
        return value
    match = re.search(r"(\d+)", str(value))
    return int(match.group(1)) if match else None


def _top_level_from_entries(entries: object, *, key: str) -> str:
    if not isinstance(entries, list) or not entries:
        return "Unknown"
    levels = [_parse_level_number(entry.get(key)) for entry in entries if isinstance(entry, dict)]
    levels = [level for level in levels if level is not None]
    if not levels:
        return "Unknown"
    return f"Level {min(levels)}"


def _safe_len_list(value: object) -> int:
    return len(value) if isinstance(value, list) else 0


def record_category_result(patient_id: str, category: str, payload: dict) -> Path:
    """
    Persist a single category agent result into analytics/data/analyze_agent_results.csv.

    Extracts common fields if present; missing values are stored as empty strings.
    """
    timestamp = datetime.now().isoformat()
    host_state = payload.get("host_state") if isinstance(payload.get("host_state"), dict) else {}
    organism_labels = payload.get("organism_labels")
    row = {
        "timestamp": timestamp,
        "patient_id": patient_id,
        "category": category,
        "rule_version": payload.get("rule_version", ""),
        "infection_likelihood": payload.get("infection_likelihood", ""),
        "probable_source": payload.get("probable_source", ""),
        "probable_pathogen_type": payload.get("probable_pathogen_type", ""),
        "support_direction": payload.get("support_direction", ""),
        "host_vulnerability_tier": host_state.get("host_vulnerability_tier", ""),
        "opportunistic_coverage_level": host_state.get("opportunistic_coverage_level", ""),
        "candidate_count": _safe_len_list(organism_labels),
        "top_candidate_level": _top_level_from_entries(organism_labels, key="causative_level"),
    }
    return _append_row(ANALYZE_DATASET, CATEGORY_FIELDNAMES, row)


def record_summary_result(patient_id: str, payload: dict) -> Path:
    """
    Persist summarize_agent output into analytics/data/summarize_agent_results.csv.
    """
    timestamp = datetime.now().isoformat()
    host_context = (
        payload.get("host_context") if isinstance(payload.get("host_context"), dict) else {}
    )
    candidates = payload.get("pathogen_candidates")
    amr_summary = (
        payload.get("amr_risk_summary") if isinstance(payload.get("amr_risk_summary"), dict) else {}
    )
    row = {
        "timestamp": timestamp,
        "patient_id": patient_id,
        "rule_version": payload.get("rule_version", ""),
        "final_infection_likelihood": payload.get("final_infection_likelihood", ""),
        "dominant_source": payload.get("dominant_source", ""),
        "dominant_pathogen_type": payload.get("dominant_pathogen_type", ""),
        "support_direction": payload.get("support_direction", ""),
        "data_quality_tier": payload.get("data_quality_tier", ""),
        "amr_risk_level": amr_summary.get("risk_level", ""),
        "host_vulnerability_tier": host_context.get("host_vulnerability_tier", ""),
        "opportunistic_coverage_level": host_context.get("opportunistic_coverage_level", ""),
        "candidate_count": _safe_len_list(candidates),
        "top_candidate_level": _top_level_from_entries(
            candidates, key="integrated_causative_level"
        ),
    }
    return _append_row(SUMMARY_DATASET, SUMMARY_FIELDNAMES, row)


def _extract_patient_id(path: Path) -> str:
    name = path.name
    if name.endswith(".json"):
        name = name[:-5]
    name = re.sub(r"_final_summary.*$", "", name)
    for suffix in AGENT_SUFFIXES.values():
        if name.endswith(suffix.replace(".json", "")):
            name = name[: -len(suffix.replace(".json", ""))]
            break
    return name


def _iter_paths(root: Path, patterns: Iterable[str]) -> Iterable[Path]:
    for pattern in patterns:
        yield from root.rglob(pattern)


def rebuild_category_dataset(inputs: Sequence[Path]) -> Path:
    """
    Re-scan patient directories for agent outputs and rebuild the category CSV.
    """
    rows: list[dict] = []
    for raw in inputs:
        base = raw.expanduser().resolve()
        if not base.exists():
            LOG.warning("Skipping non-existent path: %s", base)
            continue
        search_roots = [base]
        if base.is_dir():
            try:
                search_roots.extend(child for child in base.iterdir() if child.is_dir())
            except PermissionError:
                LOG.warning("Skipping %s (permission denied)", base)

        for root in search_roots:
            agent_dir = root / "agent_outputs"
            if not agent_dir.is_dir():
                continue
            for category, suffix in AGENT_SUFFIXES.items():
                for path in agent_dir.glob(f"*{suffix}"):
                    try:
                        payload = _load_json(path)
                    except Exception as exc:  # pragma: no cover - defensive
                        LOG.warning("Failed to load %s: %s", path, exc)
                        continue
                    patient_id = _extract_patient_id(path)
                    rows.append(
                        {
                            "timestamp": datetime.now().isoformat(),
                            "patient_id": patient_id,
                            "category": category,
                            "rule_version": payload.get("rule_version", ""),
                            "infection_likelihood": payload.get("infection_likelihood", ""),
                            "probable_source": payload.get("probable_source", ""),
                            "probable_pathogen_type": payload.get("probable_pathogen_type", ""),
                            "support_direction": payload.get("support_direction", ""),
                            "host_vulnerability_tier": (
                                payload.get("host_state", {}).get("host_vulnerability_tier", "")
                                if isinstance(payload.get("host_state"), dict)
                                else ""
                            ),
                            "opportunistic_coverage_level": (
                                payload.get("host_state", {}).get(
                                    "opportunistic_coverage_level", ""
                                )
                                if isinstance(payload.get("host_state"), dict)
                                else ""
                            ),
                            "candidate_count": _safe_len_list(payload.get("organism_labels")),
                            "top_candidate_level": _top_level_from_entries(
                                payload.get("organism_labels"), key="causative_level"
                            ),
                        }
                    )
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with ANALYZE_DATASET.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CATEGORY_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return ANALYZE_DATASET


def rebuild_summary_dataset(inputs: Sequence[Path]) -> Path:
    """
    Re-scan patient directories for *_final_summary.json files and rebuild the summary CSV.
    """
    rows: list[dict] = []
    for raw in inputs:
        base = raw.expanduser().resolve()
        if not base.exists():
            LOG.warning("Skipping non-existent path: %s", base)
            continue
        search_roots = [base]
        if base.is_dir():
            try:
                search_roots.extend(child for child in base.iterdir() if child.is_dir())
            except PermissionError:
                LOG.warning("Skipping %s (permission denied)", base)

        for root in search_roots:
            for path in _iter_paths(root, ["*_final_summary*.json"]):
                try:
                    payload = _load_json(path)
                except Exception as exc:  # pragma: no cover - defensive
                    LOG.warning("Failed to load %s: %s", path, exc)
                    continue
                patient_id = _extract_patient_id(path)
                rows.append(
                    {
                        "timestamp": datetime.now().isoformat(),
                        "patient_id": patient_id,
                        "rule_version": payload.get("rule_version", ""),
                        "final_infection_likelihood": payload.get("final_infection_likelihood", ""),
                        "dominant_source": payload.get("dominant_source", ""),
                        "dominant_pathogen_type": payload.get("dominant_pathogen_type", ""),
                        "support_direction": payload.get("support_direction", ""),
                        "data_quality_tier": payload.get("data_quality_tier", ""),
                        "amr_risk_level": (
                            payload.get("amr_risk_summary", {}).get("risk_level", "")
                            if isinstance(payload.get("amr_risk_summary"), dict)
                            else ""
                        ),
                        "host_vulnerability_tier": (
                            payload.get("host_context", {}).get("host_vulnerability_tier", "")
                            if isinstance(payload.get("host_context"), dict)
                            else ""
                        ),
                        "opportunistic_coverage_level": (
                            payload.get("host_context", {}).get("opportunistic_coverage_level", "")
                            if isinstance(payload.get("host_context"), dict)
                            else ""
                        ),
                        "candidate_count": _safe_len_list(payload.get("pathogen_candidates")),
                        "top_candidate_level": _top_level_from_entries(
                            payload.get("pathogen_candidates"), key="integrated_causative_level"
                        ),
                    }
                )
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with SUMMARY_DATASET.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return SUMMARY_DATASET


def _load_json(path: Path) -> dict:
    import json

    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)

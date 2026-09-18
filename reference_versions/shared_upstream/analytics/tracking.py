from __future__ import annotations

import csv
import logging
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


def record_category_result(patient_id: str, category: str, payload: dict) -> Path:
    """
    Persist a single category agent result into analytics/data/analyze_agent_results.csv.

    Extracts common fields if present; missing values are stored as empty strings.
    """
    timestamp = datetime.now().isoformat()
    row = {
        "timestamp": timestamp,
        "patient_id": patient_id,
        "category": category,
        "infection_likelihood": payload.get("infection_likelihood", ""),
        "probable_source": payload.get("probable_source", ""),
        "probable_pathogen_type": payload.get("probable_pathogen_type", ""),
        "confidence_score": payload.get("confidence_score", ""),
    }
    return _append_row(
        ANALYZE_DATASET,
        ["timestamp", "patient_id", "category", "infection_likelihood", "probable_source", "probable_pathogen_type", "confidence_score"],
        row,
    )


def record_summary_result(patient_id: str, payload: dict) -> Path:
    """
    Persist summarize_agent output into analytics/data/summarize_agent_results.csv.
    """
    timestamp = datetime.now().isoformat()
    row = {
        "timestamp": timestamp,
        "patient_id": patient_id,
        "final_infection_likelihood": payload.get("final_infection_likelihood", ""),
        "dominant_source": payload.get("dominant_source", ""),
        "dominant_pathogen_type": payload.get("dominant_pathogen_type", ""),
        "overall_confidence": payload.get("overall_confidence", ""),
    }
    return _append_row(
        SUMMARY_DATASET,
        ["timestamp", "patient_id", "final_infection_likelihood", "dominant_source", "dominant_pathogen_type", "overall_confidence"],
        row,
    )


def _extract_patient_id(path: Path) -> str:
    name = path.name
    if name.endswith(".json"):
        name = name[:-5]
    if name.endswith("_final_summary"):
        name = name[: -len("_final_summary")]
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
                            "infection_likelihood": payload.get("infection_likelihood", ""),
                            "probable_source": payload.get("probable_source", ""),
                            "probable_pathogen_type": payload.get("probable_pathogen_type", ""),
                            "confidence_score": payload.get("confidence_score", ""),
                        }
                    )
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with ANALYZE_DATASET.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "timestamp",
                "patient_id",
                "category",
                "infection_likelihood",
                "probable_source",
                "probable_pathogen_type",
                "confidence_score",
            ],
        )
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
            for path in _iter_paths(root, ["*_final_summary.json"]):
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
                        "final_infection_likelihood": payload.get("final_infection_likelihood", ""),
                        "dominant_source": payload.get("dominant_source", ""),
                        "dominant_pathogen_type": payload.get("dominant_pathogen_type", ""),
                        "overall_confidence": payload.get("overall_confidence", ""),
                    }
                )
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with SUMMARY_DATASET.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "timestamp",
                "patient_id",
                "final_infection_likelihood",
                "dominant_source",
                "dominant_pathogen_type",
                "overall_confidence",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return SUMMARY_DATASET


def _load_json(path: Path) -> dict:
    import json

    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)

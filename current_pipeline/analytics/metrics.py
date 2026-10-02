from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, TypedDict


class CultureThreshold(TypedDict, total=False):
    confidence_min: float | None
    confidence_max: float | None
    level_min: int | None
    level_max: int | None


# Placeholders to be filled later.
CULTURE_HC_THRESHOLD: CultureThreshold | None = {
    "confidence_min": 0.9,
    "level_max": 2,
}
CULTURE_LC_THRESHOLD: CultureThreshold | None = {
    "confidence_max": 0.15,
    "level_min": 4,
}
# Example quantile bins: [("q1", None, 100), ("q2", 100, 500), ...]
READ_QUANTILE_BINS: list[tuple[str, float | None, float | None]] = [
    ("q1", None, 20),   # <= P25
    ("q2", 20, 99),     # P25–P50
    ("q3", 99, 770),    # P50–P75
    ("q4", 770, None),  # > P75
]


@dataclass
class Candidate:
    patient_id: str
    organism_name: str
    level: int | None
    reads: float | None
    source: str  # "model" or "culture"
    culture_confidence: str | None = None  # "HC" / "LC" / None

    @property
    def organism_norm(self) -> str:
        return normalize_organism(self.organism_name)


def normalize_organism(name: str) -> str:
    """
    Lowercase, strip whitespace/punctuation, drop parentheses content.
    """
    cleaned = re.sub(r"\([^)]*\)", " ", name)
    cleaned = re.sub(r"[^a-z0-9]+", " ", cleaned.lower())
    return " ".join(cleaned.split())


def _parse_level(text: str | int | None) -> int | None:
    if text is None:
        return None
    if isinstance(text, int):
        return text
    match = re.search(r"(\d+)", str(text))
    return int(match.group(1)) if match else None


def load_mngs_max_agent(path: Path, patient_id: str | None = None) -> list[Candidate]:
    """
    Parse *_mNGS_max_agent.json pathogen candidates into Candidate rows.
    """
    import json

    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    resolved_patient = patient_id or _extract_patient_id(path)
    candidates: list[Candidate] = []
    for entry in data.get("pathogen_candidates", []) or []:
        level = _parse_level(entry.get("integrated_causative_level"))
        reads = entry.get("reads") or entry.get("read_count")
        try:
            reads = float(reads) if reads is not None else None
        except (TypeError, ValueError):
            reads = None
        candidates.append(
            Candidate(
                patient_id=resolved_patient,
                organism_name=entry.get("organism_name", ""),
                level=level,
                reads=reads,
                source="model",
            )
        )
    return candidates


def _extract_patient_id(path: Path) -> str:
    name = path.name
    if name.endswith(".json"):
        name = name[:-5]
    # NGS_patient_<id>_mNGS_max_agent.json
    match = re.search(r"NGS_patient_([^_]+)", name)
    return match.group(1) if match else name


def _names_at_level(candidates: Iterable[Candidate], predicate) -> set[str]:
    return {c.organism_norm for c in candidates if predicate(c)}


def _by_patient(candidates: Iterable[Candidate]) -> dict[str, list[Candidate]]:
    buckets: dict[str, list[Candidate]] = {}
    for c in candidates:
        buckets.setdefault(c.patient_id, []).append(c)
    return buckets


def high_conf_recall(culture_hc: Sequence[Candidate], model_preds: Sequence[Candidate]) -> float:
    """
    Organism-level: fraction of HC culture organisms captured by model (Level<=3).
    """
    gold = list(culture_hc)
    if not gold:
        return math.nan
    hit_count = 0
    preds_by_patient = _by_patient(model_preds)
    for entry in gold:
        preds = preds_by_patient.get(entry.patient_id, [])
        if any(
            p.level is not None and p.level <= 3 and p.organism_norm == entry.organism_norm
            for p in preds
        ):
            hit_count += 1
    return _safe_divide(hit_count, len(gold))


def case_level_capture(culture_hc: Sequence[Candidate], model_preds: Sequence[Candidate]) -> float:
    """
    Patient-level: 1 if any HC organism is captured (Level<=3), else 0; averaged across patients with HC.
    """
    hc_by_patient = _by_patient(culture_hc)
    pred_by_patient = _by_patient(model_preds)
    scores: list[int] = []
    for pid, hc_list in hc_by_patient.items():
        hc_names = {c.organism_norm for c in hc_list}
        if not hc_names:
            continue
        preds = pred_by_patient.get(pid, [])
        hit = any(
            (p.level is not None and p.level <= 3 and p.organism_norm in hc_names) for p in preds
        )
        scores.append(1 if hit else 0)
    return _mean(scores)


def top_k_hit_rate(
    culture_hc: Sequence[Candidate],
    model_preds: Sequence[Candidate],
    k: int,
) -> float:
    """
    Patient-level: 1 if any HC organism falls within top-k levels (ties by level allowed), averaged across patients with HC.
    """
    if k <= 0:
        return math.nan
    hc_by_patient = _by_patient(culture_hc)
    pred_by_patient = _by_patient(model_preds)
    scores: list[int] = []
    for pid, hc_list in hc_by_patient.items():
        hc_names = {c.organism_norm for c in hc_list}
        if not hc_names:
            continue
        preds = pred_by_patient.get(pid, [])
        levels = sorted({c.level for c in preds if c.level is not None})
        if not levels:
            scores.append(0)
            continue
        boundary_level = levels[min(k - 1, len(levels) - 1)]
        top_names = {
            c.organism_norm for c in preds if c.level is not None and c.level <= boundary_level
        }
        scores.append(1 if top_names & hc_names else 0)
    return _mean(scores)


def filter_rate(mngs_candidates: Sequence[Candidate]) -> float:
    """
    Organism-level: proportion of mNGS candidates labeled Level>=4.
    """
    if not mngs_candidates:
        return math.nan
    num = len([c for c in mngs_candidates if c.level is not None and c.level >= 4])
    return _safe_divide(num, len(mngs_candidates))


def deprioritization_rate(culture_lc: Sequence[Candidate], model_preds: Sequence[Candidate]) -> float:
    """
    Organism-level: fraction of LC culture organisms mapped to Level>=4.
    """
    if not culture_lc:
        return math.nan
    pred_by_patient = _by_patient(model_preds)
    hit = 0
    for entry in culture_lc:
        preds = pred_by_patient.get(entry.patient_id, [])
        if any(
            p.level is not None and p.level >= 4 and p.organism_norm == entry.organism_norm
            for p in preds
        ):
            hit += 1
    return _safe_divide(hit, len(culture_lc))


def hc_supported_rate(culture_hc: Sequence[Candidate], model_preds: Sequence[Candidate]) -> float:
    """
    Organism-level: among accepted (Level<=3) predictions, fraction supported by culture HC.
    """
    accepted = [p for p in model_preds if p.level is not None and p.level <= 3]
    if not accepted:
        return math.nan
    hc_by_patient = _by_patient(culture_hc)
    supported = 0
    for pred in accepted:
        hc_names = {c.organism_norm for c in hc_by_patient.get(pred.patient_id, [])}
        if pred.organism_norm in hc_names:
            supported += 1
    return _safe_divide(supported, len(accepted))


def acceptance_rate_by_reads_quantile(
    mngs_candidates: Sequence[Candidate],
    bins: Sequence[tuple[str, float | None, float | None]] | None = None,
) -> dict[str, float]:
    bins = list(bins) if bins is not None else READ_QUANTILE_BINS
    rates: dict[str, float] = {}
    for label, low, high in bins:
        sub = [c for c in mngs_candidates if _in_range(c.reads, low, high)]
        if not sub:
            rates[label] = math.nan
            continue
        accepted = len([c for c in sub if c.level is not None and c.level <= 3])
        rates[label] = _safe_divide(accepted, len(sub))
    return rates


def _in_range(value: float | None, low: float | None, high: float | None) -> bool:
    if value is None:
        return False
    if low is not None and value < low:
        return False
    if high is not None and value >= high:
        return False
    return True


def _safe_divide(numerator: int | float, denominator: int | float) -> float:
    return float(numerator) / float(denominator) if denominator else math.nan


def _mean(values: Sequence[float]) -> float:
    vals = [v for v in values if not math.isnan(v)]
    return sum(vals) / len(vals) if vals else math.nan


def summarize_metrics(
    culture_hc: Sequence[Candidate],
    culture_lc: Sequence[Candidate],
    model_preds: Sequence[Candidate],
    mngs_candidates: Sequence[Candidate] | None = None,
    top_k: int = 3,
    quantile_bins: Sequence[tuple[str, float | None, float | None]] | None = None,
) -> dict[str, float | dict[str, float]]:
    """
    Convenience wrapper returning all seven metrics in one dict.
    """
    mngs_candidates = list(mngs_candidates) if mngs_candidates is not None else list(model_preds)
    return {
        "high_conf_recall": high_conf_recall(culture_hc, model_preds),
        "case_level_capture": case_level_capture(culture_hc, model_preds),
        "top_k_hit_rate": top_k_hit_rate(culture_hc, model_preds, k=top_k),
        "filter_rate": filter_rate(mngs_candidates),
        "deprioritization_rate": deprioritization_rate(culture_lc, model_preds),
        "hc_supported_rate": hc_supported_rate(culture_hc, model_preds),
        "acceptance_rate_by_reads_quantile": acceptance_rate_by_reads_quantile(
            mngs_candidates, bins=quantile_bins
        ),
    }


def load_candidates_from_folder(folder: Path) -> dict[str, list[Candidate]]:
    """
    Utility to load candidates from a patient folder.
    - model: *_mNGS_max_agent.json (Level authority).
    - culture: placeholder lists; caller can filter into HC/LC using thresholds.
    """
    model_candidates: list[Candidate] = []
    culture_candidates: list[Candidate] = []
    for path in folder.glob("*_mNGS_max_agent.json"):
        model_candidates.extend(load_mngs_max_agent(path))
    # Culture loading intentionally left as a placeholder; inject when thresholds are defined.
    return {
        "model": model_candidates,
        "culture": culture_candidates,
    }

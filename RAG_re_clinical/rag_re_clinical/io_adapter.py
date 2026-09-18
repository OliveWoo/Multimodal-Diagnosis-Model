from __future__ import annotations

import copy
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


PathLike = str | Path
JSONDict = dict[str, Any]

MERGE_REVIEW_TIERS = ("review_high_priority", "review_context_needed")
MERGE_CONTEXT_EXCLUDED_KEYS = {
    "pathogen_candidates",
    "excluded_candidates",
    "best_available_summary",
}
REVIEW_CONTEXT_EXCLUDED_KEYS = {
    "review_high_priority",
    "review_context_needed",
    "review_low_specificity",
    "omitted_with_reason",
    "review_omitted_with_reason",
}


class CandidateMappingError(ValueError):
    """Base class for patient/candidate mapping failures."""


class DuplicateCandidateKeyError(CandidateMappingError):
    """Raised when a supposedly unique patient-scoped key occurs twice."""


class AmbiguousCandidateMatchError(CandidateMappingError):
    """Raised instead of choosing between multiple exact/alnum matches."""

    match_scope = "ambiguous"

    def __init__(self, message: str, *, candidate_keys: Sequence[str] = ()) -> None:
        super().__init__(message)
        self.context = {
            "match_scope": self.match_scope,
            "candidate_keys": list(candidate_keys),
        }


class PatientMismatchError(CandidateMappingError):
    """Raised when a frozen artifact and merge JSON belong to different patients."""


@dataclass(frozen=True)
class _ReviewMatch:
    tier: str
    row: JSONDict


@dataclass(frozen=True)
class _MergeIndexes:
    deterministic: dict[str, JSONDict]
    picked: dict[str, JSONDict]
    review: dict[str, _ReviewMatch]


def load_json(path: PathLike) -> JSONDict:
    """Load a UTF-8 JSON object, accepting the BOM used by the frozen merges."""

    json_path = Path(path)
    try:
        with json_path.open("r", encoding="utf-8-sig") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not load JSON object from {json_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"JSON root must be an object: {json_path}")
    return payload


def sha256_file(path: PathLike) -> str:
    """Return a streaming SHA-256 digest without changing the source file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_jsons(source: PathLike | Iterable[PathLike]) -> list[Path]:
    """Discover JSON files deterministically from paths or directories."""

    if isinstance(source, (str, Path)):
        requested: Iterable[PathLike] = (source,)
    else:
        requested = source

    discovered: dict[str, Path] = {}
    for value in requested:
        path = Path(value)
        if not path.exists():
            raise FileNotFoundError(path)
        if path.is_dir():
            candidates = path.rglob("*.json")
        elif path.is_file() and path.suffix.lower() == ".json":
            candidates = (path,)
        elif path.is_file():
            raise ValueError(f"Expected a .json file or directory: {path}")
        else:
            raise ValueError(f"Unsupported JSON source: {path}")
        for candidate in candidates:
            if candidate.is_file():
                absolute = candidate.resolve()
                discovered[str(absolute).casefold()] = absolute
    return sorted(discovered.values(), key=lambda item: str(item).casefold())


def patient_id(payload: Mapping[str, Any], source_path: PathLike | None = None) -> str:
    """Extract one patient ID, rejecting conflicting explicit identifiers."""

    explicit = {
        str(payload[key]).strip()
        for key in ("patient_id", "case_id")
        if payload.get(key) not in (None, "") and str(payload[key]).strip()
    }
    if len(explicit) > 1:
        raise PatientMismatchError(
            f"Payload has conflicting patient/case IDs: {sorted(explicit)}"
        )
    if explicit:
        return next(iter(explicit))

    if source_path is not None:
        match = re.search(
            r"(?:NGS_)?patient[_ -]?(\d+)", Path(source_path).name, re.IGNORECASE
        )
        if match:
            return match.group(1)
    raise ValueError("No patient_id/case_id was present and the filename had no patient ID")


def _validate_merge_payload(payload: Mapping[str, Any], *, source: str = "merge JSON") -> None:
    deterministic = payload.get("deterministic_max")
    review = payload.get("llm_missed_candidate_review")
    if not isinstance(deterministic, dict) or not isinstance(
        deterministic.get("pathogen_candidates"), list
    ):
        raise ValueError(
            f"{source} is not a formal merge: deterministic_max.pathogen_candidates is missing"
        )
    if not isinstance(review, dict):
        raise ValueError(
            f"{source} is not a formal merge: llm_missed_candidate_review is missing"
        )
    for tier in MERGE_REVIEW_TIERS:
        if not isinstance(review.get(tier), list):
            raise ValueError(
                f"{source} is not a formal merge: {tier} must be a list"
            )


def index_merge_files(source: PathLike | Iterable[PathLike]) -> dict[str, Path]:
    """Index formal merge JSON files by patient, rejecting duplicate patients."""

    indexed: dict[str, Path] = {}
    for path in discover_jsons(source):
        payload = load_json(path)
        _validate_merge_payload(payload, source=str(path))
        pid = patient_id(payload, path)
        if pid in indexed:
            raise PatientMismatchError(
                f"Duplicate merge JSONs for patient {pid}: {indexed[pid]} and {path}"
            )
        indexed[pid] = path
    return dict(sorted(indexed.items(), key=lambda item: item[0]))


def _candidate_name(value: Mapping[str, Any]) -> str:
    for field in (
        "organism_name",
        "pathogen",
        "pathogen_name",
        "name",
        "organism",
        "canonical_organism_name",
    ):
        text = " ".join(str(value.get(field) or "").strip().split())
        if text:
            return text
    return ""


def _alnum_key(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").casefold())


def _row_key(row: Mapping[str, Any], *, layer: str) -> str:
    name = _candidate_name(row)
    key = _alnum_key(name)
    if not key:
        raise CandidateMappingError(f"{layer} row has no usable organism name: {row!r}")
    return key


def _index_unique_rows(rows: Iterable[Any], *, layer: str) -> dict[str, JSONDict]:
    indexed: dict[str, JSONDict] = {}
    for position, raw_row in enumerate(rows):
        if not isinstance(raw_row, dict):
            raise CandidateMappingError(f"{layer}[{position}] must be an object")
        key = _row_key(raw_row, layer=f"{layer}[{position}]")
        if key in indexed:
            raise DuplicateCandidateKeyError(
                f"Duplicate patient-scoped candidate key {key!r} in {layer}"
            )
        indexed[key] = raw_row
    return indexed


def _picked_rows(deterministic: Mapping[str, Any]) -> list[Any]:
    summary = deterministic.get("best_available_summary")
    if isinstance(summary, dict) and isinstance(summary.get("picked_pathogens"), list):
        return list(summary["picked_pathogens"])
    values = deterministic.get("picked_pathogens")
    return list(values) if isinstance(values, list) else []


def _build_merge_indexes(merge_payload: Mapping[str, Any]) -> _MergeIndexes:
    _validate_merge_payload(merge_payload)
    deterministic_bundle = merge_payload["deterministic_max"]
    deterministic = _index_unique_rows(
        deterministic_bundle["pathogen_candidates"],
        layer="deterministic_max.pathogen_candidates",
    )
    picked = _index_unique_rows(
        _picked_rows(deterministic_bundle),
        layer="deterministic_max.best_available_summary.picked_pathogens",
    )
    missing_picked = sorted(set(picked) - set(deterministic))
    if missing_picked:
        raise CandidateMappingError(
            "Formal picked rows must map to deterministic candidates; missing keys: "
            + ", ".join(missing_picked)
        )

    review_payload = merge_payload["llm_missed_candidate_review"]
    review: dict[str, _ReviewMatch] = {}
    for tier in MERGE_REVIEW_TIERS:
        rows = _index_unique_rows(
            review_payload[tier], layer=f"llm_missed_candidate_review.{tier}"
        )
        overlap = sorted(set(review) & set(rows))
        if overlap:
            raise AmbiguousCandidateMatchError(
                "Review candidates occur in more than one active tier: "
                + ", ".join(overlap),
                candidate_keys=overlap,
            )
        review.update(
            {key: _ReviewMatch(tier=tier, row=row) for key, row in rows.items()}
        )
    return _MergeIndexes(deterministic=deterministic, picked=picked, review=review)


def _frozen_lookup_keys(candidate: Mapping[str, Any]) -> list[str]:
    names: list[str] = []
    for field in ("organism_name", "canonical_organism_name"):
        name = " ".join(str(candidate.get(field) or "").strip().split())
        if name and name not in names:
            names.append(name)
    if not names:
        fallback = _candidate_name(candidate)
        if fallback:
            names.append(fallback)
    keys = list(dict.fromkeys(_alnum_key(name) for name in names if _alnum_key(name)))
    if not keys:
        raise CandidateMappingError(f"Frozen candidate has no usable name: {candidate!r}")
    return keys


def _choose_match_key(
    candidate: Mapping[str, Any], indexes: _MergeIndexes
) -> tuple[str, str]:
    lookup_keys = _frozen_lookup_keys(candidate)
    deterministic_matches = [key for key in lookup_keys if key in indexes.deterministic]
    if len(deterministic_matches) > 1:
        raise AmbiguousCandidateMatchError(
            "Frozen candidate names match multiple deterministic candidates: "
            + ", ".join(deterministic_matches),
            candidate_keys=deterministic_matches,
        )
    if deterministic_matches:
        return deterministic_matches[0], "exact_species"

    metadata_matches = [
        key for key in lookup_keys if key in indexes.review or key in indexes.picked
    ]
    if len(metadata_matches) > 1:
        raise AmbiguousCandidateMatchError(
            "Frozen candidate names match multiple merge metadata rows: "
            + ", ".join(metadata_matches),
            candidate_keys=metadata_matches,
        )
    return (
        metadata_matches[0] if metadata_matches else lookup_keys[0],
        "synthetic_unmatched",
    )


def _synthetic_candidate(
    frozen_candidate: Mapping[str, Any], review_match: _ReviewMatch | None
) -> JSONDict:
    row = review_match.row if review_match else {}
    snapshot = row.get("evidence_snapshot")
    snapshot = snapshot if isinstance(snapshot, dict) else {}
    effective: JSONDict = {
        "organism_name": _candidate_name(row) or _candidate_name(frozen_candidate),
    }
    for field in (
        "classification",
        "integrated_causative_level",
        "observed_level",
        "source_category",
    ):
        if row.get(field) not in (None, ""):
            effective[field] = copy.deepcopy(row[field])
    for field, value in snapshot.items():
        if value not in (None, ""):
            effective[field] = copy.deepcopy(value)
    if review_match:
        effective["rag_re_review_tier"] = review_match.tier
    effective["rag_re_candidate_provenance"] = "merge_review_snapshot_fallback"
    return {
        "source": "review_evidence_snapshot" if review_match else "frozen_candidate_only",
        "candidate_fields": effective,
        "evidence_snapshot": copy.deepcopy(snapshot),
    }


def _merge_one(
    frozen_candidate: Mapping[str, Any],
    indexes: _MergeIndexes,
    *,
    candidate_index: int | None = None,
) -> JSONDict:
    key, match_scope = _choose_match_key(frozen_candidate, indexes)
    deterministic = indexes.deterministic.get(key)
    picked = indexes.picked.get(key)
    review_match = indexes.review.get(key)
    review_metadata = (
        {
            "tier": review_match.tier,
            "review_tier": review_match.tier,
            "row": copy.deepcopy(review_match.row),
        }
        if review_match
        else None
    )
    synthetic = (
        _synthetic_candidate(frozen_candidate, review_match)
        if match_scope == "synthetic_unmatched"
        else None
    )
    effective = (
        copy.deepcopy(deterministic)
        if deterministic is not None
        else copy.deepcopy(synthetic["candidate_fields"] if synthetic else {})
    )
    source_layers = {
        "frozen_candidate": copy.deepcopy(dict(frozen_candidate)),
        "deterministic_candidate": copy.deepcopy(deterministic),
        "picked_metadata": copy.deepcopy(picked),
        "review_metadata": copy.deepcopy(review_metadata),
        "synthetic_fallback": copy.deepcopy(synthetic),
    }
    return {
        "candidate_index": candidate_index,
        "organism_name": _candidate_name(frozen_candidate),
        "canonical_organism_name": frozen_candidate.get("canonical_organism_name"),
        "candidate_key": key,
        "match_scope": match_scope,
        "frozen_candidate": copy.deepcopy(dict(frozen_candidate)),
        "deterministic_candidate": copy.deepcopy(deterministic),
        "picked_metadata": copy.deepcopy(picked),
        "review_metadata": review_metadata,
        "synthetic_fallback": synthetic,
        "effective_candidate": effective,
        "source_layers": source_layers,
        "provenance": {
            "deterministic_exact": deterministic is not None,
            "formal_picked": picked is not None,
            "review_tier": review_match.tier if review_match else None,
            "fallback_source": synthetic.get("source") if synthetic else None,
        },
    }


def _merge_patient_context(merge_payload: Mapping[str, Any]) -> JSONDict:
    deterministic = merge_payload["deterministic_max"]
    summary = deterministic.get("best_available_summary")
    review = merge_payload["llm_missed_candidate_review"]
    return {
        "root_metadata": copy.deepcopy(
            {
                key: value
                for key, value in merge_payload.items()
                if key not in {"deterministic_max", "llm_missed_candidate_review"}
            }
        ),
        "deterministic_patient_context": copy.deepcopy(
            {
                key: value
                for key, value in deterministic.items()
                if key not in MERGE_CONTEXT_EXCLUDED_KEYS
            }
        ),
        "best_available_summary_metadata": copy.deepcopy(
            {
                key: value
                for key, value in (summary.items() if isinstance(summary, dict) else ())
                if key != "picked_pathogens"
            }
        ),
        "review_policy_metadata": copy.deepcopy(
            {
                key: value
                for key, value in review.items()
                if key not in REVIEW_CONTEXT_EXCLUDED_KEYS
            }
        ),
    }


def _coerce_payload(
    value: Mapping[str, Any] | PathLike, explicit_path: PathLike | None
) -> tuple[JSONDict, Path | None]:
    if isinstance(value, Mapping):
        return dict(value), Path(explicit_path) if explicit_path is not None else None
    path = Path(value)
    if explicit_path is not None and path.resolve() != Path(explicit_path).resolve():
        raise ValueError(f"Conflicting source paths: {path} and {explicit_path}")
    return load_json(path), path


def merge_one_candidate_context(
    frozen_candidate: Mapping[str, Any],
    merge_payload: Mapping[str, Any] | PathLike,
    *,
    candidate_index: int | None = None,
    merge_path: PathLike | None = None,
) -> JSONDict:
    """Merge one frozen candidate without genus or fuzzy evidence transfer."""

    merge, _ = _coerce_payload(merge_payload, merge_path)
    return _merge_one(
        frozen_candidate,
        _build_merge_indexes(merge),
        candidate_index=candidate_index,
    )


def merge_candidate_context(
    frozen_artifact: Mapping[str, Any] | PathLike,
    merge_payload: Mapping[str, Any] | PathLike,
    *,
    frozen_path: PathLike | None = None,
    merge_path: PathLike | None = None,
) -> JSONDict:
    """Layer frozen, deterministic, picked, and review data for one patient.

    The function is deliberately label-blind: it accepts only a frozen RAG_re
    artifact and its formal merge JSON.  It performs patient-scoped exact/alnum
    matching and never reads clinical gold or applies genus/fuzzy matching.
    """

    frozen, frozen_source = _coerce_payload(frozen_artifact, frozen_path)
    merge, merge_source = _coerce_payload(merge_payload, merge_path)
    _validate_merge_payload(merge)
    frozen_pid = patient_id(frozen, frozen_source)
    merge_pid = patient_id(merge, merge_source)
    if frozen_pid != merge_pid:
        raise PatientMismatchError(
            f"Frozen artifact patient {frozen_pid} does not match merge patient {merge_pid}"
        )
    candidates = frozen.get("candidates")
    if not isinstance(candidates, list):
        raise ValueError("Frozen RAG_re artifact must contain a candidates list")

    indexes = _build_merge_indexes(merge)
    merged_rows: list[JSONDict] = []
    seen_candidate_keys: set[str] = set()
    for index, candidate in enumerate(candidates):
        if not isinstance(candidate, dict):
            raise CandidateMappingError(f"frozen candidates[{index}] must be an object")
        row = _merge_one(candidate, indexes, candidate_index=index)
        key = str(row["candidate_key"])
        if key in seen_candidate_keys:
            raise DuplicateCandidateKeyError(
                f"Frozen artifact has duplicate patient-scoped candidate key {key!r}"
            )
        seen_candidate_keys.add(key)
        merged_rows.append(row)

    def source_metadata(path: Path | None) -> JSONDict:
        return {
            "path": str(path) if path is not None else None,
            "sha256": sha256_file(path) if path is not None else None,
        }

    frozen_source_metadata = source_metadata(frozen_source)
    merge_source_metadata = source_metadata(merge_source)
    return {
        "schema_version": "rag_re_clinical.merged_context.v1",
        "patient_id": frozen_pid,
        "sources": {
            "frozen_rag_re": frozen_source_metadata,
            "formal_merge": merge_source_metadata,
        },
        # Flat aliases keep the engine-facing envelope convenient while the
        # structured sources object remains the canonical provenance record.
        "frozen_artifact_sha256": frozen_source_metadata["sha256"],
        "merge_sha256": merge_source_metadata["sha256"],
        "matching_policy": {
            "patient_scoped": True,
            "name_matching": "exact_alnum_only",
            "genus_matching": False,
            "fuzzy_matching": False,
            "ambiguous_match_action": "raise",
        },
        "merge_context": _merge_patient_context(merge),
        "candidates": merged_rows,
        "counts": {
            "candidate_count": len(merged_rows),
            "exact_species": sum(
                row["match_scope"] == "exact_species" for row in merged_rows
            ),
            "synthetic_unmatched": sum(
                row["match_scope"] == "synthetic_unmatched" for row in merged_rows
            ),
            "ambiguous": 0,
        },
    }

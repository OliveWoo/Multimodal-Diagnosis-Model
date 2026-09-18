from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable, Mapping


PathLike = str | Path
JSONDict = dict[str, Any]

SOURCE_SCHEMA_VERSION = "rag_re_clinical.output.v1"
SOURCE_MANIFEST_SCHEMA_VERSION = "rag_re_clinical.batch_manifest.v1"

ALLOWED_D_STATES = {"D_PASS", "D_GUARDED", "D_BLOCK", "D_UNKNOWN"}
ALLOWED_C_GRADES = {"C0", "C1", "C2", "C3", "CNEG"}
ALLOWED_F_STATES = {"F_COHERENT", "F_MISMATCH", "F_UNKNOWN"}
ALLOWED_REVIEW_TIERS = {None, "review_high_priority", "review_context_needed"}


class SourceContractError(ValueError):
    """Raised when a locked clinical artifact violates the source contract."""


def load_json(path: PathLike) -> JSONDict:
    source = Path(path)
    try:
        with source.open("r", encoding="utf-8-sig") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise SourceContractError(f"Could not read JSON object from {source}: {exc}") from exc
    if not isinstance(value, dict):
        raise SourceContractError(f"JSON root must be an object: {source}")
    return value


def sha256_file(path: PathLike) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_canonical_json(value: Mapping[str, Any]) -> str:
    try:
        encoded = json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise SourceContractError(
            f"In-memory source must be JSON-serializable: {exc}"
        ) from exc
    return hashlib.sha256(encoded).hexdigest()


def discover_source_jsons(source: PathLike | Iterable[PathLike]) -> list[Path]:
    requested = (source,) if isinstance(source, (str, Path)) else source
    found: dict[str, Path] = {}
    for raw in requested:
        path = Path(raw)
        if not path.exists():
            raise FileNotFoundError(path)
        candidates = path.glob("*.json") if path.is_dir() else (path,)
        for candidate in candidates:
            if candidate.suffix.lower() != ".json" or not candidate.is_file():
                if path.is_file():
                    raise SourceContractError(
                        f"Expected a .json file or directory: {path}"
                    )
                continue
            absolute = candidate.resolve()
            found[str(absolute).casefold()] = absolute
    return sorted(found.values(), key=lambda item: str(item).casefold())


def _strict_bool_or_none(value: Any, *, where: str) -> bool | None:
    if value is True or value is False or value is None:
        return value
    raise SourceContractError(f"{where} must be true, false, or null")


def _required_mapping(value: Any, *, where: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise SourceContractError(f"{where} must be an object")
    return value


def _required_text(value: Any, *, where: str) -> str:
    text = " ".join(str(value or "").strip().split())
    if not text:
        raise SourceContractError(f"{where} must be a non-empty string")
    return text


def _name_key(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip()).casefold()


def _project_candidate(candidate: Mapping[str, Any], *, index: int) -> JSONDict:
    where = f"candidates[{index}]"
    organism_name = _required_text(
        candidate.get("organism_name"), where=f"{where}.organism_name"
    )
    canonical_name = _required_text(
        candidate.get("canonical_organism_name") or organism_name,
        where=f"{where}.canonical_organism_name",
    )
    baseline = candidate.get("baseline_selected")
    if not isinstance(baseline, bool):
        raise SourceContractError(f"{where}.baseline_selected must be boolean")

    modules = _required_mapping(candidate.get("modules"), where=f"{where}.modules")
    d = _required_mapping(modules.get("D"), where=f"{where}.modules.D")
    b = _required_mapping(
        modules.get("B_STRICT"), where=f"{where}.modules.B_STRICT"
    )
    c = _required_mapping(modules.get("C"), where=f"{where}.modules.C")
    e = _required_mapping(modules.get("E"), where=f"{where}.modules.E")
    f = _required_mapping(modules.get("F"), where=f"{where}.modules.F")

    d_state = _required_text(d.get("state"), where=f"{where}.modules.D.state")
    if d_state not in ALLOWED_D_STATES:
        raise SourceContractError(f"Unsupported D state at {where}: {d_state}")
    c_grade = _required_text(c.get("grade"), where=f"{where}.modules.C.grade")
    if c_grade not in ALLOWED_C_GRADES:
        raise SourceContractError(f"Unsupported C grade at {where}: {c_grade}")
    f_state = _required_text(f.get("state"), where=f"{where}.modules.F.state")
    if f_state not in ALLOWED_F_STATES:
        raise SourceContractError(f"Unsupported F state at {where}: {f_state}")
    e_inputs_raw = e.get("inputs")
    if e_inputs_raw is None:
        e_inputs: Mapping[str, Any] = {}
    elif isinstance(e_inputs_raw, Mapping):
        e_inputs = e_inputs_raw
    else:
        raise SourceContractError(f"{where}.modules.E.inputs must be an object or null")
    review_tier = candidate.get("review_tier")
    if review_tier is not None:
        review_tier = _required_text(review_tier, where=f"{where}.review_tier")
    if review_tier not in ALLOWED_REVIEW_TIERS:
        raise SourceContractError(
            f"Unsupported review_tier at {where}: {review_tier!r}"
        )

    # This is an intentional allowlist.  Free-text rationale, answer-like keys,
    # source diagnoses and all other fields cannot influence pruning decisions.
    return {
        "organism_name": organism_name,
        "canonical_organism_name": canonical_name,
        "baseline_selected": baseline,
        "review_tier": review_tier,
        "source_features": {
            "D": {
                "state": d_state,
            },
            "B_STRICT": {
                "positive": _strict_bool_or_none(
                    b.get("positive"), where=f"{where}.modules.B_STRICT.positive"
                ),
                "absolute_strength_axis": _strict_bool_or_none(
                    b.get("absolute_strength_axis"),
                    where=f"{where}.modules.B_STRICT.absolute_strength_axis",
                ),
                "relative_strength_axis": _strict_bool_or_none(
                    b.get("relative_strength_axis"),
                    where=f"{where}.modules.B_STRICT.relative_strength_axis",
                ),
            },
            "C": {
                "grade": c_grade,
            },
            "E": {
                "A_positive": _strict_bool_or_none(
                    e_inputs.get("A"), where=f"{where}.modules.E.inputs.A"
                ),
            },
            "F": {
                "state": f_state,
            },
        },
    }


def validate_source_artifact(
    payload: Mapping[str, Any], *, source_path: PathLike | None = None
) -> JSONDict:
    if payload.get("schema_version") != SOURCE_SCHEMA_VERSION:
        raise SourceContractError(
            f"Expected {SOURCE_SCHEMA_VERSION}, got {payload.get('schema_version')!r}"
        )
    if payload.get("run_status") != "complete":
        raise SourceContractError("Clinical source must have run_status='complete'")
    patient_id = _required_text(payload.get("patient_id"), where="patient_id")
    raw_candidates = payload.get("candidates")
    if not isinstance(raw_candidates, list) or not raw_candidates:
        raise SourceContractError("candidates must be a non-empty list")
    candidates = [
        _project_candidate(candidate, index=index)
        if isinstance(candidate, Mapping)
        else (_ for _ in ()).throw(
            SourceContractError(f"candidates[{index}] must be an object")
        )
        for index, candidate in enumerate(raw_candidates)
    ]
    name_keys = [_name_key(row["organism_name"]) for row in candidates]
    if len(set(name_keys)) != len(name_keys):
        raise SourceContractError(f"Patient {patient_id} has duplicate organism names")

    arms = _required_mapping(payload.get("arms"), where="arms")
    baseline_arm = _required_mapping(arms.get("CL0_BASELINE"), where="arms.CL0_BASELINE")
    raw_predictions = baseline_arm.get("predicted_pathogens")
    if not isinstance(raw_predictions, list):
        raise SourceContractError("arms.CL0_BASELINE.predicted_pathogens must be a list")
    baseline_names = [
        row["organism_name"] for row in candidates if row["baseline_selected"]
    ]
    if raw_predictions != baseline_names:
        raise SourceContractError(
            f"Patient {patient_id} CL0_BASELINE must exactly match baseline_selected order"
        )
    if baseline_arm.get("predicted_count") != len(baseline_names):
        raise SourceContractError(
            f"Patient {patient_id} CL0_BASELINE predicted_count is inconsistent"
        )

    counts = {
        "candidate_count": len(candidates),
        "e0_count": len(baseline_names),
        "non_e0_count": len(candidates) - len(baseline_names),
    }
    source = Path(source_path) if source_path is not None else None
    source_hash = (
        sha256_file(source)
        if source is not None
        else _sha256_canonical_json(
            {
                "patient_id": patient_id,
                "candidates": candidates,
                "counts": counts,
            }
        )
    )
    return {
        "patient_id": patient_id,
        "source": {
            "path": str(source.resolve()) if source is not None else "<memory>",
            "sha256": source_hash,
            "sha256_basis": (
                "source_file_bytes_v1"
                if source is not None
                else "canonical_projected_inputs_v1"
            ),
            "schema_version": SOURCE_SCHEMA_VERSION,
            "pipeline_version": payload.get("pipeline_version"),
            "artifact_sha256_basis": payload.get("artifact_sha256_basis"),
        },
        "candidates": candidates,
        "counts": counts,
    }


def load_source_artifact(path: PathLike) -> JSONDict:
    return validate_source_artifact(load_json(path), source_path=path)


def index_source_directory(source: PathLike) -> dict[str, tuple[Path, JSONDict]]:
    indexed: dict[str, tuple[Path, JSONDict]] = {}
    for path in discover_source_jsons(source):
        raw = load_json(path)
        if raw.get("schema_version") == SOURCE_MANIFEST_SCHEMA_VERSION:
            continue
        projected = validate_source_artifact(raw, source_path=path)
        patient_id = projected["patient_id"]
        if patient_id in indexed:
            raise SourceContractError(
                f"Duplicate clinical artifacts for patient {patient_id}: "
                f"{indexed[patient_id][0]} and {path}"
            )
        indexed[patient_id] = (path, projected)
    if not indexed:
        raise SourceContractError(f"No {SOURCE_SCHEMA_VERSION} artifacts found in {source}")
    return dict(
        sorted(indexed.items(), key=lambda item: (len(item[0]), item[0]))
    )


__all__ = [
    "SOURCE_SCHEMA_VERSION",
    "SourceContractError",
    "discover_source_jsons",
    "index_source_directory",
    "load_json",
    "load_source_artifact",
    "sha256_file",
    "validate_source_artifact",
]

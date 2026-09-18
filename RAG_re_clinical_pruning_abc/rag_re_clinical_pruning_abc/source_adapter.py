"""Strict, label-blind adapter joining frozen L5, clinical, and raw A/B/C artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from .states import C_GRADES, D_STATES, E0_STATES, NON_E0_STATES

PRUNING_SCHEMA = "rag_re_clinical_pruning.output.v1"
PRUNING_MANIFEST_SCHEMA = "rag_re_clinical_pruning.batch_manifest.v1"
CLINICAL_SCHEMA = "rag_re_clinical.output.v1"
RAW_SCHEMA = "rag_re.output.v1"
L5_SOURCE_ARM = "L5_PLUS_DIRECT_CONVERGENT_RESCUE"


class SourceValidationError(ValueError):
    """Raised when frozen source identity or schema invariants do not hold."""


def load_json(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise SourceValidationError(f"Cannot read JSON {source}: {exc}") from exc
    if not isinstance(value, dict):
        raise SourceValidationError(f"JSON root must be an object: {source}")
    return value


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_pruning_manifest_anchors(
    pruning_directory: str | Path,
    pruning_index: dict[str, Path],
    clinical_index: dict[str, Path],
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """Validate the locked pruning manifest and return per-patient anchors."""

    manifest_path = Path(pruning_directory) / "manifest.json"
    manifest = load_json(manifest_path)
    if manifest.get("schema_version") != PRUNING_MANIFEST_SCHEMA:
        raise SourceValidationError(
            f"Pruning manifest schema must be {PRUNING_MANIFEST_SCHEMA}"
        )
    if manifest.get("run_status") != "complete":
        raise SourceValidationError("Pruning manifest must have run_status=complete")
    if manifest.get("output_schema_version") != PRUNING_SCHEMA:
        raise SourceValidationError("Pruning manifest output schema mismatch")
    if manifest.get("primary_arm") != L5_SOURCE_ARM:
        raise SourceValidationError("Pruning manifest does not anchor frozen L5")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise SourceValidationError("Pruning manifest artifacts must be a list")
    if manifest.get("patient_count") != len(pruning_index) or len(artifacts) != len(pruning_index):
        raise SourceValidationError("Pruning manifest patient/artifact counts do not reconcile")
    if set(pruning_index) != set(clinical_index):
        raise SourceValidationError("Pruning and clinical patient sets differ before anchoring")

    manifest_sha = sha256_file(manifest_path)
    anchors: dict[str, dict[str, Any]] = {}
    for offset, entry in enumerate(artifacts):
        if not isinstance(entry, dict):
            raise SourceValidationError(f"Pruning manifest artifact {offset} must be an object")
        patient_id = _patient_id(entry.get("patient_id"), f"manifest.artifacts[{offset}]")
        if patient_id in anchors:
            raise SourceValidationError(f"Duplicate patient {patient_id} in pruning manifest")
        if patient_id not in pruning_index:
            raise SourceValidationError(f"Manifest patient {patient_id} has no supplied pruning file")
        pruning_path = pruning_index[patient_id].resolve()
        clinical_path = clinical_index[patient_id].resolve()
        recorded_pruning_path = entry.get("path")
        recorded_clinical_path = entry.get("source_path")
        if not isinstance(recorded_pruning_path, str) or Path(recorded_pruning_path).resolve() != pruning_path:
            raise SourceValidationError(f"Manifest pruning path mismatch for patient {patient_id}")
        if not isinstance(recorded_clinical_path, str) or Path(recorded_clinical_path).resolve() != clinical_path:
            raise SourceValidationError(f"Manifest clinical path mismatch for patient {patient_id}")
        pruning_sha = sha256_file(pruning_path)
        clinical_sha = sha256_file(clinical_path)
        if entry.get("sha256") != pruning_sha:
            raise SourceValidationError(f"Manifest pruning SHA mismatch for patient {patient_id}")
        if entry.get("source_sha256") != clinical_sha:
            raise SourceValidationError(f"Manifest clinical SHA mismatch for patient {patient_id}")
        pruning_artifact = load_json(pruning_path)
        if entry.get("decision_sha256") != pruning_artifact.get("decision_sha256"):
            raise SourceValidationError(f"Manifest decision SHA mismatch for patient {patient_id}")
        counts = pruning_artifact.get("counts")
        if not isinstance(counts, dict) or any(
            entry.get(key) != counts.get(key)
            for key in ("candidate_count", "e0_count", "non_e0_count")
        ):
            raise SourceValidationError(f"Manifest candidate counts mismatch for patient {patient_id}")
        anchors[patient_id] = {
            "manifest_path": str(manifest_path.resolve()),
            "manifest_sha256": manifest_sha,
            "pruning_sha256": pruning_sha,
            "pruning_decision_sha256": entry.get("decision_sha256"),
            "clinical_sha256": clinical_sha,
        }
    if set(anchors) != set(pruning_index):
        raise SourceValidationError("Pruning manifest and supplied artifact patient sets differ")
    return anchors, {
        "path": str(manifest_path.resolve()),
        "sha256": manifest_sha,
        "schema_version": PRUNING_MANIFEST_SCHEMA,
    }


def _patient_id(value: Any, label: str) -> str:
    if value is None or isinstance(value, bool):
        raise SourceValidationError(f"{label}.patient_id is missing or invalid")
    text = str(value).strip()
    if not text:
        raise SourceValidationError(f"{label}.patient_id is empty")
    return text


def _list(value: Any, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise SourceValidationError(f"{label} must be a list")
    return value


def _dict(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise SourceValidationError(f"{label} must be an object")
    return value


def _tri(value: Any, label: str) -> bool | None:
    if value is None or type(value) is bool:
        return value
    raise SourceValidationError(f"{label} must be true, false, or null")


def _canonical_name(candidate: dict[str, Any], label: str) -> str:
    value = candidate.get("canonical_organism_name")
    if not isinstance(value, str) or not value.strip():
        raise SourceValidationError(f"{label}.canonical_organism_name is required")
    return value.strip()


def _validate_unique_names(candidates: Iterable[dict[str, Any]], label: str) -> None:
    seen: dict[str, str] = {}
    for index, candidate in enumerate(candidates):
        name = _canonical_name(candidate, f"{label}[{index}]")
        key = name.casefold()
        if key in seen:
            raise SourceValidationError(
                f"Duplicate canonical candidate in {label}: {seen[key]!r} and {name!r}"
            )
        seen[key] = name


def _validate_raw_a(module: dict[str, Any], label: str) -> dict[str, Any]:
    positive = _tri(module.get("positive"), f"{label}.positive")
    status = module.get("status")
    if status not in {"ok", "not_requested", "insufficient_articles"}:
        raise SourceValidationError(f"{label}.status is invalid: {status!r}")
    support = module.get("support_count")
    judgeable = module.get("judgeable_count")
    ratio = module.get("support_ratio")
    if type(support) is not int or type(judgeable) is not int:
        raise SourceValidationError(f"{label} counts must be integers")
    if not (0 <= support <= judgeable <= 10):
        raise SourceValidationError(f"{label} counts violate 0<=support<=judgeable<=10")
    if ratio is not None and (
        isinstance(ratio, bool) or not isinstance(ratio, (int, float)) or not 0 <= ratio <= 1
    ):
        raise SourceValidationError(f"{label}.support_ratio must be in [0,1] or null")
    expected_ratio = support / judgeable if judgeable else None
    if ratio is None and expected_ratio is not None:
        raise SourceValidationError(f"{label}.support_ratio is missing despite judgeable articles")
    if ratio is not None and expected_ratio is None:
        raise SourceValidationError(f"{label}.support_ratio must be null when judgeable_count=0")
    if ratio is not None and abs(float(ratio) - expected_ratio) > 1e-6:
        raise SourceValidationError(f"{label}.support_ratio does not equal support/judgeable")
    if status == "not_requested":
        if positive is not None or support != 0 or judgeable != 0 or ratio is not None:
            raise SourceValidationError(f"{label} violates not_requested invariants")
    elif status == "insufficient_articles":
        if positive is not None or judgeable >= 5:
            raise SourceValidationError(f"{label} violates insufficient_articles invariants")
    elif positive is None or judgeable < 5:
        raise SourceValidationError(f"{label} violates ok invariants")
    return {
        "positive": positive,
        "status": status,
        "support_count": support,
        "judgeable_count": judgeable,
        "support_ratio": ratio,
    }


def _validate_raw_b(module: dict[str, Any], label: str) -> dict[str, Any]:
    positive = _tri(module.get("positive"), f"{label}.positive")
    status = module.get("status")
    if status != "ok" or positive is None:
        raise SourceValidationError(f"{label} must be an observed boolean with status=ok")
    site_aligned = _tri(module.get("site_aligned"), f"{label}.site_aligned")
    triggers = module.get("triggers")
    if not isinstance(triggers, list) or not all(isinstance(item, str) for item in triggers):
        raise SourceValidationError(f"{label}.triggers must be a string list")
    if positive is True and (site_aligned is not True or not triggers):
        raise SourceValidationError(f"{label} positive requires aligned site and trigger(s)")
    if positive is False and triggers:
        raise SourceValidationError(f"{label} false cannot retain positive triggers")
    return {
        "positive": positive,
        "status": status,
        "site_aligned": site_aligned,
        "trigger_count": len(triggers),
    }


def _validate_raw_c(module: dict[str, Any], label: str) -> dict[str, Any]:
    positive = _tri(module.get("positive"), f"{label}.positive")
    status = module.get("status")
    if status not in {"ok", "insufficient_input"}:
        raise SourceValidationError(f"{label}.status is invalid: {status!r}")
    site_aligned = _tri(module.get("site_aligned"), f"{label}.site_aligned")
    supports = module.get("supports")
    if not isinstance(supports, list) or not all(isinstance(item, dict) for item in supports):
        raise SourceValidationError(f"{label}.supports must be an object list")
    if positive is True and (status != "ok" or site_aligned is not True or not supports):
        raise SourceValidationError(f"{label} violates positive support/site invariants")
    if positive is None and (
        status != "insufficient_input" or site_aligned is not None or not supports
    ):
        raise SourceValidationError(f"{label} violates abstention invariants")
    if positive is False and status != "ok":
        raise SourceValidationError(f"{label} false result must have status=ok")
    if positive is False and supports and site_aligned is not False:
        raise SourceValidationError(
            f"{label} support-bearing false requires an explicit site mismatch"
        )
    return {
        "positive": positive,
        "status": status,
        "site_aligned": site_aligned,
        "support_record_count": len(supports),
    }


def _project_candidate(
    pruning_candidate: dict[str, Any],
    clinical_candidate: dict[str, Any],
    raw_candidate: dict[str, Any],
    index: int,
) -> dict[str, Any]:
    label = f"candidates[{index}]"
    names = (
        _canonical_name(pruning_candidate, f"pruning.{label}"),
        _canonical_name(clinical_candidate, f"clinical.{label}"),
        _canonical_name(raw_candidate, f"raw.{label}"),
    )
    if not names[0] == names[1] == names[2]:
        raise SourceValidationError(
            f"Exact candidate order/name mismatch at {index}: {names!r}"
        )
    baseline = pruning_candidate.get("baseline_selected")
    if type(baseline) is not bool:
        raise SourceValidationError(f"pruning.{label}.baseline_selected must be boolean")
    for source_name, candidate in (("clinical", clinical_candidate), ("raw", raw_candidate)):
        if candidate.get("baseline_selected") is not baseline:
            raise SourceValidationError(f"{source_name}.{label}.baseline_selected mismatch")

    input_features = _dict(pruning_candidate.get("input_features"), f"pruning.{label}.input_features")
    d_state = _dict(input_features.get("D"), f"pruning.{label}.input_features.D").get("state")
    if d_state not in D_STATES:
        raise SourceValidationError(f"pruning.{label} has invalid D state {d_state!r}")
    f_state = _dict(input_features.get("F"), f"pruning.{label}.input_features.F").get("state")
    if f_state not in {"F_COHERENT", "F_MISMATCH", "F_UNKNOWN"}:
        raise SourceValidationError(f"pruning.{label} has invalid F state {f_state!r}")
    corrected_a = _tri(
        _dict(input_features.get("E"), f"pruning.{label}.input_features.E").get("A_positive"),
        f"pruning.{label}.input_features.E.A_positive",
    )
    corrected_b = _tri(
        _dict(input_features.get("B_STRICT"), f"pruning.{label}.input_features.B_STRICT").get("positive"),
        f"pruning.{label}.input_features.B_STRICT.positive",
    )
    c_grade = _dict(input_features.get("C"), f"pruning.{label}.input_features.C").get("grade")
    if c_grade not in C_GRADES:
        raise SourceValidationError(f"pruning.{label} has invalid C grade {c_grade!r}")
    corrected_c = True if c_grade in {"C2", "C3"} else False if c_grade == "CNEG" else None

    # The pruning projection is only an optimization artifact. Reconcile every
    # decision-bearing corrected field against its locked clinical source so a
    # locally edited pruning JSON cannot silently redefine CORR A/B/C or D/F.
    clinical_modules = _dict(clinical_candidate.get("modules"), f"clinical.{label}.modules")
    clinical_d = _dict(clinical_modules.get("D"), f"clinical.{label}.modules.D").get("state")
    clinical_b = _tri(
        _dict(clinical_modules.get("B_STRICT"), f"clinical.{label}.modules.B_STRICT").get("positive"),
        f"clinical.{label}.modules.B_STRICT.positive",
    )
    clinical_c_grade = _dict(
        clinical_modules.get("C"), f"clinical.{label}.modules.C"
    ).get("grade")
    clinical_e_inputs = _dict(
        _dict(clinical_modules.get("E"), f"clinical.{label}.modules.E").get("inputs"),
        f"clinical.{label}.modules.E.inputs",
    )
    clinical_a = _tri(clinical_e_inputs.get("A"), f"clinical.{label}.modules.E.inputs.A")
    clinical_f = _dict(clinical_modules.get("F"), f"clinical.{label}.modules.F").get("state")
    if (d_state, corrected_b, c_grade, corrected_a, f_state) != (
        clinical_d,
        clinical_b,
        clinical_c_grade,
        clinical_a,
        clinical_f,
    ):
        raise SourceValidationError(
            f"{label} pruning corrected features differ from locked clinical modules"
        )

    arm_states = _dict(pruning_candidate.get("arm_states"), f"pruning.{label}.arm_states")
    l5 = _dict(arm_states.get(L5_SOURCE_ARM), f"pruning.{label}.arm_states.{L5_SOURCE_ARM}")
    l5_state = l5.get("state")
    allowed_states = E0_STATES if baseline else NON_E0_STATES
    if l5_state not in allowed_states:
        raise SourceValidationError(
            f"pruning.{label} L5 state {l5_state!r} invalid for baseline={baseline}"
        )

    modules = _dict(raw_candidate.get("modules"), f"raw.{label}.modules")
    raw_a = _validate_raw_a(_dict(modules.get("A"), f"raw.{label}.modules.A"), f"raw.{label}.modules.A")
    raw_b = _validate_raw_b(_dict(modules.get("B"), f"raw.{label}.modules.B"), f"raw.{label}.modules.B")
    raw_c = _validate_raw_c(_dict(modules.get("C"), f"raw.{label}.modules.C"), f"raw.{label}.modules.C")
    if corrected_a is not raw_a["positive"]:
        raise SourceValidationError(f"{label} corrected A does not equal frozen raw A")

    organism_name = pruning_candidate.get("organism_name")
    if not isinstance(organism_name, str) or not organism_name.strip():
        organism_name = names[0]
    review_tier = pruning_candidate.get("review_tier")
    if review_tier is not None and not isinstance(review_tier, str):
        raise SourceValidationError(f"pruning.{label}.review_tier must be string or null")

    return {
        "organism_name": organism_name.strip(),
        "canonical_organism_name": names[0],
        "baseline_selected": baseline,
        "review_tier": review_tier,
        "l5_state": l5_state,
        "context": {"D_state": d_state, "F_state": f_state, "clinical_C_grade": c_grade},
        "signals": {
            "CORR": {"A": corrected_a, "B": corrected_b, "C": corrected_c},
            "RAW": {
                "A": raw_a["positive"],
                "B": raw_b["positive"],
                "C": raw_c["positive"],
            },
        },
        "signal_audit": {"RAW_A": raw_a, "RAW_B": raw_b, "RAW_C": raw_c},
    }


def load_and_join_patient(
    pruning_path: str | Path,
    clinical_path: str | Path,
    raw_path: str | Path,
    *,
    anchor: dict[str, Any],
) -> dict[str, Any]:
    pruning_path = Path(pruning_path)
    clinical_path = Path(clinical_path)
    raw_path = Path(raw_path)
    pruning = load_json(pruning_path)
    clinical = load_json(clinical_path)
    raw = load_json(raw_path)
    expected = (
        (pruning, PRUNING_SCHEMA, "pruning"),
        (clinical, CLINICAL_SCHEMA, "clinical"),
        (raw, RAW_SCHEMA, "raw"),
    )
    for artifact, schema, label in expected:
        if artifact.get("schema_version") != schema:
            raise SourceValidationError(
                f"{label} schema must be {schema!r}, got {artifact.get('schema_version')!r}"
            )
        if artifact.get("run_status") != "complete":
            raise SourceValidationError(f"{label} artifact must have run_status=complete")
    ids = tuple(_patient_id(artifact.get("patient_id"), label) for artifact, _, label in expected)
    if not ids[0] == ids[1] == ids[2]:
        raise SourceValidationError(f"Patient ID mismatch: {ids!r}")

    required_anchor_keys = {
        "manifest_path",
        "manifest_sha256",
        "pruning_sha256",
        "pruning_decision_sha256",
        "clinical_sha256",
    }
    if not isinstance(anchor, dict) or set(anchor) != required_anchor_keys:
        raise SourceValidationError("A complete locked pruning-manifest anchor is required")
    manifest_path = Path(anchor["manifest_path"])
    if not manifest_path.is_file() or sha256_file(manifest_path) != anchor["manifest_sha256"]:
        raise SourceValidationError("Pruning manifest anchor path/SHA is stale or invalid")
    manifest = load_json(manifest_path)
    if (
        manifest.get("schema_version") != PRUNING_MANIFEST_SCHEMA
        or manifest.get("run_status") != "complete"
    ):
        raise SourceValidationError("Pruning manifest anchor has invalid schema/status")
    matching_entries = [
        entry
        for entry in manifest.get("artifacts", [])
        if isinstance(entry, dict)
        and _patient_id(entry.get("patient_id"), "manifest artifact") == ids[0]
    ]
    if len(matching_entries) != 1:
        raise SourceValidationError("Pruning manifest must contain exactly one patient anchor")
    manifest_entry = matching_entries[0]
    if (
        manifest_entry.get("sha256") != anchor["pruning_sha256"]
        or manifest_entry.get("decision_sha256") != anchor["pruning_decision_sha256"]
        or manifest_entry.get("source_sha256") != anchor["clinical_sha256"]
    ):
        raise SourceValidationError("Provided anchor values are not recorded in pruning manifest")
    pruning_sha = sha256_file(pruning_path)
    clinical_sha = sha256_file(clinical_path)
    raw_sha = sha256_file(raw_path)
    if pruning_sha != anchor["pruning_sha256"]:
        raise SourceValidationError("Supplied pruning artifact differs from locked manifest")
    if pruning.get("decision_sha256") != anchor["pruning_decision_sha256"]:
        raise SourceValidationError("Pruning decision SHA differs from locked manifest")
    if clinical_sha != anchor["clinical_sha256"]:
        raise SourceValidationError("Supplied clinical artifact differs from locked manifest")
    pruning_source = _dict(pruning.get("source"), "pruning.source")
    if pruning_source.get("sha256") != clinical_sha:
        raise SourceValidationError("Pruning source SHA does not match supplied clinical artifact")
    clinical_provenance = _dict(clinical.get("provenance"), "clinical.provenance")
    if clinical_provenance.get("frozen_rag_re_sha256") != raw_sha:
        raise SourceValidationError("Clinical frozen RAG_re SHA does not match supplied raw artifact")

    candidate_lists = (
        _list(pruning.get("candidates"), "pruning.candidates"),
        _list(clinical.get("candidates"), "clinical.candidates"),
        _list(raw.get("candidates"), "raw.candidates"),
    )
    if len({len(items) for items in candidate_lists}) != 1 or not candidate_lists[0]:
        raise SourceValidationError("Candidate lists must be nonempty and have identical lengths")
    for items, label in zip(candidate_lists, ("pruning.candidates", "clinical.candidates", "raw.candidates")):
        if not all(isinstance(item, dict) for item in items):
            raise SourceValidationError(f"{label} entries must be objects")
        _validate_unique_names(items, label)

    candidates = [
        _project_candidate(pruning_candidate, clinical_candidate, raw_candidate, index)
        for index, (pruning_candidate, clinical_candidate, raw_candidate) in enumerate(
            zip(*candidate_lists)
        )
    ]
    e0_count = sum(candidate["baseline_selected"] for candidate in candidates)
    source_counts = _dict(pruning.get("counts"), "pruning.counts")
    if source_counts.get("candidate_count") != len(candidates) or source_counts.get("e0_count") != e0_count:
        raise SourceValidationError("Pruning root counts do not reconcile with candidates")

    raw_input_meta = _dict(raw.get("input_meta"), "raw.input_meta")
    return {
        "patient_id": ids[0],
        "sources": {
            "pruning": {
                "path": str(pruning_path.resolve()),
                "sha256": pruning_sha,
                "schema_version": PRUNING_SCHEMA,
                "decision_sha256": pruning.get("decision_sha256"),
                "batch_manifest_path": str(manifest_path.resolve()),
                "batch_manifest_sha256": anchor["manifest_sha256"],
            },
            "clinical": {
                "path": str(clinical_path.resolve()),
                "sha256": clinical_sha,
                "schema_version": CLINICAL_SCHEMA,
            },
            "raw": {
                "path": str(raw_path.resolve()),
                "sha256": raw_sha,
                "schema_version": RAW_SCHEMA,
                "pipeline_fingerprint": raw.get("pipeline_fingerprint"),
                "config_hash": raw.get("config_hash"),
                "candidate_pool_sha256": raw_input_meta.get("candidate_pool_sha256"),
            },
        },
        "counts": {
            "candidate_count": len(candidates),
            "e0_count": e0_count,
            "non_e0_count": len(candidates) - e0_count,
        },
        "candidates": candidates,
    }


def revalidate_joined_envelope(joined: dict[str, Any]) -> None:
    """Bind a projected envelope back to its three files and locked manifest."""

    if not isinstance(joined, dict) or not isinstance(joined.get("sources"), dict):
        raise SourceValidationError("Joined envelope has no file-backed sources")
    sources = joined["sources"]
    if set(sources) != {"pruning", "clinical", "raw"}:
        raise SourceValidationError("Joined envelope source layers are invalid")
    pruning = _dict(sources["pruning"], "joined.sources.pruning")
    clinical = _dict(sources["clinical"], "joined.sources.clinical")
    raw = _dict(sources["raw"], "joined.sources.raw")
    anchor = {
        "manifest_path": pruning.get("batch_manifest_path"),
        "manifest_sha256": pruning.get("batch_manifest_sha256"),
        "pruning_sha256": pruning.get("sha256"),
        "pruning_decision_sha256": pruning.get("decision_sha256"),
        "clinical_sha256": clinical.get("sha256"),
    }
    rebuilt = load_and_join_patient(
        pruning.get("path"),
        clinical.get("path"),
        raw.get("path"),
        anchor=anchor,
    )
    if canonical_json_sha256(rebuilt) != canonical_json_sha256(joined):
        raise SourceValidationError(
            "Joined envelope differs from a fresh manifest-anchored source projection"
        )


def discover_patient_artifacts(directory: str | Path, schema: str) -> dict[str, Path]:
    root = Path(directory)
    if not root.is_dir():
        raise SourceValidationError(f"Input directory does not exist: {root}")
    index: dict[str, Path] = {}
    for path in sorted(root.glob("*.json")):
        artifact = load_json(path)
        if artifact.get("schema_version") != schema:
            continue
        pid = _patient_id(artifact.get("patient_id"), str(path))
        if pid in index:
            raise SourceValidationError(f"Duplicate patient {pid} in {root}")
        index[pid] = path
    if not index:
        raise SourceValidationError(f"No {schema} patient artifacts found in {root}")
    return index

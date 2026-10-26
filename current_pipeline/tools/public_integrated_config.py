"""Explicit private inputs for the public integrated CMV replay.

Only routing/provenance configuration lives here. No cohort roster, clinical
source row, reference answer or patient-specific inference rule is distributed.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

SCHEMA = "public_integrated_run.v1"
REQUIRED = {"schema_version", "source_manifests", "taxonomy_completion_manifest",
            "episode_identity_file", "source_clarifications", "reference_root", "development_freeze"}
ALLOWED = REQUIRED | {"partitions"}


def cohort_runs(config: dict) -> dict[str, str]:
    result = {}
    partitions = config.get("partitions", {})
    for run_name in config["source_manifests"]:
        names = list(partitions.get(run_name, {})) or [run_name]
        for name in names:
            if name in result:
                raise ValueError(f"Duplicate output cohort: {name}")
            result[name] = run_name
    return result


def cohort_for(config: dict, run_name: str, patient_id: str) -> str:
    groups = config.get("partitions", {}).get(run_name)
    if not groups:
        return run_name
    matches = [cohort for cohort, members in groups.items() if str(patient_id) in members]
    if len(matches) != 1:
        raise ValueError("A configured patient must resolve exactly one output cohort")
    return matches[0]


def load_config(path: Path) -> dict:
    path = path.resolve()
    config = json.loads(path.read_text(encoding="utf-8-sig"))
    if set(config) - ALLOWED or not REQUIRED <= set(config):
        raise ValueError("Invalid inference configuration keys; answer/label files belong only to evaluation")
    if config["schema_version"] != SCHEMA:
        raise ValueError("Unsupported run configuration schema")
    if not isinstance(config["source_manifests"], dict) or not config["source_manifests"]:
        raise ValueError("At least one explicit source manifest is required")
    hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()}

    def resolve_file(value):
        if not isinstance(value, str) or not value:
            raise ValueError("Input paths must be nonempty strings")
        resolved = (path.parent / value).resolve()
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        hashes[str(resolved)] = hashlib.sha256(resolved.read_bytes()).hexdigest()
        return resolved

    def valid_name(name):
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", name):
            raise ValueError("Cohort names must be safe directory identifiers")

    manifests = {}
    expected_ids = {}
    for name, value in config["source_manifests"].items():
        valid_name(name)
        source = resolve_file(value)
        payload = json.loads(source.read_text(encoding="utf-8-sig"))
        if not payload.get("answer_blind_decision_stages"):
            raise ValueError("Source manifest must declare answer-blind decision stages")
        inputs = payload.get("input", {})
        for role in ("inventory_root", "patient_root", "phenotype_root"):
            location = Path(inputs.get(role, ""))
            if not inputs.get(role) or not location.is_absolute() or not location.is_dir():
                raise ValueError(f"Source manifest needs an existing absolute {role}")
        if not inputs.get("source_summary_suffix"):
            raise ValueError("Source summary suffix must be explicit")
        inventory_files = sorted((Path(inputs["inventory_root"]) / "patients").glob("NGS_patient_*_multi_assay_mngs_evidence.json"))
        ids = {p.name.split("_")[2] for p in inventory_files}
        if not ids:
            raise ValueError("Source inventory has no patients")
        manifests[name] = source
        expected_ids[name] = ids
    config["source_manifests"] = manifests
    for role in ("taxonomy_completion_manifest", "episode_identity_file", "source_clarifications", "development_freeze"):
        config[role] = resolve_file(config[role])
    completion = json.loads(config["taxonomy_completion_manifest"].read_text(encoding="utf-8-sig"))
    overlay = Path(completion.get("staged_overlay", ""))
    if not overlay.is_absolute() or not overlay.is_file():
        raise ValueError("Taxonomy completion requires an existing absolute staged_overlay")
    hashes[str(overlay)] = hashlib.sha256(overlay.read_bytes()).hexdigest()
    if hashes[str(overlay)] != completion.get("staged_overlay_sha256"):
        raise ValueError("Taxonomy completion overlay hash mismatch")
    notes = json.loads(config["source_clarifications"].read_text(encoding="utf-8-sig"))
    if not isinstance(notes.get("clarifications"), list):
        raise ValueError("Source clarifications must contain an explicit list, possibly empty")
    config["reference_root"] = (path.parent / config["reference_root"]).resolve()
    if not config["reference_root"].is_dir():
        raise FileNotFoundError(config["reference_root"])
    partitions = config.setdefault("partitions", {})
    if not isinstance(partitions, dict) or set(partitions) - set(manifests):
        raise ValueError("Partitions must refer only to configured input runs")
    for run_name, groups in partitions.items():
        if not isinstance(groups, dict) or not groups:
            raise ValueError("A declared partition must be nonempty")
        seen = set()
        for name, values in groups.items():
            valid_name(name)
            if not isinstance(values, list) or not values:
                raise ValueError("Partition membership must be an explicit nonempty list")
            values = [str(v) for v in values]
            if len(set(values)) != len(values) or seen.intersection(values):
                raise ValueError("Partition patients must not overlap or repeat")
            seen.update(values)
            groups[name] = values
        if seen != expected_ids[run_name]:
            raise ValueError("Partition must cover the complete inventory exactly")
    cohort_runs(config)
    for run_name, ids in expected_ids.items():
        for patient in ids:
            baseline = config["reference_root"] / run_name / "patient_results" / f"patient_{patient}_integrated_decision.json"
            resolve_file(str(baseline))
    config["_input_hashes"] = hashes
    return config

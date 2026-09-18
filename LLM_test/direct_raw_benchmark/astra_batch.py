"""Offline preparation/collection for the frozen Astra Direct-Raw experiment.

Paid requests are made only by the unchanged direct_raw_runner.py. This helper
never reads Gold, changes the prompt, or sends clinical data over the network.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
import shutil

import direct_raw_runner as r
from input_identity import sha256_file

RUN = r.DEFAULT_OUTPUT_ROOT / "direct_raw_41_astra_idfixed_20260906"
MODEL = "gpt-6-astra"
SOURCES = ("direct_raw_runner.py", "input_identity.py", "prompt_v1.txt",
           "input_mapping_v2.json", "cohort_manifest.json", "evaluate_direct_raw.py")


def context():
    for name in SOURCES:
        if sha256_file(RUN / "source_snapshots" / name) != sha256_file(r.SCRIPT_DIR / name):
            raise ValueError(f"Frozen source changed: {name}")
    _, specs = r.load_manifest(r.DEFAULT_MANIFEST)
    audits = {s.patient_id: r.build_patient_payload(s)[1] for s in specs}
    if len(specs) != 41 or any(a['remaining_exact_date_hits'] or not a['index_date_available'] for a in audits.values()):
        raise ValueError("Frozen cohort/input validation failed")
    prompt_hash = r.sha256_text(r.DEFAULT_PROMPT.read_text(encoding="utf-8-sig").strip())
    return specs, audits, prompt_hash


def validate_output(path, spec, audit, prompt_hash):
    result = r.load_json(path)
    r.validate_cached_output(result, spec, audit, prompt_hash, MODEL, "medium", 4000)
    if result.get("input_identity") != audit["input_identity"]:
        raise ValueError(f"Identity mismatch: {path.name}")
    if result.get("model_returned") != MODEL or not result.get("response_id"):
        raise ValueError(f"Wrong/missing returned model or response ID: {path.name}")
    return result


def prepare():
    specs, audits, prompt_hash = context()
    plan_file = RUN / "execution_plan.json"
    if plan_file.exists():
        raise ValueError("Execution plan already exists; use it to resume, do not overwrite")
    pilot = specs[0]
    validate_output(RUN / f"{pilot.patient_id}_direct_raw.json", pilot, audits[pilot.patient_id], prompt_hash)
    pilot_log = RUN / "execution_logs" / "pilot"
    pilot_log.mkdir(parents=True, exist_ok=False)
    for name in ("run_summary.json", "run_configuration.json", "input_audit_v2.json"):
        shutil.copyfile(RUN / name, pilot_log / name)
    remaining = [s.patient_id for s in specs[1:]]
    groups = {f"shard_{i+1:02d}": remaining[i::3] for i in range(3)}
    plan = {
        "schema_version": "direct_raw_benchmark.astra_plan.v1",
        "model": MODEL, "reasoning_effort": "medium", "max_output_tokens": 4000,
        "pilot": pilot.patient_id, "patient_ids": [s.patient_id for s in specs],
        "shards": groups, "prompt_sha256": prompt_hash,
        "source_sha256": {name: sha256_file(r.SCRIPT_DIR / name) for name in SOURCES},
        "payload_sha256": {pid: a["payload_sha256"] for pid, a in audits.items()},
        "notes": ["Pilot is first manifest patient, selected before looking at predictions.",
                  "Three isolated output shards; each patient gets an independent request.",
                  "All 41 are new Astra responses; no Luna/Terra/Sol predictions are reused."]}
    r.write_json(plan_file, plan)
    print(r.stable_json(plan["shards"]))


def verify_execution(folder, expected_ids):
    summary = r.load_json(folder / "run_summary.json")
    if summary["failures"] or set(summary["success_patient_ids"]) != set(expected_ids):
        raise ValueError(f"Incomplete/failed execution: {folder}")
    config = r.load_json(folder / "run_configuration.json")
    expected = {"model": MODEL, "reasoning_effort": "medium", "max_output_tokens": 4000,
                "store": False, "tools": [],
                "runner_sha256": sha256_file(r.SCRIPT_DIR / "direct_raw_runner.py"),
                "output_schema_sha256": r.sha256_text(r.stable_json(r.OUTPUT_SCHEMA))}
    for field, value in expected.items():
        if config[field] != value:
            raise ValueError(f"Different request configuration: {folder}: {field}")


def collect():
    specs, audits, prompt_hash = context()
    plan = r.load_json(RUN / "execution_plan.json")
    expected_ids = [s.patient_id for s in specs]
    if plan["patient_ids"] != expected_ids or plan["prompt_sha256"] != prompt_hash:
        raise ValueError("Plan differs from current frozen cohort")
    if plan["source_sha256"] != {name: sha256_file(r.SCRIPT_DIR / name) for name in SOURCES}:
        raise ValueError("Source differs from execution plan")
    if plan["payload_sha256"] != {pid: a["payload_sha256"] for pid, a in audits.items()}:
        raise ValueError("Payload differs from execution plan")
    source_dirs = {plan["pilot"]: RUN}
    verify_execution(RUN / "execution_logs" / "pilot", [plan["pilot"]])
    for shard, ids in plan["shards"].items():
        folder = RUN / "shards" / shard
        verify_execution(folder, ids)
        files = {p.stem.removesuffix("_direct_raw") for p in folder.glob("P*_direct_raw.json")}
        if files != set(ids):
            raise ValueError(f"Shard contains missing/extra patients: {shard}")
        for pid in ids:
            if pid in source_dirs:
                raise ValueError(f"Duplicate patient assignment: {pid}")
            source_dirs[pid] = folder
    if set(source_dirs) != set(expected_ids):
        raise ValueError("Shard union differs from frozen cohort")
    existing = {p.stem.removesuffix("_direct_raw") for p in RUN.glob("P*_direct_raw.json")}
    if not existing <= set(expected_ids):
        raise ValueError("Unexpected patient file in final folder")
    records = []
    provenance = []
    for spec in specs:
        source = source_dirs[spec.patient_id] / f"{spec.patient_id}_direct_raw.json"
        record = validate_output(source, spec, audits[spec.patient_id], prompt_hash)
        destination = RUN / source.name
        if destination.exists() and sha256_file(destination) != sha256_file(source):
            raise ValueError(f"Refusing to overwrite a different output: {destination.name}")
        records.append(record)
        provenance.append({"patient_id": spec.patient_id, "source_file": str(source),
                           "output_sha256": sha256_file(source), "response_id": record["response_id"]})
    if len({x["response_id"] for x in records}) != 41:
        raise ValueError("Response IDs are not unique")
    # No copying/summary writes until every shard/output passes validation.
    for item in provenance:
        destination = RUN / f"{item['patient_id']}_direct_raw.json"
        if not destination.exists():
            shutil.copyfile(item["source_file"], destination)
    timestamp = datetime.now(timezone.utc).isoformat()
    r.write_json(RUN / "input_audit_v2.json", {
        "schema_version": "direct_raw_benchmark.run_inputs.v2",
        "manifest_sha256": sha256_file(r.DEFAULT_MANIFEST),
        "mapping_sha256": specs[0].mapping_sha256, "records": list(audits.values())})
    r.write_json(RUN / "run_summary.json", {
        "schema_version": "direct_raw_benchmark.run_summary.v1",
        "created_at_utc": timestamp, "run_name": RUN.name, "model": MODEL,
        "reasoning_effort": "medium", "prompt_sha256": prompt_hash,
        "requested_patient_ids": expected_ids, "success_patient_ids": expected_ids,
        "reused_patient_ids": [], "new_success_patient_ids": expected_ids,
        "input_mapping_sha256": specs[0].mapping_sha256, "failures": [],
        "note": "Consolidated pilot plus three shards. All 41 new Astra responses; raw execution summaries retained."})
    r.write_json(RUN / "batch_completion.json", {
        "schema_version": "direct_raw_benchmark.astra_completion.v1", "status": "complete",
        "verified_at_utc": timestamp, "verified_patient_count": 41,
        "new_response_count": 41, "legacy_reused_count": 0, "provenance": provenance,
        "usage": {k: sum(x["usage"].get(k, 0) for x in records) for k in
                  ("input_tokens", "output_tokens", "total_tokens")},
        "notes": ["Usage covers saved successful responses, not possible failed attempts.",
                  "No Gold read during preparation, requests, or collection."]})
    print(f"Verified and collected {len(records)}/41 Astra responses: {RUN}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "collect"))
    args = parser.parse_args()
    {"prepare": prepare, "collect": collect}[args.command]()

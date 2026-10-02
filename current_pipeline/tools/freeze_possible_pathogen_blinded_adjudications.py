"""Validate and freeze complete blind adjudications before identity unblinding."""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from tools.adjudicate_rag_v2_packets import validate_adjudication
from tools.build_index_event_clinical_timeline_shadow import sha256_file


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip()
    ]


def run(packets_path: Path, adjudication_dir: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    packets = read_jsonl(packets_path)
    packet_by_id = {str(packet.get("case_id") or ""): packet for packet in packets}
    if len(packet_by_id) != len(packets) or "" in packet_by_id:
        raise ValueError("Packets contain missing or duplicate case_id values")

    adjudications_path = adjudication_dir / "adjudications.jsonl"
    failures_path = adjudication_dir / "failures.jsonl"
    run_manifest_path = adjudication_dir / "run_manifest.json"
    if not adjudications_path.is_file() or not run_manifest_path.is_file():
        raise FileNotFoundError("Adjudication output is incomplete: missing adjudications or run manifest")
    failures = read_jsonl(failures_path) if failures_path.is_file() else []
    if failures:
        raise ValueError(f"Cannot freeze a run with {len(failures)} failures")

    records = read_jsonl(adjudications_path)
    record_by_id = {str(record.get("case_id") or ""): record for record in records}
    if len(record_by_id) != len(records) or "" in record_by_id:
        raise ValueError("Adjudications contain missing or duplicate case_id values")
    missing = sorted(set(packet_by_id) - set(record_by_id))
    extra = sorted(set(record_by_id) - set(packet_by_id))
    if missing or extra:
        raise ValueError(f"Adjudication coverage mismatch: missing={missing}, extra={extra}")

    for case_id, record in record_by_id.items():
        if record.get("answer_source_read") is not False:
            raise ValueError(f"{case_id} does not explicitly record answer_source_read=false")
        validate_adjudication(record.get("adjudication") or {}, packet_by_id[case_id])

    run_manifest = read_json(run_manifest_path)
    if int(run_manifest.get("completed_count") or 0) != len(packets):
        raise ValueError("Run manifest completed_count does not match packet count")
    if int(run_manifest.get("failure_count") or 0) != 0:
        raise ValueError("Run manifest reports failures")
    if run_manifest.get("answer_source_read") is not False:
        raise ValueError("Run manifest does not explicitly record answer_source_read=false")

    output_dir.mkdir(parents=True, exist_ok=True)
    frozen_adjudications = output_dir / "adjudications.frozen.jsonl"
    frozen_run_manifest = output_dir / "run_manifest.frozen.json"
    shutil.copy2(adjudications_path, frozen_adjudications)
    shutil.copy2(run_manifest_path, frozen_run_manifest)
    manifest = {
        "schema_version": "possible_pathogen_blinded_adjudication_freeze.v1",
        "frozen_at": datetime.now(timezone.utc).isoformat(),
        "answer_source_read": False,
        "private_case_map_read": False,
        "packet_count": len(packets),
        "adjudication_count": len(records),
        "packet_file": str(packets_path.resolve()),
        "packet_sha256": sha256_file(packets_path),
        "frozen_adjudications_file": str(frozen_adjudications.resolve()),
        "frozen_adjudications_sha256": sha256_file(frozen_adjudications),
        "frozen_run_manifest_file": str(frozen_run_manifest.resolve()),
        "frozen_run_manifest_sha256": sha256_file(frozen_run_manifest),
        "backend": run_manifest.get("backend"),
        "model": run_manifest.get("model"),
        "reasoning_effort": run_manifest.get("reasoning_effort"),
        "ready_for_posthoc_unblinding": True,
    }
    (output_dir / "freeze_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packets", type=Path)
    parser.add_argument("adjudication_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.packets, args.adjudication_dir, args.output_dir),
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

"""Merge Luna OBER output with selected Sol re-reviews and reapply guardrails."""

from __future__ import annotations

import argparse
import copy
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from tools.adjudicate_rag_v2_packets import (
    apply_patient_fit_visibility_boundary,
    validate_adjudication,
)


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


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows),
        encoding="utf-8",
    )


def case_records(root: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for path in sorted((root / "cases").glob("*.json")):
        record = read_json(path)
        case_id = str(record.get("case_id") or "")
        if not case_id or case_id in records:
            raise ValueError(f"Missing or duplicate case_id in {path}")
        records[case_id] = record
    return records


def run(
    packets_path: Path,
    luna_dir: Path,
    sol_dir: Path,
    escalation_path: Path,
    output_dir: Path,
    boundary_profile: str,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_dir}")
    packets = read_jsonl(packets_path)
    packet_by_id = {str(packet.get("case_id") or ""): packet for packet in packets}
    luna = case_records(luna_dir)
    sol = case_records(sol_dir)
    escalations = read_jsonl(escalation_path)
    selected = {str(row.get("case_id") or "") for row in escalations}
    if set(luna) != set(packet_by_id):
        raise ValueError("Luna arm does not cover every packet")
    if set(sol) != selected:
        raise ValueError(
            f"Sol arm coverage mismatch: missing={sorted(selected - set(sol))}, "
            f"extra={sorted(set(sol) - selected)}"
        )

    final_records: list[dict[str, Any]] = []
    override_count = 0
    for case_id in sorted(packet_by_id):
        packet = packet_by_id[case_id]
        source = sol[case_id] if case_id in selected else luna[case_id]
        if source.get("answer_source_read") is not False:
            raise ValueError(f"{case_id} does not record answer_source_read=false")
        adjudication = copy.deepcopy(source.get("adjudication") or {})
        validate_adjudication(adjudication, packet)
        original_decision = copy.deepcopy(adjudication.get("decision") or {})
        adjudication, boundary = apply_patient_fit_visibility_boundary(
            adjudication, packet, boundary_profile
        )
        validate_adjudication(adjudication, packet)
        if adjudication.get("decision") != original_decision:
            override_count += 1
        record = {
            **source,
            "adjudication": adjudication,
            "boundary_guardrail": {"profile": boundary_profile, **boundary},
            "answer_source_read": False,
            "two_stage_provenance": {
                "selected_arm": "sol_escalation" if case_id in selected else "luna_primary",
                "luna_model": read_json(luna_dir / "run_manifest.json").get("model"),
                "sol_model": read_json(sol_dir / "run_manifest.json").get("model") if selected else None,
                "deterministic_boundary_reapplied": True,
            },
        }
        output_path = output_dir / "cases" / f"{re.sub(r'[^A-Za-z0-9_.-]+', '_', case_id)}.json"
        write_json(output_path, record)
        final_records.append(record)

    write_jsonl(output_dir / "adjudications.jsonl", final_records)
    write_jsonl(output_dir / "failures.jsonl", [])
    manifest = {
        "schema_version": "ober_two_stage_adjudication_manifest_v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "backend": "codex-cli",
        "model": "luna_primary_with_sol_escalation",
        "reasoning_effort": "luna=medium;sol=high",
        "boundary_profile": boundary_profile,
        "packet_count": len(packets),
        "completed_count": len(final_records),
        "failure_count": 0,
        "luna_primary_count": len(packets) - len(selected),
        "sol_escalation_count": len(selected),
        "deterministic_boundary_override_count": override_count,
        "answer_source_read": False,
        "fixed_acceptance_actions": ["primary_pathogen", "secondary_pathogen"],
    }
    write_json(output_dir / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--packets", type=Path, required=True)
    parser.add_argument("--luna-dir", type=Path, required=True)
    parser.add_argument("--sol-dir", type=Path, required=True)
    parser.add_argument("--escalations", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--boundary-profile", default="patient-fit-v1")
    args = parser.parse_args()
    print(
        json.dumps(
            run(
                args.packets,
                args.luna_dir,
                args.sol_dir,
                args.escalations,
                args.output_dir,
                args.boundary_profile,
            ),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

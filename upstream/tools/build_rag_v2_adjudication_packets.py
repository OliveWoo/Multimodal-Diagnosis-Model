"""Join answer-blind RAG cases with PubMed retrieval results for adjudication."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

from tools.build_rag_v2_inputs import FORBIDDEN_PRODUCTION_KEYS


SCHEMA_VERSION = "rag_adjudication_packet_v2.0"


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip()
    ]


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def forbidden_paths(value: Any, prefix: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            path = f"{prefix}.{key}"
            if str(key).lower() in FORBIDDEN_PRODUCTION_KEYS:
                found.append(path)
            found.extend(forbidden_paths(child, path))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.extend(forbidden_paths(child, f"{prefix}[{index}]"))
    return found


def make_packet(case: dict[str, Any], retrieval: dict[str, Any]) -> dict[str, Any]:
    rag_task = case.get("rag_task") if isinstance(case.get("rag_task"), dict) else {}
    packet = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case.get("case_id"),
        "dataset_label": case.get("dataset_label"),
        "patient_id": case.get("patient_id"),
        "organism": case.get("organism"),
        "current_model_state": case.get("current_model_state"),
        "local_evidence": case.get("local_evidence"),
        "literature_evidence": {
            "knowledge_card_key": retrieval.get("knowledge_card_key"),
            "retrieved_at": retrieval.get("retrieved_at"),
            "retrieval_evidence_volume": retrieval.get("retrieval_evidence_volume"),
            "selected_article_count": retrieval.get("selected_article_count"),
            "selected_track_counts": retrieval.get("selected_track_counts"),
            "missing_evidence_tracks": retrieval.get("missing_evidence_tracks"),
            "retrieval_errors": retrieval.get("retrieval_errors"),
            "articles": retrieval.get("articles"),
            "adjudication_note": retrieval.get("adjudication_note"),
        },
        "adjudication_task": {
            "fixed_questions": rag_task.get("fixed_questions", []),
            "required_output_schema": rag_task.get("required_output_schema", "rag_adjudication_v2.0"),
            "constraints": [
                "Use the literature and patient evidence separately before assessing patient fit.",
                "Do not treat retrieval relevance score as evidence direction or medical quality.",
                "Do not infer support from article count alone; assess study quality and case similarity.",
                "If evidence is sparse or conflicting, preserve uncertainty instead of upgrading.",
                "Return one adjudication object for this case_id and follow the required output schema.",
            ],
        },
        "provenance": {
            **(case.get("provenance") if isinstance(case.get("provenance"), dict) else {}),
            "production_case_schema": case.get("schema_version"),
            "retrieval_schema": retrieval.get("schema_version"),
            "answer_source_read": False,
        },
    }
    leaked = forbidden_paths(packet)
    if leaked:
        raise ValueError(f"Forbidden benchmark fields in {case.get('case_id')}: {', '.join(leaked)}")
    return packet


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("production_queue", type=Path)
    parser.add_argument("retrieval_results", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--case-ids", nargs="*", help="Optional case_id allowlist.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    cases = read_jsonl(args.production_queue)
    retrieval_rows = read_jsonl(args.retrieval_results)
    retrieval_by_key = {
        str(row.get("knowledge_card_key")): row
        for row in retrieval_rows
        if row.get("knowledge_card_key")
    }
    allowlist = set(args.case_ids or [])
    if allowlist:
        cases = [case for case in cases if str(case.get("case_id")) in allowlist]

    packets: list[dict[str, Any]] = []
    missing: list[dict[str, str]] = []
    for case in cases:
        rag_task = case.get("rag_task") if isinstance(case.get("rag_task"), dict) else {}
        key = str(rag_task.get("knowledge_card_key") or "")
        retrieval = retrieval_by_key.get(key)
        if retrieval is None:
            missing.append({"case_id": str(case.get("case_id") or ""), "knowledge_card_key": key})
            continue
        packets.append(make_packet(case, retrieval))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_jsonl(args.output_dir / "adjudication_packets.jsonl", packets)
    with (args.output_dir / "packet_summary.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        fields = [
            "case_id",
            "patient_id",
            "organism_name",
            "source_tier",
            "knowledge_card_key",
            "retrieval_evidence_volume",
            "selected_article_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for packet in packets:
            writer.writerow(
                {
                    "case_id": packet.get("case_id"),
                    "patient_id": packet.get("patient_id"),
                    "organism_name": (packet.get("organism") or {}).get("display_name"),
                    "source_tier": (packet.get("current_model_state") or {}).get("source_tier"),
                    "knowledge_card_key": (packet.get("literature_evidence") or {}).get("knowledge_card_key"),
                    "retrieval_evidence_volume": (packet.get("literature_evidence") or {}).get("retrieval_evidence_volume"),
                    "selected_article_count": (packet.get("literature_evidence") or {}).get("selected_article_count"),
                }
            )
    manifest = {
        "schema_version": "rag_adjudication_packet_manifest_v2.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "production_queue": str(args.production_queue.resolve()),
        "retrieval_results": str(args.retrieval_results.resolve()),
        "input_case_count": len(cases),
        "packet_count": len(packets),
        "missing_retrieval_count": len(missing),
        "missing_retrieval": missing,
        "forbidden_benchmark_field_count": sum(len(forbidden_paths(packet)) for packet in packets),
        "answer_source_read": False,
    }
    (args.output_dir / "packet_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"packets={len(packets)}")
    print(f"missing_retrieval={len(missing)}")
    print("forbidden_benchmark_fields=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

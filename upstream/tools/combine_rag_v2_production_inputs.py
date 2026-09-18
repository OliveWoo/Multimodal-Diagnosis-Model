"""Combine multiple RAG v2 production input folders."""

from __future__ import annotations

import argparse
import json
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Sequence


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return rows


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def knowledge_key(row: dict[str, Any]) -> str:
    value = row.get("knowledge_card_key")
    if value:
        return str(value)
    organism = row.get("organism")
    if isinstance(organism, dict) and organism.get("canonical_key"):
        return str(organism["canonical_key"])
    return json.dumps(row, ensure_ascii=False, sort_keys=True)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Combine RAG v2 production input folders.")
    parser.add_argument("--input-dir", type=Path, action="append", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    cases: list[dict[str, Any]] = []
    knowledge_requests: dict[str, dict[str, Any]] = {}
    manifests: list[dict[str, Any]] = []
    source_manifest_paths: list[str] = []

    for input_dir in args.input_dir:
        manifest_path = input_dir / "production_manifest.json"
        if manifest_path.exists():
            manifests.append(read_json(manifest_path))
            source_manifest_paths.append(str(manifest_path))
        cases.extend(read_jsonl(input_dir / "production_queue.jsonl"))
        for row in read_jsonl(input_dir / "knowledge_card_requests.jsonl"):
            knowledge_requests.setdefault(knowledge_key(row), row)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_jsonl(args.output_dir / "production_queue.jsonl", cases)
    write_jsonl(
        args.output_dir / "knowledge_card_requests.jsonl",
        [knowledge_requests[key] for key in sorted(knowledge_requests)],
    )

    schema_source = next((path / "schemas" for path in args.input_dir if (path / "schemas").exists()), None)
    if schema_source is not None:
        shutil.copytree(schema_source, args.output_dir / "schemas", dirs_exist_ok=True)

    tier_counts = Counter(
        str(((case.get("current_model_state") or {}).get("source_tier") or ""))
        for case in cases
    )
    category_counts = Counter(
        str(((case.get("organism") or {}).get("pathogen_category") or ""))
        for case in cases
    )
    manifest = {
        "schema_version": "rag_production_manifest_v2.0_combined",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "dataset_label": args.label,
        "source_manifests": source_manifest_paths,
        "production_case_count": len(cases),
        "knowledge_request_count": len(knowledge_requests),
        "source_case_counts": [manifest.get("production_case_count") for manifest in manifests],
        "source_knowledge_request_counts": [
            manifest.get("knowledge_request_count") for manifest in manifests
        ],
        "tier_counts": dict(tier_counts),
        "category_counts": dict(category_counts),
        "production_validation": {
            "forbidden_benchmark_key_count": sum(
                int((manifest.get("production_validation") or {}).get("forbidden_benchmark_key_count") or 0)
                for manifest in manifests
            ),
            "answer_source_read": False,
        },
    }
    write_json(args.output_dir / "production_manifest.json", manifest)
    (args.output_dir / "README_繁中.md").write_text(
        (
            "# Combined RAG v2 production input\n\n"
            "這個資料夾合併多個 cohort 的 RAG production input。\n\n"
            "- `production_queue.jsonl`：RAG 要逐案判斷的 high/context review candidates。\n"
            "- `knowledge_card_requests.jsonl`：依菌種去重後的文獻知識卡需求。\n"
            "- `schemas/`：RAG 輸入與輸出 schema。\n\n"
            "不含 benchmark answer 或 answer-hit 標記。RAG 應判斷每個 candidate 是 "
            "`primary_pathogen`、`secondary_pathogen`、`context_only`、`reject` 或 "
            "`insufficient_evidence`。\n"
        ),
        encoding="utf-8",
    )
    print(f"combined_cases={len(cases)}")
    print(f"combined_knowledge_requests={len(knowledge_requests)}")
    print(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

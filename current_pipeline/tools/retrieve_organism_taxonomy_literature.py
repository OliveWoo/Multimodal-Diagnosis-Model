"""Retrieve answer-blind PubMed abstracts for organism-level taxonomy review."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

from tools.draft_organism_taxonomy_reviews import load_queue
from tools.retrieve_rag_v2_pubmed import PubMedClient


def retrieve(queue_path: Path, output_path: Path, *, cache_dir: Path, limit: int | None = None) -> dict[str, Any]:
    queue = load_queue(queue_path)
    client = PubMedClient(
        cache_dir=cache_dir,
        email=os.environ.get("NCBI_EMAIL"),
        api_key=os.environ.get("NCBI_API_KEY"),
        timeout=30,
        retries=3,
        refresh_cache=False,
    )
    existing: dict[str, Any] = {}
    if output_path.exists():
        existing = json.loads(output_path.read_text(encoding="utf-8"))
        if not isinstance(existing, dict):
            raise ValueError("Existing literature pack must be a JSON object")
    new_count = 0
    for item in queue[:limit]:
        key = item["canonical_key"]
        if key in existing:
            continue
        names = [item["organism_name"]]
        scientific = str(item.get("ncbi", {}).get("scientific_name") or "")
        if scientific and scientific.casefold() != names[0].casefold():
            names.append(scientific)
        articles: dict[str, dict[str, Any]] = {}
        search_counts = {}
        for name in names:
            query = f'"{name}"[Title/Abstract]'
            total, pmids = client.search(query, 10)
            search_counts[name] = total
            for article in client.fetch(pmids):
                title_and_abstract = f"{article.get('title', '')} {article.get('abstract', '')}".casefold()
                if not any(candidate.casefold() in title_and_abstract for candidate in names):
                    continue
                pmid = str(article.get("pmid") or "")
                if not pmid:
                    continue
                articles[pmid] = {
                    "pmid": pmid,
                    "title": article.get("title", ""),
                    "abstract": article.get("abstract", ""),
                    "year": article.get("year"),
                    "publication_types": article.get("publication_types", []),
                    "url": article.get("url", ""),
                }
        existing[key] = {
            "organism_name": item["organism_name"],
            "canonical_key": key,
            "search_counts": search_counts,
            "articles": list(articles.values()),
            "source": "PubMed E-utilities; abstracts are search hits, not validated clinical evidence",
        }
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(existing, ensure_ascii=False, indent=2), encoding="utf-8")
        new_count += 1
    return {
        "queued": len(queue),
        "retrieved_this_run": new_count,
        "pack_entries": len(existing),
        "entries_with_articles": sum(bool(item.get("articles")) for item in existing.values()),
        "output": str(output_path.resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch PubMed title/abstract evidence for organism-only review.")
    parser.add_argument("queue", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    print(json.dumps(retrieve(args.queue, args.output, cache_dir=args.cache_dir, limit=args.limit), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
import csv
import json
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _request_json(endpoint: str, params: dict[str, Any]) -> dict[str, Any]:
    query = dict(params)
    query["tool"] = "csie_project_taxonomy_audit"
    if os.environ.get("NCBI_EMAIL"):
        query["email"] = os.environ["NCBI_EMAIL"]
    if os.environ.get("NCBI_API_KEY"):
        query["api_key"] = os.environ["NCBI_API_KEY"]
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(query)}"
    request = urllib.request.Request(url, headers={"User-Agent": "CSIE-project-taxonomy-audit/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.loads(response.read().decode("utf-8"))


def _request_xml(endpoint: str, params: dict[str, Any]) -> ET.Element:
    query = dict(params)
    query["tool"] = "csie_project_taxonomy_audit"
    if os.environ.get("NCBI_EMAIL"):
        query["email"] = os.environ["NCBI_EMAIL"]
    if os.environ.get("NCBI_API_KEY"):
        query["api_key"] = os.environ["NCBI_API_KEY"]
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(query)}"
    request = urllib.request.Request(url, headers={"User-Agent": "CSIE-project-taxonomy-audit/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return ET.fromstring(response.read())


def _direct_text(node: ET.Element, tag: str) -> str:
    child = node.find(tag)
    return (child.text or "").strip() if child is not None else ""


def parse_taxa(root: ET.Element) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for taxon in root.findall("./Taxon"):
        taxid = _direct_text(taxon, "TaxId")
        if not taxid:
            continue
        lineage_ex = []
        for item in taxon.findall("./LineageEx/Taxon"):
            lineage_ex.append(
                {
                    "taxid": _direct_text(item, "TaxId"),
                    "scientific_name": _direct_text(item, "ScientificName"),
                    "rank": _direct_text(item, "Rank"),
                }
            )
        records[taxid] = {
            "taxid": taxid,
            "scientific_name": _direct_text(taxon, "ScientificName"),
            "rank": _direct_text(taxon, "Rank"),
            "division": _direct_text(taxon, "Division"),
            "lineage": _direct_text(taxon, "Lineage"),
            "lineage_ex": lineage_ex,
        }
    return records


def _read_gap_names(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    unique: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row["organism_name"].strip().casefold()
        unique.setdefault(
            key,
            {
                "input_name": row["organism_name"].strip(),
                "category": row.get("category", "").strip(),
            },
        )
    return sorted(unique.values(), key=lambda row: row["input_name"].casefold())


def resolve(gap_csv: Path, output_dir: Path, delay_seconds: float = 0.36) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_path = output_dir / "ncbi_taxonomy_cache.json"
    cache: dict[str, Any] = (
        json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {"searches": {}, "taxa": {}}
    )
    cache.setdefault("searches", {})
    cache.setdefault("taxa", {})

    names = _read_gap_names(gap_csv)
    for index, item in enumerate(names, start=1):
        name = item["input_name"]
        if name in cache["searches"]:
            continue
        response = _request_json(
            "esearch.fcgi",
            {
                "db": "taxonomy",
                "term": f'"{name}"[All Names]',
                "retmode": "json",
                "retmax": 20,
            },
        )
        result = response.get("esearchresult", {})
        cache["searches"][name] = {
            "count": int(result.get("count", 0)),
            "idlist": list(result.get("idlist", [])),
            "querytranslation": result.get("querytranslation", ""),
        }
        cache_path.write_text(json.dumps(cache, ensure_ascii=False, indent=2), encoding="utf-8")
        if index < len(names):
            time.sleep(delay_seconds)

    ids = sorted(
        {
            taxid
            for search in cache["searches"].values()
            for taxid in search.get("idlist", [])
            if taxid and taxid not in cache["taxa"]
        },
        key=lambda value: int(value),
    )
    for start in range(0, len(ids), 100):
        batch = ids[start : start + 100]
        root = _request_xml(
            "efetch.fcgi",
            {"db": "taxonomy", "id": ",".join(batch), "retmode": "xml"},
        )
        cache["taxa"].update(parse_taxa(root))
        cache_path.write_text(json.dumps(cache, ensure_ascii=False, indent=2), encoding="utf-8")
        if start + 100 < len(ids):
            time.sleep(delay_seconds)

    rows: list[dict[str, Any]] = []
    for item in names:
        name = item["input_name"]
        search = cache["searches"].get(name, {})
        candidates = [cache["taxa"].get(taxid, {}) for taxid in search.get("idlist", [])]
        exact = [
            record
            for record in candidates
            if str(record.get("scientific_name", "")).casefold() == name.casefold()
        ]
        if len(exact) == 1:
            chosen = exact[0]
            status = "exact_scientific_name"
        elif len(candidates) == 1:
            chosen = candidates[0]
            status = "single_ncbi_match"
        elif candidates:
            chosen = {}
            status = "ambiguous_multiple_matches"
        else:
            chosen = {}
            status = "not_found"
        rows.append(
            {
                **item,
                "resolution_status": status,
                "search_count": search.get("count", 0),
                "ncbi_taxid": chosen.get("taxid", ""),
                "ncbi_scientific_name": chosen.get("scientific_name", ""),
                "ncbi_rank": chosen.get("rank", ""),
                "ncbi_division": chosen.get("division", ""),
                "ncbi_lineage": chosen.get("lineage", ""),
                "source": "NCBI Taxonomy via E-utilities",
            }
        )

    csv_path = output_dir / "ncbi_taxonomy_resolution.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "input_unique_names": len(rows),
        "resolution_status": {},
        "output_csv": str(csv_path.resolve()),
        "cache": str(cache_path.resolve()),
    }
    for row in rows:
        status = row["resolution_status"]
        summary["resolution_status"][status] = summary["resolution_status"].get(status, 0) + 1
    (output_dir / "ncbi_taxonomy_resolution_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Resolve taxonomy-gap names against NCBI Taxonomy.")
    parser.add_argument("gap_csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--delay-seconds", type=float, default=0.36)
    args = parser.parse_args()
    print(json.dumps(resolve(args.gap_csv, args.output_dir, args.delay_seconds), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

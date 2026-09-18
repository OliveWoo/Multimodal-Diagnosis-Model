from __future__ import annotations

import json
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List
from xml.etree import ElementTree as ET

import requests

from .io_utils import stable_hash, write_json_atomic


NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class PubMedError(RuntimeError):
    pass


def _element_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return " ".join("".join(element.itertext()).split())


def _query_phrase(value: str) -> str:
    return " ".join(value.replace('"', " ").split())


def _site_query(query_terms: str) -> str:
    terms = [" ".join(term.strip().split()) for term in query_terms.split(" OR ")]
    quoted = [f'"{term}"[Title/Abstract]' for term in terms if term]
    return " OR ".join(quoted) or '"human infection"[Title/Abstract]'


def build_query(organism: str, context: Dict[str, str], cutoff: str | None) -> str:
    organism = _query_phrase(organism)
    query = f'"{organism}"[Title/Abstract] AND ({_site_query(context["query_terms"])})'
    if cutoff:
        normalized = cutoff.replace("-", "/")
        query += f" AND (1900/01/01:{normalized}[dp])"
    return query


class PubMedClient:
    def __init__(self, config: Dict[str, Any], cache_dir: Path | None = None):
        self.config = config
        self.cache_dir = cache_dir
        self.session = requests.Session()
        self.api_key = os.getenv("PUBMED_API_KEY") or ""

    def _common_params(self) -> Dict[str, Any]:
        params: Dict[str, Any] = {
            "db": "pubmed",
            "tool": self.config.get("tool") or "csie_rag_re",
            "email": self.config.get("email") or "",
        }
        if self.api_key:
            params["api_key"] = self.api_key
        return params

    def _get(self, endpoint: str, params: Dict[str, Any]) -> requests.Response:
        retries = max(1, int(self.config.get("retries", 3)))
        timeout = float(self.config.get("timeout_seconds", 30))
        interval = float(self.config.get("request_interval_seconds", 0.34))
        last_error: Exception | None = None
        for attempt in range(1, retries + 1):
            try:
                response = self.session.get(
                    f"{NCBI_BASE}/{endpoint}", params=params, timeout=timeout
                )
                if response.status_code == 429:
                    retry_after = float(response.headers.get("Retry-After", attempt))
                    time.sleep(min(max(retry_after, interval), 10.0))
                    continue
                response.raise_for_status()
                if interval > 0:
                    time.sleep(interval)
                return response
            except (requests.RequestException, ValueError) as exc:
                last_error = exc
                if attempt < retries:
                    time.sleep(min(0.8 * attempt, 5.0))
        raise PubMedError(f"PubMed request failed after {retries} attempts: {last_error}")

    def _cache_path(self, query: str, max_articles: int) -> Path | None:
        if self.cache_dir is None:
            return None
        key = stable_hash({"query": query, "max_articles": max_articles})
        return self.cache_dir / "pubmed" / f"{key}.json"

    def search(
        self, organism: str, context: Dict[str, str], max_articles: int | None = None
    ) -> Dict[str, Any]:
        max_articles = int(max_articles or self.config.get("max_articles", 10))
        query = build_query(organism, context, self.config.get("publication_cutoff"))
        cache_path = self._cache_path(query, max_articles)
        if cache_path is not None and cache_path.exists():
            cached = json.loads(cache_path.read_text(encoding="utf-8-sig"))
            cached["cache_hit"] = True
            return cached

        esearch = self._common_params()
        esearch.update(
            {
                "term": query,
                "retmode": "json",
                "retmax": max_articles,
                "sort": "relevance",
            }
        )
        response = self._get("esearch.fcgi", esearch)
        try:
            result = response.json().get("esearchresult", {})
        except (ValueError, AttributeError) as exc:
            raise PubMedError("PubMed ESearch returned invalid JSON") from exc
        ids = [str(value) for value in result.get("idlist", []) if str(value)]
        try:
            total_hits = int(result.get("count", 0) or 0)
        except (TypeError, ValueError):
            total_hits = 0

        articles: List[Dict[str, Any]] = []
        if ids:
            efetch = self._common_params()
            efetch.update({"id": ",".join(ids), "retmode": "xml"})
            fetched = self._get("efetch.fcgi", efetch)
            try:
                root = ET.fromstring(fetched.content)
            except ET.ParseError as exc:
                raise PubMedError("PubMed EFetch returned invalid XML") from exc
            articles = self._parse_articles(root)

        payload = {
            "status": "ok",
            "query": query,
            "organism": organism,
            "target_site": context.get("target_site", ""),
            "publication_cutoff": self.config.get("publication_cutoff"),
            "total_hits": total_hits,
            "requested_articles": max_articles,
            "returned_articles": len(articles),
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "cache_hit": False,
            "articles": articles,
        }
        if cache_path is not None:
            write_json_atomic(cache_path, payload)
        return payload

    @staticmethod
    def _parse_articles(root: ET.Element) -> List[Dict[str, Any]]:
        articles: List[Dict[str, Any]] = []
        for node in root.findall(".//PubmedArticle"):
            pmid = _element_text(node.find(".//MedlineCitation/PMID"))
            title = _element_text(node.find(".//Article/ArticleTitle"))
            abstract_parts = []
            for abstract in node.findall(".//Article/Abstract/AbstractText"):
                text = _element_text(abstract)
                label = (abstract.attrib.get("Label") or "").strip()
                if text:
                    abstract_parts.append(f"{label}: {text}" if label else text)
            abstract = " ".join(abstract_parts)
            journal = _element_text(node.find(".//Article/Journal/Title"))

            year = _element_text(node.find(".//Article/Journal/JournalIssue/PubDate/Year"))
            if not year:
                medline_date = _element_text(
                    node.find(".//Article/Journal/JournalIssue/PubDate/MedlineDate")
                )
                match = re.search(r"\b(19|20)\d{2}\b", medline_date)
                year = match.group(0) if match else ""

            doi = ""
            for article_id in node.findall(".//PubmedData/ArticleIdList/ArticleId"):
                if article_id.attrib.get("IdType") == "doi":
                    doi = _element_text(article_id)
                    break
            publication_types = [
                _element_text(item)
                for item in node.findall(".//Article/PublicationTypeList/PublicationType")
                if _element_text(item)
            ]
            articles.append(
                {
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "journal": journal or None,
                    "year": int(year) if year.isdigit() else None,
                    "doi": doi or None,
                    "publication_types": publication_types,
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None,
                }
            )
        return articles

"""Retrieve and cache PubMed evidence for answer-blind RAG v2 knowledge requests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_NAME = "csie_mngs_rag_v2"
PULMONARY_TERMS = (
    "pneumonia",
    "pulmonary",
    "lower respiratory",
    "bronchoalveolar",
    "balf",
    "ventilator-associated",
    "ventilator associated",
    "lung abscess",
    "empyema",
    "pneumonitis",
)
BACKGROUND_TERMS = (
    "colonization",
    "colonisation",
    "contaminant",
    "contamination",
    "commensal",
    "microbiome",
    "bystander",
    "reactivation",
    "shedding",
)
HIGHER_EVIDENCE_PUBLICATION_TYPES = {
    "Practice Guideline": 7.0,
    "Guideline": 7.0,
    "Systematic Review": 6.0,
    "Meta-Analysis": 6.0,
    "Review": 3.0,
    "Multicenter Study": 3.0,
    "Observational Study": 2.0,
    "Comparative Study": 2.0,
    "Clinical Trial": 3.0,
    "Case Reports": -2.0,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_jsonl(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_text(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def xml_text(element: ET.Element | None) -> str:
    return clean_text("".join(element.itertext())) if element is not None else ""


def dedupe(values: Sequence[str]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        cleaned = clean_text(value)
        key = cleaned.casefold()
        if cleaned and key not in seen:
            result.append(cleaned)
            seen.add(key)
    return result


def quoted_tiab(term: str) -> str:
    escaped = clean_text(term).replace('"', "")
    return f'"{escaped}"[Title/Abstract]'


def quoted_title(term: str) -> str:
    escaped = clean_text(term).replace('"', "")
    return f'"{escaped}"[Title]'


def query_tracks(request: dict[str, Any]) -> list[dict[str, str]]:
    organism = request.get("organism") if isinstance(request.get("organism"), dict) else {}
    name = clean_text(organism.get("display_name"))
    category = clean_text(organism.get("pathogen_category"))
    species = quoted_tiab(name)
    species_title = quoted_title(name)
    pulmonary = (
        '(pneumonia[Title/Abstract] OR pulmonary[Title/Abstract] OR '
        '"lower respiratory"[Title/Abstract] OR bronchoalveolar[Title/Abstract] OR '
        '"ventilator associated"[Title/Abstract] OR empyema[Title/Abstract] OR '
        '"lung abscess"[Title/Abstract])'
    )
    tracks = [
        {
            "track": "pulmonary_pathogenicity",
            "scope": "species_specific",
            "query": (
                f"(({species_title} AND {pulmonary}) OR "
                f"({species} AND (pneumonia[Title] OR pulmonary[Title] OR "
                '"lower respiratory"[Title] OR bronchoalveolar[Title] OR '
                '"ventilator associated"[Title] OR empyema[Title] OR "lung abscess"[Title])))'
            ),
        },
        {
            "track": "colonization_or_background",
            "scope": "species_specific",
            "query": (
                f"{species} AND (respiratory[Title/Abstract] OR bronchoalveolar[Title/Abstract] OR "
                "pneumonia[Title/Abstract] OR airway[Title/Abstract]) AND "
                "(colonization[Title/Abstract] OR colonisation[Title/Abstract] OR "
                "contaminant[Title/Abstract] OR contamination[Title/Abstract] OR "
                "commensal[Title/Abstract] OR microbiome[Title/Abstract] OR "
                "reactivation[Title/Abstract] OR shedding[Title/Abstract])"
            ),
        },
        {
            "track": "guideline_or_review",
            "scope": "species_specific",
            "query": f"{species} AND {pulmonary} AND (systematic review[Publication Type] OR review[Publication Type])",
        },
    ]
    if category == "Candida/yeast":
        tracks[2] = {
            "track": "guideline_or_review",
            "scope": "category_context",
            "query": (
                "Candida[Title/Abstract] AND (\"respiratory tract\"[Title/Abstract] OR "
                "\"respiratory secretions\"[Title/Abstract] OR bronchoalveolar[Title/Abstract] OR "
                "airway[Title/Abstract]) AND (colonization[Title/Abstract] OR "
                "colonisation[Title/Abstract] OR pneumonia[Title/Abstract]) AND "
                "(guideline[Publication Type] OR systematic review[Publication Type] OR review[Publication Type])"
            ),
        }
    elif category == "Herpesvirus/reactivation":
        tracks[2] = {
            "track": "guideline_or_review",
            "scope": "category_context",
            "query": (
                '("herpes simplex virus"[Title/Abstract] OR HSV[Title/Abstract]) AND '
                '(pneumonia[Title/Abstract] OR pneumonitis[Title/Abstract] OR '
                '"lower respiratory"[Title/Abstract] OR bronchoalveolar[Title/Abstract]) AND '
                '(reactivation[Title/Abstract] OR shedding[Title/Abstract] OR '
                'tracheobronchitis[Title/Abstract]) AND '
                '(systematic review[Publication Type] OR review[Publication Type])'
            ),
        }
    elif category == "Strict anaerobe/aspiration flora":
        tracks[2] = {
            "track": "guideline_or_review",
            "scope": "category_context",
            "query": (
                '(anaerobic bacteria[Title/Abstract] OR anaerobes[Title/Abstract] OR '
                'Bacteroides[Title/Abstract]) AND (aspiration pneumonia[Title/Abstract] OR '
                'lung abscess[Title/Abstract] OR empyema[Title/Abstract] OR '
                'necrotizing pneumonia[Title/Abstract]) AND '
                '(guideline[Publication Type] OR systematic review[Publication Type] OR review[Publication Type])'
            ),
        }
    if category in {"Candida/yeast", "Herpesvirus/reactivation", "Skin/airway colonizer Gram-positive"}:
        tracks.append(
            {
                "track": "high_risk_host_or_icu",
                "scope": "species_specific",
                "query": (
                    f"(({species_title} AND {pulmonary}) OR "
                    f"({species} AND (pneumonia[Title] OR pulmonary[Title] OR bronchoalveolar[Title]))) AND "
                    "(immunocompromised[Title/Abstract] OR ICU[Title/Abstract] OR "
                    "critical illness[Title/Abstract] OR ventilator[Title/Abstract])"
                ),
            }
        )
    elif category == "Strict anaerobe/aspiration flora":
        tracks.append(
            {
                "track": "syndrome_specific",
                "scope": "species_specific",
                "query": (
                    f"(({species_title} AND (aspiration pneumonia[Title/Abstract] OR "
                    "lung abscess[Title/Abstract] OR empyema[Title/Abstract] OR "
                    "necrotizing pneumonia[Title/Abstract])) OR "
                    f"({species} AND (aspiration pneumonia[Title] OR lung abscess[Title] OR "
                    "empyema[Title] OR necrotizing pneumonia[Title])))"
                ),
            }
        )
    return tracks


@dataclass
class RateLimiter:
    minimum_interval: float
    last_request_started: float = 0.0

    def wait(self) -> None:
        elapsed = time.monotonic() - self.last_request_started
        if elapsed < self.minimum_interval:
            time.sleep(self.minimum_interval - elapsed)
        self.last_request_started = time.monotonic()


class PubMedClient:
    def __init__(
        self,
        *,
        cache_dir: Path,
        email: str | None,
        api_key: str | None,
        timeout: float,
        retries: int,
        refresh_cache: bool,
    ) -> None:
        self.cache_dir = cache_dir
        self.email = clean_text(email)
        self.api_key = clean_text(api_key)
        self.timeout = timeout
        self.retries = max(retries, 1)
        self.refresh_cache = refresh_cache
        self.rate_limiter = RateLimiter(0.11 if self.api_key else 0.36)
        self.network_request_count = 0
        self.cache_hit_count = 0

    def _params(self, values: dict[str, Any]) -> dict[str, Any]:
        params: dict[str, Any] = {"db": "pubmed", "tool": TOOL_NAME, **values}
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        return params

    def _get(self, endpoint: str, params: dict[str, Any], extension: str) -> bytes:
        public_params = {key: value for key, value in params.items() if key != "api_key"}
        fingerprint = hashlib.sha256(
            json.dumps([endpoint, public_params], sort_keys=True, ensure_ascii=True).encode("utf-8")
        ).hexdigest()
        cache_path = self.cache_dir / endpoint.replace(".fcgi", "") / f"{fingerprint}.{extension}"
        if cache_path.exists() and not self.refresh_cache:
            self.cache_hit_count += 1
            return cache_path.read_bytes()

        url = f"{NCBI_BASE}/{endpoint}?{urlencode(params)}"
        last_error: Exception | None = None
        for attempt in range(1, self.retries + 1):
            try:
                self.rate_limiter.wait()
                request = Request(url, headers={"User-Agent": f"{TOOL_NAME}/2.0"})
                with urlopen(request, timeout=self.timeout) as response:  # noqa: S310
                    payload = response.read()
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                cache_path.write_bytes(payload)
                self.network_request_count += 1
                return payload
            except Exception as exc:  # noqa: BLE001
                last_error = exc
                if attempt < self.retries:
                    time.sleep(0.8 * attempt)
        raise RuntimeError(f"NCBI {endpoint} failed after {self.retries} attempts: {last_error}")

    def search(self, query: str, max_results: int) -> tuple[int, list[str]]:
        payload = self._get(
            "esearch.fcgi",
            self._params(
                {
                    "term": query,
                    "retmode": "json",
                    "retmax": max_results,
                    "sort": "relevance",
                }
            ),
            "json",
        )
        data = json.loads(payload.decode("utf-8"))
        result = data.get("esearchresult") or {}
        return int(result.get("count") or 0), [str(value) for value in result.get("idlist") or []]

    def fetch(self, pmids: Sequence[str]) -> list[dict[str, Any]]:
        if not pmids:
            return []
        payload = self._get(
            "efetch.fcgi",
            self._params({"id": ",".join(pmids), "retmode": "xml"}),
            "xml",
        )
        return parse_pubmed_xml(payload)


def parse_pubmed_xml(payload: bytes | str) -> list[dict[str, Any]]:
    root = ET.fromstring(payload)
    records: list[dict[str, Any]] = []
    for article_node in root.findall(".//PubmedArticle"):
        citation = article_node.find("MedlineCitation")
        article = citation.find("Article") if citation is not None else None
        if citation is None or article is None:
            continue
        pmid = xml_text(citation.find("PMID"))
        title = xml_text(article.find("ArticleTitle"))
        abstract_parts: list[str] = []
        for abstract_node in article.findall(".//Abstract/AbstractText"):
            label = clean_text(abstract_node.attrib.get("Label"))
            text = xml_text(abstract_node)
            if text:
                abstract_parts.append(f"{label}: {text}" if label else text)
        journal = xml_text(article.find("Journal/Title"))
        year_text = (
            xml_text(article.find("Journal/JournalIssue/PubDate/Year"))
            or xml_text(article.find("ArticleDate/Year"))
        )
        year_match = re.search(r"\b(19|20)\d{2}\b", year_text)
        year = int(year_match.group(0)) if year_match else None
        publication_types = dedupe(
            [xml_text(node) for node in article.findall("PublicationTypeList/PublicationType")]
        )
        mesh_terms = dedupe(
            [xml_text(node.find("DescriptorName")) for node in citation.findall("MeshHeadingList/MeshHeading")]
        )
        doi = ""
        pmc_id = ""
        for id_node in article_node.findall(".//PubmedData/ArticleIdList/ArticleId"):
            id_type = clean_text(id_node.attrib.get("IdType")).lower()
            value = xml_text(id_node)
            if id_type == "doi":
                doi = value
            elif id_type == "pmc":
                pmc_id = value
        records.append(
            {
                "pmid": pmid,
                "title": title,
                "abstract": " ".join(abstract_parts),
                "journal": journal,
                "year": year,
                "doi": doi,
                "pmc_id": pmc_id,
                "publication_types": publication_types,
                "mesh_terms": mesh_terms,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
            }
        )
    return records


def sentences(text: str) -> list[str]:
    return [piece.strip() for piece in re.split(r"(?<=[.!?])\s+|[;\n]", text) if piece.strip()]


def article_score(
    article: dict[str, Any],
    organism_name: str,
    category: str,
    tracks: set[str],
) -> tuple[float, list[str]]:
    title = clean_text(article.get("title")).lower()
    abstract = clean_text(article.get("abstract")).lower()
    text = title + "\n" + abstract
    name = clean_text(organism_name).lower()
    genus = name.split(" ", 1)[0] if name else ""
    score = 0.0
    reasons: list[str] = []
    exact_in_title = bool(name and name in title)
    exact_in_abstract = bool(name and name in abstract)
    genus_in_title = bool(genus and genus in title)
    genus_in_abstract = bool(genus and genus in abstract)
    pulmonary_in_title = any(term in title for term in PULMONARY_TERMS)
    pulmonary_in_abstract = any(term in abstract for term in PULMONARY_TERMS)
    background_in_title = any(term in title for term in BACKGROUND_TERMS)
    background_in_abstract = any(term in abstract for term in BACKGROUND_TERMS)
    same_sentence_pulmonary = any(
        (name in sentence or genus in sentence) and any(term in sentence for term in PULMONARY_TERMS)
        for sentence in sentences(text)
        if name or genus
    )
    same_sentence_background = any(
        (name in sentence or genus in sentence) and any(term in sentence for term in BACKGROUND_TERMS)
        for sentence in sentences(text)
        if name or genus
    )

    if exact_in_title:
        score += 14.0
        reasons.append("species_name_in_title")
    elif exact_in_abstract:
        score += 7.0
        reasons.append("species_name_present")
    elif genus_in_title:
        score += 5.0
        reasons.append("genus_name_in_title")
    elif genus_in_abstract:
        score += 3.0
        reasons.append("genus_name_present")
    if pulmonary_in_title:
        score += 9.0
        reasons.append("pulmonary_context_in_title")
    elif pulmonary_in_abstract:
        score += 4.0
        reasons.append("pulmonary_context_present")
    if background_in_title:
        score += 5.0
        reasons.append("colonization_or_background_in_title")
    elif background_in_abstract:
        score += 2.0
        reasons.append("colonization_or_background_context_present")
    if same_sentence_pulmonary:
        score += 10.0
        reasons.append("organism_and_pulmonary_same_sentence")
    if same_sentence_background:
        score += 7.0
        reasons.append("organism_and_background_same_sentence")
    for publication_type in article.get("publication_types") or []:
        weight = HIGHER_EVIDENCE_PUBLICATION_TYPES.get(str(publication_type), 0.0)
        if weight:
            score += weight
            reasons.append(f"publication_type:{publication_type}")
    if article.get("abstract"):
        score += 1.0
    if "pulmonary_pathogenicity" in tracks and (pulmonary_in_title or same_sentence_pulmonary):
        score += 2.0
    if "colonization_or_background" in tracks and (background_in_title or same_sentence_background):
        score += 2.0
    title_mismatch_terms = (
        "alzheimer",
        "dementia",
        "idiopathic pulmonary fibrosis",
        "oral microbiota",
        "periodontitis",
        "gingivitis",
        "intestinal",
        "gut microbi",
        "phage",
        "oxidative stress",
        "drug review",
        "ceftaroline",
        "doripenem",
        "fournier",
        "gangrene",
        "ocular",
        "corneal",
        "depression",
        "metabolic syndrome",
        "cancer immunity",
        "endocarditis",
        "meningitis",
        "denture",
        "transplantation for the treatment of malnutrition",
        "vaccination",
        "pulmonary embolism",
        "pulmonary fibrosis",
        "post-covid",
    )
    if any(term in title for term in title_mismatch_terms):
        score -= 16.0
        reasons.append("topic_mismatch_penalty")
    if any(term in title for term in ("children", "pediatric", "paediatric", "neonatal")):
        score -= 7.0
        reasons.append("non_adult_population_penalty")
    if any(term in title for term in ("animal model", "experimental", "rabbit", "mice", "mouse", "rat")):
        score -= 16.0
        reasons.append("nonhuman_or_experimental_penalty")
    if not same_sentence_pulmonary and not pulmonary_in_title:
        score -= 5.0
        reasons.append("weak_organism_pulmonary_link_penalty")
    if "SARS-CoV-2" in organism_name and "ventilator-associated pneumonia" in title:
        score -= 8.0
        reasons.append("covid_host_vap_not_viral_causality_penalty")
    if category == "Strict anaerobe/aspiration flora" and not any(
        term in text for term in ("aspiration", "abscess", "empyema", "necrotizing", "anaerobic pneumonia")
    ):
        score -= 5.0
        reasons.append("missing_anaerobic_pulmonary_syndrome_penalty")
    return score, dedupe(reasons)


def article_matches_track(
    article: dict[str, Any],
    organism_name: str,
    category: str,
    track: str,
    scope: str,
) -> bool:
    title = clean_text(article.get("title")).lower()
    abstract = clean_text(article.get("abstract")).lower()
    text = title + "\n" + abstract
    name = clean_text(organism_name).lower()
    genus = name.split(" ", 1)[0] if name else ""
    organism_in_title = bool(name and name in title)
    organism_in_abstract = bool(name and name in abstract)
    organism_present = organism_in_title or organism_in_abstract
    pulmonary_in_title = any(term in title for term in PULMONARY_TERMS)
    pulmonary_present = pulmonary_in_title or any(term in abstract for term in PULMONARY_TERMS)
    background_in_title = any(term in title for term in BACKGROUND_TERMS)
    background_present = background_in_title or any(term in abstract for term in BACKGROUND_TERMS)
    same_sentence_pulmonary = any(
        (name in sentence or genus in sentence) and any(term in sentence for term in PULMONARY_TERMS)
        for sentence in sentences(text)
        if name or genus
    )
    same_sentence_background = any(
        (name in sentence or genus in sentence) and any(term in sentence for term in BACKGROUND_TERMS)
        for sentence in sentences(text)
        if name or genus
    )
    publication_types = {str(value) for value in article.get("publication_types") or []}
    mesh_terms = {str(value).lower() for value in article.get("mesh_terms") or []}
    is_review = bool(
        publication_types.intersection({"Practice Guideline", "Guideline", "Systematic Review", "Meta-Analysis", "Review"})
    )
    category_present = {
        "Candida/yeast": "candida" in text,
        "Herpesvirus/reactivation": "herpes simplex" in text or "hsv" in text,
        "Strict anaerobe/aspiration flora": "anaerob" in text or "bacteroides" in text,
    }.get(category, bool(genus and genus in text))
    topic_mismatch = any(
        term in title
        for term in (
            "alzheimer",
            "dementia",
            "oral microbiota",
            "periodontitis",
            "gingivitis",
            "denture",
            "oral care",
            "endocarditis",
            "meningitis",
            "fournier",
            "gangrene",
            "vaccination",
            "pulmonary embolism",
        )
    )
    if topic_mismatch:
        return False
    if any(term in title for term in ("children", "pediatric", "paediatric", "neonatal")):
        return False
    if mesh_terms.intersection({"animals", "mice", "rats", "rabbits"}) and "humans" not in mesh_terms:
        return False
    if track == "pulmonary_pathogenicity":
        return organism_present and (pulmonary_in_title or same_sentence_pulmonary)
    if track == "colonization_or_background":
        return (
            organism_present
            and (pulmonary_in_title or same_sentence_pulmonary)
            and (background_in_title or same_sentence_background)
        )
    if track == "guideline_or_review":
        scope_match = category_present if scope == "category_context" else organism_present
        return scope_match and pulmonary_present and is_review
    if track == "high_risk_host_or_icu":
        high_risk_present = any(
            term in text for term in ("immunocompromised", "critical illness", "critically ill", "icu", "ventilator")
        )
        return organism_present and pulmonary_present and high_risk_present
    if track == "syndrome_specific":
        syndrome_present = any(
            term in text for term in ("aspiration", "lung abscess", "empyema", "necrotizing pneumonia")
        )
        return organism_present and syndrome_present and (pulmonary_in_title or same_sentence_pulmonary)
    return False


def select_diverse_articles(
    articles: Sequence[dict[str, Any]],
    *,
    max_articles: int,
) -> list[dict[str, Any]]:
    quotas = {
        "pulmonary_pathogenicity": 5,
        "colonization_or_background": 4,
        "guideline_or_review": 3,
        "high_risk_host_or_icu": 2,
        "syndrome_specific": 3,
    }
    minimum_scores = {
        "pulmonary_pathogenicity": 20.0,
        "colonization_or_background": 18.0,
        "guideline_or_review": 16.0,
        "high_risk_host_or_icu": 18.0,
        "syndrome_specific": 18.0,
    }

    def qualifies(article: dict[str, Any], track: str) -> bool:
        return float(article.get("retrieval_relevance_score") or 0) >= minimum_scores[track]

    selected: list[dict[str, Any]] = []
    selected_pmids: set[str] = set()
    for track, quota in quotas.items():
        candidates = [article for article in articles if track in article.get("retrieval_tracks", [])]
        for article in candidates:
            pmid = str(article.get("pmid") or "")
            if not pmid or pmid in selected_pmids:
                continue
            if not qualifies(article, track):
                continue
            selected.append(article)
            selected_pmids.add(pmid)
            if sum(track in item.get("retrieval_tracks", []) for item in selected) >= quota:
                break
            if len(selected) >= max_articles:
                return selected
    for article in articles:
        pmid = str(article.get("pmid") or "")
        if not pmid or pmid in selected_pmids:
            continue
        if not any(
            qualifies(article, track)
            for track in article.get("retrieval_tracks", [])
            if track in minimum_scores
        ):
            continue
        selected.append(article)
        selected_pmids.add(pmid)
        if len(selected) >= max_articles:
            break
    return selected


def retrieve_request(
    request: dict[str, Any],
    client: PubMedClient,
    *,
    max_results_per_query: int,
    max_articles: int,
) -> dict[str, Any]:
    tracks = query_tracks(request)
    pmid_tracks: dict[str, set[str]] = {}
    track_scopes = {track["track"]: track["scope"] for track in tracks}
    track_results: list[dict[str, Any]] = []
    fetched: dict[str, dict[str, Any]] = {}
    errors: list[dict[str, str]] = []
    for track in tracks:
        try:
            count, pmids = client.search(track["query"], max_results_per_query)
            articles = client.fetch(pmids)
            track_results.append(
                {
                    "track": track["track"],
                    "scope": track["scope"],
                    "query": track["query"],
                    "pubmed_total_count": count,
                    "retrieved_pmids": pmids,
                }
            )
            for article in articles:
                pmid = clean_text(article.get("pmid"))
                if not pmid:
                    continue
                fetched.setdefault(pmid, article)
                pmid_tracks.setdefault(pmid, set()).add(track["track"])
        except Exception as exc:  # noqa: BLE001
            errors.append({"track": track["track"], "error": str(exc)})

    organism = request.get("organism") if isinstance(request.get("organism"), dict) else {}
    organism_name = clean_text(organism.get("display_name"))
    category = clean_text(organism.get("pathogen_category"))
    ranked: list[dict[str, Any]] = []
    for pmid, article in fetched.items():
        article_tracks = {
            track
            for track in pmid_tracks.get(pmid, set())
            if article_matches_track(
                article,
                organism_name,
                category,
                track,
                track_scopes.get(track, "species_specific"),
            )
        }
        if not article_tracks:
            continue
        score, reasons = article_score(article, organism_name, category, article_tracks)
        ranked.append(
            {
                **article,
                "retrieval_tracks": sorted(article_tracks),
                "retrieval_relevance_score": score,
                "retrieval_relevance_reasons": reasons,
                "evidence_direction": "unassessed",
                "quality": "unassessed",
                "case_similarity": "unassessed",
            }
        )
    ranked.sort(
        key=lambda article: (
            -float(article["retrieval_relevance_score"]),
            -(int(article.get("year") or 0)),
            str(article.get("pmid") or ""),
        )
    )
    selected = select_diverse_articles(ranked, max_articles=max_articles)
    selected_track_counts = {
        track: sum(track in article.get("retrieval_tracks", []) for article in selected)
        for track in track_scopes
    }
    missing_tracks = [track for track, count in selected_track_counts.items() if count == 0]
    evidence_volume = "adequate" if len(selected) >= 10 else "limited" if len(selected) >= 5 else "sparse"
    return {
        "schema_version": "rag_pubmed_retrieval_v2.1",
        "knowledge_card_key": request.get("knowledge_card_key"),
        "organism": organism,
        "clinical_domain": request.get("clinical_domain"),
        "retrieved_at": utc_now(),
        "query_tracks": track_results,
        "retrieval_errors": errors,
        "unique_article_count_before_limit": len(ranked),
        "selected_article_count": len(selected),
        "selected_track_counts": selected_track_counts,
        "retrieval_evidence_volume": evidence_volume,
        "missing_evidence_tracks": missing_tracks,
        "articles": selected,
        "selection_policy": (
            "Sentence-level organism/pulmonary or organism/background linkage, adult pulmonary topic penalties, "
            "human/adult eligibility gates, track-specific minimum scores, and per-track diversity quotas. "
            "The article count is an upper bound; sparse evidence is reported instead of padding to the limit."
        ),
        "adjudication_note": (
            "Retrieval relevance is not a medical conclusion. evidence_direction, quality, and "
            "case_similarity remain unassessed until the adjudication stage."
        ),
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("knowledge_requests", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, help="Shared NCBI response cache; defaults to output-dir/cache.")
    parser.add_argument("--keys", nargs="*", help="Optional knowledge_card_key allowlist.")
    parser.add_argument("--max-results-per-query", type=int, default=12)
    parser.add_argument("--max-articles", type=int, default=15)
    parser.add_argument("--email", default=os.getenv("PUBMED_EMAIL"))
    parser.add_argument("--api-key", default=os.getenv("PUBMED_API_KEY"))
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    requests = read_jsonl(args.knowledge_requests)
    allowlist = {clean_text(value) for value in (args.keys or []) if clean_text(value)}
    if allowlist:
        requests = [
            request for request in requests if clean_text(request.get("knowledge_card_key")) in allowlist
        ]
        missing = sorted(allowlist - {clean_text(request.get("knowledge_card_key")) for request in requests})
        if missing:
            raise ValueError("Unknown knowledge_card_key values: " + ", ".join(missing))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plans = [
        {
            "knowledge_card_key": request.get("knowledge_card_key"),
            "organism": request.get("organism"),
            "query_tracks": query_tracks(request),
        }
        for request in requests
    ]
    write_jsonl(args.output_dir / "retrieval_plans.jsonl", plans)
    if args.dry_run:
        print(f"dry_run_plans={len(plans)}")
        return 0

    cache_dir = args.cache_dir or (args.output_dir / "cache")
    client = PubMedClient(
        cache_dir=cache_dir,
        email=args.email,
        api_key=args.api_key,
        timeout=args.timeout,
        retries=args.retries,
        refresh_cache=args.refresh_cache,
    )
    results: list[dict[str, Any]] = []
    for index, request in enumerate(requests, 1):
        result = retrieve_request(
            request,
            client,
            max_results_per_query=max(args.max_results_per_query, 1),
            max_articles=max(args.max_articles, 1),
        )
        results.append(result)
        print(
            f"[{index}/{len(requests)}] {result['knowledge_card_key']} "
            f"articles={result['selected_article_count']} errors={len(result['retrieval_errors'])}"
        )

    write_jsonl(args.output_dir / "pubmed_retrieval_results.jsonl", results)
    summary_rows = [
        {
            "knowledge_card_key": result["knowledge_card_key"],
            "organism_name": (result.get("organism") or {}).get("display_name"),
            "selected_article_count": result["selected_article_count"],
            "unique_article_count_before_limit": result["unique_article_count_before_limit"],
            "query_track_count": len(result["query_tracks"]),
            "retrieval_error_count": len(result["retrieval_errors"]),
            "top_pmids": "; ".join(article["pmid"] for article in result["articles"][:5]),
        }
        for result in results
    ]
    write_csv(
        args.output_dir / "pubmed_retrieval_summary.csv",
        summary_rows,
        (
            "knowledge_card_key",
            "organism_name",
            "selected_article_count",
            "unique_article_count_before_limit",
            "query_track_count",
            "retrieval_error_count",
            "top_pmids",
        ),
    )
    manifest = {
        "schema_version": "rag_pubmed_retrieval_manifest_v2.1",
        "created_at": utc_now(),
        "knowledge_requests": str(args.knowledge_requests),
        "cache_dir": str(cache_dir),
        "request_count": len(requests),
        "completed_count": len(results),
        "requests_with_errors": sum(bool(result["retrieval_errors"]) for result in results),
        "network_request_count": client.network_request_count,
        "cache_hit_count": client.cache_hit_count,
        "api_key_used": bool(args.api_key),
        "email_supplied": bool(args.email),
        "answer_source_read": False,
        "medical_adjudication_performed": False,
        "ncbi_documentation": [
            "https://www.ncbi.nlm.nih.gov/books/NBK25497/",
            "https://www.ncbi.nlm.nih.gov/books/NBK25499/",
        ],
    }
    write_json(args.output_dir / "retrieval_manifest.json", manifest)
    print(f"completed={len(results)}")
    print(f"requests_with_errors={manifest['requests_with_errors']}")
    print(f"network_requests={client.network_request_count}")
    print(f"cache_hits={client.cache_hit_count}")
    return 0 if not manifest["requests_with_errors"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

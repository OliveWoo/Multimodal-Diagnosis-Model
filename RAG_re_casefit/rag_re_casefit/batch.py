from __future__ import annotations

import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

from .case_card import build_case_card, find_candidate
from .judge import CaseFitJudge, SYSTEM_PROMPT_SHA256
from .rules import aggregate_judgments


WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
RAG_RE_ROOT = WORKSPACE_ROOT / "RAG_re"
if str(RAG_RE_ROOT) not in sys.path:
    sys.path.insert(0, str(RAG_RE_ROOT))

from rag_re.judge import ArticleJudge  # noqa: E402
from rag_re.pubmed import PubMedClient  # noqa: E402
from rag_re.rules import module_a_from_judgments  # noqa: E402


OUTPUT_SCHEMA = "rag_re_casefit.candidate_output.v1"
MANIFEST_SCHEMA = "rag_re_casefit.batch_manifest.v1"
TIERS = ("review_high_priority", "review_context_needed")


def _compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def _normal(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", _compact(value).casefold())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")
    temporary.replace(path)


def _patient_id(root: Dict[str, Any], path: Path) -> str:
    value = _compact(root.get("patient_id"))
    if value:
        return value
    match = re.search(r"patient[_-]?(\d+)", path.name, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"cannot identify patient from {path}")
    return match.group(1)


def _index_by_patient(directory: Path) -> Dict[str, Path]:
    result: Dict[str, Path] = {}
    for path in sorted(directory.glob("*.json")):
        if "manifest" in path.name.casefold():
            continue
        root = _read(path)
        patient_id = _patient_id(root, path)
        if patient_id in result:
            raise ValueError(f"duplicate patient {patient_id} in {directory}")
        result[patient_id] = path
    return result


def _review_rows(root: Dict[str, Any]) -> Iterable[Tuple[str, Dict[str, Any]]]:
    review = root.get("llm_missed_candidate_review")
    if not isinstance(review, dict):
        raise ValueError("latest merge is missing llm_missed_candidate_review")
    for tier in TIERS:
        rows = review.get(tier) or []
        if not isinstance(rows, list):
            raise ValueError(f"{tier} must be a list")
        for row in rows:
            if not isinstance(row, dict) or not _compact(row.get("organism_name")):
                raise ValueError(f"invalid {tier} row")
            yield tier, row


def _old_candidate(root: Dict[str, Any], organism: str) -> Dict[str, Any]:
    candidates = root.get("candidates")
    if not isinstance(candidates, list):
        raise ValueError("frozen RAG output has no candidates")
    aliases = {
        "cmv": "humancytomegalovirus",
        "hsv1": "herpessimplexvirustype1",
        "vzv": "varicellazostervirus",
    }
    target = aliases.get(_normal(organism), _normal(organism))
    matches = []
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        names = {
            aliases.get(_normal(candidate.get("organism_name")), _normal(candidate.get("organism_name"))),
            aliases.get(_normal(candidate.get("canonical_organism_name")), _normal(candidate.get("canonical_organism_name"))),
        }
        if target in names:
            matches.append(candidate)
    if len(matches) != 1:
        raise ValueError(f"expected one frozen candidate for {organism}, found {len(matches)}")
    return matches[0]


def _canonical_organism(organism: str, rag_re_config: Dict[str, Any]) -> str:
    aliases = rag_re_config.get("organism_normalization", {}).get("accepted_synonyms", {})
    normalized = {_normal(key): _compact(value) for key, value in aliases.items()}
    return normalized.get(_normal(organism), organism)


def _literature_context(card: Dict[str, Any]) -> Dict[str, str]:
    site = card.get("clinical_site")
    values = {
        "lower_respiratory": (
            "lower respiratory tract infection or pneumonia",
            "pneumonia OR pulmonary infection OR lower respiratory tract infection",
        ),
        "upper_respiratory": (
            "upper respiratory tract infection",
            "upper respiratory tract infection",
        ),
        "bloodstream": (
            "bloodstream infection, bacteremia, or sepsis",
            "bloodstream infection OR bacteremia OR sepsis",
        ),
        "cns": (
            "central nervous system infection or meningitis",
            "central nervous system infection OR meningitis",
        ),
        "urinary": ("urinary tract infection", "urinary tract infection"),
        "intra_abdominal": ("intra-abdominal infection", "intra-abdominal infection"),
    }
    target, query = values.get(site, (
        "human infection at the reported specimen or clinical site",
        "human infection",
    ))
    return {"target_site": target, "query_terms": query}


def _articles(candidate: Dict[str, Any]) -> List[Dict[str, Any]]:
    evidence = candidate.get("literature_evidence")
    retrieval = evidence.get("retrieval") if isinstance(evidence, dict) else None
    articles = retrieval.get("articles") if isinstance(retrieval, dict) else None
    if not isinstance(articles, list) or not articles or any(not isinstance(row, dict) for row in articles):
        raise ValueError("frozen candidate has no usable retrieved articles")
    pmids = [_compact(row.get("pmid")) for row in articles]
    if len([value for value in pmids if value]) != len(set(value for value in pmids if value)):
        raise ValueError("frozen retrieval contains duplicate PMIDs")
    return articles


def _old_a(candidate: Dict[str, Any]) -> Dict[str, Any]:
    modules = candidate.get("modules")
    module = modules.get("A") if isinstance(modules, dict) else None
    if not isinstance(module, dict):
        raise ValueError("frozen candidate has no old Module A")
    verdicts = module.get("verdict_counts")
    return {
        "positive": module.get("positive"),
        "status": module.get("status"),
        "support_count": module.get("support_count"),
        "judgeable_count": module.get("judgeable_count"),
        "support_ratio": module.get("support_ratio"),
        "verdict_counts": dict(verdicts) if isinstance(verdicts, dict) else {},
    }


def _slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")
    return slug or "organism"


def run_batch(
    *,
    merge_dir: Path,
    frozen_rag_dir: Path,
    output_dir: Path,
    config: Dict[str, Any],
    cache_dir: Path,
    rag_re_config: Dict[str, Any],
    rag_re_cache_dir: Path,
) -> Dict[str, Any]:
    latest = _index_by_patient(merge_dir)
    frozen = _index_by_patient(frozen_rag_dir)
    llm_config = config.get("llm")
    aggregation = config.get("aggregation")
    if not isinstance(llm_config, dict) or not isinstance(aggregation, dict):
        raise ValueError("config requires llm and aggregation objects")
    judge = CaseFitJudge(llm_config, cache_dir=cache_dir)
    baseline_pubmed = PubMedClient(rag_re_config["pubmed"], cache_dir=rag_re_cache_dir)
    baseline_judge = ArticleJudge(rag_re_config["llm"], cache_dir=rag_re_cache_dir)
    outputs = []
    errors = []
    seen_keys = set()
    output_dir.mkdir(parents=True, exist_ok=True)

    for patient_id, merge_path in sorted(latest.items(), key=lambda item: int(item[0])):
        if patient_id not in frozen:
            errors.append({"patient_id": patient_id, "error": "missing frozen RAG output"})
            continue
        merge_bytes = merge_path.read_bytes()
        merge_root = json.loads(merge_bytes.decode("utf-8-sig"))
        frozen_path = frozen[patient_id]
        frozen_root = _read(frozen_path)
        for review_tier, review_row in _review_rows(merge_root):
            organism = _compact(review_row.get("organism_name"))
            key = (patient_id, _normal(organism))
            if key in seen_keys:
                raise ValueError(f"duplicate selected candidate: patient {patient_id} {organism}")
            seen_keys.add(key)
            try:
                candidate = find_candidate(merge_root, organism)
                card = build_case_card(merge_root, candidate, source_bytes=merge_bytes)
                baseline_source = "reused_frozen_rag_output"
                try:
                    old_candidate = _old_candidate(frozen_root, organism)
                    canonical = _compact(old_candidate.get("canonical_organism_name")) or organism
                    articles = _articles(old_candidate)
                    old_a = _old_a(old_candidate)
                except ValueError:
                    baseline_source = "narrow_fallback_same_frozen_query_and_old_A_judge"
                    canonical = _canonical_organism(organism, rag_re_config)
                    retrieval = baseline_pubmed.search(
                        canonical,
                        _literature_context(card),
                        max_articles=int(rag_re_config["pubmed"].get("max_articles", 10)),
                    )
                    articles = retrieval["articles"]
                    old_judgments = baseline_judge.judge_many(
                        canonical,
                        _literature_context(card)["target_site"],
                        articles,
                    )
                    old_a = module_a_from_judgments(
                        retrieval,
                        old_judgments,
                        rag_re_config["module_a_literature"],
                    )
                judgments = judge.judge_many(canonical, card, articles)
                aggregate = aggregate_judgments(judgments, aggregation)
                artifact = {
                    "schema_version": OUTPUT_SCHEMA,
                    "patient_id": patient_id,
                    "candidate_organism": organism,
                    "canonical_organism": canonical,
                    "review_tier": review_tier,
                    "case_card": card,
                    "old_A": old_a,
                    "A1_judgments": judgments,
                    "A1_aggregate": aggregate,
                    "provenance": {
                        "latest_merge_file": merge_path.name,
                        "latest_merge_sha256": hashlib.sha256(merge_bytes).hexdigest(),
                        "frozen_rag_file": frozen_path.name,
                        "frozen_rag_sha256": _sha256(frozen_path),
                        "prompt_sha256": SYSTEM_PROMPT_SHA256,
                        "model": llm_config.get("model"),
                        "baseline_article_source": baseline_source,
                    },
                }
                output_path = output_dir / f"p{patient_id}_{_slug(canonical)}_casefit.json"
                _write(output_path, artifact)
                outputs.append({
                    "patient_id": patient_id,
                    "organism": canonical,
                    "output_file": output_path.name,
                    "output_sha256": _sha256(output_path),
                    "tier": aggregate["casefit_evidence_tier"],
                })
            except Exception as exc:
                errors.append({"patient_id": patient_id, "organism": organism, "error": str(exc)})
                raise

    manifest = {
        "schema_version": MANIFEST_SCHEMA,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "merge_dir": str(merge_dir.resolve()),
        "frozen_rag_dir": str(frozen_rag_dir.resolve()),
        "candidate_scope": list(TIERS),
        "generation_is_gold_free": True,
        "output_count": len(outputs),
        "error_count": len(errors),
        "outputs": outputs,
        "errors": errors,
    }
    _write(output_dir / "manifest.json", manifest)
    return manifest

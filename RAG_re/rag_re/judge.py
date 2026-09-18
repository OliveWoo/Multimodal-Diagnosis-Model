from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any, Dict, List

from openai import OpenAI

from .io_utils import stable_hash, write_json_atomic


PROMPT_PATH = Path(__file__).resolve().parent.parent / "prompts" / "article_support_v1.txt"
SYSTEM_PROMPT = PROMPT_PATH.read_text(encoding="utf-8").strip()
SYSTEM_PROMPT_SHA256 = hashlib.sha256(SYSTEM_PROMPT.encode("utf-8")).hexdigest()

RESPONSE_FORMAT = {
    "type": "json_schema",
    "json_schema": {
        "name": "article_support_judgment",
        "strict": True,
        "schema": {
            "type": "object",
            "additionalProperties": False,
            "properties": {
                "verdict": {
                    "type": "string",
                    "enum": ["support", "against", "irrelevant", "unclear"],
                },
                "organism_scope": {
                    "type": "string",
                    "enum": ["exact_species", "genus_only", "other", "unclear"],
                },
                "human_clinical_evidence": {"type": "boolean"},
                "target_site_match": {"type": "boolean"},
                "evidence_span": {"type": "string"},
                "rationale": {"type": "string"},
            },
            "required": [
                "verdict",
                "organism_scope",
                "human_clinical_evidence",
                "target_site_match",
                "evidence_span",
                "rationale",
            ],
        },
    },
}
RESPONSE_FORMAT_SHA256 = stable_hash(RESPONSE_FORMAT)


def _compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def _bool(value: Any) -> bool:
    return value is True


def _usage_dict(response: Any) -> Dict[str, int] | None:
    usage = getattr(response, "usage", None)
    if usage is None:
        return None
    prompt = getattr(usage, "prompt_tokens", None)
    completion = getattr(usage, "completion_tokens", None)
    total = getattr(usage, "total_tokens", None)
    details = getattr(usage, "completion_tokens_details", None)
    reasoning = getattr(details, "reasoning_tokens", None)
    return {
        "prompt_tokens": int(prompt or 0),
        "completion_tokens": int(completion or 0),
        "reasoning_tokens": int(reasoning or 0),
        "total_tokens": int(total or 0),
    }


class ArticleJudge:
    def __init__(self, config: Dict[str, Any], cache_dir: Path | None = None):
        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            raise RuntimeError(
                "OPENAI_API_KEY is required for Module A. Use --skip-literature for "
                "B/C-only runs, or replay an existing evidence artifact."
            )
        self.config = config
        self.cache_dir = cache_dir
        self.client = OpenAI(api_key=api_key)

    def _cache_path(
        self, organism: str, target_site: str, article: Dict[str, Any]
    ) -> Path | None:
        if self.cache_dir is None:
            return None
        key = stable_hash(
            {
                "prompt_version": self.config.get("prompt_version"),
                "prompt_sha256": SYSTEM_PROMPT_SHA256,
                "model": self.config.get("model"),
                "reasoning_effort": self.config.get("reasoning_effort"),
                "temperature": self.config.get("temperature"),
                "max_output_tokens": self.config.get("max_output_tokens"),
                "response_format_sha256": RESPONSE_FORMAT_SHA256,
                "organism": organism,
                "target_site": target_site,
                "pmid": article.get("pmid"),
                "title": article.get("title"),
                "abstract": article.get("abstract"),
            }
        )
        return self.cache_dir / "judgments" / f"{key}.json"

    def judge_one(
        self, organism: str, target_site: str, article: Dict[str, Any]
    ) -> Dict[str, Any]:
        title = _compact(article.get("title"))
        abstract = _compact(article.get("abstract"))
        base = {
            "pmid": str(article.get("pmid") or ""),
            "article_url": article.get("url"),
            "organism": organism,
            "target_site": target_site,
            "model": self.config.get("model"),
            "reasoning_effort": self.config.get("reasoning_effort"),
            "max_completion_tokens": int(self.config.get("max_output_tokens", 1200)),
            "prompt_version": self.config.get("prompt_version"),
            "prompt_sha256": SYSTEM_PROMPT_SHA256,
            "response_format_sha256": RESPONSE_FORMAT_SHA256,
        }
        if not abstract:
            return {
                **base,
                "status": "no_abstract",
                "verdict": "unclear",
                "organism_scope": "unclear",
                "human_clinical_evidence": False,
                "target_site_match": False,
                "evidence_span": "",
                "evidence_span_valid": False,
                "rationale": "No abstract was available for a reproducible judgment.",
                "usage": None,
                "cache_hit": False,
            }

        cache_path = self._cache_path(organism, target_site, article)
        if cache_path is not None and cache_path.exists():
            cached = json.loads(cache_path.read_text(encoding="utf-8-sig"))
            cached["cache_hit"] = True
            return cached

        user_message = json.dumps(
            {
                "organism": organism,
                "target_site": target_site,
                "article": {"pmid": base["pmid"], "title": title, "abstract": abstract},
            },
            ensure_ascii=False,
        )
        last_error: Exception | None = None
        response = None
        for attempt in range(1, 3):
            try:
                model = str(self.config.get("model"))
                request: Dict[str, Any] = {
                    "model": model,
                    "max_completion_tokens": int(
                        self.config.get("max_output_tokens", 1200)
                    ),
                    "response_format": RESPONSE_FORMAT,
                    "messages": [
                        {"role": "developer", "content": SYSTEM_PROMPT},
                        {"role": "user", "content": user_message},
                    ],
                }
                temperature = self.config.get("temperature")
                if temperature is not None and not model.lower().startswith("gpt-5"):
                    request["temperature"] = float(temperature)
                reasoning_effort = self.config.get("reasoning_effort")
                if reasoning_effort:
                    request["reasoning_effort"] = str(reasoning_effort)
                response = self.client.chat.completions.create(**request)
                break
            except Exception as exc:  # SDK exception types vary across versions.
                last_error = exc
                if attempt < 2:
                    time.sleep(1.0)
        if response is None:
            raise RuntimeError(f"Article judge failed after retries: {last_error}")

        if not getattr(response, "choices", None):
            raise RuntimeError(f"Article judge returned no choice for PMID {base['pmid']}")
        choice = response.choices[0]
        message = choice.message
        refusal = getattr(message, "refusal", None)
        if refusal:
            raise RuntimeError(
                f"Article judge refused PMID {base['pmid']}: {_compact(refusal)[:200]}"
            )
        finish_reason = getattr(choice, "finish_reason", None)
        if finish_reason not in (None, "stop"):
            raise RuntimeError(
                f"Article judge did not finish cleanly for PMID {base['pmid']}: "
                f"finish_reason={finish_reason}"
            )
        content = message.content
        if not content:
            raise RuntimeError(
                f"Article judge returned empty content for PMID {base['pmid']}"
            )
        try:
            raw = json.loads(content)
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                f"Article judge returned invalid JSON for PMID {base['pmid']}"
            ) from exc
        judgment = self._normalize(raw, title, abstract)
        payload = {
            **base,
            **judgment,
            "status": "ok",
            "usage": _usage_dict(response),
            "cache_hit": False,
        }
        if cache_path is not None:
            write_json_atomic(cache_path, payload)
        return payload

    @staticmethod
    def _normalize(raw: Dict[str, Any], title: str, abstract: str) -> Dict[str, Any]:
        verdict = _compact(raw.get("verdict")).lower()
        if verdict not in {"support", "against", "irrelevant", "unclear"}:
            verdict = "unclear"
        scope = _compact(raw.get("organism_scope")).lower()
        if scope not in {"exact_species", "genus_only", "other", "unclear"}:
            scope = "unclear"
        evidence_span = _compact(raw.get("evidence_span"))
        if len(evidence_span.split()) > 20:
            evidence_span = ""
        article_blob = f"{title} {abstract}".casefold()
        span_valid = bool(evidence_span and evidence_span.casefold() in article_blob)
        return {
            "verdict": verdict,
            "organism_scope": scope,
            "human_clinical_evidence": _bool(raw.get("human_clinical_evidence")),
            "target_site_match": _bool(raw.get("target_site_match")),
            "evidence_span": evidence_span,
            "evidence_span_valid": span_valid,
            "rationale": _compact(raw.get("rationale"))[:500],
        }

    def judge_many(
        self, organism: str, target_site: str, articles: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        # Deliberately one call per article: each paper is an independent vote and is
        # cached independently for threshold sweeps/replay.
        return [self.judge_one(organism, target_site, article) for article in articles]

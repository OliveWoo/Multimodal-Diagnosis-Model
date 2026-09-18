from __future__ import annotations

import hashlib
import json
import os
import re
import time
from pathlib import Path
from typing import Any, Dict, List

from openai import OpenAI

from .case_card import prompt_case_card
from .rules import apply_traceability_gate, derive_article_verdict


PROMPT_PATH = Path(__file__).resolve().parent.parent / "prompts" / "article_casefit_v1.txt"
SYSTEM_PROMPT = PROMPT_PATH.read_text(encoding="utf-8").strip()
SYSTEM_PROMPT_SHA256 = hashlib.sha256(SYSTEM_PROMPT.encode("utf-8")).hexdigest()
REPAIR_PROMPT_PATH = Path(__file__).resolve().parent.parent / "prompts" / "evidence_span_repair_v1.txt"
REPAIR_SYSTEM_PROMPT = REPAIR_PROMPT_PATH.read_text(encoding="utf-8").strip()
REPAIR_PROMPT_SHA256 = hashlib.sha256(REPAIR_SYSTEM_PROMPT.encode("utf-8")).hexdigest()
ALIGNMENT_SCHEMA = {"type": "string", "enum": ["match", "partial", "mismatch", "unknown"]}
RESPONSE_FORMAT = {
    "type": "json_schema",
    "json_schema": {
        "name": "article_patient_casefit_extraction",
        "strict": True,
        "schema": {
            "type": "object",
            "additionalProperties": False,
            "properties": {
                "organism_scope": {
                    "type": "string",
                    "enum": ["exact_species", "genus_only", "other", "unclear"],
                },
                "human_clinical_context": {
                    "type": "string",
                    "enum": [
                        "causal_infection",
                        "colonization_or_contamination",
                        "reactivation_or_bystander",
                        "detection_only",
                        "nonhuman",
                        "unclear",
                    ],
                },
                "site_alignment": ALIGNMENT_SCHEMA,
                "syndrome_alignment": ALIGNMENT_SCHEMA,
                "specimen_alignment": ALIGNMENT_SCHEMA,
                "host_alignment": ALIGNMENT_SCHEMA,
                "phenotype_alignment": ALIGNMENT_SCHEMA,
                "temporal_alignment": ALIGNMENT_SCHEMA,
                "support_span": {"type": "string"},
                "counter_span": {"type": "string"},
                "rationale": {"type": "string"},
            },
            "required": [
                "organism_scope",
                "human_clinical_context",
                "site_alignment",
                "syndrome_alignment",
                "specimen_alignment",
                "host_alignment",
                "phenotype_alignment",
                "temporal_alignment",
                "support_span",
                "counter_span",
                "rationale",
            ],
        },
    },
}
RESPONSE_FORMAT_SHA256 = hashlib.sha256(
    json.dumps(RESPONSE_FORMAT, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()
REPAIR_RESPONSE_FORMAT = {
    "type": "json_schema",
    "json_schema": {
        "name": "casefit_evidence_span_repair",
        "strict": True,
        "schema": {
            "type": "object",
            "additionalProperties": False,
            "properties": {
                "support_span": {"type": "string"},
                "counter_span": {"type": "string"},
            },
            "required": ["support_span", "counter_span"],
        },
    },
}


def _compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def _stable_hash(value: Any) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _span(value: Any, article_blob: str) -> tuple[str, bool]:
    text = _compact(value)
    if len(text.split()) > 20:
        return "", False
    valid = bool(text and text.casefold() in article_blob.casefold())
    return (text if valid else ""), valid


def _trigger_window(title: str, abstract: str, terms: List[str]) -> tuple[str, bool]:
    """Extract a deterministic <=20-word trace around an already-frozen concept."""

    for source in (abstract, title):
        tokens = list(re.finditer(r"\S+", source))
        for index, token in enumerate(tokens):
            normalized = re.sub(r"[^a-z]+", "", token.group(0).casefold())
            if not any(term in normalized for term in terms):
                continue
            start_index = max(0, index - 5)
            end_index = min(len(tokens), index + 10)
            if end_index - start_index > 20:
                end_index = start_index + 20
            excerpt = source[tokens[start_index].start():tokens[end_index - 1].end()]
            return _span(excerpt, f"{title} {abstract}")
    return "", False


class CaseFitJudge:
    def __init__(self, config: Dict[str, Any], cache_dir: Path | None = None):
        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            raise RuntimeError("OPENAI_API_KEY is required for case-fit judgments")
        self.config = dict(config)
        self.cache_dir = cache_dir
        self.client = OpenAI(api_key=api_key)

    def _cache_path(self, organism: str, card: Dict[str, Any], article: Dict[str, Any]) -> Path | None:
        if self.cache_dir is None:
            return None
        key = _stable_hash(
            {
                "prompt_sha256": SYSTEM_PROMPT_SHA256,
                "response_format_sha256": RESPONSE_FORMAT_SHA256,
                "model": self.config.get("model"),
                "reasoning_effort": self.config.get("reasoning_effort"),
                "max_output_tokens": self.config.get("max_output_tokens"),
                "organism": organism,
                "case_card": prompt_case_card(card),
                "pmid": article.get("pmid"),
                "title": article.get("title"),
                "abstract": article.get("abstract"),
            }
        )
        return self.cache_dir / "casefit_judgments" / f"{key}.json"

    def _request_json(
        self,
        *,
        system_prompt: str,
        payload: Dict[str, Any],
        response_format: Dict[str, Any],
        max_output_tokens: int,
    ) -> Dict[str, Any]:
        model = str(self.config.get("model") or "gpt-5.6-luna")
        request: Dict[str, Any] = {
            "model": model,
            "max_completion_tokens": max_output_tokens,
            "response_format": response_format,
            "messages": [
                {"role": "developer", "content": system_prompt},
                {"role": "user", "content": json.dumps(payload, ensure_ascii=False)},
            ],
        }
        reasoning_effort = self.config.get("reasoning_effort")
        if reasoning_effort:
            request["reasoning_effort"] = str(reasoning_effort)
        temperature = self.config.get("temperature")
        if temperature is not None and not model.casefold().startswith("gpt-5"):
            request["temperature"] = float(temperature)
        response = None
        last_error: Exception | None = None
        for attempt in range(2):
            try:
                response = self.client.chat.completions.create(**request)
                break
            except Exception as exc:  # SDK exception classes vary by release.
                last_error = exc
                if attempt == 0:
                    time.sleep(1.0)
        if response is None:
            raise RuntimeError(f"case-fit judge failed after retries: {last_error}")
        if not getattr(response, "choices", None):
            raise RuntimeError("case-fit judge returned no choice")
        choice = response.choices[0]
        message = choice.message
        if getattr(message, "refusal", None):
            raise RuntimeError(f"case-fit judge refused: {_compact(message.refusal)[:200]}")
        if getattr(choice, "finish_reason", None) not in (None, "stop"):
            raise RuntimeError(f"case-fit judge finish_reason={choice.finish_reason}")
        if not message.content:
            raise RuntimeError("case-fit judge returned empty content")
        return json.loads(message.content)

    def _repair_spans_if_needed(
        self,
        *,
        organism: str,
        title: str,
        abstract: str,
        extraction: Dict[str, Any],
        decision: Dict[str, Any],
        support_span: str,
        support_valid: bool,
        counter_span: str,
        counter_valid: bool,
    ) -> Dict[str, Any]:
        need_support = decision.get("verdict") in {"strong_match", "partial_match"} and not support_valid
        need_counter = decision.get("counterevidence") is True and not counter_valid
        result = {
            "support_span": support_span,
            "support_span_valid": support_valid,
            "counter_span": counter_span,
            "counter_span_valid": counter_valid,
            "span_repair_attempted": False,
            "span_repair_method": "not_needed",
            "span_repair_prompt_sha256": REPAIR_PROMPT_SHA256,
        }
        if not need_support and not need_counter:
            return result
        if need_support:
            fallback, valid = _trigger_window(
                title,
                abstract,
                ["pneumonia", "infection", "infectious", "bacteremia", "sepsis", "meningitis"],
            )
            if valid:
                result["support_span"] = fallback
                result["support_span_valid"] = True
        if need_counter:
            human_context = extraction.get("human_clinical_context")
            terms = (
                ["colonization", "colonisation", "contamination", "contaminant"]
                if human_context == "colonization_or_contamination"
                else ["reactivation", "bystander", "shedding", "latent"]
            )
            fallback, valid = _trigger_window(title, abstract, terms)
            if valid:
                result["counter_span"] = fallback
                result["counter_span_valid"] = True
        result["span_repair_attempted"] = True
        result["span_repair_method"] = "deterministic_trigger_window"
        need_support = need_support and not result["support_span_valid"]
        need_counter = need_counter and not result["counter_span_valid"]
        if not need_support and not need_counter:
            return result
        repaired = self._request_json(
            system_prompt=REPAIR_SYSTEM_PROMPT,
            payload={
                "organism": organism,
                "need_support_span": need_support,
                "need_counter_span": need_counter,
                "frozen_extraction": extraction,
                "article": {"title": title, "abstract": abstract},
            },
            response_format=REPAIR_RESPONSE_FORMAT,
            max_output_tokens=500,
        )
        article_blob = f"{title} {abstract}"
        if need_support:
            result["support_span"], result["support_span_valid"] = _span(
                repaired.get("support_span"), article_blob
            )
        if need_counter:
            result["counter_span"], result["counter_span_valid"] = _span(
                repaired.get("counter_span"), article_blob
            )
        result["span_repair_method"] = "bounded_llm_repair"
        return result

    def judge_one(self, organism: str, card: Dict[str, Any], article: Dict[str, Any]) -> Dict[str, Any]:
        organism = _compact(organism)
        if not organism:
            raise ValueError("organism is required")
        prompt_card = prompt_case_card(card)
        title = _compact(article.get("title"))
        abstract = _compact(article.get("abstract"))
        pmid = _compact(article.get("pmid"))
        base = {
            "schema_version": "rag_re_casefit.article_judgment.v1",
            "pmid": pmid,
            "article_url": article.get("url"),
            "organism": organism,
            "case_card_hash": _stable_hash(prompt_card),
            "model": self.config.get("model"),
            "reasoning_effort": self.config.get("reasoning_effort"),
            "prompt_sha256": SYSTEM_PROMPT_SHA256,
            "response_format_sha256": RESPONSE_FORMAT_SHA256,
        }
        if not abstract:
            return {
                **base,
                "status": "no_abstract",
                "organism_scope": "unclear",
                "human_clinical_context": "unclear",
                **{key: "unknown" for key in (
                    "site_alignment", "syndrome_alignment", "specimen_alignment",
                    "host_alignment", "phenotype_alignment", "temporal_alignment"
                )},
                "support_span": "",
                "support_span_valid": False,
                "counter_span": "",
                "counter_span_valid": False,
                "rationale": "No abstract was available.",
                **derive_article_verdict({
                    "organism_scope": "unclear",
                    "human_clinical_context": "unclear",
                    **{key: "unknown" for key in (
                        "site_alignment", "syndrome_alignment", "specimen_alignment",
                        "host_alignment", "phenotype_alignment", "temporal_alignment"
                    )},
                }),
                "cache_hit": False,
                "usage": None,
            }

        cache_path = self._cache_path(organism, card, article)
        if cache_path is not None and cache_path.exists():
            cached = json.loads(cache_path.read_text(encoding="utf-8-sig"))
            extraction = {key: cached[key] for key in (
                "organism_scope", "human_clinical_context", "site_alignment",
                "syndrome_alignment", "specimen_alignment", "host_alignment",
                "phenotype_alignment", "temporal_alignment"
            )}
            # Cache freezes the expensive LLM extraction, not the deterministic
            # research rule. This lets a reviewed truth-table fix be replayed without
            # sending patient context to the API again.
            decision = derive_article_verdict(extraction)
            repaired = self._repair_spans_if_needed(
                organism=organism,
                title=title,
                abstract=abstract,
                extraction=extraction,
                decision=decision,
                support_span=_compact(cached.get("support_span")),
                support_valid=cached.get("support_span_valid") is True,
                counter_span=_compact(cached.get("counter_span")),
                counter_valid=cached.get("counter_span_valid") is True,
            )
            cached.update(repaired)
            cached.update(apply_traceability_gate(
                decision,
                support_span_valid=repaired["support_span_valid"],
                counter_span_valid=repaired["counter_span_valid"],
            ))
            cached["cache_hit"] = True
            if repaired["span_repair_attempted"]:
                cache_path.write_text(json.dumps(cached, ensure_ascii=False, indent=2), encoding="utf-8")
            return cached

        payload = {
            "organism": organism,
            "patient_case_card": prompt_card,
            "article": {"pmid": pmid, "title": title, "abstract": abstract},
        }
        raw = self._request_json(
            system_prompt=SYSTEM_PROMPT,
            payload=payload,
            response_format=RESPONSE_FORMAT,
            max_output_tokens=int(self.config.get("max_output_tokens", 1400)),
        )
        extraction = {key: raw[key] for key in (
            "organism_scope", "human_clinical_context", "site_alignment",
            "syndrome_alignment", "specimen_alignment", "host_alignment",
            "phenotype_alignment", "temporal_alignment"
        )}
        decision = derive_article_verdict(extraction)
        article_blob = f"{title} {abstract}"
        support_span, support_valid = _span(raw.get("support_span"), article_blob)
        counter_span, counter_valid = _span(raw.get("counter_span"), article_blob)
        repaired = self._repair_spans_if_needed(
            organism=organism,
            title=title,
            abstract=abstract,
            extraction=extraction,
            decision=decision,
            support_span=support_span,
            support_valid=support_valid,
            counter_span=counter_span,
            counter_valid=counter_valid,
        )
        decision = apply_traceability_gate(
            decision,
            support_span_valid=repaired["support_span_valid"],
            counter_span_valid=repaired["counter_span_valid"],
        )
        result = {
            **base,
            "status": "ok",
            **extraction,
            **repaired,
            "rationale": _compact(raw.get("rationale"))[:800],
            **decision,
            "cache_hit": False,
            "usage": None,
        }
        if cache_path is not None:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        return result

    def judge_many(self, organism: str, card: Dict[str, Any], articles: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        return [self.judge_one(organism, card, article) for article in articles]

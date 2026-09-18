from __future__ import annotations

import base64
from dataclasses import dataclass
import json
import mimetypes
import os
from pathlib import Path
import re
from typing import Protocol


DEFAULT_VISION_PROMPT = (
    "You are extracting visible content from a single medical screenshot that came from a PowerPoint slide. "
    "Return strict JSON with keys raw_text"  #, findings, confidence. "
    "raw_text must be an array of verbatim text blocks visible in the image. "
    #"findings must be an array of short clinical or report-style statements inferred directly from the visible content. "
    #"confidence must be a number between 0 and 1, or null if uncertain. "
    "Do not add commentary outside JSON."
)


@dataclass
class VisionExtractionResult:
    raw_text: list[str]
    findings: list[str]
    confidence: float | None = None
    model_used: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "raw_text": self.raw_text,
            "findings": self.findings,
            "confidence": self.confidence,
            "model_used": self.model_used,
        }


class VisionTextExtractor(Protocol):
    """Hook point for future OCR / LLM vision backends."""

    def extract(self, image_path: Path) -> VisionExtractionResult:
        ...


def _dedupe_preserve_order(values: list[str]) -> list[str]:
    deduped: list[str] = []
    seen: set[str] = set()
    for value in values:
        normalized = value.strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(normalized)
    return deduped


def _coerce_string_list(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return [str(item).strip() for item in value if str(item).strip()]


JSON_STRING_LITERAL_PATTERN = re.compile(r'"((?:\\.|[^"\\])*)"')


def _strip_json_markers(text: str) -> str:
    stripped = text.strip()
    if not stripped.startswith("```"):
        return stripped
    lines = stripped.splitlines()
    if len(lines) >= 3 and lines[-1].strip() == "```":
        body = "\n".join(lines[1:-1]).strip()
        return body if body else stripped
    return stripped


def _extract_output_text(response) -> str:
    if getattr(response, "output_text", None):
        text = response.output_text.strip()
        if text:
            return text

    segments: list[str] = []
    for item in getattr(response, "output", []) or []:
        for content in getattr(item, "content", []) or []:
            text = getattr(content, "text", None)
            if text:
                segments.append(str(text))
    return "".join(segments).strip()


def _has_refusal(response) -> bool:
    for item in getattr(response, "output", []) or []:
        for content in getattr(item, "content", []) or []:
            refusal = getattr(content, "refusal", None)
            if refusal:
                return True
    return False


def _try_parse_embedded_payload(text: str) -> dict[str, object] | None:
    candidate = _strip_json_markers(text)
    if not candidate or candidate[:1] not in "{[":
        return None
    try:
        payload = json.loads(candidate)
    except json.JSONDecodeError:
        return None
    if isinstance(payload, dict) and any(key in payload for key in ("raw_text", "findings", "confidence")):
        return payload
    return None


def _looks_like_embedded_payload(text: str) -> bool:
    candidate = _strip_json_markers(text)
    if not candidate.startswith("{"):
        return False
    return '"raw_text"' in candidate or '"findings"' in candidate or '"confidence"' in candidate


def _extract_strings_from_json_like_payload(text: str) -> list[str]:
    candidate = _strip_json_markers(text)
    if not _looks_like_embedded_payload(candidate):
        return []

    start_index = candidate.find('"raw_text"')
    if start_index >= 0:
        candidate = candidate[start_index:]

    extracted: list[str] = []
    for match in JSON_STRING_LITERAL_PATTERN.findall(candidate):
        try:
            value = json.loads(f'"{match}"')
        except json.JSONDecodeError:
            continue
        normalized = str(value).strip()
        if not normalized or normalized in {"raw_text", "findings", "confidence", "model_used"}:
            continue
        extracted.append(normalized)
    return _dedupe_preserve_order(extracted)


def _normalize_nested_payload(
    *,
    raw_text: list[str],
    findings: list[str],
    confidence: float | None,
) -> tuple[list[str], list[str], float | None]:
    normalized_raw_text: list[str] = []
    normalized_findings: list[str] = []
    pending_payloads: list[dict[str, object]] = []
    resolved_confidence = confidence

    def _consume(values: list[str], target: list[str]) -> None:
        for value in values:
            embedded_payload = _try_parse_embedded_payload(value)
            if embedded_payload is not None:
                pending_payloads.append(embedded_payload)
                continue
            loose_payload_values = _extract_strings_from_json_like_payload(value)
            if loose_payload_values:
                target.extend(loose_payload_values)
                continue
            target.append(value)

    _consume(raw_text, normalized_raw_text)
    _consume(findings, normalized_findings)

    while pending_payloads:
        payload = pending_payloads.pop(0)
        nested_confidence = payload.get("confidence")
        if resolved_confidence is None and isinstance(nested_confidence, (int, float)):
            resolved_confidence = float(nested_confidence)
        _consume(_coerce_string_list(payload.get("raw_text")), normalized_raw_text)
        _consume(_coerce_string_list(payload.get("findings")), normalized_findings)

    return (
        _dedupe_preserve_order(normalized_raw_text),
        _dedupe_preserve_order(normalized_findings),
        resolved_confidence,
    )


def _parse_vision_response(output_text: str) -> VisionExtractionResult:
    candidate = _strip_json_markers(output_text)
    try:
        payload = json.loads(candidate)
    except json.JSONDecodeError:
        normalized_text = candidate.strip()
        return VisionExtractionResult(
            raw_text=[normalized_text] if normalized_text else [],
            findings=[],
            confidence=None,
        )

    if isinstance(payload, str):
        nested_payload = _try_parse_embedded_payload(payload)
        if nested_payload is not None:
            payload = nested_payload
        else:
            normalized_text = payload.strip()
            return VisionExtractionResult(
                raw_text=[normalized_text] if normalized_text else [],
                findings=[],
                confidence=None,
            )

    if not isinstance(payload, dict):
        normalized_text = candidate.strip()
        return VisionExtractionResult(
            raw_text=[normalized_text] if normalized_text else [],
            findings=[],
            confidence=None,
        )

    confidence = payload.get("confidence")
    if isinstance(confidence, (int, float)):
        confidence_value: float | None = float(confidence)
    else:
        confidence_value = None

    raw_text, findings, confidence_value = _normalize_nested_payload(
        raw_text=_coerce_string_list(payload.get("raw_text")),
        findings=_coerce_string_list(payload.get("findings")),
        confidence=confidence_value,
    )

    return VisionExtractionResult(
        raw_text=raw_text,
        findings=findings,
        confidence=confidence_value,
    )


def _image_to_data_url(image_path: Path) -> str:
    mime_type = mimetypes.guess_type(image_path.name)[0] or "image/png"
    encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
    return f"data:{mime_type};base64,{encoded}"


class OpenAIVisionTextExtractor:
    """OpenAI-backed image text extractor used as an OCR/vision fallback."""

    def __init__(
        self,
        *,
        model: str,
        fallback_models: tuple[str, ...] = (),
        prompt: str = DEFAULT_VISION_PROMPT,
        detail: str = "high",
        max_output_tokens: int = 1200,
    ) -> None:
        self.model = model
        self.fallback_models = tuple(m for m in fallback_models if m and m != model)
        self.prompt = prompt
        self.detail = detail
        self.max_output_tokens = max_output_tokens
        self._client = None

    def _get_client(self):
        if self._client is not None:
            return self._client

        api_key = os.environ.get("OPENAI_API_KEY")
        if not api_key:
            raise RuntimeError("OPENAI_API_KEY environment variable is not set.")

        from openai import OpenAI  # type: ignore

        try:
            self._client = OpenAI(api_key=api_key)
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(
                "Failed to initialize OpenAI client. "
                "Check Python SSL/certificate setup if the traceback mentions ssl or certifi."
            ) from exc
        return self._client

    def _single_model_extract(self, *, client, model: str, image_url: str) -> VisionExtractionResult:
        response = client.responses.create(
            model=model,
            max_output_tokens=self.max_output_tokens,
            input=[
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "input_text",
                            "text": self.prompt,
                        },
                        {
                            "type": "input_image",
                            "image_url": image_url,
                            "detail": self.detail,
                        },
                    ],
                }
            ],
        )
        output_text = _extract_output_text(response)
        if output_text:
            parsed = _parse_vision_response(output_text)
            parsed.model_used = model
            return parsed

        retry_response = client.responses.create(
            model=model,
            max_output_tokens=self.max_output_tokens,
            input=[
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "input_text",
                            "text": (
                                "Extract visible text from this medical image. "
                                "Return strict JSON: {\"raw_text\": [string], \"findings\": [], \"confidence\": null}."
                            ),
                        },
                        {
                            "type": "input_image",
                            "image_url": image_url,
                            "detail": "auto",
                        },
                    ],
                }
            ],
        )
        retry_text = _extract_output_text(retry_response)
        if retry_text:
            parsed = _parse_vision_response(retry_text)
            parsed.model_used = model
            return parsed

        if _has_refusal(response) or _has_refusal(retry_response):
            return VisionExtractionResult(raw_text=[], findings=[], confidence=None, model_used=model)
        return VisionExtractionResult(raw_text=[], findings=[], confidence=None, model_used=model)

    def extract(self, image_path: Path) -> VisionExtractionResult:
        client = self._get_client()
        image_url = _image_to_data_url(image_path)
        candidates = (self.model,) + self.fallback_models
        for candidate_model in candidates:
            result = self._single_model_extract(client=client, model=candidate_model, image_url=image_url)
            if result.raw_text or result.findings:
                return result
        return VisionExtractionResult(raw_text=[], findings=[], confidence=None, model_used=self.model)

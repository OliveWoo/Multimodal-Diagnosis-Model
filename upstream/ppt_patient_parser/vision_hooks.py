from __future__ import annotations

import base64
from dataclasses import dataclass
import json
import mimetypes
import os
from pathlib import Path
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


def _coerce_string_list(value: object) -> list[str]:
    if not isinstance(value, list):
        return []
    return [str(item).strip() for item in value if str(item).strip()]


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

    confidence = payload.get("confidence")
    if isinstance(confidence, (int, float)):
        confidence_value: float | None = float(confidence)
    else:
        confidence_value = None

    return VisionExtractionResult(
        raw_text=_coerce_string_list(payload.get("raw_text")),
        findings=_coerce_string_list(payload.get("findings")),
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

from __future__ import annotations

import base64
import os
import unittest
from unittest.mock import patch
from pathlib import Path

from ppt_patient_parser.modality_classifier import classify_modalities
from ppt_patient_parser.parsers.slide_parser import parse_slide
from ppt_patient_parser.vision_hooks import OpenAIVisionTextExtractor


def _geometry(*, left: int, top: int, width: int, height: int, shape_order: int) -> dict[str, int]:
    return {
        "left": left,
        "top": top,
        "width": width,
        "height": height,
        "right": left + width,
        "bottom": top + height,
        "center_x": left + (width // 2),
        "center_y": top + (height // 2),
        "shape_order": shape_order,
    }


class ModalityClassifierTests(unittest.TestCase):
    def test_radiology_text_does_not_trigger_abg_or_smac(self) -> None:
        text_blocks = [
            (
                "Imaging findings: Multifocal mass-like and ground glass opacities over bilateral lungs. "
                "Normal heart size is noted. Please correlate with radiography follow up."
            )
        ]

        modalities, content_guess = classify_modalities(text_blocks)

        self.assertEqual(modalities, ["radiology"])
        self.assertFalse(content_guess["contains_abg_table"])
        self.assertFalse(content_guess["contains_smac_table"])


class SlideParserTests(unittest.TestCase):
    def test_parse_slide_uses_image_text_and_normalizes_dates(self) -> None:
        text_records = [
            {
                "text": "20240529",
                "shape_name": "caption",
                "geometry": _geometry(left=120, top=70, width=180, height=24, shape_order=1),
            }
        ]
        image_records = [
            {
                "image_index": 1,
                "raw_path": "raw_01.png",
                "fixed_path": "fixed_01.png",
                "width": 900,
                "height": 700,
                "preprocessing": {"suspicious_distortion": False},
                "display_geometry": _geometry(left=100, top=100, width=600, height=420, shape_order=2),
                "vision": {
                    "status": "success",
                    "raw_text": [
                        "TS-16/9041132",
                        "20240529 Imaging findings: ground glass opacity over bilateral lungs.",
                    ],
                    "findings": ["Ground glass opacity over bilateral lungs."],
                    "confidence": 0.91,
                },
            }
        ]

        parsed = parse_slide(
            source_file=Path("demo.pptx"),
            slide_index=1,
            text_blocks=[],
            text_records=text_records,
            image_records=image_records,
        )

        self.assertEqual(parsed["case_id"], "TS-16")
        self.assertEqual(parsed["patient_id"], "9041132")
        self.assertEqual(parsed["date"], "2024-05-29")
        self.assertIn("2024-05-29", parsed["detected_dates"])
        self.assertIn("radiology", parsed["candidate_modalities"])
        self.assertEqual(parsed["images"][0]["nearby_text_blocks"], ["20240529"])
        self.assertGreaterEqual(parsed["image_extraction_status_summary"]["success"], 1)


class VisionHookTests(unittest.TestCase):
    def test_openai_vision_extractor_reuses_client_and_omits_temperature(self) -> None:
        png_bytes = base64.b64decode(
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9sXl16sAAAAASUVORK5CYII="
        )
        captured_calls: list[dict] = []
        client_inits = 0

        class _FakeResponses:
            def create(self, **kwargs):
                captured_calls.append(kwargs)

                class _Response:
                    output_text = '{"raw_text":["TS-16/9041132"],"findings":["ground glass opacity"],"confidence":0.9}'

                return _Response()

        class _FakeClient:
            def __init__(self, api_key: str):
                nonlocal client_inits
                client_inits += 1
                self.api_key = api_key
                self.responses = _FakeResponses()

        temp_dir = Path("outputs") / "test_tmp"
        temp_dir.mkdir(parents=True, exist_ok=True)
        image_path = temp_dir / "tiny.png"
        try:
            image_path.write_bytes(png_bytes)

            with patch.dict(os.environ, {"OPENAI_API_KEY": "test-key"}):
                with patch("openai.OpenAI", _FakeClient):
                    extractor = OpenAIVisionTextExtractor(model="gpt-5")
                    first = extractor.extract(image_path)
                    second = extractor.extract(image_path)
        finally:
            if image_path.exists():
                image_path.unlink()

        self.assertEqual(first.raw_text, ["TS-16/9041132"])
        self.assertEqual(second.findings, ["ground glass opacity"])
        self.assertEqual(len(captured_calls), 2)
        self.assertNotIn("temperature", captured_calls[0])
        self.assertEqual(client_inits, 1)

    def test_openai_vision_extractor_empty_response_returns_empty_result(self) -> None:
        png_bytes = base64.b64decode(
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9sXl16sAAAAASUVORK5CYII="
        )
        calls = 0

        class _FakeResponses:
            def create(self, **kwargs):
                nonlocal calls
                calls += 1

                class _Response:
                    output_text = ""
                    output = []

                return _Response()

        class _FakeClient:
            def __init__(self, api_key: str):
                self.responses = _FakeResponses()

        temp_dir = Path("outputs") / "test_tmp"
        temp_dir.mkdir(parents=True, exist_ok=True)
        image_path = temp_dir / "tiny_empty.png"
        try:
            image_path.write_bytes(png_bytes)
            with patch.dict(os.environ, {"OPENAI_API_KEY": "test-key"}):
                with patch("openai.OpenAI", _FakeClient):
                    extractor = OpenAIVisionTextExtractor(model="gpt-5")
                    result = extractor.extract(image_path)
        finally:
            if image_path.exists():
                image_path.unlink()

        self.assertEqual(result.raw_text, [])
        self.assertEqual(result.findings, [])
        self.assertIsNone(result.confidence)
        self.assertEqual(calls, 2)

    def test_openai_vision_extractor_fallback_model_is_used(self) -> None:
        png_bytes = base64.b64decode(
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9sXl16sAAAAASUVORK5CYII="
        )
        called_models: list[str] = []

        class _FakeResponses:
            def create(self, **kwargs):
                model = kwargs["model"]
                called_models.append(model)

                class _Response:
                    output_text = ""
                    output = []

                if model == "gpt-4.1-mini":
                    _Response.output_text = '{"raw_text":["fallback hit"],"findings":[],"confidence":0.7}'
                return _Response()

        class _FakeClient:
            def __init__(self, api_key: str):
                self.responses = _FakeResponses()

        temp_dir = Path("outputs") / "test_tmp"
        temp_dir.mkdir(parents=True, exist_ok=True)
        image_path = temp_dir / "tiny_fallback.png"
        try:
            image_path.write_bytes(png_bytes)
            with patch.dict(os.environ, {"OPENAI_API_KEY": "test-key"}):
                with patch("openai.OpenAI", _FakeClient):
                    extractor = OpenAIVisionTextExtractor(
                        model="gpt-5",
                        fallback_models=("gpt-4.1-mini",),
                    )
                    result = extractor.extract(image_path)
        finally:
            if image_path.exists():
                image_path.unlink()

        self.assertEqual(result.raw_text, ["fallback hit"])
        self.assertEqual(result.model_used, "gpt-4.1-mini")
        self.assertEqual(called_models[:2], ["gpt-5", "gpt-5"])
        self.assertIn("gpt-4.1-mini", called_models)


if __name__ == "__main__":
    unittest.main()

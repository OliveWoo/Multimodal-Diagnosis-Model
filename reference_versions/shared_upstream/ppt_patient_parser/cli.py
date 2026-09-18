from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any, Sequence

from .parsers.slide_parser import parse_slide
from .ppt_reader import collect_pptx_files, open_presentation
from .slide_extractor import extract_slide_payload
from .utils.io import ensure_dir, make_presentation_slug, write_json
from .utils.logging import configure_logging
from .vision_hooks import OpenAIVisionTextExtractor


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract slide text, images, and optional image-backed text from PPTX files."
    )
    parser.add_argument("--input", type=Path, required=True, help="Input folder containing .pptx files.")
    parser.add_argument("--output", type=Path, required=True, help="Output folder for extracted artifacts.")
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING, ERROR).",
    )
    parser.add_argument(
        "--vision-model",
        default="",
        help="Optional OpenAI model used to extract text from embedded images.",
    )
    parser.add_argument(
        "--vision-fallback-model",
        action="append",
        default=[],
        help="Optional fallback OpenAI model. Repeat to provide multiple fallbacks.",
    )
    return parser.parse_args(argv)


def _write_slide_text(path: Path, text_blocks: list[str]) -> None:
    content = "\n".join(block for block in text_blocks if block.strip())
    path.write_text(content, encoding="utf-8")


def _serialize_image_record(image) -> dict[str, Any]:
    return {
        "image_index": image.image_index,
        "raw_path": str(image.raw_path),
        "width": image.width,
        "height": image.height,
        "ext": image.ext,
        "shape_name": image.shape_name,
        "shape_order": image.shape_order,
        "display_geometry": dict(image.display_geometry),
        "preprocessing": {
            "suspicious_distortion": False,
            "notes": [],
        },
        "vision": {
            "status": "not_requested",
            "raw_text": [],
            "findings": [],
            "confidence": None,
        },
    }


def _extract_image_text(image_record: dict[str, Any], vision_extractor: OpenAIVisionTextExtractor | None) -> str | None:
    if vision_extractor is None:
        return None

    image_path = Path(image_record["raw_path"])
    try:
        result = vision_extractor.extract(image_path)
    except Exception as exc:  # noqa: BLE001
        image_record["vision"] = {
            "status": "error",
            "raw_text": [],
            "findings": [],
            "confidence": None,
            "error": str(exc),
        }
        return f"image_{image_record['image_index']}_vision_error: {exc}"

    status = "success" if result.raw_text or result.findings else "empty"
    vision_payload = result.to_dict()
    vision_payload["status"] = status
    image_record["vision"] = vision_payload
    return None


def _build_debug_payload(
    *,
    pptx_path: Path,
    slide_index: int,
    extracted,
    image_records: list[dict[str, Any]],
    parsed_slide: dict[str, Any],
    vision_extractor: OpenAIVisionTextExtractor | None,
    errors: list[str],
) -> dict[str, Any]:
    return {
        "source_file": pptx_path.name,
        "source_path": str(pptx_path),
        "slide_index": slide_index,
        "errors": errors,
        "vision": {
            "enabled": vision_extractor is not None,
            "model": getattr(vision_extractor, "model", None),
            "fallback_models": list(getattr(vision_extractor, "fallback_models", ())),
        },
        "text_blocks": extracted.text_blocks,
        "text_records": extracted.text_records,
        "shape_debug": extracted.debug_shapes,
        "images": image_records,
        "parser_summary": {
            "candidate_modalities": parsed_slide["candidate_modalities"],
            "detected_dates": parsed_slide["detected_dates"],
            "image_extraction_status_summary": parsed_slide["image_extraction_status_summary"],
            "needs_review": parsed_slide["needs_review"],
            "review_reasons": parsed_slide["review_reasons"],
            "text_source_summary": parsed_slide["text_source_summary"],
        },
    }


def process_presentation(
    pptx_path: Path,
    output_root: Path,
    logger: logging.Logger,
    vision_extractor: OpenAIVisionTextExtractor | None = None,
) -> None:
    presentation_name = make_presentation_slug(pptx_path.stem)
    presentation_dir = ensure_dir(output_root / presentation_name)
    slides_dir = ensure_dir(presentation_dir / "slides")

    presentation = open_presentation(pptx_path)
    total_slides = len(presentation.slides)
    logger.info("Processing %s (%d slides)", pptx_path.name, total_slides)

    for slide_index, slide in enumerate(presentation.slides, start=1):
        logger.info(
            "[%s] Slide %03d/%03d: extracting shapes",
            pptx_path.name,
            slide_index,
            total_slides,
        )
        slide_dir = ensure_dir(slides_dir / f"slide_{slide_index:03d}")
        images_dir = ensure_dir(slide_dir / "images")
        try:
            extracted = extract_slide_payload(slide, slide_index=slide_index, images_dir=images_dir)
            errors: list[str] = []
            image_records = [_serialize_image_record(image) for image in extracted.images]
            logger.info(
                "[%s] Slide %03d/%03d: found %d text blocks and %d images",
                pptx_path.name,
                slide_index,
                total_slides,
                len(extracted.text_blocks),
                len(image_records),
            )
            for image_position, image_record in enumerate(image_records, start=1):
                if vision_extractor is not None:
                    logger.info(
                        "[%s] Slide %03d/%03d: OCR image %d/%d",
                        pptx_path.name,
                        slide_index,
                        total_slides,
                        image_position,
                        len(image_records),
                    )
                vision_error = _extract_image_text(image_record, vision_extractor)
                if vision_error:
                    errors.append(vision_error)
                    logger.warning(
                        "[%s] Slide %03d/%03d: %s",
                        pptx_path.name,
                        slide_index,
                        total_slides,
                        vision_error,
                    )

            parsed_slide = parse_slide(
                source_file=Path(pptx_path.name),
                slide_index=slide_index,
                text_blocks=extracted.text_blocks,
                text_records=extracted.text_records,
                image_records=image_records,
            )
            _write_slide_text(slide_dir / "slide_text.txt", parsed_slide["analysis_text_blocks"])
            write_json(slide_dir / "parsed.json", parsed_slide)
            write_json(
                slide_dir / "debug.json",
                _build_debug_payload(
                    pptx_path=pptx_path,
                    slide_index=slide_index,
                    extracted=extracted,
                    image_records=image_records,
                    parsed_slide=parsed_slide,
                    vision_extractor=vision_extractor,
                    errors=errors,
                ),
            )
            logger.info(
                "[%s] Slide %03d/%03d: done (review=%s, native_text=%d, image_text=%d)",
                pptx_path.name,
                slide_index,
                total_slides,
                parsed_slide["needs_review"],
                parsed_slide["text_source_summary"]["native_text_block_count"],
                parsed_slide["text_source_summary"]["image_text_block_count"],
            )
        except Exception as exc:  # noqa: BLE001
            logger.exception("Failed to process slide %s in %s", slide_index, pptx_path)
    logger.info("Processed %s slides from %s", len(presentation.slides), pptx_path)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    logger = configure_logging(args.log_level)

    input_dir = args.input.expanduser().resolve()
    output_dir = ensure_dir(args.output.expanduser().resolve())
    if not input_dir.exists():
        raise FileNotFoundError(f"Input folder not found: {input_dir}")

    pptx_files = collect_pptx_files(input_dir)
    if not pptx_files:
        logger.warning("No .pptx files found under %s", input_dir)
        return
    logger.info("Found %d PPTX files under %s", len(pptx_files), input_dir)

    vision_extractor = None
    if args.vision_model.strip():
        vision_extractor = OpenAIVisionTextExtractor(
            model=args.vision_model.strip(),
            fallback_models=tuple(model.strip() for model in args.vision_fallback_model if model.strip()),
        )
        logger.info(
            "Vision OCR enabled with model=%s fallback_models=%s",
            vision_extractor.model,
            list(vision_extractor.fallback_models),
        )
    else:
        logger.info("Vision OCR disabled")

    for pptx_path in pptx_files:
        process_presentation(pptx_path, output_dir, logger, vision_extractor=vision_extractor)

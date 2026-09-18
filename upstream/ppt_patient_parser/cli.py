from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

from .ppt_reader import collect_pptx_files, open_presentation
from .slide_extractor import extract_slide_payload
from .utils.io import ensure_dir, make_presentation_slug
from .utils.logging import configure_logging


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract slide text and images from PPTX files."
    )
    parser.add_argument("--input", type=Path, required=True, help="Input folder containing .pptx files.")
    parser.add_argument("--output", type=Path, required=True, help="Output folder for extracted artifacts.")
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING, ERROR).",
    )
    return parser.parse_args(argv)


def _write_slide_text(path: Path, text_blocks: list[str]) -> None:
    content = "\n".join(block for block in text_blocks if block.strip())
    path.write_text(content, encoding="utf-8")


def process_presentation(
    pptx_path: Path,
    output_root: Path,
    logger: logging.Logger,
) -> None:
    presentation_name = make_presentation_slug(pptx_path.stem)
    presentation_dir = ensure_dir(output_root / presentation_name)
    slides_dir = ensure_dir(presentation_dir / "slides")

    presentation = open_presentation(pptx_path)

    for slide_index, slide in enumerate(presentation.slides, start=1):
        slide_dir = ensure_dir(slides_dir / f"slide_{slide_index:03d}")
        images_dir = ensure_dir(slide_dir / "images")
        try:
            extracted = extract_slide_payload(slide, slide_index=slide_index, images_dir=images_dir)
            _write_slide_text(slide_dir / "slide_text.txt", extracted.text_blocks)
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

    for pptx_path in pptx_files:
        process_presentation(pptx_path, output_dir, logger)

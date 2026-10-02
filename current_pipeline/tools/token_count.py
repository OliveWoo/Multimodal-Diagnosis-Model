from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Iterable

import tiktoken

from core.llm_parser import RAW_FILE_SUFFIX, select_prompt

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format=LOG_FORMAT)


def load_json(filepath: Path) -> dict:
    logging.debug("Loading JSON payload from %s", filepath)
    with filepath.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_prompt(template: str, payload: dict) -> str:
    serialized_payload = json.dumps(payload, ensure_ascii=False, indent=2)
    return f"{template.rstrip()}\n\n### JSON payload ###\n{serialized_payload}"


def resolve_encoding(model: str) -> tiktoken.Encoding:
    try:
        return tiktoken.encoding_for_model(model)
    except KeyError:
        logging.warning(
            "Unknown model '%s'. Falling back to cl100k_base encoding.", model
        )
        return tiktoken.get_encoding("cl100k_base")


def iter_target_files(
    target: Path, suffix: str, recursive: bool = False
) -> Iterable[Path]:
    if target.is_file():
        return [target]

    if not target.is_dir():
        raise FileNotFoundError(f"Path does not exist or is not accessible: {target}")

    pattern = suffix if suffix.startswith("*") else f"*{suffix}"
    iterator = target.rglob(pattern) if recursive else target.glob(pattern)
    return sorted(path for path in iterator if path.is_file())


def count_tokens_for_file(
    file_path: Path, encoding: tiktoken.Encoding
) -> int:
    payload = load_json(file_path)
    template = select_prompt(file_path, payload)
    prompt = build_prompt(template, payload)
    token_count = len(encoding.encode(prompt))
    logging.debug(
        "Computed %d tokens for %s (payload length=%d characters).",
        token_count,
        file_path,
        len(prompt),
    )
    return token_count


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute token usage for prompts generated from JSON payloads."
    )
    parser.add_argument(
        "path",
        type=Path,
        help="Directory or JSON file to analyze.",
    )
    parser.add_argument(
        "--model",
        default="gpt-5",
        help="Model name used to look up the tiktoken encoding (default: gpt-5).",
    )
    parser.add_argument(
        "--suffix",
        default=RAW_FILE_SUFFIX,
        help=(
            "Filename suffix to match (default: _Raw.json). "
            "Prefix with '*' to provide a full glob pattern."
        ),
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Scan directories recursively for matching files.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )

    args = parser.parse_args()
    configure_logging(args.verbose)

    target = args.path.resolve()
    try:
        files = list(iter_target_files(target, args.suffix, args.recursive))
    except FileNotFoundError as exc:
        logging.error(str(exc))
        return

    if not files:
        logging.warning("No JSON files matched in %s", target)
        return

    encoding = resolve_encoding(args.model)
    base_dir = target if target.is_dir() else target.parent
    total_tokens = 0
    processed_files = 0

    for file_path in files:
        try:
            token_count = count_tokens_for_file(file_path, encoding)
        except FileNotFoundError as exc:
            logging.error("Missing prompt for %s: %s", file_path, exc)
            continue
        except json.JSONDecodeError as exc:
            logging.error("Invalid JSON in %s: %s", file_path, exc)
            continue
        except Exception as exc:  # pylint: disable=broad-exception-caught
            logging.exception("Failed to process %s: %s", file_path, exc)
            continue

        try:
            display_path = file_path.relative_to(base_dir)
        except ValueError:
            display_path = file_path

        print(f"{display_path}: {token_count} tokens")
        total_tokens += token_count
        processed_files += 1

    if processed_files:
        print(f"Processed {processed_files} file(s); total tokens: {total_tokens}")


if __name__ == "__main__":
    main()

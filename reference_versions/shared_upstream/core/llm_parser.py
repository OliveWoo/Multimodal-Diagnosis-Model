"""LLM-powered formatter that converts Raw JSON payloads into structured outputs."""

from __future__ import annotations

import argparse
import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional

from utils import sanitize_filename

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
FAILED_DIR = Path("outputs") / "logs" / "failed"
DEFAULT_INPUT_DIR = Path("raw_excel")
RAW_FILE_SUFFIX = "_Raw.json"
PROMPTS_DIR = Path("prompts")


def configure_logging(verbose: bool = False) -> None:
    """Configure baseline logging and honor the --verbose flag."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format=LOG_FORMAT)


def load_json(filepath: Path) -> dict:
    """Load a Raw JSON payload from disk."""
    logging.debug("Loading JSON from %s", filepath)
    with filepath.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def select_prompt(filepath: Path, payload: Optional[dict] = None) -> str:
    """Derive the prompt template name that best matches the Raw JSON file."""

    def _candidate_names() -> Iterable[str]:
        if payload:
            sheet_name = payload.get("sheet")
            if isinstance(sheet_name, str) and sheet_name:
                yield f"{sheet_name}_prompt.txt"
                normalized = sheet_name.replace("_", " ")
                if normalized != sheet_name:
                    yield f"{normalized}_prompt.txt"
                compact = sheet_name.replace(" ", "_")
                if compact != sheet_name:
                    yield f"{compact}_prompt.txt"

        stem = filepath.stem
        if stem.endswith("_Raw"):
            stem = stem[: -len("_Raw")]
        yield f"{stem}_prompt.txt"
        if "_" in stem:
            yield f"{stem.replace('_', ' ')}_prompt.txt"

    seen: set[str] = set()
    for name in _candidate_names():
        if name in seen:
            continue
        seen.add(name)
        prompt_path = PROMPTS_DIR / name
        if prompt_path.exists():
            logging.debug("Using prompt template %s for %s", prompt_path, filepath.name)
            return prompt_path.read_text(encoding="utf-8")

    tried = ", ".join(seen) if seen else "(none)"
    raise FileNotFoundError(
        f"No prompt template matched for {filepath.name}; tried: {tried}"
    )


def _build_prompt(template: str, payload: dict) -> str:
    """Combine the prompt template and payload into the final request body."""
    serialized_payload = json.dumps(payload, ensure_ascii=False, indent=2)
    return f"{template.rstrip()}\n\n### JSON payload ###\n{serialized_payload}"


def _strip_json_markers(text: str) -> str:
    """Trim surrounding ```json fences when the model responds with Markdown."""
    stripped = text.strip()
    if not stripped.startswith("```"):
        return stripped
    lines = stripped.splitlines()
    if len(lines) >= 3 and lines[-1].strip() == "```":
        body = "\n".join(lines[1:-1]).strip()
        return body if body else stripped
    return stripped


def send_to_llm(prompt: str) -> str:
    """Send the prompt to OpenAI and return the raw string response."""
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY environment variable is not set.")

    logging.debug("Dispatching prompt to OpenAI (length=%d)", len(prompt))
    from openai import OpenAI  # type: ignore

    client = OpenAI(api_key=api_key)
    response = client.chat.completions.create(
        model="gpt-5",
        temperature=1,
        messages=[{"role": "user", "content": prompt}],
    )
    return response.choices[0].message.content


def validate_json(output_text: str) -> Optional[dict]:
    """Ensure the LLM response is valid JSON before persisting it."""
    candidate = _strip_json_markers(output_text)
    if candidate != output_text.strip():
        logging.debug("Stripped Markdown JSON markers before parsing.")
    try:
        parsed = json.loads(candidate)
    except json.JSONDecodeError:
        logging.warning("Failed to parse JSON output.")
        logging.debug("Raw output that failed to parse: %s", output_text)
        return None
    return parsed


def _formatted_filename_from_raw(source_path: Path) -> str:
    """Map a Raw JSON filename to its formatted JSON counterpart."""
    stem = sanitize_filename(source_path.stem)
    if stem.endswith("_Raw"):
        stem = stem[: -len("_Raw")]
    return f"{stem}.json"


def save_result(parsed: dict, source_path: Path, output_dir: Path) -> Path:
    """Persist a successful LLM response into the formatted output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / _formatted_filename_from_raw(source_path)
    destination.write_text(json.dumps(parsed, ensure_ascii=False, indent=2), encoding="utf-8")
    logging.info("Saved structured output to %s", destination)
    return destination


def save_failed_output(raw_output: str, source_path: Path) -> Path:
    """Store the raw LLM text when JSON parsing fails."""
    FAILED_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base_name = sanitize_filename(source_path.stem)
    destination = FAILED_DIR / f"{base_name}_{timestamp}.txt"
    destination.write_text(raw_output, encoding="utf-8")
    logging.info("Stored failed output to %s", destination)
    return destination


def process_file(filepath: Path, output_dir: Path) -> None:
    """Process one Raw JSON file end-to-end."""
    logging.info("Processing %s", filepath)

    try:
        payload = load_json(filepath)
        template = select_prompt(filepath, payload)
        prompt = _build_prompt(template, payload)
        raw_output = send_to_llm(prompt)
        parsed = validate_json(raw_output)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.exception("Unhandled error while processing %s: %s", filepath, exc)
        return

    if parsed is not None:
        output_path = save_result(parsed, filepath, output_dir)
        print(f"Converted successfully: {output_path}")
    else:
        failed_path = save_failed_output(raw_output, filepath)
        print(f"Conversion failed. Raw output saved to: {failed_path}")


def _iter_raw_files(input_dir: Path) -> Iterable[Path]:
    """Yield Raw JSON files located in the provided directory."""
    return sorted(path for path in input_dir.glob(f"*{RAW_FILE_SUFFIX}") if path.is_file())


def parse_with_llm(input_dir: Path, output_dir: Path, *, verbose: bool | None = None) -> None:
    """Library entry point used by other modules."""
    if verbose is not None:
        configure_logging(verbose)

    if not input_dir.exists():
        logging.warning("LLM input directory does not exist: %s", input_dir)
        return

    files = list(_iter_raw_files(input_dir))
    if not files:
        logging.warning("No Raw JSON files were found in %s", input_dir)
        return

    for file_path in files:
        process_file(file_path.resolve(), output_dir.resolve())


def gather_files(target: Optional[Path], input_dir: Path) -> Iterable[Path]:
    """Resolve CLI arguments into a list of Raw JSON files to process."""
    if target is not None:
        if not target.is_file():
            raise FileNotFoundError(f"Specified file does not exist: {target}")
        return [target]

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    return list(_iter_raw_files(input_dir))


def main() -> None:
    """CLI entry point for manual invocation."""
    parser = argparse.ArgumentParser(
        description="Send Raw JSON files to OpenAI and persist structured responses.",
    )
    parser.add_argument(
        "target",
        nargs="?",
        type=Path,
        help="Optional single Raw JSON file to process.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory to scan when no file is provided.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory for formatted JSON files. Defaults to the inferred `_json` folder.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )

    args = parser.parse_args()
    configure_logging(args.verbose)

    try:
        files = list(gather_files(args.target, args.input_dir))
    except FileNotFoundError as exc:
        logging.error(str(exc))
        return

    if not files:
        logging.warning("No JSON files found to process.")
        return

    if args.output_dir:
        output_dir = args.output_dir
    else:
        base_dir = args.target.parent if args.target else args.input_dir
        if base_dir.name.endswith("_json_Raw"):
            output_dir = base_dir.with_name(base_dir.name[: -len("_json_Raw")] + "_json")
        else:
            output_dir = base_dir / "formatted_json"

    output_dir.mkdir(parents=True, exist_ok=True)

    for file_path in files:
        process_file(file_path.resolve(), output_dir.resolve())


if __name__ == "__main__":
    main()

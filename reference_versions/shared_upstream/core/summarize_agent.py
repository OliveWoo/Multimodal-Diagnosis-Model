from __future__ import annotations

import argparse
import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

from analytics.tracking import record_summary_result
from utils import confirm_overwrite, sanitize_filename

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
PROMPTS_DIR = Path("agents") / "prompts"
BIG_PROMPT_FILE = PROMPTS_DIR / "big_prompt.txt"
FAILED_DIR = Path("outputs") / "logs" / "failed"
SUMMARY_DIR_NAME = "summary_outputs"
AGENT_OUTPUT_DIR_NAME = "agent_outputs"
DEFAULT_MODEL = "gpt-5"


class AgentOutputMissingError(RuntimeError):
    """Raised when a required agent output file is missing."""


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format=LOG_FORMAT)


def extract_patient_identifier(input_dir: Path) -> str:
    base_name = sanitize_filename(input_dir.name)
    if base_name.endswith("_json"):
        base_name = base_name[: -len("_json")] or base_name
    return base_name


def normalize_patient_directory(path: Path) -> Path | None:
    if path.name in {AGENT_OUTPUT_DIR_NAME, SUMMARY_DIR_NAME}:
        logging.info("Skipping output directory input: %s", path)
        return None
    return path.resolve()


def looks_like_patient_dir(path: Path) -> bool:
    return (path / AGENT_OUTPUT_DIR_NAME).is_dir()


def collect_patient_directories(entries: Sequence[Path]) -> list[Path]:
    collected: list[Path] = []
    seen: set[Path] = set()

    for entry in entries:
        path = entry.expanduser()
        if not path.exists():
            logging.warning("Skipping non-existent path: %s", path)
            continue
        normalized = normalize_patient_directory(path)
        if normalized is None:
            continue
        if not normalized.is_dir():
            logging.warning("Skipping non-directory path: %s", normalized)
            continue

        candidates = [normalized]
        try:
            candidates.extend(
                child.resolve()
                for child in normalized.iterdir()
                if child.is_dir() and child.name not in {AGENT_OUTPUT_DIR_NAME, SUMMARY_DIR_NAME}
            )
        except PermissionError:
            logging.warning("Unable to enumerate subdirectories for %s", normalized)

        for candidate in candidates:
            if candidate in seen:
                continue
            if looks_like_patient_dir(candidate):
                collected.append(candidate)
                seen.add(candidate)
            else:
                logging.debug("Ignoring directory without agent outputs: %s", candidate)

    collected.sort()
    return collected


def load_json(filepath: Path) -> dict:
    with filepath.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_big_prompt() -> str:
    if not BIG_PROMPT_FILE.exists():
        raise FileNotFoundError(f"Big prompt template not found: {BIG_PROMPT_FILE}")
    return BIG_PROMPT_FILE.read_text(encoding="utf-8")


def _expected_agent_paths(
    input_dir: Path,
    *,
    include_underlying: bool,
    include_filmarray_gmtest: bool,
    filmarray_agent_suffix: str,
) -> Sequence[tuple[str, Path, bool]]:
    base_name = extract_patient_identifier(input_dir)
    agent_dir = input_dir / AGENT_OUTPUT_DIR_NAME
    if not agent_dir.exists():
        raise AgentOutputMissingError(f"Agent output directory not found: {agent_dir}")
    specs: list[tuple[str, Path, bool]] = []
    if include_underlying:
        # Keep backward compatibility: include when available, but do not fail old datasets.
        specs.append(("UnderlyingAgent", agent_dir / f"{base_name}_underlying_agent.json", False))
    specs.extend(
        [
            ("CBCOtherLabAgent", agent_dir / f"{base_name}_cbc_other_lab_agent.json", True),
            ("ImageAgent", agent_dir / f"{base_name}_image_agent.json", True),
            ("CultureAgent", agent_dir / f"{base_name}_culture_agent.json", True),
        ]
    )
    if include_filmarray_gmtest:
        specs.insert(2, ("FilmArrayGMTestAgent", agent_dir / f"{base_name}_{filmarray_agent_suffix}.json", True))
    return specs


def collect_agent_outputs(
    input_dir: Path,
    *,
    include_underlying: bool,
    include_filmarray_gmtest: bool,
    filmarray_agent_suffix: str,
) -> list[tuple[str, Path, dict]]:
    sections: list[tuple[str, Path, dict]] = []
    for label, path, required in _expected_agent_paths(
        input_dir,
        include_underlying=include_underlying,
        include_filmarray_gmtest=include_filmarray_gmtest,
        filmarray_agent_suffix=filmarray_agent_suffix,
    ):
        if not path.exists():
            if not required:
                logging.warning("Optional agent output not found: %s", path)
                continue
            raise AgentOutputMissingError(f"Required agent output not found: {path}")
        data = load_json(path)
        sections.append((label, path, data))
    return sections


def build_prompt(template: str, sections: Sequence[tuple[str, Path, dict]]) -> str:
    blocks: list[str] = [template.rstrip(), "", "### Agent outputs ###"]
    for label, path, payload in sections:
        serialized = json.dumps(payload, ensure_ascii=False, indent=2)
        blocks.extend([f"#### {label} ({path.name})", serialized])
    return "\n".join(blocks)


def _strip_json_markers(text: str) -> str:
    stripped = text.strip()
    if not stripped.startswith("```"):
        return stripped
    lines = stripped.splitlines()
    if len(lines) >= 3 and lines[-1].strip() == "```":
        return "\n".join(lines[1:-1]).strip()
    return stripped


def validate_json(output_text: str) -> dict | None:
    candidate = _strip_json_markers(output_text)
    try:
        return json.loads(candidate)
    except json.JSONDecodeError:
        logging.error("Failed to parse LLM output as JSON.")
        logging.debug("Raw output that failed to parse:\n%s", output_text)
        return None


def _extract_output_text(response: Any) -> str:
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
    combined = "".join(segments).strip()
    return combined


def send_to_llm(prompt: str, *, model: str) -> str:
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY environment variable is not set.")
    from openai import OpenAI  # type: ignore

    client = OpenAI(api_key=api_key)
    response = client.responses.create(
        model=model,
        temperature=1,
        input=prompt,
    )
    message = _extract_output_text(response)
    if not message:
        raise RuntimeError("OpenAI response did not contain any content.")
    return message


def _summary_mode_suffix(*, include_underlying: bool, include_filmarray_gmtest: bool) -> str:
    if include_underlying and include_filmarray_gmtest:
        return "with_underlying_with_filmarray"
    if include_underlying and not include_filmarray_gmtest:
        return "with_underlying_no_filmarray"
    if not include_underlying and include_filmarray_gmtest:
        return "with_filmarray"
    return "no_filmarray"


def derive_output_path(
    input_dir: Path,
    *,
    include_underlying: bool,
    include_filmarray_gmtest: bool,
    no_filmarray_naming: bool = False,
) -> Path:
    base_name = extract_patient_identifier(input_dir)
    destination_dir = input_dir / SUMMARY_DIR_NAME
    destination_dir.mkdir(parents=True, exist_ok=True)
    mode_suffix = _summary_mode_suffix(
        include_underlying=include_underlying,
        include_filmarray_gmtest=(False if no_filmarray_naming else include_filmarray_gmtest),
    )
    return destination_dir / f"{base_name}_final_summary_{mode_suffix}.json"


def save_result(payload: dict, input_dir: Path, output_path: Path | None = None) -> Path:
    if output_path is None:
        output_path = derive_output_path(
            input_dir, include_underlying=True, include_filmarray_gmtest=True
        )
    output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return output_path


def save_failed_output(raw_output: str, input_dir: Path) -> Path:
    FAILED_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    identifier = sanitize_filename(input_dir.name)
    destination = FAILED_DIR / f"{identifier}_summary_{timestamp}.txt"
    destination.write_text(raw_output, encoding="utf-8")
    return destination


def summarize_directory(
    input_dir: Path,
    *,
    model: str,
    include_underlying: bool = False,
    include_filmarray_gmtest: bool = True,
    no_filmarray_naming: bool = False,
    filmarray_agent_suffix: str = "filmarray_gmtest_agent",
) -> None:
    try:
        sections = collect_agent_outputs(
            input_dir,
            include_underlying=include_underlying,
            include_filmarray_gmtest=include_filmarray_gmtest,
            filmarray_agent_suffix=filmarray_agent_suffix,
        )
    except AgentOutputMissingError as exc:
        logging.error(str(exc))
        return

    patient_id = extract_patient_identifier(input_dir)
    output_path = derive_output_path(
        input_dir,
        include_underlying=include_underlying,
        include_filmarray_gmtest=include_filmarray_gmtest,
        no_filmarray_naming=no_filmarray_naming,
    )
    if output_path.exists() and not confirm_overwrite(
        output_path, prompt=f"檔案 {output_path} 已存在，是否覆蓋？"
    ):
        logging.info("略過 %s，因為輸出檔案已存在：%s", patient_id, output_path)
        return

    try:
        template = load_big_prompt()
        prompt = build_prompt(template, sections)
        raw_output = send_to_llm(prompt, model=model)
        parsed = validate_json(raw_output)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.exception("Failed while summarizing %s: %s", input_dir, exc)
        return

    if parsed is None:
        failed_path = save_failed_output(raw_output, input_dir)
        logging.error("Stored unparsed summary output at %s", failed_path)
        return

    output_path = save_result(parsed, input_dir, output_path=output_path)
    logging.info("Wrote final summary to %s", output_path)
    try:
        dataset_path = record_summary_result(patient_id, parsed)
        logging.debug("Recorded summary analytics at %s", dataset_path)
    except Exception as exc:  # pragma: no cover - best-effort logging
        logging.exception(
            "Failed to record summary analytics for %s: %s", patient_id, exc
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine agent outputs with the big prompt and request a final summary from the GPT API.",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help=(
            "One or more patient directories containing agent outputs, "
            "or parent directories that include multiple patients."
        ),
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Model identifier to use when calling the GPT API (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    parser.add_argument(
        "--include-underlying",
        dest="include_underlying",
        action="store_true",
        help="Include underlying agent output in summary prompt (default: excluded).",
    )
    parser.add_argument(
        "--exclude-underlying",
        dest="include_underlying",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--exclude-filmarray",
        action="store_true",
        help="Exclude filmarray/GMTest agent output from summary prompt (default: included).",
    )
    parser.add_argument(
        "--no-filmarray-naming",
        action="store_true",
        help="Use no_filmarray naming suffix even when FilmArrayGMTestAgent is included.",
    )
    parser.add_argument(
        "--filmarray-agent-suffix",
        default="filmarray_gmtest_agent",
        help="Suffix of the filmarray/GMTest agent file to read (default: filmarray_gmtest_agent).",
    )

    args = parser.parse_args()
    configure_logging(args.verbose)

    patient_dirs = collect_patient_directories(args.inputs)
    if not patient_dirs:
        logging.error("No valid patient directories with agent outputs found.")
        return

    logging.info("Found %d patient directories to summarize.", len(patient_dirs))
    for patient_dir in patient_dirs:
        logging.info("Summarizing %s", patient_dir)
        if args.exclude_filmarray:
            summarize_directory(
                patient_dir,
                model=args.model,
                include_underlying=args.include_underlying,
                include_filmarray_gmtest=False,
                no_filmarray_naming=args.no_filmarray_naming,
                filmarray_agent_suffix=args.filmarray_agent_suffix,
            )
            continue

        # Default behavior: emit both with_filmarray and no_filmarray variants.
        summarize_directory(
            patient_dir,
            model=args.model,
            include_underlying=args.include_underlying,
            include_filmarray_gmtest=True,
            no_filmarray_naming=False,
            filmarray_agent_suffix=args.filmarray_agent_suffix,
        )
        summarize_directory(
            patient_dir,
            model=args.model,
            include_underlying=args.include_underlying,
            include_filmarray_gmtest=True,
            no_filmarray_naming=True,
            filmarray_agent_suffix="filmarray_gmtest_agent_no_filmarray",
        )


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from utils import confirm_overwrite, sanitize_filename

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
FULL_PROMPT_PATH = Path("agents") / "prompts" / "mNGS_max_prompt.txt"
NO_FILMARRAY_PROMPT_PATH = Path("agents") / "prompts" / "mNGS_max_prompt_no_filmarray.txt"
PROMPT_PATH = FULL_PROMPT_PATH
SUMMARY_DIR_NAME = "summary_outputs"
SUMMARY_OUTPUT_DIR_NAME = "summary_outputs"
DEFAULT_MODEL = "gpt-5"
OUTPUT_SUFFIX = "mNGS_max_agent_retry"
FAILED_DIR = Path("outputs") / "logs" / "failed"
SUMMARY_MODE_TO_SUFFIX = {
    "full": "with_filmarray",
    "no_filmarray": "no_filmarray",
}
SUMMARY_MODE_COMPAT_SUFFIXES = {
    "full": ["with_filmarray", "with_underlying_with_filmarray"],
    "no_filmarray": ["no_filmarray", "with_underlying_no_filmarray"],
}
SUMMARY_MODE_TO_PROMPT_PATH = {
    "full": FULL_PROMPT_PATH,
    "no_filmarray": NO_FILMARRAY_PROMPT_PATH,
}

try:
    import tiktoken
except ImportError:  # pragma: no cover - optional dependency
    tiktoken = None  # type: ignore[assignment]


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format=LOG_FORMAT)


def extract_patient_identifier(input_dir: Path) -> str:
    base_name = sanitize_filename(input_dir.name)
    if base_name.endswith("_json"):
        base_name = base_name[: -len("_json")] or base_name
    return base_name


def load_prompt(path: Path | None = None) -> str:
    prompt_path = path or PROMPT_PATH
    if not prompt_path.exists():
        raise FileNotFoundError(f"Prompt template not found: {prompt_path}")
    return prompt_path.read_text(encoding="utf-8")


def load_json(filepath: Path) -> dict:
    with filepath.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def find_final_summary_file(input_dir: Path, *, summary_mode: str = "auto") -> Path | None:
    if summary_mode not in {"auto", *SUMMARY_MODE_TO_SUFFIX.keys()}:
        raise ValueError(f"Unknown summary mode: {summary_mode}")

    base_name = extract_patient_identifier(input_dir)
    summary_dir = input_dir / SUMMARY_DIR_NAME
    search_dirs: list[Path] = [summary_dir, input_dir]

    if summary_mode == "auto":
        ordered_modes = ["full", "no_filmarray"]
    else:
        ordered_modes = [summary_mode]

    for mode in ordered_modes:
        for suffix in SUMMARY_MODE_COMPAT_SUFFIXES[mode]:
            filename = f"{base_name}_final_summary_{suffix}.json"
            for directory in search_dirs:
                path = directory / filename
                if path.exists():
                    return path

    legacy_filename = f"{base_name}_final_summary.json"
    for directory in search_dirs:
        path = directory / legacy_filename
        if path.exists():
            return path

    return None


def find_mngs_name_reads_file(input_dir: Path) -> Path | None:
    base_name = extract_patient_identifier(input_dir)
    candidates = [
        input_dir / f"{base_name}_mNGS_grouped.json",
        input_dir / f"{base_name}_mngs_name_reads.json",
        input_dir / f"{base_name}_mNGS_candidates_for_llm.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def find_underlying_agent_file(input_dir: Path) -> Path | None:
    base_name = extract_patient_identifier(input_dir)
    agent_dir = input_dir / "agent_outputs"
    preferred = agent_dir / f"{base_name}_underlying_agent.json"
    if preferred.exists():
        return preferred
    fallback = sorted(agent_dir.glob("*_underlying_agent.json")) if agent_dir.is_dir() else []
    return fallback[0] if fallback else None


def detect_summary_variant(base_name: str, summary_path: Path) -> str:
    name = summary_path.name
    for mode, suffixes in SUMMARY_MODE_COMPAT_SUFFIXES.items():
        for suffix in suffixes:
            expected = f"{base_name}_final_summary_{suffix}.json"
            if name == expected:
                return mode
    if name == f"{base_name}_final_summary.json":
        return "legacy"
    return "custom"


def looks_like_patient_dir(path: Path, *, summary_mode: str, require_underlying: bool) -> bool:
    has_underlying = find_underlying_agent_file(path) is not None
    return (
        find_final_summary_file(path, summary_mode=summary_mode) is not None
        and find_mngs_name_reads_file(path) is not None
        and (has_underlying or not require_underlying)
    )


def collect_patient_directories(
    entries: Sequence[Path], *, summary_mode: str, require_underlying: bool
) -> list[Path]:
    collected: list[Path] = []
    seen: set[Path] = set()

    for entry in entries:
        path = entry.expanduser()
        if not path.exists():
            logging.warning("Skipping non-existent path: %s", path)
            continue

        resolved = path.resolve()
        if resolved.is_file():
            logging.warning("Skipping file input (expecting directories): %s", resolved)
            continue
        if not resolved.is_dir():
            logging.warning("Skipping non-directory path: %s", resolved)
            continue

        candidates = [resolved]
        try:
            candidates.extend(child.resolve() for child in resolved.iterdir() if child.is_dir())
        except PermissionError:
            logging.warning("Unable to enumerate subdirectories for %s", resolved)

        for candidate in candidates:
            if candidate in seen:
                continue
            if looks_like_patient_dir(
                candidate, summary_mode=summary_mode, require_underlying=require_underlying
            ):
                collected.append(candidate)
                seen.add(candidate)
            else:
                logging.debug("Ignoring directory without required inputs: %s", candidate)

    collected.sort()
    return collected


def build_prompt(
    template: str,
    *,
    final_summary: dict,
    mngs_grouped: dict,
    final_path: Path,
    mngs_path: Path,
    underlying_agent: dict | None = None,
    underlying_path: Path | None = None,
) -> str:
    blocks: list[str] = [
        template.rstrip(),
        "",
        "### JSON payloads ###",
        f"#### final_summary ({final_path.name})",
        json.dumps(final_summary, ensure_ascii=False, indent=2),
        "",
    ]
    if underlying_agent is not None and underlying_path is not None:
        blocks.extend(
            [
                f"#### underlying_agent ({underlying_path.name})",
                json.dumps(underlying_agent, ensure_ascii=False, indent=2),
                "",
            ]
        )
    blocks.extend(
        [
        f"#### mngs_name_reads ({mngs_path.name})",
        json.dumps(mngs_grouped, ensure_ascii=False, indent=2),
        ]
    )
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


def estimate_token_count(text: str, *, model: str) -> int | None:
    if tiktoken is None:
        return None
    try:
        encoding = tiktoken.encoding_for_model(model)
    except KeyError:
        encoding = tiktoken.get_encoding("cl100k_base")
    return len(encoding.encode(text))


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
    return "".join(segments).strip()


def send_to_llm(prompt: str, *, model: str) -> tuple[str, dict[str, int | None] | None]:
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
    usage_summary: dict[str, int | None] | None = None
    usage = getattr(response, "usage", None)
    if usage is not None:
        usage_summary = {
            "prompt_tokens": getattr(usage, "input_tokens", None),
            "completion_tokens": getattr(usage, "output_tokens", None),
            "total_tokens": getattr(usage, "total_tokens", None),
        }
    return message, usage_summary


def derive_output_path(input_dir: Path, *, summary_variant: str) -> Path:
    base_name = extract_patient_identifier(input_dir)
    filename = f"{base_name}_{OUTPUT_SUFFIX}_{summary_variant}.json"
    destination_dir = input_dir / SUMMARY_OUTPUT_DIR_NAME
    destination_dir.mkdir(parents=True, exist_ok=True)
    return destination_dir / filename


def save_result(payload: dict, destination: Path) -> Path:
    destination.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return destination


def save_failed_output(raw_output: str, input_dir: Path) -> Path:
    FAILED_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    identifier = sanitize_filename(f"{input_dir.name}_{OUTPUT_SUFFIX}")
    destination = FAILED_DIR / f"{identifier}_{timestamp}.txt"
    destination.write_text(raw_output, encoding="utf-8")
    return destination


def process_directory(
    input_dir: Path,
    *,
    model: str,
    prompt_path: Path | None = None,
    summary_mode: str = "auto",
    require_underlying: bool = True,
) -> None:
    final_summary_path = find_final_summary_file(input_dir, summary_mode=summary_mode)
    if final_summary_path is None:
        logging.warning("Skipping %s; final_summary JSON not found.", input_dir)
        return

    mngs_path = find_mngs_name_reads_file(input_dir)
    if mngs_path is None:
        logging.warning("Skipping %s; mngs_name_reads JSON not found.", input_dir)
        return
    underlying_path = find_underlying_agent_file(input_dir)
    if require_underlying and underlying_path is None:
        logging.warning("Skipping %s; underlying_agent JSON not found.", input_dir)
        return

    base_name = extract_patient_identifier(input_dir)
    summary_variant = detect_summary_variant(base_name, final_summary_path)
    output_path = derive_output_path(input_dir, summary_variant=summary_variant)
    if output_path.exists() and not confirm_overwrite(
        output_path, prompt=f"File {output_path} exists. Overwrite?"
    ):
        logging.info("Skipped %s; output already exists: %s", input_dir.name, output_path)
        return

    try:
        selected_prompt = prompt_path or SUMMARY_MODE_TO_PROMPT_PATH.get(summary_mode, PROMPT_PATH)
        template = load_prompt(selected_prompt)
        final_summary = load_json(final_summary_path)
        underlying_agent = load_json(underlying_path) if underlying_path is not None else None
        mngs_grouped = load_json(mngs_path)
        prompt = build_prompt(
            template,
            final_summary=final_summary,
            mngs_grouped=mngs_grouped,
            final_path=final_summary_path,
            mngs_path=mngs_path,
            underlying_agent=underlying_agent,
            underlying_path=underlying_path,
        )
        estimated_tokens = estimate_token_count(prompt, model=model)
        if estimated_tokens is not None:
            logging.info("Estimated input tokens for %s: %d", input_dir.name, estimated_tokens)
        else:
            logging.debug("Token estimation unavailable (tiktoken missing) for %s", input_dir.name)

        raw_output, usage = send_to_llm(prompt, model=model)
        if usage is not None:
            logging.info(
                "LLM token usage - prompt: %s, completion: %s, total: %s",
                usage.get("prompt_tokens"),
                usage.get("completion_tokens"),
                usage.get("total_tokens"),
            )
        else:
            logging.debug("Token usage data not returned for %s", input_dir.name)

        parsed = validate_json(raw_output)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.exception("Failed while processing %s: %s", input_dir, exc)
        return

    if parsed is None:
        failed_path = save_failed_output(raw_output, input_dir)
        logging.error("Stored unparsed output for %s at %s", input_dir.name, failed_path)
        return

    output_path = save_result(parsed, output_path)
    logging.info("Wrote %s", output_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine final_summary and mngs_name_reads JSON with the big prompt and send to GPT."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Directories containing patient data (patient folders or their parent directories).",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Model identifier to use when calling the GPT API (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--prompt",
        type=Path,
        help=(
            "Override prompt path. "
            f"Default depends on --summary-mode: full={FULL_PROMPT_PATH}, "
            f"no_filmarray={NO_FILMARRAY_PROMPT_PATH}."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    parser.add_argument(
        "--summary-mode",
        choices=["auto", "full", "no_filmarray"],
        default="auto",
        help=(
            "Choose which final_summary variant to load. "
            "auto runs both full and no_filmarray."
        ),
    )
    parser.set_defaults(require_underlying=False)
    parser.add_argument(
        "--require-underlying",
        dest="require_underlying",
        action="store_true",
        help="Require underlying_agent input before running (default: allow missing underlying).",
    )
    parser.add_argument(
        "--allow-missing-underlying",
        dest="require_underlying",
        action="store_false",
        help="Allow running without underlying_agent input (default).",
    )

    args = parser.parse_args()
    configure_logging(args.verbose)

    target_dirs = collect_patient_directories(
        args.inputs,
        summary_mode=args.summary_mode,
        require_underlying=args.require_underlying,
    )
    if not target_dirs:
        logging.error("No patient directories found in the provided inputs.")
        return

    logging.info("Found %d patient directories to process.", len(target_dirs))
    modes_to_run = ["full", "no_filmarray"] if args.summary_mode == "auto" else [args.summary_mode]
    for directory in target_dirs:
        for mode in modes_to_run:
            logging.info("Processing %s (summary-mode=%s)", directory, mode)
            process_directory(
                directory,
                model=args.model,
                prompt_path=args.prompt,
                summary_mode=mode,
                require_underlying=args.require_underlying,
            )


if __name__ == "__main__":
    main()


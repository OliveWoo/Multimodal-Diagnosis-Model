from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analytics.tracking import record_category_result
from core.summarize_agent import DEFAULT_MODEL as SUMMARY_DEFAULT_MODEL, summarize_directory
from utils import confirm_overwrite, sanitize_filename

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
PROMPTS_DIR = Path("agents") / "prompts"
FAILED_DIR = Path("outputs") / "logs" / "failed"
AGENT_OUTPUT_DIR_NAME = "agent_outputs"
SUMMARY_DIR_NAME = "summary_outputs"
DEFAULT_MODEL = "gpt-5"

try:
    import tiktoken
except ImportError:  # pragma: no cover - optional dependency
    tiktoken = None  # type: ignore[assignment]


@dataclass(frozen=True)
class InputSpec:
    label: str
    patterns: Sequence[str]
    allow_missing_as_null: bool = False


@dataclass(frozen=True)
class CategorySpec:
    name: str
    prompt_filename: str
    inputs: Sequence[InputSpec]
    output_suffix: str


@dataclass
class PayloadSection:
    label: str
    path: Path
    data: Any


CATEGORY_SPECS: list[CategorySpec] = [
    CategorySpec(
        name="underlying",
        prompt_filename="underlying_agent_prompt.txt",
        inputs=[InputSpec(label="Underlying", patterns=["*underlying.json"])],
        output_suffix="underlying_agent",
    ),
    CategorySpec(
        name="cbc_other_lab",
        prompt_filename="CBC_otherLab_prompt.txt",
        inputs=[
            InputSpec(label="CBC", patterns=["*CBC.json"]),
            InputSpec(label="OtherLab", patterns=["*other_lab.json", "*other lab.json"]),
        ],
        output_suffix="cbc_other_lab_agent",
    ),
    CategorySpec(
        name="filmarray_gmtest",
        prompt_filename="filmarray_GMtest_prompt.txt",
        inputs=[
            InputSpec(
                label="FilmArray",
                patterns=["*filmarray.json", "*FilmArray.json"],
                allow_missing_as_null=True,
            ),
            InputSpec(
                label="GMTest",
                patterns=["*gm_test.json", "*gm test.json", "*gmTest.json", "*GMtest.json"],
            ),
            InputSpec(
                label="MolecularMicrobiology",
                patterns=["*molecular_microbiology.json", "*molecular microbiology.json"],
                allow_missing_as_null=True,
            ),
        ],
        output_suffix="filmarray_gmtest_agent",
    ),
    CategorySpec(
        name="image",
        prompt_filename="image_agent_prompt.txt",
        inputs=[InputSpec(label="Image", patterns=["*image.json"])],
        output_suffix="image_agent",
    ),
    CategorySpec(
        name="culture",
        prompt_filename="culture_agent_prompt.txt",
        inputs=[InputSpec(label="Culture", patterns=["*culture.json"])],
        output_suffix="culture_agent",
    ),
]


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format=LOG_FORMAT)


def extract_patient_identifier(input_dir: Path) -> str:
    base_name = sanitize_filename(input_dir.name)
    if base_name.endswith("_json"):
        base_name = base_name[: -len("_json")] or base_name
    return base_name


def looks_like_patient_dir(path: Path) -> bool:
    try:
        return any(child.is_file() and child.suffix.lower() == ".json" for child in path.iterdir())
    except PermissionError:
        logging.warning("Skipping directory without permission: %s", path)
        return False


def collect_patient_directories(entries: Sequence[Path]) -> list[Path]:
    collected: list[Path] = []
    seen: set[Path] = set()

    for entry in entries:
        path = entry.expanduser()
        if not path.exists():
            logging.warning("Skipping non-existent path: %s", path)
            continue
        resolved = path.resolve()
        if not resolved.is_dir():
            logging.warning("Skipping non-directory path: %s", resolved)
            continue
        if resolved.name in {AGENT_OUTPUT_DIR_NAME, SUMMARY_DIR_NAME}:
            logging.debug("Skipping directory %s because it is an outputs folder.", resolved)
            continue
        if "raw" in resolved.name.lower():
            logging.debug("Skipping directory containing 'raw' in name: %s", resolved)
            continue

        candidates = [resolved]
        try:
            candidates.extend(
                child.resolve()
                for child in resolved.iterdir()
                if child.is_dir()
                and "raw" not in child.name.lower()
                and child.name not in {AGENT_OUTPUT_DIR_NAME, SUMMARY_DIR_NAME}
            )
        except PermissionError:
            logging.warning("Unable to enumerate subdirectories for %s", resolved)

        for candidate in candidates:
            if candidate in seen:
                continue
            if "raw" in candidate.name.lower():
                logging.debug("Skipping candidate containing 'raw' in name: %s", candidate)
                continue
            if looks_like_patient_dir(candidate):
                collected.append(candidate)
                seen.add(candidate)
            else:
                logging.debug("Ignoring directory without JSON payloads: %s", candidate)

    collected.sort()
    return collected


def load_json(filepath: Path) -> dict:
    with filepath.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_prompt(filename: str) -> str:
    prompt_path = PROMPTS_DIR / filename
    if not prompt_path.exists():
        raise FileNotFoundError(f"Prompt template not found: {prompt_path}")
    return prompt_path.read_text(encoding="utf-8")


def resolve_input_path(input_dir: Path, patterns: Sequence[str]) -> Path | None:
    for pattern in patterns:
        matches = sorted(input_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def collect_payloads(
    input_dir: Path, spec: CategorySpec, *, force_null_labels: set[str] | None = None
) -> list[PayloadSection]:
    sections: list[PayloadSection] = []
    forced = {label.lower() for label in (force_null_labels or set())}
    for requirement in spec.inputs:
        if requirement.label.lower() in forced:
            sections.append(PayloadSection(label=requirement.label, path=input_dir, data=None))
            continue
        path = resolve_input_path(input_dir, requirement.patterns)
        if path is None:
            if requirement.allow_missing_as_null:
                sections.append(PayloadSection(label=requirement.label, path=input_dir, data=None))
                continue
            logging.warning(
                "Missing %s data for %s inside %s",
                requirement.label,
                spec.name,
                input_dir,
            )
            return []
        data = load_json(path)
        sections.append(PayloadSection(label=requirement.label, path=path, data=data))
    return sections


def build_prompt(template: str, sections: Sequence[PayloadSection]) -> str:
    blocks: list[str] = [template.rstrip(), "", "### JSON payloads ###"]
    for section in sections:
        serialized = json.dumps(section.data, ensure_ascii=False, indent=2)
        blocks.extend(
            [
                f"#### {section.label} ({section.path.name})",
                serialized,
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
        try:
            encoding = tiktoken.encoding_for_model(model)
        except KeyError:
            encoding = tiktoken.get_encoding("cl100k_base")
        return len(encoding.encode(text))
    except Exception as exc:  # Token estimation must never block the LLM request.
        logging.warning("Token estimation unavailable for %s: %s", model, exc)
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


def derive_output_path(
    input_dir: Path, spec: CategorySpec, *, output_suffix_override: str | None = None
) -> Path:
    base_name = extract_patient_identifier(input_dir)
    output_suffix = output_suffix_override or spec.output_suffix
    filename = f"{base_name}_{output_suffix}.json"
    destination_dir = input_dir / AGENT_OUTPUT_DIR_NAME
    destination_dir.mkdir(parents=True, exist_ok=True)
    return destination_dir / filename


def save_result(
    payload: dict, input_dir: Path, spec: CategorySpec, output_path: Path | None = None
) -> Path:
    if output_path is None:
        output_path = derive_output_path(input_dir, spec)
    output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return output_path


def save_failed_output(raw_output: str, input_dir: Path, spec: CategorySpec) -> Path:
    FAILED_DIR.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    identifier = sanitize_filename(f"{input_dir.name}_{spec.name}")
    destination = FAILED_DIR / f"{identifier}_{timestamp}.txt"
    destination.write_text(raw_output, encoding="utf-8")
    return destination


def process_category(
    input_dir: Path,
    spec: CategorySpec,
    *,
    model: str,
    force_null_labels: set[str] | None = None,
    output_suffix_override: str | None = None,
) -> None:
    patient_id = extract_patient_identifier(input_dir)
    sections = collect_payloads(input_dir, spec, force_null_labels=force_null_labels)
    if not sections:
        logging.info("Skipping %s; required inputs not found.", spec.name)
        return

    output_path = derive_output_path(input_dir, spec, output_suffix_override=output_suffix_override)
    if output_path.exists() and not confirm_overwrite(
        output_path, prompt=f"檔案 {output_path} 已存在，是否覆蓋？"
    ):
        logging.info(
            "略過 %s (%s)，因輸出檔案已存在：%s", patient_id, spec.name, output_path
        )
        return

    try:
        template = load_prompt(spec.prompt_filename)
        prompt = build_prompt(template, sections)
        estimated_tokens = estimate_token_count(prompt, model=model)
        if estimated_tokens is not None:
            logging.info("Estimated input tokens for %s: %d", spec.name, estimated_tokens)
        else:
            logging.debug("Token estimation unavailable (tiktoken missing) for %s", spec.name)
        raw_output, usage = send_to_llm(prompt, model=model)
        if usage is not None:
            logging.info(
                "LLM token usage for %s - prompt: %s, completion: %s, total: %s",
                spec.name,
                usage.get("prompt_tokens"),
                usage.get("completion_tokens"),
                usage.get("total_tokens"),
            )
        else:
            logging.debug("Token usage data not returned for %s", spec.name)
        parsed = validate_json(raw_output)
    except Exception as exc:  # pylint: disable=broad-exception-caught
        logging.exception("Failed while processing %s: %s", spec.name, exc)
        return

    if parsed is None:
        failed_path = save_failed_output(raw_output, input_dir, spec)
        logging.error("Stored unparsed output for %s at %s", spec.name, failed_path)
        return

    output_path = save_result(parsed, input_dir, spec, output_path=output_path)
    logging.info("Wrote %s", output_path)
    try:
        dataset_path = record_category_result(patient_id, spec.name, parsed)
        logging.debug(
            "Recorded category analytics for %s (%s) at %s", patient_id, spec.name, dataset_path
        )
    except Exception as exc:  # pragma: no cover - best-effort logging
        logging.exception(
            "Failed to record category analytics for %s (%s): %s", patient_id, spec.name, exc
        )


def run_for_directory(
    input_dir: Path,
    *,
    model: str,
    specs: Sequence[CategorySpec],
    force_null_map: dict[str, set[str]] | None = None,
    output_suffix_overrides: dict[str, str] | None = None,
) -> None:
    force_null_map = force_null_map or {}
    output_suffix_overrides = output_suffix_overrides or {}
    for spec in specs:
        process_category(
            input_dir,
            spec,
            model=model,
            force_null_labels=force_null_map.get(spec.name.lower()),
            output_suffix_override=output_suffix_overrides.get(spec.name.lower()),
        )


def select_category_specs(names: Sequence[str] | None) -> list[CategorySpec]:
    if not names:
        return CATEGORY_SPECS

    normalized = {name.strip().lower() for name in names if name.strip()}
    available = {spec.name.lower(): spec for spec in CATEGORY_SPECS}

    missing = [name for name in normalized if name not in available]
    if missing:
        raise ValueError(f"Unknown categories: {', '.join(sorted(missing))}")

    # Preserve the default CATEGORY_SPECS order while filtering by requested names.
    return [spec for spec in CATEGORY_SPECS if spec.name.lower() in normalized]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate patient JSON files, combine with agent prompts, and call the GPT API.",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="One or more directories containing patient JSON files (or parent directories).",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Model identifier to use when calling the GPT API (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--summarize",
        action="store_true",
        help="Automatically run core.summarize_agent after each patient directory completes.",
    )
    parser.add_argument(
        "--summary-model",
        help=(
            "Model identifier to use for the optional summarization step "
            "(default: reuse the value provided via --model; "
            f"core.summarize_agent default: {SUMMARY_DEFAULT_MODEL})."
        ),
    )
    parser.add_argument(
        "--summary-include-underlying",
        dest="summary_include_underlying",
        action="store_true",
        help="Include underlying agent output in summary prompt (default: excluded).",
    )
    parser.add_argument(
        "--summary-exclude-underlying",
        dest="summary_include_underlying",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--summary-exclude-filmarray",
        action="store_true",
        help="Exclude filmarray/GMTest agent output from summary prompt (default: included).",
    )
    parser.add_argument(
        "--categories",
        nargs="+",
        help=(
            "One or more category names to run (e.g. filmarray_gmtest image). "
            "Defaults to all categories when omitted."
        ),
    )
    parser.add_argument(
        "--exclude-filmarray",
        "--filmarray-as-null",
        dest="exclude_filmarray",
        action="store_true",
        help=(
            "Run only the no_filmarray variant (FilmArray forced to null while keeping GMTest). "
            "By default, both full and no_filmarray variants are produced."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )

    args = parser.parse_args()
    configure_logging(args.verbose)

    try:
        category_specs = select_category_specs(args.categories)
    except ValueError as exc:
        logging.error(str(exc))
        return

    target_dirs = collect_patient_directories(args.inputs)
    if not target_dirs:
        logging.error("No patient directories found in the provided inputs.")
        return

    logging.info("Found %d patient directories to process.", len(target_dirs))
    summary_model = args.summary_model or args.model
    filmarray_spec = next((spec for spec in category_specs if spec.name == "filmarray_gmtest"), None)
    for directory in target_dirs:
        logging.info("Processing %s", directory)
        if args.exclude_filmarray:
            run_for_directory(
                directory,
                model=args.model,
                specs=category_specs,
                force_null_map={"filmarray_gmtest": {"FilmArray"}},
                output_suffix_overrides={"filmarray_gmtest": "filmarray_gmtest_agent_no_filmarray"},
            )
        else:
            # Default behavior: produce full outputs first.
            run_for_directory(
                directory,
                model=args.model,
                specs=category_specs,
            )
            # Also produce no_filmarray variant for filmarray_gmtest without removing GMTest.
            if filmarray_spec is not None:
                process_category(
                    directory,
                    filmarray_spec,
                    model=args.model,
                    force_null_labels={"FilmArray"},
                    output_suffix_override="filmarray_gmtest_agent_no_filmarray",
                )

        if args.summarize:
            logging.info("Summarizing %s via core.summarize_agent.", directory)
            if args.summary_exclude_filmarray:
                summarize_directory(
                    directory,
                    model=summary_model,
                    include_underlying=args.summary_include_underlying,
                    include_filmarray_gmtest=False,
                    no_filmarray_naming=args.exclude_filmarray,
                    filmarray_agent_suffix="filmarray_gmtest_agent",
                )
            elif args.exclude_filmarray:
                summarize_directory(
                    directory,
                    model=summary_model,
                    include_underlying=args.summary_include_underlying,
                    include_filmarray_gmtest=True,
                    no_filmarray_naming=True,
                    filmarray_agent_suffix="filmarray_gmtest_agent_no_filmarray",
                )
            else:
                # Default behavior: emit both full and no_filmarray summaries.
                summarize_directory(
                    directory,
                    model=summary_model,
                    include_underlying=args.summary_include_underlying,
                    include_filmarray_gmtest=True,
                    no_filmarray_naming=False,
                    filmarray_agent_suffix="filmarray_gmtest_agent",
                )
                if filmarray_spec is not None:
                    summarize_directory(
                        directory,
                        model=summary_model,
                        include_underlying=args.summary_include_underlying,
                        include_filmarray_gmtest=True,
                        no_filmarray_naming=True,
                        filmarray_agent_suffix="filmarray_gmtest_agent_no_filmarray",
                    )


if __name__ == "__main__":
    main()

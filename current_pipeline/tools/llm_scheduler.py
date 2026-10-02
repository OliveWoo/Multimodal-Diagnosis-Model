from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from core import summarize_agent
from tools import analyze_llm_agent, mngs_big_agent

LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
DEFAULT_MODEL = "gpt-5"


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format=LOG_FORMAT)


def _union_patient_dirs(
    inputs: Sequence[Path],
    *,
    summary_mode: str,
    mngs_source: str,
) -> list[Path]:
    dirs = set(analyze_llm_agent.collect_patient_directories(inputs))
    dirs.update(
        mngs_big_agent.collect_patient_directories(
            inputs,
            summary_mode=summary_mode,
            mngs_source=mngs_source,
        )
    )
    return sorted(dirs)


def _filter_patients(patient_dirs: Iterable[Path], patient_ids: set[str] | None) -> list[Path]:
    if not patient_ids:
        return sorted(patient_dirs)
    kept: list[Path] = []
    for path in patient_dirs:
        pid = mngs_big_agent.extract_patient_identifier(path)
        if pid in patient_ids:
            kept.append(path)
    return sorted(kept)


def _task_analyze(
    patient_dir: Path,
    *,
    model: str,
    skip_existing: bool,
    overwrite: bool,
) -> None:
    specs = list(analyze_llm_agent.CATEGORY_SPECS)
    if skip_existing:
        specs = [
            spec
            for spec in specs
            if not analyze_llm_agent.derive_output_path(patient_dir, spec).exists()
        ]
        if not specs:
            logging.info("Skip analyze for %s (all outputs exist).", patient_dir.name)
            return
        logging.info(
            "Resume analyze for %s (missing: %s).",
            patient_dir.name,
            ", ".join(spec.name for spec in specs),
        )
    analyze_llm_agent.run_for_directory(
        patient_dir,
        model=model,
        specs=specs,
        overwrite=overwrite,
    )


def _task_mngs_big(
    patient_dir: Path,
    *,
    model: str,
    prompt_path: Path | None,
    ranked_mngs_path: Path | None,
    output_suffix: str,
    skip_existing: bool,
    summary_mode: str,
    mngs_source: str,
    overwrite: bool,
) -> None:
    final_summary_path = mngs_big_agent.find_final_summary_file(
        patient_dir,
        summary_mode=summary_mode,
    )
    if final_summary_path is None:
        logging.info("Skip mngs_big for %s (final summary not found).", patient_dir.name)
        return

    output_path = mngs_big_agent.derive_output_path(
        patient_dir,
        summary_mode=summary_mode,
        output_suffix=output_suffix,
    )
    if skip_existing and output_path.exists():
        logging.info("Skip mngs_big for %s (output exists: %s).", patient_dir.name, output_path)
        return

    mngs_big_agent.process_directory(
        patient_dir,
        model=model,
        prompt_path=prompt_path,
        summary_mode=summary_mode,
        output_suffix=output_suffix,
        ranked_mngs_path=ranked_mngs_path,
        mngs_source=mngs_source,
        overwrite=overwrite,
    )


def _task_summarize(
    patient_dir: Path,
    *,
    model: str,
    skip_existing: bool,
    summary_mode: str,
    overwrite: bool,
) -> None:
    no_filmarray_naming = summary_mode == "no_filmarray"
    filmarray_agent_suffix = (
        "filmarray_gmtest_agent_no_filmarray"
        if no_filmarray_naming
        else "filmarray_gmtest_agent"
    )
    output_path = summarize_agent.derive_output_path(
        patient_dir,
        include_underlying=False,
        include_filmarray_gmtest=True,
        no_filmarray_naming=no_filmarray_naming,
    )
    if skip_existing and output_path.exists():
        logging.info("Skip summarize for %s (output exists: %s).", patient_dir.name, output_path)
        return
    summarize_agent.summarize_directory(
        patient_dir,
        model=model,
        include_underlying=False,
        include_filmarray_gmtest=True,
        no_filmarray_naming=no_filmarray_naming,
        filmarray_agent_suffix=filmarray_agent_suffix,
        overwrite=overwrite,
    )


def _run_task_sequence(
    patient_dir: Path,
    tasks: Sequence[str],
    *,
    model: str,
    prompt_path: Path | None,
    ranked_mngs_path: Path | None,
    output_suffix: str,
    skip_existing: bool,
    summary_mode: str,
    mngs_source: str,
    overwrite: bool,
) -> None:
    for task in tasks:
        if task == "analyze":
            _task_analyze(
                patient_dir,
                model=model,
                skip_existing=skip_existing,
                overwrite=overwrite,
            )
        elif task == "mngs_big":
            _task_mngs_big(
                patient_dir,
                model=model,
                prompt_path=prompt_path,
                ranked_mngs_path=ranked_mngs_path,
                output_suffix=output_suffix,
                skip_existing=skip_existing,
                summary_mode=summary_mode,
                mngs_source=mngs_source,
                overwrite=overwrite,
            )
        elif task == "summarize":
            _task_summarize(
                patient_dir,
                model=model,
                skip_existing=skip_existing,
                summary_mode=summary_mode,
                overwrite=overwrite,
            )
        else:
            logging.warning("Unknown task '%s' skipped.", task)


def _parse_patient_ids(ids: Sequence[str]) -> set[str]:
    normalized: set[str] = set()
    for raw in ids:
        text = raw.strip()
        if not text:
            continue
        if text.startswith("NGS_patient_"):
            normalized.add(text.removesuffix("_json"))
        else:
            normalized.add(f"NGS_patient_{text}")
    return normalized


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run selected LLM tasks for patient directories, optionally in parallel.",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Patient directories or their parent directory.",
    )
    parser.add_argument(
        "--patients",
        nargs="*",
        help="Only process these patients, e.g. 4 6 8 or NGS_patient_4_json.",
    )
    parser.add_argument(
        "--tasks",
        nargs="+",
        choices=["analyze", "mngs_big", "summarize"],
        default=["analyze", "mngs_big"],
        help="Tasks to run. Default: analyze mngs_big.",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"LLM model name (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--prompt",
        type=Path,
        help="Override prompt path for mngs_big.",
    )
    parser.add_argument(
        "--ranked-mngs",
        type=Path,
        default=mngs_big_agent.DEFAULT_RANKED_MNGS_PATH,
        help="Ranked mNGS JSON for mngs_big. Explicit paths override local ranked files.",
    )
    parser.add_argument(
        "--output-suffix",
        default=mngs_big_agent.OUTPUT_SUFFIX,
        help="Output filename suffix for mngs_big. Empty uses the default mode-specific name.",
    )
    parser.add_argument(
        "--summary-mode",
        choices=["auto", "full", "no_filmarray"],
        default="auto",
        help="Which final_summary variant mngs_big should use.",
    )
    parser.add_argument(
        "--mngs-source",
        choices=mngs_big_agent.MNGS_SOURCE_CHOICES,
        default=mngs_big_agent.DEFAULT_MNGS_SOURCE,
        help=(
            "mNGS supplemental source for mngs_big. "
            "grouped uses *_mNGS_grouped.json; all_rk_ntc uses *_all_RK_NTC_microbes.json."
        ),
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=4,
        help="Number of patients to process concurrently.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip outputs that already exist.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs without prompting.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )

    args = parser.parse_args()
    configure_logging(args.verbose)

    all_patient_dirs = _union_patient_dirs(
        args.inputs,
        summary_mode=args.summary_mode,
        mngs_source=args.mngs_source,
    )
    if not all_patient_dirs:
        logging.error("No patient directories found in the provided inputs.")
        return

    patient_ids = _parse_patient_ids(args.patients or [])
    target_dirs = _filter_patients(all_patient_dirs, patient_ids)
    if not target_dirs:
        logging.error("No matching patient directories were selected.")
        return

    logging.info("Processing %d patient directories with tasks: %s", len(target_dirs), ", ".join(args.tasks))

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {
            executor.submit(
                _run_task_sequence,
                patient_dir,
                args.tasks,
                model=args.model,
                prompt_path=args.prompt,
                ranked_mngs_path=args.ranked_mngs,
                output_suffix=args.output_suffix,
                skip_existing=args.skip_existing,
                summary_mode=args.summary_mode,
                mngs_source=args.mngs_source,
                overwrite=args.overwrite,
            ): patient_dir
            for patient_dir in target_dirs
        }
        for future in as_completed(futures):
            patient_dir = futures[future]
            try:
                future.result()
                logging.info("Completed %s", patient_dir.name)
            except Exception as exc:  # pylint: disable=broad-exception-caught
                logging.exception("Failed %s: %s", patient_dir.name, exc)


if __name__ == "__main__":
    main()

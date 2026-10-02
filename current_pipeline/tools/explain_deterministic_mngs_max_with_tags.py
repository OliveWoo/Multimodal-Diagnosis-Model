from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import annotate_deterministic_mngs_max_tags as annotate_tags  # noqa: E402
from tools import explain_deterministic_mngs_max as explain_llm  # noqa: E402


DEFAULT_INPUT_SUFFIX = "mNGS_max_deterministic_opt_chosen_full"
DEFAULT_TAG_SUFFIX = "mNGS_max_deterministic_opt_chosen_full_tag_explained"
DEFAULT_EXPLANATION_SUFFIX = "mNGS_max_deterministic_opt_chosen_full_explanation"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run tag expansion first, then LLM explanation for deterministic mNGS max outputs."
    )
    parser.add_argument(
        "patient_root",
        nargs="?",
        type=Path,
        default=Path("outputs") / "patient_info_rich_normalized",
    )
    parser.add_argument("--patients", nargs="*", help="Optional patient IDs to process.")
    parser.add_argument("--input-json", type=Path, help="Run on one deterministic max JSON file.")
    parser.add_argument("--tag-output", type=Path, help="Tag-expanded JSON output path for --input-json.")
    parser.add_argument("--explanation-output", type=Path, help="Markdown explanation output path for --input-json.")
    parser.add_argument("--input-suffix", default=DEFAULT_INPUT_SUFFIX)
    parser.add_argument("--tag-suffix", default=DEFAULT_TAG_SUFFIX)
    parser.add_argument("--explanation-suffix", default=DEFAULT_EXPLANATION_SUFFIX)
    parser.add_argument("--model", default=explain_llm.DEFAULT_MODEL)
    parser.add_argument("--temperature", type=float, default=None)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--include-excluded", action="store_true")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Create tag-expanded JSON and explanation prompt, but do not call the LLM.",
    )
    return parser.parse_args(argv)


def tag_and_explain_one(
    input_path: Path,
    tag_output_path: Path,
    explanation_output_path: Path,
    *,
    model: str,
    temperature: float | None,
    include_excluded: bool,
    overwrite: bool,
    skip_existing: bool,
    dry_run: bool,
) -> tuple[Path, Path]:
    tag_path = annotate_tags.annotate_file(
        input_path,
        tag_output_path,
        overwrite=overwrite,
        skip_existing=skip_existing,
    )
    explanation_path = explain_llm.explain_one(
        tag_path,
        explanation_output_path,
        model=model,
        temperature=temperature,
        include_excluded=include_excluded,
        overwrite=overwrite,
        skip_existing=skip_existing,
        dry_run=dry_run,
    )
    return tag_path, explanation_path


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.input_json:
        tag_output = args.tag_output or args.input_json.with_name(f"{args.input_json.stem}_tag_explained.json")
        explanation_output = args.explanation_output or args.input_json.with_name(
            f"{args.input_json.stem}_explanation.md"
        )
        tag_path, explanation_path = tag_and_explain_one(
            args.input_json,
            tag_output,
            explanation_output,
            model=args.model,
            temperature=args.temperature,
            include_excluded=args.include_excluded,
            overwrite=args.overwrite,
            skip_existing=args.skip_existing,
            dry_run=args.dry_run,
        )
        print(f"wrote tag_json {tag_path}")
        print(f"wrote explanation {explanation_path}")
        return

    patient_dirs = explain_llm.iter_patient_dirs(args.patient_root, args.patients)
    written = 0
    failures: list[str] = []
    for patient_dir in patient_dirs:
        try:
            patient_id = explain_llm.extract_patient_id(patient_dir)
            summary_dir = patient_dir / "summary_outputs"
            input_path = summary_dir / f"NGS_patient_{patient_id}_{args.input_suffix}.json"
            tag_output_path = summary_dir / f"NGS_patient_{patient_id}_{args.tag_suffix}.json"
            explanation_output_path = summary_dir / f"NGS_patient_{patient_id}_{args.explanation_suffix}.md"
            tag_path, explanation_path = tag_and_explain_one(
                input_path,
                tag_output_path,
                explanation_output_path,
                model=args.model,
                temperature=args.temperature,
                include_excluded=args.include_excluded,
                overwrite=args.overwrite,
                skip_existing=args.skip_existing,
                dry_run=args.dry_run,
            )
            written += 1
            print(f"wrote tag_json {tag_path}")
            print(f"wrote explanation {explanation_path}")
        except Exception as exc:  # noqa: BLE001
            failures.append(f"{patient_dir}: {exc}")
            print(f"failed {patient_dir}: {exc}", file=sys.stderr)

    print(f"written_count={written}")
    if failures:
        print("failures:")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()

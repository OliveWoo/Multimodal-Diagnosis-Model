from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


DEFAULT_MODEL = "gpt-5"
DEFAULT_INPUT_SUFFIX = "mNGS_max_deterministic_opt_chosen_full"
DEFAULT_OUTPUT_SUFFIX = "mNGS_max_deterministic_opt_chosen_full_explanation"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Use an LLM to explain deterministic mNGS max JSON outputs in Traditional Chinese. "
            "The LLM is constrained to explanation only; it must not recalculate levels."
        )
    )
    parser.add_argument(
        "patient_root",
        nargs="?",
        type=Path,
        default=Path("outputs") / "patient_info_rich_normalized",
        help="Root folder containing NGS_patient_<id>_json folders.",
    )
    parser.add_argument("--patients", nargs="*", help="Optional patient IDs to process.")
    parser.add_argument("--input-json", type=Path, help="Explain one deterministic max JSON file.")
    parser.add_argument("--output", type=Path, help="Output markdown path for --input-json mode.")
    parser.add_argument("--input-suffix", default=DEFAULT_INPUT_SUFFIX)
    parser.add_argument("--output-suffix", default=DEFAULT_OUTPUT_SUFFIX)
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument(
        "--temperature",
        type=float,
        default=None,
        help=(
            "Optional temperature. Omit by default because some GPT-5 models do not "
            "support overriding temperature."
        ),
    )
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument(
        "--include-excluded",
        action="store_true",
        help="Include excluded_candidates details in the LLM prompt.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write the compact prompt payload instead of calling the LLM.",
    )
    return parser.parse_args(argv)


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def extract_patient_id(patient_dir: Path) -> str:
    parts = patient_dir.name.split("_")
    if len(parts) >= 3 and parts[0] == "NGS" and parts[1] == "patient":
        return parts[2]
    raise ValueError(f"Cannot extract patient ID from {patient_dir}")


def iter_patient_dirs(patient_root: Path, patients: Sequence[str] | None) -> list[Path]:
    wanted = {str(pid) for pid in patients or []}
    dirs = sorted(
        (path for path in patient_root.glob("NGS_patient_*_json") if path.is_dir()),
        key=lambda path: int(extract_patient_id(path)),
    )
    if wanted:
        dirs = [path for path in dirs if extract_patient_id(path) in wanted]
    return dirs


def build_paths_for_patient(
    patient_dir: Path, *, input_suffix: str, output_suffix: str
) -> tuple[Path, Path]:
    patient_id = extract_patient_id(patient_dir)
    summary_dir = patient_dir / "summary_outputs"
    input_path = summary_dir / f"NGS_patient_{patient_id}_{input_suffix}.json"
    output_path = summary_dir / f"NGS_patient_{patient_id}_{output_suffix}.md"
    return input_path, output_path


def compact_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    return {
        "organism_name": candidate.get("organism_name"),
        "classification": candidate.get("classification"),
        "integrated_causative_level": candidate.get("integrated_causative_level"),
        "mngs_signal_tier": candidate.get("mngs_signal_tier"),
        "rank_priority": candidate.get("rank_priority"),
        "rank_rule": candidate.get("rank_rule"),
        "reads": candidate.get("reads"),
        "reads_tier": candidate.get("reads_tier"),
        "reads_percentile": candidate.get("reads_percentile"),
        "dominance_tier": candidate.get("dominance_tier"),
        "specimen_alignment": candidate.get("specimen_alignment"),
        "source_category": candidate.get("source_category"),
        "is_protected_pathogen": candidate.get("is_protected_pathogen"),
        "protected_retention_condition": candidate.get("protected_retention_condition"),
        "is_likely_colonizer_or_background": candidate.get("is_likely_colonizer_or_background"),
        "module_support_summary": candidate.get("module_support_summary"),
        "applied_rules": candidate.get("applied_rules"),
        "key_evidence": candidate.get("key_evidence"),
        "integrated_reasoning": candidate.get("integrated_reasoning"),
    }


def compact_picked(picked: dict[str, Any]) -> dict[str, Any]:
    return {
        "organism_name": picked.get("organism_name"),
        "classification": picked.get("classification"),
        "picked_role": picked.get("picked_role"),
        "basis_level": picked.get("basis_level"),
        "mngs_signal_tier": picked.get("mngs_signal_tier"),
        "rank_priority": picked.get("rank_priority"),
        "reads": picked.get("reads"),
        "why_picked": picked.get("why_picked"),
        "caution_flags": picked.get("caution_flags"),
    }


def build_compact_payload(data: dict[str, Any], *, include_excluded: bool) -> dict[str, Any]:
    best = data.get("best_available_summary", {})
    payload: dict[str, Any] = {
        "rule_version": data.get("rule_version"),
        "final_infection_likelihood": data.get("final_infection_likelihood"),
        "dominant_source": data.get("dominant_source"),
        "dominant_pathogen_type": data.get("dominant_pathogen_type"),
        "host_context": data.get("host_context"),
        "best_available_summary": {
            "selection_mode": best.get("selection_mode"),
            "picked_count": best.get("picked_count"),
            "picked_pathogens": [
                compact_picked(item)
                for item in best.get("picked_pathogens", [])
                if isinstance(item, dict)
            ],
        },
        "pathogen_candidates": [
            compact_candidate(item)
            for item in data.get("pathogen_candidates", [])
            if isinstance(item, dict)
        ],
        "data_gaps": data.get("data_gaps"),
        "source_files": data.get("source_files"),
    }
    if include_excluded:
        payload["excluded_candidates"] = data.get("excluded_candidates", [])
    else:
        payload["excluded_candidates_summary"] = {
            "count": len(data.get("excluded_candidates", []) or []),
            "organisms": [
                item.get("organism_name")
                for item in data.get("excluded_candidates", [])
                if isinstance(item, dict)
            ],
        }
    return payload


def build_prompt(input_path: Path, payload: dict[str, Any]) -> str:
    patient_label = input_path.stem
    payload_text = json.dumps(payload, ensure_ascii=False, indent=2)
    return f"""你是一位 mNGS deterministic scoring 報告解釋助手。

任務：請用繁體中文解釋以下 deterministic mNGS max JSON 的內容，讓臨床研究者可以快速理解。

硬性規則：
- 只能解釋 JSON 已經存在的結果，不可以重新計算 Level。
- 不可以新增 JSON 裡沒有的病原。
- 不可以把 Level 3/4/5 說成確診。
- 不可以提供治療建議或用藥建議。
- 如果資料顯示 Evidence_limited、Possible_reactivation、Not_recommended_as_sole_treatment_basis，必須明確說明。
- 若出現 protected pathogen，請說明「保留觀察」而不是「確認致病」。
- 若出現 background/colonizer，請說明它被視為背景或定殖的原因。
- 不要輸出固定 tag glossary；固定 tag 對照已放在 README_mngs_deterministic_tag_glossary.md。
- 若需要提到 tag，請直接用白話整合在句子中，例如「O1 需要考慮機會感染」、「V2 中度免疫低下或重要共病」。
- 輸出只用繁體中文。

請固定用以下格式輸出 Markdown：

# {patient_label} 解釋報告

## 1. 總結
用 3-5 點條列說明 final_infection_likelihood、dominant_source、dominant_pathogen_type、host_context、selection_mode。

## 2. Best Available 病原
用表格列出每個 picked pathogen：病原、類型、角色、Level、mNGS signal tier、rank_priority、reads、主要原因、注意事項。

## 3. Level 與證據解讀
依 Level 1/2/3/4/5 分段解釋本案有哪些病原，以及這些 Level 的臨床解讀限制。

## 4. 重要注意事項
列出低 reads、protected pathogen、possible reactivation、background/colonizer、資料缺口等需要注意的點。

## 5. 一句話結論
用一句話總結，但必須保留不確定性。

輸入 JSON：
```json
{payload_text}
```
"""
def extract_output_text(response: Any) -> str:
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


def send_to_llm(prompt: str, *, model: str, temperature: float | None) -> tuple[str, dict[str, int | None] | None]:
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY environment variable is not set.")
    from openai import OpenAI  # type: ignore

    client = OpenAI(api_key=api_key)
    kwargs: dict[str, Any] = {"model": model, "input": prompt}
    if temperature is not None:
        kwargs["temperature"] = temperature
    response = client.responses.create(**kwargs)
    text = extract_output_text(response)
    if not text:
        raise RuntimeError("OpenAI response did not contain any content.")
    usage = getattr(response, "usage", None)
    usage_summary = None
    if usage is not None:
        usage_summary = {
            "prompt_tokens": getattr(usage, "input_tokens", None),
            "completion_tokens": getattr(usage, "output_tokens", None),
            "total_tokens": getattr(usage, "total_tokens", None),
        }
    return text, usage_summary


def explain_one(
    input_path: Path,
    output_path: Path,
    *,
    model: str,
    temperature: float | None,
    include_excluded: bool,
    overwrite: bool,
    skip_existing: bool,
    dry_run: bool,
) -> Path:
    if output_path.exists() and skip_existing:
        return output_path
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output_path}. Use --overwrite or --skip-existing.")
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    data = load_json(input_path)
    if not isinstance(data, dict):
        raise TypeError(f"Expected dict JSON from {input_path}")
    payload = build_compact_payload(data, include_excluded=include_excluded)
    prompt = build_prompt(input_path, payload)
    if dry_run:
        write_text(output_path.with_suffix(".prompt.md"), prompt)
        return output_path.with_suffix(".prompt.md")

    explanation, usage = send_to_llm(prompt, model=model, temperature=temperature)
    if usage:
        explanation += (
            "\n\n---\n"
            f"LLM token usage: prompt={usage.get('prompt_tokens')}, "
            f"completion={usage.get('completion_tokens')}, total={usage.get('total_tokens')}\n"
        )
    write_text(output_path, explanation)
    return output_path


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.input_json:
        output = args.output or args.input_json.with_name(f"{args.input_json.stem}_explanation.md")
        written = explain_one(
            args.input_json,
            output,
            model=args.model,
            temperature=args.temperature,
            include_excluded=args.include_excluded,
            overwrite=args.overwrite,
            skip_existing=args.skip_existing,
            dry_run=args.dry_run,
        )
        print(f"wrote {written}")
        return

    patient_dirs = iter_patient_dirs(args.patient_root, args.patients)
    tasks = [
        (*build_paths_for_patient(path, input_suffix=args.input_suffix, output_suffix=args.output_suffix),)
        for path in patient_dirs
    ]
    written: list[Path] = []
    failures: list[str] = []

    with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as executor:
        futures = {
            executor.submit(
                explain_one,
                input_path,
                output_path,
                model=args.model,
                temperature=args.temperature,
                include_excluded=args.include_excluded,
                overwrite=args.overwrite,
                skip_existing=args.skip_existing,
                dry_run=args.dry_run,
            ): (input_path, output_path)
            for input_path, output_path in tasks
        }
        for future in as_completed(futures):
            input_path, _ = futures[future]
            try:
                written_path = future.result()
                written.append(written_path)
                print(f"wrote {written_path}")
            except Exception as exc:  # noqa: BLE001
                failures.append(f"{input_path}: {exc}")
                print(f"failed {input_path}: {exc}", file=sys.stderr)

    print(f"written_count={len(written)}")
    if failures:
        print("failures:")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()

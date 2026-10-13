"""Validate an answer-blind small-agent workspace before API use."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

from tools import analyze_llm_agent


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_ROOT = ROOT / "outputs" / "runs" / "2026-10-09_44patients_small_agents_stage1_v5"
def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def validate(
    run_root: Path,
    model: str,
    max_input_tokens: int,
    expected_patients: int | None = 44,
) -> dict[str, Any]:
    patients_root = run_root / "patients"
    patient_dirs = analyze_llm_agent.collect_patient_directories([patients_root])
    errors: list[str] = []
    prompt_rows: list[dict[str, Any]] = []
    if expected_patients is not None and len(patient_dirs) != expected_patients:
        errors.append(f"patient_count={len(patient_dirs)} expected={expected_patients}")

    forbidden = [
        str(path.resolve())
        for path in patients_root.rglob("*")
        if "ground_truth" in path.name.lower() or "clinical_pathogen" in path.name.lower()
    ]
    if forbidden:
        errors.append(f"forbidden answer files present: {forbidden}")

    for patient_dir in patient_dirs:
        patient_name = patient_dir.name.removesuffix("_json")
        manifest_path = patient_dir / f"{patient_name}_small_agent_input_manifest.json"
        if not manifest_path.is_file():
            errors.append(f"missing patient manifest: {manifest_path}")
            continue
        manifest = load_json(manifest_path)
        if manifest.get("ground_truth_accessed") is not False:
            errors.append(f"ground-truth firewall flag failed: {manifest_path}")
        provenance_path = Path(manifest.get("provenance_path") or "")
        if not provenance_path.is_file():
            errors.append(f"missing provenance sidecar: {provenance_path}")
            provenance_modules: dict[str, Any] = {}
        else:
            provenance_modules = (load_json(provenance_path).get("modules") or {})

        for module, entries in provenance_modules.items():
            refs = [entry.get("source_ref") for entry in entries if isinstance(entry, dict)]
            if len(refs) != len(set(refs)):
                errors.append(f"duplicate source_ref: {patient_name}/{module}")

        for spec in analyze_llm_agent.CATEGORY_SPECS:
            sections = analyze_llm_agent.collect_payloads(patient_dir, spec)
            if not sections:
                errors.append(f"missing required prompt input: {patient_name}/{spec.name}")
                continue
            prompt = analyze_llm_agent.build_prompt(
                analyze_llm_agent.load_prompt(spec.prompt_filename), sections
            )
            token_count = analyze_llm_agent.estimate_token_count(prompt, model=model)
            prompt_rows.append({
                "patient": patient_name,
                "category": spec.name,
                "characters": len(prompt),
                "estimated_input_tokens": token_count,
            })
            if token_count is not None and token_count > max_input_tokens:
                errors.append(
                    f"prompt exceeds token limit: {patient_name}/{spec.name} "
                    f"tokens={token_count} limit={max_input_tokens}"
                )
            for section in sections:
                if not isinstance(section.data, list):
                    continue
                module_name = section.path.stem.removeprefix(f"{patient_name}_")
                valid_refs = {
                    entry.get("source_ref")
                    for entry in provenance_modules.get(module_name, [])
                    if isinstance(entry, dict)
                }
                for index, record in enumerate(section.data):
                    if not isinstance(record, dict):
                        continue
                    source_ref = record.get("source_ref")
                    if source_ref not in valid_refs:
                        errors.append(
                            f"unresolved source_ref: {patient_name}/{module_name}/{index}: {source_ref}"
                        )

    prompt_rows.sort(
        key=lambda row: row.get("estimated_input_tokens") or 0,
        reverse=True,
    )
    expected_prompts = (
        expected_patients * len(analyze_llm_agent.CATEGORY_SPECS)
        if expected_patients is not None
        else len(patient_dirs) * len(analyze_llm_agent.CATEGORY_SPECS)
    )
    if len(prompt_rows) != expected_prompts:
        errors.append(f"prompt_count={len(prompt_rows)} expected={expected_prompts}")
    result = {
        "schema_version": "small_agent_workspace_validation.v1",
        "run_root": str(run_root.resolve()),
        "model_for_token_estimate": model,
        "max_input_tokens": max_input_tokens,
        "patient_count": len(patient_dirs),
        "expected_patient_count": expected_patients,
        "prompt_count": len(prompt_rows),
        "expected_prompt_count": expected_prompts,
        "ground_truth_named_file_count": len(forbidden),
        "max_prompt": prompt_rows[0] if prompt_rows else None,
        "prompts_over_limit": sum(
            (row.get("estimated_input_tokens") or 0) > max_input_tokens
            for row in prompt_rows
        ),
        "status": "passed" if not errors else "failed",
        "errors": errors,
        "prompts": prompt_rows,
    }
    output_path = run_root / "workspace_validation.json"
    output_path.write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return result


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--model", default="gpt-5.6-luna")
    parser.add_argument("--max-input-tokens", type=int, default=50_000)
    parser.add_argument("--expected-patients", type=int, default=44)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    result = validate(
        args.run_root,
        args.model,
        args.max_input_tokens,
        args.expected_patients,
    )
    print(json.dumps({key: value for key, value in result.items() if key != "prompts"}, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())

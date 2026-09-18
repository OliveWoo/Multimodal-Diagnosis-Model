#!/usr/bin/env python3
"""Audit completed rationale runs and render an evidence-linked reading copy.

No API calls. Original model outputs and frozen selection inputs remain unchanged.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from generate_llm_rationales import build_payload, load_json, sha256_text, validate_output, write_json


def link(label: str, path: Path) -> str:
    return f"[{label}](<{path.resolve().as_posix()}>)"


def render_patient(result: dict, output_file: Path) -> list[str]:
    generated = result["generated_rationale"]
    lines = [
        f"## {result['hospital']} P{result['patient_id']}",
        "",
        f"檢驗編號：{'、'.join(result['test_identifiers']) or '未提供'}",
        "",
        generated["patient_result_summary"],
        "",
        link("完整 LLM 輸出與執行紀錄", output_file) + " ｜ "
        + link("輸入病例、各候選病原與原始證據", Path(result["input_source"])),
        "",
    ]
    for item in sorted(generated["selected_pathogen_explanations"], key=lambda x: x["rank"]):
        lines += [f"### {item['rank']}. {item['organism']}", "", item["concise_reason"], ""]
        lines += ["關鍵證據：" + item["key_evidence_summary"], ""]
        if item["limitations"]:
            lines += ["限制：", ""]
            lines += [
                "- " + limitation["text"]
                + ("（" + "、".join(limitation["evidence_ids"]) + "）" if limitation["evidence_ids"] else "")
                for limitation in item["limitations"]
            ]
            lines.append("")
        if item["clinician_review_question"]:
            lines += ["待醫師確認：" + item["clinician_review_question"], ""]
        ids = list(dict.fromkeys(item["evidence_ids"] + [
            eid for limitation in item["limitations"] for eid in limitation["evidence_ids"]
        ]))
        if ids:
            lines += ["證據定位（相對原始病例資料；細節可由上方輸入病例連結核對）：", ""]
            lines += [f"- {eid}：`{result['evidence_map'][eid]}`" for eid in ids]
            lines.append("")
    if generated["no_selected_pathogen"]:
        lines += ["此病例保留『無入選病原』，並非漏跑。", "", "引用證據：", ""]
        lines += [f"- {eid}：`{result['evidence_map'][eid]}`" for eid in generated["patient_level_evidence_ids"]]
        lines.append("")
    return lines


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument(
        "--input-root",
        type=Path,
        help=(
            "Legacy parent containing <cohort>/run/rationale_views/"
            "traceable_reasoning_chains. Not required when every cohort has "
            "an explicit --input-dir mapping."
        ),
    )
    parser.add_argument(
        "--cohorts",
        nargs="+",
        default=["KH", "two_hospital"],
        help="Cohort names to audit (default preserves the original two-cohort run).",
    )
    parser.add_argument(
        "--input-dir",
        action="append",
        default=[],
        metavar="COHORT=PATH",
        help="Explicit reasoning-chain directory for one cohort; repeat as needed.",
    )
    parser.add_argument("--report-root", type=Path, required=True, help="New or empty folder")
    parser.add_argument("--max-excluded-candidates", type=int, default=5)
    parser.add_argument("--retry-run", action="append", default=[], help="COHORT=retry_output_folder; only fills absent patients")
    args = parser.parse_args()
    if args.report_root.exists() and any(args.report_root.iterdir()):
        raise SystemExit("Report directory must be new or empty")
    reports = []
    all_lines = [
        "# 全部病例：LLM 選菌理由閱讀版", "",
        "以下文字為模型根據已凍結選菌及所提供病例證據生成的事後說明，不是原模型的內部思考紀錄，也不是新的診斷。",
        "",
        "本輪不新增 PubMed 文獻；未選候選最多提供前 5 項給理由產生器。完整候選及既有規則說明仍見每位病人的輸入病例連結。",
        "",
        "結構檢查通過不代表所有醫學敘述或引用相關性均經人工確認；此文件待臨床複核。", "",
    ]
    summaries = []
    errors = []
    token_usage = Counter()
    retry_dirs: dict[str, list[Path]] = {}
    input_dirs: dict[str, Path] = {}
    for value in args.input_dir:
        if "=" not in value:
            raise SystemExit("--input-dir must be COHORT=PATH")
        cohort_name, path_value = value.split("=", 1)
        if not cohort_name or cohort_name in input_dirs:
            raise SystemExit(f"Duplicate or empty --input-dir cohort: {cohort_name!r}")
        input_dirs[cohort_name] = Path(path_value).resolve()
    cohorts = list(dict.fromkeys(args.cohorts))
    if len(cohorts) != len(args.cohorts):
        raise SystemExit("--cohorts contains duplicate names")
    unknown_inputs = sorted(set(input_dirs) - set(cohorts))
    if unknown_inputs:
        raise SystemExit(f"--input-dir names are absent from --cohorts: {unknown_inputs}")
    if args.input_root is None and set(input_dirs) != set(cohorts):
        raise SystemExit("Pass --input-root or an --input-dir mapping for every cohort")
    accepted_results = []
    for value in args.retry_run:
        cohort_name, path_value = value.split("=", 1)
        if cohort_name not in cohorts:
            raise SystemExit(f"Unknown cohort: {cohort_name}")
        retry_dirs.setdefault(cohort_name, []).append(Path(path_value))
    for cohort in cohorts:
        run_dir = args.run_root / cohort
        manifest = load_json(run_dir / "llm_rationale_manifest.json")
        input_dir = input_dirs.get(cohort)
        if input_dir is None:
            assert args.input_root is not None
            input_dir = (
                args.input_root
                / cohort
                / "run"
                / "rationale_views"
                / "traceable_reasoning_chains"
            )
        if not input_dir.is_dir():
            raise SystemExit(f"Input directory is missing for {cohort}: {input_dir}")
        inputs = {str(load_json(path)["病例識別"]["病人編號"]): path for path in input_dir.glob("*_explainable_reasoning_zh.json")}
        outputs = list((run_dir / "patient_rationales").glob("*_llm_rationale.json"))
        base_output_count = len(outputs)
        manifests_by_dir = {run_dir.resolve(): manifest}
        attempt_failures = [{"attempt_dir": str(run_dir.resolve()), **failure} for failure in manifest["failures"]]
        for retry_dir in retry_dirs.get(cohort, []):
            retry_manifest = load_json(retry_dir / "llm_rationale_manifest.json")
            manifests_by_dir[retry_dir.resolve()] = retry_manifest
            outputs.extend((retry_dir / "patient_rationales").glob("*_llm_rationale.json"))
            attempt_failures.extend({"attempt_dir": str(retry_dir.resolve()), **failure} for failure in retry_manifest["failures"])
        outputs.sort(key=lambda p: int(p.name.split("_")[2]))
        seen = []
        selected_count = 0
        no_picks = []
        cohort_usage = Counter()
        for path in outputs:
            result = load_json(path)
            patient_id = result["patient_id"]
            seen.append(patient_id)
            entry = {"cohort": cohort, "patient_id": patient_id, "output_file": str(path.resolve()), "errors": []}
            try:
                result_manifest = manifests_by_dir[path.parent.parent.resolve()]
                expected_path = inputs[patient_id].resolve()
                if Path(result["input_source"]).resolve() != expected_path:
                    raise ValueError("wrong input source path")
                payload, evidence_map = build_payload(load_json(expected_path), args.max_excluded_candidates)
                if result["evidence_map"] != evidence_map:
                    raise ValueError("evidence map changed")
                expected_hash = sha256_text(json.dumps(payload, ensure_ascii=False, sort_keys=True))
                if result["input_sha256"] != expected_hash:
                    raise ValueError("input hash mismatch")
                if result["prompt_sha256"] != result_manifest["prompt_sha256"]:
                    raise ValueError("prompt hash mismatch")
                validate_output(result["generated_rationale"], payload, evidence_map)
                generated = result["generated_rationale"]["selected_pathogen_explanations"]
                expected = Counter((x["organism"], x["frozen_rank"]) for x in payload["candidates_available_to_reason_generator"] if x["frozen_decision"] == "選入")
                actual = Counter((x["organism"], x["rank"]) for x in generated)
                if expected != actual:
                    raise ValueError("selection multiset mismatch, including duplicate items")
                if not result["response_id"] or result["model_returned"] != result_manifest["model"]:
                    raise ValueError("missing response ID or model mismatch")
                cited = result["generated_rationale"]["patient_level_evidence_ids"] + [
                    eid for x in generated for eid in x["evidence_ids"]
                ] + [eid for x in generated for lim in x["limitations"] for eid in lim["evidence_ids"]]
                entry.update({
                    "passed": True, "selected_count": len(generated),
                    "unique_cited_evidence_ids": len(set(cited)),
                    "requires_clinician_review_count": sum(x["requires_clinician_review"] for x in generated),
                    "missing_test_id": not result["test_identifiers"],
                    "hospital_pending": "待補" in result["hospital"],
                    "selected_without_own_evidence": [x["organism"] for x in payload["candidates_available_to_reason_generator"] if x["frozen_decision"] == "選入" and not x["evidence"]],
                    "prompt_sha256": result["prompt_sha256"],
                })
                selected_count += len(generated)
                if result["generated_rationale"]["no_selected_pathogen"]:
                    no_picks.append(patient_id)
                usage = result["usage"] or {}
                for key in ("input_tokens", "output_tokens", "total_tokens"):
                    cohort_usage[key] += usage.get(key, 0)
                cohort_usage["cached_input_tokens"] += (usage.get("input_tokens_details") or {}).get("cached_tokens", 0)
                cohort_usage["cache_write_input_tokens"] += (usage.get("input_tokens_details") or {}).get("cache_write_tokens", 0)
                cohort_usage["reasoning_tokens"] += (usage.get("output_tokens_details") or {}).get("reasoning_tokens", 0)
                all_lines += render_patient(result, path)
                accepted_results.append({"cohort": cohort, "output_file": str(path.resolve()), "result": result})
            except Exception as exc:
                entry["passed"] = False
                entry["errors"].append(str(exc))
                errors.append({"cohort": cohort, "patient_id": patient_id, "error": str(exc)})
            reports.append(entry)
        if Counter(seen) != Counter(inputs.keys()):
            errors.append({"cohort": cohort, "error": "patient coverage mismatch", "missing": sorted(set(inputs) - set(seen)), "extra": sorted(set(seen) - set(inputs))})
        if manifest["summary"]["successful_patients"] != base_output_count:
            errors.append({"cohort": cohort, "error": "inconsistent base manifest output count"})
        unresolved = [failure for failure in attempt_failures if failure["patient_id"] not in seen]
        if unresolved:
            errors.append({"cohort": cohort, "error": "unresolved generation failures", "failures": unresolved})
        summaries.append({"cohort": cohort, "expected_patients": len(inputs), "output_patients": len(outputs), "selected_pathogen_reasons": selected_count, "no_pick_patients": no_picks, "attempt_failures": attempt_failures, "usage": dict(cohort_usage)})
        token_usage.update(cohort_usage)
    args.report_root.mkdir(parents=True, exist_ok=True)
    audit = {"created_at_utc": datetime.now(timezone.utc).isoformat(), "audit_type": "structural_not_clinical", "all_passed": not errors, "cohorts": summaries, "usage": dict(token_usage), "usage_scope": "Accepted responses only; failed attempts without stored usage are not included", "errors": errors, "patients": reports}
    write_json(args.report_root / "batch_audit.json", audit)
    write_json(args.report_root / "all_patient_rationales.json", {"status": "structurally_validated" if not errors else "audit_failed", "clinical_review_status": "not_expert_validated", "patients": accepted_results})
    (args.report_root / "all_patients_readable_zh.md").write_text("\n".join(all_lines) + "\n", encoding="utf-8")
    print(json.dumps({"all_passed": audit["all_passed"], "cohorts": summaries, "errors": errors, "report_root": str(args.report_root.resolve())}, ensure_ascii=False, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Generate evidence-bound, clinician-readable reasons for frozen OBER picks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PROMPT = SCRIPT_DIR / "prompts" / "llm_reason_generator_zh_v1.txt"

OUTPUT_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "patient_id": {"type": "string"},
        "no_selected_pathogen": {"type": "boolean"},
        "patient_result_summary": {"type": "string"},
        "patient_level_evidence_ids": {
            "type": "array",
            "items": {"type": "string"},
        },
        "selected_pathogen_explanations": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "organism": {"type": "string"},
                    "rank": {"type": "integer"},
                    "concise_reason": {"type": "string"},
                    "key_evidence_summary": {"type": "string"},
                    "evidence_ids": {
                        "type": "array",
                        "items": {"type": "string"},
                    },
                    "limitations": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "additionalProperties": False,
                            "properties": {
                                "text": {"type": "string"},
                                "evidence_ids": {
                                    "type": "array",
                                    "items": {"type": "string"},
                                },
                            },
                            "required": ["text", "evidence_ids"],
                        },
                    },
                    "requires_clinician_review": {"type": "boolean"},
                    "clinician_review_question": {"type": "string"},
                },
                "required": [
                    "organism",
                    "rank",
                    "concise_reason",
                    "key_evidence_summary",
                    "evidence_ids",
                    "limitations",
                    "requires_clinician_review",
                    "clinician_review_question",
                ],
            },
        },
    },
    "required": [
        "patient_id",
        "no_selected_pathogen",
        "patient_result_summary",
        "patient_level_evidence_ids",
        "selected_pathogen_explanations",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use an OpenAI model to explain frozen OBER selections without changing them."
    )
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--patients", required=True, help="Comma-separated patient IDs")
    parser.add_argument("--env-file", type=Path)
    parser.add_argument("--prompt", type=Path, default=DEFAULT_PROMPT)
    parser.add_argument("--model", default="gpt-5.6-luna")
    parser.add_argument("--reasoning-effort", default="low")
    parser.add_argument("--max-output-tokens", type=int, default=4000)
    parser.add_argument("--max-excluded-candidates", type=int, default=5)
    parser.add_argument("--retries", type=int, default=1)
    parser.add_argument("--pretty", action="store_true")
    return parser.parse_args()


def compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, value: Any, pretty: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(
            value,
            handle,
            ensure_ascii=False,
            indent=2 if pretty else None,
            separators=None if pretty else (",", ":"),
        )
        handle.write("\n")


def load_env_file(path: Path) -> None:
    for raw_line in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if key.startswith("$env:"):
            key = key[5:]
        value = value.strip().strip('"').strip("'")
        if key and key not in os.environ:
            os.environ[key] = value


def requested_patient_ids(raw: str) -> list[str]:
    values = [item.strip() for item in re.split(r"[,;\s]+", raw) if item.strip()]
    if not values:
        raise SystemExit("--patients did not contain any patient ID")
    return list(dict.fromkeys(values))


def patient_file(input_root: Path, patient_id: str) -> Path:
    direct = input_root / f"NGS_patient_{patient_id}_explainable_reasoning_zh.json"
    if direct.is_file():
        return direct
    matches = list(input_root.rglob(f"NGS_patient_{patient_id}_explainable_reasoning_zh.json"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected exactly one readable reasoning chain for patient {patient_id}; found {len(matches)}"
        )
    return matches[0]


def build_payload(chain: dict[str, Any], max_excluded: int) -> tuple[dict[str, Any], dict[str, str]]:
    patient = chain.get("病例識別") or {}
    overview = chain.get("一眼看懂") or {}
    candidates = [
        item
        for item in chain.get("步驟2至5_逐候選病原判讀", [])
        if isinstance(item, dict)
    ]
    selected = [item for item in candidates if item.get("最終處置") == "選入"]
    excluded = [item for item in candidates if item.get("最終處置") != "選入"]
    included = selected + excluded[: max(0, max_excluded)]
    evidence_map: dict[str, str] = {}
    candidate_payloads: list[dict[str, Any]] = []
    for candidate_index, candidate in enumerate(included, start=1):
        evidence_payloads: list[dict[str, Any]] = []
        for evidence_index, evidence in enumerate(
            candidate.get("病例內支持或反對證據", []), start=1
        ):
            if not isinstance(evidence, dict):
                continue
            evidence_id = f"C{candidate_index:02d}-E{evidence_index:02d}"
            locator = compact(evidence.get("原始位置")) or "未提供"
            evidence_map[evidence_id] = locator
            evidence_payloads.append(
                {
                    "evidence_id": evidence_id,
                    "summary": evidence.get("證據摘要"),
                    "interpretation": evidence.get("如何解讀"),
                    "source_locator": locator,
                    "checkable_fields": evidence.get("可核對欄位"),
                }
            )
        candidate_payloads.append(
            {
                "organism": candidate.get("病原"),
                "frozen_decision": candidate.get("最終處置"),
                "frozen_rank": candidate.get("最終排序"),
                "upstream_conclusion": candidate.get("一句話結論"),
                "why_considered": candidate.get("為什麼會被考慮"),
                "why_selected_or_excluded": candidate.get("為什麼選入或排除"),
                "key_features": candidate.get("關鍵判斷資訊"),
                "known_limitations": candidate.get("重要限制"),
                "evidence": evidence_payloads,
            }
        )
    payload = {
        "task": "explain_frozen_ober_selection",
        "patient_id": str(patient.get("病人編號") or ""),
        "hospital": patient.get("醫院"),
        "test_identifiers": (chain.get("步驟1_病人資料摘要") or {}).get("檢驗編號"),
        "patient_context": chain.get("步驟1_病人資料摘要"),
        "frozen_decision": {
            "selected_pathogens": overview.get("最終入選病原"),
            "no_selected_pathogen": int(overview.get("入選病原數") or 0) == 0,
            "selected_count": int(overview.get("入選病原數") or 0),
            "warning": overview.get("重要提醒"),
        },
        "candidates_available_to_reason_generator": candidate_payloads,
    }
    return payload, evidence_map


def validate_output(
    output: dict[str, Any],
    payload: dict[str, Any],
    evidence_map: dict[str, str],
) -> None:
    expected_patient = payload["patient_id"]
    if output.get("patient_id") != expected_patient:
        raise ValueError(
            f"patient_id mismatch: expected {expected_patient}, got {output.get('patient_id')}"
        )
    expected_no_selected = bool(payload["frozen_decision"]["no_selected_pathogen"])
    if bool(output.get("no_selected_pathogen")) != expected_no_selected:
        raise ValueError("no_selected_pathogen differs from the frozen decision")
    expected_selected = {
        item["organism"]: item["frozen_rank"]
        for item in payload["candidates_available_to_reason_generator"]
        if item["frozen_decision"] == "選入"
    }
    generated = output.get("selected_pathogen_explanations")
    if not isinstance(generated, list):
        raise ValueError("selected_pathogen_explanations is not a list")
    generated_selected = {item.get("organism"): item.get("rank") for item in generated}
    if len(generated_selected) != len(generated):
        raise ValueError("selected_pathogen_explanations contains duplicate organisms")
    if generated_selected != expected_selected:
        raise ValueError(
            f"selected pathogens/ranks changed: expected {expected_selected}, got {generated_selected}"
        )
    allowed_ids = set(evidence_map)
    cited_ids = list(output.get("patient_level_evidence_ids") or [])
    for item in generated:
        cited_ids.extend(item.get("evidence_ids") or [])
        for limitation in item.get("limitations") or []:
            cited_ids.extend(limitation.get("evidence_ids") or [])
    invalid_ids = sorted(set(cited_ids) - allowed_ids)
    if invalid_ids:
        raise ValueError(f"unknown evidence IDs: {invalid_ids}")
    if expected_no_selected and allowed_ids and not output.get("patient_level_evidence_ids"):
        raise ValueError("no-pick explanation has available evidence but cited none")
    evidence_by_organism = {
        item["organism"]: {evidence["evidence_id"] for evidence in item["evidence"]}
        for item in payload["candidates_available_to_reason_generator"]
    }
    for item in generated:
        organism = item["organism"]
        own_ids = evidence_by_organism.get(organism, set())
        item_ids = list(item.get("evidence_ids") or [])
        for limitation in item.get("limitations") or []:
            item_ids.extend(limitation.get("evidence_ids") or [])
        foreign_ids = sorted(set(item_ids) - own_ids)
        if foreign_ids:
            raise ValueError(f"{organism} cited another candidate's evidence: {foreign_ids}")
        if own_ids and not item_ids:
            raise ValueError(f"{organism} has evidence but cited none")


def usage_dict(response: Any) -> dict[str, Any] | None:
    usage = getattr(response, "usage", None)
    if usage is None:
        return None
    if hasattr(usage, "model_dump"):
        return usage.model_dump()
    return usage if isinstance(usage, dict) else {"value": str(usage)}


def call_model(
    *,
    client: Any,
    prompt: str,
    payload: dict[str, Any],
    model: str,
    reasoning_effort: str,
    max_output_tokens: int,
    retries: int,
) -> Any:
    request_input = (
        "以下是單一病人的凍結 OBER 決策與可回查病例證據。"
        "請只解釋，不得改變決策。\n\n"
        + json.dumps(payload, ensure_ascii=False, indent=2)
    )
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            return client.responses.create(
                model=model,
                instructions=prompt,
                input=request_input,
                reasoning={"effort": reasoning_effort},
                text={
                    "format": {
                        "type": "json_schema",
                        "name": "ober_clinician_readable_rationale",
                        "strict": True,
                        "schema": OUTPUT_SCHEMA,
                    }
                },
                max_output_tokens=max_output_tokens,
                store=False,
            )
        except Exception as exc:
            last_error = exc
            if attempt >= retries:
                break
            time.sleep(min(2**attempt, 8))
    raise RuntimeError(f"OpenAI request failed: {last_error}") from last_error


def write_summary_csv(path: Path, results: list[dict[str, Any]], model: str) -> None:
    fieldnames = (
        "醫院",
        "病人編號",
        "檢驗編號",
        "選擇菌種",
        "理由產生器短理由",
        "關鍵證據摘要",
        "重要限制",
        "證據ID",
        "原始位置",
        "需醫師確認",
        "建議確認問題",
        "模型",
    )
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            generated = result["generated_rationale"]
            explanations = generated["selected_pathogen_explanations"]
            if not explanations:
                writer.writerow(
                    {
                        "醫院": result["hospital"],
                        "病人編號": result["patient_id"],
                        "檢驗編號": "; ".join(result["test_identifiers"]) or "未提供",
                        "選擇菌種": "無入選病原",
                        "理由產生器短理由": generated["patient_result_summary"],
                        "關鍵證據摘要": "",
                        "重要限制": "",
                        "證據ID": "; ".join(generated["patient_level_evidence_ids"]),
                        "原始位置": "; ".join(
                            result["evidence_map"].get(item, "")
                            for item in generated["patient_level_evidence_ids"]
                        ),
                        "需醫師確認": "是",
                        "建議確認問題": "請確認是否確實沒有足以列為感染源的病原證據。",
                        "模型": model,
                    }
                )
                continue
            for item in explanations:
                limitation_ids = [
                    evidence_id
                    for limitation in item["limitations"]
                    for evidence_id in limitation["evidence_ids"]
                ]
                all_evidence_ids = list(
                    dict.fromkeys(list(item["evidence_ids"]) + limitation_ids)
                )
                writer.writerow(
                    {
                        "醫院": result["hospital"],
                        "病人編號": result["patient_id"],
                        "檢驗編號": "; ".join(result["test_identifiers"]) or "未提供",
                        "選擇菌種": item["organism"],
                        "理由產生器短理由": item["concise_reason"],
                        "關鍵證據摘要": item["key_evidence_summary"],
                        "重要限制": "；".join(
                            limitation["text"] for limitation in item["limitations"]
                        ),
                        "證據ID": "; ".join(all_evidence_ids),
                        "原始位置": "; ".join(
                            result["evidence_map"].get(evidence_id, "")
                            for evidence_id in all_evidence_ids
                        ),
                        "需醫師確認": "是" if item["requires_clinician_review"] else "否",
                        "建議確認問題": item["clinician_review_question"],
                        "模型": model,
                    }
                )


def main() -> int:
    args = parse_args()
    input_root = args.input_root.resolve()
    output_root = args.output_root.resolve()
    if not input_root.is_dir():
        raise SystemExit(f"input root is not a directory: {input_root}")
    if output_root.exists() and any(output_root.iterdir()):
        raise SystemExit(f"output root must be new or empty: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    if args.env_file:
        load_env_file(args.env_file.resolve())
    if not os.getenv("OPENAI_API_KEY"):
        raise SystemExit("OPENAI_API_KEY is not set; pass --env-file")
    try:
        from openai import OpenAI
    except ImportError as exc:
        raise SystemExit("Install optional dependencies from requirements-llm.txt") from exc

    prompt_path = args.prompt.resolve()
    prompt = prompt_path.read_text(encoding="utf-8-sig").strip()
    prompt_hash = sha256_text(prompt)
    client = OpenAI()
    results: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    for patient_id in requested_patient_ids(args.patients):
        response = None
        try:
            source_file = patient_file(input_root, patient_id)
            chain = load_json(source_file)
            payload, evidence_map = build_payload(chain, args.max_excluded_candidates)
            response = call_model(
                client=client,
                prompt=prompt,
                payload=payload,
                model=args.model,
                reasoning_effort=args.reasoning_effort,
                max_output_tokens=args.max_output_tokens,
                retries=args.retries,
            )
            generated = json.loads(response.output_text)
            validate_output(generated, payload, evidence_map)
            result = {
                "schema_version": "ober_llm_rationale.v1",
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "patient_id": patient_id,
                "hospital": payload.get("hospital") or "未提供",
                "test_identifiers": [
                    compact(item)
                    for item in (payload.get("test_identifiers") or [])
                    if compact(item) and compact(item) != "未提供"
                ],
                "model_requested": args.model,
                "model_returned": getattr(response, "model", None),
                "reasoning_effort": args.reasoning_effort,
                "response_id": getattr(response, "id", None),
                "usage": usage_dict(response),
                "prompt_path": str(prompt_path),
                "prompt_sha256": prompt_hash,
                "input_source": str(source_file.resolve()),
                "input_sha256": sha256_text(json.dumps(payload, ensure_ascii=False, sort_keys=True)),
                "evidence_map": evidence_map,
                "generated_rationale": generated,
            }
            output_file = output_root / "patient_rationales" / f"NGS_patient_{patient_id}_llm_rationale.json"
            write_json(output_file, result, args.pretty)
            result["output_file"] = str(output_file)
            results.append(result)
            print(f"OK patient={patient_id} selected={len(generated['selected_pathogen_explanations'])}", flush=True)
        except Exception as exc:
            if response is not None:
                write_json(
                    output_root / "failed_responses" / f"NGS_patient_{patient_id}_failed_response.json",
                    {
                        "patient_id": patient_id,
                        "validation_error": str(exc),
                        "response_id": getattr(response, "id", None),
                        "model_returned": getattr(response, "model", None),
                        "usage": usage_dict(response),
                        "output_text": getattr(response, "output_text", None),
                        "prompt_sha256": prompt_hash,
                        "not_accepted_as_result": True,
                    },
                    True,
                )
            failures.append({"patient_id": patient_id, "error": str(exc)})
            print(f"FAILED patient={patient_id}: {exc}", flush=True)

    write_summary_csv(output_root / "llm_reasoned_patient_results.csv", results, args.model)
    manifest = {
        "schema_version": "ober_llm_rationale_manifest.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_root": str(input_root),
        "output_root": str(output_root),
        "model": args.model,
        "reasoning_effort": args.reasoning_effort,
        "prompt_path": str(prompt_path),
        "prompt_sha256": prompt_hash,
        "summary": {
            "requested_patients": len(requested_patient_ids(args.patients)),
            "successful_patients": len(results),
            "failed_patients": len(failures),
            "generated_selected_pathogen_reasons": sum(
                len(item["generated_rationale"]["selected_pathogen_explanations"])
                for item in results
            ),
        },
        "patients": [
            {
                "patient_id": item["patient_id"],
                "output_file": item["output_file"],
                "usage": item["usage"],
            }
            for item in results
        ],
        "failures": failures,
    }
    write_json(output_root / "llm_rationale_manifest.json", manifest, True)
    print(
        f"completed={len(results)} failed={len(failures)} output={output_root}",
        flush=True,
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

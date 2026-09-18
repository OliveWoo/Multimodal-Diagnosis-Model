"""Assemble clinical modules into baseline-preserving experiment arms.

This module deliberately does not import the evaluator or load gold labels.
Generation and evaluation are separate processes so predictions cannot depend on
the answer file.
"""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .clinical_rules import evaluate_modules
from .io_adapter import merge_candidate_context, sha256_file


PRIMARY_ARM = "CL2_CONVERGENT"
ARM_DEFINITIONS = {
    "CL0_BASELINE": "E0 (existing picked pathogens; immutable)",
    "CL1_DIRECT": "E0 OR eligible non-E0 with E3 direct-high evidence",
    "CL2_CONVERGENT": "CL1 OR eligible non-E0 with E2 convergent-high evidence",
    "CL3_BALANCED": "CL2 OR D-pass non-E0 with E1 literature-plus-strict-mNGS evidence",
    "CL4_TIERED_EXPLORATORY": (
        "CL2 OR D-pass non-E0 with high-priority+B-strict, or context+A+B-strict"
    ),
}


def _canonical_json(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _hash_object(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value)).hexdigest()


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8-sig") as handle:
        value = json.load(handle)
    if value.get("schema_version") != "rag_re_clinical.config.v1":
        raise ValueError("Expected rag_re_clinical.config.v1 config")
    decision = value.get("decision") or {}
    if decision.get("baseline_immutable") is not True:
        raise ValueError("Clinical replay requires decision.baseline_immutable=true")
    if decision.get("primary_arm") not in ARM_DEFINITIONS:
        raise ValueError("Unknown decision.primary_arm")
    return value


def _evidence_hint(modules: dict[str, Any], frozen: dict[str, Any]) -> bool:
    a = ((frozen.get("modules") or {}).get("A") or {}).get("positive")
    b = (modules.get("B_STRICT") or modules.get("B_strict") or {}).get("positive")
    c_grade = (modules.get("C") or {}).get("grade")
    d_state = (modules.get("D") or {}).get("state")
    return bool(
        a is True
        or b is True
        or c_grade not in {None, "C0", "CNEG", "C0_NONE_OBSERVED", "CNEG_PERFORMED_NEGATIVE"}
        or d_state in {"D_GUARDED", "D_UNKNOWN"}
    )


def _a_strong(frozen: dict[str, Any], config: dict[str, Any]) -> bool | None:
    module = ((frozen.get("modules") or {}).get("A") or {})
    positive = module.get("positive")
    if positive is not True:
        return None if positive is None else False
    rule = config.get("module_a") or {}
    support = module.get("support_count")
    judgeable = module.get("judgeable_count")
    if not isinstance(support, int) or not isinstance(judgeable, int):
        return None
    return support >= int(rule.get("strong_min_support", 3)) and judgeable >= int(
        rule.get("strong_min_judgeable", 5)
    )


def _review_tier(frozen: dict[str, Any], merged: dict[str, Any]) -> str | None:
    tier = frozen.get("review_tier")
    if tier:
        return str(tier)
    review = merged.get("review_metadata") or {}
    return review.get("tier") or review.get("review_tier")


def _arm_values(
    *,
    baseline: bool,
    modules: dict[str, Any],
    frozen: dict[str, Any],
    merged: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, bool]:
    if baseline:
        return {arm: True for arm in ARM_DEFINITIONS}

    d_state = (modules.get("D") or {}).get("state")
    e_state = (modules.get("E") or {}).get("tier") or (modules.get("E") or {}).get("state")
    b_strict = (
        modules.get("B_STRICT") or modules.get("B_strict") or {}
    ).get("positive") is True
    a = _a_strong(frozen, config)
    tier = (_review_tier(frozen, merged) or "").lower()
    d_pass = d_state == "D_PASS"

    cl1 = e_state == "E3_DIRECT_HIGH"
    cl2 = cl1 or e_state == "E2_CONVERGENT_HIGH"
    cl3 = cl2 or (
        d_pass and e_state in {"E1_LITERATURE_MNGS", "E1_BIOLOGIC_PLUS_SIGNAL"}
    )
    tiered = d_pass and b_strict and (
        tier == "review_high_priority"
        or (tier == "review_context_needed" and a is True)
    )
    cl4 = cl2 or tiered
    return {
        "CL0_BASELINE": False,
        "CL1_DIRECT": cl1,
        "CL2_CONVERGENT": cl2,
        "CL3_BALANCED": cl3,
        "CL4_TIERED_EXPLORATORY": cl4,
    }


def _disposition(
    baseline: bool,
    arms: dict[str, bool],
    modules: dict[str, Any],
    frozen: dict[str, Any],
    *,
    primary_arm: str,
) -> str:
    if baseline:
        return "BASELINE_RETAINED"
    if arms.get(primary_arm):
        return "AUTO_RESCUE"
    if _evidence_hint(modules, frozen):
        return "MANUAL_REVIEW"
    return "NOT_RESCUED"


def _validate_patient_output(result: dict[str, Any]) -> None:
    candidates = result.get("candidates") or []
    names = [row.get("organism_name") for row in candidates]
    if not names or any(not isinstance(name, str) or not name.strip() for name in names):
        raise ValueError("Every candidate must have a non-empty organism_name")
    if len(names) != len(set(names)):
        raise ValueError("Duplicate organism_name in clinical output")
    baseline_names = {
        row["organism_name"] for row in candidates if row.get("baseline_selected") is True
    }
    for arm_name in ARM_DEFINITIONS:
        arm = (result.get("arms") or {}).get(arm_name)
        if not isinstance(arm, dict):
            raise ValueError(f"Missing arm {arm_name}")
        predicted = set(arm.get("predicted_pathogens") or [])
        if not predicted.issubset(set(names)):
            raise ValueError(f"{arm_name} predicts a pathogen outside the frozen pool")
        if not baseline_names.issubset(predicted):
            raise ValueError(f"{arm_name} violates E0 immutability")
    if set((result.get("arms") or {})["CL0_BASELINE"]["predicted_pathogens"]) != baseline_names:
        raise ValueError("CL0_BASELINE differs from frozen E0")


def run_patient(
    frozen_artifact: dict[str, Any],
    merge_payload: dict[str, Any],
    config: dict[str, Any],
    *,
    frozen_path: str | Path | None = None,
    merge_path: str | Path | None = None,
) -> dict[str, Any]:
    """Run a single patient without loading or receiving gold labels."""

    if frozen_artifact.get("run_status") != "complete":
        raise ValueError("Frozen RAG_re artifact must have run_status=complete")
    if frozen_artifact.get("schema_version") != "rag_re.output.v1":
        raise ValueError("Expected rag_re.output.v1 frozen artifact")

    envelope = merge_candidate_context(
        frozen_artifact,
        merge_payload,
        frozen_path=frozen_path,
        merge_path=merge_path,
    )
    rows: list[dict[str, Any]] = []
    for merged in envelope["candidates"]:
        frozen = deepcopy(merged["frozen_candidate"])
        baseline = frozen.get("baseline_selected")
        if not isinstance(baseline, bool):
            raise ValueError("baseline_selected must be explicitly boolean")
        rule_context = deepcopy(merged)
        # Keep label-blind patient context (host/time/specimen metadata) visible
        # to F and future clinical modules.  It is deliberately separate from
        # the candidate layers and never includes the gold sidecar.
        rule_context["patient_context"] = deepcopy(envelope.get("merge_context") or {})
        modules = evaluate_modules(rule_context, frozen, config)
        arm_values = _arm_values(
            baseline=baseline,
            modules=modules,
            frozen=frozen,
            merged=merged,
            config=config,
        )
        disposition = _disposition(
            baseline,
            arm_values,
            modules,
            frozen,
            primary_arm=config["decision"]["primary_arm"],
        )
        rows.append(
            {
                "organism_name": frozen.get("organism_name"),
                "canonical_organism_name": frozen.get("canonical_organism_name")
                or frozen.get("organism_name"),
                "baseline_selected": baseline,
                "review_tier": _review_tier(frozen, merged),
                "candidate_provenance": frozen.get("candidate_provenance"),
                "match_scope": merged.get("match_scope"),
                "disposition": disposition,
                "modules": modules,
                "arms": arm_values,
                "audit": {
                    "source_layers": merged.get("provenance"),
                    "critical_unknown": disposition == "MANUAL_REVIEW"
                    and any(
                        (modules.get(key) or {}).get("status") in {
                            "unknown",
                            "insufficient_input",
                            "partial",
                            "discordant_context",
                        }
                        for key in ("D", "B_STRICT", "C", "E", "F")
                    ),
                },
            }
        )

    arm_outputs: dict[str, Any] = {}
    for arm_name, expression in ARM_DEFINITIONS.items():
        predicted = [row["organism_name"] for row in rows if row["arms"][arm_name]]
        manual = [
            row["organism_name"]
            for row in rows
            if not row["arms"][arm_name] and row["disposition"] == "MANUAL_REVIEW"
        ]
        arm_outputs[arm_name] = {
            "expression": expression,
            "predicted_pathogens": predicted,
            "predicted_count": len(predicted),
            "manual_review_pathogens": manual,
            "manual_review_count": len(manual),
        }

    result = {
        "schema_version": "rag_re_clinical.output.v1",
        "pipeline_version": "0.1.0",
        "run_status": "complete",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "patient_id": str(envelope["patient_id"]),
        "study_status": "exploratory_development_only",
        "primary_arm": config["decision"]["primary_arm"],
        "methodology": {
            "baseline_policy": "E0 is immutable in every arm",
            "unknown_policy": "unknown is routed to manual review, never coerced to negative",
            "literature_policy": "reuse frozen A only; no new API or LLM call",
            "gold_isolation": "generation code does not import or load evaluation/gold",
        },
        "provenance": {
            "frozen_rag_re_sha256": (
                (envelope.get("sources") or {}).get("frozen_rag_re") or {}
            ).get("sha256"),
            "merge_input_sha256": (
                (envelope.get("sources") or {}).get("formal_merge") or {}
            ).get("sha256"),
            "frozen_rag_re_pipeline_fingerprint": frozen_artifact.get("pipeline_fingerprint"),
            "frozen_candidate_pool_sha256": (frozen_artifact.get("input_meta") or {}).get(
                "candidate_pool_sha256"
            ),
            "config_sha256": _hash_object(config),
        },
        "candidates": rows,
        "arms": arm_outputs,
    }
    _validate_patient_output(result)
    result["artifact_sha256_basis"] = _hash_object(result)
    return result


def write_json(path: str | Path, value: dict[str, Any]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

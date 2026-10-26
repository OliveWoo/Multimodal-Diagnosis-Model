"""Source audit and annotation helpers, without private cohort configuration."""
from __future__ import annotations
import copy
import json
import re
from pathlib import Path
from tools.pathogen_normalization import canonical_key
from tools.bind_viral_evidence_shadow import read, sha
EXTRA_FIELD = "viral_evidence_shadow"
RUN_FIELD = "viral_evidence_shadow_run"



def pointer_value(payload, pointer: str):
    value = payload
    for name, index in re.findall(r"([A-Za-z_][A-Za-z_0-9]*)|\[(\d+)\]", pointer):
        value = value[int(index)] if index else value[name]
    return value


def audit_binding(binding_root: Path) -> dict:
    summary = read(binding_root / "summary.json")
    if summary.get("candidate_loss") != 0 or not summary.get("source_mapping_completed"):
        raise ValueError("Binding is incomplete or lost candidates")
    if summary.get("production_decisions_changed") is not False or summary.get("benchmark_answers_loaded") is not False:
        raise ValueError("Expected answer-blind annotation-only binding")
    files = {}
    errors = []
    count = 0

    def verify(ref):
        nonlocal count
        count += 1
        path = Path(ref["file"])
        if str(path) not in files:
            files[str(path)] = (sha(path), read(path))
        actual, value = files[str(path)]
        if actual != ref["sha256"]:
            errors.append({"type": "hash_mismatch", "file": str(path)})
        try:
            pointer_value(value, ref["json_path"])
        except (KeyError, IndexError, TypeError) as exc:
            errors.append({"type": "invalid_pointer", "file": str(path), "pointer": ref["json_path"], "error": str(exc)})

    facts = read(binding_root / "source_facts.json")
    bound = read(binding_root / "bound_candidate_packets.json")
    frozen = read(binding_root / "frozen_affected_evidence_packets.json")
    for fact in facts:
        for ref in fact["source_refs"]:
            verify(ref)
    for item in bound:
        for claim in item["claims"]:
            verify({"file": claim["source_file"], "sha256": claim["source_sha256"], "json_path": claim["source_row"]})
            if str(claim["patient_id"]) != str(item["patient_id"]):
                errors.append({"type": "patient_mismatch", "claim": claim["evidence_id"]})
            if claim["episode_id"] != item["episode_id"] or claim["organism_key"] != item["organism_key"]:
                errors.append({"type": "episode_or_organism_mismatch", "claim": claim["evidence_id"]})
        reconciliation = item["frozen_hospital_observation_reconciliation"]
        if reconciliation["missing"]:
            errors.append({"type": "lost_hospital_observation", "patient": item["patient_id"]})
    before = {(p["cohort"], str(p["patient_id"]), canonical_key(p["organism_name"])) for p in frozen}
    after = {(p["cohort"], str(p["patient_id"]), p["organism_key"]) for p in bound}
    if before != after or len(after) != len(bound):
        errors.append({"type": "candidate_set_mismatch"})
    for key, path in (
        ("binding_policy_sha256", binding_root / "viral_evidence_binding_shadow_v1.json"),
        ("axes_policy_sha256", binding_root / "viral_evidence_axes_shadow_v1.json"),
        ("evidence_packets_sha256", binding_root / "frozen_affected_evidence_packets.json"),
        ("adapter_code_sha256", binding_root / "code_snapshot/bind_viral_evidence_shadow.py"),
        ("axes_code_sha256", binding_root / "code_snapshot/viral_evidence_axes_shadow.py"),
    ):
        if not path.is_file() or sha(path) != summary[key]:
            errors.append({"type": "frozen_component_mismatch", "component": key})
    audit = {"passed": not errors, "checked_source_pointers": count,
             "checked_source_files": len(files), "candidate_count": len(bound), "errors": errors}
    if errors:
        raise ValueError(json.dumps(audit, ensure_ascii=False))
    return audit


def without_annotations(payload):
    if isinstance(payload, dict):
        return {k: without_annotations(v) for k, v in payload.items() if k not in {EXTRA_FIELD, RUN_FIELD}}
    if isinstance(payload, list):
        return [without_annotations(v) for v in payload]
    return payload


def evidence_annotation(bound: dict, index: int, binding_file: Path) -> dict:
    axes = bound["axes"]
    return {"mode": "annotation_only_shadow", "source_mapping_status": axes["source_mapping_status"],
            "detection": axes["detection"], "activity": axes["activity"], "organ_support": axes["organ_support"],
            "detection_description_zh": "另有可核對的檢驗／檢體支持檢出" if axes["detection"] == "D2" else "具有可判讀的檢出訊號",
            "activity_description_zh": "目前來源尚未建立活動／再活化支持；未知不等於陰性" if axes["activity"] == "A0" else "請依所附來源解讀",
            "pulmonary_description_zh": "有同次肺部背景，病毒歸因仍未知" if axes["organ_support"] == "O1" else "目前未建立相符肺部疾病支持",
            "evidence_gap_codes": axes["gaps"], "episode_id": bound["episode_id"],
            "claim_count": len(bound["claims"]), "unique_supporting_evidence_count": axes["unique_supporting_evidence_count"],
            "frozen_hospital_observation_reconciliation": bound["frozen_hospital_observation_reconciliation"],
            "evidence_consumption_ledger": axes["evidence_consumption_ledger"],
            "binding_source": {"file": str(binding_file.resolve()), "sha256": sha(binding_file), "json_path": f"[{index}]"},
            "production_tier_unchanged": True, "clinical_adjudication_complete": False,
            "decision_effect": "none; shadow organ-attribution policy has not replaced clinical decisions"}


def attach_annotations(payload: dict, bindings: dict, *, binding_file: Path) -> tuple[dict, int]:
    output = copy.deepcopy(payload)
    patient = str(payload["patient_id"])
    added = 0
    # Both compact results and full reporting packets have multiple views of
    # the same candidate. Every view receives exactly the same annotation.
    for values in output.values():
        if not isinstance(values, list):
            continue
        for candidate in values:
            if not isinstance(candidate, dict) or not candidate.get("organism_name"):
                continue
            key = (patient, canonical_key(candidate["organism_name"]))
            if key not in bindings:
                continue
            index, bound = bindings[key]
            if str(candidate.get("patient_id")) != patient:
                raise ValueError("Candidate/patient identity mismatch")
            tier = (candidate.get("reporting") or {}).get("tier", candidate.get("final_reporting_tier"))
            selected = (candidate.get("reporting") or {}).get("selected_for_complete_report", candidate.get("selected_for_complete_report"))
            if tier != bound["current_tier"] or bool(selected) != bound["current_selected"]:
                raise ValueError(f"Bound evidence does not match frozen decision: {key}")
            candidate[EXTRA_FIELD] = evidence_annotation(bound, index, binding_file)
            added += 1
    output[RUN_FIELD] = {"mode": "annotation_only_shadow", "decision_effect": "none",
        "binding_file": str(binding_file.resolve()), "binding_sha256": sha(binding_file)}
    if without_annotations(output) != payload:
        raise ValueError("Annotation changed an original decision or evidence field")
    return output, added

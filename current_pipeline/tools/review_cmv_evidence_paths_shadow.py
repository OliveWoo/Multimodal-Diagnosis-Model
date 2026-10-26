"""Frozen, answer-blind CMV reporting comparison; never edits production output.

The reporting policy distinguishes a low-confidence etiologic candidate from a
clinical probable/proven diagnosis. The existing three biological axes are not
upgraded simply because reporting selects Possible.
"""
from __future__ import annotations

import argparse
import copy
import json
import re
from collections import Counter
from pathlib import Path

from tools.bind_viral_evidence_shadow import read, sha, write, write_csv
from tools.project_viral_evidence_patient_outputs import audit_binding, pointer_value
from tools.viral_evidence_axes_shadow import interpret

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINDING = ROOT / "outputs/runs/2026-10-19_viral_evidence_binding_shadow_v4"
DEFAULT_POLICY = ROOT / "rules/cmv_evidence_reporting_shadow_v1.json"
DEFAULT_CLARIFICATIONS = ROOT / "rules/clinical_source_clarifications_template.json"
DEFAULT_OUTPUT = ROOT / "outputs/runs/2026-10-20_cmv_evidence_reporting_shadow_v1"


def apply_clarifications(bound: dict, clarifications: list[dict]) -> tuple[dict, list[str]]:
    """Add user provenance without rewriting the original clinical source."""
    output = copy.deepcopy(bound)
    applied = []
    for note in clarifications:
        if (note["cohort"], str(note["patient_id"]), note["organism_key"], note["episode_id"]) != (
            bound["cohort"], str(bound["patient_id"]), bound["organism_key"], bound["episode_id"]
        ):
            continue
        matches = [c for c in output["claims"] if c["evidence_id"] == note["evidence_id"]]
        if len(matches) != 1:
            raise ValueError("Clarification must resolve exactly one accepted source claim")
        claim = matches[0]
        if (claim.get("source_sha256"), claim.get("source_row")) != (
            note["original_source_sha256"], note["original_source_row"]
        ) or claim.get("claim_type") != "detection" or claim.get("assay_method") == "mNGS":
            raise ValueError("Clarification source identity changed or is not hospital detection")
        if note["assertion"] != "assay_independent_of_mngs":
            raise ValueError("Unsupported source clarification")
        claim["independence_class_before_user_clarification"] = claim.get("independence_class")
        claim["independence_class"] = "independent_clinical_assay"
        claim["independence_basis"] = "user_confirmed_different_method_not_source_document_verified"
        claim["source_metadata_clarification"] = copy.deepcopy(note)
        applied.append(note["clarification_id"])
    return output, applied


def _specimen_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def decide(bound: dict, candidate: dict, policy: dict, axes_policy: dict) -> dict:
    """No reference answers, patient exceptions, rank or test-count decisions."""
    if bound["organism_key"] not in policy["organism_keys"]:
        raise ValueError("CMV-only shadow: other viruses require their own reviewed paths")
    claims = [c for c in bound["claims"] if str(c.get("patient_id")) == str(bound["patient_id"])]
    axes = interpret(candidate, episode_id=bound["episode_id"], claims=claims, policy=axes_policy)
    rejected = {c["evidence_id"] for c in axes["rejected_claims"]}
    accepted = [c for c in claims if c["evidence_id"] not in rejected]
    respiratory = {_specimen_key(s) for s in policy["respiratory_specimens"]}
    valid_mngs = bool(
        (candidate.get("taxonomy_profile") or {}).get("mapping_status") == "exact_species"
        and any(float(s.get("reads") or 0) > 0 and s.get("qc_status") in axes_policy["detection_qc_statuses"]
                for s in (candidate.get("analytical_profile") or {}).get("per_test_signals", []))
    )
    mngs_groups = {c.get("biological_evidence_group", c["evidence_id"])
                   for c in accepted if c.get("assay_method") == "mNGS"}
    hospital = [c for c in accepted if c.get("claim_type") == "detection"
                and c.get("assay_method") != "mNGS" and c.get("identity_match") == "exact_or_alias"]
    lower_resp = [c for c in hospital if _specimen_key(c.get("specimen")) in respiratory]
    supports = [c for c in lower_resp
                if c.get("independence_class") in policy["independent_detection_classes"]
                and c.get("biological_evidence_group", c["evidence_id"]) not in mngs_groups]
    gaps = list(axes["gaps"])
    if any(c.get("independence_class") == "not_established" for c in lower_resp):
        gaps.append("respiratory_hospital_assay_independence_unknown")
    if axes["organ_support"] == "O0":
        gaps.append("aligned_pulmonary_syndrome_not_established_in_available_input")
    tier, route = "Context", "detected_significance_unknown"
    if axes["detection"] == "D0":
        tier, route = "Audit", "no_accepted_detection"
    elif axes["activity"] == "A3":
        route = "explicit_shedding_or_resolved_status_requires_review"
        gaps.append("resolved_or_shedding_interpretation_requires_etiologic_review")
    elif axes["organ_support"] in {"O2", "O3"}:
        tier = "Picked" if axes["organ_support"] == "O3" else "Possible"
        route = "source_reviewed_pulmonary_attribution"
    elif valid_mngs and supports:
        tier, route = "Possible", policy["possible_route"]
    return {
        "shadow_reporting_tier": tier,
        "selected_as_pulmonary_candidate": tier in {"Possible", "Picked"},
        "doctor_visible": axes["detection"] != "D0",
        "report_label_zh": policy["possible_label_zh"] if tier == "Possible" and route == policy["possible_route"]
            else policy["unknown_label_zh"] if tier == "Context" else tier,
        "route": route, "axes": axes,
        "reporting_selection_is_clinical_diagnosis": False,
        "clinical_probable_or_proven_label_inferred_from_reporting": False,
        "activity_or_organ_axis_upgraded_by_reporting": False,
        "independent_respiratory_support_ids": sorted({c["evidence_id"] for c in supports}),
        "respiratory_hospital_detection_ids": sorted({c["evidence_id"] for c in lower_resp}),
        "gaps": sorted(set(gaps)),
        "co_pathogen_veto_applied": False,
        "coinfection_assertion": "co-detection/candidate co-infection only; confirmed co-infection not inferred",
        "benchmark_answers_used": False, "production_changed": False,
    }


def counterfactual_independence(bound: dict, candidate: dict, policy: dict, axes_policy: dict) -> list[dict]:
    """Non-applied sensitivity: unknown method is confirmed vs remains unknown."""
    respiratory = {_specimen_key(s) for s in policy["respiratory_specimens"]}
    targets = [c for c in bound["claims"] if c.get("claim_type") == "detection"
               and c.get("assay_method") != "mNGS"
               and _specimen_key(c.get("specimen")) in respiratory
               and c.get("independence_class") == "not_established"]
    out = []
    for target in targets:
        alternate = copy.deepcopy(bound)
        for claim in alternate["claims"]:
            if claim["evidence_id"] == target["evidence_id"]:
                claim["independence_class"] = "independent_clinical_assay"
        result = decide(alternate, candidate, policy, axes_policy)
        out.append({"evidence_id": target["evidence_id"], "assumption": "different_method_confirmed",
                    "counterfactual_only_not_applied": True,
                    "shadow_tier_if_confirmed": result["shadow_reporting_tier"]})
    return out


def run(binding_root: Path, policy_path: Path, clarification_path: Path, output_root: Path) -> dict:
    if output_root.exists():
        raise FileExistsError(f"Refusing to overwrite {output_root}")
    audit = audit_binding(binding_root)
    bound = read(binding_root / "bound_candidate_packets.json")
    frozen = read(binding_root / "frozen_affected_evidence_packets.json")
    axes_policy = read(binding_root / "viral_evidence_axes_shadow_v1.json")
    policy, clarification = read(policy_path), read(clarification_path)
    original = copy.deepcopy(bound)
    result_rows, preserved, applied = [], [], []
    inputs = [binding_root / name for name in ("bound_candidate_packets.json", "frozen_affected_evidence_packets.json",
                                               "viral_evidence_axes_shadow_v1.json", "source_registry.json")]
    inputs += [policy_path, clarification_path, Path(__file__).resolve(), ROOT / "tools/viral_evidence_axes_shadow.py"]
    registry = {str(p.resolve()): sha(p) for p in inputs}
    for item in bound:
        if item["organism_key"] not in policy["organism_keys"]:
            preserved.append(copy.deepcopy(item))
            continue
        ref = item["source_candidate"]
        if sha(Path(ref["file"])) != ref["sha256"]:
            raise ValueError("Frozen candidate source hash mismatch")
        candidate = pointer_value(frozen, ref["json_path"])
        enriched, notes = apply_clarifications(item, clarification["clarifications"])
        applied.extend(notes)
        decision = decide(enriched, candidate, policy, axes_policy)
        row = {"cohort": item["cohort"], "patient_id": item["patient_id"], "organism_name": item["organism_name"],
               "episode_id": item["episode_id"], "before_tier": item["current_tier"],
               "before_selected": item["current_selected"], "selected_event": item["selected_event"],
               "source_candidate": ref, "applied_source_clarifications": notes,
               "decision": decision,
               "method_confirmation_sensitivity_not_applied": counterfactual_independence(enriched, candidate, policy, axes_policy)}
        result_rows.append(row)
        packet = copy.deepcopy(item)
        packet["cmv_reporting_shadow"] = copy.deepcopy(row)
        # Preserve the original claims, axes and current tier; enriched claims are explicitly separate.
        packet["cmv_reporting_shadow"]["claims_with_source_clarifications"] = enriched["claims"]
        preserved.append(packet)
    expected = sorted(n["clarification_id"] for n in clarification["clarifications"])
    if sorted(applied) != expected:
        raise ValueError("Not all source clarifications resolved exactly once")
    if bound != original or len(preserved) != len(original):
        raise AssertionError("Candidate loss or source mutation")
    if any({k: v for k, v in after.items() if k != "cmv_reporting_shadow"} != before
           for before, after in zip(original, preserved)):
        raise AssertionError("Original candidate field changed")
    if any(sha(Path(p)) != digest for p, digest in registry.items()):
        raise ValueError("Input changed during run")
    transitions = Counter(f"{r['before_tier']} -> {r['decision']['shadow_reporting_tier']}" for r in result_rows)
    summary = {"policy_id": policy["policy_id"], "mode": "cmv_reporting_shadow_only",
        "benchmark_answers_loaded": False, "development_has_seen_reference_answers": True,
        "independent_external_validation": False, "production_decisions_changed": False,
        "original_candidate_fields_identical": True, "candidate_count_preserved": len(preserved),
        "cmv_candidates_reviewed": len(result_rows), "other_viral_candidates_unchanged": len(preserved) - len(result_rows),
        "source_clarifications_applied": applied,
        "before_selected": sum(bool(r["before_selected"]) for r in result_rows),
        "shadow_selected": sum(r["decision"]["selected_as_pulmonary_candidate"] for r in result_rows),
        "all_detected_doctor_visible": all(r["decision"]["doctor_visible"] for r in result_rows),
        "transitions": dict(transitions),
        "shadow_activity_counts": dict(Counter(r["decision"]["axes"]["activity"] for r in result_rows)),
        "shadow_organ_counts": dict(Counter(r["decision"]["axes"]["organ_support"] for r in result_rows)),
        "scope": "CMV only; no automatic changes to HSV/EBV/PJP/other pathogen gates",
        "counterfactuals_not_applied": True}
    output_root.mkdir(parents=True)
    write(output_root / "frozen_registry.json", {"answer_access": False, "inputs": registry})
    for name, source in (("policy.json", policy_path), ("source_clarifications.json", clarification_path),
                         ("axes_policy.json", binding_root / "viral_evidence_axes_shadow_v1.json")):
        write(output_root / name, read(source))
    snapshot = output_root / "code_snapshot" / Path(__file__).name
    snapshot.parent.mkdir()
    snapshot.write_bytes(Path(__file__).read_bytes())
    axes_code = ROOT / "tools/viral_evidence_axes_shadow.py"
    (snapshot.parent / axes_code.name).write_bytes(axes_code.read_bytes())
    write(output_root / "source_audit.json", audit)
    write(output_root / "all_viral_candidates_preserved.json", preserved)
    write(output_root / "all_cmv_decisions.json", result_rows)
    write(output_root / "summary.json", summary)
    write_csv(output_root / "all_cmv_before_after.csv", [
        {"cohort": r["cohort"], "patient_id": r["patient_id"], "before": r["before_tier"],
         "after_shadow": r["decision"]["shadow_reporting_tier"], "doctor_visible": r["decision"]["doctor_visible"],
         "activity": r["decision"]["axes"]["activity"], "organ_support": r["decision"]["axes"]["organ_support"],
         "route": r["decision"]["route"], "clarification_ids": ";".join(r["applied_source_clarifications"]),
         "gaps": ";".join(r["decision"]["gaps"])} for r in result_rows])
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binding-root", type=Path, default=DEFAULT_BINDING)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--clarifications", type=Path, default=DEFAULT_CLARIFICATIONS)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run(args.binding_root, args.policy, args.clarifications, args.output_root), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

"""One answer-free entry: deterministic replay -> source-bound CMV -> reports.

Uses explicit private manifests and a frozen baseline, then applies the approved
CMV reporting policy and the existing patient-level fallback. No cohort roster,
source clarification or answer file is embedded in this public entry.
"""
from __future__ import annotations

import argparse
import copy
import json
from collections import Counter
from pathlib import Path

from tools import current_pipeline_defaults as defaults
from tools.public_integrated_config import load_config, cohort_for, cohort_runs
from tools import bind_viral_evidence_shadow as binding
from tools import run_viral_policy_consistency_review as replay
from tools import build_test_aware_possible_pathogen_shadow as reporting
from tools.pathogen_normalization import canonical_key
from tools.project_viral_evidence_patient_outputs import audit_binding
from tools.review_cmv_evidence_paths_shadow import apply_clarifications, decide
from tools.run_kh_integrated_decision_pipeline import build_integrated_results

ROOT = defaults.ROOT
DEFAULT_OUTPUT = ROOT / "outputs/runs/2026-10-20_current_integrated_cmv_v1"
REFERENCE = ROOT / "outputs/runs/2026-10-19_viral_policy_consistency_review_v5"
VIRAL_FAMILIES = {"herpes_or_reactivation_virus", "polyomavirus_persistence_or_shedding"}
POLICIES = {"cmv_reporting": defaults.CMV_POLICY, "fallback": defaults.FALLBACK_POLICY,
            "source_clarifications": defaults.SOURCE_CLARIFICATIONS,
            "binding": binding.POLICY_PATH, "axes": binding.AXES_POLICY}


ACTIVE_CONFIG = None

def remap_clarifications(items: list[dict], notes: list[dict]) -> list[dict]:
    """A snapshot path may change IDs, but never the verified row/content/event."""
    remapped = []
    for note in notes:
        matches = []
        for item in items:
            if (item["cohort"], str(item["patient_id"]), item["organism_key"], item["episode_id"]) != (
                note["cohort"], str(note["patient_id"]), note["organism_key"], note["episode_id"]
            ):
                continue
            for claim in item["claims"]:
                if (claim.get("source_sha256"), claim.get("source_row")) == (
                    note["original_source_sha256"], note["original_source_row"]
                ) and claim.get("claim_type") == "detection" and claim.get("assay_method") != "mNGS":
                    matches.append(claim)
        if len(matches) != 1:
            raise ValueError(f"Source clarification must resolve exactly one unchanged row: {note['clarification_id']}")
        resolved = copy.deepcopy(note)
        resolved["prior_evidence_id"] = note["evidence_id"]
        resolved["evidence_id"] = matches[0]["evidence_id"]
        resolved["resolution_basis"] = "same_cohort_patient_organism_episode_content_hash_and_json_row"
        remapped.append(resolved)
    return remapped


def apply_patient_policy(packet: dict, bounds: dict, clarifications: list[dict],
                         cmv_policy: dict, axes_policy: dict, base_reporting_policy: dict) -> tuple[dict, list[dict]]:
    """Final CMV policy overrides old CMV reporting gates, never clinical axes."""
    output = copy.deepcopy(packet)
    candidates = output["all_forwarded_candidates"]
    before = {canonical_key(c["organism_name"]): copy.deepcopy(c) for c in candidates}
    if len(before) != len(candidates):
        raise ValueError("Duplicate candidate identity in patient")
    for candidate in candidates:
        key = canonical_key(candidate["organism_name"])
        bound = bounds.get(key)
        if bound:
            enriched, applied = apply_clarifications(bound, clarifications)
            candidate["viral_evidence"] = {"binding": enriched, "source_clarifications_applied": applied,
                                           "decision_effect": "annotation_only_for_non_cmv"}
        # First clear prior fallback so it is selected once at the final boundary.
        if candidate.get("final_reporting_tier") == "Fallback-Possible":
            candidate["final_reporting_tier"] = "Context"
            candidate["selected_for_complete_report"] = False
            candidate["possible_reporting_role"] = reporting.NOT_SELECTED_ROLE
            candidate["possible_pathogen_gate"].update(selected_for_complete_report=False,
                reporting_role=reporting.NOT_SELECTED_ROLE, outcome="fallback_pending_final_patient_recalculation")
        if key not in cmv_policy["organism_keys"]:
            continue
        if not bound:
            raise ValueError("CMV candidate has no source-bound evidence")
        decision = decide(enriched, candidate, cmv_policy, axes_policy)
        candidate["viral_evidence"].update(decision_effect="cmv_reporting_policy_applied", cmv_decision=decision)
        tier = decision["shadow_reporting_tier"]
        role = reporting.STRICT_ROLE if tier == "Picked" else (
            "cmv_independent_respiratory_detection_possible" if tier == "Possible" else reporting.NOT_SELECTED_ROLE)
        old_gate = copy.deepcopy(candidate["possible_pathogen_gate"])
        gate = copy.deepcopy(old_gate)
        gate.update(outcome="cmv_source_bound_reporting", route=decision["route"], reporting_role=role,
                    selected_for_complete_report=decision["selected_as_pulmonary_candidate"],
                    blockers=[] if decision["selected_as_pulmonary_candidate"] else decision["gaps"],
                    cautions=[decision["report_label_zh"], "Reporting selection does not establish CMV activity or pneumonia."])
        gate["prior_v7_gate_for_audit_only"] = old_gate
        gate["axes"]["viral_evidence_axes"] = copy.deepcopy(decision["axes"])
        candidate.update(possible_pathogen_gate=gate, possible_reporting_role=role,
                         selected_for_complete_report=decision["selected_as_pulmonary_candidate"],
                         final_reporting_tier=tier, final_reporting_confidence="low" if tier == "Possible" else None)
    profile = output.get("decision_evidence_profile") or reporting.decision_evidence_profile(candidates)
    eligible = [c for c in candidates if not (
        (c.get("viral_evidence") or {}).get("cmv_decision", {}).get("axes", {}).get("activity") == "A3")]
    # Reuse v7's existing whole-patient ranker; do not create a CMV-only fallback.
    reporting.apply_best_available_fallback(eligible, base_reporting_policy, evidence_profile=profile)
    changes = []
    for candidate in candidates:
        key = canonical_key(candidate["organism_name"])
        original = before[key]
        for field in ("clinical_decision", "clinical_level", "formal_pick_allowed", "analytical_profile",
                      "history_rule_ids", "history_evidence_consumption", "host_evidence", "hospital_profile"):
            if candidate.get(field) != original.get(field):
                raise AssertionError(f"CMV reporting changed upstream evidence/clinical field: {field}")
        if "viral_evidence" in candidate:
            candidate["viral_evidence"]["final_reporting_tier"] = candidate["final_reporting_tier"]
            candidate["viral_evidence"]["counts_as_positive_prediction"] = bool(candidate["selected_for_complete_report"])
        if candidate.get("final_reporting_tier") != original.get("final_reporting_tier") or bool(
            candidate.get("selected_for_complete_report")) != bool(original.get("selected_for_complete_report")):
            changes.append({"patient_id": str(packet["patient_id"]), "organism_name": candidate["organism_name"],
                "before_tier": original.get("final_reporting_tier"), "after_tier": candidate.get("final_reporting_tier"),
                "before_selected": bool(original.get("selected_for_complete_report")),
                "after_selected": bool(candidate.get("selected_for_complete_report"))})
    complete = [c for c in candidates if c["selected_for_complete_report"]]
    ordinary = [c for c in complete if c["final_reporting_tier"] in {"Picked", "Possible"}]
    if ordinary and any(c["final_reporting_tier"] == "Fallback-Possible" for c in complete):
        raise AssertionError("Fallback conflicts with ordinary selection")
    output.update(strict_picked=[c for c in candidates if c["final_reporting_tier"] == "Picked" and c["selected_for_complete_report"]],
        possible_pathogens=[c for c in complete if c["final_reporting_tier"] in {"Possible", "Fallback-Possible"}],
        complete_report=complete,
        remaining_high=[c for c in candidates if c["clinical_decision"] == "review_high_priority" and not c["selected_for_complete_report"]],
        remaining_context=[c for c in candidates if c["clinical_decision"] == "review_context_needed" and not c["selected_for_complete_report"]],
        abstention={"abstained": not bool(complete), "reason": "no_credible_candidate" if not complete else "not_applicable"},
        policy_id="current_integrated_cmv_v1_20261020", base_reporting_policy_id=packet["policy_id"],
        cmv_policy_id=cmv_policy["policy_id"], answer_blind=True)
    if len(candidates) != len(packet["all_forwarded_candidates"]):
        raise AssertionError("Candidate loss")
    return output, changes


def collect_affected_packets(base: Path) -> list[dict]:
    packets = []
    for run_name in replay.SOURCES:
        for path in sorted((base / run_name / "10_reporting_v3_route_a/patient_outputs").glob("*.json")):
            payload = binding.read(path)
            for candidate in payload["all_forwarded_candidates"]:
                if (candidate.get("taxonomy_profile") or {}).get("primary_rule_family") not in VIRAL_FAMILIES:
                    continue
                packets.append({"cohort": cohort_for(ACTIVE_CONFIG, run_name, payload["patient_id"]),
                    "patient_id": str(payload["patient_id"]), "organism_name": candidate["organism_name"],
                    "source_packet": str(path.resolve()), "source_packet_sha256": binding.sha(path),
                    "current_candidate_evidence_unchanged": candidate})
    return packets


def write_reporting(output_root: Path, packets: list[dict], policy_id: str) -> dict:
    target = output_root / "reporting"
    (target / "patient_outputs").mkdir(parents=True)
    all_candidates = []
    for packet in packets:
        binding.write(target / "patient_outputs" / f"NGS_patient_{packet['patient_id']}_test_aware_possible_pathogen_shadow.json", packet)
        all_candidates.extend(packet["all_forwarded_candidates"])
    for name, selected in (("possible_pathogen_decisions.csv", all_candidates),
            ("complete_report.csv", [c for c in all_candidates if c["selected_for_complete_report"]]),
            ("possible_pathogens.csv", [c for c in all_candidates if c["selected_for_complete_report"] and c["final_reporting_tier"] != "Picked"])):
        reporting._write_csv(target / name, selected, include_v2=True)
    summary = {"policy_id": policy_id, "answer_blind": True, "patient_count": len(packets),
        "candidate_count": len(all_candidates), "complete_report_count": sum(bool(c["selected_for_complete_report"]) for c in all_candidates),
        "final_tiers": dict(Counter(c["final_reporting_tier"] for c in all_candidates))}
    binding.write(target / "summary.json", summary)
    summary["integrated_output"] = build_integrated_results(target, output_root)
    binding.write(output_root / "summary.json", summary)
    return summary


def run(output_root: Path, *, config_path: Path) -> dict:
    global ACTIVE_CONFIG
    config = load_config(config_path)
    ACTIVE_CONFIG = config
    reference = config["reference_root"]
    replay.SOURCES = {name: name for name in config["source_manifests"]}
    replay.SOURCE_MANIFESTS = config["source_manifests"]
    replay.TAXONOMY_COMPLETION_SOURCE = config["taxonomy_completion_manifest"]
    replay.EPISODE_IDENTITY_SOURCE = config["episode_identity_file"]
    replay.DEVELOPMENT_FREEZE_SOURCE = config["development_freeze"]
    policy_sources = {**POLICIES, "source_clarifications": config["source_clarifications"]}
    if output_root.exists():
        raise FileExistsError(f"Refusing to overwrite {output_root}")
    output_root.mkdir(parents=True)
    frozen_policies = output_root / "00_current_rules"
    frozen_policies.mkdir()
    policy_registry = {}
    for role, path in policy_sources.items():
        target = frozen_policies / path.name
        target.write_bytes(path.read_bytes())
        policy_registry[role] = {"path": str(target.resolve()), "sha256": binding.sha(target)}
    config_snapshot = {k: v for k, v in config.items() if k != "_input_hashes"}
    binding.write(output_root / "run_configuration.json", json.loads(json.dumps(config_snapshot, default=str)))
    code_registry = {str(p.resolve()): binding.sha(p) for p in (ROOT / "tools").glob("*.py")}
    binding.write(output_root / "implementation_freeze.json", {"policies": policy_registry, "code_sha256": code_registry,
        "private_input_sha256": config["_input_hashes"], "answers_loaded": False, "status": "integration_replay_not_independent_external_validation"})
    print("Rebuilding deterministic stages from explicit private manifests", flush=True)
    base = output_root / "01_v7_deterministic_replay"
    consistency = replay.run(base)
    # Compare replay with the previously audited v7 at candidate/selection level.
    replay_comparison = []
    for run_name in replay.SOURCES:
        for path in (base / run_name / "patient_results").glob("*.json"):
            now = binding.read(path)
            old = binding.read(reference / run_name / "patient_results" / path.name)
            def view(payload):
                return {canonical_key(c["organism_name"]): (c["decision_path"]["final_reporting_tier"],
                    bool(c["reporting"]["selected_for_complete_report"])) for c in payload["all_candidates"]}
            if view(now) != view(old):
                replay_comparison.append({"run": run_name, "patient_id": now["patient_id"], "old": view(old), "new": view(now)})
    binding.write(output_root / "v7_replay_comparison.json", {"matches_prior_v7": not replay_comparison, "differences": replay_comparison})
    if replay_comparison:
        raise ValueError("Baseline replay changed; inspect before integrating policy")
    binding.write(base / "review/all_affected_evidence_packets.json", collect_affected_packets(base))
    binding_root = output_root / "02_viral_evidence_binding"
    print("Binding sources to all viral candidates", flush=True)
    binding_summary = binding.run(base, binding_root, cohort_runs=cohort_runs(config))
    source_audit = audit_binding(binding_root)
    bound = binding.read(binding_root / "bound_candidate_packets.json")
    notes = remap_clarifications(bound, binding.read(Path(policy_registry["source_clarifications"]["path"]))["clarifications"])
    binding.write(output_root / "resolved_source_clarifications.json", {"clarifications": notes})
    cmv_policy = binding.read(frozen_policies / defaults.CMV_POLICY.name)
    axes_policy = binding.read(frozen_policies / binding.AXES_POLICY.name)
    # The approved fallback v2 overrides only reporting endpoint semantics;
    # the actual whole-patient fallback ranker remains v7.
    fallback_policy = binding.read(frozen_policies / defaults.FALLBACK_POLICY.name)
    if not fallback_policy["fallback"]["include_in_positive_prediction_endpoint"]:
        raise ValueError("Approved positive endpoint must retain Fallback-Possible")
    bounds = {(b["cohort"], str(b["patient_id"]), b["organism_key"]): b for b in bound}
    grouped = {co: [] for co in cohort_runs(config)}
    all_changes = []
    seen = set()
    for run_name in replay.SOURCES:
        base_policy, _ = reporting.load_policy(base / "00_frozen_rules" / defaults.REPORTING_POLICY.name)
        for path in sorted((base / run_name / "10_reporting_v3_route_a/patient_outputs").glob("*.json")):
            payload = binding.read(path)
            cohort = cohort_for(config, run_name, payload["patient_id"])
            patient_bounds = {key[2]: item for key, item in bounds.items() if key[:2] == (cohort, str(payload["patient_id"]))}
            result, changes = apply_patient_policy(payload, patient_bounds, notes, cmv_policy, axes_policy, base_policy)
            for c in result["all_forwarded_candidates"]:
                if "viral_evidence" in c:
                    seen.add((cohort, str(payload["patient_id"]), canonical_key(c["organism_name"])))
            grouped[cohort].append(result)
            all_changes.extend(dict(c, cohort=cohort) for c in changes)
    if seen != set(bounds):
        raise AssertionError("Not all viral evidence bindings reached final output")
    summaries = {co: write_reporting(output_root / co, packets, "current_integrated_cmv_v1_20261020") for co, packets in grouped.items()}
    binding.write_csv(output_root / "all_selection_changes.csv", all_changes)
    if any(binding.sha(Path(p)) != digest for p, digest in config["_input_hashes"].items()):
        raise ValueError("Private input/config changed during replay")
    if any(binding.sha(Path(p)) != digest for p, digest in code_registry.items()):
        raise ValueError("Code changed during replay")
    if any(binding.sha(Path(item["path"])) != item["sha256"] for item in policy_registry.values()):
        raise ValueError("Frozen policy changed during replay")
    artifacts = {str(p.resolve()): binding.sha(p) for co in grouped
                 for p in (output_root / co).rglob("*") if p.is_file()}
    summary = {"release_id": output_root.name, "status": "integrated_research_release_completed",
        "answers_loaded": False, "ober_included": False, "legacy_pipeline_changed": False,
        "base_reporting_policy": defaults.REPORTING_POLICY.name, "cmv_policy_id": cmv_policy["policy_id"],
        "fallback_policy_id": fallback_policy["policy_id"],
        "reference_v7_preserved": True, "baseline_replay_matches_audited_v7": not replay_comparison,
        "cohorts": summaries, "patient_count": sum(s["patient_count"] for s in summaries.values()),
        "candidate_count": sum(s["candidate_count"] for s in summaries.values()),
        "viral_candidate_count": len(bound), "cmv_candidate_count": sum(b["organism_key"] in cmv_policy["organism_keys"] for b in bound),
        "selection_change_count": len(all_changes), "candidate_loss": 0,
        "shared_base_component_hashes": consistency["shared_component_hashes"],
        "source_audit": source_audit, "binding_summary": binding_summary,
        "evaluation_separation": "Evaluation is a separate postfreeze command; no labels used here.",
        "positive_endpoint": "Picked + Possible + Fallback-Possible",
        "validation_status": "development-cohort integration stress test, not independent external validation"}
    binding.write(output_root / "summary.json", summary)
    binding.write(output_root / "answer_blind_freeze.json", {"passed": True, "answers_loaded": False,
        "patient_count": summary["patient_count"], "candidate_count": summary["candidate_count"],
        "artifacts_sha256": artifacts, "summary_sha256": binding.sha(output_root / "summary.json"),
        "configuration_sha256": binding.sha(output_root / "run_configuration.json")})
    print("Completed configured cohort decisions; answers have not been loaded", flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--config", type=Path, required=True, help="Private run configuration; never commit it")
    args = parser.parse_args()
    summary = run(args.output_root, config_path=args.config)
    print(json.dumps({k: summary[k] for k in ("release_id", "patient_count", "candidate_count", "selection_change_count", "answers_loaded")}, indent=2))


if __name__ == "__main__":
    main()

"""Validated fallback finalization helper; no private cohort replay."""
from __future__ import annotations
import copy
from pathlib import Path
from tools.pathogen_normalization import canonical_key
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_POLICY = ROOT / "rules/cmv_evidence_reporting_shadow_v2_fallback.json"
POSITIVE = {"Picked", "Possible", "Fallback-Possible"}
ORDINARY = {"Picked", "Possible"}



def finalize_patient(patient: dict, cmv_rows: list[dict], policy: dict) -> tuple[list[dict], dict]:
    """Use all pathogens to enforce ordinary-vs-fallback mutual exclusivity."""
    rows = copy.deepcopy(cmv_rows)
    originals = {canonical_key(c["organism_name"]): c for c in patient["all_candidates"]}
    if len(originals) != len(patient["all_candidates"]):
        raise ValueError("Duplicate patient candidate identities")
    if len({canonical_key(r["organism_name"]) for r in rows}) != len(rows):
        raise ValueError("Duplicate CMV rows")
    reporting = {key: {"organism_name": c["organism_name"],
                      "tier": c["reporting"]["tier"],
                      "selected": bool(c["reporting"]["selected_for_complete_report"])}
                 for key, c in originals.items()}
    for row in rows:
        key = canonical_key(row["organism_name"])
        if str(row["patient_id"]) != str(patient["patient_id"]) or key not in reporting:
            raise ValueError("Patient/candidate linkage mismatch")
        decision = row["decision"]
        reporting[key].update(tier=decision["shadow_reporting_tier"],
                              selected=decision["selected_as_pulmonary_candidate"])
    has_ordinary = any(c["selected"] and c["tier"] in ORDINARY for c in reporting.values())
    restored = []
    for row in rows:
        decision = row["decision"]
        key = canonical_key(row["organism_name"])
        original = originals[key]["reporting"]
        row["decision_before_fallback_finalization"] = copy.deepcopy(decision)
        decision["fallback_reporting_policy"] = policy["policy_id"]
        decision["fallback_preserved"] = False
        eligible = (not has_ordinary and original["tier"] == "Fallback-Possible"
                    and original["selected_for_complete_report"]
                    and decision["shadow_reporting_tier"] == "Context"
                    and decision["axes"]["detection"] != "D0"
                    and decision["axes"]["activity"] != "A3")
        if eligible:
            decision.update(shadow_reporting_tier="Fallback-Possible", selected_as_pulmonary_candidate=True,
                            doctor_visible=True, route="preserve_validated_fallback_when_no_ordinary_selection",
                            report_label_zh=policy["fallback"]["label_zh"], fallback_preserved=True)
            reporting[key].update(tier="Fallback-Possible", selected=True)
            restored.append(row["organism_name"])
        decision["included_in_positive_prediction_endpoint"] = bool(
            decision["selected_as_pulmonary_candidate"] and decision["shadow_reporting_tier"] in POSITIVE)
        decision["fallback_means_clinical_disease_confirmed"] = False
    # The frozen production output already enforces exclusivity; do not silently
    # repair a newly conflicting fixture/input or re-rank unrelated candidates.
    final = [c for c in reporting.values() if c["selected"] and c["tier"] in POSITIVE]
    if has_ordinary and any(c["tier"] == "Fallback-Possible" for c in final):
        raise ValueError("Existing patient-level fallback conflicts with ordinary selection")
    report = {"patient_id": patient["patient_id"], "mode": "cmv_reporting_shadow_with_preserved_fallback",
              "has_picked_or_ordinary_possible": has_ordinary, "fallback_restored": restored,
              "final_report_after_shadow": final,
              "positive_prediction_definition": "Picked + Possible + Fallback-Possible",
              "has_at_least_one_reported_candidate": bool(final),
              "candidate_count": len(originals), "production_changed": False}
    return rows, report

"""Source-verified, three-axis viral interpretation; no answer access."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from tools.pathogen_normalization import canonical_key

DEFAULT_POLICY = Path(__file__).resolve().parents[1] / "rules/viral_evidence_axes_shadow_v1.json"


def interpret(candidate: dict[str, Any], *, episode_id: str,
              claims: list[dict[str, Any]] | None = None,
              policy: dict[str, Any] | None = None) -> dict[str, Any]:
    policy = policy or json.loads(DEFAULT_POLICY.read_text(encoding="utf-8"))
    organism = canonical_key(candidate.get("organism_name"))
    taxonomy = candidate.get("taxonomy_profile") or {}
    signals = [s for s in (candidate.get("analytical_profile") or {}).get("per_test_signals") or []
               if float(s.get("reads") or 0) > 0
               and s.get("qc_status") in policy["detection_qc_statuses"]]
    precise = taxonomy.get("mapping_status") == "exact_species"
    accepted, rejected = [], []
    for claim in claims or []:
        failures = [f"missing:{k}" for k in policy["required_claim_provenance"] if not claim.get(k)]
        if canonical_key(claim.get("organism_key")) != organism:
            failures.append("organism_mismatch")
        if str(claim.get("episode_id")) != str(episode_id):
            failures.append("episode_mismatch")
        if claim.get("review_complete") is not True or claim.get("event_aligned") is not True:
            failures.append("source_review_or_event_alignment_incomplete")
        if claim.get("claim_type") not in {"detection", *policy["activity_claims"], *policy["organ_claims"]}:
            failures.append("unsupported_claim_type")
        if claim.get("claim_type") in {"reviewed_probable_pulmonary_disease", "reviewed_proven_pulmonary_disease"}:
            allowed_bases = policy["organ_attribution_bases"][claim["claim_type"]]
            if claim.get("organ") != "lung" or claim.get("attribution_basis") not in allowed_bases:
                failures.append("no_lung_specific_attribution")
            if not claim.get("clinical_interpretation"):
                failures.append("missing_clinical_interpretation")
        if failures:
            rejected.append({"evidence_id": claim.get("evidence_id"), "failures": failures})
        else:
            accepted.append(claim)
    kinds = {c.get("claim_type") for c in accepted}
    detection_claims = [c for c in accepted if c.get("claim_type") == "detection"]
    detection = "D1" if (signals and precise) or detection_claims else "D0"
    # Corroboration needs a verified second biological specimen or independent
    # clinical assay; seq_id/molecule counts never qualify on their own.
    independent = [c for c in detection_claims if c.get("independence_class") in
                   {"independent_biological_specimen", "independent_clinical_assay"}]
    independent_ids = {c["evidence_id"] for c in independent}
    if detection != "D0" and independent and ((signals and precise) or len(independent_ids) >= 2):
        detection = "D2"
    activity = "A0"
    if "clinician_current_activity" in kinds:
        activity = "A1"
    if "validated_activity_interpretation" in kinds:
        activity = "A2"
    if "resolved_or_asymptomatic_shedding" in kinds:
        activity = "A3"
    organ = "O0"
    if "compatible_pulmonary_syndrome" in kinds:
        organ = "O1"
    if "reviewed_probable_pulmonary_disease" in kinds:
        organ = "O2"
    if "reviewed_proven_pulmonary_disease" in kinds:
        organ = "O3"
    tier = "Audit" if detection == "D0" else "Context"
    if detection != "D0" and organ == "O2":
        tier = "Possible"
    if detection != "D0" and organ == "O3":
        tier = "Picked"
    gaps = []
    if not accepted:
        gaps.append("clinical_claim_extraction_or_source_verification_incomplete")
    if activity == "A0":
        gaps.append("activity_unknown")
    if organ in {"O0", "O1"}:
        gaps.append("pulmonary_attribution_unknown")
    ledger: dict[str, list[str]] = {}
    for c in accepted:
        axis = "detection" if c["claim_type"] == "detection" else (
            "activity" if c["claim_type"] in policy["activity_claims"] else "organ_support")
        ledger.setdefault(c["evidence_id"], [])
        if axis not in ledger[c["evidence_id"]]:
            ledger[c["evidence_id"]].append(axis)
    return {"detection": detection, "activity": activity, "organ_support": organ,
            "shadow_tier": tier, "production_changed": False,
            "report_label": policy["unknown_reporting_label"] if tier == "Context" else tier,
            "analytical_positive_branch_count": len(signals),
            "technical_branches_are_independent_clinical_evidence": False,
            "gaps": gaps, "accepted_claim_count": len(accepted), "rejected_claims": rejected,
            "evidence_consumption_ledger": ledger,
            "unique_supporting_evidence_count": len(ledger),
            "clinical_absence_inferred_from_missing_claims": False}

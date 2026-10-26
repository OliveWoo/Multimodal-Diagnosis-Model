"""Shared current deterministic defaults; historical frozen bundles stay intact.

Use run_current_integrated_pipeline for the full source-bound CMV integration.
Lower-level runners retain explicit policy overrides for historical replays.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANALYTICAL_POLICY = ROOT / "rules/test_aware_deterministic_shadow_v6_identity_resolved_20261020.json"
HISTORY_POLICY = ROOT / "rules/test_aware_history_route_v4c_combined.json"
CLINICAL_POLICY = ROOT / "rules/test_aware_clinical_deterministic_shadow_v1_20261020.json"
PROMOTION_POLICY = ROOT / "rules/test_aware_clinical_promotion_v9_pathogen_specific_convergence_20261020.json"
REPORTING_POLICY = ROOT / "rules/test_aware_unified_reporting_v7_respiratory_virus_possible_20261020.json"
CMV_POLICY = ROOT / "rules/cmv_evidence_reporting_shadow_v1.json"
FALLBACK_POLICY = ROOT / "rules/cmv_evidence_reporting_shadow_v2_fallback.json"
SOURCE_CLARIFICATIONS = ROOT / "rules/clinical_source_clarifications_template.json"

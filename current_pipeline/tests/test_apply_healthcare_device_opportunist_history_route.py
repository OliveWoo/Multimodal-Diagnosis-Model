import json
from pathlib import Path

from tools.apply_healthcare_device_opportunist_history_route import apply_patient


POLICY = json.loads(
    Path("rules/healthcare_device_opportunist_history_route_v2.json").read_text(
        encoding="utf-8"
    )
)


def patient_payload() -> dict:
    return {
        "patient_id": "synthetic",
        "hap_mdr_history_profile": {
            "support_level": "hap_setting_context",
            "matched_evidence": [
                {
                    "evidence_id": "SYNTHETIC:HAP:1",
                    "categories": ["hap_setting"],
                }
            ],
        },
        "candidates": [
            {
                "organism_name": "Synthetic healthcare-associated bacterium",
                "taxonomy_profile": {
                    "primary_rule_family": "skin_airway_colonizer_prone",
                    "clinical_traits": [
                        "healthcare_associated",
                        "pulmonary_pathogen_possible",
                        "colonization_possible",
                    ],
                },
                "history_adjusted_route": "Context",
                "history_route_changed": False,
                "history_routed_decision": "review_context_needed",
                "history_rule_ids": ["HIST-S1-PRESERVE-ANALYTIC-ROUTE"],
                "history_evidence_ids": [],
                "consumed_evidence_ids": [],
                "history_adjustment_blockers": [
                    "no implemented history route applies to this taxonomy family"
                ],
                "specimen_context": "lower_respiratory",
                "analytical_profile": {
                    "best_rank_in_retained_universe": 2,
                    "selected_positive_test_count": 1,
                    "reproducibility_axis": False,
                    "cross_molecule_selected": False,
                },
                "colonization_interpretation": {
                    "has_direct_nonsterile_match": False,
                    "has_explicit_colonization": False,
                    "has_explicit_current_not_infection": False,
                },
            }
        ],
    }


def test_hap_context_candidate_moves_to_priority() -> None:
    output, rows = apply_patient(patient_payload(), POLICY)
    candidate = output["candidates"][0]
    assert candidate["history_adjusted_route"] == "Priority"
    assert candidate["history_routed_decision"] == "review_high_priority"
    assert candidate["history_policy_family"] == "healthcare_device_opportunist"
    assert rows[0]["promoted_to_priority"] is True


def test_missing_taxonomy_trait_blocks_route_change() -> None:
    payload = patient_payload()
    payload["candidates"][0]["taxonomy_profile"]["clinical_traits"].remove(
        "pulmonary_pathogen_possible"
    )
    output, rows = apply_patient(payload, POLICY)
    assert output["candidates"][0]["history_adjusted_route"] == "Context"
    assert rows[0]["promoted_to_priority"] is False
    assert "required_taxonomy_traits" in rows[0]["blockers"]


def test_explicit_colonization_blocks_route_change() -> None:
    payload = patient_payload()
    payload["candidates"][0]["colonization_interpretation"].update(
        {
            "has_direct_nonsterile_match": True,
            "has_explicit_colonization": True,
        }
    )
    output, rows = apply_patient(payload, POLICY)
    assert output["candidates"][0]["history_adjusted_route"] == "Context"
    assert rows[0]["promoted_to_priority"] is False
    assert "no_direct_colonization_or_not_infection" in rows[0]["blockers"]


def test_policy_contains_no_patient_or_organism_allowlist() -> None:
    serialized = json.dumps(POLICY, sort_keys=True).casefold()
    assert "corynebacterium" not in serialized
    assert "patient_ids" not in POLICY
    assert "organism_names" not in POLICY


def test_patient_cap_keeps_two_and_leaves_overflow_in_context() -> None:
    payload = patient_payload()
    template = payload["candidates"][0]
    payload["candidates"] = []
    for rank in (1, 2, 3):
        value = json.loads(json.dumps(template))
        value["organism_name"] = f"Synthetic candidate {rank}"
        value["analytical_profile"]["best_rank_in_retained_universe"] = rank
        payload["candidates"].append(value)
    output, rows = apply_patient(payload, POLICY)
    assert [row["promoted_to_priority"] for row in rows] == [True, True, False]
    assert [row["after_route"] for row in rows] == ["Priority", "Priority", "Context"]
    assert "patient_representative_cap" in rows[2]["blockers"]


def test_exact_hospital_support_is_exempt_from_patient_cap() -> None:
    payload = patient_payload()
    template = payload["candidates"][0]
    payload["candidates"] = []
    for rank in (1, 2, 3):
        value = json.loads(json.dumps(template))
        value["organism_name"] = f"Synthetic candidate {rank}"
        value["analytical_profile"]["best_rank_in_retained_universe"] = rank
        if rank == 3:
            value["hospital_profile"] = {"exact_or_alias_evidence": True}
        payload["candidates"].append(value)
    output, rows = apply_patient(payload, POLICY)
    assert all(row["promoted_to_priority"] for row in rows)
    assert rows[2]["patient_cap_exempt_exact_hospital_support"] is True

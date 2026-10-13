import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "rules/kh_legacy_pipeline_finalized_v1.json"


def load_contract() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_finalized_release_has_only_legacy_labels() -> None:
    labels = set(load_contract()["output_labels"].values())
    assert labels == {"Picked", "Priority", "Context", "Audit"}
    assert "Possible" not in labels


def test_history_cannot_directly_create_picked() -> None:
    constraints = " ".join(load_contract()["constraints"])
    assert "HAP or airway-device history alone never creates Picked" in constraints
    assert "Priority to Picked requires a separate disease-specific gate" in constraints


def test_representative_cap_is_two_with_exact_support_exemption() -> None:
    policy_path = ROOT / "rules/healthcare_device_opportunist_history_route_v2.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    assert policy["max_promoted_representatives_per_patient"] == 2
    assert policy["cap_exempt_when_exact_hospital_support"] is True


def test_release_contract_contains_no_patient_or_organism_allowlist() -> None:
    contract = load_contract()
    serialized = json.dumps(contract, ensure_ascii=False).casefold()
    assert "patient_allowlist" not in serialized
    assert "organism_allowlist" not in serialized

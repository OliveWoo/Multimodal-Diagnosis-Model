import json
from pathlib import Path

from tools.run_two_hospital_22_legacy_finalized_compat_v1 import (
    build_selected_dna_payload,
    representative_dna_tests,
    selected_rank,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "rules/two_hospital_22_legacy_finalized_compat_v1.json"


def test_contract_has_two_disjoint_eleven_patient_cohorts() -> None:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    cohorts = contract["cohorts"]
    assert [len(item["expected_patient_ids"]) for item in cohorts] == [11, 11]
    first = set(cohorts[0]["expected_patient_ids"])
    second = set(cohorts[1]["expected_patient_ids"])
    assert first.isdisjoint(second)
    assert len(first | second) == 22
    assert "Possible" not in json.dumps(contract, ensure_ascii=False)
    assert cohorts[0]["summary_suffix"] == (
        "release_input_new_mngs_11_evidence_bridge_v1"
    )
    assert cohorts[0]["required_summary_fields"] == [
        "hospital_organism_evidence"
    ]
    assert cohorts[1]["required_summary_fields"] == [
        "hospital_organism_evidence"
    ]


def test_kh18_contract_requires_structured_hospital_evidence() -> None:
    contract = json.loads(
        (ROOT / "rules/kh_18_legacy_finalized_compat_v1.json").read_text(
            encoding="utf-8-sig"
        )
    )
    cohort = contract["cohorts"][0]
    assert cohort["summary_suffix"] == "final_summary_51patients_no_zero_20260930"
    assert cohort["required_summary_fields"] == ["hospital_organism_evidence"]
    assert cohort["summary_schema_role"] == "structured_hospital_evidence_bridge"


def test_representative_test_prefers_lowest_numbered_dna_and_excludes_rna() -> None:
    raw = {
        "cases": [
            {
                "case_review_id": "event-1",
                "tests": [
                    {"nucleic_type": "RNA", "test_number": 1, "seq_id": "rna"},
                    {"nucleic_type": "DNA", "test_number": 3, "seq_id": "dna-3"},
                    {"nucleic_type": "DNA", "test_number": 2, "seq_id": "dna-2"},
                ],
            }
        ]
    }
    selected = representative_dna_tests(raw)
    assert selected["event-1"]["seq_id"] == "dna-2"


def test_representative_test_skips_empty_earlier_dna() -> None:
    raw = {
        "cases": [
            {
                "case_review_id": "event-1",
                "tests": [
                    {"nucleic_type": "DNA", "test_number": 1, "seq_id": "empty"},
                    {"nucleic_type": "DNA", "test_number": 2, "seq_id": "positive"},
                ],
            }
        ],
        "organisms": [
            {
                "name": "Example organism",
                "measurements": [
                    {
                        "case_review_id": "event-1",
                        "seq_id": "positive",
                        "reads": 10,
                        "selected": True,
                    }
                ],
            }
        ],
    }
    selected = representative_dna_tests(raw)
    assert selected["event-1"]["seq_id"] == "positive"


def test_payload_keeps_only_positive_selected_measurements_from_chosen_dna() -> None:
    raw = {
        "patient_id": "58",
        "cases": [
            {
                "case_review_id": "event-1",
                "specimen_code": "BALF-1",
                "tests": [
                    {"nucleic_type": "DNA", "test_number": 1, "seq_id": "dna-1"},
                    {"nucleic_type": "RNA", "test_number": 2, "seq_id": "rna-2"},
                ],
            }
        ],
        "organisms": [
            {
                "name": "Included organism",
                "category": "bacteria",
                "measurements": [
                    {
                        "case_review_id": "event-1",
                        "seq_id": "dna-1",
                        "reads": 12,
                        "selected": True,
                        "selection_tags": ["Code_NTC=00", "Code_RK_NTC=RK0"],
                    }
                ],
            },
            {
                "name": "Zero-read organism",
                "category": "bacteria",
                "measurements": [
                    {
                        "case_review_id": "event-1",
                        "seq_id": "dna-1",
                        "reads": 0,
                        "selected": True,
                        "selection_tags": ["Code_NTC=00"],
                    }
                ],
            },
            {
                "name": "RNA organism",
                "category": "virus",
                "measurements": [
                    {
                        "case_review_id": "event-1",
                        "seq_id": "rna-2",
                        "reads": 999,
                        "selected": True,
                        "selection_tags": ["Code_NTC=00", "Code_RK_NTC=RK0"],
                    }
                ],
            },
        ],
    }
    payload, audit = build_selected_dna_payload(raw)
    names = [item["organism_name"] for item in payload["records"][0]["candidates"]]
    assert names == ["Included organism"]
    assert payload["ranking_metadata"]["rna_excluded"] is True
    assert audit["events"][0]["selected_candidate_count"] == 1


def test_selection_tag_rank_order_is_stable() -> None:
    assert selected_rank(["Code_NTC=00", "Code_RK_NTC=RK00"]) == 1
    assert selected_rank(["Code_NTC=00", "Code_RK_NTC=RK8"]) == 2
    assert selected_rank(["Code_RK_NTC=RK0"]) == 3
    assert selected_rank(["Code_RK_NTC=RK8"]) == 4
    assert selected_rank(["Code_NTC=00"]) == 5


def test_blank_measurement_case_id_links_by_selected_seq_id() -> None:
    raw = {
        "patient_id": "107",
        "cases": [
            {
                "case_review_id": None,
                "specimen_code": "240111001A",
                "tests": [
                    {
                        "nucleic_type": "DNA",
                        "test_number": 1,
                        "seq_id": "24S0000015",
                    }
                ],
            }
        ],
        "organisms": [
            {
                "name": "Pneumocystis jirovecii",
                "category": "fungal",
                "measurements": [
                    {
                        "case_review_id": None,
                        "seq_id": "24S0000015",
                        "reads": 57,
                        "selected": True,
                        "selection_tags": ["Code_RK_NTC=RK0"],
                    }
                ],
            }
        ],
    }
    payload, audit = build_selected_dna_payload(raw)
    assert payload["records"][0]["candidates"][0]["organism_name"] == (
        "Pneumocystis jirovecii"
    )
    assert audit["events"][0]["selected_candidate_count"] == 1

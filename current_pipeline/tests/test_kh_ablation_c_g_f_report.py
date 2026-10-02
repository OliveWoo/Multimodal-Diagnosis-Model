from __future__ import annotations

from tools.build_kh_ablation_c_g_f_report import compare_candidate_maps


def test_strict_comparison_reports_only_picked_changes() -> None:
    before = {
        (1, "organisma"): {
            "patient_id": "1",
            "organism_name": "Organism A",
            "_tier": "Picked",
        }
    }
    after = {
        (1, "organisma"): {
            "patient_id": "1",
            "organism_name": "Organism A",
            "_tier": "High",
        },
        (1, "organismb"): {
            "patient_id": "1",
            "organism_name": "Organism B",
            "_tier": "Context",
        },
    }
    changes = compare_candidate_maps(
        "C", before, after, {1: ["Organism A"]}, strict_only=True
    )
    assert len(changes) == 1
    assert changes[0]["change_type"] == "picked_lost"
    assert changes[0]["answer_status"] == "matched_answer"


def test_history_comparison_retains_tier_change_and_consumed_evidence() -> None:
    before = {
        (30, "bacteroidesfragilis"): {
            "patient_id": "30",
            "organism_name": "Bacteroides fragilis",
            "_tier": "Context",
        }
    }
    after = {
        (30, "bacteroidesfragilis"): {
            "patient_id": "30",
            "organism_name": "Bacteroides fragilis",
            "_tier": "High",
            "rule_ids": "HIST-ASP-01",
            "consumed_evidence_ids": "HX-30-1",
        }
    }
    changes = compare_candidate_maps(
        "F", before, after, {30: ["Bacteroides fragilis"]}
    )
    assert len(changes) == 1
    assert changes[0]["change_type"] == "report_tier_changed"
    assert changes[0]["rule_ids"] == "HIST-ASP-01"
    assert changes[0]["evidence_ids"] == "HX-30-1"

import json
from pathlib import Path

from tools.materialize_multi_assay_ranked_context import run


def test_missing_clinical_patient_is_reported(tmp_path: Path) -> None:
    compact = tmp_path / "compact"
    packets = compact / "patient_packets"
    patient_root = tmp_path / "patients"
    patient_root.mkdir(parents=True)
    packets.mkdir(parents=True)
    packet = {
        "patient_id": "7",
        "cases": [
            {
                "case_review_id": "C7",
                "specimen_code": "C7",
                "collected_time": "2026-01-01",
                "specimen_site": "BALF",
                "scorer_entries": [],
            }
        ],
    }
    (packets / "NGS_patient_7_compact_entry_shadow.json").write_text(
        json.dumps(packet), encoding="utf-8"
    )
    summary = run(compact, patient_root, tmp_path / "summary.json", replace_existing=True)
    assert summary["written_patient_count"] == 0
    assert summary["missing_patient_directory_ids"] == [7]

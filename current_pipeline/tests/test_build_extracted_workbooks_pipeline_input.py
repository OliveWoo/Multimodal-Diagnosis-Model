import json
from pathlib import Path

from tools.build_extracted_workbooks_pipeline_input import run


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_only_patients_with_current_mngs_are_forwarded(tmp_path: Path) -> None:
    source = tmp_path / "source"
    output = tmp_path / "output"
    ready = source / "extracted_patient_3_json"
    pending = source / "extracted_patient_1_json"
    ready.mkdir(parents=True)
    pending.mkdir(parents=True)
    for patient_dir, patient_id in ((ready, 3), (pending, 1)):
        write_json(patient_dir / f"extracted_patient_{patient_id}_CBC.json", [])
        write_json(patient_dir / "enrichment_manifest.json", {"patient_id": patient_id})
    write_json(
        ready / "extracted_patient_3_all_mNGS_DNA_RNA_RK_NTC_by_species.json",
        {"cases": [{"tests": [{"seq_id": "S1"}]}]},
    )
    write_json(ready / "extracted_patient_3_mNGS_ranked_candidates.json", {"records": []})

    summary = run(source, output)

    assert summary["forwarded_patient_ids"] == [3]
    assert summary["pending_mngs_patient_ids"] == [1]
    assert summary["multi_assay_case_count"] == 1
    assert summary["multi_assay_test_count"] == 1
    assert (output / "patients" / "NGS_patient_3_json" / "NGS_patient_3_CBC.json").is_file()
    assert not (output / "patients" / "NGS_patient_1_json").exists()
    pending_payload = json.loads((output / "mngs_pending_patients.json").read_text())
    assert pending_payload["patient_count"] == 1

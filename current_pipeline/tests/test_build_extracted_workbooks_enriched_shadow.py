from pathlib import Path

from tools.build_extracted_workbooks_enriched_shadow import (
    build_linkage_audit,
    normalize_value,
    patient_id_from_name,
)


def test_patient_id_and_numeric_normalization() -> None:
    assert patient_id_from_name("NGS_patient_10_json") == 10
    assert patient_id_from_name("extracted_patient_39_CBC.json") == 39
    assert normalize_value("1,200.00") == "1200"


def test_same_number_linkage_uses_content_not_identifier(tmp_path: Path) -> None:
    legacy = tmp_path / "NGS_patient_3_json"
    staging_3 = tmp_path / "extracted_patient_3_json"
    staging_4 = tmp_path / "extracted_patient_4_json"
    for directory in (legacy, staging_3, staging_4):
        directory.mkdir()
    (legacy / "NGS_patient_3_CBC.json").write_text(
        '[{"item":"WBC","value":"5.2","reported_time":"2024-01-02"}]',
        encoding="utf-8",
    )
    (staging_3 / "extracted_patient_3_CBC.json").write_text(
        '[{"item":"WBC","value":"5.20","collected_time":"2024/01/02 08:00"}]',
        encoding="utf-8",
    )
    (staging_4 / "extracted_patient_4_CBC.json").write_text(
        '[{"item":"WBC","value":"8.1","collected_time":"2024/01/02 08:00"}]',
        encoding="utf-8",
    )
    rows = build_linkage_audit({3: legacy}, {3: staging_3, 4: staging_4})
    assert rows[0]["best_match_patient_id"] == 3
    assert rows[0]["same_number_rank"] == 1
    assert rows[0]["same_number_exact_dated_lab_matches"] == 1

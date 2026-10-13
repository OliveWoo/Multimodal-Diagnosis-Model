from tools.build_44_patient_small_agent_workspace import (
    merge_records,
    parser_source_key,
    source_aware_image_prompt_records,
)


def test_parser_source_key_accepts_both_layer_shapes() -> None:
    existing = {
        "parser_reference": {
            "parser_file": "synthetic_CBC.json",
            "json_pointer": "/49",
        }
    }
    supplement = {
        "supplement_provenance": {
            "parser_file": "synthetic_CBC.json",
            "json_pointer": "/49",
        }
    }

    assert parser_source_key(existing) == parser_source_key(supplement)


def test_merge_records_deduplicates_reconciled_supplement_by_parser_source() -> None:
    existing = [{
        "test": "CBC",
        "item": "WBC",
        "value": "10.0",
        "reported_time": "2000-01-10 06:00",
        "parser_reference": {
            "parser_file": "synthetic_CBC.json",
            "json_pointer": "/49",
        },
    }]
    supplement = [{
        "test": "CBC",
        "item": "WBC",
        "value": "10.0",
        "collected_time": "2000-01-10 06:00",
        "supplement_provenance": {
            "parser_file": "synthetic_CBC.json",
            "json_pointer": "/49",
        },
    }]

    merged, duplicate_count = merge_records(existing, supplement)

    assert merged == existing
    assert duplicate_count == 1


def test_merge_records_keeps_distinct_parser_source() -> None:
    existing = [{
        "test": "CBC",
        "item": "WBC",
        "value": "10.0",
        "parser_reference": {
            "parser_file": "synthetic_CBC.json",
            "json_pointer": "/49",
        },
    }]
    supplement = [{
        "test": "CBC",
        "item": "WBC",
        "value": "20.0",
        "supplement_provenance": {
            "parser_file": "synthetic_CBC.json",
            "json_pointer": "/50",
        },
    }]

    merged, duplicate_count = merge_records(existing, supplement)

    assert len(merged) == 2
    assert duplicate_count == 0


def test_image_prompt_drops_source_time_only_copy_but_keeps_distinct_study() -> None:
    records = [
        {
            "exam_type": "Chest AP",
            "collected_time": "2000-01-10 09:00",
            "reported_time": "2000-01-10 10:00",
            "findings": ["Pneumonia in both lungs.", "Small pleural effusions."],
        },
        {
            "exam_type": "Chest AP",
            "source_time": "2000-01-10 09:00",
            "findings": ["Pneumonia in both lungs. Small pleural effusions."],
        },
        {
            "exam_type": "KUB",
            "source_time": "2000-01-10 11:00",
            "findings": ["Prominent bowel gas pattern."],
        },
    ]

    selected, omitted = source_aware_image_prompt_records(records)

    assert selected == [records[0], records[2]]
    assert len(omitted) == 1
    assert omitted[0]["reason"] == "source_time_only_duplicate_of_richer_image_record"

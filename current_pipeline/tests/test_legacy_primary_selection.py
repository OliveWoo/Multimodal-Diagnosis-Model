from tools.legacy_primary_selection import (
    FALLBACK_ROLE,
    FORMAL_ROLE,
    NO_CANDIDATE_ROLE,
    build_primary_selections,
    complete_report_rows,
    doctor_report_rows,
)


def row(patient, organism, label, rank, family="typical_respiratory_pathogen"):
    return {
        "patient_id": patient,
        "organism_name": organism,
        "final_label": label,
        "best_rank": rank,
        "taxonomy_family": family,
        "specimen_context": "lower_respiratory",
    }


def test_existing_picked_is_primary() -> None:
    rows = [
        row(1, "Context organism", "Context", 1),
        row(1, "Picked organism", "Picked", 2),
    ]
    selected = build_primary_selections(rows, [1])
    assert selected[0]["organism_name"] == "Picked organism"
    assert selected[0]["primary_selection_role"] == FORMAL_ROLE


def test_no_picked_gets_one_low_confidence_fallback() -> None:
    rows = [
        row(2, "Second", "Context", 2),
        row(2, "First", "Priority", 3),
    ]
    selected = build_primary_selections(rows, [2])
    assert selected[0]["organism_name"] == "First"
    assert selected[0]["primary_selection_role"] == FALLBACK_ROLE
    assert selected[0]["doctor_visible_label"] == "Fallback-BestAvailable"
    report = complete_report_rows(rows, selected)
    assert [item["organism_name"] for item in report] == ["First"]


def test_unmapped_does_not_beat_interpretable_candidate() -> None:
    rows = [
        row(3, "Unknown organism", "Priority", 1, "unmapped_or_uncertain"),
        row(3, "Known organism", "Context", 3),
    ]
    selected = build_primary_selections(rows, [3])
    assert selected[0]["organism_name"] == "Known organism"


def test_zero_candidate_is_explicit_not_fabricated() -> None:
    selected = build_primary_selections([], [4])
    assert selected[0]["organism_name"] == ""
    assert selected[0]["primary_selection_role"] == NO_CANDIDATE_ROLE


def test_doctor_report_labels_fallback_separately_from_picked() -> None:
    rows = [row(5, "Fallback organism", "Context", 1)]
    selected = build_primary_selections(rows, [5])
    report = doctor_report_rows(rows, selected)
    assert report[0]["doctor_visible_tier"] == "Fallback-BestAvailable"
    assert report[0]["final_label"] == "Context"
    assert report[0]["is_primary_output"] is True

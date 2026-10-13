"""Select one answer-blind primary organism for the finalized legacy flow.

The historical user-facing flow always exposed at least one best-available
organism when candidates existed.  The frozen legacy decision pipeline kept
the candidate tiers but lost that final reporting step.  This module restores
the reporting contract without relabeling a low-confidence fallback as a
formal ``Picked`` result.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Iterable


FORMAL_ROLE = "formal_picked_primary"
FALLBACK_ROLE = "fallback_best_available"
NO_CANDIDATE_ROLE = "no_candidate_available"

_LABEL_ORDER = {
    "Picked": 0,
    "Priority": 1,
    "Context": 2,
    "Audit": 3,
}

# These families can remain visible in the audit, but should not win the
# mandatory fallback while a more interpretable organism is available.
_LAST_RESORT_FAMILIES = {
    "unmapped_or_uncertain",
    "commensal_virome_or_endogenous_element",
}

_FAMILY_ORDER = {
    "typical_respiratory_pathogen": 0,
    "hospital_or_nonfermenter_gnb": 1,
    "high_consequence_opportunistic": 1,
    "other_respiratory_virus": 1,
    "mold_or_opportunistic_fungus": 2,
    "oral_aspiration_or_anaerobe": 3,
    "herpes_or_reactivation_virus": 4,
    "environmental_low_specificity": 5,
    "skin_airway_colonizer_prone": 6,
    "candida_or_yeast": 6,
    "commensal_virome_or_endogenous_element": 8,
    "unmapped_or_uncertain": 9,
}

_SPECIMEN_ORDER = {
    "lower_respiratory": 0,
    "sterile_or_systemic": 1,
    "unknown": 2,
}


def _as_rank(value: Any) -> int:
    try:
        rank = int(float(str(value)))
    except (TypeError, ValueError):
        return 999
    return rank if rank > 0 else 999


def _candidate_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        _LABEL_ORDER.get(str(row.get("final_label") or ""), 99),
        _as_rank(row.get("best_rank")),
        _FAMILY_ORDER.get(str(row.get("taxonomy_family") or ""), 7),
        _SPECIMEN_ORDER.get(str(row.get("specimen_context") or ""), 3),
        str(row.get("organism_name") or "").casefold(),
    )


def _eligible_pool(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    non_last_resort = [
        row
        for row in rows
        if str(row.get("taxonomy_family") or "") not in _LAST_RESORT_FAMILIES
    ]
    return non_last_resort or rows


def build_primary_selections(
    final_rows: list[dict[str, Any]],
    patient_ids: Iterable[int],
    cohort_by_patient: dict[int, str] | None = None,
) -> list[dict[str, Any]]:
    """Return exactly one primary row per patient when any candidate exists.

    Formal Picked candidates always take precedence.  If none exists, one
    answer-blind fallback is chosen from the existing Priority/Context/Audit
    candidates.  Unmapped and endogenous/commensal-virome families are only
    eligible when no more interpretable candidate exists.
    """

    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in final_rows:
        grouped[int(row["patient_id"])].append(row)

    output: list[dict[str, Any]] = []
    for patient in sorted({int(value) for value in patient_ids}):
        rows = grouped.get(patient, [])
        picked = [row for row in rows if row.get("final_label") == "Picked"]
        if picked:
            chosen = min(picked, key=_candidate_key)
            role = FORMAL_ROLE
            confidence = "formal_picked"
            source_label = "Picked"
            reason = "existing_formal_picked_ranked_primary"
        elif rows:
            chosen = min(_eligible_pool(rows), key=_candidate_key)
            role = FALLBACK_ROLE
            confidence = "low_confidence_best_available"
            source_label = str(chosen.get("final_label") or "")
            reason = (
                "no_formal_picked_choose_one_by_tier_rank_family_specimen"
            )
        else:
            chosen = {}
            role = NO_CANDIDATE_ROLE
            confidence = "not_available"
            source_label = ""
            reason = "no_candidate_entered_legacy_pipeline"

        record = {
            "patient_id": patient,
            "organism_name": str(chosen.get("organism_name") or ""),
            "primary_selection_role": role,
            "selection_confidence": confidence,
            "source_final_label": source_label,
            "taxonomy_family": str(chosen.get("taxonomy_family") or ""),
            "specimen_context": str(chosen.get("specimen_context") or ""),
            "best_rank": str(chosen.get("best_rank") or ""),
            "selection_reason": reason,
            "doctor_visible_label": (
                "Picked" if role == FORMAL_ROLE else
                "Fallback-BestAvailable" if role == FALLBACK_ROLE else
                "No candidate available"
            ),
        }
        if cohort_by_patient is not None:
            record["cohort_id"] = cohort_by_patient.get(patient, "")
        output.append(record)
    return output


def complete_report_rows(
    final_rows: list[dict[str, Any]],
    primary_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Return all formal Picked rows plus one fallback for no-Picked patients."""

    result = [row for row in final_rows if row.get("final_label") == "Picked"]
    by_key = {
        (int(row["patient_id"]), str(row.get("organism_name") or "")): row
        for row in final_rows
    }
    for primary in primary_rows:
        if primary.get("primary_selection_role") != FALLBACK_ROLE:
            continue
        key = (int(primary["patient_id"]), str(primary.get("organism_name") or ""))
        source = by_key.get(key)
        if source is not None:
            result.append(source)
    return result


def doctor_report_rows(
    final_rows: list[dict[str, Any]],
    primary_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build the doctor-visible legacy report with an explicit fallback tier."""

    primary_by_patient = {
        int(row["patient_id"]): row for row in primary_rows
    }
    output: list[dict[str, Any]] = []
    for source in complete_report_rows(final_rows, primary_rows):
        row = dict(source)
        patient = int(row["patient_id"])
        primary = primary_by_patient[patient]
        is_primary = (
            str(row.get("organism_name") or "")
            == str(primary.get("organism_name") or "")
        )
        is_fallback = primary.get("primary_selection_role") == FALLBACK_ROLE
        row.update(
            {
                "doctor_visible_tier": (
                    "Fallback-BestAvailable" if is_fallback else "Picked"
                ),
                "is_primary_output": is_primary,
                "selection_confidence": (
                    "low_confidence_best_available"
                    if is_fallback
                    else "formal_picked"
                ),
                "reporting_caution": (
                    "Best available only; not a formal Picked pathogen."
                    if is_fallback
                    else ""
                ),
            }
        )
        output.append(row)
    return output

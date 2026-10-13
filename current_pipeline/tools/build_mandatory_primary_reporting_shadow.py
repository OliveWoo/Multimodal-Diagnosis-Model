"""Build an answer-blind reporting shadow with one primary candidate per patient.

Existing Picked/Possible decisions are never promoted or demoted.  If a patient
has no selected report row but has a non-QC-excluded positive-read mNGS
candidate, the best generic candidate is exposed as Fallback-BestAvailable.
This is a reporting requirement, not a causal or clinical-certainty claim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


DEFAULT_POLICY = Path("rules/mandatory_primary_reporting_shadow_v1.json")
TRUTHY = {"1", "true", "yes"}
FALLBACK_ROLE = "mandatory_fallback_best_available"


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path, fieldnames: list[str], rows: list[dict[str, Any]]
) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)", path.name)
    if not match:
        raise ValueError(f"Cannot determine patient id: {path}")
    return int(match.group(1))


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().casefold() in TRUTHY


def as_int(value: Any, default: int = 999) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def canonical_name(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").casefold())


def family_priority(family: str, policy: dict[str, Any]) -> int:
    return int(
        (policy.get("family_priority_order") or {}).get(
            family, policy.get("default_family_priority", 99)
        )
    )


def specimen_penalty(context: str, policy: dict[str, Any]) -> int:
    accepted = set(policy.get("compatible_specimen_contexts") or [])
    return 0 if context in accepted else 1


def taxonomy_penalty(status: str, policy: dict[str, Any]) -> int:
    return int(
        (policy.get("taxonomy_mapping_priority") or {}).get(
            status, policy.get("default_taxonomy_mapping_priority", 2)
        )
    )


def selected_row_rank(row: dict[str, Any], policy: dict[str, Any]) -> tuple[Any, ...]:
    tiers = policy.get("formal_tier_order") or {}
    tier = str(row.get("final_reporting_tier") or "Context")
    return (
        int(tiers.get(tier, 99)),
        family_priority(str(row.get("taxonomy_family") or ""), policy),
        specimen_penalty(str(row.get("specimen_context") or ""), policy),
        taxonomy_penalty(str(row.get("identity_resolution") or ""), policy),
        0 if as_bool(row.get("reproducibility_axis")) else 1,
        as_int(row.get("best_rank")),
        -as_int(row.get("selected_positive_test_count"), default=0),
        canonical_name(row.get("organism_name")),
    )


def fallback_features(
    candidate: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any] | None:
    profile = candidate.get("analytical_profile") or {}
    source = candidate.get("source_entry") or {}
    signals = profile.get("per_test_signals") or source.get("selected_test_signals") or []
    accepted_qc = set(policy.get("accepted_qc_statuses") or [])
    valid_signals = [
        signal
        for signal in signals
        if float(signal.get("reads") or 0) > 0
        and str(signal.get("qc_status") or "") in accepted_qc
    ]
    if not valid_signals:
        return None
    family = str(candidate.get("taxonomy_family") or source.get("taxonomy_family") or "")
    context = str(candidate.get("specimen_context") or source.get("specimen_context") or "")
    mapping = str(source.get("taxonomy_mapping_status") or "")
    best_rank = min(
        as_int(signal.get("rank_in_retained_universe")) for signal in valid_signals
    )
    qc_penalty = min(
        0 if str(signal.get("qc_status")) == "evaluable" else 1
        for signal in valid_signals
    )
    reproducible = as_bool(profile.get("reproducibility_axis"))
    cross_molecule = as_bool(profile.get("cross_molecule_selected"))
    reproducibility_penalty = 0 if cross_molecule else 1 if reproducible else 2
    maximum_single_test_reads = max(float(signal.get("reads") or 0) for signal in valid_signals)
    decision_order = {
        "picked_shadow": 0,
        "review_high_priority": 1,
        "review_context_needed": 2,
        "review_low_specificity": 3,
    }
    decision_penalty = decision_order.get(str(candidate.get("decision") or ""), 4)
    rank_vector = (
        decision_penalty,
        family_priority(family, policy),
        specimen_penalty(context, policy),
        taxonomy_penalty(mapping, policy),
        best_rank,
        qc_penalty,
        reproducibility_penalty,
        -len(valid_signals),
        -maximum_single_test_reads,
        canonical_name(candidate.get("organism_name")),
    )
    return {
        "sort_key": rank_vector,
        "rank_vector": list(rank_vector[:-1]),
        "decision_penalty": decision_penalty,
        "family_priority": family_priority(family, policy),
        "specimen_penalty": specimen_penalty(context, policy),
        "taxonomy_penalty": taxonomy_penalty(mapping, policy),
        "best_rank": best_rank,
        "best_qc_status": "evaluable" if qc_penalty == 0 else "partial_evaluable",
        "reproducibility_penalty": reproducibility_penalty,
        "valid_positive_test_count": len(valid_signals),
        "maximum_single_test_reads": maximum_single_test_reads,
        "family": family,
        "specimen_context": context,
        "taxonomy_mapping_status": mapping,
    }


def synthetic_row(
    candidate: dict[str, Any], features: dict[str, Any]
) -> dict[str, Any]:
    profile = candidate.get("analytical_profile") or {}
    source = candidate.get("source_entry") or {}
    hospital = candidate.get("hospital_profile") or {}
    mapping = str(source.get("taxonomy_mapping_status") or "")
    return {
        "patient_id": candidate.get("patient_id"),
        "organism_name": candidate.get("organism_name"),
        "taxonomy_family": features["family"],
        "clinical_decision": candidate.get("decision"),
        "original_selection_role": "analytical_unforwarded_candidate",
        "possible_reporting_role": FALLBACK_ROLE,
        "selected_for_complete_report": "True",
        "possible_outcome": "selected_as_mandatory_primary_fallback",
        "possible_route": "generic_best_available_answer_blind",
        "possible_blockers": "",
        "possible_cautions": (
            "Low-confidence mandatory primary; not a formal causal claim|"
            "Clinical certainty was not increased|"
            "Missing corroborating evidence remains unknown"
        ),
        "specimen_context": features["specimen_context"],
        "precise_identity": str(mapping not in {"", "unmapped", "provisional", "rejected"}),
        "best_rank": features["best_rank"],
        "selected_positive_test_count": features["valid_positive_test_count"],
        "reproducibility_axis": profile.get("reproducibility_axis", False),
        "cross_molecule_selected": profile.get("cross_molecule_selected", False),
        "direct_hospital_level": hospital.get("best_direct_level"),
        "event_aligned_direct_positive_count": 0,
        "repeated_sterile_positive_count": 0,
        "identity_resolution": source.get("taxonomy_identity_resolution"),
        "taxonomy_rule_source": source.get("taxonomy_rule_source"),
        "resistance_target": "",
        "resistance_status": "unknown",
        "final_reporting_tier": "Fallback-BestAvailable",
        "fallback_rank": 1,
        "fallback_primary": "True",
        "final_reporting_confidence": "low",
        "patient_decision_evidence_mode": "mngs_only_or_unforwarded",
        "candidate_decision_evidence_mode": "mngs_only",
        "positive_or_contextual_sources": "mngs",
    }


def run(release_root: Path, output_dir: Path, policy_path: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty directory: {output_dir}")
    policy = read_json(policy_path)
    reporting_path = (
        release_root / "10_reporting_v3_route_a" / "possible_pathogen_decisions.csv"
    )
    analytical_root = release_root / "02_analytical_scorer" / "patient_outputs"
    if not reporting_path.is_file() or not analytical_root.is_dir():
        raise FileNotFoundError("Release is missing reporting or analytical scorer artifacts")

    rows = read_csv(reporting_path)
    if not rows:
        raise ValueError("Reporting decisions are empty")
    base_fields = list(rows[0])
    extra_fields = [
        "primary_candidate",
        "primary_selection_source",
        "primary_selection_reason",
        "primary_rank_vector",
        "clinical_certainty_changed",
    ]
    fields = base_fields + [field for field in extra_fields if field not in base_fields]
    for row in rows:
        row.update({
            "primary_candidate": "False",
            "primary_selection_source": "",
            "primary_selection_reason": "",
            "primary_rank_vector": "",
            "clinical_certainty_changed": "False",
        })

    rows_by_patient: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        rows_by_patient[int(row["patient_id"])].append(row)

    analytical_paths = sorted(
        analytical_root.glob("NGS_patient_*_test_aware_deterministic_shadow.json"),
        key=patient_number,
    )
    if not analytical_paths:
        raise ValueError("No analytical patient outputs found")

    primary_rows: list[dict[str, Any]] = []
    ranking_audit: list[dict[str, Any]] = []
    patient_payloads: list[dict[str, Any]] = []
    selected_existing_count = 0
    fallback_count = 0
    no_candidate_count = 0

    for path in analytical_paths:
        patient_id = patient_number(path)
        payload = read_json(path)
        patient_rows = rows_by_patient.get(patient_id, [])
        selected_rows = [
            row
            for row in patient_rows
            if str(row.get("selected_for_complete_report") or "").casefold() in TRUTHY
        ]
        selected_source = ""
        selection_reason = ""
        rank_vector: list[Any] = []
        chosen: dict[str, Any] | None = None

        if selected_rows:
            chosen = min(selected_rows, key=lambda row: selected_row_rank(row, policy))
            rank_vector = list(selected_row_rank(chosen, policy)[:-1])
            selected_source = "existing_complete_report"
            selection_reason = "Highest existing report tier, then generic ranking features."
            selected_existing_count += 1
        else:
            ranked: list[tuple[tuple[Any, ...], dict[str, Any], dict[str, Any]]] = []
            for candidate in payload.get("patient_organism_decisions") or []:
                features = fallback_features(candidate, policy)
                ranking_audit.append({
                    "patient_id": patient_id,
                    "organism_name": candidate.get("organism_name"),
                    "eligible": features is not None,
                    "selected": False,
                    "clinical_decision": candidate.get("decision"),
                    "taxonomy_family": candidate.get("taxonomy_family"),
                    "rank_vector": (
                        json.dumps(features["rank_vector"], ensure_ascii=False)
                        if features else ""
                    ),
                    "exclusion_reason": "" if features else "no_positive_accepted_qc_test",
                })
                if features is not None:
                    ranked.append((features["sort_key"], features, candidate))
            if ranked:
                ranked.sort(key=lambda item: item[0])
                _, features, candidate = ranked[0]
                rank_vector = features["rank_vector"]
                existing = next(
                    (
                        row
                        for row in patient_rows
                        if canonical_name(row.get("organism_name"))
                        == canonical_name(candidate.get("organism_name"))
                    ),
                    None,
                )
                fallback = synthetic_row(candidate, features)
                fallback.update({
                    "primary_candidate": "False",
                    "primary_selection_source": "",
                    "primary_selection_reason": "",
                    "primary_rank_vector": "",
                    "clinical_certainty_changed": "False",
                })
                if existing is None:
                    rows.append(fallback)
                    rows_by_patient[patient_id].append(fallback)
                    chosen = fallback
                else:
                    existing.update(fallback)
                    chosen = existing
                selected_source = "mandatory_fallback_from_analytical_candidates"
                selection_reason = (
                    "No Picked/Possible report row; selected by generic clinical, "
                    "specimen, taxonomy, per-test rank, QC, and reproducibility order."
                )
                fallback_count += 1
                for audit in ranking_audit:
                    if (
                        audit["patient_id"] == patient_id
                        and canonical_name(audit["organism_name"])
                        == canonical_name(candidate.get("organism_name"))
                    ):
                        audit["selected"] = True
                        break
            else:
                no_candidate_count += 1

        if chosen is not None:
            chosen["primary_candidate"] = "True"
            chosen["primary_selection_source"] = selected_source
            chosen["primary_selection_reason"] = selection_reason
            chosen["primary_rank_vector"] = json.dumps(rank_vector, ensure_ascii=False)
            chosen["clinical_certainty_changed"] = "False"
            primary_rows.append({field: chosen.get(field, "") for field in fields})

        patient_payloads.append({
            "patient_id": patient_id,
            "primary_candidate": (
                {field: chosen.get(field, "") for field in fields}
                if chosen is not None else None
            ),
            "primary_candidate_required": chosen is not None,
            "selection_source": selected_source or "no_valid_positive_read_candidate",
            "clinical_certainty_changed": False,
            "source_analytical_path": str(path.resolve()),
        })

    if len(primary_rows) + no_candidate_count != len(analytical_paths):
        raise AssertionError("Primary-candidate accounting mismatch")
    if any(
        sum(
            str(row.get("primary_candidate") or "").casefold() in TRUTHY
            for row in rows_by_patient.get(patient_number(path), [])
        ) > 1
        for path in analytical_paths
    ):
        raise AssertionError("A patient received more than one primary candidate")

    output_dir.mkdir(parents=True, exist_ok=True)
    patient_dir = output_dir / "patient_outputs"
    patient_dir.mkdir()
    for payload in patient_payloads:
        write_json(
            patient_dir / f"NGS_patient_{payload['patient_id']}_mandatory_primary.json",
            payload,
        )
    rows.sort(key=lambda row: (int(row["patient_id"]), canonical_name(row["organism_name"])))
    primary_rows.sort(key=lambda row: int(row["patient_id"]))
    complete_rows = [
        row
        for row in rows
        if str(row.get("selected_for_complete_report") or "").casefold() in TRUTHY
    ]
    write_csv(output_dir / "possible_pathogen_decisions.csv", fields, rows)
    write_csv(output_dir / "doctor_visible_report.csv", fields, complete_rows)
    write_csv(output_dir / "primary_candidates.csv", fields, primary_rows)
    write_csv(
        output_dir / "fallback_ranking_audit.csv",
        [
            "patient_id", "organism_name", "eligible", "selected",
            "clinical_decision", "taxonomy_family", "rank_vector",
            "exclusion_reason",
        ],
        ranking_audit,
    )
    summary = {
        "schema_version": "mandatory_primary_reporting_shadow.v1",
        "policy_id": policy["policy_id"],
        "answer_blind": True,
        "patient_count": len(analytical_paths),
        "primary_candidate_count": len(primary_rows),
        "existing_complete_report_primary_count": selected_existing_count,
        "mandatory_fallback_primary_count": fallback_count,
        "no_valid_positive_read_candidate_count": no_candidate_count,
        "doctor_visible_report_count": len(complete_rows),
        "doctor_visible_tier_counts": dict(
            Counter(row.get("final_reporting_tier") or "" for row in complete_rows)
        ),
        "source_release_root": str(release_root.resolve()),
        "source_reporting_path": str(reporting_path.resolve()),
        "source_reporting_sha256": sha256_file(reporting_path),
        "policy_path": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "constraints": policy.get("constraints") or [],
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    print(json.dumps(
        run(args.release_root, args.output_dir, args.policy),
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()

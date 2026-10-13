"""Run the finalized legacy-compatible rules on two 11-patient cohorts.

The original manually chosen DNA file does not exist for these patients.  This
adapter therefore selects one representative DNA test per infection event by
source order, runs the locked legacy scorer, and applies the same taxonomy,
history, HAP/device Priority, aspiration Picked, and four-label contract as the
KH finalized legacy release.  Answers are loaded only after decisions freeze.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools import deterministic_mngs_max_scorer as legacy
from tools.apply_healthcare_device_opportunist_history_route import (
    run as run_healthcare_device_route,
)
from tools.apply_legacy_history_picked_promotion import run as run_picked_gate
from tools.apply_test_aware_history_route_shadow import run as run_history_route
from tools.evaluate_multi_assay_candidate_entry import evaluate_set, load_answer_rows
from tools.legacy_primary_selection import (
    build_primary_selections,
    complete_report_rows,
    doctor_report_rows,
)
from tools.legacy_taxonomy_snapshot import materialize_locked_overlays
from tools.pathogen_normalization import canonical_key
from tools.run_kh_ablation_d_history_only import (
    answer_status,
    bridge_candidate,
    picked_keys,
    read_json,
    sha256_file,
    write_csv,
    write_json,
)
from tools.run_legacy_scorer_on_multi_assay_shadow import (
    RANK_RULES,
    SOURCE_CATEGORY,
    configure_taxonomy_mode,
    empirical_reads_percentile,
    possibility,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "rules/two_hospital_22_legacy_finalized_compat_v1.json"
DEFAULT_OUTPUT = (
    ROOT / "outputs/runs/2026-10-07_two_hospital_22_legacy_finalized_compat_v1"
)
SCHEMA_VERSION = "two_hospital_22_legacy_finalized_compat.v1"


def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot identify patient: {path}")
    return int(match.group(1))


def as_int(value: Any, default: int = 999) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def selected_rank(tags: list[Any]) -> int:
    normalized = {
        str(tag).replace("Code_RK_NTC=RK00", "Code_RK_NTC=RK0")
        for tag in tags
    }
    ntc = "Code_NTC=00" in normalized
    rk0 = "Code_RK_NTC=RK0" in normalized
    rk8 = "Code_RK_NTC=RK8" in normalized
    if ntc and rk0:
        return 1
    if ntc and rk8:
        return 2
    if rk0:
        return 3
    if rk8:
        return 4
    if ntc:
        return 5
    raise ValueError(f"Selected observation has no recognized selection tag: {tags}")


def representative_dna_tests(raw: dict[str, Any]) -> dict[str, dict[str, Any]]:
    positive_selected_seq_ids = {
        str(measurement.get("seq_id") or "")
        for organism in raw.get("organisms") or []
        for measurement in organism.get("measurements") or []
        if bool(measurement.get("selected"))
        and float(measurement.get("reads") or 0) > 0
        and str(measurement.get("seq_id") or "")
    }
    selected: dict[str, dict[str, Any]] = {}
    for case in raw.get("cases") or []:
        case_id = str(case.get("case_review_id") or case.get("specimen_code") or "")
        dna = [
            test
            for test in case.get("tests") or []
            if str(test.get("nucleic_type") or "").upper() == "DNA"
        ]
        if not dna:
            continue
        dna.sort(
            key=lambda test: (
                str(test.get("seq_id") or "")
                not in positive_selected_seq_ids,
                as_int(test.get("test_number")),
                str(test.get("seq_id") or ""),
            )
        )
        selected[case_id] = dna[0]
    return selected


def build_selected_dna_payload(
    raw: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    chosen = representative_dna_tests(raw)
    case_lookup = {
        str(case.get("case_review_id") or case.get("specimen_code") or ""): case
        for case in raw.get("cases") or []
    }
    records: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    for case_id, test in sorted(chosen.items()):
        seq_id = str(test.get("seq_id") or "")
        candidates: list[dict[str, Any]] = []
        for organism in raw.get("organisms") or []:
            name = str(organism.get("name") or "").strip()
            if not name:
                continue
            for measurement in organism.get("measurements") or []:
                measurement_case_id = str(
                    measurement.get("case_review_id") or ""
                )
                # Some upstream exports omit case_review_id on measurements
                # even though seq_id still identifies the selected test.  A
                # blank ID must therefore be treated as unknown, not as a
                # contradictory linkage.
                if measurement_case_id and measurement_case_id != case_id:
                    continue
                if str(measurement.get("seq_id") or "") != seq_id:
                    continue
                reads = float(measurement.get("reads") or 0)
                if not bool(measurement.get("selected")) or reads <= 0:
                    continue
                rank = selected_rank(list(measurement.get("selection_tags") or []))
                candidates.append(
                    {
                        "organism_name": name,
                        "source_category": SOURCE_CATEGORY.get(
                            str(organism.get("category") or "").lower(), "Unknown"
                        ),
                        "reads": reads,
                        "sec_hit": reads,
                        "selection_reasons": list(
                            measurement.get("selection_tags") or []
                        ),
                        "ranking": {
                            "rank_priority": rank,
                            "rank_rule": RANK_RULES[rank],
                            "reads_percentile": None,
                            "possibility_level": possibility(rank),
                        },
                    }
                )
        read_values = [float(item["reads"]) for item in candidates]
        for candidate in candidates:
            candidate["ranking"]["reads_percentile"] = empirical_reads_percentile(
                float(candidate["reads"]), read_values
            )
        candidates.sort(
            key=lambda item: (
                as_int(item["ranking"]["rank_priority"]),
                -float(item["ranking"]["reads_percentile"]),
                -float(item["reads"]),
                canonical_key(item["organism_name"]),
            )
        )
        case = case_lookup[case_id]
        records.append(
            {
                "specimen_code": case.get("specimen_code") or case_id,
                "collected_time": case.get("collected_time"),
                "specimen_site": case.get("specimen_site"),
                "seq_id": seq_id,
                "nucleic_type": "DNA",
                "status": "selected_dna_emulation_first_source_dna_test",
                "candidates": candidates,
            }
        )
        audit_rows.append(
            {
                "patient_id": raw.get("patient_id"),
                "case_review_id": case_id,
                "selected_seq_id": seq_id,
                "selected_test_number": test.get("test_number"),
                "selected_candidate_count": len(candidates),
                "available_dna_test_count": sum(
                    str(item.get("nucleic_type") or "").upper() == "DNA"
                    for item in case.get("tests") or []
                ),
                "available_rna_test_count": sum(
                    str(item.get("nucleic_type") or "").upper() == "RNA"
                    for item in case.get("tests") or []
                ),
            }
        )
    payload = {
        "patient_id": str(raw.get("patient_id")),
        "records": records,
        "ranking_metadata": {
            "source": "two_hospital_per_test_json",
            "adapter": SCHEMA_VERSION,
            "selection_rule": "one lowest-test-number DNA test per event",
                "empty_test_fallback_rule": (
                    "prefer the earliest DNA test containing at least one "
                    "selected positive-read measurement"
                ),
            "selected_observation_rule": "selected=true and reads>0",
            "reads_percentile": "empirical CDF within representative test",
            "no_cross_test_read_sum": True,
            "rna_excluded": True,
        },
    }
    return payload, {"events": audit_rows}


def level_rank(value: Any) -> int:
    match = re.search(r"([1-5])", str(value or ""))
    return int(match.group(1)) if match else 99


def baseline_decision(candidate: dict[str, Any], picked: set[str]) -> str:
    if canonical_key(candidate.get("organism_name")) in picked:
        return "picked_shadow"
    rank = level_rank(candidate.get("integrated_causative_level"))
    if rank <= 3:
        return "review_high_priority"
    if rank == 4:
        return "review_context_needed"
    return "review_low_specificity"


def score_and_bridge_cohort(
    cohort: dict[str, Any], output_root: Path
) -> tuple[Path, dict[str, Any]]:
    cohort_id = str(cohort["cohort_id"])
    source_root = resolve_path(cohort["source_root"])
    patient_root = resolve_path(cohort["patient_root"])
    expected = {int(value) for value in cohort["expected_patient_ids"]}
    all_source_paths = sorted(
        source_root.glob("NGS_patient_*_all_mNGS_DNA_RNA_RK_NTC_by_species.json"),
        key=patient_number,
    )
    observed_all = {patient_number(path) for path in all_source_paths}
    allow_source_superset = bool(cohort.get("allow_source_superset", False))
    if allow_source_superset:
        source_paths = [
            path for path in all_source_paths if patient_number(path) in expected
        ]
        observed = {patient_number(path) for path in source_paths}
    else:
        source_paths = all_source_paths
        observed = observed_all
    if observed != expected:
        raise ValueError(
            f"{cohort_id} patient mismatch: observed={sorted(observed)}, "
            f"expected={sorted(expected)}, source_all={sorted(observed_all)}"
        )

    scorer_root = output_root / "01_legacy_selected_dna_scorer" / cohort_id
    bridge_root = output_root / "02_legacy_route_bridge" / cohort_id
    ranked_root = scorer_root / "legacy_ranked_inputs"
    scored_root = scorer_root / "patient_outputs"
    bridge_output = bridge_root / "patient_outputs"
    for path in (ranked_root, scored_root, bridge_output):
        path.mkdir(parents=True, exist_ok=True)

    score_rows: list[dict[str, Any]] = []
    bridge_rows: list[dict[str, Any]] = []
    source_audit: list[dict[str, Any]] = []
    source_hashes: dict[str, str] = {}
    summary_hashes: dict[str, str] = {}
    counts: Counter[str] = Counter()
    summary_contract_rows: list[dict[str, Any]] = []
    for source_path in source_paths:
        patient = patient_number(source_path)
        raw = read_json(source_path)
        ranked, audit = build_selected_dna_payload(raw)
        patient_dir = patient_root / f"NGS_patient_{patient}_json"
        summary_path = patient_dir / "summary_outputs" / (
            f"NGS_patient_{patient}_{cohort['summary_suffix']}.json"
        )
        if not summary_path.is_file():
            raise FileNotFoundError(summary_path)
        final_summary = read_json(summary_path)
        required_summary_fields = list(cohort.get("required_summary_fields") or [])
        missing_summary_fields = [
            field for field in required_summary_fields if field not in final_summary
        ]
        if missing_summary_fields:
            raise ValueError(
                f"{cohort_id} patient {patient} summary is incompatible: "
                f"missing {missing_summary_fields} in {summary_path}"
            )
        hospital_evidence = final_summary.get("hospital_organism_evidence")
        if hospital_evidence is not None and not isinstance(hospital_evidence, list):
            raise ValueError(
                f"{cohort_id} patient {patient} hospital_organism_evidence "
                f"must be a list: {summary_path}"
            )
        summary_contract_rows.append(
            {
                "cohort_id": cohort_id,
                "patient_id": patient,
                "summary_path": str(summary_path.resolve()),
                "summary_sha256": sha256_file(summary_path),
                "summary_schema_role": cohort.get("summary_schema_role", ""),
                "required_summary_fields": "|".join(required_summary_fields),
                "hospital_organism_evidence_count": len(hospital_evidence or []),
                "schema_contract_passed": True,
            }
        )
        score = legacy.score_payload(
            patient_dir=patient_dir,
            ranked_mngs=ranked,
            final_summary=final_summary,
        )
        score.update(
            {
                "schema_version": SCHEMA_VERSION,
                "answer_blind": True,
                "cohort_id": cohort_id,
                "selected_dna_emulation": ranked["ranking_metadata"],
                "source_files": {
                    "raw_mngs": str(source_path.resolve()),
                    "raw_mngs_sha256": sha256_file(source_path),
                    "hospital_summary": str(summary_path.resolve()),
                    "hospital_summary_sha256": sha256_file(summary_path),
                },
            }
        )
        ranked_path = ranked_root / f"NGS_patient_{patient}_selected_dna_emulation.json"
        score_path = scored_root / f"NGS_patient_{patient}_legacy_scorer.json"
        write_json(ranked_path, ranked)
        write_json(score_path, score)
        source_hashes[str(source_path.resolve())] = sha256_file(source_path)
        summary_hashes[str(summary_path.resolve())] = sha256_file(summary_path)
        source_audit.extend(audit["events"])

        picked = picked_keys(score)
        patient_candidates: list[dict[str, Any]] = []
        hospital_only: list[dict[str, Any]] = []
        for candidate in score.get("pathogen_candidates") or []:
            if not isinstance(candidate, dict):
                continue
            decision = baseline_decision(candidate, picked)
            bridged = bridge_candidate(
                candidate,
                decision=decision,
                review_item=None,
                taxonomy_mode="locked_taxonomy_v2_6_external_compat",
            )
            target = (
                hospital_only
                if candidate.get("evidence_source") == "hospital_only"
                else patient_candidates
            )
            target.append(bridged)
            score_rows.append(
                {
                    "cohort_id": cohort_id,
                    "patient_id": patient,
                    "organism_name": candidate.get("organism_name"),
                    "integrated_level": candidate.get("integrated_causative_level"),
                    "legacy_picked": canonical_key(candidate.get("organism_name")) in picked,
                    "baseline_decision": decision,
                    "evidence_source": candidate.get("evidence_source"),
                }
            )
            bridge_rows.append(
                {
                    "cohort_id": cohort_id,
                    "patient_id": patient,
                    "organism_name": candidate.get("organism_name"),
                    "baseline_decision": decision,
                    "taxonomy_family": bridged.get("taxonomy_family"),
                }
            )
            counts[f"baseline:{decision}"] += 1
            counts["candidates"] += 1
        bridge_payload = {
            "schema_version": SCHEMA_VERSION,
            "answer_blind": True,
            "patient_id": str(patient),
            "patient_organism_decisions": patient_candidates,
            "hospital_only_decisions": hospital_only,
            "source": {
                "legacy_score": str(score_path.resolve()),
                "legacy_score_sha256": sha256_file(score_path),
            },
            "constraints": [
                "One representative DNA test per event; no RNA or repeat axis.",
                "Non-Picked legacy levels are deterministically mapped to Priority, Context, or Audit because no frozen patient-specific Luna review exists.",
            ],
        }
        write_json(
            bridge_output
            / f"NGS_patient_{patient}_test_aware_deterministic_shadow.json",
            bridge_payload,
        )
        counts["patients"] += 1

    write_csv(scorer_root / "legacy_candidate_decisions.csv", score_rows)
    write_csv(bridge_root / "baseline_candidate_decisions.csv", bridge_rows)
    write_csv(scorer_root / "selected_dna_event_audit.csv", source_audit)
    write_csv(scorer_root / "source_summary_contract_audit.csv", summary_contract_rows)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "cohort_id": cohort_id,
        "answer_blind": True,
        "clinical_data_mode": cohort["clinical_data_mode"],
        "counts": dict(sorted(counts.items())),
        "source_hashes": source_hashes,
        "summary_hashes": summary_hashes,
        "source_root": str(source_root.resolve()),
        "source_superset_allowed": allow_source_superset,
        "source_patient_count_before_filter": len(observed_all),
        "patient_root": str(patient_root.resolve()),
        "summary_suffix": cohort["summary_suffix"],
        "summary_schema_role": cohort.get("summary_schema_role", ""),
        "required_summary_fields": list(cohort.get("required_summary_fields") or []),
        "summary_contract_passed": all(
            bool(row["schema_contract_passed"]) for row in summary_contract_rows
        ),
        "structured_hospital_evidence_count": sum(
            int(row["hospital_organism_evidence_count"])
            for row in summary_contract_rows
        ),
    }
    write_json(scorer_root / "summary.json", summary)
    write_json(bridge_root / "summary.json", summary)
    return bridge_root, summary


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def grouped(
    rows: list[dict[str, Any]], allowed: set[str], patients: set[int] | None = None
) -> dict[int, list[str]]:
    result: dict[int, list[str]] = {}
    for row in rows:
        patient = int(row["patient_id"])
        if patients is not None and patient not in patients:
            continue
        if row["final_label"] not in allowed:
            continue
        name = str(row.get("organism_name") or "")
        if name and name not in result.setdefault(patient, []):
            result[patient].append(name)
    return result


def run(contract_path: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    contract = read_json(contract_path)
    parent_path = resolve_path(contract["parent_release_contract"])
    parent = read_json(parent_path)
    paths = parent["paths"]
    locked = parent["locked_sha256"]
    history_policy = resolve_path(paths["history_policy"])
    healthcare_policy = resolve_path(paths["healthcare_device_policy"])
    picked_policy = resolve_path(paths["picked_policy"])
    for path, expected, label in (
        (history_policy, locked["history_policy"], "history policy"),
        (
            healthcare_policy,
            locked["healthcare_device_policy"],
            "healthcare-device policy",
        ),
        (picked_policy, locked["picked_policy"], "Picked policy"),
    ):
        actual = sha256_file(path)
        if actual != expected:
            raise ValueError(f"Locked {label} changed: {actual} != {expected}")

    output_root.mkdir(parents=True)
    locked_taxonomy_overlays = materialize_locked_overlays(
        repo_root=ROOT,
        output_root=output_root,
        contract=parent,
    )
    taxonomy = configure_taxonomy_mode(
        "current-overlays", locked_taxonomy_overlays
    )
    overlay_hashes = {
        Path(path).name: digest
        for path, digest in taxonomy["loaded_overlay_sha256"].items()
    }
    expected_overlays = {
        "organism_taxonomy_gap_resolution_v1.json": locked[
            "taxonomy_gap_resolution_overlay"
        ],
        "organism_taxonomy_auto_v1.json": locked["taxonomy_auto_overlay"],
        "organism_taxonomy_reviewed_v1.json": locked[
            "taxonomy_reviewed_overlay"
        ],
    }
    if taxonomy["base_rules_sha256"] != locked["taxonomy_base_rules"]:
        raise ValueError("Locked taxonomy base changed")
    if overlay_hashes != expected_overlays:
        raise ValueError("Locked taxonomy overlays changed")

    final_rows: list[dict[str, Any]] = []
    cohort_stage_summaries: dict[str, Any] = {}
    cohort_patients: dict[str, set[int]] = {}
    for cohort in contract["cohorts"]:
        cohort_id = str(cohort["cohort_id"])
        cohort_patients[cohort_id] = {
            int(value) for value in cohort["expected_patient_ids"]
        }
        bridge_root, scorer_summary = score_and_bridge_cohort(cohort, output_root)
        patient_root = resolve_path(cohort["patient_root"])
        history_root = output_root / "03_history_route_v4c" / cohort_id
        history_summary = run_history_route(
            bridge_root,
            patient_root,
            history_root,
            policy_path=history_policy,
            phenotype_root=None,
            patients=cohort_patients[cohort_id],
        )
        if history_summary["resolved_policy_sha256"] != locked["resolved_history_policy"]:
            raise ValueError("Locked resolved history policy changed")
        healthcare_root = output_root / "04_healthcare_device_priority_route" / cohort_id
        healthcare_summary = run_healthcare_device_route(
            history_root, healthcare_root, healthcare_policy
        )
        picked_root = output_root / "05_aspiration_picked_gate" / cohort_id
        picked_summary = run_picked_gate(healthcare_root, picked_root, picked_policy)
        rows = read_csv(picked_root / "promotion_decisions.csv")
        for row in rows:
            row["cohort_id"] = cohort_id
            row["final_label"] = parent["output_labels"][row["final_decision"]]
        final_rows.extend(rows)
        cohort_stage_summaries[cohort_id] = {
            "scorer": scorer_summary["counts"],
            "history": history_summary["counts"],
            "healthcare_device": healthcare_summary["counts"],
            "picked_gate": picked_summary["counts"],
        }

    expected_ids = set().union(*cohort_patients.values())
    expected_patient_count = int(
        contract.get("expected_patient_count", len(expected_ids))
    )
    observed_ids = {int(row["patient_id"]) for row in final_rows}
    unexpected_ids = observed_ids - expected_ids
    processed_ids = set().union(
        *(
            patients
            for cohort_id, patients in cohort_patients.items()
            if cohort_stage_summaries[cohort_id]["scorer"].get("patients")
            == len(patients)
        )
    )
    if (
        unexpected_ids
        or processed_ids != expected_ids
        or len(expected_ids) != expected_patient_count
    ):
        raise ValueError(
            f"Final patient mismatch: candidates={sorted(observed_ids)}, "
            f"processed={sorted(processed_ids)}, expected={sorted(expected_ids)}, "
            f"expected_patient_count={expected_patient_count}"
        )
    if any(row["final_label"] == "Possible" for row in final_rows):
        raise ValueError("Legacy finalized compatibility run must not emit Possible")

    final_root = output_root / "06_frozen_decisions"
    final_path = final_root / "final_decisions.csv"
    write_csv(final_path, final_rows)
    cohort_by_patient = {
        patient: cohort_id
        for cohort_id, patients in cohort_patients.items()
        for patient in patients
    }
    primary_rows = build_primary_selections(
        final_rows, expected_ids, cohort_by_patient
    )
    primary_path = final_root / "patient_primary_output.csv"
    write_csv(primary_path, primary_rows)
    complete_fallback_rows = complete_report_rows(final_rows, primary_rows)
    doctor_report = doctor_report_rows(final_rows, primary_rows)
    doctor_report_path = final_root / "doctor_visible_report.csv"
    write_csv(doctor_report_path, doctor_report)
    label_counts = dict(sorted(Counter(row["final_label"] for row in final_rows).items()))
    patient_status_rows: list[dict[str, Any]] = []
    for cohort_id, patients in cohort_patients.items():
        for patient in sorted(patients):
            rows = [row for row in final_rows if int(row["patient_id"]) == patient]
            row_counts = Counter(row["final_label"] for row in rows)
            patient_status_rows.append(
                {
                    "cohort_id": cohort_id,
                    "patient_id": patient,
                    "processing_status": (
                        "complete_with_candidates" if rows else "complete_zero_candidates"
                    ),
                    "candidate_count": len(rows),
                    "picked_count": row_counts.get("Picked", 0),
                    "priority_count": row_counts.get("Priority", 0),
                    "context_count": row_counts.get("Context", 0),
                    "audit_count": row_counts.get("Audit", 0),
                    "primary_output": next(
                        (
                            item["organism_name"]
                            for item in primary_rows
                            if int(item["patient_id"]) == patient
                        ),
                        "",
                    ),
                    "primary_output_role": next(
                        (
                            item["primary_selection_role"]
                            for item in primary_rows
                            if int(item["patient_id"]) == patient
                        ),
                        "",
                    ),
                }
            )
    write_csv(final_root / "patient_processing_status.csv", patient_status_rows)
    zero_candidate_ids = sorted(expected_ids - observed_ids)
    freeze = {
        "schema_version": SCHEMA_VERSION,
        "release_id": contract["release_id"],
        "status": "decisions_frozen_before_answers",
        "answers_loaded_during_decision_generation": False,
        "independent_validation_claim": contract["independent_validation_claim"],
        "parent_release_contract": str(parent_path.resolve()),
        "parent_release_contract_sha256": sha256_file(parent_path),
        "cohort_contract": str(contract_path.resolve()),
        "cohort_contract_sha256": sha256_file(contract_path),
        "patient_count": len(processed_ids),
        "patients_with_candidates": len(observed_ids),
        "zero_candidate_patient_ids": zero_candidate_ids,
        "candidate_count": len(final_rows),
        "label_counts": label_counts,
        "primary_selection_counts": dict(
            sorted(
                Counter(
                    row["primary_selection_role"] for row in primary_rows
                ).items()
            )
        ),
        "cohort_stage_summaries": cohort_stage_summaries,
        "selected_dna_emulation": contract["selected_dna_emulation"],
        "constraints": contract["constraints"],
        "final_decisions_sha256": sha256_file(final_path),
        "patient_primary_output_sha256": sha256_file(primary_path),
        "doctor_visible_report_sha256": sha256_file(doctor_report_path),
        "runner_sha256": sha256_file(Path(__file__)),
    }
    freeze_path = final_root / "decision_freeze.json"
    write_json(freeze_path, freeze)

    # Post-freeze evaluation starts here. No answer file is opened above.
    answer_sets: dict[str, dict[int, list[str]]] = {}
    combined_answers: dict[int, list[str]] = {}
    answer_hashes: dict[str, str] = {}
    for cohort in contract["cohorts"]:
        cohort_id = str(cohort["cohort_id"])
        answer_path = resolve_path(cohort["postfreeze_answer_path"])
        answers = load_answer_rows(answer_path)
        if set(answers) - cohort_patients[cohort_id]:
            raise ValueError(f"{cohort_id} answer IDs are outside the cohort")
        if set(answers) & set(combined_answers):
            raise ValueError("Answer contracts contain overlapping patient IDs")
        answer_sets[cohort_id] = answers
        combined_answers.update(answers)
        answer_hashes[cohort_id] = sha256_file(answer_path)

    metric_labels = {
        "strict_picked": {"Picked"},
        "picked_or_priority": {"Picked", "Priority"},
        "complete_visible": {"Picked", "Priority", "Context"},
    }
    combined_metric_key = str(
        contract.get("combined_metric_key", "combined_labeled19")
    )
    metrics: dict[str, Any] = {combined_metric_key: {}}
    for endpoint, labels in metric_labels.items():
        metrics[combined_metric_key][endpoint] = evaluate_set(
            grouped(final_rows, labels), combined_answers
        )
    metrics[combined_metric_key]["picked_plus_mandatory_fallback"] = evaluate_set(
        grouped(complete_fallback_rows, {"Picked", "Priority", "Context", "Audit"}),
        combined_answers,
    )
    for cohort_id, answers in answer_sets.items():
        metrics[cohort_id] = {}
        for endpoint, labels in metric_labels.items():
            metrics[cohort_id][endpoint] = evaluate_set(
                grouped(final_rows, labels, cohort_patients[cohort_id]), answers
            )
        metrics[cohort_id]["picked_plus_mandatory_fallback"] = evaluate_set(
            grouped(
                complete_fallback_rows,
                {"Picked", "Priority", "Context", "Audit"},
                cohort_patients[cohort_id],
            ),
            answers,
        )

    audit_rows: list[dict[str, Any]] = []
    for row in final_rows:
        if row["final_label"] == "Audit":
            continue
        patient = int(row["patient_id"])
        audit_rows.append(
            {
                "cohort_id": row["cohort_id"],
                "patient_id": patient,
                "organism_name": row["organism_name"],
                "final_label": row["final_label"],
                "postfreeze_answer_status": answer_status(
                    patient, row["organism_name"], combined_answers
                ),
            }
        )
    evaluation_root = output_root / "07_postfreeze_evaluation"
    write_csv(evaluation_root / "visible_candidate_status.csv", audit_rows)
    evaluation = {
        "schema_version": SCHEMA_VERSION,
        "scope": "Post-freeze compatibility evaluation; not untouched external validation.",
        "decision_freeze_sha256": sha256_file(freeze_path),
        "labeled_patient_count": len(combined_answers),
        "answer_count": sum(len(value) for value in combined_answers.values()),
        "answer_contract_sha256": answer_hashes,
        "metrics": metrics,
    }
    write_json(evaluation_root / "metrics.json", evaluation)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "release_id": contract["release_id"],
        "status": "complete_frozen_compatibility_run",
        "decision_generation_answer_blind": True,
        "benchmark_evaluation_is_postfreeze": True,
        "independent_validation_claim": False,
        "patient_count": len(processed_ids),
        "patients_with_candidates": len(observed_ids),
        "zero_candidate_patient_ids": zero_candidate_ids,
        "candidate_count": len(final_rows),
        "label_counts": label_counts,
        "primary_selection_counts": dict(
            sorted(
                Counter(
                    row["primary_selection_role"] for row in primary_rows
                ).items()
            )
        ),
        "metrics": metrics,
        "decision_freeze": str(freeze_path.resolve()),
        "decision_freeze_sha256": sha256_file(freeze_path),
        "final_decisions": str(final_path.resolve()),
        "final_decisions_sha256": sha256_file(final_path),
        "patient_primary_output": str(primary_path.resolve()),
        "patient_primary_output_sha256": sha256_file(primary_path),
        "doctor_visible_report": str(doctor_report_path.resolve()),
        "doctor_visible_report_sha256": sha256_file(doctor_report_path),
    }
    write_json(output_root / "run_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run(args.contract, args.output), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

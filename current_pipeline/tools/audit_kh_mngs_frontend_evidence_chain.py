"""Audit the KH per-test mNGS front end from QC through scorer entry.

This audit is intentionally answer blind.  It checks structural invariants and
provenance only; physician benchmark answers are not read.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


DEFAULT_INVENTORY = Path(
    "outputs/runs/2026-09-23_KH_multi_assay_mngs_inventory_source_contract_v3"
)
DEFAULT_COMPACT = Path(
    "outputs/runs/2026-09-25_KH_compact_multi_assay_entry_taxonomy_auto_v1"
)
DEFAULT_SCORER = Path(
    "outputs/runs/2026-09-25_KH_test_aware_deterministic_taxonomy_auto_v1"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-30_KH_mngs_frontend_evidence_chain_audit_v1"
)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_float(value: Any) -> float | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def as_int(value: Any) -> int | None:
    number = as_float(value)
    return int(number) if number is not None else None


def as_bool(value: Any) -> bool:
    return str(value or "").strip().lower() == "true"


def key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (
        str(row.get("patient_id") or ""),
        str(row.get("case_review_id") or ""),
        str(row.get("organism_key") or ""),
    )


def signal_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return (*key(row), str(row.get("seq_id") or ""))


def count_values(rows: Iterable[dict[str, Any]], field: str) -> dict[str, int]:
    return dict(sorted(Counter(str(row.get(field) or "") for row in rows).items()))


def compact_profiles(compact_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    scorer_profiles: list[dict[str, Any]] = []
    clinical_profiles: list[dict[str, Any]] = []
    for path in sorted((compact_dir / "patient_packets").glob("*.json")):
        packet = read_json(path)
        patient_id = str(packet["patient_id"])
        for case in packet.get("cases") or []:
            case_id = str(case.get("case_review_id") or "")
            for disposition, target in (
                ("scorer_entry", scorer_profiles),
                ("clinical_review", clinical_profiles),
            ):
                field = "scorer_entries" if disposition == "scorer_entry" else "clinical_review"
                for profile in case.get(field) or []:
                    target.append(
                        {
                            **profile,
                            "patient_id": patient_id,
                            "case_review_id": case_id,
                            "disposition": disposition,
                        }
                    )
    return scorer_profiles, clinical_profiles


def audit(
    inventory_dir: Path,
    compact_dir: Path,
    scorer_dir: Path,
    output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    test_rows = read_csv(inventory_dir / "test_qc.csv")
    observation_rows = read_csv(inventory_dir / "observation_qc.csv")
    positive_rows = read_csv(inventory_dir / "positive_observations.csv")
    candidate_rows = read_csv(inventory_dir / "candidate_union.csv")
    compact_rows = read_csv(compact_dir / "all_decisions.csv")
    inventory_summary = read_json(inventory_dir / "inventory_summary.json")
    compact_summary = read_json(compact_dir / "summary.json")
    scorer_summary = read_json(scorer_dir / "summary.json")
    scorer_profiles, clinical_profiles = compact_profiles(compact_dir)

    issues: list[dict[str, Any]] = []

    def issue(scope: str, identifier: str, code: str, detail: str) -> None:
        issues.append(
            {"scope": scope, "identifier": identifier, "code": code, "detail": detail}
        )

    test_index: dict[tuple[str, str, str], dict[str, str]] = {}
    duplicate_test_keys: list[tuple[str, str, str]] = []
    for row in test_rows:
        test_key = (
            str(row["patient_id"]),
            str(row["case_review_id"]),
            str(row["seq_id"]),
        )
        if test_key in test_index:
            duplicate_test_keys.append(test_key)
        test_index[test_key] = row

    observations_by_test: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    observations_by_candidate: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in observation_rows:
        test_key = (
            str(row["patient_id"]),
            str(row["case_review_id"]),
            str(row["seq_id"]),
        )
        observations_by_test[test_key].append(row)
        observations_by_candidate[key(row)].append(row)

    positive_keys = {signal_key(row) for row in positive_rows}
    expected_positive_keys = {
        signal_key(row)
        for row in observation_rows
        if row.get("availability_status") == "evaluable"
        and (as_float(row.get("reads")) or 0.0) > 0
    }

    rank_mismatches = 0
    rpm_mismatches = 0
    positive_definition_mismatches = 0
    for test_key, rows in observations_by_test.items():
        test = test_index.get(test_key)
        if test is None:
            issue("observation", ":".join(test_key), "orphan_test", "No test_qc row")
            continue
        positives = [
            row
            for row in rows
            if row.get("availability_status") == "evaluable"
            and (as_float(row.get("reads")) or 0.0) > 0
        ]
        expected_rank = {
            str(row["organism_key"]): index + 1
            for index, row in enumerate(
                sorted(
                    positives,
                    key=lambda item: (
                        -(as_float(item.get("reads")) or 0.0),
                        str(item.get("organism_key") or ""),
                    ),
                )
            )
        }
        denominator = as_float(test.get("input_reads"))
        for row in rows:
            reads = as_float(row.get("reads")) or 0.0
            is_positive = row.get("availability_status") == "evaluable" and reads > 0
            if as_bool(row.get("observed_in_current_test")) != is_positive:
                positive_definition_mismatches += 1
            observed_rank = as_int(row.get("rank_by_reads_in_retained_universe"))
            wanted_rank = expected_rank.get(str(row.get("organism_key") or ""))
            if observed_rank != wanted_rank:
                rank_mismatches += 1
            observed_rpm = as_float(row.get("rpm_total"))
            wanted_rpm = reads / denominator * 1_000_000 if denominator and row.get("availability_status") == "evaluable" else None
            if wanted_rpm is None:
                if observed_rpm is not None:
                    rpm_mismatches += 1
            elif observed_rpm is None or not math.isclose(observed_rpm, wanted_rpm, rel_tol=1e-12, abs_tol=1e-12):
                rpm_mismatches += 1

    candidate_index = {key(row): row for row in candidate_rows}
    compact_index = {key(row): row for row in compact_rows}
    duplicate_candidate_keys = len(candidate_index) != len(candidate_rows)
    duplicate_compact_keys = len(compact_index) != len(compact_rows)
    candidate_mismatches = 0
    candidate_trace: list[dict[str, Any]] = []
    for candidate_key, row in candidate_index.items():
        observations = [
            item
            for item in observations_by_candidate.get(candidate_key, [])
            if item.get("availability_status") == "evaluable"
            and as_bool(item.get("test_eligible_for_union"))
            and (as_float(item.get("reads")) or 0.0) > 0
        ]
        selected = [item for item in observations if as_bool(item.get("selected"))]
        dna = [item for item in observations if item.get("nucleic_type") == "DNA"]
        rna = [item for item in observations if item.get("nucleic_type") == "RNA"]
        checks = {
            "detected_test_count": as_int(row.get("detected_test_count")) == len(observations),
            "selected_test_count": as_int(row.get("selected_test_count")) == len(selected),
            "selected_union": as_bool(row.get("selected_union")) == bool(selected),
            "cross_molecule_support": as_bool(row.get("cross_molecule_support")) == bool(dna and rna),
            "repeated_detection": as_bool(row.get("repeated_detection")) == (len(observations) >= 2),
            "no_reads_sum": not str(row.get("reads_sum_across_tests") or "").strip(),
        }
        if not all(checks.values()):
            candidate_mismatches += 1
            issue(
                "candidate",
                ":".join(candidate_key),
                "candidate_aggregation_mismatch",
                json.dumps(checks, ensure_ascii=False, sort_keys=True),
            )
        compact = compact_index.get(candidate_key)
        candidate_trace.append(
            {
                "patient_id": candidate_key[0],
                "case_review_id": candidate_key[1],
                "organism_key": candidate_key[2],
                "organism_name": row.get("organism_name"),
                "union_status": row.get("union_status"),
                "detected_test_count": len(observations),
                "selected_test_count": len(selected),
                "dna_positive_test_count": len(dna),
                "rna_positive_test_count": len(rna),
                "cross_molecule_support": bool(dna and rna),
                "repeated_detection": len(observations) >= 2,
                "best_rank_in_retained_universe": row.get("best_rank_in_retained_universe"),
                "max_reads_single_test": row.get("max_reads_single_test"),
                "compact_disposition": compact.get("disposition") if compact else "MISSING",
                "screening_tier": compact.get("screening_tier") if compact else "",
                "compact_selected_positive_test_count": compact.get("selected_positive_test_count") if compact else "",
                "trace_status": "PASS" if compact is not None and all(checks.values()) else "REVIEW",
            }
        )

    scorer_signals: list[dict[str, Any]] = []
    for profile in scorer_profiles:
        for signal in profile.get("selected_test_signals") or []:
            scorer_signals.append(
                {
                    "patient_id": profile["patient_id"],
                    "case_review_id": profile["case_review_id"],
                    "organism_key": profile["organism_key"],
                    **signal,
                }
            )
    scorer_signal_keys = {signal_key(row) for row in scorer_signals}
    scorer_signals_without_rpm = [
        row for row in scorer_signals if as_float(row.get("rpm_total")) is None
    ]
    selected_positive_keys = {
        signal_key(row)
        for row in positive_rows
        if as_bool(row.get("selected"))
    }

    source_checks: list[dict[str, Any]] = []
    source_hash_mismatches = 0
    for path in sorted((inventory_dir / "patients").glob("*.json")):
        packet = read_json(path)
        source = packet.get("source") or {}
        source_path = Path(str(source.get("path") or ""))
        exists = source_path.is_file()
        actual = sha256_file(source_path) if exists else None
        expected = source.get("sha256")
        match = exists and actual == expected
        source_hash_mismatches += not match
        source_checks.append(
            {
                "patient_id": packet.get("patient_id"),
                "inventory_packet": str(path.resolve()),
                "source_path": str(source_path),
                "source_exists": exists,
                "source_sha256_expected": expected,
                "source_sha256_actual": actual,
                "source_hash_match": match,
            }
        )

    compact_summary_hash = sha256_file(compact_dir / "summary.json")
    scorer_expected_hash = scorer_summary.get("source_compact_entry_summary_sha256")
    all_decision_scorer_entries = [
        row for row in compact_rows if row.get("disposition") == "scorer_entry"
    ]
    scorer_selected_count = sum(
        as_int(row.get("selected_positive_test_count")) or 0
        for row in all_decision_scorer_entries
    )

    assertions = {
        "patient_count_is_33": len({row["patient_id"] for row in test_rows}) == 33,
        "case_count_is_34": len({(row["patient_id"], row["case_review_id"]) for row in test_rows}) == 34,
        "test_count_is_117": len(test_rows) == 117,
        "test_keys_unique": not duplicate_test_keys,
        "nucleic_types_only_dna_rna": {row["nucleic_type"] for row in test_rows} <= {"DNA", "RNA"},
        "no_invalid_test_qc": not any(row.get("qc_status") == "invalid" for row in test_rows),
        "raw_unavailable_is_not_union_eligible": all(
            not as_bool(row.get("eligible_for_union"))
            for row in test_rows
            if row.get("qc_status") == "raw_unavailable"
        ),
        "source_patient_packets_hash_match": source_hash_mismatches == 0,
        "positive_csv_exactly_matches_evaluable_reads_gt_zero": positive_keys == expected_positive_keys,
        "observed_in_current_test_matches_reads_gt_zero": positive_definition_mismatches == 0,
        "case_manifest_only_rows_are_zero_read": all(
            (as_float(row.get("reads")) or 0.0) == 0.0
            for row in observation_rows
            if row.get("selection_status") == "case_manifest_only"
        ),
        "case_manifest_only_rows_do_not_enter_positive_csv": not any(
            signal_key(row) in positive_keys
            for row in observation_rows
            if row.get("selection_status") == "case_manifest_only"
        ),
        "rpm_formula_and_missingness_correct": rpm_mismatches == 0,
        "within_test_rank_recomputes_exactly": rank_mismatches == 0,
        "candidate_keys_unique": not duplicate_candidate_keys,
        "candidate_generation_recomputes_exactly": candidate_mismatches == 0,
        "compact_keys_unique": not duplicate_compact_keys,
        "candidate_to_compact_is_one_to_one": set(candidate_index) == set(compact_index),
        "scorer_entry_count_is_428": len(all_decision_scorer_entries) == 428,
        "scorer_entry_profiles_complete": len(scorer_profiles) == 428,
        "scorer_selected_signal_count_is_856": len(scorer_signals) == 856,
        "scorer_signal_count_matches_all_decisions": scorer_selected_count == len(scorer_signals),
        "all_scorer_signals_are_selected_positive_observations": scorer_signal_keys <= selected_positive_keys,
        "scorer_does_not_use_cross_test_read_sum": scorer_summary.get("read_handling", "").startswith("Per-test signals only"),
        "scorer_input_summary_hash_matches": compact_summary_hash == scorer_expected_hash,
        "scorer_expected_entry_count_matches": scorer_summary.get("expected_case_scorer_entries_from_input") == 428,
        "answer_blind_frontend": scorer_summary.get("answer_blind") is True and not scorer_summary.get("answer_keys_used"),
    }

    per_test_rows: list[dict[str, Any]] = []
    scorer_signal_counts = Counter(
        (str(row["patient_id"]), str(row["case_review_id"]), str(row["seq_id"]))
        for row in scorer_signals
    )
    for row in test_rows:
        test_key = (str(row["patient_id"]), str(row["case_review_id"]), str(row["seq_id"]))
        observations = observations_by_test.get(test_key, [])
        per_test_rows.append(
            {
                **row,
                "observation_count": len(observations),
                "positive_observation_count": sum(
                    item.get("availability_status") == "evaluable"
                    and (as_float(item.get("reads")) or 0.0) > 0
                    for item in observations
                ),
                "selected_positive_count": sum(
                    item.get("selection_status") == "selected" for item in observations
                ),
                "filtered_positive_count": sum(
                    item.get("selection_status") == "filtered" for item in observations
                ),
                "case_manifest_only_count": sum(
                    item.get("selection_status") == "case_manifest_only" for item in observations
                ),
                "scorer_signal_count": scorer_signal_counts[test_key],
                "audit_status": "PASS",
            }
        )

    summary = {
        "schema_version": "kh_mngs_frontend_evidence_chain_audit.v1",
        "answer_blind": True,
        "inputs": {
            "inventory_dir": str(inventory_dir.resolve()),
            "compact_dir": str(compact_dir.resolve()),
            "scorer_dir": str(scorer_dir.resolve()),
            "inventory_summary_sha256": sha256_file(inventory_dir / "inventory_summary.json"),
            "compact_summary_sha256": compact_summary_hash,
            "scorer_summary_sha256": sha256_file(scorer_dir / "summary.json"),
        },
        "counts": {
            "patients": len({row["patient_id"] for row in test_rows}),
            "cases": len({(row["patient_id"], row["case_review_id"]) for row in test_rows}),
            "tests": len(test_rows),
            "nucleic_types": count_values(test_rows, "nucleic_type"),
            "test_qc_status": count_values(test_rows, "qc_status"),
            "tests_missing_input_reads": sum(as_float(row.get("input_reads")) is None for row in test_rows),
            "tests_missing_report_time": sum(row.get("report_time_status") != "available" for row in test_rows),
            "observations": len(observation_rows),
            "positive_observations": len(positive_rows),
            "selection_status": count_values(observation_rows, "selection_status"),
            "candidate_union": len(candidate_rows),
            "candidate_union_status": count_values(candidate_rows, "union_status"),
            "compact_disposition": count_values(compact_rows, "disposition"),
            "scorer_entries": len(scorer_profiles),
            "scorer_selected_test_signals": len(scorer_signals),
            "scorer_signals_without_rpm": len(scorer_signals_without_rpm),
            "scorer_signals_without_rpm_by_molecule": count_values(
                scorer_signals_without_rpm, "nucleic_type"
            ),
            "clinical_review_profiles": len(clinical_profiles),
            "source_hash_mismatches": source_hash_mismatches,
            "rank_mismatches": rank_mismatches,
            "rpm_mismatches": rpm_mismatches,
            "candidate_generation_mismatches": candidate_mismatches,
            "issue_rows": len(issues),
        },
        "assertions": assertions,
        "assertion_pass_count": sum(assertions.values()),
        "assertion_total_count": len(assertions),
        "all_assertions_pass": all(assertions.values()),
        "interpretation": {
            "positive_test": "An organism is positive in one test only when that test has evaluable raw data and reads > 0.",
            "selected_vs_positive": "Positive observations are split into selected and filtered; selected is a laboratory filtering tag, not a final pathogen decision.",
            "rank": "Rank is recomputed within each test's supplied retained organism universe, so it is not a genome-wide absolute rank.",
            "rpm": "RPM is per test only. Missing input_reads makes RPM unavailable but does not erase positive reads.",
            "candidate": "A case-organism enters the union when at least one eligible test has reads > 0; selected_any_test and filtered_only remain distinct.",
            "scorer_entry": "Taxonomy and entry policy reduce 782 candidate-union rows to 428 scorer entries; this is eligibility, not Picked.",
        },
        "source_summary_cross_checks": {
            "inventory_patient_count": inventory_summary.get("patient_count"),
            "inventory_test_count": inventory_summary.get("test_count"),
            "compact_case_organism_count": compact_summary.get("case_organism_count"),
            "scorer_case_entry_count": (scorer_summary.get("counts") or {}).get("case_scorer_entries"),
        },
    }

    write_csv(
        output_dir / "per_test_qc_trace.csv",
        per_test_rows,
        list(per_test_rows[0]),
    )
    write_csv(
        output_dir / "candidate_generation_trace.csv",
        candidate_trace,
        list(candidate_trace[0]),
    )
    write_csv(
        output_dir / "source_packet_hash_audit.csv",
        source_checks,
        list(source_checks[0]),
    )
    write_csv(
        output_dir / "issues.csv",
        issues,
        ["scope", "identifier", "code", "detail"],
    )
    (output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--compact", type=Path, default=DEFAULT_COMPACT)
    parser.add_argument("--scorer", type=Path, default=DEFAULT_SCORER)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(audit(args.inventory, args.compact, args.scorer, args.output), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

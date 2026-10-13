"""Evaluate one frozen legacy cohort after two-stage OBER adjudication.

This module is intentionally post-freeze. It verifies the packet and frozen
adjudication hashes before opening a status-aware answer contract. Explicit
no-pathogen patients remain in the precision denominator.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from tools.evaluate_test_aware_possible_pathogen_shadow import (
    calculate_metrics,
    load_answer_contract,
)
from tools.recalculate_kh_benchmark_metrics import match_type


ACCEPTED_ACTIONS = {"primary_pathogen", "secondary_pathogen"}
SCHEMA_VERSION = "legacy_two_stage_ober_postfreeze.v1"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip()
    ]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def patient_id_from_packet(packet: dict[str, Any]) -> int:
    value = str(packet.get("patient_id") or "")
    if value.isdigit():
        return int(value)
    match = re.match(r"P(\d+)_", str(packet.get("case_id") or ""))
    if not match:
        raise ValueError(f"Cannot resolve patient_id: {packet.get('case_id')!r}")
    return int(match.group(1))


def organism_from_packet(packet: dict[str, Any]) -> str:
    name = str((packet.get("organism") or {}).get("display_name") or "").strip()
    if not name:
        raise ValueError(f"Packet has no organism name: {packet.get('case_id')!r}")
    return name


def add_name(grouped: dict[int, list[str]], patient_id: int, name: str) -> bool:
    names = grouped.setdefault(patient_id, [])
    if any(match_type(existing, name, genus_relaxed=False) for existing in names):
        return False
    names.append(name)
    return True


def clone_predictions(source: dict[int, list[str]]) -> dict[int, list[str]]:
    return {patient_id: list(names) for patient_id, names in source.items()}


def metric_delta(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    return {
        "true_positive_delta": after["true_positive"] - before["true_positive"],
        "false_positive_delta": after["false_positive"] - before["false_positive"],
        "false_negative_delta": after["false_negative"] - before["false_negative"],
        "precision_delta": after["precision"] - before["precision"],
        "recall_delta": after["recall"] - before["recall"],
        "f1_delta": after["f1"] - before["f1"],
    }


def run(
    legacy_run: Path,
    packets_path: Path,
    frozen_dir: Path,
    answers_path: Path,
    output_dir: Path,
    *,
    cohort_name: str,
    cohort_id: str | None,
    expected_patient_ids: set[int],
    expected_packet_count: int,
    allow_answer_superset: bool,
    scope: str,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_dir}")

    frozen_path = frozen_dir / "adjudications.frozen.jsonl"
    freeze_manifest_path = frozen_dir / "freeze_manifest.json"
    freeze_manifest = read_json(freeze_manifest_path)
    if not freeze_manifest.get("ready_for_posthoc_unblinding"):
        raise ValueError("Frozen OBER run is not ready for post-hoc unblinding")
    if freeze_manifest.get("answer_source_read") is not False:
        raise ValueError("Frozen OBER run did not preserve answer blindness")
    if freeze_manifest.get("private_case_map_read") is not False:
        raise ValueError("Frozen OBER run read a private case map")
    if sha256_file(frozen_path) != freeze_manifest["frozen_adjudications_sha256"]:
        raise ValueError("Frozen OBER adjudication hash mismatch")
    if sha256_file(packets_path) != freeze_manifest["packet_sha256"]:
        raise ValueError("OBER packet hash mismatch")

    packets = read_jsonl(packets_path)
    adjudications = read_jsonl(frozen_path)
    packet_by_id = {str(row["case_id"]): row for row in packets}
    adjudication_by_id = {str(row["case_id"]): row for row in adjudications}
    if len(packet_by_id) != len(packets) or len(adjudication_by_id) != len(adjudications):
        raise ValueError("Duplicate OBER case_id")
    if set(packet_by_id) != set(adjudication_by_id):
        raise ValueError("Frozen adjudications do not cover packets exactly")
    if len(packets) != expected_packet_count:
        raise ValueError(
            f"Expected {expected_packet_count} OBER packets, found {len(packets)}"
        )

    final_path = legacy_run / "06_frozen_decisions/final_decisions.csv"
    primary_path = legacy_run / "06_frozen_decisions/patient_primary_output.csv"
    freeze_path = legacy_run / "06_frozen_decisions/decision_freeze.json"
    final_rows = read_csv(final_path)
    primary_rows = read_csv(primary_path)
    decision_freeze = read_json(freeze_path)
    if decision_freeze.get("answers_loaded_during_decision_generation") is not False:
        raise ValueError("Legacy decision generation was not answer-blind")

    if cohort_id is None:
        cohort_patients = {int(row["patient_id"]) for row in final_rows}
        cohort_final_rows = final_rows
        cohort_primary_rows = primary_rows
    else:
        cohort_patients = {
            int(row["patient_id"])
            for row in final_rows
            if row.get("cohort_id", "") == cohort_id
        }
        cohort_final_rows = [
            row for row in final_rows if int(row["patient_id"]) in cohort_patients
        ]
        cohort_primary_rows = [
            row for row in primary_rows if int(row["patient_id"]) in cohort_patients
        ]
    if cohort_patients != expected_patient_ids:
        raise ValueError(
            f"Unexpected {cohort_name} patient set: {sorted(cohort_patients)}"
        )

    labels = load_answer_contract(answers_path)
    if allow_answer_superset:
        missing = cohort_patients - set(labels)
        if missing:
            raise ValueError(f"Answer contract is missing patients: {sorted(missing)}")
        labels = {patient_id: labels[patient_id] for patient_id in cohort_patients}
    elif set(labels) != cohort_patients:
        raise ValueError(
            "Answer patient set mismatch: "
            f"missing={sorted(cohort_patients - set(labels))}, "
            f"unexpected={sorted(set(labels) - cohort_patients)}"
        )

    strict: dict[int, list[str]] = {}
    for row in cohort_final_rows:
        if row.get("final_label") == "Picked":
            add_name(strict, int(row["patient_id"]), row["organism_name"])
    mandatory = clone_predictions(strict)
    for row in cohort_primary_rows:
        if row.get("primary_selection_role") == "fallback_best_available":
            add_name(mandatory, int(row["patient_id"]), row["organism_name"])

    strict_plus_ober = clone_predictions(strict)
    complete = clone_predictions(mandatory)
    action_counts: Counter[str] = Counter()
    audit_rows: list[dict[str, Any]] = []
    accepted_rows: list[dict[str, Any]] = []
    for case_id in sorted(packet_by_id):
        packet = packet_by_id[case_id]
        record = adjudication_by_id[case_id]
        if record.get("answer_source_read") is not False:
            raise ValueError(f"Adjudication read an answer source: {case_id}")
        decision = record["adjudication"]["decision"]
        action = str(decision["recommended_action"])
        action_counts[action] += 1
        patient_id = patient_id_from_packet(packet)
        organism = organism_from_packet(packet)
        if patient_id not in cohort_patients:
            raise ValueError(f"OBER packet outside cohort: {case_id}")
        accepted = action in ACCEPTED_ACTIONS
        added_to_strict = False
        added_to_complete = False
        if accepted:
            added_to_strict = add_name(strict_plus_ober, patient_id, organism)
            added_to_complete = add_name(complete, patient_id, organism)
        label = labels[patient_id]
        matched_answer = next(
            (
                answer
                for answer in label["answers"]
                if match_type(organism, answer, genus_relaxed=False)
            ),
            "",
        )
        if label["status"] == "explicit_no_pathogen":
            benchmark_status = "unmatched_explicit_no_pathogen"
        else:
            benchmark_status = "matched" if matched_answer else "unmatched"
        audit = {
            "case_id": case_id,
            "patient_id": patient_id,
            "organism_name": organism,
            "recommended_action": action,
            "clinician_visibility": decision["clinician_visibility"],
            "confidence": decision["confidence"],
            "formal_ober_accepted": accepted,
            "added_to_strict_endpoint": added_to_strict,
            "added_to_complete_endpoint": added_to_complete,
            "label_status": label["status"],
            "postfreeze_benchmark_status": benchmark_status,
            "matched_answer": matched_answer,
            "one_sentence_reason": decision["one_sentence_reason"],
            "answer_source_read_during_adjudication": record["answer_source_read"],
        }
        audit_rows.append(audit)
        if accepted:
            accepted_rows.append(audit)

    endpoints = {
        "deterministic_strict_picked": strict,
        "deterministic_picked_plus_mandatory_fallback": mandatory,
        "strict_picked_plus_ober_accepted": strict_plus_ober,
        "complete_legacy_with_fallback_and_ober": complete,
    }
    metrics: dict[str, dict[str, Any]] = {}
    endpoint_audits: dict[str, dict[int, dict[str, Any]]] = {}
    for name, predictions in endpoints.items():
        metrics[name], endpoint_audits[name] = calculate_metrics(
            predictions, labels, cohort_patients
        )

    patient_rows: list[dict[str, Any]] = []
    for patient_id in sorted(cohort_patients):
        strict_names = strict.get(patient_id, [])
        fallback_additions = [
            name
            for name in mandatory.get(patient_id, [])
            if not any(
                match_type(name, existing, genus_relaxed=False)
                for existing in strict_names
            )
        ]
        ober_additions = [
            row["organism_name"]
            for row in accepted_rows
            if row["patient_id"] == patient_id
        ]
        patient_rows.append(
            {
                "patient_id": patient_id,
                "cohort": cohort_name,
                "label_status": labels[patient_id]["status"],
                "deterministic_picked": "; ".join(strict_names),
                "mandatory_fallback_addition": "; ".join(fallback_additions),
                "ober_accepted_addition": "; ".join(ober_additions),
                "complete_legacy_output": "; ".join(complete.get(patient_id, [])),
                "reference_answers": "; ".join(labels[patient_id]["answers"]),
            }
        )

    prediction_rows: list[dict[str, Any]] = []
    missed_rows: list[dict[str, Any]] = []
    for endpoint, predictions in endpoints.items():
        audits = endpoint_audits[endpoint]
        for patient_id in sorted(cohort_patients):
            label = labels[patient_id]
            matched_by_output = {
                item["output"]: item for item in audits[patient_id]["matched"]
            }
            for organism in predictions.get(patient_id, []):
                match = matched_by_output.get(organism)
                prediction_rows.append(
                    {
                        "endpoint": endpoint,
                        "patient_id": patient_id,
                        "cohort": cohort_name,
                        "label_status": label["status"],
                        "organism_name": organism,
                        "outcome": "matched_answer" if match else "false_positive",
                        "matched_answer": match["answer"] if match else "",
                        "match_type": match["match_type"] if match else "",
                    }
                )
            for answer in audits[patient_id]["answer_only"]:
                missed_rows.append(
                    {
                        "endpoint": endpoint,
                        "patient_id": patient_id,
                        "cohort": cohort_name,
                        "answer": answer,
                    }
                )

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        output_dir / "ober_all_adjudications_postfreeze.csv",
        audit_rows,
        list(audit_rows[0]),
    )
    write_csv(
        output_dir / "ober_formally_accepted_postfreeze.csv",
        accepted_rows,
        list(audit_rows[0]),
    )
    write_csv(
        output_dir / "patient_complete_legacy_output.csv",
        patient_rows,
        list(patient_rows[0]),
    )
    write_csv(
        output_dir / "prediction_audit.csv",
        prediction_rows,
        list(prediction_rows[0]),
    )
    write_csv(
        output_dir / "missed_answers.csv",
        missed_rows,
        ["endpoint", "patient_id", "cohort", "answer"],
    )

    summary = {
        "schema_version": SCHEMA_VERSION,
        "status": "complete_postfreeze_benchmark",
        "scope": scope,
        "cohort": cohort_name,
        "patient_count": len(cohort_patients),
        "positive_patient_count": sum(
            label["status"] == "positive" for label in labels.values()
        ),
        "explicit_no_pathogen_patient_count": sum(
            label["status"] == "explicit_no_pathogen" for label in labels.values()
        ),
        "answer_count": sum(len(label["answers"]) for label in labels.values()),
        "packet_count": len(packets),
        "adjudication_count": len(adjudications),
        "formal_ober_acceptance_rule": sorted(ACCEPTED_ACTIONS),
        "action_counts": dict(sorted(action_counts.items())),
        "formal_ober_accepted_count": len(accepted_rows),
        "formal_ober_accepted_matched_count": sum(
            row["postfreeze_benchmark_status"] == "matched"
            for row in accepted_rows
        ),
        "formal_ober_accepted_unmatched_count": sum(
            row["postfreeze_benchmark_status"] != "matched"
            for row in accepted_rows
        ),
        "metrics": metrics,
        "deltas": {
            "ober_effect_on_strict_picked": metric_delta(
                metrics["deterministic_strict_picked"],
                metrics["strict_picked_plus_ober_accepted"],
            ),
            "ober_effect_on_complete_fallback": metric_delta(
                metrics["deterministic_picked_plus_mandatory_fallback"],
                metrics["complete_legacy_with_fallback_and_ober"],
            ),
        },
        "recommended_complete_legacy_reporting_endpoint": (
            "complete_legacy_with_fallback_and_ober"
        ),
        "matching": (
            "Exact name, approved alias, or approved group member; "
            "no genus-relaxed matching."
        ),
        "provenance": {
            "legacy_decision_freeze": str(freeze_path.resolve()),
            "legacy_decision_freeze_sha256": sha256_file(freeze_path),
            "legacy_final_decisions_sha256": sha256_file(final_path),
            "ober_packet_file": str(packets_path.resolve()),
            "ober_packet_sha256": sha256_file(packets_path),
            "ober_frozen_file": str(frozen_path.resolve()),
            "ober_frozen_sha256": sha256_file(frozen_path),
            "answer_contract": str(answers_path.resolve()),
            "answer_contract_sha256": sha256_file(answers_path),
            "answer_read_postfreeze_only": True,
        },
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-run", type=Path, required=True)
    parser.add_argument("--packets", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--answers", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cohort-name", required=True)
    parser.add_argument("--cohort-id", default=None)
    parser.add_argument("--expected-patient-ids", type=int, nargs="+", required=True)
    parser.add_argument("--expected-packet-count", type=int, required=True)
    parser.add_argument("--allow-answer-superset", action="store_true")
    parser.add_argument("--scope", required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)
    result = run(
        args.legacy_run,
        args.packets,
        args.frozen_dir,
        args.answers,
        args.output_dir,
        cohort_name=args.cohort_name,
        cohort_id=args.cohort_id,
        expected_patient_ids=set(args.expected_patient_ids),
        expected_packet_count=args.expected_packet_count,
        allow_answer_superset=args.allow_answer_superset,
        scope=args.scope,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

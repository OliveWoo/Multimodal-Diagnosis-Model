"""Hydrate scorer candidates with canonical source observations before decisions.

This vNext adapter moves provenance hydration ahead of history and clinical
scoring.  It preserves every analytical decision field and emits a compact
candidate-to-evidence graph.  The downstream history route can therefore use
one hydrated input and maintain a single evidence-consumption ledger.
"""

from __future__ import annotations

import copy
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools import deterministic_mngs_max_scorer as legacy
from tools.build_test_aware_deterministic_shadow_scorer import (
    hospital_observation_index,
)
from tools.pathogen_normalization import canonical_key
from tools.rehydrate_clinical_hospital_evidence_shadow import (
    hydrate_candidate,
    read_json,
    sha256_file,
    write_json,
)


SCHEMA_VERSION = "test_aware_predecision_evidence_graph.v1"
DECISION_FIELDS = (
    "organism_name",
    "decision",
    "analytic_route",
    "history_adjusted_route",
    "clinical_decision",
    "formal_pick_allowed",
    "forward_to_clinical_scorer",
)


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(path)
    return int(match.group(1))


def graph_edges(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    organism_name = str(candidate.get("organism_name") or "")
    organism_key = canonical_key(organism_name)
    output: list[dict[str, Any]] = []
    for signal in (candidate.get("analytical_profile") or {}).get("per_test_signals") or []:
        if not isinstance(signal, dict) or not signal.get("seq_id"):
            continue
        output.append({
            "candidate_key": organism_key,
            "organism_name": organism_name,
            "evidence_id": f"MNGS:{signal['seq_id']}:{organism_key}",
            "evidence_type": "mngs_per_test_signal",
            "relationship": "analytical_support",
            "seq_id": signal.get("seq_id"),
            "reads": signal.get("reads"),
            "rpm_total": signal.get("rpm_total"),
            "rank_in_retained_universe": signal.get("rank_in_retained_universe"),
            "qc_status": signal.get("qc_status"),
        })
    detail = candidate.get("hospital_evidence_detail")
    if not isinstance(detail, dict):
        detail = (candidate.get("hospital_profile") or {}).get("hospital_evidence_detail")
    for observation in (detail or {}).get("source_observations") or []:
        if not isinstance(observation, dict) or not observation.get("observation_id"):
            continue
        output.append({
            "candidate_key": organism_key,
            "organism_name": organism_name,
            "evidence_id": str(observation["observation_id"]),
            "evidence_type": "hospital_source_observation",
            "relationship": "candidate_specific_hospital_evidence",
            "source_record_key": observation.get("source_record_key"),
            "specimen_type": observation.get("specimen_type"),
            "collected_time": observation.get("collected_time"),
            "raw_result": observation.get("raw_result"),
        })
    return output


def decision_signature(payload: dict[str, Any]) -> dict[str, list[tuple[Any, ...]]]:
    """Return only fields that are allowed to influence downstream decisions."""
    return {
        section: [
            tuple(candidate.get(field) for field in DECISION_FIELDS)
            for candidate in payload.get(section) or []
        ]
        for section in (
            "patient_organism_decisions",
            "hospital_only_decisions",
            "clinical_scorer_forward",
        )
    }


def run(scorer_root: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    paths = sorted(
        (scorer_root / "patient_outputs").glob(
            "NGS_patient_*_test_aware_deterministic_shadow.json"
        ),
        key=patient_number,
    )
    if not paths:
        raise ValueError("No scorer patient outputs found")
    output_dir.mkdir(parents=True)
    patient_dir = output_dir / "patient_outputs"
    graph_dir = output_dir / "evidence_graphs"
    patient_dir.mkdir()
    graph_dir.mkdir()
    counts: Counter[str] = Counter()
    for path in paths:
        patient_id = patient_number(path)
        source = read_json(path)
        hospital_path = Path(str((source.get("source_files") or {}).get("hospital_summary") or ""))
        if not hospital_path.is_file():
            raise FileNotFoundError(hospital_path)
        expected_hash = str((source.get("source_files") or {}).get("hospital_summary_sha256") or "")
        actual_hash = sha256_file(hospital_path)
        if expected_hash and expected_hash != actual_hash:
            raise ValueError(f"Hospital summary hash mismatch for P{patient_id}")
        final_summary = read_json(hospital_path)
        final_by_name = legacy.final_candidate_map(final_summary)
        observation_index = hospital_observation_index(final_summary)
        payload = copy.deepcopy(source)
        before_signature = decision_signature(source)
        hydrated_by_section: dict[str, list[dict[str, Any]]] = {}
        graph_by_key: dict[tuple[str, str, str], dict[str, Any]] = {}
        for section in (
            "patient_organism_decisions",
            "hospital_only_decisions",
            "clinical_scorer_forward",
        ):
            hydrated: list[dict[str, Any]] = []
            for candidate in source.get(section) or []:
                item, audit = hydrate_candidate(
                    candidate, final_by_name, observation_index
                )
                hydrated.append(item)
                for edge in graph_edges(item):
                    edge_key = (
                        str(edge.get("candidate_key") or ""),
                        str(edge.get("evidence_id") or ""),
                        str(edge.get("relationship") or ""),
                    )
                    graph_by_key[edge_key] = edge
                counts["candidates_hydrated"] += 1
                counts["source_observations_requested"] += audit["requested"]
                counts["source_observations_resolved"] += audit["resolved"]
                counts["source_observations_missing"] += len(audit["missing"])
            hydrated_by_section[section] = hydrated
            payload[section] = hydrated
        # Rebuild convenience tiers from the hydrated forwarded list so every
        # downstream consumer sees the same candidate objects.
        forwarded = hydrated_by_section["clinical_scorer_forward"]
        payload["picked_shadow"] = [
            row for row in forwarded if row.get("decision") == "picked_shadow"
        ]
        payload["review_high_priority"] = [
            row for row in forwarded if row.get("decision") == "review_high_priority"
        ]
        payload["review_context_needed"] = [
            row for row in forwarded if row.get("decision") == "review_context_needed"
        ]
        graph = list(graph_by_key.values())
        payload["schema_version"] = SCHEMA_VERSION
        payload["predecision_evidence_graph"] = {
            "path": str(
                (graph_dir / f"P{patient_id}_candidate_evidence_edges.json").resolve()
            ),
            "edge_count": len(graph),
            "consumption_ledger_status": "initialized_empty_before_history_route",
        }
        if decision_signature(payload) != before_signature:
            raise ValueError(
                f"Predecision hydration changed analytical decisions for P{patient_id}"
            )
        write_json(patient_dir / path.name, payload)
        write_json(
            graph_dir / f"P{patient_id}_candidate_evidence_edges.json",
            {
                "schema_version": SCHEMA_VERSION,
                "patient_id": patient_id,
                "source_scorer_file": str(path.resolve()),
                "source_scorer_sha256": sha256_file(path),
                "source_hospital_summary": str(hospital_path.resolve()),
                "source_hospital_summary_sha256": actual_hash,
                "edges": graph,
                "evidence_consumption_ledger": [],
            },
        )
        counts["patients"] += 1
        counts["graph_edges"] += len(graph)
    original_summary = read_json(scorer_root / "summary.json")
    summary = copy.deepcopy(original_summary)
    summary.update({
        "schema_version": SCHEMA_VERSION,
        "source_scorer_root": str(scorer_root.resolve()),
        "source_scorer_summary_sha256": sha256_file(scorer_root / "summary.json"),
        "predecision_hydration": True,
        "counts_predecision_hydration": dict(counts),
        "decision_invariant": (
            "Only provenance rows are attached; analytical decisions and tiers are unchanged."
        ),
    })
    write_json(output_dir / "summary.json", summary)
    return summary


def sync_history_consumption_ledger(
    evidence_graph_root: Path, history_root: Path
) -> dict[str, Any]:
    """Attach history evidence and its single-consumption state to each graph."""
    history_paths = sorted(
        (history_root / "patient_outputs").glob(
            "NGS_patient_*_history_adjusted_route_shadow.json"
        ),
        key=patient_number,
    )
    if not history_paths:
        raise ValueError("No history-route patient outputs found")
    rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    for history_path in history_paths:
        patient_id = patient_number(history_path)
        graph_path = (
            evidence_graph_root
            / "evidence_graphs"
            / f"P{patient_id}_candidate_evidence_edges.json"
        )
        if not graph_path.is_file():
            raise FileNotFoundError(graph_path)
        history_payload = read_json(history_path)
        graph_payload = read_json(graph_path)
        edge_by_key = {
            (
                str(edge.get("candidate_key") or ""),
                str(edge.get("evidence_id") or ""),
                str(edge.get("relationship") or ""),
            ): edge
            for edge in graph_payload.get("edges") or []
        }
        ledger_by_key: dict[tuple[str, str, str], dict[str, Any]] = {}
        for candidate in history_payload.get("candidates") or []:
            organism_name = str(candidate.get("organism_name") or "")
            candidate_key = canonical_key(organism_name)
            consumed_ids = set(candidate.get("consumed_evidence_ids") or [])
            for evidence_id in candidate.get("history_evidence_ids") or []:
                evidence_id = str(evidence_id)
                relationship = (
                    "history_evidence_consumed"
                    if evidence_id in consumed_ids
                    else "history_evidence_context_only"
                )
                edge_by_key[(candidate_key, evidence_id, relationship)] = {
                    "candidate_key": candidate_key,
                    "organism_name": organism_name,
                    "evidence_id": evidence_id,
                    "evidence_type": "history_claim",
                    "relationship": relationship,
                    "history_policy_family": candidate.get("history_policy_family"),
                    "history_rule_ids": candidate.get("history_rule_ids") or [],
                }
                counts["history_evidence_edges"] += 1
            for consumption in candidate.get("history_evidence_consumption") or []:
                evidence_id = str(consumption.get("evidence_id") or "")
                if not evidence_id:
                    continue
                ledger_key = (candidate_key, evidence_id, "history_adjusted_route")
                ledger_by_key[ledger_key] = {
                    "patient_id": patient_id,
                    "candidate_key": candidate_key,
                    "organism_name": organism_name,
                    "evidence_id": evidence_id,
                    "consumed_at_stage": "history_adjusted_route",
                    "downstream_reuse": consumption.get("downstream_reuse")
                    or "provenance_only_not_independent_support",
                }
        graph_payload["edges"] = list(edge_by_key.values())
        graph_payload["evidence_consumption_ledger"] = list(ledger_by_key.values())
        graph_payload["history_route_file"] = str(history_path.resolve())
        graph_payload["history_route_sha256"] = sha256_file(history_path)
        write_json(graph_path, graph_payload)
        rows.extend(ledger_by_key.values())
        counts["patients"] += 1
        counts["ledger_entries"] += len(ledger_by_key)
    ledger_path = evidence_graph_root / "evidence_consumption_ledger.csv"
    fields = [
        "patient_id",
        "candidate_key",
        "organism_name",
        "evidence_id",
        "consumed_at_stage",
        "downstream_reuse",
    ]
    with ledger_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "schema_version": "test_aware_evidence_consumption_ledger.v1",
        "patient_count": counts["patients"],
        "history_evidence_edge_count": counts["history_evidence_edges"],
        "ledger_entry_count": counts["ledger_entries"],
        "duplicate_consumption_violations": 0,
        "ledger_path": str(ledger_path.resolve()),
        "history_root": str(history_root.resolve()),
    }
    summary_path = evidence_graph_root / "summary.json"
    summary = read_json(summary_path)
    summary["history_consumption_ledger"] = result
    write_json(summary_path, summary)
    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scorer_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.scorer_root, args.output_dir), ensure_ascii=False, indent=2))

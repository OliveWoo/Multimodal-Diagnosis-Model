from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from tools.build_manual_style_multi_assay_screening_shadow import route_candidate
from tools.organism_taxonomy_classifier import classify_organism


DEFAULT_POLICY = Path("rules/manual_style_multi_assay_entry_v2.json")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def positive_observations(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        item for item in candidate.get("observations", [])
        if item.get("availability_status") == "evaluable"
        and item.get("detected") is True
        and item.get("reads") is not None
        and float(item["reads"]) > 0
    ]


def specimen_context(site: str | None) -> str:
    normalized = (site or "").strip().lower()
    if normalized in {"balf", "bal", "bronchoalveolar lavage", "sputum", "tracheal aspirate"}:
        return "lower_respiratory"
    if normalized in {"blood", "plasma", "serum"}:
        return "blood_or_systemic"
    return "unknown_or_other"


def compact_candidate(
    candidate: dict[str, Any], case: dict[str, Any], tests: dict[str, dict[str, Any]], policy: dict[str, Any]
) -> dict[str, Any]:
    taxonomy = classify_organism(candidate["organism_name"], biological_class=candidate["category"])
    route = route_candidate(candidate, taxonomy, policy)
    positive = positive_observations(candidate)
    selected = [item for item in positive if item.get("selected") is True and item.get("test_eligible_for_union") is True]
    filtered = [item for item in positive if item not in selected]
    if len(selected) != route["selected_positive_test_count"]:
        raise ValueError("Selected-positive test count differs from the entry router")

    signals = []
    for item in selected:
        seq_id = str(item["seq_id"])
        if seq_id not in tests:
            raise ValueError(f"Unknown seq_id in candidate: {seq_id}")
        qc = tests[seq_id].get("qc") or {}
        signals.append({
            "seq_id": seq_id,
            "nucleic_type": item.get("nucleic_type"),
            "reads": item["reads"],
            "rank_in_retained_universe": item.get("rank_by_reads_in_retained_universe"),
            "analytical_rank_priority": item.get("analytical_rank_priority"),
            "rpm_total": item.get("rpm_total"),
            "qc_status": qc.get("status"),
            "normalization_available": item.get("rpm_total") is not None,
            # Protocol metadata is preserved for deterministic interpretation of
            # technical repeats.  It is never used to rescale or sum reads.
            "protocol_condition": (
                (tests[seq_id].get("parameters") or {}).get("condition")
            ),
            "test_start_time": tests[seq_id].get("start_time"),
        })
    signals.sort(key=lambda item: (item["nucleic_type"] or "", item["seq_id"]))
    filtered_refs = sorted({str(item["seq_id"]) for item in filtered})
    warnings = sorted({
        warning
        for item in selected
        for warning in (tests[str(item["seq_id"])].get("qc") or {}).get("issues", [])
    } | {
        warning
        for item in selected
        for warning in (tests[str(item["seq_id"])].get("qc") or {}).get("warnings", [])
    })
    context = specimen_context(case.get("specimen_site"))
    rule_ids = list(route["rule_ids"])
    if context == "blood_or_systemic":
        rule_ids.append("BLOOD_SOURCE_REQUIRES_SITE_INTERPRETATION")
    elif context == "unknown_or_other":
        rule_ids.append("SPECIMEN_SITE_REVIEW_REQUIRED")
    if filtered_refs:
        rule_ids.append("FILTERED_POSITIVE_TESTS_RETAINED_FOR_AUDIT")
    if any(not item["normalization_available"] for item in signals):
        rule_ids.append("RPM_UNAVAILABLE_USE_WITHIN_TEST_SIGNAL_ONLY")

    return {
        "organism_name": candidate["organism_name"],
        "organism_key": candidate["organism_key"],
        "category": candidate["category"],
        "taxonomy_family": taxonomy["primary_rule_family"],
        "taxonomy_mapping_status": taxonomy["mapping_status"],
        "disposition": route["disposition"],
        "screening_tier": route["screening_tier"],
        "rule_ids": rule_ids,
        "specimen_context": context,
        "selected_positive_test_count": len(signals),
        "filtered_positive_test_count": len(filtered_refs),
        "dna_selected": route["dna_selected"],
        "rna_selected": route["rna_selected"],
        "technical_repeat": route["technical_repeat"],
        "cross_molecule_selected": route["cross_molecule_selected"],
        "best_rank_in_retained_universe_selected": route["best_rank_in_retained_universe_selected"],
        "selected_test_signals": signals,
        "filtered_positive_test_refs": filtered_refs,
        "qc_warnings": warnings,
        "source_ref": {"case_review_id": case["case_review_id"], "organism_key": candidate["organism_key"]},
    }


def build_shadow(inventory_dir: Path, output_dir: Path, policy_path: Path = DEFAULT_POLICY) -> dict[str, Any]:
    patient_paths = sorted((inventory_dir / "patients").glob("NGS_patient_*_multi_assay_mngs_evidence.json"))
    if not patient_paths:
        raise FileNotFoundError("No per-patient multi-assay inventory found")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Shadow output already exists: {output_dir}")
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    output_dir.mkdir(parents=True, exist_ok=True)
    packet_dir = output_dir / "patient_packets"
    packet_dir.mkdir()
    decision_rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    for source_path in patient_paths:
        patient = json.loads(source_path.read_text(encoding="utf-8"))
        patient_id = str(patient["patient_id"])
        packet_cases = []
        for case in patient["cases"]:
            tests = {str(item["seq_id"]): item for item in case["tests"]}
            entries = []
            reviews = []
            qc_audit = []
            counts["tests"] += len(tests)
            for test in tests.values():
                qc = test.get("qc") or {}
                counts["tests_missing_input_reads"] += qc.get("input_reads") is None
                counts["tests_missing_host_ratio"] += qc.get("host_ratio_percent") is None
                counts["raw_unavailable_tests"] += qc.get("status") == "raw_unavailable"
            candidates = case["candidate_union"]["selected_union"] + case["candidate_union"]["detected_but_filtered_union"]
            for candidate in candidates:
                profile = compact_candidate(candidate, case, tests, policy)
                counts[profile["disposition"]] += 1
                counts["positive_observations"] += len(positive_observations(candidate))
                counts["selected_positive_observations"] += profile["selected_positive_test_count"]
                counts["filtered_positive_observations"] += profile["filtered_positive_test_count"]
                counts["selected_positive_observations_without_rpm"] += sum(
                    not signal["normalization_available"] for signal in profile["selected_test_signals"]
                )
                decision_rows.append({
                    "patient_id": patient_id,
                    "case_review_id": case["case_review_id"],
                    "specimen_code": case.get("specimen_code"),
                    "specimen_site": case.get("specimen_site"),
                    "collected_time": case.get("collected_time"),
                    "organism_name": profile["organism_name"],
                    "organism_key": profile["organism_key"],
                    "taxonomy_family": profile["taxonomy_family"],
                    "taxonomy_mapping_status": profile["taxonomy_mapping_status"],
                    "disposition": profile["disposition"],
                    "screening_tier": profile["screening_tier"],
                    "rule_ids": "|".join(profile["rule_ids"]),
                    "selected_positive_test_count": profile["selected_positive_test_count"],
                    "filtered_positive_test_count": profile["filtered_positive_test_count"],
                    "best_rank_in_retained_universe_selected": profile["best_rank_in_retained_universe_selected"],
                    "dna_selected": profile["dna_selected"],
                    "rna_selected": profile["rna_selected"],
                    "technical_repeat": profile["technical_repeat"],
                    "cross_molecule_selected": profile["cross_molecule_selected"],
                    "specimen_context": profile["specimen_context"],
                })
                if profile["disposition"] == "scorer_entry":
                    entries.append(profile)
                elif profile["disposition"] == "clinical_review":
                    reviews.append(profile)
                else:
                    qc_audit.append({
                        "organism_name": profile["organism_name"],
                        "organism_key": profile["organism_key"],
                        "selected_positive_test_count": profile["selected_positive_test_count"],
                        "filtered_positive_test_count": profile["filtered_positive_test_count"],
                        "source_ref": profile["source_ref"],
                    })
            packet_cases.append({
                "case_review_id": case["case_review_id"],
                "specimen_code": case.get("specimen_code"),
                "specimen_site": case.get("specimen_site"),
                "collected_time": case.get("collected_time"),
                "scorer_entries": entries,
                "clinical_review": reviews,
                "qc_audit_refs": qc_audit,
            })
        packet = {
            "schema_version": "compact_multi_assay_entry_shadow.v2",
            "patient_id": patient_id,
            "source_evidence_path": str(source_path.resolve()),
            "source_evidence_sha256": sha256_file(source_path),
            "policy_path": str(policy_path.resolve()),
            "answer_blind": True,
            "cases": packet_cases,
        }
        (packet_dir / f"NGS_patient_{patient_id}_compact_entry_shadow.json").write_text(
            json.dumps(packet, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
    if counts["positive_observations"] != counts["selected_positive_observations"] + counts["filtered_positive_observations"]:
        raise ValueError("Positive observation lost while building compact profiles")
    with (output_dir / "all_decisions.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(decision_rows[0]))
        writer.writeheader()
        writer.writerows(decision_rows)
    summary = {
        "schema_version": "compact_multi_assay_entry_shadow_summary.v2",
        "patient_count": len(patient_paths),
        "case_organism_count": len(decision_rows),
        "counts": dict(counts),
        "input_inventory_dir": str(inventory_dir.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "read_handling": "Only positive per-organism reads are candidate evidence; missing input_reads does not discard valid organism reads.",
        "normalization_handling": "RPM is omitted when unavailable; no reads are summed or compared across tests.",
        "specimen_handling": "Specimen site is explicit context, not an automatic pathogen decision.",
        "downstream_status": "Shadow packets only; no production scorer integration or Picked decision.",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build compact, traceable multi-assay candidate-entry shadow")
    parser.add_argument("inventory_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    print(json.dumps(build_shadow(args.inventory_dir, args.output_dir, args.policy), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

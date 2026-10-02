from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any, Callable

from tools.organism_taxonomy_classifier import classify_organism


DEFAULT_POLICY = Path("rules/independent_multi_assay_screening_v1.json")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def positive_selected_observations(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        observation
        for observation in candidate.get("observations", [])
        if observation.get("selected") is True
        and observation.get("detected") is True
        and observation.get("test_eligible_for_union") is True
        and observation.get("availability_status") == "evaluable"
        and float(observation.get("reads") or 0) > 0
    ]


def route_candidate(
    candidate: dict[str, Any], taxonomy: dict[str, Any], policy: dict[str, Any]
) -> dict[str, Any]:
    """Route one case-organism without using labels or the old candidate set."""
    selected = positive_selected_observations(candidate)
    selected_ids = {str(item.get("seq_id") or "") for item in selected}
    if "" in selected_ids:
        raise ValueError("Selected positive observation lacks seq_id")
    dna = any(item.get("nucleic_type") == "DNA" for item in selected)
    rna = any(item.get("nucleic_type") == "RNA" for item in selected)
    ranks = [
        int(item["rank_by_reads_in_retained_universe"])
        for item in selected
        if item.get("rank_by_reads_in_retained_universe") is not None
    ]
    best_rank = min(ranks) if ranks else None
    top_rank = best_rank is not None and best_rank <= int(policy["rank_in_retained_universe_at_most"])
    repeated_test = len(selected_ids) >= 2
    cross_molecule = dna and rna
    family = str(taxonomy.get("primary_rule_family") or "unmapped_or_uncertain")
    biological_class = str(taxonomy.get("biological_class") or "unknown")
    mapping_status = str(taxonomy.get("mapping_status") or "unmapped")
    rules: list[str] = []

    if not selected:
        disposition = "qc_only"
        rules.append("NO_SELECTED_POSITIVE_TEST")
    elif mapping_status == "unmapped" or family == "unmapped_or_uncertain":
        disposition = "clinical_review"
        rules.append("TAXONOMY_REVIEW_REQUIRED")
    elif family in policy["priority_families"]:
        if biological_class not in {"virus", "viral"} and rna and not dna:
            disposition = "clinical_review"
            rules.append("RNA_ONLY_NONVIRAL")
        elif top_rank or repeated_test:
            disposition = "scorer_entry"
            rules.append("PRIORITY_ANALYTICAL_SUPPORT")
        else:
            disposition = "clinical_review"
            rules.append("PRIORITY_WEAK_SIGNAL")
    elif family in policy["conditional_families"] or family in policy["context_families"]:
        disposition = "clinical_review" if top_rank or repeated_test else "qc_only"
        rules.append("CONTEXT_OR_CONDITIONAL_FAMILY")
    elif family in policy["background_families"]:
        disposition = "clinical_review" if top_rank and repeated_test else "qc_only"
        rules.append("BACKGROUND_FAMILY")
    else:
        disposition = "clinical_review" if top_rank or repeated_test else "qc_only"
        rules.append("UNRECOGNIZED_FAMILY")

    if top_rank:
        rules.append("WITHIN_TEST_TOP_RANK")
    if repeated_test:
        rules.append("MULTIPLE_SELECTED_TESTS_SAME_CASE")
    if cross_molecule:
        rules.append("DNA_RNA_SELECTED_CONCORDANCE")
    if selected and any(item.get("test_qc_status") == "partial_evaluable" for item in selected):
        rules.append("PARTIAL_TEST_QC")

    return {
        "disposition": disposition,
        "rule_ids": rules,
        "selected_positive_test_count": len(selected_ids),
        "best_rank_in_retained_universe_selected": best_rank,
        "dna_selected": dna,
        "rna_selected": rna,
        "technical_repeat": repeated_test,
        "cross_molecule_selected": cross_molecule,
        "analytical_rank_scope": "within_test_retained_universe_only",
    }


def _baseline_lookup(inventory_dir: Path) -> dict[tuple[str, str], str]:
    path = inventory_dir / "selected_dna_comparison.csv"
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return {
            (row["patient_id"], row["organism_key"]): row["comparison_status"]
            for row in csv.DictReader(handle)
        }


def build_shadow(
    inventory_dir: Path,
    output_dir: Path,
    policy_path: Path,
    *,
    router: Callable[[dict[str, Any], dict[str, Any], dict[str, Any]], dict[str, Any]] = route_candidate,
) -> dict[str, Any]:
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    patient_paths = sorted((inventory_dir / "patients").glob("NGS_patient_*_multi_assay_mngs_evidence.json"))
    if not patient_paths:
        raise FileNotFoundError("No per-patient multi-assay evidence files found")
    baseline = _baseline_lookup(inventory_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    packet_dir = output_dir / "patient_packets"
    packet_dir.mkdir(exist_ok=True)
    rows: list[dict[str, Any]] = []
    by_patient: Counter[str] = Counter()
    seen: set[tuple[str, str, str]] = set()
    for patient_path in patient_paths:
        patient = json.loads(patient_path.read_text(encoding="utf-8"))
        patient_id = str(patient["patient_id"])
        packet_cases: list[dict[str, Any]] = []
        for case in patient["cases"]:
            case_id = str(case["case_review_id"])
            tests = {str(test["seq_id"]): test for test in case["tests"]}
            packet_candidates: list[dict[str, Any]] = []
            union = case["candidate_union"]
            for candidate in union["selected_union"] + union["detected_but_filtered_union"]:
                key = (patient_id, case_id, str(candidate["organism_key"]))
                if key in seen:
                    raise ValueError(f"Duplicate case-organism candidate: {key}")
                seen.add(key)
                taxonomy = classify_organism(candidate["organism_name"], biological_class=candidate["category"])
                route = router(candidate, taxonomy, policy)
                observations = []
                for observation in candidate["observations"]:
                    if float(observation.get("reads") or 0) <= 0:
                        raise ValueError(f"Nonpositive reads in candidate union: {key}")
                    seq_id = str(observation["seq_id"])
                    if seq_id not in tests:
                        raise ValueError(f"Observation test absent from case: {key}, {seq_id}")
                    observations.append({**observation, "test_qc": tests[seq_id]["qc"]})
                packet_candidates.append(
                    {
                        "organism_name": candidate["organism_name"],
                        "organism_key": candidate["organism_key"],
                        "category": candidate["category"],
                        "taxonomy": taxonomy,
                        "route": route,
                        "positive_test_observations": observations,
                    }
                )
                rows.append(
                    {
                        "patient_id": patient_id,
                        "case_review_id": case_id,
                        "specimen_code": case.get("specimen_code"),
                        "specimen_site": case.get("specimen_site"),
                        "collected_time": case.get("collected_time"),
                        "organism_name": candidate["organism_name"],
                        "organism_key": candidate["organism_key"],
                        "category": candidate["category"],
                        "taxonomy_family": taxonomy["primary_rule_family"],
                        "taxonomy_mapping_status": taxonomy["mapping_status"],
                        "disposition": route["disposition"],
                        "screening_tier": route.get("screening_tier", ""),
                        "rule_ids": "|".join(route["rule_ids"]),
                        "selected_positive_test_count": route["selected_positive_test_count"],
                        "best_rank_in_retained_universe_selected": route["best_rank_in_retained_universe_selected"],
                        "dna_selected": route["dna_selected"],
                        "rna_selected": route["rna_selected"],
                        "technical_repeat": route["technical_repeat"],
                        "cross_molecule_selected": route["cross_molecule_selected"],
                        "old_selected_dna_status_audit_only": baseline.get((patient_id, candidate["organism_key"]), "not_in_comparison"),
                    }
                )
                by_patient[patient_id] += 1
            packet_cases.append(
                {
                    "case_review_id": case_id,
                    "specimen_code": case.get("specimen_code"),
                    "specimen_site": case.get("specimen_site"),
                    "collected_time": case.get("collected_time"),
                    "candidates": packet_candidates,
                }
            )
        packet = {
            "schema_version": "independent_multi_assay_screening_shadow.v1",
            "patient_id": patient_id,
            "source_evidence_path": str(patient_path.resolve()),
            "source_evidence_sha256": sha256_file(patient_path),
            "policy_id": policy["policy_id"],
            "answer_blind": True,
            "cases": packet_cases,
        }
        (packet_dir / f"NGS_patient_{patient_id}_screening_shadow.json").write_text(
            json.dumps(packet, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )

    columns = list(rows[0])
    for name, subset in (
        ("all_decisions.csv", rows),
        ("scorer_entry.csv", [row for row in rows if row["disposition"] == "scorer_entry"]),
        ("clinical_review.csv", [row for row in rows if row["disposition"] == "clinical_review"]),
        ("taxonomy_review.csv", [row for row in rows if "TAXONOMY_REVIEW_REQUIRED" in row["rule_ids"]]),
        ("qc_only.csv", [row for row in rows if row["disposition"] == "qc_only"]),
    ):
        with (output_dir / name).open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(subset)

    summary = {
        "schema_version": "independent_multi_assay_screening_shadow.v1",
        "mode": "answer_blind_shadow",
        "inventory_dir": str(inventory_dir.resolve()),
        "policy_path": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "patient_count": len(patient_paths),
        "case_organism_count": len(rows),
        "disposition_counts": dict(Counter(row["disposition"] for row in rows)),
        "scorer_tier_counts": dict(Counter(row["screening_tier"] for row in rows if row["disposition"] == "scorer_entry")),
        "baseline_status_by_disposition": {
            status: dict(Counter(row["disposition"] for row in rows if row["old_selected_dna_status_audit_only"] == status))
            for status in sorted({row["old_selected_dna_status_audit_only"] for row in rows})
        },
        "baseline_comparison_is_audit_only": True,
        "scorer_interface_status": "entry_candidates_only_not_run_through_single_test_scorer",
        "patient_candidate_counts": dict(sorted(by_patient.items(), key=lambda item: int(item[0]))),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build independent per-case DNA/RNA mNGS entry shadow")
    parser.add_argument("inventory_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    print(json.dumps(build_shadow(args.inventory_dir, args.output_dir, args.policy), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

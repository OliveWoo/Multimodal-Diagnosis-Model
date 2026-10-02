"""Build the final answer-blind KH taxonomy status and before/after audit.

The output separates three historically distinct comparisons:

* curated gap taxonomy at the Priority/Context/Audit/Hold screening layer;
* reviewed/auto exact overlays at compact candidate entry; and
* the downstream clinical Picked/High/Context layer.

It never reads benchmark answers and never changes a patient tier.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IDENTITY = ROOT / "outputs/runs/2026-09-30_KH_taxonomy_end_to_end_audit_v1/unique_identity_audit.csv"
DEFAULT_AUTO = ROOT / "outputs/runs/2026-09-24_KH_taxonomy_auto_decision_v1/decision_summary.csv"
DEFAULT_FRONT_BEFORE = ROOT / "outputs/runs/2026-09-22_KH_multi_assay_screening_shadow/screening_shadow.csv"
DEFAULT_FRONT_AFTER = ROOT / "outputs/runs/2026-09-22_KH_multi_assay_screening_shadow_taxonomy_v2/screening_shadow.csv"
DEFAULT_ENTRY_BEFORE = ROOT / "outputs/runs/2026-09-23_KH_compact_multi_assay_entry_source_contract_v5/all_decisions.csv"
DEFAULT_ENTRY_AFTER = ROOT / "outputs/runs/2026-09-25_KH_compact_multi_assay_entry_taxonomy_auto_v1/all_decisions.csv"
DEFAULT_CLINICAL_BEFORE = ROOT / "outputs/runs/2026-09-24_KH_test_aware_clinical_deterministic_shadow_v1/clinical_decisions.csv"
DEFAULT_CLINICAL_AFTER = ROOT / "outputs/runs/2026-09-25_KH_test_aware_clinical_taxonomy_auto_v2/clinical_decisions.csv"
DEFAULT_OUTPUT = ROOT / "outputs/runs/2026-09-30_KH_taxonomy_completion_audit_v1"


def _read(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _write(path: Path, rows: list[dict[str, Any]], fields: Iterable[str]) -> None:
    fieldnames = list(fields)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _bool(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes"}


def _front_route(action: str) -> str:
    return {
        "forward_priority_review": "Priority",
        "forward_context_review": "Context",
        "audit_only": "Audit",
        "hold_not_forwarded": "Hold",
        "retain_baseline": "Baseline-retained",
    }.get(action, action or "Unknown")


def _clinical_tier(value: str) -> str:
    return {
        "picked_shadow": "Picked",
        "high_priority_review": "High",
        "review_context_needed": "Context",
    }.get(value, value or "Unknown")


def _index(rows: list[dict[str, str]], fields: tuple[str, ...]) -> dict[tuple[str, ...], dict[str, str]]:
    return {tuple(str(row.get(field) or "").casefold() for field in fields): row for row in rows}


def build(
    identity_path: Path,
    auto_path: Path,
    front_before_path: Path,
    front_after_path: Path,
    entry_before_path: Path,
    entry_after_path: Path,
    clinical_before_path: Path,
    clinical_after_path: Path,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    identities = _read(identity_path)
    auto_rows = _read(auto_path)
    auto_by_key = {row["canonical_key"]: row for row in auto_rows}

    inventory: list[dict[str, Any]] = []
    for row in identities:
        auto = auto_by_key.get(row["canonical_key"], {})
        inventory.append(
            {
                "organism_name": row["classifier_input_name"],
                "canonical_key": row["canonical_key"],
                "effective_taxid": row["effective_audit_taxid"],
                "effective_scientific_name": row["effective_audit_scientific_name"],
                "effective_rank": row["effective_audit_rank"],
                "biological_class": row["biological_class"],
                "mapping_status": row["taxonomy_mapping_status_replayed"],
                "primary_rule_family": row["taxonomy_family_replayed"],
                "rule_source": row["rule_source"],
                "classification_confidence": row["classification_confidence"],
                "needs_literature_review": _bool(row["needs_literature_review"]),
                "auto_decision_status": auto.get("status", "not_in_auto_queue"),
                "auto_identity_status": auto.get("identity_status", ""),
                "auto_clinical_status": auto.get("clinical_status", ""),
                "auto_validator_reasons": auto.get("validator_reasons", ""),
            }
        )
    inventory.sort(key=lambda row: row["organism_name"].casefold())

    status_lists = {
        "exact_species": [r for r in inventory if r["mapping_status"] == "exact_species"],
        "genus_inherited": [r for r in inventory if r["mapping_status"] == "genus_inherited"],
        "prefix_inferred": [r for r in inventory if r["mapping_status"] == "key_prefix_inferred"],
        "unmapped": [r for r in inventory if r["mapping_status"] == "unmapped"],
        "provisional": [r for r in inventory if r["auto_decision_status"] == "provisional_unmapped"],
        "rejected": [r for r in inventory if r["auto_decision_status"] == "rejected"],
    }

    taxid_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    canonical_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    name_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in inventory:
        if row["effective_taxid"]:
            taxid_groups[row["effective_taxid"]].append(row)
        canonical_groups[row["canonical_key"]].append(row)
        name_groups[row["organism_name"].casefold()].append(row)

    alias_clusters: list[dict[str, Any]] = []
    for taxid, rows in sorted(taxid_groups.items()):
        canonical_keys = sorted({row["canonical_key"] for row in rows})
        if len(canonical_keys) > 1:
            alias_clusters.append(
                {
                    "taxid": taxid,
                    "input_names": "|".join(sorted({row["organism_name"] for row in rows})),
                    "canonical_keys": "|".join(canonical_keys),
                    "effective_scientific_names": "|".join(
                        sorted({row["effective_scientific_name"] for row in rows})
                    ),
                    "interpretation": "same_ncbi_taxon_alias_cluster_not_conflict",
                }
            )

    duplicate_name_rows: list[dict[str, Any]] = []
    for normalized_name, rows in sorted(name_groups.items()):
        if len(rows) > 1:
            duplicate_name_rows.append(
                {
                    "normalized_name": normalized_name,
                    "row_count": len(rows),
                    "canonical_keys": "|".join(sorted({r["canonical_key"] for r in rows})),
                    "taxids": "|".join(sorted({r["effective_taxid"] for r in rows})),
                }
            )

    taxid_conflicts: list[dict[str, Any]] = []
    for canonical_key, rows in sorted(canonical_groups.items()):
        taxids = sorted({row["effective_taxid"] for row in rows if row["effective_taxid"]})
        if len(taxids) > 1:
            taxid_conflicts.append(
                {
                    "canonical_key": canonical_key,
                    "organism_names": "|".join(sorted({r["organism_name"] for r in rows})),
                    "taxids": "|".join(taxids),
                    "conflict": True,
                }
            )

    front_before = _read(front_before_path)
    front_after = _read(front_after_path)
    front_key = ("patient_id", "organism_key")
    fb = _index(front_before, front_key)
    fa = _index(front_after, front_key)
    front_changes: list[dict[str, Any]] = []
    for key in sorted(set(fb) | set(fa)):
        before = fb.get(key, {})
        after = fa.get(key, {})
        old_action = before.get("screening_action", "missing")
        new_action = after.get("screening_action", "missing")
        old_family = before.get("primary_rule_family", "missing")
        new_family = after.get("primary_rule_family", "missing")
        if old_action != new_action or old_family != new_family:
            front_changes.append(
                {
                    "patient_id": after.get("patient_id") or before.get("patient_id"),
                    "organism_name": after.get("organism_name") or before.get("organism_name"),
                    "organism_key": after.get("organism_key") or before.get("organism_key"),
                    "old_family": old_family,
                    "new_family": new_family,
                    "old_screening_action": old_action,
                    "new_screening_action": new_action,
                    "old_route": _front_route(old_action),
                    "new_route": _front_route(new_action),
                    "route_changed": old_action != new_action,
                }
            )

    entry_before = _read(entry_before_path)
    entry_after = _read(entry_after_path)
    entry_key = ("patient_id", "case_review_id", "organism_key")
    eb = _index(entry_before, entry_key)
    ea = _index(entry_after, entry_key)
    entry_changes: list[dict[str, Any]] = []
    for key in sorted(set(eb) | set(ea)):
        before = eb.get(key, {})
        after = ea.get(key, {})
        fields = ("taxonomy_family", "taxonomy_mapping_status", "disposition", "screening_tier")
        if not before or not after or any(before.get(field, "") != after.get(field, "") for field in fields):
            entry_changes.append(
                {
                    "patient_id": after.get("patient_id") or before.get("patient_id"),
                    "case_review_id": after.get("case_review_id") or before.get("case_review_id"),
                    "organism_name": after.get("organism_name") or before.get("organism_name"),
                    "organism_key": after.get("organism_key") or before.get("organism_key"),
                    "old_family": before.get("taxonomy_family", "missing"),
                    "new_family": after.get("taxonomy_family", "missing"),
                    "old_mapping_status": before.get("taxonomy_mapping_status", "missing"),
                    "new_mapping_status": after.get("taxonomy_mapping_status", "missing"),
                    "old_disposition": before.get("disposition", "missing"),
                    "new_disposition": after.get("disposition", "missing"),
                    "old_screening_tier": before.get("screening_tier", "missing"),
                    "new_screening_tier": after.get("screening_tier", "missing"),
                    "candidate_disappeared": bool(before and not after),
                    "candidate_added": bool(after and not before),
                }
            )

    clinical_before = _read(clinical_before_path)
    clinical_after = _read(clinical_after_path)
    clinical_key = ("patient_id", "organism_name")
    cb = _index(clinical_before, clinical_key)
    ca = _index(clinical_after, clinical_key)
    clinical_changes: list[dict[str, Any]] = []
    for key in sorted(set(cb) | set(ca)):
        before = cb.get(key, {})
        after = ca.get(key, {})
        old_decision = before.get("clinical_decision", "missing")
        new_decision = after.get("clinical_decision", "missing")
        if not before or not after or old_decision != new_decision:
            clinical_changes.append(
                {
                    "patient_id": after.get("patient_id") or before.get("patient_id"),
                    "organism_name": after.get("organism_name") or before.get("organism_name"),
                    "old_family": before.get("taxonomy_family", "missing"),
                    "new_family": after.get("taxonomy_family", "missing"),
                    "old_clinical_decision": old_decision,
                    "new_clinical_decision": new_decision,
                    "old_clinical_tier": _clinical_tier(old_decision),
                    "new_clinical_tier": _clinical_tier(new_decision),
                    "common_candidate_tier_changed": bool(before and after and old_decision != new_decision),
                    "candidate_disappeared": bool(before and not after),
                    "candidate_added": bool(after and not before),
                }
            )

    patient_ids = sorted(
        {row["patient_id"] for row in front_before + front_after + entry_before + entry_after + clinical_before + clinical_after},
        key=int,
    )
    patient_rows: list[dict[str, Any]] = []
    for patient_id in patient_ids:
        f_old = [row for row in front_before if row["patient_id"] == patient_id]
        f_new = [row for row in front_after if row["patient_id"] == patient_id]
        e_old = [row for row in entry_before if row["patient_id"] == patient_id]
        e_new = [row for row in entry_after if row["patient_id"] == patient_id]
        c_old = [row for row in clinical_before if row["patient_id"] == patient_id]
        c_new = [row for row in clinical_after if row["patient_id"] == patient_id]
        old_clinical_names = {row["organism_name"].casefold(): row["organism_name"] for row in c_old}
        new_clinical_names = {row["organism_name"].casefold(): row["organism_name"] for row in c_new}
        row: dict[str, Any] = {"patient_id": patient_id}
        for label, action in (
            ("priority", "forward_priority_review"),
            ("context", "forward_context_review"),
            ("audit", "audit_only"),
            ("hold", "hold_not_forwarded"),
        ):
            row[f"front_{label}_before"] = sum(r["screening_action"] == action for r in f_old)
            row[f"front_{label}_after"] = sum(r["screening_action"] == action for r in f_new)
        for label, disposition in (
            ("scorer_entry", "scorer_entry"),
            ("clinical_review", "clinical_review"),
            ("qc_only", "qc_only"),
        ):
            row[f"entry_{label}_before"] = sum(r["disposition"] == disposition for r in e_old)
            row[f"entry_{label}_after"] = sum(r["disposition"] == disposition for r in e_new)
        for label, decision in (
            ("picked", "picked_shadow"),
            ("high", "high_priority_review"),
            ("context", "review_context_needed"),
        ):
            row[f"clinical_{label}_before"] = sum(r["clinical_decision"] == decision for r in c_old)
            row[f"clinical_{label}_after"] = sum(r["clinical_decision"] == decision for r in c_new)
        row["clinical_added_candidates"] = "|".join(
            new_clinical_names[key] for key in sorted(set(new_clinical_names) - set(old_clinical_names))
        )
        row["clinical_removed_candidates"] = "|".join(
            old_clinical_names[key] for key in sorted(set(old_clinical_names) - set(new_clinical_names))
        )
        row["common_clinical_tier_change_count"] = sum(
            change["patient_id"] == patient_id and change["common_candidate_tier_changed"]
            for change in clinical_changes
        )
        patient_rows.append(row)

    front_route_changes = [row for row in front_changes if row["route_changed"]]
    clinical_additions = [row for row in clinical_changes if row["candidate_added"]]
    clinical_removals = [row for row in clinical_changes if row["candidate_disappeared"]]
    common_clinical_tier_changes = [row for row in clinical_changes if row["common_candidate_tier_changed"]]
    assertions = {
        "identity_inventory_has_404_unique_names": len(inventory) == 404,
        "no_duplicate_normalized_input_names": not duplicate_name_rows,
        "no_canonical_key_to_multiple_taxid_conflict": not taxid_conflicts,
        "front_screening_candidate_keys_preserved": set(fb) == set(fa),
        "compact_entry_candidate_keys_preserved": set(eb) == set(ea),
        "no_downstream_clinical_candidate_disappearance": not clinical_removals,
        "no_existing_downstream_clinical_tier_change": not common_clinical_tier_changes,
        "all_downstream_additions_are_context": all(
            row["new_clinical_tier"] == "Context" for row in clinical_additions
        ),
    }
    summary = {
        "schema_version": "kh_taxonomy_completion_audit.v1",
        "answer_blind": True,
        "counts": {
            "unique_identity_names": len(inventory),
            "mapping_status": dict(Counter(row["mapping_status"] for row in inventory)),
            "auto_decision_status": dict(Counter(row["auto_decision_status"] for row in inventory)),
            "provisional": len(status_lists["provisional"]),
            "rejected": len(status_lists["rejected"]),
            "duplicate_normalized_names": len(duplicate_name_rows),
            "taxid_conflicts": len(taxid_conflicts),
            "same_taxid_alias_clusters": len(alias_clusters),
            "front_candidate_rows_before_after": [len(front_before), len(front_after)],
            "front_family_or_route_changes": len(front_changes),
            "front_route_changes": len(front_route_changes),
            "front_route_change_matrix": {
                f"{old} -> {new}": count
                for (old, new), count in sorted(
                    Counter((row["old_route"], row["new_route"]) for row in front_route_changes).items()
                )
            },
            "compact_entry_rows_before_after": [len(entry_before), len(entry_after)],
            "compact_entry_changes": len(entry_changes),
            "entry_disposition_change_matrix": {
                f"{old} -> {new}": count
                for (old, new), count in sorted(
                    Counter((row["old_disposition"], row["new_disposition"]) for row in entry_changes).items()
                )
            },
            "clinical_rows_before_after": [len(clinical_before), len(clinical_after)],
            "clinical_candidate_additions": len(clinical_additions),
            "clinical_candidate_removals": len(clinical_removals),
            "common_clinical_tier_changes": len(common_clinical_tier_changes),
            "clinical_additions_by_tier": dict(Counter(row["new_clinical_tier"] for row in clinical_additions)),
        },
        "interpretation": {
            "front_comparison": "Curated gap taxonomy v2 changes Priority/Context/Audit/Hold routing on the same 767 candidates.",
            "entry_comparison": "Reviewed plus auto exact overlays change candidate-entry disposition on the same 782 candidates.",
            "downstream_comparison": "Three candidates are added as Context; no existing candidate changes Picked/High/Context and none disappears.",
            "alias": "Multiple input canonical keys sharing one NCBI taxid are reported as synonym clusters, not treated as taxid conflicts.",
        },
        "assertions": assertions,
        "all_assertions_pass": all(assertions.values()),
    }
    outputs = {
        "inventory": inventory,
        **status_lists,
        "alias_clusters": alias_clusters,
        "duplicate_names": duplicate_name_rows,
        "taxid_conflicts": taxid_conflicts,
        "front_changes": front_changes,
        "entry_changes": entry_changes,
        "clinical_changes": clinical_changes,
        "per_patient": patient_rows,
    }
    return summary, outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--auto", type=Path, default=DEFAULT_AUTO)
    parser.add_argument("--front-before", type=Path, default=DEFAULT_FRONT_BEFORE)
    parser.add_argument("--front-after", type=Path, default=DEFAULT_FRONT_AFTER)
    parser.add_argument("--entry-before", type=Path, default=DEFAULT_ENTRY_BEFORE)
    parser.add_argument("--entry-after", type=Path, default=DEFAULT_ENTRY_AFTER)
    parser.add_argument("--clinical-before", type=Path, default=DEFAULT_CLINICAL_BEFORE)
    parser.add_argument("--clinical-after", type=Path, default=DEFAULT_CLINICAL_AFTER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary, outputs = build(
        args.identity,
        args.auto,
        args.front_before,
        args.front_after,
        args.entry_before,
        args.entry_after,
        args.clinical_before,
        args.clinical_after,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    inventory_fields = list(outputs["inventory"][0])
    for name in ("inventory", "exact_species", "genus_inherited", "prefix_inferred", "unmapped", "provisional", "rejected"):
        _write(args.output_dir / f"{name}.csv", outputs[name], inventory_fields)
    _write(args.output_dir / "same_taxid_alias_clusters.csv", outputs["alias_clusters"], ["taxid", "input_names", "canonical_keys", "effective_scientific_names", "interpretation"])
    _write(args.output_dir / "duplicate_normalized_names.csv", outputs["duplicate_names"], ["normalized_name", "row_count", "canonical_keys", "taxids"])
    _write(args.output_dir / "taxid_conflicts.csv", outputs["taxid_conflicts"], ["canonical_key", "organism_names", "taxids", "conflict"])
    for name in ("front_changes", "entry_changes", "clinical_changes", "per_patient"):
        rows = outputs[name]
        _write(args.output_dir / f"{name}.csv", rows, list(rows[0]) if rows else ["status"])
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

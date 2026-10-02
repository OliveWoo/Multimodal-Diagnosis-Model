"""Audit the KH raw-name -> canonical identity -> taxid -> family -> route chain.

This audit is answer blind and non-mutating.  It replays the current taxonomy
classifier against every frozen case-organism candidate, checks downstream
metadata parity, and emits an identity-resolution queue without applying any
new taxonomy profile.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools import organism_taxonomy_classifier as taxonomy


DEFAULT_UNION = Path(
    "outputs/runs/2026-09-23_KH_multi_assay_mngs_inventory_source_contract_v3"
    "/candidate_union.csv"
)
DEFAULT_COMPACT = Path(
    "outputs/runs/2026-09-30_KH_integrated_release_v2_5_taxonomy_auto_review"
    "/01_compact_taxonomy"
    "/all_decisions.csv"
)
DEFAULT_SCORER = Path(
    "outputs/runs/2026-09-30_KH_integrated_release_v2_5_taxonomy_auto_review"
    "/02_analytical_scorer"
    "/case_candidate_decisions.csv"
)
DEFAULT_INTEGRATED = Path(
    "outputs/runs/2026-09-30_KH_integrated_release_v2_5_taxonomy_auto_review"
    "/integrated_decisions.csv"
)
DEFAULT_AUTO_DECISIONS = Path(
    "outputs/runs/2026-09-30_KH_taxonomy_auto_decision_v2_4/decision_summary.csv"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-30_KH_taxonomy_end_to_end_audit_v1"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def row_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (
        str(row.get("patient_id") or ""),
        str(row.get("case_review_id") or ""),
        str(row.get("organism_key") or ""),
    )


def patient_name_key(row: dict[str, Any]) -> tuple[str, str]:
    return str(row.get("patient_id") or ""), str(row.get("organism_name") or "")


def normalized_class(category: str) -> str:
    return taxonomy.normalize_biological_class(category)


def load_resolution(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.is_file():
        return {}
    return {row["input_name"].casefold(): row for row in read_csv(path)}


def audit(
    union_path: Path,
    compact_path: Path,
    scorer_path: Path,
    integrated_path: Path,
    auto_decisions_path: Path,
    output_dir: Path,
    *,
    resolution_path: Path | None = None,
    manual_resolution_path: Path | None = None,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    union_rows = read_csv(union_path)
    compact_rows = read_csv(compact_path)
    scorer_rows = read_csv(scorer_path)
    integrated_rows = read_csv(integrated_path)
    auto_rows = read_csv(auto_decisions_path)
    resolution = load_resolution(resolution_path)
    manual_resolution = load_resolution(manual_resolution_path)
    resolution.update(manual_resolution)

    union_index = {row_key(row): row for row in union_rows}
    compact_index = {row_key(row): row for row in compact_rows}
    scorer_index = {row_key(row): row for row in scorer_rows}
    integrated_index = {patient_name_key(row): row for row in integrated_rows}

    trace_rows: list[dict[str, Any]] = []
    unique_profiles: dict[str, dict[str, Any]] = {}
    normalized_groups: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for compact in compact_rows:
        source = union_index.get(row_key(compact), {})
        profile = taxonomy.classify_organism(
            compact["organism_name"], biological_class=source.get("category")
        )
        resolved = resolution.get(compact["organism_name"].casefold(), {})
        effective_taxid = profile.get("taxid") or resolved.get("ncbi_taxid") or ""
        effective_name = (
            profile.get("display_name")
            if profile.get("taxid")
            else resolved.get("ncbi_scientific_name") or profile.get("display_name")
        )
        effective_rank = (
            profile.get("taxonomic_rank")
            if profile.get("taxid")
            else resolved.get("ncbi_rank") or profile.get("taxonomic_rank")
        )
        scorer = scorer_index.get(row_key(compact), {})
        integrated = integrated_index.get(patient_name_key(compact), {})
        family_matches = compact.get("taxonomy_family") == profile.get(
            "primary_rule_family"
        )
        status_matches = compact.get("taxonomy_mapping_status") == profile.get(
            "mapping_status"
        )
        trace = {
            "patient_id": compact.get("patient_id"),
            "case_review_id": compact.get("case_review_id"),
            "specimen_code": compact.get("specimen_code"),
            "raw_organism_name": source.get("organism_name"),
            "source_organism_key": source.get("organism_key"),
            "source_category": source.get("category"),
            "classifier_input_name": compact.get("organism_name"),
            "canonical_display_name": profile.get("display_name"),
            "canonical_key": profile.get("canonical_key"),
            "canonical_key_changed_from_source": (
                source.get("organism_key") != profile.get("canonical_key")
            ),
            "classifier_taxid": profile.get("taxid") or "",
            "resolver_status": resolved.get("resolution_status") or "not_run",
            "resolver_taxid": resolved.get("ncbi_taxid") or "",
            "resolver_scientific_name": resolved.get("ncbi_scientific_name") or "",
            "resolver_source": resolved.get("source") or "",
            "effective_audit_taxid": effective_taxid,
            "effective_audit_scientific_name": effective_name,
            "effective_audit_rank": effective_rank,
            "taxonomy_mapping_status_frozen": compact.get("taxonomy_mapping_status"),
            "taxonomy_mapping_status_replayed": profile.get("mapping_status"),
            "taxonomy_family_frozen": compact.get("taxonomy_family"),
            "taxonomy_family_replayed": profile.get("primary_rule_family"),
            "family_replay_matches": family_matches,
            "mapping_status_replay_matches": status_matches,
            "biological_class": profile.get("biological_class"),
            "biological_class_conflict": profile.get("biological_class_conflict"),
            "rule_source": profile.get("rule_source"),
            "classification_confidence": profile.get("classification_confidence"),
            "needs_literature_review": profile.get("needs_literature_review"),
            "compact_disposition": compact.get("disposition"),
            "reached_scorer": bool(scorer),
            "scorer_taxonomy_family": scorer.get("taxonomy_family") or "",
            "scorer_family_matches_compact": (
                not scorer
                or scorer.get("taxonomy_family") == compact.get("taxonomy_family")
            ),
            "reached_integrated_patient_candidate": bool(integrated),
            "integrated_taxonomy_family": integrated.get("taxonomy_family") or "",
            "integrated_family_matches_compact": (
                not integrated
                or integrated.get("taxonomy_family") == compact.get("taxonomy_family")
            ),
            "routing_only": profile.get("routing_only"),
            "answer_blind": True,
        }
        trace_rows.append(trace)
        key = str(compact.get("organism_key") or "")
        unique_profiles.setdefault(key, trace)
        normalized_groups[str(profile.get("canonical_key") or "")].add(
            (str(source.get("organism_key") or ""), str(compact.get("organism_name") or ""))
        )

    unique_rows = sorted(
        unique_profiles.values(), key=lambda row: str(row["classifier_input_name"]).casefold()
    )
    resolution_queue = [
        {
            "organism_name": row["classifier_input_name"],
            "category": row["source_category"],
            "current_mapping_status": row["taxonomy_mapping_status_replayed"],
            "current_family": row["taxonomy_family_replayed"],
            "current_taxid": row["classifier_taxid"],
            "reason": "species_identity_taxid_missing",
        }
        for row in unique_rows
        if not row["classifier_taxid"]
    ]
    drift_rows = [
        row
        for row in trace_rows
        if not row["family_replay_matches"] or not row["mapping_status_replay_matches"]
    ]
    collisions = [
        {
            "canonical_key": canonical,
            "source_variants": json.dumps(sorted(values), ensure_ascii=False),
            "variant_count": len(values),
        }
        for canonical, values in sorted(normalized_groups.items())
        if len({value[0] for value in values}) > 1
    ]
    resolver_counts = Counter(
        row["resolver_status"] for row in unique_rows if row["resolver_status"] != "not_run"
    )
    auto_counts = Counter(row.get("status") or "unknown" for row in auto_rows)

    write_csv(output_dir / "taxonomy_trace.csv", trace_rows, list(trace_rows[0]))
    write_csv(
        output_dir / "unique_identity_audit.csv", unique_rows, list(unique_rows[0])
    )
    write_csv(
        output_dir / "taxid_resolution_queue.csv",
        resolution_queue,
        list(resolution_queue[0]) if resolution_queue else ["organism_name"],
    )
    write_csv(
        output_dir / "taxonomy_replay_drift.csv",
        drift_rows,
        list(trace_rows[0]),
    )
    write_csv(
        output_dir / "canonical_collision_audit.csv",
        collisions,
        ["canonical_key", "source_variants", "variant_count"],
    )

    assertions = {
        "inventory_and_compact_have_782_rows": len(union_rows)
        == len(compact_rows)
        == 782,
        "inventory_to_compact_keys_are_one_to_one": set(union_index)
        == set(compact_index),
        "all_scorer_entries_trace_to_compact": (
            len(scorer_rows)
            == sum(row.get("disposition") == "scorer_entry" for row in compact_rows)
            and len(scorer_index) == len(scorer_rows)
            and set(scorer_index).issubset(set(compact_index))
        ),
        "no_biological_class_conflicts": not any(
            row["biological_class_conflict"] for row in trace_rows
        ),
        "no_taxonomy_family_replay_drift": not any(
            not row["family_replay_matches"] for row in trace_rows
        ),
        "scorer_family_matches_compact": all(
            row["scorer_family_matches_compact"] for row in trace_rows
        ),
        "integrated_family_matches_compact_when_present": all(
            row["integrated_family_matches_compact"] for row in trace_rows
        ),
        "taxonomy_is_routing_only": all(row["routing_only"] for row in trace_rows),
        "no_cross_source_canonical_collisions": not collisions,
    }
    result = {
        "schema_version": "kh_taxonomy_end_to_end_audit.v1",
        "answer_blind": True,
        "scope": (
            "Raw candidate name to canonical identity, optional NCBI taxid, central "
            "family routing, scorer, and integrated metadata. No clinical tier changes."
        ),
        "inputs": {
            "candidate_union": str(union_path.resolve()),
            "compact_decisions": str(compact_path.resolve()),
            "scorer_decisions": str(scorer_path.resolve()),
            "integrated_decisions": str(integrated_path.resolve()),
            "auto_decisions": str(auto_decisions_path.resolve()),
            "ncbi_resolution": str(resolution_path.resolve())
            if resolution_path and resolution_path.exists()
            else None,
            "manual_ncbi_alias_resolution": str(manual_resolution_path.resolve())
            if manual_resolution_path and manual_resolution_path.exists()
            else None,
        },
        "hashes": {
            "candidate_union": sha256_file(union_path),
            "compact_decisions": sha256_file(compact_path),
            "scorer_decisions": sha256_file(scorer_path),
            "integrated_decisions": sha256_file(integrated_path),
        },
        "counts": {
            "case_organism_rows": len(trace_rows),
            "unique_source_organism_keys": len(unique_rows),
            "mapping_status_unique": dict(
                sorted(Counter(row["taxonomy_mapping_status_replayed"] for row in unique_rows).items())
            ),
            "family_unique": dict(
                sorted(Counter(row["taxonomy_family_replayed"] for row in unique_rows).items())
            ),
            "classifier_taxid_present_unique": sum(
                bool(row["classifier_taxid"]) for row in unique_rows
            ),
            "classifier_taxid_missing_unique": len(resolution_queue),
            "effective_taxid_present_unique": sum(
                bool(row["effective_audit_taxid"]) for row in unique_rows
            ),
            "resolver_status_unique": dict(sorted(resolver_counts.items())),
            "auto_review_status": dict(sorted(auto_counts.items())),
            "family_replay_drift_rows": sum(
                not row["family_replay_matches"] for row in trace_rows
            ),
            "mapping_status_replay_drift_rows": sum(
                not row["mapping_status_replay_matches"] for row in trace_rows
            ),
            "canonical_key_changed_from_source_rows": sum(
                row["canonical_key_changed_from_source"] for row in trace_rows
            ),
            "canonical_collision_count": len(collisions),
            "scorer_rows": len(scorer_rows),
            "integrated_rows": len(integrated_rows),
        },
        "interpretation": {
            "family_routing": (
                "Current family routing is replayed independently from frozen compact "
                "metadata; mapping-status-only drift is reported separately."
            ),
            "taxid": (
                "A missing classifier taxid does not mean identity was disproven. It is "
                "queued for NCBI resolution and must not change clinical tier by itself."
            ),
            "provisional": (
                "Provisional/unmapped organisms remain visible for review and are not "
                "silently dropped or auto-promoted."
            ),
        },
        "assertions": assertions,
        "assertion_pass_count": sum(assertions.values()),
        "assertion_total_count": len(assertions),
        "all_required_assertions_pass": all(assertions.values()),
    }
    write_json(output_dir / "summary.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--union", type=Path, default=DEFAULT_UNION)
    parser.add_argument("--compact", type=Path, default=DEFAULT_COMPACT)
    parser.add_argument("--scorer", type=Path, default=DEFAULT_SCORER)
    parser.add_argument("--integrated", type=Path, default=DEFAULT_INTEGRATED)
    parser.add_argument("--auto-decisions", type=Path, default=DEFAULT_AUTO_DECISIONS)
    parser.add_argument("--ncbi-resolution", type=Path)
    parser.add_argument("--manual-ncbi-alias-resolution", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = audit(
        args.union,
        args.compact,
        args.scorer,
        args.integrated,
        args.auto_decisions,
        args.output_dir,
        resolution_path=args.ncbi_resolution,
        manual_resolution_path=args.manual_ncbi_alias_resolution,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

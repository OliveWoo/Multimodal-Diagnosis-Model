"""Run the legacy deterministic scorer on test-aware multi-assay mNGS data.

The adapter creates one legacy ranked record per real seq_id test. It uses the
frozen selected-DNA hospital summary and deliberately omits phenotype. The
legacy scorer then deduplicates repeated organism observations as it normally
does, exposing information lost when a non-test-aware model consumes the new
data. Existing patient outputs are never modified.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools import deterministic_mngs_max_scorer as legacy
from tools.build_index_event_clinical_timeline_shadow import sha256_file
from tools.pathogen_normalization import canonical_key


SCHEMA_VERSION = "legacy_scorer_multi_assay_compatibility_shadow.v1"
DEFAULT_SUMMARY_SUFFIX = "final_summary_selected_dna_20260918"
TAXONOMY_MODES = ("base-only", "current-overlays")
DEFAULT_TAXONOMY_OVERLAY_PATHS = tuple(legacy.organism_taxonomy.RULE_OVERLAY_PATHS)
RANK_RULES = {
    1: "Code_NTC=00 + Code_RK_NTC=RK00",
    2: "Code_RK_NTC=RK8 + Code_NTC=00",
    3: "Code_RK_NTC=RK00",
    4: "Code_RK_NTC=RK8",
    5: "Code_NTC=00",
}
SOURCE_CATEGORY = {
    "bacterial": "1.Bac",
    "fungal": "2.Fungi",
    "viral": "3.Virus",
    "parasitic": "4.Parasite",
}


def configure_taxonomy_mode(mode: str) -> dict[str, Any]:
    """Make the legacy scorer's taxonomy dependency explicit and reproducible.

    The historical A comparator predates the later gap-resolution, automatic,
    and reviewed overlays.  A base-only run therefore has to opt out of those
    overlays even when the files now exist in the working tree.
    """
    if mode not in TAXONOMY_MODES:
        raise ValueError(f"Unknown taxonomy mode: {mode}")
    taxonomy = legacy.organism_taxonomy
    taxonomy.RULE_OVERLAY_PATHS = (
        () if mode == "base-only" else DEFAULT_TAXONOMY_OVERLAY_PATHS
    )
    taxonomy.rules_payload.cache_clear()
    payload = taxonomy.rules_payload()
    taxonomy.PROFILE_VERSION = str(
        payload.get("profile_version") or "organism_taxonomy_profile_v1"
    )
    taxonomy.RULE_SCHEMA_VERSION = str(
        payload.get("schema_version") or "organism_taxonomy_rules_v1"
    )
    loaded_overlays = [str(path) for path in taxonomy.RULE_OVERLAY_PATHS if path.exists()]
    return {
        "mode": mode,
        "base_rules_path": str(taxonomy.RULES_PATH.resolve()),
        "base_rules_sha256": sha256_file(taxonomy.RULES_PATH),
        "profile_version": taxonomy.PROFILE_VERSION,
        "rule_schema_version": taxonomy.RULE_SCHEMA_VERSION,
        "loaded_overlay_paths": loaded_overlays,
        "loaded_overlay_sha256": {
            str(path.resolve()): sha256_file(path)
            for path in taxonomy.RULE_OVERLAY_PATHS
            if path.exists()
        },
    }


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot identify patient from {path}")
    return int(match.group(1))


def as_int(value: Any, default: int = 999) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def possibility(rank_priority: int) -> str:
    if rank_priority == 1:
        return "high"
    if rank_priority in {2, 3}:
        return "medium"
    return "low"


def empirical_reads_percentile(reads: float, all_reads: list[float]) -> float:
    if not all_reads:
        return 0.0
    return round(sum(value <= reads for value in all_reads) / len(all_reads), 6)


def build_legacy_ranked_payload(packet: dict[str, Any]) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    records: dict[tuple[str, str], dict[str, Any]] = {}
    provenance: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for case in packet.get("cases") or []:
        case_id = str(case.get("case_review_id") or case.get("specimen_code") or "")
        for entry in case.get("scorer_entries") or []:
            organism = str(entry.get("organism_name") or "").strip()
            if not organism:
                continue
            for signal in entry.get("selected_test_signals") or []:
                reads = float(signal.get("reads") or 0)
                if reads <= 0:
                    raise ValueError(f"Nonpositive reads in compact scorer entry: {organism}")
                seq_id = str(signal.get("seq_id") or "").strip()
                if not seq_id:
                    raise ValueError(f"Missing seq_id for {organism}")
                key = (case_id, seq_id)
                record = records.setdefault(key, {
                    "specimen_code": case.get("specimen_code") or case_id,
                    "collected_time": case.get("collected_time"),
                    "specimen_site": case.get("specimen_site"),
                    "seq_id": seq_id,
                    "nucleic_type": signal.get("nucleic_type"),
                    "status": "multi_assay_compatibility_adapter",
                    "candidates": [],
                })
                rank_priority = as_int(signal.get("analytical_rank_priority"), 5)
                candidate = {
                    "organism_name": organism,
                    "source_category": SOURCE_CATEGORY.get(str(entry.get("category") or "").lower(), "Unknown"),
                    "reads": reads,
                    "sec_hit": reads,
                    "selection_reasons": list(entry.get("rule_ids") or []),
                    "ranking": {
                        "rank_priority": rank_priority,
                        "rank_rule": RANK_RULES.get(rank_priority, "multi_assay_rank_unavailable"),
                        "reads_percentile": None,
                        "possibility_level": possibility(rank_priority),
                    },
                    "compatibility_source": {
                        "case_review_id": case_id,
                        "seq_id": seq_id,
                        "nucleic_type": signal.get("nucleic_type"),
                        "rank_in_retained_universe": signal.get("rank_in_retained_universe"),
                        "rpm_total": signal.get("rpm_total"),
                        "qc_status": signal.get("qc_status"),
                        "normalization_available": signal.get("normalization_available"),
                    },
                }
                record["candidates"].append(candidate)
                provenance[canonical_key(organism)].append(candidate["compatibility_source"])

    ordered_records = []
    for key in sorted(records):
        record = records[key]
        reads_values = [float(item["reads"]) for item in record["candidates"]]
        for candidate in record["candidates"]:
            candidate["ranking"]["reads_percentile"] = empirical_reads_percentile(
                float(candidate["reads"]), reads_values
            )
        record["candidates"].sort(key=lambda item: (
            as_int(item["ranking"]["rank_priority"]),
            -float(item["ranking"]["reads_percentile"]),
            -float(item["reads"]),
            canonical_key(item["organism_name"]),
        ))
        ordered_records.append(record)
    payload = {
        "patient_id": str(packet.get("patient_id")),
        "records": ordered_records,
        "ranking_metadata": {
            "source": "compact_multi_assay_entry_shadow_v2",
            "adapter": SCHEMA_VERSION,
            "reads_percentile": "empirical CDF within each real seq_id record",
            "no_cross_test_read_sum": True,
            "known_information_loss": [
                "The legacy scorer has no DNA/RNA concordance or technical-repeat feature.",
                "The legacy scorer deduplicates repeated organism observations after per-test scoring.",
                "The legacy scorer does not consume RPM, normalization availability, or phenotype.",
            ],
        },
    }
    return payload, dict(provenance)


def picked_names(score: dict[str, Any]) -> set[str]:
    return {
        canonical_key(item.get("organism_name"))
        for item in (score.get("best_available_summary") or {}).get("picked_pathogens") or []
        if isinstance(item, dict)
    }


def score_rows(patient: int, score: dict[str, Any]) -> list[dict[str, Any]]:
    picked = picked_names(score)
    rows = []
    for candidate in score.get("pathogen_candidates") or []:
        if not isinstance(candidate, dict):
            continue
        rows.append({
            "patient_id": patient,
            "organism_name": candidate.get("organism_name"),
            "evidence_source": candidate.get("evidence_source", "mNGS_ranked"),
            "integrated_level": candidate.get("integrated_causative_level"),
            "picked": canonical_key(candidate.get("organism_name")) in picked,
            "rank_priority": candidate.get("rank_priority"),
            "reads": candidate.get("reads"),
            "reads_percentile": candidate.get("reads_percentile"),
            "dominance_tier": candidate.get("dominance_tier"),
            "specimen_class": candidate.get("specimen_class"),
            "formal_pick_exclusion_rule": (
                (candidate.get("formal_pick_exclusion_rule") or {}).get("rule_id")
                if isinstance(candidate.get("formal_pick_exclusion_rule"), dict) else ""
            ),
            "source_test_count": len(candidate.get("compatibility_source_tests") or []),
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0]) if rows else ["patient_id", "organism_name"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(
    compact_root: Path, patient_root: Path, output_dir: Path,
    *, summary_suffix: str = DEFAULT_SUMMARY_SUFFIX,
    taxonomy_mode: str = "base-only",
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    taxonomy_snapshot = configure_taxonomy_mode(taxonomy_mode)
    packet_paths = sorted((compact_root / "patient_packets").glob(
        "NGS_patient_*_compact_entry_shadow.json"), key=patient_number)
    if not packet_paths:
        raise ValueError(f"No compact multi-assay packets in {compact_root}")
    output_dir.mkdir(parents=True, exist_ok=True)
    ranked_dir = output_dir / "legacy_ranked_inputs"
    score_dir = output_dir / "patient_outputs"
    ranked_dir.mkdir()
    score_dir.mkdir()
    counts: Counter[str] = Counter()
    rows: list[dict[str, Any]] = []
    for packet_path in packet_paths:
        patient = patient_number(packet_path)
        patient_dir = patient_root / f"NGS_patient_{patient}_json"
        summary_path = patient_dir / "summary_outputs" / f"NGS_patient_{patient}_{summary_suffix}.json"
        if not summary_path.is_file():
            raise FileNotFoundError(summary_path)
        packet = read_json(packet_path)
        ranked, provenance = build_legacy_ranked_payload(packet)
        final_summary = read_json(summary_path)
        if any("phenotype" in str(key).casefold() for key in final_summary):
            raise ValueError(f"Phenotype key found in frozen hospital summary: {summary_path}")
        score = legacy.score_payload(
            patient_dir=patient_dir,
            ranked_mngs=ranked,
            final_summary=final_summary,
        )
        for candidate in score.get("pathogen_candidates") or []:
            if not isinstance(candidate, dict):
                continue
            candidate["compatibility_source_tests"] = provenance.get(
                canonical_key(candidate.get("organism_name")), []
            )
        score.update({
            "compatibility_schema_version": SCHEMA_VERSION,
            "experiment": "new_multi_assay_mngs_plus_frozen_old_hospital_data_through_legacy_scorer",
            "phenotype_used": False,
            "source_files": {
                "compact_packet": str(packet_path.resolve()),
                "compact_packet_sha256": sha256_file(packet_path),
                "frozen_hospital_summary": str(summary_path.resolve()),
                "frozen_hospital_summary_sha256": sha256_file(summary_path),
            },
            "adapter_constraints": ranked["ranking_metadata"]["known_information_loss"],
        })
        ranked_path = ranked_dir / f"NGS_patient_{patient}_legacy_ranked_multi_assay_compat.json"
        score_path = score_dir / f"NGS_patient_{patient}_legacy_scorer_multi_assay_shadow.json"
        write_json(ranked_path, ranked)
        write_json(score_path, score)
        patient_rows = score_rows(patient, score)
        rows.extend(patient_rows)
        counts["patients"] += 1
        counts["legacy_test_records"] += len(ranked["records"])
        counts["positive_test_signals"] += sum(len(record["candidates"]) for record in ranked["records"])
        counts["scored_patient_organisms"] += len(patient_rows)
        counts["picked"] += len(picked_names(score))
        counts["hospital_only"] += sum(row["evidence_source"] == "hospital_only" for row in patient_rows)
        for row in patient_rows:
            counts[f"level:{row['integrated_level']}"] += 1

    write_csv(output_dir / "legacy_multi_assay_candidate_decisions.csv", rows)
    summary = {
        "schema_version": SCHEMA_VERSION,
        "experiment": "new multi-assay mNGS + frozen old hospital data + legacy deterministic scorer",
        "answer_blind": True,
        "phenotype_used": False,
        "legacy_rule_version": legacy.RANKED_RULE_VERSION,
        "taxonomy_snapshot": taxonomy_snapshot,
        "code_snapshot": {
            "legacy_scorer_path": str(Path(legacy.__file__).resolve()),
            "legacy_scorer_sha256": sha256_file(Path(legacy.__file__)),
            "adapter_path": str(Path(__file__).resolve()),
            "adapter_sha256": sha256_file(Path(__file__)),
        },
        "counts": dict(sorted(counts.items())),
        "compact_root": str(compact_root.resolve()),
        "patient_root": str(patient_root.resolve()),
        "hospital_summary_suffix": summary_suffix,
        "constraints": [
            "Each seq_id remains a separate legacy ranked record; reads are never summed across tests.",
            "The legacy scorer receives no phenotype fields.",
            "The adapter reconstructs per-test reads percentile because the compact source does not preserve the legacy percentile.",
            "The legacy scorer cannot use DNA/RNA concordance, repeated-test support, RPM, or normalization availability.",
            "Base-only taxonomy mode excludes all later taxonomy overlays so the scorer-side taxonomy factor stays old.",
            "This is a compatibility shadow, not a production or historically exact replay."
        ],
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("compact_root", type=Path)
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--summary-suffix", default=DEFAULT_SUMMARY_SUFFIX)
    parser.add_argument("--taxonomy-mode", choices=TAXONOMY_MODES, default="base-only")
    args = parser.parse_args()
    print(json.dumps(run(
        args.compact_root, args.patient_root, args.output_dir,
        summary_suffix=args.summary_suffix,
        taxonomy_mode=args.taxonomy_mode,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

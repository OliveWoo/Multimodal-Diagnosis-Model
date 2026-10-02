from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from tools.organism_taxonomy_classifier import classify_organism


DEFAULT_POLICY = Path("rules/multi_assay_mngs_screening_policy.json")


def _as_int(value: Any, default: int | None = None) -> int | None:
    if value in (None, ""):
        return default
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def _as_float(value: Any, default: float | None = None) -> float | None:
    if value in (None, ""):
        return default
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return default


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _min_present(values: Iterable[int | None]) -> int | None:
    present = [value for value in values if value is not None]
    return min(present) if present else None


def _max_present(values: Iterable[float | int | None]) -> float | int | None:
    present = [value for value in values if value is not None]
    return max(present) if present else None


@dataclass(frozen=True)
class EvidenceAxes:
    cross_molecule_support: bool
    repeated_detection: bool
    priority_rank_support: bool
    reproducibility_support: bool
    independent_dimension_count: int


def _evidence_axes(candidate: dict[str, Any], policy: dict[str, Any]) -> EvidenceAxes:
    rank_threshold = int(policy["strong_analytical_axes"]["best_priority_rank_at_most"])
    priority_rank = candidate.get("best_analytical_rank_priority")
    cross = bool(candidate.get("cross_molecule_support"))
    repeated = int(candidate.get("detected_test_count") or 0) >= 2
    rank_support = priority_rank is not None and priority_rank <= rank_threshold
    reproducibility = cross or repeated
    return EvidenceAxes(
        cross_molecule_support=cross,
        repeated_detection=repeated,
        priority_rank_support=rank_support,
        reproducibility_support=reproducibility,
        independent_dimension_count=sum((reproducibility, rank_support)),
    )


def _molecule_support(candidate: dict[str, Any]) -> str:
    dna = int(candidate.get("dna_detected_test_count") or 0)
    rna = int(candidate.get("rna_detected_test_count") or 0)
    if dna and rna:
        return "DNA+RNA"
    if dna:
        return "DNA-only"
    if rna:
        return "RNA-only"
    return "none"


def route_candidate(
    candidate: dict[str, Any],
    taxonomy: dict[str, Any],
    comparison_status: str,
    policy: dict[str, Any],
) -> tuple[str, list[str], list[str]]:
    """Route one candidate without using answers or benchmark hit status."""

    if comparison_status == "retained_from_selected_dna":
        return (
            "retain_baseline",
            ["BASELINE_001"],
            ["Candidate already exists in the frozen selected-DNA baseline."],
        )

    if str(candidate.get("union_status") or "") != "selected_any_test":
        return (
            "hold_not_forwarded",
            ["QC_001"],
            ["Candidate is filtered-only and has no selected observation in an evaluable test."],
        )

    family = str(taxonomy.get("primary_rule_family") or "unmapped_or_uncertain")
    biological_class = str(taxonomy.get("biological_class") or candidate.get("category") or "unknown")
    axes = _evidence_axes(candidate, policy)
    molecule = _molecule_support(candidate)
    rules: list[str] = []
    reasons: list[str] = []

    if axes.cross_molecule_support:
        rules.append("AXIS_DNA_RNA")
        reasons.append("Detected by both DNA and RNA tests.")
    if axes.repeated_detection:
        rules.append("AXIS_REPEAT")
        reasons.append("Detected in at least two distinct tests.")
    if axes.priority_rank_support:
        rules.append("AXIS_RANK")
        reasons.append("Best RK/NTC priority rank is within the top two of its test.")

    priority_families = set(policy["priority_families"])
    conditional_families = set(policy["conditional_families"])
    context_families = set(policy["context_families"])
    background_families = set(policy["background_families"])
    uncertain_families = set(policy["uncertain_families"])

    rna_only_nonviral = molecule == "RNA-only" and biological_class not in {"viral", "virus"}
    trace_like = (
        int(candidate.get("detected_test_count") or 0) <= 1
        and (candidate.get("max_reads_single_test") or 0) <= 3
        and not axes.priority_rank_support
    )

    if family in priority_families:
        rules.append("FAMILY_PRIORITY")
        if axes.independent_dimension_count > 0 and not rna_only_nonviral:
            return "forward_priority_review", rules, reasons + ["Priority clinical family with analytical support."]
        if axes.independent_dimension_count > 1:
            return "forward_priority_review", rules, reasons + ["RNA-only nonviral signal has more than one independent analytical support axis."]
        if axes.independent_dimension_count > 0:
            return "forward_context_review", rules + ["RNA_ONLY_CAUTION"], reasons + ["RNA-only bacterial/fungal evidence is retained as context, not promoted directly."]
        return "audit_only", rules + ["WEAK_PRIORITY_AUDIT"], reasons + ["Clinically important family, but current analytical support is weak."]

    if family in conditional_families:
        rules.append("FAMILY_CONDITIONAL")
        if axes.independent_dimension_count >= 2 and not trace_like:
            return "forward_priority_review", rules, reasons + ["Opportunistic mold/fungus has at least two analytical support axes."]
        if axes.independent_dimension_count >= 1 and not trace_like:
            return "forward_context_review", rules, reasons + ["Opportunistic mold/fungus has limited analytical support and requires host/specimen review."]
        return "audit_only", rules + ["MOLD_WEAK_AUDIT"], reasons + ["Weak opportunistic fungal signal retained for audit only."]

    if family in context_families:
        rules.append("FAMILY_CONTEXT")
        if axes.independent_dimension_count >= 1 and not trace_like:
            return "forward_context_review", rules, reasons + ["Colonization/reactivation-prone family is forwarded only as context despite analytical support."]
        return "audit_only", rules + ["CONTEXT_WEAK_AUDIT"], reasons + ["Colonization/reactivation-prone family lacks sufficient analytical support for automatic forwarding."]

    if family in background_families:
        rules.append("FAMILY_BACKGROUND")
        if axes.reproducibility_support and axes.priority_rank_support and not trace_like:
            return "audit_only", rules + ["BACKGROUND_STRONG_AUDIT"], reasons + ["Background/nonpulmonary-prone family has reproducible signal but still requires manual audit."]
        return "hold_not_forwarded", rules + ["BACKGROUND_HOLD"], reasons + ["Background/nonpulmonary-prone family is not forwarded automatically."]

    if family in uncertain_families or taxonomy.get("needs_literature_review"):
        rules.append("FAMILY_UNCERTAIN")
        if axes.reproducibility_support and not trace_like:
            return "audit_only", rules + ["UNCERTAIN_STRONG_AUDIT"], reasons + ["Unmapped organism has reproducible analytical support and needs taxonomy/literature audit."]
        return "hold_not_forwarded", rules + ["UNCERTAIN_HOLD"], reasons + ["Unmapped or uncertain organism is weakly supported and is not forwarded."]

    return "audit_only", ["FALLBACK_AUDIT"], reasons + ["No explicit family route matched; retained for audit without automatic forwarding."]


def aggregate_candidates(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["patient_id"], row["organism_key"])].append(row)

    result: list[dict[str, Any]] = []
    for (patient_id, organism_key), items in sorted(grouped.items()):
        first = items[0]
        detected_count = sum(_as_int(item.get("detected_test_count"), 0) or 0 for item in items)
        selected_count = sum(_as_int(item.get("selected_test_count"), 0) or 0 for item in items)
        dna_count = sum(_as_int(item.get("dna_detected_test_count"), 0) or 0 for item in items)
        rna_count = sum(_as_int(item.get("rna_detected_test_count"), 0) or 0 for item in items)
        result.append(
            {
                "patient_id": patient_id,
                "case_review_ids": "|".join(sorted({item["case_review_id"] for item in items})),
                "specimen_codes": "|".join(sorted({item["specimen_code"] for item in items})),
                "specimen_sites": "|".join(sorted({item["specimen_site"] for item in items})),
                "collected_times": "|".join(sorted({item["collected_time"] for item in items})),
                "category": first["category"],
                "organism_name": first["organism_name"],
                "organism_key": organism_key,
                "union_status": "selected_any_test" if any(item["union_status"] == "selected_any_test" for item in items) else "filtered_only",
                "detected_test_count": detected_count,
                "selected_test_count": selected_count,
                "dna_detected_test_count": dna_count,
                "rna_detected_test_count": rna_count,
                "cross_molecule_support": bool(dna_count and rna_count),
                "repeated_detection": detected_count >= 2,
                "distinct_protocol_conditions": sum(_as_int(item.get("distinct_protocol_conditions"), 0) or 0 for item in items),
                "max_reads_single_test": _max_present(_as_int(item.get("max_reads_single_test")) for item in items),
                "max_rpm_total_single_test": _max_present(_as_float(item.get("max_rpm_total_single_test")) for item in items),
                "max_rpm_nonhost_estimated_single_test": _max_present(_as_float(item.get("max_rpm_nonhost_estimated_single_test")) for item in items),
                "best_analytical_rank_priority": _min_present(_as_int(item.get("best_analytical_rank_priority")) for item in items),
                "best_rank_in_retained_universe": _min_present(_as_int(item.get("best_rank_in_retained_universe")) for item in items),
                "source_case_row_count": len(items),
                "reads_sum_across_tests": None,
                "aggregation_warning": "Raw reads are not summed or directly compared across tests.",
            }
        )
    return result


def build_shadow(inventory_dir: Path, output_dir: Path, policy_path: Path) -> dict[str, Any]:
    candidate_path = inventory_dir / "candidate_union.csv"
    comparison_path = inventory_dir / "selected_dna_comparison.csv"
    if not candidate_path.exists() or not comparison_path.exists():
        raise FileNotFoundError("inventory_dir must contain candidate_union.csv and selected_dna_comparison.csv")

    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    candidates = aggregate_candidates(_read_csv(candidate_path))
    comparison = _read_csv(comparison_path)
    comparison_lookup = {
        (row["patient_id"], row["organism_key"]): row["comparison_status"]
        for row in comparison
    }

    output_rows: list[dict[str, Any]] = []
    for candidate in candidates:
        status = comparison_lookup.get(
            (candidate["patient_id"], candidate["organism_key"]),
            "new_from_multi_assay_union",
        )
        taxonomy = classify_organism(candidate["organism_name"], biological_class=candidate["category"])
        action, rule_ids, reasons = route_candidate(candidate, taxonomy, status, policy)
        axes = _evidence_axes(candidate, policy)
        output_rows.append(
            {
                **candidate,
                "comparison_status": status,
                "molecule_support": _molecule_support(candidate),
                "primary_rule_family": taxonomy["primary_rule_family"],
                "taxonomy_mapping_status": taxonomy["mapping_status"],
                "taxonomy_confidence": taxonomy["classification_confidence"],
                "taxonomy_needs_literature_review": taxonomy["needs_literature_review"],
                "clinical_traits": "|".join(taxonomy.get("clinical_traits", [])),
                "reproducibility_support": axes.reproducibility_support,
                "independent_dimension_count": axes.independent_dimension_count,
                "screening_action": action,
                "rule_ids": "|".join(rule_ids),
                "screening_reasons": " | ".join(reasons),
                "answer_blind": True,
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    new_rows = [row for row in output_rows if row["comparison_status"] == "new_from_multi_assay_union"]
    forwarded = [
        row
        for row in output_rows
        if row["screening_action"] in {"retain_baseline", "forward_priority_review", "forward_context_review"}
    ]
    fieldnames = list(output_rows[0].keys()) if output_rows else []
    _write_csv(output_dir / "screening_shadow.csv", output_rows, fieldnames)
    _write_csv(output_dir / "new_candidate_analysis.csv", new_rows, fieldnames)
    _write_csv(output_dir / "shadow_candidate_pool.csv", forwarded, fieldnames)

    route_counts = Counter(row["screening_action"] for row in output_rows)
    new_route_counts = Counter(row["screening_action"] for row in new_rows)
    new_family_counts = Counter(row["primary_rule_family"] for row in new_rows)
    new_molecule_counts = Counter(row["molecule_support"] for row in new_rows)
    new_category_counts = Counter(row["category"] for row in new_rows)
    new_mapping_counts = Counter(row["taxonomy_mapping_status"] for row in new_rows)
    new_axis_counts = Counter(str(row["independent_dimension_count"]) for row in new_rows)
    rna_only_rows = [row for row in new_rows if row["molecule_support"] == "RNA-only"]
    rna_only_route_counts = Counter(row["screening_action"] for row in rna_only_rows)
    rna_only_category_counts = Counter(row["category"] for row in rna_only_rows)
    family_route: dict[str, Counter[str]] = defaultdict(Counter)
    for row in new_rows:
        family_route[row["primary_rule_family"]][row["screening_action"]] += 1

    summary = {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "policy_id": policy["policy_id"],
        "answer_blind": True,
        "answer_keys_used": [],
        "inventory_dir": str(inventory_dir.resolve()),
        "policy_path": str(policy_path.resolve()),
        "source_hashes": {
            "candidate_union.csv": _sha256(candidate_path),
            "selected_dna_comparison.csv": _sha256(comparison_path),
            "policy": _sha256(policy_path),
        },
        "counts": {
            "all_distinct_patient_organism": len(output_rows),
            "retained_baseline": sum(row["comparison_status"] == "retained_from_selected_dna" for row in output_rows),
            "new_from_multi_assay_union": len(new_rows),
            "shadow_candidate_pool": len(forwarded),
        },
        "routes_all": dict(sorted(route_counts.items())),
        "routes_new": dict(sorted(new_route_counts.items())),
        "new_by_family": dict(sorted(new_family_counts.items())),
        "new_by_molecule_support": dict(sorted(new_molecule_counts.items())),
        "new_by_category": dict(sorted(new_category_counts.items())),
        "new_by_taxonomy_mapping": dict(sorted(new_mapping_counts.items())),
        "new_by_independent_dimension_count": dict(sorted(new_axis_counts.items())),
        "rna_only_routes": dict(sorted(rna_only_route_counts.items())),
        "rna_only_categories": dict(sorted(rna_only_category_counts.items())),
        "new_family_route_matrix": {
            family: dict(sorted(counts.items())) for family, counts in sorted(family_route.items())
        },
        "new_priority_candidates": [
            {
                "patient_id": row["patient_id"],
                "organism_name": row["organism_name"],
                "category": row["category"],
                "molecule_support": row["molecule_support"],
                "detected_test_count": row["detected_test_count"],
                "best_analytical_rank_priority": row["best_analytical_rank_priority"],
                "primary_rule_family": row["primary_rule_family"],
            }
            for row in new_rows
            if row["screening_action"] == "forward_priority_review"
        ],
        "interpretation_constraints": policy["safety_constraints"],
    }
    taxonomy_gaps = [
        row
        for row in new_rows
        if row["primary_rule_family"] == "unmapped_or_uncertain"
        and bool(row["reproducibility_support"])
    ]
    taxonomy_gaps.sort(
        key=lambda row: (
            -int(row["independent_dimension_count"] or 0),
            -int(row["max_reads_single_test"] or 0),
            row["patient_id"],
            row["organism_name"],
        )
    )
    _write_csv(output_dir / "taxonomy_gap_candidates.csv", taxonomy_gaps, fieldnames)
    summary["counts"]["taxonomy_gaps_with_reproducibility"] = len(taxonomy_gaps)
    (output_dir / "screening_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (output_dir / "SCREENING_REPORT.md").write_text(_render_report(summary), encoding="utf-8")
    return summary


def _render_report(summary: dict[str, Any]) -> str:
    counts = summary["counts"]
    lines = [
        "# Multi-assay mNGS test-aware screening shadow",
        "",
        f"- Policy: `{summary['policy_id']}`",
        "- Mode: answer-blind shadow",
        f"- Distinct patient-organism candidates: {counts['all_distinct_patient_organism']}",
        f"- Frozen selected-DNA baseline retained: {counts['retained_baseline']}",
        f"- New multi-assay candidates analyzed: {counts['new_from_multi_assay_union']}",
        f"- Shadow candidate pool (baseline + forwarded new): {counts['shadow_candidate_pool']}",
        "",
        "## Routing of the new candidates",
        "",
        "| Action | Count |",
        "|---|---:|",
    ]
    for action, count in summary["routes_new"].items():
        lines.append(f"| `{action}` | {count} |")
    lines.extend(["", "## New candidates by molecule support", "", "| Support | Count |", "|---|---:|"])
    for support, count in summary["new_by_molecule_support"].items():
        lines.append(f"| {support} | {count} |")
    lines.extend(["", "## New candidates by taxonomy family", "", "| Family | Count |", "|---|---:|"])
    for family, count in summary["new_by_family"].items():
        lines.append(f"| `{family}` | {count} |")
    lines.extend(["", "## Taxonomy coverage", "", "| Mapping status | Count |", "|---|---:|"])
    for status, count in summary["new_by_taxonomy_mapping"].items():
        lines.append(f"| `{status}` | {count} |")
    lines.extend(
        [
            "",
            f"Unmapped candidates with cross-test reproducibility: {counts['taxonomy_gaps_with_reproducibility']}. These are exported to `taxonomy_gap_candidates.csv`; they are not treated as rare pathogens automatically.",
            "",
            "## Family-by-route matrix",
            "",
            "| Family | Priority | Context | Audit | Hold |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for family, routes in summary["new_family_route_matrix"].items():
        lines.append(
            f"| `{family}` | {routes.get('forward_priority_review', 0)} | "
            f"{routes.get('forward_context_review', 0)} | {routes.get('audit_only', 0)} | "
            f"{routes.get('hold_not_forwarded', 0)} |"
        )
    lines.extend(
        [
            "",
            "## Priority-review candidates",
            "",
            "| Patient | Organism | Class | DNA/RNA support | Tests | Best within-test priority rank | Family |",
            "|---|---|---|---|---:|---:|---|",
        ]
    )
    for row in summary["new_priority_candidates"]:
        lines.append(
            f"| P{row['patient_id']} | {row['organism_name']} | {row['category']} | "
            f"{row['molecule_support']} | {row['detected_test_count']} | "
            f"{row['best_analytical_rank_priority'] or ''} | `{row['primary_rule_family']}` |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This shadow preserves every baseline candidate and does not assign clinical pathogenicity from mNGS alone. New candidates are routed by within-test rank, repeat detection, DNA/RNA concordance, taxonomy family, and QC state. Raw reads are retained for provenance but are not summed or compared directly across assays.",
            "",
            "`forward_priority_review` and `forward_context_review` mean the candidate should be visible to the downstream reviewer. They do not mean that the organism is a final pathogen or should enter Picked.",
            "",
            "## Safety constraints",
            "",
        ]
    )
    for constraint in summary["interpretation_constraints"]:
        lines.append(f"- {constraint}")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build answer-blind test-aware mNGS screening shadow.")
    parser.add_argument("inventory_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    summary = build_shadow(args.inventory_dir, args.output_dir, args.policy)
    print(json.dumps(summary["counts"], ensure_ascii=False, indent=2))
    print(json.dumps(summary["routes_new"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

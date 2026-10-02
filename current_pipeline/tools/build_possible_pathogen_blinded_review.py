"""Build answer-blind OBER review inputs for every added possible pathogen.

The production queue excludes benchmark labels, post-hoc outcomes, patient
numbers, and answer files. A private case map restores provenance only after
adjudication has been frozen.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from tools import build_rag_prototype_cases as rag_prototype
from tools.build_index_event_clinical_timeline_shadow import sha256_file
from tools.build_rag_v2_inputs import FORBIDDEN_PRODUCTION_KEYS
from tools.pathogen_normalization import canonical_key


SCHEMA_VERSION = "possible_pathogen_blinded_review_queue.v1"
PRODUCTION_SCHEMA_VERSION = "rag_production_case_v2.0"
POSSIBLE_ROLES = {
    "analytical_or_direct_possible_pathogen",
    "reproducible_low_specificity_possible_pathogen",
    "possible_concurrent_systemic_pathogen",
}


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")


def forbidden_paths(value: Any, prefix: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            path = f"{prefix}.{key}"
            if str(key).casefold() in FORBIDDEN_PRODUCTION_KEYS:
                found.append(path)
            found.extend(forbidden_paths(child, path))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.extend(forbidden_paths(child, f"{prefix}[{index}]"))
    return found


def _compact_rows(rows: Any, keys: tuple[str, ...], limit: int = 40) -> list[dict[str, Any]]:
    output = []
    for row in rows or []:
        if not isinstance(row, dict):
            continue
        compact = {key: row.get(key) for key in keys if row.get(key) not in (None, "", [], {})}
        if compact:
            output.append(compact)
        if len(output) >= limit:
            break
    return output


def _hospital_observations(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    if candidate.get("candidate_source") == "hospital_only":
        detail = candidate.get("hospital_evidence_detail") or {}
    else:
        detail = (candidate.get("hospital_profile") or {}).get("hospital_evidence_detail") or {}
    rows = detail.get("source_observations") or []
    return _compact_rows(
        rows,
        (
            "observation_id", "test_type", "panel_or_assay", "specimen_type",
            "specimen_category", "collected_time", "reported_time", "target",
            "organism_name", "raw_result", "detection_status", "numeric_value",
            "unit", "reference_or_cutoff", "semiquant_bin", "quantitation_status",
            "source_text", "contaminant_flag",
        ),
    )


def _host_context(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    return _compact_rows(
        candidate.get("host_evidence"),
        (
            "evidence_id", "evidence_role", "label", "assertion", "source_text",
            "observed_at", "reported_at", "episode_evidence_status",
            "relative_to_treatment",
        ),
        limit=30,
    )


def _image_context(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    return _compact_rows(
        candidate.get("image_evidence"),
        (
            "evidence_id", "label", "assertion", "source_text", "observed_at",
            "reported_at", "episode_evidence_status", "microbiologic_confirmation",
        ),
        limit=20,
    )


def _phenotype_context(candidate: dict[str, Any]) -> dict[str, Any]:
    source = candidate.get("phenotype_evidence") or {}
    return {
        "production_link_verified": bool(source.get("production_link_verified")),
        "shadow_reasoning_allowed": bool(source.get("shadow_reasoning_allowed")),
        "linkage_status": source.get("linkage_status"),
        "routed_context": _compact_rows(
            source.get("routed_phenotype_context"),
            (
                "source_id", "phenotype", "status", "confidence", "time_relation",
                "reason", "qc_status",
            ),
            limit=30,
        ),
        "constraints": [
            "Phenotype is host/syndrome context, not direct organism confirmation.",
            "NO or missing phenotype data is not proven absence.",
        ],
    }


def _event_context(candidate: dict[str, Any]) -> dict[str, Any]:
    source = candidate.get("source_candidate") or {}
    return {
        "specimen_context": source.get("specimen_context"),
        "collected_time": source.get("collected_time"),
        "case_review_id_hash": (
            hashlib.sha256(str(source.get("case_review_id")).encode("utf-8")).hexdigest()[:12]
            if source.get("case_review_id") else None
        ),
        "linked_events": _compact_rows(
            candidate.get("linked_events"),
            ("event_id", "specimen_site", "collected_time", "event_relation"),
            limit=10,
        ),
    }


def _analytical_evidence(candidate: dict[str, Any]) -> dict[str, Any]:
    profile = candidate.get("analytical_profile") or {}
    return {
        "best_rank_in_retained_universe": profile.get("best_rank_in_retained_universe"),
        "rank_band": profile.get("rank_band"),
        "selected_positive_test_count": profile.get("selected_positive_test_count"),
        "dna_selected": profile.get("dna_selected"),
        "rna_selected": profile.get("rna_selected"),
        "cross_molecule_selected": profile.get("cross_molecule_selected"),
        "technical_repeat": profile.get("technical_repeat"),
        "reproducibility_axis": profile.get("reproducibility_axis"),
        "reproducibility_note": profile.get("reproducibility_note"),
        "normalization_available_test_count": profile.get("normalization_available_test_count"),
        "fully_evaluable_test_count": profile.get("fully_evaluable_test_count"),
        "per_test_signals": _compact_rows(
            profile.get("per_test_signals"),
            (
                "seq_id", "nucleic_type", "reads", "rank_in_retained_universe",
                "analytical_rank_priority", "rpm_total", "qc_status",
                "normalization_available",
            ),
        ),
        "reads_sum_across_tests": None,
    }


def _direct_timing(candidate: dict[str, Any]) -> dict[str, Any]:
    gate = candidate.get("possible_pathogen_gate") or {}
    axes = gate.get("axes") or {}
    profile = ((candidate.get("promotion_gate") or {}).get("axes") or {}).get(
        "direct_evidence_timing_profile"
    ) or {}
    return {
        "direct_hospital_level": axes.get("direct_hospital_level"),
        "event_window_hours": profile.get("event_window_hours"),
        "event_aligned_positive_count": profile.get("event_aligned_positive_count"),
        "rows": _compact_rows(
            profile.get("rows"),
            (
                "observation_id", "positive", "timing", "time_basis",
                "collected_time", "reported_time", "specimen_type", "target",
                "raw_result",
            ),
        ),
    }


def make_case(
    blind_id: str,
    candidate: dict[str, Any],
    strict_competitors: list[dict[str, Any]],
    source_sha256: str,
) -> dict[str, Any]:
    name = str(candidate.get("organism_name") or "")
    taxonomy = candidate.get("taxonomy_profile") or {}
    biological_class = str(taxonomy.get("biological_class") or "Other")
    category = rag_prototype.pathogen_category(name, biological_class)
    possible_gate = candidate.get("possible_pathogen_gate") or {}
    case = {
        "schema_version": PRODUCTION_SCHEMA_VERSION,
        "case_id": blind_id,
        "dataset_label": "KH_multi_assay_possible_pathogen_blinded_review_v1",
        "patient_id": blind_id,
        "organism": {
            "display_name": name,
            "canonical_key": canonical_key(name),
            "classification": biological_class,
            "pathogen_category": category,
            "taxonomy_unknown": taxonomy.get("mapping_status") == "unmapped",
            "taxonomy_profile": taxonomy,
        },
        "current_model_state": {
            "source_tier": candidate.get("clinical_decision"),
            "formal_picked": False,
            "possible_reporting_role": candidate.get("possible_reporting_role"),
            "possible_route": possible_gate.get("route"),
            "possible_cautions": possible_gate.get("cautions") or [],
            "strict_competing_pathogens": [
                {
                    "organism_name": item.get("organism_name"),
                    "taxonomy_family": (item.get("taxonomy_profile") or {}).get(
                        "primary_rule_family"
                    ),
                    "selection_role": item.get("selection_role"),
                }
                for item in strict_competitors
                if canonical_key(item.get("organism_name")) != canonical_key(name)
            ],
        },
        "local_evidence": {
            "index_event": _event_context(candidate),
            "mngs_per_test_evidence": _analytical_evidence(candidate),
            "same_organism_hospital_observations": _hospital_observations(candidate),
            "direct_evidence_timing": _direct_timing(candidate),
            "host_context": _host_context(candidate),
            "image_context": _image_context(candidate),
            "phenotype_context": _phenotype_context(candidate),
            "rule_trace": {
                "clinical_reasons": candidate.get("reasons") or [],
                "clinical_rule_ids": candidate.get("rule_ids") or [],
                "possible_gate": possible_gate,
            },
            "interpretation_constraints": [
                "Do not sum reads across tests.",
                "DNA/RNA concordance and technical repeats are one analytical axis, not independent clinical confirmation.",
                "A possible-pathogen role is not a proven causal or primary pathogen claim.",
                "Missing timing or normalization is unknown, not negative evidence.",
                "Evaluate the candidate against competing pathogens and specimen-specific colonization risk.",
            ],
        },
        "rag_task": {
            "knowledge_card_key": f"adult_pulmonary:{canonical_key(name)}",
            "retrieval_priority": "high",
            "fixed_questions": rag_prototype.rag_questions(),
            "suggested_search_queries": rag_prototype.search_queries(name, category),
            "required_output_schema": "rag_adjudication_v2.0",
        },
        "provenance": {
            "possible_pathogen_policy_id": "test_aware_possible_pathogen_v1_answer_blind_generic",
            "source_candidate_sha256": source_sha256,
            "private_identity_map_required": True,
            "answer_source_read": False,
        },
    }
    leaked = forbidden_paths(case)
    if leaked:
        raise ValueError(f"Forbidden benchmark fields in {blind_id}: {', '.join(leaked)}")
    return case


def make_knowledge_request(case: dict[str, Any]) -> dict[str, Any]:
    task = case["rag_task"]
    return {
        "schema_version": "rag_knowledge_request_v2.0",
        "knowledge_card_key": task["knowledge_card_key"],
        "organism": case["organism"],
        "retrieval_priority": task["retrieval_priority"],
        "fixed_questions": task["fixed_questions"],
        "suggested_search_queries": task["suggested_search_queries"],
        "answer_source_read": False,
    }


def run(input_root: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    paths = sorted((input_root / "patient_outputs").glob(
        "NGS_patient_*_test_aware_possible_pathogen_shadow.json"
    ))
    if not paths:
        raise ValueError("No possible-pathogen patient outputs found")

    source_rows: list[tuple[Path, dict[str, Any], list[dict[str, Any]]]] = []
    for path in paths:
        payload = read_json(path)
        strict = payload.get("strict_picked") or []
        for candidate in payload.get("possible_pathogens") or []:
            if candidate.get("possible_reporting_role") not in POSSIBLE_ROLES:
                continue
            source_rows.append((path, candidate, strict))
    source_rows.sort(key=lambda row: (
        canonical_key(row[1].get("organism_name")),
        str(row[1].get("patient_id") or ""),
    ))
    if not source_rows:
        raise ValueError("No added possible-pathogen candidates found")

    cases = []
    private_rows = []
    for index, (path, candidate, strict) in enumerate(source_rows, 1):
        blind_id = f"BLIND-{index:03d}"
        candidate_fingerprint = hashlib.sha256(
            json.dumps(candidate, ensure_ascii=False, sort_keys=True).encode("utf-8")
        ).hexdigest()
        cases.append(make_case(blind_id, candidate, strict, candidate_fingerprint))
        private_rows.append({
            "blind_case_id": blind_id,
            "patient_id": candidate.get("patient_id"),
            "organism_name": candidate.get("organism_name"),
            "possible_reporting_role": candidate.get("possible_reporting_role"),
            "source_file": str(path.resolve()),
            "source_file_sha256": sha256_file(path),
            "source_candidate_sha256": candidate_fingerprint,
        })

    requests_by_key = {}
    for case in cases:
        request = make_knowledge_request(case)
        requests_by_key.setdefault(request["knowledge_card_key"], request)
    requests = [requests_by_key[key] for key in sorted(requests_by_key)]

    output_dir.mkdir(parents=True, exist_ok=True)
    write_jsonl(output_dir / "production_queue.jsonl", cases)
    write_jsonl(output_dir / "knowledge_requests.jsonl", requests)
    private_dir = output_dir / "private_do_not_send_to_model"
    private_dir.mkdir()
    write_json(private_dir / "case_identity_map.json", {
        "schema_version": "possible_pathogen_private_case_map.v1",
        "private": True,
        "rows": private_rows,
    })
    with (private_dir / "case_identity_map.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(private_rows[0]))
        writer.writeheader()
        writer.writerows(private_rows)

    role_counts = Counter(
        case["current_model_state"]["possible_reporting_role"] for case in cases
    )
    family_counts = Counter(
        case["organism"]["taxonomy_profile"].get("primary_rule_family") for case in cases
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "answer_blind": True,
        "patient_numbers_exposed_to_model": False,
        "input_root": str(input_root.resolve()),
        "input_summary_sha256": (
            sha256_file(input_root / "summary.json")
            if (input_root / "summary.json").is_file() else None
        ),
        "case_count": len(cases),
        "knowledge_request_count": len(requests),
        "forbidden_benchmark_field_count": sum(len(forbidden_paths(case)) for case in cases),
        "role_counts": dict(sorted(role_counts.items())),
        "family_counts": dict(sorted(family_counts.items(), key=lambda item: str(item[0]))),
        "private_case_map": str((private_dir / "case_identity_map.json").resolve()),
        "constraints": [
            "Do not provide the private case map or benchmark evaluator outputs to the adjudication model.",
            "Freeze adjudications before joining the private case map or benchmark answers.",
            "This cohort was used to design the possible-pathogen thresholds; it is an internal blinded review, not held-out validation.",
        ],
    }
    write_json(output_dir / "manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.input_root, args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

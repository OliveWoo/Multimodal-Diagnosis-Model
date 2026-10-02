"""Build answer-blind phenotype context for patient-organism shadow review.

All 22 phenotype statuses remain visible in a compact overview. Organism
taxonomy only controls which phenotype rows are expanded with source evidence.
This module never changes deterministic tiers or Picked.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

from tools import organism_taxonomy_classifier as taxonomy
from tools import pathogen_normalization as pathogen_names


REPO_ROOT = Path(__file__).resolve().parent.parent
POLICY_PATH = REPO_ROOT / "rules" / "phenotype_candidate_policy.json"
SCHEMA_VERSION = "candidate_phenotype_context_v1.0"
LINKAGE_MODES = {"disabled", "same_number_shadow", "verified"}
FORBIDDEN_KEYS = {"answer", "answers", "benchmark_answer", "answer_hit", "is_answer"}


def as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def find_forbidden_keys(value: Any, path: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, dict):
        for raw_key, child in value.items():
            key = str(raw_key)
            child_path = f"{path}.{key}"
            if key.lower() in FORBIDDEN_KEYS:
                found.append(child_path)
            found.extend(find_forbidden_keys(child, child_path))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.extend(find_forbidden_keys(child, f"{path}[{index}]"))
    return found


@lru_cache(maxsize=1)
def policy() -> dict[str, Any]:
    payload = read_json(POLICY_PATH)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected object in {POLICY_PATH}")
    forbidden = find_forbidden_keys(payload)
    if forbidden:
        raise ValueError(f"Answer-derived keys are forbidden in phenotype policy: {forbidden}")
    taxonomy_families = set(as_dict(taxonomy.rules_payload().get("families")))
    configured_families = set(as_dict(payload.get("families")))
    if configured_families != taxonomy_families:
        missing = sorted(taxonomy_families - configured_families)
        extra = sorted(configured_families - taxonomy_families)
        raise ValueError(f"Phenotype policy family mismatch: missing={missing}, extra={extra}")
    return payload


POLICY_VERSION = str(policy().get("policy_version") or "phenotype_candidate_shadow_v1")


def patient_number(value: Any) -> int | None:
    match = re.search(r"(?:patient[_\s-]*|^)(\d+)", str(value or ""), re.I)
    return int(match.group(1)) if match else None


def ledger_path(phenotype_root: Path, patient_id: Any) -> Path | None:
    number = patient_number(patient_id)
    if number is None:
        return None
    path = phenotype_root / f"patient_{number}_phenotype_decision_evidence_v1.json"
    return path if path.is_file() else None


def load_ledger(phenotype_root: Path, patient_id: Any) -> tuple[dict[str, Any] | None, Path | None]:
    path = ledger_path(phenotype_root, patient_id)
    if path is None:
        return None, None
    payload = read_json(path)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected object in {path}")
    return payload, path


def linkage(mode: str, patient_id: Any, ledger: dict[str, Any] | None) -> dict[str, Any]:
    if mode not in LINKAGE_MODES:
        raise ValueError(f"Unknown phenotype linkage mode: {mode}")
    requested = patient_number(patient_id)
    source = patient_number(ledger.get("patient_id")) if ledger else None
    same_number = bool(requested is not None and requested == source)
    if ledger is None:
        status = "phenotype_file_missing"
    elif mode == "disabled":
        status = "available_but_disabled"
    elif mode == "verified":
        status = "verified_patient_link"
    else:
        status = "provisional_same_number_link"
    return {
        "mode": mode,
        "status": status,
        "pipeline_patient_number": requested,
        "phenotype_patient_number": source,
        "same_number": same_number,
        "shadow_reasoning_allowed": bool(
            ledger is not None and same_number and mode in {"same_number_shadow", "verified"}
        ),
        "production_link_verified": bool(ledger is not None and same_number and mode == "verified"),
        "note": (
            "same_number_shadow is a provisional research link and must not be represented as verified provenance"
            if mode == "same_number_shadow"
            else None
        ),
    }


def record_rows(ledger: dict[str, Any], sheet: str) -> list[dict[str, Any]]:
    return [item for item in as_list(as_dict(ledger.get("records")).get(sheet)) if isinstance(item, dict)]


def compact_phenotype(item: dict[str, Any]) -> dict[str, Any]:
    row = as_dict(item.get("row"))
    return {
        "source_id": item.get("id"),
        "phenotype": row.get("phenotype"),
        "status": row.get("status"),
        "confidence": row.get("confidence"),
        "temporal_relation": row.get("temporal_relation"),
        "phenotype_qc_status": row.get("phenotype_qc_status"),
        "model_role": phenotype_model_role(row),
    }


def phenotype_model_role(row: dict[str, Any]) -> str:
    status = str(row.get("status") or "").upper()
    confidence = str(row.get("confidence") or "").upper()
    temporal = str(row.get("temporal_relation") or "").upper()
    qc = str(row.get("phenotype_qc_status") or "").upper()
    if status != "YES":
        return "not_established"
    if qc != "PASS" or confidence not in {"HIGH", "MEDIUM"}:
        return "manual_review"
    if temporal in {"PRE_PNEUMONIA", "AT_ONSET"}:
        return "etiologic_prior"
    if temporal == "AFTER_ONSET":
        return "outcome_or_context"
    return "context_only"


def overview(ledger: dict[str, Any]) -> list[dict[str, Any]]:
    return [compact_phenotype(item) for item in record_rows(ledger, "Phenotypes_Long")]


def family_policy(profile: dict[str, Any]) -> dict[str, str]:
    families = [profile.get("primary_rule_family"), *as_list(profile.get("secondary_rule_families"))]
    output: dict[str, str] = {}
    configured = as_dict(policy().get("families"))
    for family in families:
        family_rules = as_dict(configured.get(str(family)))
        for phenotype, relation in as_dict(family_rules.get("phenotypes")).items():
            output.setdefault(str(phenotype), str(relation))
    return output


def interpretation(role: str, relation: str) -> str:
    if role == "not_established":
        return "not_established_no_negative_effect"
    if role == "manual_review":
        return "manual_review_only"
    if role == "outcome_or_context":
        return "outcome_context_not_etiologic_prior"
    if role == "context_only":
        return "timing_unclear_context_only"
    return relation


def evidence_by_phenotype(ledger: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    output: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for item in record_rows(ledger, "Evidence"):
        name = str(as_dict(item.get("row")).get("phenotype") or "").strip()
        if name:
            output[name].append(item)
    return output


def review_by_phenotype(ledger: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    output: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for item in record_rows(ledger, "Review_Only"):
        name = str(as_dict(item.get("row")).get("phenotype") or "__SYSTEM__").strip()
        output[name].append(item)
    return output


def compact_source_record(item: dict[str, Any]) -> dict[str, Any]:
    row = as_dict(item.get("row"))
    return {
        "source_id": item.get("id"),
        "candidate_status": row.get("candidate_status"),
        "confidence": row.get("confidence"),
        "temporal_relation": row.get("temporal_relation"),
        "reason_codes": row.get("reason_codes"),
        "reason": row.get("reason"),
        "evidence": row.get("evidence"),
        "evidence_page": row.get("evidence_page"),
        "source_type": row.get("source_type"),
        "source_date": row.get("source_date"),
        "chunk_id": row.get("chunk_id"),
    }


def compact_review_record(item: dict[str, Any]) -> dict[str, Any]:
    row = as_dict(item.get("row"))
    return {
        "source_id": item.get("id"),
        "status": row.get("status"),
        "review_reason": row.get("review_reason"),
        "evidence": row.get("evidence"),
        "evidence_page": row.get("evidence_page"),
        "source_type": row.get("source_type"),
        "source_date": row.get("source_date"),
    }


def relevant_context(
    ledger: dict[str, Any],
    *,
    organism_name: str,
    classification: Any = None,
    taxonomy_profile: dict[str, Any] | None = None,
) -> dict[str, Any]:
    profile = taxonomy_profile or taxonomy.classify_organism(
        organism_name, biological_class=classification
    )
    routes = family_policy(profile)
    source_evidence = evidence_by_phenotype(ledger)
    review_evidence = review_by_phenotype(ledger)
    expanded: list[dict[str, Any]] = []
    for item in record_rows(ledger, "Phenotypes_Long"):
        row = as_dict(item.get("row"))
        name = str(row.get("phenotype") or "")
        relation = routes.get(name)
        if relation is None:
            continue
        role = phenotype_model_role(row)
        expanded.append(
            {
                **compact_phenotype(item),
                "policy_relation": relation,
                "shadow_interpretation": interpretation(role, relation),
                "eligible_etiologic_prior": role == "etiologic_prior",
                "reason": row.get("reason"),
                "evidence": row.get("evidence"),
                "evidence_page": row.get("evidence_page"),
                "source_type": row.get("source_type"),
                "source_date": row.get("source_date"),
                "corroborating_evidence": [
                    compact_source_record(record) for record in source_evidence.get(name, [])
                ],
                "phenotype_review_records": [
                    compact_review_record(record) for record in review_evidence.get(name, [])
                ],
            }
        )
    counts = Counter(item["shadow_interpretation"] for item in expanded)
    return {
        "organism_name": organism_name,
        "canonical_key": pathogen_names.canonical_key(organism_name),
        "taxonomy_profile": profile,
        "expanded_phenotypes": expanded,
        "summary": {
            "configured_relevant_phenotype_count": len(routes),
            "expanded_phenotype_count": len(expanded),
            "eligible_etiologic_prior_count": sum(
                1 for item in expanded if item["eligible_etiologic_prior"]
            ),
            "interpretation_counts": dict(sorted(counts.items())),
            "system_review_record_count": len(review_evidence.get("__SYSTEM__", [])),
            "allowed_effect": "shadow_reasoning_only_no_deterministic_change",
        },
    }


def candidate_values(candidates: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        name = str(candidate.get("organism_name") or candidate.get("name") or "").strip()
        canonical = pathogen_names.canonical_key(name)
        if not name or not canonical or canonical in seen:
            continue
        seen.add(canonical)
        output.append(candidate)
    return output


def unavailable_candidate(candidate: dict[str, Any], reason: str) -> dict[str, Any]:
    name = str(candidate.get("organism_name") or candidate.get("name") or "").strip()
    profile = taxonomy.classify_organism(name, biological_class=candidate.get("classification"))
    return {
        "organism_name": name,
        "canonical_key": pathogen_names.canonical_key(name),
        "taxonomy_profile": profile,
        "expanded_phenotypes": [],
        "summary": {
            "configured_relevant_phenotype_count": len(family_policy(profile)),
            "expanded_phenotype_count": 0,
            "eligible_etiologic_prior_count": 0,
            "interpretation_counts": {},
            "system_review_record_count": 0,
            "allowed_effect": "none",
            "unavailable_reason": reason,
        },
    }


def build_patient_shadow(
    phenotype_root: Path,
    *,
    patient_id: Any,
    candidates: Iterable[dict[str, Any]],
    linkage_mode: str,
) -> dict[str, Any]:
    ledger, path = load_ledger(phenotype_root, patient_id)
    link = linkage(linkage_mode, patient_id, ledger)
    candidate_list = candidate_values(candidates)
    if ledger is None:
        contexts = [unavailable_candidate(item, "phenotype_file_missing") for item in candidate_list]
        all_phenotypes: list[dict[str, Any]] = []
        source: dict[str, Any] = {"ledger_file": None}
        system_review_ids: list[str] = []
    else:
        contexts = [
            relevant_context(
                ledger,
                organism_name=str(item.get("organism_name") or item.get("name")),
                classification=item.get("classification"),
                taxonomy_profile=(
                    item.get("taxonomy_profile")
                    if isinstance(item.get("taxonomy_profile"), dict)
                    else None
                ),
            )
            for item in candidate_list
        ]
        all_phenotypes = overview(ledger)
        source = {
            "ledger_file": str(path),
            "ledger_schema_version": ledger.get("schema_version"),
            "source_workbook": as_dict(ledger.get("source")).get("workbook_file"),
            "source_workbook_sha256": as_dict(ledger.get("source")).get("workbook_sha256"),
            "patient_qc_status": as_dict(ledger.get("source")).get("patient_qc_status"),
        }
        system_review_ids = [
            str(item.get("id"))
            for item in record_rows(ledger, "Review_Only")
            if str(as_dict(item.get("row")).get("phenotype") or "__SYSTEM__") == "__SYSTEM__"
        ]
    output = {
        "schema_version": SCHEMA_VERSION,
        "policy_version": POLICY_VERSION,
        "patient_id": str(patient_id),
        "source": source,
        "linkage": link,
        "all_phenotype_overview": all_phenotypes,
        "candidate_contexts": contexts,
        "qc_summary": {
            "patient_qc_status": source.get("patient_qc_status"),
            "system_review_record_count": len(system_review_ids),
            "system_review_source_ids": system_review_ids,
        },
        "constraints": as_list(policy().get("global_constraints")),
    }
    forbidden = find_forbidden_keys(output)
    if forbidden:
        raise ValueError(f"Forbidden answer-derived keys in phenotype shadow: {forbidden}")
    return output


def candidate_context_from_patient_shadow(
    patient_shadow: dict[str, Any], organism_name: Any
) -> dict[str, Any]:
    canonical = pathogen_names.canonical_key(organism_name)
    selected = next(
        (
            item
            for item in as_list(patient_shadow.get("candidate_contexts"))
            if as_dict(item).get("canonical_key") == canonical
        ),
        None,
    )
    return {
        "schema_version": patient_shadow.get("schema_version"),
        "policy_version": patient_shadow.get("policy_version"),
        "patient_id": patient_shadow.get("patient_id"),
        "source": patient_shadow.get("source"),
        "linkage": patient_shadow.get("linkage"),
        "all_phenotype_overview": patient_shadow.get("all_phenotype_overview"),
        "candidate_relevant_context": selected,
        "qc_summary": patient_shadow.get("qc_summary"),
        "constraints": patient_shadow.get("constraints"),
    }


def queue_path(patient_dir: Path, suffix: str) -> Path:
    match = re.search(r"NGS_patient_(\d+)_json$", patient_dir.name)
    if not match:
        raise ValueError(f"Unexpected patient directory name: {patient_dir.name}")
    return patient_dir / "summary_outputs" / f"NGS_patient_{match.group(1)}_{suffix}.json"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--phenotype-root", type=Path, required=True)
    parser.add_argument("--queue-suffix", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--linkage-mode", choices=sorted(LINKAGE_MODES), default="disabled")
    args = parser.parse_args()

    patient_dirs = sorted(
        args.patient_root.glob("NGS_patient_*_json"),
        key=lambda path: patient_number(path.name) or 10**9,
    )
    records: list[dict[str, Any]] = []
    for patient_dir in patient_dirs:
        number = patient_number(patient_dir.name)
        qpath = queue_path(patient_dir, args.queue_suffix)
        if not qpath.is_file():
            records.append({"patient_id": number, "status": "queue_missing", "queue": str(qpath)})
            continue
        queue = read_json(qpath)
        candidates = as_list(as_dict(queue).get("review_queue"))
        shadow = build_patient_shadow(
            args.phenotype_root,
            patient_id=number,
            candidates=candidates,
            linkage_mode=args.linkage_mode,
        )
        output_path = args.output / f"patient_{number}_candidate_phenotype_shadow_v1.json"
        write_json(output_path, shadow)
        records.append(
            {
                "patient_id": number,
                "status": shadow["linkage"]["status"],
                "shadow_reasoning_allowed": shadow["linkage"]["shadow_reasoning_allowed"],
                "phenotype_count": len(shadow["all_phenotype_overview"]),
                "candidate_count": len(shadow["candidate_contexts"]),
                "output_file": output_path.name,
            }
        )
    audit = {
        "schema_version": SCHEMA_VERSION,
        "policy_version": POLICY_VERSION,
        "linkage_mode": args.linkage_mode,
        "patient_count": len(patient_dirs),
        "status_counts": dict(Counter(item["status"] for item in records)),
        "candidate_count": sum(int(item.get("candidate_count") or 0) for item in records),
        "patients": records,
    }
    write_json(args.output / "candidate_phenotype_shadow_audit.json", audit)
    print(json.dumps({key: value for key, value in audit.items() if key != "patients"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

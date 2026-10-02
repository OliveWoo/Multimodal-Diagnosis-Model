"""Build candidate-specific phenotype context for test-aware scorer inputs.

This adapter keeps the high-recall forwarding set intact, adds reusable
phenotype context, and records which multi-assay events each mNGS candidate
belongs to. Phenotype remains provisional and cannot change deterministic
decisions in this stage.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools.build_index_event_clinical_timeline_shadow import sha256_file
from tools.candidate_phenotype_context import build_patient_shadow, find_forbidden_keys
from tools.pathogen_normalization import canonical_key


SCHEMA_VERSION = "test_aware_candidate_phenotype_shadow.v1"


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return payload


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot identify patient from {path}")
    return int(match.group(1))


def event_refs(candidate: dict[str, Any]) -> list[str]:
    refs = {
        str(item.get("case_review_id") or "").strip()
        for item in candidate.get("case_decision_refs") or []
        if isinstance(item, dict) and str(item.get("case_review_id") or "").strip()
    }
    direct = str(candidate.get("case_review_id") or "").strip()
    if direct:
        refs.add(direct)
    return sorted(refs)


def enrich_contexts(
    shadow: dict[str, Any], candidates: list[dict[str, Any]]
) -> dict[str, Any]:
    by_key = {
        canonical_key(item.get("organism_name")): item
        for item in candidates
        if isinstance(item, dict)
    }
    for context in shadow.get("candidate_contexts") or []:
        source = by_key.get(str(context.get("canonical_key") or ""))
        if source is None:
            raise ValueError(f"No forwarded source for {context.get('organism_name')}")
        refs = event_refs(source)
        hospital_only = str(source.get("evidence_source") or "") == "hospital_only"
        context.update({
            "candidate_event_refs": refs,
            "event_scope": (
                "patient_episode_all_events"
                if hospital_only else
                "listed_multi_assay_cases"
            ),
            "candidate_source": "hospital_only" if hospital_only else "mngs_positive_reads",
            "incoming_decision": source.get("decision"),
            "incoming_level": source.get("integrated_level"),
            "incoming_rule_ids": source.get("rule_ids") or [],
            "phenotype_deterministic_effect": "none_provisional_link",
        })
    shadow.update({
        "schema_version": SCHEMA_VERSION,
        "source_candidate_count": len(candidates),
        "event_link_policy": (
            "mNGS candidates are limited to listed case_review_id values; "
            "hospital-only candidates remain patient-episode context"
        ),
        "decision_effect": "none_shadow_only",
    })
    forbidden = find_forbidden_keys(shadow)
    if forbidden:
        raise ValueError(f"Benchmark answer keys found: {forbidden[:5]}")
    return shadow


def run(
    scorer_root: Path,
    phenotype_root: Path,
    output_dir: Path,
    *,
    linkage_mode: str = "same_number_shadow",
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    patient_root = scorer_root / "patient_outputs"
    paths = sorted(patient_root.glob("NGS_patient_*_test_aware_deterministic_shadow.json"),
                   key=patient_number)
    if not paths:
        raise ValueError(f"No test-aware patient outputs found in {patient_root}")

    output_dir.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    for path in paths:
        patient = patient_number(path)
        source = read_json(path)
        candidates = source.get("clinical_scorer_forward") or []
        shadow = build_patient_shadow(
            phenotype_root,
            patient_id=patient,
            candidates=candidates,
            linkage_mode=linkage_mode,
        )
        shadow = enrich_contexts(shadow, candidates)
        shadow["source_test_aware_file"] = str(path.resolve())
        shadow["source_test_aware_sha256"] = sha256_file(path)
        target = output_dir / f"patient_{patient}_candidate_phenotype_shadow_v1.json"
        target.write_text(json.dumps(shadow, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        counts["patients"] += 1
        counts["candidates"] += len(shadow.get("candidate_contexts") or [])
        counts[f"phenotype_link:{shadow.get('linkage', {}).get('status')}"] += 1
        for item in shadow.get("candidate_contexts") or []:
            counts[f"candidate_source:{item.get('candidate_source')}"] += 1
            counts[f"incoming:{item.get('incoming_decision')}"] += 1

    summary = {
        "schema_version": SCHEMA_VERSION,
        "answer_blind": True,
        "decision_effect": "none_shadow_only",
        "counts": dict(sorted(counts.items())),
        "source_scorer_root": str(scorer_root.resolve()),
        "phenotype_root": str(phenotype_root.resolve()),
        "linkage_mode": linkage_mode,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scorer_root", type=Path)
    parser.add_argument("phenotype_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--linkage-mode",
        choices=("disabled", "same_number_shadow", "verified"),
        default="same_number_shadow",
    )
    args = parser.parse_args()
    print(json.dumps(run(
        args.scorer_root,
        args.phenotype_root,
        args.output_dir,
        linkage_mode=args.linkage_mode,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

"""Rehydrate source hospital observations without changing frozen clinical tiers.

This answer-blind adapter fills the source-observation rows referenced by an
existing clinical shadow.  It never reruns analytical, history, clinical, or
promotion decisions.  Semantic hospital levels are checked against the frozen
candidate before the richer provenance is attached.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools import deterministic_mngs_max_scorer as legacy
from tools.build_test_aware_deterministic_shadow_scorer import (
    enrich_hospital_evidence_detail,
    hospital_observation_index,
    hospital_profile,
)


SCHEMA_VERSION = "test_aware_clinical_evidence_rehydrated_shadow.v1"
SEMANTIC_PROFILE_KEYS = (
    "exact_or_alias_evidence",
    "matched_hospital_name",
    "best_direct_level",
    "direct_support_modules_level_1_2",
    "direct_support_modules_level_1_3",
    "module_level_summary",
)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot identify patient from {path}")
    return int(match.group(1))


def _semantic_profile(profile: dict[str, Any]) -> dict[str, Any]:
    return {key: profile.get(key) for key in SEMANTIC_PROFILE_KEYS}


def _resolution(detail: dict[str, Any] | None) -> dict[str, Any]:
    if not isinstance(detail, dict):
        return {"requested": 0, "resolved": 0, "missing": []}
    value = detail.get("source_observation_resolution") or {}
    return {
        "requested": int(value.get("requested") or 0),
        "resolved": int(value.get("resolved") or 0),
        "missing": list(value.get("missing") or []),
    }


def hydrate_candidate(
    candidate: dict[str, Any],
    final_by_name: dict[str, dict[str, Any]],
    observation_index: dict[str, dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    output = copy.deepcopy(candidate)
    name = str(candidate.get("organism_name") or "")
    source = str(candidate.get("candidate_source") or candidate.get("evidence_source") or "")
    location = "none"
    detail: dict[str, Any] | None = None

    if source == "hospital_only":
        location = "hospital_evidence_detail"
        detail = enrich_hospital_evidence_detail(
            candidate.get("hospital_evidence_detail"), observation_index
        )
        output["hospital_evidence_detail"] = detail
    elif isinstance(candidate.get("hospital_profile"), dict):
        location = "hospital_profile.hospital_evidence_detail"
        frozen_profile = candidate["hospital_profile"]
        enriched_profile = hospital_profile(name, final_by_name, observation_index)
        frozen_semantics = _semantic_profile(frozen_profile)
        enriched_semantics = _semantic_profile(enriched_profile)
        if frozen_semantics != enriched_semantics:
            raise ValueError(
                f"Hospital profile semantics changed for {name}: "
                f"frozen={frozen_semantics!r}, enriched={enriched_semantics!r}"
            )
        output["hospital_profile"] = copy.deepcopy(frozen_profile)
        detail = enriched_profile.get("hospital_evidence_detail")
        output["hospital_profile"]["hospital_evidence_detail"] = detail
    elif isinstance(candidate.get("hospital_evidence_detail"), dict):
        location = "hospital_evidence_detail"
        detail = enrich_hospital_evidence_detail(
            candidate.get("hospital_evidence_detail"), observation_index
        )
        output["hospital_evidence_detail"] = detail

    resolution = _resolution(detail)
    output["evidence_rehydration"] = {
        "stage": SCHEMA_VERSION,
        "location": location,
        "semantic_decision_fields_changed": False,
        "source_observation_resolution": resolution,
    }
    audit = {
        "organism_name": name,
        "location": location,
        **resolution,
    }
    return output, audit


def _decision_signature(candidate: dict[str, Any]) -> dict[str, Any]:
    return {
        key: copy.deepcopy(candidate.get(key))
        for key in (
            "organism_name",
            "organism_key",
            "clinical_decision",
            "clinical_level",
            "formal_pick_allowed",
            "selection_role",
            "analytic_route",
            "history_adjusted_route",
            "history_route_changed",
            "history_rule_ids",
            "history_evidence_ids",
            "consumed_evidence_ids",
        )
    }


def run(
    clinical_root: Path,
    scorer_root: Path,
    output_dir: Path,
    *,
    patients: set[int] | None = None,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    clinical_summary_path = clinical_root / "summary.json"
    clinical_summary = read_json(clinical_summary_path)
    paths = sorted(
        (clinical_root / "patient_outputs").glob(
            "NGS_patient_*_test_aware_clinical_shadow.json"
        ),
        key=patient_number,
    )
    if patients is not None:
        paths = [path for path in paths if patient_number(path) in patients]
    if not paths:
        raise ValueError("No selected clinical shadow patient outputs found")

    output_dir.mkdir(parents=True, exist_ok=True)
    patient_dir = output_dir / "patient_outputs"
    patient_dir.mkdir()
    counts: Counter[str] = Counter()
    audit_rows: list[dict[str, Any]] = []

    for clinical_path in paths:
        patient = patient_number(clinical_path)
        scorer_path = (
            scorer_root
            / "patient_outputs"
            / f"NGS_patient_{patient}_test_aware_deterministic_shadow.json"
        )
        if not scorer_path.is_file():
            raise FileNotFoundError(scorer_path)
        scorer = read_json(scorer_path)
        hospital_summary_path = Path(
            str((scorer.get("source_files") or {}).get("hospital_summary") or "")
        )
        if not hospital_summary_path.is_file():
            raise FileNotFoundError(hospital_summary_path)
        expected_hospital_hash = str(
            (scorer.get("source_files") or {}).get("hospital_summary_sha256") or ""
        )
        actual_hospital_hash = sha256_file(hospital_summary_path)
        if expected_hospital_hash and expected_hospital_hash != actual_hospital_hash:
            raise ValueError(
                f"Hospital summary hash mismatch for P{patient}: {hospital_summary_path}"
            )

        final_summary = read_json(hospital_summary_path)
        final_by_name = legacy.final_candidate_map(final_summary)
        observation_index = hospital_observation_index(final_summary)
        source = read_json(clinical_path)
        original = list(source.get("all_forwarded_candidates") or [])
        hydrated: list[dict[str, Any]] = []
        for candidate in original:
            before = _decision_signature(candidate)
            item, audit = hydrate_candidate(candidate, final_by_name, observation_index)
            if _decision_signature(item) != before:
                raise ValueError(
                    f"Evidence rehydration changed a frozen decision for P{patient} "
                    f"{candidate.get('organism_name')}"
                )
            hydrated.append(item)
            audit_rows.append({"patient_id": str(patient), **audit})
            counts["candidates"] += 1
            counts[f"location:{audit['location']}"] += 1
            counts["source_observations_requested"] += audit["requested"]
            counts["source_observations_resolved"] += audit["resolved"]
            counts["source_observations_missing"] += len(audit["missing"])
            if audit["resolved"]:
                counts["candidates_with_resolved_source_observations"] += 1

        payload = copy.deepcopy(source)
        payload["schema_version"] = SCHEMA_VERSION
        payload["all_forwarded_candidates"] = hydrated
        payload["picked_shadow"] = [
            item for item in hydrated if item.get("clinical_decision") == "picked_shadow"
        ]
        payload["review_high_priority"] = [
            item
            for item in hydrated
            if item.get("clinical_decision") == "review_high_priority"
        ]
        payload["review_context_needed"] = [
            item
            for item in hydrated
            if item.get("clinical_decision") == "review_context_needed"
        ]
        payload["source_clinical_shadow"] = {
            "file": str(clinical_path.resolve()),
            "sha256": sha256_file(clinical_path),
        }
        payload["source_scorer_shadow"] = {
            "file": str(scorer_path.resolve()),
            "sha256": sha256_file(scorer_path),
        }
        payload["source_hospital_summary"] = {
            "file": str(hospital_summary_path.resolve()),
            "sha256": actual_hospital_hash,
        }
        write_json(patient_dir / clinical_path.name, payload)
        counts["patients"] += 1

    csv_path = output_dir / "evidence_rehydration_audit.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        fields = [
            "patient_id", "organism_name", "location", "requested", "resolved", "missing"
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in audit_rows:
            writer.writerow({**row, "missing": "|".join(row.get("missing") or [])})

    summary = {
        "schema_version": SCHEMA_VERSION,
        "answer_blind": True,
        "patient_count": len(paths),
        "candidate_count": counts["candidates"],
        "counts": dict(sorted(counts.items())),
        "history_route_root": clinical_summary.get("history_route_root"),
        "source_scorer_root": clinical_summary.get("source_scorer_root"),
        "timeline_root": clinical_summary.get("timeline_root"),
        "phenotype_packet_root": clinical_summary.get("phenotype_packet_root"),
        "policy_file": clinical_summary.get("policy_file"),
        "policy_sha256": clinical_summary.get("policy_sha256"),
        "source_clinical_root": str(clinical_root.resolve()),
        "source_clinical_summary_sha256": sha256_file(clinical_summary_path),
        "source_scorer_shadow_root": str(scorer_root.resolve()),
        "decision_invariant": (
            "Only source hospital observations were attached; all analytical, history, "
            "clinical, selection-role, and formal-pick fields are unchanged."
        ),
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("clinical_root", type=Path)
    parser.add_argument("scorer_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--patients", type=int, nargs="*")
    args = parser.parse_args()
    selected = set(args.patients) if args.patients else None
    print(json.dumps(run(
        args.clinical_root,
        args.scorer_root,
        args.output_dir,
        patients=selected,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

"""Build an answer-blind OBER/RAG bridge from a frozen external legacy run.

The legacy compatibility runner emits deterministic scorer JSON plus frozen
Priority/Context decisions, while the shared OBER builder expects the historic
``deterministic_max`` + ``llm_missed_candidate_review`` merged shape.  This
adapter changes no tier and reads no benchmark answer.  It only serializes the
already-frozen review candidates into that shared interface.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from tools import pathogen_normalization as pathogen_names


BRIDGE_SUFFIX = "legacy_external_ober_bridge_v1"
REVIEW_POLICY_VERSION = "legacy_external_frozen_priority_context_bridge_v1"
TIER_MAP = {
    "Priority": "review_high_priority",
    "Context": "review_context_needed",
}


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def patient_id_from_dir(path: Path) -> str:
    name = path.name
    prefix = "NGS_patient_"
    suffix = "_json"
    if not (name.startswith(prefix) and name.endswith(suffix)):
        raise ValueError(f"Unexpected patient directory: {path}")
    return name[len(prefix) : -len(suffix)]


def patient_id_from_scorer(path: Path) -> str:
    match = re.fullmatch(r"NGS_patient_(\d+)_.*\.json", path.name)
    if match is None:
        raise ValueError(f"Unexpected scorer filename: {path}")
    return match.group(1)


def candidate_index(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    candidates = payload.get("pathogen_candidates")
    if candidates is None:
        candidates = payload.get("candidates")
    for candidate in candidates or []:
        if not isinstance(candidate, dict):
            continue
        key = pathogen_names.canonical_key(str(candidate.get("organism_name") or ""))
        if key and key not in result:
            result[key] = candidate
    return result


def review_item(row: dict[str, str], candidate: dict[str, Any]) -> dict[str, Any]:
    passed = str(row.get("passed_checks") or "")
    blockers = str(row.get("blockers") or "")
    return {
        "organism_name": str(row.get("organism_name") or "").strip(),
        "classification": candidate.get("classification"),
        "taxonomy_profile": candidate.get("taxonomy_profile"),
        "evidence_snapshot": {
            "rank_priority": candidate.get("rank_priority") or row.get("best_rank"),
            "reads": candidate.get("reads"),
            "reads_tier": candidate.get("reads_tier"),
            "reads_percentile": candidate.get("reads_percentile"),
            "dominance_tier": candidate.get("dominance_tier"),
            "specimen_class": candidate.get("specimen_class"),
            "rank_rule": candidate.get("rank_rule"),
            "support_modules": (candidate.get("key_evidence") or {}).get("support_modules") or [],
        },
        "review_tier_reason": (
            f"Frozen legacy route={row.get('final_label')}; "
            f"passed_checks={passed or 'none'}; blockers={blockers or 'none'}"
        ),
        "risk_if_ignored": "Requires answer-blind OBER review before any optional legacy addition.",
        "why_not_auto_upgrade": "Not in the frozen deterministic Picked tier.",
        "suggested_next_check": "Review organism plausibility, local evidence, competing pathogens, and literature.",
    }


def run(
    *,
    compatibility_run: Path,
    cohort_id: str,
    source_patient_root: Path,
    summary_suffix: str,
    output_dir: Path,
    scorer_root_override: Path | None = None,
) -> dict[str, Any]:
    decisions_path = compatibility_run / "06_frozen_decisions" / "final_decisions.csv"
    scorer_root = scorer_root_override or (
        compatibility_run
        / "01_legacy_selected_dna_scorer"
        / cohort_id
        / "patient_outputs"
    )
    if not decisions_path.is_file():
        raise FileNotFoundError(decisions_path)
    if not scorer_root.is_dir():
        raise FileNotFoundError(scorer_root)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_dir}")

    rows_by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    with decisions_path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            row_cohort_id = str(row.get("cohort_id") or "").strip()
            if row_cohort_id and row_cohort_id != cohort_id:
                continue
            if row.get("final_label") not in TIER_MAP:
                continue
            rows_by_patient[str(row.get("patient_id") or "")].append(row)

    source_dirs_by_id = {
        patient_id_from_dir(path): path
        for path in source_patient_root.glob("NGS_patient_*_json")
        if path.is_dir()
    }
    scorer_paths_by_id: dict[str, Path] = {}
    for path in scorer_root.glob("NGS_patient_*_*.json"):
        if not path.is_file():
            continue
        patient_id = patient_id_from_scorer(path)
        if patient_id in scorer_paths_by_id:
            raise ValueError(
                f"Multiple scorer payloads for patient {patient_id}: "
                f"{scorer_paths_by_id[patient_id]}, {path}"
            )
        scorer_paths_by_id[patient_id] = path
    if not scorer_paths_by_id:
        raise ValueError(f"No patient scorer files under {scorer_root}")
    missing_source_ids = sorted(
        set(scorer_paths_by_id) - set(source_dirs_by_id), key=int
    )
    if missing_source_ids:
        raise ValueError(
            "Source patient directories are missing for frozen cohort patients: "
            f"{missing_source_ids}"
        )
    patient_ids = sorted(scorer_paths_by_id, key=int)

    bridge_root = output_dir / "patients"
    audit_rows: list[dict[str, Any]] = []
    tier_counts: Counter[str] = Counter()
    source_hashes: list[dict[str, str]] = []
    unresolved: list[dict[str, str]] = []

    for patient_id in patient_ids:
        source_patient_dir = source_dirs_by_id[patient_id]
        scorer_path = scorer_paths_by_id[patient_id]
        summary_path = (
            source_patient_dir
            / "summary_outputs"
            / f"NGS_patient_{patient_id}_{summary_suffix}.json"
        )
        if not scorer_path.is_file():
            raise FileNotFoundError(scorer_path)
        if not summary_path.is_file():
            raise FileNotFoundError(summary_path)
        scorer = read_json(scorer_path)
        summary = read_json(summary_path)
        candidates = candidate_index(scorer)
        review: dict[str, Any] = {
            "review_high_priority": [],
            "review_context_needed": [],
            "review_tiering_policy": {
                "version": REVIEW_POLICY_VERSION,
                "answer_blind": True,
                "source": "frozen legacy Priority/Context decisions",
            },
        }
        for row in rows_by_patient.get(patient_id, []):
            name = str(row.get("organism_name") or "").strip()
            key = pathogen_names.canonical_key(name)
            candidate = candidates.get(key)
            if candidate is None:
                unresolved.append({"patient_id": patient_id, "organism_name": name})
                continue
            tier = TIER_MAP[str(row.get("final_label"))]
            review[tier].append(review_item(row, candidate))
            tier_counts[tier] += 1
            audit_rows.append(
                {
                    "patient_id": patient_id,
                    "organism_name": name,
                    "frozen_label": row.get("final_label"),
                    "ober_source_tier": tier,
                    "candidate_found": True,
                }
            )

        target_summary_dir = bridge_root / source_patient_dir.name / "summary_outputs"
        merged_path = target_summary_dir / f"NGS_patient_{patient_id}_{BRIDGE_SUFFIX}.json"
        copied_summary_path = target_summary_dir / f"NGS_patient_{patient_id}_{summary_suffix}.json"
        write_json(
            merged_path,
            {
                "schema_version": "legacy_external_ober_bridge.v1",
                "answer_blind": True,
                "patient_id": patient_id,
                "deterministic_max": scorer,
                "llm_missed_candidate_review": review,
                "provenance": {
                    "compatibility_run": str(compatibility_run.resolve()),
                    "frozen_decisions": str(decisions_path.resolve()),
                    "source_scorer": str(scorer_path.resolve()),
                },
            },
        )
        write_json(copied_summary_path, summary)
        source_hashes.append(
            {
                "patient_id": patient_id,
                "scorer_sha256": sha256(scorer_path),
                "summary_sha256": sha256(summary_path),
                "bridge_sha256": sha256(merged_path),
            }
        )

    if unresolved:
        raise RuntimeError(f"Frozen candidates missing from scorer payload: {unresolved}")

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "candidate_bridge_audit.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        fields = [
            "patient_id",
            "organism_name",
            "frozen_label",
            "ober_source_tier",
            "candidate_found",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(audit_rows)

    manifest = {
        "schema_version": "legacy_external_ober_bridge_manifest.v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "answer_source_read": False,
        "answer_blind": True,
        "cohort_id": cohort_id,
        "compatibility_run": str(compatibility_run.resolve()),
        "scorer_root": str(scorer_root.resolve()),
        "source_patient_root": str(source_patient_root.resolve()),
        "summary_suffix": summary_suffix,
        "bridge_suffix": BRIDGE_SUFFIX,
        "patient_count": len(patient_ids),
        "candidate_count": len(audit_rows),
        "tier_counts": dict(tier_counts),
        "unresolved_candidate_count": 0,
        "source_hashes": source_hashes,
        "constraints": [
            "No benchmark answer file is read.",
            "Frozen Picked, Priority, and Context decisions are not changed.",
            "Only frozen Priority and Context candidates enter the OBER queue.",
            "OBER output remains separate and cannot overwrite deterministic Picked.",
        ],
    }
    write_json(output_dir / "bridge_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compatibility-run", type=Path, required=True)
    parser.add_argument("--cohort-id", default="clinical11")
    parser.add_argument("--source-patient-root", type=Path, required=True)
    parser.add_argument("--summary-suffix", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--scorer-root",
        type=Path,
        help=(
            "Optional frozen scorer patient_outputs directory. Defaults to the "
            "two-hospital selected-DNA compatibility layout."
        ),
    )
    args = parser.parse_args()
    manifest = run(
        compatibility_run=args.compatibility_run,
        cohort_id=args.cohort_id,
        source_patient_root=args.source_patient_root,
        summary_suffix=args.summary_suffix,
        output_dir=args.output_dir,
        scorer_root_override=args.scorer_root,
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

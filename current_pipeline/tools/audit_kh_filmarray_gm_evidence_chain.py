"""Audit KH FilmArray/GM raw -> observation -> summary -> candidate preservation."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools import build_deterministic_summary as summary_builder
from tools import normalized_agent_fallback
from tools.pathogen_normalization import (
    approved_group_member_match,
    canonical_key,
)


DEFAULT_ROOT = Path("outputs/patient_info_KH_0728_2Days")
DEFAULT_INTEGRATED_ROOT = Path(
    "outputs/runs/2026-09-29_KH_integrated_decision_pipeline_v4_filmarray_gm_shadow"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-29_KH_filmarray_gm_evidence_chain_audit_v2"
)
AGENT_TAG = "evidence_v2_full33_20260917"


def read_json_any(path: Path) -> Any:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    return [value] if isinstance(value, dict) else []


def patient_number(path: Path) -> int:
    match = re.fullmatch(r"NGS_patient_(\d+)_json", path.name)
    if not match:
        raise ValueError(path.name)
    return int(match.group(1))


def patient_dirs(root: Path) -> list[Path]:
    output = [
        path for path in root.iterdir()
        if path.is_dir() and re.fullmatch(r"NGS_patient_\d+_json", path.name)
    ]
    return sorted(output, key=patient_number)


def source_path(patient_dir: Path, suffix: str) -> Path:
    patient = patient_number(patient_dir)
    return patient_dir / f"NGS_patient_{patient}_{suffix}.json"


def agent_path(patient_dir: Path) -> Path:
    patient = patient_number(patient_dir)
    return (
        patient_dir
        / "agent_outputs"
        / f"NGS_patient_{patient}_filmarray_gmtest_agent_{AGENT_TAG}.json"
    )


def integrated_candidates(root: Path, patient: int) -> list[dict[str, Any]]:
    path = root / "patient_results" / f"patient_{patient}_integrated_decision.json"
    payload = read_json_any(path)
    return list(payload.get("all_candidates") or []) if isinstance(payload, dict) else []


def candidate_matches(label_name: str, candidates: list[dict[str, Any]]) -> list[str]:
    output: list[str] = []
    for candidate in candidates:
        candidate_name = str(candidate.get("organism_name") or "")
        if (
            canonical_key(label_name) == canonical_key(candidate_name)
            or approved_group_member_match(label_name, candidate_name)
        ):
            output.append(candidate_name)
    return list(dict.fromkeys(output))


def linked_observation_ids(payload: dict[str, Any]) -> tuple[set[str], set[str]]:
    label_ids: set[str] = set()
    resistance_ids: set[str] = set()
    for item in payload.get("organism_labels") or []:
        label_ids.update(str(value) for value in item.get("evidence_observation_ids") or [])
    for item in payload.get("resistance_findings") or []:
        resistance_ids.update(str(value) for value in item.get("evidence_observation_ids") or [])
    return label_ids, resistance_ids


def audit(
    root: Path,
    integrated_root: Path,
    output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    chain_rows: list[dict[str, Any]] = []
    patient_rows: list[dict[str, Any]] = []
    summary_only_positive_labels: list[dict[str, Any]] = []
    totals: Counter[str] = Counter()

    for patient_dir in patient_dirs(root):
        patient = patient_number(patient_dir)
        filmarray = rows(read_json_any(source_path(patient_dir, "filmarray")))
        gm = rows(read_json_any(source_path(patient_dir, "gm_test")))
        source_rows = [("filmarray", index, row) for index, row in enumerate(filmarray)]
        source_rows.extend(("gm_test", index, row) for index, row in enumerate(gm))
        agent = read_json_any(agent_path(patient_dir)) or {}
        fallback = normalized_agent_fallback.filmarray_gm_agent_from_sources(
            filmarray, gm
        )
        agent_observations = {
            str(item.get("observation_id")): item
            for item in agent.get("assay_observations") or []
            if isinstance(item, dict) and item.get("observation_id")
        }
        fallback_observations = {
            str(item.get("observation_id")): item
            for item in fallback.get("assay_observations") or []
            if isinstance(item, dict) and item.get("observation_id")
        }
        agent_label_ids, agent_resistance_ids = linked_observation_ids(agent)
        fallback_label_ids, fallback_resistance_ids = linked_observation_ids(fallback)
        summary = summary_builder.build_summary(
            patient_dir, agent_suffix_tag=AGENT_TAG
        )
        summary_names = [
            str(item.get("organism_name") or "")
            for item in summary.get("hospital_organism_evidence") or []
        ]
        candidates = integrated_candidates(integrated_root, patient)
        fallback_label_by_id: dict[str, str] = {}
        for item in fallback.get("organism_labels") or []:
            for observation_id in item.get("evidence_observation_ids") or []:
                fallback_label_by_id[str(observation_id)] = str(
                    item.get("organism_name") or ""
                )

        patient_errors: list[str] = []
        expected_count = len(source_rows)
        if len(agent_observations) != expected_count:
            patient_errors.append("agent_observation_count_mismatch")
        if len(fallback_observations) != expected_count:
            patient_errors.append("fallback_observation_count_mismatch")

        for combined_index, (source_type, source_index, row) in enumerate(source_rows):
            observation_id = f"ASSAY-{combined_index + 1:03d}"
            agent_observation = agent_observations.get(observation_id) or {}
            fallback_observation = fallback_observations.get(observation_id) or {}
            raw_target = str(row.get("target") or row.get("test") or "")
            label_name = fallback_label_by_id.get(observation_id, "")
            summary_matches = [
                name for name in summary_names
                if canonical_key(name) == canonical_key(label_name)
            ] if label_name else []
            matches = candidate_matches(label_name, candidates) if label_name else []
            record = {
                "patient_id": patient,
                "source_type": source_type,
                "source_record_index": source_index,
                "observation_id": observation_id,
                "raw_target": raw_target,
                "raw_result": row.get("result") or row.get("value") or "",
                "agent_observation_present": bool(agent_observation),
                "agent_target_matches": bool(
                    agent_observation
                    and canonical_key(agent_observation.get("target"))
                    == canonical_key(raw_target)
                ),
                "fallback_observation_present": bool(fallback_observation),
                "fallback_target_matches": bool(
                    fallback_observation
                    and canonical_key(fallback_observation.get("target"))
                    == canonical_key(raw_target)
                ),
                "agent_detection_status": agent_observation.get("detection_status"),
                "fallback_detection_status": fallback_observation.get("detection_status"),
                "agent_label_link": observation_id in agent_label_ids,
                "agent_resistance_link": observation_id in agent_resistance_ids,
                "fallback_label_link": observation_id in fallback_label_ids,
                "fallback_resistance_link": observation_id in fallback_resistance_ids,
                "normalized_label_name": label_name,
                "summary_match": "|".join(summary_matches),
                "candidate_match": "|".join(matches),
            }
            chain_rows.append(record)
            totals["agent_target_mismatches"] += int(
                not record["agent_target_matches"]
            )
            totals["fallback_target_mismatches"] += int(
                not record["fallback_target_matches"]
            )
            totals["positive_label_observations"] += int(
                record["fallback_label_link"]
            )
            totals["summary_matched_positive_label_observations"] += int(
                bool(record["summary_match"])
            )
            totals["candidate_matched_positive_label_observations"] += int(
                bool(record["candidate_match"])
            )
            if (
                record["fallback_label_link"]
                and record["summary_match"]
                and not record["candidate_match"]
            ):
                summary_only_positive_labels.append({
                    "patient_id": patient,
                    "observation_id": observation_id,
                    "organism_name": label_name,
                    "raw_result": record["raw_result"],
                    "interpretation": (
                        "preserved_in_hospital_summary_but_not_forwarded_as_"
                        "a_scorer_candidate"
                    ),
                })
            if not record["agent_observation_present"]:
                patient_errors.append(f"missing_agent_observation:{observation_id}")
            if not record["fallback_observation_present"]:
                patient_errors.append(f"missing_fallback_observation:{observation_id}")
            if not record["agent_target_matches"]:
                patient_errors.append(f"agent_target_mismatch:{observation_id}")
            if not record["fallback_target_matches"]:
                patient_errors.append(f"fallback_target_mismatch:{observation_id}")

        unresolved_agent_ids = sorted(
            (agent_label_ids | agent_resistance_ids) - set(agent_observations)
        )
        unresolved_fallback_ids = sorted(
            (fallback_label_ids | fallback_resistance_ids) - set(fallback_observations)
        )
        if unresolved_agent_ids:
            patient_errors.append("unresolved_agent_links")
        if unresolved_fallback_ids:
            patient_errors.append("unresolved_fallback_links")
        patient_rows.append({
            "patient_id": patient,
            "filmarray_raw_count": len(filmarray),
            "gm_raw_count": len(gm),
            "raw_total": expected_count,
            "agent_observation_count": len(agent_observations),
            "fallback_observation_count": len(fallback_observations),
            "agent_label_count": len(agent.get("organism_labels") or []),
            "fallback_label_count": len(fallback.get("organism_labels") or []),
            "unresolved_agent_link_count": len(unresolved_agent_ids),
            "unresolved_fallback_link_count": len(unresolved_fallback_ids),
            "complete": not patient_errors,
            "errors": "|".join(dict.fromkeys(patient_errors)),
        })
        totals["patients"] += 1
        totals["filmarray_raw_records"] += len(filmarray)
        totals["gm_raw_records"] += len(gm)
        totals["raw_records"] += expected_count
        totals["agent_observations"] += len(agent_observations)
        totals["fallback_observations"] += len(fallback_observations)
        totals["complete_patients"] += int(not patient_errors)
        totals["incomplete_patients"] += int(bool(patient_errors))

    def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
        fields = list(records[0]) if records else []
        with path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)

    write_csv(output_dir / "filmarray_gm_source_chain.csv", chain_rows)
    write_csv(output_dir / "filmarray_gm_patient_completeness.csv", patient_rows)
    summary = {
        "schema_version": "kh_filmarray_gm_evidence_chain_audit.v1",
        "answer_blind": True,
        "source_root": str(root.resolve()),
        "integrated_root": str(integrated_root.resolve()),
        "agent_tag": AGENT_TAG,
        "totals": dict(totals),
        "source_complete_agent": (
            totals["agent_observations"] == totals["raw_records"]
            and totals["incomplete_patients"] == 0
        ),
        "source_complete_fallback": (
            totals["fallback_observations"] == totals["raw_records"]
            and all(row["unresolved_fallback_link_count"] == 0 for row in patient_rows)
        ),
        "patient_failures": [
            row for row in patient_rows if not row["complete"]
        ],
        "summary_only_positive_labels": summary_only_positive_labels,
        "constraints": [
            "Every raw FilmArray, GM, serology, or antigen row must have one observation.",
            "Negative, equivocal, and invalid rows remain observations but do not create organism labels.",
            "AMR targets remain resistance findings and never become organism candidates.",
            "Candidate matching permits exact/alias or explicitly approved group/member relationships only.",
            "This audit does not use benchmark answers and does not modify decision tiers.",
        ],
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument(
        "--integrated-root", type=Path, default=DEFAULT_INTEGRATED_ROOT
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(
        audit(args.root, args.integrated_root, args.output_dir),
        ensure_ascii=False,
        indent=2,
    ))


if __name__ == "__main__":
    main()

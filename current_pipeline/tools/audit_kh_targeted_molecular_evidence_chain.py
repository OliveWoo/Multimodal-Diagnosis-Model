"""Audit KH targeted-molecular availability and assay-source semantics.

The KH cohort contains multiplex FilmArray rows but currently no dedicated
``*_molecular_microbiology.json`` source.  This audit keeps those two evidence
axes distinct and verifies that any targeted-respiratory promotion is traceable
to its real module, panel, observation, and event window.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any

from tools import normalized_agent_fallback
from tools.pathogen_normalization import approved_group_member_match, canonical_key


DEFAULT_ROOT = Path("outputs/patient_info_KH_0728_2Days")
DEFAULT_INTEGRATED_ROOT = Path(
    "outputs/runs/2026-09-29_KH_integrated_decision_pipeline_v4_filmarray_gm_shadow"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-29_KH_targeted_molecular_evidence_chain_audit_v1"
)
PROMOTION_ROUTE = "bacterial_targeted_respiratory_assay_bridge"


def read_json(path: Path) -> Any:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    return [value] if isinstance(value, dict) else []


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)", path.name)
    if not match:
        raise ValueError(path.name)
    return int(match.group(1))


def patient_dirs(root: Path) -> list[Path]:
    return sorted(
        [
            path
            for path in root.iterdir()
            if path.is_dir() and re.fullmatch(r"NGS_patient_\d+_json", path.name)
        ],
        key=patient_number,
    )


def write_csv(path: Path, records: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def matched_candidate(target: str, organism: str) -> bool:
    return (
        canonical_key(target) == canonical_key(organism)
        or approved_group_member_match(target, organism)
        or approved_group_member_match(organism, target)
    )


def audit(root: Path, integrated_root: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    completeness: list[dict[str, Any]] = []
    panel_targets: list[dict[str, Any]] = []
    promotions: list[dict[str, Any]] = []

    for patient_dir in patient_dirs(root):
        patient = patient_number(patient_dir)
        raw_molecular_path = patient_dir / f"NGS_patient_{patient}_molecular_microbiology.json"
        molecular_agent_paths = sorted(
            (patient_dir / "agent_outputs").glob(
                f"NGS_patient_{patient}_molecular_microbiology_agent*.json"
            )
        ) if (patient_dir / "agent_outputs").is_dir() else []
        raw_molecular = read_json(raw_molecular_path)
        molecular_rows = rows(raw_molecular)
        molecular_fallback = (
            normalized_agent_fallback.molecular_microbiology_agent_from_source(
                raw_molecular
            )
            if raw_molecular is not None
            else None
        )
        completeness.append(
            {
                "patient_id": patient,
                "dedicated_molecular_source_present": raw_molecular_path.is_file(),
                "dedicated_molecular_raw_row_count": len(molecular_rows),
                "dedicated_molecular_agent_file_count": len(molecular_agent_paths),
                "fallback_observation_count": len(
                    (molecular_fallback or {}).get("molecular_observations") or []
                ),
                "fallback_positive_label_count": len(
                    (molecular_fallback or {}).get("organism_labels") or []
                ),
                "availability_interpretation": (
                    "available"
                    if raw_molecular_path.is_file()
                    else "not_available_not_negative"
                ),
            }
        )

        filmarray_path = patient_dir / f"NGS_patient_{patient}_filmarray.json"
        filmarray_payload = normalized_agent_fallback.filmarray_agent_from_source(
            read_json(filmarray_path)
        )
        observations = {
            str(item.get("observation_id")): item
            for item in filmarray_payload.get("assay_observations") or []
            if isinstance(item, dict)
        }
        for label in filmarray_payload.get("organism_labels") or []:
            for observation_id in label.get("evidence_observation_ids") or []:
                observation = observations.get(str(observation_id)) or {}
                panel_targets.append(
                    {
                        "patient_id": patient,
                        "observation_id": observation_id,
                        "target": observation.get("target") or label.get("organism_name"),
                        "detection_status": observation.get("detection_status"),
                        "panel_or_assay": observation.get("panel_or_assay"),
                        "test_type": observation.get("test_type"),
                        "specimen_type": observation.get("specimen_type"),
                        "reported_time": observation.get("reported_time"),
                        "causative_level": label.get("causative_level"),
                        "evidence_axis": "multiplex_filmarray_target",
                    }
                )

    decision_dir = integrated_root / "decision_engine" / "patient_outputs"
    for path in sorted(decision_dir.glob("NGS_patient_*_test_aware_clinical_shadow.json"), key=patient_number):
        patient = patient_number(path)
        payload = read_json(path) or {}
        for candidate in payload.get("all_forwarded_candidates") or []:
            gate = candidate.get("promotion_gate") or {}
            if gate.get("route") != PROMOTION_ROUTE:
                continue
            hospital = candidate.get("hospital_profile") or {}
            detail = hospital.get("hospital_evidence_detail") or {}
            source_observations = detail.get("source_observations") or []
            axes = gate.get("axes") or {}
            timing = axes.get("direct_evidence_timing_profile") or {}
            matching_observations = [
                item
                for item in source_observations
                if matched_candidate(
                    str(item.get("target") or ""),
                    str(candidate.get("organism_name") or ""),
                )
            ]
            promotions.append(
                {
                    "patient_id": patient,
                    "organism_name": candidate.get("organism_name"),
                    "promotion_route": gate.get("route"),
                    "outcome": gate.get("outcome"),
                    "evidence_modules": "|".join(detail.get("evidence_modules") or []),
                    "molecular_microbiology_level": (
                        hospital.get("module_level_summary") or {}
                    ).get("molecular_microbiology"),
                    "filmarray_gmtest_level": (
                        hospital.get("module_level_summary") or {}
                    ).get("filmarray_gmtest"),
                    "matching_observation_count": len(matching_observations),
                    "observation_ids": "|".join(
                        str(item.get("observation_id") or "")
                        for item in matching_observations
                    ),
                    "test_types": "|".join(
                        dict.fromkeys(
                            str(item.get("test_type") or "")
                            for item in matching_observations
                        )
                    ),
                    "panels_or_assays": "|".join(
                        dict.fromkeys(
                            str(item.get("panel_or_assay") or "")
                            for item in matching_observations
                        )
                    ),
                    "targets": "|".join(
                        str(item.get("target") or "")
                        for item in matching_observations
                    ),
                    "event_window_hours": timing.get("event_window_hours"),
                    "event_aligned_positive_count": timing.get(
                        "event_aligned_positive_count"
                    ),
                    "semantic_interpretation": (
                        "organism_specific_target_inside_multiplex_filmarray_panel"
                        if matching_observations
                        and all(
                            item.get("test_type") == "filmarray"
                            for item in matching_observations
                        )
                        else "dedicated_or_unresolved_targeted_molecular_assay"
                    ),
                }
            )

    source_present = sum(
        bool(item["dedicated_molecular_source_present"]) for item in completeness
    )
    raw_row_count = sum(
        int(item["dedicated_molecular_raw_row_count"]) for item in completeness
    )
    promotion_patient_ids = sorted(int(item["patient_id"]) for item in promotions)
    all_panel_semantics = all(
        item["semantic_interpretation"]
        == "organism_specific_target_inside_multiplex_filmarray_panel"
        and item["evidence_modules"] == "filmarray_gmtest"
        and item["molecular_microbiology_level"] == "Not_available"
        and int(item["matching_observation_count"]) >= 1
        for item in promotions
    )
    assertions = {
        "all_33_patients_audited": len(completeness) == 33,
        "no_dedicated_molecular_source_files_in_kh": source_present == 0,
        "no_dedicated_molecular_rows_in_kh": raw_row_count == 0,
        "targeted_bridge_is_two_expected_cases": promotion_patient_ids == [12, 16],
        "targeted_bridge_sources_are_filmarray_not_dedicated_pcr": all_panel_semantics,
        "multiplex_positive_targets_remain_traceable": bool(panel_targets),
    }
    summary = {
        "schema_version": "kh_targeted_molecular_evidence_chain_audit.v1",
        "answer_blind": True,
        "patient_root": str(root.resolve()),
        "integrated_root": str(integrated_root.resolve()),
        "totals": {
            "patients": len(completeness),
            "dedicated_molecular_source_files": source_present,
            "dedicated_molecular_raw_rows": raw_row_count,
            "dedicated_molecular_agent_files": sum(
                int(item["dedicated_molecular_agent_file_count"])
                for item in completeness
            ),
            "multiplex_filmarray_positive_organism_targets": len(panel_targets),
            "targeted_respiratory_bridge_promotions": len(promotions),
        },
        "promotion_patient_ids": promotion_patient_ids,
        "assertions": assertions,
        "pass": all(assertions.values()),
        "interpretation": {
            "dedicated_targeted_molecular": (
                "not available in the KH 33-patient source; absence is not a negative test"
            ),
            "p12_p16": (
                "the supporting target is inside a BioFire FilmArray Pneumonia Panel "
                "and must not be counted again as an independent PCR module"
            ),
        },
    }

    write_csv(
        output_dir / "patient_availability.csv",
        completeness,
        list(completeness[0]) if completeness else [],
    )
    write_csv(
        output_dir / "multiplex_panel_positive_targets.csv",
        panel_targets,
        list(panel_targets[0]) if panel_targets else ["patient_id"],
    )
    write_csv(
        output_dir / "promotion_assay_semantics.csv",
        promotions,
        list(promotions[0]) if promotions else ["patient_id"],
    )
    (output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--integrated-root", type=Path, default=DEFAULT_INTEGRATED_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = audit(args.root, args.integrated_root, args.output_dir)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

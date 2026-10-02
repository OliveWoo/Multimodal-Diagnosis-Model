"""Audit KH underlying disease, medication timing, and deterministic H tiers."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools import build_deterministic_summary as summary_builder
from tools import evidence_preservation
from tools.audit_kh_cbc_other_lab_evidence_chain import promotion_guardrail_rows


DEFAULT_ROOT = Path("outputs/patient_info_KH_0728_2Days")
DEFAULT_INTEGRATED_ROOT = Path(
    "outputs/runs/2026-09-29_KH_integrated_decision_pipeline_v4_filmarray_gm_shadow"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-30_KH_underlying_medication_evidence_chain_audit_v2"
)
AGENT_TAG = "evidence_v2_full33_20260917"


def read_json(path: Path) -> Any:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8-sig"))


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


def source_path(patient_dir: Path, suffix: str) -> Path:
    patient = patient_number(patient_dir)
    return patient_dir / f"NGS_patient_{patient}_{suffix}.json"


def agent_path(patient_dir: Path) -> Path:
    patient = patient_number(patient_dir)
    return (
        patient_dir
        / "agent_outputs"
        / f"NGS_patient_{patient}_cbc_other_lab_agent_{AGENT_TAG}.json"
    )


def rows(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, list):
        return [item for item in value if isinstance(item, dict)]
    return [value] if isinstance(value, dict) else []


def write_csv(path: Path, records: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def audit(root: Path, integrated_root: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    patients: list[dict[str, Any]] = []
    medications: list[dict[str, Any]] = []
    tier_differences: list[dict[str, Any]] = []
    high_risk_basis: list[dict[str, Any]] = []
    uncertain_mentions: list[dict[str, Any]] = []
    json_errors: list[dict[str, Any]] = []
    agent_h_counts: Counter[str] = Counter()
    deterministic_h_counts: Counter[str] = Counter()
    support_status_counts: Counter[str] = Counter()

    for patient_dir in patient_dirs(root):
        patient = patient_number(patient_dir)
        try:
            underlying = read_json(source_path(patient_dir, "underlying"))
            admission = read_json(source_path(patient_dir, "admission_diagnosis"))
            ranked = read_json(source_path(patient_dir, "mNGS_ranked_candidates"))
            agent = read_json(agent_path(patient_dir)) or {}
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            json_errors.append({"patient_id": patient, "error": str(exc)})
            continue

        sample_time = evidence_preservation.mngs_sample_time_from_ranked(ranked)
        medication = evidence_preservation.extract_host_risk_medications(
            {"underlying": underlying, "admission_diagnosis": admission},
            sample_time=sample_time,
        )
        profile = evidence_preservation.deterministic_host_risk_profile(
            underlying,
            admission,
            sample_time=sample_time,
            medication_profile=medication,
        )
        summary = summary_builder.build_summary(
            patient_dir, agent_suffix_tag=AGENT_TAG
        )
        host = summary.get("host_context") or {}
        agent_state = agent.get("host_state") or {}
        agent_h = str(agent_state.get("immunocompromise_tier") or "Unknown")
        deterministic_h = str(profile.get("immunocompromise_tier") or "Unknown")
        agent_h_counts[agent_h] += 1
        deterministic_h_counts[deterministic_h] += 1
        support_status_counts[str(profile["opportunistic_host_support_status"])] += 1

        for item in medication["medications"]:
            medications.append(
                {
                    "patient_id": patient,
                    "sample_time": sample_time,
                    "medication_name": item["medication_name"],
                    "medication_class": item["medication_class"],
                    "role": item["role"],
                    "start_date": item["start_date"],
                    "end_date": item["end_date"],
                    "ongoing_at_source": item["ongoing_at_source"],
                    "temporal_status": item["temporal_status"],
                    "days_from_start_to_sample": item[
                        "days_from_start_to_sample"
                    ],
                    "days_from_end_to_sample": item["days_from_end_to_sample"],
                    "dose": item["dose"],
                    "source_path": item["source_path"],
                    "source_text": item["source_text"],
                }
            )

        for item in profile["tier_triggers"]:
            if deterministic_h in {"H2", "H3"}:
                high_risk_basis.append(
                    {
                        "patient_id": patient,
                        "deterministic_H": deterministic_h,
                        "host_support_eligible": profile["host_support_eligible"],
                        "support_status": profile[
                            "opportunistic_host_support_status"
                        ],
                        "rule": item["rule"],
                        "assertion": item["assertion"],
                        "matched_text": item["matched_text"],
                        "source_section": item["source_section"],
                        "source_path": item["source_path"],
                        "source_text": item["source_text"],
                    }
                )

        for item in profile["uncertain_or_negated_mentions"]:
            uncertain_mentions.append(
                {
                    "patient_id": patient,
                    "tier_if_confirmed": item["tier"],
                    "rule": item["rule"],
                    "assertion": item["assertion"],
                    "matched_text": item["matched_text"],
                    "source_path": item["source_path"],
                    "source_text": item["source_text"],
                }
            )

        if agent_h != deterministic_h:
            tier_differences.append(
                {
                    "patient_id": patient,
                    "agent_H": agent_h,
                    "deterministic_H": deterministic_h,
                    "agent_V": agent_state.get("host_vulnerability_tier"),
                    "deterministic_V": host.get("host_vulnerability_tier"),
                    "opportunistic_host_support": host.get(
                        "opportunistic_host_support"
                    ),
                    "support_status": host.get(
                        "opportunistic_host_support_status"
                    ),
                    "tier_triggers": json.dumps(
                        profile["tier_triggers"],
                        ensure_ascii=False,
                        separators=(",", ":"),
                    ),
                }
            )

        patients.append(
            {
                "patient_id": patient,
                "underlying_source_available": underlying is not None,
                "admission_source_available": admission is not None,
                "underlying_row_count": len(rows(underlying)),
                "admission_row_count": len(rows(admission)),
                "host_history_observation_count": len(
                    profile["host_history_observations"]
                ),
                "sample_time": sample_time,
                "agent_H": agent_h,
                "deterministic_H": deterministic_h,
                "agent_A": agent_state.get("acute_instability_tier"),
                "deterministic_A": host.get("acute_instability_tier"),
                "agent_V": agent_state.get("host_vulnerability_tier"),
                "deterministic_V": host.get("host_vulnerability_tier"),
                "host_support_eligible": profile["host_support_eligible"],
                "summary_host_support": host.get("opportunistic_host_support"),
                "support_status": profile["opportunistic_host_support_status"],
                "medication_observation_count": medication[
                    "medication_observation_count"
                ],
                "distinct_risk_medication_count": medication[
                    "risk_medication_count"
                ],
                "distinct_risk_medications": "|".join(
                    medication["distinct_risk_medications"]
                ),
                "medication_evidence_strength": medication[
                    "host_risk_evidence_strength"
                ],
                "medication_data_gaps": "|".join(medication["data_gaps"]),
                "tier_trigger_count": len(profile["tier_triggers"]),
                "uncertain_mention_count": len(
                    profile["uncertain_or_negated_mentions"]
                ),
            }
        )

    promotion_rows = promotion_guardrail_rows(integrated_root)
    host_only_picks = [row for row in promotion_rows if row["host_only_pick"]]
    by_patient = {int(item["patient_id"]): item for item in patients}
    med_names_by_patient: dict[int, set[str]] = {}
    for row in medications:
        med_names_by_patient.setdefault(int(row["patient_id"]), set()).add(
            str(row["medication_name"])
        )
    p11_methotrexate = [
        row
        for row in medications
        if row["patient_id"] == 11 and row["medication_name"] == "methotrexate"
    ]
    assertions = {
        "all_33_patients_audited": len(patients) == 33,
        "no_json_read_errors": not json_errors,
        "underlying_and_admission_sources_available_for_all_patients": all(
            row["underlying_source_available"]
            and row["admission_source_available"]
            for row in patients
        ),
        "all_h2_h3_have_traceable_tier_basis": all(
            any(
                basis["patient_id"] == row["patient_id"]
                for basis in high_risk_basis
            )
            for row in patients
            if row["deterministic_H"] in {"H2", "H3"}
        ),
        "summary_host_support_matches_deterministic_eligibility": all(
            bool(row["summary_host_support"]) == bool(row["host_support_eligible"])
            for row in patients
        ),
        "provisional_medication_only_h2_does_not_enable_host_support": bool(
            by_patient.get(16)
            and by_patient[16]["deterministic_H"] == "H2"
            and not by_patient[16]["summary_host_support"]
        ),
        "parathyroid_autotransplant_does_not_trigger_h3": bool(
            by_patient.get(6) and by_patient[6]["deterministic_H"] != "H3"
        ),
        "copd_is_preserved_as_h1": bool(
            by_patient.get(37) and by_patient[37]["deterministic_H"] == "H1"
        ),
        "decitabine_is_preserved": "decitabine" in med_names_by_patient.get(10, set()),
        "pemetrexed_and_cisplatin_are_preserved": {
            "pemetrexed",
            "cisplatin",
        }.issubset(med_names_by_patient.get(17, set())),
        "methotrexate_does_not_inherit_prednisolone_dose": bool(
            p11_methotrexate
            and all(row["dose"] == "Unknown" for row in p11_methotrexate)
        ),
        "no_picked_candidate_created_from_host_context_alone": not host_only_picks,
    }
    report = {
        "schema_version": "kh_underlying_medication_evidence_chain_audit.v1",
        "answer_blind": True,
        "patient_root": str(root.resolve()),
        "integrated_root": str(integrated_root.resolve()),
        "totals": {
            "patients": len(patients),
            "underlying_rows": sum(row["underlying_row_count"] for row in patients),
            "admission_rows": sum(row["admission_row_count"] for row in patients),
            "host_history_observations": sum(
                row["host_history_observation_count"] for row in patients
            ),
            "medication_observations": len(medications),
            "patients_with_risk_medications": sum(
                row["distinct_risk_medication_count"] > 0 for row in patients
            ),
            "h_tier_differences": len(tier_differences),
            "h2_h3_patients": sum(
                row["deterministic_H"] in {"H2", "H3"} for row in patients
            ),
            "picked_candidates_audited": len(promotion_rows),
            "host_only_picks": len(host_only_picks),
        },
        "agent_H_counts": dict(agent_h_counts),
        "deterministic_H_counts": dict(deterministic_h_counts),
        "support_status_counts": dict(support_status_counts),
        "assertions": assertions,
        "pass": all(assertions.values()),
        "json_errors": json_errors,
    }

    write_csv(
        output_dir / "patient_host_risk.csv",
        patients,
        list(patients[0]) if patients else ["patient_id"],
    )
    write_csv(
        output_dir / "medication_observations.csv",
        medications,
        list(medications[0]) if medications else ["patient_id"],
    )
    write_csv(
        output_dir / "h_tier_differences.csv",
        tier_differences,
        list(tier_differences[0]) if tier_differences else ["patient_id"],
    )
    write_csv(
        output_dir / "high_risk_basis.csv",
        high_risk_basis,
        list(high_risk_basis[0]) if high_risk_basis else ["patient_id"],
    )
    write_csv(
        output_dir / "uncertain_host_mentions.csv",
        uncertain_mentions,
        list(uncertain_mentions[0]) if uncertain_mentions else ["patient_id"],
    )
    (output_dir / "summary.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--integrated-root", type=Path, default=DEFAULT_INTEGRATED_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(
        json.dumps(
            audit(args.root, args.integrated_root, args.output_dir),
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

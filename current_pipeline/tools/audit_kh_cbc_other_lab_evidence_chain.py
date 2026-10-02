"""Audit KH CBC/other-lab preservation and non-pathogen-specific guardrails."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from tools import build_deterministic_summary as summary_builder
from tools import evidence_preservation


DEFAULT_ROOT = Path("outputs/patient_info_KH_0728_2Days")
DEFAULT_INTEGRATED_ROOT = Path(
    "outputs/runs/2026-09-29_KH_integrated_decision_pipeline_v4_filmarray_gm_shadow"
)
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-30_KH_cbc_other_lab_evidence_chain_audit_v3"
)
AGENT_TAG = "evidence_v2_full33_20260917"


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
        match = re.search(r"patient_(\d+)", path.name)
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


def item_key(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").casefold())


def compatible_item(raw_item: str, agent_item: str) -> bool:
    left = item_key(raw_item)
    right = item_key(agent_item)
    if left == right:
        return True
    # Mojibake suffixes occur in the normalized CRP label; the biological
    # analyte remains explicit at the beginning of both strings.
    return bool(left and right and (left.startswith(right) or right.startswith(left)))


def host_tier_trace(agent: dict[str, Any]) -> dict[str, Any]:
    """Return the explicit S1 trace that supports the agent's H tier."""

    state = agent.get("host_state") or {}
    tier = str(state.get("immunocompromise_tier") or "Unknown")
    traces = [
        item
        for item in (agent.get("decision_trace") or [])
        if isinstance(item, dict) and item.get("step") == "S1_underlying_risk"
    ]
    evidence = [
        item
        for trace in traces
        for item in (trace.get("evidence") or [])
        if isinstance(item, dict)
    ]
    explicit_tier = any(
        str(item.get("field") or "").casefold() == "immunocompromise_tier"
        and tier in str(item.get("value") or "")
        for item in evidence
    )
    # Some valid H3 traces encode the tier in the rule id and keep the evidence
    # field focused on the exposure (for example P9 long-term immunosuppression).
    rule_tier = any(tier in str(trace.get("rule_id") or "") for trace in traces)
    substantive_evidence = [
        item
        for item in evidence
        if str(item.get("field") or "").casefold()
        not in {"immunocompromise_tier", "age"}
    ]
    return {
        "tier": tier,
        "trace_present": bool(traces),
        "tier_explicit_in_trace": explicit_tier or rule_tier,
        "substantive_evidence_present": bool(substantive_evidence),
        "rule_ids": "|".join(str(item.get("rule_id") or "") for item in traces),
        "evidence": substantive_evidence,
        "reasoning": [str(item) for item in (agent.get("reasoning") or [])],
    }


def parse_datetime(value: Any) -> datetime | None:
    text = str(value or "").strip()
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None


def write_csv(path: Path, records: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def agent_lab_comparison(
    raw_cbc: Any, raw_other: Any, agent: dict[str, Any]
) -> dict[str, Any]:
    raw_by_section = {"CBC": rows(raw_cbc), "OtherLab": rows(raw_other)}
    agent_rows = [
        item
        for item in ((agent.get("structured_raw_evidence") or {}).get("lab_observations") or [])
        if isinstance(item, dict)
    ]
    agent_by_section = {
        section: [
            item for item in agent_rows if str(item.get("source_section") or "") == section
        ]
        for section in raw_by_section
    }
    comparison: dict[str, Any] = {
        "agent_total_observation_count": len(agent_rows),
        "agent_cbc_observation_count": len(agent_by_section["CBC"]),
        "agent_other_lab_observation_count": len(agent_by_section["OtherLab"]),
        "agent_count_match": True,
        "agent_item_mismatch_count": 0,
        "agent_value_mismatch_count": 0,
        "agent_reported_time_mismatch_count": 0,
    }
    for section, raw_section in raw_by_section.items():
        agent_section = agent_by_section[section]
        if len(raw_section) != len(agent_section):
            comparison["agent_count_match"] = False
        for raw, normalized in zip(raw_section, agent_section):
            raw_item = str(raw.get("item") or raw.get("item_full") or raw.get("test") or "")
            if not compatible_item(raw_item, str(normalized.get("test_or_item") or "")):
                comparison["agent_item_mismatch_count"] += 1
            if str(raw.get("value") or "").strip() != str(normalized.get("value") or "").strip():
                comparison["agent_value_mismatch_count"] += 1
            if str(raw.get("reported_time") or "").strip() != str(
                normalized.get("reported_time") or ""
            ).strip():
                comparison["agent_reported_time_mismatch_count"] += 1
    return comparison


def promotion_guardrail_rows(integrated_root: Path) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    decision_dir = integrated_root / "decision_engine" / "patient_outputs"
    for path in sorted(
        decision_dir.glob("NGS_patient_*_test_aware_clinical_shadow.json"),
        key=patient_number,
    ):
        patient = patient_number(path)
        payload = read_json(path) or {}
        for candidate in payload.get("all_forwarded_candidates") or []:
            if str(candidate.get("clinical_decision")) != "picked_shadow":
                continue
            gate = candidate.get("promotion_gate") or {}
            axes = gate.get("axes") or {}
            timing = axes.get("direct_evidence_timing_profile") or {}
            analytical = candidate.get("analytical_profile") or {}
            route = str(gate.get("route") or "")
            independent_non_host_axis = bool(
                candidate.get("pre_promotion_decision") == "picked_shadow"
                or axes.get("direct_hospital_level_1_2")
                or timing.get("event_aligned_positive")
                or axes.get("exact_image_context")
                or analytical.get("reproducibility_axis")
                or int(analytical.get("selected_positive_test_count") or 0) > 0
            )
            output.append(
                {
                    "patient_id": patient,
                    "organism_name": candidate.get("organism_name"),
                    "candidate_source": candidate.get("candidate_source"),
                    "pre_promotion_decision": candidate.get("pre_promotion_decision"),
                    "promotion_route": route,
                    "host_support_flag": bool(axes.get("host_support")),
                    "verified_host_support_eligible": bool(
                        (axes.get("verified_host_support_profile") or {}).get("eligible")
                    ),
                    "selected_positive_test_count": analytical.get(
                        "selected_positive_test_count"
                    ),
                    "reproducibility_axis": analytical.get("reproducibility_axis"),
                    "direct_hospital_level_1_2": axes.get(
                        "direct_hospital_level_1_2"
                    ),
                    "event_aligned_positive": timing.get("event_aligned_positive"),
                    "exact_image_context": axes.get("exact_image_context"),
                    "independent_non_host_axis_present": independent_non_host_axis,
                    "host_only_pick": not independent_non_host_axis,
                }
            )
    return output


def audit(root: Path, integrated_root: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    patient_records: list[dict[str, Any]] = []
    qc_records: list[dict[str, Any]] = []
    tier_differences: list[dict[str, Any]] = []
    high_host_risk_records: list[dict[str, Any]] = []
    json_errors: list[dict[str, Any]] = []
    cbc_status_counts: Counter[str] = Counter()
    agent_a_counts: Counter[str] = Counter()
    deterministic_a_counts: Counter[str] = Counter()
    summary_h_counts: Counter[str] = Counter()
    summary_v_counts: Counter[str] = Counter()

    for patient_dir in patient_dirs(root):
        patient = patient_number(patient_dir)
        try:
            raw_cbc = read_json(source_path(patient_dir, "CBC"))
            raw_other = read_json(source_path(patient_dir, "other_lab"))
            agent = read_json(agent_path(patient_dir)) or {}
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            json_errors.append({"patient_id": patient, "error": str(exc)})
            continue
        raw_cbc_rows = rows(raw_cbc)
        raw_other_rows = rows(raw_other)
        observations = evidence_preservation.lab_observations_from_sources(
            raw_cbc, raw_other
        )
        acute = evidence_preservation.acute_instability_profile(raw_cbc, raw_other)
        cbc_status = evidence_preservation.cbc_data_status(raw_cbc)
        comparison = agent_lab_comparison(raw_cbc, raw_other, agent)
        summary = summary_builder.build_summary(
            patient_dir, agent_suffix_tag=AGENT_TAG
        )
        host = summary.get("host_context") or {}
        agent_state = agent.get("host_state") or {}
        host_trace = host_tier_trace(agent)
        agent_a = str(agent_state.get("acute_instability_tier") or "Unknown")
        deterministic_a = str(acute.get("acute_instability_tier") or "Unknown")
        cbc_status_counts[cbc_status["status"]] += 1
        agent_a_counts[agent_a] += 1
        deterministic_a_counts[deterministic_a] += 1
        summary_h_counts[str(host.get("immunocompromise_tier") or "Unknown")] += 1
        summary_v_counts[str(host.get("host_vulnerability_tier") or "Unknown")] += 1

        for observation in observations:
            if observation.get("qc_flags"):
                qc_records.append(
                    {
                        "patient_id": patient,
                        "observation_id": observation["observation_id"],
                        "source_section": observation["source_section"],
                        "source_record_index": observation["source_record_index"],
                        "test_or_item": observation["test_or_item"],
                        "value": observation["value"],
                        "unit": observation["unit"],
                        "reported_time": observation["reported_time"],
                        "qc_flags": "|".join(observation["qc_flags"]),
                        "eligible_for_acute_tiering": observation[
                            "eligible_for_acute_tiering"
                        ],
                    }
                )

        if agent_a != deterministic_a:
            tier_differences.append(
                {
                    "patient_id": patient,
                    "agent_acute_instability_tier": agent_a,
                    "deterministic_acute_instability_tier": deterministic_a,
                    "deterministic_trigger_count": acute["trigger_count"],
                    "deterministic_triggers": json.dumps(
                        acute["triggers"], ensure_ascii=False, separators=(",", ":")
                    ),
                    "interpretation": "review_required_no_automatic_tier_change",
                }
            )

        sample_time = evidence_preservation.mngs_sample_time_from_ranked(
            read_json(source_path(patient_dir, "mNGS_ranked_candidates"))
        )
        sample_dt = parse_datetime(sample_time)
        with_collected = sum(bool(item.get("collected_time")) for item in observations)
        within_report_window = 0
        outside_report_window = 0
        unknown_report_timing = 0
        for observation in observations:
            reported = parse_datetime(observation.get("reported_time"))
            if sample_dt is None or reported is None:
                unknown_report_timing += 1
            elif abs((reported - sample_dt).total_seconds()) <= 48 * 3600:
                within_report_window += 1
            else:
                outside_report_window += 1

        medication = host.get("risk_medications") or []
        patient_records.append(
            {
                "patient_id": patient,
                "cbc_raw_count": len(raw_cbc_rows),
                "other_lab_raw_count": len(raw_other_rows),
                "raw_total": len(raw_cbc_rows) + len(raw_other_rows),
                "deterministic_observation_count": len(observations),
                "source_complete": len(observations)
                == len(raw_cbc_rows) + len(raw_other_rows),
                **comparison,
                "cbc_status": cbc_status["status"],
                "collection_time_present_count": with_collected,
                "reported_within_48h_count": within_report_window,
                "reported_outside_48h_count": outside_report_window,
                "reported_timing_unknown_count": unknown_report_timing,
                "qc_flagged_observation_count": sum(
                    bool(item.get("qc_flags")) for item in observations
                ),
                "agent_H": agent_state.get("immunocompromise_tier"),
                "agent_A": agent_a,
                "deterministic_A": deterministic_a,
                "agent_V": agent_state.get("host_vulnerability_tier"),
                "summary_H": host.get("immunocompromise_tier"),
                "summary_A": host.get("acute_instability_tier"),
                "summary_V": host.get("host_vulnerability_tier"),
                "summary_opportunistic_host_support": host.get(
                    "opportunistic_host_support"
                ),
                "host_risk_evidence_strength": host.get(
                    "host_risk_evidence_strength"
                ),
                "risk_medications": "|".join(
                    str(item.get("medication_name") or "") for item in medication
                ),
                "host_tier_trace_present": host_trace["trace_present"],
                "host_tier_explicit_in_trace": host_trace["tier_explicit_in_trace"],
                "host_tier_substantive_evidence_present": host_trace[
                    "substantive_evidence_present"
                ],
            }
        )
        if str(host.get("immunocompromise_tier") or "") in {"H2", "H3"}:
            high_host_risk_records.append(
                {
                    "patient_id": patient,
                    "immunocompromise_tier": host.get("immunocompromise_tier"),
                    "rule_ids": host_trace["rule_ids"],
                    "trace_present": host_trace["trace_present"],
                    "tier_explicit_in_trace": host_trace["tier_explicit_in_trace"],
                    "substantive_evidence_present": host_trace[
                        "substantive_evidence_present"
                    ],
                    "supporting_evidence": json.dumps(
                        host_trace["evidence"],
                        ensure_ascii=False,
                        separators=(",", ":"),
                    ),
                    "reasoning": " | ".join(host_trace["reasoning"]),
                    "opportunistic_host_support": host.get(
                        "opportunistic_host_support"
                    ),
                }
            )

    promotion_rows = promotion_guardrail_rows(integrated_root)
    host_only_picks = [item for item in promotion_rows if item["host_only_pick"]]
    summary_opportunistic_valid = all(
        not item["summary_opportunistic_host_support"]
        or item["summary_H"] in {"H2", "H3"}
        or item["host_risk_evidence_strength"] == "strong"
        for item in patient_records
    )
    totals = {
        "patients": len(patient_records),
        "cbc_raw_rows": sum(int(item["cbc_raw_count"]) for item in patient_records),
        "other_lab_raw_rows": sum(
            int(item["other_lab_raw_count"]) for item in patient_records
        ),
        "raw_lab_rows": sum(int(item["raw_total"]) for item in patient_records),
        "deterministic_lab_observations": sum(
            int(item["deterministic_observation_count"]) for item in patient_records
        ),
        "agent_lab_observations": sum(
            int(item["agent_total_observation_count"]) for item in patient_records
        ),
        "patients_with_agent_count_mismatch": sum(
            not bool(item["agent_count_match"]) for item in patient_records
        ),
        "agent_item_mismatches": sum(
            int(item["agent_item_mismatch_count"]) for item in patient_records
        ),
        "agent_value_mismatches": sum(
            int(item["agent_value_mismatch_count"]) for item in patient_records
        ),
        "agent_reported_time_mismatches": sum(
            int(item["agent_reported_time_mismatch_count"])
            for item in patient_records
        ),
        "qc_flagged_observations": len(qc_records),
        "acute_tier_differences": len(tier_differences),
        "picked_candidates_audited": len(promotion_rows),
        "host_only_picks": len(host_only_picks),
        "high_host_risk_patients": len(high_host_risk_records),
    }
    assertions = {
        "all_33_patients_audited": len(patient_records) == 33,
        "no_json_read_errors": not json_errors,
        "deterministic_one_observation_per_raw_row": totals["raw_lab_rows"]
        == totals["deterministic_lab_observations"],
        "deterministic_source_complete_for_every_patient": all(
            item["source_complete"] for item in patient_records
        ),
        "summary_opportunistic_support_requires_immune_risk_not_acute_severity": summary_opportunistic_valid,
        "all_high_host_risk_tiers_have_explicit_s1_trace": all(
            item["trace_present"]
            and item["tier_explicit_in_trace"]
            and item["substantive_evidence_present"]
            for item in high_host_risk_records
        ),
        "no_picked_candidate_created_from_host_labs_alone": not host_only_picks,
        "cbc_other_lab_never_enters_hospital_organism_identity": all(
            str(item.get("candidate_source")) != "cbc_other_lab"
            for item in promotion_rows
        ),
    }
    report = {
        "schema_version": "kh_cbc_other_lab_evidence_chain_audit.v2",
        "answer_blind": True,
        "patient_root": str(root.resolve()),
        "integrated_root": str(integrated_root.resolve()),
        "totals": totals,
        "cbc_status_counts": dict(cbc_status_counts),
        "agent_acute_tier_counts": dict(agent_a_counts),
        "deterministic_acute_tier_counts": dict(deterministic_a_counts),
        "summary_immunocompromise_tier_counts": dict(summary_h_counts),
        "summary_vulnerability_tier_counts": dict(summary_v_counts),
        "assertions": assertions,
        "pass": all(assertions.values()),
        "review_items": {
            "agent_row_mapping_mismatches_require_review": bool(
                totals["patients_with_agent_count_mismatch"]
                or totals["agent_item_mismatches"]
                or totals["agent_value_mismatches"]
                or totals["agent_reported_time_mismatches"]
            ),
            "qc_flagged_normalized_values_require_source_review": bool(qc_records),
            "acute_tier_differences_are_shadow_only": bool(tier_differences),
        },
        "json_errors": json_errors,
    }

    write_csv(
        output_dir / "patient_completeness.csv",
        patient_records,
        list(patient_records[0]) if patient_records else ["patient_id"],
    )
    write_csv(
        output_dir / "lab_qc_flags.csv",
        qc_records,
        list(qc_records[0]) if qc_records else ["patient_id"],
    )
    write_csv(
        output_dir / "acute_tier_differences.csv",
        tier_differences,
        list(tier_differences[0]) if tier_differences else ["patient_id"],
    )
    write_csv(
        output_dir / "picked_host_guardrail_audit.csv",
        promotion_rows,
        list(promotion_rows[0]) if promotion_rows else ["patient_id"],
    )
    write_csv(
        output_dir / "high_host_risk_review.csv",
        high_host_risk_records,
        list(high_host_risk_records[0])
        if high_host_risk_records
        else ["patient_id"],
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

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from tools import deterministic_mngs_max_scorer as legacy
from tools import normalized_agent_fallback
from tools import pathogen_rule_lists


DEFAULT_POLICY = Path("rules/test_aware_deterministic_shadow_v2.json")
DECISION_ORDER = {
    "picked_shadow": 0,
    "review_high_priority": 1,
    "review_context_needed": 2,
    "review_low_specificity": 3,
}


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_int(value: Any, default: int = 999) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def technical_repeat_group_type(signals: list[dict[str, Any]]) -> str:
    """Classify selected branches without claiming biological independence."""
    if len(signals) < 2:
        return "single_test"
    conditions = [
        str(signal.get("protocol_condition") or "-").strip()
        for signal in signals
    ]
    starts = {
        str(signal.get("test_start_time") or "").strip()
        for signal in signals
        if str(signal.get("test_start_time") or "").strip()
    }
    if any(
        "原library重新上機" in value
        or "original library rerun" in value.casefold()
        for value in conditions
    ):
        return "explicit_original_library_rerun"
    if len(starts) > 1:
        return "multiple_run_times_unclassified"
    if any("12000" in value for value in conditions) and any(
        value in {"", "-"} for value in conditions
    ):
        return "same_run_parallel_12000g_and_unannotated"
    return "same_run_parallel_protocols"


def rank_band(rank: int, policy: dict[str, Any]) -> str:
    bands = policy["rank_bands"]
    if rank <= int(bands["top"]):
        return "top_1"
    if rank <= int(bands["strong"]):
        return "top_3"
    if rank <= int(bands["moderate"]):
        return "top_10"
    return "below_top_10"


def analytical_profile(entry: dict[str, Any], policy: dict[str, Any]) -> dict[str, Any]:
    signals = entry.get("selected_test_signals") or []
    if not signals:
        raise ValueError(f"Scorer entry has no selected test signal: {entry.get('organism_name')}")
    if any(float(signal.get("reads") or 0) <= 0 for signal in signals):
        raise ValueError(f"Zero/nonpositive reads reached scorer: {entry.get('organism_name')}")
    seq_ids = [str(signal.get("seq_id") or "") for signal in signals]
    if len(seq_ids) != len(set(seq_ids)):
        raise ValueError(f"Duplicate seq_id reached scorer: {entry.get('organism_name')}")
    best_rank = min(as_int(signal.get("rank_in_retained_universe")) for signal in signals)
    selected_count = len(signals)
    cross_molecule = bool(entry.get("cross_molecule_selected"))
    repeated_test = selected_count >= 2
    reproducible = cross_molecule or repeated_test
    return {
        "best_rank_in_retained_universe": best_rank,
        "rank_band": rank_band(best_rank, policy),
        "selected_positive_test_count": selected_count,
        "dna_selected": bool(entry.get("dna_selected")),
        "rna_selected": bool(entry.get("rna_selected")),
        "cross_molecule_selected": cross_molecule,
        "technical_repeat": repeated_test,
        "technical_repeat_group_type": technical_repeat_group_type(signals),
        "reproducibility_axis": reproducible,
        "reproducibility_note": (
            "DNA/RNA concordance or repeated seq_id support within one specimen; counted once and not treated as independent clinical confirmation."
            if reproducible else
            "Single-test analytical signal only."
        ),
        "normalization_available_test_count": sum(bool(signal.get("normalization_available")) for signal in signals),
        "fully_evaluable_test_count": sum(signal.get("qc_status") == "evaluable" for signal in signals),
        "per_test_signals": signals,
        "reads_sum_across_tests": None,
    }


def low_specificity_repeat_high_eligible(
    entry: dict[str, Any], analytical: dict[str, Any], policy: dict[str, Any]
) -> tuple[bool, dict[str, Any]]:
    """Evaluate the answer-blind Route A analytical High gate."""
    config = policy.get("low_specificity_repeat_high") or {}
    signals = analytical.get("per_test_signals") or []
    families = set(config.get("families") or [])
    accepted_qc = set(
        config.get("accepted_qc_statuses") or ["evaluable", "partial_evaluable"]
    )
    ranks = [as_int(signal.get("rank_in_retained_universe")) for signal in signals]
    rpms = [signal.get("rpm_total") for signal in signals]
    all_normalized = bool(signals) and all(value is not None for value in rpms)
    minimum_rpm = min(float(value) for value in rpms) if all_normalized else None
    group_type = str(analytical.get("technical_repeat_group_type") or "")
    axes = {
        "enabled": bool(config.get("enabled")),
        "taxonomy_family": entry.get("taxonomy_family"),
        "taxonomy_mapping_status": entry.get("taxonomy_mapping_status"),
        "specimen_context": entry.get("specimen_context"),
        "selected_branch_count": len(signals),
        "all_selected_branches_rank_one": bool(signals) and all(rank == 1 for rank in ranks),
        "all_selected_branches_normalized": all_normalized,
        "minimum_selected_branch_rpm": minimum_rpm,
        "technical_repeat_group_type": group_type,
        "all_selected_branches_qc_accepted": bool(signals) and all(
            str(signal.get("qc_status") or "") in accepted_qc for signal in signals
        ),
        "cross_molecule_selected": bool(
            analytical.get("cross_molecule_selected")
        ),
    }
    eligible = bool(
        axes["enabled"]
        and entry.get("taxonomy_family") in families
        and entry.get("taxonomy_mapping_status") == config.get(
            "required_mapping_status", "exact_species"
        )
        and entry.get("specimen_context") == config.get(
            "required_specimen_context", "lower_respiratory"
        )
        and len(signals) >= int(config.get("minimum_selected_branches", 2))
        and axes["all_selected_branches_rank_one"]
        and all_normalized
        and minimum_rpm is not None
        and minimum_rpm >= float(config.get("minimum_per_branch_rpm", 10.0))
        and axes["all_selected_branches_qc_accepted"]
        and (
            not config.get("require_no_cross_molecule", False)
            or not axes["cross_molecule_selected"]
        )
        and group_type not in set(config.get("excluded_repeat_group_types") or [])
    )
    axes["eligible"] = eligible
    return eligible, axes


def context_rescue_eligible(analytical: dict[str, Any], policy: dict[str, Any]) -> bool:
    """Keep analytically credible signals visible without treating them as pathogens."""
    rescue = policy["context_rescue"]
    best_rank = int(analytical["best_rank_in_retained_universe"])
    if best_rank <= int(rescue["single_test_rank_at_most"]):
        return True
    return bool(
        analytical["reproducibility_axis"]
        and best_rank <= int(rescue["reproducible_rank_at_most"])
    )


def forward_to_clinical_scorer(decision: str, policy: dict[str, Any]) -> bool:
    return decision in set(policy["clinical_scorer_forwarding"]["include_decisions"])


def _collect_observation_rows(value: Any, output: dict[str, dict[str, Any]]) -> None:
    if isinstance(value, dict):
        observation_id = str(value.get("observation_id") or "").strip()
        if observation_id:
            output[observation_id] = value
        for child in value.values():
            _collect_observation_rows(child, output)
    elif isinstance(value, list):
        for child in value:
            _collect_observation_rows(child, output)


def hospital_observation_index(final_summary: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Index source observations without making their availability a scoring gate."""
    output: dict[str, dict[str, Any]] = {}
    source_files = final_summary.get("source_files") or {}
    for module in ("culture", "molecular_microbiology", "filmarray_gmtest"):
        source = str(source_files.get(module) or "").strip()
        if not source:
            continue
        path = Path(source)
        if not path.exists():
            continue
        try:
            _collect_observation_rows(read_json(path), output)
        except (OSError, ValueError, json.JSONDecodeError):
            continue

    # Older normalized agent files can predate the source-complete observation
    # contract.  Reconstruct the same deterministic IDs from raw rows so
    # evidence_observation_ids remain resolvable without changing any tier.
    def raw_payload(key: str) -> Any:
        source = str(source_files.get(key) or "").strip()
        if not source:
            return None
        path = Path(source)
        if not path.exists():
            return None
        try:
            return json.loads(path.read_text(encoding="utf-8-sig"))
        except (OSError, ValueError, json.JSONDecodeError):
            return None

    raw_culture = raw_payload("raw_culture")
    if raw_culture is not None:
        _collect_observation_rows(
            normalized_agent_fallback.culture_agent_from_source(raw_culture),
            output,
        )
    raw_filmarray = raw_payload("raw_filmarray")
    raw_gm = raw_payload("raw_gm_test")
    if raw_filmarray is not None or raw_gm is not None:
        _collect_observation_rows(
            normalized_agent_fallback.filmarray_gm_agent_from_sources(
                raw_filmarray, raw_gm
            ),
            output,
        )
    raw_molecular = raw_payload("raw_molecular_microbiology")
    if raw_molecular is not None:
        _collect_observation_rows(
            normalized_agent_fallback.molecular_microbiology_agent_from_source(
                raw_molecular
            ),
            output,
        )
    return output


def _observation_ids(value: Any) -> list[str]:
    output: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            if key == "evidence_observation_ids" and isinstance(child, list):
                output.extend(str(item) for item in child if str(item).strip())
            else:
                output.extend(_observation_ids(child))
    elif isinstance(value, list):
        for child in value:
            output.extend(_observation_ids(child))
    return list(dict.fromkeys(output))


def enrich_hospital_evidence_detail(
    detail: dict[str, Any] | None,
    observation_index: dict[str, dict[str, Any]] | None,
) -> dict[str, Any] | None:
    if not isinstance(detail, dict):
        return None
    output = dict(detail)
    ids = _observation_ids(detail)
    index = observation_index or {}
    output["source_observations"] = [index[item] for item in ids if item in index]
    output["source_observation_ids"] = ids
    output["source_observation_resolution"] = {
        "requested": len(ids),
        "resolved": sum(item in index for item in ids),
        "missing": [item for item in ids if item not in index],
    }
    return output


def _exact_hospital_evidence_detail(
    evidence: dict[str, Any],
    observation_index: dict[str, dict[str, Any]] | None,
) -> dict[str, Any]:
    detail = {
        "organism_name": evidence.get("organism_name") or evidence.get("name"),
        "classification": evidence.get("classification"),
        "best_hospital_level": evidence.get("best_hospital_level"),
        "evidence_modules": evidence.get("evidence_modules") or [],
        "module_level_summary": evidence.get("module_level_summary") or {},
        "module_evidence": evidence.get("module_evidence") or {},
    }
    return enrich_hospital_evidence_detail(detail, observation_index) or detail


def _related_hospital_context(
    final_by_name: dict[str, dict[str, Any]],
    name: str,
    observation_index: dict[str, dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for item in legacy.related_representative_hospital_support(
        final_by_name, name
    ):
        output.append(
            enrich_hospital_evidence_detail(item, observation_index) or item
        )
    return output


def hospital_profile(
    name: str,
    final_by_name: dict[str, dict[str, Any]],
    observation_index: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    evidence = legacy.final_candidate_for_name(final_by_name, name)
    if not evidence:
        return {
            "exact_or_alias_evidence": False,
            "best_direct_level": None,
            "direct_support_modules_level_1_2": [],
            "direct_support_modules_level_1_3": [],
            "hospital_evidence_detail": None,
            "related_representative_context": _related_hospital_context(
                final_by_name, name, observation_index
            ),
        }
    level_1_2 = legacy.direct_hospital_support_modules(evidence, max_level_rank=2)
    level_1_3 = legacy.direct_hospital_support_modules(evidence, max_level_rank=3)
    levels = evidence.get("module_level_summary") if isinstance(evidence.get("module_level_summary"), dict) else {}
    direct_levels = [legacy.level_rank(levels.get(module)) for module in level_1_3]
    best_direct = min(direct_levels) if direct_levels else None
    return {
        "exact_or_alias_evidence": True,
        "matched_hospital_name": evidence.get("organism_name") or evidence.get("name"),
        "best_direct_level": f"Level {best_direct}" if best_direct is not None else None,
        "direct_support_modules_level_1_2": level_1_2,
        "direct_support_modules_level_1_3": level_1_3,
        "module_level_summary": levels,
        "hospital_evidence_detail": _exact_hospital_evidence_detail(evidence, observation_index),
        "related_representative_context": _related_hospital_context(
            final_by_name, name, observation_index
        ),
    }


def decide_entry(
    entry: dict[str, Any],
    *,
    final_by_name: dict[str, dict[str, Any]],
    host: dict[str, Any],
    policy: dict[str, Any],
    observation_index: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    analytical = analytical_profile(entry, policy)
    hospital = hospital_profile(
        str(entry.get("organism_name") or ""), final_by_name, observation_index
    )
    family = str(entry.get("taxonomy_family") or "unmapped_or_uncertain")
    specimen = str(entry.get("specimen_context") or "unknown_or_other")
    entry_tier = str(entry.get("screening_tier") or "low_colonizer")
    host_support = legacy.host_support(host)
    direct_l12 = bool(hospital["direct_support_modules_level_1_2"])
    direct_l13 = bool(hospital["direct_support_modules_level_1_3"])
    best_rank = analytical["best_rank_in_retained_universe"]
    top = best_rank <= int(policy["rank_bands"]["top"])
    strong = best_rank <= int(policy["rank_bands"]["strong"])
    reproducible = analytical["reproducibility_axis"]
    rescue_for_context = context_rescue_eligible(analytical, policy)
    lower_respiratory = specimen == "lower_respiratory"
    impossible = pathogen_rule_lists.impossible_infection_source_rule(entry.get("organism_name"))
    analytical_pick_families = set(policy["analytical_pick_families"])
    direct_pick_families = set(policy["direct_pick_families"])
    opportunistic_families = set(policy["opportunistic_families"])
    context_only_families = set(policy["context_only_families"])
    low_specificity_families = set(policy["low_specificity_families"])
    rules: list[str] = ["TA-S1-PER-TEST-SIGNALS", "TA-S1-NO-READ-SUM"]
    reasons: list[str] = []

    if reproducible:
        rules.append("TA-S2-REPRODUCIBILITY-ONE-AXIS")
        reasons.append("Signal is reproduced across selected tests and/or DNA/RNA, counted as one technical reproducibility axis.")
    if strong:
        rules.append("TA-S2-WITHIN-TEST-TOP3")
        reasons.append("At least one selected test places the organism within the retained-universe top three.")
    if hospital["related_representative_context"]:
        rules.append("TA-S3-RELATED-HOSPITAL-CONTEXT")
        reasons.append("A related representative/group has hospital-side evidence; this is context, not exact-organism confirmation.")
    if impossible:
        return {
            "decision": "review_low_specificity",
            "integrated_level": "Level 5",
            "formal_pick_allowed": False,
            "rule_ids": rules + ["TA-S4-IMPOSSIBLE-SOURCE"],
            "reasons": reasons + ["A deterministic impossible/nonpulmonary-source rule prevents formal pulmonary selection."],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if direct_l12 and family in direct_pick_families:
        level = hospital["best_direct_level"] or "Level 2"
        return {
            "decision": "picked_shadow",
            "integrated_level": level,
            "formal_pick_allowed": True,
            "rule_ids": rules + ["TA-S4-DIRECT-HOSPITAL-L1L2-PICK"],
            "reasons": reasons + ["Exact/approved-alias organism has Level 1/2 direct hospital evidence in a family eligible for direct selection."],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if direct_l12 and family in opportunistic_families:
        decision = "picked_shadow" if host_support else "review_high_priority"
        return {
            "decision": decision,
            "integrated_level": hospital["best_direct_level"] or "Level 2",
            "formal_pick_allowed": decision == "picked_shadow",
            "rule_ids": rules + ["TA-S4-OPPORTUNISTIC-DIRECT-HOST" if host_support else "TA-S4-OPPORTUNISTIC-DIRECT-NO-HOST"],
            "reasons": reasons + [
                "Direct opportunistic-pathogen evidence is accompanied by host support."
                if host_support else
                "Direct opportunistic-pathogen evidence lacks verified host support and remains high-priority review."
            ],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if direct_l12 and (family in context_only_families or family in low_specificity_families):
        return {
            "decision": "review_high_priority",
            "integrated_level": hospital["best_direct_level"] or "Level 2",
            "formal_pick_allowed": False,
            "rule_ids": rules + ["TA-S4-CONTEXT-FAMILY-DIRECT-REVIEW"],
            "reasons": reasons + ["Direct detection is important, but this colonization/reactivation/background-prone family requires clinical adjudication before Picked."],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if direct_l13 and family in analytical_pick_families and strong and reproducible and lower_respiratory:
        return {
            "decision": "picked_shadow",
            "integrated_level": "Level 3",
            "formal_pick_allowed": True,
            "rule_ids": rules + ["TA-S4-DIRECT-L3-ANALYTICAL-CONVERGENCE"],
            "reasons": reasons + ["Level 3 direct hospital context converges with reproducible top-three lower-respiratory mNGS evidence."],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if family in analytical_pick_families and top and reproducible and lower_respiratory and entry_tier in {"high", "medium"}:
        return {
            "decision": "picked_shadow",
            "integrated_level": "Level 3",
            "formal_pick_allowed": True,
            "rule_ids": rules + ["TA-S4-ANALYTICAL-TOP1-REPRODUCIBLE-PICK"],
            "reasons": reasons + ["Priority bacterial family is top-ranked in at least one lower-respiratory test and technically reproducible; no reads were combined."],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if family in analytical_pick_families and strong and lower_respiratory:
        return {
            "decision": "review_high_priority",
            "integrated_level": "Level 3",
            "formal_pick_allowed": False,
            "rule_ids": rules + ["TA-S4-PRIORITY-FAMILY-HIGH-REVIEW"],
            "reasons": reasons + ["Priority bacterial family has a top-three lower-respiratory signal but does not meet the conservative analytical-only Picked rule."],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if family == "other_respiratory_virus" and strong:
        return {
            "decision": "review_high_priority",
            "integrated_level": "Level 3" if reproducible else "Level 4",
            "formal_pick_allowed": False,
            "rule_ids": rules + ["TA-S4-RESP-VIRUS-REVIEW"],
            "reasons": reasons + ["Respiratory-virus mNGS signal remains high-priority review without matching direct hospital confirmation."],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if family in opportunistic_families:
        high = strong and reproducible and host_support
        context = strong or host_support or rescue_for_context
        return {
            "decision": "review_high_priority" if high else "review_context_needed" if context else "review_low_specificity",
            "integrated_level": "Level 3" if high else "Level 4" if context else "Level 5",
            "formal_pick_allowed": False,
            "rule_ids": rules + [
                "TA-S4-OPPORTUNISTIC-HOST-REVIEW"
                if high else
                "TA-S4-OPPORTUNISTIC-CONTEXT-RESCUE"
                if rescue_for_context else
                "TA-S4-OPPORTUNISTIC-CONTEXT"
            ],
            "reasons": reasons + ["Opportunistic candidates require host and/or pathogen-specific confirmation; analytical evidence alone cannot make them Picked."],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    repeat_high, repeat_high_axes = low_specificity_repeat_high_eligible(
        entry, analytical, policy
    )
    analytical["low_specificity_repeat_high_gate"] = repeat_high_axes
    if repeat_high:
        return {
            "decision": "review_high_priority",
            "integrated_level": "Level 3",
            "formal_pick_allowed": False,
            "rule_ids": rules + ["TA-S4-LOW-SPECIFICITY-REPEAT-HIGH"],
            "reasons": reasons + [
                "An exact-species low-specificity organism is rank 1 with normalized burden above the fixed threshold in every accepted technical branch; it enters High review but cannot become Picked from technical repetition alone."
            ],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if family in context_only_families:
        context = strong or reproducible or direct_l13 or bool(hospital["related_representative_context"])
        return {
            "decision": "review_context_needed" if context else "review_low_specificity",
            "integrated_level": "Level 4" if context else "Level 5",
            "formal_pick_allowed": False,
            "rule_ids": rules + ["TA-S4-CONTEXT-FAMILY-NO-AUTOPICK"],
            "reasons": reasons + ["Colonization/reactivation-prone family is retained for context but cannot become analytical-only Picked."],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    if rescue_for_context:
        return {
            "decision": "review_context_needed",
            "integrated_level": "Level 4",
            "formal_pick_allowed": False,
            "rule_ids": rules + ["TA-S4-ANALYTICAL-CONTEXT-RESCUE"],
            "reasons": reasons + [
                "A top-three single-test signal or reproducible top-ten signal is forwarded for clinical scoring; taxonomy limits automatic picking but does not erase the candidate."
            ],
            "analytical_profile": analytical,
            "hospital_profile": hospital,
            "host_support": host_support,
        }

    return {
        "decision": "review_low_specificity",
        "integrated_level": "Level 5",
        "formal_pick_allowed": False,
        "rule_ids": rules + ["TA-S4-LOW-SPECIFICITY"],
        "reasons": reasons + ["Low-specificity, nonpulmonary-prone, or uncertain family remains audit-visible but is not automatically promoted."],
        "analytical_profile": analytical,
        "hospital_profile": hospital,
        "host_support": host_support,
    }


def aggregate_patient_candidates(case_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in case_rows:
        grouped.setdefault(str(row["organism_key"]), []).append(row)
    output = []
    for organism_key, rows in grouped.items():
        rows.sort(key=lambda row: (DECISION_ORDER[row["decision"]], legacy.level_rank(row["integrated_level"])))
        representative = dict(rows[0])
        representative["case_decision_refs"] = [
            {
                "case_review_id": row["case_review_id"],
                "decision": row["decision"],
                "integrated_level": row["integrated_level"],
            }
            for row in rows
        ]
        representative["patient_case_count"] = len(rows)
        output.append(representative)
    output.sort(key=lambda row: (DECISION_ORDER[row["decision"]], legacy.level_rank(row["integrated_level"]), row["organism_name"]))
    return output


def hospital_only_shadow(
    final_summary: dict[str, Any],
    patient_candidates: list[dict[str, Any]],
    host: dict[str, Any],
    observation_index: dict[str, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    existing = [{"organism_name": item["organism_name"]} for item in patient_candidates]
    candidates = legacy.hospital_only_candidates(final_summary, existing, host)
    legacy.apply_formal_pick_guardrails(candidates, host)
    output = []
    for item in candidates:
        decision = "picked_shadow" if legacy.pickable(item, host) else (
            "review_high_priority" if legacy.level_rank(item.get("integrated_causative_level")) <= 3 else "review_context_needed"
        )
        output.append({
            "organism_name": item.get("organism_name"),
            "organism_key": legacy.mngs.normalize_organism_name(item.get("organism_name")),
            "decision": decision,
            "integrated_level": item.get("integrated_causative_level"),
            "formal_pick_allowed": decision == "picked_shadow",
            "forward_to_clinical_scorer": True,
            "evidence_source": item.get("evidence_source"),
            "rule_ids": item.get("applied_rules") or [],
            "reasons": item.get("integrated_reasoning") or [],
            "hospital_evidence_detail": enrich_hospital_evidence_detail(
                (item.get("key_evidence") or {}).get("hospital_evidence_detail"),
                observation_index,
            ),
        })
    return output


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def expected_scorer_entry_count(compact_summary: dict[str, Any]) -> int:
    value = (compact_summary.get("counts") or {}).get("scorer_entry")
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("Compact entry summary lacks a valid counts.scorer_entry")
    return value


def build_shadow(
    compact_entry_dir: Path,
    patient_root: Path,
    output_dir: Path,
    *,
    summary_suffix: str,
    summary_root: Path | None = None,
    policy_path: Path = DEFAULT_POLICY,
) -> dict[str, Any]:
    packet_paths = sorted((compact_entry_dir / "patient_packets").glob("NGS_patient_*_compact_entry_shadow.json"))
    if not packet_paths:
        raise FileNotFoundError("No compact entry patient packets found")
    compact_summary_path = compact_entry_dir / "summary.json"
    if not compact_summary_path.exists():
        raise FileNotFoundError(compact_summary_path)
    compact_summary = read_json(compact_summary_path)
    expected_entry_count = expected_scorer_entry_count(compact_summary)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Shadow output already exists: {output_dir}")
    policy = read_json(policy_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    patient_output_dir = output_dir / "patient_outputs"
    patient_output_dir.mkdir()
    case_csv_rows: list[dict[str, Any]] = []
    patient_csv_rows: list[dict[str, Any]] = []
    hospital_only_csv_rows: list[dict[str, Any]] = []
    forward_csv_rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()

    for packet_path in packet_paths:
        packet = read_json(packet_path)
        patient_id = str(packet["patient_id"])
        patient_dir = patient_root / f"NGS_patient_{patient_id}_json"
        if summary_root is None:
            summary_path = (
                patient_dir / "summary_outputs"
                / f"NGS_patient_{patient_id}_{summary_suffix}.json"
            )
        else:
            summary_path = summary_root / f"NGS_patient_{patient_id}_{summary_suffix}.json"
        if not summary_path.exists():
            raise FileNotFoundError(summary_path)
        final_summary = read_json(summary_path)
        final_by_name = legacy.final_candidate_map(final_summary)
        observation_index = hospital_observation_index(final_summary)
        host = legacy.host_context(final_summary)
        case_rows: list[dict[str, Any]] = []
        case_outputs = []
        for case in packet.get("cases") or []:
            decisions = []
            for entry in case.get("scorer_entries") or []:
                result = decide_entry(
                    entry,
                    final_by_name=final_by_name,
                    host=host,
                    policy=policy,
                    observation_index=observation_index,
                )
                result["forward_to_clinical_scorer"] = forward_to_clinical_scorer(result["decision"], policy)
                row = {
                    "patient_id": patient_id,
                    "case_review_id": case.get("case_review_id"),
                    "specimen_code": case.get("specimen_code"),
                    "specimen_site": case.get("specimen_site"),
                    "collected_time": case.get("collected_time"),
                    "organism_name": entry.get("organism_name"),
                    "organism_key": entry.get("organism_key"),
                    "category": entry.get("category"),
                    "taxonomy_family": entry.get("taxonomy_family"),
                    "screening_tier": entry.get("screening_tier"),
                    "specimen_context": entry.get("specimen_context"),
                    **result,
                    "source_entry": entry,
                }
                decisions.append(row)
                case_rows.append(row)
                counts[f"case_{row['decision']}"] += 1
                counts["case_scorer_entries"] += 1
                counts["selected_test_signals"] += len(result["analytical_profile"]["per_test_signals"])
            case_outputs.append({
                "case_review_id": case.get("case_review_id"),
                "specimen_code": case.get("specimen_code"),
                "specimen_site": case.get("specimen_site"),
                "collected_time": case.get("collected_time"),
                "candidate_decisions": decisions,
                "clinical_review_count_not_scored": len(case.get("clinical_review") or []),
                "qc_audit_count_not_scored": len(case.get("qc_audit_refs") or []),
            })
        patient_candidates = aggregate_patient_candidates(case_rows)
        hospital_only = hospital_only_shadow(
            final_summary, patient_candidates, host, observation_index
        )
        counts["hospital_only_candidates"] += len(hospital_only)
        for row in patient_candidates:
            counts[f"patient_{row['decision']}"] += 1
            counts["patient_forward_to_clinical_scorer"] += row["forward_to_clinical_scorer"]
            patient_csv_rows.append({
                "patient_id": patient_id,
                "organism_name": row["organism_name"],
                "organism_key": row["organism_key"],
                "taxonomy_family": row["taxonomy_family"],
                "screening_tier": row["screening_tier"],
                "decision": row["decision"],
                "integrated_level": row["integrated_level"],
                "formal_pick_allowed": row["formal_pick_allowed"],
                "forward_to_clinical_scorer": row["forward_to_clinical_scorer"],
                "best_rank_in_retained_universe": row["analytical_profile"]["best_rank_in_retained_universe"],
                "selected_positive_test_count": row["analytical_profile"]["selected_positive_test_count"],
                "cross_molecule_selected": row["analytical_profile"]["cross_molecule_selected"],
                "direct_hospital_level": row["hospital_profile"]["best_direct_level"],
                "direct_hospital_modules": "|".join(row["hospital_profile"]["direct_support_modules_level_1_3"]),
                "rule_ids": "|".join(row["rule_ids"]),
            })
            if row["forward_to_clinical_scorer"]:
                forward_csv_rows.append({
                    "patient_id": patient_id,
                    "organism_name": row["organism_name"],
                    "organism_key": row["organism_key"],
                    "evidence_source": "multi_assay_mngs",
                    "decision": row["decision"],
                    "integrated_level": row["integrated_level"],
                    "taxonomy_family": row["taxonomy_family"],
                    "screening_tier": row["screening_tier"],
                    "best_rank_in_retained_universe": row["analytical_profile"]["best_rank_in_retained_universe"],
                    "selected_positive_test_count": row["analytical_profile"]["selected_positive_test_count"],
                    "rule_ids": "|".join(row["rule_ids"]),
                })
        for row in hospital_only:
            counts[f"hospital_only_{row['decision']}"] += 1
            counts["hospital_only_forward_to_clinical_scorer"] += row["forward_to_clinical_scorer"]
            hospital_only_csv_rows.append({"patient_id": patient_id, **row, "rule_ids": "|".join(row["rule_ids"])})
            forward_csv_rows.append({
                "patient_id": patient_id,
                "organism_name": row["organism_name"],
                "organism_key": row["organism_key"],
                "evidence_source": row["evidence_source"],
                "decision": row["decision"],
                "integrated_level": row["integrated_level"],
                "taxonomy_family": "hospital_only",
                "screening_tier": "hospital_only",
                "best_rank_in_retained_universe": "",
                "selected_positive_test_count": "",
                "rule_ids": "|".join(row["rule_ids"]),
            })
        for row in case_rows:
            case_csv_rows.append({
                "patient_id": patient_id,
                "case_review_id": row["case_review_id"],
                "organism_name": row["organism_name"],
                "organism_key": row["organism_key"],
                "taxonomy_family": row["taxonomy_family"],
                "screening_tier": row["screening_tier"],
                "decision": row["decision"],
                "integrated_level": row["integrated_level"],
                "formal_pick_allowed": row["formal_pick_allowed"],
                "forward_to_clinical_scorer": row["forward_to_clinical_scorer"],
                "best_rank_in_retained_universe": row["analytical_profile"]["best_rank_in_retained_universe"],
                "selected_positive_test_count": row["analytical_profile"]["selected_positive_test_count"],
                "cross_molecule_selected": row["analytical_profile"]["cross_molecule_selected"],
                "direct_hospital_level": row["hospital_profile"]["best_direct_level"],
                "direct_hospital_modules": "|".join(row["hospital_profile"]["direct_support_modules_level_1_3"]),
                "rule_ids": "|".join(row["rule_ids"]),
            })
        patient_output = {
            "schema_version": "test_aware_deterministic_shadow.v2",
            "policy_id": policy["policy_id"],
            "answer_blind": True,
            "patient_id": patient_id,
            "host_context": host,
            "cases": case_outputs,
            "patient_organism_decisions": patient_candidates,
            "hospital_only_decisions": hospital_only,
            "picked_shadow": [row for row in patient_candidates + hospital_only if row["decision"] == "picked_shadow"],
            "review_high_priority": [row for row in patient_candidates + hospital_only if row["decision"] == "review_high_priority"],
            "review_context_needed": [row for row in patient_candidates + hospital_only if row["decision"] == "review_context_needed"],
            "review_low_specificity": [row for row in patient_candidates if row["decision"] == "review_low_specificity"],
            "clinical_scorer_forward": [
                row for row in patient_candidates + hospital_only
                if row.get("forward_to_clinical_scorer", row["decision"] != "review_low_specificity")
            ],
            "source_files": {
                "compact_entry_packet": str(packet_path.resolve()),
                "compact_entry_packet_sha256": sha256_file(packet_path),
                "hospital_summary": str(summary_path.resolve()),
                "hospital_summary_sha256": sha256_file(summary_path),
                "policy": str(policy_path.resolve()),
            },
            "constraints": policy["constraints"],
        }
        write_json(patient_output_dir / f"NGS_patient_{patient_id}_test_aware_deterministic_shadow.json", patient_output)

    if counts["case_scorer_entries"] != expected_entry_count:
        raise ValueError(
            f"Expected {expected_entry_count} compact scorer entries from input summary, "
            f"got {counts['case_scorer_entries']}"
        )
    if any(float(signal.get("reads") or 0) <= 0 for path in patient_output_dir.glob("*.json") for case in read_json(path)["cases"] for row in case["candidate_decisions"] for signal in row["analytical_profile"]["per_test_signals"]):
        raise ValueError("Zero-read signal reached output")
    write_csv(output_dir / "case_candidate_decisions.csv", case_csv_rows, list(case_csv_rows[0]))
    write_csv(output_dir / "patient_organism_decisions.csv", patient_csv_rows, list(patient_csv_rows[0]))
    hospital_fields = ["patient_id", "organism_name", "organism_key", "decision", "integrated_level", "formal_pick_allowed", "forward_to_clinical_scorer", "evidence_source", "rule_ids"]
    write_csv(output_dir / "hospital_only_decisions.csv", hospital_only_csv_rows, hospital_fields)
    write_csv(output_dir / "clinical_scorer_forward.csv", forward_csv_rows, list(forward_csv_rows[0]))
    summary = {
        "schema_version": "test_aware_deterministic_shadow_summary.v2",
        "policy_id": policy["policy_id"],
        "answer_blind": True,
        "answer_keys_used": [],
        "patient_count": len(packet_paths),
        "counts": dict(counts),
        "expected_case_scorer_entries_from_input": expected_entry_count,
        "source_compact_entry_dir": str(compact_entry_dir.resolve()),
        "source_compact_entry_summary_sha256": sha256_file(compact_entry_dir / "summary.json"),
        "policy_path": str(policy_path.resolve()),
        "policy_sha256": sha256_file(policy_path),
        "read_handling": "Per-test signals only; reads_sum_across_tests is always null and is never used for a decision.",
        "reproducibility_handling": "DNA/RNA concordance and multiple seq_id tests from one specimen form one technical reproducibility axis.",
        "forwarding_handling": "Candidate forwarding is separate from Picked; top-three single-test or reproducible top-ten signals remain available to the next clinical deterministic scorer.",
        "phenotype_handling": "Not used in deterministic scoring while patient/event linkage remains provisional.",
        "forced_pick_policy": "No patient is forced to have a Picked candidate.",
        "downstream_status": "Shadow only; no production patient summary or Picked output was overwritten.",
    }
    write_json(output_dir / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build an answer-blind test-aware deterministic mNGS shadow scorer")
    parser.add_argument("compact_entry_dir", type=Path)
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--summary-suffix", required=True)
    parser.add_argument(
        "--summary-root",
        type=Path,
        help="Optional flat directory containing the selected patient summaries.",
    )
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    args = parser.parse_args()
    print(json.dumps(build_shadow(
        args.compact_entry_dir,
        args.patient_root,
        args.output_dir,
        summary_suffix=args.summary_suffix,
        summary_root=args.summary_root,
        policy_path=args.policy,
    ), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

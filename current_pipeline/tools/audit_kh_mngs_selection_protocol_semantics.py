"""Audit KH mNGS source-selection tags and protocol/repeat semantics.

The audit is answer blind.  It does not read physician benchmark answers and
does not alter any frozen scorer output.  It separates three questions:

1. Which positive observations pass the supplied Code_NTC / Code_RK_NTC gate?
2. Does the test-aware scorer use those codes or protocol labels as weights?
3. What would change *inside the scorer stage only* if same-molecule technical
   repeats stopped contributing to its one reproducibility axis?
"""

from __future__ import annotations

import argparse
import copy
import csv
import inspect
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable, Iterable

from tools import build_test_aware_deterministic_shadow_scorer as scorer


DEFAULT_INVENTORY = Path(
    "outputs/runs/2026-09-23_KH_multi_assay_mngs_inventory_source_contract_v3"
)
DEFAULT_COMPACT = Path(
    "outputs/runs/2026-09-25_KH_compact_multi_assay_entry_taxonomy_auto_v1"
)
DEFAULT_SCORER = Path(
    "outputs/runs/2026-09-25_KH_test_aware_deterministic_taxonomy_auto_v1"
)
DEFAULT_PATIENT_ROOT = Path("outputs/patient_info_KH_0728_2Days")
DEFAULT_POLICY = Path("rules/test_aware_deterministic_shadow_v2.json")
DEFAULT_OUTPUT = Path(
    "outputs/runs/2026-09-30_KH_mngs_selection_protocol_semantics_audit_v1"
)
DEFAULT_SUMMARY_SUFFIX = "final_summary_selected_dna_20260918"


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: Any) -> float | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def as_int(value: Any) -> int | None:
    number = as_float(value)
    return int(number) if number is not None else None


def as_bool(value: Any) -> bool:
    return str(value or "").strip().lower() == "true"


def parse_tags(value: Any) -> set[str]:
    if isinstance(value, list):
        return {str(item).strip() for item in value if str(item).strip()}
    text = str(value or "").strip()
    if not text:
        return set()
    parsed = json.loads(text)
    if not isinstance(parsed, list):
        raise ValueError(f"selection_tags must be a list: {value!r}")
    return {str(item).strip() for item in parsed if str(item).strip()}


def normalized_tags(value: Any) -> set[str]:
    output = set()
    for tag in parse_tags(value):
        tag = tag.replace("Code_RK_NTC=RK00", "Code_RK_NTC=RK0")
        tag = tag.replace("Code_RK_NTC=RK08", "Code_RK_NTC=RK8")
        output.add(tag)
    return output


def tag_combo(value: Any) -> str:
    tags = normalized_tags(value)
    ordered = [
        tag
        for tag in ("Code_NTC=00", "Code_RK_NTC=RK0", "Code_RK_NTC=RK8")
        if tag in tags
    ]
    extra = sorted(tags - set(ordered))
    return "+".join(ordered + extra) if ordered or extra else "none"


def expected_priority(value: Any) -> tuple[int | None, str]:
    tags = normalized_tags(value)
    ntc = "Code_NTC=00" in tags
    rk0 = "Code_RK_NTC=RK0" in tags
    rk8 = "Code_RK_NTC=RK8" in tags
    if ntc and rk0:
        return 1, "NTC00+RK00"
    if ntc and rk8:
        return 2, "NTC00+RK8"
    if rk0:
        return 3, "RK00"
    if rk8:
        return 4, "RK8"
    if ntc:
        return 5, "NTC00"
    return None, "unranked"


def observation_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(row.get("patient_id") or ""),
        str(row.get("case_review_id") or ""),
        str(row.get("organism_key") or ""),
        str(row.get("seq_id") or ""),
    )


def candidate_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return observation_key(row)[:3]


def is_positive(row: dict[str, Any]) -> bool:
    return (
        row.get("availability_status") == "evaluable"
        and (as_float(row.get("reads")) or 0.0) > 0
    )


def condition_flags(condition: str) -> list[str]:
    text = str(condition or "").strip()
    flags = []
    if not text or text == "-":
        return ["unannotated"]
    if "12000" in text:
        flags.append("centrifuge_12000g_10min")
    if "2倍稀釋" in text:
        flags.append("adapter_2x_dilution")
    if "稀釋一半" in text:
        flags.append("adapter_half_wording")
    if "10倍稀釋" in text:
        flags.append("adapter_10x_dilution")
    if "片段化" in text:
        flags.append("fragmentation_94c_5min")
    if "原library重新上機" in text:
        flags.append("original_library_rerun")
    if "Input 500ng去建庫" in text:
        flags.append("input_500ng_library_pcr10")
    if "考核" in text or "教育訓練" in text:
        flags.append("operator_training_or_assessment")
    return flags or ["other_annotated_protocol"]


def repeat_group_type(test_rows: list[dict[str, Any]]) -> str:
    conditions = [str(row.get("condition") or "-") for row in test_rows]
    starts = {str(row.get("start_time") or "") for row in test_rows}
    if any("原library重新上機" in value for value in conditions):
        return "explicit_original_library_rerun"
    if len(starts) > 1:
        return "multiple_run_times_unclassified"
    if any("12000" in value for value in conditions) and any(
        value == "-" for value in conditions
    ):
        return "same_run_parallel_12000g_and_unannotated"
    return "same_run_parallel_protocols"


def compact_entries(compact_dir: Path) -> list[dict[str, Any]]:
    output = []
    for path in sorted((compact_dir / "patient_packets").glob("*.json")):
        packet = read_json(path)
        patient_id = str(packet["patient_id"])
        for case in packet.get("cases") or []:
            case_id = str(case.get("case_review_id") or "")
            for entry in case.get("scorer_entries") or []:
                output.append(
                    {
                        **entry,
                        "patient_id": patient_id,
                        "case_review_id": case_id,
                    }
                )
    return output


def _policy_current(tags: set[str]) -> bool:
    return bool(tags & {"Code_NTC=00", "Code_RK_NTC=RK0", "Code_RK_NTC=RK8"})


def _policy_exclude_rk8_only(tags: set[str]) -> bool:
    return _policy_current(tags) and tags != {"Code_RK_NTC=RK8"}


SELECTION_POLICIES: list[tuple[str, str, Callable[[set[str]], bool]]] = [
    (
        "current_source_or",
        "Code_NTC=00 OR Code_RK_NTC in {RK00,RK8}",
        _policy_current,
    ),
    (
        "exclude_RK8_only",
        "Current source gate except observations carrying RK8 alone",
        _policy_exclude_rk8_only,
    ),
    ("NTC00_any", "Require Code_NTC=00", lambda tags: "Code_NTC=00" in tags),
    (
        "RK00_any",
        "Require Code_RK_NTC=RK00/RK0",
        lambda tags: "Code_RK_NTC=RK0" in tags,
    ),
    (
        "NTC00_and_RK00",
        "Require both Code_NTC=00 and Code_RK_NTC=RK00/RK0",
        lambda tags: {"Code_NTC=00", "Code_RK_NTC=RK0"} <= tags,
    ),
]


def audit(
    inventory_dir: Path,
    compact_dir: Path,
    scorer_dir: Path,
    patient_root: Path,
    policy_path: Path,
    output_dir: Path,
    *,
    summary_suffix: str = DEFAULT_SUMMARY_SUFFIX,
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    test_rows = read_csv(inventory_dir / "test_qc.csv")
    observation_rows = read_csv(inventory_dir / "observation_qc.csv")
    positive_rows = [row for row in observation_rows if is_positive(row)]
    entries = compact_entries(compact_dir)
    frozen_decisions = read_csv(scorer_dir / "case_candidate_decisions.csv")
    policy = read_json(policy_path)

    observations = {observation_key(row): row for row in observation_rows}
    tests = {
        (str(row["patient_id"]), str(row["case_review_id"]), str(row["seq_id"])): row
        for row in test_rows
    }
    frozen_index = {candidate_key(row): row for row in frozen_decisions}

    selected_formula_mismatches = []
    priority_mismatches = []
    for row in positive_rows:
        priority, rule = expected_priority(row.get("selection_tags"))
        wanted_selected = priority is not None
        if as_bool(row.get("selected")) != wanted_selected:
            selected_formula_mismatches.append(observation_key(row))
        if as_int(row.get("analytical_rank_priority")) != priority:
            priority_mismatches.append((*observation_key(row), "priority"))
        observed_rule = str(row.get("analytical_rank_rule") or "")
        if observed_rule != rule:
            priority_mismatches.append((*observation_key(row), "rule"))

    scorer_signal_keys: set[tuple[str, str, str, str]] = set()
    signals_by_entry: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    scorer_signal_mismatches = []
    for entry in entries:
        entry_key = candidate_key(entry)
        matched = []
        for signal in entry.get("selected_test_signals") or []:
            signal_key = (*entry_key, str(signal.get("seq_id") or ""))
            row = observations.get(signal_key)
            if row is None or not is_positive(row) or not as_bool(row.get("selected")):
                scorer_signal_mismatches.append(signal_key)
                continue
            matched.append(row)
            scorer_signal_keys.add(signal_key)
        signals_by_entry[entry_key] = matched

    selection_profile = []
    profile_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in positive_rows:
        profile_groups[tag_combo(row.get("selection_tags"))].append(row)
    for combo, rows in sorted(profile_groups.items()):
        selection_profile.append(
            {
                "tag_combination": combo,
                "positive_observation_count": len(rows),
                "selected_observation_count": sum(as_bool(row.get("selected")) for row in rows),
                "filtered_observation_count": sum(not as_bool(row.get("selected")) for row in rows),
                "unique_test_count": len(
                    {
                        (row["patient_id"], row["case_review_id"], row["seq_id"])
                        for row in rows
                    }
                ),
                "unique_case_organism_count": len({candidate_key(row) for row in rows}),
                "scorer_signal_count": sum(
                    observation_key(row) in scorer_signal_keys for row in rows
                ),
            }
        )

    selection_sensitivity = []
    for policy_name, definition, selector in SELECTION_POLICIES:
        kept_signal_count = 0
        kept_candidate_count = 0
        lost_candidates = []
        for entry_key, rows in signals_by_entry.items():
            kept = [row for row in rows if selector(normalized_tags(row["selection_tags"]))]
            kept_signal_count += len(kept)
            if kept:
                kept_candidate_count += 1
            else:
                lost_candidates.append(entry_key)
        selection_sensitivity.append(
            {
                "counterfactual_policy": policy_name,
                "definition": definition,
                "kept_scorer_signal_count": kept_signal_count,
                "scorer_entries_with_at_least_one_signal": kept_candidate_count,
                "scorer_entries_losing_all_signals": len(lost_candidates),
                "lost_entry_examples": json.dumps(
                    ["|".join(value) for value in lost_candidates[:10]],
                    ensure_ascii=False,
                ),
                "interpretation": "eligibility sensitivity only; clinical performance was not evaluated",
            }
        )

    condition_groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in test_rows:
        condition_groups[str(row.get("condition") or "-")].append(row)
    protocol_profile = []
    for condition, rows in sorted(condition_groups.items()):
        protocol_profile.append(
            {
                "condition": condition,
                "flags": "|".join(condition_flags(condition)),
                "test_count": len(rows),
                "patient_count": len({row["patient_id"] for row in rows}),
                "case_molecule_count": len(
                    {
                        (row["patient_id"], row["case_review_id"], row["nucleic_type"])
                        for row in rows
                    }
                ),
                "dna_test_count": sum(row.get("nucleic_type") == "DNA" for row in rows),
                "rna_test_count": sum(row.get("nucleic_type") == "RNA" for row in rows),
                "evaluable_count": sum(row.get("qc_status") == "evaluable" for row in rows),
                "partial_evaluable_count": sum(
                    row.get("qc_status") == "partial_evaluable" for row in rows
                ),
                "raw_unavailable_count": sum(
                    row.get("qc_status") == "raw_unavailable" for row in rows
                ),
            }
        )

    repeat_groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in test_rows:
        repeat_groups[
            (str(row["patient_id"]), str(row["case_review_id"]), str(row["nucleic_type"]))
        ].append(row)
    repeat_profile = []
    for group_key, rows in sorted(repeat_groups.items()):
        if len(rows) < 2:
            continue
        seq_ids = {str(row["seq_id"]) for row in rows}
        group_observations = [
            row
            for row in positive_rows
            if str(row["patient_id"]) == group_key[0]
            and str(row["case_review_id"]) == group_key[1]
            and str(row["nucleic_type"]) == group_key[2]
            and str(row["seq_id"]) in seq_ids
        ]
        repeat_profile.append(
            {
                "patient_id": group_key[0],
                "case_review_id": group_key[1],
                "nucleic_type": group_key[2],
                "repeat_group_type": repeat_group_type(rows),
                "test_count": len(rows),
                "seq_ids": "|".join(str(row["seq_id"]) for row in rows),
                "start_times": "|".join(sorted({str(row.get("start_time") or "") for row in rows})),
                "conditions": " || ".join(str(row.get("condition") or "-") for row in rows),
                "qc_statuses": "|".join(str(row.get("qc_status") or "") for row in rows),
                "positive_observation_count": len(group_observations),
                "selected_positive_count": sum(
                    as_bool(row.get("selected")) for row in group_observations
                ),
                "independence_status": "requires_laboratory_definition",
            }
        )

    replay_mismatches = []
    repeat_counterfactual = []
    for packet_path in sorted((compact_dir / "patient_packets").glob("*.json")):
        packet = read_json(packet_path)
        patient_id = str(packet["patient_id"])
        summary_path = (
            patient_root
            / f"NGS_patient_{patient_id}_json"
            / "summary_outputs"
            / f"NGS_patient_{patient_id}_{summary_suffix}.json"
        )
        final_summary = scorer.read_json(summary_path)
        final_by_name = scorer.legacy.final_candidate_map(final_summary)
        observation_index = scorer.hospital_observation_index(final_summary)
        host = scorer.legacy.host_context(final_summary)
        for case in packet.get("cases") or []:
            case_id = str(case.get("case_review_id") or "")
            for entry in case.get("scorer_entries") or []:
                entry_key = (patient_id, case_id, str(entry.get("organism_key") or ""))
                current = scorer.decide_entry(
                    entry,
                    final_by_name=final_by_name,
                    host=host,
                    policy=policy,
                    observation_index=observation_index,
                )
                frozen = frozen_index.get(entry_key)
                if frozen is None or (
                    current["decision"], current["integrated_level"]
                ) != (frozen.get("decision"), frozen.get("integrated_level")):
                    replay_mismatches.append(entry_key)

                signals = entry.get("selected_test_signals") or []
                if len(signals) < 2 or bool(entry.get("cross_molecule_selected")):
                    continue
                reduced = copy.deepcopy(entry)
                reduced["selected_test_signals"] = [
                    sorted(
                        signals,
                        key=lambda item: (
                            as_int(item.get("rank_in_retained_universe")) or 10**9,
                            -(as_float(item.get("reads")) or 0.0),
                            str(item.get("seq_id") or ""),
                        ),
                    )[0]
                ]
                reduced["selected_positive_test_count"] = 1
                reduced["technical_repeat"] = False
                counterfactual = scorer.decide_entry(
                    reduced,
                    final_by_name=final_by_name,
                    host=host,
                    policy=policy,
                    observation_index=observation_index,
                )
                signal_tests = [
                    tests[(patient_id, case_id, str(signal.get("seq_id") or ""))]
                    for signal in signals
                ]
                changed = (
                    current["decision"], current["integrated_level"]
                ) != (
                    counterfactual["decision"],
                    counterfactual["integrated_level"],
                )
                repeat_counterfactual.append(
                    {
                        "patient_id": patient_id,
                        "case_review_id": case_id,
                        "organism_name": entry.get("organism_name"),
                        "organism_key": entry.get("organism_key"),
                        "repeat_group_type": repeat_group_type(signal_tests),
                        "selected_test_count": len(signals),
                        "best_rank_in_retained_universe": current["analytical_profile"][
                            "best_rank_in_retained_universe"
                        ],
                        "current_decision": current["decision"],
                        "current_level": current["integrated_level"],
                        "no_technical_repeat_decision": counterfactual["decision"],
                        "no_technical_repeat_level": counterfactual["integrated_level"],
                        "decision_or_level_changed": changed,
                        "current_rule_ids": "|".join(current["rule_ids"]),
                        "counterfactual_rule_ids": "|".join(counterfactual["rule_ids"]),
                        "scope_note": "scorer-stage-only; candidate entry was held fixed",
                    }
                )

    changed_repeat_rows = [
        row for row in repeat_counterfactual if row["decision_or_level_changed"]
    ]
    analytical_source = inspect.getsource(scorer.analytical_profile)
    decision_source = inspect.getsource(scorer.decide_entry)
    priority_not_consumed = "analytical_rank_priority" not in analytical_source + decision_source
    protocol_not_consumed = "protocol_condition" not in analytical_source + decision_source

    write_csv(
        output_dir / "selection_tag_profile.csv",
        selection_profile,
        [
            "tag_combination",
            "positive_observation_count",
            "selected_observation_count",
            "filtered_observation_count",
            "unique_test_count",
            "unique_case_organism_count",
            "scorer_signal_count",
        ],
    )
    write_csv(
        output_dir / "selection_gate_sensitivity.csv",
        selection_sensitivity,
        list(selection_sensitivity[0]),
    )
    write_csv(
        output_dir / "protocol_condition_profile.csv",
        protocol_profile,
        list(protocol_profile[0]),
    )
    write_csv(
        output_dir / "multi_test_group_profile.csv",
        repeat_profile,
        list(repeat_profile[0]),
    )
    write_csv(
        output_dir / "technical_repeat_scorer_counterfactual.csv",
        repeat_counterfactual,
        list(repeat_counterfactual[0]),
    )

    selection_sensitivity_map = {
        row["counterfactual_policy"]: row for row in selection_sensitivity
    }
    assertions = {
        "answer_blind": True,
        "positive_observation_count_is_1898": len(positive_rows) == 1898,
        "source_selected_count_is_1267": sum(
            as_bool(row.get("selected")) for row in positive_rows
        )
        == 1267,
        "source_filtered_count_is_631": sum(
            not as_bool(row.get("selected")) for row in positive_rows
        )
        == 631,
        "selection_formula_matches_all_positive_observations": not selected_formula_mismatches,
        "priority_mapping_matches_all_positive_observations": not priority_mismatches,
        "all_856_scorer_signals_map_to_selected_positive_observations": (
            len(scorer_signal_keys) == 856 and not scorer_signal_mismatches
        ),
        "frozen_scorer_replay_matches_428_entries": (
            len(entries) == 428 and not replay_mismatches
        ),
        "analytical_rank_priority_is_not_consumed_as_scorer_weight": priority_not_consumed,
        "protocol_condition_is_not_consumed_as_scorer_weight": protocol_not_consumed,
        "technical_repeat_only_entry_count_is_48": len(repeat_counterfactual) == 48,
        "technical_repeat_counterfactual_changes_28_routes": len(changed_repeat_rows) == 28,
        "technical_repeat_counterfactual_changes_no_picked_route": not any(
            "picked_shadow"
            in {row["current_decision"], row["no_technical_repeat_decision"]}
            and row["current_decision"] != row["no_technical_repeat_decision"]
            for row in changed_repeat_rows
        ),
        "exclude_RK8_only_loses_26_scorer_entries": (
            selection_sensitivity_map["exclude_RK8_only"][
                "scorer_entries_losing_all_signals"
            ]
            == 26
        ),
    }
    result = {
        "schema_version": "kh_mngs_selection_protocol_semantics_audit.v1",
        "answer_blind": True,
        "inputs": {
            "inventory_dir": str(inventory_dir.resolve()),
            "compact_dir": str(compact_dir.resolve()),
            "scorer_dir": str(scorer_dir.resolve()),
            "patient_root": str(patient_root.resolve()),
            "policy_path": str(policy_path.resolve()),
            "summary_suffix": summary_suffix,
        },
        "counts": {
            "tests": len(test_rows),
            "positive_observations": len(positive_rows),
            "selected_positive_observations": sum(
                as_bool(row.get("selected")) for row in positive_rows
            ),
            "filtered_positive_observations": sum(
                not as_bool(row.get("selected")) for row in positive_rows
            ),
            "scorer_entries": len(entries),
            "scorer_signals": len(scorer_signal_keys),
            "distinct_protocol_conditions": len(condition_groups),
            "annotated_protocol_tests": sum(
                str(row.get("condition") or "-") != "-" for row in test_rows
            ),
            "multi_test_case_molecule_groups": len(repeat_profile),
            "multi_test_group_types": dict(
                sorted(Counter(row["repeat_group_type"] for row in repeat_profile).items())
            ),
            "technical_repeat_only_scorer_entries": len(repeat_counterfactual),
            "technical_repeat_scorer_route_changes": len(changed_repeat_rows),
            "technical_repeat_route_change_types": dict(
                sorted(
                    Counter(
                        f"{row['current_decision']}->{row['no_technical_repeat_decision']}"
                        for row in changed_repeat_rows
                    ).items()
                )
            ),
            "technical_repeat_changed_group_types": dict(
                sorted(Counter(row["repeat_group_type"] for row in changed_repeat_rows).items())
            ),
        },
        "selection_gate_sensitivity": {
            row["counterfactual_policy"]: {
                "kept_scorer_signal_count": row["kept_scorer_signal_count"],
                "scorer_entries_with_at_least_one_signal": row[
                    "scorer_entries_with_at_least_one_signal"
                ],
                "scorer_entries_losing_all_signals": row[
                    "scorer_entries_losing_all_signals"
                ],
            }
            for row in selection_sensitivity
        },
        "current_semantics": {
            "Code_NTC_RK": (
                "Upstream eligibility gate. Any NTC00, RK00/RK0, or RK8 tag selects a positive observation. "
                "The supplied priority remains attached for audit, but the current test-aware scorer does not consume it as a weight."
            ),
            "adapter_and_centrifugation": (
                "Protocol metadata only. Reads/RPM are not rescaled and the labels are not direct scorer weights."
            ),
            "multiple_seq_ids": (
                "Indirect scorer effect: same-molecule repeated selected tests and DNA/RNA concordance currently share one reproducibility axis."
            ),
        },
        "recommended_disposition_before_lab_confirmation": {
            "Code_NTC_RK": "keep_as_source_filter_and_audit_metadata_only",
            "adapter_dilution": "metadata_warning_only_do_not_rescale_reads",
            "centrifugation": "metadata_only_until_fraction_and_aliquot_relationship_are_known",
            "rerun": "separate_same_library_resequence_from_new_extraction_before_weighting",
            "reproducibility": (
                "keep DNA/RNA concordance distinct; do not promote same-run fraction pairs or same-library reruns to Picked without a confirmed independence model"
            ),
        },
        "laboratory_questions": [
            "What exactly generate Code_NTC=00, RK00/RK0, and RK8, including thresholds and whether RK8 is weaker or a different state?",
            "For 12000g/10min paired tests, was pellet, supernatant, or another fraction sequenced, and did paired seq_ids come from the same original aliquot?",
            "For later seq_ids and 原library重新上機, was this the same library, a new library from the same extract, or a new extraction?",
            "Does adapter 2x/10x/half wording mean adapter-reagent dilution or sample/library dilution, and is any abundance correction validated?",
        ],
        "assertions": assertions,
        "assertion_pass_count": sum(assertions.values()),
        "assertion_total_count": len(assertions),
        "all_assertions_pass": all(assertions.values()),
    }
    write_json(output_dir / "summary.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit mNGS selection-code and protocol/repeat semantics"
    )
    parser.add_argument("--inventory-dir", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--compact-dir", type=Path, default=DEFAULT_COMPACT)
    parser.add_argument("--scorer-dir", type=Path, default=DEFAULT_SCORER)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-suffix", default=DEFAULT_SUMMARY_SUFFIX)
    args = parser.parse_args()
    result = audit(
        args.inventory_dir,
        args.compact_dir,
        args.scorer_dir,
        args.patient_root,
        args.policy,
        args.output_dir,
        summary_suffix=args.summary_suffix,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

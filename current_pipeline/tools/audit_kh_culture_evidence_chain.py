"""Audit KH culture evidence from normalized rows through scorer candidates.

The audit is answer blind.  It never changes the frozen agent, summary, or
scorer outputs.  It also emits a deterministic provenance-repair shadow that
preserves the current organism labels and causative levels while rebuilding
one culture observation per normalized source row.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import re
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any

from tools import normalized_agent_fallback as fallback
from tools import build_deterministic_summary as summary_builder
from tools import build_test_aware_deterministic_shadow_scorer as scorer_builder
from tools.pathogen_normalization import canonical_key


SCHEMA_VERSION = "kh_culture_evidence_chain_audit.v1"
DEFAULT_SCORER_ROOT = Path(
    "outputs/runs/2026-09-28_KH_ablation_G_new_scorer_taxonomy_v1/scorer"
)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def patient_number(path: Path) -> int:
    match = re.search(r"NGS_patient_(\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot identify patient number from {path}")
    return int(match.group(1))


def _norm(value: Any) -> str:
    return re.sub(r"\s+", "", str(value or "").strip().casefold())


def _source_signature(row: dict[str, Any]) -> tuple[str, ...]:
    return tuple(_norm(value) for value in (
        row.get("test"),
        row.get("sample") or row.get("specimen"),
        row.get("collected_time"),
        row.get("reported_time"),
        row.get("organism"),
        row.get("status"),
        row.get("colony_count") or row.get("quantity") or row.get("value"),
    ))


def _observation_signature(row: dict[str, Any]) -> tuple[str, ...]:
    return tuple(_norm(value) for value in (
        row.get("test"),
        row.get("specimen_type"),
        row.get("collected_time"),
        row.get("reported_time"),
        row.get("organism_name"),
        row.get("raw_result"),
        row.get("raw_quantity"),
    ))


def _is_positive_source(row: dict[str, Any]) -> bool:
    name = str(row.get("organism") or "").strip()
    status = str(row.get("status") or "").casefold()
    return bool(name) and not any(
        term in status
        for term in ("not detected", "negative", "no growth", "not isolated")
    )


def _labels_by_observation(
    labels: list[dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for label in labels:
        for observation_id in label.get("evidence_observation_ids") or []:
            result[str(observation_id)].append(label)
    return result


def _culture_summary_map(summary: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        canonical_key(item.get("organism_name")): item
        for item in summary.get("hospital_organism_evidence") or []
        if "culture" in (item.get("evidence_modules") or [])
    }


def _candidate_keys(scorer: dict[str, Any]) -> set[str]:
    rows = [
        *(scorer.get("patient_organism_decisions") or []),
        *(scorer.get("hospital_only_decisions") or []),
    ]
    return {canonical_key(row.get("organism_name")) for row in rows}


def repair_agent_provenance(
    current: dict[str, Any], source: list[dict[str, Any]]
) -> dict[str, Any]:
    """Repair row provenance without changing current labels or levels."""
    generated = fallback.culture_agent_from_source(source)
    generated_labels = {
        canonical_key(item.get("organism_name")): item
        for item in generated.get("organism_labels") or []
    }
    repaired = copy.deepcopy(current)
    repaired["culture_observations"] = generated["culture_observations"]
    repaired_labels: list[dict[str, Any]] = []
    unmatched_current: list[str] = []
    for current_label in current.get("organism_labels") or []:
        item = copy.deepcopy(current_label)
        generated_label = generated_labels.get(
            canonical_key(current_label.get("organism_name"))
        )
        if generated_label is None:
            unmatched_current.append(str(current_label.get("organism_name") or ""))
            repaired_labels.append(item)
            continue
        old_ids = list(item.get("evidence_observation_ids") or [])
        new_ids = list(generated_label.get("evidence_observation_ids") or [])
        item["evidence_observation_ids"] = new_ids
        item["source_complete_evidence_profile"] = copy.deepcopy(
            generated_label.get("key_evidence") or {}
        )
        if len(new_ids) > len(old_ids):
            preserved_level = item.get("causative_level")
            preserved_rules = copy.deepcopy(item.get("applied_rules") or [])
            item["key_evidence"] = copy.deepcopy(
                generated_label.get("key_evidence") or item.get("key_evidence") or {}
            )
            item["causative_level"] = preserved_level
            item["applied_rules"] = preserved_rules
            item["provenance_repair"] = {
                "old_evidence_observation_ids": old_ids,
                "new_evidence_observation_ids": new_ids,
                "semantic_level_changed": False,
            }
        repaired_labels.append(item)
    current_keys = {
        canonical_key(item.get("organism_name"))
        for item in current.get("organism_labels") or []
    }
    generated_only = [
        item.get("organism_name")
        for key, item in generated_labels.items()
        if key not in current_keys
    ]
    repaired["organism_labels"] = repaired_labels
    repaired["provenance_repair"] = {
        "schema_version": "culture_source_complete_provenance_shadow.v1",
        "answer_blind": True,
        "source_record_count": len(source),
        "observation_count_before": len(current.get("culture_observations") or []),
        "observation_count_after": len(generated.get("culture_observations") or []),
        "organism_label_count_before": len(current.get("organism_labels") or []),
        "organism_label_count_after": len(repaired_labels),
        "causative_levels_preserved": True,
        "unmatched_current_labels": unmatched_current,
        "generated_only_labels_not_auto_added": generated_only,
    }
    repaired["source_record_contract"] = generated["source_record_contract"]
    return repaired


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _decision_rows(payload: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    output: dict[tuple[str, str], dict[str, Any]] = {}
    for section in ("patient_organism_decisions", "hospital_only_decisions"):
        for row in payload.get(section) or []:
            key = (section, canonical_key(row.get("organism_name")))
            if key in output:
                raise ValueError(f"Duplicate scorer decision key: {key}")
            output[key] = row
    return output


def _compare_scorer_semantics(
    current_root: Path, shadow_root: Path
) -> dict[str, Any]:
    fields = (
        "decision", "integrated_level", "formal_pick_allowed",
        "forward_to_clinical_scorer",
    )
    differences: list[dict[str, Any]] = []
    current_count = 0
    shadow_count = 0
    for current_path in sorted(
        (current_root / "patient_outputs").glob(
            "NGS_patient_*_test_aware_deterministic_shadow.json"
        ),
        key=patient_number,
    ):
        patient = patient_number(current_path)
        shadow_path = (
            shadow_root / "patient_outputs"
            / f"NGS_patient_{patient}_test_aware_deterministic_shadow.json"
        )
        current = _decision_rows(read_json(current_path))
        shadow = _decision_rows(read_json(shadow_path))
        current_count += len(current)
        shadow_count += len(shadow)
        for key in sorted(set(current) | set(shadow)):
            if key not in current or key not in shadow:
                differences.append({
                    "patient_id": patient,
                    "decision_section": key[0],
                    "organism_key": key[1],
                    "field": "__candidate__",
                    "current_value": "present" if key in current else "missing",
                    "shadow_value": "present" if key in shadow else "missing",
                })
                continue
            for field in fields:
                current_value = current[key].get(field)
                shadow_value = shadow[key].get(field)
                if current_value != shadow_value:
                    differences.append({
                        "patient_id": patient,
                        "decision_section": key[0],
                        "organism_key": key[1],
                        "field": field,
                        "current_value": current_value,
                        "shadow_value": shadow_value,
                    })
    return {
        "passed": not differences,
        "current_candidate_count": current_count,
        "shadow_candidate_count": shadow_count,
        "difference_count": len(differences),
        "differences": differences,
    }


def run(
    scorer_root: Path,
    output_dir: Path,
    *,
    agent_suffix_tag: str = "evidence_v2_full33_20260917",
) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    paths = sorted(
        (scorer_root / "patient_outputs").glob(
            "NGS_patient_*_test_aware_deterministic_shadow.json"
        ),
        key=patient_number,
    )
    if not paths:
        raise ValueError(f"No scorer patient outputs found under {scorer_root}")
    output_dir.mkdir(parents=True, exist_ok=True)
    shadow_dir = output_dir / "repaired_culture_agent_shadow"
    shadow_dir.mkdir()
    summary_shadow_dir = output_dir / "repaired_hospital_summary_shadow"
    summary_shadow_dir.mkdir()
    summary_suffix = "culture_repaired_shadow_v1"
    chain_rows: list[dict[str, Any]] = []
    issue_rows: list[dict[str, Any]] = []
    label_rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    patient_root: Path | None = None

    for scorer_path in paths:
        patient = patient_number(scorer_path)
        scorer = read_json(scorer_path)
        summary_path = Path(str((scorer.get("source_files") or {}).get(
            "hospital_summary"
        ) or ""))
        summary = read_json(summary_path)
        source_files = summary.get("source_files") or {}
        source_path = Path(str(source_files.get("raw_culture") or ""))
        agent_path = Path(str(source_files.get("culture") or ""))
        if not source_path.is_file() or not agent_path.is_file():
            raise FileNotFoundError(
                f"P{patient}: culture source={source_path}, agent={agent_path}"
            )
        source = read_json(source_path)
        if not isinstance(source, list):
            raise ValueError(f"P{patient}: culture source must be a list")
        current = read_json(agent_path)
        shadow = fallback.culture_agent_from_source(source)
        repaired = repair_agent_provenance(current, source)
        repaired_path = (
            shadow_dir / f"NGS_patient_{patient}_culture_agent_repaired_shadow.json"
        )
        write_json(repaired_path, repaired)
        patient_dir = source_path.parent
        patient_root = patient_dir.parent
        shadow_summary = summary_builder.build_summary(
            patient_dir,
            agent_suffix_tag=agent_suffix_tag,
            agent_payload_overrides={"culture": repaired},
            agent_path_overrides={"culture": repaired_path.resolve()},
        )
        write_json(
            summary_shadow_dir / f"NGS_patient_{patient}_{summary_suffix}.json",
            shadow_summary,
        )

        current_by_signature: dict[tuple[str, ...], deque[dict[str, Any]]] = (
            defaultdict(deque)
        )
        for observation in current.get("culture_observations") or []:
            current_by_signature[_observation_signature(observation)].append(observation)
        shadow_by_index = {
            int(observation["source_record_index"]): observation
            for observation in shadow.get("culture_observations") or []
        }
        current_links = _labels_by_observation(current.get("organism_labels") or [])
        shadow_links = _labels_by_observation(shadow.get("organism_labels") or [])
        summary_map = _culture_summary_map(summary)
        candidate_keys = _candidate_keys(scorer)

        counts["patients"] += 1
        counts["source_rows"] += len(source)
        counts["current_observations"] += len(
            current.get("culture_observations") or []
        )
        counts["shadow_observations"] += len(
            shadow.get("culture_observations") or []
        )
        counts["current_labels"] += len(current.get("organism_labels") or [])
        counts["shadow_labels"] += len(shadow.get("organism_labels") or [])

        for index, row in enumerate(source):
            current_match = None
            matches = current_by_signature.get(_source_signature(row))
            if matches:
                current_match = matches.popleft()
            shadow_match = shadow_by_index.get(index)
            positive = _is_positive_source(row)
            current_id = (
                str(current_match.get("observation_id")) if current_match else ""
            )
            shadow_id = (
                str(shadow_match.get("observation_id")) if shadow_match else ""
            )
            current_label_items = current_links.get(current_id) or []
            shadow_label_items = shadow_links.get(shadow_id) or []
            organism = str(row.get("organism") or "")
            key = canonical_key(organism)
            summary_item = summary_map.get(key)
            candidate_present = bool(key and key in candidate_keys)
            summary_level = str(
                (summary_item or {}).get("best_hospital_level") or ""
            )
            if not positive:
                candidate_status = "not_expected_for_negative_or_no_growth"
            elif candidate_present:
                candidate_status = "present"
            elif summary_item is None:
                candidate_status = "not_in_hospital_summary_review_required"
            elif summary_level in {"Level 3", "Level 4", "Level 5"}:
                candidate_status = "hospital_only_below_candidate_entry_threshold"
            else:
                candidate_status = "candidate_gap_review_required"

            record_key = str(
                (shadow_match or {}).get("source_record_key")
                or f"culture[{index}]"
            )
            chain_rows.append({
                "patient_id": patient,
                "source_record_index": index,
                "source_record_key": record_key,
                "test": row.get("test"),
                "specimen": row.get("sample") or row.get("specimen"),
                "collected_time": row.get("collected_time"),
                "reported_time": row.get("reported_time"),
                "organism_name": organism,
                "status": row.get("status"),
                "quantity": row.get("colony_count") or row.get("quantity") or row.get("value"),
                "positive_or_provisional": positive,
                "current_observation_id": current_id,
                "current_observation_present": current_match is not None,
                "current_label_linked": bool(current_label_items),
                "current_label_names": "|".join(
                    str(item.get("organism_name") or "")
                    for item in current_label_items
                ),
                "hospital_summary_present": summary_item is not None,
                "hospital_summary_level": summary_level,
                "scorer_candidate_present": candidate_present,
                "candidate_status": candidate_status,
                "shadow_observation_id": shadow_id,
                "shadow_observation_present": shadow_match is not None,
                "shadow_label_linked": bool(shadow_label_items),
            })

            if current_match is None:
                counts["current_missing_source_rows"] += 1
                issue_rows.append({
                    "severity": "error",
                    "patient_id": patient,
                    "issue_type": "current_observation_missing_source_row",
                    "source_record_key": record_key,
                    "organism_name": organism,
                    "detail": (
                        f"{row.get('sample')} | {row.get('collected_time')} | "
                        f"{row.get('status')} | {row.get('colony_count') or ''}"
                    ),
                })
            if shadow_match is None:
                counts["shadow_missing_source_rows"] += 1
                issue_rows.append({
                    "severity": "error",
                    "patient_id": patient,
                    "issue_type": "shadow_observation_missing_source_row",
                    "source_record_key": record_key,
                    "organism_name": organism,
                    "detail": "deterministic source-complete shadow failed",
                })
            if positive and current_match is not None and not current_label_items:
                counts["positive_current_observations_without_label"] += 1
                issue_rows.append({
                    "severity": "error",
                    "patient_id": patient,
                    "issue_type": "positive_observation_not_linked_to_label",
                    "source_record_key": record_key,
                    "organism_name": organism,
                    "detail": current_id,
                })
            if positive and summary_item is None:
                counts["positive_source_rows_without_summary"] += 1
                issue_rows.append({
                    "severity": "error",
                    "patient_id": patient,
                    "issue_type": "positive_source_organism_missing_from_summary",
                    "source_record_key": record_key,
                    "organism_name": organism,
                    "detail": "culture label did not reach hospital summary",
                })
            if positive and not candidate_present and summary_item is not None:
                counts[f"candidate_policy_gap:{summary_level}"] += 1

        current_label_map = {
            canonical_key(item.get("organism_name")): item
            for item in current.get("organism_labels") or []
        }
        shadow_label_map = {
            canonical_key(item.get("organism_name")): item
            for item in shadow.get("organism_labels") or []
        }
        for key in sorted(set(current_label_map) | set(shadow_label_map)):
            current_label = current_label_map.get(key) or {}
            shadow_label = shadow_label_map.get(key) or {}
            current_ids = current_label.get("evidence_observation_ids") or []
            shadow_ids = shadow_label.get("evidence_observation_ids") or []
            current_level = current_label.get("causative_level")
            shadow_level = shadow_label.get("causative_level")
            label_rows.append({
                "patient_id": patient,
                "organism_name": current_label.get("organism_name") or shadow_label.get("organism_name"),
                "current_present": bool(current_label),
                "shadow_present": bool(shadow_label),
                "current_level": current_level,
                "shadow_rule_level_not_applied": shadow_level,
                "current_observation_ids": "|".join(current_ids),
                "shadow_complete_observation_ids": "|".join(shadow_ids),
                "observation_link_count_change": len(shadow_ids) - len(current_ids),
                "current_quantity_tier": (current_label.get("key_evidence") or {}).get("quantity_tier"),
                "shadow_quantity_tier": (shadow_label.get("key_evidence") or {}).get("quantity_tier"),
                "semantic_level_preserved_in_repaired_shadow": True,
            })
            if current_level != shadow_level:
                counts["deterministic_rule_level_differences_not_applied"] += 1

        repaired_levels = {
            canonical_key(item.get("organism_name")): item.get("causative_level")
            for item in repaired.get("organism_labels") or []
        }
        current_levels = {
            canonical_key(item.get("organism_name")): item.get("causative_level")
            for item in current.get("organism_labels") or []
        }
        if repaired_levels != current_levels:
            issue_rows.append({
                "severity": "error",
                "patient_id": patient,
                "issue_type": "repaired_shadow_changed_semantic_levels",
                "source_record_key": "",
                "organism_name": "",
                "detail": "current and repaired label-level maps differ",
            })

    chain_fields = [
        "patient_id", "source_record_index", "source_record_key", "test",
        "specimen", "collected_time", "reported_time", "organism_name",
        "status", "quantity", "positive_or_provisional",
        "current_observation_id", "current_observation_present",
        "current_label_linked", "current_label_names",
        "hospital_summary_present", "hospital_summary_level",
        "scorer_candidate_present", "candidate_status",
        "shadow_observation_id", "shadow_observation_present",
        "shadow_label_linked",
    ]
    issue_fields = [
        "severity", "patient_id", "issue_type", "source_record_key",
        "organism_name", "detail",
    ]
    label_fields = [
        "patient_id", "organism_name", "current_present", "shadow_present",
        "current_level", "shadow_rule_level_not_applied",
        "current_observation_ids", "shadow_complete_observation_ids",
        "observation_link_count_change", "current_quantity_tier",
        "shadow_quantity_tier", "semantic_level_preserved_in_repaired_shadow",
    ]
    _write_csv(output_dir / "culture_chain_rows.csv", chain_rows, chain_fields)
    _write_csv(output_dir / "culture_chain_issues.csv", issue_rows, issue_fields)
    _write_csv(output_dir / "culture_label_comparison.csv", label_rows, label_fields)

    if patient_root is None:
        raise ValueError("No patient root resolved")
    current_scorer_summary = read_json(scorer_root / "summary.json")
    compact_entry_dir = Path(
        str(current_scorer_summary.get("source_compact_entry_dir") or "")
    )
    policy_path = Path(str(current_scorer_summary.get("policy_path") or ""))
    scorer_shadow_dir = output_dir / "repaired_scorer_shadow"
    scorer_builder.build_shadow(
        compact_entry_dir,
        patient_root,
        scorer_shadow_dir,
        summary_suffix=summary_suffix,
        summary_root=summary_shadow_dir,
        policy_path=policy_path,
    )
    scorer_parity = _compare_scorer_semantics(scorer_root, scorer_shadow_dir)
    scorer_difference_fields = [
        "patient_id", "decision_section", "organism_key", "field",
        "current_value", "shadow_value",
    ]
    _write_csv(
        output_dir / "scorer_semantic_differences.csv",
        scorer_parity["differences"],
        scorer_difference_fields,
    )

    p24_summary_path = (
        summary_shadow_dir / f"NGS_patient_24_{summary_suffix}.json"
    )
    p24_trichosporon_ids: list[str] = []
    p24_trichosporon_quantity = None
    if p24_summary_path.is_file():
        p24_summary = read_json(p24_summary_path)
        for item in p24_summary.get("hospital_organism_evidence") or []:
            if canonical_key(item.get("organism_name")) != canonical_key(
                "Trichosporon asahii"
            ):
                continue
            culture_rows = (item.get("module_evidence") or {}).get("culture") or []
            for culture_row in culture_rows:
                p24_trichosporon_ids.extend(
                    culture_row.get("evidence_observation_ids") or []
                )
                p24_trichosporon_quantity = culture_row.get("quantity_tier")

    error_count = sum(row["severity"] == "error" for row in issue_rows)
    repaired_shadow_error_count = (
        counts["shadow_missing_source_rows"]
        + scorer_parity["difference_count"]
        + (0 if (
            len(set(p24_trichosporon_ids)) == 4
            and p24_trichosporon_quantity == "Q4"
        ) else 1)
    )
    repaired_shadow_error_count = (
        counts["shadow_missing_source_rows"]
        + scorer_parity["difference_count"]
        + (0 if (
            len(set(p24_trichosporon_ids)) == 4
            and p24_trichosporon_quantity == "Q4"
        ) else 1)
    )
    summary = {
        "schema_version": SCHEMA_VERSION,
        "answer_blind": True,
        "scorer_root": str(scorer_root.resolve()),
        "patient_count": counts["patients"],
        "source_record_count": counts["source_rows"],
        "current_observation_count": counts["current_observations"],
        "source_rows_missing_from_current_observations": counts[
            "current_missing_source_rows"
        ],
        "deterministic_shadow_observation_count": counts["shadow_observations"],
        "source_rows_missing_from_deterministic_shadow": counts[
            "shadow_missing_source_rows"
        ],
        "current_organism_label_count": counts["current_labels"],
        "deterministic_shadow_organism_label_count": counts["shadow_labels"],
        "positive_current_observations_without_label": counts[
            "positive_current_observations_without_label"
        ],
        "positive_source_rows_without_hospital_summary": counts[
            "positive_source_rows_without_summary"
        ],
        "deterministic_rule_level_differences_not_applied": counts[
            "deterministic_rule_level_differences_not_applied"
        ],
        "candidate_policy_gap_counts": {
            key.removeprefix("candidate_policy_gap:"): value
            for key, value in sorted(counts.items())
            if key.startswith("candidate_policy_gap:")
        },
        "error_count": error_count,
        "current_frozen_chain_error_count": error_count,
        "repaired_shadow_error_count": repaired_shadow_error_count,
        "repair_ready_for_rule_review": repaired_shadow_error_count == 0,
        "current_frozen_chain_error_count": error_count,
        "repaired_shadow_error_count": repaired_shadow_error_count,
        "repair_ready_for_rule_review": repaired_shadow_error_count == 0,
        "repair_shadow": {
            "root": str(shadow_dir.resolve()),
            "semantic_levels_preserved": True,
            "production_outputs_overwritten": False,
        },
        "repaired_summary_shadow_root": str(summary_shadow_dir.resolve()),
        "repaired_scorer_shadow_root": str(scorer_shadow_dir.resolve()),
        "scorer_semantic_parity": {
            key: value for key, value in scorer_parity.items()
            if key != "differences"
        },
        "p24_trichosporon_repair": {
            "hospital_summary_observation_ids": list(dict.fromkeys(
                p24_trichosporon_ids
            )),
            "quantity_tier": p24_trichosporon_quantity,
            "expected_observation_count": 4,
            "passed": (
                len(set(p24_trichosporon_ids)) == 4
                and p24_trichosporon_quantity == "Q4"
            ),
        },
        "limitations": [
            "This run starts from normalized culture JSON; original workbook-to-normalized completeness remains unverified where source workbooks are unavailable.",
            "Deterministic shadow rule-level differences are reported but are not applied to the frozen pipeline.",
            "Hospital-only Level 3-5 records may remain outside the scorer candidate universe by current policy; they remain visible in the chain audit.",
        ],
    }
    write_json(output_dir / "summary.json", summary)
    report = [
        "# KH Culture evidence-chain audit v1",
        "",
        "This is an answer-blind audit. No frozen output was overwritten.",
        "",
        f"- Patients: `{summary['patient_count']}`",
        f"- Normalized culture records: `{summary['source_record_count']}`",
        f"- Current culture observations: `{summary['current_observation_count']}`",
        f"- Missing from current observations: `{summary['source_rows_missing_from_current_observations']}`",
        f"- Deterministic shadow observations: `{summary['deterministic_shadow_observation_count']}`",
        f"- Missing from deterministic shadow: `{summary['source_rows_missing_from_deterministic_shadow']}`",
        f"- Current / shadow organism labels: `{summary['current_organism_label_count']}` / `{summary['deterministic_shadow_organism_label_count']}`",
        f"- Positive observations without a label: `{summary['positive_current_observations_without_label']}`",
        f"- Positive source rows absent from hospital summary: `{summary['positive_source_rows_without_hospital_summary']}`",
        f"- Audit errors: `{summary['error_count']}`",
        f"- Repaired-shadow errors: `{summary['repaired_shadow_error_count']}`",
        f"- Repaired-shadow errors: `{summary['repaired_shadow_error_count']}`",
        f"- Repaired scorer semantic parity: `{summary['scorer_semantic_parity']['passed']}`",
        f"- Repaired scorer semantic differences: `{summary['scorer_semantic_parity']['difference_count']}`",
        f"- P24 Trichosporon repair: `{summary['p24_trichosporon_repair']['passed']}`",
        "",
        "The repaired shadow rebuilds provenance from source rows while preserving the current organism labels and causative levels. Rule-level differences from the independent deterministic classifier are audit-only and were not applied.",
    ]
    (output_dir / "README.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--scorer-root", type=Path, default=DEFAULT_SCORER_ROOT)
    parser.add_argument(
        "--agent-suffix-tag", default="evidence_v2_full33_20260917"
    )
    args = parser.parse_args()
    print(json.dumps(
        run(
            args.scorer_root,
            args.output_dir,
            agent_suffix_tag=args.agent_suffix_tag,
        ), ensure_ascii=False, indent=2
    ))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build a conservative v14-to-current evidence compatibility adapter.

For KH33, the current reviewed hospital evidence is the non-regression base.
New v14 Culture/targeted-molecular evidence is additive only when it has a
mapped canonical organism identity, explicit source observations, a preserved
specimen, and a collection/report timestamp within +/-48 hours of a fixed mNGS
index event.  Unmapped or non-aligned evidence remains audit-only.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import re
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

from tools.clinical_event_window import nearest_calendar_event

from tools.organism_taxonomy_classifier import classify_organism


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PATIENT_ROOT = ROOT / "outputs/kh_v14_linkage_gated_shadow_input_v1/patients"
DEFAULT_BRIDGE_ROOT = ROOT / "outputs/runs/2026-10-07_v14_hospital_evidence_bridge_v2"
DEFAULT_BASELINE_RELEASE = ROOT / "outputs/runs/2026-10-07_KH33_opportunistic_family_fallback_shadow_v2"
DEFAULT_CURRENT_PATIENT_ROOT = ROOT / "outputs/patient_info_KH_0728_2Days"
DEFAULT_INVENTORY = ROOT / "outputs/runs/2026-10-07_51patients_inventory_kh33_frozen_mngs_v1"
DEFAULT_LINKAGE = ROOT / "outputs/kh_v14_linkage_gated_shadow_input_v1/input_manifest.json"
DEFAULT_OUTPUT = ROOT / "outputs/kh_v14_compatibility_adapter_v4"
SOURCE_V14_SUFFIX = "final_summary_v14_linkage_gated_20261007"
OUTPUT_SUMMARY_SUFFIX = "final_summary_v14_compatibility_v4_20261007"
ALLOWED_NEW_MODULES = {"culture", "molecular_microbiology"}
ALLOWED_SPECIMEN_CATEGORIES = {"Lower_Respiratory", "Sterile_Site", "CNS"}
LEVEL_RANK = {"Level 1": 1, "Level 2": 2, "Level 3": 3, "Level 4": 4, "Not_available": 99, "Unknown": 99, "": 99}


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_time(value: Any) -> datetime | None:
    text = str(value or "").strip()
    if not text or text.casefold() in {"unknown", "na", "n/a", "none", "null"}:
        return None
    text = text.replace("T", " ").replace("/", "-")
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None


def patient_id_from_name(value: str) -> int:
    match = re.search(r"(?:NGS_patient_|patient_)(\d+)", value)
    if not match:
        raise ValueError(f"Cannot parse patient ID from {value!r}")
    return int(match.group(1))


def best_level(values: list[str]) -> str:
    return min(values, key=lambda value: LEVEL_RANK.get(value, 99), default="Not_available")


def event_times(inventory_root: Path, patient_id: int) -> list[tuple[str, datetime]]:
    path = inventory_root / "patients" / f"NGS_patient_{patient_id}_multi_assay_mngs_evidence.json"
    payload = read_json(path)
    result: list[tuple[str, datetime]] = []
    for index, case in enumerate(payload.get("cases") or []):
        value = parse_time(case.get("collected_time"))
        if value is not None:
            result.append((str(case.get("case_review_id") or f"case_{index + 1}"), value))
    if not result:
        raise ValueError(f"P{patient_id} has no parseable mNGS collection time")
    return result


def alignment(row: dict[str, Any], events: list[tuple[str, datetime]]) -> dict[str, Any]:
    collected = parse_time(row.get("collected_time"))
    reported = parse_time(row.get("reported_time"))
    evidence_time = collected or reported
    basis = "collected_time" if collected is not None else "reported_time_proxy" if reported is not None else "missing"
    if evidence_time is None:
        return {"passed": False, "reason": "time_unavailable", "time_basis": basis, "event_id": "", "delta_hours": None}
    nearest = nearest_calendar_event(evidence_time, (event_time for _, event_time in events))
    assert nearest is not None
    event_time = nearest["event"]
    event_id = next(event_id for event_id, value in events if value == event_time)
    delta = float(nearest["signed_delta_hours"])
    passed = bool(nearest["within_window"])
    return {
        "passed": passed,
        "reason": "within_two_calendar_dates" if passed else "outside_two_calendar_dates",
        "time_basis": basis,
        "event_id": event_id,
        "delta_hours": round(delta, 3),
        "calendar_day_delta": nearest["signed_calendar_day_delta"],
        "event_window_mode": nearest["window_mode"],
    }


def raw_source_path(patient_dir: Path, patient_id: int, module: str) -> Path:
    suffix = "culture" if module == "culture" else "molecular_microbiology"
    return patient_dir / f"NGS_patient_{patient_id}_{suffix}.json"


def source_references(
    patient_dir: Path,
    patient_id: int,
    module: str,
    observation_ids: list[str],
) -> list[dict[str, Any]]:
    path = raw_source_path(patient_dir, patient_id, module)
    if not path.is_file():
        return []
    payload = read_json(path)
    rows = payload if isinstance(payload, list) else [payload] if isinstance(payload, dict) else []
    prefix = "CUL" if module == "culture" else "MOL"
    result: list[dict[str, Any]] = []
    for observation_id in observation_ids:
        match = re.fullmatch(fr"{prefix}-(\d+)", str(observation_id))
        if not match:
            continue
        index = int(match.group(1)) - 1
        if not 0 <= index < len(rows):
            continue
        row = rows[index] if isinstance(rows[index], dict) else {}
        source = row.get("source") if isinstance(row.get("source"), dict) else {}
        result.append({
            "observation_id": observation_id,
            "source_file": str(path.resolve()),
            "source_file_sha256": sha256(path),
            "source_record_index": index,
            "record_id": row.get("record_id"),
            "source_workbook": source.get("source_file"),
            "source_sheet": source.get("source_sheet"),
            "source_row": source.get("row"),
        })
    return result


def evidence_key(item: dict[str, Any]) -> str:
    profile = classify_organism(item.get("organism_name"), biological_class=item.get("classification"))
    return str(profile.get("canonical_key") or "")


def taxonomy_profile(item: dict[str, Any]) -> dict[str, Any]:
    return classify_organism(item.get("organism_name"), biological_class=item.get("classification"))


def is_mapped(profile: dict[str, Any]) -> bool:
    return bool(
        profile.get("canonical_key")
        and profile.get("mapping_status") != "unmapped"
        and profile.get("primary_rule_family") != "unmapped_or_uncertain"
    )


def is_explicit_positive(row: dict[str, Any]) -> bool:
    """Require a normalized positive result before evidence can leave Audit."""
    values = [
        str(row.get("result") or "").strip().casefold(),
        str(row.get("detection_status") or "").strip().casefold(),
        str(row.get("quantitation_status") or "").strip().casefold(),
    ]
    negative_tokens = ("negative", "not detected", "non-reactive", "nonreactive", "no growth")
    if any(any(token in value for token in negative_tokens) for value in values if value):
        return False
    positive_tokens = ("positive", "detected", "reactive", "growth")
    return any(any(token in value for token in positive_tokens) for value in values if value)


def recalculate(item: dict[str, Any]) -> None:
    levels = [str(value) for value in (item.get("module_level_summary") or {}).values()]
    item["best_hospital_level"] = best_level(levels)


def merge_evidence(base: dict[str, Any], incoming: dict[str, Any]) -> None:
    for module in incoming.get("evidence_modules") or []:
        if module not in base.setdefault("evidence_modules", []):
            base["evidence_modules"].append(module)
    base_levels = base.setdefault("module_level_summary", {})
    for module, level in (incoming.get("module_level_summary") or {}).items():
        if module == "cbc_other_lab":
            base_levels.setdefault(module, level)
            continue
        current = str(base_levels.get(module) or "Not_available")
        base_levels[module] = best_level([current, str(level or "Not_available")])
    base_rows = base.setdefault("module_evidence", {})
    for module, rows in (incoming.get("module_evidence") or {}).items():
        destination = base_rows.setdefault(module, [])
        seen = {json.dumps(row, ensure_ascii=False, sort_keys=True) for row in destination if isinstance(row, dict)}
        for row in rows if isinstance(rows, list) else []:
            signature = json.dumps(row, ensure_ascii=False, sort_keys=True)
            if signature not in seen:
                destination.append(row)
                seen.add(signature)
    base.setdefault("evidence_trace", []).extend(copy.deepcopy(incoming.get("evidence_trace") or []))
    recalculate(base)


def audit_row(
    *,
    patient_id: int,
    input_name: str,
    canonical_name: str,
    module: str,
    source_basis: str,
    disposition: str,
    reason: str,
    profile: dict[str, Any],
    row: dict[str, Any] | None = None,
    aligned: dict[str, Any] | None = None,
    refs: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    row = row or {}
    aligned = aligned or {}
    return {
        "patient_id": patient_id,
        "input_organism_name": input_name,
        "canonical_organism_name": canonical_name,
        "module": module,
        "source_basis": source_basis,
        "disposition": disposition,
        "reason": reason,
        "taxonomy_mapping_status": profile.get("mapping_status"),
        "taxonomy_family": profile.get("primary_rule_family"),
        "taxid": profile.get("taxid"),
        "specimen_type": row.get("specimen_type"),
        "specimen_category": row.get("specimen_category"),
        "positive_result": row.get("result") or row.get("detection_status") or row.get("quantitation_status"),
        "collected_time": row.get("collected_time"),
        "reported_time": row.get("reported_time"),
        "time_basis": aligned.get("time_basis"),
        "matched_event_id": aligned.get("event_id"),
        "delta_hours": aligned.get("delta_hours"),
        "source_observation_ids": "|".join(str(value) for value in row.get("evidence_observation_ids") or []),
        "source_rows": "|".join(
            f"{ref.get('source_workbook')}:{ref.get('source_sheet')}:{ref.get('source_row')}"
            for ref in refs or []
        ),
    }


def baseline_summary_path(release_root: Path, patient_id: int) -> Path:
    candidates = [
        path for path in (release_root / "00_current_source_summary").glob(f"NGS_patient_{patient_id}_*.json")
        if path.name != "summary.json"
    ]
    if len(candidates) != 1:
        raise ValueError(f"Expected one baseline summary for P{patient_id}; found {candidates}")
    return candidates[0]


def source_v14_summary_path(patient_dir: Path, patient_id: int) -> Path:
    return patient_dir / "summary_outputs" / f"NGS_patient_{patient_id}_{SOURCE_V14_SUFFIX}.json"


def adapt_patient(
    *,
    patient_id: int,
    patient_dir: Path,
    bridge_packet: dict[str, Any],
    events: list[tuple[str, datetime]],
    baseline_summary: dict[str, Any] | None,
    v14_summary: dict[str, Any],
    allow_v14: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    output = copy.deepcopy(baseline_summary if baseline_summary is not None else v14_summary)
    evidence_by_key: dict[str, dict[str, Any]] = {}
    audits: list[dict[str, Any]] = []
    preserved_keys: set[str] = set()
    preserved_best_rank: dict[str, int] = {}

    for item in copy.deepcopy((baseline_summary or {}).get("hospital_organism_evidence") or []):
        profile = taxonomy_profile(item)
        input_name = str(item.get("organism_name") or "")
        if not is_mapped(profile):
            audits.append(audit_row(
                patient_id=patient_id,
                input_name=input_name,
                canonical_name=str(profile.get("display_name") or input_name),
                module="|".join(item.get("evidence_modules") or []),
                source_basis="current_KH_reviewed_summary",
                disposition="audit_only",
                reason="taxonomy_unmapped_or_uncertain",
                profile=profile,
            ))
            continue
        item["organism_name"] = profile["display_name"]
        item["taxonomy_profile"] = profile
        item["compatibility_adapter"] = {"basis": "preserved_current_verified_hospital_evidence", "non_regression_protected": True}
        key = str(profile["canonical_key"])
        preserved_keys.add(key)
        preserved_best_rank[key] = min(
            preserved_best_rank.get(key, 99),
            LEVEL_RANK.get(str(item.get("best_hospital_level") or "Not_available"), 99),
        )
        if key in evidence_by_key:
            merge_evidence(evidence_by_key[key], item)
        else:
            evidence_by_key[key] = item
        audits.append(audit_row(
            patient_id=patient_id,
            input_name=input_name,
            canonical_name=str(profile["display_name"]),
            module="|".join(item.get("evidence_modules") or []),
            source_basis="current_KH_reviewed_summary",
            disposition="preserved",
            reason="mapped_current_evidence_non_regression_base",
            profile=profile,
        ))

    for item in bridge_packet.get("hospital_organism_evidence") or []:
        input_name = str(item.get("organism_name") or "")
        profile = taxonomy_profile(item)
        canonical_name = str(profile.get("display_name") or input_name)
        for module, rows in (item.get("module_evidence") or {}).items():
            if module not in ALLOWED_NEW_MODULES:
                continue
            for source_row in rows if isinstance(rows, list) else []:
                observation_ids = [str(value) for value in source_row.get("evidence_observation_ids") or []]
                refs = source_references(patient_dir, patient_id, module, observation_ids)
                aligned = alignment(source_row, events)
                reason = "accepted"
                if not allow_v14:
                    reason = "patient_linkage_not_approved_for_v14"
                elif not is_mapped(profile):
                    reason = "taxonomy_unmapped_or_uncertain"
                elif not observation_ids or len(refs) != len(observation_ids):
                    reason = "source_row_reference_incomplete"
                elif not str(source_row.get("specimen_type") or "").strip():
                    reason = "specimen_missing"
                elif not is_explicit_positive(source_row):
                    reason = "positive_result_not_explicit"
                elif str(source_row.get("specimen_category") or "") not in ALLOWED_SPECIMEN_CATEGORIES:
                    reason = "specimen_not_pulmonary_or_invasive"
                elif not aligned["passed"]:
                    reason = aligned["reason"]
                disposition = "accepted" if reason == "accepted" else "audit_only"
                audits.append(audit_row(
                    patient_id=patient_id,
                    input_name=input_name,
                    canonical_name=canonical_name,
                    module=module,
                    source_basis="v14_normalized_raw_bridge",
                    disposition=disposition,
                    reason=reason,
                    profile=profile,
                    row=source_row,
                    aligned=aligned,
                    refs=refs,
                ))
                if disposition != "accepted":
                    continue
                accepted_row = copy.deepcopy(source_row)
                accepted_row["source_references"] = refs
                accepted_row["event_alignment"] = aligned
                accepted = {
                    "organism_name": canonical_name,
                    "classification": item.get("classification"),
                    "evidence_modules": [module],
                    "module_level_summary": {
                        "culture": "Not_available",
                        "filmarray_gmtest": "Not_available",
                        "image": "Not_available",
                        "molecular_microbiology": "Not_available",
                        "cbc_other_lab": "Not_pathogen_specific",
                        module: (item.get("module_level_summary") or {}).get(module, "Unknown"),
                    },
                    "module_evidence": {module: [accepted_row]},
                    "evidence_trace": [{
                        "step": "v14_compatibility_adapter",
                        "rule_id": "V14-COMPAT-CANONICAL-EVENT-ALIGNED",
                        "result": "pass",
                        "evidence": [{"field": "event_alignment.delta_hours", "value": aligned["delta_hours"]}],
                    }],
                    "taxonomy_profile": profile,
                    "compatibility_adapter": {"basis": "v14_normalized_event_aligned_addition", "non_regression_protected": False},
                }
                recalculate(accepted)
                key = str(profile["canonical_key"])
                if key in evidence_by_key:
                    merge_evidence(evidence_by_key[key], accepted)
                else:
                    evidence_by_key[key] = accepted

    output["hospital_organism_evidence"] = sorted(
        evidence_by_key.values(), key=lambda item: str(item.get("organism_name") or "").casefold()
    )
    output.setdefault("evidence_preservation", {})["v14_compatibility_adapter"] = {
        "version": "v3_canonical_positive_event_aligned_non_regression",
        "baseline_current_evidence_present": baseline_summary is not None,
        "baseline_mapped_key_count": len(preserved_keys),
        "final_mapped_key_count": len(evidence_by_key),
        "accepted_v14_row_count": sum(row["disposition"] == "accepted" for row in audits),
        "audit_only_row_count": sum(row["disposition"] == "audit_only" for row in audits),
        "rules": [
            "taxonomy_unmapped_or_uncertain_is_audit_only",
            "new_culture_or_targeted_molecular_requires_source_row_specimen_positive_normalization_and_plus_minus_48h_alignment",
            "mapped_current_KH_hospital_evidence_is_non_regression_protected",
        ],
    }
    final_keys = {evidence_key(item) for item in output["hospital_organism_evidence"]}
    if not preserved_keys.issubset(final_keys):
        raise AssertionError(f"P{patient_id} lost mapped current hospital evidence: {sorted(preserved_keys-final_keys)}")
    for key, baseline_rank in preserved_best_rank.items():
        final_rank = LEVEL_RANK.get(str(evidence_by_key[key].get("best_hospital_level") or "Not_available"), 99)
        if final_rank > baseline_rank:
            raise AssertionError(
                f"P{patient_id} weakened current hospital evidence {key}: "
                f"baseline_rank={baseline_rank}, final_rank={final_rank}"
            )
    if any(not is_mapped(taxonomy_profile(item)) for item in output["hospital_organism_evidence"]):
        raise AssertionError(f"P{patient_id} output contains unmapped hospital evidence")
    patient_summary = {
        "patient_id": patient_id,
        "baseline_mapped_key_count": len(preserved_keys),
        "final_mapped_key_count": len(evidence_by_key),
        "accepted_v14_row_count": sum(row["disposition"] == "accepted" for row in audits),
        "audit_only_row_count": sum(row["disposition"] == "audit_only" for row in audits),
        "non_regression_level_passed": True,
    }
    return output, audits, patient_summary


def copy_patient_root_files(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for path in source.glob("*.json"):
        shutil.copy2(path, destination / path.name)


def build(
    output_root: Path,
    *,
    patient_root: Path,
    bridge_root: Path,
    baseline_release: Path,
    current_patient_root: Path,
    inventory_root: Path,
    linkage_manifest: Path,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    linkage = read_json(linkage_manifest)
    linkage_by_patient = {int(item["patient_id"]): item for item in linkage["patients"]}
    patient_dirs = sorted(patient_root.glob("NGS_patient_*_json"), key=lambda path: patient_id_from_name(path.name))
    audits: list[dict[str, Any]] = []
    patients: list[dict[str, Any]] = []
    for source_dir in patient_dirs:
        patient_id = patient_id_from_name(source_dir.name)
        link = linkage_by_patient[patient_id]
        is_kh33 = link["cohort_role"] == "KH33_development"
        allow_v14 = (not is_kh33) or link["linkage_status"] == "supported_by_local_fingerprint"
        v14_summary_path = source_v14_summary_path(source_dir, patient_id)
        baseline_path = baseline_summary_path(baseline_release, patient_id) if is_kh33 else None
        baseline = read_json(baseline_path) if baseline_path is not None else None
        v14_summary = read_json(v14_summary_path)
        bridge_path = bridge_root / "patient_packets" / f"NGS_patient_{patient_id}_hospital_evidence_bridge.json"
        bridge = read_json(bridge_path)
        adapted, patient_audits, patient_summary = adapt_patient(
            patient_id=patient_id,
            patient_dir=source_dir,
            bridge_packet=bridge,
            events=event_times(inventory_root, patient_id),
            baseline_summary=baseline,
            v14_summary=v14_summary,
            allow_v14=allow_v14,
        )
        destination = output_root / "patients" / source_dir.name
        if is_kh33:
            raw_base_dir = current_patient_root / f"NGS_patient_{patient_id}_json"
            if not raw_base_dir.is_dir():
                raise FileNotFoundError(f"Missing current verified KH raw patient directory: {raw_base_dir}")
        else:
            raw_base_dir = source_dir
        copy_patient_root_files(raw_base_dir, destination)
        summary_path = destination / "summary_outputs" / f"NGS_patient_{patient_id}_{OUTPUT_SUMMARY_SUFFIX}.json"
        write_json(summary_path, adapted)
        audits.extend(patient_audits)
        patients.append({
            **patient_summary,
            "cohort_role": link["cohort_role"],
            "linkage_status": link["linkage_status"],
            "v14_additions_allowed": allow_v14,
            "clinical_raw_non_regression_base": str(raw_base_dir.resolve()),
            "baseline_summary": str(baseline_path.resolve()) if baseline_path else None,
            "v14_summary": str(v14_summary_path.resolve()),
            "bridge_packet": str(bridge_path.resolve()),
            "output_summary": str(summary_path.resolve()),
            "output_summary_sha256": sha256(summary_path),
        })

    output_root.mkdir(parents=True, exist_ok=True)
    audit_path = output_root / "evidence_audit.csv"
    fields = list(audits[0]) if audits else []
    with audit_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(audits)
    counts = Counter(f"{row['disposition']}:{row['reason']}" for row in audits)
    result = {
        "schema_version": "v14_compatibility_adapter.v4",
        "status": "shadow_only_not_production",
        "output_summary_suffix": OUTPUT_SUMMARY_SUFFIX,
        "patient_count": len(patients),
        "kh33_patient_count": sum(item["cohort_role"] == "KH33_development" for item in patients),
        "accepted_v14_row_count": sum(row["disposition"] == "accepted" for row in audits),
        "audit_only_row_count": sum(row["disposition"] == "audit_only" for row in audits),
        "audit_counts": dict(sorted(counts.items())),
        "non_regression_passed": True,
        "non_regression_level_passed": all(item["non_regression_level_passed"] for item in patients),
        "unmapped_output_count": 0,
        "sources": {
            "patient_root": str(patient_root.resolve()),
            "bridge_root": str(bridge_root.resolve()),
            "baseline_release": str(baseline_release.resolve()),
            "current_patient_root": str(current_patient_root.resolve()),
            "inventory_root": str(inventory_root.resolve()),
            "linkage_manifest": str(linkage_manifest.resolve()),
        },
        "evidence_audit": str(audit_path.resolve()),
        "patients": patients,
    }
    write_json(output_root / "adapter_manifest.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--bridge-root", type=Path, default=DEFAULT_BRIDGE_ROOT)
    parser.add_argument("--baseline-release", type=Path, default=DEFAULT_BASELINE_RELEASE)
    parser.add_argument("--current-patient-root", type=Path, default=DEFAULT_CURRENT_PATIENT_ROOT)
    parser.add_argument("--inventory-root", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--linkage-manifest", type=Path, default=DEFAULT_LINKAGE)
    args = parser.parse_args()
    result = build(
        args.output_root,
        patient_root=args.patient_root,
        bridge_root=args.bridge_root,
        baseline_release=args.baseline_release,
        current_patient_root=args.current_patient_root,
        inventory_root=args.inventory_root,
        linkage_manifest=args.linkage_manifest,
    )
    print(json.dumps({key: result[key] for key in (
        "patient_count", "kh33_patient_count", "accepted_v14_row_count",
        "audit_only_row_count", "audit_counts", "non_regression_passed", "non_regression_level_passed",
        "unmapped_output_count",
    )}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

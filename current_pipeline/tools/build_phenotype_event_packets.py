"""Align physician phenotype evidence to fixed mNGS events and candidates.

This is a provenance-preserving shadow export, not a pathogen scorer.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from tools.audit_phenotype_patient_linkage import lab_feature
from tools.build_index_event_clinical_timeline_shadow import episode_window_relation, sha256_file
from tools.pathogen_normalization import canonical_key


SCHEMA_VERSION = "phenotype_event_candidate_packets_shadow.v2"
WINDOW_STATUSES = {"within_window", "boundary_uncertain", "outside_window", "unknown"}
SHEETS = ("Phenotypes_Long", "Evidence", "Clinical_Facts", "Clinical_Timeline", "Microbiology", "Lab_Results")
FORBIDDEN_KEYS = {"answer", "answers", "benchmark_answer", "answer_hit", "is_answer"}


def find_forbidden_keys(value: Any, path: str = "$") -> list[str]:
    if isinstance(value, dict):
        found = []
        for key, child in value.items():
            if str(key).lower() in FORBIDDEN_KEYS:
                found.append(f"{path}.{key}")
            found.extend(find_forbidden_keys(child, f"{path}.{key}"))
        return found
    if isinstance(value, list):
        return [key for index, child in enumerate(value)
                for key in find_forbidden_keys(child, f"{path}[{index}]")]
    return []


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def patient_number(value: Any) -> int:
    match = re.match(r"\s*(?:patient[\s_]*|P)(\d+)(?=$|[\s_-])", str(value), re.I)
    if not match:
        raise ValueError(f"No patient number in {value!r}")
    return int(match.group(1))


def normalize_text(value: Any) -> str:
    return " ".join(str(value or "").casefold().split())


def overlap_ids(phenotype_text: Any, clinical_rows: list[dict[str, Any]]) -> list[str]:
    text = normalize_text(phenotype_text)
    if len(text) < 24:
        return []
    return [str(item["evidence_id"]) for item in clinical_rows
            if (source := normalize_text(item.get("source_text"))) and len(source) >= 24
            and (source in text or text in source)]


def lab_overlap_ids(item: dict[str, Any], clinical_rows: list[dict[str, Any]]) -> list[str]:
    row = item.get("row") or {}
    key = lab_feature(row.get("date") or row.get("datetime"), row.get("test_name"), row.get("value"))
    if key is None:
        return []
    matches = []
    for clinical in clinical_rows:
        if clinical.get("evidence_role") != "illness_lab":
            continue
        clinical_key = lab_feature(clinical.get("reported_at"), clinical.get("label"), clinical.get("value"))
        if key == clinical_key:
            matches.append(str(clinical["evidence_id"]))
    return matches


def record_date(sheet: str, row: dict[str, Any]) -> Any:
    if sheet in {"Phenotypes_Long", "Evidence"}:
        return row.get("source_date")
    return row.get("date") or row.get("datetime")


def record_text(sheet: str, row: dict[str, Any]) -> Any:
    return row.get("evidence") or row.get("description") or row.get("reason")


def compact_record(
    sheet: str, item: dict[str, Any], event_time: Any,
    clinical_rows: list[dict[str, Any]], ledger_path: Path,
    ledger_sha256: str,
) -> dict[str, Any]:
    row = item.get("row") or {}
    date = record_date(sheet, row)
    overlap = overlap_ids(record_text(sheet, row), clinical_rows)
    if sheet == "Lab_Results":
        overlap = lab_overlap_ids(item, clinical_rows)
    output = {
        "source_id": item.get("id"),
        "source_ref": {
            "ledger_file": str(ledger_path.resolve()),
            "ledger_sha256": ledger_sha256,
            "sheet": sheet,
            "source_row": row.get("_source_row"),
            "source_page": row.get("evidence_page") or row.get("pdf_page"),
        },
        "date": date,
        "date_basis": "lab_report_or_parser_date_unverified" if sheet == "Lab_Results" else "source_record_date",
        "event_window_relation": episode_window_relation(date, event_time),
        "source_type": row.get("source_type"),
        "possible_same_kh_evidence_ids": overlap,
        "independent_corroboration": False,
        "clinical_role": "phenotype_derived_context_not_independent_microbiology",
        "source_qc": "workbook_system_review_required_source_pdf_unavailable",
    }
    if sheet == "Phenotypes_Long":
        output.update(phenotype=row.get("phenotype"), status=row.get("status"),
                      confidence=row.get("confidence"), qc=row.get("phenotype_qc_status"),
                      temporal_relation_to_pneumonia=row.get("temporal_relation"),
                      reason=row.get("reason"), source_text=row.get("evidence"))
    elif sheet == "Evidence":
        output.update(phenotype=row.get("phenotype"), candidate_status=row.get("candidate_status"),
                      confidence=row.get("confidence"), qc="patient_system_review_required",
                      temporal_relation_to_pneumonia=row.get("temporal_relation"),
                      source_text=row.get("evidence"))
    elif sheet == "Clinical_Facts":
        output.update(fact_type=row.get("fact_type"), description=row.get("description"),
                      temporal_relation_to_pneumonia=row.get("temporal_relation"),
                      source_text=row.get("evidence"))
    elif sheet == "Clinical_Timeline":
        output.update(event_type=row.get("event_type"), description=row.get("description"),
                      source_text=row.get("evidence"))
    elif sheet == "Microbiology":
        output.update(organism=row.get("organism"), specimen=row.get("specimen"),
                      method=row.get("method"), resistance_marker=row.get("resistance_marker"),
                      source_text=row.get("evidence"),
                      clinical_role="physician_workbook_microbiology_requires_source_validation")
    elif sheet == "Lab_Results":
        output.update(test_name=row.get("test_name"), value=row.get("value"), unit=row.get("unit"),
                      parser=row.get("parser"), confidence=row.get("confidence"),
                      clinical_role="parser_lab_audit_only_not_independent_support")
    return output


def build_packet(
    timeline: dict[str, Any], ledger: dict[str, Any] | None,
    link: dict[str, Any], candidate_shadow: dict[str, Any],
    *, timeline_path: Path, ledger_path: Path | None, candidate_path: Path,
) -> dict[str, Any]:
    patient = int(timeline["patient_id"])
    if int(candidate_shadow["patient_id"]) != patient:
        raise ValueError("Candidate shadow belongs to a different patient")
    if int(link["pipeline_patient_number"]) != patient:
        raise ValueError("Linkage audit belongs to a different patient")
    if ledger is not None and patient_number(ledger.get("patient_id")) != patient:
        raise ValueError("Provisional same-number ledger does not match this patient")
    if ledger is None and ledger_path is not None:
        raise ValueError("Ledger path supplied without ledger")
    if ledger is not None and ledger_path is None:
        raise ValueError("Ledger missing source path")
    if not timeline.get("index_events"):
        raise ValueError("Timeline lacks index event")
    if candidate_shadow.get("linkage", {}).get("production_link_verified"):
        raise ValueError("Provisional packet must not assert a verified link")

    events = []
    clinical_rows = timeline.get("clinical_evidence") or []
    candidates = candidate_shadow.get("candidate_contexts") or []
    ledger_sha256 = sha256_file(ledger_path) if ledger_path else None
    for event in timeline["index_events"]:
        event_id = event["event_id"]
        event_time = event.get("collected_time")
        records: dict[str, list[dict[str, Any]]] = {sheet: [] for sheet in SHEETS}
        lab_total = 0
        if ledger is not None:
            for sheet in SHEETS:
                for item in ledger.get("records", {}).get(sheet, []):
                    if sheet == "Lab_Results":
                        lab_total += 1
                    compact = compact_record(
                        sheet,
                        item,
                        event_time,
                        clinical_rows,
                        ledger_path,
                        ledger_sha256,
                    )
                    if sheet != "Lab_Results" or compact["event_window_relation"] in {"within_window", "boundary_uncertain"}:
                        records[sheet].append(compact)
        by_id = {row["source_id"]: row for sheet in ("Phenotypes_Long", "Evidence") for row in records[sheet]}
        candidate_packets = []
        for candidate in candidates:
            candidate_event_refs = {
                str(value).strip()
                for value in candidate.get("candidate_event_refs") or []
                if str(value).strip()
            }
            event_case = str(event.get("case_review_id") or event.get("specimen_code") or "").strip()
            if candidate_event_refs and event_case not in candidate_event_refs:
                continue
            expanded = []
            for item in (candidate.get("expanded_phenotypes") or []) if ledger is not None else []:
                phenotype_row = by_id.get(item.get("source_id"))
                if phenotype_row is None:
                    raise ValueError(f"Missing routed phenotype source {item.get('source_id')}")
                evidence_refs = [by_id[record.get("source_id")] for record in item.get("corroborating_evidence") or []
                                 if record.get("source_id") in by_id]
                expanded.append({**item, "event_alignment": phenotype_row,
                                 "source_evidence_alignment": evidence_refs,
                                 "independent_corroboration": False})
            exact_micro = [row for row in records["Microbiology"]
                           if canonical_key(row.get("organism")) == candidate.get("canonical_key")
                           and row["event_window_relation"] in {"within_window", "boundary_uncertain"}]
            candidate_packets.append({
                "organism_name": candidate.get("organism_name"),
                "canonical_key": candidate.get("canonical_key"),
                "taxonomy_profile": candidate.get("taxonomy_profile"),
                "candidate_event_refs": sorted(candidate_event_refs),
                "event_scope": candidate.get("event_scope") or "patient_episode_all_events",
                "candidate_source": candidate.get("candidate_source"),
                "incoming_decision": candidate.get("incoming_decision"),
                "incoming_level": candidate.get("incoming_level"),
                "phenotype_deterministic_effect": "none_provisional_link",
                "phenotype_context": expanded,
                "exact_organism_workbook_microbiology": exact_micro,
                "model_use": "provisional_shadow_review_only_event_unverified" if ledger is not None else "unavailable",
                "model_constraint": "phenotype context is not independent organism proof or a Picked trigger",
            })
        events.append({
            "index_event": event,
            "linkage": {
                "patient_identity_assessment": link.get("patient_identity_assessment"),
                "infection_episode_assessment": link.get("infection_episode_assessment"),
                "duplicate_export_cluster": link.get("phenotype_duplicate_export_cluster") or [],
                "production_link_verified": False,
                "external_confirmation_required": True,
            },
            "source_records": records,
            "lab_source_total_in_ledger": lab_total,
            "lab_packet_rows_temporally_retrieved": len(records["Lab_Results"]),
            "candidate_packets": candidate_packets,
            "candidate_count": len(candidate_packets),
        })
    output = {
        "schema_version": SCHEMA_VERSION,
        "patient_id": patient,
        "primary_analysis": "retrospective_infection_episode_window",
        "episode_window_hours_each_side": 48,
        "decision_effect": "none_shadow_only",
        "answer_blind": True,
        "source": {
            "timeline_file": str(timeline_path.resolve()),
            "timeline_sha256": sha256_file(timeline_path),
            "candidate_shadow_file": str(candidate_path.resolve()),
            "candidate_shadow_sha256": sha256_file(candidate_path),
            "phenotype_ledger_file": str(ledger_path.resolve()) if ledger_path else None,
            "phenotype_ledger_sha256": ledger_sha256,
            "source_workbook": (ledger or {}).get("source", {}).get("workbook_file"),
            "source_workbook_sha256": (ledger or {}).get("source", {}).get("workbook_sha256"),
            "source_pdf_available": (ledger or {}).get("source", {}).get("source_pdf_available", False),
            "patient_qc_status": (ledger or {}).get("source", {}).get("patient_qc_status"),
        },
        "phenotype_status": "provisional_link" if ledger is not None else "phenotype_file_missing",
        "all_phenotype_overview": candidate_shadow.get("all_phenotype_overview") or [],
        "events": events,
    }
    forbidden = find_forbidden_keys(output)
    if forbidden:
        raise ValueError(f"Benchmark answer keys found: {forbidden[:5]}")
    return output


def run(timeline_root: Path, phenotype_root: Path, audit_path: Path,
        candidate_root: Path, output_dir: Path) -> dict[str, Any]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_dir}")
    audit = read_json(audit_path)
    links = {int(item["pipeline_patient_number"]): item for item in audit["patients"]}
    packets = []
    for timeline_path in sorted(timeline_root.glob("NGS_patient_*_clinical_timeline_shadow.json")):
        timeline = read_json(timeline_path)
        patient = int(timeline["patient_id"])
        if patient not in links:
            raise ValueError(f"No linkage audit entry for P{patient}")
        candidate_path = candidate_root / f"patient_{patient}_candidate_phenotype_shadow_v1.json"
        if not candidate_path.is_file():
            raise FileNotFoundError(candidate_path)
        ledger_path = phenotype_root / f"patient_{patient}_phenotype_decision_evidence_v1.json"
        ledger = read_json(ledger_path) if ledger_path.is_file() else None
        packet = build_packet(timeline, ledger, links[patient], read_json(candidate_path),
                              timeline_path=timeline_path, ledger_path=ledger_path if ledger else None,
                              candidate_path=candidate_path)
        packets.append(packet)
    if not packets:
        raise ValueError("No clinical timeline packets found")
    output_dir.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    for packet in packets:
        counts[packet["phenotype_status"]] += 1
        for event in packet["events"]:
            counts["events"] += 1
            counts["candidates"] += event["candidate_count"]
            for sheet, rows in event["source_records"].items():
                counts[f"rows:{sheet}"] += len(rows)
                counts[f"window:{sheet}:within"] += sum(
                    row["event_window_relation"] == "within_window" for row in rows)
        path = output_dir / f"patient_{packet['patient_id']}_phenotype_event_shadow.json"
        path.write_text(json.dumps(packet, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    summary = {"schema_version": SCHEMA_VERSION, "patient_count": len(packets),
               "counts": dict(sorted(counts.items())), "decision_effect": "none_shadow_only",
               "linkage_audit_file": str(audit_path.resolve()),
               "linkage_audit_sha256": sha256_file(audit_path)}
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("timeline_root", type=Path)
    parser.add_argument("phenotype_root", type=Path)
    parser.add_argument("linkage_audit", type=Path)
    parser.add_argument("candidate_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.timeline_root, args.phenotype_root, args.linkage_audit,
                         args.candidate_root, args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

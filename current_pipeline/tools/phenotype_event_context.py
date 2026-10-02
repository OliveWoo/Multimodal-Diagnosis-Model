"""Compact candidate-specific view of an audited phenotype event packet."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from tools.pathogen_normalization import canonical_key


def read_patient_packet(root: Path, patient_id: Any) -> dict[str, Any]:
    number = int(str(patient_id).removeprefix("P"))
    path = root / f"patient_{number}_phenotype_event_shadow.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    packet = json.loads(path.read_text(encoding="utf-8-sig"))
    if int(packet["patient_id"]) != number:
        raise ValueError(f"Phenotype event packet patient mismatch: {path}")
    if packet.get("decision_effect") != "none_shadow_only":
        raise ValueError("Expected a shadow-only phenotype event packet")
    return packet


def candidate_event_view(packet: dict[str, Any], organism_name: Any) -> dict[str, Any]:
    key = canonical_key(organism_name)
    views = []
    for event in packet.get("events") or []:
        candidate = next((item for item in event.get("candidate_packets") or []
                          if item.get("canonical_key") == key), None)
        if candidate is None:
            continue
        records = event.get("source_records") or {}
        facts = [item for item in records.get("Clinical_Facts") or []
                 if item.get("event_window_relation") in {"within_window", "boundary_uncertain"}
                 or item.get("temporal_relation_to_pneumonia") == "PRE_PNEUMONIA"]
        timeline = [item for item in records.get("Clinical_Timeline") or []
                    if item.get("event_window_relation") in {"within_window", "boundary_uncertain"}]
        views.append({
            "index_event": event.get("index_event"),
            "linkage": event.get("linkage"),
            "phenotype_context": candidate.get("phenotype_context") or [],
            "exact_organism_workbook_microbiology": candidate.get("exact_organism_workbook_microbiology") or [],
            "dated_clinical_facts": facts,
            "dated_clinical_timeline": timeline,
            "lab_parser_rows_in_packet": event.get("lab_packet_rows_temporally_retrieved"),
            "lab_values_supplied_to_model": False,
            "model_use": candidate.get("model_use"),
        })
    return {
        "schema_version": "phenotype_event_model_context_shadow.v1",
        "patient_id": packet.get("patient_id"),
        "organism_name": str(organism_name),
        "canonical_key": key,
        "phenotype_status": packet.get("phenotype_status"),
        "source": packet.get("source"),
        "event_views": views,
        "constraints": [
            "Retrospective mNGS collection +/-48h view; not an at-collection prediction.",
            "Patient/episode links are provisional, not production verified.",
            "Phenotype-derived claims are context, not independent organism confirmation.",
            "Workbook microbiology requires source-PDF, specimen, method, and patient/event validation.",
            "NO is not proven absence; uncertain time and QC remain uncertain.",
            "No change to deterministic Picked is permitted.",
        ],
    }

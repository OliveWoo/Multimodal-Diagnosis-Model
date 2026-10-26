"""Bind source facts to the viral three-axis shadow without benchmark access.

The adapter preserves all frozen candidates. It verifies source fields, creates
traceable detection/context claims, and carries unresolved clinical assertions
for review. It does not adjudicate viral pneumonia or modify production tiers.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from tools.pathogen_normalization import canonical_key
from tools.build_index_event_clinical_timeline_shadow import time_interval
from tools.clinical_event_window import PRIMARY_WINDOW_MODE
from tools.viral_evidence_axes_shadow import interpret, DEFAULT_POLICY as AXES_POLICY

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "outputs/runs/2026-10-19_viral_policy_consistency_review_v5"
DEFAULT_OUTPUT = ROOT / "outputs/runs/2026-10-19_viral_evidence_binding_shadow_v4"
POLICY_PATH = ROOT / "rules/viral_evidence_binding_shadow_v1.json"


def read(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, values: list[dict]) -> None:
    fields = list(values[0]) if values else ["patient_id"]
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(values)


def text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def scalar_text(row: dict, keys: tuple[str, ...]) -> str:
    return " ".join(text(row.get(k)) for k in keys if row.get(k) is not None)


def target_matches(value: str, policy: dict) -> dict[str, str]:
    """Use bounded assay/name terms, keeping group and species distinct."""
    found = {}
    for key, aliases in policy["organism_targets"].items():
        if any(re.search(r"(?<![a-z0-9])" + re.escape(a) + r"(?![a-z0-9])", value, re.I) for a in aliases):
            found[key] = "exact_or_alias"
    for alias, members in policy["group_targets"].items():
        if re.search(r"(?<![a-z0-9])" + re.escape(alias) + r"(?![a-z0-9-])", value, re.I):
            for key in members:
                found.setdefault(key, "group_only")
    # A combined HSV-1/2 target identifies the panel group, not either species.
    if re.search(r"HSV[- ]?1\s*[/&+]\s*(?:HSV[- ]?)?2\b|herpes simplex(?: virus)?\s*1\s*[/&+]\s*2", value, re.I):
        for key in ("humanalphaherpesvirus1", "humanalphaherpesvirus2"):
            found[key] = "group_only"
    return found


def align(collected: Any, index: Any, *, usable: bool = True) -> dict:
    observed, anchor = time_interval(collected), time_interval(index)
    if not usable:
        return {"status": "source_time_unusable", "signed_calendar_days": None, "signed_hours": None}
    if not observed or not anchor or observed[2] == "month" or anchor[2] == "month":
        return {"status": "unknown", "signed_calendar_days": None, "signed_hours": None}
    delta = (observed[0].date() - anchor[0].date()).days
    return {"status": "within_window" if abs(delta) <= 2 else "outside_window",
            "signed_calendar_days": delta,
            "signed_hours": (observed[0] - anchor[0]).total_seconds() / 3600,
            "window_mode": PRIMARY_WINDOW_MODE,
            "precision": observed[2]}


def specimen_context(value: Any) -> str:
    value = text(value).casefold()
    if re.search(r"\bblood\b|\bplasma\b|\bserum\b|血", value):
        return "blood"
    if re.search(r"\bbalf?\b|bronch|sputum|endo|trache|lung|肺|痰", value):
        return "respiratory"
    if re.search(r"urine|尿", value):
        return "urine"
    if re.search(r"lesion|vesicle|skin", value):
        return "skin_lesion"
    return "unknown"


def result_status(value: Any) -> str:
    value = re.sub(r"\s+", " ", text(value)).casefold()
    if not value:
        return "unknown"
    negative = re.search(r"not\s*detected|negative|nonreactive|no growth|not isolated|未檢出|陰性", value)
    residual = re.sub(r"not\s*detected|negative|nonreactive|no growth|not isolated|未檢出|陰性", " ", value)
    positive = re.search(r"\bdetected\b|\bpositive\b|\breactive\b|\bisolated\b", residual)
    if negative and positive:
        return "conflicting_or_composite_result"
    if negative:
        return "negative"
    if re.search(r"\bdetected\b|\bpositive\b|\breactive\b|\bisolated\b", value):
        return "positive"
    if re.fullmatch(r"[<>≤≥]?\s*[\d,.]+(?:[eE][+-]?\d+)?", value):
        return "below_limit" if value.startswith(("<", "≤")) else "numeric_measurement"
    return "unresolved"


def source_ref(path: Path, pointer: str, row: dict | None = None) -> dict:
    origin = (row or {}).get("source") or {}
    return {"file": str(path.resolve()), "sha256": sha(path), "json_path": pointer,
            "upstream_file": origin.get("source_file"), "upstream_sheet": origin.get("source_sheet"),
            "upstream_row": origin.get("row"), "upstream_page": (origin.get("original_row") or {}).get("source_page")}


def fact_id(run: str, patient: str, ref: dict) -> str:
    # Original source row identifies copied records when available; otherwise
    # the full JSON source pointer is the only established origin.
    identity = [run, patient, ref.get("upstream_file"), ref.get("upstream_sheet"), ref.get("upstream_row")]
    if not all(identity[2:]):
        identity = [run, patient, ref["file"], ref["json_path"]]
    digest = hashlib.sha256(json.dumps(identity, ensure_ascii=False).encode()).hexdigest()[:20]
    return f"{run}:P{patient}:VFACT:{digest}"


def row_facts(run: str, patient: str, module: str, row: dict, index: int,
              path: Path, policy: dict) -> list[dict]:
    facts = []
    ref = source_ref(path, f"[{index}]", row)
    if module in {"image", "admission_diagnosis", "underlying"}:
        if module == "image":
            parts = [(f"findings[{i}]", v) for i, v in enumerate(row.get("findings") or [])]
            parts += [(k, row[k]) for k in ("impression", "raw_report") if row.get(k)]
        elif module == "admission_diagnosis":
            parts = [("diagnosis", row.get("diagnosis"))]
            parts += [(f"notes[{i}]", v) for i, v in enumerate(row.get("notes") or [])]
        else:
            parts = [(f"underlying_diseases[{i}]", v) for i, v in enumerate(row.get("underlying_diseases") or [])]
        for field, value in parts:
            if not isinstance(value, str) or not value.strip():
                continue
            matches = target_matches(value, policy)
            pulmonary = bool(any(re.search(p, value, re.I) for p in policy["pulmonary_context_patterns"]))
            if not matches and not pulmonary:
                continue
            pointer = source_ref(path, f"[{index}].{field}", row)
            blockers = []
            # source_time has explicitly unspecified semantics in KH v14.
            collected = row.get("collected_time") or row.get("recorded_time")
            if row.get("usable_for_temporal_inference") is False:
                blockers.append("source_time_semantics_unverified")
            facts.append({"evidence_id": fact_id(run, patient, ref), "patient_id": patient,
                "module": module, "fact_kind": "narrative", "target_matches": matches,
                "pulmonary_context_mentioned": pulmonary, "source_text": value,
                "collected_time": collected, "reported_time": row.get("reported_time"),
                "unspecified_source_time": row.get("source_time"),
                "specimen": row.get("specimen_type") or row.get("sample"),
                "source_refs": [pointer], "blockers": blockers,
                "result_status": "narrative_only", "raw_record": row,
                "temporal_usable": row.get("usable_for_temporal_inference") is not False,
                "source_field_verified": True})
        return facts
    label = scalar_text(row, ("target", "organism_name", "organism", "item", "item_full", "test"))
    matches = target_matches(label, policy)
    if not matches:
        return []
    result = row.get("result", row.get("value", row.get("raw_result")))
    source = row.get("source") or {}
    normalization = (source.get("field_normalization") or {}).get("fields") or {}
    unit = row.get("unit") or normalization.get("unit")
    # Recover only a literal unit in the same original result, not a neighboring
    # composite raw-evidence row and never an assumed unit for mNGS reads.
    original_result = text(source.get("raw_result"))
    if not unit:
        match = re.search(r"\b(IU/mL|copies/mL|copy/mL|copies/µL|copies/uL)\b", original_result, re.I)
        if match:
            unit = match.group(1)
    serology = bool(re.search(r"\bIg[GM]\b|antibod|serolog|S/CO|AU/mL", label + " " + text(unit), re.I))
    assay_method = row.get("assay_method") or row.get("method")
    if not assay_method:
        method = re.search(r"\b(qPCR|RT-PCR|PCR|GeneXpert|IHC|immunohistochemistry)\b", label, re.I)
        assay_method = method.group(1) if method else None
    blockers = []
    for field in ("usable_for_result_inference", "usable_for_value_inference"):
        if row.get(field) is False:
            blockers.append(field + "_false")
    if row.get("review_required") is True or row.get("result_status") == "unresolved":
        blockers.append("source_record_requires_review")
    kind = "serology" if serology else ("molecular" if module == "molecular_microbiology" or assay_method else module)
    numeric_input = next((row[k] for k in ("numeric_value", "viral_load_or_copies", "value")
                          if row.get(k) not in (None, "")), result)
    numeric_text = text(numeric_input)
    number = re.fullmatch(r"\s*([<>≤≥]?)\s*([\d,.]+(?:[eE][+-]?\d+)?)\s*", numeric_text)
    numeric = float(number.group(2).replace(",", "")) if number else None
    facts.append({"evidence_id": fact_id(run, patient, ref), "patient_id": patient,
        "module": module, "fact_kind": kind, "target_matches": matches,
        "test": row.get("test") or row.get("item_full") or row.get("item"),
        "target": row.get("target") or row.get("organism_name") or row.get("organism"),
        "assay_method": assay_method, "specimen": row.get("sample") or row.get("specimen_type") or row.get("specimen"),
        "collected_time": row.get("collected_time"), "reported_time": row.get("reported_time"),
        "result": result, "result_status": result_status(result), "numeric_value": numeric,
        "numeric_comparator": number.group(1) if number else None, "unit": unit,
        "numeric_source_field": next((k for k in ("numeric_value", "viral_load_or_copies", "value")
                                      if row.get(k) not in (None, "")), "result"),
        "ct_value": row.get("ct_value"), "reference_or_cutoff": row.get("reference_or_cutoff"),
        "temporal_usable": row.get("usable_for_temporal_inference") is not False,
        "source_refs": [ref], "source_text": source.get("raw_evidence") or text(result),
        "blockers": blockers, "raw_record": row, "source_field_verified": not blockers})
    return facts


def choose_event(candidate: dict) -> tuple[dict | None, str]:
    events = candidate.get("linked_events") or []
    selected = (candidate.get("source_candidate") or {}).get("case_review_id")
    if selected:
        matches = [e for e in events if e.get("case_review_id") == selected]
        if len(matches) == 1:
            return matches[0], "selected_candidate_case"
    if len(events) == 1:
        return events[0], "single_linked_event"
    return None, "missing_or_ambiguous_selected_event"


def current_activity_assertion(value: str, key: str, policy: dict) -> bool:
    """A dated explicit clinician assertion supports A1, not proven activity."""
    for clause in re.split(r"[;\n。]", value):
        if target_matches(clause, policy).get(key) != "exact_or_alias":
            continue
        aliases = policy["organism_targets"].get(key, [])
        named = r"(?<![a-z0-9])(?:" + "|".join(re.escape(a) for a in aliases) + r")(?![a-z0-9])"
        activity = r"(?:reactivat(?:ion|ed)|再活化)"
        connector = r"\s*(?:\([^)]{0,50}\)\s*)?(?:(?:viral|virus|current|currently|infection|has|is|with)\s+){0,2}"
        if not (re.search(named + connector + activity, clause, re.I)
                or re.search(activity + r"\s+(?:of\s+)?" + named, clause, re.I)):
            continue
        if re.search(r"\b(?:no|not|suspect(?:ed)?|possible|prior|previous|history|resolved|exclude|risk)\b|排除|疑似|曾經|已緩解", clause, re.I):
            continue
        return True
    return False


def deduplicate_claims(claims: list[dict]) -> list[dict]:
    unique = {}
    for claim in claims:
        key = (claim["evidence_id"], claim["claim_type"], claim["organism_key"], claim["episode_id"])
        if key not in unique:
            unique[key] = copy.deepcopy(claim)
            unique[key]["source_texts"] = [claim["source_text"]] if claim.get("source_text") else []
            continue
        existing = unique[key]
        for ref in claim.get("source_refs") or []:
            if ref not in existing.setdefault("source_refs", []):
                existing["source_refs"].append(ref)
        if claim.get("source_text") and claim["source_text"] not in existing["source_texts"]:
            existing["source_texts"].append(claim["source_text"])
    return list(unique.values())


def reconcile_frozen_observations(candidate: dict, facts: list[dict], policy: dict) -> dict:
    """Audit that previously accepted hospital rows are still represented."""
    detail = candidate.get("hospital_evidence_detail") or (
        candidate.get("hospital_profile") or {}).get("hospital_evidence_detail") or {}
    observations = detail.get("source_observations") or []
    resolution = []
    for observation in observations:
        targets = target_matches(text(observation.get("target") or detail.get("organism_name")), policy)
        matches = []
        for fact in facts:
            if not set(targets).intersection(fact.get("target_matches") or {}):
                continue
            left, right = time_interval(observation.get("collected_time")), time_interval(fact.get("collected_time"))
            if not left or not right or left[0] != right[0]:
                continue
            if text(observation.get("specimen_type")).casefold() != text(fact.get("specimen")).casefold():
                continue
            if result_status(observation.get("raw_result")) != fact.get("result_status"):
                continue
            matches.append(fact["evidence_id"])
        resolution.append({"observation_id": observation.get("observation_id"),
            "status": "resolved_to_source_fact" if matches else "unresolved_preserved_in_frozen_profile",
            "source_fact_ids": sorted(set(matches)), "frozen_observation": observation})
    return {"requested": len(observations), "resolved": sum(bool(r["source_fact_ids"]) for r in resolution),
            "missing": [r["observation_id"] for r in resolution if not r["source_fact_ids"]],
            "observations": resolution}


def claim_from_fact(fact: dict, key: str, episode_id: str, event: dict) -> tuple[dict | None, list[str]]:
    temporal = align(fact.get("collected_time"), event.get("collected_time"), usable=fact.get("temporal_usable", True))
    blockers = list(fact.get("blockers") or [])
    if fact.get("case_review_id") and fact["case_review_id"] != event.get("case_review_id"):
        blockers.append("different_mngs_case")
    if temporal["status"] != "within_window":
        blockers.append(temporal["status"])
    match = (fact.get("target_matches") or {}).get(key)
    if match == "group_only":
        blockers.append("identity_group_not_exact_species")
    is_narrative = fact["fact_kind"] in {"narrative", "index_symptom"}
    if is_narrative:
        # Generic pulmonary context does not constitute a virus-specific
        # diagnosis; historical/current viral assertions are carried for review.
        if fact.get("explicit_current_activity") and match == "exact_or_alias":
            claim_type = "clinician_current_activity"
        elif not fact.get("pulmonary_context_mentioned"):
            return None, blockers + ["viral_narrative_requires_clinical_attribution_review"]
        else:
            value = fact.get("source_text", "")
            if re.search(r"\b(?:no|without|ruled out|negative for)\b.{0,30}\b(?:pneumonia|infiltrat|consolidation)", value, re.I):
                blockers.append("negated_pulmonary_context")
            claim_type = "compatible_pulmonary_syndrome"
    else:
        if not match:
            return None, blockers + ["different_organism"]
        if fact["fact_kind"] == "serology":
            return None, blockers + ["serology_not_active_viral_detection"]
        quantitative_positive = (
            fact.get("result_status") == "numeric_measurement"
            and (fact.get("numeric_value") or 0) > 0
            and not fact.get("numeric_comparator")
            and bool(re.fullmatch(r"IU/mL|copies/mL|copy/mL|copies/µL|copies/uL", text(fact.get("unit")), re.I))
            and bool(re.search(r"\b(?:qPCR|PCR|viral load)\b|CMV-VL", text(fact.get("assay_method")) + " " + text(fact.get("target")), re.I))
        )
        if fact.get("result_status") != "positive" and not quantitative_positive:
            # Numbers, negatives and below-quantification results are preserved
            # as facts; their clinical interpretation is not guessed here.
            return None, blockers + ["result_is_" + fact.get("result_status", "unknown")]
        if not fact.get("specimen"):
            blockers.append("specimen_missing")
        claim_type = "detection"
    if blockers:
        return None, sorted(set(blockers))
    ref = fact["source_refs"][0]
    assay = text(fact.get("assay_method")).casefold()
    compartment = specimen_context(fact.get("specimen"))
    mngs_compartment = specimen_context(event.get("specimen_site"))
    if not is_narrative and compartment != "unknown" and mngs_compartment != "unknown" and compartment != mngs_compartment:
        independence = "independent_biological_specimen"
        independence_basis = "explicit_different_specimen_compartments"
    elif not is_narrative and re.search(r"pcr|genexpert|ihc|culture", assay):
        independence = "independent_clinical_assay"
        independence_basis = "explicit_assay_method_distinct_from_mngs"
    else:
        independence = "not_established"
        independence_basis = "method_or_biological_independence_not_documented"
    return {"claim_type": claim_type, "evidence_id": fact["evidence_id"], "organism_key": key,
        "patient_id": fact["patient_id"], "episode_id": episode_id,
        "specimen": fact.get("specimen") if not is_narrative else "not_applicable_clinical_context",
        "index_specimen": event.get("specimen_site"),
        "collected_time": fact["collected_time"], "source_file": ref["file"], "source_row": ref["json_path"],
        "source_sha256": ref["sha256"], "source_refs": fact["source_refs"],
        "review_complete": True, "verification_kind": "deterministic_source_field_validation",
        "human_clinical_adjudication": False, "event_aligned": True, "event_alignment": temporal,
        "identity_match": match or "patient_episode_context_not_pathogen_specific",
        "pathogen_specific": claim_type != "compatible_pulmonary_syndrome", "assay_method": fact.get("assay_method"),
        "independence_class": independence, "independence_basis": independence_basis,
        "biological_evidence_group": fact["evidence_id"],
        "source_text": fact.get("source_text"), "result": fact.get("result"),
        "numeric_value": fact.get("numeric_value"), "unit": fact.get("unit")}, []


def bind_candidate(candidate: dict, facts: list[dict], *, run: str,
                   provenance: dict, policy: dict, axes_policy: dict) -> dict:
    patient, key = str(candidate["patient_id"]), canonical_key(candidate["organism_name"])
    event, selection = choose_event(candidate)
    episode = f"{run}:{event['event_id']}" if event else f"{run}:P{patient}:episode_unresolved"
    claims, links = [], []
    if event:
        signals = [s for s in (candidate.get("analytical_profile") or {}).get("per_test_signals") or []
                   if float(s.get("reads") or 0) > 0 and s.get("qc_status") in axes_policy["detection_qc_statuses"]]
        # All branches of the selected case become one biological detection fact.
        if signals and event.get("collected_time") and event.get("specimen_site"):
            claims.append({"claim_type": "detection", "evidence_id": f"{episode}:{key}:mngs",
                "patient_id": patient, "organism_key": key, "episode_id": episode,
                "specimen": event["specimen_site"], "collected_time": event["collected_time"],
                "source_file": provenance["file"], "source_row": provenance["json_path"] + ".analytical_profile.per_test_signals",
                "source_sha256": provenance["sha256"], "review_complete": True, "event_aligned": True,
                "verification_kind": "deterministic_source_field_validation", "human_clinical_adjudication": False,
                "independence_class": "technical_branches_one_specimen", "biological_evidence_group": episode + ":mngs_specimen",
                "analytical_branches": signals, "reads_sum": None, "assay_method": "mNGS"})
    for fact in facts:
        matching = key in (fact.get("target_matches") or {})
        generic_context = fact.get("pulmonary_context_mentioned", False)
        if not matching and not generic_context:
            continue
        fact_view = dict(fact)
        fact_view["explicit_current_activity"] = current_activity_assertion(fact.get("source_text", ""), key, policy) if fact["fact_kind"] == "narrative" else False
        claim, blockers = claim_from_fact(fact_view, key, episode, event) if event else (None, [selection])
        if claim:
            claims.append(claim)
        links.append({"evidence_id": fact["evidence_id"], "module": fact["module"], "fact_kind": fact["fact_kind"],
            "identity_match": (fact.get("target_matches") or {}).get(key) or "patient_episode_context",
            "accepted_claim_type": claim["claim_type"] if claim else None,
            "binding_status": "connected_claim" if claim else "connected_context_or_audit",
            "blockers": blockers,
            "event_alignment": align(fact.get("collected_time"), event.get("collected_time") if event else None,
                                     usable=fact.get("temporal_usable", True)),
            "specimen": fact.get("specimen"), "collected_time": fact.get("collected_time"),
            "reported_time": fact.get("reported_time"), "assay_method": fact.get("assay_method"),
            "result": fact.get("result"), "result_status": fact.get("result_status"),
            "numeric_value": fact.get("numeric_value"), "unit": fact.get("unit"),
            "source_text": fact.get("source_text"), "source_refs": fact["source_refs"]})
    claims = deduplicate_claims(claims)
    result = interpret(candidate, episode_id=episode, claims=claims, policy=axes_policy)
    # Mapping is now complete for the supplied source modules. Remaining unknown
    # axes express documented evidence gaps, rather than unexecuted extraction.
    result["source_mapping_status"] = "completed_for_available_sources"
    result["gaps"] = [g for g in result["gaps"] if g != "clinical_claim_extraction_or_source_verification_incomplete"]
    result["gaps"].extend(sorted({b for link in links for b in link["blockers"]
                                  if b in {"specimen_missing", "source_record_requires_review", "unknown", "source_time_unusable"}}))
    if any(c["claim_type"] == "detection" and c.get("assay_method") is None for c in claims):
        result["gaps"].append("hospital_assay_method_not_documented")
    reconciliation = reconcile_frozen_observations(candidate, facts, policy)
    return {"patient_id": patient, "organism_name": candidate["organism_name"], "organism_key": key,
        "episode_id": episode, "selected_event": event, "event_selection_basis": selection,
        "claims": claims, "evidence_links": links, "axes": result,
        "current_tier": candidate.get("final_reporting_tier"),
        "current_selected": bool(candidate.get("selected_for_complete_report")),
        "clinical_adjudication_complete": False, "production_changed": False,
        "frozen_hospital_observation_reconciliation": reconciliation,
        "host_context_preserved": candidate.get("host_evidence") or [],
        "prior_history_consumption_preserved": candidate.get("history_evidence_consumption") or [],
        "phenotype_context_preserved": candidate.get("phenotype_evidence") or {},
        "source_candidate": provenance}


def run(input_root: Path, output_root: Path, *, cohort_runs: dict[str, str] | None = None) -> dict:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output_root}")
    source_packets_path = input_root / "review/all_affected_evidence_packets.json"
    packets = read(source_packets_path)
    policy, axes_policy = read(POLICY_PATH), read(AXES_POLICY)
    output_root.mkdir(parents=True)
    shutil.copy2(POLICY_PATH, output_root / POLICY_PATH.name)
    shutil.copy2(AXES_POLICY, output_root / AXES_POLICY.name)
    for code_path in (Path(__file__), ROOT / "tools/viral_evidence_axes_shadow.py"):
        code_snapshot = output_root / "code_snapshot" / code_path.name
        code_snapshot.parent.mkdir(exist_ok=True)
        shutil.copy2(code_path, code_snapshot)
    shutil.copy2(source_packets_path, output_root / "frozen_affected_evidence_packets.json")
    packets = read(output_root / "frozen_affected_evidence_packets.json")
    source_registry = [{"source": str(source_packets_path.resolve()), "sha256": sha(source_packets_path),
        "snapshot": str((output_root / "frozen_affected_evidence_packets.json").resolve()),
        "snapshot_sha256": sha(output_root / "frozen_affected_evidence_packets.json"), "status": "present"}]
    facts_by_patient, inventories = {}, {}
    all_facts, all_results, link_rows = [], [], []
    run_names = cohort_runs or {}
    for packet in packets:
        run_name, patient = run_names.get(packet["cohort"], packet["cohort"]), str(packet["patient_id"])
        if (run_name, patient) in facts_by_patient:
            continue
        manifest_path = input_root / run_name / "run_manifest.json"
        manifest = read(manifest_path)
        source_registry.append({"source": str(manifest_path.resolve()), "sha256": sha(manifest_path)})
        patient_root = Path(manifest["input"]["patient_root"])
        folder = patient_root / f"NGS_patient_{patient}_json"
        facts = []
        for module in policy["source_modules"]:
            path = folder / f"NGS_patient_{patient}_{module}.json"
            if not path.is_file():
                source_registry.append({"patient_id": patient, "run": run_name, "module": module, "status": "missing"})
                continue
            snapshot = output_root / "00_source_snapshot" / run_name / patient / path.name
            snapshot.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, snapshot)
            source_registry.append({"patient_id": patient, "run": run_name, "module": module,
                "source": str(path.resolve()), "sha256": sha(path), "snapshot": str(snapshot.resolve()),
                "snapshot_sha256": sha(snapshot), "status": "present"})
            payload = read(snapshot)
            if not isinstance(payload, list):
                raise TypeError(f"Unsupported source module: {path}")
            for i, row in enumerate(payload):
                if isinstance(row, dict):
                    facts.extend(row_facts(run_name, patient, module, row, i, snapshot, policy))
        inv = Path(manifest["input"]["inventory_root"]) / "patients" / f"NGS_patient_{patient}_multi_assay_mngs_evidence.json"
        snapshot = output_root / "00_source_snapshot" / run_name / patient / inv.name
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(inv, snapshot)
        source_registry.append({"patient_id": patient, "run": run_name, "module": "mngs_inventory",
            "source": str(inv.resolve()), "sha256": sha(inv), "snapshot": str(snapshot.resolve()), "snapshot_sha256": sha(snapshot), "status": "present"})
        inventory = read(snapshot)
        inventories[(run_name, patient)] = inventory
        # Only the recorded symptom field is consumed; expected-pathogen and
        # outcome/clinical-result fields are deliberately not inspected.
        for i, case in enumerate(inventory.get("cases") or []):
            symptom = text((case.get("case_review") or {}).get("symptom_description"))
            if not symptom or not any(re.search(p, symptom, re.I) for p in policy["pulmonary_context_patterns"]):
                continue
            ref = source_ref(snapshot, f"cases[{i}].case_review.symptom_description")
            facts.append({"evidence_id": fact_id(run_name, patient, ref), "patient_id": patient,
                "module": "mngs_case_symptom", "fact_kind": "index_symptom", "target_matches": {},
                "case_review_id": case.get("case_review_id"), "pulmonary_context_mentioned": True,
                "source_text": symptom, "collected_time": case.get("collected_time"),
                "reported_time": None, "specimen": case.get("specimen_site"),
                "source_refs": [ref], "blockers": [], "temporal_usable": True,
                "source_field_verified": True, "result_status": "context_only"})
        facts_by_patient[(run_name, patient)] = facts
        all_facts.extend(dict(f, run=run_name) for f in facts)
    for packet in packets:
        candidate = packet["current_candidate_evidence_unchanged"]
        run_name, patient = run_names.get(packet["cohort"], packet["cohort"]), str(packet["patient_id"])
        # The frozen copy is the authority for mNGS claims; no mutable original
        # packet is silently read as a replacement.
        provenance = {"file": str((output_root / "frozen_affected_evidence_packets.json").resolve()),
            "sha256": sha(output_root / "frozen_affected_evidence_packets.json"),
            "json_path": f"[{len(all_results)}].current_candidate_evidence_unchanged"}
        result = bind_candidate(candidate, facts_by_patient[(run_name, patient)], run=run_name,
                                provenance=provenance, policy=policy, axes_policy=axes_policy)
        result["cohort"] = packet["cohort"]
        all_results.append(result)
        for link in result["evidence_links"]:
            link_rows.append({"cohort": packet["cohort"], "patient_id": patient,
                "organism_name": candidate["organism_name"], **{k: v for k, v in link.items()
                    if k not in {"source_refs", "event_alignment"}},
                "blockers": "; ".join(link["blockers"]),
                "window_status": link["event_alignment"]["status"],
                "source_file": link["source_refs"][0]["file"],
                "source_row": link["source_refs"][0]["json_path"]})
    before_keys = {(p["cohort"], p["patient_id"], p["organism_name"]) for p in packets}
    after_keys = {(p["cohort"], p["patient_id"], p["organism_name"]) for p in all_results}
    if before_keys != after_keys or len(all_results) != len(packets):
        raise ValueError("Candidate preservation failed")
    source_issues = [r["source"] for r in source_registry if r.get("sha256") and (
        sha(Path(r["source"])) != r["sha256"]
        or r.get("snapshot_sha256", r["sha256"]) != r["sha256"])]
    if source_issues:
        raise ValueError(f"Sources changed during binding: {source_issues}")
    summaries = []
    for p in all_results:
        a = p["axes"]
        summaries.append({"cohort": p["cohort"], "patient_id": p["patient_id"], "organism_name": p["organism_name"],
            "episode_id": p["episode_id"], "current_tier": p["current_tier"], "current_selected": p["current_selected"],
            "detection": a["detection"], "activity": a["activity"], "organ_support": a["organ_support"],
            "shadow_tier": a["shadow_tier"], "accepted_claim_count": a["accepted_claim_count"],
            "clinical_claims": sum(c["claim_type"] != "detection" for c in p["claims"]),
            "hospital_positive_claims": sum(c["claim_type"] == "detection" and c.get("assay_method") != "mNGS" for c in p["claims"]),
            "linked_source_count": len(p["evidence_links"]), "gaps": "; ".join(a["gaps"])})
    write(output_root / "source_registry.json", source_registry)
    write(output_root / "source_facts.json", all_facts)
    write(output_root / "bound_candidate_packets.json", all_results)
    write_csv(output_root / "candidate_binding_summary.csv", summaries)
    write_csv(output_root / "source_binding_audit.csv", link_rows)
    linked_ids = {link["evidence_id"] for p in all_results for link in p["evidence_links"]}
    unassigned = [f for f in all_facts if f["target_matches"] and f["evidence_id"] not in linked_ids]
    write(output_root / "viral_facts_without_matching_candidate.json", unassigned)
    summary = {"schema_version": "viral_evidence_binding_shadow.v1", "answer_blind": True,
        "benchmark_answers_loaded": False, "production_decisions_changed": False,
        "candidate_count": len(all_results), "candidate_loss": len(before_keys - after_keys),
        "patient_count": len(facts_by_patient), "source_fact_count": len(all_facts),
        "unique_source_evidence_group_count": len({f["evidence_id"] for f in all_facts}),
        "viral_facts_without_existing_candidate": len(unassigned),
        "cohort_counts": dict(Counter(p["cohort"] for p in all_results)),
        "detection_counts": dict(Counter(p["axes"]["detection"] for p in all_results)),
        "activity_counts": dict(Counter(p["axes"]["activity"] for p in all_results)),
        "organ_support_counts": dict(Counter(p["axes"]["organ_support"] for p in all_results)),
        "candidates_with_connected_hospital_detection": sum(r["hospital_positive_claims"] > 0 for r in summaries),
        "candidates_with_pulmonary_context": sum(p["axes"]["organ_support"] == "O1" for p in all_results),
        "candidates_with_no_relevant_non_mngs_fact": sum(not p["evidence_links"] for p in all_results),
        "rejected_claim_count": sum(len(p["axes"]["rejected_claims"]) for p in all_results),
        "frozen_hospital_observations_requested": sum(p["frozen_hospital_observation_reconciliation"]["requested"] for p in all_results),
        "frozen_hospital_observations_resolved": sum(p["frozen_hospital_observation_reconciliation"]["resolved"] for p in all_results),
        "frozen_hospital_observations_unresolved": sum(len(p["frozen_hospital_observation_reconciliation"]["missing"]) for p in all_results),
        "source_mapping_completed": True, "clinical_adjudication_complete": False,
        "binding_policy_sha256": sha(output_root / POLICY_PATH.name),
        "axes_policy_sha256": sha(output_root / AXES_POLICY.name),
        "evidence_packets_sha256": sha(output_root / "frozen_affected_evidence_packets.json"),
        "adapter_code_sha256": sha(Path(__file__)),
        "axes_code_sha256": sha(ROOT / "tools/viral_evidence_axes_shadow.py")}
    write(output_root / "summary.json", summary)
    lines = ["# 病毒三軸證據接線結果", "", "本輪完成可用來源的接線與 shadow replay；沒有讀取 benchmark 答案，也沒有改正式分級。", "",
        f"共 {summary['candidate_count']} 筆候選、{summary['patient_count']} 位病人；候選遺失 {summary['candidate_loss']} 筆。",
        f"院端陽性檢出接入 {summary['candidates_with_connected_hospital_detection']} 筆候選；相符事件肺部背景接入 {summary['candidates_with_pulmonary_context']} 筆。",
        f"檢出狀態：{summary['detection_counts']}；活動狀態：{summary['activity_counts']}；肺部支持：{summary['organ_support_counts']}。",
        f"原有院端觀察 {summary['frozen_hospital_observations_requested']} 列，已核回原始輸入 {summary['frozen_hospital_observations_resolved']} 列，未解決 {summary['frozen_hospital_observations_unresolved']} 列。", "",
        "## 重要病例與來源限制", "",
        "| 病人 | 接入的證據 | 時間／限制 | 本輪解讀 |", "|---|---|---|---|"]
    for p in all_results:
        if p["organism_key"] != "humanbetaherpesvirus5":
            continue
        clinical_links = [l for l in p["evidence_links"] if l["fact_kind"] == "molecular"]
        if not clinical_links:
            continue
        for link in clinical_links:
            raw = text(link.get("result"))
            if link.get("unit"):
                raw += " " + link["unit"]
            site = text(link.get("specimen")) or "檢體未記載"
            delta = link["event_alignment"].get("signed_calendar_days")
            when = text(link.get("collected_time"))
            meaning = "接入檢出支持；不直接判定活動或肺炎" if link["accepted_claim_type"] else "僅保留來源與限制，未成為本次事件的陽性支持"
            lines.append(f"| {p['cohort']} P{p['patient_id']} | CMV：{site}；{raw} | {when}；相差{delta}曆日；{'; '.join(link['blockers']) or '事件內，方法未記載' } | {meaning} |")
    lines += ["", "## 如何解讀", "",
        "- 已接線表示每筆來源、檢體、時間、結果與限制都進入可追溯資料結構；不代表已由醫師確認病毒肺炎。",
        "- D2 只代表檢出有可確認的另一檢驗或生物檢體支持；血液或呼吸道陽性本身不證明肺炎。",
        "- O1 只代表同次事件具有相符肺部背景，未建立病毒特異歸因；沒有直接升成 Possible／Picked。",
        "- A0 是目前已接入資料仍未建立病毒活動／再活化支持，不是陰性。",
        "- O1的統計單位是病例－病毒候選，並非獨立病人或獨立證據；來源片段與重複副本已按evidence_id去重。",
        "- 本輪不自動產生probable/proven病毒肺炎；原型shadow_tier不能當作已定案的降級建議。正式Picked／Possible／Fallback與既有統計維持原值。",
        "- source_time 未指定採檢語意、未知日期、未定種、陰性、抗體、窗口外及 parser 未釐清資料都有獨立狀態，不會消失。",
        "- 主要事件窗口沿用前兩天至後兩天的曆日定義；報告時間不替代採檢時間，另保留精確小時差供稽核。", "",
        "## 每筆候選", "", "| Cohort | 病人 | 病毒 | 現行層級 | 接線後三軸 | 院端陽性檢出筆數 | 來源連結數 |", "|---|---|---|---|---|---:|---:|"]
    for r in summaries:
        lines.append(f"| {r['cohort']} | P{r['patient_id']} | {r['organism_name']} | {r['current_tier']} | {r['detection']}／{r['activity']}／{r['organ_support']} | {r['hospital_positive_claims']} | {r['linked_source_count']} |")
    lines += ["", "完整來源與排除理由：source_binding_audit.csv；逐候選packet：bound_candidate_packets.json；沒有對應候選的病毒資料：viral_facts_without_matching_candidate.json。", "",
        "下一階段再確認病毒特異的活動與肺炎歸因規則；本次沒有用已知答案、自動fallback或mNGS排名建立probable/proven疾病。", ""]
    (output_root / "BINDING_REPORT_zh.md").write_text("\n".join(lines), encoding="utf-8")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run(args.input_root, args.output_root), ensure_ascii=False, indent=2))

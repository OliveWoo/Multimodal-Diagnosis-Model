"""Build a loss-aware DNA/RNA mNGS inventory without changing model outputs.

The source files are already grouped by patient, infection/specimen event, test,
and organism.  This tool adds explicit QC/missingness states, within-test
normalization, and per-case organism unions.  Raw reads from different tests
are intentionally never summed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SOURCE_ROOT = REPO_ROOT / "33_patients_all_mNGS_DNA_RNA_RK_NTC_by_species"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "outputs" / "runs" / "2026-09-22_KH_multi_assay_mngs_inventory"
DEFAULT_SELECTED_DNA_ROOT = (
    REPO_ROOT / "outputs" / "patient_info_KH_0728_2Days" / "33_patients_selected_DNA_all_RK_NTC_microbes"
)
SOURCE_RE = re.compile(r"^NGS_patient_(\d+)_all_mNGS_DNA_RNA_RK_NTC_by_species\.json$")
SELECTED_DNA_RE = re.compile(r"^NGS_patient_(\d+)_DNA_.+_all_RK_NTC_microbes\.json$")
CATEGORIES = ("bacterial", "viral", "fungal", "others")
MISSING_MARKERS = {"", "-", "na", "n/a", "none", "null", "unknown"}

SOURCE_FIELD_CONTRACT = {
    "input_reads": {
        "source_field": "檢出總序列數",
        "unit": "reads",
        "missingness": "source_missing_when_null",
    },
    "host_ratio": {
        "source_field": "宿主序列佔比",
        "unit": "percent",
        "example": "99.66 means 99.66%",
        "missingness": "source_missing_when_null",
    },
    "report_time": {
        "missingness": "source_missing_when_blank",
        "decision_time_anchor": "specimen_collected_time",
    },
    "detected": {
        "scope": "organism_detected_in_another_test_of_the_same_case_or_in_this_test",
        "current_test_positive_definition": "reads > 0",
    },
    "raw_result": {
        "missingness": "no_corresponding_raw_file",
        "interpretation": "unevaluable_not_negative",
    },
}


def load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected object JSON: {path}")
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


def to_float(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def normalize_name(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def normalize_tags(values: Any) -> list[str]:
    if not isinstance(values, list):
        return []
    tags: list[str] = []
    for value in values:
        tag = str(value or "").strip().replace("Code_RK_NTC=RK08", "Code_RK_NTC=RK8")
        tag = tag.replace("Code_RK_NTC=RK00", "Code_RK_NTC=RK0")
        if tag and tag not in tags:
            tags.append(tag)
    return tags


def analytical_rank(tags: Iterable[str]) -> tuple[int | None, str]:
    values = set(tags)
    ntc = "Code_NTC=00" in values
    rk0 = "Code_RK_NTC=RK0" in values
    rk8 = "Code_RK_NTC=RK8" in values
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


def raw_category_state(raw: Any) -> str:
    value = str(raw or "").strip().lower()
    if value == "available":
        return "available"
    if value in {"missing", "unavailable", "not_available"}:
        return "raw_unavailable"
    return "unknown"


def is_missing(value: Any) -> bool:
    return str(value or "").strip().lower() in MISSING_MARKERS


def derive_test_qc(test: dict[str, Any]) -> dict[str, Any]:
    issues: list[str] = []
    warnings: list[str] = []
    seq_id = str(test.get("seq_id") or "").strip()
    nucleic_type = str(test.get("nucleic_type") or "").strip().upper()
    params = test.get("parameters") if isinstance(test.get("parameters"), dict) else {}
    input_reads = to_float(params.get("input_reads"))
    host_ratio = to_float(params.get("host_ratio"))
    raw_status = test.get("raw_result_status") if isinstance(test.get("raw_result_status"), dict) else {}
    category_status = {category: raw_category_state(raw_status.get(category)) for category in CATEGORIES}

    if not seq_id:
        issues.append("missing_seq_id")
    if nucleic_type not in {"DNA", "RNA"}:
        issues.append("invalid_or_missing_nucleic_type")
    if input_reads is None:
        issues.append("source_missing_input_reads")
    elif input_reads <= 0:
        issues.append("invalid_input_reads")
    if host_ratio is None:
        issues.append("source_missing_host_ratio")
    elif not 0 <= host_ratio <= 100:
        issues.append("invalid_host_ratio")
    if all(value != "available" for value in category_status.values()):
        issues.append("no_available_raw_category")
    elif any(value != "available" for value in category_status.values()):
        warnings.append("partial_raw_category_availability")

    specimen_status = str(params.get("specimen_status") or "").strip()
    if specimen_status.lower() not in MISSING_MARKERS:
        warnings.append(f"specimen_status:{specimen_status}")

    identity_invalid = any(issue in {"missing_seq_id", "invalid_or_missing_nucleic_type"} for issue in issues)
    no_raw = "no_available_raw_category" in issues
    numeric_problem = any(
        issue in {
            "source_missing_input_reads",
            "invalid_input_reads",
            "source_missing_host_ratio",
            "invalid_host_ratio",
        }
        for issue in issues
    )
    if identity_invalid:
        status = "invalid"
    elif no_raw:
        status = "raw_unavailable"
    elif numeric_problem or warnings:
        status = "partial_evaluable"
    else:
        status = "evaluable"

    estimated_nonhost_reads = None
    if input_reads and host_ratio is not None and 0 <= host_ratio < 100:
        estimated_nonhost_reads = input_reads * (1.0 - host_ratio / 100.0)
        if estimated_nonhost_reads <= 0:
            estimated_nonhost_reads = None

    normalization_status = "available"
    if input_reads is None and host_ratio is None:
        normalization_status = "source_denominators_missing"
    elif input_reads is None:
        normalization_status = "source_input_reads_missing"
    elif host_ratio is None:
        normalization_status = "source_host_ratio_missing"
    elif input_reads <= 0 or not 0 <= host_ratio <= 100:
        normalization_status = "invalid_source_denominator"

    report_time_status = "source_missing_confirmed" if is_missing(test.get("report_time")) else "available"
    raw_data_interpretation = (
        "no_corresponding_raw_file_unevaluable_not_negative"
        if no_raw
        else "partial_raw_files_available" if any(value != "available" for value in category_status.values())
        else "raw_files_available"
    )

    return {
        "status": status,
        "eligible_for_union": status in {"evaluable", "partial_evaluable"},
        "issues": issues,
        "warnings": warnings,
        "category_status": category_status,
        "input_reads": input_reads,
        "host_ratio_percent": host_ratio,
        "estimated_nonhost_reads": estimated_nonhost_reads,
        "normalization_status": normalization_status,
        "normalization_missingness": (
            "confirmed_absent_in_source"
            if normalization_status.startswith("source_")
            else None
        ),
        "report_time_status": report_time_status,
        "decision_time_anchor": "specimen_collected_time",
        "raw_data_interpretation": raw_data_interpretation,
    }


def observation_status(measurement: dict[str, Any], category_state: str) -> tuple[str, str, list[str]]:
    warnings: list[str] = []
    raw_available = measurement.get("raw_available") is True and category_state == "available"
    source_case_detected = measurement.get("detected") is True
    selected = measurement.get("selected") is True

    if not raw_available:
        return "raw_unavailable", "not_applicable", warnings
    reads = to_float(measurement.get("reads"))
    if reads is None or reads < 0:
        return "invalid", "not_applicable", ["invalid_reads"]
    observed_in_current_test = reads > 0
    if source_case_detected and not observed_in_current_test:
        warnings.append("case_manifest_presence_without_current_test_reads")
    if not source_case_detected and observed_in_current_test:
        warnings.append("positive_reads_but_case_detected_false")
    if selected and not observed_in_current_test:
        warnings.append("selected_without_current_test_reads")

    if observed_in_current_test and selected:
        selection = "selected"
    elif observed_in_current_test:
        selection = "filtered"
    elif source_case_detected:
        selection = "case_manifest_only"
    else:
        selection = "not_detected"
    return "evaluable", selection, warnings


def dominance_map(observations: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    positive = [item for item in observations if item["availability_status"] == "evaluable" and item["reads"] > 0]
    positive.sort(key=lambda item: (-item["reads"], item["organism_key"]))
    total_reads = sum(item["reads"] for item in positive)
    selected = [item for item in positive if item["selection_status"] == "selected"]
    selected.sort(key=lambda item: (-item["reads"], item["organism_key"]))

    result: dict[str, dict[str, Any]] = {}
    top_reads = positive[0]["reads"] if positive else 0.0
    second_reads = positive[1]["reads"] if len(positive) > 1 else 0.0
    top_ratio = top_reads / max(second_reads, 1.0) if positive else 0.0
    top_tier = "D3_dominant" if top_ratio >= 10 else "D2_moderate" if top_ratio >= 3 else "D1_low"
    overall_rank = {item["organism_key"]: index + 1 for index, item in enumerate(positive)}
    selected_rank = {item["organism_key"]: index + 1 for index, item in enumerate(selected)}
    for item in observations:
        key = item["organism_key"]
        result[key] = {
            "rank_by_reads_in_retained_universe": overall_rank.get(key),
            "rank_selected_by_reads": selected_rank.get(key),
            "read_fraction_in_retained_universe": (item["reads"] / total_reads) if total_reads > 0 else None,
            "dominance_tier_in_retained_universe": top_tier if overall_rank.get(key) == 1 else "D0_not_top",
            "dominance_scope_warning": "Computed only among organisms retained in the supplied by-species files.",
        }
    return result


def build_patient(source_path: Path) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    source = load_json(source_path)
    patient_id = str(source.get("patient_id") or "").strip()
    cases = source.get("cases") if isinstance(source.get("cases"), list) else []
    organisms = source.get("organisms") if isinstance(source.get("organisms"), list) else []

    test_lookup: dict[tuple[str, str], dict[str, Any]] = {}
    normalized_cases: list[dict[str, Any]] = []
    test_rows: list[dict[str, Any]] = []
    duplicate_test_keys: list[str] = []
    for case in cases:
        if not isinstance(case, dict):
            continue
        case_id = str(case.get("case_review_id") or case.get("specimen_code") or "").strip()
        normalized_tests: list[dict[str, Any]] = []
        for test in case.get("tests") or []:
            if not isinstance(test, dict):
                continue
            seq_id = str(test.get("seq_id") or "").strip()
            key = (case_id, seq_id)
            if key in test_lookup:
                duplicate_test_keys.append(f"{case_id}:{seq_id}")
            qc = derive_test_qc(test)
            enriched_test = {**test, "qc": qc}
            test_lookup[key] = enriched_test
            normalized_tests.append(enriched_test)
            test_rows.append(
                {
                    "patient_id": patient_id,
                    "case_review_id": case_id,
                    "specimen_code": case.get("specimen_code"),
                    "collected_time": case.get("collected_time"),
                    "specimen_site": case.get("specimen_site"),
                    "seq_id": seq_id,
                    "nucleic_type": test.get("nucleic_type"),
                    "start_time": test.get("start_time"),
                    "report_time": None if is_missing(test.get("report_time")) else test.get("report_time"),
                    "report_time_status": qc["report_time_status"],
                    "decision_time_anchor": qc["decision_time_anchor"],
                    "condition": (test.get("parameters") or {}).get("condition"),
                    "specimen_status": (test.get("parameters") or {}).get("specimen_status"),
                    "input_reads": qc["input_reads"],
                    "host_ratio_percent": qc["host_ratio_percent"],
                    "estimated_nonhost_reads": qc["estimated_nonhost_reads"],
                    "normalization_status": qc["normalization_status"],
                    "normalization_missingness": qc["normalization_missingness"],
                    "raw_data_interpretation": qc["raw_data_interpretation"],
                    "qc_status": qc["status"],
                    "eligible_for_union": qc["eligible_for_union"],
                    "qc_issues": ";".join(qc["issues"]),
                    "qc_warnings": ";".join(qc["warnings"]),
                    "category_status": json.dumps(qc["category_status"], ensure_ascii=False, sort_keys=True),
                }
            )
        normalized_cases.append({**case, "tests": normalized_tests})

    observations_by_test: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    orphan_measurements: list[dict[str, Any]] = []
    for organism in organisms:
        if not isinstance(organism, dict):
            continue
        category = str(organism.get("category") or "").strip().lower()
        name = str(organism.get("name") or "").strip()
        organism_key = normalize_name(name)
        for measurement in organism.get("measurements") or []:
            if not isinstance(measurement, dict):
                continue
            case_id = str(measurement.get("case_review_id") or "").strip()
            seq_id = str(measurement.get("seq_id") or "").strip()
            test = test_lookup.get((case_id, seq_id))
            if test is None:
                orphan_measurements.append({"organism": name, "case_review_id": case_id, "seq_id": seq_id})
                continue
            qc = test["qc"]
            category_state = qc["category_status"].get(category, "unknown")
            availability, selection, warnings = observation_status(measurement, category_state)
            reads = to_float(measurement.get("reads")) or 0.0
            rpm_total = reads / qc["input_reads"] * 1_000_000 if availability == "evaluable" and qc["input_reads"] else None
            rpm_nonhost = (
                reads / qc["estimated_nonhost_reads"] * 1_000_000
                if availability == "evaluable" and qc["estimated_nonhost_reads"]
                else None
            )
            tags = normalize_tags(measurement.get("selection_tags"))
            rank_priority, rank_rule = analytical_rank(tags)
            source_case_detected = measurement.get("detected") is True
            observed_in_current_test = availability == "evaluable" and reads > 0
            observations_by_test[(case_id, seq_id)].append(
                {
                    "patient_id": patient_id,
                    "case_review_id": case_id,
                    "seq_id": seq_id,
                    "nucleic_type": str(test.get("nucleic_type") or "").upper(),
                    "category": category,
                    "organism_name": name,
                    "organism_key": organism_key,
                    "reads": reads,
                    "detected": source_case_detected,
                    "source_case_detected": source_case_detected,
                    "source_detected_scope": "same_case_any_test",
                    "observed_in_current_test": observed_in_current_test,
                    "current_test_positive_basis": "reads > 0",
                    "selected": measurement.get("selected") is True,
                    "availability_status": availability,
                    "selection_status": selection,
                    "selection_tags": tags,
                    "analytical_rank_priority": rank_priority,
                    "analytical_rank_rule": rank_rule,
                    "rpm_total": rpm_total,
                    "rpm_nonhost_estimated": rpm_nonhost,
                    "normalization_note": (
                        "input_reads is total detected sequence count; host_ratio is a percentage. "
                        "RPM_nonhost remains an estimate and is not compared across protocol variants."
                    ),
                    "observation_warnings": warnings,
                    "test_qc_status": qc["status"],
                    "test_eligible_for_union": qc["eligible_for_union"],
                    "protocol_condition": (test.get("parameters") or {}).get("condition"),
                }
            )

    all_observations: list[dict[str, Any]] = []
    for test_key, observations in observations_by_test.items():
        ranks = dominance_map(observations)
        for observation in observations:
            observation.update(ranks.get(observation["organism_key"], {}))
            all_observations.append(observation)

    case_lookup = {
        str(case.get("case_review_id") or case.get("specimen_code") or "").strip(): case
        for case in normalized_cases
    }
    aggregate_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for observation in all_observations:
        aggregate_groups[(observation["case_review_id"], observation["organism_key"])].append(observation)

    candidate_rows: list[dict[str, Any]] = []
    case_unions: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for (case_id, _), values in sorted(aggregate_groups.items()):
        evaluable = [value for value in values if value["availability_status"] == "evaluable" and value["test_eligible_for_union"]]
        detected = [value for value in evaluable if value["observed_in_current_test"]]
        selected = [value for value in detected if value["selected"]]
        if not detected:
            continue
        dna_detected = [value for value in detected if value["nucleic_type"] == "DNA"]
        rna_detected = [value for value in detected if value["nucleic_type"] == "RNA"]
        selected_union = bool(selected)
        row = {
            "patient_id": patient_id,
            "case_review_id": case_id,
            "specimen_code": case_lookup.get(case_id, {}).get("specimen_code"),
            "collected_time": case_lookup.get(case_id, {}).get("collected_time"),
            "specimen_site": case_lookup.get(case_id, {}).get("specimen_site"),
            "category": values[0]["category"],
            "organism_name": values[0]["organism_name"],
            "organism_key": values[0]["organism_key"],
            "detected_union": True,
            "selected_union": selected_union,
            "union_status": "selected_any_test" if selected_union else "filtered_only",
            "evaluable_test_count": len(evaluable),
            "detected_test_count": len(detected),
            "selected_test_count": len(selected),
            "dna_detected_test_count": len(dna_detected),
            "rna_detected_test_count": len(rna_detected),
            "cross_molecule_support": bool(dna_detected and rna_detected),
            "repeated_detection": len(detected) >= 2,
            "distinct_protocol_conditions": sorted({str(value.get("protocol_condition") or "-") for value in detected}),
            "max_reads_single_test": max(value["reads"] for value in detected),
            "max_rpm_total_single_test": max((value["rpm_total"] for value in detected if value["rpm_total"] is not None), default=None),
            "max_rpm_nonhost_estimated_single_test": max(
                (value["rpm_nonhost_estimated"] for value in detected if value["rpm_nonhost_estimated"] is not None),
                default=None,
            ),
            "best_analytical_rank_priority": min(
                (value["analytical_rank_priority"] for value in selected if value["analytical_rank_priority"] is not None),
                default=None,
            ),
            "best_rank_in_retained_universe": min(
                (value["rank_by_reads_in_retained_universe"] for value in detected if value["rank_by_reads_in_retained_universe"] is not None),
                default=None,
            ),
            "reads_sum_across_tests": None,
            "aggregation_warning": "Raw reads are not summed across DNA/RNA, library preparations, or reruns.",
            "observations": detected,
        }
        case_unions[case_id].append(row)
        candidate_rows.append({key: value for key, value in row.items() if key != "observations"})

    normalized_case_output: list[dict[str, Any]] = []
    for case in normalized_cases:
        case_id = str(case.get("case_review_id") or case.get("specimen_code") or "").strip()
        candidates = sorted(
            case_unions.get(case_id, []),
            key=lambda item: (
                not item["selected_union"],
                item["best_analytical_rank_priority"] or 99,
                -item["max_reads_single_test"],
                item["organism_name"],
            ),
        )
        normalized_case_output.append(
            {
                **case,
                "candidate_union": {
                    "selected_union": [item for item in candidates if item["selected_union"]],
                    "detected_but_filtered_union": [item for item in candidates if not item["selected_union"]],
                },
            }
        )

    patient_output = {
        "schema_version": "multi_assay_mngs_evidence.v2",
        "patient_id": patient_id,
        "source": {"path": str(source_path.resolve()), "sha256": sha256_file(source_path)},
        "source_field_contract": SOURCE_FIELD_CONTRACT,
        "scope_notes": [
            "All tests are retained separately by case/specimen and seq_id.",
            "Current-test positivity is defined by reads > 0. The source detected flag has same-case scope and is not current-test evidence by itself.",
            "Candidate observations contain only evaluable, union-eligible measurements with reads > 0; all source observations remain in the QC audit.",
            "Raw reads are never summed across tests.",
            "Missing raw files are unevaluable and are never interpreted as negative tests.",
            "Missing report time is source missingness; specimen collection time remains the decision anchor.",
            "Dominance is provisional because the supplied files contain a retained organism universe, not every raw taxon.",
            "Adapter dilution is preserved as protocol metadata and is not used to rescale reads.",
        ],
        "cases": normalized_case_output,
        "observations": all_observations,
        "audit": {
            "duplicate_test_keys": duplicate_test_keys,
            "orphan_measurements": orphan_measurements,
        },
    }
    return patient_output, test_rows, candidate_rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8-sig")
        return
    fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict)) else value
                    for key, value in row.items()
                }
            )


def load_selected_dna_index(root: Path) -> dict[tuple[str, str], dict[str, Any]]:
    index: dict[tuple[str, str], dict[str, Any]] = {}
    if not root.is_dir():
        return index
    for path in sorted(root.glob("NGS_patient_*_DNA_*_all_RK_NTC_microbes.json")):
        match = SELECTED_DNA_RE.fullmatch(path.name)
        if match is None:
            continue
        patient_id = match.group(1)
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        if not isinstance(payload, list):
            continue
        for record in payload:
            pathogens = record.get("pathogens") if isinstance(record, dict) else None
            if not isinstance(pathogens, dict):
                continue
            for category in CATEGORIES:
                for item in pathogens.get(category) or []:
                    if not isinstance(item, dict):
                        continue
                    name = str(item.get("name") or "").strip()
                    key = normalize_name(name)
                    if not key:
                        continue
                    index[(patient_id, key)] = {
                        "patient_id": patient_id,
                        "organism_key": key,
                        "organism_name": name,
                        "category": category,
                        "source_path": str(path.resolve()),
                    }
    return index


def compare_selected_dna(
    candidate_rows: list[dict[str, Any]], selected_dna_root: Path
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    baseline = load_selected_dna_index(selected_dna_root)
    current: dict[tuple[str, str], dict[str, Any]] = {}
    for row in candidate_rows:
        if not row["selected_union"]:
            continue
        key = (str(row["patient_id"]), str(row["organism_key"]))
        existing = current.get(key)
        if existing is None or (bool(row["cross_molecule_support"]), row["detected_test_count"]) > (
            bool(existing["cross_molecule_support"]),
            existing["detected_test_count"],
        ):
            current[key] = row

    rows: list[dict[str, Any]] = []
    for key in sorted(set(baseline) | set(current), key=lambda item: (int(item[0]), item[1])):
        old = baseline.get(key)
        new = current.get(key)
        if old and new:
            status = "retained_from_selected_dna"
        elif new:
            status = "new_from_multi_assay_union"
        else:
            status = "missing_from_multi_assay_union"
        rows.append(
            {
                "patient_id": key[0],
                "organism_key": key[1],
                "organism_name": (new or old or {}).get("organism_name"),
                "category": (new or old or {}).get("category"),
                "comparison_status": status,
                "cross_molecule_support": new.get("cross_molecule_support") if new else None,
                "dna_detected_test_count": new.get("dna_detected_test_count") if new else None,
                "rna_detected_test_count": new.get("rna_detected_test_count") if new else None,
                "max_reads_single_test": new.get("max_reads_single_test") if new else None,
            }
        )

    new_only = [row for row in rows if row["comparison_status"] == "new_from_multi_assay_union"]
    molecule_counts: Counter[str] = Counter()
    for row in new_only:
        if row["cross_molecule_support"]:
            molecule_counts["DNA+RNA"] += 1
        elif row["dna_detected_test_count"]:
            molecule_counts["DNA-only"] += 1
        elif row["rna_detected_test_count"]:
            molecule_counts["RNA-only"] += 1
    summary = {
        "selected_dna_root": str(selected_dna_root.resolve()),
        "selected_dna_distinct_patient_organism_count": len(baseline),
        "multi_assay_selected_union_distinct_patient_organism_count": len(current),
        "retained_from_selected_dna_count": sum(row["comparison_status"] == "retained_from_selected_dna" for row in rows),
        "new_from_multi_assay_union_count": len(new_only),
        "missing_from_multi_assay_union_count": sum(
            row["comparison_status"] == "missing_from_multi_assay_union" for row in rows
        ),
        "new_from_multi_assay_by_molecule_support": dict(sorted(molecule_counts.items())),
        "new_from_multi_assay_by_category": dict(sorted(Counter(str(row["category"]) for row in new_only).items())),
    }
    return rows, summary


def build_report(summary: dict[str, Any]) -> str:
    qc = summary["test_qc_status_counts"]
    union = summary["candidate_union_counts"]
    return f"""# KH multi-assay mNGS inventory

This is a shadow inventory. It does not replace the selected-DNA baseline or change Picked/OBER outputs.

## Coverage

| Item | Count |
| --- | ---: |
| Patients | {summary['patient_count']} |
| Infection/specimen events | {summary['case_count']} |
| Sequencing tests | {summary['test_count']} |
| DNA tests | {summary['nucleic_type_counts'].get('DNA', 0)} |
| RNA tests | {summary['nucleic_type_counts'].get('RNA', 0)} |
| BALF tests | {summary['specimen_site_counts'].get('BALF', 0)} |
| Blood tests | {summary['specimen_site_counts'].get('Blood', 0)} |
| Organism-test observations | {summary['observation_count']} |

## QC states

| State | Tests |
| --- | ---: |
| evaluable | {qc.get('evaluable', 0)} |
| partial_evaluable | {qc.get('partial_evaluable', 0)} |
| raw_unavailable | {qc.get('raw_unavailable', 0)} |
| invalid | {qc.get('invalid', 0)} |

`partial_evaluable` remains eligible for the union but carries explicit warnings. Hemolysis and other specimen notes are warnings, not automatic invalidation.

## Candidate unions

| Union | Patient-case-organism rows |
| --- | ---: |
| selected in at least one evaluable test | {union.get('selected_any_test', 0)} |
| detected but filtered in every evaluable test | {union.get('filtered_only', 0)} |
| DNA and RNA both detected | {summary['cross_molecule_candidate_count']} |
| detected in at least two tests | {summary['repeated_detection_candidate_count']} |

Selected-union molecular support:

| Support | Candidates |
| --- | ---: |
| DNA + RNA | {summary['selected_union_molecule_support_counts'].get('DNA+RNA', 0)} |
| DNA only | {summary['selected_union_molecule_support_counts'].get('DNA-only', 0)} |
| RNA only | {summary['selected_union_molecule_support_counts'].get('RNA-only', 0)} |

Selected-union organism categories:

| Category | Candidates |
| --- | ---: |
| bacterial | {summary['selected_union_category_counts'].get('bacterial', 0)} |
| viral | {summary['selected_union_category_counts'].get('viral', 0)} |
| fungal | {summary['selected_union_category_counts'].get('fungal', 0)} |
| others | {summary['selected_union_category_counts'].get('others', 0)} |

Observation selection states:

| State | Observations |
| --- | ---: |
| selected | {summary['observation_selection_status_counts'].get('selected', 0)} |
| filtered | {summary['observation_selection_status_counts'].get('filtered', 0)} |
| case manifest only (zero reads in this test) | {summary['observation_selection_status_counts'].get('case_manifest_only', 0)} |
| not_detected | {summary['observation_selection_status_counts'].get('not_detected', 0)} |
| not_applicable | {summary['observation_selection_status_counts'].get('not_applicable', 0)} |

## Comparison with the selected-DNA source

| Item | Distinct patient-organism pairs |
| --- | ---: |
| Old selected-DNA source | {summary['selected_dna_comparison'].get('selected_dna_distinct_patient_organism_count', 0)} |
| Multi-assay selected union | {summary['selected_dna_comparison'].get('multi_assay_selected_union_distinct_patient_organism_count', 0)} |
| Old candidates retained | {summary['selected_dna_comparison'].get('retained_from_selected_dna_count', 0)} |
| Newly exposed by all tests | {summary['selected_dna_comparison'].get('new_from_multi_assay_union_count', 0)} |
| Old candidates missing | {summary['selected_dna_comparison'].get('missing_from_multi_assay_union_count', 0)} |

Newly exposed candidates include {summary['selected_dna_comparison'].get('new_from_multi_assay_by_molecule_support', {}).get('RNA-only', 0)} RNA-only, {summary['selected_dna_comparison'].get('new_from_multi_assay_by_molecule_support', {}).get('DNA-only', 0)} alternate-DNA-only, and {summary['selected_dna_comparison'].get('new_from_multi_assay_by_molecule_support', {}).get('DNA+RNA', 0)} DNA+RNA-supported patient-organism pairs. These are discovery candidates, not automatic pathogens.

## Interpretation boundary

- `selected_union` is the closest answer-blind replacement for the old single selected-DNA candidate source.
- `detected_but_filtered_union` is retained for audit and future high-consequence rescue; it is not automatically promoted into the model.
- RPM is calculated per test. Raw reads and RPM are never added across DNA/RNA tests or protocol variants.
- `input_reads` is the Excel `檢出總序列數` in reads; `host_ratio` is the Excel `宿主序列佔比` stored as percent (for example, 99.66 means 99.66%).
- `rpm_nonhost_estimated` remains a derived estimate and is not compared across protocol variants.
- Rank and dominance are explicitly labelled as being within the supplied retained organism universe.
- Adapter dilution, centrifugation, PCR cycles, and reruns remain technical protocol metadata.
- {summary['tests_missing_normalization_denominator_count']} tests lack `input_reads` and/or `host_ratio`, so RPM normalization is unavailable for them.
- {summary['tests_missing_report_time_count']} tests lack report time in the source; specimen collection time remains the decision anchor.
- {summary['raw_unavailable_test_count']} tests have no corresponding raw file; they are unevaluable, not negative.
- {summary['observation_warning_counts'].get('case_manifest_presence_without_current_test_reads', 0)} observations carry the source case-level `detected=true` flag but have zero reads in this test. They mean another test in the same case detected the organism and do not enter the positive candidate union.

## Downstream dependency that must change

The current selected-DNA pipeline still uses the simulated manual-screening branch:

1. ranked candidates from one selected DNA file;
2. `run_mngs_to_specimen_v4_deterministic.py`;
3. `optimize_mngs_rank_read_precision.py`;
4. `build_chosen_ranked_mngs_from_mngs_to_specimen.py`;
5. deterministic scorer and Luna review.

That branch must be made test-aware before multi-assay candidates become a formal model input. Until then, this inventory stays shadow-only.
"""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--selected-dna-root", type=Path, default=DEFAULT_SELECTED_DNA_ROOT)
    parser.add_argument("--expected-patients", type=int, default=33)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    source_files = sorted(args.source_root.glob("NGS_patient_*_all_mNGS_DNA_RNA_RK_NTC_by_species.json"))
    if len(source_files) != args.expected_patients:
        raise ValueError(f"Expected {args.expected_patients} patient files, found {len(source_files)}")

    patient_outputs: list[dict[str, Any]] = []
    test_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    observation_rows: list[dict[str, Any]] = []
    seen_patients: set[str] = set()
    patient_dir = args.output_root / "patients"
    for source_path in source_files:
        match = SOURCE_RE.fullmatch(source_path.name)
        if match is None:
            raise ValueError(f"Unexpected source filename: {source_path.name}")
        expected_id = match.group(1)
        patient_output, patient_tests, patient_candidates = build_patient(source_path)
        if patient_output["patient_id"] != expected_id:
            raise ValueError(f"Filename/payload patient mismatch: {source_path}")
        if expected_id in seen_patients:
            raise ValueError(f"Duplicate patient id: {expected_id}")
        seen_patients.add(expected_id)
        patient_outputs.append(patient_output)
        test_rows.extend(patient_tests)
        candidate_rows.extend(patient_candidates)
        observation_rows.extend(patient_output["observations"])
        write_json(patient_dir / f"NGS_patient_{expected_id}_multi_assay_mngs_evidence.json", patient_output)

    qc_counts = Counter(row["qc_status"] for row in test_rows)
    nucleic_counts = Counter(str(row["nucleic_type"]) for row in test_rows)
    site_counts = Counter(str(row["specimen_site"]) for row in test_rows)
    union_counts = Counter(row["union_status"] for row in candidate_rows)
    observation_availability = Counter(row["availability_status"] for row in observation_rows)
    observation_selection = Counter(row["selection_status"] for row in observation_rows)
    observation_warnings: Counter[str] = Counter()
    for row in observation_rows:
        observation_warnings.update(row["observation_warnings"])
    positive_observation_rows = [
        row for row in observation_rows
        if row["availability_status"] == "evaluable"
        and row["test_eligible_for_union"]
        and row["observed_in_current_test"]
    ]
    selected_support: Counter[str] = Counter()
    selected_categories: Counter[str] = Counter()
    for row in candidate_rows:
        if not row["selected_union"]:
            continue
        selected_categories[str(row["category"])] += 1
        if row["cross_molecule_support"]:
            selected_support["DNA+RNA"] += 1
        elif row["dna_detected_test_count"]:
            selected_support["DNA-only"] += 1
        elif row["rna_detected_test_count"]:
            selected_support["RNA-only"] += 1
        else:
            selected_support["neither"] += 1
    comparison_rows, comparison_summary = compare_selected_dna(candidate_rows, args.selected_dna_root)
    summary = {
        "schema_version": "multi_assay_mngs_inventory_summary.v2",
        "source_field_contract": SOURCE_FIELD_CONTRACT,
        "source_root": str(args.source_root.resolve()),
        "output_root": str(args.output_root.resolve()),
        "patient_count": len(patient_outputs),
        "case_count": sum(len(patient["cases"]) for patient in patient_outputs),
        "test_count": len(test_rows),
        "observation_count": len(observation_rows),
        "positive_observation_count": len(positive_observation_rows),
        "nucleic_type_counts": dict(sorted(nucleic_counts.items())),
        "specimen_site_counts": dict(sorted(site_counts.items())),
        "test_qc_status_counts": dict(sorted(qc_counts.items())),
        "candidate_union_counts": dict(sorted(union_counts.items())),
        "selected_union_molecule_support_counts": dict(sorted(selected_support.items())),
        "selected_union_category_counts": dict(sorted(selected_categories.items())),
        "observation_availability_status_counts": dict(sorted(observation_availability.items())),
        "observation_selection_status_counts": dict(sorted(observation_selection.items())),
        "observation_warning_counts": dict(sorted(observation_warnings.items())),
        "tests_missing_normalization_denominator_count": sum(
            row["input_reads"] in {None, ""} or row["host_ratio_percent"] in {None, ""} for row in test_rows
        ),
        "tests_missing_report_time_count": sum(row["report_time_status"] == "source_missing_confirmed" for row in test_rows),
        "raw_unavailable_test_count": sum(row["qc_status"] == "raw_unavailable" for row in test_rows),
        "cross_molecule_candidate_count": sum(bool(row["cross_molecule_support"]) for row in candidate_rows),
        "repeated_detection_candidate_count": sum(bool(row["repeated_detection"]) for row in candidate_rows),
        "duplicate_test_key_count": sum(len(patient["audit"]["duplicate_test_keys"]) for patient in patient_outputs),
        "orphan_measurement_count": sum(len(patient["audit"]["orphan_measurements"]) for patient in patient_outputs),
        "selected_dna_comparison": comparison_summary,
    }
    args.output_root.mkdir(parents=True, exist_ok=True)
    write_json(args.output_root / "inventory_summary.json", summary)
    write_csv(args.output_root / "test_qc.csv", test_rows)
    write_csv(args.output_root / "observation_qc.csv", observation_rows)
    write_csv(args.output_root / "positive_observations.csv", positive_observation_rows)
    write_csv(args.output_root / "candidate_union.csv", candidate_rows)
    write_csv(args.output_root / "selected_dna_comparison.csv", comparison_rows)
    (args.output_root / "INVENTORY_REPORT.md").write_text(build_report(summary), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

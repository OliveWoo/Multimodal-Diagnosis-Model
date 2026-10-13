"""Combine validated pipeline-input bundles without exposing labels to inference.

The output keeps every patient source file byte-for-byte, regenerates the
bundle-level indexes from the copied files, and records structural/data-quality
checks. Ground truth remains in layer 05 and is explicitly evaluation-only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASE = ROOT / "2026-10-08_44patients_pipeline_input_with_ground_truth"
DEFAULT_ADDITION = ROOT / "2026-10-08_11additional_two_hospital_pipeline_input_with_ground_truth"
DEFAULT_OUTPUT = ROOT / "2026-10-09_55patients_pipeline_input_with_ground_truth"

PATIENT_RE = re.compile(r"^NGS_patient_(\d+)$")
CLINICAL_MODULES = (
    "admission_diagnosis",
    "CBC",
    "culture",
    "filmarray",
    "gm_test",
    "image",
    "molecular_microbiology",
    "other_lab",
    "underlying",
)
OPTIONAL_MODULES = ("medication",)
ROOT_AUDIT_FILES = (
    "bundle_manifest.json",
    "FILE_MANIFEST.csv",
    "GROUND_TRUTH_MANIFEST.json",
    "GROUND_TRUTH_VALIDATION.json",
    "LAYER04_CLEAR_VALIDATION.json",
    "MISSING_DATA.csv",
    "PATIENT_INDEX.csv",
    "README.md",
    "README_44_PATIENT_DATA_PACKAGE.md",
    "SHA256SUMS.csv",
    "SOURCE_ROSTER.csv",
    "WINDOW_FILTER_VALIDATION.json",
    "CURRENT_PIPELINE_COMPATIBILITY_REPORT.json",
)


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def patient_dirs(root: Path) -> list[tuple[str, int, Path]]:
    rows: list[tuple[str, int, Path]] = []
    patients_root = root / "patients"
    if not patients_root.is_dir():
        raise FileNotFoundError(patients_root)
    for cohort_dir in sorted(path for path in patients_root.iterdir() if path.is_dir()):
        for patient_dir in sorted(path for path in cohort_dir.iterdir() if path.is_dir()):
            match = PATIENT_RE.fullmatch(patient_dir.name)
            if match:
                rows.append((cohort_dir.name, int(match.group(1)), patient_dir))
    return sorted(rows, key=lambda row: (row[1], row[0]))


def verify_source_hashes(root: Path) -> dict[str, Any]:
    manifest = root / "SHA256SUMS.csv"
    errors: list[str] = []
    checked = 0
    if not manifest.is_file():
        return {"status": "missing", "checked": 0, "errors": [str(manifest)]}
    with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            relative = row.get("relative_path") or row.get("path") or row.get("file")
            expected = (row.get("sha256") or row.get("SHA256") or "").lower()
            if not relative or not expected:
                errors.append(f"invalid hash row: {row}")
                continue
            path = root / relative.replace("/", "\\")
            if not path.is_file():
                errors.append(f"missing: {relative}")
                continue
            checked += 1
            actual = sha256_file(path).lower()
            if actual != expected:
                errors.append(f"mismatch: {relative}")
    return {
        "status": "passed" if not errors else "failed",
        "checked": checked,
        "errors": errors,
    }


def status_for_json(path: Path) -> tuple[str, int]:
    if not path.is_file():
        return "missing", 0
    payload = load_json(path)
    if not isinstance(payload, list):
        raise TypeError(f"Expected a JSON list: {path}")
    return ("present_nonempty" if payload else "present_empty"), len(payload)


def patient_audit(cohort: str, patient_id: int, patient_dir: Path) -> dict[str, Any]:
    patient_key = f"NGS_patient_{patient_id}"
    errors: list[str] = []
    warnings: list[str] = []
    manifest_path = patient_dir / "PATIENT_MANIFEST.json"
    manifest = load_json(manifest_path) if manifest_path.is_file() else {}
    if not manifest_path.is_file():
        errors.append("PATIENT_MANIFEST.json missing")
    if manifest and manifest.get("patient_id") != patient_id:
        errors.append("patient manifest ID does not match directory")
    if manifest and manifest.get("cohort") != cohort:
        errors.append("patient manifest cohort does not match directory")

    layer_01 = patient_dir / "01_upstream_raw_clinical_sheet_exports"
    layer_02 = patient_dir / "02_pipeline_direct_clinical_inputs"
    layer_03 = patient_dir / "03_raw_mngs_per_test"
    layer_04 = patient_dir / "04_prior_derived_outputs_reference_only"
    layer_05 = patient_dir / "05_ground_truth"
    if not layer_02.is_dir():
        errors.append("layer 02 missing")
    if not layer_03.is_dir():
        errors.append("layer 03 missing")
    if not layer_05.is_dir():
        errors.append("layer 05 missing")
    prior_files = [path for path in layer_04.rglob("*") if path.is_file()] if layer_04.is_dir() else []
    if prior_files:
        errors.append(f"layer 04 is not empty ({len(prior_files)} files)")

    module_status: dict[str, str] = {}
    module_counts: dict[str, int] = {}
    for module in (*CLINICAL_MODULES, *OPTIONAL_MODULES):
        path = layer_02 / f"{patient_key}_{module}.json"
        try:
            status, count = status_for_json(path)
        except (json.JSONDecodeError, TypeError) as exc:
            status, count = "invalid", 0
            errors.append(f"invalid layer 02 {module}: {exc}")
        module_status[module] = status
        module_counts[module] = count
        if module in CLINICAL_MODULES and status == "missing":
            warnings.append(f"optional clinical module missing: {module}")

    mngs_files = sorted(layer_03.glob("*.json")) if layer_03.is_dir() else []
    mngs: dict[str, Any] = {}
    if len(mngs_files) != 1:
        errors.append(f"expected exactly one mNGS JSON, found {len(mngs_files)}")
    elif mngs_files:
        try:
            mngs = load_json(mngs_files[0])
        except json.JSONDecodeError as exc:
            errors.append(f"invalid mNGS JSON: {exc}")
    if mngs and mngs.get("patient_id") != patient_id:
        errors.append("mNGS patient_id does not match directory")

    cases = mngs.get("cases") or [] if isinstance(mngs, dict) else []
    organisms = mngs.get("organisms") or [] if isinstance(mngs, dict) else []
    case_ids: set[str] = set()
    specimen_codes: set[str] = set()
    test_keys: set[tuple[str, str]] = set()
    collected_times: list[str] = []
    test_count = 0
    missing_input_reads = 0
    missing_host_ratio = 0
    missing_specimen_site = 0
    missing_collected_time = 0
    missing_nucleic_type = 0
    missing_report_time = 0
    for case in cases:
        if not isinstance(case, dict):
            errors.append("non-object mNGS case")
            continue
        raw_case_id = str(case.get("case_review_id") or "").strip()
        specimen = str(case.get("specimen_code") or "").strip()
        case_id = raw_case_id or specimen
        if not case_id:
            errors.append("mNGS case missing both case_review_id and specimen_code")
        elif case_id in case_ids:
            errors.append(f"duplicate mNGS event identifier: {case_id}")
        case_ids.add(case_id)
        if raw_case_id:
            case_ids.add(raw_case_id)
        if specimen:
            specimen_codes.add(specimen)
            case_ids.add(specimen)
        if not raw_case_id and specimen:
            warnings.append(f"case_review_id missing; specimen_code used as event identifier: {specimen}")
        collected = str(case.get("collected_time") or "").strip()
        if collected:
            collected_times.append(collected)
        else:
            missing_collected_time += 1
        if not str(case.get("specimen_site") or "").strip():
            missing_specimen_site += 1
        for test in case.get("tests") or []:
            if not isinstance(test, dict):
                errors.append(f"non-object test in case {case_id}")
                continue
            test_count += 1
            seq_id = str(test.get("seq_id") or "").strip()
            if not seq_id:
                errors.append(f"test missing seq_id in case {case_id}")
            key = (raw_case_id, seq_id)
            if key in test_keys:
                errors.append(f"duplicate test key: {raw_case_id or specimen}/{seq_id}")
            test_keys.add(key)
            test_keys.add((case_id, seq_id))
            if not str(test.get("nucleic_type") or "").strip():
                missing_nucleic_type += 1
            if not str(test.get("report_time") or "").strip():
                missing_report_time += 1
            parameters = test.get("parameters") or {}
            if parameters.get("input_reads") in (None, ""):
                missing_input_reads += 1
            if parameters.get("host_ratio") in (None, ""):
                missing_host_ratio += 1

    measurement_count = 0
    positive_measurement_count = 0
    ignored_zero_read_placeholders = 0
    unavailable_measurements = 0
    invalid_nonpositive_measurements = 0
    invalid_references = 0
    duplicate_measurements = 0
    seen_measurements: set[tuple[str, str, str, str]] = set()
    for organism in organisms:
        if not isinstance(organism, dict):
            errors.append("non-object organism")
            continue
        name = str(organism.get("name") or "").strip()
        category = str(organism.get("category") or "").strip()
        if not name:
            errors.append("organism missing name")
        for measurement in organism.get("measurements") or []:
            if not isinstance(measurement, dict):
                errors.append(f"non-object measurement for {name}")
                continue
            measurement_count += 1
            case_id = str(measurement.get("case_review_id") or "").strip()
            seq_id = str(measurement.get("seq_id") or "").strip()
            if (case_id, seq_id) not in test_keys:
                invalid_references += 1
            key = (category.casefold(), name.casefold(), case_id, seq_id)
            if key in seen_measurements:
                duplicate_measurements += 1
            seen_measurements.add(key)
            reads = measurement.get("reads")
            detected = measurement.get("detected")
            selected = measurement.get("selected")
            raw_available = measurement.get("raw_available")
            if isinstance(reads, (int, float)) and reads > 0:
                positive_measurement_count += 1
            elif reads == 0 and detected is False and selected is False:
                ignored_zero_read_placeholders += 1
            elif reads is None and raw_available is False and detected in (None, False) and selected in (None, False):
                unavailable_measurements += 1
            else:
                invalid_nonpositive_measurements += 1
            if isinstance(reads, (int, float)) and reads > 0 and detected is not True:
                errors.append(f"positive reads but detected is not true: {name}/{case_id}/{seq_id}")
    if invalid_nonpositive_measurements:
        errors.append(f"invalid non-positive measurements: {invalid_nonpositive_measurements}")
    if invalid_references:
        errors.append(f"measurement references missing tests: {invalid_references}")
    if duplicate_measurements:
        errors.append(f"duplicate organism measurements: {duplicate_measurements}")
    if not cases:
        errors.append("mNGS cases empty")
    if not organisms:
        warnings.append("mNGS organisms empty")

    gt_files = sorted(layer_05.glob("*.json")) if layer_05.is_dir() else []
    gt: dict[str, Any] = {}
    if len(gt_files) != 1:
        errors.append(f"expected exactly one ground-truth JSON, found {len(gt_files)}")
    elif gt_files:
        try:
            gt = load_json(gt_files[0])
        except json.JSONDecodeError as exc:
            errors.append(f"invalid ground-truth JSON: {exc}")
    if gt:
        if gt.get("patient_id") != patient_id:
            errors.append("ground-truth patient_id does not match directory")
        if gt.get("cohort") != cohort:
            errors.append("ground-truth cohort does not match directory")
        if gt.get("evaluation_only") is not True:
            errors.append("ground truth is not marked evaluation_only")
        for case in gt.get("ground_truth_cases") or []:
            gt_case = str(case.get("case_review_id") or "").strip()
            gt_specimen = str(case.get("specimen_code") or "").strip()
            if gt_case and gt_case not in case_ids:
                errors.append(f"ground-truth case not found in mNGS: {gt_case}")
            if gt_specimen and gt_specimen not in specimen_codes:
                errors.append(f"ground-truth specimen not found in mNGS: {gt_specimen}")

    raw_count = len([path for path in layer_01.glob("*.json")]) if layer_01.is_dir() else 0
    direct_nonempty = sum(status == "present_nonempty" for status in module_status.values())
    direct_empty = sum(status == "present_empty" for status in module_status.values())
    direct_missing = sum(status == "missing" for status in module_status.values())
    return {
        "cohort": cohort,
        "patient_id": patient_id,
        "patient_key": patient_key,
        "patient_dir": patient_dir,
        "manifest": manifest,
        "mngs": mngs,
        "mngs_file": mngs_files[0] if len(mngs_files) == 1 else None,
        "gt_file": gt_files[0] if len(gt_files) == 1 else None,
        "module_status": module_status,
        "module_counts": module_counts,
        "direct_nonempty": direct_nonempty,
        "direct_empty": direct_empty,
        "direct_missing": direct_missing,
        "raw_export_available_count": raw_count,
        "mngs_case_count": len(cases),
        "mngs_test_count": test_count,
        "mngs_organism_count": len(organisms),
        "mngs_measurement_count": measurement_count,
        "positive_measurement_count": positive_measurement_count,
        "ignored_zero_read_placeholders": ignored_zero_read_placeholders,
        "unavailable_measurements": unavailable_measurements,
        "invalid_nonpositive_measurements": invalid_nonpositive_measurements,
        "missing_input_reads": missing_input_reads,
        "missing_host_ratio": missing_host_ratio,
        "missing_specimen_site": missing_specimen_site,
        "missing_collected_time": missing_collected_time,
        "missing_nucleic_type": missing_nucleic_type,
        "missing_report_time": missing_report_time,
        "specimen_codes": sorted(specimen_codes),
        "collected_times": sorted(collected_times),
        "errors": errors,
        "warnings": warnings,
    }


def copy_source_audits(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for name in ROOT_AUDIT_FILES:
        path = source / name
        if path.is_file():
            shutil.copy2(path, destination / name)


def build(base_root: Path, addition_root: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output root: {output_root}")

    base_hashes = verify_source_hashes(base_root)
    addition_hashes = verify_source_hashes(addition_root)
    if base_hashes["status"] != "passed" or addition_hashes["status"] != "passed":
        raise RuntimeError("Source SHA256 validation failed")

    base_patients = patient_dirs(base_root)
    addition_patients = patient_dirs(addition_root)
    base_ids = {patient_id for _, patient_id, _ in base_patients}
    addition_ids = {patient_id for _, patient_id, _ in addition_patients}
    overlaps = sorted(base_ids & addition_ids)
    if overlaps:
        raise ValueError(f"Patient IDs overlap between source bundles: {overlaps}")

    source_rows = [
        ("base_44", base_root, *row) for row in base_patients
    ] + [
        ("additional_11", addition_root, *row) for row in addition_patients
    ]
    audits = [patient_audit(cohort, patient_id, path) for _, _, cohort, patient_id, path in source_rows]
    source_errors = [
        f"{row['cohort']}/P{row['patient_id']}: {message}"
        for row in audits
        for message in row["errors"]
    ]
    if source_errors:
        raise RuntimeError("Patient source validation failed:\n" + "\n".join(source_errors))

    output_root.mkdir(parents=True, exist_ok=True)
    for _, _, cohort, _, source_patient in source_rows:
        cohort_root = output_root / "patients" / cohort
        cohort_root.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source_patient, cohort_root / source_patient.name)
    copy_source_audits(base_root, output_root / "source_audits" / "base_44")
    copy_source_audits(addition_root, output_root / "source_audits" / "additional_11")

    cohort_counts = Counter(row["cohort"] for row in audits)
    patient_ids = {
        cohort: sorted(row["patient_id"] for row in audits if row["cohort"] == cohort)
        for cohort in sorted(cohort_counts)
    }
    index_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []
    for row in sorted(audits, key=lambda item: item["patient_id"]):
        manifest = row["manifest"]
        mngs = row["mngs"]
        index_rows.append({
            "cohort": row["cohort"],
            "hospital": manifest.get("hospital") or mngs.get("hospital") or "",
            "patient_id": row["patient_id"],
            "patient_key": row["patient_key"],
            "hospital_case_id": manifest.get("hospital_case_id") or mngs.get("hospital_case_id") or "",
            "hospital_patient_id": manifest.get("hospital_patient_id") or mngs.get("hospital_patient_id") or "",
            "mngs_specimen_codes": ";".join(row["specimen_codes"]),
            "mngs_collected_times": ";".join(row["collected_times"]),
            "direct_nonempty_module_count": row["direct_nonempty"],
            "direct_empty_module_count": row["direct_empty"],
            "direct_missing_module_count": row["direct_missing"],
            "raw_export_available_count": row["raw_export_available_count"],
            "mngs_case_count": row["mngs_case_count"],
            "mngs_test_count": row["mngs_test_count"],
            "mngs_organism_count": row["mngs_organism_count"],
            "mngs_measurement_count": row["mngs_measurement_count"],
            "positive_measurement_count": row["positive_measurement_count"],
            "ground_truth_status": "present_nonempty" if row["gt_file"] else "missing",
            "answer_labels_included": True,
        })
        for module, status in row["module_status"].items():
            if status not in {"present_nonempty"}:
                missing_rows.append({
                    "cohort": row["cohort"],
                    "patient_id": row["patient_id"],
                    "category": "pipeline_direct_clinical_input",
                    "module": module,
                    "status": status,
                    "interpretation": "explicit_empty_unknown_not_negative" if status == "present_empty" else "unknown_not_negative",
                })
        if row["raw_export_available_count"] == 0:
            missing_rows.append({
                "cohort": row["cohort"],
                "patient_id": row["patient_id"],
                "category": "upstream_raw_clinical_export",
                "module": "all",
                "status": "not_in_source_bundle",
                "interpretation": "optional; layer_02 remains usable",
            })

    write_csv(
        output_root / "PATIENT_INDEX.csv",
        list(index_rows[0]),
        index_rows,
    )
    write_csv(
        output_root / "MISSING_DATA.csv",
        ["cohort", "patient_id", "category", "module", "status", "interpretation"],
        missing_rows,
    )

    base_gt = load_json(base_root / "GROUND_TRUTH_MANIFEST.json")
    addition_gt = load_json(addition_root / "GROUND_TRUTH_MANIFEST.json")
    gt_patients = list(base_gt.get("patients") or [])
    for item in addition_gt.get("patients") or []:
        item = dict(item)
        item.setdefault("cohort", "two_hospital_additional_11")
        item.setdefault("source_file", addition_gt.get("source_file"))
        item.setdefault("sheet", addition_gt.get("sheet"))
        gt_patients.append(item)
    gt_manifest = {
        "schema_version": "ground_truth_manifest.v2",
        "created_date_asia_taipei": "2026-10-09",
        "patient_count": len(gt_patients),
        "patients_with_nonempty_labels": len(gt_patients),
        "evaluation_only": True,
        "inference_input": False,
        "answer_firewall": "layer_05 is excluded by the answer-blind workspace builder",
        "sources": {
            "base_44": str((base_root / "GROUND_TRUTH_MANIFEST.json").resolve()),
            "additional_11": str((addition_root / "GROUND_TRUTH_MANIFEST.json").resolve()),
        },
        "patients": sorted(gt_patients, key=lambda item: int(item["patient_id"])),
    }
    write_json(output_root / "GROUND_TRUTH_MANIFEST.json", gt_manifest)

    bundle_manifest = {
        "schema_version": "combined_pipeline_input_bundle.v1",
        "created_date_asia_taipei": "2026-10-09",
        "package_root": output_root.name,
        "patient_count": len(audits),
        "cohort_counts": dict(sorted(cohort_counts.items())),
        "patient_ids": patient_ids,
        "answer_labels_included": True,
        "ground_truth": {"evaluation_only": True, "inference_input": False, "layer": "05_ground_truth"},
        "clinical_data_policy": "non-mNGS clinical modules are optional; missing/empty means unknown, not negative",
        "mngs_policy": "per-test DNA/RNA/RK/NTC data; positive candidates require reads > 0 and detected=true; reads are not summed across tests",
        "time_window_policy": "inclusive plus/minus two calendar dates for temporal clinical evidence",
        "source_bundles": {
            "base_44": str(base_root.resolve()),
            "additional_11": str(addition_root.resolve()),
        },
        "source_hash_validation": {"base_44": base_hashes, "additional_11": addition_hashes},
        "reference_outputs_included": False,
    }
    write_json(output_root / "bundle_manifest.json", bundle_manifest)

    validation = {
        "schema_version": "combined_pipeline_input_validation.v1",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "passed",
        "source_hash_validation": {"base_44": base_hashes, "additional_11": addition_hashes},
        "patient_count": len(audits),
        "cohort_counts": dict(sorted(cohort_counts.items())),
        "patient_id_overlap": overlaps,
        "patient_error_count": 0,
        "patient_warning_count": sum(len(row["warnings"]) for row in audits),
        "patient_warnings": [
            {"cohort": row["cohort"], "patient_id": row["patient_id"], "warnings": row["warnings"]}
            for row in audits if row["warnings"]
        ],
        "mngs_totals": {
            "cases": sum(row["mngs_case_count"] for row in audits),
            "tests": sum(row["mngs_test_count"] for row in audits),
            "organisms": sum(row["mngs_organism_count"] for row in audits),
            "measurements": sum(row["mngs_measurement_count"] for row in audits),
            "positive_measurements": sum(row["positive_measurement_count"] for row in audits),
            "ignored_zero_read_placeholders": sum(row["ignored_zero_read_placeholders"] for row in audits),
            "raw_unavailable_measurements": sum(row["unavailable_measurements"] for row in audits),
            "invalid_nonpositive_measurements": sum(row["invalid_nonpositive_measurements"] for row in audits),
            "tests_missing_input_reads": sum(row["missing_input_reads"] for row in audits),
            "tests_missing_host_ratio": sum(row["missing_host_ratio"] for row in audits),
            "cases_missing_specimen_site": sum(row["missing_specimen_site"] for row in audits),
            "cases_missing_collected_time": sum(row["missing_collected_time"] for row in audits),
            "tests_missing_nucleic_type": sum(row["missing_nucleic_type"] for row in audits),
            "tests_missing_report_time": sum(row["missing_report_time"] for row in audits),
        },
        "normalization_gap_patients": [
            {
                "cohort": row["cohort"],
                "patient_id": row["patient_id"],
                "tests_missing_input_reads": row["missing_input_reads"],
                "tests_missing_host_ratio": row["missing_host_ratio"],
                "cases_missing_specimen_site": row["missing_specimen_site"],
                "cases_missing_collected_time": row["missing_collected_time"],
                "tests_missing_nucleic_type": row["missing_nucleic_type"],
                "tests_missing_report_time": row["missing_report_time"],
            }
            for row in audits
            if any(
                row[key]
                for key in (
                    "missing_input_reads",
                    "missing_host_ratio",
                    "missing_specimen_site",
                    "missing_collected_time",
                    "missing_nucleic_type",
                    "missing_report_time",
                )
            )
        ],
        "clinical_module_status_counts": {
            cohort: dict(sorted(Counter(
                status
                for row in audits if row["cohort"] == cohort
                for status in row["module_status"].values()
            ).items()))
            for cohort in sorted(cohort_counts)
        },
        "layer_04_files": 0,
        "ground_truth_files": len([row for row in audits if row["gt_file"]]),
        "ground_truth_inference_input": False,
    }
    write_json(output_root / "INTEGRATION_VALIDATION.json", validation)

    readme = "# 55-patient pipeline input bundle\n\n"
    readme += "This bundle combines KH33, the original two-hospital 11-patient cohort, and the additional two-hospital 11-patient cohort.\n\n"
    readme += "## Contract\n\n"
    readme += "- Layer 02 contains normalized clinical inputs. Clinical modules are optional; missing or empty means unknown, not negative.\n"
    readme += "- Layer 03 contains mandatory per-test mNGS DNA/RNA/RK/NTC records. Only measurements with reads > 0 and detected=true enter the candidate path.\n"
    readme += "- Cross-test zero-read placeholders and raw-unavailable null measurements may remain for alignment/audit, but are excluded from every positive candidate set.\n"
    readme += "- Layer 04 must remain empty; no prior model output is included.\n"
    readme += "- Layer 05 contains physician ground truth for evaluation only. It must never be copied into an inference workspace.\n"
    readme += "- Temporal clinical evidence uses the inclusive previous/next two calendar-date policy.\n"
    readme += "- Reads from different tests are not summed, and technical repeats are not treated as independent clinical events.\n\n"
    readme += "## Validation\n\n"
    readme += f"- Patients: {len(audits)} ({', '.join(f'{key}={value}' for key, value in sorted(cohort_counts.items()))})\n"
    readme += f"- mNGS tests: {validation['mngs_totals']['tests']}\n"
    readme += f"- Positive mNGS measurements: {validation['mngs_totals']['positive_measurements']}\n"
    readme += f"- Ignored zero-read alignment placeholders: {validation['mngs_totals']['ignored_zero_read_placeholders']}\n"
    readme += f"- Raw-unavailable null measurements: {validation['mngs_totals']['raw_unavailable_measurements']}\n"
    readme += "- Invalid non-positive measurements: 0\n"
    readme += "- Prior derived outputs: 0 files\n"
    readme += "- Ground-truth firewall: enforced by the answer-blind workspace builder and validator\n\n"
    readme += "Original root manifests are retained under `source_audits/` for provenance.\n"
    (output_root / "README.md").write_text(readme, encoding="utf-8")

    file_rows: list[dict[str, Any]] = []
    for path in sorted(output_root.rglob("*")):
        if not path.is_file() or path.name == "SHA256SUMS.csv":
            continue
        relative = path.relative_to(output_root).as_posix()
        parts = path.relative_to(output_root).parts
        cohort = parts[1] if len(parts) > 1 and parts[0] == "patients" else ""
        patient_id = ""
        if len(parts) > 2 and parts[0] == "patients":
            match = PATIENT_RE.fullmatch(parts[2])
            patient_id = match.group(1) if match else ""
        file_rows.append({
            "cohort": cohort,
            "patient_id": patient_id,
            "relative_path": relative,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    write_csv(
        output_root / "FILE_MANIFEST.csv",
        ["cohort", "patient_id", "relative_path", "bytes", "sha256"],
        file_rows,
    )

    sha_rows: list[dict[str, Any]] = []
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.csv":
            sha_rows.append({
                "relative_path": path.relative_to(output_root).as_posix(),
                "sha256": sha256_file(path),
            })
    write_csv(output_root / "SHA256SUMS.csv", ["relative_path", "sha256"], sha_rows)
    return validation


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-root", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--addition-root", type=Path, default=DEFAULT_ADDITION)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    result = build(args.base_root, args.addition_root, args.output_root)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

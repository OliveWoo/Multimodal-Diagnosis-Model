"""Build an auditable enriched shadow for extracted workbook JSON.

The source staging tree is never modified.  For the 33 KH patients with the
current multi-assay mNGS bundle, this tool:

* audits same-number identity with date/item/value fingerprints;
* copies the complete per-test DNA/RNA mNGS bundle and ranked candidates;
* fills only empty compatible clinical modules from the legacy KH folder;
* preserves all legacy clinical modules under ``legacy_raw_sources``;
* writes patient-level and cohort-level provenance and missing-data reports.

Non-empty v14 clinical modules are not concatenated with legacy rows because
the two exports may contain the same observation under different structures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STAGING = ROOT / "extracted_workbooks_staging_v14"
DEFAULT_LEGACY = ROOT / "outputs" / "patient_info_KH_0728_2Days"
DEFAULT_MULTI_ASSAY = ROOT / "52_patients_all_mNGS_DNA_RNA_RK_NTC_by_species_no_zero_reads"
DEFAULT_OUTPUT = ROOT / "outputs" / "extracted_workbooks_enriched_shadow_v4_52mngs"

CLINICAL_MODULES = (
    "CBC",
    "other_lab",
    "culture",
    "filmarray",
    "gm_test",
    "molecular_microbiology",
    "image",
    "underlying",
    "admission_diagnosis",
)
FILL_IF_EMPTY_MODULES = CLINICAL_MODULES
FINGERPRINT_MODULES = ("CBC", "other_lab")
MICROBIOLOGY_MODULES = ("culture", "filmarray", "gm_test", "molecular_microbiology")


@dataclass(frozen=True)
class PatientSources:
    patient_id: int
    staging_dir: Path
    legacy_dir: Path | None
    multi_assay_file: Path | None


def patient_id_from_name(name: str) -> int | None:
    match = re.search(r"(?:patient|extracted_patient|NGS_patient)_(\d+)", name, re.I)
    return int(match.group(1)) if match else None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path | None) -> Any:
    if path is None or not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def as_records(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [row for row in payload if isinstance(row, dict)]
    if isinstance(payload, dict):
        records = payload.get("records")
        if isinstance(records, list):
            return [row for row in records if isinstance(row, dict)]
        return [payload]
    return []


def is_empty_json(path: Path) -> bool:
    payload = load_json(path)
    return payload in (None, [], {})


def normalize_text(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def normalize_value(value: Any) -> str:
    text = str(value or "").strip().replace(",", "")
    if not text:
        return ""
    try:
        number = Decimal(text)
    except InvalidOperation:
        return normalize_text(text)
    normalized = format(number.normalize(), "f")
    return "0" if normalized in {"-0", ""} else normalized


def date_part(value: Any) -> str:
    text = str(value or "").strip()
    match = re.search(r"(20\d{2})[-/](\d{1,2})[-/](\d{1,2})", text)
    if not match:
        return ""
    return f"{int(match.group(1)):04d}-{int(match.group(2)):02d}-{int(match.group(3)):02d}"


def candidate_dates(row: dict[str, Any]) -> set[str]:
    dates = {
        date_part(row.get(key))
        for key in ("collected_time", "reported_time", "received_time", "source_time")
    }
    source = row.get("source") if isinstance(row.get("source"), dict) else {}
    original = source.get("original_date_fields") if isinstance(source, dict) else {}
    if isinstance(original, dict):
        dates.update(date_part(value) for value in original.values())
    return {value for value in dates if value}


def find_module_file(directory: Path | None, patient_id: int, module: str) -> Path | None:
    if directory is None:
        return None
    exact_names = (
        f"extracted_patient_{patient_id}_{module}.json",
        f"NGS_patient_{patient_id}_{module}.json",
    )
    for name in exact_names:
        path = directory / name
        if path.exists():
            return path
    matches = sorted(directory.glob(f"*_{module}.json"))
    return matches[0] if matches else None


def lab_fingerprint(directory: Path | None, patient_id: int) -> dict[str, Counter[str]]:
    dated: Counter[str] = Counter()
    undated: Counter[str] = Counter()
    for module in FINGERPRINT_MODULES:
        path = find_module_file(directory, patient_id, module)
        for row in as_records(load_json(path)):
            item = normalize_text(row.get("item") or row.get("item_full") or row.get("test"))
            value = normalize_value(row.get("value") if "value" in row else row.get("result"))
            if not item or not value:
                continue
            token = f"{item}|{value}"
            undated[token] += 1
            for date in candidate_dates(row):
                dated[f"{date}|{token}"] += 1
    return {"dated": dated, "undated": undated}


def microbiology_fingerprint(directory: Path | None, patient_id: int) -> Counter[str]:
    output: Counter[str] = Counter()
    for module in MICROBIOLOGY_MODULES:
        path = find_module_file(directory, patient_id, module)
        for row in as_records(load_json(path)):
            identity = normalize_text(
                row.get("organism") or row.get("target") or row.get("test") or row.get("item")
            )
            result = normalize_text(row.get("status") or row.get("result") or row.get("value"))
            if not identity:
                continue
            base = f"{identity}|{result}"
            dates = candidate_dates(row)
            if dates:
                for date in dates:
                    output[f"{date}|{base}"] += 1
            else:
                output[f"undated|{base}"] += 1
    return output


def intersection_count(left: Counter[str], right: Counter[str]) -> int:
    return sum((left & right).values())


def linkage_score(
    legacy_fp: dict[str, Counter[str]], staging_fp: dict[str, Counter[str]],
    legacy_micro: Counter[str], staging_micro: Counter[str],
) -> tuple[float, int, int, int]:
    dated = intersection_count(legacy_fp["dated"], staging_fp["dated"])
    undated = intersection_count(legacy_fp["undated"], staging_fp["undated"])
    micro = intersection_count(legacy_micro, staging_micro)
    score = dated * 4.0 + micro * 3.0 + undated * 0.15
    return score, dated, undated, micro


def build_linkage_audit(
    legacy_dirs: dict[int, Path], staging_dirs: dict[int, Path]
) -> list[dict[str, Any]]:
    legacy_labs = {pid: lab_fingerprint(path, pid) for pid, path in legacy_dirs.items()}
    staging_labs = {pid: lab_fingerprint(path, pid) for pid, path in staging_dirs.items()}
    legacy_micro = {pid: microbiology_fingerprint(path, pid) for pid, path in legacy_dirs.items()}
    staging_micro = {pid: microbiology_fingerprint(path, pid) for pid, path in staging_dirs.items()}
    rows: list[dict[str, Any]] = []
    for pid in sorted(legacy_dirs):
        comparisons: list[dict[str, Any]] = []
        for staging_pid in sorted(staging_dirs):
            score, dated, undated, micro = linkage_score(
                legacy_labs[pid], staging_labs[staging_pid],
                legacy_micro[pid], staging_micro[staging_pid],
            )
            comparisons.append(
                {
                    "staging_patient_id": staging_pid,
                    "score": round(score, 3),
                    "exact_dated_lab_matches": dated,
                    "item_value_matches": undated,
                    "dated_microbiology_matches": micro,
                }
            )
        comparisons.sort(key=lambda row: (-row["score"], row["staging_patient_id"]))
        same = next(row for row in comparisons if row["staging_patient_id"] == pid)
        same_rank = next(index + 1 for index, row in enumerate(comparisons) if row["staging_patient_id"] == pid)
        best = comparisons[0]
        second = comparisons[1] if len(comparisons) > 1 else None
        locally_supported = same_rank == 1 and (
            same["exact_dated_lab_matches"] >= 3 or same["dated_microbiology_matches"] >= 1
        )
        rows.append(
            {
                "patient_id": pid,
                "user_asserted_same_number": True,
                "same_number_local_status": (
                    "supported_by_local_fingerprint" if locally_supported
                    else "user_asserted_requires_manual_confirmation"
                ),
                "same_number_rank": same_rank,
                "same_number_score": same["score"],
                "same_number_exact_dated_lab_matches": same["exact_dated_lab_matches"],
                "same_number_item_value_matches": same["item_value_matches"],
                "same_number_dated_microbiology_matches": same["dated_microbiology_matches"],
                "best_match_patient_id": best["staging_patient_id"],
                "best_match_score": best["score"],
                "second_match_patient_id": second["staging_patient_id"] if second else None,
                "second_match_score": second["score"] if second else None,
                "same_number_margin_over_best_other": round(
                    same["score"] - max(
                        (row["score"] for row in comparisons if row["staging_patient_id"] != pid),
                        default=0.0,
                    ),
                    3,
                ),
            }
        )
    return rows


def discover_sources(staging_root: Path, legacy_root: Path, multi_assay_root: Path) -> list[PatientSources]:
    staging_dirs = {
        pid: path
        for path in staging_root.glob("extracted_patient_*_json")
        if path.is_dir() and (pid := patient_id_from_name(path.name)) is not None
    }
    legacy_dirs = {
        pid: path
        for path in legacy_root.glob("NGS_patient_*_json")
        if path.is_dir() and (pid := patient_id_from_name(path.name)) is not None
    }
    multi_files = {
        pid: path
        for path in multi_assay_root.glob("NGS_patient_*_all_mNGS_DNA_RNA_RK_NTC_by_species.json")
        if (pid := patient_id_from_name(path.name)) is not None
    }
    return [
        PatientSources(pid, staging_dirs[pid], legacy_dirs.get(pid), multi_files.get(pid))
        for pid in sorted(staging_dirs)
    ]


def copy_with_provenance(source: Path, destination: Path) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return {
        "source": str(source.resolve()),
        "destination": str(destination.resolve()),
        "sha256": sha256_file(source),
        "bytes": source.stat().st_size,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def update_source_mapping(patient_dir: Path, additions: dict[str, list[str]]) -> None:
    mapping_path = patient_dir / "source_mapping.json"
    mapping = load_json(mapping_path)
    if not isinstance(mapping, dict):
        mapping = {}
    mapping.update(additions)
    write_json(mapping_path, mapping)


def build_patient_shadow(
    source: PatientSources,
    output_root: Path,
    linkage_row: dict[str, Any] | None = None,
) -> dict[str, Any]:
    pid = source.patient_id
    destination = output_root / source.staging_dir.name
    shutil.copytree(source.staging_dir, destination)
    copied: list[dict[str, Any]] = []
    filled_modules: list[str] = []
    retained_v14_modules: list[str] = []
    legacy_available: list[str] = []

    if source.legacy_dir is not None:
        # The analyzer skips directories containing "raw".  This prevents the
        # preserved audit copy from being mistaken for another patient input.
        legacy_dir = destination / "legacy_raw_sources"
        for module in CLINICAL_MODULES:
            old_file = find_module_file(source.legacy_dir, pid, module)
            target_file = find_module_file(destination, pid, module)
            if old_file is None:
                continue
            legacy_available.append(module)
            copied.append(copy_with_provenance(old_file, legacy_dir / old_file.name))
            if (
                module in FILL_IF_EMPTY_MODULES
                and target_file is not None
                and is_empty_json(target_file)
                and not is_empty_json(old_file)
            ):
                payload = load_json(old_file)
                write_json(target_file, payload)
                filled_modules.append(module)
            elif target_file is not None and not is_empty_json(target_file):
                retained_v14_modules.append(module)

    mapping_additions: dict[str, list[str]] = {}
    if source.multi_assay_file is not None:
        name = f"extracted_patient_{pid}_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
        copied.append(copy_with_provenance(source.multi_assay_file, destination / name))
        mapping_additions["mNGS_multi_assay"] = [name]

    ranked_source = (
        source.legacy_dir / f"NGS_patient_{pid}_mNGS_ranked_candidates.json"
        if source.legacy_dir is not None else None
    )
    if ranked_source is not None and ranked_source.exists():
        name = f"extracted_patient_{pid}_mNGS_ranked_candidates.json"
        copied.append(copy_with_provenance(ranked_source, destination / name))
        mapping_additions["mNGS_ranked_candidates"] = [name]

    if source.legacy_dir is not None:
        mapping_additions["legacy_raw_sources"] = ["legacy_raw_sources"]
    if mapping_additions:
        update_source_mapping(destination, mapping_additions)

    manifest = {
        "schema_version": "extracted_workbooks_enriched_shadow_v4",
        "patient_id": pid,
        "source_staging_dir": str(source.staging_dir.resolve()),
        "source_staging_manifest_sha256": (
            sha256_file(source.staging_dir / "extraction_manifest.json")
            if (source.staging_dir / "extraction_manifest.json").exists() else None
        ),
        "same_number_linkage_policy": "user_asserted_with_local_fingerprint_audit",
        "same_number_linkage_audit": linkage_row,
        "multi_assay_available": source.multi_assay_file is not None,
        "legacy_patient_available": source.legacy_dir is not None,
        "filled_only_when_v14_empty": filled_modules,
        "retained_nonempty_v14_modules": sorted(set(retained_v14_modules)),
        "legacy_modules_preserved_separately": sorted(set(legacy_available)),
        "copied_sources": copied,
        "mngs_grouped_remains_unmodified": True,
        "ready_for_full_pipeline": False,
        "blocking_reason": (
            None if source.multi_assay_file is not None and source.legacy_dir is not None
            else "missing_current_multi_assay_or_legacy_patient_source"
        ),
    }
    write_json(destination / "enrichment_manifest.json", manifest)
    return manifest


def missing_report_rows(sources: Iterable[PatientSources], output_root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for source in sources:
        pid = source.patient_id
        patient_dir = output_root / source.staging_dir.name
        missing_modules = []
        for module in CLINICAL_MODULES:
            path = find_module_file(patient_dir, pid, module)
            if path is None or is_empty_json(path):
                missing_modules.append(module)
        rows.append(
            {
                "patient_id": pid,
                "has_current_multi_assay_mngs": source.multi_assay_file is not None,
                "has_legacy_kh_patient": source.legacy_dir is not None,
                "has_ranked_candidates": (
                    patient_dir / f"extracted_patient_{pid}_mNGS_ranked_candidates.json"
                ).exists(),
                "missing_clinical_modules_after_fill": "|".join(missing_modules),
                "full_pipeline_source_status": (
                    "source_bundle_present_needs_adapter_validation"
                    if source.multi_assay_file is not None and source.legacy_dir is not None
                    else "missing_source_bundle"
                ),
            }
        )
    return rows


def write_markdown_report(
    path: Path, linkage_rows: list[dict[str, Any]], coverage_rows: list[dict[str, Any]]
) -> None:
    supported = [r["patient_id"] for r in linkage_rows if r["same_number_local_status"] == "supported_by_local_fingerprint"]
    review = [r["patient_id"] for r in linkage_rows if r["same_number_local_status"] != "supported_by_local_fingerprint"]
    current = [r["patient_id"] for r in coverage_rows if r["has_current_multi_assay_mngs"]]
    missing = [r["patient_id"] for r in coverage_rows if not r["has_current_multi_assay_mngs"]]
    current_gaps = [
        r for r in coverage_rows
        if r["has_current_multi_assay_mngs"] and r["missing_clinical_modules_after_fill"]
    ]
    missing_source_rows = [r for r in coverage_rows if not r["has_current_multi_assay_mngs"]]
    lines = [
        "# Extracted workbooks enriched shadow v4 (52-patient mNGS source)",
        "",
        "This output does not modify `extracted_workbooks_staging_v14`.",
        "",
        "## Same-number linkage",
        "",
        f"- User-confirmed same-number KH patients: {len(linkage_rows)}.",
        f"- Locally supported by dated lab or microbiology fingerprints: {len(supported)} ({', '.join('P'+str(x) for x in supported) or 'none'}).",
        f"- Still requiring manual linkage confirmation: {len(review)} ({', '.join('P'+str(x) for x in review) or 'none'}).",
        "- Local fingerprints support patient identity only; they do not prove the same infection episode.",
        "",
        "## Current multi-assay mNGS coverage",
        "",
        f"- Present: {len(current)} ({', '.join('P'+str(x) for x in current)}).",
        f"- Missing: {len(missing)} ({', '.join('P'+str(x) for x in missing)}).",
        "",
        "## Remaining clinical gaps in the 33-patient cohort",
        "",
        *(
            [
                f"- P{row['patient_id']}: {row['missing_clinical_modules_after_fill'].replace('|', ', ')}"
                for row in current_gaps
            ]
            or ["- None."]
        ),
        "- Missing means no usable record was found; it must not be interpreted as a negative result.",
        "",
        "## Staging patients without the KH source bundle",
        "",
        *[
            (
                f"- P{row['patient_id']}: current multi-assay mNGS, ranked candidates; "
                f"clinical gaps: {row['missing_clinical_modules_after_fill'].replace('|', ', ') or 'none'}"
            )
            for row in missing_source_rows
        ],
        "",
        "## Merge policy",
        "",
        "- Complete per-test DNA/RNA mNGS and ranked candidates are copied as separate traced files.",
        "- Empty v14 clinical modules may be filled from the same-number legacy KH source.",
        "- Non-empty v14 clinical modules are retained; legacy versions are stored under `legacy_raw_sources` and are not concatenated.",
        "- The preserved legacy directory contains `raw`, so the clinical analyzer will not mistake it for a patient input.",
        "- `mNGS_grouped.json` is left unchanged because its old schema differs from the current per-test multi-assay schema.",
        "- All enriched patients remain shadow inputs until the compatibility adapter and event-linkage checks pass.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(staging_root: Path, legacy_root: Path, multi_assay_root: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output root: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    sources = discover_sources(staging_root, legacy_root, multi_assay_root)
    staging_dirs = {source.patient_id: source.staging_dir for source in sources}
    all_multi_assay_ids = {
        pid
        for path in multi_assay_root.glob(
            "NGS_patient_*_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
        )
        if (pid := patient_id_from_name(path.name)) is not None
    }
    multi_assay_without_staging = sorted(all_multi_assay_ids - set(staging_dirs))
    legacy_dirs = {source.patient_id: source.legacy_dir for source in sources if source.legacy_dir is not None}
    linkage_rows = build_linkage_audit(legacy_dirs, staging_dirs)
    linkage_by_patient = {row["patient_id"]: row for row in linkage_rows}
    manifests = [
        build_patient_shadow(source, output_root, linkage_by_patient.get(source.patient_id))
        for source in sources
    ]
    coverage_rows = missing_report_rows(sources, output_root)
    write_csv(output_root / "identity_linkage_audit.csv", linkage_rows)
    write_csv(output_root / "patient_coverage.csv", coverage_rows)
    write_markdown_report(output_root / "MISSING_DATA_REPORT.md", linkage_rows, coverage_rows)

    summary = {
        "schema_version": "extracted_workbooks_enriched_shadow_summary_v4",
        "source_staging_root": str(staging_root.resolve()),
        "legacy_patient_root": str(legacy_root.resolve()),
        "multi_assay_root": str(multi_assay_root.resolve()),
        "output_root": str(output_root.resolve()),
        "staging_patients": len(sources),
        "legacy_patients": sum(bool(source.legacy_dir) for source in sources),
        "multi_assay_patients": sum(bool(source.multi_assay_file) for source in sources),
        "multi_assay_files_found": len(all_multi_assay_ids),
        "multi_assay_without_staging_patient": multi_assay_without_staging,
        "locally_supported_same_number": sum(
            row["same_number_local_status"] == "supported_by_local_fingerprint"
            for row in linkage_rows
        ),
        "same_number_manual_confirmation_required": [
            row["patient_id"] for row in linkage_rows
            if row["same_number_local_status"] != "supported_by_local_fingerprint"
        ],
        "patients_missing_multi_assay": [
            row["patient_id"] for row in coverage_rows if not row["has_current_multi_assay_mngs"]
        ],
        "current_cohort_clinical_gaps": {
            f"P{row['patient_id']}": row["missing_clinical_modules_after_fill"].split("|")
            for row in coverage_rows
            if row["has_current_multi_assay_mngs"] and row["missing_clinical_modules_after_fill"]
        },
        "filled_module_counts": dict(
            Counter(module for manifest in manifests for module in manifest["filled_only_when_v14_empty"])
        ),
        "source_staging_modified": False,
        "ready_for_full_pipeline": False,
        "next_required_step": "validate compatibility adapter and infection-event linkage",
    }
    write_json(output_root / "enrichment_summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging-root", type=Path, default=DEFAULT_STAGING)
    parser.add_argument("--legacy-root", type=Path, default=DEFAULT_LEGACY)
    parser.add_argument("--multi-assay-root", type=Path, default=DEFAULT_MULTI_ASSAY)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = run(args.staging_root, args.legacy_root, args.multi_assay_root, args.output_root)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

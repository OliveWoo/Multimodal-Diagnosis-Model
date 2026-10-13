"""Adapt the legacy two-hospital artifacts to the current v2.6 input contract.

The adapter is deliberately loss-aware.  It never invents input-read
denominators, specimen site, collection time, molecule identity, or technical
repeat independence.  Candidate rows retain their original seq_id, reads,
taxid and NTC/RK tags.  Tests with unknown molecule identity are explicitly
marked ``source_unknown`` so the current pipeline can keep the signal visible
without assigning DNA/RNA concordance or repeat credit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "outputs/patient_info_two_hospitals_0611_rich_moderate"
DEFAULT_OUTPUT = (
    ROOT
    / "outputs/runs/2026-10-06_two_hospital_v2_6_compatibility_input_v1"
)
DEFAULT_SUMMARY_SUFFIX = "final_summary_with_underlying_with_filmarray"
PATIENT_RE = re.compile(r"^NGS_patient_(\d+)_json$")
CATEGORIES = ("bacterial", "viral", "fungal", "others")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def infer_molecule(source_file: Any) -> tuple[str, str]:
    text = str(source_file or "").upper()
    has_dna = "DNA" in text
    has_rna = "RNA" in text
    if has_dna and not has_rna:
        return "DNA", "inferred_from_source_path"
    if has_rna and not has_dna:
        return "RNA", "inferred_from_source_path"
    return "UNKNOWN", "source_unknown"


def category_name(value: Any) -> str:
    text = str(value or "").lower()
    if "bac" in text:
        return "bacterial"
    if "virus" in text or "viral" in text:
        return "viral"
    if "fung" in text:
        return "fungal"
    return "others"


def numeric_reads(row: dict[str, Any]) -> float:
    source = row.get("row") if isinstance(row.get("row"), dict) else {}
    try:
        value = float(str(source.get("Sec.hit") or "").strip())
    except ValueError as error:
        raise ValueError(
            f"Missing/invalid Sec.hit for {row.get('organism_name')!r}"
        ) from error
    if value <= 0:
        raise ValueError(
            f"Nonpositive retained Sec.hit for {row.get('organism_name')!r}: {value}"
        )
    return value


def selection_tags(row: dict[str, Any]) -> list[str]:
    tags: list[str] = []
    for value in row.get("selection_reasons") or []:
        tag = str(value or "").strip().replace("RK00", "RK0").replace("RK08", "RK8")
        if tag and tag not in tags:
            tags.append(tag)
    source = row.get("row") if isinstance(row.get("row"), dict) else {}
    if str(source.get("Code_NTC") or "").strip() == "00" and "Code_NTC=00" not in tags:
        tags.append("Code_NTC=00")
    rk_ntc = str(source.get("Code_RK_NTC") or "").strip().upper().replace("-", "")
    if rk_ntc in {"RK0", "RK00"} and "Code_RK_NTC=RK0" not in tags:
        tags.append("Code_RK_NTC=RK0")
    elif rk_ntc in {"RK8", "RK08"} and "Code_RK_NTC=RK8" not in tags:
        tags.append("Code_RK_NTC=RK8")
    if not tags:
        tags.append("TWO_HOSPITAL_LEGACY_RETAINED_ROW")
    return tags


def ranked_metadata(patient_dir: Path, patient_id: int) -> dict[str, dict[str, str]]:
    path = patient_dir / f"NGS_patient_{patient_id}_mNGS_ranked_candidates.json"
    if not path.is_file():
        return {}
    payload = read_json(path)
    result: dict[str, dict[str, str]] = {}
    for record in payload.get("records") or []:
        if not isinstance(record, dict):
            continue
        seq_id = str(record.get("seq_id") or record.get("specimen_code") or "").strip()
        if not seq_id:
            continue
        result[seq_id] = {
            "collected_time": str(record.get("collected_time") or "").strip(),
            "specimen_site": str(record.get("specimen_site") or "").strip(),
            "specimen_code": str(record.get("specimen_code") or seq_id).strip(),
        }
    return result


def build_multi_assay_document(
    patient_dir: Path, patient_id: int, raw_rows: list[dict[str, Any]]
) -> dict[str, Any]:
    metadata = ranked_metadata(patient_dir, patient_id)
    rows_by_case: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        case_id = str(row.get("specimen_id") or row.get("seq_id") or "").strip()
        if not case_id:
            raise ValueError(f"Patient {patient_id} has a raw row without specimen/seq identity")
        rows_by_case[case_id].append(row)

    cases: list[dict[str, Any]] = []
    for case_id, case_rows in sorted(rows_by_case.items()):
        seq_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in case_rows:
            seq_id = str(row.get("seq_id") or "").strip()
            if not seq_id:
                raise ValueError(f"Patient {patient_id} has a raw row without seq_id")
            seq_rows[seq_id].append(row)
        tests: list[dict[str, Any]] = []
        case_times: list[str] = []
        case_sites: list[str] = []
        specimen_codes: list[str] = []
        for index, (seq_id, values) in enumerate(sorted(seq_rows.items()), start=1):
            molecule, molecule_status = infer_molecule(values[0].get("source_file"))
            present_categories = {category_name(value.get("source_category")) for value in values}
            meta = metadata.get(seq_id, {})
            if meta.get("collected_time"):
                case_times.append(meta["collected_time"])
            if meta.get("specimen_site"):
                case_sites.append(meta["specimen_site"])
            if meta.get("specimen_code"):
                specimen_codes.append(meta["specimen_code"])
            tests.append({
                "test_number": index,
                "seq_id": seq_id,
                "chip_id": "",
                "nucleic_type": molecule,
                "molecule_status": molecule_status,
                "start_time": "",
                "report_time": "",
                "parameters": {
                    "condition": "-",
                    "specimen_status": "",
                    "specimen_type": meta.get("specimen_site") or "",
                    "input_reads": None,
                    "host_ratio": None,
                    "source_contract_note": "two_hospital_legacy_denominators_absent",
                },
                "raw_result_status": {
                    category: "available" if category in present_categories else "unknown"
                    for category in CATEGORIES
                },
            })
        cases.append({
            "case_review_id": case_id,
            "specimen_code": specimen_codes[0] if specimen_codes else case_id,
            "specimen_code_source": "legacy_ranked_artifact_or_raw_specimen_id",
            "collected_time": case_times[0] if case_times else "",
            "specimen_site": case_sites[0] if case_sites else "",
            "case_review": {
                "source": "two_hospital_legacy_flattened_candidate_artifact",
                "missingness_note": "Blank collection/site fields are source unknowns, not negative evidence.",
            },
            "tests": tests,
        })

    organism_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        name = str(row.get("organism_name") or "").strip()
        if not name:
            continue
        organism_groups[(category_name(row.get("source_category")), name)].append(row)
    organisms: list[dict[str, Any]] = []
    for (category, name), values in sorted(
        organism_groups.items(), key=lambda item: (item[0][0], item[0][1].lower())
    ):
        measurements = []
        # The legacy export can contain the same physical seq_id/result more
        # than once because the same AJed archive was discovered through
        # different source folders.  A pipeline observation is defined at the
        # case + seq_id + organism level, so those export duplicates must not
        # masquerade as independent tests.  Preserve all source paths for
        # audit, and use the maximum reported reads if duplicated exports
        # disagree (the disagreement remains explicit below).
        observation_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
        for row in values:
            observation_groups[(
                str(row.get("specimen_id") or row.get("seq_id") or "").strip(),
                str(row.get("seq_id") or "").strip(),
            )].append(row)
        for (case_review_id, seq_id), duplicates in sorted(observation_groups.items()):
            ranked = sorted(
                duplicates,
                key=lambda row: (
                    -numeric_reads(row),
                    str(row.get("source_file") or ""),
                    int(row.get("row_index") or 0),
                ),
            )
            row = ranked[0]
            source = row.get("row") if isinstance(row.get("row"), dict) else {}
            read_values = sorted({numeric_reads(item) for item in duplicates})
            taxid_values = sorted({
                str((item.get("row") or {}).get("taxid") or "").strip()
                for item in duplicates
                if isinstance(item.get("row"), dict)
            })
            merged_tags: list[str] = []
            for item in duplicates:
                for tag in selection_tags(item):
                    if tag not in merged_tags:
                        merged_tags.append(tag)
            measurements.append({
                "case_review_id": case_review_id,
                "seq_id": seq_id,
                "reads": max(read_values),
                "detected": True,
                "selected": True,
                "selection_tags": merged_tags,
                "raw_available": True,
                "source_taxid": str(source.get("taxid") or "").strip(),
                "source_file": str(row.get("source_file") or ""),
                "source_row_index": row.get("row_index"),
                "duplicate_export_row_count": len(duplicates),
                "duplicate_export_source_files": sorted({
                    str(item.get("source_file") or "") for item in duplicates
                }),
                "duplicate_export_reads_conflict": len(read_values) > 1,
                "duplicate_export_taxid_conflict": len(taxid_values) > 1,
            })
        organisms.append({"category": category, "name": name, "measurements": measurements})

    return {
        "schema_version": "two_hospital_multi_assay_compatibility.v1",
        "patient_id": patient_id,
        "source_contract": {
            "kind": "legacy_flattened_retained_candidate_rows",
            "input_reads": "source_missing",
            "host_ratio": "source_missing",
            "technical_repeat_relation": "source_unknown",
            "molecule_identity": "inferred_from_source_path_when_unambiguous_else_unknown",
            "duplicate_export_handling": (
                "collapse case+seq_id+organism to one observation; preserve all source paths"
            ),
        },
        "cases": cases,
        "organisms": organisms,
    }


def copy_patient_view(
    source_dir: Path, destination: Path, patient_id: int, summary_suffix: str
) -> list[dict[str, Any]]:
    destination.mkdir(parents=True, exist_ok=False)
    copied: list[dict[str, Any]] = []
    for source in sorted(source_dir.glob("*.json")):
        target = destination / source.name
        shutil.copy2(source, target)
        copied.append({
            "source": str(source.resolve()),
            "destination": str(target.resolve()),
            "sha256": sha256_file(source),
        })
    summary_source = (
        source_dir
        / "summary_outputs"
        / f"NGS_patient_{patient_id}_{summary_suffix}.json"
    )
    if not summary_source.is_file():
        raise FileNotFoundError(summary_source)
    summary_target = destination / "summary_outputs" / summary_source.name
    summary_target.parent.mkdir()
    shutil.copy2(summary_source, summary_target)
    copied.append({
        "source": str(summary_source.resolve()),
        "destination": str(summary_target.resolve()),
        "sha256": sha256_file(summary_source),
    })
    (destination / "agent_outputs").mkdir()
    return copied


def supplied_multi_assay_document(
    source_root: Path, patient_id: int
) -> tuple[dict[str, Any], Path]:
    path = (
        source_root
        / f"NGS_patient_{patient_id}_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
    )
    if not path.is_file():
        raise FileNotFoundError(path)
    payload = read_json(path)
    if int(payload.get("patient_id") or -1) != patient_id:
        raise ValueError(
            f"Filename/payload patient mismatch for supplied mNGS: {path}"
        )
    if not payload.get("cases") or not payload.get("organisms"):
        raise ValueError(f"Supplied mNGS lacks cases or organisms: {path}")
    return payload, path


def run(
    source_root: Path,
    output_root: Path,
    summary_suffix: str,
    patient_ids: set[int] | None = None,
    multi_assay_source_root: Path | None = None,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_root}")
    patients_root = output_root / "patients"
    multi_assay_root = output_root / "multi_assay_source"
    phenotype_root = output_root / "phenotype_unavailable"
    patients_root.mkdir(parents=True, exist_ok=True)
    multi_assay_root.mkdir(parents=True, exist_ok=True)
    phenotype_root.mkdir(parents=True, exist_ok=True)

    included: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    molecule_counts: Counter[str] = Counter()
    for source_dir in sorted(
        (path for path in source_root.iterdir() if path.is_dir()),
        key=lambda path: int(PATIENT_RE.fullmatch(path.name).group(1))
        if PATIENT_RE.fullmatch(path.name) else 10**9,
    ):
        match = PATIENT_RE.fullmatch(source_dir.name)
        if not match:
            continue
        patient_id = int(match.group(1))
        if patient_ids is not None and patient_id not in patient_ids:
            continue
        supplied_path: Path | None = None
        raw_rows: list[dict[str, Any]] = []
        if multi_assay_source_root is not None:
            document, supplied_path = supplied_multi_assay_document(
                multi_assay_source_root, patient_id
            )
            source_kind = "supplied_full_per_test_multi_assay"
        else:
            raw_path = source_dir / f"NGS_patient_{patient_id}_raw_RK_NTC_rows.json"
            raw_rows = read_json(raw_path) if raw_path.is_file() else []
            if not isinstance(raw_rows, list) or not raw_rows:
                excluded.append({
                    "patient_id": patient_id,
                    "reason": "no_retained_mngs_candidate_rows",
                    "source": str(raw_path.resolve()),
                })
                continue
            document = build_multi_assay_document(source_dir, patient_id, raw_rows)
            source_kind = "adapted_legacy_flattened_rows"
        for case in document["cases"]:
            molecule_counts.update(test["nucleic_type"] for test in case["tests"])
        output_source = (
            multi_assay_root
            / f"NGS_patient_{patient_id}_all_mNGS_DNA_RNA_RK_NTC_by_species.json"
        )
        if supplied_path is not None:
            shutil.copy2(supplied_path, output_source)
        else:
            write_json(output_source, document)
        copied = copy_patient_view(
            source_dir,
            patients_root / f"NGS_patient_{patient_id}_json",
            patient_id,
            summary_suffix,
        )
        positive_observation_count = sum(
            1
            for organism in document["organisms"]
            for measurement in organism.get("measurements") or []
            if float(measurement.get("reads") or 0) > 0
        )
        adapted_observation_count = (
            positive_observation_count
            if supplied_path is not None
            else sum(
                len(organism["measurements"])
                for organism in document["organisms"]
            )
        )
        included.append({
            "patient_id": patient_id,
            "mngs_source_kind": source_kind,
            "raw_candidate_row_count": (
                len(raw_rows) if raw_rows else positive_observation_count
            ),
            "adapted_observation_count": adapted_observation_count,
            "collapsed_duplicate_export_row_count": (
                len(raw_rows) - adapted_observation_count
                if supplied_path is None else 0
            ),
            "case_count": len(document["cases"]),
            "test_count": sum(len(case["tests"]) for case in document["cases"]),
            "organism_count": len(document["organisms"]),
            "adapted_mngs_source": str(output_source.resolve()),
            "adapted_mngs_sha256": sha256_file(output_source),
            "supplied_mngs_source": (
                str(supplied_path.resolve()) if supplied_path is not None else None
            ),
            "supplied_mngs_sha256": (
                sha256_file(supplied_path) if supplied_path is not None else None
            ),
            "supplied_sha256_preserved": (
                sha256_file(supplied_path) == sha256_file(output_source)
                if supplied_path is not None else None
            ),
            "copied_patient_files": copied,
        })

    if patient_ids is not None:
        included_ids = {int(row["patient_id"]) for row in included}
        missing_ids = sorted(patient_ids - included_ids)
        if missing_ids:
            raise ValueError(
                f"Requested patients were not integrated: {missing_ids}"
            )

    summary = {
        "schema_version": "two_hospital_v2_6_compatibility_input.v1",
        "source_root": str(source_root.resolve()),
        "output_root": str(output_root.resolve()),
        "summary_suffix": summary_suffix,
        "supplied_multi_assay_source_root": (
            str(multi_assay_source_root.resolve())
            if multi_assay_source_root is not None else None
        ),
        "requested_patient_ids": sorted(patient_ids) if patient_ids is not None else None,
        "included_patient_count": len(included),
        "included_patient_ids": [row["patient_id"] for row in included],
        "excluded_patient_count": len(excluded),
        "excluded_patients": excluded,
        "test_count": sum(row["test_count"] for row in included),
        "candidate_row_count": sum(row["raw_candidate_row_count"] for row in included),
        "adapted_observation_count": sum(
            row["adapted_observation_count"] for row in included
        ),
        "collapsed_duplicate_export_row_count": sum(
            row["collapsed_duplicate_export_row_count"] for row in included
        ),
        "molecule_counts": dict(sorted(molecule_counts.items())),
        "patient_root": str(patients_root.resolve()),
        "multi_assay_source_root": str(multi_assay_root.resolve()),
        "phenotype_root": str(phenotype_root.resolve()),
        "patients": included,
        "constraints": [
            "This adapter does not invent RPM denominators or non-host reads.",
            "When a supplied full per-test multi-assay file is provided, it is copied byte-for-byte and its SHA-256 is preserved.",
            "Unknown molecule tests stay visible but cannot create DNA/RNA concordance or repeat credit.",
            "Missing collection time and specimen site stay explicit unknowns.",
            "Technical-repeat independence is unknown and cannot independently promote a candidate.",
            "Repeated legacy export rows for the same case/seq_id/organism are collapsed and retained as provenance, not counted as repeated tests.",
            "A requested patient without a supplied/retained mNGS source causes the run to fail rather than being treated as a negative test.",
        ],
    }
    write_json(output_root / "adapter_summary.json", summary)
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-suffix", default=DEFAULT_SUMMARY_SUFFIX)
    parser.add_argument(
        "--multi-assay-source-root",
        type=Path,
        help=(
            "Optional full per-test DNA/RNA mNGS source. When supplied, these "
            "files replace legacy flattened-row reconstruction."
        ),
    )
    parser.add_argument(
        "--patient-ids",
        help="Optional comma-separated patient IDs to adapt (for a frozen subset run).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    patient_ids = None
    if args.patient_ids:
        patient_ids = {
            int(value.strip())
            for value in args.patient_ids.split(",")
            if value.strip()
        }
        if not patient_ids:
            raise ValueError("--patient-ids did not contain any patient IDs")
    print(json.dumps(
        run(
            args.source_root,
            args.output_root,
            args.summary_suffix,
            patient_ids=patient_ids,
            multi_assay_source_root=args.multi_assay_source_root,
        ),
        ensure_ascii=False,
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

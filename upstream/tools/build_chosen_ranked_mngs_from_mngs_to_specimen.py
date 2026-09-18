from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.build_chosen_ranked_mngs_from_excel_mapping import expanded_name_keys, normalize_name  # noqa: E402
from tools.mngs_common import deduplicate_ranked_mngs_payload  # noqa: E402


DEFAULT_PATIENT_ROOT = Path("outputs") / "patient_info_rich_normalized"
DEFAULT_SOURCE_SUFFIX = "mngs_to_specimen_all_rk_ntc_opt"
DEFAULT_COMBINED_OUTPUT = Path("outputs") / "mngs_candidate_microbes_rich_all_rk_ntc_opt_chosen_ranked.json"
DEFAULT_REPORT = Path("outputs") / "mngs_candidate_microbes_rich_all_rk_ntc_opt_chosen_ranked_report.json"
DEFAULT_CHOSEN_SUFFIX = "mNGS_ranked_candidates_chosen"
PROTECTED_KEYWORDS = (
    "pneumocystis",
    "pjp",
    "nocardia",
    "legionella",
    "mycobacterium",
    "mycobacteroides",
    "mycolicibacterium",
    "mycolicibacter",
    "aspergillus",
    "mucor",
    "rhizopus",
    "rhizomucor",
    "lichtheimia",
    "cunninghamella",
    "saksenaea",
    "apophysomyces",
    "histoplasma",
    "cryptococcus",
    "cmv",
    "cytomegalovirus",
    "humanbetaherpesvirus5",
    "hsv1",
    "hsv2",
    "herpessimplex",
    "humanalphaherpesvirus1",
    "humanalphaherpesvirus2",
    "toxoplasma",
)


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def patient_id_from_dir(path: Path) -> str:
    return path.name.removesuffix("_json").removeprefix("NGS_patient_")


def is_protected_pathogen_name(value: Any) -> bool:
    keys = expanded_name_keys(value)
    keys.add(normalize_name(value))
    return any(any(keyword in key for keyword in PROTECTED_KEYWORDS) for key in keys)


def add_selected_name(
    *,
    name: str,
    specimen_code: str,
    by_specimen: dict[str, set[str]],
    all_keys: set[str],
    source_names: set[str],
) -> bool:
    keys = expanded_name_keys(name)
    if not keys:
        return False
    source_names.add(name)
    all_keys.update(keys)
    if specimen_code:
        by_specimen.setdefault(specimen_code, set()).update(keys)
    return True


def collect_organism_names_from_group(group: Any) -> list[str]:
    names: list[str] = []
    if not isinstance(group, dict):
        return names
    for organisms in group.values():
        if not isinstance(organisms, list):
            continue
        for organism in organisms:
            if not isinstance(organism, dict):
                continue
            name = str(organism.get("name") or organism.get("organism_name") or "").strip()
            if name:
                names.append(name)
    return names


def collect_selected_keys(specimen_outputs: Any) -> tuple[dict[str, set[str]], set[str], set[str]]:
    by_specimen: dict[str, set[str]] = {}
    all_keys: set[str] = set()
    source_names: set[str] = set()
    if not isinstance(specimen_outputs, list):
        return by_specimen, all_keys, source_names

    for record in specimen_outputs:
        if not isinstance(record, dict):
            continue
        specimen_code = str(record.get("specimen_code") or "").strip()
        pathogens = record.get("pathogens")
        for name in collect_organism_names_from_group(pathogens):
            add_selected_name(
                name=name,
                specimen_code=specimen_code,
                by_specimen=by_specimen,
                all_keys=all_keys,
                source_names=source_names,
            )

        likelihood = record.get("pathogens_by_likelihood")
        if not isinstance(likelihood, dict):
            continue
        for likelihood_group in ("high", "medium", "low_colonizer"):
            group = likelihood.get(likelihood_group)
            if not isinstance(group, dict):
                continue
            for name in collect_organism_names_from_group(group):
                add_selected_name(
                    name=name,
                    specimen_code=specimen_code,
                    by_specimen=by_specimen,
                    all_keys=all_keys,
                    source_names=source_names,
                )
    return by_specimen, all_keys, source_names


def filter_ranked_payload(
    ranked_payload: dict[str, Any],
    *,
    selected_by_specimen: dict[str, set[str]],
    selected_all_keys: set[str],
    selected_source_names: set[str],
    patient_id: str,
    source_summary_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    input_count = 0
    output_count = 0
    matched_keys: set[str] = set()
    matched_names: set[str] = set()
    filtered_records: list[dict[str, Any]] = []

    for record in ranked_payload.get("records") or []:
        if not isinstance(record, dict):
            continue
        specimen_code = str(record.get("specimen_code") or "").strip()
        wanted_keys = selected_by_specimen.get(specimen_code) or selected_all_keys
        filtered_candidates: list[dict[str, Any]] = []
        for candidate in record.get("candidates") or []:
            if not isinstance(candidate, dict):
                continue
            input_count += 1
            organism_name = candidate.get("organism_name") or candidate.get("name")
            candidate_keys = expanded_name_keys(organism_name)
            if candidate_keys & wanted_keys:
                filtered_candidates.append(candidate)
                matched_keys.update(candidate_keys & wanted_keys)
                matched_names.add(str(organism_name or "").strip())
        if filtered_candidates:
            new_record = dict(record)
            new_record["candidates"] = filtered_candidates
            new_record["status"] = "parsed"
            filtered_records.append(new_record)
            output_count += len(filtered_candidates)

    unmatched_source_names = [
        name
        for name in sorted(selected_source_names, key=lambda item: item.lower())
        if not (expanded_name_keys(name) & matched_keys)
    ]

    metadata = dict(ranked_payload.get("ranking_metadata") or {})
    metadata["chosen_filter"] = {
        "enabled": True,
        "source": "mngs_to_specimen_selected_pathogens",
        "source_summary_path": str(source_summary_path),
        "source_species_count": len(selected_source_names),
        "matched_ranked_species_count": len(matched_names),
        "input_candidate_rows": input_count,
        "output_candidate_rows": output_count,
        "unmatched_source_species": unmatched_source_names,
    }

    chosen_payload = deduplicate_ranked_mngs_payload(
        {
            "patient_id": str(ranked_payload.get("patient_id") or patient_id),
            "records": filtered_records,
            "ranking_metadata": metadata,
        }
    )
    dedup = (chosen_payload.get("ranking_metadata") or {}).get("deduplication") or {}
    report_entry = {
        "patient_id": patient_id,
        "source_species_count": len(selected_source_names),
        "matched_ranked_species_count": len(matched_names),
        "input_candidate_rows": input_count,
        "output_candidate_rows": output_count,
        "output_candidate_rows_after_dedup": dedup.get("output_candidate_rows", output_count),
        "matched_ranked_species": sorted(matched_names, key=lambda item: item.lower()),
        "unmatched_source_species": unmatched_source_names,
        "source_summary_path": str(source_summary_path),
    }
    return chosen_payload, report_entry


def build_chosen_ranked(
    patient_root: Path,
    source_suffix: str,
    combined_output: Path,
    report_output: Path,
    chosen_suffix: str = DEFAULT_CHOSEN_SUFFIX,
) -> dict[str, Any]:
    patients: dict[str, dict[str, Any]] = {}
    report: dict[str, Any] = {
        "patient_root": str(patient_root),
        "source_suffix": source_suffix,
        "chosen_suffix": chosen_suffix,
        "combined_output": str(combined_output),
        "written_patient_count": 0,
        "skipped_missing_source": 0,
        "skipped_missing_ranked": 0,
        "patients": [],
        "skipped": [],
    }

    for patient_dir in sorted(patient_root.glob("NGS_patient_*_json")):
        if not patient_dir.is_dir():
            continue
        patient_id = patient_id_from_dir(patient_dir)
        base = f"NGS_patient_{patient_id}"
        source_path = patient_dir / "summary_outputs" / f"{base}_{source_suffix}.json"
        ranked_path = patient_dir / f"{base}_mNGS_ranked_candidates.json"
        if not source_path.exists():
            report["skipped_missing_source"] += 1
            report["skipped"].append({"patient_id": patient_id, "reason": "missing_source", "path": str(source_path)})
            continue
        if not ranked_path.exists():
            report["skipped_missing_ranked"] += 1
            report["skipped"].append({"patient_id": patient_id, "reason": "missing_ranked", "path": str(ranked_path)})
            continue

        specimen_outputs = load_json(source_path)
        selected_by_specimen, selected_all_keys, selected_source_names = collect_selected_keys(specimen_outputs)
        ranked_payload = load_json(ranked_path)
        chosen_payload, report_entry = filter_ranked_payload(
            ranked_payload,
            selected_by_specimen=selected_by_specimen,
            selected_all_keys=selected_all_keys,
            selected_source_names=selected_source_names,
            patient_id=patient_id,
            source_summary_path=source_path,
        )
        chosen_path = patient_dir / f"{base}_{chosen_suffix}.json"
        write_json(chosen_path, chosen_payload)
        patients[patient_id] = chosen_payload
        report_entry["ranked_path"] = str(ranked_path)
        report_entry["chosen_path"] = str(chosen_path)
        report["patients"].append(report_entry)
        report["written_patient_count"] += 1

    combined = {
        "patients": patients,
        "ranking_metadata": {
            "source": "chosen_from_mngs_to_specimen_outputs",
            "patient_root": str(patient_root),
            "source_suffix": source_suffix,
            "chosen_suffix": chosen_suffix,
            "patient_count": len(patients),
        },
    }
    write_json(combined_output, combined)
    write_json(report_output, report)
    return report


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a combined ranked mNGS JSON by keeping only organisms selected in "
            "per-patient mNGS-to-specimen outputs."
        )
    )
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--source-suffix", default=DEFAULT_SOURCE_SUFFIX)
    parser.add_argument("--combined-output", type=Path, default=DEFAULT_COMBINED_OUTPUT)
    parser.add_argument("--report-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--chosen-suffix", default=DEFAULT_CHOSEN_SUFFIX)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report = build_chosen_ranked(
        args.patient_root,
        args.source_suffix,
        args.combined_output,
        args.report_output,
        args.chosen_suffix,
    )
    print(f"written_patient_count={report['written_patient_count']}")
    print(f"skipped_missing_source={report['skipped_missing_source']}")
    print(f"skipped_missing_ranked={report['skipped_missing_ranked']}")
    print(f"combined_output={args.combined_output}")
    print(f"report={args.report_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

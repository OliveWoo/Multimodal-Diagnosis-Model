"""Audit critical evidence preservation across normalized sources and agents.

This audit is answer-blind. It reports field-class loss and missing modules but
does not score organisms or alter any patient output.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Sequence

from tools import evidence_preservation


AUDIT_VERSION = "evidence_preservation_audit_v1"

MODULES: dict[str, dict[str, Any]] = {
    "mngs_ranked": {
        "source_patterns": ("*mNGS_ranked_candidates.json",),
        "agent_patterns": (),
        "downstream": "deterministic_scorer",
    },
    "culture": {
        "source_patterns": ("*culture.json",),
        "agent_patterns": ("*culture_agent.json",),
        "downstream": "small_agent",
    },
    "filmarray": {
        "source_patterns": ("*filmarray.json", "*FilmArray.json"),
        "agent_patterns": ("*filmarray_gmtest_agent.json",),
        "downstream": "small_agent",
    },
    "gm_test": {
        "source_patterns": ("*gm_test.json", "*gm test.json", "*GMtest.json"),
        "agent_patterns": ("*filmarray_gmtest_agent.json",),
        "downstream": "small_agent_plus_raw_fallback",
    },
    "molecular_microbiology": {
        "source_patterns": ("*molecular_microbiology.json", "*molecular microbiology.json"),
        "agent_patterns": ("*molecular_microbiology_agent.json",),
        "downstream": "small_agent",
    },
    "image": {
        "source_patterns": ("*image.json",),
        "agent_patterns": ("*image_agent.json",),
        "downstream": "small_agent_plus_raw_preservation",
    },
    "admission_diagnosis": {
        "source_patterns": ("*admission_diagnosis.json", "*admission diagnosis.json"),
        "agent_patterns": (),
        "downstream": "host_evidence_preservation",
    },
    "underlying": {
        "source_patterns": ("*underlying.json",),
        "agent_patterns": ("*cbc_other_lab_agent.json",),
        "downstream": "small_agent_plus_raw_preservation",
    },
    "cbc": {
        "source_patterns": ("*CBC.json", "*cbc.json"),
        "agent_patterns": ("*cbc_other_lab_agent.json",),
        "downstream": "small_agent_plus_raw_preservation",
    },
    "other_lab": {
        "source_patterns": ("*other_lab.json", "*other lab.json"),
        "agent_patterns": ("*cbc_other_lab_agent.json",),
        "downstream": "small_agent",
    },
}

def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def patient_id(patient_dir: Path) -> str:
    return patient_dir.name.removeprefix("NGS_patient_").removesuffix("_json")


def patient_sort_key(value: str) -> tuple[bool, int | str]:
    return (not value.isdigit(), int(value) if value.isdigit() else value)


def find_first(parent: Path, patterns: Iterable[str]) -> Path | None:
    for pattern in patterns:
        matches = sorted(path for path in parent.glob(pattern) if path.is_file())
        if matches:
            return matches[0]
    return None


def iter_keys(value: Any) -> Iterable[str]:
    if isinstance(value, dict):
        for key, item in value.items():
            yield str(key).lower()
            yield from iter_keys(item)
    elif isinstance(value, list):
        for item in value:
            yield from iter_keys(item)


def semantic_field_counts(value: Any) -> dict[str, int]:
    counts = Counter()
    for key in iter_keys(value):
        field_class = evidence_preservation.source_field_class(key)
        if field_class:
            counts[field_class] += 1
    return {
        field: counts.get(field, 0)
        for field in evidence_preservation.SOURCE_FIELD_CLASS_TERMS
    }


def row_count(value: Any) -> int:
    if isinstance(value, list):
        return len(value)
    if isinstance(value, dict):
        records = value.get("records")
        if isinstance(records, list):
            return len(records)
        return 1
    return 0


def load_optional(path: Path | None, errors: list[dict[str, str]]) -> Any:
    if path is None:
        return None
    try:
        return read_json(path)
    except Exception as exc:  # pragma: no cover - CLI diagnostics
        errors.append({"path": str(path), "error": str(exc)})
        return None


def module_record(patient_dir: Path, module: str, errors: list[dict[str, str]]) -> dict[str, Any]:
    spec = MODULES[module]
    source_path = find_first(patient_dir, spec["source_patterns"])
    agent_dir = patient_dir / "agent_outputs"
    agent_path = find_first(agent_dir, spec["agent_patterns"]) if spec["agent_patterns"] else None
    source = load_optional(source_path, errors)
    agent = load_optional(agent_path, errors)
    source_counts = semantic_field_counts(source)
    agent_counts = semantic_field_counts(agent)
    source_index = evidence_preservation.build_source_evidence_index({module: source})
    indexed_counts = source_index["sections"][0]["field_counts"]
    return {
        "module": module,
        "downstream": spec["downstream"],
        "source_path": str(source_path) if source_path else None,
        "source_available": source_path is not None,
        "source_row_count": row_count(source),
        "source_field_classes": source_counts,
        "agent_path": str(agent_path) if agent_path else None,
        "agent_available": agent_path is not None,
        "agent_field_classes": agent_counts,
        "preserved_index_field_classes": indexed_counts,
        "field_classes_absent_from_agent": [
            field
            for field, count in source_counts.items()
            if count > 0 and agent_path is not None and agent_counts.get(field, 0) == 0
        ],
        "field_classes_absent_from_preserved_index": [
            field
            for field, count in source_counts.items()
            if count > 0 and indexed_counts.get(field, 0) == 0
        ],
        "source_payload": source,
        "agent_payload": agent,
    }


def audit_patient(patient_dir: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    errors: list[dict[str, str]] = []
    modules = {module: module_record(patient_dir, module, errors) for module in MODULES}
    ranked = modules["mngs_ranked"]["source_payload"]
    sample_time = evidence_preservation.mngs_sample_time_from_ranked(ranked)

    image = modules["image"]
    enriched_image = evidence_preservation.enrich_image_agent_output(
        image["agent_payload"], image["source_payload"]
    )
    raw_mentions = evidence_preservation.extract_imaging_diagnostic_mentions(image["source_payload"])
    original_mentions = (
        image["agent_payload"].get("diagnostic_mentions")
        if isinstance(image["agent_payload"], dict)
        else []
    )

    host_agent = modules["cbc"]["agent_payload"] or modules["underlying"]["agent_payload"]
    enriched_host = evidence_preservation.enrich_host_agent_output(
        host_agent,
        raw_cbc=modules["cbc"]["source_payload"],
        raw_underlying=modules["underlying"]["source_payload"],
        raw_admission=modules["admission_diagnosis"]["source_payload"],
        sample_time=sample_time,
    )
    medication_profile = enriched_host["structured_host_evidence"]["medication_profile"]
    cbc_status = enriched_host["structured_host_evidence"]["cbc_data_status"]

    issues: list[dict[str, Any]] = []
    for module, record in modules.items():
        expects_agent = bool(MODULES[module]["agent_patterns"])
        if record["source_available"] and expects_agent and not record["agent_available"]:
            issues.append(
                {
                    "issue_type": "source_available_agent_missing",
                    "module": module,
                    "detail": record["source_path"],
                }
            )
        if record["field_classes_absent_from_agent"]:
            issues.append(
                {
                    "issue_type": "semantic_field_class_absent_from_agent",
                    "module": module,
                    "detail": record["field_classes_absent_from_agent"],
                }
            )
        if record["field_classes_absent_from_preserved_index"]:
            issues.append(
                {
                    "issue_type": "semantic_field_class_absent_from_preserved_index",
                    "module": module,
                    "detail": record["field_classes_absent_from_preserved_index"],
                }
            )
    if raw_mentions and not original_mentions:
        issues.append(
            {
                "issue_type": "imaging_diagnostic_mentions_missing_from_original_agent",
                "module": "image",
                "detail": [item["concept_name"] for item in raw_mentions],
            }
        )
    if medication_profile["risk_medication_count"] and not isinstance(host_agent, dict):
        issues.append(
            {
                "issue_type": "host_risk_medication_source_present_agent_missing",
                "module": "cbc_other_lab",
                "detail": [item["medication_name"] for item in medication_profile["medications"]],
            }
        )

    public_modules = []
    for module, record in modules.items():
        public_modules.append({key: value for key, value in record.items() if not key.endswith("_payload")})
    return (
        {
            "patient_id": patient_id(patient_dir),
            "patient_dir": str(patient_dir),
            "mngs_sample_time": sample_time,
            "modules": public_modules,
            "imaging_diagnostic_mentions": raw_mentions,
            "original_agent_imaging_mention_count": len(original_mentions or []),
            "preserved_imaging_mention_count": len(enriched_image.get("diagnostic_mentions") or []),
            "host_medication_profile": medication_profile,
            "cbc_data_status": cbc_status,
            "issues": issues,
        },
        errors,
    )


def audit_roots(roots: Sequence[Path]) -> dict[str, Any]:
    patients: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []
    seen: set[Path] = set()
    for root in roots:
        candidates = [root] if root.name.startswith("NGS_patient_") else sorted(root.glob("NGS_patient_*_json"))
        for candidate in candidates:
            if not candidate.is_dir():
                continue
            resolved = candidate.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            patient, patient_errors = audit_patient(candidate)
            patients.append(patient)
            errors.extend(patient_errors)
    patients.sort(key=lambda item: patient_sort_key(item["patient_id"]))

    source_available = Counter()
    agent_available = Counter()
    issue_counts = Counter()
    cbc_status_counts = Counter()
    imaging_patients = 0
    medication_patients = 0
    for patient in patients:
        for module in patient["modules"]:
            if module["source_available"]:
                source_available[module["module"]] += 1
            if module["agent_available"]:
                agent_available[module["module"]] += 1
        for issue in patient["issues"]:
            issue_counts[issue["issue_type"]] += 1
        cbc_status_counts[patient["cbc_data_status"]["status"]] += 1
        imaging_patients += int(bool(patient["imaging_diagnostic_mentions"]))
        medication_patients += int(patient["host_medication_profile"]["risk_medication_count"] > 0)

    return {
        "audit_version": AUDIT_VERSION,
        "answer_blind": True,
        "roots": [str(root) for root in roots],
        "patient_count": len(patients),
        "source_available_counts": dict(source_available),
        "agent_available_counts": dict(agent_available),
        "issue_counts": dict(issue_counts),
        "cbc_status_counts": dict(cbc_status_counts),
        "patients_with_imaging_diagnostic_mentions": imaging_patients,
        "patients_with_host_risk_medications": medication_patients,
        "json_read_errors": errors,
        "patients": patients,
    }


def markdown_summary(report: dict[str, Any]) -> str:
    preserved_gap_count = report["issue_counts"].get(
        "semantic_field_class_absent_from_preserved_index", 0
    )
    agent_gap_count = report["issue_counts"].get("semantic_field_class_absent_from_agent", 0)
    lines = [
        "# Evidence Preservation Audit",
        "",
        f"Audit version: `{report['audit_version']}`",
        "",
        "This report is answer-blind and does not change model tiers or outputs.",
        "",
        "## Summary",
        "",
        f"- Patients: {report['patient_count']}",
        f"- JSON read errors: {len(report['json_read_errors'])}",
        f"- Patients with imaging diagnostic mentions: {report['patients_with_imaging_diagnostic_mentions']}",
        f"- Patients with host-risk medication evidence: {report['patients_with_host_risk_medications']}",
        f"- Semantic field classes absent from original agents: {agent_gap_count}",
        f"- Semantic field classes absent from deterministic source index: {preserved_gap_count}",
        "",
        "## Module Availability",
        "",
        "| Module | Normalized source | Agent output |",
        "| --- | ---: | ---: |",
    ]
    for module in MODULES:
        lines.append(
            f"| {module} | {report['source_available_counts'].get(module, 0)} | "
            f"{report['agent_available_counts'].get(module, 0)} |"
        )
    lines.extend(["", "## Issue Counts", ""])
    for issue, count in sorted(report["issue_counts"].items()):
        lines.append(f"- `{issue}`: {count}")
    lines.extend(["", "## CBC Status", ""])
    for status, count in sorted(report["cbc_status_counts"].items()):
        lines.append(f"- `{status}`: {count}")
    lines.extend(["", "## Missing Agent Outputs", ""])
    missing_rows = []
    for patient in report["patients"]:
        modules = [
            issue["module"]
            for issue in patient["issues"]
            if issue["issue_type"] == "source_available_agent_missing"
        ]
        if modules:
            missing_rows.append((patient["patient_id"], modules))
    if missing_rows:
        for patient, modules in missing_rows:
            lines.append(f"- P{patient}: {', '.join(modules)}")
    else:
        lines.append("- None")
    lines.extend(["", "## Priority Acceptance Cases", ""])
    by_id = {patient["patient_id"]: patient for patient in report["patients"]}
    p14 = by_id.get("14", {})
    p9 = by_id.get("9", {})
    p14_pjp = [
        item
        for item in p14.get("imaging_diagnostic_mentions", [])
        if item.get("concept_name") == "Pneumocystis jirovecii"
    ]
    p9_meds = p9.get("host_medication_profile", {}).get("medications", [])
    lines.append(
        "- P14 imaging PJP differential preserved as context-only: "
        + ("PASS" if p14_pjp and p14_pjp[0].get("assertion") == "differential" else "NOT FOUND")
    )
    lines.append(
        "- P9 Prednisolone and Methotrexate preserved as structured medication evidence: "
        + (
            "PASS"
            if {item.get("medication_name") for item in p9_meds} >= {"prednisolone", "methotrexate"}
            else "NOT FOUND"
        )
    )
    lines.append(
        "- P9 CBC distinguished from missing differential/ALC: "
        + ("PASS" if p9.get("cbc_data_status", {}).get("status") == "cbc_available_differential_missing" else "NOT FOUND")
    )
    lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            "- Agent omissions remain visible for audit, but indexed critical source fields are preserved deterministically.",
            "- Imaging diagnostic mentions are context only and are never converted into culture, PCR, or other microbiologic confirmation.",
            "- Medication evidence adds structured host context but does not automatically change host tier or pathogen tier.",
            "- This audit begins at normalized patient JSON. Raw Excel/PPT-to-normalized extraction still requires a separate source-level contract audit.",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Audit evidence preservation without reading benchmark answers.")
    parser.add_argument("roots", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    report = audit_roots(args.roots)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8-sig")
    print(f"wrote {args.output}")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(markdown_summary(report), encoding="utf-8-sig")
        print(f"wrote {args.markdown}")
    print(
        f"patients={report['patient_count']} errors={len(report['json_read_errors'])} "
        f"imaging_mentions={report['patients_with_imaging_diagnostic_mentions']} "
        f"host_medications={report['patients_with_host_risk_medications']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Audit Raw JSON to normalized JSON traceability without benchmark answers."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence

from tools.source_provenance import (
    SPECIMEN_PATTERNS,
    build_source_provenance,
    extract_source_anchors,
)


AUDIT_VERSION = "source_to_normalized_contract_v2_text_integrity"

MOJIBAKE_MARKERS = frozenset(
    "嚗銝雿蝝憭摰璆閮霈撘芣迤銵瘨蝚璊隞靘詻頛蝯撠甇敶"
)
NORMALIZED_MODULE_SUFFIXES = (
    "_admission_diagnosis.json",
    "_CBC.json",
    "_culture.json",
    "_filmarray.json",
    "_gm_test.json",
    "_image.json",
    "_other_lab.json",
    "_underlying.json",
)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def flatten_text(value: Any) -> str:
    values: list[str] = []
    if isinstance(value, dict):
        for key, item in value.items():
            if key == "_source_provenance":
                continue
            values.append(str(key))
            values.append(flatten_text(item))
    elif isinstance(value, list):
        values.extend(flatten_text(item) for item in value)
    elif value is not None:
        values.append(str(value))
    return "\n".join(item for item in values if item)


def iter_text_nodes(value: Any, path: str = "$") -> Iterable[tuple[str, str]]:
    if isinstance(value, dict):
        for key, item in value.items():
            yield from iter_text_nodes(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from iter_text_nodes(item, f"{path}[{index}]")
    elif isinstance(value, str):
        yield path, value


def text_integrity_issues(text: str) -> list[str]:
    issues: list[str] = []
    if "\ufffd" in text:
        issues.append("unicode_replacement_character")
    if any(0xE000 <= ord(char) <= 0xF8FF for char in text):
        issues.append("private_use_character")
    if any(0xD800 <= ord(char) <= 0xDFFF for char in text):
        issues.append("unpaired_surrogate")
    if any(ord(char) < 32 and char not in "\t\n\r" for char in text):
        issues.append("unexpected_control_character")

    marker_positions = [index for index, char in enumerate(text) if char in MOJIBAKE_MARKERS]
    clustered_markers = any(
        marker_positions[index + 1] - marker_positions[index] <= 12
        for index in range(len(marker_positions) - 1)
    )
    if len(marker_positions) >= 3 or clustered_markers:
        issues.append("mojibake_marker_cluster")
    if text.count("?") >= 3 and any(ord(char) > 127 for char in text):
        issues.append("possible_lossy_question_marks")
    return issues


def text_integrity_summary(value: Any, *, max_examples: int = 5) -> dict[str, Any]:
    issue_counts: Counter[str] = Counter()
    examples: list[dict[str, Any]] = []
    text_node_count = 0
    suspicious_node_count = 0
    for path, text in iter_text_nodes(value):
        text_node_count += 1
        issues = text_integrity_issues(text)
        if not issues:
            continue
        suspicious_node_count += 1
        issue_counts.update(issues)
        if len(examples) < max_examples:
            suspicious_chars = sorted(
                {
                    char
                    for char in text
                    if char == "\ufffd"
                    or 0xE000 <= ord(char) <= 0xF8FF
                    or 0xD800 <= ord(char) <= 0xDFFF
                    or char in MOJIBAKE_MARKERS
                }
            )
            examples.append(
                {
                    "path": path,
                    "issues": issues,
                    "snippet": text[:240],
                    "suspicious_codepoints": [
                        f"U+{ord(char):04X}" for char in suspicious_chars[:20]
                    ],
                }
            )
    return {
        "text_node_count": text_node_count,
        "suspicious_node_count": suspicious_node_count,
        "has_suspicious_text": suspicious_node_count > 0,
        "issue_counts": dict(issue_counts),
        "examples": examples,
    }


def infer_text_issue_stage(
    raw_summary: dict[str, Any],
    normalized_summary: dict[str, Any] | None,
) -> str:
    raw_suspicious = bool(raw_summary.get("has_suspicious_text"))
    if normalized_summary is None:
        return "raw_only_normalized_unavailable" if raw_suspicious else "not_assessable"
    normalized_suspicious = bool(normalized_summary.get("has_suspicious_text"))
    if raw_suspicious and normalized_suspicious:
        return "present_in_raw_and_normalized"
    if raw_suspicious:
        return "raw_only_or_filtered_during_normalization"
    if normalized_suspicious:
        return "introduced_or_exposed_during_normalization"
    return "no_high_confidence_text_issue"


def normalize_loose(value: str) -> str:
    return re.sub(r"[^a-z0-9\u4e00-\u9fff]+", "", value.lower())


def counterpart_path(raw_path: Path, normalized_dir: Path) -> Path:
    name = raw_path.name
    if name.endswith("_Raw.json"):
        name = name[: -len("_Raw.json")] + ".json"
    return normalized_dir / name


def anchor_retained(anchor_type: str, anchor: dict[str, Any], normalized_text: str) -> bool:
    if anchor_type == "dates":
        normalized = anchor.get("normalized")
        variants = {
            str(normalized),
            str(normalized).replace("-", "/"),
            str(normalized).replace("-", ""),
        }
        return any(value and value in normalized_text for value in variants)
    if anchor_type == "specimens":
        canonical = anchor.get("canonical")
        for name, pattern in SPECIMEN_PATTERNS:
            if name == canonical:
                return bool(pattern.search(normalized_text))
        return False
    raw = str(anchor.get("raw") or "")
    return bool(raw) and normalize_loose(raw) in normalize_loose(normalized_text)


def audit_pair(raw_path: Path, normalized_path: Path) -> dict[str, Any]:
    errors: list[str] = []
    try:
        raw = read_json(raw_path)
    except Exception as exc:
        return {
            "raw_path": str(raw_path),
            "normalized_path": str(normalized_path),
            "status": "raw_json_error",
            "errors": [str(exc)],
        }
    normalized = None
    if normalized_path.exists():
        try:
            normalized = read_json(normalized_path)
        except Exception as exc:
            errors.append(f"normalized_json_error: {exc}")
    provenance = build_source_provenance(raw, raw_path)
    raw_text_integrity = text_integrity_summary(raw)
    normalized_text_integrity = (
        text_integrity_summary(normalized) if normalized is not None else None
    )
    anchors = extract_source_anchors(raw)
    normalized_text = flatten_text(normalized) if normalized is not None else ""
    anchor_summary: dict[str, dict[str, Any]] = {}
    for anchor_type, values in anchors.items():
        retained = [item for item in values if anchor_retained(anchor_type, item, normalized_text)]
        missing = [item for item in values if item not in retained]
        anchor_summary[anchor_type] = {
            "raw_count": len(values),
            "retained_count": len(retained),
            "missing_count": len(missing),
            "missing_examples": missing[:10],
        }
    embedded = normalized.get("_source_provenance") if isinstance(normalized, dict) else None
    embedded_hash_matches = bool(
        isinstance(embedded, dict)
        and embedded.get("raw_payload_sha256") == provenance["raw_payload_sha256"]
    )
    status = "paired"
    if not normalized_path.exists():
        status = "normalized_missing"
    elif errors:
        status = "normalized_json_error"
    return {
        "raw_path": str(raw_path),
        "normalized_path": str(normalized_path),
        "source_file": provenance["source_file"],
        "sheet": provenance["sheet"],
        "status": status,
        "raw_payload_sha256": provenance["raw_payload_sha256"],
        "raw_row_count": provenance["raw_row_count"],
        "raw_nonempty_cell_count": provenance["raw_nonempty_cell_count"],
        "embedded_provenance": isinstance(embedded, dict),
        "embedded_hash_matches": embedded_hash_matches,
        "raw_text_integrity": raw_text_integrity,
        "normalized_text_integrity": normalized_text_integrity,
        "text_issue_stage": infer_text_issue_stage(
            raw_text_integrity, normalized_text_integrity
        ),
        "anchor_summary": anchor_summary,
        "errors": errors,
    }


def patient_id_from_dir(path: Path) -> str:
    match = re.search(r"NGS_patient_(.+?)_json(?:_Raw)?$", path.name)
    return match.group(1) if match else path.name


def iter_raw_dirs(root: Path) -> Iterable[Path]:
    if root.name.endswith("_json_Raw"):
        yield root
        return
    yield from sorted(path for path in root.glob("NGS_patient_*_json_Raw") if path.is_dir())


def iter_normalized_dirs(root: Path) -> Iterable[Path]:
    if root.name.startswith("NGS_patient_") and root.name.endswith("_json"):
        yield root
        return
    yield from sorted(
        path
        for path in root.glob("NGS_patient_*_json")
        if path.is_dir() and not path.name.endswith("_json_Raw")
    )


def audit_normalized_only_dir(normalized_dir: Path) -> dict[str, Any]:
    files: list[dict[str, Any]] = []
    for path in sorted(normalized_dir.glob("*.json")):
        if not path.name.endswith(NORMALIZED_MODULE_SUFFIXES):
            continue
        try:
            payload = read_json(path)
            integrity = text_integrity_summary(payload)
            files.append(
                {
                    "path": str(path),
                    "status": "readable",
                    "text_integrity": integrity,
                }
            )
        except Exception as exc:
            files.append(
                {
                    "path": str(path),
                    "status": "json_error",
                    "error": str(exc),
                }
            )
    suspicious_files = [
        item
        for item in files
        if (item.get("text_integrity") or {}).get("has_suspicious_text")
    ]
    return {
        "module_file_count": len(files),
        "suspicious_module_file_count": len(suspicious_files),
        "text_issue_stage": (
            "normalized_only_unknown_origin"
            if suspicious_files
            else "no_high_confidence_text_issue"
        ),
        "files": files,
    }


def index_workbooks(roots: Sequence[Path]) -> dict[str, list[str]]:
    indexed: dict[str, list[str]] = defaultdict(list)
    for root in roots:
        if root.is_file() and root.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
            indexed[root.name.lower()].append(str(root))
            continue
        if not root.is_dir():
            continue
        for suffix in ("*.xlsx", "*.xlsm", "*.xls"):
            for path in root.rglob(suffix):
                indexed[path.name.lower()].append(str(path))
    return dict(indexed)


def audit_roots(
    roots: Sequence[Path],
    workbook_roots: Sequence[Path] = (),
) -> dict[str, Any]:
    pairs: list[dict[str, Any]] = []
    normalized_inventory: dict[str, dict[str, Any]] = {}
    paired_normalized_dirs: set[Path] = set()
    workbook_index = index_workbooks(workbook_roots)
    for root in roots:
        for raw_dir in iter_raw_dirs(root):
            normalized_dir = raw_dir.with_name(raw_dir.name[: -len("_Raw")])
            paired_normalized_dirs.add(normalized_dir.resolve())
            patient_id = patient_id_from_dir(raw_dir)
            raw_files = sorted(raw_dir.glob("*_Raw.json"))
            if not raw_files:
                pairs.append(
                    {
                        "patient_id": patient_id,
                        "status": "raw_directory_empty",
                        "raw_path": str(raw_dir),
                    }
                )
            for raw_path in raw_files:
                item = audit_pair(raw_path, counterpart_path(raw_path, normalized_dir))
                item["patient_id"] = patient_id
                source_file = str(item.get("source_file") or "").lower()
                item["source_workbook_matches"] = workbook_index.get(source_file, [])
                item["source_workbook_resolved"] = bool(item["source_workbook_matches"])
                pairs.append(item)
        for normalized_dir in iter_normalized_dirs(root):
            resolved = normalized_dir.resolve()
            if resolved in paired_normalized_dirs:
                continue
            patient_id = patient_id_from_dir(normalized_dir)
            raw_dir = normalized_dir.with_name(normalized_dir.name + "_Raw")
            normalized_inventory[str(resolved)] = {
                "patient_id": patient_id,
                "normalized_dir": str(normalized_dir),
                "raw_dir_expected": str(raw_dir),
                "raw_dir_available": raw_dir.exists(),
                "status": "raw_pair_available" if raw_dir.exists() else "normalized_only_provenance_gap",
            }
            if not raw_dir.exists():
                normalized_inventory[str(resolved)]["normalized_only_text_audit"] = (
                    audit_normalized_only_dir(normalized_dir)
                )

    status_counts = Counter(item.get("status") for item in pairs)
    text_stage_counts = Counter(item.get("text_issue_stage") for item in pairs)
    sheet_stats: dict[str, Counter[str]] = defaultdict(Counter)
    for item in pairs:
        sheet = str(item.get("sheet") or "unknown")
        stats = sheet_stats[sheet]
        stats["files"] += 1
        stats["paired"] += int(item.get("status") == "paired")
        stats["embedded_provenance"] += int(bool(item.get("embedded_provenance")))
        for anchor_type, summary in (item.get("anchor_summary") or {}).items():
            stats[f"{anchor_type}_raw"] += int(summary.get("raw_count", 0))
            stats[f"{anchor_type}_retained"] += int(summary.get("retained_count", 0))

    inventory = list(normalized_inventory.values())
    normalized_only_suspicious_patients = [
        item
        for item in inventory
        if (item.get("normalized_only_text_audit") or {}).get(
            "suspicious_module_file_count", 0
        )
    ]
    unique_source_files = sorted(
        {str(item.get("source_file")) for item in pairs if item.get("source_file")}
    )
    resolved_source_files = sorted(
        {
            str(item.get("source_file"))
            for item in pairs
            if item.get("source_file") and item.get("source_workbook_resolved")
        }
    )
    return {
        "audit_version": AUDIT_VERSION,
        "answer_blind": True,
        "tier_effect": "none_audit_only",
        "roots": [str(root) for root in roots],
        "raw_file_count": len(pairs),
        "paired_file_count": status_counts.get("paired", 0),
        "status_counts": dict(status_counts),
        "text_integrity_summary": {
            "paired_stage_counts": dict(text_stage_counts),
            "paired_raw_suspicious_file_count": sum(
                bool((item.get("raw_text_integrity") or {}).get("has_suspicious_text"))
                for item in pairs
            ),
            "paired_normalized_suspicious_file_count": sum(
                bool(
                    (item.get("normalized_text_integrity") or {}).get(
                        "has_suspicious_text"
                    )
                )
                for item in pairs
            ),
            "normalized_only_suspicious_patient_count": len(
                normalized_only_suspicious_patients
            ),
            "normalized_only_suspicious_patients": [
                item["patient_id"] for item in normalized_only_suspicious_patients
            ],
        },
        "unique_source_workbook_names": unique_source_files,
        "resolved_source_workbook_names": resolved_source_files,
        "unresolved_source_workbook_names": sorted(set(unique_source_files) - set(resolved_source_files)),
        "normalized_only_patient_count": sum(
            item["status"] == "normalized_only_provenance_gap" for item in inventory
        ),
        "normalized_only_patients": inventory,
        "sheet_stats": {sheet: dict(values) for sheet, values in sorted(sheet_stats.items())},
        "pairs": pairs,
    }


def ratio(numerator: int, denominator: int) -> str:
    return "n/a" if denominator == 0 else f"{numerator}/{denominator} = {numerator / denominator:.3f}"


def markdown_summary(report: dict[str, Any]) -> str:
    lines = [
        "# Source-to-Normalized Contract Audit",
        "",
        f"Audit version: {report['audit_version']}",
        "",
        "This audit is answer-blind and does not change any tier or picked output.",
        "",
        "## Traceability",
        "",
        f"- Raw module files inspected: {report['raw_file_count']}",
        f"- Raw/normalized files paired: {report['paired_file_count']}",
        f"- Source workbook names referenced by Raw JSON: {len(report['unique_source_workbook_names'])}",
        f"- Source workbooks still resolvable: {len(report['resolved_source_workbook_names'])}",
        f"- Normalized-only patients without Raw pair: {report['normalized_only_patient_count']}",
        "",
        "## Text Integrity",
        "",
        f"- Paired Raw files with high-confidence suspicious text: {report['text_integrity_summary']['paired_raw_suspicious_file_count']}",
        f"- Paired normalized files with high-confidence suspicious text: {report['text_integrity_summary']['paired_normalized_suspicious_file_count']}",
        f"- Normalized-only patients with suspicious module text: {report['text_integrity_summary']['normalized_only_suspicious_patient_count']}",
        "- Paired issue stages: "
        + ", ".join(
            f"{key}={value}"
            for key, value in sorted(
                report["text_integrity_summary"]["paired_stage_counts"].items()
            )
        ),
        "",
        "## Anchor Retention by Sheet",
        "",
        "| Sheet | Paired files | Dates | Quantities | Specimens | Embedded provenance |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for sheet, stats in report["sheet_stats"].items():
        lines.append(
            f"| {sheet} | {stats.get('paired', 0)}/{stats.get('files', 0)} | "
            f"{ratio(stats.get('dates_retained', 0), stats.get('dates_raw', 0))} | "
            f"{ratio(stats.get('quantities_retained', 0), stats.get('quantities_raw', 0))} | "
            f"{ratio(stats.get('specimens_retained', 0), stats.get('specimens_raw', 0))} | "
            f"{stats.get('embedded_provenance', 0)}/{stats.get('files', 0)} |"
        )
    lines.extend(["", "## Provenance Gaps", ""])
    if report["unresolved_source_workbook_names"]:
        lines.append(
            "- Raw JSON references source workbooks that are not present under the searched workbook roots: "
            + ", ".join(report["unresolved_source_workbook_names"])
        )
    gaps = [
        item for item in report["normalized_only_patients"]
        if item["status"] == "normalized_only_provenance_gap"
    ]
    if gaps:
        lines.append(
            "- The following normalized patient folders have no paired _json_Raw directory, "
            "so Excel/PPT-to-normalized completeness cannot be verified directly: "
            + ", ".join(f"P{item['patient_id']}" for item in gaps)
        )
    else:
        lines.append("- None")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Anchor retention is a deterministic traceability check, not a clinical accuracy score.",
            "- Missing anchors identify records that require source review; they do not prove that the normalized clinical interpretation is wrong.",
            "- Existing normalized files predate embedded provenance. Future core.llm_parser outputs include source path, SHA-256, row/cell counts, and a critical-anchor index.",
            "- Text-integrity flags are conservative triage signals, not proof of semantic corruption. Raw-stage flags cannot be attributed to Excel versus Excel-to-Raw conversion while the source workbooks are unavailable.",
            "- The original Raw JSON remains the authoritative source; LLM formatting must never be treated as lossless extraction by itself.",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument(
        "--workbook-root",
        action="append",
        type=Path,
        default=[],
        help="Optional file or directory roots used to resolve source_file workbook names.",
    )
    args = parser.parse_args(argv)
    report = audit_roots(args.roots, args.workbook_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8-sig",
    )
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(markdown_summary(report), encoding="utf-8-sig")
    print(
        f"raw_files={report['raw_file_count']} paired={report['paired_file_count']} "
        f"normalized_only_patients={report['normalized_only_patient_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

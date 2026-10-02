"""Audit phenotype source lineage, content hashes, and duplicate exports.

This audit deliberately separates four claims that were previously easy to
conflate:

1. a normalized workbook exists;
2. the downstream JSON was built from that exact workbook;
3. the workbook is linked to the intended pipeline patient/episode; and
4. the original clinical PDF is actually available for verification.

The script never treats an Excel/JSON hash match as proof of (3) or (4).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKBOOK_ROOT = ROOT / "phenotype for students" / "patient phenotype_send"
DEFAULT_CONTEXT_ROOT = ROOT / "outputs" / "phenotype_context_v1"
DEFAULT_DECISION_ROOT = ROOT / "outputs" / "phenotype_decision_evidence_v1"
DEFAULT_LINKAGE = (
    ROOT
    / "outputs"
    / "runs"
    / "2026-09-22_KH_phenotype_patient_linkage_audit_v1"
    / "phenotype_patient_linkage_audit.json"
)
DEFAULT_EPISODE_POLICY = ROOT / "rules" / "phenotype_episode_identity_v1.json"
DEFAULT_OUTPUT = (
    ROOT
    / "outputs"
    / "runs"
    / "2026-09-30_KH_phenotype_source_linkage_dedup_v2"
)

PATIENT_RE = re.compile(r"patient[ _](\d+)", re.IGNORECASE)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _patient_number(path_or_name: str) -> int | None:
    match = PATIENT_RE.search(path_or_name)
    return int(match.group(1)) if match else None


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: Iterable[str]) -> None:
    fieldnames = list(fields)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _cluster_members(cluster: Any) -> list[int]:
    if isinstance(cluster, dict):
        values = cluster.get("value", [])
    else:
        values = cluster
    return [int(value) for value in values]


def _cluster_map(raw_clusters: list[Any]) -> dict[int, str]:
    result: dict[int, str] = {}
    for index, cluster in enumerate(raw_clusters, start=1):
        values = _cluster_members(cluster)
        label = f"semantic_cluster_{index}:" + "-".join(str(v) for v in values)
        for value in values:
            result[int(value)] = label
    return result


def _confirmed_episode_map(path: Path) -> tuple[dict[int, dict[str, Any]], dict[tuple[int, ...], dict[str, Any]]]:
    if not path.exists():
        return {}, {}
    payload = _read_json(path)
    by_member: dict[int, dict[str, Any]] = {}
    by_cluster: dict[tuple[int, ...], dict[str, Any]] = {}
    for raw in payload.get("clusters", []):
        if not isinstance(raw, dict):
            continue
        members = tuple(sorted(int(value) for value in raw.get("members", [])))
        if not members:
            continue
        record = {**raw, "members": list(members)}
        by_cluster[members] = record
        for member in members:
            by_member[member] = record
    return by_member, by_cluster


def audit(
    workbook_root: Path,
    context_root: Path,
    decision_root: Path,
    linkage_path: Path,
    workspace_root: Path,
    episode_policy_path: Path = DEFAULT_EPISODE_POLICY,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    linkage = _read_json(linkage_path)
    linkage_by_patient = {
        int(row["pipeline_patient_number"]): row for row in linkage.get("patients", [])
    }
    kh_patients = set(linkage_by_patient)
    semantic_cluster_by_patient = _cluster_map(
        linkage.get("phenotype_duplicate_export_clusters", [])
    )
    confirmed_episode_by_patient, confirmed_episode_by_cluster = _confirmed_episode_map(
        episode_policy_path
    )

    workbooks: dict[int, Path] = {}
    for path in sorted(workbook_root.glob("patient *_phenotype.xlsx")):
        number = _patient_number(path.name)
        if number is not None:
            workbooks[number] = path

    contexts: dict[int, Path] = {}
    for path in sorted(context_root.glob("patient_*_phenotype_context.json")):
        number = _patient_number(path.name)
        if number is not None:
            contexts[number] = path

    decisions: dict[int, Path] = {}
    for path in sorted(decision_root.glob("patient_*_phenotype_decision_evidence_v1.json")):
        number = _patient_number(path.name)
        if number is not None:
            decisions[number] = path

    all_patients = sorted(set(workbooks) | set(contexts) | set(decisions) | kh_patients)
    actual_hash_by_patient = {number: _sha256(path) for number, path in workbooks.items()}
    hash_groups: dict[str, list[int]] = defaultdict(list)
    for number, digest in actual_hash_by_patient.items():
        hash_groups[digest].append(number)
    exact_duplicate_by_patient: dict[int, str] = {}
    for digest, members in hash_groups.items():
        if len(members) > 1:
            label = f"sha256:{digest[:12]}:" + "-".join(map(str, members))
            for member in members:
                exact_duplicate_by_patient[member] = label

    rows: list[dict[str, Any]] = []
    for number in all_patients:
        workbook = workbooks.get(number)
        context_path = contexts.get(number)
        decision_path = decisions.get(number)
        context = _read_json(context_path) if context_path else {}
        decision = _read_json(decision_path) if decision_path else {}
        context_source = context.get("source", {})
        decision_source = decision.get("source", {})
        actual_hash = actual_hash_by_patient.get(number, "")
        context_hash = str(context_source.get("workbook_sha256") or "")
        decision_hash = str(decision_source.get("workbook_sha256") or "")
        full_evidence_name = str(context_source.get("full_evidence_json") or "")
        full_evidence_path = context_root / full_evidence_name if full_evidence_name else None
        full_evidence_exists = bool(full_evidence_path and full_evidence_path.exists())
        full_evidence_hash = _sha256(full_evidence_path) if full_evidence_exists else ""
        expected_full_hash = str(decision_source.get("full_evidence_sha256") or "")
        link = linkage_by_patient.get(number, {})
        confirmed_episode = confirmed_episode_by_patient.get(number, {})
        identity_status = str(
            confirmed_episode.get("patient_identity")
            or link.get("patient_identity_assessment")
            or "not_in_kh_cohort"
        )
        episode_status = str(
            confirmed_episode.get("episode_relationship")
            or link.get("infection_episode_assessment")
            or "not_in_kh_cohort"
        )
        source_pdf_available = bool(decision_source.get("source_pdf_available", False))
        lineage_hash_match = bool(
            workbook
            and context_path
            and decision_path
            and actual_hash
            and actual_hash == context_hash == decision_hash
            and full_evidence_exists
            and full_evidence_hash == expected_full_hash
        )
        rows.append(
            {
                "patient_number": number,
                "in_kh_33_cohort": number in kh_patients,
                "workbook_file": workbook.name if workbook else "",
                "workbook_exists": bool(workbook),
                "context_exists": bool(context_path),
                "decision_evidence_exists": bool(decision_path),
                "full_evidence_exists": full_evidence_exists,
                "workbook_sha256": actual_hash,
                "context_workbook_hash_match": bool(actual_hash and actual_hash == context_hash),
                "decision_workbook_hash_match": bool(actual_hash and actual_hash == decision_hash),
                "full_evidence_hash_match": bool(
                    full_evidence_hash and full_evidence_hash == expected_full_hash
                ),
                "complete_normalized_lineage": lineage_hash_match,
                "patient_identity_status": identity_status,
                "infection_episode_status": episode_status,
                "semantic_duplicate_export_cluster": semantic_cluster_by_patient.get(number, ""),
                "exact_workbook_duplicate_cluster": exact_duplicate_by_patient.get(number, ""),
                "source_pdf_available": source_pdf_available,
                "production_link_verified": bool(link.get("production_link_verified", False))
                and source_pdf_available,
                "review_required": (
                    number in kh_patients
                    and (
                        not lineage_hash_match
                        or not source_pdf_available
                        or identity_status
                        in {
                            "phenotype_file_missing",
                            "local_same_patient_duplicate_export_cluster",
                        }
                        or episode_status
                        in {
                            "phenotype_file_missing",
                            "index_infection_episode_unresolved",
                            "multiple_exports_share_index_date_anchor",
                        }
                    )
                ),
            }
        )

    duplicate_rows: list[dict[str, Any]] = []
    for cluster in linkage.get("phenotype_duplicate_export_clusters", []):
        members = _cluster_members(cluster)
        confirmed = confirmed_episode_by_cluster.get(tuple(sorted(members)))
        duplicate_rows.append(
            {
                "duplicate_type": "semantic_export_cluster",
                "members": "|".join(map(str, members)),
                "member_count": len(members),
                "all_workbooks_present": all(member in workbooks for member in members),
                "action": (
                    "keep_separate_confirmed_distinct_pneumonia_episodes"
                    if confirmed
                    else "keep_separate_until_episode_or_identity_confirmation"
                ),
            }
        )
    for digest, members in sorted(hash_groups.items()):
        if len(members) > 1:
            duplicate_rows.append(
                {
                    "duplicate_type": "byte_identical_workbook",
                    "members": "|".join(map(str, members)),
                    "member_count": len(members),
                    "all_workbooks_present": True,
                    "action": f"safe_content_duplicate_sha256_{digest[:12]}",
                }
            )

    pdf_paths = sorted(
        path for path in workspace_root.rglob("*.pdf") if path.is_file()
    )
    pdf_rows = [
        {
            "pdf_path": str(path.relative_to(workspace_root)),
            "patient_number_from_name": _patient_number(path.name) or "",
            "sha256": _sha256(path),
        }
        for path in pdf_paths
    ]

    kh_rows = [row for row in rows if row["in_kh_33_cohort"]]
    assertions = {
        "all_49_workbooks_have_context": len(workbooks) == 49
        and all(number in contexts for number in workbooks),
        "all_49_workbooks_have_decision_evidence": len(workbooks) == 49
        and all(number in decisions for number in workbooks),
        "all_present_workbooks_have_hash_verified_lineage": all(
            row["complete_normalized_lineage"] for row in rows if row["workbook_exists"]
        ),
        "kh_cohort_size_is_33": len(kh_rows) == 33,
        "kh_missing_phenotype_is_explicit": sorted(
            row["patient_number"] for row in kh_rows if not row["workbook_exists"]
        )
        == [14, 16],
        "semantic_duplicate_clusters_are_preserved_not_merged": len(
            linkage.get("phenotype_duplicate_export_clusters", [])
        )
        == 4,
        "confirmed_recurrent_pneumonia_clusters_are_kept_separate": (
            len(confirmed_episode_by_cluster) == 4
            and all(
                str(row.get("action"))
                == "keep_separate_confirmed_distinct_pneumonia_episodes"
                for row in duplicate_rows
                if row.get("duplicate_type") == "semantic_export_cluster"
            )
        ),
        "no_unverified_pdf_is_claimed_available": all(
            not row["source_pdf_available"] for row in rows
        )
        if not pdf_paths
        else True,
    }
    summary = {
        "schema_version": "phenotype_source_linkage_dedup_audit.v2",
        "answer_blind": True,
        "scope": (
            "File lineage and duplicate/export audit only. It does not approve clinical "
            "phenotypes or patient/episode identity."
        ),
        "counts": {
            "workbooks": len(workbooks),
            "context_files": len(contexts),
            "decision_evidence_files": len(decisions),
            "kh_cohort_patients": len(kh_rows),
            "kh_complete_normalized_lineage": sum(
                bool(row["complete_normalized_lineage"]) for row in kh_rows
            ),
            "kh_missing_workbook": sum(not row["workbook_exists"] for row in kh_rows),
            "kh_source_pdf_available": sum(row["source_pdf_available"] for row in kh_rows),
            "kh_production_link_verified": sum(
                row["production_link_verified"] for row in kh_rows
            ),
            "kh_review_required": sum(row["review_required"] for row in kh_rows),
            "workspace_pdf_files": len(pdf_paths),
            "semantic_duplicate_clusters": len(
                linkage.get("phenotype_duplicate_export_clusters", [])
            ),
            "byte_identical_workbook_clusters": sum(
                len(members) > 1 for members in hash_groups.values()
            ),
            "identity_status": dict(
                Counter(row["patient_identity_status"] for row in kh_rows)
            ),
            "episode_status": dict(
                Counter(row["infection_episode_status"] for row in kh_rows)
            ),
        },
        "assertions": assertions,
        "all_local_lineage_assertions_pass": all(assertions.values()),
        "production_source_linkage_complete": all(
            row["production_link_verified"] for row in kh_rows
        ),
        "remaining_external_requirements": [
            "Original de-identified source PDFs or a hospital source-record identifier per phenotype export.",
            "Physician/source-holder confirmation for P14 and P16 missing phenotype exports.",
            "Episode disambiguation only for future export clusters not already covered by the confirmed recurrent-pneumonia policy.",
        ],
        "dedup_policy": (
            "The four confirmed clusters belong to the same patients but distinct pneumonia "
            "episodes and must remain separate. Future unresolved clusters must not be collapsed "
            "without explicit episode confirmation."
        ),
    }
    return summary, rows, duplicate_rows, pdf_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook-root", type=Path, default=DEFAULT_WORKBOOK_ROOT)
    parser.add_argument("--context-root", type=Path, default=DEFAULT_CONTEXT_ROOT)
    parser.add_argument("--decision-root", type=Path, default=DEFAULT_DECISION_ROOT)
    parser.add_argument("--linkage", type=Path, default=DEFAULT_LINKAGE)
    parser.add_argument(
        "--episode-policy", type=Path, default=DEFAULT_EPISODE_POLICY
    )
    parser.add_argument("--workspace-root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    summary, rows, duplicate_rows, pdf_rows = audit(
        args.workbook_root,
        args.context_root,
        args.decision_root,
        args.linkage,
        args.workspace_root,
        args.episode_policy,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "phenotype_source_linkage_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    _write_csv(
        args.output_dir / "phenotype_source_linkage.csv",
        rows,
        [
            "patient_number",
            "in_kh_33_cohort",
            "workbook_file",
            "workbook_exists",
            "context_exists",
            "decision_evidence_exists",
            "full_evidence_exists",
            "workbook_sha256",
            "context_workbook_hash_match",
            "decision_workbook_hash_match",
            "full_evidence_hash_match",
            "complete_normalized_lineage",
            "patient_identity_status",
            "infection_episode_status",
            "semantic_duplicate_export_cluster",
            "exact_workbook_duplicate_cluster",
            "source_pdf_available",
            "production_link_verified",
            "review_required",
        ],
    )
    _write_csv(
        args.output_dir / "duplicate_export_audit.csv",
        duplicate_rows,
        ["duplicate_type", "members", "member_count", "all_workbooks_present", "action"],
    )
    _write_csv(
        args.output_dir / "workspace_pdf_inventory.csv",
        pdf_rows,
        ["pdf_path", "patient_number_from_name", "sha256"],
    )
    _write_csv(
        args.output_dir / "physician_source_confirmation_remaining.csv",
        [row for row in rows if row["review_required"]],
        [
            "patient_number",
            "workbook_exists",
            "complete_normalized_lineage",
            "patient_identity_status",
            "infection_episode_status",
            "semantic_duplicate_export_cluster",
            "source_pdf_available",
            "production_link_verified",
        ],
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

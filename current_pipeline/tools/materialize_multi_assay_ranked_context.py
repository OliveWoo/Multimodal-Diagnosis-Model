"""Materialize per-test ranked context from compact multi-assay packets.

The generated file supplies the index collection time and traceable per-test
context expected by small-agent evidence preservation.  It is not the source
of candidate decisions; the test-aware scorer continues to read compact
multi-assay packets directly.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Sequence

from tools.build_index_event_clinical_timeline_shadow import sha256_file
from tools.run_legacy_scorer_on_multi_assay_shadow import build_legacy_ranked_payload


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_COMPACT = ROOT / "outputs/runs/2026-09-30_52patients_compact_entry_no_zero_v1"
DEFAULT_PATIENT_ROOT = ROOT / "outputs/extracted_workbooks_pipeline_input_v3_51patients/patients"
DEFAULT_REPORT = ROOT / "outputs/runs/2026-09-30_52patients_ranked_context_materialization_v1/summary.json"
PATIENT_RE = re.compile(r"^NGS_patient_(\d+)_compact_entry_shadow\.json$")


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected object JSON: {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def patient_id(path: Path) -> int:
    match = PATIENT_RE.fullmatch(path.name)
    if match is None:
        raise ValueError(f"Unexpected compact packet name: {path.name}")
    return int(match.group(1))


def run(
    compact_root: Path,
    patient_root: Path,
    report_path: Path,
    *,
    replace_existing: bool,
) -> dict[str, Any]:
    packets = sorted(
        (compact_root / "patient_packets").glob("NGS_patient_*_compact_entry_shadow.json"),
        key=patient_id,
    )
    if not packets:
        raise FileNotFoundError(f"No compact packets found: {compact_root}")

    written: list[int] = []
    skipped_existing: list[int] = []
    missing_patient_dirs: list[int] = []
    record_count = 0
    candidate_count = 0
    for packet_path in packets:
        pid = patient_id(packet_path)
        patient_dir = patient_root / f"NGS_patient_{pid}_json"
        if not patient_dir.is_dir():
            missing_patient_dirs.append(pid)
            continue
        destination = patient_dir / f"NGS_patient_{pid}_mNGS_ranked_candidates.json"
        if destination.exists() and not replace_existing:
            skipped_existing.append(pid)
            continue
        packet = read_json(packet_path)
        ranked, _ = build_legacy_ranked_payload(packet)
        ranked["ranking_metadata"].update(
            {
                "schema_version": "multi_assay_ranked_context.v1",
                "decision_input_role": (
                    "collection-time and per-test context only; candidate decisions use the compact test-aware packet"
                ),
                "source_compact_packet": str(packet_path.resolve()),
                "source_compact_packet_sha256": sha256_file(packet_path),
                "zero_read_candidates_allowed": False,
            }
        )
        if any(
            float(candidate.get("reads") or 0) <= 0
            for record in ranked.get("records") or []
            for candidate in record.get("candidates") or []
        ):
            raise ValueError(f"Zero-read candidate reached ranked context for P{pid}")
        write_json(destination, ranked)
        written.append(pid)
        record_count += len(ranked.get("records") or [])
        candidate_count += sum(
            len(record.get("candidates") or []) for record in ranked.get("records") or []
        )

    summary = {
        "schema_version": "multi_assay_ranked_context_materialization_summary.v1",
        "compact_root": str(compact_root.resolve()),
        "patient_root": str(patient_root.resolve()),
        "packet_count": len(packets),
        "written_patient_count": len(written),
        "written_patient_ids": written,
        "skipped_existing_patient_ids": skipped_existing,
        "missing_patient_directory_ids": missing_patient_dirs,
        "ranked_record_count": record_count,
        "ranked_candidate_observation_count": candidate_count,
        "replace_existing": replace_existing,
        "clinical_decision_effect": "none_context_adapter_only",
    }
    write_json(report_path, summary)
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compact-root", type=Path, default=DEFAULT_COMPACT)
    parser.add_argument("--patient-root", type=Path, default=DEFAULT_PATIENT_ROOT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--replace-existing", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    print(json.dumps(run(
        args.compact_root,
        args.patient_root,
        args.report,
        replace_existing=args.replace_existing,
    ), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

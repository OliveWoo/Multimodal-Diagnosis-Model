#!/usr/bin/env python3
"""Create final-decision-aligned reasoning chains for rationale generation.

The source reasoning chains are never changed.  Each output chain keeps the
case evidence from the source and applies an explicit ``ober.final_decision``
file to the disposition/rank fields that the rationale generator freezes.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compact_unique(values: list[Any]) -> list[Any]:
    output: list[Any] = []
    for value in values:
        if value not in output:
            output.append(value)
    return output


def parse_patient_ids(value: str | None) -> set[str] | None:
    if not value:
        return None
    output = {item.strip().removeprefix("P") for item in value.split(",") if item.strip()}
    if not output or any(not item.isdigit() for item in output):
        raise ValueError("--patients must be comma-separated numeric patient IDs")
    return output


def chain_root_for_selected_root(selected_root: Path) -> Path:
    root = selected_root.parent / "run" / "rationale_views" / "traceable_reasoning_chains"
    if not root.is_dir():
        raise FileNotFoundError(f"reasoning-chain root not found: {root}")
    return root


def apply_final_decision(
    chain: dict[str, Any], decision: dict[str, Any], decision_file: Path
) -> dict[str, Any]:
    output = copy.deepcopy(chain)
    patient_id = str((output.get("病例識別") or {}).get("病人編號") or "")
    if patient_id != str(decision.get("patient_id") or ""):
        raise ValueError(f"patient mismatch for {decision_file}")
    if decision.get("stage") != "workflow_final":
        raise ValueError(f"not a workflow-final decision: {decision_file}")

    candidates = output.get("步驟2至5_逐候選病原判讀")
    if not isinstance(candidates, list):
        raise ValueError(f"missing candidate list for P{patient_id}")
    chain_by_name = {item.get("病原"): item for item in candidates if isinstance(item, dict)}
    decisions = decision.get("candidates") or []
    decision_by_name = {item.get("organism"): item for item in decisions}
    missing_in_decision = sorted(set(chain_by_name) - set(decision_by_name))
    if missing_in_decision:
        raise ValueError(
            f"candidate mismatch P{patient_id}; missing_in_decision={missing_in_decision}"
        )
    # The R5 finalizer can append a review-tier candidate that the older
    # presentation-oriented chain did not contain.  Preserve it explicitly
    # instead of dropping it or pretending that case evidence was available.
    for organism in sorted(set(decision_by_name) - set(chain_by_name)):
        item = decision_by_name[organism]
        signals = item.get("R5_signals") or {}
        candidate = {
            "病原": organism,
            "最終處置": "未選入",
            "最終排序": None,
            "一句話結論": "由 OBER review 候選池補入完整候選清單；病例鏈未提供可定位證據。",
            "為什麼會被考慮": [
                f"由 {signals.get('review_tier') or 'OBER review'} 送入 A1／病例適配判讀。"
            ],
            "為什麼選入或排除": list(item.get("reasons") or []),
            "關鍵判斷資訊": {"OBER R5 訊號": signals},
            "重要限制": [
                "原始 reasoning chain 未列此候選，因此目前沒有可定位的病例證據；不代表醫院端沒有原始資料。"
            ],
            "病例內支持或反對證據": [],
            "外部文獻病例或網站": {
                "狀態": "已有保存的 A1 判讀；由最終交付流程另行串接",
                "說明": "此 reasoning chain 僅供理由產生器凍結選擇與病例證據引用。",
                "證據": [],
            },
        }
        candidates.append(candidate)
        chain_by_name[organism] = candidate

    selected = sorted(decision.get("selected") or [], key=lambda item: item["rank"])
    selected_by_name = {item["organism"]: item for item in selected}
    expected_ranks = list(range(1, len(selected) + 1))
    if [item["rank"] for item in selected] != expected_ranks:
        raise ValueError(f"non-contiguous final ranks for P{patient_id}")

    for organism, candidate in chain_by_name.items():
        item = decision_by_name[organism]
        is_selected = organism in selected_by_name
        candidate["最終處置"] = "選入" if is_selected else "未選入"
        candidate["最終排序"] = selected_by_name.get(organism, {}).get("rank")
        decision_reasons = [str(value) for value in (item.get("reasons") or []) if str(value).strip()]
        candidate["為什麼選入或排除"] = compact_unique(
            list(candidate.get("為什麼選入或排除") or []) + decision_reasons
        )
        if item.get("selection_origin") == "OBER_R5_rescue":
            candidate["一句話結論"] = "選入（OBER R5 補漏）：" + " ".join(decision_reasons)
            candidate["為什麼會被考慮"] = compact_unique(
                list(candidate.get("為什麼會被考慮") or [])
                + ["此病原未被上游 picked 選入，進入 OBER 候選補漏判讀。"]
            )
            candidate["重要限制"] = compact_unique(
                list(candidate.get("重要限制") or [])
                + [
                    "此項是 R5 規則的計算補入，尚未經臨床專家逐例確認；文獻支持不等同此病人已確診。",
                    "補入順位接在上游 picked 之後，並非臨床信心排序。",
                ]
            )
        elif is_selected:
            candidate["一句話結論"] = candidate.get("一句話結論") or (
                "選入：" + " ".join(decision_reasons)
            )
        else:
            candidate["一句話結論"] = "未選入（OBER R5）：" + " ".join(decision_reasons)
        features = candidate.setdefault("關鍵判斷資訊", {})
        features["OBER R5 最終處置"] = candidate["最終處置"]
        features["OBER R5 規則路徑"] = item.get("R5_signals", {}).get("rule_path") if item.get("R5_signals") else None
        if item.get("R5_signals"):
            features["OBER R5 訊號"] = item["R5_signals"]

    overview = output.setdefault("一眼看懂", {})
    overview["候選病原數"] = len(candidates)
    overview["入選病原數"] = len(selected)
    overview["最終入選病原"] = [item["organism"] for item in selected]
    overview["尚待補充"] = "臨床專家逐例複核"
    overview["重要提醒"] = (
        "此處是 OBER R5 計算最終選擇；R5 在本次執行未讀取 gold answer，"
        "但規則曾在開發資料上查看答案，且尚未完成外部測試與臨床專家複核。"
    )
    output["步驟6_最終鑑別診斷排序"] = [
        {
            "排序": item["rank"],
            "病原": item["organism"],
            "選入理由": " ".join(decision_by_name[item["organism"]].get("reasons") or []),
            "選擇來源": item.get("selection_origin"),
            "規則路徑": item.get("rule_path"),
        }
        for item in selected
    ]
    output["OBER_R5_最終決策"] = {
        "決策狀態": decision.get("status"),
        "決策來源": decision.get("source"),
        "規則": decision.get("rule"),
        "驗證聲明": decision.get("validation"),
        "來源決策檔": decision_file.name,
        "來源決策檔SHA256": sha256_file(decision_file),
        "說明": "此副本只為最終理由產生器對齊；原始 reasoning chain 未被改寫。",
    }
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--final-decisions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--patients", help="optional comma-separated numeric patient IDs")
    parser.add_argument(
        "--cohorts",
        nargs="+",
        help="optional cohort names; useful when patient IDs overlap across cohorts",
    )
    args = parser.parse_args()

    config_file = args.config.resolve()
    config = read_json(config_file)
    config_root = config_file.parent
    decisions_root = args.final_decisions.resolve()
    output_root = args.output.resolve()
    patient_filter = parse_patient_ids(args.patients)
    cohort_filter = set(args.cohorts or [])
    unknown_cohorts = sorted(cohort_filter - set(config.get("cohorts", {})))
    if unknown_cohorts:
        raise ValueError(f"unknown cohorts: {unknown_cohorts}")
    if output_root.exists() and any(output_root.iterdir()):
        raise SystemExit(f"output must be new or empty: {output_root}")

    manifest_rows: list[dict[str, Any]] = []
    for cohort, cohort_config in config["cohorts"].items():
        if cohort_filter and cohort not in cohort_filter:
            continue
        selected_root = (config_root / cohort_config["selected_root"]).resolve()
        chain_root = chain_root_for_selected_root(selected_root)
        decision_dir = decisions_root / cohort
        for decision_file in sorted(decision_dir.glob("P*_final_decision.json")):
            decision = read_json(decision_file)
            patient_id = str(decision["patient_id"])
            if patient_filter is not None and patient_id not in patient_filter:
                continue
            source_file = chain_root / f"NGS_patient_{patient_id}_explainable_reasoning_zh.json"
            if not source_file.is_file():
                raise FileNotFoundError(source_file)
            expected_hash = (decision.get("source_snapshot") or {}).get("reasoning_chain_sha256")
            actual_hash = sha256_file(source_file)
            if expected_hash and expected_hash != actual_hash:
                raise ValueError(f"source reasoning-chain hash changed: {cohort}/P{patient_id}")
            output_chain = apply_final_decision(read_json(source_file), decision, decision_file)
            destination = output_root / cohort / f"NGS_patient_{patient_id}_explainable_reasoning_zh.json"
            write_json(destination, output_chain)
            manifest_rows.append(
                {
                    "cohort": cohort,
                    "patient_id": patient_id,
                    "source_file": str(source_file),
                    "source_sha256": actual_hash,
                    "decision_file": str(decision_file),
                    "decision_sha256": sha256_file(decision_file),
                    "output_file": str(destination),
                    "output_sha256": sha256_file(destination),
                    "selected": len(decision.get("selected") or []),
                }
            )

    if patient_filter is not None:
        found = {row["patient_id"] for row in manifest_rows}
        if found != patient_filter:
            raise ValueError(f"patient filter mismatch; requested={sorted(patient_filter)}, found={sorted(found)}")
    manifest = {
        "schema_version": "ober.final_rationale_input_manifest.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "config": str(config_file),
        "final_decisions_root": str(decisions_root),
        "output_root": str(output_root),
        "cohort_filter": sorted(cohort_filter) if cohort_filter else None,
        "summary": {
            "patients": len(manifest_rows),
            "selected": sum(row["selected"] for row in manifest_rows),
        },
        "patients": manifest_rows,
    }
    write_json(output_root / "manifest.json", manifest)
    print(json.dumps(manifest["summary"], ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

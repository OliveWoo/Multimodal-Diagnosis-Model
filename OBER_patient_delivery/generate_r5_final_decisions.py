from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import unicodedata
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping


SCHEMA_VERSION = "ober.final_decision.v1"
MANIFEST_VERSION = "ober.final_decision_manifest.v1"
RULE_ID = "R5_A1_UNION_AND_B_LOOSE"
RULE_FORMULA = (
    "picked OR C3 OR ((A1_strong >= 2 OR A1_match >= 7) "
    "AND site_aligned AND (B_absolute_axis OR B_relative_axis))"
)
VALID_REVIEW_TIERS = {"review_high_priority", "review_context_needed"}


class FinalizerError(ValueError):
    pass


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FinalizerError(f"could not read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise FinalizerError(f"expected JSON object: {path}")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def normal(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value or "")).strip().lower()
    return re.sub(r"\s+", " ", text)


def candidate_name(row: Mapping[str, Any]) -> str:
    for key in (
        "canonical_organism_name",
        "organism_name",
        "organism",
        "pathogen_name",
        "pathogen",
        "name",
    ):
        if row.get(key) not in (None, ""):
            return str(row[key])
    return ""


def canonical(value: Any, aliases: Mapping[str, str]) -> str:
    current = normal(value)
    if not current:
        raise FinalizerError("pathogen name must not be empty")
    seen: set[str] = set()
    while current in aliases:
        if current in seen:
            raise FinalizerError(f"alias cycle for {value!r}")
        seen.add(current)
        current = aliases[current]
    return current


def resolved_aliases(raw: Mapping[str, str]) -> dict[str, str]:
    aliases = {normal(key): normal(value) for key, value in raw.items()}
    for start in list(aliases):
        aliases[start] = canonical(start, aliases)
    return aliases


def strict_bool_or_none(value: Any, where: str) -> bool | None:
    if value is True or value is False or value is None:
        return value
    raise FinalizerError(f"{where} must be boolean or null")


def evaluate_r5(signals: Mapping[str, Any]) -> dict[str, Any]:
    c3 = signals.get("C_grade") == "C3"
    a1 = signals.get("A1_strong", 0) >= 2 or signals.get("A1_match", 0) >= 7
    site = signals.get("B_site_aligned") is True
    b_axis = (
        signals.get("B_absolute_axis") is True
        or signals.get("B_relative_axis") is True
    )
    selected = c3 or (a1 and site and b_axis)
    return {
        "selected": selected,
        "rule_path": "C3_direct" if c3 else "A1_site_B" if selected else None,
        "gates": {
            "C3_direct": c3,
            "A1_threshold": a1,
            "site_aligned": site,
            "B_axis_any": b_axis,
        },
    }


def signal_reason(signal: Mapping[str, Any], evaluation: Mapping[str, Any]) -> list[str]:
    gates = evaluation["gates"]
    values = (
        f"A1 strong={signal['A1_strong']}、partial={signal['A1_partial']}、"
        f"match={signal['A1_match']}；site_aligned={signal['B_site_aligned']}；"
        f"B absolute={signal['B_absolute_axis']}、relative={signal['B_relative_axis']}；"
        f"C={signal['C_grade']}。"
    )
    if evaluation["selected"]:
        if evaluation["rule_path"] == "C3_direct":
            return [f"R5 由 C3 直接補入。{values}"]
        return [
            "R5 補入：A1 文獻門檻、病例部位對齊及至少一個 B 強度軸均通過。"
            + values
        ]
    failed = [
        label
        for key, label in (
            ("A1_threshold", "A1 門檻"),
            ("site_aligned", "site align"),
            ("B_axis_any", "B absolute／relative 軸"),
        )
        if not gates[key]
    ]
    return [
        "未由 R5 補入：不是 C3，且 A1＋site＋B 路徑未全部通過"
        f"（未通過：{'、'.join(failed) or '無'}）。{values}",
        "這只表示未通過本次補漏規則，不等同臨床排除感染。",
    ]


def compose_decision(
    *,
    cohort: str,
    patient_id: str,
    candidate_names: Iterable[str],
    picked: list[Mapping[str, Any]],
    signals: list[Mapping[str, Any]],
    source_snapshot: Mapping[str, Any],
) -> dict[str, Any]:
    names = list(candidate_names)
    if any(not isinstance(name, str) or not name.strip() for name in names):
        raise FinalizerError(f"{cohort}/P{patient_id}: invalid candidate name")
    keys = [normal(name) for name in names]
    if len(keys) != len(set(keys)):
        raise FinalizerError(f"{cohort}/P{patient_id}: duplicate candidate")

    picked_sorted = sorted(picked, key=lambda item: item["rank"])
    ranks = [item["rank"] for item in picked_sorted]
    if ranks != list(range(1, len(ranks) + 1)):
        raise FinalizerError(f"{cohort}/P{patient_id}: picked ranks are not contiguous")
    picked_by_key = {normal(item["canonical_name"]): item for item in picked_sorted}
    if len(picked_by_key) != len(picked_sorted):
        raise FinalizerError(f"{cohort}/P{patient_id}: duplicate picked organism")
    missing_picked = sorted(set(picked_by_key) - set(keys))
    if missing_picked:
        raise FinalizerError(
            f"{cohort}/P{patient_id}: picked organism absent from reasoning chain: {missing_picked}"
        )

    signal_by_key: dict[str, Mapping[str, Any]] = {}
    for signal in signals:
        key = normal(signal["candidate_organism"])
        if key in signal_by_key:
            raise FinalizerError(f"{cohort}/P{patient_id}: duplicate casefit signal {key}")
        signal_by_key[key] = signal
        if key not in keys:
            names.append(signal["candidate_organism"])
            keys.append(key)
    overlap = set(picked_by_key) & set(signal_by_key)
    if overlap:
        raise FinalizerError(
            f"{cohort}/P{patient_id}: casefit review candidate overlaps picked: {sorted(overlap)}"
        )

    rescued = sorted(
        [
            signal
            for signal in signal_by_key.values()
            if evaluate_r5(signal)["selected"]
        ],
        key=lambda signal: normal(signal["candidate_organism"]),
    )
    final_selected = [
        {
            "organism": item["canonical_name"],
            "rank": item["rank"],
            "selection_origin": "upstream_picked",
            "rule_path": "picked_immutable",
        }
        for item in picked_sorted
    ]
    for offset, signal in enumerate(rescued, start=len(final_selected) + 1):
        final_selected.append(
            {
                "organism": signal["candidate_organism"],
                "rank": offset,
                "selection_origin": "OBER_R5_rescue",
                "rule_path": evaluate_r5(signal)["rule_path"],
            }
        )

    selected_by_key = {normal(item["organism"]): item for item in final_selected}
    candidate_decisions: list[dict[str, Any]] = []
    for name in names:
        key = normal(name)
        signal = signal_by_key.get(key)
        if key in picked_by_key:
            reasons = [
                "保留上游 picked；R5 是只補漏、不剪枝的規則，不移除既有選擇。"
            ]
            origin = "upstream_picked"
        elif signal is not None:
            reasons = signal_reason(signal, evaluate_r5(signal))
            origin = "OBER_R5_rescue" if key in selected_by_key else "not_selected"
        else:
            reasons = [
                "未在上游 picked，且不在本批 review_high_priority／review_context_needed "
                "的 A1 casefit 輸入中；本次 R5 沒有可用的 A1＋B＋C3 訊號，因此不補入。",
                "這只表示不在本次補漏判讀範圍，不等同臨床排除感染。",
            ]
            origin = "not_selected"
        candidate_decisions.append(
            {
                "organism": name,
                "decision": "selected" if key in selected_by_key else "not_selected",
                "rank": selected_by_key.get(key, {}).get("rank"),
                "selection_origin": origin,
                "reasons": reasons,
                "R5_signals": signal,
            }
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "patient_id": patient_id,
        "cohort": cohort,
        "stage": "workflow_final",
        "status": "computational_final_pending_rationale_and_clinical_review",
        "source": f"OBER {RULE_ID} frozen replay; gold answers not read",
        "rule": {
            "id": RULE_ID,
            "formula": RULE_FORMULA,
            "picked_policy": "immutable_no_pruning",
            "rescued_rank_policy": (
                "append after picked; alphabetical only for simultaneous rescues, "
                "not a clinical-confidence ranking"
            ),
            "D_policy": "reported for audit but not used as a hidden eligibility gate",
        },
        "validation": {
            "gold_answer_used_in_this_run": False,
            "development_rule_answers_previously_inspected": True,
            "untouched_external_test_validated": False,
            "clinical_expert_reviewed": False,
        },
        "selected": final_selected,
        "candidates": candidate_decisions,
        "counts": {
            "upstream_picked": len(picked_sorted),
            "R5_rescued": len(rescued),
            "final_selected": len(final_selected),
            "all_candidates": len(candidate_decisions),
        },
        "source_snapshot": dict(source_snapshot),
    }


def runtime_modules(workspace: Path):
    casefit_root = workspace / "RAG_re_casefit"
    clinical_root = workspace / "RAG_re_clinical"
    for root in (casefit_root, clinical_root):
        if not root.is_dir():
            raise FinalizerError(f"missing required package: {root}")
        if str(root) not in sys.path:
            sys.path.insert(0, str(root))
    from rag_re_casefit.rescue_evaluation import FROZEN_SYNONYMS
    from rag_re_clinical.clinical_rules import evaluate_modules

    return resolved_aliases(FROZEN_SYNONYMS), evaluate_modules


def load_casefit_signals(
    casefit_root: Path,
    clinical_config_path: Path,
    aliases: Mapping[str, str],
    evaluate_modules,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    manifest_path = casefit_root / "manifest.json"
    manifest = load_json(manifest_path)
    if manifest.get("schema_version") != "rag_re_casefit.batch_manifest.v1":
        raise FinalizerError(f"unexpected casefit manifest: {manifest_path}")
    if manifest.get("generation_is_gold_free") is not True:
        raise FinalizerError(f"casefit manifest is not marked gold-free: {manifest_path}")
    if manifest.get("errors"):
        raise FinalizerError(f"casefit manifest contains errors: {manifest_path}")
    outputs = manifest.get("outputs")
    if not isinstance(outputs, list) or len(outputs) != manifest.get("output_count"):
        raise FinalizerError(f"casefit output count mismatch: {manifest_path}")
    merge_dir = Path(manifest["merge_dir"]).resolve()
    clinical_config = load_json(clinical_config_path)
    if clinical_config.get("schema_version") != "rag_re_clinical.config.v1":
        raise FinalizerError(f"unexpected clinical config: {clinical_config_path}")

    result: dict[str, list[dict[str, Any]]] = defaultdict(list)
    merge_cache: dict[Path, dict[str, Any]] = {}
    for output in outputs:
        casefit_file = (casefit_root / output["output_file"]).resolve()
        if sha256_file(casefit_file) != output.get("output_sha256"):
            raise FinalizerError(f"casefit hash mismatch: {casefit_file}")
        artifact = load_json(casefit_file)
        if artifact.get("schema_version") != "rag_re_casefit.candidate_output.v1":
            raise FinalizerError(f"unexpected casefit schema: {casefit_file}")
        patient_id = str(artifact.get("patient_id") or "")
        review_tier = str(artifact.get("review_tier") or "")
        if review_tier not in VALID_REVIEW_TIERS:
            raise FinalizerError(f"candidate outside R5 review scope: {casefit_file}")
        provenance = artifact.get("provenance") or {}
        merge_file = (merge_dir / provenance.get("latest_merge_file", "")).resolve()
        if merge_file.parent != merge_dir:
            raise FinalizerError(f"unsafe merge provenance path: {merge_file}")
        if sha256_file(merge_file) != provenance.get("latest_merge_sha256"):
            raise FinalizerError(f"latest merge hash mismatch: {merge_file}")
        payload = merge_cache.setdefault(merge_file, load_json(merge_file))
        if str(payload.get("patient_id") or patient_id) != patient_id:
            raise FinalizerError(f"casefit/merge patient mismatch: {casefit_file}")

        display_name = str(artifact.get("candidate_organism") or "").strip()
        canonical_name = str(artifact.get("canonical_organism") or display_name).strip()
        organism_key = canonical(canonical_name, aliases)
        review = payload.get("llm_missed_candidate_review") or {}
        review_matches: list[tuple[str, Mapping[str, Any]]] = []
        for tier in sorted(VALID_REVIEW_TIERS):
            rows = review.get(tier)
            if not isinstance(rows, list):
                raise FinalizerError(f"P{patient_id} missing review tier {tier}")
            for row in rows:
                if isinstance(row, dict) and canonical(candidate_name(row), aliases) == organism_key:
                    review_matches.append((tier, row))
        if len(review_matches) != 1 or review_matches[0][0] != review_tier:
            raise FinalizerError(
                f"P{patient_id} {canonical_name}: expected one matching latest review candidate"
            )
        review_row = review_matches[0][1]

        deterministic_rows = (
            (payload.get("deterministic_max") or {}).get("pathogen_candidates")
        )
        if not isinstance(deterministic_rows, list):
            raise FinalizerError(f"P{patient_id} missing deterministic candidates")
        deterministic_matches = [
            row
            for row in deterministic_rows
            if isinstance(row, dict)
            and canonical(candidate_name(row), aliases) == organism_key
        ]
        exact_matches = [
            row
            for row in deterministic_matches
            if normal(candidate_name(row)) == normal(canonical_name)
        ]
        if len(exact_matches) == 1:
            deterministic = exact_matches[0]
        elif len(deterministic_matches) == 1:
            deterministic = deterministic_matches[0]
        elif not deterministic_matches:
            deterministic = None
        else:
            raise FinalizerError(
                f"P{patient_id} {canonical_name}: ambiguous deterministic candidates"
            )
        snapshot = review_row.get("evidence_snapshot")
        snapshot = snapshot if isinstance(snapshot, dict) else {}
        effective = dict(deterministic or snapshot)
        effective["organism_name"] = canonical_name
        effective["canonical_organism_name"] = canonical_name
        case_card = artifact.get("case_card") or {}
        source_context = {
            "site_code": case_card.get("clinical_site"),
            "target_site": case_card.get("clinical_site"),
            "specimen_site_code": case_card.get("specimen_site"),
            "specimen_class": effective.get("specimen_class"),
            "specimen_alignment": effective.get("specimen_alignment"),
            "dominant_source": case_card.get("target_syndrome"),
        }
        merged_context = {
            "source_context": source_context,
            "effective_candidate": effective,
            "deterministic_candidate": deterministic,
            "review_metadata": {"tier": review_tier, "row": review_row},
        }
        modules = evaluate_modules(merged_context, effective, clinical_config)
        b = modules.get("B_STRICT") or {}
        c = modules.get("C") or {}
        d = modules.get("D") or {}
        c_grade = str(c.get("grade") or "").strip()
        if c_grade not in {"C0", "C1", "C2", "C3", "CNEG", "CU", "CX"}:
            raise FinalizerError(f"invalid C grade for P{patient_id} {canonical_name}")
        a1 = artifact.get("A1_aggregate") or {}
        strong = a1.get("strong_match_count")
        partial = a1.get("partial_match_count")
        if type(strong) is not int or type(partial) is not int or strong < 0 or partial < 0:
            raise FinalizerError(f"invalid A1 counts: {casefit_file}")
        signal = {
            "candidate_organism": display_name,
            "canonical_organism": canonical_name,
            "review_tier": review_tier,
            "A1_strong": strong,
            "A1_partial": partial,
            "A1_match": strong + partial,
            "B_strict": strict_bool_or_none(
                b.get("positive"), f"P{patient_id}.{canonical_name}.B_strict"
            ),
            "B_absolute_axis": strict_bool_or_none(
                b.get("absolute_strength_axis"),
                f"P{patient_id}.{canonical_name}.B_absolute_axis",
            ),
            "B_relative_axis": strict_bool_or_none(
                b.get("relative_strength_axis"),
                f"P{patient_id}.{canonical_name}.B_relative_axis",
            ),
            "B_site_aligned": strict_bool_or_none(
                b.get("site_aligned"), f"P{patient_id}.{canonical_name}.B_site_aligned"
            ),
            "C_grade": c_grade,
            "D_state": str(d.get("state") or "").strip(),
            "clinical_features_source": (
                "latest_merge_deterministic_candidate"
                if deterministic is not None
                else "latest_merge_review_evidence_snapshot"
            ),
            "casefit_file": casefit_file.name,
            "casefit_sha256": sha256_file(casefit_file),
            "latest_merge_file": merge_file.name,
            "latest_merge_sha256": sha256_file(merge_file),
        }
        signal.update(evaluate_r5(signal))
        result[patient_id].append(signal)
    return dict(result), {
        "manifest_file": manifest_path.name,
        "manifest_sha256": sha256_file(manifest_path),
        "merge_dir_name": merge_dir.name,
        "clinical_config_file": clinical_config_path.name,
        "clinical_config_sha256": sha256_file(clinical_config_path),
        "candidate_count": len(outputs),
    }


def resolve_from(base: Path, value: str) -> Path:
    return (base / value).resolve()


def make_all_decisions(
    workspace: Path, config_path: Path, clinical_config_path: Path
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    aliases, evaluate_modules = runtime_modules(workspace)
    config = load_json(config_path)
    config_base = config_path.parent
    accepted_path = resolve_from(config_base, config["accepted_rationales"])
    accepted = load_json(accepted_path)
    if accepted.get("status") != "structurally_validated":
        raise FinalizerError("accepted rationale collection is not structurally validated")
    entries: dict[tuple[str, str], dict[str, Any]] = {}
    for entry in accepted.get("patients", []):
        cohort = str(entry.get("cohort") or "")
        patient_id = str((entry.get("result") or {}).get("patient_id") or "")
        key = (cohort, patient_id)
        if key in entries:
            raise FinalizerError(f"duplicate patient entry: {key}")
        entries[key] = entry

    all_decisions: list[dict[str, Any]] = []
    cohort_meta: dict[str, Any] = {}
    for cohort, cohort_config in config.get("cohorts", {}).items():
        casefit_root = resolve_from(config_base, cohort_config["casefit_root"])
        signals_by_patient, signal_meta = load_casefit_signals(
            casefit_root,
            clinical_config_path,
            aliases,
            evaluate_modules,
        )
        cohort_meta[cohort] = signal_meta
        selected_root = resolve_from(config_base, cohort_config["selected_root"])
        cohort_entries = sorted(
            ((key, value) for key, value in entries.items() if key[0] == cohort),
            key=lambda pair: int(pair[0][1]),
        )
        for (_, patient_id), entry in cohort_entries:
            result = entry["result"]
            override = (config.get("rationale_overrides") or {}).get(
                f"{cohort}/{patient_id}"
            )
            if override:
                result = load_json(resolve_from(config_base, override))
            if str(result.get("patient_id")) != patient_id:
                raise FinalizerError(f"rationale patient mismatch: {cohort}/P{patient_id}")
            chain_path = Path(result["input_source"]).resolve()
            chain = load_json(chain_path)
            chain_patient = str((chain.get("病例識別") or {}).get("病人編號") or "")
            if chain_patient != patient_id:
                raise FinalizerError(f"chain patient mismatch: {cohort}/P{patient_id}")
            candidate_rows = chain.get("步驟2至5_逐候選病原判讀")
            if not isinstance(candidate_rows, list):
                raise FinalizerError(f"missing candidate rows: {cohort}/P{patient_id}")
            candidate_names = [str(row.get("病原") or "") for row in candidate_rows]
            selected_path = selected_root / f"NGS_patient_{patient_id}_selected_pathogens.json"
            selected = load_json(selected_path)
            if str(selected.get("patient_id")) != patient_id:
                raise FinalizerError(f"selected patient mismatch: {cohort}/P{patient_id}")
            source = selected.get("source") or {}
            source_file = Path(source.get("file", "")).resolve()
            if sha256_file(source_file) != source.get("sha256"):
                raise FinalizerError(f"changed upstream selection source: {cohort}/P{patient_id}")
            signals = signals_by_patient.get(patient_id, [])
            for signal in signals:
                if signal["latest_merge_sha256"] != source.get("sha256"):
                    raise FinalizerError(
                        f"casefit/selection source mismatch: {cohort}/P{patient_id}"
                    )
            decision = compose_decision(
                cohort=cohort,
                patient_id=patient_id,
                candidate_names=candidate_names,
                picked=selected.get("selected_pathogens") or [],
                signals=signals,
                source_snapshot={
                    "reasoning_chain_file": chain_path.name,
                    "reasoning_chain_sha256": sha256_file(chain_path),
                    "selected_pathogens_file": selected_path.name,
                    "selected_pathogens_sha256": sha256_file(selected_path),
                    "latest_merge_file": source_file.name,
                    "latest_merge_sha256": source.get("sha256"),
                    "casefit_files": [
                        {
                            "file": signal["casefit_file"],
                            "sha256": signal["casefit_sha256"],
                        }
                        for signal in signals
                    ],
                    "clinical_config_file": clinical_config_path.name,
                    "clinical_config_sha256": sha256_file(clinical_config_path),
                },
            )
            all_decisions.append(decision)
        unused = sorted(set(signals_by_patient) - {key[1] for key, _ in cohort_entries})
        if unused:
            raise FinalizerError(f"casefit patients absent from delivery cohort {cohort}: {unused}")
    if set(entries) != {(d["cohort"], d["patient_id"]) for d in all_decisions}:
        raise FinalizerError("not all accepted patients were assigned to a configured cohort")
    return all_decisions, {
        "delivery_config_file": config_path.name,
        "delivery_config_sha256": sha256_file(config_path),
        "accepted_rationales_file": accepted_path.name,
        "accepted_rationales_sha256": sha256_file(accepted_path),
        "cohorts": cohort_meta,
    }


def render_summary(decisions: list[Mapping[str, Any]]) -> str:
    picked = sum(item["counts"]["upstream_picked"] for item in decisions)
    rescued = sum(item["counts"]["R5_rescued"] for item in decisions)
    lines = [
        "# OBER R5 計算最終選擇清單",
        "",
        f"規則：`{RULE_FORMULA}`",
        "",
        f"共 {len(decisions)} 位病人；picked {picked} 筆，R5 新增 {rescued} 筆，合計 {picked + rescued} 筆。",
        "",
        "> 本次執行沒有讀 gold answer；但 R5 是在先前看過開發答案後提出，不能把這批結果當成 untouched test performance。",
        "",
        "## R5 新增",
        "",
        "| 資料組 | 病人 | 新順位 | 病原 | 通過路徑 | A1 strong | A1 match | site | B abs | B rel | C |",
        "|---|---:|---:|---|---|---:|---:|---|---|---|---|",
    ]
    rescue_count = 0
    for decision in decisions:
        candidate_by_key = {
            normal(item["organism"]): item for item in decision["candidates"]
        }
        for selected in decision["selected"]:
            if selected["selection_origin"] != "OBER_R5_rescue":
                continue
            rescue_count += 1
            signal = candidate_by_key[normal(selected["organism"])]["R5_signals"]
            lines.append(
                f"| {decision['cohort']} | P{decision['patient_id']} | {selected['rank']} | "
                f"{selected['organism']} | {selected['rule_path']} | {signal['A1_strong']} | "
                f"{signal['A1_match']} | {signal['B_site_aligned']} | "
                f"{signal['B_absolute_axis']} | {signal['B_relative_axis']} | {signal['C_grade']} |"
            )
    if not rescue_count:
        lines.append("| — | — | — | 本批無新增 | — | — | — | — | — | — | — |")
    lines.extend(
        [
            "",
            "## 每位病人最終清單",
            "",
            "| 資料組 | 病人 | picked | R5 新增 | 最終數 | 最終病原（依順位） |",
            "|---|---:|---:|---:|---:|---|",
        ]
    )
    for decision in decisions:
        selected_text = "；".join(
            f"{item['rank']}. {item['organism']}" for item in decision["selected"]
        ) or "無入選病原"
        lines.append(
            f"| {decision['cohort']} | P{decision['patient_id']} | "
            f"{decision['counts']['upstream_picked']} | {decision['counts']['R5_rescued']} | "
            f"{decision['counts']['final_selected']} | {selected_text} |"
        )
    lines.extend(
        [
            "",
            "## 使用限制",
            "",
            "- picked 完全保留；R5 不負責剪枝。",
            "- 新增候選接在 picked 後；若同病人同時新增多個，字母排序只為固定輸出，不代表臨床信心高低。",
            "- D 狀態保留供查核，但不是隱藏門檻。",
            "- 這是計算最終清單，尚未針對新增病原重跑理由產生器，也未經醫師複核。",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(
    output: Path,
    decisions: list[dict[str, Any]],
    input_snapshot: Mapping[str, Any],
) -> dict[str, Any]:
    if output.exists() and any(output.iterdir()):
        raise FinalizerError("output directory must be new or empty")
    output.mkdir(parents=True, exist_ok=True)
    file_rows: list[dict[str, Any]] = []
    selection_rows: list[dict[str, Any]] = []
    for decision in decisions:
        cohort_dir = output / decision["cohort"]
        cohort_dir.mkdir(parents=True, exist_ok=True)
        path = cohort_dir / f"P{decision['patient_id']}_final_decision.json"
        write_json(path, decision)
        file_rows.append(
            {
                "cohort": decision["cohort"],
                "patient_id": decision["patient_id"],
                "file": path.relative_to(output).as_posix(),
                "sha256": sha256_file(path),
                **decision["counts"],
            }
        )
        selection_rows.append(
            {
                "cohort": decision["cohort"],
                "patient_id": decision["patient_id"],
                **decision["counts"],
                "selected": decision["selected"],
            }
        )
    summary_json = output / "final_selection_list.json"
    write_json(
        summary_json,
        {
            "schema_version": "ober.final_selection_list.v1",
            "rule_id": RULE_ID,
            "formula": RULE_FORMULA,
            "gold_answer_used_in_this_run": False,
            "patients": selection_rows,
        },
    )
    summary_md = output / "final_selection_list.markdown"
    summary_md.write_text(render_summary(decisions), encoding="utf-8")
    totals = {
        "patients": len(decisions),
        "upstream_picked": sum(d["counts"]["upstream_picked"] for d in decisions),
        "R5_rescued": sum(d["counts"]["R5_rescued"] for d in decisions),
        "final_selected": sum(d["counts"]["final_selected"] for d in decisions),
        "all_candidates": sum(d["counts"]["all_candidates"] for d in decisions),
    }
    manifest = {
        "schema_version": MANIFEST_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "stage": "workflow_final",
        "decision_status": "computational_final_pending_rationale_and_clinical_review",
        "rule_id": RULE_ID,
        "formula": RULE_FORMULA,
        "gold_answer_used_in_this_run": False,
        "development_rule_answers_previously_inspected": True,
        "totals": totals,
        "input_snapshot": dict(input_snapshot),
        "decision_files": file_rows,
        "summary_files": {
            "json": {
                "file": summary_json.name,
                "sha256": sha256_file(summary_json),
            },
            "markdown": {
                "file": summary_md.name,
                "sha256": sha256_file(summary_md),
            },
        },
    }
    write_json(output / "manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate gold-free OBER R5 final-decision files from frozen picked and casefit inputs"
    )
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--delivery-config", required=True, type=Path)
    parser.add_argument("--clinical-config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    config = args.delivery_config.resolve()
    clinical = args.clinical_config.resolve()
    decisions, snapshot = make_all_decisions(workspace, config, clinical)
    manifest = write_outputs(args.output.resolve(), decisions, snapshot)
    print(json.dumps({"status": "complete", "output": str(args.output.resolve()), "totals": manifest["totals"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build a human-readable view and a traceable reasoning-chain view."""

from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


IDENTIFIER_FIELDS = (
    "specimen_code",
    "test_id",
    "accession_number",
    "accession_no",
    "sample_id",
    "order_id",
    "lab_no",
    "檢驗編號",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create simple selected-pathogen summaries and traceable reasoning chains."
    )
    parser.add_argument("--candidate-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--pretty", action="store_true")
    return parser.parse_args()


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, value: Any, pretty: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        if pretty:
            json.dump(value, handle, ensure_ascii=False, indent=2)
        else:
            json.dump(value, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")


def compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def strings(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    return [text for item in value if (text := compact(item))]


def iter_evidence(candidate: dict[str, Any]) -> Iterable[dict[str, Any]]:
    for item in candidate.get("candidate_sources", []):
        if isinstance(item, dict):
            yield item
    grouped = candidate.get("all_time_evidence")
    if not isinstance(grouped, dict):
        return
    for items in grouped.values():
        if not isinstance(items, list):
            continue
        for item in items:
            if isinstance(item, dict):
                yield item


def record_identifiers(record: Any) -> list[str]:
    if not isinstance(record, dict):
        return []
    values: list[str] = []
    for field in IDENTIFIER_FIELDS:
        value = compact(record.get(field))
        if value:
            values.append(value)
    return values


def test_identifiers(candidate: dict[str, Any]) -> list[str]:
    values: list[str] = []
    seen: set[str] = set()
    for evidence in iter_evidence(candidate):
        for value in record_identifiers(evidence.get("record")):
            if value not in seen:
                seen.add(value)
                values.append(value)
    return values


def evidence_locator(evidence: dict[str, Any]) -> dict[str, Any]:
    source_file = evidence.get("source_file")
    record_index = evidence.get("record_index")
    locator = (
        f"{source_file}#{record_index}"
        if source_file is not None and record_index is not None
        else compact(source_file)
    )
    return {
        "locator": locator or None,
        "record_type": evidence.get("record_type"),
        "source_file": source_file,
        "record_index": record_index,
        "time": evidence.get("time"),
        "time_field": evidence.get("time_field"),
        "microbe_name": evidence.get("microbe_name"),
        "result": evidence.get("result"),
        "record": evidence.get("record"),
    }


def unique_evidence(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw in iter_evidence(candidate):
        item = evidence_locator(raw)
        key = json.dumps(item, ensure_ascii=False, sort_keys=True)
        if key in seen:
            continue
        seen.add(key)
        output.append(item)
    return output


def evidence_summary(candidate: dict[str, Any]) -> list[str]:
    summaries: list[str] = []
    for item in unique_evidence(candidate):
        record_type = compact(item.get("record_type"))
        record_type_label = {
            "all_RK_NTC_microbes": "mNGS 完整菌表",
            "mNGS_grouped": "mNGS 分組結果",
            "culture": "培養",
            "molecular_microbiology": "分子微生物檢驗",
            "filmarray": "FilmArray",
            "gm_test": "GM／抗原檢驗",
        }.get(record_type, record_type)
        parts = [record_type_label]
        if compact(item.get("microbe_name")):
            parts.append(compact(item.get("microbe_name")))
        if compact(item.get("result")):
            parts.append(compact(item.get("result")))
        if compact(item.get("time")):
            parts.append(compact(item.get("time")))
        summaries.append(" / ".join(part for part in parts if part))
    return summaries


def provenance(candidate: dict[str, Any]) -> dict[str, Any]:
    value = candidate.get("selection_provenance")
    return value if isinstance(value, dict) else {}


def consideration_reasons(candidate: dict[str, Any]) -> list[str]:
    source = provenance(candidate)
    reasons = strings(source.get("consideration_reasons"))
    if reasons:
        return reasons
    summaries = evidence_summary(candidate)
    if summaries:
        return [f"病例資料中存在與此病原相符的檢驗證據：{text}" for text in summaries]
    return ["上游將此病原列入候選，但目前病例資料包沒有找到可直接定位的同名紀錄。"]


def decision_reasons(candidate: dict[str, Any]) -> list[str]:
    source = provenance(candidate)
    reasons = strings(source.get("decision_reasons"))
    if not reasons:
        reasons = strings(source.get("selection_reasons"))
    return reasons


REASON_TRANSLATIONS = {
    "Related representative-group hospital evidence is present but is not treated as same-organism direct support.":
        "醫院端有同一代表群的相關證據，但不視為同一菌種的直接支持。",
    "Typical pneumonia bacterium has lower-respiratory culture/FilmArray support plus clean mNGS signal; promote to M2_moderate so it is not lost because the hospital label and mNGS species label differ.":
        "此為典型肺炎病原；下呼吸道培養／FilmArray 有支持，且 mNGS 訊號乾淨，因此提升為 M2（中度支持），避免因醫院標籤與 mNGS 物種名稱不同而漏選。",
    "Respiratory Candida/yeast meets context-retention criteria and is capped at M3_protected / Level 3. Respiratory Candida alone is not formal picked unless same-species blood, sterile-site, tissue, pleural, abscess, or histopathology evidence is present.":
        "呼吸道 Candida／酵母菌符合保留脈絡的條件，但最高只列為 M3（保護性保留）／Level 3；若只有呼吸道 Candida，除非另有同種血液、無菌部位、組織、胸水、膿瘍或病理證據，否則不單獨視為正式感染源。",
    "Same-organism culture evidence is present only in a nonpulmonary source, so it is retained as patient context but does not support pulmonary causality.":
        "同菌培養證據只來自非肺部檢體，因此僅保留為病例脈絡，不能支持它是肺部感染源。",
    "Background-risk organism such as Corynebacterium has same-organism hospital support plus high-burden mNGS signal; retain as a secondary Level 3 candidate without requiring repeated culture.":
        "此類背景風險菌（如 Corynebacterium）有醫院端同菌支持，且 mNGS 菌量訊號高，因此保留為次要 Level 3 候選，不強制要求重複培養。",
    "Oral/aspiration flora has same-organism lower-respiratory hospital support or D2/D3 dominant mNGS signal, so it is retained but capped at M3_protected / Level 3.":
        "此口腔／吸入相關菌有下呼吸道同菌支持，或在 mNGS 中呈 D2／D3 優勢訊號，因此保留，但最高列為 M3（保護性保留）／Level 3。",
    "Rare opportunistic yeast has same-organism lower-respiratory hospital support, host vulnerability, and high-burden or high-percentile clean mNGS signal. It is retained as a Level 3 formal candidate instead of being handled by the respiratory Candida colonization guardrail.":
        "此罕見伺機性酵母菌同時具有下呼吸道同菌支持、宿主易感性，以及高菌量或高百分位的乾淨 mNGS 訊號，因此保留為正式 Level 3 候選，不套用一般呼吸道 Candida 定植規則。",
    "SARS-CoV-2/COVID-19 has high-clean mNGS signal; retain as picked or strong viral candidate.":
        "SARS-CoV-2／COVID-19 具有高且乾淨的 mNGS 訊號，因此保留為入選或強病毒候選。",
    "Candida/yeast has blood or sterile-site hospital evidence, so it is not treated as respiratory Candida colonization. This preserves invasive/systemic Candida context but does not by itself prove Candida pneumonia.":
        "Candida／酵母菌有血液或無菌部位證據，因此不當成單純呼吸道定植；這支持保留侵襲性／全身性 Candida 的可能性，但不能單獨證明 Candida 肺炎。",
    "mNGS-only D0/D1 non-dominant signal lacks same-organism hospital support and is outside protected, common hospital pneumonia, core respiratory virus, deferred Enterococcus/C. difficile, and Acinetobacter answer-sensitive groups. Keep out of formal picked_pathogens; use review tiering by signal strength.":
        "只有 mNGS 訊號，且為 D0／D1 非優勢型；沒有醫院端同菌支持，也不屬於需保護保留的常見肺炎病原或核心呼吸道病毒，因此不列入正式結果，只按訊號強度保留供人工複核。",
    "Level 3 mNGS-only signal has no pathogen-specific hospital support, R0/R1 reads, and D0/D1 non-dominant pattern; keep for audit/review context rather than formal picked_pathogens.":
        "此為 Level 3 的 mNGS-only 訊號：沒有醫院端同病原支持、reads 僅 R0／R1，且為 D0／D1 非優勢型，因此只保留供稽核／複核，不列入正式結果。",
    "Common hospital pneumonia pathogen has only R0/R1 mNGS signal, D0/D1 non-dominant pattern, and no direct same-organism culture, FilmArray, or molecular support. Keep it out of formal picked_pathogens; retain as review/audit context depending on related hospital context.":
        "雖屬常見院內肺炎病原，但目前只有 R0／R1 的低量 mNGS 訊號、呈 D0／D1 非優勢型，且沒有同菌培養、FilmArray 或分子檢驗支持，因此不列入正式結果，只保留供複核。",
    "Respiratory Candida/yeast is retained as Level 3 context but excluded from formal picked_pathogens because there is no same-species blood, sterile-site, tissue, pleural, abscess, or histopathology evidence.":
        "呼吸道 Candida／酵母菌僅保留為 Level 3 脈絡；因缺乏同種血液、無菌部位、組織、胸水、膿瘍或病理證據，所以不列入正式結果。",
    "mNGS-only R0/R1 host-context signal is D0/D1 non-dominant and is not in protected, answer-sensitive respiratory virus, or common hospital pneumonia pathogen groups. Keep it out of formal picked_pathogens.":
        "只有 R0／R1 的低量 mNGS 宿主脈絡訊號，且為 D0／D1 非優勢型；也不屬於需保護的呼吸道病毒或常見院內肺炎病原，因此不列入正式結果。",
    "Non-core respiratory virus detection is clinically relevant but often reflects recent infection, shedding, coinfection, or syndrome context rather than a formal single-causative picked pathogen. Keep it RAG-visible/context instead of formal picked.":
        "此非核心呼吸道病毒的檢出仍有臨床意義，但可能代表近期感染、持續排毒、共感染或症候群脈絡，不能直接當成單一感染源，因此只保留供後續文獻／人工判讀。",
    "Mold/opportunistic fungal signal is host-context only, R0/R1/R2, and D0/D1 non-dominant without culture, GM, molecular, sterile-site, or related-representative support. Keep it out of formal picked_pathogens and leave for review/RAG context if needed.":
        "黴菌／伺機性真菌目前只有宿主脈絡下的 mNGS 訊號，reads 為 R0–R2、呈 D0／D1 非優勢型，且沒有培養、GM、分子檢驗、無菌部位或相關代表菌支持，因此不列入正式結果，必要時交由文獻／人工複核。",
    "HSV/CMV/VZV has rank <=2, R2/R3 reads, high reads percentile, host vulnerability, and D0/D1 non-dominant host-only support. Keep RAG-visible as Level 3 context, but exclude from formal picked_pathogens without PCR/viral-load/direct hospital support.":
        "HSV／CMV／VZV 雖具有前二順位、R2／R3 reads、高百分位與宿主易感性，但仍是 D0／D1 非優勢且只有宿主脈絡支持；在沒有 PCR、病毒量或醫院端直接支持前，只列為 Level 3 複核對象，不列入正式結果。",
    "Generic unidentified hospital-side label is kept as culture/microbiology evidence but excluded from formal picked_pathogens until species-level identification is available.":
        "醫院端只有未鑑定到種的廣義標籤，因此保留為培養／微生物證據；在取得 species-level 鑑定前，不獨立列入正式病原。",
    "Water/environmental low-pulmonary-specificity organism lacks same-organism culture, FilmArray, molecular, or related-representative hospital support. Even with a strong mNGS burden, it should be reviewed as context rather than formal picked without local support.":
        "此水生／環境菌對肺部感染的特異性較低，且缺乏同菌培養、FilmArray、分子檢驗或相關代表菌支持；即使 mNGS 菌量高，在沒有本地檢驗支持前仍只保留作複核脈絡。",
}


def plain_language_reason(reason: str) -> str:
    text = compact(reason)
    for source, translated in REASON_TRANSLATIONS.items():
        text = text.replace(source, translated)
    text = text.replace("deterministic scorer：", "確定性規則評分：")
    text = text.replace("mNGS ranked/chosen 候選池", "mNGS 排名／入選候選池")
    text = text.replace(
        "以 hospital_only 候選納入 final pick",
        "以僅由醫院檢驗支持的候選納入最終結果",
    )
    text = text.replace("hospital_only 候選", "僅由醫院檢驗支持的候選")
    text = text.replace("final pick", "最終結果")
    text = text.replace("molecular_microbiology", "分子微生物檢驗")
    text = text.replace("filmarray_gmtest", "FilmArray／GM test")
    text = re.sub(r"\bculture\b", "培養", text)
    text = text.replace("rank_priority=", "候選順位=")
    text = text.replace("reads_tier=", "reads 分層=")
    text = text.replace("reads_percentile=", "同病例 reads 百分位=")
    text = text.replace(" -> ", "，判為 ")
    text = text.replace("mNGS-only", "僅有 mNGS")
    text = text.replace("species-level", "物種層級")
    for source, translated in {
        "R0_trace": "R0（微量）",
        "R1_low": "R1（低）",
        "R2_medium": "R2（中）",
        "R3_high": "R3（高）",
        "R4_very_high": "R4（很高）",
        "M1_strong": "M1（強支持）",
        "M2_moderate": "M2（中度支持）",
        "M3_protected": "M3（保護性保留）",
        "M4_weak": "M4（弱支持）",
    }.items():
        text = text.replace(source, translated)
    representative_match = re.fullmatch(
        r"Same representative group \((.+?)\) already has a stronger formal picked representative: "
        r"(.+?)\. Keep (.+?) as related/audit evidence rather than a separate picked pathogen\.",
        text,
    )
    if representative_match:
        group, stronger, retained = representative_match.groups()
        return (
            f"同一代表群（{group}）已有更強且正式入選的代表病原：{stronger}；"
            f"因此 {retained} 只保留為相關／稽核證據，不重複列為另一個入選病原。"
        )
    if text == "Insufficient evidence":
        return "現有證據不足。"
    return text


def selected_rank(candidate: dict[str, Any]) -> int | None:
    source = provenance(candidate)
    value = source.get("selected_rank", source.get("rank"))
    return value if isinstance(value, int) else None


def simple_reason(candidate: dict[str, Any]) -> str:
    reasons = decision_reasons(candidate)
    if reasons:
        return "；".join(plain_language_reason(reason) for reason in reasons)
    return "上游結果將此病原列為 picked，但未提供可直接使用的文字理由。"


def simple_row(bundle: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    selected_input = bundle.get("selected_input") or {}
    identifiers = test_identifiers(candidate)
    source = provenance(candidate)
    return {
        "hospital": selected_input.get("hospital") or "未提供",
        "patient_id": str(bundle.get("patient_id") or ""),
        "test_identifiers": identifiers,
        "selected_organism": candidate.get("name"),
        "reason": simple_reason(candidate),
        "picked_role": source.get("picked_role"),
        "evidence_source": source.get("evidence_source"),
        "selected_rank": selected_rank(candidate),
        "selection_status": "picked",
    }


def no_selection_row(bundle: dict[str, Any]) -> dict[str, Any]:
    selected_input = bundle.get("selected_input") or {}
    identifiers: list[str] = []
    seen: set[str] = set()
    for candidate in bundle.get("candidates", []):
        if not isinstance(candidate, dict):
            continue
        for identifier in test_identifiers(candidate):
            if identifier not in seen:
                seen.add(identifier)
                identifiers.append(identifier)
    return {
        "hospital": selected_input.get("hospital") or "未提供",
        "patient_id": str(bundle.get("patient_id") or ""),
        "test_identifiers": identifiers,
        "selected_organism": "無入選病原",
        "reason": "上游的 picked_pathogens 為空；候選病原與排除理由請查看同病人的可追溯推理鏈。",
        "picked_role": None,
        "evidence_source": None,
        "selected_rank": None,
        "selection_status": "none_picked",
    }


def candidate_chain(candidate: dict[str, Any]) -> dict[str, Any]:
    status = candidate.get("selection_status")
    reasons = decision_reasons(candidate)
    return {
        "organism": candidate.get("name"),
        "selection_status": status,
        "candidate_rank": provenance(candidate).get("candidate_rank"),
        "selected_rank": selected_rank(candidate),
        "why_considered": consideration_reasons(candidate),
        "why_selected_or_excluded": reasons,
        "plain_language_decision_summary": [plain_language_reason(reason) for reason in reasons],
        "decision_reason_status": "available" if reasons else "not_provided_by_upstream",
        "upstream_features": provenance(candidate).get("upstream_features", {}),
        "internal_case_evidence": unique_evidence(candidate),
        "external_evidence": {
            "status": "pending_ober",
            "items": [],
            "note": "No literature, case-report, or website evidence is fabricated by this tool.",
        },
    }


def first_present(record: dict[str, Any], fields: tuple[str, ...]) -> Any:
    for field in fields:
        value = record.get(field)
        if value not in (None, "", [], {}):
            return value
    return None


def evidence_interpretation(record_type: str, record: dict[str, Any]) -> str:
    specimen = compact(
        first_present(
            record,
            ("specimen_site", "sample", "specimen", "sample_type", "source"),
        )
    ).casefold()
    if record_type in {"all_RK_NTC_microbes", "mNGS_grouped", "mngs"}:
        if any(term in specimen for term in ("balf", "bal", "sputum", "trache", "bronch")):
            return "下呼吸道 mNGS 檢出此病原；它可直接支持肺部相關性，但仍需結合菌量、優勢程度與其他檢驗。"
        if "blood" in specimen:
            return "血液 mNGS 檢出此病原；可支持系統性感染脈絡，但不能單獨證明它是肺部感染源。"
        return "mNGS 檢出此病原；是否為感染源仍需結合檢體部位、菌量、優勢程度與其他檢驗。"
    if record_type == "culture":
        if any(term in specimen for term in ("urine", "foley")):
            return "尿液培養檢出同名病原；支持泌尿道／全身感染脈絡，但不直接支持肺部感染源。"
        if "blood" in specimen:
            return "血液培養檢出同名病原；支持侵襲性／系統性感染，但仍需判斷是否與肺部事件相關。"
        if any(term in specimen for term in ("balf", "bal", "sputum", "trache", "bronch")):
            return "下呼吸道培養檢出同名病原，是支持肺部感染來源的病例內直接證據。"
        return "培養檢出同名病原；需依檢體部位判斷它對肺部感染的支持程度。"
    if record_type in {"molecular_microbiology", "filmarray", "gm_test"}:
        return "醫院端直接病原／抗原檢驗結果；需依檢體部位與檢驗方法判斷其肺部特異性。"
    if record_type == "image":
        return "影像可支持肺部感染事件存在，但不能單獨確認是哪一個病原。"
    return "這是病例內可回查的相關紀錄；其臨床意義需和其他證據一起判讀。"


def readable_evidence(candidate: dict[str, Any]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for evidence in unique_evidence(candidate):
        record = evidence.get("record") if isinstance(evidence.get("record"), dict) else {}
        record_type = compact(evidence.get("record_type"))
        specimen = first_present(
            record,
            ("specimen_site", "sample", "specimen", "sample_type"),
        )
        microbe = evidence.get("microbe_name") or first_present(
            record,
            ("organism", "name", "target", "pathogen", "microbe"),
        )
        result = evidence.get("result") or first_present(
            record,
            ("result", "status", "interpretation", "value"),
        )
        identifiers = record_identifiers(record)
        reads = record.get("reads")
        output.append(
            {
                "證據摘要": "｜".join(
                    part
                    for part in (
                        compact(evidence.get("time")),
                        record_type,
                        compact(specimen),
                        compact(microbe),
                        f"結果={compact(result)}"
                        if compact(result)
                        else "",
                        f"reads={reads}" if reads not in (None, "") else "",
                    )
                    if part
                ),
                "如何解讀": evidence_interpretation(record_type, record),
                "檢驗編號": identifiers or ["未提供"],
                "原始位置": evidence.get("locator") or "未提供",
                "可核對欄位": {
                    "時間": evidence.get("time") or "未提供",
                    "檢體": specimen or "未提供",
                    "病原／檢驗目標": microbe or "未提供",
                    "結果": result or "未提供",
                    "reads": reads if reads not in (None, "") else "未提供",
                },
            }
        )
    return output


def concise_features(candidate: dict[str, Any]) -> dict[str, Any]:
    features = provenance(candidate).get("upstream_features", {})
    if not isinstance(features, dict):
        features = {}
    labels = (
        ("integrated_causative_level", "綜合病原等級"),
        ("mngs_signal_tier", "mNGS 支持層級"),
        ("reads", "mNGS reads"),
        ("reads_tier", "reads 分層"),
        ("reads_percentile", "同病例 reads 百分位"),
        ("dominance_tier", "優勢程度"),
        ("specimen_alignment", "檢體部位對齊"),
        ("specimen_class", "檢體分類"),
        ("rank_priority", "候選順位"),
        ("possibility_level", "上游可能性"),
    )
    value_labels = {
        "M1_strong": "M1_strong（強支持）",
        "M2_moderate": "M2_moderate（中度支持）",
        "M3_protected": "M3_protected（保護性保留）",
        "M4_weak": "M4_weak（弱支持）",
        "R0_trace": "R0_trace（微量）",
        "R1_low": "R1_low（低）",
        "R2_medium": "R2_medium（中）",
        "R3_high": "R3_high（高）",
        "R4_very_high": "R4_very_high（很高）",
        "D0_not_top": "D0_not_top（非優勢）",
        "D1_top_not_dominant": "D1_top_not_dominant（最高但未形成優勢）",
        "D2_dominant": "D2_dominant（優勢）",
        "D3_highly_dominant": "D3_highly_dominant（高度優勢）",
        "Sterile_or_Systemic": "Sterile_or_Systemic（無菌部位或系統性）",
        "Lower_Respiratory": "Lower_Respiratory（下呼吸道）",
        "S1_sterile_systemic": "S1_sterile_systemic（無菌部位／系統性檢體）",
        "S2_lower_respiratory": "S2_lower_respiratory（下呼吸道檢體）",
        "high": "high（高）",
        "medium": "medium（中）",
        "low": "low（低）",
    }
    output: dict[str, Any] = {}
    for field, label in labels:
        value = features.get(field)
        if value in (None, "", [], {}):
            continue
        output[label] = value_labels.get(str(value), value)
    return output


def readable_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    status = candidate.get("selection_status")
    considered = list(dict.fromkeys(plain_language_reason(item) for item in consideration_reasons(candidate)))
    decisions = list(dict.fromkeys(plain_language_reason(item) for item in decision_reasons(candidate)))
    if not decisions:
        decisions = [
            "上游只記錄了最終未入選，但沒有提供明確排除理由；此處不能自行補造，需人工複核。"
            if status != "picked"
            else "上游將此病原列為入選，但沒有提供可直接閱讀的文字理由。"
        ]
    if status == "picked":
        positive = next(
            (
                item
                for item in decisions
                if any(
                    marker in item
                    for marker in ("判為", "納入最終結果", "提升為", "正式 Level")
                )
            ),
            decisions[-1],
        )
        caveats = [
            item
            for item in decisions
            if item != positive
            and any(
                marker in item
                for marker in ("不能支持", "不直接支持", "不能單獨證明", "只來自非肺部")
            )
        ]
        conclusion = f"選入：{positive}"
        if caveats:
            conclusion += f" 注意限制：{caveats[0]}"
    else:
        conclusion = f"未選入：{decisions[0]}"
    evidence_items = readable_evidence(candidate)
    evidence_limitations = list(
        dict.fromkeys(
            item["如何解讀"]
            for item in evidence_items
            if any(
                marker in item["如何解讀"]
                for marker in ("不能單獨證明", "不直接支持", "需依檢體部位", "仍需結合")
            )
        )
    )
    return {
        "病原": candidate.get("name"),
        "最終處置": "選入" if status == "picked" else "未選入",
        "最終排序": selected_rank(candidate) if status == "picked" else None,
        "一句話結論": conclusion,
        "為什麼會被考慮": considered,
        "為什麼選入或排除": decisions,
        "關鍵判斷資訊": concise_features(candidate),
        "重要限制": evidence_limitations or ["目前病例內證據沒有額外標示可自動摘要的限制。"],
        "病例內支持或反對證據": evidence_items,
        "外部文獻病例或網站": {
            "狀態": "尚未執行 OBER 文獻查核",
            "說明": "目前留空不代表沒有文獻支持；待 OBER 執行後填入標題、來源、URL 與支持／反對結論。",
            "證據": [],
        },
    }


def context_records(bundle: dict[str, Any], category: str) -> list[dict[str, Any]]:
    context = bundle.get("all_time_context")
    if not isinstance(context, dict):
        return []
    category_data = context.get(category)
    if not isinstance(category_data, dict):
        return []
    return [item for item in category_data.get("records", []) if isinstance(item, dict)]


def readable_patient_summary(bundle: dict[str, Any]) -> dict[str, Any]:
    diagnoses = []
    for record in context_records(bundle, "admission_diagnosis"):
        value = first_present(record, ("diagnosis", "admission_diagnosis", "value"))
        if value not in (None, ""):
            text = compact(value)
            text = text.replace("Pneumonia type:", "肺炎類型：")
            text = text.replace("Ward / unit:", "病房／單位：")
            diagnoses.append(text)
    image_findings: list[str] = []
    for record in context_records(bundle, "image"):
        findings = record.get("findings")
        if isinstance(findings, list):
            image_findings.extend(compact(item) for item in findings if compact(item))
        elif compact(findings):
            image_findings.append(compact(findings))
    pulmonary_terms = (
        "pneumonia",
        "lung",
        "infiltrat",
        "opacity",
        "consolidat",
        "ground-glass",
        "ggo",
        "airspace",
        "pleural effusion",
        "cp angle",
    )
    pulmonary_findings = [
        finding for finding in image_findings if any(term in finding.casefold() for term in pulmonary_terms)
    ]
    demographics: dict[str, Any] = {}
    underlying = context_records(bundle, "underlying")
    if underlying:
        demographics = {
            "年齡": underlying[0].get("age", "未提供"),
            "性別": underlying[0].get("gender", "未提供"),
        }
    candidates = [item for item in bundle.get("candidates", []) if isinstance(item, dict)]
    identifiers = list(
        dict.fromkeys(
            identifier
            for candidate in candidates
            for identifier in test_identifiers(candidate)
        )
    )
    return {
        "基本資料": demographics or "未提供",
        "住院診斷或照護脈絡": list(dict.fromkeys(diagnoses)) or ["未提供"],
        "影像重點（保留原文）": list(dict.fromkeys(pulmonary_findings or image_findings))[:6]
        or ["未提供"],
        "檢驗編號": identifiers or ["未提供"],
        "證據中心時間": (bundle.get("anchor") or {}).get("time") or "未提供",
        "納入證據時間範圍": bundle.get("candidate_window") or "未提供",
        "原始病例資料夾": bundle.get("input_patient_dir") or "未提供",
    }


def readable_reasoning_chain(
    bundle: dict[str, Any],
    technical_output_file: Path,
) -> dict[str, Any]:
    selected_input = bundle.get("selected_input") or {}
    candidates = [item for item in bundle.get("candidates", []) if isinstance(item, dict)]
    candidates.sort(
        key=lambda item: (
            item.get("selection_status") != "picked",
            selected_rank(item) is None,
            selected_rank(item) or 10**9,
            provenance(item).get("candidate_rank") or 10**9,
            compact(item.get("name")).casefold(),
        )
    )
    selected = [item for item in candidates if item.get("selection_status") == "picked"]
    final_ranking = [
        {
            "排序": selected_rank(item),
            "病原": item.get("name"),
            "選入理由": simple_reason(item),
        }
        for item in selected
    ]
    return {
        "文件名稱": "OBER 可解釋病原判讀鏈",
        "文件用途": (
            "讓醫師或研究者能從最終結果往回查到每個候選的判斷理由與原始病例證據；"
            "這是決策支援與研究稽核資料，不等同臨床確診。"
        ),
        "閱讀順序": [
            "步驟1：先看病人資料摘要與檢體脈絡",
            "步驟2：確認本病例有哪些候選病原",
            "步驟3：逐一閱讀為什麼考慮、選入或排除",
            "步驟4：查看支持或反對判斷的病例內證據",
            "步驟5：由原始位置回查來源；外部文獻待 OBER 補入 URL",
            "步驟6：查看最終鑑別診斷排序",
        ],
        "病例識別": {
            "醫院": selected_input.get("hospital") or "未提供",
            "病人編號": str(bundle.get("patient_id") or ""),
            "資料集": selected_input.get("dataset") or "未提供",
            "流程版本": selected_input.get("pipeline_version") or "未提供",
        },
        "一眼看懂": {
            "候選病原數": len(candidates),
            "入選病原數": len(selected),
            "最終入選病原": [item.get("name") for item in selected] or ["無入選病原"],
            "尚待補充": "OBER 外部文獻、病例報告與網站查核",
            "重要提醒": "此處的『選入』表示通過目前上游規則，不等同已由臨床醫師確診為感染源。",
        },
        "步驟1_病人資料摘要": readable_patient_summary(bundle),
        "步驟2至5_逐候選病原判讀": [readable_candidate(item) for item in candidates],
        "步驟6_最終鑑別診斷排序": final_ranking or [
            {"排序": None, "病原": "無入選病原", "選入理由": "所有候選均未通過目前規則。"}
        ],
        "完整技術追溯": {
            "用途": "需要查看完整原始欄位、上游 features 或程式稽核時再開啟。",
            "檔案": str(technical_output_file.resolve()),
        },
    }


def traceable_chain(bundle: dict[str, Any], source_file: Path) -> dict[str, Any]:
    selected_input = bundle.get("selected_input") or {}
    candidates = [item for item in bundle.get("candidates", []) if isinstance(item, dict)]
    chains = [candidate_chain(candidate) for candidate in candidates]
    final_ranking = [
        {
            "rank": chain["selected_rank"],
            "organism": chain["organism"],
            "selection_status": "picked",
        }
        for chain in sorted(
            (item for item in chains if item["selection_status"] == "picked"),
            key=lambda item: (item["selected_rank"] is None, item["selected_rank"] or 10**9),
        )
    ]
    return {
        "schema_version": "traceable_pathogen_reasoning_chain.v1",
        "hospital": selected_input.get("hospital") or "未提供",
        "patient_id": str(bundle.get("patient_id") or ""),
        "dataset": selected_input.get("dataset"),
        "selection_source": selected_input.get("selection_source"),
        "pipeline_version": selected_input.get("pipeline_version"),
        "patient_data": {
            "input_patient_dir": bundle.get("input_patient_dir"),
            "candidate_bundle": str(source_file.resolve()),
            "anchor": bundle.get("anchor"),
            "candidate_window": bundle.get("candidate_window"),
            "record_type_counts": bundle.get("record_type_counts"),
            "clinical_context": bundle.get("all_time_context"),
        },
        "candidate_diagnoses": chains,
        "final_differential_diagnosis_ranking": final_ranking,
        "traceability_status": {
            "internal_case_evidence": "available",
            "external_literature_evidence": "pending_ober",
            "final_ranking_source": "frozen_upstream_picked_pathogens",
        },
        "reasoning_chain_flow": [
            {"stage": 1, "name": "病人資料", "field": "patient_data"},
            {"stage": 2, "name": "候選診斷", "field": "candidate_diagnoses[].organism"},
            {
                "stage": 3,
                "name": "考慮／排除理由",
                "field": "candidate_diagnoses[].why_considered / why_selected_or_excluded",
            },
            {
                "stage": 4,
                "name": "支持來源",
                "field": "candidate_diagnoses[].internal_case_evidence / external_evidence",
            },
            {
                "stage": 5,
                "name": "可回查證據",
                "field": "candidate_diagnoses[].internal_case_evidence[].locator / external_evidence.items[].url",
            },
            {
                "stage": 6,
                "name": "最終鑑別診斷排序",
                "field": "final_differential_diagnosis_ranking",
            },
        ],
    }


def patient_file_sort_key(path: Path) -> tuple[int, str]:
    match = re.search(r"NGS_patient_(\d+)", path.name)
    return (int(match.group(1)) if match else 10**9, path.name.casefold())


def write_simple_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        fieldnames = (
            "醫院",
            "病人編號",
            "檢驗編號",
            "選擇菌種",
            "原因",
            "角色",
            "證據來源",
            "排序",
            "選擇狀態",
        )
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in rows:
            writer.writerow(
                {
                    "醫院": item["hospital"],
                    "病人編號": item["patient_id"],
                    "檢驗編號": "; ".join(item["test_identifiers"]) or "未提供",
                    "選擇菌種": item["selected_organism"],
                    "原因": item["reason"],
                    "角色": item["picked_role"] or "未提供",
                    "證據來源": item["evidence_source"] or "未提供",
                    "排序": item["selected_rank"] or "",
                    "選擇狀態": item["selection_status"],
                }
            )


def main() -> int:
    args = parse_args()
    candidate_root = args.candidate_root.resolve()
    output_root = args.output_root.resolve()
    if not candidate_root.is_dir():
        raise SystemExit(f"candidate root is not a directory: {candidate_root}")
    if output_root.exists() and any(output_root.iterdir()):
        raise SystemExit(f"output root must be new or empty: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    simple_rows: list[dict[str, Any]] = []
    patient_result_rows: list[dict[str, Any]] = []
    patient_outputs: list[dict[str, Any]] = []
    for source_file in sorted(
        candidate_root.glob("NGS_patient_*_candidate_window_*_evidence.json"),
        key=patient_file_sort_key,
    ):
        bundle = load_json(source_file)
        if not isinstance(bundle, dict):
            continue
        selected = [
            candidate
            for candidate in bundle.get("candidates", [])
            if isinstance(candidate, dict) and candidate.get("selection_status") == "picked"
        ]
        selected.sort(
            key=lambda candidate: (
                selected_rank(candidate) is None,
                selected_rank(candidate) or 10**9,
                compact(candidate.get("name")).casefold(),
            )
        )
        selected_rows = [simple_row(bundle, candidate) for candidate in selected]
        simple_rows.extend(selected_rows)
        patient_result_rows.extend(selected_rows or [no_selection_row(bundle)])
        chain = traceable_chain(bundle, source_file)
        patient_id = chain["patient_id"]
        technical_output_file = (
            output_root
            / "technical_traces"
            / f"NGS_patient_{patient_id}_technical_trace.json"
        )
        readable_output_file = (
            output_root
            / "traceable_reasoning_chains"
            / f"NGS_patient_{patient_id}_explainable_reasoning_zh.json"
        )
        write_json(technical_output_file, chain, args.pretty)
        write_json(
            readable_output_file,
            readable_reasoning_chain(bundle, technical_output_file),
            args.pretty,
        )
        patient_outputs.append(
            {
                "patient_id": patient_id,
                "hospital": chain["hospital"],
                "candidate_count": len(chain["candidate_diagnoses"]),
                "picked_count": len(chain["final_differential_diagnosis_ranking"]),
                "readable_output_file": str(readable_output_file),
                "technical_output_file": str(technical_output_file),
            }
        )

    write_json(output_root / "human_readable_selected_pathogens.json", simple_rows, args.pretty)
    write_simple_csv(output_root / "human_readable_selected_pathogens.csv", simple_rows)
    write_json(output_root / "human_readable_patient_results.json", patient_result_rows, args.pretty)
    write_simple_csv(output_root / "human_readable_patient_results.csv", patient_result_rows)

    manifest = {
        "schema_version": "rationale_views.manifest.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_root": str(candidate_root),
        "output_root": str(output_root),
        "summary": {
            "patients": len(patient_outputs),
            "selected_pathogen_rows": len(simple_rows),
            "human_readable_patient_result_rows": len(patient_result_rows),
            "patients_with_no_picked_pathogen": sum(
                item["picked_count"] == 0 for item in patient_outputs
            ),
            "traceable_reasoning_chains": len(patient_outputs),
            "technical_traces": len(patient_outputs),
        },
        "patients": patient_outputs,
    }
    write_json(output_root / "rationale_views_manifest.json", manifest, True)
    print(
        f"patients={len(patient_outputs)} selected_rows={len(simple_rows)} "
        f"trace_chains={len(patient_outputs)} output={output_root}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

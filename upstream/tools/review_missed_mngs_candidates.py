from __future__ import annotations

import argparse
import json
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import candidate_evidence_profile as evidence_profiles  # noqa: E402
from tools import mngs_common as mngs  # noqa: E402
from tools import pathogen_normalization as pathogen_names  # noqa: E402
from tools import pathogen_rule_lists  # noqa: E402


DEFAULT_MODEL = "gpt-5"
DEFAULT_QUEUE_SUFFIX = "mNGS_missed_candidate_review_queue"
DEFAULT_MAX_SUFFIX = "mNGS_max_deterministic_resp_commensal_dominance_guardrail_opt_chosen_full"
DEFAULT_MNGS_TO_SPECIMEN_SUFFIX = "mngs_to_specimen_all_rk_ntc_opt"
DEFAULT_OUTPUT_SUFFIX = "mNGS_missed_candidate_review"


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def extract_patient_id(patient_dir: Path) -> str:
    return mngs.extract_patient_identifier(patient_dir).replace("NGS_patient_", "")


def summary_dir_for(patient_dir: Path) -> Path:
    return patient_dir / mngs.SUMMARY_OUTPUT_DIR_NAME


def path_for_suffix(
    patient_dir: Path,
    suffix: str,
    extension: str = ".json",
    artifact_root: Path | None = None,
) -> Path:
    base = mngs.extract_patient_identifier(patient_dir)
    if artifact_root is not None:
        return artifact_root / patient_dir.name / mngs.SUMMARY_OUTPUT_DIR_NAME / f"{base}_{suffix}{extension}"
    return summary_dir_for(patient_dir) / f"{base}_{suffix}{extension}"


def parse_patient_ids(values: Sequence[str] | None) -> set[str]:
    output: set[str] = set()
    for value in values or []:
        text = str(value).strip()
        if not text:
            continue
        output.add(text.removesuffix("_json") if text.startswith("NGS_patient_") else f"NGS_patient_{text}")
    return output


def compact_picked(item: dict[str, Any]) -> dict[str, Any]:
    return {
        "organism_name": item.get("organism_name"),
        "picked_role": item.get("picked_role"),
        "basis_level": item.get("basis_level"),
        "mngs_signal_tier": item.get("mngs_signal_tier"),
        "rank_priority": item.get("rank_priority"),
        "reads": item.get("reads"),
        "caution_flags": item.get("caution_flags"),
    }


def compact_candidate(item: dict[str, Any]) -> dict[str, Any]:
    return {
        "organism_name": item.get("organism_name"),
        "classification": item.get("classification"),
        "integrated_causative_level": item.get("integrated_causative_level"),
        "mngs_signal_tier": item.get("mngs_signal_tier"),
        "rank_priority": item.get("rank_priority"),
        "rank_rule": item.get("rank_rule"),
        "reads": item.get("reads"),
        "reads_tier": item.get("reads_tier"),
        "reads_percentile": item.get("reads_percentile"),
        "dominance_tier": item.get("dominance_tier"),
        "specimen_class": item.get("specimen_class"),
        "specimen_alignment": item.get("specimen_alignment"),
        "source_category": item.get("source_category"),
        "is_protected_pathogen": item.get("is_protected_pathogen"),
        "is_likely_colonizer_or_background": item.get("is_likely_colonizer_or_background"),
        "support_modules": item.get("support_modules") or item.get("module_support_summary"),
        "non_host_support_modules": item.get("non_host_support_modules"),
        "guardrail_rule": item.get("guardrail_rule") or (item.get("key_evidence") or {}).get("guardrail_rule"),
        "related_representative_hospital_support": item.get("related_representative_hospital_support"),
        "candida_invasive_hospital_support": item.get("candida_invasive_hospital_support"),
        "candida_related_invasive_hospital_context": item.get("candida_related_invasive_hospital_context"),
        "candida_evidence_strength_profile": item.get("candida_evidence_strength_profile"),
        "candida_related_evidence_strength_profile": item.get("candida_related_evidence_strength_profile"),
        "review_reasons": item.get("review_reasons"),
        "applied_rules": item.get("applied_rules"),
    }


def compact_excluded(item: dict[str, Any]) -> dict[str, Any]:
    return {
        "organism_name": item.get("organism_name"),
        "classification": item.get("classification"),
        "observed_level": item.get("observed_level"),
        "rank_priority": item.get("rank_priority"),
        "reads": item.get("reads"),
        "exclusion_reason_code": item.get("exclusion_reason_code"),
    }


def compact_hospital_evidence(item: dict[str, Any]) -> dict[str, Any]:
    return {
        "organism_name": item.get("organism_name") or item.get("name"),
        "classification": item.get("classification"),
        "evidence_modules": item.get("evidence_modules"),
        "module_level_summary": item.get("module_level_summary"),
        "module_evidence": item.get("module_evidence"),
        "best_hospital_level": item.get("best_hospital_level"),
        "reasoning": item.get("reasoning") or item.get("rationale"),
    }


def compact_summary(summary: dict[str, Any]) -> dict[str, Any]:
    module_summaries = summary.get("module_summaries")
    if not isinstance(module_summaries, dict):
        module_summaries = {}
    return {
        "rule_version": summary.get("rule_version"),
        "final_infection_likelihood": summary.get("final_infection_likelihood"),
        "dominant_source": summary.get("dominant_source"),
        "dominant_pathogen_type": summary.get("dominant_pathogen_type"),
        "support_direction": summary.get("support_direction"),
        "data_quality_tier": summary.get("data_quality_tier"),
        "module_availability": summary.get("module_availability"),
        "amr_risk_summary": summary.get("amr_risk_summary"),
        "host_context": summary.get("host_context"),
        "hospital_organism_evidence": [
            compact_hospital_evidence(item)
            for item in summary.get("hospital_organism_evidence") or []
            if isinstance(item, dict)
        ],
        "priority_flags": summary.get("priority_flags"),
        "data_gaps": summary.get("data_gaps"),
        "cross_module_reasoning": summary.get("cross_module_reasoning"),
        "module_summaries": {
            key: module_summaries.get(key)
            for key in (
                "cbc_other_lab",
                "culture",
                "filmarray_gmtest",
                "molecular_microbiology",
                "image",
            )
            if key in module_summaries
        },
    }


def compact_max_output(payload: dict[str, Any]) -> dict[str, Any]:
    best = payload.get("best_available_summary") if isinstance(payload.get("best_available_summary"), dict) else {}
    return {
        "rule_version": payload.get("rule_version"),
        "final_infection_likelihood": payload.get("final_infection_likelihood"),
        "dominant_source": payload.get("dominant_source"),
        "dominant_pathogen_type": payload.get("dominant_pathogen_type"),
        "host_context": payload.get("host_context"),
        "best_available_summary": {
            "selection_mode": best.get("selection_mode"),
            "picked_count": best.get("picked_count"),
            "picked_pathogens": [
                compact_picked(item)
                for item in best.get("picked_pathogens") or []
                if isinstance(item, dict)
            ],
        },
        "excluded_candidates": [
            compact_excluded(item)
            for item in payload.get("excluded_candidates") or []
            if isinstance(item, dict)
        ],
        "data_gaps": payload.get("data_gaps"),
    }


def compact_mngs_to_specimen(payload: dict[str, Any] | None) -> dict[str, Any] | None:
    if not isinstance(payload, dict):
        return None
    by_likelihood = payload.get("pathogens_by_likelihood")
    if not isinstance(by_likelihood, dict):
        return {
            key: payload.get(key)
            for key in ("patient_id", "specimen_summary", "data_gaps")
            if key in payload
        }
    compact: dict[str, Any] = {
        "patient_id": payload.get("patient_id"),
        "specimen_summary": payload.get("specimen_summary"),
        "pathogens_by_likelihood": {},
        "data_gaps": payload.get("data_gaps"),
    }
    for group in ("high", "medium", "low_colonizer"):
        values = by_likelihood.get(group) or []
        compact["pathogens_by_likelihood"][group] = [
            {
                "organism_name": item.get("organism_name") or item.get("name"),
                "classification": item.get("classification"),
                "specimen_id": item.get("specimen_id") or item.get("Specimen_ID"),
                "likelihood": item.get("likelihood"),
                "reasoning": item.get("reasoning") or item.get("rationale"),
            }
            for item in values
            if isinstance(item, dict)
        ]
    return compact


def build_review_payload(
    *,
    patient_dir: Path,
    queue_path: Path,
    max_path: Path,
    summary_path: Path,
    mngs_to_specimen_path: Path | None,
) -> dict[str, Any]:
    queue_payload = load_json(queue_path)
    max_payload = load_json(max_path)
    summary_payload = load_json(summary_path)
    mngs_to_specimen_payload = load_json(mngs_to_specimen_path) if mngs_to_specimen_path and mngs_to_specimen_path.exists() else None
    return {
        "patient_id": extract_patient_id(patient_dir),
        "review_purpose": "LLM safety review only; do not overwrite deterministic picked_pathogens.",
        "source_files": {
            "queue": str(queue_path),
            "deterministic_max": str(max_path),
            "deterministic_summary": str(summary_path),
            "mngs_to_specimen": str(mngs_to_specimen_path) if mngs_to_specimen_path else None,
        },
        "review_queue": {
            "version": queue_payload.get("review_queue_version"),
            "host_context": queue_payload.get("host_context"),
            "formal_picked_count": queue_payload.get("formal_picked_count"),
            "full_ranked_scored_count": queue_payload.get("full_ranked_scored_count"),
            "review_queue_count": queue_payload.get("review_queue_count"),
            "review_queue_policy": queue_payload.get("review_queue_policy"),
            "candidates": [
                compact_candidate(item)
                for item in queue_payload.get("review_queue") or []
                if isinstance(item, dict)
            ],
        },
        "deterministic_formal_output": compact_max_output(max_payload),
        "hospital_side_summary": compact_summary(summary_payload),
        "mngs_to_specimen": compact_mngs_to_specimen(mngs_to_specimen_payload),
    }


def build_prompt(payload: dict[str, Any]) -> str:
    payload_text = json.dumps(payload, ensure_ascii=False, indent=2)
    return f"""You are a clinical mNGS missed-candidate safety-review assistant for suspected pulmonary infection.
Your task is NOT to replace deterministic picked_pathogens. Your task is to review candidates that deterministic rules did not formally pick and identify organisms that still deserve clinician attention.

Input sources:
- review_queue.candidates: organisms selected for safety review, including rank, reads, reads_tier, reads_percentile, dominance, deterministic level, and guardrail reasons.
- deterministic_formal_output: organisms already formally picked by the deterministic max scorer.
- hospital_side_summary: culture, FilmArray/GM, molecular microbiology, CBC/underlying, host vulnerability, and hospital organism evidence.
- mngs_to_specimen: high / medium / low_colonizer organisms from the deterministic mNGS-to-specimen stage.

Output categories:
- review_high_priority: RAG-visible organisms with the strongest missed-review concern. Use for actionable/common respiratory pathogens, respiratory panel targets, or very strong high-burden signals where missing the organism could plausibly change pulmonary infection interpretation.
- review_context_needed: RAG-visible organisms with incomplete but clinically plausible pulmonary evidence. Use when literature, specimen context, imaging, device context, repeated culture, or host context is needed before deciding whether to retain or discard.
- review_low_specificity: audit-only organisms. Use for weak, low-specificity, colonizer/background/environmental, non-pulmonary culture-only, respiratory Candida/yeast colonization, skin flora, or weak oral flora signals that should not enter RAG by default.
- omitted_with_reason: audit-only items that should not be treated as candidate organisms, including impossible pulmonary infection sources, same-pathogen aliases already retained, same-genus duplicate representatives, or clearly non-pulmonary findings.

Important constraints:
- Do not simply upgrade because a patient is immunocompromised.
- Do not simply upgrade because an organism is a protected pathogen.
- Be cautious with respiratory Candida, oral flora, gut flora, skin flora, water/environmental organisms, and reactivation viruses.
- Respiratory Candida or generic yeast evidence alone should be review_low_specificity, not RAG-visible by default. Different Candida species must not be treated as aliases. Same-species blood, sterile-site, tissue, pleural, abscess, or histopathology evidence is strong invasive evidence and should be review_high_priority when not already formal picked, but it still does not automatically prove Candida pneumonia. BAL/lower-respiratory mNGS Candida with rank <= 2 or high percentile in a highly vulnerable host plus same-species or same-genus respiratory hospital context should be review_high_priority for RAG, while still not formal picked without same-species strong invasive evidence. High-rank/high-percentile respiratory Candida with host context but without convergent hospital evidence should be review_context_needed.
- Trichosporon and other rare opportunistic yeasts should not be collapsed into the respiratory Candida colonization rule. Same-organism lower-respiratory culture plus high mNGS burden in an immunocompromised host should be review_context_needed or higher depending on strength.
- CMV/HSV/VZV often represent reactivation or shedding; flag caution unless evidence strongly supports primary pulmonary disease. EBV should not be treated as a pulmonary infection-source candidate.
- Do not ignore very high herpesvirus reads. However, R4 very high reads with D0/D1 non-dominant pattern and host-only support should be review_context_needed rather than formal picked; R4 with D2/D3 dominance or same-organism PCR/viral-load/molecular/hospital support can remain higher concern.
- Rare mold/environmental fungi require host, GM/culture/molecular support, compatible imaging, or meaningful mNGS signal.
- Multiple oral/aspiration flora should not be listed one by one; keep only 1-2 representative organisms and describe the pattern using oral_aspiration_flora_pattern.
- Split oral/aspiration organisms into three clinical ecology groups:
  1. strict aspiration anaerobes: Bacteroides/Phocaeicola/Parabacteroides, Prevotella, Fusobacterium, Segatella, Veillonella, Porphyromonas, Parvimonas, Peptostreptococcus/Finegoldia, and Actinomyces-like aspiration anaerobes. These can be review_context_needed when lower-respiratory clean mNGS signal is strong.
  2. broad oral/upper-airway commensals: non-pathogenic Neisseria-like oral flora, Moraxella osloensis, Rothia, Gemella, Lactobacillus, and similar common mouth/nasopharynx residents. These should stay review_low_specificity unless there is direct same-organism lower-respiratory/sterile support, D2/D3 dominance, or a strong compatible syndrome.
  3. atypical oral-associated organisms: Capnocytophaga, Leptotrichia, Mycoplasma salivarium/orale, Kingella oralis, Campylobacter concisus, Abiotrophia, and similar opportunistic oral organisms. These need stronger host/specimen/cluster or hospital support than strict anaerobes before entering review_context_needed.
- Do not upgrade broad oral/upper-airway commensals or atypical oral-associated organisms solely because rank, reads, percentile, or host vulnerability is high. High percentile with mNGS-only, D0/D1, and no same-organism hospital support should remain review_low_specificity.
- Do not demote a true oral/aspiration pattern solely because deterministic output already has another plausible primary pathogen. If multiple oral/anaerobe organisms appear together with compatible aspiration/abscess/empyema/necrotizing pneumonia context, keep 1-2 representative organisms in review_context_needed and mark them as non-primary.

Pulmonary-specific safety rules:
- Lower-respiratory specimen + single strict aspiration anaerobe + clean mNGS signal + rank <= 3 + (reads_tier >= R3 or reads_percentile >= 0.85): keep as review_context_needed, not picked. This is a RAG-visible safety-net signal because strict anaerobes can matter in aspiration, abscess, empyema, or necrotizing pneumonia, even when they are D0/D1 and lack same-organism hospital support.
- Promote strict anaerobes or oral-associated organisms to review_high_priority only with unusually strong convergence: same-organism lower-respiratory culture/molecular support, pus/pleural/abscess/sterile-site support, D3 dominance with high burden, or explicit abscess/empyema/necrotizing pneumonia context. Otherwise keep high-signal single strict anaerobes in review_context_needed.
- Water/environmental GNB with mNGS-only, no dominance, and no same-organism hospital support should not become formal picked.
- Water/environmental GNB with clean mNGS signal and high reads or high rank should usually be review_low_specificity. High host vulnerability alone is not enough to make these RAG-visible.
- Water/environmental GNB should be review_context_needed only with blood culture, effective BALF/respiratory culture, repeated same organism, D2/D3 dominance, sterile/systemic evidence, or a clear device/source context. mNGS-only high reads/percentile with D0/D1 and no same-organism support should remain review_low_specificity.
- GI-only pathogens without pulmonary/systemic support should generally be omitted_with_reason.
- Corynebacterium/Cutibacterium/non-aureus Staphylococcus are high colonization/contamination-risk organisms. Do not treat percentile alone as enough for formal picked. Lower-respiratory clean NTC/RK_NTC + rank <= 2 + reads_tier >= R3 + percentile >= 0.95 should be review_low_specificity or review_context_needed; move to review_high_priority only if there is D3 dominance with very high reads, or another unusually strong reason to prioritize RAG.

Return only valid JSON. Do not use markdown or code fences.
All fields ending with _zh, plus risk_if_ignored, why_not_auto_upgrade, suggested_next_check, and interpretation_zh, must be written in Traditional Chinese.

Required JSON schema:
{{
  "patient_id": "string",
  "review_purpose": "safety_review_only",
  "overall_assessment_zh": "string",
  "review_high_priority": [
    {{
      "organism_name": "string",
      "confidence": "low|moderate|high",
      "rationale_zh": "string",
      "risk_if_ignored": "string",
      "why_not_auto_upgrade": "string",
      "suggested_next_check": "string"
    }}
  ],
  "review_context_needed": [
    {{
      "organism_name": "string",
      "confidence": "low|moderate|high",
      "rationale_zh": "string",
      "risk_if_ignored": "string",
      "why_not_auto_upgrade": "string",
      "suggested_next_check": "string"
    }}
  ],
  "review_low_specificity": [
    {{
      "organism_name": "string",
      "confidence": "low|moderate|high",
      "rationale_zh": "string",
      "risk_if_ignored": "string",
      "why_not_auto_upgrade": "string",
      "suggested_next_check": "string"
    }}
  ],
  "omitted_with_reason": [
    {{
      "organism_name": "string",
      "rationale_zh": "string"
    }}
  ],
  "oral_aspiration_flora_pattern": {{
    "enabled": true,
    "pattern_detected": true,
    "representative_organisms": ["string"],
    "grouped_non_representative_organisms": ["string"],
    "interpretation_zh": "string"
  }},
  "llm_limitations": [
    "This LLM review is advisory only and does not override deterministic picked_pathogens.",
    "Clinician review is required before changing the final interpretation."
  ]
}}

Input payload:
{payload_text}
"""


def extract_output_text(response: Any) -> str:
    if getattr(response, "output_text", None):
        text = response.output_text.strip()
        if text:
            return text
    segments: list[str] = []
    for item in getattr(response, "output", []) or []:
        for content in getattr(item, "content", []) or []:
            text = getattr(content, "text", None)
            if text:
                segments.append(str(text))
    return "".join(segments).strip()


def model_accepts_temperature(model: str) -> bool:
    normalized = (model or "").strip().lower()
    return not normalized.startswith("gpt-5")


def send_to_llm(
    prompt: str,
    *,
    model: str,
    temperature: float | None,
    reasoning_effort: str | None = None,
) -> tuple[str, dict[str, int | None] | None]:
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY environment variable is not set.")
    from openai import OpenAI  # type: ignore

    client = OpenAI(api_key=api_key)
    kwargs: dict[str, Any] = {"model": model, "input": prompt}
    if temperature is not None and model_accepts_temperature(model):
        kwargs["temperature"] = temperature
    if reasoning_effort:
        kwargs["reasoning"] = {"effort": reasoning_effort}
    response = client.responses.create(**kwargs)
    text = extract_output_text(response)
    if not text:
        raise RuntimeError("LLM response did not contain text.")
    usage = getattr(response, "usage", None)
    usage_summary = None
    if usage is not None:
        usage_summary = {
            "prompt_tokens": getattr(usage, "input_tokens", None),
            "completion_tokens": getattr(usage, "output_tokens", None),
            "total_tokens": getattr(usage, "total_tokens", None),
        }
    return text, usage_summary


def close_unterminated_json_containers(text: str) -> str | None:
    """Close trailing JSON containers without changing response content."""
    stack: list[str] = []
    in_string = False
    escaped = False
    matching = {"}": "{", "]": "["}

    for char in text:
        if in_string:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                in_string = False
            continue

        if char == '"':
            in_string = True
        elif char in "[{":
            stack.append(char)
        elif char in "]}":
            if not stack or stack[-1] != matching[char]:
                return None
            stack.pop()

    if in_string or not stack:
        return None
    closers = {"{": "}", "[": "]"}
    return text + "".join(closers[opening] for opening in reversed(stack))


def remove_trailing_json_commas(text: str) -> str:
    """Remove commas immediately before JSON container closers, outside strings."""
    output: list[str] = []
    in_string = False
    escaped = False
    index = 0

    while index < len(text):
        char = text[index]
        if in_string:
            output.append(char)
            if escaped:
                escaped = False
            elif ord(char) == 92:
                escaped = True
            elif char == '"':
                in_string = False
            index += 1
            continue

        if char == '"':
            in_string = True
            output.append(char)
            index += 1
            continue

        if char == ",":
            next_index = index + 1
            while next_index < len(text) and text[next_index].isspace():
                next_index += 1
            if next_index < len(text) and text[next_index] in "]}":
                index += 1
                continue

        output.append(char)
        index += 1

    return "".join(output)


def parse_json_response(text: str) -> Any:
    stripped = text.strip()
    if stripped.startswith("```"):
        stripped = re.sub(r"^```(?:json)?\s*", "", stripped, flags=re.IGNORECASE)
        stripped = re.sub(r"\s*```$", "", stripped)
    candidates = [stripped]
    match = re.search(r"\{.*\}", stripped, flags=re.DOTALL)
    if match and match.group(0) != stripped:
        candidates.append(match.group(0))

    last_error: json.JSONDecodeError | None = None
    for candidate in candidates:
        attempts = [candidate]
        closed = close_unterminated_json_containers(candidate)
        if closed is not None:
            attempts.append(closed)
        for attempt in tuple(attempts):
            without_trailing_commas = remove_trailing_json_commas(attempt)
            if without_trailing_commas != attempt:
                attempts.append(without_trailing_commas)

        for attempt in attempts:
            try:
                return json.loads(attempt)
            except json.JSONDecodeError as exc:
                last_error = exc

    if last_error is not None:
        raise last_error
    raise json.JSONDecodeError("No JSON object found", stripped, 0)


GENUS_SKIP_TOKENS = {
    "cmv",
    "ebv",
    "hsv",
    "epec",
    "esbl",
    "covid",
    "coronavirus",
    "human",
    "rhino",
    "enterovirus",
    "rhinovirus",
    "severe",
    "sars",
    "aspergillosis",
}


def normalized_name(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def genus_key(value: Any) -> str:
    text = str(value or "").strip()
    if not text or "/" in text:
        return ""
    text = re.sub(r"\([^)]*\)", " ", text)
    tokens = [token for token in re.split(r"\s+", text) if token]
    if len(tokens) < 2:
        return ""
    first = re.sub(r"[^A-Za-z]+", "", tokens[0]).lower()
    second = re.sub(r"[^A-Za-z]+", "", tokens[1]).lower()
    if not first or not second or first in GENUS_SKIP_TOKENS:
        return ""
    if first.endswith("virus") or second in {"virus", "viruses", "dna", "rna"}:
        return ""
    return first


def numeric_value(value: Any, default: float = 0.0) -> float:
    try:
        if value in (None, ""):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def level_score(value: Any) -> float:
    text = str(value or "")
    match = re.search(r"(?:Level|L)([1-5])", text, flags=re.IGNORECASE)
    if not match:
        return 0.0
    return {1: 500.0, 2: 400.0, 3: 250.0, 4: 80.0, 5: 0.0}.get(int(match.group(1)), 0.0)


def tier_score(value: Any) -> float:
    text = str(value or "")
    if "M1" in text or "R4" in text:
        return 220.0
    if "M2" in text or "R3" in text:
        return 170.0
    if "M3" in text or "R2" in text:
        return 110.0
    if "M4" in text or "R1" in text:
        return 45.0
    return 0.0


def dominance_score(value: Any) -> float:
    text = str(value or "")
    if "D3" in text:
        return 150.0
    if "D2" in text:
        return 100.0
    if "D1" in text:
        return 35.0
    return 0.0


def confidence_score(value: Any) -> float:
    text = str(value or "").strip().lower()
    return {"high": 300.0, "moderate": 180.0, "low": 60.0}.get(text, 0.0)


def support_score(value: Any) -> float:
    if not isinstance(value, list):
        return 0.0
    score = 0.0
    for item in value:
        text = str(item or "").lower()
        if any(key in text for key in ("culture", "filmarray", "molecular", "gm", "serology")):
            score += 90.0
        elif text and text != "host":
            score += 25.0
    return score


def tier_rank(value: Any) -> int:
    text = str(value or "")
    if "R4" in text:
        return 4
    if "R3" in text:
        return 3
    if "R2" in text:
        return 2
    if "R1" in text:
        return 1
    return 0


def build_queue_evidence_index(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    queue = payload.get("review_queue") or []
    if isinstance(queue, dict):
        candidates: list[Any] = list(queue.get("candidates") or [])
    elif isinstance(queue, list):
        candidates = list(queue)
    else:
        candidates = []

    for key in (
        "supplemental_full_ranked_review_candidates",
        "hospital_detected_review_candidates",
        "hospital_detected_merged_candidates",
    ):
        values = payload.get(key) or []
        if isinstance(values, list):
            candidates.extend(values)

    index: dict[str, dict[str, Any]] = {}
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        for key in {
            normalized_name(candidate.get("organism_name")),
            pathogen_names.canonical_key(candidate.get("organism_name")),
            *pathogen_names.expanded_name_keys(candidate.get("organism_name")),
        }:
            if key:
                index[key] = candidate
    return index


LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_PREFIXES = (
    "arthrobacter",
    "agrobacterium",
    "brachybacterium",
    "janibacter",
    "kocuria",
    "microbacterium",
    "micrococcus",
    "methylobacterium",
    "novosphingobium",
    "paracoccus",
    "polynucleobacter",
    "sphingobium",
    "sphingomonas",
    "undibacterium",
)

WATER_ENVIRONMENTAL_GNB_PREFIXES = (
    "achromobacter",
    "acidovorax",
    "brevundimonas",
    "chryseobacterium",
    "cloacibacterium",
    "comamonas",
    "delftia",
    "elizabethkingia",
    "flavobacterium",
    "paraburkholderia",
    "phytobacter",
    "ralstonia",
    "pseudoxanthomonas",
    "xanthomonas",
)

LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_EXACT_NAMES = {
    "bacillusmycoides",
}

WATER_ENVIRONMENTAL_GNB_EXACT_NAMES = {
    "acinetobacterjohnsonii",
    "brucellaintermedia",
    "enterobactersoli",
    "pseudomonasalcaligenes",
    "pseudomonasceruminis",
    "pseudomonasjaponica",
    "pseudomonaskhazarica",
    "pseudomonaslactis",
    "pseudomonassoli",
    "pseudomonaszhaodongensis",
}

ORAL_ASPIRATION_FLORA_PREFIXES = (
    "prevotella",
    "fusobacterium",
    "bacteroides",
    "phocaeicola",
    "segatella",
    "parabacteroides",
    "veillonella",
    "parvimonas",
    "porphyromonas",
    "peptostreptococcus",
    "finegoldia",
    "actinomyces",
    "capnocytophaga",
    "leptotrichia",
    "gemella",
    "rothia",
    "lactobacillus",
    "mycoplasmasalivarium",
    "mycoplasmaorale",
    "abiotrophia",
)

STRICT_ASPIRATION_ANAEROBE_PREFIXES = (
    "prevotella",
    "fusobacterium",
    "bacteroides",
    "phocaeicola",
    "segatella",
    "parabacteroides",
    "veillonella",
    "parvimonas",
    "porphyromonas",
    "peptostreptococcus",
    "finegoldia",
    "actinomyces",
)

BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_PREFIXES = (
    "gemella",
    "lacticaseibacillus",
    "ligilactobacillus",
    "limosilactobacillus",
    "rothia",
    "lactobacillus",
)

ATYPICAL_ORAL_ASSOCIATED_PREFIXES = (
    "capnocytophaga",
    "leptotrichia",
    "mycoplasmasalivarium",
    "mycoplasmaorale",
    "abiotrophia",
)

ORAL_ASPIRATION_FLORA_EXACT_NAMES = {
    "kingellaoralis",
    "campylobacterconcisus",
    "moraxellaosloensis",
    "mycoplasmasalivarium",
    "neisseriasubflava",
    "neisseriaflavescens",
    "neisseriamucosa",
    "neisseriasicca",
    "neisseriacinerea",
    "neisseriaperflava",
    "neisseriaelongata",
    "streptococcustaonis",
}

STRICT_ASPIRATION_ANAEROBE_EXACT_NAMES: set[str] = set()

BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_EXACT_NAMES = {
    "moraxellaosloensis",
    "neisseriasubflava",
    "neisseriaflavescens",
    "neisseriamucosa",
    "neisseriasicca",
    "neisseriacinerea",
    "neisseriaperflava",
    "neisseriaelongata",
    "streptococcustaonis",
}

ATYPICAL_ORAL_ASSOCIATED_EXACT_NAMES = {
    "kingellaoralis",
    "campylobacterconcisus",
    "mycoplasmasalivarium",
}

ASPIRATION_PATTERN_CONTEXT_TERMS = (
    "aspiration",
    "aspirated",
    "aspirating",
    "aspiration pneumonia",
    "lung abscess",
    "abscess",
    "empyema",
    "necrotizing pneumonia",
    "necrotising pneumonia",
    "膿胸",
    "肺膿瘍",
    "吸入",
    "吸入性",
    "壞死性肺炎",
)

GI_PATHOGEN_EXACT_NAMES = {
    "clostridioidesdifficile",
    "clostridiumdifficile",
    "epec",
    "enteropathogenicescherichiacoli",
    "enteropathogenicecoli",
    "etec",
    "ehec",
    "stec",
    "eaec",
    "shigella",
    "norovirus",
    "rotavirus",
    "sapovirus",
    "astrovirus",
}

GI_CONTEXT_TERMS = (
    "stool",
    "feces",
    "faeces",
    "rectal",
    "gastrointestinal",
    "gastroenteritis",
    "diarrhea",
    "diarrhoea",
    "colon",
    "intestine",
)

PULMONARY_OR_SYSTEMIC_CONTEXT_TERMS = (
    "blood",
    "balf",
    "bal",
    "sputum",
    "endotracheal",
    "eta",
    "respiratory",
    "lung",
    "pleural",
    "tissue",
    "csf",
    "abscess",
    "pus",
)

HIGH_BURDEN_CONTAMINANT_RISK_PREFIXES = (
    "corynebacterium",
    "cutibacterium",
)


def starts_with_any_normalized(name: Any, prefixes: tuple[str, ...]) -> bool:
    normalized = normalized_name(name)
    return any(normalized.startswith(prefix) for prefix in prefixes)


def is_non_aureus_staphylococcus_name(name: Any) -> bool:
    normalized = normalized_name(name)
    return normalized.startswith("staphylococcus") and not normalized.startswith("staphylococcusaureus")


def is_high_burden_contaminant_risk_name(name: Any) -> bool:
    return is_non_aureus_staphylococcus_name(name) or starts_with_any_normalized(
        name, HIGH_BURDEN_CONTAMINANT_RISK_PREFIXES
    )


def is_water_or_low_specificity_environmental_name(name: Any) -> bool:
    normalized = normalized_name(name)
    return (
        normalized in LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_EXACT_NAMES
        or
        normalized in WATER_ENVIRONMENTAL_GNB_EXACT_NAMES
        or starts_with_any_normalized(name, WATER_ENVIRONMENTAL_GNB_PREFIXES)
        or starts_with_any_normalized(name, LOW_PULMONARY_SPECIFICITY_ENVIRONMENTAL_PREFIXES)
    )


def is_oral_aspiration_flora_name(name: Any) -> bool:
    normalized = normalized_name(name)
    return normalized in ORAL_ASPIRATION_FLORA_EXACT_NAMES or starts_with_any_normalized(
        name, ORAL_ASPIRATION_FLORA_PREFIXES
    )


def oral_aspiration_flora_category(name: Any) -> str:
    normalized = normalized_name(name)
    if normalized in STRICT_ASPIRATION_ANAEROBE_EXACT_NAMES or starts_with_any_normalized(
        name, STRICT_ASPIRATION_ANAEROBE_PREFIXES
    ):
        return "strict_aspiration_anaerobe"
    if normalized in BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_EXACT_NAMES or starts_with_any_normalized(
        name, BROAD_ORAL_UPPER_AIRWAY_COMMENSAL_PREFIXES
    ):
        return "broad_oral_upper_airway_commensal"
    if normalized in ATYPICAL_ORAL_ASSOCIATED_EXACT_NAMES or starts_with_any_normalized(
        name, ATYPICAL_ORAL_ASSOCIATED_PREFIXES
    ):
        return "atypical_oral_associated"
    if is_oral_aspiration_flora_name(name):
        return "other_oral_aspiration_flora"
    return "not_oral_aspiration_flora"


def is_strict_aspiration_anaerobe_name(name: Any) -> bool:
    return oral_aspiration_flora_category(name) == "strict_aspiration_anaerobe"


def recursive_text(value: Any) -> str:
    if isinstance(value, dict):
        return " ".join(recursive_text(item) for item in value.values())
    if isinstance(value, list):
        return " ".join(recursive_text(item) for item in value)
    return str(value or "")


def has_non_host_support(evidence: dict[str, Any]) -> bool:
    modules = evidence.get("non_host_support_modules")
    if not isinstance(modules, list):
        modules = evidence.get("support_modules")
    if not isinstance(modules, list):
        return False
    for module in modules:
        text = str(module or "").strip().lower()
        if text and text not in {"host", "legacy_final_summary"}:
            return True
    return False


def has_related_representative_context(evidence: dict[str, Any]) -> bool:
    related = evidence.get("related_representative_hospital_support")
    if isinstance(related, list) and related:
        return True
    key_evidence = evidence.get("key_evidence")
    if isinstance(key_evidence, dict):
        related = key_evidence.get("related_representative_hospital_support")
        return isinstance(related, list) and bool(related)
    return False


def evidence_rule_ids(evidence: dict[str, Any]) -> set[str]:
    rule_ids: set[str] = set()
    applied_rules = evidence.get("applied_rules")
    if isinstance(applied_rules, list):
        rule_ids.update(str(rule) for rule in applied_rules if rule)
    guardrail_rule = evidence.get("guardrail_rule")
    if guardrail_rule:
        rule_ids.add(str(guardrail_rule))
    key_evidence = evidence.get("key_evidence")
    if isinstance(key_evidence, dict) and key_evidence.get("guardrail_rule"):
        rule_ids.add(str(key_evidence.get("guardrail_rule")))
    formal_rule = evidence.get("formal_pick_exclusion_rule")
    if isinstance(formal_rule, dict) and formal_rule.get("rule_id"):
        rule_ids.add(str(formal_rule.get("rule_id")))
    review_reasons = evidence.get("review_reasons")
    if isinstance(review_reasons, list):
        rule_ids.update(str(reason) for reason in review_reasons if reason)
    return rule_ids


def is_reactivation_herpesvirus_name(name: Any) -> bool:
    return canonical_name_key(name) in REACTIVATION_HERPESVIRUS_KEYS


def herpes_r4_nondominant_host_only_review(name: Any, evidence: dict[str, Any]) -> bool:
    if HERPES_R4_NONDOMINANT_HOST_ONLY_REVIEW_RULE in evidence_rule_ids(evidence):
        return True
    return (
        is_reactivation_herpesvirus_name(name)
        and tier_rank(evidence.get("reads_tier")) >= 4
        and str(evidence.get("dominance_tier") or "") in {"D0_not_top", "D1_low"}
        and not has_non_host_support(evidence)
    )


def herpes_ranked_host_context_review(name: Any, evidence: dict[str, Any], payload: dict[str, Any]) -> bool:
    if HERPES_HOST_R2R3_NONDOMINANT_REVIEW_RULE in evidence_rule_ids(evidence):
        return True
    return (
        is_reactivation_herpesvirus_name(name)
        and host_is_highly_vulnerable(payload)
        and numeric_value(evidence.get("rank_priority"), default=99.0) <= 2
        and tier_rank(evidence.get("reads_tier")) in {2, 3}
        and numeric_value(evidence.get("reads_percentile")) >= 0.85
        and str(evidence.get("dominance_tier") or "") in {"D0_not_top", "D1_low"}
        and not has_non_host_support(evidence)
    )


def is_candida_or_generic_yeast_name(name: Any) -> bool:
    return pathogen_names.is_candida_or_generic_yeast(name)


def is_rare_opportunistic_yeast_name(name: Any) -> bool:
    return direct_genus_key(name) in {"trichosporon"}


def candida_has_invasive_context(evidence: dict[str, Any]) -> bool:
    return evidence_profiles.candida_has_invasive_context(evidence)


def candida_evidence_profile(evidence: dict[str, Any]) -> dict[str, Any]:
    return evidence_profiles.candida_evidence_profile(evidence)


def candida_best_evidence_strength(evidence: dict[str, Any]) -> str:
    return evidence_profiles.candida_best_evidence_strength(evidence)


def candida_has_intermediate_context(evidence: dict[str, Any]) -> bool:
    return evidence_profiles.candida_has_intermediate_context(evidence)


def evidence_row_is_lower_respiratory(row: Any) -> bool:
    if not isinstance(row, dict):
        return False
    category = canonical_name_key(row.get("specimen_category"))
    if category in {"lowerrespiratory", "lower_respiratory"}:
        return True
    specimen = recursive_text(
        {
            "specimen_type": row.get("specimen_type"),
            "specimen_category": row.get("specimen_category"),
            "sample": row.get("sample"),
            "sample_type": row.get("sample_type"),
            "source": row.get("source"),
            "site": row.get("site"),
            "body_site": row.get("body_site"),
        }
    ).lower()
    return any(
        term in specimen
        for term in ("bal", "balf", "bronchoalveolar", "lower respiratory", "sputum", "eta", "endotracheal", "tracheal")
    )


def candida_hospital_evidence_has_lower_respiratory_context(item: dict[str, Any]) -> bool:
    module_evidence = item.get("module_evidence")
    if not isinstance(module_evidence, dict):
        module_evidence = item.get("hospital_module_evidence")
    if not isinstance(module_evidence, dict):
        return False
    for rows in module_evidence.values():
        if not isinstance(rows, list):
            continue
        if any(evidence_row_is_lower_respiratory(row) for row in rows):
            return True
    return False


def candida_related_respiratory_hospital_context(name: Any, evidence: dict[str, Any], payload: dict[str, Any]) -> bool:
    if not is_candida_or_generic_yeast_name(name):
        return False
    related = evidence.get("related_representative_hospital_support")
    if isinstance(related, list) and related:
        return True
    evidence_hospital = evidence.get("hospital_evidence")
    if isinstance(evidence_hospital, list) and any(
        is_candida_or_generic_yeast_name(item.get("organism_name") or item.get("name"))
        and candida_hospital_evidence_has_lower_respiratory_context(item)
        for item in evidence_hospital
        if isinstance(item, dict)
    ):
        return True
    target_group = pathogen_names.representative_group_key(name)
    summary = payload.get("hospital_side_summary")
    if not isinstance(summary, dict) or not target_group:
        return False
    for item in summary.get("hospital_organism_evidence") or []:
        if not isinstance(item, dict):
            continue
        item_name = item.get("organism_name") or item.get("name")
        if pathogen_names.representative_group_key(item_name) != target_group:
            continue
        if candida_hospital_evidence_has_lower_respiratory_context(item):
            return True
    return False


def candida_respiratory_high_priority_review(name: Any, evidence: dict[str, Any], payload: dict[str, Any]) -> bool:
    if not is_candida_or_generic_yeast_name(name):
        return False
    if str(evidence.get("specimen_class") or "") != "S2_lower_respiratory":
        return False
    if not host_is_highly_vulnerable(payload):
        return False
    rank = numeric_value(evidence.get("rank_priority"), default=99.0)
    reads_rank = tier_rank(evidence.get("reads_tier"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    has_mngs_signal = rank <= 2 and (reads_rank >= 2 or percentile >= 0.90)
    has_hospital_context = (
        has_non_host_support(evidence)
        or candida_has_invasive_context(evidence)
        or bool(evidence.get("candida_related_invasive_hospital_context"))
        or candida_related_respiratory_hospital_context(name, evidence, payload)
    )
    return has_mngs_signal and has_hospital_context


def candida_respiratory_context_review(name: Any, evidence: dict[str, Any], payload: dict[str, Any]) -> bool:
    if not is_candida_or_generic_yeast_name(name):
        return False
    if candida_has_invasive_context(evidence):
        return True
    if str(evidence.get("specimen_class") or "") != "S2_lower_respiratory":
        return False
    rank = numeric_value(evidence.get("rank_priority"), default=99.0)
    reads_rank = tier_rank(evidence.get("reads_tier"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    has_contextual_support = (
        host_is_highly_vulnerable(payload)
        or has_non_host_support(evidence)
        or bool(evidence.get("candida_related_invasive_hospital_context"))
        or candida_related_respiratory_hospital_context(name, evidence, payload)
    )
    has_mngs_signal = rank <= 2 and (reads_rank >= 2 or percentile >= 0.85)
    return has_contextual_support and has_mngs_signal


def weak_l3_host_only_nondominant_context(evidence: dict[str, Any]) -> bool:
    rule_ids = evidence_rule_ids(evidence)
    return (
        (
            "R-S4-L3-MNGS-ONLY-HOST-LOW-NONDOMINANT" in rule_ids
            or "borderline_level_3_not_in_formal_pick" in rule_ids
        )
        and tier_rank(evidence.get("reads_tier")) <= 1
        and not has_non_host_support(evidence)
        and str(evidence.get("dominance_tier") or "") in {"D0_not_top", "D1_low"}
    )


def weak_trace_protected_context(evidence: dict[str, Any]) -> bool:
    rule_ids = evidence_rule_ids(evidence)
    return (
        (
            "R-S4-TRACE-PROTECTED-WATCHLIST-PREFERRED" in rule_ids
            or "protected_pathogen" in rule_ids
        )
        and tier_rank(evidence.get("reads_tier")) <= 1
        and not has_non_host_support(evidence)
        and str(evidence.get("dominance_tier") or "") in {"D0_not_top", "D1_low"}
    )


def weak_high_risk_mold_context(evidence: dict[str, Any]) -> bool:
    return (
        "full_ranked_high_risk_mold_or_dimorphic_fungus" in evidence_rule_ids(evidence)
        and tier_rank(evidence.get("reads_tier")) <= 2
        and not has_non_host_support(evidence)
        and str(evidence.get("dominance_tier") or "") in {"D0_not_top", "D1_low"}
    )


def has_strong_environmental_signal(evidence: dict[str, Any]) -> bool:
    rank = numeric_value(evidence.get("rank_priority"), default=99.0)
    reads_rank = tier_rank(evidence.get("reads_tier"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    dominance = str(evidence.get("dominance_tier") or "")
    specimen_cls = evidence.get("specimen_class")
    if has_non_host_support(evidence) or dominance in {"D2_moderate", "D3_dominant"}:
        return True
    if specimen_cls == "S1_sterile_systemic" and rank <= 2 and (reads_rank >= 2 or percentile >= 0.90):
        return True
    return False


def has_minimal_environmental_watch_signal(evidence: dict[str, Any]) -> bool:
    rank = numeric_value(evidence.get("rank_priority"), default=99.0)
    reads = numeric_value(evidence.get("reads"))
    reads_rank = tier_rank(evidence.get("reads_tier"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    specimen_cls = evidence.get("specimen_class")
    if has_non_host_support(evidence):
        return True
    if reads <= 0:
        return False
    if specimen_cls == "S1_sterile_systemic" and rank <= 3 and reads_rank >= 2:
        return True
    if specimen_cls == "S2_lower_respiratory" and rank <= 3 and (reads_rank >= 2 or percentile >= 0.85):
        return True
    return False


def has_strong_aspiration_anaerobe_watch_signal(evidence: dict[str, Any], name: Any = None) -> bool:
    name = name if name is not None else evidence.get("organism_name")
    if not is_strict_aspiration_anaerobe_name(name):
        return False
    rank = numeric_value(evidence.get("rank_priority"), default=99.0)
    reads = numeric_value(evidence.get("reads"))
    reads_rank = tier_rank(evidence.get("reads_tier"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    specimen_cls = evidence.get("specimen_class")
    if reads <= 0 or specimen_cls != "S2_lower_respiratory" or rank > 3:
        return False
    return reads_rank >= 3 or percentile >= 0.85 or has_non_host_support(evidence)


def has_aspiration_clinical_context(payload: dict[str, Any]) -> bool:
    text = recursive_text(
        {
            "hospital_side_summary": payload.get("hospital_side_summary"),
            "mngs_to_specimen": payload.get("mngs_to_specimen"),
        }
    ).lower()
    return any(term in text for term in ASPIRATION_PATTERN_CONTEXT_TERMS)


def oral_aspiration_cluster_count(payload: dict[str, Any]) -> int:
    review_queue = payload.get("review_queue")
    if isinstance(review_queue, dict):
        candidates = review_queue.get("candidates") or []
    elif isinstance(review_queue, list):
        candidates = review_queue
    else:
        candidates = []
    seen: set[str] = set()
    for candidate in candidates:
        if not isinstance(candidate, dict):
            continue
        name = candidate.get("organism_name")
        if not is_oral_aspiration_flora_name(name):
            continue
        if str(candidate.get("specimen_class") or "") != "S2_lower_respiratory":
            continue
        if numeric_value(candidate.get("reads")) <= 0:
            continue
        if numeric_value(candidate.get("rank_priority"), default=99.0) > 3:
            continue
        key = canonical_name_key(name)
        if key:
            seen.add(key)
    return len(seen)


def has_oral_aspiration_context_signal(name: Any, evidence: dict[str, Any], payload: dict[str, Any]) -> bool:
    if not is_oral_aspiration_flora_name(name):
        return False
    category = oral_aspiration_flora_category(name)
    specimen_cls = evidence.get("specimen_class")
    dominance = str(evidence.get("dominance_tier") or "")
    reads = numeric_value(evidence.get("reads"))
    if specimen_cls == "S1_sterile_systemic" and reads > 0:
        return True
    if specimen_cls != "S2_lower_respiratory" or reads <= 0:
        return False
    if has_non_host_support(evidence) or dominance in {"D2_moderate", "D3_dominant"}:
        return True
    if category == "strict_aspiration_anaerobe" and has_strong_aspiration_anaerobe_watch_signal(evidence, name):
        return True
    if category == "atypical_oral_associated":
        rank = numeric_value(evidence.get("rank_priority"), default=99.0)
        reads_rank = tier_rank(evidence.get("reads_tier"))
        percentile = numeric_value(evidence.get("reads_percentile"))
        has_minimal_signal = rank <= 3 and (reads_rank >= 2 or percentile >= 0.85)
        return (
            has_minimal_signal
            and has_aspiration_clinical_context(payload)
            and oral_aspiration_cluster_count(payload) >= 2
        )
    return False


def has_lower_respiratory_hospital_support(evidence: dict[str, Any]) -> bool:
    module_evidence = evidence.get("hospital_module_evidence")
    if not isinstance(module_evidence, dict):
        module_evidence = evidence.get("module_evidence")
    if not isinstance(module_evidence, dict):
        return False
    for rows in module_evidence.values():
        if not isinstance(rows, list):
            continue
        if any(evidence_row_is_lower_respiratory(row) for row in rows):
            return True
    return False


def has_oral_aspiration_high_priority_signal(name: Any, evidence: dict[str, Any], payload: dict[str, Any]) -> bool:
    if not is_oral_aspiration_flora_name(name):
        return False
    category = oral_aspiration_flora_category(name)
    reads = numeric_value(evidence.get("reads"))
    if reads <= 0:
        return False
    specimen_cls = evidence.get("specimen_class")
    dominance = str(evidence.get("dominance_tier") or "")
    reads_rank = tier_rank(evidence.get("reads_tier"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    if specimen_cls == "S1_sterile_systemic":
        return True
    if specimen_cls != "S2_lower_respiratory":
        return False
    strong_signal = category == "strict_aspiration_anaerobe" and has_strong_aspiration_anaerobe_watch_signal(evidence, name)
    if not strong_signal:
        return False
    if has_lower_respiratory_hospital_support(evidence):
        return True
    if dominance == "D3_dominant" and (reads_rank >= 3 or percentile >= 0.85):
        return True
    return False


DIRECT_REVIEW_TIERS = (
    "review_high_priority",
    "review_context_needed",
    "review_low_specificity",
    "omitted_with_reason",
)

DIRECT_ACTIVE_TIERS = (
    "review_high_priority",
    "review_context_needed",
    "review_low_specificity",
)

LEGACY_REVIEW_SECTIONS = (
    "possible_missed_pathogens",
    "watchlist_candidates",
    "not_recommended_for_upgrade",
)

ACTIONABLE_RESPIRATORY_REVIEW_KEYS = {
    "haemophilusinfluenzae",
    "enterobactercloacaecomplex",
}

CONTEXT_REVIEW_KEYS = {
    "acinetobacterbaumannii",
    "cronobactersakazakii",
    "enterobactersoli",
    "escherichiacoli",
    "mycobacteriumtuberculosis",
    "nocardiathailandica",
    "rhodococcusqingshengii",
    "schizophyllumcommune",
    "severeacuterespiratorysyndromerelatedcoronavirus",
    "stenotrophomonasmaltophilia",
    "streptococcuspneumoniae",
    "talaromycespinophilus",
}

REACTIVATION_HERPESVIRUS_KEYS = {
    "humanbetaherpesvirus5",  # CMV
    "humanalphaherpesvirus1",  # HSV-1
    "humanalphaherpesvirus2",  # HSV-2
    "humanalphaherpesvirus3",  # VZV
}

HERPES_R4_NONDOMINANT_HOST_ONLY_REVIEW_RULE = "R-S4-HERPES-R4-NONDOMINANT-HOST-ONLY-REVIEW"
HERPES_HOST_R2R3_NONDOMINANT_REVIEW_RULE = "R-S4-HERPES-HOST-R2R3-NONDOMINANT-REVIEW"
NON_CORE_RESPIRATORY_VIRUS_FORMAL_PICK_GUARDRAIL = "R-S4-NON-CORE-RESP-VIRUS-PICKED-OMIT"
MNGS_ONLY_NONDOMINANT_NON_TYPICAL_FORMAL_PICK_GUARDRAIL = "R-S4-MNGS-ONLY-NONDOMINANT-NON-TYPICAL-PICKED-OMIT"
COMMON_HOSPITAL_LOW_READ_NODOM_FORMAL_PICK_GUARDRAIL = "R-S4-COMMON-HOSPITAL-R0R1-NONDOMINANT-PICKED-OMIT"
CANDIDA_EVIDENCE_STRONG_INVASIVE = evidence_profiles.CANDIDA_EVIDENCE_STRONG_INVASIVE
CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC = evidence_profiles.CANDIDA_EVIDENCE_INTERMEDIATE_SYSTEMIC
CANDIDA_EVIDENCE_WEAK_NONINVASIVE = evidence_profiles.CANDIDA_EVIDENCE_WEAK_NONINVASIVE

ASPIRATION_CONTEXT_KEYS = {
    "bacteroidesfragilis",
    "fusobacteriumnucleatum",
    "porphyromonasendodontalis",
}

TIER_PRIORITY = {
    "review_high_priority": 4,
    "review_context_needed": 3,
    "review_low_specificity": 2,
    "omitted_with_reason": 1,
}


def canonical_name_key(value: Any) -> str:
    return pathogen_names.canonical_key(value) or normalized_name(value)


def direct_genus_key(value: Any) -> str:
    return pathogen_names.genus_name(value) or genus_key(value)


def evidence_for_review_item(
    item: dict[str, Any],
    evidence_index: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    snapshot = item.get("evidence_snapshot")
    if isinstance(snapshot, dict) and snapshot:
        return snapshot
    for key in {
        normalized_name(item.get("organism_name")),
        canonical_name_key(item.get("organism_name")),
        *pathogen_names.expanded_name_keys(item.get("organism_name")),
    }:
        evidence = evidence_index.get(key)
        if evidence:
            return evidence
    return {}


def item_found_in_review_evidence(item: dict[str, Any], evidence: dict[str, Any]) -> bool:
    """Return True when the review item can be tied back to an input queue organism."""
    if evidence:
        return True
    name = clean_review_organism_name(item.get("organism_name") or item.get("name"))
    if not name:
        return False
    key = canonical_name_key(name)
    text = normalized_name(name)
    if not key and not text:
        return False
    # Allow only explicit no-organism placeholders to fall through to omitted.
    return False


ORGANISM_NAME_LABEL_PREFIX_TERMS = (
    "代表",
    "菌群",
    "相關菌",
    "oral",
    "aspiration",
    "flora",
    "representative",
    "pattern",
)


def clean_review_organism_name(value: Any) -> str:
    """Remove explanatory labels that some smaller models put into organism_name."""
    text = pathogen_names.clean_display_text(value)
    if not text:
        return ""
    text = re.sub(r"^\s*[-*•\d.)、\s]+", "", text).strip()
    for delimiter in ("：", ":"):
        if delimiter not in text:
            continue
        prefix, suffix = text.rsplit(delimiter, 1)
        if suffix.strip() and any(term in prefix.lower() for term in ORGANISM_NAME_LABEL_PREFIX_TERMS):
            text = suffix.strip()
            break
    return pathogen_names.clean_display_text(text)


def clean_review_item_organism_name(item: dict[str, Any]) -> dict[str, Any]:
    cleaned = clean_review_organism_name(item.get("organism_name") or item.get("name"))
    if not cleaned:
        return item
    output = dict(item)
    if output.get("organism_name") != cleaned:
        output["raw_organism_name_before_cleanup"] = output.get("organism_name") or output.get("name")
        output["organism_name_cleanup_applied"] = True
    output["organism_name"] = cleaned
    return output


def picked_name_sets_from_payload(payload: dict[str, Any]) -> tuple[set[str], set[str]]:
    picked = (
        ((payload.get("deterministic_formal_output") or {}).get("best_available_summary") or {}).get("picked_pathogens")
        or []
    )
    keys: set[str] = set()
    representative_groups: set[str] = set()
    if not isinstance(picked, list):
        return keys, representative_groups
    for item in picked:
        if not isinstance(item, dict):
            continue
        name = item.get("organism_name") or item.get("name")
        key = canonical_name_key(name)
        group = pathogen_names.representative_group_key(name)
        if key:
            keys.add(key)
        if group and not group.startswith("name:"):
            representative_groups.add(group)
    return keys, representative_groups


def direct_review_tier_for_item(
    item: dict[str, Any],
    *,
    source_section: str,
    evidence: dict[str, Any],
    payload: dict[str, Any],
    picked_keys: set[str],
    picked_genera: set[str],
) -> tuple[str, str]:
    name = item.get("organism_name") or evidence.get("organism_name")
    key = canonical_name_key(name)
    group = pathogen_names.representative_group_key(name)
    reads = numeric_value(evidence.get("reads"))
    dominance = str(evidence.get("dominance_tier") or "")
    confidence = str(item.get("confidence") or "").strip().lower()
    source_category = str(evidence.get("source_category") or "")

    impossible_rule = pathogen_rule_lists.impossible_infection_source_rule(name)
    if impossible_rule:
        return "omitted_with_reason", "impossible pulmonary infection source rule"
    if key and key in picked_keys:
        return "omitted_with_reason", "same pathogen or alias is already retained in picked_pathogens"
    if (
        group
        and group in picked_genera
        and pathogen_names.representative_group_suppression_allowed(name)
    ):
        return "omitted_with_reason", "same representative alias/generic label already has a picked representative"

    hospital_summary = payload.get("hospital_side_summary")
    hospital_items = (
        hospital_summary.get("hospital_organism_evidence")
        if isinstance(hospital_summary, dict)
        else []
    )
    exact_hospital_items = []
    name_keys = pathogen_names.expanded_name_keys(name)
    for hospital_item in hospital_items or []:
        if not isinstance(hospital_item, dict):
            continue
        hospital_keys = pathogen_names.expanded_name_keys(
            hospital_item.get("organism_name") or hospital_item.get("name")
        )
        if name_keys & hospital_keys:
            exact_hospital_items.append(hospital_item)
    serology_profile = evidence_profiles.atypical_respiratory_serology_profile(
        name,
        exact_hospital_items,
    )
    source_is_case_fit = source_section in {
        "review_high_priority",
        "review_context_needed",
        "possible_missed_pathogens",
    }
    pulmonary_context = evidence_profiles.has_compatible_pulmonary_context(
        {
            "hospital_side_summary": payload.get("hospital_side_summary"),
            "mngs_to_specimen": payload.get("mngs_to_specimen"),
        }
    )
    if serology_profile["context_rescue_match"] and source_is_case_fit and pulmonary_context:
        return (
            "review_context_needed",
            f"{serology_profile['rule_id']}: {serology_profile['reason']}",
        )

    if key == "corynebacteriumstriatum" and dominance == "D3_dominant" and reads >= 10000:
        return "review_high_priority", "D3 dominant, very high-burden lower-respiratory C. striatum safety signal"
    if key in ACTIONABLE_RESPIRATORY_REVIEW_KEYS and source_section in {
        "review_high_priority",
        "review_context_needed",
        "possible_missed_pathogens",
    }:
        return "review_high_priority", "actionable respiratory pathogen or respiratory panel target with non-picked evidence"

    if key == "corynebacteriumstriatum" and reads >= 1000:
        return "review_context_needed", "high-burden C. striatum needs airway-device, culture, imaging, and treatment-context review"
    if herpes_r4_nondominant_host_only_review(name, evidence):
        return (
            "review_context_needed",
            "R4 very high herpesvirus reads with D0/D1 non-dominant host-only support needs PCR/viral-load and clinical context review",
        )
    if herpes_ranked_host_context_review(name, evidence, payload):
        return (
            "review_context_needed",
            "HSV/CMV/VZV rank/reads/percentile signal in a highly vulnerable host needs PCR/viral-load and clinical context review; not formal picked without direct support",
        )
    rule_ids = evidence_rule_ids(evidence)
    if (
        NON_CORE_RESPIRATORY_VIRUS_FORMAL_PICK_GUARDRAIL in rule_ids
        or "non_core_respiratory_virus_context_not_formal_picked" in rule_ids
    ):
        return (
            "review_context_needed",
            "non-core respiratory virus was moved out of formal picked but remains RAG-visible for symptom timing, panel/PCR, imaging, and coinfection context",
        )
    if (
        MNGS_ONLY_NONDOMINANT_NON_TYPICAL_FORMAL_PICK_GUARDRAIL in rule_ids
        or "mngs_only_nondominant_non_typical_strong_context" in rule_ids
        or "mngs_only_nondominant_non_typical_low_specificity" in rule_ids
        or "broad_or_atypical_oral_mngs_only_low_specificity" in rule_ids
    ):
        if (
            "broad_or_atypical_oral_mngs_only_low_specificity" in rule_ids
            or (
                is_oral_aspiration_flora_name(name)
                and oral_aspiration_flora_category(name)
                in {"broad_oral_upper_airway_commensal", "atypical_oral_associated"}
                and not has_oral_aspiration_context_signal(name, evidence, payload)
            )
        ):
            return (
                "review_low_specificity",
                "broad oral/upper-airway or atypical oral-associated organism has mNGS-only D0/D1 high-percentile signal but lacks same-organism support, D2/D3 dominance, sterile/systemic signal, or coherent aspiration pattern",
            )
        if (
            "mngs_only_nondominant_non_typical_strong_context" in rule_ids
            or tier_rank(evidence.get("reads_tier")) >= 3
            or numeric_value(evidence.get("reads_percentile")) >= 0.85
        ):
            return (
                "review_context_needed",
                "mNGS-only D0/D1 non-dominant non-typical organism has high reads tier or percentile; review as context rather than formal picked",
            )
        return (
            "review_low_specificity",
            "mNGS-only D0/D1 non-dominant non-typical organism lacks same-organism hospital support and does not have strong enough signal for RAG-visible context",
        )
    if (
        COMMON_HOSPITAL_LOW_READ_NODOM_FORMAL_PICK_GUARDRAIL in rule_ids
        or "common_hospital_low_read_nondominant_related_context" in rule_ids
        or "common_hospital_low_read_nondominant_low_specificity" in rule_ids
    ):
        if (
            "common_hospital_low_read_nondominant_related_context" in rule_ids
            or has_related_representative_context(evidence)
        ):
            return (
                "review_context_needed",
                "common hospital pneumonia pathogen has only R0/R1 non-dominant mNGS signal, but related representative hospital context merits RAG review",
            )
        return (
            "review_low_specificity",
            "common hospital pneumonia pathogen has only R0/R1 non-dominant mNGS signal without direct same-organism or related hospital support; retain as audit rather than RAG-visible context",
        )
    if is_rare_opportunistic_yeast_name(name):
        if has_non_host_support(evidence) or tier_rank(evidence.get("reads_tier")) >= 3:
            return (
                "review_context_needed",
                "rare opportunistic yeast needs host, lower-respiratory culture, imaging, and treatment-context review",
            )
        return (
            "review_low_specificity",
            "rare opportunistic yeast has weak or mNGS-only evidence; do not treat as formal pulmonary attribution without context",
        )
    if is_candida_or_generic_yeast_name(name):
        if candida_respiratory_high_priority_review(name, evidence, payload):
            return (
                "review_high_priority",
                "BAL/lower-respiratory Candida mNGS signal plus highly vulnerable host and same-species or same-genus respiratory hospital context; prioritize RAG while keeping it out of formal picked without same-species strong invasive proof",
            )
        if candida_has_invasive_context(evidence):
            return (
                "review_high_priority",
                "same-species Candida has strong invasive evidence such as blood, sterile-site, tissue, pleural, abscess, or histopathology; high-priority review is needed but this is not automatic Candida pneumonia attribution",
            )
        if candida_has_intermediate_context(evidence):
            return (
                "review_context_needed",
                "Candida has intermediate systemic or invasive-candidiasis clues; review clinical context, source, and pulmonary linkage before attribution",
            )
        if candida_respiratory_context_review(name, evidence, payload):
            return (
                "review_context_needed",
                "respiratory Candida/generic yeast has lower-respiratory mNGS or vulnerable-host context, but lacks same-species strong invasive evidence; RAG should assess colonization versus infection",
            )
        return (
            "review_low_specificity",
            "respiratory or non-invasive Candida/generic yeast is low-specificity for pneumonia and should not enter RAG by default",
        )
    if is_oral_aspiration_flora_name(name):
        if has_oral_aspiration_high_priority_signal(name, evidence, payload):
            return (
                "review_high_priority",
                "oral/aspiration anaerobe has unusually strong convergence such as lower-respiratory hospital support, sterile/systemic signal, D3 high-burden dominance, or compatible abscess/empyema/necrotizing aspiration pattern",
            )
        if has_oral_aspiration_context_signal(name, evidence, payload):
            return (
                "review_context_needed",
                "oral/aspiration anaerobe or flora has lower-respiratory clean high mNGS signal, same-organism support, D2/D3 dominance, sterile/systemic signal, or a coherent aspiration pattern; RAG context is needed but it is not a formal picked pathogen",
            )
        return (
            "review_low_specificity",
            "single oral/upper-airway/aspiration flora signal lacks clean high lower-respiratory mNGS signal, same-organism lower-respiratory support, D2/D3 dominance, abscess/empyema/necrotizing pneumonia, or coherent multi-organism aspiration pattern with compatible clinical context",
        )
    if is_water_or_low_specificity_environmental_name(name):
        if weak_hospital_only_environmental_signal(name, evidence, payload):
            return (
                "review_low_specificity",
                "weak hospital-only water/environmental signal lacks repeated/effective respiratory culture, sterile-site evidence, or mNGS support",
            )
        if has_strong_environmental_signal(evidence):
            return (
                "review_context_needed",
                "water/environmental low-pulmonary-specificity organism has strong enough burden, dominance, non-host support, sterile/systemic context, or source/device context for RAG review",
            )
        return (
            "review_low_specificity",
            "water/environmental low-pulmonary-specificity organism lacks same-organism hospital support, D2/D3 dominance, very high mNGS burden, or source/device context",
        )
    if weak_l3_host_only_nondominant_context(evidence):
        return (
            "review_low_specificity",
            "Level 3 mNGS-only host-context signal is low/non-dominant and lacks non-host support; retain as audit rather than RAG-visible",
        )
    if weak_trace_protected_context(evidence):
        return (
            "review_low_specificity",
            "trace protected-pathogen signal is non-dominant and lacks non-host support; retain as audit unless stronger microbiologic evidence appears",
        )
    if weak_high_risk_mold_context(evidence):
        return (
            "review_low_specificity",
            "high-risk mold/dimorphic fungus signal is low-to-medium burden, non-dominant, and lacks non-host support; retain as audit rather than RAG-visible",
        )
    if key in CONTEXT_REVIEW_KEYS:
        return "review_context_needed", "clinically plausible or high-consequence organism with incomplete pulmonary evidence"
    if key == "chryseobacteriumindologenes" and reads >= 10000:
        if has_non_host_support(evidence) or dominance in {"D2_moderate", "D3_dominant"}:
            return "review_context_needed", "water/environment-associated GNB has unusually strong local convergence"
        return (
            "review_low_specificity",
            "high-read Chryseobacterium without same-organism support or D2/D3 dominance remains low pulmonary specificity",
        )
    if key == "elizabethkingiaanophelis" and confidence == "moderate":
        if has_non_host_support(evidence) or dominance in {"D2_moderate", "D3_dominant"}:
            return "review_context_needed", "potential nosocomial water-associated GNB has convergent support"
        return (
            "review_low_specificity",
            "Elizabethkingia signal lacks same-organism support or D2/D3 dominance; retain as low-specificity audit",
        )

    if source_section in {"omitted_with_reason", "not_recommended_for_upgrade"} and source_category == "hospital_detected_only":
        return "omitted_with_reason", "hospital-only finding has insufficient pulmonary relevance"
    return (
        "review_low_specificity",
        "low-specificity, colonizer/background/environmental, non-pulmonary culture-only, or weak mNGS signal; not sent to RAG",
    )


def ensure_review_item_fields(
    item: dict[str, Any],
    *,
    tier: str,
    reason: str,
    source_section: str,
    evidence: dict[str, Any],
) -> dict[str, Any]:
    output = dict(item)
    name = output.get("organism_name") or evidence.get("organism_name")
    if name:
        output["organism_name"] = pathogen_names.display_name(name)
        output["canonical_key"] = canonical_name_key(name)
    output["review_tier"] = tier
    output["review_tier_reason"] = reason
    output["original_review_section"] = source_section
    output.setdefault("confidence", "low")
    if tier != "omitted_with_reason":
        output.setdefault("rationale_zh", reason)
        output.setdefault(
            "risk_if_ignored",
            "若此項目代表真實共同病原，可能低估混合感染或特殊宿主風險；目前證據仍不足以作為正式病原。",
        )
        output.setdefault(
            "why_not_auto_upgrade",
            "仍需臨床脈絡、有效檢體、重複檢出、優勢度、培養或其他佐證，不能只憑單一低特異性訊號自動升級。",
        )
        output.setdefault(
            "suggested_next_check",
            "核對檢體來源與品質、影像、同菌培養/分子證據、重複檢出、治療反應及宿主風險。",
        )
    else:
        output.setdefault("rationale_zh", reason)
    if evidence and not isinstance(output.get("evidence_snapshot"), dict):
        output["evidence_snapshot"] = evidence_snapshot_for_review(evidence)
    return output


def review_item_direct_score(
    item: dict[str, Any],
    tier: str,
    evidence_index: dict[str, dict[str, Any]],
) -> float:
    evidence = evidence_for_review_item(item, evidence_index)
    rank = numeric_value(evidence.get("rank_priority"), default=99.0)
    reads = numeric_value(evidence.get("reads"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    return (
        TIER_PRIORITY.get(tier, 0) * 10000.0
        + confidence_score(item.get("confidence"))
        + level_score(evidence.get("integrated_causative_level"))
        + tier_score(evidence.get("mngs_signal_tier"))
        + tier_score(evidence.get("reads_tier"))
        + dominance_score(evidence.get("dominance_tier"))
        + support_score(evidence.get("support_modules"))
        + min(reads, 100000.0) / 1000.0
        + percentile * 100.0
        + max(0.0, 80.0 - rank * 10.0)
    )


def consolidate_direct_review_items(review: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    evidence_index = build_queue_evidence_index(payload)
    tracked: list[dict[str, Any]] = []
    for tier in DIRECT_ACTIVE_TIERS:
        values = review.get(tier) or []
        if not isinstance(values, list):
            review[tier] = []
            continue
        for index, item in enumerate(values):
            if not isinstance(item, dict):
                continue
            tracked.append(
                {
                    "tier": tier,
                    "index": index,
                    "item": item,
                    "key": canonical_name_key(item.get("organism_name")),
                    "genus": direct_genus_key(item.get("organism_name")),
                    "representative_group": pathogen_names.representative_group_key(item.get("organism_name")),
                    "score": review_item_direct_score(item, tier, evidence_index),
                }
            )

    winners_by_key: dict[str, dict[str, Any]] = {}
    removed: list[dict[str, Any]] = []
    for entry in tracked:
        key = entry["key"]
        if not key:
            continue
        current = winners_by_key.get(key)
        if current is None or entry["score"] > current["score"]:
            if current is not None:
                removed.append(current)
            winners_by_key[key] = entry
        else:
            removed.append(entry)

    active_after_key = [entry for entry in tracked if entry not in removed]
    grouped_by_representative: dict[str, list[dict[str, Any]]] = {}
    for entry in active_after_key:
        group = entry.get("representative_group")
        if not group or str(group).startswith("name:"):
            continue
        grouped_by_representative.setdefault(str(group), []).append(entry)

    winners_by_representative: dict[str, dict[str, Any]] = {}
    for group, entries in grouped_by_representative.items():
        suppressible = [
            entry
            for entry in entries
            if pathogen_names.representative_group_suppression_allowed(
                (entry.get("item") or {}).get("organism_name")
            )
        ]
        if not suppressible:
            continue
        retained_pool = [
            entry
            for entry in entries
            if not pathogen_names.representative_group_suppression_allowed(
                (entry.get("item") or {}).get("organism_name")
            )
        ] or entries
        winner = max(retained_pool, key=lambda entry: entry["score"])
        winners_by_representative[group] = winner
        for entry in suppressible:
            if entry is not winner:
                removed.append(entry)

    removed_ids = {(entry["tier"], entry["index"]) for entry in removed}
    for tier in DIRECT_ACTIVE_TIERS:
        values = review.get(tier) or []
        review[tier] = [
            item
            for index, item in enumerate(values)
            if isinstance(item, dict) and (tier, index) not in removed_ids
        ]

    omitted = review.setdefault("omitted_with_reason", [])
    if not isinstance(omitted, list):
        omitted = []
        review["omitted_with_reason"] = omitted
    for entry in removed:
        winner = (
            winners_by_key.get(entry["key"])
            or winners_by_representative.get(entry.get("representative_group"))
            or {}
        )
        winner_name = ((winner or {}).get("item") or {}).get("organism_name")
        item = dict(entry["item"] or {})
        item["review_tier"] = "omitted_with_reason"
        item["original_review_section"] = entry["tier"]
        item["same_representative_or_alias_representative"] = winner_name
        item["rationale_zh"] = "同名 alias 或同屬代表已保留；此項目只作為 related/audit evidence，不送 RAG。"
        omitted.append(item)

    review["same_genus_consolidation"] = {
        "enabled": True,
        "policy": "Direct review tiers consolidate exact aliases and broad generic labels only; distinct species in the same genus remain visible.",
        "removed_count": len(removed),
    }
    return review


def consolidate_direct_oral_aspiration_pattern(review: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    evidence_index = build_queue_evidence_index(payload)
    tracked: list[dict[str, Any]] = []
    for tier in DIRECT_ACTIVE_TIERS:
        values = review.get(tier) or []
        if not isinstance(values, list):
            continue
        for index, item in enumerate(values):
            if not isinstance(item, dict) or not is_oral_aspiration_flora_name(item.get("organism_name")):
                continue
            tracked.append(
                {
                    "tier": tier,
                    "index": index,
                    "item": item,
                    "score": review_item_direct_score(item, tier, evidence_index),
                }
            )

    if len(tracked) <= 2:
        for entry in tracked:
            entry["item"]["oral_aspiration_flora_pattern_representative"] = True
        review["oral_aspiration_flora_pattern"] = {
            "enabled": True,
            "pattern_detected": bool(tracked),
            "policy": "Keep only 1-2 representative oral/aspiration flora organisms in direct review tiers.",
            "representative_organisms": [entry["item"].get("organism_name") for entry in tracked],
            "grouped_non_representative_organisms": [],
            "interpretation_zh": "偵測到口腔/吸入相關菌叢代表項目，需結合吸入風險、影像與厭氧感染證據判讀。" if tracked else "",
        }
        return review

    ranked = sorted(tracked, key=lambda entry: entry["score"], reverse=True)
    keep = ranked[:2]
    remove = ranked[2:]
    remove_ids = {(entry["tier"], entry["index"]) for entry in remove}
    keep_ids = {(entry["tier"], entry["index"]) for entry in keep}
    grouped = []
    for tier in DIRECT_ACTIVE_TIERS:
        values = review.get(tier) or []
        new_values = []
        for index, item in enumerate(values):
            if not isinstance(item, dict):
                continue
            item_id = (tier, index)
            if item_id in remove_ids:
                grouped.append({"organism_name": item.get("organism_name"), "original_section": tier})
                continue
            if item_id in keep_ids:
                item = dict(item)
                item["oral_aspiration_flora_pattern_representative"] = True
            new_values.append(item)
        review[tier] = new_values

    review["oral_aspiration_flora_pattern"] = {
        "enabled": True,
        "pattern_detected": True,
        "policy": "Keep only 1-2 representative oral/aspiration flora organisms in direct review tiers.",
        "representative_organisms": [entry["item"].get("organism_name") for entry in keep],
        "grouped_non_representative_organisms": grouped,
        "interpretation_zh": "多個口腔/吸入相關菌同時出現時，輸出以代表菌與菌叢型態呈現，避免逐隻放大低特異性訊號。",
    }
    return review


def normalize_direct_review_output(review: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    evidence_index = build_queue_evidence_index(payload)
    picked_keys, picked_genera = picked_name_sets_from_payload(payload)
    output = dict(review)
    buckets: dict[str, list[dict[str, Any]]] = {tier: [] for tier in DIRECT_REVIEW_TIERS}
    input_counts: dict[str, int] = {}
    non_queue_items: list[dict[str, Any]] = []

    for source_section in (*DIRECT_REVIEW_TIERS, *LEGACY_REVIEW_SECTIONS):
        values = review.get(source_section) or []
        if not isinstance(values, list):
            input_counts[source_section] = 0
            continue
        input_counts[source_section] = len(values)
        for item in values:
            if not isinstance(item, dict):
                continue
            item = clean_review_item_organism_name(item)
            evidence = evidence_for_review_item(item, evidence_index)
            if not item_found_in_review_evidence(item, evidence):
                omitted_item = dict(item)
                omitted_item["postprocess_guardrail"] = "REVIEW_ITEM_NOT_IN_INPUT_QUEUE"
                omitted_item["review_tier"] = "omitted_with_reason"
                omitted_item["review_tier_reason"] = (
                    "LLM output organism could not be matched to the review queue or evidence index"
                )
                omitted_item.setdefault(
                    "rationale_zh",
                    "LLM 輸出的菌名無法對回本次 review queue 或 evidence index，因此只保留在 audit，不送 RAG。",
                )
                buckets["omitted_with_reason"].append(omitted_item)
                non_queue_items.append(
                    {
                        "organism_name": omitted_item.get("organism_name") or omitted_item.get("name"),
                        "source_section": source_section,
                    }
                )
                continue
            tier, reason = direct_review_tier_for_item(
                item,
                source_section=source_section,
                evidence=evidence,
                payload=payload,
                picked_keys=picked_keys,
                picked_genera=picked_genera,
            )
            buckets[tier].append(
                ensure_review_item_fields(
                    item,
                    tier=tier,
                    reason=reason,
                    source_section=source_section,
                    evidence=evidence,
                )
            )

    for key in (*DIRECT_REVIEW_TIERS, *LEGACY_REVIEW_SECTIONS):
        output.pop(key, None)
    output.update(buckets)
    output["review_tiering_policy"] = {
        "enabled": True,
        "version": "direct_review_triage_high_context_only_v11_luna_rescue_cleanup",
        "input_counts": input_counts,
        "rag_visible_tiers": ["review_high_priority", "review_context_needed"],
        "audit_only_tiers": ["review_low_specificity", "omitted_with_reason"],
        "counts_before_consolidation": {key: len(value) for key, value in buckets.items()},
        "non_queue_item_guardrail": {
            "enabled": True,
            "removed_from_rag_visible_count": len(non_queue_items),
            "items": non_queue_items,
        },
        "policy_zh": (
            "LLM review 直接輸出 review_high_priority、review_context_needed、review_low_specificity、omitted_with_reason。"
            "只有前兩層送 RAG；低特異性與省略項目只保留 audit。"
        ),
    }
    output = consolidate_direct_review_items(output, payload)
    output = consolidate_direct_oral_aspiration_pattern(output, payload)
    output = apply_strict_aspiration_anaerobe_rescue(output, payload)
    policy = output.get("review_tiering_policy")
    if isinstance(policy, dict):
        policy["counts_after_consolidation"] = {
            key: len(output.get(key) or []) if isinstance(output.get(key), list) else 0
            for key in DIRECT_REVIEW_TIERS
        }
    return output


def move_not_recommended_signal_to_watchlist(
    review: dict[str, Any],
    payload: dict[str, Any],
    *,
    predicate,
    guardrail_name: str,
    reason_zh: str,
) -> tuple[int, list[str]]:
    evidence_index = build_queue_evidence_index(payload)
    values = review.get("not_recommended_for_upgrade") or []
    if not isinstance(values, list):
        review["not_recommended_for_upgrade"] = []
        return 0, []
    kept: list[dict[str, Any]] = []
    moved = 0
    moved_names: list[str] = []
    for item in values:
        if not isinstance(item, dict):
            continue
        name = item.get("organism_name")
        evidence = evidence_index.get(normalized_name(name), {})
        if predicate(evidence):
            watch_item = dict(item)
            watch_item["confidence"] = watch_item.get("confidence") or "low"
            watch_item["postprocess_guardrail"] = guardrail_name
            watch_item["original_not_recommended_rationale_zh"] = item.get("rationale_zh")
            watch_item["rationale_zh"] = reason_zh
            watch_item.setdefault(
                "risk_if_ignored",
                "May underestimate a true co-pathogen, aspiration-pattern signal, or low-specificity environmental signal; evidence is still insufficient for primary attribution.",
            )
            watch_item.setdefault(
                "why_not_auto_upgrade",
                "Watchlist only: possible_missed/picked still requires dominance, same-organism hospital support, or a clear compatible clinical syndrome.",
            )
            watch_item.setdefault(
                "suggested_next_check",
                "Review specimen quality, imaging for abscess/empyema/necrotizing pneumonia/aspiration, culture quantity/same organism, and treatment response.",
            )
            watch_item.setdefault(
                "evidence_snapshot",
                {
                    "rank_priority": evidence.get("rank_priority"),
                    "rank_rule": evidence.get("rank_rule"),
                    "reads": evidence.get("reads"),
                    "reads_tier": evidence.get("reads_tier"),
                    "reads_percentile": evidence.get("reads_percentile"),
                    "dominance_tier": evidence.get("dominance_tier"),
                    "specimen_class": evidence.get("specimen_class"),
                    "support_modules": evidence.get("support_modules"),
                    "non_host_support_modules": evidence.get("non_host_support_modules"),
                    "review_reasons": evidence.get("review_reasons"),
                }
                if evidence
                else {},
            )
            watchlist = review.setdefault("watchlist_candidates", [])
            if isinstance(watchlist, list):
                watchlist.append(watch_item)
            moved += 1
            moved_names.append(str(name or evidence.get("organism_name") or ""))
            continue
        kept.append(item)
    review["not_recommended_for_upgrade"] = kept
    return moved, moved_names


def is_gi_only_pathogen(name: Any, evidence: dict[str, Any]) -> bool:
    normalized = normalized_name(name)
    if normalized not in GI_PATHOGEN_EXACT_NAMES and not any(
        normalized.startswith(prefix) for prefix in ("clostridioidesdifficile", "clostridiumdifficile")
    ):
        return False
    haystack = recursive_text(evidence).lower()
    has_gi_context = any(term in haystack for term in GI_CONTEXT_TERMS)
    has_pulmonary_or_systemic_context = any(term in haystack for term in PULMONARY_OR_SYSTEMIC_CONTEXT_TERMS)
    has_mngs_signal = numeric_value(evidence.get("reads")) > 0 and str(evidence.get("rank_priority") or "") != "Not_available"
    if has_gi_context and not has_pulmonary_or_systemic_context:
        return True
    if normalized in GI_PATHOGEN_EXACT_NAMES and evidence.get("source_category") == "hospital_detected_only" and not has_mngs_signal:
        return True
    return False


def host_is_highly_vulnerable(payload: dict[str, Any]) -> bool:
    host = payload.get("host_context")
    if not isinstance(host, dict):
        review_queue = payload.get("review_queue")
        if isinstance(review_queue, dict):
            host = review_queue.get("host_context")
    if not isinstance(host, dict):
        return False
    return bool(
        str(host.get("host_vulnerability_tier") or "").upper() == "V3"
        or str(host.get("opportunistic_coverage_level") or "").upper() in {"O2", "O3"}
    )


def has_formal_picked_pathogen(payload: dict[str, Any]) -> bool:
    review_queue = payload.get("review_queue")
    if isinstance(review_queue, dict):
        try:
            return int(review_queue.get("formal_picked_count") or 0) > 0
        except (TypeError, ValueError):
            return False
    return False


def weak_hospital_only_environmental_signal(name: Any, evidence: dict[str, Any], payload: dict[str, Any]) -> bool:
    """Generic rule for water/environment taxa detected only in weak non-sterile respiratory culture."""
    if not is_water_or_low_specificity_environmental_name(name):
        return False
    if evidence.get("source_category") != "hospital_detected_only":
        return False
    if numeric_value(evidence.get("reads")) > 0 or str(evidence.get("rank_priority") or "") != "Not_available":
        return False
    haystack = recursive_text(evidence).lower()
    has_sterile_or_systemic = any(
        term in haystack
        for term in ("blood", "csf", "tissue", "abscess", "pus", "pleural", "sterile")
    )
    if has_sterile_or_systemic:
        return False
    is_nonsterile_respiratory = any(term in haystack for term in ("sputum", "respiratory", "tracheal", "eta"))
    if not is_nonsterile_respiratory:
        return False
    quantified_or_repeated = any(
        term in haystack
        for term in (
            "heavy",
            "many",
            "moderate",
            "3+",
            "4+",
            "high",
            "abundant",
            "predominant",
            "pure",
            "repeated",
            "repeat",
        )
    )
    if "balf" in haystack and quantified_or_repeated:
        return False
    if quantified_or_repeated:
        return False
    # Do not demote solely because another formal pathogen is already picked.
    # A weak hospital-only environmental signal is excluded because the signal itself is weak;
    # highly vulnerable hosts are kept for manual review even with weak evidence.
    if host_is_highly_vulnerable(payload):
        return False
    return True


def add_not_recommended(review: dict[str, Any], item: dict[str, Any], rationale: str, extra: dict[str, Any] | None = None) -> None:
    not_recommended = review.setdefault("not_recommended_for_upgrade", [])
    if not isinstance(not_recommended, list):
        not_recommended = []
        review["not_recommended_for_upgrade"] = not_recommended
    payload = {
        "organism_name": item.get("organism_name"),
        "rationale_zh": rationale,
    }
    if extra:
        payload.update(extra)
    existing_keys = {
        (normalized_name(entry.get("organism_name")), entry.get("rationale_zh"))
        for entry in not_recommended
        if isinstance(entry, dict)
    }
    key = (normalized_name(payload.get("organism_name")), payload.get("rationale_zh"))
    if key not in existing_keys:
        not_recommended.append(payload)


def review_reasons_text(evidence: dict[str, Any]) -> str:
    reasons = evidence.get("review_reasons")
    if isinstance(reasons, list):
        return " ".join(str(reason or "") for reason in reasons)
    return str(reasons or "")


def has_clean_ntc_rank_rule(evidence: dict[str, Any]) -> bool:
    rank_rule = str(evidence.get("rank_rule") or "").lower()
    return "ntc=00" in rank_rule or "rk_ntc=00" in rank_rule or "rk00" in rank_rule



def high_burden_contaminant_possible_missed_reason(evidence: dict[str, Any], payload: dict[str, Any]) -> str | None:
    name = evidence.get("organism_name")
    if not is_high_burden_contaminant_risk_name(name):
        return None

    rank = numeric_value(evidence.get("rank_priority"), default=99.0)
    reads = numeric_value(evidence.get("reads"))
    reads_rank = tier_rank(evidence.get("reads_tier"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    dominance = str(evidence.get("dominance_tier") or "")
    specimen_cls = evidence.get("specimen_class")
    has_aux_support = has_non_host_support(evidence)
    clean_ntc = has_clean_ntc_rank_rule(evidence)

    if specimen_cls != "S2_lower_respiratory" or reads <= 0:
        return None
    if rank > 2 or reads_rank < 3 or percentile < 0.95 or not clean_ntc:
        return None

    strong_mngs_burden = reads_rank >= 4 or reads >= 50000
    strong_context = (
        strong_mngs_burden
        or dominance in {"D2_moderate", "D3_dominant"}
        or has_aux_support
        or host_is_highly_vulnerable(payload)
    )
    if not strong_context:
        return None

    return (
        "High-burden contaminant-risk organism kept as possible_missed: lower-respiratory specimen, "
        "clean NTC/RK_NTC, front-rank signal, reads_tier >= R3, percentile >= 0.95, and at least one stronger support "
        "feature such as R4/reads>=50000, D2/D3 dominance, same-organism hospital support, or high host vulnerability. "
        "This remains non-primary and does not become formal picked without culture/Gram stain/clinical confirmation."
    )

def high_burden_contaminant_watchlist_reason(evidence: dict[str, Any]) -> str | None:
    name = evidence.get("organism_name")
    if not is_high_burden_contaminant_risk_name(name):
        return None

    rank = numeric_value(evidence.get("rank_priority"), default=99.0)
    reads = numeric_value(evidence.get("reads"))
    reads_rank = tier_rank(evidence.get("reads_tier"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    dominance = str(evidence.get("dominance_tier") or "")
    specimen_cls = evidence.get("specimen_class")
    reasons_text = review_reasons_text(evidence)
    has_aux_support = has_non_host_support(evidence)
    clean_ntc = has_clean_ntc_rank_rule(evidence)

    if reads <= 0:
        return None

    if specimen_cls == "S2_lower_respiratory":
        very_high_burden = reads_rank >= 4 or reads >= 50000 or percentile >= 0.90
        if rank <= 2 and (very_high_burden or dominance in {"D2_moderate", "D3_dominant"} or has_aux_support):
            return (
                "高 mNGS 負荷但低特異性菌種 watchlist：此菌屬 CoNS/skin flora 或 Corynebacterium 類群，"
                "在呼吸道檢體可能是污染或定植，但因 rank 靠前且 reads/percentile 很高，仍應列入人工覆核。"
            )
        if "full_ranked_common_background_high_burden_review" in reasons_text and very_high_burden:
            return (
                "高負荷背景風險菌 watchlist：原始 full-ranked queue 已標記為 common/background high-burden review，"
                "不建議直接升級為主病原，但應保留供人工確認。"
            )

    if specimen_cls == "S1_sterile_systemic" and is_non_aureus_staphylococcus_name(name):
        sterile_signal = rank <= 2 and (percentile >= 0.75 or reads >= 50 or reads_rank >= 1)
        if sterile_signal and (clean_ntc or "pulmonary_opt_removed_sterile_systemic_non_aureus_staphylococcus" in reasons_text):
            return (
                "無菌/系統性檢體 non-aureus Staphylococcus watchlist：此類菌仍有污染風險，"
                "但 rank 靠前且 NTC/percentile 條件支持 mNGS 訊號，應列入人工覆核而非完全排除。"
            )

    return None



def evidence_snapshot_for_review(evidence: dict[str, Any]) -> dict[str, Any]:
    return {
        "rank_priority": evidence.get("rank_priority"),
        "reads": evidence.get("reads"),
        "reads_tier": evidence.get("reads_tier"),
        "reads_percentile": evidence.get("reads_percentile"),
        "dominance_tier": evidence.get("dominance_tier"),
        "specimen_class": evidence.get("specimen_class"),
        "rank_rule": evidence.get("rank_rule"),
        "support_modules": evidence.get("support_modules"),
        "non_host_support_modules": evidence.get("non_host_support_modules"),
        "related_representative_hospital_support": evidence.get("related_representative_hospital_support"),
        "candida_invasive_hospital_support": evidence.get("candida_invasive_hospital_support"),
        "candida_related_invasive_hospital_context": evidence.get("candida_related_invasive_hospital_context"),
        "candida_evidence_strength_profile": evidence.get("candida_evidence_strength_profile"),
        "candida_related_evidence_strength_profile": evidence.get("candida_related_evidence_strength_profile"),
        "review_reasons": evidence.get("review_reasons"),
    }


def move_item_to_possible_missed(
    review: dict[str, Any],
    item: dict[str, Any],
    evidence: dict[str, Any],
    reason: str,
) -> bool:
    name_key = normalized_name(item.get("organism_name") or evidence.get("organism_name"))
    if not name_key:
        return False
    possible = review.get("possible_missed_pathogens")
    watchlist = review.get("watchlist_candidates")
    if not isinstance(possible, list):
        possible = []
        review["possible_missed_pathogens"] = possible
    if not isinstance(watchlist, list):
        watchlist = []
        review["watchlist_candidates"] = watchlist

    # If the organism is already in possible_missed, do not duplicate it.
    if any(normalized_name(entry.get("organism_name")) == name_key for entry in possible if isinstance(entry, dict)):
        return False

    # Remove the same organism from watchlist when upgrading it.
    review["watchlist_candidates"] = [
        entry
        for entry in watchlist
        if not (isinstance(entry, dict) and normalized_name(entry.get("organism_name")) == name_key)
    ]

    possible_item = dict(item)
    possible_item["organism_name"] = item.get("organism_name") or evidence.get("organism_name")
    possible_item["confidence"] = possible_item.get("confidence") or "low"
    possible_item["rationale_zh"] = reason
    possible_item["postprocess_guardrail"] = "HIGH_BURDEN_CONTAMINANT_RISK_POSSIBLE_MISSED"
    possible_item["why_not_auto_upgrade"] = (
        "Possible missed only: this organism has high colonization/contamination risk and still requires culture, Gram stain, "
        "D2/D3 dominance, sterile-site evidence, or clinical response support before formal picked attribution."
    )
    possible_item["evidence_snapshot"] = evidence_snapshot_for_review(evidence)
    possible.append(possible_item)
    return True

def move_item_to_watchlist(
    review: dict[str, Any],
    item: dict[str, Any],
    evidence: dict[str, Any],
    reason: str,
) -> bool:
    name_key = normalized_name(item.get("organism_name") or evidence.get("organism_name"))
    if not name_key:
        return False
    possible = review.get("possible_missed_pathogens")
    watchlist = review.get("watchlist_candidates")
    if not isinstance(possible, list):
        possible = []
        review["possible_missed_pathogens"] = possible
    if not isinstance(watchlist, list):
        watchlist = []
        review["watchlist_candidates"] = watchlist
    existing = {
        normalized_name(entry.get("organism_name"))
        for entry in [*possible, *watchlist]
        if isinstance(entry, dict)
    }
    if name_key in existing:
        return False

    watch_item = dict(item)
    watch_item["organism_name"] = item.get("organism_name") or evidence.get("organism_name")
    watch_item["rationale_zh"] = reason
    watch_item["original_not_recommended_rationale_zh"] = item.get("rationale_zh")
    watch_item["postprocess_guardrail"] = "HIGH_BURDEN_CONTAMINANT_RISK_WATCHLIST"
    watch_item["why_not_auto_upgrade"] = (
        "此菌屬污染/定植風險較高的低特異性菌，現階段只列 watchlist；"
        "Upgrade requires effective single culture support such as high/semquant BALF or ETA, acceptable-quality sputum moderate/heavy growth, repeated culture, Gram stain, specimen quality, device/ventilator context, or treatment-response support."
    )
    watch_item["evidence_snapshot"] = evidence_snapshot_for_review(evidence)
    watchlist.append(watch_item)
    return True


def apply_high_burden_contaminant_watchlist_guardrail(review: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    """Recover high-burden contaminant-risk taxa as advisory review items, not formal picks."""
    evidence_index = build_queue_evidence_index(payload)
    possible_moved = 0
    possible_names: list[str] = []
    watch_moved = 0
    watch_names: list[str] = []

    # Upgrade existing watchlist items to possible_missed only when their own signal is very strong.
    watchlist = review.get("watchlist_candidates") or []
    if not isinstance(watchlist, list):
        review["watchlist_candidates"] = []
        watchlist = []
    for item in list(watchlist):
        if not isinstance(item, dict):
            continue
        name = item.get("organism_name")
        evidence = evidence_index.get(normalized_name(name), {})
        reason = high_burden_contaminant_possible_missed_reason(evidence, payload)
        if reason and move_item_to_possible_missed(review, item, evidence, reason):
            possible_moved += 1
            possible_names.append(str(name or evidence.get("organism_name") or ""))

    not_recommended = review.get("not_recommended_for_upgrade") or []
    if not isinstance(not_recommended, list):
        review["not_recommended_for_upgrade"] = []
        not_recommended = []

    kept_not_recommended: list[dict[str, Any]] = []
    for item in not_recommended:
        if not isinstance(item, dict):
            continue
        name = item.get("organism_name")
        evidence = evidence_index.get(normalized_name(name), {})
        possible_reason = high_burden_contaminant_possible_missed_reason(evidence, payload)
        if possible_reason and move_item_to_possible_missed(review, item, evidence, possible_reason):
            possible_moved += 1
            possible_names.append(str(name or evidence.get("organism_name") or ""))
            continue
        watch_reason = high_burden_contaminant_watchlist_reason(evidence)
        if watch_reason and move_item_to_watchlist(review, item, evidence, watch_reason):
            watch_moved += 1
            watch_names.append(str(name or evidence.get("organism_name") or ""))
            continue
        kept_not_recommended.append(item)

    review["not_recommended_for_upgrade"] = kept_not_recommended
    review["high_burden_contaminant_risk_watchlist_guardrail"] = {
        "enabled": True,
        "policy": (
            "High-burden contaminant-risk organisms such as Corynebacterium/Cutibacterium/non-aureus Staphylococcus "
            "remain non-formal-pick. Very strong clean lower-respiratory mNGS signal can move them to possible_missed; "
            "moderate high-burden signal stays watchlist. Formal picked still requires culture, Gram stain, dominance, "
            "sterile-site evidence, or compatible clinical confirmation."
        ),
        "moved_to_possible_missed_count": possible_moved,
        "moved_to_possible_missed": possible_names,
        "moved_to_watchlist_count": watch_moved,
        "moved_to_watchlist": watch_names,
    }
    return review


def apply_impossible_infection_source_guardrail(review: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    """Remove user-curated impossible pulmonary infection-source organisms from advisory output."""
    removed: list[str] = []
    for section in ("possible_missed_pathogens", "watchlist_candidates"):
        values = review.get(section) or []
        if not isinstance(values, list):
            review[section] = []
            continue
        kept: list[dict[str, Any]] = []
        for item in values:
            if not isinstance(item, dict):
                continue
            rule = pathogen_rule_lists.impossible_infection_source_rule(item.get("organism_name"))
            if rule:
                add_not_recommended(
                    review,
                    item,
                    str(rule.get("rationale_zh") or "???????????????????? possible_missed ? watchlist?"),
                    {
                        "postprocess_guardrail": "IMPOSSIBLE_INFECTION_SOURCE_NOT_RECOMMENDED",
                        "rule_file": str(pathogen_rule_lists.DEFAULT_IMPOSSIBLE_RULE_PATH),
                        "rule_organism_name": rule.get("organism_name"),
                    },
                )
                removed.append(str(item.get("organism_name") or rule.get("organism_name") or ""))
                continue
            kept.append(item)
        review[section] = kept
    review["impossible_infection_source_guardrail"] = {
        "enabled": True,
        "rule_file": str(pathogen_rule_lists.DEFAULT_IMPOSSIBLE_RULE_PATH),
        "removed_count": len([name for name in removed if name]),
        "removed_from_possible_or_watchlist": [name for name in removed if name],
        "policy_zh": "?? rules/impossible_infection_sources.json ?????????? possible_missed ? watchlist?",
    }
    return review


def apply_pulmonary_review_guardrails(review: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    """Constrain advisory LLM review output to pulmonary-relevant candidates."""
    evidence_index = build_queue_evidence_index(payload)
    anaerobe_moved, anaerobe_names = move_not_recommended_signal_to_watchlist(
        review,
        payload,
        predicate=has_strong_aspiration_anaerobe_watch_signal,
        guardrail_name="ASPIRATION_ANAEROBE_CLEAN_HIGH_SIGNAL_WATCHLIST",
        reason_zh=(
            "Lower-respiratory aspiration anaerobe/oral flora has clean high mNGS signal. "
            "Do not exclude only because another primary pathogen is more plausible; keep as non-primary watchlist."
        ),
    )
    water_moved, water_names = move_not_recommended_signal_to_watchlist(
        review,
        payload,
        predicate=has_minimal_environmental_watch_signal,
        guardrail_name="WATER_ENVIRONMENTAL_CLEAN_SIGNAL_WATCHLIST",
        reason_zh=(
            "Water/environmental or low-pulmonary-specificity organism has clean high-rank or moderate/high mNGS signal. "
            "Not enough for primary attribution, but keep as watchlist when device/water-source clues are unavailable."
        ),
    )
    removed_counts = {
        "gi_only_to_not_recommended": 0,
        "water_environment_possible_to_watchlist": 0,
        "water_environment_watchlist_to_not_recommended": 0,
        "aspiration_anaerobe_not_recommended_to_watchlist": anaerobe_moved,
        "water_environment_not_recommended_to_watchlist": water_moved,
    }

    for section in ("possible_missed_pathogens", "watchlist_candidates"):
        values = review.get(section) or []
        if not isinstance(values, list):
            review[section] = []
            continue
        kept: list[dict[str, Any]] = []
        demoted_to_watchlist: list[dict[str, Any]] = []
        for item in values:
            if not isinstance(item, dict):
                continue
            name = item.get("organism_name")
            evidence = evidence_index.get(normalized_name(name), {})
            if is_gi_only_pathogen(name, evidence):
                add_not_recommended(
                    review,
                    item,
                    "GI-only guardrail：目前證據只支持腸胃道檢出或腸道感染，沒有肺部、血液、胸腔或系統性感染支持，因此不列為肺部感染候選。",
                    {"postprocess_guardrail": "GI_ONLY_NOT_PULMONARY"},
                )
                removed_counts["gi_only_to_not_recommended"] += 1
                continue

            if weak_hospital_only_environmental_signal(name, evidence, payload):
                add_not_recommended(
                    review,
                    item,
                    "Low-specificity water/environmental guardrail: only non-sterile respiratory or old culture detected, no quantity/repeated culture, and mNGS=0; therefore it is not kept in watchlist. Reassess if BALF/ETA high quantity, blood/sterile-site evidence, effective single culture, repeated culture, or other clinical support appears.",
                    {"postprocess_guardrail": "WEAK_HOSPITAL_ONLY_ENVIRONMENTAL_NOT_RECOMMENDED"},
                )
                removed_counts.setdefault("weak_hospital_only_environment_to_not_recommended", 0)
                removed_counts["weak_hospital_only_environment_to_not_recommended"] += 1
                continue

            if is_water_or_low_specificity_environmental_name(name):
                strong_signal = has_strong_environmental_signal(evidence)
                minimal_signal = has_minimal_environmental_watch_signal(evidence)
                if section == "possible_missed_pathogens" and not strong_signal:
                    item = dict(item)
                    item["postprocess_guardrail"] = "WATER_ENVIRONMENTAL_GNB_DEMOTED_TO_WATCHLIST"
                    item["why_not_auto_upgrade"] = (
                        str(item.get("why_not_auto_upgrade") or "").strip()
                        + " 水源/環境相關菌若缺乏同菌醫院端支持、D2/D3 dominance 或明確強訊號，最多列為 watchlist。"
                    ).strip()
                    demoted_to_watchlist.append(item)
                    removed_counts["water_environment_possible_to_watchlist"] += 1
                    continue
                if section == "watchlist_candidates" and not minimal_signal:
                    add_not_recommended(
                        review,
                        item,
                        "水源/環境或低肺炎特異性菌種 guardrail：目前缺乏同菌醫院端支持、dominance 或足夠 mNGS 訊號，較可能是背景、污染或定殖，不列入 watchlist。",
                        {"postprocess_guardrail": "WATER_ENVIRONMENTAL_WEAK_SIGNAL_NOT_RECOMMENDED"},
                    )
                    removed_counts["water_environment_watchlist_to_not_recommended"] += 1
                    continue

            kept.append(item)
        review[section] = kept
        if section == "possible_missed_pathogens" and demoted_to_watchlist:
            existing = review.get("watchlist_candidates")
            if not isinstance(existing, list):
                existing = []
                review["watchlist_candidates"] = existing
            existing.extend(demoted_to_watchlist)

    review["pulmonary_review_guardrails"] = {
        "enabled": True,
        "policy": (
            "GI-only 病原不列入肺部感染候選；水源/環境或低肺炎特異性 GNB 若缺乏同菌支持、"
            "dominance 或強 mNGS 訊號，possible missed 會降到 watchlist，弱訊號 watchlist 會移到 not_recommended。"
        ),
        "aspiration_anaerobe_recovered_to_watchlist": anaerobe_names,
        "water_environment_recovered_to_watchlist": water_names,
        **removed_counts,
    }
    return review


def review_item_score(item: dict[str, Any], section: str, evidence_index: dict[str, dict[str, Any]]) -> float:
    evidence = evidence_index.get(normalized_name(item.get("organism_name")), {})
    rank = numeric_value(evidence.get("rank_priority"), default=99.0)
    reads = numeric_value(evidence.get("reads"))
    percentile = numeric_value(evidence.get("reads_percentile"))
    section_bonus = 1000.0 if section == "possible_missed_pathogens" else 0.0
    return (
        section_bonus
        + confidence_score(item.get("confidence"))
        + level_score(evidence.get("integrated_causative_level"))
        + tier_score(evidence.get("mngs_signal_tier"))
        + tier_score(evidence.get("reads_tier"))
        + dominance_score(evidence.get("dominance_tier"))
        + support_score(evidence.get("support_modules"))
        + min(reads, 100000.0) / 1000.0
        + percentile * 100.0
        + max(0.0, 80.0 - rank * 10.0)
    )


def consolidate_same_genus_review_items(review: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    """Keep one most plausible species per genus across possible/watchlist review output."""
    evidence_index = build_queue_evidence_index(payload)
    tracked: list[dict[str, Any]] = []
    for section in ("possible_missed_pathogens", "watchlist_candidates"):
        values = review.get(section) or []
        if not isinstance(values, list):
            review[section] = []
            continue
        for index, item in enumerate(values):
            if not isinstance(item, dict):
                continue
            genus = genus_key(item.get("organism_name"))
            tracked.append(
                {
                    "section": section,
                    "index": index,
                    "item": item,
                    "genus": genus,
                    "representative_group": pathogen_names.representative_group_key(item.get("organism_name")),
                    "score": review_item_score(item, section, evidence_index),
                }
            )

    winners: dict[str, dict[str, Any]] = {}
    removed: list[dict[str, Any]] = []
    grouped_by_representative: dict[str, list[dict[str, Any]]] = {}
    for entry in tracked:
        group = entry.get("representative_group")
        if not group or str(group).startswith("name:"):
            continue
        grouped_by_representative.setdefault(str(group), []).append(entry)

    for group, entries in grouped_by_representative.items():
        suppressible = [
            entry
            for entry in entries
            if pathogen_names.representative_group_suppression_allowed(
                (entry.get("item") or {}).get("organism_name")
            )
        ]
        if not suppressible:
            continue
        retained_pool = [
            entry
            for entry in entries
            if not pathogen_names.representative_group_suppression_allowed(
                (entry.get("item") or {}).get("organism_name")
            )
        ] or entries
        winner = max(retained_pool, key=lambda entry: entry["score"])
        winners[group] = winner
        for entry in suppressible:
            if entry is not winner:
                removed.append(entry)

    removed_ids = {(entry["section"], entry["index"]) for entry in removed}
    for section in ("possible_missed_pathogens", "watchlist_candidates"):
        values = review.get(section) or []
        review[section] = [
            item
            for index, item in enumerate(values)
            if isinstance(item, dict) and (section, index) not in removed_ids
        ]

    if removed:
        not_recommended = review.setdefault("not_recommended_for_upgrade", [])
        if not isinstance(not_recommended, list):
            not_recommended = []
            review["not_recommended_for_upgrade"] = not_recommended
        for entry in removed:
            winner = winners.get(entry.get("representative_group"))
            winner_name = ((winner or {}).get("item") or {}).get("organism_name") or "同屬較強代表"
            organism_name = (entry["item"] or {}).get("organism_name")
            not_recommended.append(
                {
                    "organism_name": organism_name,
                    "rationale_zh": (
                        f"同屬合併後未保留；同一屬已有較可能代表菌種「{winner_name}」。"
                        "此菌保留為同屬背景/輔助訊號，不另外列為 possible missed 或 watchlist。"
                    ),
                    "same_genus_representative": winner_name,
                }
            )
        review["same_genus_consolidation"] = {
            "enabled": True,
            "policy": "possible_missed_pathogens 與 watchlist_candidates 內同屬候選只保留證據分數最高的一隻，不改成 spp.",
            "removed_count": len(removed),
        }
    else:
        review.setdefault(
            "same_genus_consolidation",
            {
                "enabled": True,
                "policy": "possible_missed_pathogens 與 watchlist_candidates 內同屬候選只保留證據分數最高的一隻，不改成 spp.",
                "removed_count": 0,
            },
        )
    return review


def consolidate_oral_aspiration_flora_pattern(review: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    """Represent oral/aspiration flora as a pattern plus 1-2 representative organisms."""
    evidence_index = build_queue_evidence_index(payload)
    tracked: list[dict[str, Any]] = []
    for section in ("possible_missed_pathogens", "watchlist_candidates"):
        values = review.get(section) or []
        if not isinstance(values, list):
            review[section] = []
            continue
        for index, item in enumerate(values):
            if not isinstance(item, dict):
                continue
            name = item.get("organism_name")
            if not is_oral_aspiration_flora_name(name):
                continue
            tracked.append(
                {
                    "section": section,
                    "index": index,
                    "item": item,
                    "score": review_item_score(item, section, evidence_index),
                    "evidence": evidence_index.get(normalized_name(name), {}),
                }
            )

    policy = "口腔/aspiration flora 不逐隻全部列為候選；同時出現多隻時，只保留 1-2 隻代表菌，其餘作為群聚背景訊號。"
    if len(tracked) <= 2:
        if tracked:
            for entry in tracked:
                entry["item"]["oral_aspiration_flora_pattern_representative"] = True
        review["oral_aspiration_flora_pattern"] = {
            "enabled": True,
            "pattern_detected": bool(tracked),
            "policy": policy,
            "representative_organisms": [entry["item"].get("organism_name") for entry in tracked],
            "grouped_non_representative_organisms": [],
            "interpretation_zh": "偵測到口腔/吸入相關菌群訊號；目前候選數不多，保留代表菌即可。" if tracked else "",
        }
        return review

    ranked = sorted(tracked, key=lambda entry: entry["score"], reverse=True)
    keep = ranked[:2]
    remove = ranked[2:]
    keep_ids = {(entry["section"], entry["index"]) for entry in keep}
    remove_ids = {(entry["section"], entry["index"]) for entry in remove}

    for section in ("possible_missed_pathogens", "watchlist_candidates"):
        values = review.get(section) or []
        new_values: list[dict[str, Any]] = []
        for index, item in enumerate(values):
            if not isinstance(item, dict):
                continue
            item_id = (section, index)
            if item_id in remove_ids:
                continue
            if item_id in keep_ids:
                item = dict(item)
                item["oral_aspiration_flora_pattern_representative"] = True
            new_values.append(item)
        review[section] = new_values

    representatives = [entry["item"].get("organism_name") for entry in keep]
    grouped = []
    for entry in remove:
        evidence = entry.get("evidence") or {}
        grouped.append(
            {
                "organism_name": entry["item"].get("organism_name"),
                "original_section": entry["section"],
                "rank_priority": evidence.get("rank_priority"),
                "reads": evidence.get("reads"),
                "reads_tier": evidence.get("reads_tier"),
                "dominance_tier": evidence.get("dominance_tier"),
                "reason_zh": "同屬於口腔/吸入相關菌群，為避免逐隻列出造成臨床解讀過度複雜，收進 oral/aspiration flora pattern。",
            }
        )

    review["oral_aspiration_flora_pattern"] = {
        "enabled": True,
        "pattern_detected": True,
        "policy": policy,
        "representative_organisms": representatives,
        "grouped_non_representative_organisms": grouped,
        "interpretation_zh": (
            "偵測到多個口腔/吸入相關菌群。臨床上較適合解讀為 oral/aspiration flora pattern，"
            "代表可能有吸入、口腔菌群混入或多菌型下呼吸道訊號；是否需要 anaerobic coverage "
            "仍需依 aspiration risk、肺膿瘍、empyema、壞死性肺炎或培養證據判斷。"
        ),
    }
    return review


def active_direct_review_keys(review: dict[str, Any]) -> set[str]:
    keys: set[str] = set()
    for tier in DIRECT_ACTIVE_TIERS:
        values = review.get(tier) or []
        if not isinstance(values, list):
            continue
        for item in values:
            if not isinstance(item, dict):
                continue
            key = canonical_name_key(item.get("organism_name"))
            if key:
                keys.add(key)
    return keys


def remove_direct_review_key_from_audit_tiers(review: dict[str, Any], key: str) -> None:
    for tier in ("review_low_specificity", "omitted_with_reason"):
        values = review.get(tier) or []
        if not isinstance(values, list):
            continue
        review[tier] = [
            item
            for item in values
            if not isinstance(item, dict) or canonical_name_key(item.get("organism_name")) != key
        ]


def apply_strict_aspiration_anaerobe_rescue(review: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    evidence_index = build_queue_evidence_index(payload)
    picked_keys, _ = picked_name_sets_from_payload(payload)
    active_keys = active_direct_review_keys(review)
    rescued: list[dict[str, Any]] = []
    context_items = review.get("review_context_needed")
    if not isinstance(context_items, list):
        context_items = []
        review["review_context_needed"] = context_items

    for evidence in evidence_index.values():
        if not isinstance(evidence, dict):
            continue
        name = clean_review_organism_name(evidence.get("organism_name"))
        key = canonical_name_key(name)
        if not key or key in picked_keys or key in active_keys:
            continue
        if not has_strong_aspiration_anaerobe_watch_signal(evidence, name):
            continue
        item = ensure_review_item_fields(
            {"organism_name": name, "confidence": "low"},
            tier="review_context_needed",
            reason=(
                "strict aspiration anaerobe rescue: lower-respiratory clean mNGS signal is strong enough "
                "for RAG-visible syndrome/context review even if model grouped or demoted it"
            ),
            source_section="postprocess_strict_aspiration_anaerobe_rescue",
            evidence=evidence,
        )
        item["postprocess_guardrail"] = "STRICT_ASPIRATION_ANAEROBE_RESCUE"
        remove_direct_review_key_from_audit_tiers(review, key)
        context_items.append(item)
        active_keys.add(key)
        rescued.append(
            {
                "organism_name": item.get("organism_name"),
                "rank_priority": evidence.get("rank_priority"),
                "reads": evidence.get("reads"),
                "reads_tier": evidence.get("reads_tier"),
                "reads_percentile": evidence.get("reads_percentile"),
                "dominance_tier": evidence.get("dominance_tier"),
            }
        )

    review["strict_aspiration_anaerobe_rescue"] = {
        "enabled": True,
        "policy": (
            "After direct oral-flora consolidation, restore strict aspiration anaerobes with lower-respiratory "
            "rank<=3 and R3/R4 or high-percentile clean mNGS signal to review_context_needed."
        ),
        "rescued_count": len(rescued),
        "rescued_items": rescued,
    }
    return review


def process_patient(
    patient_dir: Path,
    *,
    queue_suffix: str,
    deterministic_max_suffix: str,
    mngs_to_specimen_suffix: str,
    output_suffix: str,
    artifact_root: Path | None,
    summary_mode: str,
    model: str,
    temperature: float | None,
    reasoning_effort: str | None,
    overwrite: bool,
    skip_existing: bool,
    dry_run: bool,
    keep_prompt: bool,
    keep_raw: bool,
    reuse_existing_raw: bool,
) -> tuple[Path, Path]:
    summary_dir = summary_dir_for(patient_dir)
    queue_path = path_for_suffix(patient_dir, queue_suffix, artifact_root=artifact_root)
    max_path = path_for_suffix(patient_dir, deterministic_max_suffix)
    summary_path = mngs.find_final_summary_file(patient_dir, summary_mode=summary_mode)
    if summary_path is None:
        raise FileNotFoundError(f"final summary not found for {patient_dir}")
    mngs_to_specimen_path = path_for_suffix(patient_dir, mngs_to_specimen_suffix)
    if not mngs_to_specimen_path.exists():
        mngs_to_specimen_path = None

    output_path = path_for_suffix(patient_dir, output_suffix, artifact_root=artifact_root)
    prompt_base_dir = output_path.parent if artifact_root is not None else summary_dir
    prompt_path = prompt_base_dir / f"{mngs.extract_patient_identifier(patient_dir)}_{output_suffix}.prompt.md"
    raw_output_path = prompt_base_dir / f"{mngs.extract_patient_identifier(patient_dir)}_{output_suffix}.raw.txt"

    if output_path.exists() and skip_existing:
        return prompt_path, output_path
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {output_path}. Use --overwrite or --skip-existing.")
    if not queue_path.exists():
        raise FileNotFoundError(queue_path)
    if not max_path.exists():
        raise FileNotFoundError(max_path)

    payload = build_review_payload(
        patient_dir=patient_dir,
        queue_path=queue_path,
        max_path=max_path,
        summary_path=summary_path,
        mngs_to_specimen_path=mngs_to_specimen_path,
    )
    prompt = build_prompt(payload)
    if dry_run or keep_prompt:
        write_text(prompt_path, prompt)

    if dry_run:
        dry_payload = {
            "dry_run": True,
            "message": "Prompt generated; LLM was not called.",
            "prompt_path": str(prompt_path),
            "review_payload_summary": {
                "patient_id": payload.get("patient_id"),
                "review_queue_count": payload.get("review_queue", {}).get("review_queue_count"),
                "formal_picked_count": payload.get("review_queue", {}).get("formal_picked_count"),
            },
        }
        write_json(output_path, dry_payload)
        return prompt_path, output_path

    if reuse_existing_raw:
        if not raw_output_path.exists():
            raise FileNotFoundError(f"Raw LLM response not found: {raw_output_path}")
        raw_text = raw_output_path.read_text(encoding="utf-8")
        usage = None
    else:
        raw_text, usage = send_to_llm(
            prompt,
            model=model,
            temperature=temperature,
            reasoning_effort=reasoning_effort,
        )
        if keep_raw:
            write_text(raw_output_path, raw_text)
    parsed = parse_json_response(raw_text)
    if not isinstance(parsed, dict):
        raise RuntimeError("LLM output JSON is not an object.")
    parsed.setdefault("review_purpose", "safety_review_only")
    parsed["source_files"] = payload.get("source_files")
    parsed["llm_request"] = {
        "model": model,
        "temperature_requested": temperature,
        "temperature_sent": temperature if model_accepts_temperature(model) else None,
        "reasoning_effort": reasoning_effort,
    }
    parsed["usage"] = usage
    parsed["reused_existing_raw"] = reuse_existing_raw
    parsed = normalize_direct_review_output(parsed, payload)
    if keep_raw:
        parsed["raw_output_path"] = str(raw_output_path)
    if keep_prompt:
        parsed["prompt_path"] = str(prompt_path)
    write_json(output_path, parsed)
    return prompt_path, output_path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run LLM safety review for deterministic mNGS missed-candidate queues.")
    parser.add_argument("patient_root", type=Path)
    parser.add_argument("--patients", nargs="*")
    parser.add_argument("--queue-suffix", default=DEFAULT_QUEUE_SUFFIX)
    parser.add_argument("--deterministic-max-suffix", default=DEFAULT_MAX_SUFFIX)
    parser.add_argument("--mngs-to-specimen-suffix", default=DEFAULT_MNGS_TO_SPECIMEN_SUFFIX)
    parser.add_argument("--output-suffix", default=DEFAULT_OUTPUT_SUFFIX)
    parser.add_argument(
        "--artifact-root",
        type=Path,
        help="Root where queue/review/prompt artifacts are stored. Reads formal patient inputs from patient_root.",
    )
    parser.add_argument("--summary-mode", choices=["deterministic", "full"], default="deterministic")
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--temperature", type=float, default=None)
    parser.add_argument("--reasoning-effort", choices=["none", "low", "medium", "high", "xhigh", "max"], default=None)
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--keep-prompt", action="store_true", help="Keep generated prompt markdown files.")
    parser.add_argument("--keep-raw", action="store_true", help="Keep raw LLM text responses.")
    parser.add_argument(
        "--reuse-existing-raw",
        action="store_true",
        help="Re-parse and postprocess an existing raw response without calling the model again.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    patient_dirs = mngs.collect_patient_directories(
        [args.patient_root],
        summary_mode=args.summary_mode,
        mngs_source="ranked_only",
    )
    requested = parse_patient_ids(args.patients)
    if requested:
        patient_dirs = [
            path for path in patient_dirs if mngs.extract_patient_identifier(path) in requested
        ]

    written = 0
    failures: list[str] = []

    def run_one(path: Path) -> tuple[Path, Path]:
        return process_patient(
            path,
            queue_suffix=args.queue_suffix,
            deterministic_max_suffix=args.deterministic_max_suffix,
            mngs_to_specimen_suffix=args.mngs_to_specimen_suffix,
            output_suffix=args.output_suffix,
            artifact_root=args.artifact_root,
            summary_mode=args.summary_mode,
            model=args.model,
            temperature=args.temperature,
            reasoning_effort=args.reasoning_effort,
            overwrite=args.overwrite,
            skip_existing=args.skip_existing,
            dry_run=args.dry_run,
            keep_prompt=args.keep_prompt,
            keep_raw=args.keep_raw,
            reuse_existing_raw=args.reuse_existing_raw,
        )

    if args.max_workers <= 1:
        for patient_dir in patient_dirs:
            try:
                prompt_path, output_path = run_one(patient_dir)
                print(f"wrote prompt {prompt_path}")
                print(f"wrote review {output_path}")
                written += 1
            except Exception as exc:  # noqa: BLE001
                failures.append(f"{patient_dir}: {exc}")
                print(f"failed {patient_dir}: {exc}", file=sys.stderr)
    else:
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            future_map = {executor.submit(run_one, patient_dir): patient_dir for patient_dir in patient_dirs}
            for future in as_completed(future_map):
                patient_dir = future_map[future]
                try:
                    prompt_path, output_path = future.result()
                    print(f"wrote prompt {prompt_path}")
                    print(f"wrote review {output_path}")
                    written += 1
                except Exception as exc:  # noqa: BLE001
                    failures.append(f"{patient_dir}: {exc}")
                    print(f"failed {patient_dir}: {exc}", file=sys.stderr)

    print(f"written_count={written}")
    if failures:
        print("failures:")
        for failure in failures:
            print(f"- {failure}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = next(
    (
        path
        for path in (ROOT, ROOT.parent, ROOT / ".artifact_work" / "multimodal-publish-20261002")
        if (path / "RAG_re_casefit").is_dir() and (path / "RAG_re_clinical").is_dir()
    ),
    ROOT / ".artifact_work" / "multimodal-publish-20261002",
)
CASEFIT_ROOT = ARCHIVE_ROOT / "RAG_re_casefit"
CLINICAL_ROOT = ARCHIVE_ROOT / "RAG_re_clinical"
HISTORICAL_CASEFIT_DIR = ROOT / "full_high_context_casefit_v1_20260819"
DEFAULT_OUTPUT_DIR = (
    ROOT / "outputs" / "runs" / "2026-10-13_73patients_historical_ober_r5_replay_v1"
)

for package_root in (CASEFIT_ROOT, CLINICAL_ROOT):
    package_text = str(package_root)
    if package_text not in sys.path:
        sys.path.insert(0, package_text)

from rag_re_casefit.case_card import build_case_card, prompt_case_card  # noqa: E402
from rag_re_casefit.rules import (  # noqa: E402
    aggregate_judgments,
    apply_traceability_gate,
    derive_article_verdict,
)
from rag_re_clinical.clinical_rules import evaluate_modules  # noqa: E402
from tools.adjudicate_rag_v2_packets import (  # noqa: E402
    codex_cli_environment,
    ensure_codex_cli_authenticated,
)


SYSTEM_PROMPT_PATH = CASEFIT_ROOT / "prompts" / "article_casefit_v1.txt"
REPAIR_PROMPT_PATH = CASEFIT_ROOT / "prompts" / "evidence_span_repair_v1.txt"
SYSTEM_PROMPT = SYSTEM_PROMPT_PATH.read_text(encoding="utf-8").strip()
REPAIR_SYSTEM_PROMPT = REPAIR_PROMPT_PATH.read_text(encoding="utf-8").strip()
SYSTEM_PROMPT_SHA256 = hashlib.sha256(SYSTEM_PROMPT.encode("utf-8")).hexdigest()
REPAIR_PROMPT_SHA256 = hashlib.sha256(REPAIR_SYSTEM_PROMPT.encode("utf-8")).hexdigest()
MODEL_DEFAULT = "gpt-5.6-luna"
REASONING_DEFAULT = "low"

ALIGNMENT_SCHEMA = {"type": "string", "enum": ["match", "partial", "mismatch", "unknown"]}
EXTRACTION_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "organism_scope": {
            "type": "string",
            "enum": ["exact_species", "genus_only", "other", "unclear"],
        },
        "human_clinical_context": {
            "type": "string",
            "enum": [
                "causal_infection",
                "colonization_or_contamination",
                "reactivation_or_bystander",
                "detection_only",
                "nonhuman",
                "unclear",
            ],
        },
        "site_alignment": ALIGNMENT_SCHEMA,
        "syndrome_alignment": ALIGNMENT_SCHEMA,
        "specimen_alignment": ALIGNMENT_SCHEMA,
        "host_alignment": ALIGNMENT_SCHEMA,
        "phenotype_alignment": ALIGNMENT_SCHEMA,
        "temporal_alignment": ALIGNMENT_SCHEMA,
        "support_span": {"type": "string"},
        "counter_span": {"type": "string"},
        "rationale": {"type": "string"},
    },
    "required": [
        "organism_scope",
        "human_clinical_context",
        "site_alignment",
        "syndrome_alignment",
        "specimen_alignment",
        "host_alignment",
        "phenotype_alignment",
        "temporal_alignment",
        "support_span",
        "counter_span",
        "rationale",
    ],
}
REPAIR_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "properties": {
        "support_span": {"type": "string"},
        "counter_span": {"type": "string"},
    },
    "required": ["support_span", "counter_span"],
}
EXTRACTION_KEYS = (
    "organism_scope",
    "human_clinical_context",
    "site_alignment",
    "syndrome_alignment",
    "specimen_alignment",
    "host_alignment",
    "phenotype_alignment",
    "temporal_alignment",
)
AGGREGATION_CONFIG = {
    "min_judgeable_articles": 5,
    "min_strong_for_high": 2,
    "min_strong_for_moderate": 1,
    "min_partial_for_moderate": 3,
}


@dataclass(frozen=True)
class CohortSpec:
    label: str
    packets: Path
    manifest: Path


@dataclass
class CandidateCase:
    cohort: str
    case_id: str
    patient_id: str
    organism: str
    review_tier: str
    packet: dict[str, Any]
    merged_path: Path
    merged: dict[str, Any]
    deterministic: dict[str, Any] | None
    review_row: dict[str, Any]
    case_card: dict[str, Any]
    articles: list[dict[str, Any]]
    historical_artifact: dict[str, Any] | None
    artifact_path: Path


COHORTS = (
    CohortSpec(
        "KH33",
        ROOT / "outputs" / "runs" / "2026-10-12_KH33_legacy_ober_packets_two_stage_v1" / "adjudication_packets.jsonl",
        ROOT / "outputs" / "runs" / "2026-10-12_KH33_legacy_ober_production_two_stage_v1" / "production_manifest.json",
    ),
    CohortSpec(
        "KH18",
        ROOT / "outputs" / "runs" / "2026-10-12_KH18_legacy_ober_packets_v1" / "adjudication_packets.jsonl",
        ROOT / "outputs" / "runs" / "2026-10-12_KH18_legacy_ober_production_v1" / "production_manifest.json",
    ),
    CohortSpec(
        "A11",
        ROOT / "outputs" / "runs" / "2026-10-08_two_hospital_A11_legacy_ober_packets_v1" / "adjudication_packets.jsonl",
        ROOT / "outputs" / "runs" / "2026-10-08_two_hospital_A11_legacy_ober_production_v1" / "production_manifest.json",
    ),
    CohortSpec(
        "B11",
        ROOT / "outputs" / "runs" / "2026-10-12_two_hospital_B11_legacy_ober_packets_v1" / "adjudication_packets.jsonl",
        ROOT / "outputs" / "runs" / "2026-10-12_two_hospital_B11_legacy_ober_production_v1" / "production_manifest.json",
    ),
)


def compact(value: Any) -> str:
    return " ".join(str(value or "").split())


def normal(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", compact(value).casefold())


def stable_hash(value: Any) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")
    temporary.replace(path)


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line.strip():
            continue
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError(f"Expected an object at {path}:{line_number}")
        rows.append(value)
    return rows


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "item"


def candidate_names(row: dict[str, Any]) -> set[str]:
    values = {
        normal(row.get(key))
        for key in (
            "canonical_organism_name",
            "organism_name",
            "organism",
            "pathogen_name",
            "pathogen",
            "name",
        )
        if compact(row.get(key))
    }
    taxonomy = row.get("taxonomy_profile")
    if isinstance(taxonomy, dict):
        for key in ("canonical_key", "display_name", "input_name"):
            if compact(taxonomy.get(key)):
                values.add(normal(taxonomy.get(key)))
    return values


def packet_names(packet: dict[str, Any]) -> set[str]:
    organism = packet.get("organism") or {}
    values = {
        normal(organism.get(key))
        for key in ("display_name", "canonical_key")
        if compact(organism.get(key))
    }
    profile = organism.get("taxonomy_profile")
    if isinstance(profile, dict):
        for key in ("canonical_key", "display_name", "input_name"):
            if compact(profile.get(key)):
                values.add(normal(profile.get(key)))
    return values


def select_matching_row(rows: list[dict[str, Any]], packet: dict[str, Any]) -> dict[str, Any] | None:
    targets = packet_names(packet)
    matches = [row for row in rows if candidate_names(row).intersection(targets)]
    if len(matches) == 1:
        return matches[0]
    display = compact((packet.get("organism") or {}).get("display_name"))
    exact = [row for row in matches if normal(row.get("organism_name")) == normal(display)]
    if len(exact) == 1:
        return exact[0]
    if not matches:
        return None
    raise ValueError(f"Ambiguous candidate match for {packet.get('case_id')}: {len(matches)} rows")


def merged_path_for(manifest: dict[str, Any], patient_id: str) -> Path:
    patient_root = Path(str(manifest["patient_root"]))
    if not patient_root.is_absolute():
        patient_root = ROOT / patient_root
    suffix = compact(manifest["merged_suffix"])
    return (
        patient_root
        / f"NGS_patient_{patient_id}_json"
        / "summary_outputs"
        / f"NGS_patient_{patient_id}_{suffix}.json"
    )


def historical_artifacts(directory: Path = HISTORICAL_CASEFIT_DIR) -> dict[tuple[str, str], dict[str, Any]]:
    result: dict[tuple[str, str], dict[str, Any]] = {}
    if not directory.is_dir():
        return result
    for path in sorted(directory.glob("*_casefit.json")):
        artifact = read_json(path)
        key = (
            compact(artifact.get("patient_id")),
            normal(artifact.get("canonical_organism") or artifact.get("candidate_organism")),
        )
        artifact["_historical_source_path"] = str(path.resolve())
        result[key] = artifact
    return result


def load_cohort_specs(path: Path) -> tuple[CohortSpec, ...]:
    config = read_json(path)
    values = config.get("cohorts")
    if not isinstance(values, list) or not values:
        raise ValueError("Cohort config requires a non-empty cohorts list")
    specs: list[CohortSpec] = []
    seen: set[str] = set()
    for row in values:
        if not isinstance(row, dict):
            raise ValueError("Each cohort config entry must be an object")
        label = compact(row.get("label"))
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", label) or label in seen:
            raise ValueError(f"Invalid or duplicate cohort label: {label!r}")
        seen.add(label)
        paths = []
        for key in ("packets", "manifest"):
            value = compact(row.get(key))
            if not value:
                raise ValueError(f"Missing {key} for cohort {label}")
            target = Path(value)
            paths.append(target if target.is_absolute() else path.resolve().parent / target)
        specs.append(CohortSpec(label, paths[0], paths[1]))
    return tuple(specs)


def prepare_cases(
    output_dir: Path,
    selected_cohorts: set[str],
    cohort_specs: tuple[CohortSpec, ...] = COHORTS,
    historical_directory: Path = HISTORICAL_CASEFIT_DIR,
) -> list[CandidateCase]:
    historical = historical_artifacts(historical_directory)
    cases: list[CandidateCase] = []
    seen: set[tuple[str, str]] = set()
    for spec in cohort_specs:
        if spec.label not in selected_cohorts:
            continue
        if not spec.packets.is_file() or not spec.manifest.is_file():
            raise FileNotFoundError(f"Missing {spec.label} packet or manifest input")
        manifest = read_json(spec.manifest)
        for packet in read_jsonl(spec.packets):
            patient_id = compact(packet.get("patient_id"))
            case_id = compact(packet.get("case_id"))
            organism = compact((packet.get("organism") or {}).get("display_name"))
            review_tier = compact((packet.get("current_model_state") or {}).get("source_tier"))
            if not patient_id or not case_id or not organism:
                raise ValueError(f"Incomplete packet identity in {spec.packets}")
            if review_tier not in {"review_high_priority", "review_context_needed"}:
                raise ValueError(f"Unsupported review tier for {case_id}: {review_tier}")
            cohort_key = (spec.label, case_id)
            if cohort_key in seen:
                raise ValueError(f"Duplicate case id in {spec.label}: {case_id}")
            seen.add(cohort_key)
            merged_path = merged_path_for(manifest, patient_id)
            if not merged_path.is_file():
                raise FileNotFoundError(merged_path)
            merged_bytes = merged_path.read_bytes()
            merged = json.loads(merged_bytes.decode("utf-8-sig"))
            deterministic_root = merged.get("deterministic_max") or {}
            deterministic_rows = deterministic_root.get("candidates")
            if not isinstance(deterministic_rows, list):
                deterministic_rows = deterministic_root.get("pathogen_candidates") or []
            deterministic = select_matching_row(
                [row for row in deterministic_rows if isinstance(row, dict)], packet
            )
            review = merged.get("llm_missed_candidate_review") or {}
            review_rows: list[dict[str, Any]] = []
            for tier in ("review_high_priority", "review_context_needed"):
                values = review.get(tier) or []
                if isinstance(values, list):
                    review_rows.extend(row for row in values if isinstance(row, dict))
            review_row = select_matching_row(review_rows, packet)
            if review_row is None:
                raise ValueError(f"No review row for {spec.label} {case_id}")
            card_root = copy.deepcopy(merged)
            card_bundle = card_root.setdefault("deterministic_max", {})
            case_context = (packet.get("local_evidence") or {}).get("case_context") or {}
            if not compact(card_bundle.get("dominant_source")):
                card_bundle["dominant_source"] = compact(case_context.get("dominant_source"))
            effective_for_card = deterministic or dict(review_row.get("evidence_snapshot") or {})
            effective_for_card = dict(effective_for_card)
            effective_for_card.setdefault("organism_name", organism)
            case_card = build_case_card(card_root, effective_for_card, source_bytes=merged_bytes)
            literature = packet.get("literature_evidence") or {}
            articles = literature.get("articles") or []
            if not isinstance(articles, list) or any(not isinstance(row, dict) for row in articles):
                raise ValueError(f"Invalid article list for {spec.label} {case_id}")
            historical_match = None
            if spec.label == "KH33":
                historical_match = historical.get((patient_id, normal(organism)))
            artifact_path = output_dir / spec.label / "casefit" / f"{safe_name(case_id)}_casefit.json"
            cases.append(
                CandidateCase(
                    cohort=spec.label,
                    case_id=case_id,
                    patient_id=patient_id,
                    organism=organism,
                    review_tier=review_tier,
                    packet=packet,
                    merged_path=merged_path,
                    merged=merged,
                    deterministic=deterministic,
                    review_row=review_row,
                    case_card=case_card,
                    articles=articles,
                    historical_artifact=historical_match,
                    artifact_path=artifact_path,
                )
            )
    return cases


def no_abstract_judgment(
    organism: str, card: dict[str, Any], article: dict[str, Any], model: str, reasoning: str
) -> dict[str, Any]:
    extraction = {
        "organism_scope": "unclear",
        "human_clinical_context": "unclear",
        **{key: "unknown" for key in EXTRACTION_KEYS[2:]},
    }
    decision = derive_article_verdict(extraction)
    return {
        "schema_version": "rag_re_casefit.article_judgment.v1",
        "pmid": compact(article.get("pmid")),
        "article_url": article.get("url"),
        "organism": organism,
        "case_card_hash": stable_hash(prompt_case_card(card)),
        "model": model,
        "reasoning_effort": reasoning,
        "prompt_sha256": SYSTEM_PROMPT_SHA256,
        "status": "no_abstract",
        **extraction,
        "support_span": "",
        "support_span_valid": False,
        "counter_span": "",
        "counter_span_valid": False,
        "span_repair_attempted": False,
        "span_repair_method": "not_applicable",
        "span_repair_prompt_sha256": REPAIR_PROMPT_SHA256,
        "rationale": "No abstract was available.",
        **decision,
        "traceability_status": "valid",
        "cache_hit": False,
        "usage": None,
    }


def span(value: Any, article_blob: str) -> tuple[str, bool]:
    text = compact(value)
    if len(text.split()) > 20:
        return "", False
    valid = bool(text and text.casefold() in article_blob.casefold())
    return (text if valid else ""), valid


def trigger_window(title: str, abstract: str, terms: list[str]) -> tuple[str, bool]:
    for source in (abstract, title):
        tokens = list(re.finditer(r"\S+", source))
        for index, token in enumerate(tokens):
            normalized = re.sub(r"[^a-z]+", "", token.group(0).casefold())
            if not any(term in normalized for term in terms):
                continue
            start = max(0, index - 5)
            end = min(len(tokens), index + 10)
            if end - start > 20:
                end = start + 20
            excerpt = source[tokens[start].start() : tokens[end - 1].end()]
            return span(excerpt, f"{title} {abstract}")
    return "", False


def cache_key(
    organism: str,
    card: dict[str, Any],
    article: dict[str, Any],
    model: str,
    reasoning: str,
) -> str:
    return stable_hash(
        {
            "prompt_sha256": SYSTEM_PROMPT_SHA256,
            "model": model,
            "reasoning_effort": reasoning,
            "organism": organism,
            "case_card": prompt_case_card(card),
            "pmid": article.get("pmid"),
            "title": article.get("title"),
            "abstract": article.get("abstract"),
        }
    )


def call_codex_json(
    *,
    executable: str,
    model: str,
    reasoning: str,
    system_prompt: str,
    payload: dict[str, Any],
    schema: dict[str, Any],
    runtime_dir: Path,
    result_key: str,
    timeout_seconds: int,
) -> dict[str, Any]:
    schema_dir = runtime_dir / "schemas"
    result_dir = runtime_dir / "results"
    schema_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    schema_name = "repair.schema.json" if schema is REPAIR_SCHEMA else "extraction.schema.json"
    schema_path = schema_dir / schema_name
    if not schema_path.is_file():
        write_json(schema_path, schema)
    result_path = result_dir / f"{safe_name(result_key)}.{threading.get_ident()}.json"
    result_path.unlink(missing_ok=True)
    prompt = (
        "Do not browse, use tools, inspect files, or modify files. Use only the supplied policy "
        "and payload. Return only the JSON object required by the output schema.\n\n"
        f"SYSTEM POLICY:\n{system_prompt}\n\n"
        f"INPUT PAYLOAD:\n{json.dumps(payload, ensure_ascii=False)}"
    )
    command = [
        executable,
        "exec",
        "--model",
        model,
        "--sandbox",
        "read-only",
        "--ephemeral",
        "--ignore-rules",
        "--skip-git-repo-check",
        "--output-schema",
        str(schema_path.resolve()),
        "--output-last-message",
        str(result_path.resolve()),
        "--cd",
        str(runtime_dir.resolve()),
        "--config",
        f'model_reasoning_effort="{reasoning}"',
        "-",
    ]
    completed = subprocess.run(
        command,
        input=prompt,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        env=codex_cli_environment(),
        timeout=timeout_seconds,
        check=False,
    )
    if completed.returncode != 0:
        diagnostic = "\n".join((completed.stderr or completed.stdout or "").splitlines()[-20:])
        raise RuntimeError(f"Codex CLI failed ({completed.returncode}): {diagnostic}")
    if not result_path.is_file():
        raise RuntimeError("Codex CLI completed without a structured response")
    raw = result_path.read_text(encoding="utf-8").strip()
    result_path.unlink(missing_ok=True)
    if not raw:
        raise RuntimeError("Codex CLI wrote an empty response")
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError("Codex CLI response is not a JSON object")
    return value


def judge_article(
    *,
    executable: str,
    organism: str,
    card: dict[str, Any],
    article: dict[str, Any],
    model: str,
    reasoning: str,
    cache_dir: Path,
    runtime_dir: Path,
    timeout_seconds: int,
) -> dict[str, Any]:
    key = cache_key(organism, card, article, model, reasoning)
    cache_path = cache_dir / "judgments" / f"{key}.json"
    if cache_path.is_file():
        cached = read_json(cache_path)
        extraction = {name: cached[name] for name in EXTRACTION_KEYS}
        decision = derive_article_verdict(extraction)
        decision = apply_traceability_gate(
            decision,
            support_span_valid=cached.get("support_span_valid") is True,
            counter_span_valid=cached.get("counter_span_valid") is True,
        )
        cached.update(decision)
        cached["cache_hit"] = True
        return cached

    title = compact(article.get("title"))
    abstract = compact(article.get("abstract"))
    payload = {
        "organism": organism,
        "patient_case_card": prompt_case_card(card),
        "article": {"pmid": compact(article.get("pmid")), "title": title, "abstract": abstract},
    }
    raw = call_codex_json(
        executable=executable,
        model=model,
        reasoning=reasoning,
        system_prompt=SYSTEM_PROMPT,
        payload=payload,
        schema=EXTRACTION_SCHEMA,
        runtime_dir=runtime_dir,
        result_key=key,
        timeout_seconds=timeout_seconds,
    )
    extraction = {name: raw[name] for name in EXTRACTION_KEYS}
    decision = derive_article_verdict(extraction)
    article_blob = f"{title} {abstract}"
    support_span, support_valid = span(raw.get("support_span"), article_blob)
    counter_span, counter_valid = span(raw.get("counter_span"), article_blob)
    need_support = decision.get("verdict") in {"strong_match", "partial_match"} and not support_valid
    need_counter = decision.get("counterevidence") is True and not counter_valid
    repair_attempted = False
    repair_method = "not_needed"
    if need_support:
        support_span, support_valid = trigger_window(
            title,
            abstract,
            ["pneumonia", "infection", "infectious", "bacteremia", "sepsis", "meningitis"],
        )
    if need_counter:
        terms = (
            ["colonization", "colonisation", "contamination", "contaminant"]
            if extraction["human_clinical_context"] == "colonization_or_contamination"
            else ["reactivation", "bystander", "shedding", "latent"]
        )
        counter_span, counter_valid = trigger_window(title, abstract, terms)
    if need_support or need_counter:
        repair_attempted = True
        repair_method = "deterministic_trigger_window"
    need_support = need_support and not support_valid
    need_counter = need_counter and not counter_valid
    if need_support or need_counter:
        repair_key = f"repair_{key}"
        repair_cache = cache_dir / "repairs" / f"{repair_key}.json"
        if repair_cache.is_file():
            repaired = read_json(repair_cache)
        else:
            repaired = call_codex_json(
                executable=executable,
                model=model,
                reasoning=reasoning,
                system_prompt=REPAIR_SYSTEM_PROMPT,
                payload={
                    "organism": organism,
                    "need_support_span": need_support,
                    "need_counter_span": need_counter,
                    "frozen_extraction": extraction,
                    "article": {"title": title, "abstract": abstract},
                },
                schema=REPAIR_SCHEMA,
                runtime_dir=runtime_dir,
                result_key=repair_key,
                timeout_seconds=timeout_seconds,
            )
            write_json(repair_cache, repaired)
        if need_support:
            support_span, support_valid = span(repaired.get("support_span"), article_blob)
        if need_counter:
            counter_span, counter_valid = span(repaired.get("counter_span"), article_blob)
        repair_method = "bounded_llm_repair"
    decision = apply_traceability_gate(
        decision,
        support_span_valid=support_valid,
        counter_span_valid=counter_valid,
    )
    result = {
        "schema_version": "rag_re_casefit.article_judgment.v1",
        "pmid": compact(article.get("pmid")),
        "article_url": article.get("url"),
        "organism": organism,
        "case_card_hash": stable_hash(prompt_case_card(card)),
        "model": model,
        "reasoning_effort": reasoning,
        "prompt_sha256": SYSTEM_PROMPT_SHA256,
        "status": "ok",
        **extraction,
        "support_span": support_span,
        "support_span_valid": support_valid,
        "counter_span": counter_span,
        "counter_span_valid": counter_valid,
        "span_repair_attempted": repair_attempted,
        "span_repair_method": repair_method,
        "span_repair_prompt_sha256": REPAIR_PROMPT_SHA256,
        "rationale": compact(raw.get("rationale"))[:800],
        **decision,
        "cache_hit": False,
        "usage": None,
    }
    write_json(cache_path, result)
    return result


def judgment_ref(case: CandidateCase, article: dict[str, Any], model: str, reasoning: str) -> str:
    return cache_key(case.organism, case.case_card, article, model, reasoning)


def build_artifact(
    case: CandidateCase,
    judgments_by_key: dict[str, dict[str, Any]],
    model: str,
    reasoning: str,
) -> dict[str, Any]:
    if case.historical_artifact is not None:
        artifact = copy.deepcopy(case.historical_artifact)
        source_path = artifact.pop("_historical_source_path", None)
        artifact["review_tier"] = case.review_tier
        artifact.setdefault("provenance", {})["historical_artifact_reused"] = True
        artifact["provenance"]["historical_source_path"] = source_path
        artifact["provenance"]["replay_cohort"] = case.cohort
        artifact["provenance"]["replay_case_id"] = case.case_id
        return artifact
    judgments: list[dict[str, Any]] = []
    for article in case.articles:
        if not compact(article.get("abstract")):
            judgments.append(no_abstract_judgment(case.organism, case.case_card, article, model, reasoning))
        else:
            judgments.append(judgments_by_key[judgment_ref(case, article, model, reasoning)])
    aggregate = aggregate_judgments(judgments, AGGREGATION_CONFIG)
    return {
        "schema_version": "rag_re_casefit.candidate_output.v1",
        "patient_id": case.patient_id,
        "candidate_organism": case.organism,
        "canonical_organism": case.organism,
        "review_tier": case.review_tier,
        "case_card": case.case_card,
        "old_A": {
            "module": "A_LITERATURE",
            "positive": None,
            "status": "not_replayed_not_used_by_R5",
            "reason": "Historical R5 uses A1, B, and C; the older Module A is outside the R5 formula.",
        },
        "A1_judgments": judgments,
        "A1_aggregate": aggregate,
        "provenance": {
            "replay_cohort": case.cohort,
            "replay_case_id": case.case_id,
            "merged_file": str(case.merged_path.resolve()),
            "merged_sha256": file_sha256(case.merged_path),
            "packet_case_id": case.case_id,
            "prompt_sha256": SYSTEM_PROMPT_SHA256,
            "model": model,
            "reasoning_effort": reasoning,
            "historical_artifact_reused": False,
            "generation_is_answer_blind": True,
        },
    }


def strict_bool_or_none(value: Any) -> bool | None:
    if value is True or value is False or value is None:
        return value
    raise ValueError(f"Expected boolean or null, got {value!r}")


def evaluate_r5(case: CandidateCase, artifact: dict[str, Any], clinical_config: dict[str, Any]) -> dict[str, Any]:
    snapshot = case.review_row.get("evidence_snapshot")
    snapshot = snapshot if isinstance(snapshot, dict) else {}
    effective = dict(case.deterministic or snapshot)
    effective["organism_name"] = case.organism
    effective["canonical_organism_name"] = case.organism
    source_context = {
        "site_code": case.case_card.get("clinical_site"),
        "target_site": case.case_card.get("clinical_site"),
        "specimen_site_code": case.case_card.get("specimen_site"),
        "specimen_class": effective.get("specimen_class"),
        "specimen_alignment": effective.get("specimen_alignment"),
        "dominant_source": case.case_card.get("target_syndrome"),
    }
    merged_context = {
        "source_context": source_context,
        "effective_candidate": effective,
        "deterministic_candidate": case.deterministic,
        "review_metadata": {"tier": case.review_tier, "row": case.review_row},
    }
    modules = evaluate_modules(merged_context, effective, clinical_config)
    b = modules.get("B_STRICT") or {}
    c = modules.get("C") or {}
    d = modules.get("D") or {}
    a1 = artifact.get("A1_aggregate") or {}
    strong = int(a1.get("strong_match_count") or 0)
    partial = int(a1.get("partial_match_count") or 0)
    match_count = strong + partial
    b_absolute = strict_bool_or_none(b.get("absolute_strength_axis"))
    b_relative = strict_bool_or_none(b.get("relative_strength_axis"))
    b_site = strict_bool_or_none(b.get("site_aligned"))
    c_grade = compact(c.get("grade"))
    selected = c_grade == "C3" or (
        (strong >= 2 or match_count >= 7)
        and b_site is True
        and (b_absolute is True or b_relative is True)
    )
    if c_grade == "C3":
        selection_basis = "C3_same_organism_direct_evidence"
    elif selected:
        selection_basis = "A1_casefit_plus_site_aligned_B_axis"
    else:
        blockers: list[str] = []
        if strong < 2 and match_count < 7:
            blockers.append("A1_threshold_not_met")
        if b_site is not True:
            blockers.append("B_site_not_aligned_or_unknown")
        if b_absolute is not True and b_relative is not True:
            blockers.append("B_strength_axis_not_met")
        if c_grade != "C3":
            blockers.append("no_C3_direct_evidence")
        selection_basis = "|".join(blockers)
    return {
        "schema_version": "historical_ober_r5_candidate_replay_v1",
        "cohort": case.cohort,
        "case_id": case.case_id,
        "patient_id": case.patient_id,
        "organism": case.organism,
        "review_tier": case.review_tier,
        "historical_casefit_reused": case.historical_artifact is not None,
        "article_count": int(a1.get("article_count") or len(artifact.get("A1_judgments") or [])),
        "A1_strong": strong,
        "A1_partial": partial,
        "A1_match": match_count,
        "A1_tier": compact(a1.get("casefit_evidence_tier")),
        "B_strict": strict_bool_or_none(b.get("positive")),
        "B_absolute_axis": b_absolute,
        "B_relative_axis": b_relative,
        "B_site_aligned": b_site,
        "C_grade": c_grade,
        "D_state": compact(d.get("state")),
        "selected_R5": selected,
        "selection_basis": selection_basis,
        "casefit_file": str(case.artifact_path.resolve()),
        "answer_source_read": False,
    }


def write_results(output_dir: Path, rows: list[dict[str, Any]], manifest: dict[str, Any]) -> None:
    rows = sorted(rows, key=lambda row: (row["cohort"], int(row["patient_id"]), row["organism"]))
    write_json(output_dir / "candidate_results.json", rows)
    with (output_dir / "candidate_results.jsonl").open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")
    fields = list(rows[0]) if rows else []
    with (output_dir / "candidate_results.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    selected = [row for row in rows if row["selected_R5"]]
    write_json(output_dir / "r5_selected_additions.json", selected)
    summary: dict[str, Any] = {
        "schema_version": "historical_ober_r5_replay_summary_v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_count": len(rows),
        "scope_patient_count": sum(
            int(value) for value in (manifest.get("cohort_patient_counts") or {}).values()
        ),
        "candidate_patient_count": len({(row["cohort"], row["patient_id"]) for row in rows}),
        "selected_R5_count": len(selected),
        "historical_casefit_reused_count": sum(row["historical_casefit_reused"] for row in rows),
        "cohorts": {},
        "answer_source_read": False,
        "formula": (
            "C3 OR ((A1_strong >= 2 OR A1_strong + A1_partial >= 7) "
            "AND B_site_aligned AND (B_absolute_axis OR B_relative_axis))"
        ),
    }
    for cohort in sorted({row["cohort"] for row in rows}):
        cohort_rows = [row for row in rows if row["cohort"] == cohort]
        summary["cohorts"][cohort] = {
            "scope_patient_count": int((manifest.get("cohort_patient_counts") or {}).get(cohort, 0)),
            "candidate_count": len(cohort_rows),
            "selected_R5_count": sum(row["selected_R5"] for row in cohort_rows),
            "candidate_patient_count": len({row["patient_id"] for row in cohort_rows}),
        }
    historical_selected = ROOT / "outputs" / "patient_database_20260826_r5" / "r5_selected_additions.csv"
    historical_manifest = HISTORICAL_CASEFIT_DIR / "manifest.json"
    if historical_selected.is_file() and historical_manifest.is_file():
        with historical_selected.open("r", encoding="utf-8-sig", newline="") as handle:
            historical_selected_rows = list(csv.DictReader(handle))
        historical_casefit = read_json(historical_manifest)
        summary["exact_historical_KH33_snapshot"] = {
            "status": "archived_original_not_recomputed",
            "candidate_count": int(historical_casefit.get("output_count") or 0),
            "selected_R5_count": len(historical_selected_rows),
            "selected_path": str(historical_selected.resolve()),
            "casefit_manifest_path": str(historical_manifest.resolve()),
            "interpretation": (
                "This is the original KH33 development snapshot. The main replay metrics above "
                "apply the same historical rule to the currently frozen candidate inputs."
            ),
        }
    write_json(output_dir / "summary.json", summary)
    manifest.update(
        {
            "completion_status": "complete",
            "candidate_count": len(rows),
            "selected_R5_count": len(selected),
            "summary_path": str((output_dir / "summary.json").resolve()),
        }
    )
    write_json(output_dir / "run_manifest.json", manifest)
    (output_dir / "failures.json").unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replay the historical answer-blind OBER A1/B/C/R5 workflow across all cohorts."
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--model", default=MODEL_DEFAULT)
    parser.add_argument("--reasoning-effort", default=REASONING_DEFAULT)
    parser.add_argument("--max-workers", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=int, default=900)
    parser.add_argument("--cohort", action="append")
    parser.add_argument("--cohort-config", type=Path, help="Private JSON containing cohort labels and packet/manifest paths")
    parser.add_argument("--historical-casefit-dir", type=Path, default=HISTORICAL_CASEFIT_DIR)
    parser.add_argument("--case-id", action="append")
    parser.add_argument("--prepare-only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cohort_specs = load_cohort_specs(args.cohort_config) if args.cohort_config else COHORTS
    selected_cohorts = set(args.cohort or [spec.label for spec in cohort_specs])
    unknown_cohorts = selected_cohorts - {spec.label for spec in cohort_specs}
    if unknown_cohorts:
        raise ValueError(f"Unknown cohorts: {sorted(unknown_cohorts)}")
    cases = prepare_cases(output_dir, selected_cohorts, cohort_specs, args.historical_casefit_dir)
    if args.case_id:
        requested = set(args.case_id)
        cases = [case for case in cases if case.case_id in requested]
        missing = requested - {case.case_id for case in cases}
        if missing:
            raise ValueError(f"Unknown case ids in selected cohorts: {sorted(missing)}")
    clinical_config = read_json(CLINICAL_ROOT / "config" / "clinical_v1.json")
    cache_dir = output_dir / "_cache"
    runtime_dir = output_dir / "_codex_cli_runtime"
    manifest: dict[str, Any] = {
        "schema_version": "historical_ober_r5_replay_manifest_v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "generation_is_answer_blind": True,
        "ground_truth_read": False,
        "cohorts": sorted(selected_cohorts),
        "model": args.model,
        "reasoning_effort": args.reasoning_effort,
        "historical_prompt_sha256": SYSTEM_PROMPT_SHA256,
        "historical_clinical_config": str((CLINICAL_ROOT / "config" / "clinical_v1.json").resolve()),
        "candidate_count": len(cases),
        "historical_candidate_reuse_count": sum(case.historical_artifact is not None for case in cases),
        "source_packets": {
            spec.label: str(spec.packets.resolve())
            for spec in cohort_specs
            if spec.label in selected_cohorts
        },
        "cohort_patient_counts": {
            spec.label: int(read_json(spec.manifest).get("patient_directory_count") or 0)
            for spec in cohort_specs
            if spec.label in selected_cohorts
        },
        "completion_status": "prepared",
    }

    unique_tasks: dict[str, tuple[CandidateCase, dict[str, Any]]] = {}
    no_abstract_count = 0
    for case in cases:
        if case.historical_artifact is not None:
            continue
        for article in case.articles:
            if not compact(article.get("abstract")):
                no_abstract_count += 1
                continue
            key = judgment_ref(case, article, args.model, args.reasoning_effort)
            unique_tasks.setdefault(key, (case, article))
    cached_count = sum((cache_dir / "judgments" / f"{key}.json").is_file() for key in unique_tasks)
    manifest["article_slot_count"] = sum(len(case.articles) for case in cases)
    manifest["unique_llm_extraction_count"] = len(unique_tasks)
    manifest["cached_extraction_count_at_start"] = cached_count
    manifest["no_abstract_slot_count"] = no_abstract_count
    write_json(output_dir / "run_manifest.json", manifest)
    print(
        f"Prepared {len(cases)} candidates; historical reuse="
        f"{manifest['historical_candidate_reuse_count']}; unique LLM extractions="
        f"{len(unique_tasks)}; cached={cached_count}; no-abstract slots={no_abstract_count}.",
        flush=True,
    )
    if args.prepare_only:
        return 0

    executable = ensure_codex_cli_authenticated() if unique_tasks else ""
    judgments: dict[str, dict[str, Any]] = {}
    failures: list[dict[str, Any]] = []
    pending = [key for key in unique_tasks if not (cache_dir / "judgments" / f"{key}.json").is_file()]
    completed_count = 0
    if pending:
        print(
            f"Running {len(pending)} historical A1 extractions with {args.model}/"
            f"{args.reasoning_effort}, workers={args.max_workers}.",
            flush=True,
        )
        with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as executor:
            futures = {}
            for key in pending:
                case, article = unique_tasks[key]
                future = executor.submit(
                    judge_article,
                    executable=executable,
                    organism=case.organism,
                    card=case.case_card,
                    article=article,
                    model=args.model,
                    reasoning=args.reasoning_effort,
                    cache_dir=cache_dir,
                    runtime_dir=runtime_dir,
                    timeout_seconds=args.timeout_seconds,
                )
                futures[future] = (key, case, article)
            for future in as_completed(futures):
                key, case, article = futures[future]
                try:
                    judgments[key] = future.result()
                except Exception as exc:
                    failures.append(
                        {
                            "key": key,
                            "cohort": case.cohort,
                            "case_id": case.case_id,
                            "patient_id": case.patient_id,
                            "organism": case.organism,
                            "pmid": compact(article.get("pmid")),
                            "error": str(exc),
                        }
                    )
                completed_count += 1
                if completed_count == 1 or completed_count % 10 == 0 or completed_count == len(pending):
                    print(
                        f"A1 progress {completed_count}/{len(pending)}; failures={len(failures)}",
                        flush=True,
                    )
        if failures:
            write_json(output_dir / "failures.json", failures)
            manifest["completion_status"] = "incomplete_llm_failures"
            manifest["failure_count"] = len(failures)
            write_json(output_dir / "run_manifest.json", manifest)
            raise RuntimeError(
                f"{len(failures)} A1 extraction(s) failed; rerun the same command to resume"
            )

    for key, (case, article) in unique_tasks.items():
        if key in judgments:
            continue
        judgments[key] = judge_article(
            executable=executable,
            organism=case.organism,
            card=case.case_card,
            article=article,
            model=args.model,
            reasoning=args.reasoning_effort,
            cache_dir=cache_dir,
            runtime_dir=runtime_dir,
            timeout_seconds=args.timeout_seconds,
        )

    rows: list[dict[str, Any]] = []
    for index, case in enumerate(cases, 1):
        artifact = build_artifact(case, judgments, args.model, args.reasoning_effort)
        write_json(case.artifact_path, artifact)
        rows.append(evaluate_r5(case, artifact, clinical_config))
        if index % 25 == 0 or index == len(cases):
            print(f"R5 evaluation {index}/{len(cases)}", flush=True)
    write_results(output_dir, rows, manifest)
    print(
        f"Complete: candidates={len(rows)}, R5 additions="
        f"{sum(row['selected_R5'] for row in rows)}. Output={output_dir}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

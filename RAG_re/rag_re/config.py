from __future__ import annotations

import copy
import hashlib
import json
import os
from datetime import date
from pathlib import Path
from typing import Any, Dict


DEFAULT_CONFIG: Dict[str, Any] = {
    "schema_version": "rag_re.config.v1",
    "pubmed": {
        "max_articles": 10,
        "publication_cutoff": "2026-06-30",
        "timeout_seconds": 30,
        "retries": 3,
        "request_interval_seconds": 0.34,
        "email": "",
        "tool": "csie_rag_re",
    },
    "llm": {
        "model": "gpt-5.6-luna",
        "reasoning_effort": "low",
        # GPT-5 snapshots do not need an explicit temperature.  Leaving this null also
        # prevents unsupported-parameter failures in the Chat Completions endpoint.
        "temperature": None,
        "max_output_tokens": 1200,
        "prompt_version": "article_support_v1",
    },
    "organism_normalization": {
        # This map is part of the frozen experiment config.  Add aliases only on the
        # development split, then rerun retrieval instead of changing them on holdout.
        "accepted_synonyms": {
            "pneumocystis jiroveci": "Pneumocystis jirovecii",
            "pjp": "Pneumocystis jirovecii",
            "hsv": "Herpes simplex virus type 1",
            "hsv-1": "Herpes simplex virus type 1",
            "herpes simplex virus 1": "Herpes simplex virus type 1",
            "cmv": "Human cytomegalovirus",
            "crkp": "Klebsiella pneumoniae",
            "candida albican": "Candida albicans"
        }
    },
    "candidate_pool": {
        # Prior LLM missed-candidate expansion is a separate, secondary experiment.
        "allow_prior_llm_expansion": False,
        # Formal 2026-08 merge envelopes freeze these as the only rescue tiers.
        "merge_review_tiers": [
            "review_high_priority",
            "review_context_needed",
        ],
        # False follows READ: keep picked, but send only rescue tiers to RAG.
        "literature_on_baseline_picked": False,
    },
    "module_a_literature": {
        "min_support": 3,
        "min_support_ratio": 0.0,
        "min_judgeable_articles": 5,
        "require_exact_species": True,
        "require_human_clinical_evidence": True,
        "require_target_site_match": True,
        "require_evidence_span": True,
    },
    "module_b_mngs": {
        "require_site_alignment": True,
        "positive_specimen_alignments": [
            "Aligned",
            "Matched",
            "Site_aligned",
            "Sterile_or_Systemic"
        ],
        "negative_specimen_alignments": ["Not_aligned", "Mismatched", "Site_mismatch"],
        "positive_specimen_classes": [
            "S1_sterile_systemic",
            "S2_lower_respiratory"
        ],
        "max_rank_priority": 1,
        "min_reads_percentile": 0.80,
        "positive_reads_tiers": ["R3_HIGH", "R4_VERY_HIGH"],
        "min_dominance_tier": 2,
        "positive_signal_tiers": ["M1_STRONG", "M2_MODERATE"],
    },
    "module_c_direct_support": {
        "require_site_alignment": True,
        "positive_specimen_alignments": [
            "Aligned",
            "Matched",
            "Site_aligned",
            "Sterile_or_Systemic"
        ],
        "negative_specimen_alignments": ["Not_aligned", "Mismatched", "Site_mismatch"],
        "positive_specimen_classes": [
            "S1_sterile_systemic",
            "S2_lower_respiratory"
        ],
        "allow_direct_evidence_site_inference": True,
        "positive_modules": [
            "culture",
            "filmarray",
            "filmarray_gmtest",
            "targeted_molecular",
            "molecular_microbiology",
            "pcr",
            "antigen",
            "gm_test",
        ]
    },
    "decision": {
        "primary_experiment": "E0_OR_C_OR_A_AND_B",
        "upgrade_steps": 1,
        "best_level_after_upgrade": 2,
    },
}


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def load_config(
    path: Path | None = None, *, base: Dict[str, Any] | None = None
) -> Dict[str, Any]:
    config = copy.deepcopy(base if isinstance(base, dict) else DEFAULT_CONFIG)
    if path is not None:
        override = json.loads(path.read_text(encoding="utf-8-sig"))
        if not isinstance(override, dict):
            raise ValueError(f"Config root must be an object: {path}")
        config = _deep_merge(config, override)

    # Environment variables are runtime secrets/settings, not part of config files.
    config["pubmed"]["email"] = (
        os.getenv("PUBMED_EMAIL") or config["pubmed"].get("email") or ""
    )
    config["llm"]["model"] = (
        os.getenv("RAG_RE_MODEL")
        or os.getenv("RAG_MODEL")
        or config["llm"]["model"]
    )
    validate_config(config)
    return config


def validate_config(config: Dict[str, Any]) -> None:
    merge_tiers = config.get("candidate_pool", {}).get("merge_review_tiers")
    allowed_merge_tiers = {"review_high_priority", "review_context_needed"}
    if not isinstance(merge_tiers, list) or not merge_tiers:
        raise ValueError("candidate_pool.merge_review_tiers must be a non-empty list")
    if any(str(tier) not in allowed_merge_tiers for tier in merge_tiers):
        raise ValueError(
            "candidate_pool.merge_review_tiers may contain only review_high_priority "
            "and review_context_needed"
        )
    if len(set(map(str, merge_tiers))) != len(merge_tiers):
        raise ValueError("candidate_pool.merge_review_tiers must not contain duplicates")
    if not isinstance(
        config.get("candidate_pool", {}).get("literature_on_baseline_picked"), bool
    ):
        raise ValueError("candidate_pool.literature_on_baseline_picked must be boolean")
    max_articles = int(config["pubmed"]["max_articles"])
    if max_articles < 1:
        raise ValueError("pubmed.max_articles must be >= 1")
    threshold = int(config["module_a_literature"]["min_support"])
    if threshold < 1 or threshold > max_articles:
        raise ValueError("module_a_literature.min_support must be within 1..max_articles")
    ratio = float(config["module_a_literature"]["min_support_ratio"])
    if not 0 <= ratio <= 1:
        raise ValueError("module_a_literature.min_support_ratio must be within 0..1")
    min_judgeable = int(config["module_a_literature"]["min_judgeable_articles"])
    if min_judgeable < 1 or min_judgeable > max_articles:
        raise ValueError(
            "module_a_literature.min_judgeable_articles must be within 1..max_articles"
        )
    cutoff = str(config["pubmed"].get("publication_cutoff") or "")
    if cutoff:
        try:
            date.fromisoformat(cutoff)
        except ValueError as exc:
            raise ValueError("pubmed.publication_cutoff must be YYYY-MM-DD") from exc
    if float(config["pubmed"].get("request_interval_seconds", 0)) < 0:
        raise ValueError("pubmed.request_interval_seconds must be >= 0")
    temperature = config["llm"].get("temperature")
    if temperature is not None and not 0 <= float(temperature) <= 2:
        raise ValueError("llm.temperature must be null or within 0..2")
    if int(config["llm"].get("max_output_tokens", 0)) < 1:
        raise ValueError("llm.max_output_tokens must be >= 1")
    if config["llm"].get("reasoning_effort") not in {
        "none",
        "low",
        "medium",
        "high",
        "xhigh",
        "max",
    }:
        raise ValueError(
            "llm.reasoning_effort must be one of none/low/medium/high/xhigh/max"
        )
    for module_name in ("module_b_mngs", "module_c_direct_support"):
        module = config[module_name]
        if not isinstance(module.get("positive_specimen_alignments"), list):
            raise ValueError(f"{module_name}.positive_specimen_alignments must be a list")
        if not isinstance(module.get("positive_specimen_classes"), list):
            raise ValueError(f"{module_name}.positive_specimen_classes must be a list")
        if not isinstance(module.get("negative_specimen_alignments"), list):
            raise ValueError(f"{module_name}.negative_specimen_alignments must be a list")
    module_b = config["module_b_mngs"]
    if int(module_b.get("max_rank_priority", 0)) < 1:
        raise ValueError("module_b_mngs.max_rank_priority must be >= 1")
    percentile = float(module_b.get("min_reads_percentile", 0))
    if not 0 <= percentile <= 1:
        raise ValueError("module_b_mngs.min_reads_percentile must be within 0..1")
    if int(module_b.get("min_dominance_tier", 0)) < 0:
        raise ValueError("module_b_mngs.min_dominance_tier must be >= 0")
    decision = config["decision"]
    if not 0 <= int(decision.get("upgrade_steps", 0)) <= 4:
        raise ValueError("decision.upgrade_steps must be within 0..4")
    if not 1 <= int(decision.get("best_level_after_upgrade", 0)) <= 5:
        raise ValueError("decision.best_level_after_upgrade must be within 1..5")
    aliases = config.get("organism_normalization", {}).get("accepted_synonyms", {})
    if not isinstance(aliases, dict) or any(
        not str(alias).strip() or not str(canonical).strip()
        for alias, canonical in aliases.items()
    ):
        raise ValueError("organism_normalization.accepted_synonyms must map nonempty names")
    if config["llm"].get("prompt_version") != "article_support_v1":
        raise ValueError(
            "This release implements only llm.prompt_version=article_support_v1"
        )


def public_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """Return a JSON-safe config snapshot without API keys."""
    snapshot = copy.deepcopy(config)
    # The contact email changes neither retrieval nor decisions and may be personal.
    if isinstance(snapshot.get("pubmed"), dict):
        snapshot["pubmed"].pop("email", None)
    return snapshot


def config_hash(config: Dict[str, Any]) -> str:
    payload = json.dumps(public_config(config), ensure_ascii=False, sort_keys=True).encode(
        "utf-8"
    )
    return hashlib.sha256(payload).hexdigest()


def load_env_file(path: Path) -> None:
    """Load NAME=VALUE pairs without overwriting variables already in the process."""
    for raw_line in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        name, value = line.split("=", 1)
        name = name.strip()
        if name.startswith("$env:"):
            name = name[5:]
        value = value.strip().strip('"').strip("'")
        if name and name not in os.environ:
            os.environ[name] = value

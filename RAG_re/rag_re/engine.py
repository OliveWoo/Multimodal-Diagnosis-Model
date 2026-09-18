from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

from .config import config_hash, public_config
from .input_adapter import (
    baseline_picked_names,
    candidate_literature_eligible,
    candidate_provenance,
    candidate_review_tier,
    candidate_level,
    candidate_name,
    canonical_organism_name,
    load_bundle,
    normalize_name,
    patient_id_from,
    select_candidate_pool,
    source_context,
)
from .io_utils import sha256_file, stable_hash
from .judge import ArticleJudge, SYSTEM_PROMPT_SHA256
from .pubmed import PubMedClient
from .rules import (
    EXPERIMENT_META,
    build_experiment_summary,
    experiment_values,
    module_a_from_judgments,
    module_a_not_requested,
    module_a_skipped,
    module_b_mngs,
    module_c_direct_support,
)
from .versioning import (
    PIPELINE_FINGERPRINT,
    PIPELINE_VERSION,
    REPLAY_COMPATIBLE_PIPELINE_FINGERPRINTS,
)


RULE_INPUT_FIELDS = {
    "organism_name",
    "pathogen",
    "pathogen_name",
    "name",
    "classification",
    "organism",
    "integrated_causative_level",
    "recommended_level",
    "basis_level",
    "level",
    "original_level",
    "mngs_signal_tier",
    "rank_priority",
    "reads",
    "reads_tier",
    "reads_percentile",
    "dominance_tier",
    "specimen_alignment",
    "specimen_class",
    "module_support_summary",
    "key_evidence",
    "culture",
    "filmarray",
    "filmarray_gmtest",
    "targeted_molecular",
    "molecular_microbiology",
    "pcr",
    "antigen",
    "gm_test",
    "source_category",
    "rank_rule",
    "applied_rules",
    "rag_re_candidate_provenance",
    "rag_re_review_tier",
    "rag_re_literature_eligible",
}


def _rule_input(candidate: Dict[str, Any]) -> Dict[str, Any]:
    return {key: value for key, value in candidate.items() if key in RULE_INPUT_FIELDS}


def _level_text(value: int | None) -> str:
    return f"Level {value}" if value is not None else "Level unknown"


def _upgraded_level(level: int | None, config: Dict[str, Any]) -> int | None:
    if level is None:
        return None
    return min(
        level,
        max(
            int(config["decision"].get("best_level_after_upgrade", 2)),
            level - int(config["decision"].get("upgrade_steps", 1)),
        ),
    )


def _compatibility_row(
    row: Dict[str, Any], primary_experiment: str, config: Dict[str, Any]
) -> Dict[str, Any]:
    selected = row["experiments"].get(primary_experiment)
    original = row.get("original_level_number")
    primary_semantics = EXPERIMENT_META[primary_experiment]["semantics"]
    is_precision_pruning = primary_semantics == "precision_pruning_policy"
    evidence_gate = {
        "E0_OR_A": "A_LITERATURE",
        "E0_OR_A_OR_B": "AB_OR",
        "E0_OR_A_OR_B_OR_C": "ABC_OR",
        "E0_OR_C": "C_DIRECT_SUPPORT",
        "E0_OR_C_OR_A_AND_B": "GATED_C_OR_A_AND_B",
    }.get(primary_experiment, primary_experiment)
    new_evidence_positive = row["experiments"].get(evidence_gate) is True
    upgrade_selected = (
        selected is True
        and primary_experiment != "E0_BASELINE"
        and new_evidence_positive
        and not is_precision_pruning
    )
    recommended = _upgraded_level(original, config) if upgrade_selected else original
    if is_precision_pruning and selected is True:
        decision = "保留"
        category = "validated_candidate"
        status = "positive"
    elif is_precision_pruning and selected is False:
        decision = "排除"
        category = "pruned_insufficient_confirmation"
        status = "negative"
    elif is_precision_pruning:
        decision = "人工複核"
        category = "insufficient_evidence_monitor"
        status = "abstain"
    elif selected is True:
        decision = "升級" if recommended != original else "維持"
        category = "reasonable_coinfection_or_treat_worthy"
        status = "positive"
    elif selected is False:
        decision = "維持"
        category = "insufficient_evidence_monitor"
        status = "negative"
    else:
        decision = "無法判定"
        category = "insufficient_evidence_monitor"
        status = "abstain"

    modules = row["modules"]
    a, b, c = modules["A"], modules["B"], modules["C"]
    rationale = (
        f"預先指定策略 {primary_experiment}={status}。"
        f"A 文獻票決={a.get('positive')}（有效支持 {a.get('support_count', 0)}/"
        f"{a.get('judgeable_count', 0)}）；"
        f"B mNGS 訊號={b.get('positive')}；C 非 mNGS 直接支持={c.get('positive')}。"
        "文獻模組只代表一般可致病性，不能單獨證明本病人的感染來源。"
    )
    return {
        "pathogen": row["organism_name"],
        "original_level": _level_text(original),
        "recommended_level": _level_text(recommended),
        "decision": decision,
        "decision_status": status,
        "selected_experiment": primary_experiment,
        "clinical_reasonableness_category": category,
        "is_infection_source": selected if selected in (True, False) else None,
        "rationale": rationale,
        "supporting_references": [
            f"PMID: {pmid}" for pmid in a.get("valid_support_pmids", []) if pmid
        ],
    }


def _candidate_row(
    *,
    candidate: Dict[str, Any],
    context: Dict[str, str],
    baseline_names: set[str],
    module_a: Dict[str, Any],
    retrieval: Dict[str, Any] | None,
    judgments: List[Dict[str, Any]],
    config: Dict[str, Any],
) -> Dict[str, Any]:
    name = candidate_name(candidate)
    canonical_name = canonical_organism_name(name, config)
    baseline = normalize_name(canonical_name) in baseline_names
    module_b = module_b_mngs(candidate, config["module_b_mngs"])
    module_c = module_c_direct_support(
        candidate, config["module_c_direct_support"], context=context
    )
    experiments = experiment_values(
        baseline=baseline,
        a=module_a.get("positive"),
        b=module_b.get("positive"),
        c=module_c.get("positive"),
    )
    original = candidate_level(candidate)
    literature_level = (
        _upgraded_level(original, config) if module_a.get("positive") is True else original
    )
    return {
        "organism_name": name,
        "canonical_organism_name": canonical_name,
        "candidate_provenance": candidate_provenance(candidate),
        "review_tier": candidate_review_tier(candidate),
        "literature_eligible": candidate_literature_eligible(candidate),
        "classification": candidate.get("classification"),
        "original_level": _level_text(original),
        "original_level_number": original,
        "source_context": context,
        "baseline_selected": baseline,
        "rule_input": _rule_input(candidate),
        "modules": {
            "A": {
                **module_a,
                "literature_recommended_level": _level_text(literature_level),
                "upgrade_is_patient_causality": False,
            },
            "B": module_b,
            "C": module_c,
        },
        "literature_evidence": {"retrieval": retrieval, "judgments": judgments},
        "experiments": experiments,
    }


class RagReEngine:
    def __init__(
        self,
        config: Dict[str, Any],
        *,
        cache_dir: Path | None = None,
        pubmed_client: PubMedClient | None = None,
        article_judge: ArticleJudge | None = None,
        skip_literature: bool = False,
    ):
        self.config = config
        self.skip_literature = skip_literature
        primary = str(config.get("decision", {}).get("primary_experiment") or "")
        if primary not in EXPERIMENT_META:
            raise ValueError(f"Unknown decision.primary_experiment: {primary}")
        self.pubmed = pubmed_client or PubMedClient(config["pubmed"], cache_dir)
        self.judge = article_judge
        if not skip_literature and self.judge is None:
            self.judge = ArticleJudge(config["llm"], cache_dir)

    def run_file(self, input_path: Path) -> Dict[str, Any]:
        root, bundle = load_bundle(input_path)
        patient_id = patient_id_from(input_path, root, bundle)
        baseline_names_list = baseline_picked_names(bundle)
        baseline_names = {
            normalize_name(canonical_organism_name(name, self.config))
            for name in baseline_names_list
        }
        rows: List[Dict[str, Any]] = []
        errors: List[Dict[str, Any]] = []
        raw_candidates = bundle.get("pathogen_candidates")
        if not isinstance(raw_candidates, list):
            raise ValueError("pathogen_candidates must be a list")
        invalid_indexes = [
            index
            for index, candidate in enumerate(raw_candidates)
            if not isinstance(candidate, dict) or not candidate_name(candidate)
        ]
        if invalid_indexes:
            raise ValueError(
                "Every pathogen_candidates entry must be an object with an organism name; "
                f"invalid indexes: {invalid_indexes}"
            )
        all_candidates, pool_meta = select_candidate_pool(
            root,
            bundle,
            self.config.get("candidate_pool", {}),
        )
        excluded_prior_llm = [
            candidate
            for candidate in all_candidates
            if candidate_provenance(candidate) == "prior_llm_missed_candidate_expansion"
        ]
        allow_prior_llm = bool(
            self.config.get("candidate_pool", {}).get("allow_prior_llm_expansion", False)
        )
        candidates = (
            all_candidates
            if allow_prior_llm
            else [
                candidate
                for candidate in all_candidates
                if candidate_provenance(candidate) != "prior_llm_missed_candidate_expansion"
            ]
        )
        seen_names = set()
        duplicate_names = set()
        for candidate in candidates:
            key = normalize_name(
                canonical_organism_name(candidate_name(candidate), self.config)
            )
            if key in seen_names:
                duplicate_names.add(candidate_name(candidate))
            seen_names.add(key)
        if duplicate_names:
            raise ValueError(
                "Duplicate patient-pathogen candidates are not a valid experimental unit: "
                + ", ".join(sorted(duplicate_names))
            )
        candidate_name_set = {
            normalize_name(canonical_organism_name(candidate_name(row), self.config))
            for row in candidates
        }
        missing_baseline = sorted(baseline_names - candidate_name_set)
        if missing_baseline:
            raise ValueError(
                "E0 baseline contains pathogens outside the frozen candidate pool: "
                + ", ".join(missing_baseline)
            )

        for candidate in candidates:
            name = candidate_name(candidate)
            context = source_context(bundle, candidate)
            retrieval: Dict[str, Any] | None = None
            judgments: List[Dict[str, Any]] = []
            if self.skip_literature:
                module_a = module_a_skipped()
            elif not candidate_literature_eligible(candidate):
                retrieval = {
                    "status": "not_requested",
                    "organism": name,
                    "articles": [],
                    "reason": (
                        "Formal merge baseline picked pathogens are retained without "
                        "a new literature call."
                    ),
                }
                module_a = module_a_not_requested(
                    "Baseline picked pathogen retained per formal merge input policy."
                )
            else:
                try:
                    canonical_name = canonical_organism_name(name, self.config)
                    retrieval = self.pubmed.search(
                        canonical_name,
                        context,
                        int(self.config["pubmed"]["max_articles"]),
                    )
                    for article in retrieval.get("articles", []):
                        try:
                            assert self.judge is not None
                            judgments.append(
                                self.judge.judge_one(
                                    canonical_name, context["target_site"], article
                                )
                            )
                        except Exception as exc:
                            judgments.append(
                                {
                                    "pmid": str(article.get("pmid") or ""),
                                    "status": "judge_error",
                                    "verdict": "unclear",
                                    "error": f"{type(exc).__name__}: {exc}",
                                }
                            )
                            errors.append(
                                {
                                    "organism_name": name,
                                    "stage": "article_judge",
                                    "pmid": str(article.get("pmid") or ""),
                                    "error": f"{type(exc).__name__}: {exc}",
                                }
                            )
                    module_a = module_a_from_judgments(
                        retrieval, judgments, self.config["module_a_literature"]
                    )
                except Exception as exc:
                    retrieval = {
                        "status": "error",
                        "organism": name,
                        "error": f"{type(exc).__name__}: {exc}",
                        "articles": [],
                    }
                    module_a = module_a_from_judgments(
                        retrieval, [], self.config["module_a_literature"]
                    )
                    errors.append(
                        {
                            "organism_name": name,
                            "stage": "pubmed_or_judge",
                            "error": f"{type(exc).__name__}: {exc}",
                        }
                    )

            rows.append(
                _candidate_row(
                    candidate=candidate,
                    context=context,
                    baseline_names=baseline_names,
                    module_a=module_a,
                    retrieval=retrieval,
                    judgments=judgments,
                    config=self.config,
                )
            )

        is_formal_merge = pool_meta.get("input_candidate_mode") == "formal_merge_review_tiered"
        candidate_pool_provenance = (
            str(pool_meta["candidate_pool_provenance"])
            if is_formal_merge
            else (
                "expanded_with_prior_llm_missed_candidate_review"
                if allow_prior_llm and excluded_prior_llm
                else "original_candidate_pool"
            )
        )
        return self._assemble(
            patient_id=patient_id,
            rows=rows,
            baseline_names=baseline_names_list,
            errors=errors,
            input_meta={
                **pool_meta,
                "input_file": input_path.name,
                "input_sha256": sha256_file(input_path),
                "candidate_count": len(rows),
                "original_candidate_count": len(all_candidates) - len(excluded_prior_llm),
                "prior_llm_expansion_count": len(excluded_prior_llm),
                "excluded_prior_llm_expansion_count": (
                    0 if allow_prior_llm else len(excluded_prior_llm)
                ),
                "excluded_prior_llm_expansion_names": (
                    []
                    if allow_prior_llm
                    else [candidate_name(row) for row in excluded_prior_llm]
                ),
                "literature_eligible_candidate_count": sum(
                    candidate_literature_eligible(row) for row in candidates
                ),
                "candidate_pool_provenance": candidate_pool_provenance,
                "candidate_pool_sha256": stable_hash(
                    [
                        {
                            "organism_name": candidate_name(row),
                            "provenance": candidate_provenance(row),
                        }
                        for row in candidates
                    ]
                ),
            },
        )

    def _assemble(
        self,
        *,
        patient_id: str,
        rows: List[Dict[str, Any]],
        baseline_names: List[str],
        errors: List[Dict[str, Any]],
        input_meta: Dict[str, Any],
        replayed_from: str | None = None,
        evidence_provenance: Dict[str, Any] | None = None,
        literature_attempted: bool | None = None,
    ) -> Dict[str, Any]:
        primary = str(self.config["decision"]["primary_experiment"])
        if primary not in EXPERIMENT_META:
            raise ValueError(f"Unknown decision.primary_experiment: {primary}")
        conclusions = [_compatibility_row(row, primary, self.config) for row in rows]
        pmids = sorted(
            {
                pmid
                for row in rows
                for pmid in row["modules"]["A"].get("valid_support_pmids", [])
                if pmid
            }
        )
        status = "complete" if not errors else "partial_with_errors"
        payload = {
            "schema_version": "rag_re.output.v1",
            "pipeline_version": PIPELINE_VERSION,
            "pipeline_fingerprint": PIPELINE_FINGERPRINT,
            "run_status": status,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "patient_id": str(patient_id),
            "input_meta": input_meta,
            "config_hash": config_hash(self.config),
            "config_snapshot": public_config(self.config),
            "prompt_sha256": (evidence_provenance or {}).get(
                "prompt_sha256", SYSTEM_PROMPT_SHA256
            ),
            "evidence_provenance": evidence_provenance
            or {
                "source_config_hash": config_hash(self.config),
                "pubmed": public_config(self.config)["pubmed"],
                "llm": public_config(self.config)["llm"],
                "prompt_sha256": SYSTEM_PROMPT_SHA256,
            },
            "methodology": {
                "E0": "Existing best_available_summary.picked_pathogens baseline.",
                "A": "Per-pathogen PubMed top-10, one independent LLM verdict per article.",
                "B": "mNGS-only strength features; no literature or non-mNGS evidence.",
                "C": "Same-organism non-mNGS direct support only.",
                "primary_policy": EXPERIMENT_META[primary]["expression"],
                "primary_semantics": EXPERIMENT_META[primary]["semantics"],
                "warning": (
                    "A is biological plausibility, not patient-level causality. "
                    "Screening arms must not be interpreted as deployed clinical decisions. "
                    "Precision-pruning arms retain only E0 candidates confirmed by B/C."
                ),
            },
            "execution": {
                "mode": "replay" if replayed_from else "run",
                "skip_literature": self.skip_literature,
                "literature_eligible_candidate_count": sum(
                    row.get("literature_eligible") is True for row in rows
                ),
                "literature_not_requested_candidate_count": sum(
                    (row.get("modules") or {}).get("A", {}).get("status")
                    == "not_requested"
                    for row in rows
                ),
                "literature_attempted_candidate_count": sum(
                    (row.get("modules") or {}).get("A", {}).get("status")
                    not in {"skipped", "not_requested"}
                    for row in rows
                ),
                "literature_evidence_available": any(
                    (row.get("literature_evidence") or {}).get("judgments") for row in rows
                ),
                "literature_attempted": (
                    not self.skip_literature
                    if literature_attempted is None
                    else literature_attempted
                ),
            },
            "baseline_picked_pathogens": baseline_names,
            "candidates": rows,
            "experiments": build_experiment_summary(rows),
            "rag_answer": {
                "selected_experiment": primary,
                "summary": (
                    f"RAG_re 以 {primary} 作為相容輸出的預先指定策略；"
                    "A/B/C 與所有消融結果另存於 experiments。"
                ),
                "level_reassessment_conclusion": conclusions,
                "recommendations": [],
                "limitations": [
                    "研究用輸出，未經前瞻性臨床驗證。",
                    "PubMed 文獻支持不能單獨證明個別病人的感染來源。",
                ],
                "cited_pmids": pmids,
            },
            "errors": errors,
        }
        if replayed_from:
            payload["replayed_from_sha256"] = replayed_from
        payload["artifact_sha256_basis"] = stable_hash(
            {
                "patient_id": patient_id,
                "input_meta": input_meta,
                "config_hash": payload["config_hash"],
                "pipeline_fingerprint": payload["pipeline_fingerprint"],
                "prompt_sha256": payload["prompt_sha256"],
                "candidates": rows,
            }
        )
        return payload

    def replay(self, source: Dict[str, Any]) -> Dict[str, Any]:
        if source.get("schema_version") != "rag_re.output.v1":
            raise ValueError("Replay requires a rag_re.output.v1 artifact")
        source_fingerprint = str(source.get("pipeline_fingerprint") or "")
        if source_fingerprint not in {
            PIPELINE_FINGERPRINT,
            *REPLAY_COMPATIBLE_PIPELINE_FINGERPRINTS,
        }:
            raise ValueError(
                "Replay source was generated by a different pipeline implementation; "
                "rerun it with this version before changing thresholds"
            )
        source_config = source.get("config_snapshot")
        if not isinstance(source_config, dict):
            raise ValueError("Replay source is missing its frozen config_snapshot")
        current_public = public_config(self.config)
        for section in ("pubmed", "llm", "organism_normalization", "candidate_pool"):
            if source_config.get(section) != current_public.get(section):
                raise ValueError(
                    f"Replay cannot change evidence-generation section {section}; "
                    "rerun retrieval/judgment instead"
                )
        source_prompt = str(source.get("prompt_sha256") or "")
        if source_prompt != SYSTEM_PROMPT_SHA256:
            raise ValueError(
                "Replay evidence prompt hash differs from this code version; "
                "use the original version or rerun Module A"
            )
        rows: List[Dict[str, Any]] = []
        baseline_names_list = list(source.get("baseline_picked_pathogens") or [])
        baseline_names = {
            normalize_name(canonical_organism_name(name, self.config))
            for name in baseline_names_list
        }
        for old in source.get("candidates", []):
            candidate = old.get("rule_input") or {}
            evidence = old.get("literature_evidence") or {}
            retrieval = evidence.get("retrieval") or {"status": "skipped"}
            judgments = evidence.get("judgments") or []
            for judgment in judgments:
                if not isinstance(judgment, dict):
                    raise ValueError("Replay judgments must be objects")
                judgment_prompt = judgment.get("prompt_sha256")
                judgment_model = judgment.get("model")
                if judgment_prompt not in (None, "", source_prompt):
                    raise ValueError("Replay contains mixed article-judgment prompt hashes")
                if judgment_model not in (None, "", source_config["llm"].get("model")):
                    raise ValueError("Replay contains mixed article-judgment models")
            module_a = module_a_from_judgments(
                retrieval, judgments, self.config["module_a_literature"]
            )
            rows.append(
                _candidate_row(
                    candidate=candidate,
                    context=old.get("source_context") or {},
                    baseline_names=baseline_names,
                    module_a=module_a,
                    retrieval=retrieval,
                    judgments=judgments,
                    config=self.config,
                )
            )
        return self._assemble(
            patient_id=str(source.get("patient_id") or ""),
            rows=rows,
            baseline_names=baseline_names_list,
            errors=list(source.get("errors") or []),
            input_meta=dict(source.get("input_meta") or {}),
            replayed_from=stable_hash(
                {key: value for key, value in source.items() if not str(key).startswith("_")}
            ),
            evidence_provenance=dict(source.get("evidence_provenance") or {
                "source_config_hash": source.get("config_hash"),
                "pubmed": source_config.get("pubmed"),
                "llm": source_config.get("llm"),
                "prompt_sha256": source_prompt,
            }),
            literature_attempted=bool(
                (source.get("execution") or {}).get("literature_attempted")
                or (
                    (source.get("execution") or {}).get("skip_literature") is False
                )
            ),
        )

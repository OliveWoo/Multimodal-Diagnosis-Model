import copy
import json
import tempfile
import unittest
from pathlib import Path

from rag_re.config import DEFAULT_CONFIG
from rag_re.engine import RagReEngine


class FakePubMed:
    def search(self, organism, context, max_articles):
        return {
            "status": "ok",
            "query": f"fixed query for {organism}",
            "articles": [
                {
                    "pmid": str(index),
                    "title": f"{organism} pneumonia study {index}",
                    "abstract": "Human pneumonia evidence at the target site.",
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{index}/",
                }
                for index in range(1, 6)
            ],
        }


class FakeJudge:
    def __init__(self):
        self.calls = []

    def judge_one(self, organism, target_site, article):
        self.calls.append((organism, target_site, set(article)))
        support = int(article["pmid"]) <= 3
        return {
            "pmid": article["pmid"],
            "status": "ok",
            "verdict": "support" if support else "unclear",
            "organism_scope": "exact_species",
            "human_clinical_evidence": support,
            "target_site_match": support,
            "evidence_span": "Human pneumonia evidence",
            "evidence_span_valid": support,
        }


class FailingPubMed:
    def search(self, organism, context, max_articles):
        raise RuntimeError("network unavailable")


class EngineTests(unittest.TestCase):
    def setUp(self):
        self.config = copy.deepcopy(DEFAULT_CONFIG)

    def _input(self, root: Path) -> Path:
        payload = {
            "agent_bundle": {
                "dominant_source": "Respiratory",
                "best_available_summary": {
                    "picked_pathogens": [{"organism_name": "Organism alpha"}]
                },
                "pathogen_candidates": [
                    {
                        "organism_name": "Organism alpha",
                        "classification": "Bacterial",
                        "integrated_causative_level": "Level 4",
                        "mngs_signal_tier": "M2_moderate",
                        "rank_priority": 1,
                        "reads_percentile": 0.9,
                        "specimen_alignment": "Aligned",
                        "specimen_class": "S2_lower_respiratory",
                        "module_support_summary": {
                            "culture": "Not_available",
                            "filmarray_gmtest": "Support",
                            "host": "Support",
                        },
                        "integrated_reasoning": ["must never reach article judge"],
                        "llm_missed_candidate_review": ["must never reach article judge"],
                    }
                ],
            }
        }
        path = root / "NGS_patient_900.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_engine_runs_atomic_arms_and_allowlists_judge_payload(self):
        with tempfile.TemporaryDirectory() as temp:
            input_path = self._input(Path(temp))
            judge = FakeJudge()
            engine = RagReEngine(
                self.config,
                pubmed_client=FakePubMed(),
                article_judge=judge,
            )
            result = engine.run_file(input_path)
        row = result["candidates"][0]
        self.assertEqual(result["patient_id"], "900")
        self.assertIs(row["modules"]["A"]["positive"], True)
        self.assertIs(row["modules"]["B"]["positive"], True)
        self.assertIs(row["modules"]["C"]["positive"], True)
        self.assertIs(row["experiments"]["E0_BASELINE"], True)
        self.assertIs(row["experiments"]["GATED_C_OR_A_AND_B"], True)
        self.assertEqual(
            result["rag_answer"]["level_reassessment_conclusion"][0]["recommended_level"],
            "Level 3",
        )
        self.assertEqual(len(judge.calls), 5)
        self.assertTrue(all("integrated_reasoning" not in keys for _, _, keys in judge.calls))

    def test_precision_pruning_primary_retains_or_prunes_without_level_upgrade(self):
        with tempfile.TemporaryDirectory() as temp:
            input_path = self._input(Path(temp))
            config = copy.deepcopy(self.config)
            config["decision"]["primary_experiment"] = "E0_AND_B_OR_C"
            result = RagReEngine(
                config,
                pubmed_client=FakePubMed(),
                article_judge=FakeJudge(),
            ).run_file(input_path)
        row = result["candidates"][0]
        conclusion = result["rag_answer"]["level_reassessment_conclusion"][0]
        self.assertIs(row["experiments"]["E0_AND_B_OR_C"], True)
        self.assertEqual(conclusion["decision"], "保留")
        self.assertEqual(conclusion["recommended_level"], conclusion["original_level"])
        self.assertEqual(result["methodology"]["primary_semantics"], "precision_pruning_policy")

    def test_formal_merge_calls_literature_only_for_required_review_tiers(self):
        payload = {
            "patient_id": "901",
            "deterministic_max": {
                "dominant_source": "Respiratory",
                "best_available_summary": {
                    "picked_pathogens": [{"organism_name": "Picked species"}]
                },
                "pathogen_candidates": [
                    {
                        "organism_name": "Picked species",
                        "integrated_causative_level": "Level 2",
                        "specimen_alignment": "Aligned",
                        "specimen_class": "S2_lower_respiratory",
                    },
                    {
                        "organism_name": "Review species",
                        "integrated_causative_level": "Level 4",
                        "rank_priority": 1,
                        "reads_percentile": 0.9,
                        "specimen_alignment": "Aligned",
                        "specimen_class": "S2_lower_respiratory",
                    },
                    {"organism_name": "Low species", "integrated_causative_level": "Level 5"},
                ],
            },
            "llm_missed_candidate_review": {
                "review_high_priority": [{"organism_name": "Review species"}],
                "review_context_needed": [],
                "review_low_specificity": [{"organism_name": "Low species"}],
                "omitted_with_reason": [],
                "review_omitted_with_reason": [],
            },
        }
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "NGS_patient_901.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            judge = FakeJudge()
            result = RagReEngine(
                self.config,
                pubmed_client=FakePubMed(),
                article_judge=judge,
            ).run_file(path)
        self.assertEqual(
            [row["organism_name"] for row in result["candidates"]],
            ["Picked species", "Review species"],
        )
        self.assertEqual(result["candidates"][0]["modules"]["A"]["status"], "not_requested")
        self.assertEqual(result["candidates"][1]["modules"]["A"]["status"], "ok")
        self.assertEqual(len(judge.calls), 5)
        self.assertEqual(result["input_meta"]["literature_eligible_candidate_count"], 1)
        self.assertEqual(result["input_meta"]["excluded_review_tier_counts"]["review_low_specificity"], 1)

    def test_pubmed_failure_produces_abstention_and_partial_status(self):
        with tempfile.TemporaryDirectory() as temp:
            input_path = self._input(Path(temp))
            engine = RagReEngine(
                self.config,
                pubmed_client=FailingPubMed(),
                article_judge=FakeJudge(),
            )
            result = engine.run_file(input_path)
        self.assertEqual(result["run_status"], "partial_with_errors")
        self.assertIsNone(result["candidates"][0]["modules"]["A"]["positive"])
        replayed = RagReEngine(self.config, skip_literature=True).replay(result)
        self.assertEqual(replayed["run_status"], "partial_with_errors")
        self.assertTrue(replayed["errors"])

    def test_replay_changes_threshold_without_provider_calls(self):
        with tempfile.TemporaryDirectory() as temp:
            input_path = self._input(Path(temp))
            engine = RagReEngine(
                self.config,
                pubmed_client=FakePubMed(),
                article_judge=FakeJudge(),
            )
            original = engine.run_file(input_path)
        stricter = copy.deepcopy(self.config)
        stricter["module_a_literature"]["min_support"] = 4
        replayed = RagReEngine(stricter, skip_literature=True).replay(original)
        self.assertIs(replayed["candidates"][0]["modules"]["A"]["positive"], False)
        self.assertNotEqual(original["config_hash"], replayed["config_hash"])

    def test_duplicate_patient_pathogen_is_rejected_before_provider_calls(self):
        with tempfile.TemporaryDirectory() as temp:
            input_path = self._input(Path(temp))
            payload = json.loads(input_path.read_text(encoding="utf-8"))
            payload["agent_bundle"]["pathogen_candidates"].append(
                dict(payload["agent_bundle"]["pathogen_candidates"][0])
            )
            input_path.write_text(json.dumps(payload), encoding="utf-8")
            judge = FakeJudge()
            engine = RagReEngine(
                self.config,
                pubmed_client=FakePubMed(),
                article_judge=judge,
            )
            with self.assertRaisesRegex(ValueError, "Duplicate patient-pathogen"):
                engine.run_file(input_path)
            self.assertEqual(judge.calls, [])

    def test_replay_preserves_organism_alias_key(self):
        with tempfile.TemporaryDirectory() as temp:
            input_path = self._input(Path(temp))
            payload = json.loads(input_path.read_text(encoding="utf-8"))
            candidate = payload["agent_bundle"]["pathogen_candidates"][0]
            candidate["organism"] = candidate.pop("organism_name")
            payload["agent_bundle"]["best_available_summary"]["picked_pathogens"] = []
            input_path.write_text(json.dumps(payload), encoding="utf-8")
            original = RagReEngine(
                self.config,
                pubmed_client=FakePubMed(),
                article_judge=FakeJudge(),
            ).run_file(input_path)
        replayed = RagReEngine(self.config, skip_literature=True).replay(original)
        self.assertEqual(replayed["candidates"][0]["organism_name"], "Organism alpha")

    def test_malformed_candidate_and_baseline_outside_pool_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            input_path = self._input(Path(temp))
            payload = json.loads(input_path.read_text(encoding="utf-8"))
            payload["agent_bundle"]["pathogen_candidates"].append({})
            input_path.write_text(json.dumps(payload), encoding="utf-8")
            engine = RagReEngine(
                self.config,
                pubmed_client=FakePubMed(),
                article_judge=FakeJudge(),
            )
            with self.assertRaisesRegex(ValueError, "invalid indexes"):
                engine.run_file(input_path)

            payload["agent_bundle"]["pathogen_candidates"] = []
            input_path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "outside the frozen candidate pool"):
                engine.run_file(input_path)

    def test_level_one_is_never_downgraded_and_abstention_stays_null(self):
        with tempfile.TemporaryDirectory() as temp:
            input_path = self._input(Path(temp))
            payload = json.loads(input_path.read_text(encoding="utf-8"))
            candidate = payload["agent_bundle"]["pathogen_candidates"][0]
            candidate["integrated_causative_level"] = "Level 1"
            input_path.write_text(json.dumps(payload), encoding="utf-8")
            result = RagReEngine(
                self.config,
                pubmed_client=FakePubMed(),
                article_judge=FakeJudge(),
            ).run_file(input_path)
        conclusion = result["rag_answer"]["level_reassessment_conclusion"][0]
        self.assertEqual(conclusion["recommended_level"], "Level 1")
        self.assertEqual(conclusion["decision"], "維持")

        row = result["candidates"][0]
        row["experiments"]["E0_OR_C_OR_A_AND_B"] = None
        rebuilt = RagReEngine(self.config, skip_literature=True)._assemble(
            patient_id="900",
            rows=[row],
            baseline_names=["Organism alpha"],
            errors=[],
            input_meta=result["input_meta"],
        )
        self.assertIsNone(
            rebuilt["rag_answer"]["level_reassessment_conclusion"][0][
                "is_infection_source"
            ]
        )

    def test_replay_rejects_evidence_generation_changes(self):
        with tempfile.TemporaryDirectory() as temp:
            original = RagReEngine(
                self.config,
                pubmed_client=FakePubMed(),
                article_judge=FakeJudge(),
            ).run_file(self._input(Path(temp)))
        changed = copy.deepcopy(self.config)
        changed["pubmed"]["publication_cutoff"] = "2025-12-31"
        with self.assertRaisesRegex(ValueError, "cannot change evidence-generation"):
            RagReEngine(changed, skip_literature=True).replay(original)

        wrong_pipeline = copy.deepcopy(original)
        wrong_pipeline["pipeline_fingerprint"] = "different"
        with self.assertRaisesRegex(ValueError, "different pipeline implementation"):
            RagReEngine(self.config, skip_literature=True).replay(wrong_pipeline)

    def test_replay_source_hash_ignores_runtime_file_location(self):
        with tempfile.TemporaryDirectory() as temp:
            original = RagReEngine(
                self.config,
                pubmed_client=FakePubMed(),
                article_judge=FakeJudge(),
            ).run_file(self._input(Path(temp)))
        first = copy.deepcopy(original)
        second = copy.deepcopy(original)
        first["_source_file"] = "C:/first/location.json"
        second["_source_file"] = "D:/second/location.json"
        engine = RagReEngine(self.config, skip_literature=True)
        self.assertEqual(
            engine.replay(first)["replayed_from_sha256"],
            engine.replay(second)["replayed_from_sha256"],
        )

    def test_prior_llm_appended_candidate_is_excluded_from_primary_pool(self):
        with tempfile.TemporaryDirectory() as temp:
            input_path = self._input(Path(temp))
            payload = json.loads(input_path.read_text(encoding="utf-8"))
            payload["agent_bundle"]["pathogen_candidates"].append(
                {
                    "organism_name": "Prior review organism",
                    "integrated_causative_level": "Level 4",
                    "rank_rule": "LLM_MISSED_CANDIDATE_REVIEW_APPENDED",
                    "applied_rules": ["LLM_MISSED_CANDIDATE_REVIEW_APPENDED"],
                    "key_evidence": {"manual_review_expansion": True},
                }
            )
            input_path.write_text(json.dumps(payload), encoding="utf-8")
            result = RagReEngine(self.config, skip_literature=True).run_file(input_path)
        self.assertEqual(result["input_meta"]["candidate_count"], 1)
        self.assertEqual(result["input_meta"]["excluded_prior_llm_expansion_count"], 1)
        self.assertNotIn(
            "Prior review organism",
            [row["organism_name"] for row in result["candidates"]],
        )


if __name__ == "__main__":
    unittest.main()

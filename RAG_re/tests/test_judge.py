import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from rag_re.judge import ArticleJudge


class FakeCompletions:
    def __init__(self, *, refusal=None):
        self.calls = []
        self.refusal = refusal

    def create(self, **kwargs):
        self.calls.append(kwargs)
        content = json.dumps(
            {
                "verdict": "support",
                "organism_scope": "exact_species",
                "human_clinical_evidence": True,
                "target_site_match": True,
                "evidence_span": "Human pneumonia evidence",
                "rationale": "Direct case evidence.",
            }
        )
        message = SimpleNamespace(content=content, refusal=self.refusal)
        return SimpleNamespace(
            choices=[SimpleNamespace(message=message, finish_reason="stop")],
            usage=None,
        )


def fake_judge(completions, cache_dir=None, max_tokens=450):
    judge = ArticleJudge.__new__(ArticleJudge)
    judge.config = {
        "model": "gpt-5.6-luna",
        "reasoning_effort": "low",
        "temperature": None,
        "max_output_tokens": max_tokens,
        "prompt_version": "article_support_v1",
    }
    judge.cache_dir = cache_dir
    judge.client = SimpleNamespace(chat=SimpleNamespace(completions=completions))
    return judge


class JudgeTests(unittest.TestCase):
    ARTICLE = {
        "pmid": "1",
        "title": "Example",
        "abstract": "Human pneumonia evidence in a clinical case.",
        "url": "https://pubmed.ncbi.nlm.nih.gov/1/",
    }

    def test_gpt5_request_omits_temperature_and_refusal_is_error(self):
        completions = FakeCompletions()
        result = fake_judge(completions).judge_one(
            "Example species", "pneumonia", self.ARTICLE
        )
        self.assertEqual(result["status"], "ok")
        self.assertEqual(completions.calls[0]["model"], "gpt-5.6-luna")
        self.assertEqual(completions.calls[0]["reasoning_effort"], "low")
        self.assertNotIn("temperature", completions.calls[0])

        with self.assertRaisesRegex(RuntimeError, "refused"):
            fake_judge(FakeCompletions(refusal="cannot comply")).judge_one(
                "Example species", "pneumonia", self.ARTICLE
            )

    def test_cache_key_changes_with_judgment_configuration(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first = fake_judge(FakeCompletions(), root, max_tokens=300)._cache_path(
                "Example species", "pneumonia", self.ARTICLE
            )
            second = fake_judge(FakeCompletions(), root, max_tokens=450)._cache_path(
                "Example species", "pneumonia", self.ARTICLE
            )
        self.assertNotEqual(first, second)


if __name__ == "__main__":
    unittest.main()

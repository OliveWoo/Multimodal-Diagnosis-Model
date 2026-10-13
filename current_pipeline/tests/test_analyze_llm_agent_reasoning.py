import os
import sys
import types
import unittest
from unittest.mock import patch

from tools import analyze_llm_agent


class _Usage:
    input_tokens = 10
    output_tokens = 5
    total_tokens = 15


class _Response:
    output_text = '{"status":"ok"}'
    output = []
    usage = _Usage()


class _Responses:
    def __init__(self, captured):
        self.captured = captured

    def create(self, **kwargs):
        self.captured.update(kwargs)
        return _Response()


class _OpenAI:
    captured = {}
    init_options = {}

    def __init__(self, *, api_key, **kwargs):
        self.api_key = api_key
        self.init_options = kwargs
        type(self).init_options = kwargs
        self.responses = _Responses(self.captured)


class AnalyzeReasoningTests(unittest.TestCase):
    def setUp(self) -> None:
        _OpenAI.captured = {}
        _OpenAI.init_options = {}

    def test_reasoning_effort_is_sent_to_responses_api(self) -> None:
        fake_module = types.SimpleNamespace(OpenAI=_OpenAI)
        with patch.dict(os.environ, {"OPENAI_API_KEY": "sk-test"}, clear=False), patch.dict(
            sys.modules, {"openai": fake_module}
        ):
            text, usage = analyze_llm_agent.send_to_llm(
                "prompt",
                model="gpt-6.1-sol",
                reasoning_effort="high",
            )

        self.assertEqual(text, '{"status":"ok"}')
        self.assertEqual(_OpenAI.captured["model"], "gpt-6.1-sol")
        self.assertEqual(_OpenAI.captured["reasoning"], {"effort": "high"})
        self.assertEqual(_OpenAI.init_options["max_retries"], 5)
        self.assertEqual(_OpenAI.init_options["timeout"], 900.0)
        self.assertEqual(usage["total_tokens"], 15)

    def test_reasoning_is_omitted_when_unspecified(self) -> None:
        fake_module = types.SimpleNamespace(OpenAI=_OpenAI)
        with patch.dict(os.environ, {"OPENAI_API_KEY": "sk-test"}, clear=False), patch.dict(
            sys.modules, {"openai": fake_module}
        ):
            analyze_llm_agent.send_to_llm("prompt", model="gpt-5.6-luna")

        self.assertNotIn("reasoning", _OpenAI.captured)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MODEL = "gpt-5"

# 預設搜尋 prompt 的目錄，順序越前優先度越高
PROMPT_ROOTS: tuple[Path, ...] = (
    REPO_ROOT / "prompts",
    REPO_ROOT / "agents" / "prompts",
)

# 明確指定任務代號與檔案路徑的對應關係（優先於目錄自動搜尋）
PROMPT_REGISTRY: dict[str, Path] = {
    # prompts/
    "admission_diagnosis": REPO_ROOT / "prompts" / "admission diagnosis_prompt.txt",
    "cbc": REPO_ROOT / "prompts" / "CBC_prompt.txt",
    "culture": REPO_ROOT / "prompts" / "culture_prompt.txt",
    "filmarray": REPO_ROOT / "prompts" / "filmarray_prompt.txt",
    "gm_test": REPO_ROOT / "prompts" / "gm test_prompt.txt",
    "image": REPO_ROOT / "prompts" / "image_prompt.txt",
    "other_lab": REPO_ROOT / "prompts" / "other lab_prompt.txt",
    # agents/prompts/
    "big_summary": REPO_ROOT / "agents" / "prompts" / "big_prompt.txt",
    "cbc_other_lab_agent": REPO_ROOT / "agents" / "prompts" / "CBC_Lab_underlying_prompt.txt",
    "filmarray_gmtest_agent": REPO_ROOT / "agents" / "prompts" / "filmarray_GMtest_prompt.txt",
    "image_agent": REPO_ROOT / "agents" / "prompts" / "image_agent_prompt.txt",
    "culture_agent": REPO_ROOT / "agents" / "prompts" / "culture_agent_prompt.txt",
    "mngs_to_specimen_agent": REPO_ROOT / "agents" / "prompts" / "mNGS_to_specimen_prompt.txt",
}


def _candidate_paths(task: str) -> Iterable[Path]:
    """依任務代號產生可能的檔案路徑，供後續檢查是否存在。"""
    for root in PROMPT_ROOTS:
        yield root / f"{task}.txt"
        yield root / f"{task}_prompt.txt"


def resolve_prompt_path(task: str) -> Path:
    """根據任務代號找到實際的 prompt 檔案路徑。"""
    if task in PROMPT_REGISTRY:
        path = PROMPT_REGISTRY[task]
        if path.exists():
            return path
        logging.warning("PROMPT_REGISTRY 配置的檔案不存在: %s", path)

    for candidate in _candidate_paths(task):
        if candidate.exists():
            return candidate

    raise FileNotFoundError(f"找不到符合任務 '{task}' 的 prompt 檔案。")


@lru_cache(maxsize=None)
def load_prompt_template(task: str) -> str:
    """讀取並快取指定任務的 prompt 內容。"""
    path = resolve_prompt_path(task)
    return path.read_text(encoding="utf-8")


def render_prompt(task: str, variables: Mapping[str, Any] | None = None) -> str:
    """用變數填入模板（若無變數則直接返回模板）。"""
    template = load_prompt_template(task)
    if not variables:
        return template
    try:
        return template.format(**variables)
    except KeyError as exc:  # 明確指出缺少的變數
        missing = exc.args[0]
        raise KeyError(f"prompt '{task}' 缺少變數: {missing}") from exc


def _extract_output_text(response: Any) -> str:
    """兼容 Responses API 或其他回傳格式，提取實際文本內容。"""
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


def send_to_llm(prompt: str, *, model: str = DEFAULT_MODEL) -> tuple[str, dict[str, int | None] | None]:
    """實際送出 LLM 請求並回傳結果與 token 使用量。"""
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY 環境變數未設定。")

    from openai import OpenAI  # type: ignore

    client = OpenAI(api_key=api_key)
    response = client.responses.create(
        model=model,
        temperature=1,
        input=prompt,
    )
    message = _extract_output_text(response)
    if not message:
        raise RuntimeError("LLM 回應沒有內容。")

    usage = getattr(response, "usage", None)
    usage_summary: dict[str, int | None] | None = None
    if usage is not None:
        usage_summary = {
            "prompt_tokens": getattr(usage, "input_tokens", None),
            "completion_tokens": getattr(usage, "output_tokens", None),
            "total_tokens": getattr(usage, "total_tokens", None),
        }
    return message, usage_summary


def run_llm_request(
    task: str,
    *,
    variables: Mapping[str, Any] | None = None,
    model: str = DEFAULT_MODEL,
) -> tuple[str, dict[str, int | None] | None]:
    """給定任務代號與變數，組合 prompt 並送出 LLM 請求。"""
    prompt = render_prompt(task, variables)
    return send_to_llm(prompt, model=model)


def run_llm_requests_parallel(
    requests: Sequence[tuple[str, Mapping[str, Any] | None]],
    *,
    model: str = DEFAULT_MODEL,
    max_workers: int = 4,
) -> list[tuple[str, dict[str, int | None] | None]]:
    """以 ThreadPoolExecutor 平行處理多組任務，結果順序與輸入一致。"""
    results: list[tuple[str, dict[str, int | None] | None]] = [("", None)] * len(requests)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {
            executor.submit(run_llm_request, task, variables=variables, model=model): idx
            for idx, (task, variables) in enumerate(requests)
        }
        for future in as_completed(future_map):
            idx = future_map[future]
            results[idx] = future.result()
    return results

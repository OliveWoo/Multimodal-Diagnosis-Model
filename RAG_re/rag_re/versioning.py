from __future__ import annotations

import hashlib
from pathlib import Path


PIPELINE_VERSION = "0.4.0"
EVALUATOR_VERSION = "0.3.0"

# Version 0.3.0 produced the frozen GPT-5.6-luna evidence used by the precision
# pruning experiment.  Version 0.4.0 only adds decision combinations and output
# semantics; retrieval, prompt, judge schema, B, and C evidence rules are unchanged.
REPLAY_COMPATIBLE_PIPELINE_FINGERPRINTS = {
    "609f5b957780e45c537d3ff3457bc7b536d4bc99fa2ad0084c2ca3d566ea2c7e",
}
_FINGERPRINT_FILES = (
    "versioning.py",
    "config.py",
    "input_adapter.py",
    "pubmed.py",
    "judge.py",
    "rules.py",
    "engine.py",
)


def pipeline_fingerprint() -> str:
    root = Path(__file__).resolve().parent
    digest = hashlib.sha256()
    digest.update(PIPELINE_VERSION.encode("utf-8"))
    for name in _FINGERPRINT_FILES:
        path = root / name
        digest.update(name.encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


PIPELINE_FINGERPRINT = pipeline_fingerprint()


def evaluator_fingerprint() -> str:
    root = Path(__file__).resolve().parent
    digest = hashlib.sha256()
    digest.update(EVALUATOR_VERSION.encode("utf-8"))
    for name in ("versioning.py", "evaluation.py", "cli.py"):
        digest.update(name.encode("utf-8"))
        digest.update((root / name).read_bytes())
    return digest.hexdigest()


EVALUATOR_FINGERPRINT = evaluator_fingerprint()

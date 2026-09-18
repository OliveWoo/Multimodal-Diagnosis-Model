from __future__ import annotations

import json
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools import mngs_common as mngs  # noqa: E402

DEFAULT_IMPOSSIBLE_RULE_PATH = REPO_ROOT / "rules" / "impossible_infection_sources.json"


@lru_cache(maxsize=8)
def load_impossible_infection_sources(path_text: str | None = None) -> dict[str, Any]:
    path = Path(path_text) if path_text else DEFAULT_IMPOSSIBLE_RULE_PATH
    if not path.exists():
        return {"version": "missing", "items": []}
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        return {"version": "invalid", "items": []}
    payload.setdefault("items", [])
    return payload


def _normalized_aliases(item: dict[str, Any]) -> set[str]:
    aliases = set()
    for value in [item.get("organism_name"), *(item.get("normalized_names") or []), *(item.get("aliases") or [])]:
        norm = mngs.normalize_organism_name(value)
        if norm:
            aliases.add(norm)
        text = str(value or "").strip().lower()
        if text:
            aliases.add(text.replace(" ", ""))
    return aliases


@lru_cache(maxsize=8)
def impossible_infection_source_index(path_text: str | None = None) -> dict[str, dict[str, Any]]:
    payload = load_impossible_infection_sources(path_text)
    index: dict[str, dict[str, Any]] = {}
    for item in payload.get("items") or []:
        if not isinstance(item, dict):
            continue
        for alias in _normalized_aliases(item):
            index[alias] = item
    return index


def impossible_infection_source_rule(value: Any, *, path_text: str | None = None) -> dict[str, Any] | None:
    norm = mngs.normalize_organism_name(value)
    if not norm:
        return None
    return impossible_infection_source_index(path_text).get(norm)


def is_impossible_infection_source(value: Any, *, path_text: str | None = None) -> bool:
    return impossible_infection_source_rule(value, path_text=path_text) is not None


def impossible_infection_source_names(path_text: str | None = None) -> list[str]:
    payload = load_impossible_infection_sources(path_text)
    names = []
    for item in payload.get("items") or []:
        if isinstance(item, dict) and item.get("organism_name"):
            names.append(str(item.get("organism_name")))
    return names

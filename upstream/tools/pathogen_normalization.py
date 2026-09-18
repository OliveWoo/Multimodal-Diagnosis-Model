"""Shared organism-name normalization for the mNGS pipeline.

The goal is to keep matching keys stable while allowing output names to stay
clinically readable. Rules live in rules/pathogen_aliases.json so ranked-input
selection, max scoring, merge/reconcile, and evaluation scripts do not drift.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
RULES_PATH = REPO_ROOT / "rules" / "pathogen_aliases.json"
NAME_SANITIZER = re.compile(r"[^a-z0-9]+")
NUMERIC_TRAILING_PAREN = re.compile(r"\s*\(\s*\d+(?:\.\d+)?\s*\)\s*$")
PARENTHETICAL = re.compile(r"\([^)]*\)")
RESISTANCE_LABELS = {
    "crkp": "CRKP phenotype",
    "cre": "CRE phenotype",
    "esbl": "ESBL phenotype",
    "mrsa": "MRSA phenotype",
    "vre": "VRE phenotype",
    "vrefm": "VRE phenotype",
}
GENERIC_TAXON_WORDS = {
    "complex",
    "evidence",
    "fungal",
    "group",
    "organism",
    "species",
    "sp",
    "spp",
}


@lru_cache(maxsize=1)
def rules_payload() -> dict[str, Any]:
    if not RULES_PATH.exists():
        return {
            "aliases": {},
            "approved_group_member_matches": {},
            "canonical_display_names": {},
            "representative_genus_groups": {},
        }
    payload = json.loads(RULES_PATH.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict JSON from {RULES_PATH}")
    return payload


def _alias_record(alias_key: str) -> dict[str, Any]:
    aliases = rules_payload().get("aliases")
    if not isinstance(aliases, dict):
        return {}
    record = aliases.get(alias_key)
    if isinstance(record, str):
        return {"canonical_key": record}
    if isinstance(record, dict):
        return record
    return {}


def alias_mapping() -> dict[str, str]:
    aliases = rules_payload().get("aliases")
    if not isinstance(aliases, dict):
        return {}
    result: dict[str, str] = {}
    for alias, record in aliases.items():
        if isinstance(record, str):
            result[alias] = record
        elif isinstance(record, dict):
            canonical = record.get("canonical_key")
            if canonical:
                result[alias] = str(canonical)
    return result


def display_mapping() -> dict[str, str]:
    displays = rules_payload().get("canonical_display_names")
    return dict(displays) if isinstance(displays, dict) else {}


def alias_expansions() -> dict[str, set[str]]:
    result: dict[str, set[str]] = {}
    mapping = alias_mapping()
    for alias, canonical in mapping.items():
        result.setdefault(alias, set()).add(canonical)
        result.setdefault(canonical, set()).add(alias)
    return result


def clean_display_text(value: Any) -> str:
    text = str(value or "").strip()
    text = NUMERIC_TRAILING_PAREN.sub("", text)
    return re.sub(r"\s+", " ", text).strip()


def _strip_qualifiers(text: str) -> str:
    text = NUMERIC_TRAILING_PAREN.sub(" ", text)
    text = PARENTHETICAL.sub(" ", text)
    text = re.sub(r"\b(subsp|ssp)\.?\b", "subsp", text, flags=re.I)
    for label in RESISTANCE_LABELS:
        candidate = re.sub(rf"\b{re.escape(label)}\b", " ", text, flags=re.I)
        if NAME_SANITIZER.sub("", candidate):
            text = candidate
    return text


def compact_key(value: Any, *, strip_qualifiers: bool = True) -> str:
    text = str(value or "").strip().lower()
    text = text.replace("β", "beta")
    if strip_qualifiers:
        text = _strip_qualifiers(text)
    else:
        text = NUMERIC_TRAILING_PAREN.sub(" ", text)
        text = re.sub(r"\b(subsp|ssp)\.?\b", "subsp", text, flags=re.I)
    return NAME_SANITIZER.sub("", text)


def raw_key(value: Any) -> str:
    return compact_key(value, strip_qualifiers=True)


def canonical_key(value: Any) -> str:
    key = raw_key(value)
    if not key:
        return ""
    record = _alias_record(key)
    canonical = record.get("canonical_key")
    return str(canonical) if canonical else key


def expanded_name_keys(value: Any) -> set[str]:
    key = raw_key(value)
    if not key:
        return set()
    canonical = canonical_key(value)
    expansions = alias_expansions()
    return {key, canonical, *expansions.get(key, set()), *expansions.get(canonical, set())}


def approved_group_member_mapping() -> dict[str, set[str]]:
    """Return explicit group/complex labels and their approved member species.

    These relationships are formal evaluation matches, not aliases. Keeping
    them separate preserves the original group-level laboratory label while
    allowing a validated member-species answer to count as a true positive.
    """
    rules = rules_payload().get("approved_group_member_matches")
    if not isinstance(rules, dict):
        return {}
    output: dict[str, set[str]] = {}
    for group, members in rules.items():
        group_key = raw_key(group)
        if not group_key or not isinstance(members, list):
            continue
        normalized_members = {canonical_key(member) for member in members if canonical_key(member)}
        if normalized_members:
            output[group_key] = normalized_members
    return output


def approved_group_member_match(left: Any, right: Any) -> bool:
    """Whether one name is an approved group/complex label for the other."""
    mapping = approved_group_member_mapping()
    left_raw = raw_key(left)
    right_raw = raw_key(right)
    left_canonical = canonical_key(left)
    right_canonical = canonical_key(right)
    return (
        bool(left_raw and right_canonical and right_canonical in mapping.get(left_raw, set()))
        or bool(right_raw and left_canonical and left_canonical in mapping.get(right_raw, set()))
    )


def approved_group_members(value: Any) -> set[str]:
    """Return explicitly approved member species for a group/complex label."""
    return set(approved_group_member_mapping().get(raw_key(value), set()))


def name_match_type(left: Any, right: Any, *, allow_genus_relaxed: bool = True) -> str:
    """Return the centralized evaluation relationship between two names."""
    if approved_group_member_match(left, right):
        return "approved_group_member_match"
    left_canonical = canonical_key(left)
    right_canonical = canonical_key(right)
    if left_canonical and left_canonical == right_canonical:
        return "exact_or_alias"
    if allow_genus_relaxed:
        left_genus = genus_name(left)
        right_genus = genus_name(right)
        if left_genus and left_genus == right_genus:
            return "genus_relaxed"
    return ""


def display_name(value: Any) -> str:
    canonical = canonical_key(value)
    if not canonical:
        return ""
    display = display_mapping().get(canonical)
    if display:
        return display
    return clean_display_text(value)


def genus_name(value: Any) -> str:
    text = display_name(value) or clean_display_text(value)
    if not text:
        return ""
    token = re.split(r"\s+", text.strip())[0]
    return compact_key(token, strip_qualifiers=True)


def resistance_note(value: Any) -> str | None:
    compact = compact_key(value, strip_qualifiers=False)
    lowered = str(value or "").lower()
    for label, note in RESISTANCE_LABELS.items():
        if re.search(rf"\b{re.escape(label)}\b", lowered, flags=re.I) or compact.endswith(label):
            return note
    return None


def alias_info(value: Any) -> tuple[str, str | None, str | None]:
    key = raw_key(value)
    if not key:
        return "", None, None
    record = _alias_record(key)
    canonical = str(record.get("canonical_key") or key)
    relationship = record.get("relationship")
    note = record.get("note")
    resistance = resistance_note(value)
    if resistance:
        relationship = "resistance_phenotype"
        note = resistance
    return canonical, str(relationship) if relationship else None, str(note) if note else None


def is_candida_or_generic_yeast(value: Any) -> bool:
    key = canonical_key(value)
    genus = genus_name(value)
    raw = raw_key(value)
    return key == "yeastfungalevidence" or genus in {"candida", "nakaseomyces"} or raw.startswith("candida")


def representative_genus_key(value: Any) -> str:
    genus = genus_name(value)
    groups = rules_payload().get("representative_genus_groups")
    if isinstance(groups, dict):
        mapped = groups.get(genus)
        if mapped:
            return str(mapped)
    return ""


def representative_group_key(value: Any) -> str:
    key, relationship, _ = alias_info(value)
    if not key:
        return "name:"
    genus_group = representative_genus_key(value)
    if key == "aspergillusspp" or genus_group == "aspergillus":
        return "genus:aspergillus"
    if is_candida_or_generic_yeast(value):
        return "yeast:candida_or_generic_yeast"
    if genus_group:
        return f"genus:{genus_group}"
    if relationship:
        return f"alias:{key}"
    return f"name:{key}"


def is_generic_representative_label(value: Any) -> bool:
    """Return True for broad labels that may be folded into a species.

    This is intentionally narrower than representative_group_key(). Different
    species in the same genus can have different clinical meaning and treatment
    implications, so they should not disappear solely because one genus-level
    representative was already retained.
    """
    key, relationship, _ = alias_info(value)
    raw = raw_key(value)
    text = clean_display_text(value).lower()
    tokens = [token.strip(".") for token in re.split(r"\s+", text) if token.strip(".")]
    if key == "yeastfungalevidence":
        return True
    if relationship in {"broad_morphology_label", "disease_or_test_label"}:
        return True
    if any(token in GENERIC_TAXON_WORDS for token in tokens):
        return True
    if raw.endswith("spp") or " species" in text or " spp" in text:
        return True
    if " group" in text or " complex" in text or " evidence" in text:
        return True
    return False


def representative_group_suppression_allowed(value: Any) -> bool:
    """Whether a same-representative-group hit may suppress this item.

    True aliases remain suppressible. Genus-level grouping is only suppressive
    for generic labels such as spp./group/complex/Yeast, not for a distinct
    species like Candida albicans, Aspergillus terreus, or Klebsiella variicola.
    """
    group = representative_group_key(value)
    if not group or group.startswith("name:"):
        return False
    if group.startswith("alias:"):
        return True
    return is_generic_representative_label(value)


def representative_specificity_score(value: Any) -> float:
    text = clean_display_text(value).lower()
    tokens = [token.strip(".") for token in re.split(r"\s+", text) if token.strip(".")]
    score = 0.0
    if any(token in GENERIC_TAXON_WORDS for token in tokens):
        score -= 150.0
    if any(word in text for word in (" group", " complex", " spp", "species", " evidence")):
        score -= 150.0
    if len(tokens) >= 2 and tokens[1] not in GENERIC_TAXON_WORDS:
        score += 120.0
    return score

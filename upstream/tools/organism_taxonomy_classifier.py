"""Answer-blind organism identity, traits, and clinical rule routing.

This module deliberately does not assign picked/review tiers. It centralizes
the stable organism-side metadata that downstream scorer, review, merge, and
RAG stages may consume together with patient-level evidence.
"""

from __future__ import annotations

import argparse
import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

from tools import pathogen_normalization as pathogen_names


REPO_ROOT = Path(__file__).resolve().parent.parent
RULES_PATH = REPO_ROOT / "rules" / "organism_taxonomy_rules.json"
FORBIDDEN_ANSWER_KEYS = {
    "answer",
    "answers",
    "answer_hit",
    "benchmark_answer",
    "ground_truth",
    "is_answer",
    "match_type",
}
GENERIC_LABEL_TERMS = {"sp", "spp", "species", "group", "complex", "evidence", "organism"}
BIOLOGICAL_CLASS_ALIASES = {
    "1bac": "bacterium",
    "bac": "bacterium",
    "bacteria": "bacterium",
    "bacterial": "bacterium",
    "bacterium": "bacterium",
    "2fungi": "fungus",
    "fungal": "fungus",
    "fungi": "fungus",
    "fungus": "fungus",
    "3virus": "virus",
    "viral": "virus",
    "virus": "virus",
    "parasite": "parasite",
    "parasitic": "parasite",
}


def _compact(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").strip().lower())


def _as_dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _unique_strings(*groups: Iterable[Any]) -> list[str]:
    output: list[str] = []
    seen: set[str] = set()
    for group in groups:
        for value in group:
            text = str(value or "").strip()
            if text and text not in seen:
                seen.add(text)
                output.append(text)
    return output


def _find_forbidden_rule_keys(value: Any, path: str = "$") -> list[str]:
    found: list[str] = []
    if isinstance(value, dict):
        for key, item in value.items():
            key_text = str(key)
            child_path = f"{path}.{key_text}"
            if key_text.lower() in FORBIDDEN_ANSWER_KEYS:
                found.append(child_path)
            found.extend(_find_forbidden_rule_keys(item, child_path))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            found.extend(_find_forbidden_rule_keys(item, f"{path}[{index}]"))
    return found


def _validate_rules(payload: dict[str, Any]) -> None:
    families = _as_dict(payload.get("families"))
    if "unmapped_or_uncertain" not in families:
        raise ValueError("organism taxonomy rules require unmapped_or_uncertain family")
    forbidden = _find_forbidden_rule_keys(payload)
    if forbidden:
        raise ValueError(f"answer-derived keys are forbidden in taxonomy rules: {forbidden}")
    profiles: list[dict[str, Any]] = []
    profiles.extend(value for value in _as_dict(payload.get("exact_profiles")).values() if isinstance(value, dict))
    profiles.extend(value for value in _as_dict(payload.get("genus_profiles")).values() if isinstance(value, dict))
    profiles.extend(value for value in _as_list(payload.get("key_prefix_profiles")) if isinstance(value, dict))
    for profile in profiles:
        family = str(profile.get("primary_rule_family") or "")
        if family not in families:
            raise ValueError(f"unknown primary_rule_family in taxonomy rules: {family!r}")
        unknown_secondary = [
            item for item in _as_list(profile.get("secondary_rule_families")) if str(item) not in families
        ]
        if unknown_secondary:
            raise ValueError(f"unknown secondary_rule_families: {unknown_secondary}")


@lru_cache(maxsize=1)
def rules_payload() -> dict[str, Any]:
    payload = json.loads(RULES_PATH.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise TypeError(f"Expected object in {RULES_PATH}")
    _validate_rules(payload)
    return payload


PROFILE_VERSION = str(rules_payload().get("profile_version") or "organism_taxonomy_profile_v1")
RULE_SCHEMA_VERSION = str(rules_payload().get("schema_version") or "organism_taxonomy_rules_v1")


def normalize_biological_class(value: Any) -> str:
    return BIOLOGICAL_CLASS_ALIASES.get(_compact(value), "unknown")


def inferred_taxonomic_rank(name: Any) -> str:
    text = pathogen_names.clean_display_text(name)
    if not text:
        return "unclassified"
    tokens = [token.strip(".,()[]").lower() for token in re.split(r"\s+", text) if token.strip(".,()[]")]
    if any(token in GENERIC_LABEL_TERMS for token in tokens):
        return "genus_or_group_label"
    if len(tokens) >= 2 and re.fullmatch(r"[a-z][a-z-]+", tokens[0]) and re.fullmatch(r"[a-z0-9][a-z0-9.-]+", tokens[1]):
        return "species_candidate"
    return "unclassified"


def _match_profile(name: Any) -> tuple[dict[str, Any], str, str]:
    payload = rules_payload()
    canonical = pathogen_names.canonical_key(name)
    exact = _as_dict(payload.get("exact_profiles"))
    profile = exact.get(canonical)
    if isinstance(profile, dict):
        return profile, "exact_species", "curated_exact_species"

    genus = pathogen_names.genus_name(name)
    genus_profiles = _as_dict(payload.get("genus_profiles"))
    profile = genus_profiles.get(genus)
    if isinstance(profile, dict):
        return profile, "genus_inherited", "curated_genus_profile"

    for raw_profile in _as_list(payload.get("key_prefix_profiles")):
        if not isinstance(raw_profile, dict):
            continue
        prefixes = [str(value) for value in _as_list(raw_profile.get("prefixes")) if value]
        if canonical and any(canonical.startswith(prefix) for prefix in prefixes):
            return raw_profile, "key_prefix_inferred", "curated_key_prefix_profile"
    return {}, "unmapped", "fallback"


def classify_organism(name: Any, *, biological_class: Any = None) -> dict[str, Any]:
    """Return routing metadata without deciding whether the organism caused disease."""
    input_name = str(name or "").strip()
    canonical = pathogen_names.canonical_key(input_name)
    display = pathogen_names.display_name(input_name) or input_name
    matched, mapping_status, rule_source = _match_profile(input_name)
    families = _as_dict(rules_payload().get("families"))
    primary_family = str(matched.get("primary_rule_family") or "unmapped_or_uncertain")
    family = _as_dict(families.get(primary_family))

    supplied_class = normalize_biological_class(biological_class)
    rule_class = normalize_biological_class(matched.get("biological_class"))
    resolved_class = rule_class if rule_class != "unknown" else supplied_class
    biological_class_conflict = bool(
        supplied_class != "unknown"
        and rule_class != "unknown"
        and supplied_class != rule_class
    )
    traits = _unique_strings(
        _as_list(family.get("default_traits")),
        _as_list(matched.get("clinical_traits")),
    )
    if primary_family == "unmapped_or_uncertain" and "taxonomy_review_required" not in traits:
        traits.append("taxonomy_review_required")
    confidence = {
        "exact_species": "high",
        "genus_inherited": "medium",
        "key_prefix_inferred": "medium",
        "unmapped": "low",
    }[mapping_status]
    rank = str(matched.get("taxonomic_rank") or inferred_taxonomic_rank(display))
    review_recommended = bool(
        matched.get("literature_review_recommended")
        or mapping_status == "unmapped"
        or resolved_class == "unknown"
    )
    return {
        "profile_version": PROFILE_VERSION,
        "input_name": input_name,
        "display_name": display,
        "canonical_key": canonical,
        "taxid": matched.get("taxid"),
        "taxonomic_rank": rank,
        "biological_class": resolved_class,
        "supplied_biological_class": supplied_class,
        "biological_class_conflict": biological_class_conflict,
        "primary_rule_family": primary_family,
        "secondary_rule_families": _unique_strings(_as_list(matched.get("secondary_rule_families"))),
        "clinical_traits": traits,
        "mapping_status": mapping_status,
        "classification_confidence": confidence,
        "rule_source": rule_source,
        "matched_rule_ids": _unique_strings([matched.get("rule_id")]),
        "needs_literature_review": review_recommended,
        "routing_only": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Inspect answer-blind organism taxonomy routing profiles.")
    parser.add_argument("organisms", nargs="+")
    parser.add_argument("--biological-class", default=None)
    args = parser.parse_args()
    output = [
        classify_organism(name, biological_class=args.biological_class)
        for name in args.organisms
    ]
    print(json.dumps(output, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

from dataclasses import dataclass
import re


@dataclass(frozen=True)
class Rule:
    modality: str
    contains_other_label: str
    keywords: tuple[str, ...]
    min_hits: int = 1


RULES = (
    Rule("abg", "abg_keywords", ("ph", "po2", "pco2", "hco3", "o2sat", "fio2"), min_hits=2),
    Rule("cbc", "cbc_keywords", ("wbc", "hb", "hct", "plt", "neutrophil"), min_hits=2),
    Rule(
        "smac",
        "smac_keywords",
        ("bun", "creatinine", "ast", "alt", "albumin", "crp", "sodium", "potassium", "na"),
        min_hits=2,
    ),
    Rule(
        "radiology",
        "radiology_keywords",
        (
            "imaging findings",
            "pleural effusion",
            "opacity",
            "ground glass",
            "ground-glass",
            "consolidation",
            "radiography",
        ),
    ),
)


def _keyword_pattern(keyword: str) -> re.Pattern[str]:
    escaped = re.escape(keyword.strip().lower())
    escaped = escaped.replace(r"\ ", r"[-\s]+").replace(r"\-", r"[-\s]+")
    prefix = r"(?<![a-z0-9])" if keyword[:1].isalnum() else ""
    suffix = r"(?![a-z0-9])" if keyword[-1:].isalnum() else ""
    return re.compile(prefix + escaped + suffix, re.IGNORECASE)


def _count_hits(text: str, keywords: tuple[str, ...]) -> int:
    hit_count = 0
    for keyword in keywords:
        if _keyword_pattern(keyword).search(text):
            hit_count += 1
    return hit_count


def classify_modalities(text_blocks: list[str]) -> tuple[list[str], dict]:
    text = "\n".join(text_blocks).lower()
    candidate_modalities: list[str] = []
    content_type_guess = {
        "contains_lab_table": False,
        "contains_abg_table": False,
        "contains_radiology_text": False,
        "contains_cbc_table": False,
        "contains_smac_table": False,
        "contains_other": [],
    }

    for rule in RULES:
        hits = _count_hits(text, rule.keywords)
        if hits >= rule.min_hits:
            candidate_modalities.append(rule.modality)
            if rule.modality == "abg":
                content_type_guess["contains_abg_table"] = True
                content_type_guess["contains_lab_table"] = True
            elif rule.modality == "cbc":
                content_type_guess["contains_cbc_table"] = True
                content_type_guess["contains_lab_table"] = True
            elif rule.modality == "smac":
                content_type_guess["contains_smac_table"] = True
                content_type_guess["contains_lab_table"] = True
            elif rule.modality == "radiology":
                content_type_guess["contains_radiology_text"] = True
            content_type_guess["contains_other"].append(rule.contains_other_label)

    deduped_modalities = sorted(set(candidate_modalities))
    content_type_guess["contains_other"] = sorted(set(content_type_guess["contains_other"]))
    return deduped_modalities, content_type_guess

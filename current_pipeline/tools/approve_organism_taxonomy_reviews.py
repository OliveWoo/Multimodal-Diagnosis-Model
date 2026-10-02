"""Compile human-approved organism reviews into an exact-name rule overlay."""

from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from tools import pathogen_normalization
from tools.organism_taxonomy_classifier import REPO_ROOT, _validate_rules, classify_organism, rules_payload


APPROVED_OVERLAY = REPO_ROOT / "rules" / "organism_taxonomy_reviewed_v1.json"


def _is_url(value: Any) -> bool:
    parsed = urlparse(str(value or ""))
    return parsed.scheme == "https" and bool(parsed.netloc)


def validate_approval(row: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    if row.get("decision") != "approved" or row.get("evidence_verified") is not True:
        raise ValueError("Only explicit, evidence-verified human approvals may become rules")
    review = row.get("manual_review") or {}
    if not isinstance(review, dict) or not str(review.get("reviewer") or "").strip():
        raise ValueError("manual_review.reviewer is required")
    try:
        date.fromisoformat(str(review.get("reviewed_at") or ""))
    except ValueError as exc:
        raise ValueError("manual_review.reviewed_at must be YYYY-MM-DD") from exc
    if not str(review.get("rationale") or "").strip():
        raise ValueError("manual_review.rationale is required")
    name = str(row.get("organism_name") or "").strip()
    key = pathogen_normalization.canonical_key(name)
    if not name or key != row.get("canonical_key"):
        raise ValueError("organism_name and canonical_key disagree")
    family = str(row.get("primary_rule_family") or "")
    if family not in rules_payload()["families"] or family == "unmapped_or_uncertain":
        raise ValueError("Approved profile needs an existing non-unmapped clinical family")
    if classify_organism(name)["mapping_status"] == "exact_species":
        raise ValueError("Organism already has an exact profile; change it through normal code review")
    refs = row.get("source_references")
    if not isinstance(refs, list) or not refs:
        raise ValueError("At least one verified clinical source is required")
    if not any(isinstance(ref, dict) and ref.get("role") == "clinical" and _is_url(ref.get("url")) for ref in refs):
        raise ValueError("A verified HTTPS clinical source is required; NCBI taxonomy alone is not enough")
    if any(
        not isinstance(ref, dict)
        or set(ref) != {"role", "url"}
        or ref["role"] not in {"identity", "clinical"}
        or not _is_url(ref["url"])
        for ref in refs
    ):
        raise ValueError("References may contain only role and a verified HTTPS URL")
    taxid = row.get("taxid")
    if isinstance(taxid, bool) or not str(taxid or "").isdigit() or int(taxid) <= 0:
        raise ValueError("A verified positive NCBI taxid is required")
    biological_class = str(row.get("biological_class") or "")
    if biological_class not in {"bacterium", "fungus", "virus", "parasite"}:
        raise ValueError("Unrecognized biological class")
    traits = row.get("clinical_traits")
    if not isinstance(traits, list) or any(not isinstance(item, str) or not item.strip() for item in traits):
        raise ValueError("clinical_traits must be a list of nonempty strings")
    rank = str(row.get("taxonomic_rank") or "")
    if rank not in {"species", "genus", "subspecies"}:
        raise ValueError("taxonomic_rank must be verified")
    return key, {
        "rule_id": f"TAX-REVIEWED-{key.upper()}",
        "taxid": int(taxid),
        "taxonomic_rank": rank,
        "biological_class": biological_class,
        "primary_rule_family": family,
        "clinical_traits": sorted(set(traits)),
        "literature_review_recommended": False,
        "review_provenance": {
            "reviewer": str(review["reviewer"]).strip(),
            "reviewed_at": str(review["reviewed_at"]),
            "rationale": str(review["rationale"]).strip(),
            "source_references": refs,
        },
    }


def compile_approvals(approvals_path: Path, output_path: Path, *, apply: bool = False) -> dict[str, Any]:
    rows = json.loads(approvals_path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError("Approvals file must be a JSON array")
    profiles: dict[str, dict[str, Any]] = {}
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("Every approval must be an object")
        key, profile = validate_approval(row)
        if key in profiles:
            raise ValueError(f"Duplicate approval for {key}")
        profiles[key] = profile
    overlay = {
        "schema_version": "organism_taxonomy_reviewed_v1",
        "purpose": "Human-approved organism metadata only; not patient-level causality or a picked decision.",
        "exact_profiles": dict(sorted(profiles.items())),
    }
    if apply and output_path.resolve() != APPROVED_OVERLAY.resolve():
        raise ValueError("--apply must target the central reviewed overlay")
    if output_path.exists():
        existing = json.loads(output_path.read_text(encoding="utf-8"))
        existing_profiles = existing.get("exact_profiles", {})
        for key, profile in existing_profiles.items():
            if key in profiles and profiles[key] != profile:
                raise ValueError(f"Existing approved rule differs: {key}")
            profiles.setdefault(key, profile)
        overlay["exact_profiles"] = dict(sorted(profiles.items()))
    combined = rules_payload().copy()
    combined["exact_profiles"] = {**combined["exact_profiles"], **overlay["exact_profiles"]}
    _validate_rules(combined)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(overlay, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"mode": "applied" if apply else "staged", "approved_profiles": len(profiles), "output": str(output_path.resolve())}


def main() -> None:
    parser = argparse.ArgumentParser(description="Compile manually verified organism reviews. Does not use LLM drafts as approvals.")
    parser.add_argument("approvals_json", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--apply", action="store_true", help="Allow writing rules/organism_taxonomy_reviewed_v1.json")
    args = parser.parse_args()
    if not args.apply and args.output.resolve() == APPROVED_OVERLAY.resolve():
        parser.error("Writing the live reviewed overlay requires --apply")
    print(json.dumps(compile_approvals(args.approvals_json, args.output, apply=args.apply), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

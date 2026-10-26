"""Freeze one policy bundle and replay KH51/A11/B11 without loading answers.

KH51 is partitioned into KH33 and KH18 in the review report. Existing releases
are never overwritten. The CMV evidence-axis proposal is a separate shadow.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "outputs/runs/2026-10-19_viral_policy_consistency_review_v5"
SOURCES = {}
SOURCE_MANIFESTS = {}
TAXONOMY_COMPLETION_SOURCE = None
EPISODE_IDENTITY_SOURCE = None
DEVELOPMENT_FREEZE_SOURCE = None

RULES = {
    "COMPACT_POLICY": "manual_style_multi_assay_entry_v2_20261020.json",
    "ANALYTICAL_POLICY": "test_aware_deterministic_shadow_v6_identity_resolved_20261020.json",
    "HISTORY_POLICY": "test_aware_history_route_v4c_combined.json",
    "CLINICAL_POLICY": "test_aware_clinical_deterministic_shadow_v1_20261020.json",
    "PROMOTION_POLICY": "test_aware_clinical_promotion_v9_pathogen_specific_convergence_20261020.json",
    "REPORTING_POLICY": "test_aware_unified_reporting_v7_respiratory_virus_possible_20261020.json",
    "IDENTITY_ALIAS_FILE": "pathogen_aliases.json",
    "TAXONOMY_BASE_FILE": "organism_taxonomy_rules_20261020.json",
    "TAXONOMY_GAP_FILE": "organism_taxonomy_gap_resolution_v1_20261020.json",
    "TAXONOMY_REVIEWED_FILE": "organism_taxonomy_reviewed_v1_20261020.json",
    "EPISODE_IDENTITY_FILE": "phenotype_episode_identity_v1.json",
    "PHENOTYPE_POLICY": "phenotype_candidate_policy_20261020.json",
}


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output: Path):
    if not SOURCE_MANIFESTS or TAXONOMY_COMPLETION_SOURCE is None or EPISODE_IDENTITY_SOURCE is None or DEVELOPMENT_FREEZE_SOURCE is None:
        raise ValueError("Use run_current_integrated_pipeline with an explicit private run config")
    output.mkdir(parents=True, exist_ok=True)
    bundle = output / "00_frozen_rules"
    registry_path = output / "frozen_registry.json"
    if not registry_path.exists():
        bundle.mkdir(exist_ok=True)
        paths = {}
        for role, name in RULES.items():
            target = bundle / name
            shutil.copy2(EPISODE_IDENTITY_SOURCE if role == "EPISODE_IDENTITY_FILE" else ROOT / "rules" / name, target)
            paths[role] = {"path": str(target.resolve()), "sha256": digest(target)}
        dependencies = {}
        pending = list(RULES.values())
        visited = set()
        while pending:
            name = pending.pop()
            if name in visited:
                continue
            visited.add(name)
            source = EPISODE_IDENTITY_SOURCE if name == RULES["EPISODE_IDENTITY_FILE"] else ROOT / "rules" / name
            payload = read(source)
            parent = payload.get("extends") or payload.get("base_policy")
            if parent:
                parent_source = (source.parent / parent).resolve()
                relative = parent_source.relative_to(ROOT / "rules")
                target = bundle / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(parent_source, target)
                dependencies[str(relative)] = {"path": str(target.resolve()), "sha256": digest(target)}
                pending.append(str(relative))
        completion_source = TAXONOMY_COMPLETION_SOURCE
        completion = read(completion_source)
        overlay = bundle / "staged_organism_taxonomy_auto_v1.json"
        shutil.copy2(completion["staged_overlay"], overlay)
        completion.update({"staged_overlay": str(overlay.resolve()), "staged_overlay_sha256": digest(overlay),
                           "frozen_from": str(completion_source), "frozen_from_sha256": digest(completion_source)})
        write(bundle / "taxonomy_completion_manifest.json", completion)
        # Source snapshots are recorded before any decisions or benchmark access.
        code = {str(p.relative_to(ROOT)): digest(p) for p in sorted((ROOT / "tools").glob("*.py"))}
        write(registry_path, {"schema_version": "viral_policy_consistency_review.v1", "answers_loaded": False,
              "rule_files": paths, "inherited_rule_files": dependencies, "code_sha256": code,
              "taxonomy_completion": str((bundle / "taxonomy_completion_manifest.json").resolve()),
              "taxonomy_completion_sha256": digest(bundle / "taxonomy_completion_manifest.json"),
              "taxonomy_overlay_sha256": digest(overlay), "source_manifests": {c: str(p) for c, p in SOURCE_MANIFESTS.items()},
              "cmv_policy_changed": False, "pjp_threshold_changed": False})
    registry = read(registry_path)
    for item in registry["rule_files"].values():
        if digest(Path(item["path"])) != item["sha256"]:
            raise ValueError("Frozen rule was modified")
    for item in registry.get("inherited_rule_files", {}).values():
        if digest(Path(item["path"])) != item["sha256"]:
            raise ValueError("Frozen inherited rule was modified")
    completion_path = Path(registry["taxonomy_completion"])
    completion = read(completion_path)
    os.environ["MNGS_TAXONOMY_AUTO_OVERLAY"] = completion["staged_overlay"]
    # Import after the overlay is selected, then bind every central taxonomy layer.
    from tools import run_external_integrated_release_v2 as runner
    from tools import organism_taxonomy_classifier as taxonomy
    from tools import candidate_phenotype_context as phenotype
    for role, item in registry["rule_files"].items():
        setattr(runner.kh_release, role, Path(item["path"]))
    taxonomy.RULES_PATH = Path(registry["rule_files"]["TAXONOMY_BASE_FILE"]["path"])
    taxonomy.RULE_OVERLAY_PATHS = (
        Path(registry["rule_files"]["TAXONOMY_GAP_FILE"]["path"]),
        Path(completion["staged_overlay"]),
        Path(registry["rule_files"]["TAXONOMY_REVIEWED_FILE"]["path"]),
    )
    taxonomy.rules_payload.cache_clear()
    phenotype.POLICY_PATH = Path(registry["rule_files"]["PHENOTYPE_POLICY"]["path"])
    phenotype.policy.cache_clear()
    # Resolve both inheritance formats before the expensive replay starts.
    runner.history.read_json(runner.kh_release.HISTORY_POLICY)
    runner.reporting.load_policy(runner.kh_release.REPORTING_POLICY)
    phenotype.policy()
    manifests = {}
    for cohort, source_name in SOURCES.items():
        target = output / cohort
        old = read(SOURCE_MANIFESTS[cohort])
        data = old["input"]
        if not (target / "run_manifest.json").exists():
            if target.exists() and any(target.iterdir()):
                raise FileExistsError(f"Incomplete run requires inspection: {target}")
            print(f"START {cohort}", flush=True)
            runner.run(target, inventory_root=Path(data["inventory_root"]),
                patient_root=Path(data["patient_root"]), phenotype_root=Path(data["phenotype_root"]),
                source_summary_suffix=data["source_summary_suffix"],
                development_freeze=DEVELOPMENT_FREEZE_SOURCE,
                upstream_extraction_model=old.get("upstream_extraction_model", "not_declared"),
                compatibility_target="source-bound CMV public replay, explicit private inputs",
                analytical_policy=runner.kh_release.ANALYTICAL_POLICY,
                promotion_policy=runner.kh_release.PROMOTION_POLICY,
                reporting_policy=runner.kh_release.REPORTING_POLICY,
                taxonomy_completion_manifest=completion_path,
                mngs_only=data.get("evidence_mode") == "mngs_only_optional_clinical_absent",
                predecision_evidence_graph=bool(data.get("predecision_evidence_graph", False)))
        manifest = read(target / "run_manifest.json")
        manifests[cohort] = manifest
        print(f"DONE {cohort}: {manifest['patient_count']} patients, {manifest['candidate_count']} candidates", flush=True)
    component_hashes = {}
    for cohort, manifest in manifests.items():
        component_hashes[cohort] = {x["role"]: x["sha256"] for x in manifest["components"]}
    roles = set.intersection(*(set(v) for v in component_hashes.values()))
    mismatches = {role: {c: h[role] for c, h in component_hashes.items()}
                  for role in roles if len({h[role] for h in component_hashes.values()}) != 1}
    for item in registry["rule_files"].values():
        if digest(Path(item["path"])) != item["sha256"]:
            raise ValueError("Frozen rule changed during replay")
    summary = {"answers_loaded": False, "passed": not mismatches, "component_hash_mismatches": mismatches,
               "run_roots": {c: str((output / c).resolve()) for c in SOURCES},
               "patient_counts": {c: m["patient_count"] for c, m in manifests.items()},
               "candidate_counts": {c: m["candidate_count"] for c, m in manifests.items()},
               "shared_component_hashes": next(iter(component_hashes.values()))}
    write(output / "consistency_summary.json", summary)
    if mismatches:
        raise ValueError(mismatches)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run(args.output_root), ensure_ascii=False, indent=2))

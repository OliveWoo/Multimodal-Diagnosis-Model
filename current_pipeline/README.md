# Current pipeline: publication snapshot 2026-10-13

This directory contains the source, policies, and offline tests for the multi-assay
DNA/RNA workflow and the legacy compatibility workflow. It contains no patient
inputs, reference answers, run outputs, or credentials.

## Install

From the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r current_pipeline\requirements-dev.txt
cd current_pipeline
python -m pytest tests -q
```

The PowerShell launchers use the Python interpreter on PATH. Activate the virtual
environment first. Paid model calls need an existing API credential or an
authenticated Codex CLI, depending on the selected backend; scripts never prompt
for an API key.

## Workflow

1. Standardize clinical inputs; preserve source references in audit sidecars.
2. Run five clinical evidence agents with explicit missing-data handling.
3. Build the per-test mNGS inventory without pooling reads across tests.
4. Detect unmapped names and resolve NCBI identity plus PubMed literature.
5. Run identity and clinical-family reviewers, then a deterministic validator.
6. Verify every reviewed name is Auto-approved, Provisional, or Rejected.
7. Use a run-scoped taxonomy overlay for deterministic scoring and reporting.
8. Keep OBER review and historical replay as explicitly selected workflows.
9. Freeze outputs before running the separate evaluator on private answers.

Clinical evidence, history, phenotype, and resistance data are optional. Missing
data remain unknown; they do not become negative results. The mNGS-only mode
explicitly creates unknown clinical inputs. Stronger claims still require the
appropriate supporting evidence.

## Public entry points

| Entry | Purpose |
| --- | --- |
| `scripts/run_small_agents.ps1` | Rebuild an answer-blind workspace, validate prompts, preflight, and resume five agents |
| `tools.build_multi_assay_mngs_evidence` | Build per-event and per-test DNA/RNA inventory |
| `scripts/run_taxonomy_complete_pipeline.ps1` | Resolve all unmapped names before allowing the integrated scorer |
| `tools.run_kh_legacy_pipeline_finalized_v1` | Legacy scorer with a caller-supplied release contract |
| `tools.run_two_hospital_22_legacy_finalized_compat_v1` | Deterministic representative-DNA legacy compatibility |
| `tools.run_historical_ober_r5_replay` | Historical OBER evidence review plus fixed acceptance conditions |
| `tools.select_ober_sol_escalations` | Select difficult Luna reviews for a second Sol review |
| `tools.merge_ober_two_stage_adjudications` | Revalidate and merge both reviewers under deterministic boundaries |
| `tools.build_mandatory_primary_reporting_shadow` | Optional best-available primary reporting; not a certainty upgrade |

See [the update guide](../docs/PIPELINE_UPDATE_20261013.md) for commands and workflow
boundaries. The repository-root `run_full_pipeline.py` remains the historical full
workflow; it is not silently replaced by the new DNA/RNA scorer.

## Reproducibility

The public launcher uses scorer v6, history v4C, clinical scorer v1, calendar-window
promotion v9, and reporting v7. Additional reporting, candidate-completeness, and
mandatory-primary policies are explicit research alternatives, not simultaneous
defaults. No new clinical performance claim is made by this source publication.

`rules/current_pipeline_policy_bundle_v6_v9_calendar_v7_frozen_20261013.json`
pins the public default component bytes and IDs. The promotion policy has a new ID
because its inclusive +/-2-calendar-date window is not the old +/-48-hour
reference. The older development manifest has not been silently rewritten.

```powershell
python -m tools.validate_frozen_policy_bundle
```

`PUBLIC_SOURCE_MANIFEST.csv` records the published source files and SHA-256 hashes.
A source update is not proof that historical private results can be reproduced
without their matching input and answer versions.

## Data safety

Do not commit patient JSON, Excel workbooks, phenotype tables, benchmark answers,
prompts containing patient records, raw model responses, logs, or API keys.
Cohort contracts include dataset layout and anonymous selection metadata only;
their default local run paths are examples, not shipped input files.

This is a research pipeline for candidate review, not an autonomous clinical
diagnosis system.

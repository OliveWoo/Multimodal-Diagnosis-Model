param(
    [Parameter(Mandatory = $true)][string]$InventoryRoot,
    [Parameter(Mandatory = $true)][string]$TaxonomyRoot,
    [Parameter(Mandatory = $true)][string]$OutputRoot,
    [string]$PatientRoot,
    [string]$PhenotypeRoot,
    [string]$SourceSummarySuffix,
    [string]$UpstreamExtractionModel = "not_declared",
    [ValidateSet("codex-cli", "openai-api")]
    [string]$Backend = "codex-cli",
    [string]$IdentityModel = "gpt-6-luna",
    [string]$ClinicalModel = "gpt-6.1-sol",
    [ValidateSet("low", "medium", "high", "xhigh", "max")]
    [string]$TaxonomyReasoningEffort = "high",
    [int]$CallTimeoutSeconds = 900,
    [switch]$TaxonomyOnly,
    [switch]$MngsOnly
)

$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent $PSScriptRoot)
$InventoryRoot = [IO.Path]::GetFullPath($InventoryRoot)
$TaxonomyRoot = [IO.Path]::GetFullPath($TaxonomyRoot)
$OutputRoot = [IO.Path]::GetFullPath($OutputRoot)

python -B -X utf8 -m tools.validate_frozen_policy_bundle
if ($LASTEXITCODE -ne 0) { throw "Publication policy fingerprint validation failed." }

$baselineRoot = Join-Path $TaxonomyRoot "00_baseline_compact"
$initialQueueRoot = Join-Path $TaxonomyRoot "01_queue_initial"
$ncbiRoot = Join-Path $TaxonomyRoot "02_ncbi"
$queueRoot = Join-Path $TaxonomyRoot "03_queue_resolved"
$literatureRoot = Join-Path $TaxonomyRoot "04_literature"
$decisionRoot = Join-Path $TaxonomyRoot "05_auto_decision"
$auditRoot = Join-Path $TaxonomyRoot "06_route_audit"
$completionManifest = Join-Path $TaxonomyRoot "taxonomy_completion_manifest.json"
$baselineCsv = Join-Path $baselineRoot "all_decisions.csv"
$initialQueue = Join-Path $initialQueueRoot "organism_review_queue.jsonl"
$ncbiInput = Join-Path $initialQueueRoot "ncbi_lookup_input.csv"
$ncbiCsv = Join-Path $ncbiRoot "ncbi_taxonomy_resolution.csv"
$queue = Join-Path $queueRoot "organism_review_queue.jsonl"
$literature = Join-Path $literatureRoot "pubmed_literature_pack.json"
$overlay = Join-Path $decisionRoot "staged_organism_taxonomy_auto_v1.json"
$decisionSummary = Join-Path $decisionRoot "run_summary.json"
$auditSummary = Join-Path $auditRoot "audit_summary.json"

foreach ($path in @($InventoryRoot)) {
    if (-not (Test-Path -LiteralPath $path -PathType Container)) {
        throw "Required input is missing: $path"
    }
}
if (-not $TaxonomyOnly -and -not $MngsOnly) {
    if (-not $PatientRoot -or -not $SourceSummarySuffix) {
        throw "Provide PatientRoot and SourceSummarySuffix, or use MngsOnly for explicitly missing clinical data."
    }
    if (-not (Test-Path -LiteralPath $PatientRoot -PathType Container)) {
        throw "PatientRoot is missing: $PatientRoot"
    }
    $PatientRoot = [IO.Path]::GetFullPath($PatientRoot)
}
if ($PhenotypeRoot) {
    if (-not (Test-Path -LiteralPath $PhenotypeRoot -PathType Container)) {
        throw "PhenotypeRoot is missing: $PhenotypeRoot"
    }
    $PhenotypeRoot = [IO.Path]::GetFullPath($PhenotypeRoot)
}
# Refuse to reuse taxonomy caches for different inventory bytes or model settings.
$fingerprints = @(Get-ChildItem -LiteralPath $InventoryRoot -Recurse -File |
    Where-Object { $_.Extension -in @(".csv", ".json", ".jsonl") } |
    Sort-Object FullName | ForEach-Object {
        [ordered]@{ path = $_.FullName; sha256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash }
    })
if ($fingerprints.Count -eq 0) { throw "InventoryRoot contains no inventory records." }
$runConfig = [ordered]@{
    inventory = $fingerprints; backend = $Backend; identity_model = $IdentityModel
    clinical_model = $ClinicalModel; reasoning_effort = $TaxonomyReasoningEffort
}
$configJson = $runConfig | ConvertTo-Json -Depth 6 -Compress
$configPath = Join-Path $TaxonomyRoot "public_run_config.json"
if (Test-Path -LiteralPath $configPath) {
    $saved = Get-Content -Raw -Encoding utf8 -LiteralPath $configPath | ConvertFrom-Json
    if (($saved | ConvertTo-Json -Depth 6 -Compress) -ne $configJson) {
        throw "TaxonomyRoot belongs to different inputs/models. Use a new TaxonomyRoot."
    }
} elseif ((Test-Path -LiteralPath $TaxonomyRoot) -and (Get-ChildItem -LiteralPath $TaxonomyRoot -Force | Select-Object -First 1)) {
    throw "TaxonomyRoot has no matching run configuration. Use a fresh directory."
}
New-Item -ItemType Directory -Force -Path $TaxonomyRoot | Out-Null
$configJson | Set-Content -Encoding utf8 -LiteralPath $configPath

if (-not (Test-Path -LiteralPath $baselineCsv)) {
    Write-Host "[1/9] Building the preliminary taxonomy mapping from the inventory."
    python -B -X utf8 -m tools.build_compact_multi_assay_entry_shadow `
        $InventoryRoot $baselineRoot `
        --policy "rules\manual_style_multi_assay_entry_v2.json"
    if ($LASTEXITCODE -ne 0) { throw "Preliminary taxonomy mapping failed." }
} else {
    Write-Host "[1/9] Reusing preliminary taxonomy mapping: $baselineCsv"
}

if (-not (Test-Path -LiteralPath $initialQueue)) {
    Write-Host "[2/9] Detecting the unique unmapped organism names."
    python -B -X utf8 -m tools.build_organism_taxonomy_review_queue `
        $baselineCsv --output-dir $initialQueueRoot
    if ($LASTEXITCODE -ne 0) { throw "Initial taxonomy queue build failed." }
} else {
    Write-Host "[2/9] Reusing initial taxonomy queue: $initialQueue"
}

if (-not (Test-Path -LiteralPath $ncbiCsv)) {
    Write-Host "[3/9] Resolving canonical name, taxid, and lineage against NCBI Taxonomy."
    python -B -X utf8 -m tools.resolve_taxonomy_gaps_ncbi `
        $ncbiInput --output-dir $ncbiRoot
    if ($LASTEXITCODE -ne 0) { throw "NCBI taxonomy resolution failed. Rerun to resume from cache." }
} else {
    Write-Host "[3/9] Reusing NCBI taxonomy resolution: $ncbiCsv"
}

if (-not (Test-Path -LiteralPath $queue)) {
    Write-Host "[4/9] Rebuilding the answer-blind queue with NCBI identity fields."
    python -B -X utf8 -m tools.build_organism_taxonomy_review_queue `
        $baselineCsv --output-dir $queueRoot --ncbi-resolution-csv $ncbiCsv
    if ($LASTEXITCODE -ne 0) { throw "Resolved taxonomy queue build failed." }
} else {
    Write-Host "[4/9] Reusing resolved taxonomy queue: $queue"
}

if (-not (Test-Path -LiteralPath $literature)) {
    Write-Host "[5/9] Retrieving organism-level PubMed evidence."
    New-Item -ItemType Directory -Force -Path $literatureRoot | Out-Null
    python -B -X utf8 -m tools.retrieve_organism_taxonomy_literature `
        $queue --output $literature --cache-dir (Join-Path $literatureRoot "cache")
    if ($LASTEXITCODE -ne 0) { throw "Literature retrieval failed. Rerun to resume from cache." }
} else {
    Write-Host "[5/9] Reusing literature pack: $literature"
}

$decisionComplete = $false
if (Test-Path -LiteralPath $decisionSummary) {
    $savedDecisionSummary = Get-Content -Raw -Encoding utf8 -LiteralPath $decisionSummary | ConvertFrom-Json
    $decisionComplete = (
        $savedDecisionSummary.completion_status -eq "complete" -and
        [int]$savedDecisionSummary.failures -eq 0 -and
        [int]$savedDecisionSummary.remaining_eligible_not_selected -eq 0
    )
}
if (-not $decisionComplete) {
    Write-Host "[6/9] Running/resuming LLM identity and clinical-family reviewers."
    $decisionArgs = @(
        "-B", "-X", "utf8", "-m", "tools.auto_decide_organism_taxonomy",
        "decide", $queue,
        "--literature-json", $literature,
        "--output-dir", $decisionRoot,
        "--backend", $Backend,
        "--identity-model", $IdentityModel,
        "--clinical-model", $ClinicalModel,
        "--reasoning-effort", $TaxonomyReasoningEffort,
        "--call-timeout-seconds", $CallTimeoutSeconds,
        "--execute"
    )
    & python @decisionArgs
    if ($LASTEXITCODE -ne 0) { throw "LLM taxonomy review failed. Rerun to resume completed organisms only." }
} else {
    Write-Host "[6/9] Reusing complete LLM taxonomy decisions: $decisionSummary"
}

$savedDecisionSummary = Get-Content -Raw -Encoding utf8 -LiteralPath $decisionSummary | ConvertFrom-Json
if (
    $savedDecisionSummary.completion_status -ne "complete" -or
    [int]$savedDecisionSummary.failures -ne 0 -or
    [int]$savedDecisionSummary.remaining_eligible_not_selected -ne 0
) {
    throw "Taxonomy review is incomplete. Rerun this script; scorer execution is blocked."
}

Write-Host "[7/9] Running deterministic route-impact validation."
python -B -X utf8 -m tools.audit_reviewed_taxonomy_overlay `
    $overlay $InventoryRoot $baselineRoot `
    --output-dir $auditRoot `
    --policy "rules\manual_style_multi_assay_entry_v2.json"
if ($LASTEXITCODE -ne 0) { throw "Taxonomy route-impact audit failed." }

Write-Host "[8/9] Verifying complete Auto-approved / Provisional / Rejected coverage."
python -B -X utf8 -m tools.validate_taxonomy_review_completion `
    --baseline-csv $baselineCsv `
    --queue $queue `
    --decision-dir $decisionRoot `
    --staged-overlay $overlay `
    --audit-summary $auditSummary `
    --output $completionManifest
if ($LASTEXITCODE -ne 0) { throw "Taxonomy completion gate failed; scorer execution is blocked." }

if ($TaxonomyOnly) {
    Write-Host "Taxonomy-only run complete: $completionManifest"
    return
}

if ((Test-Path -LiteralPath $OutputRoot) -and (Get-ChildItem -LiteralPath $OutputRoot -Force | Select-Object -First 1)) {
    $freeze = Join-Path $OutputRoot "answer_blind_freeze.json"
    if (Test-Path -LiteralPath $freeze) {
        $freezePayload = Get-Content -Raw -Encoding utf8 -LiteralPath $freeze | ConvertFrom-Json
        if ($freezePayload.passed -eq $true) {
            Write-Host "[9/9] Integrated release already completed: $OutputRoot"
            return
        }
    }
    throw "OutputRoot is non-empty but not a completed release: $OutputRoot"
}

Write-Host "[9/9] Running the deterministic scorer with the validated run-scoped taxonomy overlay."
$previousOverlay = $env:MNGS_TAXONOMY_AUTO_OVERLAY
try {
    $env:MNGS_TAXONOMY_AUTO_OVERLAY = $overlay
    if (-not $PhenotypeRoot) {
        $PhenotypeRoot = Join-Path $TaxonomyRoot "empty_optional_phenotype"
        New-Item -ItemType Directory -Force -Path $PhenotypeRoot | Out-Null
    }
    $releaseArgs = @(
        "-B", "-X", "utf8", "-m", "tools.run_external_integrated_release_v2",
        "--output-root", $OutputRoot, "--inventory-root", $InventoryRoot,
        "--phenotype-root", $PhenotypeRoot,
        "--upstream-extraction-model", $UpstreamExtractionModel,
        "--compatibility-target", "public current v6 + history v4C + clinical v1 + calendar promotion v9 + reporting v7",
        "--analytical-policy", "rules\test_aware_deterministic_shadow_v6_identity_resolved.json",
        "--promotion-policy", "rules\test_aware_clinical_promotion_v9_calendar_dates_v1.json",
        "--reporting-policy", "rules\test_aware_unified_reporting_v7_respiratory_virus_possible.json",
        "--taxonomy-completion-manifest", $completionManifest
    )
    if ($MngsOnly) {
        $releaseArgs += "--mngs-only"
    } else {
        $releaseArgs += @("--patient-root", $PatientRoot, "--source-summary-suffix", $SourceSummarySuffix)
    }
    & python @releaseArgs
    if ($LASTEXITCODE -ne 0) { throw "Integrated deterministic pipeline failed." }
}
finally {
    $env:MNGS_TAXONOMY_AUTO_OVERLAY = $previousOverlay
}

Write-Host "Complete taxonomy-aware pipeline finished: $OutputRoot"

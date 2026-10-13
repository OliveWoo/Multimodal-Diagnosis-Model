param(
    [Parameter(Mandatory = $true)][string]$SourceRoot,
    [Parameter(Mandatory = $true)][string]$RunRoot,
    [string]$Model = "gpt-5.6-luna",
    [ValidateSet("low", "medium", "high", "xhigh", "max")]
    [string]$ReasoningEffort = "high",
    [int]$MaxWorkers = 4,
    [int]$PreflightPatient = 0,
    [switch]$PreflightOnly,
    [switch]$SkipPreflight
)

$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent $PSScriptRoot)
$SourceRoot = [IO.Path]::GetFullPath($SourceRoot)
$RunRoot = [IO.Path]::GetFullPath($RunRoot)

if ($PreflightOnly -and $SkipPreflight) {
    throw "PreflightOnly and SkipPreflight cannot be used together."
}

$primaryApiKey = if ($env:OPENAI_API_KEY) { $env:OPENAI_API_KEY.Trim() } else { $null }
$legacyApiKey = if ($env:OPENAI_APIKEY) { $env:OPENAI_APIKEY.Trim() } else { $null }
if ($primaryApiKey -and $legacyApiKey -and $primaryApiKey -ne $legacyApiKey) {
    throw "OPENAI_API_KEY and OPENAI_APIKEY contain different values. Clear the stale variable before running."
}
$apiKey = if ($primaryApiKey) { $primaryApiKey } else { $legacyApiKey }
if (-not $apiKey -or $apiKey -notmatch '^sk-\S+$') {
    throw "OPENAI_API_KEY is missing or malformed in this PowerShell session. This script never prompts for a key."
}
$env:OPENAI_API_KEY = $apiKey
Write-Host "Using the API credential already present in this session."

$summaryPath = Join-Path $RunRoot "workspace_summary.json"
if (-not (Test-Path -LiteralPath $summaryPath)) {
    if ((Test-Path -LiteralPath $RunRoot) -and (Get-ChildItem -LiteralPath $RunRoot -Force | Select-Object -First 1)) {
        throw "RunRoot is non-empty but workspace_summary.json is missing: $RunRoot"
    }
    Write-Host "Building a fresh answer-blind workspace with source-aware image deduplication."
    python -B -X utf8 -m tools.build_44_patient_small_agent_workspace `
        --source-root $SourceRoot `
        --output-root $RunRoot
    if ($LASTEXITCODE -ne 0) {
        throw "Small-agent workspace build failed."
    }
}

$summary = Get-Content -Raw -Encoding utf8 -LiteralPath $summaryPath | ConvertFrom-Json
if ($summary.patient_count -lt 1 -or $summary.ground_truth_firewall -ne "passed") {
    throw "Workspace validation failed: expected patients and a passed ground-truth firewall."
}
$runConfig = [ordered]@{ source_root = $SourceRoot; model = $Model; reasoning_effort = $ReasoningEffort }
$configJson = $runConfig | ConvertTo-Json -Compress
$configPath = Join-Path $RunRoot "public_small_agent_run_config.json"
if (Test-Path -LiteralPath $configPath) {
    $saved = Get-Content -Raw -Encoding utf8 -LiteralPath $configPath | ConvertFrom-Json
    if (($saved | ConvertTo-Json -Compress) -ne $configJson) {
        throw "This workspace belongs to a different model or source. Use a fresh RunRoot."
    }
} elseif (Get-ChildItem -LiteralPath (Join-Path $RunRoot "patients") -Recurse -File -Filter "*_agent.json" | Select-Object -First 1) {
    throw "Existing model outputs lack a matching run configuration. Use a fresh RunRoot."
}
$configJson | Set-Content -Encoding utf8 -LiteralPath $configPath
$patientCount = [int]$summary.patient_count
$expectedOutputs = $patientCount * 5
if ($PreflightPatient -eq 0) { $PreflightPatient = [int]$summary.patient_ids[0] }

python -B -X utf8 -m tools.validate_44_patient_small_agent_workspace `
    --run-root $RunRoot `
    --model $Model `
    --max-input-tokens 50000 `
    --expected-patients $patientCount
if ($LASTEXITCODE -ne 0) {
    throw "Small-agent workspace validation failed. See $RunRoot\workspace_validation.json"
}

$InputRoot = Join-Path $RunRoot "patients"
$preflightDir = Join-Path $InputRoot "NGS_patient_${PreflightPatient}_json"
$preflightOutput = Join-Path $preflightDir "agent_outputs\NGS_patient_${PreflightPatient}_cbc_other_lab_agent.json"
if (-not $SkipPreflight) {
    if (-not (Test-Path -LiteralPath $preflightDir)) {
        throw "Preflight patient directory is missing: $preflightDir"
    }
    if (-not (Test-Path -LiteralPath $preflightOutput)) {
        Write-Host "Preflight: patient $PreflightPatient, cbc_other_lab, model=$Model, reasoning=$ReasoningEffort"
        python -B -X utf8 -m tools.analyze_llm_agent $preflightDir `
            --model $Model `
            --reasoning-effort $ReasoningEffort `
            --categories cbc_other_lab
        if ($LASTEXITCODE -ne 0) {
            throw "Small-agent API preflight process failed."
        }
    }
    Get-Content -Raw -Encoding utf8 -LiteralPath $preflightOutput | ConvertFrom-Json | Out-Null
    Write-Host "Preflight passed: $preflightOutput"
    if ($PreflightOnly) {
        return
    }
}

$existingCount = @(
    Get-ChildItem -LiteralPath $InputRoot -Recurse -File -Filter "*_agent.json" -ErrorAction SilentlyContinue
).Count
Write-Host "Running/resuming $patientCount patients with $Model ($ReasoningEffort). Existing output files: $existingCount/$expectedOutputs"
python -B -X utf8 -m tools.llm_scheduler $InputRoot `
    --tasks analyze `
    --model $Model `
    --reasoning-effort $ReasoningEffort `
    --max-workers $MaxWorkers `
    --skip-existing
if ($LASTEXITCODE -ne 0) {
    throw "Small-agent scheduler process failed. Rerun this script to resume missing modules only."
}

$suffixes = @(
    "cbc_other_lab_agent",
    "filmarray_gmtest_agent",
    "image_agent",
    "culture_agent",
    "molecular_microbiology_agent"
)
$missing = @()
$invalid = @()
$valid = 0
foreach ($patientDir in Get-ChildItem -LiteralPath $InputRoot -Directory -Filter "NGS_patient_*_json") {
    $prefix = $patientDir.Name -replace "_json$", ""
    foreach ($suffix in $suffixes) {
        $path = Join-Path $patientDir.FullName "agent_outputs\${prefix}_${suffix}.json"
        if (-not (Test-Path -LiteralPath $path)) {
            $missing += $path
            continue
        }
        try {
            Get-Content -Raw -Encoding utf8 -LiteralPath $path | ConvertFrom-Json | Out-Null
            $valid++
        }
        catch {
            $invalid += $path
        }
    }
}

$status = [ordered]@{
    schema_version = "small_agent_stage1_status.v2"
    model = $Model
    reasoning_effort = $ReasoningEffort
    source_root = $SourceRoot
    run_root = $RunRoot
    expected_patients = $patientCount
    expected_outputs = $expectedOutputs
    valid_outputs = $valid
    missing_count = $missing.Count
    invalid_json_count = $invalid.Count
    complete = ($valid -eq $expectedOutputs -and $missing.Count -eq 0 -and $invalid.Count -eq 0)
    missing_outputs = $missing
    invalid_json_outputs = $invalid
    completed_at = (Get-Date).ToString("o")
}
$statusPath = Join-Path $RunRoot "small_agent_stage1_status.json"
$status | ConvertTo-Json -Depth 8 | Set-Content -Encoding utf8 -LiteralPath $statusPath
if (-not $status.complete) {
    throw "Small-agent stage is incomplete. See $statusPath"
}
Write-Host "Small-agent stage complete: $valid/$expectedOutputs valid JSON outputs. Status: $statusPath"

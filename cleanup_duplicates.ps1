# PowerShell Script to Clean Up Duplicate Files
# Identifies and removes files with "1" suffix that are duplicates

param(
    [Parameter(Position=0)]
    [ValidateSet('dryrun', 'execute', 'force')]
    [string]$Mode = 'dryrun'
)

$ErrorActionPreference = "Stop"

Write-Host "`n=======================================================" -ForegroundColor Cyan
Write-Host "  DUPLICATE FILE CLEANUP SCRIPT" -ForegroundColor Cyan
Write-Host "=======================================================" -ForegroundColor Cyan
Write-Host "Mode: $($Mode.ToUpper())" -ForegroundColor Yellow
Write-Host "Start time: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
Write-Host "=======================================================`n" -ForegroundColor Cyan

# Get script directory (repository root)
$repoRoot = $PSScriptRoot

# Define file patterns to search for duplicates
$patterns = @('*1.mat', '*1.fig', '*1.dll', '*1.m', '*1.ecp', '*1.pdf',
              '*1.pptx', '*1.txt', '*1.html', '*1.png', '*1.jpg')

# Exclusion list (legitimate numbered files)
$exclusions = @(
    'fig1_1.m', 'fig3_1.m', 'fig4_1.m', 'fig5_1.m', 'fig6_1.m',
    'fig7_1.m', 'fig8_1.m', 'fig9_1.m', 'fig10_1.m', 'fig11_1.m', 'fig13_1.m',
    'lfex1.m', 'spence21.m', 'Exsnip1.m', 'sub2ind1.m', 'RMAOV31.m',
    'VnCArtifactRemoval11.m', 'VnCArtifactRemoval111.m',
    'GetMUAeLFPwithInterpolationSubroutine31.m',
    'GetMUAewithInterpolationSubroutine1.m',
    'VnCArtifactRemovalIndie1.m', 'VnCGetTrialData21.m',
    'VnCGetTrialDataCell1.m', 'VnCMicroStimTrialAlignTTL1.m',
    'VnCReadLogFile1.m', 'VnCgetNeuraldata31.m',
    'VnCRemoveArtifactsSubroutine11.m', 'VnCRemoveArtifactsSubroutine91.m',
    'GetFilepath1.m', 'openNSxFast_v1.m'
)

# Initialize arrays
$confirmed = @()
$skipped = @()

Write-Host "Scanning for duplicate files..." -ForegroundColor Green

# Find all potential duplicate files
foreach ($pattern in $patterns) {
    $files = Get-ChildItem -Path $repoRoot -Filter $pattern -Recurse -File -ErrorAction SilentlyContinue

    foreach ($file in $files) {
        # Skip if in exclusion list
        if ($exclusions -contains $file.Name) {
            $skipped += [PSCustomObject]@{
                Path = $file.FullName
                Reason = "Exclusion list (legitimate numbered file)"
            }
            continue
        }

        # Check if filename ends with "1" before extension
        $nameWithoutExt = [System.IO.Path]::GetFileNameWithoutExtension($file.Name)
        if (-not $nameWithoutExt.EndsWith('1')) {
            $skipped += [PSCustomObject]@{
                Path = $file.FullName
                Reason = "Does not match duplicate pattern"
            }
            continue
        }

        # Construct original filename
        $originalName = $nameWithoutExt.Substring(0, $nameWithoutExt.Length - 1)
        $originalFile = Join-Path $file.DirectoryName ($originalName + $file.Extension)

        # Check if original exists
        if (-not (Test-Path $originalFile)) {
            $skipped += [PSCustomObject]@{
                Path = $file.FullName
                Reason = "Original file not found"
            }
            continue
        }

        # Compare file sizes
        $originalInfo = Get-Item $originalFile
        if ($file.Length -eq $originalInfo.Length) {
            $confirmed += [PSCustomObject]@{
                Path = $file.FullName
                Original = $originalFile
                Size = $file.Length
                Reason = "Same size as original (likely identical)"
            }
        } else {
            $skipped += [PSCustomObject]@{
                Path = $file.FullName
                Reason = "Different size (dup: $($file.Length), orig: $($originalInfo.Length))"
            }
        }
    }
}

# Display results
Write-Host "`n=======================================================" -ForegroundColor Cyan
Write-Host "ANALYSIS RESULTS" -ForegroundColor Cyan
Write-Host "=======================================================" -ForegroundColor Cyan
Write-Host "Total files scanned:       $($confirmed.Count + $skipped.Count)"
Write-Host "Confirmed duplicates:      $($confirmed.Count)" -ForegroundColor Yellow
Write-Host "Skipped (not duplicates):  $($skipped.Count)"
Write-Host "=======================================================`n" -ForegroundColor Cyan

# Format bytes helper
function Format-Bytes {
    param([long]$Bytes)
    if ($Bytes -lt 1KB) { return "$Bytes B" }
    elseif ($Bytes -lt 1MB) { return "{0:N1} KB" -f ($Bytes / 1KB) }
    elseif ($Bytes -lt 1GB) { return "{0:N1} MB" -f ($Bytes / 1MB) }
    else { return "{0:N1} GB" -f ($Bytes / 1GB) }
}

# Show confirmed duplicates
if ($confirmed.Count -gt 0) {
    Write-Host "CONFIRMED DUPLICATES (will be deleted):" -ForegroundColor Yellow
    Write-Host "-------------------------------------------------------"

    $totalSize = 0
    for ($i = 0; $i -lt $confirmed.Count; $i++) {
        $relPath = $confirmed[$i].Path.Replace("$repoRoot\", "")
        Write-Host ("{0,3}. {1}" -f ($i+1), $relPath)
        Write-Host ("     Size: {0} | Reason: {1}" -f (Format-Bytes $confirmed[$i].Size), $confirmed[$i].Reason)
        $totalSize += $confirmed[$i].Size
    }

    Write-Host "-------------------------------------------------------"
    Write-Host "Total space to reclaim: $(Format-Bytes $totalSize)`n" -ForegroundColor Green
} else {
    Write-Host "No confirmed duplicates found.`n" -ForegroundColor Green
}

# Show skipped files (first 20)
if ($skipped.Count -gt 0) {
    Write-Host "SKIPPED FILES (first 20):" -ForegroundColor Gray
    Write-Host "-------------------------------------------------------"

    $displayCount = [Math]::Min(20, $skipped.Count)
    for ($i = 0; $i -lt $displayCount; $i++) {
        $relPath = $skipped[$i].Path.Replace("$repoRoot\", "")
        Write-Host ("{0,3}. {1}" -f ($i+1), $relPath)
        Write-Host ("     Reason: {0}" -f $skipped[$i].Reason)
    }

    if ($skipped.Count -gt 20) {
        Write-Host "... and $($skipped.Count - 20) more (see cleanup_report.txt for full list)"
    }
    Write-Host "-------------------------------------------------------`n"
}

# Generate report file
$reportFile = Join-Path $repoRoot "cleanup_report.txt"
$reportContent = @"
DUPLICATE FILE CLEANUP REPORT
Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')
Mode: $($Mode.ToUpper())

=======================================================
CONFIRMED DUPLICATES ($($confirmed.Count) files)
=======================================================
$($confirmed | ForEach-Object {
    "$($_.Path)`n  Original: $($_.Original)`n  Size: $($_.Size) bytes`n  Reason: $($_.Reason)`n"
} | Out-String)

=======================================================
SKIPPED FILES ($($skipped.Count) files)
=======================================================
$($skipped | ForEach-Object {
    "$($_.Path)`n  Reason: $($_.Reason)`n"
} | Out-String)
"@

$reportContent | Out-File -FilePath $reportFile -Encoding UTF8
Write-Host "Full report saved to: cleanup_report.txt`n" -ForegroundColor Green

# Exit if no duplicates found
if ($confirmed.Count -eq 0) {
    Write-Host "Nothing to delete." -ForegroundColor Green
    exit 0
}

# Handle dry run mode
if ($Mode -eq 'dryrun') {
    Write-Host "=======================================================" -ForegroundColor Cyan
    Write-Host "DRY RUN MODE - No files deleted" -ForegroundColor Yellow
    Write-Host "=======================================================" -ForegroundColor Cyan
    Write-Host "To actually delete these files, run:"
    Write-Host "  .\cleanup_duplicates.ps1 execute`n" -ForegroundColor Green
    Write-Host "To delete without confirmation prompts:"
    Write-Host "  .\cleanup_duplicates.ps1 force`n" -ForegroundColor Green
    exit 0
}

# Confirm deletion
if ($Mode -eq 'execute') {
    Write-Host "=======================================================" -ForegroundColor Red
    Write-Host "WARNING: About to delete $($confirmed.Count) files" -ForegroundColor Red
    Write-Host "=======================================================" -ForegroundColor Red
    $response = Read-Host "Continue? (yes/no)"
    if ($response -ne 'yes') {
        Write-Host "Deletion cancelled by user." -ForegroundColor Yellow
        exit 0
    }
}

# Create backup list
$timestamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$backupFile = Join-Path $repoRoot "deleted_files_$timestamp.txt"
"Files deleted on $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')`n" | Out-File -FilePath $backupFile -Encoding UTF8

# Delete files
Write-Host "`nDeleting files..." -ForegroundColor Yellow
$deleteCount = 0
$errorCount = 0

for ($i = 0; $i -lt $confirmed.Count; $i++) {
    try {
        $confirmed[$i].Path | Out-File -FilePath $backupFile -Append -Encoding UTF8
        Remove-Item -Path $confirmed[$i].Path -Force
        $deleteCount++
        $relPath = $confirmed[$i].Path.Replace("$repoRoot\", "")
        Write-Host ("  [{0,3}/{1,3}] Deleted: {2}" -f ($i+1), $confirmed.Count, $relPath) -ForegroundColor Green
    } catch {
        $errorCount++
        $relPath = $confirmed[$i].Path.Replace("$repoRoot\", "")
        Write-Host ("  [{0,3}/{1,3}] ERROR: {2}" -f ($i+1), $confirmed.Count, $relPath) -ForegroundColor Red
        Write-Host "    $($_.Exception.Message)" -ForegroundColor Red
    }
}

# Final summary
Write-Host "`n=======================================================" -ForegroundColor Cyan
Write-Host "CLEANUP COMPLETE" -ForegroundColor Cyan
Write-Host "=======================================================" -ForegroundColor Cyan
Write-Host "Successfully deleted: $deleteCount files" -ForegroundColor Green
Write-Host "Errors:               $errorCount files" -ForegroundColor $(if ($errorCount -gt 0) { 'Red' } else { 'Green' })
Write-Host "Backup list saved:    $backupFile"
Write-Host "Completion time:      $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
Write-Host "=======================================================`n" -ForegroundColor Cyan

Write-Host "Next steps:" -ForegroundColor Yellow
Write-Host "1. Review cleanup_report.txt to verify deletions"
Write-Host "2. Test your code to ensure nothing broke"
Write-Host "3. Commit changes:"
Write-Host "   git add -A"
Write-Host '   git commit -m "Remove duplicate files identified by cleanup script"'
Write-Host ""

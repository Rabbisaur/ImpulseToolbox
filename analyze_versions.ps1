# Analyze File Versions - Check if *1 files are newer or duplicates
# This script compares timestamps and content before recommending deletion

Write-Host "`n=======================================================" -ForegroundColor Cyan
Write-Host "  VERSION ANALYSIS TOOL" -ForegroundColor Cyan
Write-Host "=======================================================" -ForegroundColor Cyan
Write-Host "Checking timestamps and content...`n" -ForegroundColor Yellow

$repoRoot = $PSScriptRoot

# File patterns
$patterns = @('*1.mat', '*1.fig', '*1.dll', '*1.m')

# Exclusions
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

$results = @{
    Identical = @()
    NewerVersion = @()
    OlderVersion = @()
    Different = @()
}

# Scan files
foreach ($pattern in $patterns) {
    Get-ChildItem -Path $repoRoot -Filter $pattern -Recurse -File -ErrorAction SilentlyContinue | ForEach-Object {
        if ($exclusions -contains $_.Name) { return }

        $nameWithoutExt = [System.IO.Path]::GetFileNameWithoutExtension($_.Name)
        if (-not $nameWithoutExt.EndsWith('1')) { return }

        $originalName = $nameWithoutExt.Substring(0, $nameWithoutExt.Length - 1)
        $originalFile = Join-Path $_.DirectoryName ($originalName + $_.Extension)

        if (-not (Test-Path $originalFile)) { return }

        $origInfo = Get-Item $originalFile
        $versInfo = Get-Item $_.FullName

        # Compare content
        $hash1 = (Get-FileHash $originalFile -Algorithm MD5 -ErrorAction SilentlyContinue).Hash
        $hash2 = (Get-FileHash $_.FullName -Algorithm MD5 -ErrorAction SilentlyContinue).Hash
        $identical = ($hash1 -eq $hash2)

        $timeDiff = ($versInfo.LastWriteTime - $origInfo.LastWriteTime).TotalDays

        $item = [PSCustomObject]@{
            Original = $originalFile
            Versioned = $_.FullName
            OriginalDate = $origInfo.LastWriteTime
            VersionedDate = $versInfo.LastWriteTime
            DaysDiff = [Math]::Round($timeDiff, 1)
            OriginalSize = $origInfo.Length
            VersionedSize = $versInfo.Length
            Identical = $identical
        }

        if ($identical) {
            $results.Identical += $item
        } elseif ($timeDiff -gt 0.01) {
            $results.NewerVersion += $item
        } elseif ($timeDiff -lt -0.01) {
            $results.OlderVersion += $item
        } else {
            $results.Different += $item
        }
    }
}

# Display results
Write-Host "=======================================================" -ForegroundColor Cyan
Write-Host "ANALYSIS COMPLETE" -ForegroundColor Cyan
Write-Host "=======================================================" -ForegroundColor Cyan
Write-Host "Identical files:        $($results.Identical.Count)" -ForegroundColor Green
Write-Host "Newer *1 versions:      $($results.NewerVersion.Count)" -ForegroundColor Yellow
Write-Host "Older *1 versions:      $($results.OlderVersion.Count)" -ForegroundColor Gray
Write-Host "Different (same time):  $($results.Different.Count)" -ForegroundColor Magenta
Write-Host "=======================================================`n" -ForegroundColor Cyan

# Show identical files (safe to delete)
if ($results.Identical.Count -gt 0) {
    Write-Host "✓ IDENTICAL FILES - Safe to delete *1 versions:" -ForegroundColor Green
    Write-Host "-------------------------------------------------------"
    $results.Identical | Select-Object -First 15 | ForEach-Object {
        $relPath = $_.Versioned.Replace("$repoRoot\", "")
        Write-Host "  $relPath"
    }
    if ($results.Identical.Count -gt 15) {
        Write-Host "  ... and $($results.Identical.Count - 15) more`n"
    }
    Write-Host ""
}

# Show NEWER versions (CAUTION!)
if ($results.NewerVersion.Count -gt 0) {
    Write-Host "⚠️  NEWER VERSIONS - *1 file is NEWER than original!" -ForegroundColor Red
    Write-Host "These MAY be improved versions - DO NOT auto-delete!" -ForegroundColor Yellow
    Write-Host "-------------------------------------------------------"
    $results.NewerVersion | Select-Object -First 20 | ForEach-Object {
        $relPath = $_.Versioned.Replace("$repoRoot\", "")
        Write-Host "  $relPath" -ForegroundColor Yellow
        Write-Host "    Original: $($_.OriginalDate.ToString('yyyy-MM-dd HH:mm'))"
        Write-Host "    Versioned: $($_.VersionedDate.ToString('yyyy-MM-dd HH:mm')) [+$($_.DaysDiff) days]" -ForegroundColor Green
        if ($_.OriginalSize -ne $_.VersionedSize) {
            Write-Host "    Size changed: $($_.OriginalSize) -> $($_.VersionedSize) bytes"
        }
        Write-Host ""
    }
    if ($results.NewerVersion.Count -gt 20) {
        Write-Host "  ... and $($results.NewerVersion.Count - 20) more`n"
    }
}

# Show older versions
if ($results.OlderVersion.Count -gt 0) {
    Write-Host "OLDER VERSIONS - *1 file is OLDER (likely old backup):" -ForegroundColor Gray
    Write-Host "-------------------------------------------------------"
    $results.OlderVersion | Select-Object -First 10 | ForEach-Object {
        $relPath = $_.Versioned.Replace("$repoRoot\", "")
        Write-Host "  $relPath"
        Write-Host "    Older by: $([Math]::Abs($_.DaysDiff)) days"
    }
    if ($results.OlderVersion.Count -gt 10) {
        Write-Host "  ... and $($results.OlderVersion.Count - 10) more`n"
    }
    Write-Host ""
}

# Save detailed report
$reportFile = Join-Path $repoRoot "version_analysis_report.csv"
$allItems = @()
$allItems += $results.Identical | ForEach-Object { $_ | Add-Member -NotePropertyName Category -NotePropertyValue "Identical" -PassThru }
$allItems += $results.NewerVersion | ForEach-Object { $_ | Add-Member -NotePropertyName Category -NotePropertyValue "NewerVersion" -PassThru }
$allItems += $results.OlderVersion | ForEach-Object { $_ | Add-Member -NotePropertyName Category -NotePropertyValue "OlderVersion" -PassThru }
$allItems += $results.Different | ForEach-Object { $_ | Add-Member -NotePropertyName Category -NotePropertyValue "Different" -PassThru }

$allItems | Select-Object Category, Versioned, Original, DaysDiff, OriginalSize, VersionedSize, Identical |
    Export-Csv -Path $reportFile -NoTypeInformation -Encoding UTF8

Write-Host "Detailed report saved to: version_analysis_report.csv`n" -ForegroundColor Green

# Recommendations
Write-Host "=======================================================" -ForegroundColor Cyan
Write-Host "RECOMMENDATIONS" -ForegroundColor Cyan
Write-Host "=======================================================`n" -ForegroundColor Cyan

if ($results.NewerVersion.Count -gt 0) {
    Write-Host "⚠️  CRITICAL: $($results.NewerVersion.Count) files have NEWER *1 versions!" -ForegroundColor Red
    Write-Host ""
    Write-Host "ACTION REQUIRED for newer versions:" -ForegroundColor Yellow
    Write-Host "1. Review each file manually to determine which version to keep"
    Write-Host "2. For .m files, use MATLAB or a diff tool to compare"
    Write-Host "3. Options:"
    Write-Host "   - Keep *1 version: Rename file1.m -> file.m (delete old file.m)"
    Write-Host "   - Keep original: Delete file1.m"
    Write-Host "   - Keep both: Add comment explaining why`n"
}

Write-Host "SAFE ACTIONS:" -ForegroundColor Green
Write-Host "1. Delete $($results.Identical.Count) identical files (duplicates)"
Write-Host "2. Review $($results.OlderVersion.Count) older *1 versions (likely old backups)`n"

Write-Host "To review a specific file pair:"
Write-Host "  fc original.m file1.m        # Windows file compare"
Write-Host "  git diff original.m file1.m  # If both in git`n"

Write-Host "=======================================================" -ForegroundColor Cyan

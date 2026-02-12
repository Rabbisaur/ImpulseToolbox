# Smart Duplicate File Cleanup - Compares versions before deletion
# Checks timestamps and content to identify which version is newer

param(
    [Parameter(Position=0)]
    [ValidateSet('analyze', 'interactive', 'auto-keep-newer')]
    [string]$Mode = 'analyze'
)

$ErrorActionPreference = "Stop"

Write-Host "`n=======================================================" -ForegroundColor Cyan
Write-Host "  SMART DUPLICATE FILE CLEANUP" -ForegroundColor Cyan
Write-Host "=======================================================" -ForegroundColor Cyan
Write-Host "Mode: $($Mode.ToUpper())" -ForegroundColor Yellow
Write-Host "Start time: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
Write-Host "=======================================================`n" -ForegroundColor Cyan

# Get script directory
$repoRoot = $PSScriptRoot

# File patterns to check
$patterns = @('*1.mat', '*1.fig', '*1.dll', '*1.m', '*1.ecp', '*1.pdf',
              '*1.pptx', '*1.txt', '*1.html', '*1.png', '*1.jpg')

# Exclusion list
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

# Results categories
$trueDuplicates = @()      # Same content, same time
$newerVersions = @()       # Different content, file1 is newer
$olderVersions = @()       # Different content, file1 is older
$sameSizeDiffContent = @() # Same size but different content
$skipped = @()

Write-Host "Scanning and analyzing files..." -ForegroundColor Green

# Function to compare file content
function Compare-FileContent {
    param($file1, $file2)

    try {
        $hash1 = (Get-FileHash $file1 -Algorithm MD5).Hash
        $hash2 = (Get-FileHash $file2 -Algorithm MD5).Hash
        return $hash1 -eq $hash2
    } catch {
        return $false
    }
}

# Function to show file comparison
function Show-FileComparison {
    param($original, $versioned)

    $origInfo = Get-Item $original
    $versInfo = Get-Item $versioned

    $timeDiff = ($versInfo.LastWriteTime - $origInfo.LastWriteTime).TotalDays

    Write-Host "`n  Original: $($origInfo.Name)" -ForegroundColor White
    Write-Host "    Size:     $($origInfo.Length) bytes"
    Write-Host "    Modified: $($origInfo.LastWriteTime.ToString('yyyy-MM-dd HH:mm:ss'))"

    Write-Host "`n  Versioned: $($versInfo.Name)" -ForegroundColor Yellow
    Write-Host "    Size:     $($versInfo.Length) bytes"
    Write-Host "    Modified: $($versInfo.LastWriteTime.ToString('yyyy-MM-dd HH:mm:ss'))"

    if ($timeDiff -gt 0) {
        Write-Host "    Age:      NEWER by $([Math]::Abs($timeDiff).ToString('F1')) days" -ForegroundColor Green
    } elseif ($timeDiff -lt 0) {
        Write-Host "    Age:      OLDER by $([Math]::Abs($timeDiff).ToString('F1')) days" -ForegroundColor Red
    } else {
        Write-Host "    Age:      SAME timestamp" -ForegroundColor Gray
    }
}

# Scan all files
foreach ($pattern in $patterns) {
    $files = Get-ChildItem -Path $repoRoot -Filter $pattern -Recurse -File -ErrorAction SilentlyContinue

    foreach ($file in $files) {
        # Skip exclusions
        if ($exclusions -contains $file.Name) {
            $skipped += [PSCustomObject]@{
                Path = $file.FullName
                Reason = "Exclusion list"
            }
            continue
        }

        # Check naming pattern
        $nameWithoutExt = [System.IO.Path]::GetFileNameWithoutExtension($file.Name)
        if (-not $nameWithoutExt.EndsWith('1')) {
            $skipped += [PSCustomObject]@{
                Path = $file.FullName
                Reason = "Does not match pattern"
            }
            continue
        }

        # Find original file
        $originalName = $nameWithoutExt.Substring(0, $nameWithoutExt.Length - 1)
        $originalFile = Join-Path $file.DirectoryName ($originalName + $file.Extension)

        if (-not (Test-Path $originalFile)) {
            $skipped += [PSCustomObject]@{
                Path = $file.FullName
                Reason = "No original file found"
            }
            continue
        }

        $origInfo = Get-Item $originalFile
        $fileInfo = Get-Item $file.FullName

        # Compare content
        $sameContent = Compare-FileContent $originalFile $file.FullName
        $sameSize = ($origInfo.Length -eq $fileInfo.Length)
        $timeDiff = ($fileInfo.LastWriteTime - $origInfo.LastWriteTime).TotalSeconds

        if ($sameContent) {
            # True duplicate - identical content
            $trueDuplicates += [PSCustomObject]@{
                Versioned = $file.FullName
                Original = $originalFile
                Size = $fileInfo.Length
                TimeDiff = $timeDiff
                Status = "IDENTICAL"
            }
        } elseif ($timeDiff -gt 60) {
            # File1 is newer (more than 1 minute)
            $newerVersions += [PSCustomObject]@{
                Versioned = $file.FullName
                Original = $originalFile
                Size = $fileInfo.Length
                OrigSize = $origInfo.Length
                TimeDiff = $timeDiff
                Status = "NEWER_VERSION"
            }
        } elseif ($timeDiff -lt -60) {
            # File1 is older
            $olderVersions += [PSCustomObject]@{
                Versioned = $file.FullName
                Original = $originalFile
                Size = $fileInfo.Length
                OrigSize = $origInfo.Length
                TimeDiff = $timeDiff
                Status = "OLDER_VERSION"
            }
        } else {
            # Same time but different content
            $sameSizeDiffContent += [PSCustomObject]@{
                Versioned = $file.FullName
                Original = $originalFile
                Size = $fileInfo.Length
                TimeDiff = $timeDiff
                Status = "DIFFERENT"
            }
        }
    }
}

# Format bytes
function Format-Bytes {
    param([long]$Bytes)
    if ($Bytes -lt 1KB) { return "$Bytes B" }
    elseif ($Bytes -lt 1MB) { return "{0:N1} KB" -f ($Bytes / 1KB) }
    elseif ($Bytes -lt 1GB) { return "{0:N1} MB" -f ($Bytes / 1MB) }
    else { return "{0:N1} GB" -f ($Bytes / 1GB) }
}

# Display results
Write-Host "`n=======================================================" -ForegroundColor Cyan
Write-Host "ANALYSIS RESULTS" -ForegroundColor Cyan
Write-Host "=======================================================" -ForegroundColor Cyan
Write-Host "True duplicates (identical):     $($trueDuplicates.Count)" -ForegroundColor Green
Write-Host "Newer versions (file1 newer):    $($newerVersions.Count)" -ForegroundColor Yellow
Write-Host "Older versions (file1 older):    $($olderVersions.Count)" -ForegroundColor Red
Write-Host "Different content (same time):   $($sameSizeDiffContent.Count)" -ForegroundColor Magenta
Write-Host "Skipped files:                   $($skipped.Count)" -ForegroundColor Gray
Write-Host "=======================================================`n" -ForegroundColor Cyan

# Show true duplicates (safe to delete)
if ($trueDuplicates.Count -gt 0) {
    Write-Host "TRUE DUPLICATES - Safe to delete (identical content):" -ForegroundColor Green
    Write-Host "-------------------------------------------------------"
    for ($i = 0; $i -lt [Math]::Min(10, $trueDuplicates.Count); $i++) {
        $item = $trueDuplicates[$i]
        $relPath = $item.Versioned.Replace("$repoRoot\", "")
        Write-Host ("{0,3}. {1}" -f ($i+1), $relPath)
        Write-Host ("     Size: {0}" -f (Format-Bytes $item.Size))
    }
    if ($trueDuplicates.Count -gt 10) {
        Write-Host "... and $($trueDuplicates.Count - 10) more"
    }
    Write-Host ""
}

# Show newer versions (CAUTION - might want to keep these!)
if ($newerVersions.Count -gt 0) {
    Write-Host "NEWER VERSIONS - *1 file is NEWER than original!" -ForegroundColor Yellow
    Write-Host "⚠️  CAUTION: These may be improved versions!" -ForegroundColor Red
    Write-Host "-------------------------------------------------------"
    for ($i = 0; $i -lt [Math]::Min(10, $newerVersions.Count); $i++) {
        $item = $newerVersions[$i]
        $relPath = $item.Versioned.Replace("$repoRoot\", "")
        $days = [Math]::Abs($item.TimeDiff / 86400)
        Write-Host ("{0,3}. {1}" -f ($i+1), $relPath) -ForegroundColor Yellow
        Write-Host ("     Newer by: {0:F1} days | Size: {1} vs {2}" -f $days, (Format-Bytes $item.Size), (Format-Bytes $item.OrigSize))
    }
    if ($newerVersions.Count -gt 10) {
        Write-Host "... and $($newerVersions.Count - 10) more"
    }
    Write-Host ""
}

# Show older versions
if ($olderVersions.Count -gt 0) {
    Write-Host "OLDER VERSIONS - *1 file is OLDER than original (likely safe to delete):" -ForegroundColor Red
    Write-Host "-------------------------------------------------------"
    for ($i = 0; $i -lt [Math]::Min(10, $olderVersions.Count); $i++) {
        $item = $olderVersions[$i]
        $relPath = $item.Versioned.Replace("$repoRoot\", "")
        $days = [Math]::Abs($item.TimeDiff / 86400)
        Write-Host ("{0,3}. {1}" -f ($i+1), $relPath)
        Write-Host ("     Older by: {0:F1} days" -f $days)
    }
    if ($olderVersions.Count -gt 10) {
        Write-Host "... and $($olderVersions.Count - 10) more"
    }
    Write-Host ""
}

# Generate detailed report
$reportFile = Join-Path $repoRoot "cleanup_analysis_report.txt"
$report = @"
SMART DUPLICATE FILE CLEANUP - DETAILED ANALYSIS
Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')

=======================================================
SUMMARY
=======================================================
True duplicates (identical):     $($trueDuplicates.Count)
Newer versions (file1 newer):    $($newerVersions.Count)
Older versions (file1 older):    $($olderVersions.Count)
Different content (same time):   $($sameSizeDiffContent.Count)
Skipped files:                   $($skipped.Count)

=======================================================
TRUE DUPLICATES (Safe to delete)
=======================================================
$($trueDuplicates | ForEach-Object {
    "$($_.Versioned)`n  Original: $($_.Original)`n  Size: $($_.Size) bytes`n  Status: IDENTICAL CONTENT`n"
} | Out-String)

=======================================================
NEWER VERSIONS (⚠️ CAUTION - *1 file is NEWER!)
=======================================================
$($newerVersions | ForEach-Object {
    $days = [Math]::Abs($_.TimeDiff / 86400)
    "$($_.Versioned)`n  Original: $($_.Original)`n  Versioned size: $($_.Size) bytes`n  Original size: $($_.OrigSize) bytes`n  Newer by: $($days.ToString('F1')) days`n  ⚠️ REVIEW BEFORE DELETING`n"
} | Out-String)

=======================================================
OLDER VERSIONS (file1 is older - likely safe to delete)
=======================================================
$($olderVersions | ForEach-Object {
    $days = [Math]::Abs($_.TimeDiff / 86400)
    "$($_.Versioned)`n  Original: $($_.Original)`n  Older by: $($days.ToString('F1')) days`n"
} | Out-String)

=======================================================
DIFFERENT CONTENT (same timestamp)
=======================================================
$($sameSizeDiffContent | ForEach-Object {
    "$($_.Versioned)`n  Original: $($_.Original)`n  Status: Different content but same timestamp`n"
} | Out-String)
"@

$report | Out-File -FilePath $reportFile -Encoding UTF8
Write-Host "Detailed report saved to: cleanup_analysis_report.txt`n" -ForegroundColor Green

# Interactive mode
if ($Mode -eq 'interactive' -and ($trueDuplicates.Count -gt 0 -or $olderVersions.Count -gt 0 -or $newerVersions.Count -gt 0)) {
    Write-Host "`n=======================================================" -ForegroundColor Cyan
    Write-Host "INTERACTIVE CLEANUP MODE" -ForegroundColor Cyan
    Write-Host "=======================================================`n" -ForegroundColor Cyan

    $deletionList = @()

    # Auto-add true duplicates
    if ($trueDuplicates.Count -gt 0) {
        Write-Host "Adding $($trueDuplicates.Count) true duplicates to deletion list..." -ForegroundColor Green
        $deletionList += $trueDuplicates | ForEach-Object { $_.Versioned }
    }

    # Ask about older versions
    if ($olderVersions.Count -gt 0) {
        Write-Host "`nFound $($olderVersions.Count) files where *1 version is OLDER."
        $response = Read-Host "Delete older versions? (yes/no/review)"

        if ($response -eq 'yes') {
            $deletionList += $olderVersions | ForEach-Object { $_.Versioned }
        } elseif ($response -eq 'review') {
            foreach ($item in $olderVersions) {
                Write-Host "`n-------------------------------------------------------"
                Show-FileComparison $item.Original $item.Versioned
                $decision = Read-Host "`nDelete $($item.Versioned.Split('\')[-1])? (y/n)"
                if ($decision -eq 'y') {
                    $deletionList += $item.Versioned
                }
            }
        }
    }

    # Warn about newer versions
    if ($newerVersions.Count -gt 0) {
        Write-Host "`n⚠️  WARNING: Found $($newerVersions.Count) files where *1 version is NEWER!" -ForegroundColor Red
        Write-Host "These files may contain improved code. Review carefully!" -ForegroundColor Yellow
        Write-Host "`nOptions:"
        Write-Host "  1. Review each file individually (recommended)"
        Write-Host "  2. Keep all newer versions (skip deletion)"
        Write-Host "  3. Promote newer versions (rename file1.m to file.m)"
        $choice = Read-Host "Choose option (1/2/3)"

        if ($choice -eq '1') {
            foreach ($item in $newerVersions) {
                Write-Host "`n-------------------------------------------------------"
                Show-FileComparison $item.Original $item.Versioned

                if ($item.Versioned -match '\.m$') {
                    Write-Host "`nThis is a MATLAB file. Options:"
                    Write-Host "  k - Keep versioned file, delete original"
                    Write-Host "  d - Delete versioned file, keep original"
                    Write-Host "  p - Promote versioned (rename to original name)"
                    Write-Host "  s - Skip (keep both)"
                    Write-Host "  v - View diff (requires MATLAB or diff tool)"
                    $decision = Read-Host "Decision"
                } else {
                    $decision = Read-Host "Delete versioned file? (y/n/s=skip)"
                }

                switch ($decision) {
                    'k' { $deletionList += $item.Original }
                    'd' { $deletionList += $item.Versioned }
                    'y' { $deletionList += $item.Versioned }
                    'p' {
                        Write-Host "  Will rename $($item.Versioned) to $($item.Original)" -ForegroundColor Yellow
                        # Store for promotion
                    }
                }
            }
        } elseif ($choice -eq '3') {
            Write-Host "`nPromotion mode selected. This will:" -ForegroundColor Yellow
            Write-Host "  1. Delete the original file"
            Write-Host "  2. Rename file1.m to file.m"
            Write-Host "`nThis is a destructive operation!" -ForegroundColor Red
            $confirm = Read-Host "Are you sure? (yes/no)"
            if ($confirm -eq 'yes') {
                Write-Host "Promotion mode not yet implemented. Use manual renaming." -ForegroundColor Yellow
            }
        }
    }

    # Execute deletions
    if ($deletionList.Count -gt 0) {
        Write-Host "`n======================================================="
        Write-Host "Ready to delete $($deletionList.Count) files"
        Write-Host "=======================================================`n"
        $finalConfirm = Read-Host "Proceed with deletion? (yes/no)"

        if ($finalConfirm -eq 'yes') {
            $timestamp = Get-Date -Format 'yyyyMMdd_HHmmss'
            $backupFile = Join-Path $repoRoot "deleted_files_$timestamp.txt"

            foreach ($file in $deletionList) {
                try {
                    $file | Out-File -FilePath $backupFile -Append -Encoding UTF8
                    Remove-Item -Path $file -Force
                    Write-Host "✓ Deleted: $($file.Replace("$repoRoot\", ""))" -ForegroundColor Green
                } catch {
                    Write-Host "✗ Error: $($file.Replace("$repoRoot\", ""))" -ForegroundColor Red
                }
            }

            Write-Host "`n✓ Deletion complete. Backup list: $backupFile" -ForegroundColor Green
        } else {
            Write-Host "Deletion cancelled." -ForegroundColor Yellow
        }
    } else {
        Write-Host "No files selected for deletion." -ForegroundColor Yellow
    }
}

# Final recommendations
Write-Host "`n=======================================================" -ForegroundColor Cyan
Write-Host "RECOMMENDATIONS" -ForegroundColor Cyan
Write-Host "=======================================================`n" -ForegroundColor Cyan

if ($Mode -eq 'analyze') {
    Write-Host "1. Review the detailed report: cleanup_analysis_report.txt"
    Write-Host "2. For interactive cleanup, run:"
    Write-Host "   .\cleanup_duplicates_smart.ps1 interactive`n" -ForegroundColor Green

    if ($newerVersions.Count -gt 0) {
        Write-Host "⚠️  WARNING: $($newerVersions.Count) files have NEWER *1 versions!" -ForegroundColor Red
        Write-Host "   Review these carefully before deleting!`n"
    }
}

Write-Host "=======================================================" -ForegroundColor Cyan

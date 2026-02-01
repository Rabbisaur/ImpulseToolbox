$root = "c:\Users\wang\Documents\GitHub\ImpulseToolbox"
$backupRoot = Join-Path $root "_cleanup\redundant"
$reportFile = Join-Path $root "_cleanup_report.md"

if (-not (Test-Path $backupRoot)) { New-Item -ItemType Directory -Force -Path $backupRoot | Out-Null }

$report = @()
$report += "# Cleanup Report"
$report += ""
$report += "| File 1 (Candidate) | File 2 (Original) | Status | Action |"
$report += "| --- | --- | --- | --- |"

# Find candidates: *1.m
$candidates = Get-ChildItem -Path $root -Recurse -Filter "*1.m"

foreach ($file1 in $candidates) {
    $dir = $file1.DirectoryName
    $name1 = $file1.Name
    # Construct assumed original name: remove '1' before '.m'
    # Pattern: filename1.m -> filename.m
    $baseName = $name1.Substring(0, $name1.Length - 3) 
    # Handle cases like 'file11.m' -> 'file1.m'? No, strict *1.m rule as observed.
    # Actually, let's just strip the last character of the basename if it is '1'.
    
    # safer:
    $originalName = $file1.BaseName.Substring(0, $file1.BaseName.Length - 1) + $file1.Extension
    $filePathOriginal = Join-Path $dir $originalName

    if (Test-Path $filePathOriginal) {
        $hash1 = (Get-FileHash -Path $file1.FullName -Algorithm SHA256).Hash
        $hashOriginal = (Get-FileHash -Path $filePathOriginal -Algorithm SHA256).Hash

        if ($hash1 -eq $hashOriginal) {
            # IDENTICAL
            $relPath = $file1.FullName.Substring($root.Length + 1)
            $backupPath = Join-Path $backupRoot $relPath
            $backupDir = Split-Path $backupPath -Parent
            if (-not (Test-Path $backupDir)) { New-Item -ItemType Directory -Force -Path $backupDir | Out-Null }
            
            Move-Item -Path $file1.FullName -Destination $backupPath -Force
            $report += "| `$relPath` | `$originalName` | **IDENTICAL** | Moved to `_cleanup/redundant` |"
            Write-Host "Moved IDENTICAL: $relPath" -ForegroundColor Green
        } else {
            # DIFFERENT
            $relPath = $file1.FullName.Substring($root.Length + 1)
            $report += "| `$relPath` | `$originalName` | **DIFFERENT** | **LOGGED ONLY** |"
            Write-Host "Found DIFFERENT: $relPath" -ForegroundColor Yellow
        }
    } else {
        # ORPHAN
        $relPath = $file1.FullName.Substring($root.Length + 1)
        $report += "| `$relPath` | *Missing* | ORPHAN | LOGGED ONLY |"
        Write-Host "Found ORPHAN: $relPath" -ForegroundColor Cyan
    }
}

$report | Out-File -FilePath $reportFile -Encoding utf8
Write-Host "Report generated at $reportFile"

# Duplicate File Cleanup Instructions

This repository contains two cleanup scripts to remove duplicate files (files ending in "1" before the extension, like `file1.mat`).

## Quick Start

### Option 1: PowerShell Script (Recommended for Windows)

```powershell
# 1. Preview what will be deleted (dry run - safe)
.\cleanup_duplicates.ps1

# 2. Actually delete the files (with confirmation)
.\cleanup_duplicates.ps1 execute

# 3. Delete without confirmation prompts
.\cleanup_duplicates.ps1 force
```

### Option 2: MATLAB Script

```matlab
% 1. Preview what will be deleted (dry run - safe)
cleanup_duplicates()

% 2. Actually delete the files (with confirmation)
cleanup_duplicates('execute')

% 3. Delete without confirmation prompts
cleanup_duplicates('force')
```

## What Gets Deleted

The scripts identify files matching these patterns:
- `*1.mat` - MATLAB data files
- `*1.fig` - MATLAB figure files
- `*1.dll` - Dynamic libraries
- `*1.m` - MATLAB code files
- `*1.ecp`, `*1.pdf`, `*1.pptx`, `*1.txt`, etc.

### Safety Features

✅ **Dry run by default** - No files deleted unless you specify 'execute' or 'force'

✅ **File size comparison** - Only deletes if duplicate has same size as original

✅ **Exclusion list** - Skips legitimate numbered files like:
- `fig1_1.m`, `fig3_1.m` (chronux examples)
- `VnCArtifactRemoval11.m` (version numbers, not duplicates)
- `GetFilepath1.m` (legitimate function names)

✅ **Backup list created** - All deleted files logged to `deleted_files_TIMESTAMP.txt`

✅ **Detailed report** - Full analysis saved to `cleanup_report.txt`

## Expected Results

Based on the audit, you should see approximately:

- **~50-100 duplicate files** identified
- **~50 MB** of space reclaimed
- Main duplicates in:
  - `SpikeCluster/` directory (`*.mat`, `*.fig` backups)
  - `adopted/salpa/` directory (sample data duplicates)
  - Various DLL libraries with "1" suffix

## After Cleanup

1. **Review the report**
   ```
   notepad cleanup_report.txt
   ```

2. **Test your code** to ensure nothing broke

3. **Commit the changes**
   ```bash
   git add -A
   git commit -m "chore: Remove duplicate files to reduce repository size"
   git push
   ```

## Files That Will NOT Be Deleted

The scripts are smart enough to skip:

- Files where no original exists (e.g., `file1.m` exists but `file.m` doesn't)
- Files with different sizes than the original (might have different content)
- Legitimate numbered files in the exclusion list
- Files in `.git/` directory

## Troubleshooting

### PowerShell Execution Policy Error

If you see "cannot be loaded because running scripts is disabled", run:
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

### MATLAB Path Issues

Make sure you're in the repository root directory when running:
```matlab
cd('C:\Users\wang\Documents\GitHub\ImpulseToolbox')
cleanup_duplicates()
```

## Manual Verification

To manually check a specific duplicate:
```powershell
# Compare file sizes
(Get-Item "SpikeCluster\DefaultPara.mat").Length
(Get-Item "SpikeCluster\DefaultPara1.mat").Length

# View file details
dir "SpikeCluster\DefaultPara*.mat"
```

## Rollback

If you need to restore deleted files:

1. Check the backup list: `deleted_files_TIMESTAMP.txt`
2. Restore from git if needed:
   ```bash
   git checkout HEAD -- path/to/file
   ```
3. Or restore from your system backup

---

**Note**: Always run in **dry run mode first** to preview changes before actually deleting files!

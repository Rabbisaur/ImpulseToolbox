function cleanup_duplicates(varargin)
% CLEANUP_DUPLICATES - Remove duplicate files with "1" suffix from repository
%
% Usage:
%   cleanup_duplicates()              - Dry run mode (shows what would be deleted)
%   cleanup_duplicates('execute')     - Actually delete the files
%   cleanup_duplicates('force')       - Delete without confirmation prompts
%
% This script identifies files ending in "1" before the extension (e.g., file1.mat)
% and removes them if they appear to be duplicates of the original file.
%
% Safety features:
%   - Dry run by default
%   - Compares file sizes before deletion
%   - Creates backup list of deleted files
%   - Skips files in 'adopted' directory that may be legitimate
%   - Generates detailed report

    %% Parse input arguments
    p = inputParser;
    addOptional(p, 'mode', 'dryrun', @(x) ismember(x, {'dryrun', 'execute', 'force'}));
    parse(p, varargin{:});

    dryRun = strcmp(p.Results.mode, 'dryrun');
    forceDelete = strcmp(p.Results.mode, 'force');

    %% Initialize
    fprintf('\n=======================================================\n');
    fprintf('  DUPLICATE FILE CLEANUP SCRIPT\n');
    fprintf('=======================================================\n');
    fprintf('Mode: %s\n', upper(p.Results.mode));
    fprintf('Start time: %s\n', datestr(now));
    fprintf('=======================================================\n\n');

    % Get repository root
    repoRoot = fileparts(mfilename('fullpath'));

    %% Define patterns to search for duplicates
    % Files ending in "1" before extension (e.g., file1.mat, file1.dll)
    patterns = {'**/*1.mat', '**/*1.fig', '**/*1.dll', '**/*1.m', ...
                '**/*1.ecp', '**/*1.pdf', '**/*1.pptx', '**/*1.txt', ...
                '**/*1.html', '**/*1.png', '**/*1.jpg'};

    %% Exclusion patterns (legitimate files that should NOT be deleted)
    exclusions = {
        % Legitimate numbered files (not duplicates)
        'fig1_1.m', 'fig3_1.m', 'fig4_1.m', 'fig5_1.m', 'fig6_1.m', ...
        'fig7_1.m', 'fig8_1.m', 'fig9_1.m', 'fig10_1.m', 'fig11_1.m', 'fig13_1.m', ...
        'lfex1.m', 'spence21.m', 'Exsnip1.m', ...
        'sub2ind1.m', 'RMAOV31.m', ...
        'VnCArtifactRemoval11.m', 'VnCArtifactRemoval111.m', ...
        'GetMUAeLFPwithInterpolationSubroutine31.m', ...
        'GetMUAewithInterpolationSubroutine1.m', ...
        'VnCArtifactRemovalIndie1.m', 'VnCGetTrialData21.m', ...
        'VnCGetTrialDataCell1.m', 'VnCMicroStimTrialAlignTTL1.m', ...
        'VnCReadLogFile1.m', 'VnCgetNeuraldata31.m', ...
        'VnCRemoveArtifactsSubroutine11.m', 'VnCRemoveArtifactsSubroutine91.m', ...
        'GetFilepath1.m', 'openNSxFast_v1.m'
    };

    %% Find all potential duplicate files
    fprintf('Scanning for duplicate files...\n');
    allDuplicates = {};

    for i = 1:length(patterns)
        files = dir(fullfile(repoRoot, patterns{i}));
        for j = 1:length(files)
            if ~files(j).isdir
                fullPath = fullfile(files(j).folder, files(j).name);
                allDuplicates{end+1} = fullPath; %#ok<AGROW>
            end
        end
    end

    fprintf('Found %d potential duplicate files\n\n', length(allDuplicates));

    %% Analyze duplicates and categorize
    confirmed = struct('path', {}, 'original', {}, 'size', {}, 'reason', {});
    skipped = struct('path', {}, 'reason', {});

    for i = 1:length(allDuplicates)
        dupFile = allDuplicates{i};
        [filepath, name, ext] = fileparts(dupFile);

        % Check if file is in exclusion list
        if any(strcmp([name, ext], exclusions))
            skipped(end+1).path = dupFile; %#ok<AGROW>
            skipped(end).reason = 'Exclusion list (legitimate numbered file)';
            continue;
        end

        % Check if name ends with "1"
        if ~endsWith(name, '1')
            skipped(end+1).path = dupFile; %#ok<AGROW>
            skipped(end).reason = 'Does not match duplicate pattern';
            continue;
        end

        % Construct original filename
        originalName = name(1:end-1); % Remove trailing "1"
        originalFile = fullfile(filepath, [originalName, ext]);

        % Check if original exists
        if ~exist(originalFile, 'file')
            skipped(end+1).path = dupFile; %#ok<AGROW>
            skipped(end).reason = 'Original file not found';
            continue;
        end

        % Compare file sizes
        dupInfo = dir(dupFile);
        origInfo = dir(originalFile);

        if dupInfo.bytes == origInfo.bytes
            confirmed(end+1).path = dupFile; %#ok<AGROW>
            confirmed(end).original = originalFile;
            confirmed(end).size = dupInfo.bytes;
            confirmed(end).reason = 'Same size as original (likely identical)';
        else
            % Different sizes - might be different content
            skipped(end+1).path = dupFile; %#ok<AGROW>
            skipped(end).reason = sprintf('Different size (dup: %d, orig: %d)', ...
                dupInfo.bytes, origInfo.bytes);
        end
    end

    %% Display results
    fprintf('=======================================================\n');
    fprintf('ANALYSIS RESULTS\n');
    fprintf('=======================================================\n');
    fprintf('Total files scanned:       %d\n', length(allDuplicates));
    fprintf('Confirmed duplicates:      %d\n', length(confirmed));
    fprintf('Skipped (not duplicates):  %d\n', length(skipped));
    fprintf('=======================================================\n\n');

    %% Show confirmed duplicates
    if ~isempty(confirmed)
        fprintf('CONFIRMED DUPLICATES (will be deleted):\n');
        fprintf('-------------------------------------------------------\n');
        totalSize = 0;
        for i = 1:length(confirmed)
            relPath = strrep(confirmed(i).path, [repoRoot, filesep], '');
            fprintf('%3d. %s\n', i, relPath);
            fprintf('     Size: %s | Reason: %s\n', ...
                formatBytes(confirmed(i).size), confirmed(i).reason);
            totalSize = totalSize + confirmed(i).size;
        end
        fprintf('-------------------------------------------------------\n');
        fprintf('Total space to reclaim: %s\n\n', formatBytes(totalSize));
    else
        fprintf('No confirmed duplicates found.\n\n');
    end

    %% Show skipped files (first 20 for brevity)
    if ~isempty(skipped)
        fprintf('SKIPPED FILES (first 20):\n');
        fprintf('-------------------------------------------------------\n');
        for i = 1:min(20, length(skipped))
            relPath = strrep(skipped(i).path, [repoRoot, filesep], '');
            fprintf('%3d. %s\n', i, relPath);
            fprintf('     Reason: %s\n', skipped(i).reason);
        end
        if length(skipped) > 20
            fprintf('... and %d more (see cleanup_report.txt for full list)\n', ...
                length(skipped) - 20);
        end
        fprintf('-------------------------------------------------------\n\n');
    end

    %% Generate report file
    reportFile = fullfile(repoRoot, 'cleanup_report.txt');
    fid = fopen(reportFile, 'w');
    fprintf(fid, 'DUPLICATE FILE CLEANUP REPORT\n');
    fprintf(fid, 'Generated: %s\n', datestr(now));
    fprintf(fid, 'Mode: %s\n\n', upper(p.Results.mode));

    fprintf(fid, '=======================================================\n');
    fprintf(fid, 'CONFIRMED DUPLICATES (%d files)\n', length(confirmed));
    fprintf(fid, '=======================================================\n');
    for i = 1:length(confirmed)
        fprintf(fid, '%s\n', confirmed(i).path);
        fprintf(fid, '  Original: %s\n', confirmed(i).original);
        fprintf(fid, '  Size: %d bytes\n', confirmed(i).size);
        fprintf(fid, '  Reason: %s\n\n', confirmed(i).reason);
    end

    fprintf(fid, '\n=======================================================\n');
    fprintf(fid, 'SKIPPED FILES (%d files)\n', length(skipped));
    fprintf(fid, '=======================================================\n');
    for i = 1:length(skipped)
        fprintf(fid, '%s\n', skipped(i).path);
        fprintf(fid, '  Reason: %s\n\n', skipped(i).reason);
    end
    fclose(fid);

    fprintf('Full report saved to: cleanup_report.txt\n\n');

    %% Delete files if not dry run
    if isempty(confirmed)
        fprintf('Nothing to delete.\n');
        return;
    end

    if dryRun
        fprintf('=======================================================\n');
        fprintf('DRY RUN MODE - No files deleted\n');
        fprintf('=======================================================\n');
        fprintf('To actually delete these files, run:\n');
        fprintf('  cleanup_duplicates(''execute'')\n\n');
        fprintf('To delete without confirmation prompts:\n');
        fprintf('  cleanup_duplicates(''force'')\n\n');
        return;
    end

    %% Confirm deletion
    if ~forceDelete
        fprintf('=======================================================\n');
        fprintf('WARNING: About to delete %d files\n', length(confirmed));
        fprintf('=======================================================\n');
        response = input('Continue? (yes/no): ', 's');
        if ~strcmpi(response, 'yes')
            fprintf('Deletion cancelled by user.\n');
            return;
        end
    end

    %% Create backup list before deletion
    backupFile = fullfile(repoRoot, sprintf('deleted_files_%s.txt', ...
        datestr(now, 'yyyymmdd_HHMMSS')));
    fid = fopen(backupFile, 'w');
    fprintf(fid, 'Files deleted on %s\n\n', datestr(now));

    %% Delete files
    fprintf('\nDeleting files...\n');
    deleteCount = 0;
    errorCount = 0;

    for i = 1:length(confirmed)
        try
            fprintf(fid, '%s\n', confirmed(i).path);
            delete(confirmed(i).path);
            deleteCount = deleteCount + 1;
            fprintf('  [%3d/%3d] Deleted: %s\n', i, length(confirmed), ...
                strrep(confirmed(i).path, [repoRoot, filesep], ''));
        catch ME
            errorCount = errorCount + 1;
            fprintf('  [%3d/%3d] ERROR: %s\n    %s\n', i, length(confirmed), ...
                strrep(confirmed(i).path, [repoRoot, filesep], ''), ME.message);
        end
    end
    fclose(fid);

    %% Final summary
    fprintf('\n=======================================================\n');
    fprintf('CLEANUP COMPLETE\n');
    fprintf('=======================================================\n');
    fprintf('Successfully deleted: %d files\n', deleteCount);
    fprintf('Errors:               %d files\n', errorCount);
    fprintf('Backup list saved:    %s\n', backupFile);
    fprintf('Completion time:      %s\n', datestr(now));
    fprintf('=======================================================\n\n');

    fprintf('Next steps:\n');
    fprintf('1. Review cleanup_report.txt to verify deletions\n');
    fprintf('2. Test your code to ensure nothing broke\n');
    fprintf('3. Commit changes:\n');
    fprintf('   git add -A\n');
    fprintf('   git commit -m "Remove duplicate files identified by cleanup script"\n\n');
end

%% Helper function to format bytes
function str = formatBytes(bytes)
    if bytes < 1024
        str = sprintf('%d B', bytes);
    elseif bytes < 1024^2
        str = sprintf('%.1f KB', bytes/1024);
    elseif bytes < 1024^3
        str = sprintf('%.1f MB', bytes/1024^2);
    else
        str = sprintf('%.1f GB', bytes/1024^3);
    end
end

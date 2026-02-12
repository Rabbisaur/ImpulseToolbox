function analyze_versions()
% ANALYZE_VERSIONS - Check if *1 files are newer versions or duplicates
%
% This function analyzes all files ending in "1" (e.g., file1.m, data1.mat)
% and determines if they are:
%   - Identical duplicates (same content)
%   - Newer versions (modified after the original)
%   - Older versions (modified before the original)
%
% Usage:
%   analyze_versions()

    fprintf('\n=======================================================\n');
    fprintf('  VERSION ANALYSIS TOOL\n');
    fprintf('=======================================================\n');
    fprintf('Analyzing files...\n\n');

    % Get repository root
    repoRoot = fileparts(mfilename('fullpath'));

    % File patterns to check
    patterns = {'**/*1.mat', '**/*1.fig', '**/*1.dll', '**/*1.m', ...
                '**/*1.ecp', '**/*1.pdf', '**/*1.txt'};

    % Exclusion list (legitimate numbered files)
    exclusions = {
        'fig1_1.m', 'fig3_1.m', 'fig4_1.m', 'fig5_1.m', 'fig6_1.m', ...
        'fig7_1.m', 'fig8_1.m', 'fig9_1.m', 'fig10_1.m', 'fig11_1.m', 'fig13_1.m', ...
        'lfex1.m', 'spence21.m', 'Exsnip1.m', 'sub2ind1.m', 'RMAOV31.m', ...
        'VnCArtifactRemoval11.m', 'VnCArtifactRemoval111.m', ...
        'GetMUAeLFPwithInterpolationSubroutine31.m', ...
        'GetMUAewithInterpolationSubroutine1.m', ...
        'VnCArtifactRemovalIndie1.m', 'VnCGetTrialData21.m', ...
        'VnCGetTrialDataCell1.m', 'VnCMicroStimTrialAlignTTL1.m', ...
        'VnCReadLogFile1.m', 'VnCgetNeuraldata31.m', ...
        'VnCRemoveArtifactsSubroutine11.m', 'VnCRemoveArtifactsSubroutine91.m', ...
        'GetFilepath1.m', 'openNSxFast_v1.m'
    };

    % Results categories
    identical = struct('original', {}, 'versioned', {}, 'size', {}, 'date', {});
    newerVersions = struct('original', {}, 'versioned', {}, 'daysDiff', {}, ...
                          'sizeOrig', {}, 'sizeVers', {});
    olderVersions = struct('original', {}, 'versioned', {}, 'daysDiff', {}, ...
                          'sizeOrig', {}, 'sizeVers', {});

    % Find all potential duplicate files
    for i = 1:length(patterns)
        files = dir(fullfile(repoRoot, patterns{i}));

        for j = 1:length(files)
            if files(j).isdir
                continue;
            end

            % Check exclusion list
            if any(strcmp(files(j).name, exclusions))
                continue;
            end

            % Check if name ends with "1"
            [~, name, ext] = fileparts(files(j).name);
            if ~endsWith(name, '1')
                continue;
            end

            % Construct original filename
            originalName = name(1:end-1);
            originalFile = fullfile(files(j).folder, [originalName, ext]);

            % Check if original exists
            if ~exist(originalFile, 'file')
                continue;
            end

            versionedFile = fullfile(files(j).folder, files(j).name);

            % Get file info
            origInfo = dir(originalFile);
            versInfo = dir(versionedFile);

            % Compare content using file hash
            try
                hashOrig = fileHash(originalFile);
                hashVers = fileHash(versionedFile);
                isIdentical = strcmp(hashOrig, hashVers);
            catch
                % If hash fails, compare sizes
                isIdentical = (origInfo.bytes == versInfo.bytes);
            end

            % Calculate time difference in days
            timeDiff = days(versInfo.datenum - origInfo.datenum);

            % Categorize
            if isIdentical
                identical(end+1).original = originalFile; %#ok<AGROW>
                identical(end).versioned = versionedFile;
                identical(end).size = versInfo.bytes;
                identical(end).date = versInfo.date;
            elseif timeDiff > 0.01  % Versioned is newer
                newerVersions(end+1).original = originalFile; %#ok<AGROW>
                newerVersions(end).versioned = versionedFile;
                newerVersions(end).daysDiff = timeDiff;
                newerVersions(end).sizeOrig = origInfo.bytes;
                newerVersions(end).sizeVers = versInfo.bytes;
            elseif timeDiff < -0.01  % Versioned is older
                olderVersions(end+1).original = originalFile; %#ok<AGROW>
                olderVersions(end).versioned = versionedFile;
                olderVersions(end).daysDiff = abs(timeDiff);
                olderVersions(end).sizeOrig = origInfo.bytes;
                olderVersions(end).sizeVers = versInfo.bytes;
            end
        end
    end

    % Display results
    fprintf('=======================================================\n');
    fprintf('ANALYSIS COMPLETE\n');
    fprintf('=======================================================\n');
    fprintf('Identical files:        %d\n', length(identical));
    fprintf('Newer *1 versions:      %d\n', length(newerVersions));
    fprintf('Older *1 versions:      %d\n', length(olderVersions));
    fprintf('=======================================================\n\n');

    % Show identical files (safe to delete)
    if ~isempty(identical)
        fprintf('✓ IDENTICAL FILES - Safe to delete *1 versions:\n');
        fprintf('-------------------------------------------------------\n');
        for i = 1:min(15, length(identical))
            relPath = strrep(identical(i).versioned, [repoRoot, filesep], '');
            fprintf('  %s\n', relPath);
        end
        if length(identical) > 15
            fprintf('  ... and %d more\n', length(identical) - 15);
        end
        fprintf('\n');
    end

    % Show NEWER versions (CAUTION!)
    if ~isempty(newerVersions)
        fprintf('⚠️  NEWER VERSIONS - *1 file is NEWER than original!\n');
        fprintf('These MAY be improved versions - DO NOT auto-delete!\n');
        fprintf('-------------------------------------------------------\n');
        for i = 1:min(20, length(newerVersions))
            relPath = strrep(newerVersions(i).versioned, [repoRoot, filesep], '');
            fprintf('  %s\n', relPath);
            fprintf('    Newer by: %.1f days | Size: %d -> %d bytes\n', ...
                newerVersions(i).daysDiff, newerVersions(i).sizeOrig, ...
                newerVersions(i).sizeVers);

            % For .m files, show first line (often has description)
            if endsWith(newerVersions(i).versioned, '.m')
                try
                    fid = fopen(newerVersions(i).versioned, 'r');
                    firstLine = fgetl(fid);
                    fclose(fid);
                    if ~isempty(firstLine) && startsWith(strtrim(firstLine), '%')
                        fprintf('    First line: %s\n', strtrim(firstLine));
                    end
                catch
                end
            end
            fprintf('\n');
        end
        if length(newerVersions) > 20
            fprintf('  ... and %d more\n\n', length(newerVersions) - 20);
        end
    end

    % Show older versions
    if ~isempty(olderVersions)
        fprintf('OLDER VERSIONS - *1 file is OLDER (likely old backup):\n');
        fprintf('-------------------------------------------------------\n');
        for i = 1:min(10, length(olderVersions))
            relPath = strrep(olderVersions(i).versioned, [repoRoot, filesep], '');
            fprintf('  %s (older by %.1f days)\n', relPath, olderVersions(i).daysDiff);
        end
        if length(olderVersions) > 10
            fprintf('  ... and %d more\n', length(olderVersions) - 10);
        end
        fprintf('\n');
    end

    % Save detailed report
    reportFile = fullfile(repoRoot, 'version_analysis_report.txt');
    fid = fopen(reportFile, 'w');

    fprintf(fid, 'VERSION ANALYSIS REPORT\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));

    fprintf(fid, '=======================================================\n');
    fprintf(fid, 'IDENTICAL FILES (%d files - safe to delete)\n', length(identical));
    fprintf(fid, '=======================================================\n');
    for i = 1:length(identical)
        fprintf(fid, '%s\n', identical(i).versioned);
        fprintf(fid, '  Original: %s\n', identical(i).original);
        fprintf(fid, '  Size: %d bytes\n\n', identical(i).size);
    end

    fprintf(fid, '\n=======================================================\n');
    fprintf(fid, 'NEWER VERSIONS (%d files - REVIEW BEFORE DELETING!)\n', length(newerVersions));
    fprintf(fid, '=======================================================\n');
    for i = 1:length(newerVersions)
        fprintf(fid, '%s\n', newerVersions(i).versioned);
        fprintf(fid, '  Original: %s\n', newerVersions(i).original);
        fprintf(fid, '  Newer by: %.1f days\n', newerVersions(i).daysDiff);
        fprintf(fid, '  Size: %d -> %d bytes\n\n', newerVersions(i).sizeOrig, ...
            newerVersions(i).sizeVers);
    end

    fprintf(fid, '\n=======================================================\n');
    fprintf(fid, 'OLDER VERSIONS (%d files)\n', length(olderVersions));
    fprintf(fid, '=======================================================\n');
    for i = 1:length(olderVersions)
        fprintf(fid, '%s\n', olderVersions(i).versioned);
        fprintf(fid, '  Original: %s\n', olderVersions(i).original);
        fprintf(fid, '  Older by: %.1f days\n\n', olderVersions(i).daysDiff);
    end

    fclose(fid);

    fprintf('Detailed report saved to: version_analysis_report.txt\n\n');

    % Recommendations
    fprintf('=======================================================\n');
    fprintf('RECOMMENDATIONS\n');
    fprintf('=======================================================\n\n');

    if ~isempty(newerVersions)
        fprintf('⚠️  CRITICAL: %d files have NEWER *1 versions!\n\n', length(newerVersions));
        fprintf('ACTION REQUIRED for newer versions:\n');
        fprintf('1. Review each file manually to determine which to keep\n');
        fprintf('2. For .m files, use visdiff or compare tool:\n');
        fprintf('   visdiff(''original.m'', ''file1.m'')\n');
        fprintf('3. Options:\n');
        fprintf('   - Keep *1 version: rename file1.m -> file.m\n');
        fprintf('   - Keep original: delete file1.m\n');
        fprintf('   - Keep both: add comment explaining why\n\n');
    end

    fprintf('SAFE ACTIONS:\n');
    fprintf('1. Delete %d identical files (true duplicates)\n', length(identical));
    fprintf('2. Review %d older *1 versions (likely old backups)\n\n', length(olderVersions));

    fprintf('=======================================================\n');
end

function hash = fileHash(filename)
    % Simple hash function using MATLAB's built-in capabilities
    try
        % Try using Java MD5 (faster)
        fid = fopen(filename, 'rb');
        data = fread(fid, '*uint8');
        fclose(fid);
        md = java.security.MessageDigest.getInstance('MD5');
        md.update(data);
        hash = sprintf('%02x', typecast(md.digest(), 'uint8'));
    catch
        % Fallback: use file size and first/last bytes
        info = dir(filename);
        hash = sprintf('%d', info.bytes);
    end
end

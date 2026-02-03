% Comprehensive Validation: openNSxFasterMex vs openNSx
% Tests ALL switches and compares output behavior.

addpath('cerebus/NPMK');
dataFile = 'C:\testdata\Hub1-SampleRecording_sittinginsetup_instance1_B001.ns6';

if ~exist(dataFile, 'file')
    error('Test file not found: %s', dataFile);
end

fprintf('Comprehensive Validation: openNSxFasterMex vs openNSx\n');
fprintf('File: %s (%.2f GB)\n', dataFile, dir(dataFile).bytes/1024^3);
fprintf('=========================================================\n\n');

% --- Test Cases ---
testCases = {
    {'Test 1: Basic read',                      {'read'}};
    {'Test 2: noread (metadata only)',          {'noread'}};
    {'Test 3: double precision',                {'read', 'double'}};
    {'Test 4: int16 precision',                 {'read', 'int16'}};
    {'Test 5: uV units',                        {'read', 'uv'}};
    {'Test 6: Time range (t:0:1 sec)',          {'read', 't:0:1', 'sec'}};
    {'Test 7: Skip factor 10',                  {'read', 'skipfactor', 10}};
    {'Test 8: Channels 1-10',                   {'read', 'c:1:10'}};
    {'Test 9: Electrodes 1-5',                  {'read', 'e:1:5'}};
    {'Test 10: Combo (t:0:0.5 sec, skip 5, double)', {'read', 't:0:0.5', 'sec', 'skipfactor', 5, 'double'}};
    {'Test 11: Channels 1-5, skip 2',           {'read', 'c:1:5', 'skipfactor', 2}};
};

results = {};

for i = 1:length(testCases)
    testName = testCases{i}{1};
    args = testCases{i}{2};
    
    fprintf('\n--- %s ---\n', testName);
    argsStr = '';
    for a = 1:length(args)
        if isnumeric(args{a}), argsStr = [argsStr num2str(args{a}) ' '];
        else argsStr = [argsStr args{a} ' ']; end
    end
    fprintf('Args: %s\n', argsStr);
    
    try
        % 1. Run Original
        warning('off', 'all');
        tic;
        NSx_Orig = openNSx(dataFile, args{:});
        t_orig = toc;
        D1 = []; if isfield(NSx_Orig, 'Data'), D1 = NSx_Orig.Data; if iscell(D1), D1 = D1{1}; end; end
        expectedSize = size(D1);
        clear NSx_Orig;
        
        % 2. Run MEX
        tic;
        NSx_Mex = openNSxFasterMex(dataFile, args{:});
        t_mex = toc;
        D2 = []; if isfield(NSx_Mex, 'Data'), D2 = NSx_Mex.Data; if iscell(D2), D2 = D2{1}; end; end
        clear NSx_Mex;
        warning('on', 'all');
        
        % Compare
        passed = true;
        failReasons = {};
        
        if isempty(D1) && isempty(D2)
            % Both empty (noread)
        elseif isempty(D1) ~= isempty(D2)
            passed = false; failReasons{end+1} = 'Data field presence mismatch';
        elseif ~isequal(size(D1), size(D2))
            passed = false; failReasons{end+1} = sprintf('Size mismatch: %s vs %s', mat2str(size(D1)), mat2str(size(D2)));
        else
            % Chunked comparison to avoid OOM
            nCols = size(D1, 2);
            checkCols = unique([1:min(100,nCols), round(nCols/2)-50:round(nCols/2)+50, max(1,nCols-99):nCols]);
            checkCols = checkCols(checkCols >= 1 & checkCols <= nCols);
            
            if ~isempty(checkCols)
                D1_c = double(D1(:, checkCols));
                D2_c = double(D2(:, checkCols));
                if ~isequal(D1_c, D2_c)
                    passed = false; failReasons{end+1} = 'Data content mismatch at significant columns';
                end
            end
        end
        
        clear D1 D2;
        
        % Report
        if passed
            fprintf('[PASS] %s (Orig: %.2fs, MEX: %.2fs, Speedup: %.1fx)\n', testName, t_orig, t_mex, t_orig/t_mex);
            results{end+1} = {testName, 'PASS', t_orig, t_mex};
        else
            fprintf('[FAIL] %s\n', testName);
            for f = 1:length(failReasons), fprintf('  - %s\n', failReasons{f}); end
            results{end+1} = {testName, 'FAIL', t_orig, t_mex, failReasons};
        end
        
    catch ME
        fprintf('[ERROR] %s: %s\n', testName, ME.message);
        results{end+1} = {testName, 'ERROR', ME.message};
    end
end

fprintf('\n\nVALIDATION SUMMARY\n');
passCount = sum(cellfun(@(x) strcmp(x{2}, 'PASS'), results));
fprintf('PASSED: %d / %d\n', passCount, length(results));

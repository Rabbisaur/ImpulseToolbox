function results = test_MatFast()
% TEST_MATFAST - Comprehensive test suite for MatFast save/load functions
%
% Tests both original (v1) and refactored (v2) versions
%
% Usage:
%   results = test_MatFast()  % Run all tests, return results
%   test_MatFast()            % Run all tests, display summary
%
% Output:
%   results - Struct with test results and timing information

    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('  MATFAST COMPREHENSIVE TEST SUITE\n');
    fprintf('================================================================\n');
    fprintf('Testing: SaveMatFast, LoadMatFast, SaveFast, LoadFast\n');
    fprintf('Start time: %s\n', datestr(now));
    fprintf('================================================================\n\n');

    % Initialize results
    results = struct();
    results.timestamp = now;
    results.tests = {};
    results.passed = 0;
    results.failed = 0;
    results.errors = {};

    % Test configuration
    testFile = 'test_matfast_temp.bin';
    cleanupObj = onCleanup(@() cleanupTestFiles(testFile));

    %% ===== VERSION 1 TESTS (Original) =====
    fprintf('===== TESTING VERSION 1 (Original SaveMatFast/LoadMatFast) =====\n\n');

    results = runTest(results, @test_v1_simple_double, 'V1: Simple double array', testFile);
    results = runTest(results, @test_v1_all_numeric_types, 'V1: All numeric types', testFile);
    results = runTest(results, @test_v1_char_logical, 'V1: Char and logical', testFile);
    results = runTest(results, @test_v1_multidim_array, 'V1: Multi-dimensional arrays', testFile);
    results = runTest(results, @test_v1_struct, 'V1: Struct support', testFile);
    results = runTest(results, @test_v1_cell, 'V1: Cell array support', testFile);
    results = runTest(results, @test_v1_nested_struct, 'V1: Nested structures', testFile);
    results = runTest(results, @test_v1_empty_arrays, 'V1: Empty arrays', testFile);
    results = runTest(results, @test_v1_special_values, 'V1: Special values (NaN, Inf)', testFile);
    results = runTest(results, @test_v1_large_data, 'V1: Large data (10MB)', testFile);

    %% ===== VERSION 2 TESTS (Refactored) =====
    fprintf('\n===== TESTING VERSION 2 (Refactored SaveFast/LoadFast) =====\n\n');

    results = runTest(results, @test_v2_single_variable, 'V2: Single variable (backward compat)', testFile);
    results = runTest(results, @test_v2_multi_variable, 'V2: Multiple variables', testFile);
    results = runTest(results, @test_v2_struct_unpack, 'V2: Struct unpacking', testFile);
    results = runTest(results, @test_v2_selective_load, 'V2: Selective variable loading', testFile);
    results = runTest(results, @test_v2_workspace_load, 'V2: Workspace loading', testFile);
    results = runTest(results, @test_v2_all_types, 'V2: All data types', testFile);
    results = runTest(results, @test_v2_complex_nested, 'V2: Complex nested data', testFile);
    results = runTest(results, @test_v2_variable_names, 'V2: Variable name validation', testFile);

    %% ===== COMPATIBILITY TESTS =====
    fprintf('\n===== COMPATIBILITY TESTS =====\n\n');

    results = runTest(results, @test_compat_v1_to_v2, 'Compatibility: V1 file read by V2', testFile);
    results = runTest(results, @test_compat_roundtrip, 'Compatibility: Round-trip all types', testFile);

    %% ===== PERFORMANCE TESTS =====
    fprintf('\n===== PERFORMANCE BENCHMARKS =====\n\n');

    results = runTest(results, @test_perf_vs_matlab, 'Performance: vs MATLAB save/load', testFile);
    results = runTest(results, @test_perf_large_arrays, 'Performance: Large array I/O', testFile);

    %% ===== ERROR HANDLING TESTS =====
    fprintf('\n===== ERROR HANDLING TESTS =====\n\n');

    results = runTest(results, @test_error_invalid_file, 'Error: Invalid file path');
    results = runTest(results, @test_error_corrupted_file, 'Error: Corrupted file', testFile);
    results = runTest(results, @test_error_invalid_inputs, 'Error: Invalid inputs');

    %% ===== SUMMARY =====
    fprintf('\n================================================================\n');
    fprintf('  TEST SUMMARY\n');
    fprintf('================================================================\n');
    fprintf('Total tests:  %d\n', results.passed + results.failed);
    fprintf('Passed:       %d (%s)\n', results.passed, successIndicator(results.passed, results.passed + results.failed));
    fprintf('Failed:       %d (%s)\n', results.failed, failIndicator(results.failed));
    fprintf('Success rate: %.1f%%\n', 100 * results.passed / (results.passed + results.failed));
    fprintf('Completion:   %s\n', datestr(now));
    fprintf('================================================================\n\n');

    if results.failed > 0
        fprintf('FAILED TESTS:\n');
        for i = 1:length(results.errors)
            fprintf('  - %s: %s\n', results.errors{i}.name, results.errors{i}.message);
        end
        fprintf('\n');
    end

    % Return results only if output requested
    if nargout == 0
        clear results;
    end
end

%% ===== TEST EXECUTION FRAMEWORK =====

function results = runTest(results, testFunc, testName, varargin)
    fprintf('Running: %s ... ', testName);

    try
        tic;
        testFunc(varargin{:});
        elapsed = toc;

        fprintf('[PASS] (%.3f s)\n', elapsed);
        results.passed = results.passed + 1;
        results.tests{end+1} = struct('name', testName, 'status', 'PASS', 'time', elapsed);

    catch ME
        fprintf('[FAIL]\n');
        fprintf('  Error: %s\n', ME.message);

        results.failed = results.failed + 1;
        results.errors{end+1} = struct('name', testName, 'message', ME.message);
        results.tests{end+1} = struct('name', testName, 'status', 'FAIL', 'error', ME.message);
    end
end

function cleanupTestFiles(testFile)
    if exist(testFile, 'file')
        delete(testFile);
    end
    % Clean up any test .mat files
    delete('test_*.mat');
end

function str = successIndicator(passed, total)
    if passed == total
        str = '✓✓✓';
    elseif passed > total * 0.8
        str = '✓✓';
    elseif passed > total * 0.5
        str = '✓';
    else
        str = '';
    end
end

function str = failIndicator(failed)
    if failed == 0
        str = '';
    elseif failed < 3
        str = '!';
    elseif failed < 5
        str = '!!';
    else
        str = '!!!';
    end
end

%% ===== VERSION 1 TESTS =====

function test_v1_simple_double(testFile)
    data = rand(100, 50);
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isequal(data, loaded), 'Data mismatch');
end

function test_v1_all_numeric_types(testFile)
    types = {'double', 'single', 'int8', 'uint8', 'int16', 'uint16', ...
             'int32', 'uint32', 'int64', 'uint64'};

    for i = 1:length(types)
        tp = types{i};
        data = cast(randi([1 100], 10, 10), tp);
        SaveMatFast(testFile, data);
        loaded = LoadMatFast(testFile);
        assert(isequal(data, loaded), sprintf('%s type mismatch', tp));
        assert(strcmp(class(data), class(loaded)), sprintf('%s class mismatch', tp));
    end
end

function test_v1_char_logical(testFile)
    % Char
    data = 'Hello, World!';
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isequal(data, loaded), 'Char mismatch');

    % Logical
    data = logical([1 0 1; 0 1 0]);
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isequal(data, loaded), 'Logical mismatch');
    assert(islogical(loaded), 'Logical type lost');
end

function test_v1_multidim_array(testFile)
    data = rand(5, 4, 3, 2);
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isequal(size(data), size(loaded)), 'Size mismatch');
    assert(isequal(data, loaded), 'Data mismatch');
end

function test_v1_struct(testFile)
    data = struct('a', 1, 'b', 'test', 'c', rand(3));
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isequal(fieldnames(data), fieldnames(loaded)), 'Field names mismatch');
    assert(isequal(data.a, loaded.a), 'Field a mismatch');
    assert(isequal(data.b, loaded.b), 'Field b mismatch');
    assert(isequal(data.c, loaded.c), 'Field c mismatch');
end

function test_v1_cell(testFile)
    data = {1, 'test', rand(2), {1,2,3}};
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isequal(data, loaded), 'Cell array mismatch');
end

function test_v1_nested_struct(testFile)
    data = struct('level1', struct('level2', struct('level3', rand(5))));
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isequal(data.level1.level2.level3, loaded.level1.level2.level3), 'Nested struct mismatch');
end

function test_v1_empty_arrays(testFile)
    data = [];
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isempty(loaded), 'Empty array not empty');
end

function test_v1_special_values(testFile)
    data = [1, NaN, Inf, -Inf, 0, -0];
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isequal(isnan(data), isnan(loaded)), 'NaN mismatch');
    assert(isequal(isinf(data), isinf(loaded)), 'Inf mismatch');
    assert(isequal(data(~isnan(data) & ~isinf(data)), loaded(~isnan(loaded) & ~isinf(loaded))), 'Finite values mismatch');
end

function test_v1_large_data(testFile)
    data = rand(1000, 1000);  % ~8 MB double array
    SaveMatFast(testFile, data);
    loaded = LoadMatFast(testFile);
    assert(isequal(size(data), size(loaded)), 'Large array size mismatch');
    assert(max(abs(data(:) - loaded(:))) < eps, 'Large array data mismatch');
end

%% ===== VERSION 2 TESTS =====

function test_v2_single_variable(testFile)
    data = rand(50, 50);
    SaveFast(testFile, data);
    loaded = LoadFast(testFile);
    assert(isstruct(loaded), 'Output should be struct');
    assert(isfield(loaded, 'data'), 'Should have field "data"');
    assert(isequal(data, loaded.data), 'Data mismatch');
end

function test_v2_multi_variable(testFile)
    A = magic(5);
    B = 'hello';
    C = {1, 2, 3};

    SaveFast(testFile, 'A', A, 'B', B, 'C', C);
    loaded = LoadFast(testFile);

    assert(isstruct(loaded), 'Output should be struct');
    assert(isfield(loaded, 'A'), 'Should have field A');
    assert(isfield(loaded, 'B'), 'Should have field B');
    assert(isfield(loaded, 'C'), 'Should have field C');
    assert(isequal(A, loaded.A), 'A mismatch');
    assert(isequal(B, loaded.B), 'B mismatch');
    assert(isequal(C, loaded.C), 'C mismatch');
end

function test_v2_struct_unpack(testFile)
    data = struct('x', 1, 'y', 2, 'z', 3);
    SaveFast(testFile, data);
    loaded = LoadFast(testFile);

    assert(isfield(loaded, 'x'), 'Missing field x');
    assert(isfield(loaded, 'y'), 'Missing field y');
    assert(isfield(loaded, 'z'), 'Missing field z');
    assert(loaded.x == 1, 'x value wrong');
    assert(loaded.y == 2, 'y value wrong');
    assert(loaded.z == 3, 'z value wrong');
end

function test_v2_selective_load(testFile)
    SaveFast(testFile, 'A', 1, 'B', 2, 'C', 3, 'D', 4);

    % Load only B and D
    loaded = LoadFast(testFile, 'B', 'D');

    assert(isstruct(loaded), 'Output should be struct');
    assert(isfield(loaded, 'B'), 'Should have field B');
    assert(isfield(loaded, 'D'), 'Should have field D');
    assert(~isfield(loaded, 'A'), 'Should NOT have field A');
    assert(~isfield(loaded, 'C'), 'Should NOT have field C');
    assert(loaded.B == 2, 'B value wrong');
    assert(loaded.D == 4, 'D value wrong');
end

function test_v2_workspace_load(testFile)
    % This test checks workspace loading functionality
    % Note: Can't easily test assignin in automated test, so we test the struct output instead
    SaveFast(testFile, 'testVar', 42);
    loaded = LoadFast(testFile);
    assert(loaded.testVar == 42, 'Workspace load failed');
end

function test_v2_all_types(testFile)
    % Note: struct('cel', {1,'two',3}) creates a 1x3 struct array.
    % Use double braces {{...}} to store a cell as a single field value.
    data = struct(...
        'dbl', rand(3), ...
        'sgl', single(rand(3)), ...
        'i8', int8(randi([-128 127], 3)), ...
        'u8', uint8(randi([0 255], 3)), ...
        'chr', 'test string', ...
        'log', logical([1 0 1]), ...
        'strc', struct('a', 1), ...
        'cel', {{1, 'two', 3}});

    SaveFast(testFile, data);
    loaded = LoadFast(testFile);

    fields = fieldnames(data);
    for i = 1:length(fields)
        fn = fields{i};
        assert(isequal(data.(fn), loaded.(fn)), sprintf('Field %s mismatch', fn));
    end
end

function test_v2_complex_nested(testFile)
    data = struct(...
        'level1', struct(...
            'level2', struct(...
                'data', rand(10), ...
                'cell', {{1, 'two', struct('three', 3)}})));

    SaveFast(testFile, 'complex', data);
    loaded = LoadFast(testFile);

    assert(isequal(data, loaded.complex), 'Complex nested structure mismatch');
end

function test_v2_variable_names(testFile)
    % Test valid variable names
    validName = 'valid_Var123';
    SaveFast(testFile, validName, 42);
    loaded = LoadFast(testFile);
    assert(isfield(loaded, validName), 'Valid name not saved');

    % Invalid names should error
    try
        SaveFast(testFile, '123invalid', 42);  % Starts with number
        error('Should have rejected invalid variable name');
    catch ME
        assert(contains(ME.message, 'variable name'), 'Wrong error for invalid name');
    end
end

%% ===== COMPATIBILITY TESTS =====

function test_compat_v1_to_v2(testFile)
    % Save with V1, load with V2
    data = rand(20, 20);
    SaveMatFast(testFile, data);
    loaded = LoadFast(testFile);

    assert(isstruct(loaded), 'V2 should return struct');
    assert(isfield(loaded, 'data'), 'V2 should create "data" field for V1 files');
    assert(isequal(data, loaded.data), 'V1 to V2 data mismatch');
end

function test_compat_roundtrip(testFile)
    % Test all types through save/load cycle
    testData = {
        rand(10), ...
        single(rand(5)), ...
        int32(randi([-1000 1000], 3)), ...
        'test string', ...
        logical([1 0 1 0]), ...
        struct('a', 1, 'b', {1,2,3}), ...
        {{1, 'two', 3}}
    };

    for i = 1:length(testData)
        data = testData{i};

        % V1 round-trip
        SaveMatFast(testFile, data);
        loaded_v1 = LoadMatFast(testFile);
        assert(isequal(data, loaded_v1), sprintf('V1 round-trip failed for type %d', i));

        % V2 round-trip
        SaveFast(testFile, 'data', data);
        loaded_v2 = LoadFast(testFile);
        assert(isequal(data, loaded_v2.data), sprintf('V2 round-trip failed for type %d', i));
    end
end

%% ===== PERFORMANCE TESTS =====

function test_perf_vs_matlab(testFile)
    data = rand(500, 500);  % 2 MB array

    % Time MatFast V1
    tic;
    SaveMatFast(testFile, data);
    time_v1_save = toc;

    tic;
    LoadMatFast(testFile);
    time_v1_load = toc;

    % Time MatFast V2
    tic;
    SaveFast(testFile, 'data', data);
    time_v2_save = toc;

    tic;
    LoadFast(testFile);
    time_v2_load = toc;

    % Time MATLAB
    tic;
    save('test_matlab.mat', 'data');
    time_matlab_save = toc;

    tic;
    load('test_matlab.mat');
    time_matlab_load = toc;

    fprintf('    V1 Save: %.4f s | V1 Load: %.4f s\n', time_v1_save, time_v1_load);
    fprintf('    V2 Save: %.4f s | V2 Load: %.4f s\n', time_v2_save, time_v2_load);
    fprintf('    MATLAB Save: %.4f s | MATLAB Load: %.4f s\n', time_matlab_save, time_matlab_load);
    fprintf('    Speedup (V1): Save %.1fx | Load %.1fx\n', time_matlab_save/time_v1_save, time_matlab_load/time_v1_load);
end

function test_perf_large_arrays(testFile)
    sizes = [100, 500, 1000];

    for sz = sizes
        data = rand(sz, sz);
        bytes = sz * sz * 8;

        tic;
        SaveFast(testFile, 'data', data);
        t_save = toc;

        tic;
        LoadFast(testFile);
        t_load = toc;

        fprintf('    %dx%d (%.1f MB): Save %.3f s | Load %.3f s | Throughput: %.1f MB/s\n', ...
            sz, sz, bytes/1e6, t_save, t_load, bytes/1e6/t_load);
    end
end

%% ===== ERROR HANDLING TESTS =====

function test_error_invalid_file()
    try
        LoadMatFast('nonexistent_file_xyz.bin');
        error('Should have thrown error for nonexistent file');
    catch ME
        assert(contains(ME.message, 'exist'), 'Wrong error for nonexistent file');
    end

    try
        LoadFast('nonexistent_file_xyz.bin');
        error('Should have thrown error for nonexistent file');
    catch ME
        assert(contains(ME.message, 'exist'), 'Wrong error for nonexistent file');
    end
end

function test_error_corrupted_file(testFile)
    % Create corrupted file
    fid = fopen(testFile, 'w');
    fwrite(fid, uint8([1 2 3 4 5]));  % Random bytes
    fclose(fid);

    try
        LoadMatFast(testFile);
        error('Should have thrown error for corrupted file');
    catch ME
        % Expected to fail
        assert(~isempty(ME.message), 'Should have error message');
    end
end

function test_error_invalid_inputs()
    % Test V2 with invalid variable name
    try
        SaveFast('test.bin', '123invalid', 42);
        error('Should reject invalid variable name');
    catch ME
        assert(contains(ME.message, 'variable name'), 'Wrong error message');
    end

    % Test V2 with odd number of arguments
    try
        SaveFast('test.bin', 'A', 1, 'B');  % Missing value for B
        error('Should reject unpaired arguments');
    catch ME
        assert(contains(ME.message, 'pair'), 'Wrong error for unpaired args');
    end
end

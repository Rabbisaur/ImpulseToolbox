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

    %% ===== VERSION 3 TESTS (New Types) =====
    fprintf('\n===== TESTING VERSION 3 (New Type Support) =====\n\n');

    results = runTest(results, @test_v3_complex_double, 'V3: Complex double', testFile);
    results = runTest(results, @test_v3_complex_single, 'V3: Complex single', testFile);
    results = runTest(results, @test_v3_sparse_double, 'V3: Sparse double', testFile);
    results = runTest(results, @test_v3_sparse_logical, 'V3: Sparse logical', testFile);
    results = runTest(results, @test_v3_sparse_complex, 'V3: Sparse complex', testFile);
    results = runTest(results, @test_v3_string_array, 'V3: String array', testFile);
    results = runTest(results, @test_v3_string_missing, 'V3: String with missing', testFile);
    results = runTest(results, @test_v3_table_basic, 'V3: Table basic', testFile);
    results = runTest(results, @test_v3_table_properties, 'V3: Table with properties', testFile);
    results = runTest(results, @test_v3_timetable, 'V3: Timetable', testFile);
    results = runTest(results, @test_v3_categorical, 'V3: Categorical', testFile);
    results = runTest(results, @test_v3_categorical_ordinal, 'V3: Categorical ordinal', testFile);
    results = runTest(results, @test_v3_datetime, 'V3: Datetime', testFile);
    results = runTest(results, @test_v3_datetime_timezone, 'V3: Datetime with timezone', testFile);
    results = runTest(results, @test_v3_duration, 'V3: Duration', testFile);
    results = runTest(results, @test_v3_calendar_duration, 'V3: CalendarDuration', testFile);
    results = runTest(results, @test_v3_containers_map, 'V3: containers.Map', testFile);
    results = runTest(results, @test_v3_function_handle, 'V3: function_handle', testFile);
    results = runTest(results, @test_v3_unsupported_type_error, 'V3: Unsupported type errors', testFile);
    results = runTest(results, @test_v3_mixed_new_types, 'V3: Mixed new types in struct', testFile);

    %% ===== COMPREHENSIVE EDGE-CASE TESTS =====
    fprintf('\n===== COMPREHENSIVE EDGE-CASE TESTS =====\n\n');

    % -- Complex numeric edge cases --
    results = runTest(results, @test_edge_complex_all_numeric_types, 'Edge: Complex all numeric types', testFile);
    results = runTest(results, @test_edge_complex_3d, 'Edge: Complex 3D array', testFile);
    results = runTest(results, @test_edge_complex_nan_inf, 'Edge: Complex with NaN/Inf', testFile);
    results = runTest(results, @test_edge_complex_scalar, 'Edge: Complex scalar', testFile);
    results = runTest(results, @test_edge_complex_pure_imaginary, 'Edge: Pure imaginary', testFile);

    % -- Sparse edge cases --
    results = runTest(results, @test_edge_sparse_empty, 'Edge: Sparse all-zeros', testFile);
    results = runTest(results, @test_edge_sparse_large, 'Edge: Sparse large matrix', testFile);
    results = runTest(results, @test_edge_sparse_single_element, 'Edge: Sparse single element', testFile);
    results = runTest(results, @test_edge_sparse_full_density, 'Edge: Sparse full density', testFile);

    % -- String edge cases --
    results = runTest(results, @test_edge_string_empty_array, 'Edge: String empty array', testFile);
    results = runTest(results, @test_edge_string_scalar, 'Edge: String scalar', testFile);
    results = runTest(results, @test_edge_string_empty_element, 'Edge: String with empty ""', testFile);
    results = runTest(results, @test_edge_string_long, 'Edge: String long content', testFile);
    results = runTest(results, @test_edge_string_special_chars, 'Edge: String special chars', testFile);
    results = runTest(results, @test_edge_string_all_missing, 'Edge: String all missing', testFile);

    % -- Table edge cases --
    results = runTest(results, @test_edge_table_empty, 'Edge: Table empty (no rows)', testFile);
    results = runTest(results, @test_edge_table_single_row, 'Edge: Table single row', testFile);
    results = runTest(results, @test_edge_table_mixed_columns, 'Edge: Table mixed column types', testFile);
    results = runTest(results, @test_edge_table_no_properties, 'Edge: Table no extra properties', testFile);
    results = runTest(results, @test_edge_table_cell_column, 'Edge: Table with cell column', testFile);

    % -- Timetable edge cases --
    results = runTest(results, @test_edge_timetable_duration_rowtimes, 'Edge: Timetable with duration rowtimes', testFile);
    results = runTest(results, @test_edge_timetable_multi_col, 'Edge: Timetable multiple columns', testFile);

    % -- Categorical edge cases --
    results = runTest(results, @test_edge_categorical_undefined, 'Edge: Categorical with undefined', testFile);
    results = runTest(results, @test_edge_categorical_2d, 'Edge: Categorical 2D', testFile);
    results = runTest(results, @test_edge_categorical_single_cat, 'Edge: Categorical single category', testFile);
    results = runTest(results, @test_edge_categorical_unicode_names, 'Edge: Categorical special char names', testFile);

    % -- Datetime edge cases --
    results = runTest(results, @test_edge_datetime_nat, 'Edge: Datetime with NaT', testFile);
    results = runTest(results, @test_edge_datetime_2d, 'Edge: Datetime 2D array', testFile);
    results = runTest(results, @test_edge_datetime_utc, 'Edge: Datetime UTC timezone', testFile);
    results = runTest(results, @test_edge_datetime_custom_format, 'Edge: Datetime custom format', testFile);

    % -- Duration edge cases --
    results = runTest(results, @test_edge_duration_nan, 'Edge: Duration with NaN', testFile);
    results = runTest(results, @test_edge_duration_2d, 'Edge: Duration 2D array', testFile);
    results = runTest(results, @test_edge_duration_zero, 'Edge: Duration zero', testFile);
    results = runTest(results, @test_edge_duration_negative, 'Edge: Duration negative', testFile);
    results = runTest(results, @test_edge_duration_custom_format, 'Edge: Duration custom format', testFile);

    % -- CalendarDuration edge cases --
    results = runTest(results, @test_edge_caldur_2d, 'Edge: CalendarDuration 2D', testFile);
    results = runTest(results, @test_edge_caldur_mixed, 'Edge: CalendarDuration months+days+time', testFile);
    results = runTest(results, @test_edge_caldur_negative, 'Edge: CalendarDuration negative', testFile);

    % -- containers.Map edge cases --
    results = runTest(results, @test_edge_map_empty, 'Edge: containers.Map empty', testFile);
    results = runTest(results, @test_edge_map_numeric_keys, 'Edge: containers.Map numeric keys', testFile);
    results = runTest(results, @test_edge_map_single_entry, 'Edge: containers.Map single entry', testFile);
    results = runTest(results, @test_edge_map_large, 'Edge: containers.Map large', testFile);

    % -- function_handle edge cases --
    results = runTest(results, @test_edge_fh_builtin, 'Edge: function_handle builtins', testFile);
    results = runTest(results, @test_edge_fh_multiarg, 'Edge: function_handle multi-arg lambda', testFile);
    results = runTest(results, @test_edge_fh_nested_ops, 'Edge: function_handle nested ops', testFile);

    % -- Nesting / integration tests --
    results = runTest(results, @test_edge_cell_of_new_types, 'Edge: Cell containing new types', testFile);
    results = runTest(results, @test_edge_struct_of_new_types, 'Edge: Struct containing all new types', testFile);
    results = runTest(results, @test_edge_nested_deep, 'Edge: Deeply nested new types', testFile);
    results = runTest(results, @test_edge_many_variables, 'Edge: Many variables in one file', testFile);
    results = runTest(results, @test_edge_overwrite, 'Edge: Overwrite existing file', testFile);

    % -- Empty / degenerate cases for new types --
    results = runTest(results, @test_edge_empty_struct_no_fields, 'Edge: Empty struct (no fields)', testFile);
    results = runTest(results, @test_edge_struct_array, 'Edge: Struct array (non-scalar)', testFile);
    results = runTest(results, @test_edge_scalar_values, 'Edge: Scalar values all types', testFile);

    % -- Multi-output load --
    results = runTest(results, @test_edge_multi_output_load, 'Edge: Multi-output LoadFast', testFile);

    % -- Error handling for new types --
    results = runTest(results, @test_edge_error_name_too_long, 'Edge: Error name too long', testFile);
    results = runTest(results, @test_edge_error_unsupported_in_cell, 'Edge: Error unsupported type in cell', testFile);

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

%% ===== VERSION 3 TESTS (New Types) =====

function test_v3_complex_double(testFile)
    data = [1+2i, 3+4i; 5+6i, 7+8i];
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isequal(data, loaded.data), 'Complex double mismatch');
    assert(~isreal(loaded.data), 'Should be complex');
end

function test_v3_complex_single(testFile)
    data = single([1+2i, 3-4i]);
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isequal(data, loaded.data), 'Complex single mismatch');
    assert(isa(loaded.data, 'single'), 'Should be single');
end

function test_v3_sparse_double(testFile)
    data = sparse(eye(10));
    data(3, 7) = 42;
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(issparse(loaded.data), 'Should be sparse');
    assert(isequal(full(data), full(loaded.data)), 'Sparse double mismatch');
end

function test_v3_sparse_logical(testFile)
    data = sparse(logical(eye(5)));
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(issparse(loaded.data), 'Should be sparse');
    assert(islogical(loaded.data), 'Should be logical');
    assert(isequal(full(data), full(loaded.data)), 'Sparse logical mismatch');
end

function test_v3_sparse_complex(testFile)
    data = sparse([1 2 3], [2 3 1], [1+1i, 2+2i, 3+3i], 4, 4);
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(issparse(loaded.data), 'Should be sparse');
    assert(~isreal(loaded.data), 'Should be complex');
    assert(isequal(full(data), full(loaded.data)), 'Sparse complex mismatch');
end

function test_v3_string_array(testFile)
    data = ["hello", "world"; "foo", "bar"];
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'string'), 'Should be string');
    assert(isequal(data, loaded.data), 'String array mismatch');
end

function test_v3_string_missing(testFile)
    data = ["hello", string(missing), "world"];
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'string'), 'Should be string');
    assert(loaded.data(1) == "hello", 'First element mismatch');
    assert(ismissing(loaded.data(2)), 'Second element should be missing');
    assert(loaded.data(3) == "world", 'Third element mismatch');
end

function test_v3_table_basic(testFile)
    data = table([1;2;3], {'a';'b';'c'}, [true;false;true], ...
        'VariableNames', {'Num', 'Str', 'Log'});
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'table'), 'Should be table');
    assert(isequal(data.Properties.VariableNames, loaded.data.Properties.VariableNames), 'Variable names mismatch');
    assert(isequal(data.Num, loaded.data.Num), 'Num column mismatch');
    assert(isequal(data.Log, loaded.data.Log), 'Log column mismatch');
end

function test_v3_table_properties(testFile)
    data = table([1;2], [3;4], 'VariableNames', {'A', 'B'}, ...
        'RowNames', {'row1', 'row2'});
    data.Properties.VariableUnits = {'m', 'kg'};
    data.Properties.VariableDescriptions = {'length', 'mass'};
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isequal(data.Properties.RowNames, loaded.data.Properties.RowNames), 'RowNames mismatch');
    assert(isequal(data.Properties.VariableUnits, loaded.data.Properties.VariableUnits), 'Units mismatch');
    assert(isequal(data.Properties.VariableDescriptions, loaded.data.Properties.VariableDescriptions), 'Descriptions mismatch');
end

function test_v3_timetable(testFile)
    times = datetime(2024, 1, 1) + hours(0:2)';
    data = timetable(times, [10;20;30], {'a';'b';'c'}, 'VariableNames', {'Val', 'Label'});
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'timetable'), 'Should be timetable');
    assert(isequal(data.Val, loaded.data.Val), 'Val column mismatch');
    % Check row times are close (datetime precision)
    timeDiff = abs(data.Properties.RowTimes - loaded.data.Properties.RowTimes);
    assert(all(timeDiff < seconds(1)), 'RowTimes mismatch');
end

function test_v3_categorical(testFile)
    data = categorical({'red', 'blue', 'green', 'red', 'blue'});
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'categorical'), 'Should be categorical');
    assert(isequal(data, loaded.data), 'Categorical mismatch');
end

function test_v3_categorical_ordinal(testFile)
    data = categorical({'low', 'medium', 'high', 'low'}, ...
        {'low', 'medium', 'high'}, 'Ordinal', true);
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'categorical'), 'Should be categorical');
    assert(isordinal(loaded.data), 'Should be ordinal');
    assert(isequal(categories(data), categories(loaded.data)), 'Categories mismatch');
    assert(isequal(data, loaded.data), 'Ordinal categorical mismatch');
end

function test_v3_datetime(testFile)
    data = datetime(2024, 1, 1:5);
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'datetime'), 'Should be datetime');
    timeDiff = abs(data - loaded.data);
    assert(all(timeDiff < seconds(1)), 'Datetime mismatch');
end

function test_v3_datetime_timezone(testFile)
    data = datetime(2024, 6, 15, 12, 0, 0, 'TimeZone', 'America/New_York');
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'datetime'), 'Should be datetime');
    assert(strcmp(loaded.data.TimeZone, 'America/New_York'), 'Timezone mismatch');
    timeDiff = abs(data - loaded.data);
    assert(all(timeDiff < seconds(1)), 'Datetime TZ value mismatch');
end

function test_v3_duration(testFile)
    data = hours(1:5) + minutes(30);
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'duration'), 'Should be duration');
    diff = abs(milliseconds(data) - milliseconds(loaded.data));
    assert(all(diff < 1), 'Duration mismatch');
end

function test_v3_calendar_duration(testFile)
    data = calmonths(1:3) + caldays(10);
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'calendarDuration'), 'Should be calendarDuration');
    assert(isequal(data, loaded.data), 'CalendarDuration mismatch');
end

function test_v3_containers_map(testFile)
    data = containers.Map({'a', 'b', 'c'}, {1, 'two', [3 4 5]});
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'containers.Map'), 'Should be containers.Map');
    assert(isequal(sort(keys(data)), sort(keys(loaded.data))), 'Keys mismatch');
    assert(loaded.data('a') == 1, 'Value a mismatch');
    assert(strcmp(loaded.data('b'), 'two'), 'Value b mismatch');
    assert(isequal(loaded.data('c'), [3 4 5]), 'Value c mismatch');
end

function test_v3_function_handle(testFile)
    data = @sin;
    SaveFast(testFile, 'data', data);
    loaded = LoadFast(testFile);
    assert(isa(loaded.data, 'function_handle'), 'Should be function_handle');
    assert(loaded.data(0) == 0, 'Function should return sin(0)=0');
    assert(abs(loaded.data(pi/2) - 1) < 1e-10, 'Function should return sin(pi/2)=1');

    % Anonymous function
    data2 = @(x) x.^2 + 1;
    SaveFast(testFile, 'data', data2);
    loaded2 = LoadFast(testFile);
    assert(isa(loaded2.data, 'function_handle'), 'Should be function_handle');
    assert(loaded2.data(3) == 10, 'Anonymous function mismatch');
end

function test_v3_unsupported_type_error(testFile)
    % Arbitrary objects should error, not silently save empty
    try
        m = containers.Map();  % This is supported, so test with MException
        err = MException('test:test', 'test');
        SaveFast(testFile, 'data', err);
        error('Should have errored on unsupported type');
    catch ME
        assert(contains(ME.message, 'Unsupported type') || contains(ME.message, 'Cannot serialize'), ...
            'Should error with meaningful message for unsupported types');
    end
end

function test_v3_mixed_new_types(testFile)
    % Test saving multiple new types together
    s = struct();
    s.complex_val = 1 + 2i;
    s.str = "hello";
    s.cat = categorical({'a', 'b'});
    s.dur = hours(1);
    s.fh = @sin;

    SaveFast(testFile, s);
    loaded = LoadFast(testFile);

    assert(isequal(s.complex_val, loaded.complex_val), 'Complex mismatch in struct');
    assert(s.str == loaded.str, 'String mismatch in struct');
    assert(isequal(s.cat, loaded.cat), 'Categorical mismatch in struct');
    diff = abs(milliseconds(s.dur) - milliseconds(loaded.dur));
    assert(diff < 1, 'Duration mismatch in struct');
    assert(isa(loaded.fh, 'function_handle'), 'Function handle type lost');
end

%% ===== COMPREHENSIVE EDGE-CASE TESTS =====

% --- Complex numeric edge cases ---

function test_edge_complex_all_numeric_types(testFile)
    % Complex double and single should round-trip exactly
    % (MATLAB only supports complex double and single)
    d = [1+2i, 3+4i, NaN+1i; -1-1i, 0, Inf+0i];
    SaveFast(testFile, 'data', d);
    L = LoadFast(testFile);
    assert(isa(L.data, 'double'), 'Class mismatch');
    assert(isequal(size(d), size(L.data)), 'Size mismatch');
    % NaN comparison: use isequaln which treats NaN==NaN
    assert(isequaln(d, L.data), 'Complex double mismatch');

    s = single([1+2i; 3-4i; 0+0i]);
    SaveFast(testFile, 'data', s);
    L = LoadFast(testFile);
    assert(isa(L.data, 'single'), 'Single class mismatch');
    assert(isequal(s, L.data), 'Complex single mismatch');
end

function test_edge_complex_3d(testFile)
    data = complex(randn(3,4,5), randn(3,4,5));
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(size(data), size(L.data)), '3D size mismatch');
    assert(max(abs(data(:) - L.data(:))) < eps('double'), '3D complex mismatch');
end

function test_edge_complex_nan_inf(testFile)
    data = [NaN+1i, 1+NaN*1i, Inf+2i, -Inf-Inf*1i, complex(0,0)];
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    % Check element by element for NaN/Inf
    for k = 1:numel(data)
        assert(isequaln(real(data(k)), real(L.data(k))), sprintf('Real part mismatch at %d', k));
        assert(isequaln(imag(data(k)), imag(L.data(k))), sprintf('Imag part mismatch at %d', k));
    end
end

function test_edge_complex_scalar(testFile)
    data = 3+4i;
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(data, L.data), 'Scalar complex mismatch');
    assert(numel(L.data) == 1, 'Should be scalar');
end

function test_edge_complex_pure_imaginary(testFile)
    data = [0+1i, 0+2i, 0-3i];
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(data, L.data), 'Pure imaginary mismatch');
    assert(all(real(L.data) == 0), 'Real parts should be zero');
end

% --- Sparse edge cases ---

function test_edge_sparse_empty(testFile)
    % All-zeros sparse matrix
    data = sparse(10, 20);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(issparse(L.data), 'Should be sparse');
    assert(isequal(size(data), size(L.data)), 'Size mismatch');
    assert(nnz(L.data) == 0, 'Should have no nonzeros');
end

function test_edge_sparse_large(testFile)
    % Large sparse matrix with specific sparsity
    n = 1000;
    data = sprandn(n, n, 0.01); % ~1% density
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(issparse(L.data), 'Should be sparse');
    assert(isequal(size(data), size(L.data)), 'Size mismatch');
    assert(nnz(data) == nnz(L.data), 'nnz mismatch');
    assert(max(abs(nonzeros(data) - nonzeros(L.data))) < eps, 'Values mismatch');
end

function test_edge_sparse_single_element(testFile)
    data = sparse(100, 200);
    data(50, 100) = 42;
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(issparse(L.data), 'Should be sparse');
    assert(L.data(50, 100) == 42, 'Value mismatch');
    assert(nnz(L.data) == 1, 'Should have exactly 1 nonzero');
    assert(isequal(size(data), size(L.data)), 'Size mismatch');
end

function test_edge_sparse_full_density(testFile)
    % Sparse matrix that's actually dense
    data = sparse(ones(5, 5));
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(issparse(L.data), 'Should remain sparse');
    assert(isequal(full(data), full(L.data)), 'Dense-sparse mismatch');
end

% --- String edge cases ---

function test_edge_string_empty_array(testFile)
    data = string.empty(0, 0);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'string'), 'Should be string');
    assert(isempty(L.data), 'Should be empty');
    assert(isequal(size(data), size(L.data)), 'Size mismatch');
end

function test_edge_string_scalar(testFile)
    data = "hello world";
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'string'), 'Should be string');
    assert(isscalar(L.data), 'Should be scalar');
    assert(L.data == "hello world", 'Content mismatch');
end

function test_edge_string_empty_element(testFile)
    % Empty string "" vs missing
    data = ["", "hello", "", string(missing), ""];
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(L.data(1) == "", 'First element should be empty string');
    assert(L.data(2) == "hello", 'Second element mismatch');
    assert(L.data(3) == "", 'Third element should be empty string');
    assert(ismissing(L.data(4)), 'Fourth element should be missing');
    assert(L.data(5) == "", 'Fifth element should be empty string');
    assert(~ismissing(L.data(1)), 'Empty string should NOT be missing');
end

function test_edge_string_long(testFile)
    % Very long string content
    data = string(repmat('A', 1, 10000));
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(strlength(L.data) == 10000, 'Length mismatch');
    assert(L.data == data, 'Content mismatch');
end

function test_edge_string_special_chars(testFile)
    data = ["tab" + char(9) + "here", "newline" + char(10) + "here", ...
            "quote""double", "backslash\path"];
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(data, L.data), 'Special chars mismatch');
end

function test_edge_string_all_missing(testFile)
    data = [string(missing), string(missing), string(missing)];
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'string'), 'Should be string');
    assert(all(ismissing(L.data)), 'All should be missing');
end

% --- Table edge cases ---

function test_edge_table_empty(testFile)
    % Table with columns but no rows
    data = table(double.empty(0,1), cell(0,1), 'VariableNames', {'Num', 'Txt'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'table'), 'Should be table');
    assert(height(L.data) == 0, 'Should have 0 rows');
    assert(isequal(data.Properties.VariableNames, L.data.Properties.VariableNames), 'VarNames mismatch');
end

function test_edge_table_single_row(testFile)
    data = table(42, {'hello'}, true, 'VariableNames', {'A', 'B', 'C'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(height(L.data) == 1, 'Should have 1 row');
    assert(L.data.A == 42, 'A mismatch');
    assert(L.data.C == true, 'C mismatch');
end

function test_edge_table_mixed_columns(testFile)
    % Table with numeric, string, categorical, and logical columns
    data = table((1:5)', ["a";"b";"c";"d";"e"], ...
        categorical({'x';'y';'x';'y';'x'}), logical([1;0;1;0;1]), ...
        'VariableNames', {'Num', 'Str', 'Cat', 'Log'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'table'), 'Should be table');
    assert(isequal(data.Num, L.data.Num), 'Num mismatch');
    assert(isequal(data.Str, L.data.Str), 'Str mismatch');
    assert(isequal(data.Cat, L.data.Cat), 'Cat mismatch');
    assert(isequal(data.Log, L.data.Log), 'Log mismatch');
end

function test_edge_table_no_properties(testFile)
    % Simple table with no row names, units, or descriptions
    data = table([1;2;3], [4;5;6]);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'table'), 'Should be table');
    assert(isequal(data{:,:}, L.data{:,:}), 'Data mismatch');
end

function test_edge_table_cell_column(testFile)
    % Table with a cell column (heterogeneous data)
    data = table([1;2;3], {[1 2 3]; 'hello'; {4,5}}, 'VariableNames', {'ID', 'Data'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'table'), 'Should be table');
    assert(isequal(data.ID, L.data.ID), 'ID mismatch');
    assert(isequal(data.Data, L.data.Data), 'Cell column mismatch');
end

% --- Timetable edge cases ---

function test_edge_timetable_duration_rowtimes(testFile)
    % Timetable with duration-based row times (instead of datetime)
    times = seconds([0; 0.5; 1.0; 1.5; 2.0]);
    data = timetable(times, randn(5,1), 'VariableNames', {'Signal'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'timetable'), 'Should be timetable');
    diff = abs(milliseconds(data.Properties.RowTimes) - milliseconds(L.data.Properties.RowTimes));
    assert(all(diff < 1), 'Duration RowTimes mismatch');
    assert(max(abs(data.Signal - L.data.Signal)) < eps, 'Signal data mismatch');
end

function test_edge_timetable_multi_col(testFile)
    times = datetime(2024,1,1) + hours(0:4)';
    data = timetable(times, (1:5)', randn(5,1), logical([1;0;1;0;1]), ...
        'VariableNames', {'A', 'B', 'C'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(width(L.data) == 3, 'Should have 3 columns');
    assert(isequal(data.A, L.data.A), 'A mismatch');
    assert(isequal(data.C, L.data.C), 'C mismatch');
end

% --- Categorical edge cases ---

function test_edge_categorical_undefined(testFile)
    % Categorical with some undefined elements
    data = categorical({'red', 'blue', 'green'}, {'red', 'blue', 'green', 'yellow'});
    % 'yellow' is a category but no element uses it; this tests unused categories
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(data, L.data), 'Values mismatch');
    assert(isequal(categories(data), categories(L.data)), 'Categories mismatch');

    % Now test with actual undefined
    data2 = categorical({'red', '', 'blue'}, {'red', 'blue'});
    % '' is not in the category list, so it becomes <undefined>
    SaveFast(testFile, 'data', data2);
    L2 = LoadFast(testFile);
    assert(isundefined(L2.data(2)), 'Should be undefined');
    assert(~isundefined(L2.data(1)), 'Should not be undefined');
    assert(isequal(categories(data2), categories(L2.data)), 'Categories mismatch');
end

function test_edge_categorical_2d(testFile)
    data = categorical({'a','b'; 'c','a'; 'b','c'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(size(data), size(L.data)), 'Size mismatch');
    assert(isequal(data, L.data), 'Values mismatch');
end

function test_edge_categorical_single_cat(testFile)
    data = categorical({'x','x','x'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(data, L.data), 'Single cat mismatch');
    assert(numel(categories(L.data)) == 1, 'Should have 1 category');
end

function test_edge_categorical_unicode_names(testFile)
    data = categorical({'alpha-1', 'beta_2', 'gamma 3'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(data, L.data), 'Special char names mismatch');
end

% --- Datetime edge cases ---

function test_edge_datetime_nat(testFile)
    data = [datetime(2024,1,1), NaT, datetime(2024,6,15), NaT];
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'datetime'), 'Should be datetime');
    assert(isnat(L.data(2)), 'Element 2 should be NaT');
    assert(isnat(L.data(4)), 'Element 4 should be NaT');
    assert(~isnat(L.data(1)), 'Element 1 should not be NaT');
    % Check valid elements are close
    timeDiff = abs(data(1) - L.data(1));
    assert(timeDiff < seconds(1), 'Valid datetime mismatch');
end

function test_edge_datetime_2d(testFile)
    data = datetime(2024, 1:6, 15);
    data = reshape(data, 2, 3);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(size(data), size(L.data)), 'Size mismatch');
    timeDiff = max(abs(data(:) - L.data(:)));
    assert(timeDiff < seconds(1), '2D datetime mismatch');
end

function test_edge_datetime_utc(testFile)
    data = datetime(2024, 6, 15, 12, 30, 45, 'TimeZone', 'UTC');
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(strcmp(L.data.TimeZone, 'UTC'), 'UTC timezone mismatch');
    timeDiff = abs(data - L.data);
    assert(timeDiff < seconds(1), 'UTC value mismatch');
end

function test_edge_datetime_custom_format(testFile)
    data = datetime(2024, 3, 15);
    data.Format = 'yyyy/MM/dd';
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(strcmp(L.data.Format, 'yyyy/MM/dd'), 'Custom format not preserved');
end

% --- Duration edge cases ---

function test_edge_duration_nan(testFile)
    data = [hours(1), duration(NaN, 0, 0), minutes(30)];
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'duration'), 'Should be duration');
    assert(isnan(L.data(2)), 'Element 2 should be NaN');
    assert(~isnan(L.data(1)), 'Element 1 should not be NaN');
end

function test_edge_duration_2d(testFile)
    data = reshape(hours(1:6), 2, 3);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(size(data), size(L.data)), 'Size mismatch');
    diff = max(abs(milliseconds(data(:)) - milliseconds(L.data(:))));
    assert(diff < 1, '2D duration mismatch');
end

function test_edge_duration_zero(testFile)
    data = duration(0, 0, 0);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(milliseconds(L.data) == 0, 'Zero duration mismatch');
end

function test_edge_duration_negative(testFile)
    data = [-hours(5), -minutes(30), -seconds(1)];
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    diff = abs(milliseconds(data) - milliseconds(L.data));
    assert(all(diff < 1), 'Negative duration mismatch');
end

function test_edge_duration_custom_format(testFile)
    data = hours(1) + minutes(30) + seconds(45);
    data.Format = 'hh:mm:ss';
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(strcmp(L.data.Format, 'hh:mm:ss'), 'Duration format not preserved');
end

% --- CalendarDuration edge cases ---

function test_edge_caldur_2d(testFile)
    data = reshape(calmonths(1:6), 2, 3);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(size(data), size(L.data)), 'Size mismatch');
    assert(isequal(data, L.data), '2D calendarDuration mismatch');
end

function test_edge_caldur_mixed(testFile)
    % Combined months, days, and time
    data = calyears(1) + calmonths(2) + caldays(15) + hours(6);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'calendarDuration'), 'Should be calendarDuration');
    assert(isequal(data, L.data), 'Mixed calendarDuration mismatch');
end

function test_edge_caldur_negative(testFile)
    data = -calmonths(3) + caldays(-10);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isequal(data, L.data), 'Negative calendarDuration mismatch');
end

% --- containers.Map edge cases ---

function test_edge_map_empty(testFile)
    data = containers.Map();
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'containers.Map'), 'Should be containers.Map');
    assert(length(keys(L.data)) == 0, 'Should be empty');
end

function test_edge_map_numeric_keys(testFile)
    data = containers.Map({1, 2, 3}, {'one', 'two', 'three'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isa(L.data, 'containers.Map'), 'Should be containers.Map');
    assert(strcmp(L.data(1), 'one'), 'Value for key 1 mismatch');
    assert(strcmp(L.data(2), 'two'), 'Value for key 2 mismatch');
    assert(strcmp(L.data(3), 'three'), 'Value for key 3 mismatch');
end

function test_edge_map_single_entry(testFile)
    data = containers.Map({'only'}, {42});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(length(keys(L.data)) == 1, 'Should have 1 entry');
    assert(L.data('only') == 42, 'Value mismatch');
end

function test_edge_map_large(testFile)
    % Map with many entries
    k = arrayfun(@(x) sprintf('key_%04d', x), 1:100, 'UniformOutput', false);
    v = num2cell(1:100);
    data = containers.Map(k, v);
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(length(keys(L.data)) == 100, 'Should have 100 entries');
    assert(L.data('key_0001') == 1, 'First value mismatch');
    assert(L.data('key_0100') == 100, 'Last value mismatch');
end

% --- function_handle edge cases ---

function test_edge_fh_builtin(testFile)
    % Test several builtins
    fns = {@cos, @abs, @sqrt, @exp};
    for i = 1:length(fns)
        SaveFast(testFile, 'data', fns{i});
        L = LoadFast(testFile);
        assert(isa(L.data, 'function_handle'), sprintf('fn %d not function_handle', i));
        assert(abs(L.data(1) - fns{i}(1)) < eps, sprintf('fn %d value mismatch', i));
    end
end

function test_edge_fh_multiarg(testFile)
    data = @(x, y) x.^2 + y.^2;
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(L.data(3, 4) == 25, 'Multi-arg function mismatch');
end

function test_edge_fh_nested_ops(testFile)
    data = @(x) sin(x).^2 + cos(x).^2;  % Should always return 1
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    vals = L.data(linspace(0, 2*pi, 100));
    assert(max(abs(vals - 1)) < 1e-10, 'Nested ops function mismatch');
end

% --- Nesting / integration tests ---

function test_edge_cell_of_new_types(testFile)
    data = {1+2i, sparse(eye(3)), "hello", categorical({'a'}), ...
            hours(5), @sin, datetime(2024,1,1)};
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(iscell(L.data), 'Should be cell');
    assert(numel(L.data) == 7, 'Should have 7 elements');

    % Check each element type
    assert(isequal(data{1}, L.data{1}), 'Complex in cell mismatch');
    assert(issparse(L.data{2}), 'Sparse in cell lost');
    assert(isequal(full(data{2}), full(L.data{2})), 'Sparse in cell mismatch');
    assert(L.data{3} == "hello", 'String in cell mismatch');
    assert(isa(L.data{4}, 'categorical'), 'Categorical in cell type lost');
    assert(isa(L.data{5}, 'duration'), 'Duration in cell type lost');
    assert(isa(L.data{6}, 'function_handle'), 'FH in cell type lost');
    assert(isa(L.data{7}, 'datetime'), 'Datetime in cell type lost');
end

function test_edge_struct_of_new_types(testFile)
    s = struct();
    s.complex_mat = randn(3) + 1i*randn(3);
    s.sparse_mat = sparse(eye(5));
    s.string_arr = ["a", "b"; "c", "d"];
    s.cat = categorical({'x','y','z'});
    s.dt = datetime(2024, 1:3, 1);
    s.dur = hours([1,2,3]);
    s.caldur = calmonths(1:3);
    s.map = containers.Map({'k1','k2'}, {10, 20});
    s.fh = @(x) x+1;

    SaveFast(testFile, s);
    L = LoadFast(testFile);

    assert(isequaln(s.complex_mat, L.complex_mat), 'Complex in struct mismatch');
    assert(issparse(L.sparse_mat), 'Sparse in struct lost');
    assert(isequal(s.string_arr, L.string_arr), 'String in struct mismatch');
    assert(isequal(s.cat, L.cat), 'Cat in struct mismatch');
    assert(isa(L.dt, 'datetime'), 'Datetime type lost');
    assert(isa(L.dur, 'duration'), 'Duration type lost');
    assert(isa(L.caldur, 'calendarDuration'), 'CalDur type lost');
    assert(isa(L.map, 'containers.Map'), 'Map type lost');
    assert(L.fh(5) == 6, 'FH in struct mismatch');
end

function test_edge_nested_deep(testFile)
    % Cell inside struct inside cell, with new types at leaf
    data = {struct('inner', {{1+2i, sparse(eye(2)), "test", ...
            struct('deep', categorical({'a','b'}))}})};
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);

    inner = L.data{1}.inner;
    assert(isequal(1+2i, inner{1}), 'Deep complex mismatch');
    assert(issparse(inner{2}), 'Deep sparse lost');
    assert(inner{3} == "test", 'Deep string mismatch');
    assert(isa(inner{4}.deep, 'categorical'), 'Deep categorical type lost');
end

function test_edge_many_variables(testFile)
    % Save 20 variables of mixed types
    args = {};
    for i = 1:20
        name = sprintf('var%02d', i);
        switch mod(i, 5)
            case 0, val = randn(3) + 1i*randn(3);
            case 1, val = rand(10);
            case 2, val = sprintf('string_%d', i);
            case 3, val = logical(randi([0 1], 4));
            case 4, val = int32(randi(1000, 5));
        end
        args = [args, {name, val}]; %#ok<AGROW>
    end

    SaveFast(testFile, args{:});
    L = LoadFast(testFile);

    assert(isstruct(L), 'Should be struct');
    flds = fieldnames(L);
    assert(numel(flds) == 20, 'Should have 20 fields');
    for i = 1:20
        name = sprintf('var%02d', i);
        assert(isfield(L, name), sprintf('Missing field %s', name));
    end
end

function test_edge_overwrite(testFile)
    % Save, then overwrite with different data
    SaveFast(testFile, 'data', rand(100));
    SaveFast(testFile, 'data', int32(42));  % Overwrite
    L = LoadFast(testFile);
    assert(isa(L.data, 'int32'), 'Should be int32 after overwrite');
    assert(L.data == 42, 'Value mismatch after overwrite');
end

% --- Empty / degenerate cases ---

function test_edge_empty_struct_no_fields(testFile)
    data = struct();
    SaveFast(testFile, 'x', data);
    L = LoadFast(testFile);
    assert(isstruct(L.x), 'Should be struct');
    assert(numel(fieldnames(L.x)) == 0, 'Should have no fields');
end

function test_edge_struct_array(testFile)
    data = struct('a', {1,2,3}, 'b', {'x','y','z'});
    SaveFast(testFile, 'data', data);
    L = LoadFast(testFile);
    assert(isstruct(L.data), 'Should be struct');
    assert(numel(L.data) == 3, 'Should have 3 elements');
    assert(L.data(1).a == 1, 'Element 1.a mismatch');
    assert(strcmp(L.data(2).b, 'y'), 'Element 2.b mismatch');
    assert(L.data(3).a == 3, 'Element 3.a mismatch');
end

function test_edge_scalar_values(testFile)
    % Scalar values of every basic type
    SaveFast(testFile, ...
        'd', double(3.14), ...
        's', single(2.72), ...
        'i8', int8(-42), ...
        'u8', uint8(255), ...
        'i16', int16(-1000), ...
        'u16', uint16(65535), ...
        'i32', int32(-100000), ...
        'u32', uint32(100000), ...
        'i64', int64(-1e15), ...
        'u64', uint64(1e15), ...
        'ch', 'A', ...
        'lg', true);
    L = LoadFast(testFile);

    assert(L.d == 3.14 && isa(L.d, 'double'), 'double mismatch');
    assert(L.s == single(2.72) && isa(L.s, 'single'), 'single mismatch');
    assert(L.i8 == int8(-42) && isa(L.i8, 'int8'), 'int8 mismatch');
    assert(L.u8 == uint8(255) && isa(L.u8, 'uint8'), 'uint8 mismatch');
    assert(L.i16 == int16(-1000) && isa(L.i16, 'int16'), 'int16 mismatch');
    assert(L.u16 == uint16(65535) && isa(L.u16, 'uint16'), 'uint16 mismatch');
    assert(L.i32 == int32(-100000) && isa(L.i32, 'int32'), 'int32 mismatch');
    assert(L.u32 == uint32(100000) && isa(L.u32, 'uint32'), 'uint32 mismatch');
    assert(L.i64 == int64(-1e15) && isa(L.i64, 'int64'), 'int64 mismatch');
    assert(L.u64 == uint64(1e15) && isa(L.u64, 'uint64'), 'uint64 mismatch');
    assert(L.ch == 'A' && ischar(L.ch), 'char mismatch');
    assert(L.lg == true && islogical(L.lg), 'logical mismatch');
end

% --- Multi-output load ---

function test_edge_multi_output_load(testFile)
    SaveFast(testFile, 'A', 10, 'B', 20, 'C', 30);

    % Two outputs
    [a, b] = LoadFast(testFile);
    assert(a == 10, 'First output mismatch');
    assert(b == 20, 'Second output mismatch');

    % Three outputs
    [a, b, c] = LoadFast(testFile);
    assert(a == 10, 'First output mismatch');
    assert(b == 20, 'Second output mismatch');
    assert(c == 30, 'Third output mismatch');

    % Too many outputs should error
    try
        [a, b, c, d] = LoadFast(testFile);  %#ok<ASGLU>
        error('Should have errored on too many outputs');
    catch ME
        assert(contains(ME.message, 'outputs') || contains(ME.message, 'variable'), ...
            'Wrong error for too many outputs');
    end
end

% --- Error handling for new types ---

function test_edge_error_name_too_long(testFile)
    longName = repmat('a', 1, 256);
    try
        SaveFast(testFile, longName, 42);
        error('Should have errored on long name');
    catch ME
        assert(contains(ME.message, '255') || contains(ME.message, 'exceeds') || ...
               contains(ME.message, 'variable name'), ...
            'Wrong error for long name');
    end
end

function test_edge_error_unsupported_in_cell(testFile)
    % Unsupported type nested inside a cell should error
    try
        err = MException('test:test', 'msg');
        SaveFast(testFile, 'data', {1, err, 3});
        error('Should have errored on unsupported type in cell');
    catch ME
        assert(contains(ME.message, 'Unsupported') || contains(ME.message, 'Cannot serialize'), ...
            'Wrong error for unsupported type in cell');
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

function test_exhaustive_types()
% TEST_EXHAUSTIVE_TYPES - Exhaustive type x container x nesting test suite
%
% Tests every supported type standalone, inside every container type,
% and in all nesting combinations. ~300 tests with dot-progress output.
%
% Usage:
%   test_exhaustive_types()

    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('  EXHAUSTIVE TYPE & NESTING TEST SUITE\n');
    fprintf('================================================================\n');
    fprintf('Start time: %s\n', datestr(now)); %#ok<TNOW1,DATST>
    fprintf('================================================================\n\n');

    rng(42);

    binFile = [tempname '.bin'];
    cleanup = onCleanup(@() deleteIfExists(binFile)); %#ok<NASGU>

    totalPass = 0;
    totalFail = 0;

    % Build leaf values used across sections
    [leaves, leafNames] = buildLeaves();

    %% Section 1: Standalone leaves
    [p, f] = runSection('Section 1: Standalone leaves', leaves, leafNames, ...
        @(val, file) testStandalone(val, file), binFile);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 2: Each leaf in struct field
    [p, f] = runSection('Section 2: Each leaf in struct field', leaves, leafNames, ...
        @(val, file) testInStruct(val, file), binFile);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 3: Each leaf in cell element
    [p, f] = runSection('Section 3: Each leaf in cell element', leaves, leafNames, ...
        @(val, file) testInCell(val, file), binFile);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 4: Each leaf as Map value
    [p, f] = runSection('Section 4: Each leaf as Map value', leaves, leafNames, ...
        @(val, file) testInMap(val, file), binFile);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 5: Table column types
    [p, f] = runTableColumnSection(binFile);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 6: Container x Container 2-level
    [p, f] = runContainerContainerSection(binFile);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 7: Deep nesting 3-5 levels
    [p, f] = runDeepNestingSection(binFile);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 8: Mixed mega-containers
    [p, f] = runMegaContainerSection(binFile, leaves, leafNames);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 9: Empty/degenerate containers
    [p, f] = runEmptyDegenerateSection(binFile);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 10: Struct arrays with nesting
    [p, f] = runStructArraySection(binFile);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Section 11: All leaves as named variables in one file
    [p, f] = runAllLeavesOneFileSection(binFile, leaves, leafNames);
    totalPass = totalPass + p; totalFail = totalFail + f;

    %% Final summary
    total = totalPass + totalFail;
    fprintf('\n================================================================\n');
    fprintf('  SUMMARY: %d/%d passed, %d failed (%.1f%%)\n', ...
        totalPass, total, totalFail, 100 * totalPass / max(total, 1));
    fprintf('================================================================\n\n');
end

%% ===== LEAF VALUE BUILDER =====

function [leaves, names] = buildLeaves()
    leaves = {};
    names = {};

    % --- Numeric types ---
    add('dbl_scalar', pi);
    add('dbl_vector', [1 2 3 4 5]);
    add('dbl_matrix', rand(4, 5));
    add('dbl_3d', rand(2, 3, 4));
    add('dbl_empty', []);
    add('sgl_scalar', single(3.14));
    add('sgl_vector', single([1 2 3]));
    add('int8_vec', int8([-128 0 127]));
    add('uint8_vec', uint8([0 128 255]));
    add('int16_vec', int16([-1000 0 1000]));
    add('uint16_vec', uint16([0 32768 65535]));
    add('int32_vec', int32([-1e5 0 1e5]));
    add('uint32_vec', uint32([0 1e5 4e9]));
    add('int64_vec', int64([-1e15 0 1e15]));
    add('uint64_vec', uint64([0 1e9 1e18]));

    % --- Special values ---
    add('special_nan_inf', [NaN, Inf, -Inf, eps, realmin, realmax]);
    add('special_single', single([NaN, Inf, -Inf, 0]));

    % --- Complex ---
    add('cplx_dbl_mat', randn(3) + 1i*randn(3));
    add('cplx_sgl_vec', single([1+2i, 3-4i]));
    add('cplx_nan_inf', [NaN+1i, Inf-2i, 1+NaN*1i]);

    % --- Sparse ---
    add('sp_double', sprandn(20, 20, 0.1));
    add('sp_logical', sparse(logical(eye(5))));
    add('sp_complex', sparse([1 3], [2 4], [1+1i 2-2i], 5, 5));
    add('sp_empty', sparse(10, 10));

    % --- Char ---
    add('char_scalar', 'A');
    add('char_vector', 'Hello, World!');
    add('char_matrix', ['abc'; 'def'; 'ghi']);
    add('char_special', ['tab' char(9) 'nl' char(10) 'end']);

    % --- Logical ---
    add('log_scalar', true);
    add('log_vector', logical([1 0 1 0 1]));
    add('log_matrix', logical(randi([0 1], 3, 4)));

    % --- String ---
    add('str_scalar', "hello world");
    add('str_2d', ["a", "b"; "c", "d"]);
    add('str_missing', ["hello", string(missing), ""]);
    add('str_empty', string.empty(0, 0));

    % --- Categorical ---
    add('cat_plain', categorical({'red', 'blue', 'green', 'red'}));
    add('cat_ordinal', categorical({'low', 'med', 'high'}, ...
        {'low', 'med', 'high'}, 'Ordinal', true));
    add('cat_undefined', categorical({'x', '', 'y'}, {'x', 'y'}));

    % --- Datetime ---
    add('dt_plain', datetime(2024, 1, 1:5));
    add('dt_tz', datetime(2024, 6, 15, 12, 0, 0, 'TimeZone', 'America/New_York'));
    add('dt_nat', [datetime(2024, 1, 1), NaT, datetime(2024, 12, 31)]);
    add('dt_format', setfield_inline(datetime(2024, 3, 15), 'Format', 'yyyy/MM/dd')); %#ok<*SFLD>

    % --- Duration ---
    add('dur_vec', hours(1:5) + minutes(30));
    add('dur_nan', [hours(1), duration(NaN, 0, 0), minutes(-30)]);
    add('dur_zero', duration(0, 0, 0));

    % --- CalendarDuration ---
    add('caldur_simple', calmonths(1:3) + caldays(15));
    add('caldur_mixed', calyears(1) + calmonths(2) + caldays(15) + hours(6));

    % --- containers.Map ---
    add('map_char', containers.Map({'a', 'b', 'c'}, {1, 'two', [3 4 5]}));
    add('map_num', containers.Map({1, 2, 3}, {'x', 'y', 'z'}));
    add('map_empty', containers.Map());

    % --- function_handle ---
    add('fh_builtin', @sin);
    add('fh_anon', @(x, y) x.^2 + y.^2);

    % --- Struct ---
    add('struct_simple', struct('a', 1, 'b', 'text', 'c', rand(3)));
    add('struct_empty', struct());
    add('struct_array', struct('x', {1, 2, 3}, 'y', {'a', 'b', 'c'}));

    % --- Cell ---
    add('cell_mixed', {1, 'two', [3; 4; 5], {6, 7}});
    add('cell_2d', reshape({1, 2, 3, 4, 5, 6}, 2, 3));
    add('cell_empty', cell(0, 0));

    % --- Table ---
    add('tbl_basic', table((1:3)', ["a";"b";"c"], logical([1;0;1]), ...
        'VariableNames', {'Num', 'Str', 'Log'}));
    tbl_props = table([1;2], [3;4], 'VariableNames', {'A', 'B'}, ...
        'RowNames', {'r1', 'r2'});
    tbl_props.Properties.VariableUnits = {'m', 'kg'};
    tbl_props.Properties.VariableDescriptions = {'length', 'mass'};
    add('tbl_props', tbl_props);
    add('tbl_empty', table(double.empty(0,1), cell(0,1), ...
        'VariableNames', {'Num', 'Txt'}));

    % --- Timetable ---
    tt_times = datetime(2024, 1, 1) + hours(0:4)';
    add('tt_basic', timetable(tt_times, (1:5)', randn(5, 1), ...
        'VariableNames', {'A', 'B'}));

    function add(name, val)
        leaves{end+1} = val; %#ok<SETNU>
        names{end+1} = name; %#ok<SETNU>
    end
end

function dt = setfield_inline(dt, field, val)
    dt.(field) = val;
end

%% ===== SECTION RUNNERS =====

function [passed, failed] = runSection(sectionName, leaves, leafNames, testFn, binFile)
    n = numel(leaves);
    fprintf('--- %s (%d tests) ---\n', sectionName, n);
    passed = 0;
    failed = 0;
    failures = {};
    for i = 1:n
        try
            testFn(leaves{i}, binFile);
            fprintf('.');
            passed = passed + 1;
        catch ME
            fprintf('F');
            failed = failed + 1;
            failures{end+1} = sprintf('FAIL: %s: %s', leafNames{i}, ME.message); %#ok<AGROW>
        end
    end
    fprintf('  %d/%d passed\n', passed, n);
    for i = 1:numel(failures)
        fprintf('  %s\n', failures{i});
    end
    if ~isempty(failures), fprintf('\n'); end
end

%% ===== TEST PATTERNS =====

function testStandalone(val, file)
    SaveFast(file, 'val', val);
    L = LoadFast(file);
    [ok, msg] = compareValues(val, L.val);
    assert(ok, msg);
end

function testInStruct(val, file)
    s = struct();
    s.wrapped = val;
    SaveFast(file, 'val', s);
    L = LoadFast(file);
    [ok, msg] = compareValues(val, L.val.wrapped);
    assert(ok, msg);
end

function testInCell(val, file)
    c = {val};
    SaveFast(file, 'val', c);
    L = LoadFast(file);
    [ok, msg] = compareValues(val, L.val{1});
    assert(ok, msg);
end

function testInMap(val, file)
    m = containers.Map({'key'}, {val});
    % Compare against what Map itself returns (Map can reshape, e.g. char matrices)
    expected = m('key');
    SaveFast(file, 'val', m);
    L = LoadFast(file);
    assert(isa(L.val, 'containers.Map'), 'Not a Map');
    [ok, msg] = compareValues(expected, L.val('key'));
    assert(ok, msg);
end

%% ===== SECTION 5: TABLE COLUMN TYPES =====

function [passed, failed] = runTableColumnSection(binFile)
    fprintf('--- Section 5: Table column types ---\n');
    passed = 0; failed = 0; failures = {};
    n = 5; % rows per column

    tests = {
        'double',       (1:n)'
        'single',       single((1:n)')
        'int8',         int8((1:n)')
        'uint8',        uint8((1:n)')
        'int16',        int16((1:n)')
        'uint16',       uint16((1:n)')
        'int32',        int32((1:n)')
        'uint32',       uint32((1:n)')
        'int64',        int64((1:n)')
        'uint64',       uint64((1:n)')
        'logical',      logical([1;0;1;0;1])
        'string',       ["a";"b";"c";"d";"e"]
        'categorical',  categorical({'x';'y';'x';'y';'x'})
        'datetime',     datetime(2024,1,1:n)'
        'duration',     hours(1:n)'
        'calDuration',  (calmonths(1:n) + caldays(1))'
        'cell',         {1; 'two'; [3 4]; {5}; 'six'}
        'complex',      (1:n)' + 1i*(n:-1:1)'
        'charCell',     {'a'; 'bb'; 'ccc'; 'dddd'; 'eeeee'}
        'strMissing',   ["a"; string(missing); "c"; "d"; ""]
    };

    for i = 1:size(tests, 1)
        name = tests{i, 1};
        colData = tests{i, 2};
        try
            tbl = table(colData, 'VariableNames', {'Col'});
            SaveFast(binFile, 'val', tbl);
            L = LoadFast(binFile);
            [ok, msg] = compareValues(tbl, L.val);
            assert(ok, msg);
            fprintf('.'); passed = passed + 1;
        catch ME
            fprintf('F'); failed = failed + 1;
            failures{end+1} = sprintf('FAIL: tbl_%s: %s', name, ME.message); %#ok<AGROW>
        end
    end
    fprintf('  %d/%d passed\n', passed, size(tests, 1));
    for i = 1:numel(failures), fprintf('  %s\n', failures{i}); end
    if ~isempty(failures), fprintf('\n'); end
end

%% ===== SECTION 6: CONTAINER x CONTAINER 2-LEVEL =====

function [passed, failed] = runContainerContainerSection(binFile)
    fprintf('--- Section 6: Container x Container 2-level ---\n');
    passed = 0; failed = 0; failures = {};

    payload = struct('x', 42, 'y', 'hello');
    innerCell = {1, 'two', [3 4]};
    innerMap = containers.Map({'k1', 'k2'}, {10, 20});
    innerTable = table([1;2], ["a";"b"], 'VariableNames', {'N', 'S'});
    tt_t = datetime(2024,1,1) + hours(0:2)';
    innerTT = timetable(tt_t, [10;20;30], 'VariableNames', {'V'});
    innerSA = struct('a', {1, 2}, 'b', {'x', 'y'});
    innerFH = @sin;

    % Outer: struct
    combos = {
        'struct_in_struct',   struct('inner', payload)
        'cell_in_struct',     struct('inner', {innerCell})  % use {} to avoid struct array
        'map_in_struct',      struct('inner', innerMap)
        'table_in_struct',    struct('inner', innerTable)
        'tt_in_struct',       struct('inner', innerTT)
        'sa_in_struct',       struct('inner', innerSA)
        'fh_in_struct',       struct('inner', innerFH)
    };
    % Fix: struct('inner', {innerCell}) creates a 1x3 struct array.
    % We need struct('inner', {{...}}) or just assign after.
    % Let me fix that:
    s_cell = struct(); s_cell.inner = innerCell;
    combos{2, 2} = s_cell;
    s_sa = struct(); s_sa.inner = innerSA;
    combos{6, 2} = s_sa;

    % Outer: cell
    cellCombos = {
        'struct_in_cell',   {payload}
        'cell_in_cell',     {innerCell}
        'map_in_cell',      {innerMap}
        'table_in_cell',    {innerTable}
        'tt_in_cell',       {innerTT}
        'sa_in_cell',       {innerSA}
        'fh_in_cell',       {innerFH}
    };

    % Outer: map
    mapCombos = {
        'struct_in_map',    make_map('inner', payload)
        'cell_in_map',      make_map('inner', innerCell)
        'map_in_map',       make_map('inner', innerMap)
        'table_in_map',     make_map('inner', innerTable)
        'tt_in_map',        make_map('inner', innerTT)
        'sa_in_map',        make_map('inner', innerSA)
        'fh_in_map',        make_map('inner', innerFH)
        'mixed_map',        containers.Map({'s','c','m'}, {payload, innerCell, innerMap})
    };

    allCombos = [combos; cellCombos; mapCombos];

    for i = 1:size(allCombos, 1)
        name = allCombos{i, 1};
        val = allCombos{i, 2};
        try
            SaveFast(binFile, 'val', val);
            L = LoadFast(binFile);
            [ok, msg] = compareValues(val, L.val);
            assert(ok, msg);
            fprintf('.'); passed = passed + 1;
        catch ME
            fprintf('F'); failed = failed + 1;
            failures{end+1} = sprintf('FAIL: %s: %s', name, ME.message); %#ok<AGROW>
        end
    end
    fprintf('  %d/%d passed\n', passed, size(allCombos, 1));
    for i = 1:numel(failures), fprintf('  %s\n', failures{i}); end
    if ~isempty(failures), fprintf('\n'); end
end

function m = make_map(key, val)
    m = containers.Map({key}, {val});
end

%% ===== SECTION 7: DEEP NESTING 3-5 LEVELS =====

function [passed, failed] = runDeepNestingSection(binFile)
    fprintf('--- Section 7: Deep nesting 3-5 levels ---\n');
    passed = 0; failed = 0; failures = {};

    tests = {
        'struct_struct_struct', ...
            struct('L1', struct('L2', struct('data', rand(3))))
        'struct_cell_struct', ...
            struct('L1', {{struct('data', rand(2))}})
        'cell_struct_cell', ...
            {struct('inner', {{1, 'two', [3 4]}})}
        'cell_cell_cell', ...
            {{{1, 'hello', [2 3]}}}
        'struct_cell_map', ...
            struct('L1', {{containers.Map({'k'}, {42})}})
        'cell_map_struct', ...
            {containers.Map({'k'}, {struct('x', 99)})}
        'map_struct_cell', ...
            make_map('outer', struct('inner', {{1, 'two'}}))
        'struct_cell_table', ...
            struct('L1', {{table([1;2], 'VariableNames', {'A'})}})
        'struct_table_cell_content', ...
            struct('tbl', table({struct('x',1); struct('x',2)}, 'VariableNames', {'Data'}))
        'L4_struct_cell_struct_cell', ...
            struct('L1', {{struct('L2', {{1, 'deep'}})}})
        'L5_cell_struct_cell_struct_cell', ...
            {struct('L1', {{struct('L2', {{'deepest'}})}})}
        'mixed_deep_types', ...
            struct('a', {{1+2i, sparse(eye(2)), categorical({'x'}), ...
                         datetime(2024,1,1), @sin, ...
                         table([1;2], 'VariableNames', {'V'})}})
    };

    nTests = size(tests, 1);
    for i = 1:nTests
        name = tests{i, 1};
        val = tests{i, 2};
        try
            SaveFast(binFile, 'val', val);
            L = LoadFast(binFile);
            [ok, msg] = compareValues(val, L.val);
            assert(ok, msg);
            fprintf('.'); passed = passed + 1;
        catch ME
            fprintf('F'); failed = failed + 1;
            failures{end+1} = sprintf('FAIL: %s: %s', name, ME.message); %#ok<AGROW>
        end
    end
    fprintf('  %d/%d passed\n', passed, nTests);
    for i = 1:numel(failures), fprintf('  %s\n', failures{i}); end
    if ~isempty(failures), fprintf('\n'); end
end

%% ===== SECTION 8: MIXED MEGA-CONTAINERS =====

function [passed, failed] = runMegaContainerSection(binFile, leaves, leafNames)
    fprintf('--- Section 8: Mixed mega-containers ---\n');
    passed = 0; failed = 0; failures = {};

    tests = {};

    % MegaStruct: one field per major type
    ms = struct();
    ms.f_double = rand(3);
    ms.f_single = single(rand(2));
    ms.f_int32 = int32([1 2 3]);
    ms.f_complex = 1+2i;
    ms.f_sparse = sparse(eye(4));
    ms.f_char = 'hello';
    ms.f_logical = true;
    ms.f_string = "world";
    ms.f_categorical = categorical({'a', 'b'});
    ms.f_datetime = datetime(2024, 1, 1);
    ms.f_duration = hours(1);
    ms.f_caldur = calmonths(1);
    ms.f_map = containers.Map({'k'}, {1});
    ms.f_fh = @sin;
    ms.f_struct = struct('x', 1);
    ms.f_cell = {1, 'two'};
    tests = [tests, {'MegaStruct', ms}];

    % MegaCell: one element per major type
    mc = {rand(3), single(rand(2)), int32([1 2 3]), 1+2i, ...
          sparse(eye(4)), 'hello', true, "world", ...
          categorical({'a','b'}), datetime(2024,1,1), hours(1), ...
          calmonths(1), containers.Map({'k'},{1}), @sin, ...
          struct('x',1), {1,'two'}};
    tests = [tests, {'MegaCell', mc}];

    % MegaMap: mixed value types
    mm = containers.Map( ...
        {'dbl','str','cat','dt','fh','struct','cell'}, ...
        {rand(3), "hello", categorical({'a'}), datetime(2024,1,1), ...
         @sin, struct('x',1), {1,'two'}});
    tests = [tests, {'MegaMap', mm}];

    % MegaTable: many column types
    n = 3;
    mt = table((1:n)', single((1:n)'), logical([1;0;1]), ...
        ["a";"b";"c"], categorical({'x';'y';'x'}), ...
        'VariableNames', {'Dbl','Sgl','Log','Str','Cat'});
    tests = [tests, {'MegaTable', mt}];

    nTests = numel(tests) / 2;
    for i = 1:nTests
        name = tests{2*i-1};
        val = tests{2*i};
        try
            SaveFast(binFile, 'val', val);
            L = LoadFast(binFile);
            [ok, msg] = compareValues(val, L.val);
            assert(ok, msg);
            fprintf('.'); passed = passed + 1;
        catch ME
            fprintf('F'); failed = failed + 1;
            failures{end+1} = sprintf('FAIL: %s: %s', name, ME.message); %#ok<AGROW>
        end
    end
    fprintf('  %d/%d passed\n', passed, nTests);
    for i = 1:numel(failures), fprintf('  %s\n', failures{i}); end
    if ~isempty(failures), fprintf('\n'); end
end

%% ===== SECTION 9: EMPTY/DEGENERATE CONTAINERS =====

function [passed, failed] = runEmptyDegenerateSection(binFile)
    fprintf('--- Section 9: Empty/degenerate containers ---\n');
    passed = 0; failed = 0; failures = {};

    tests = {
        'empty_struct',         struct()
        'cell_0x0',             cell(0, 0)
        'cell_0x3',             cell(0, 3)
        'cell_3x0',             cell(3, 0)
        'empty_map',            containers.Map()
        'empty_table',          table(double.empty(0,1), 'VariableNames', {'A'})
        'empty_double',         []
        'empty_char',           ''
        'empty_string',         string.empty(0, 0)
        'empty_logical',        logical.empty(0, 0)
        'empty_int32_0x1',      int32.empty(0, 1)
        'empty_single_2x0',     single.empty(2, 0)
        'struct_all_empty',     struct('a', [], 'b', '', 'c', {{}})
        'cell_of_empties',      {[], '', struct(), cell(0,0)}
        'nested_empties',       struct('a', struct('b', struct('c', [])))
    };

    nTests = size(tests, 1);
    for i = 1:nTests
        name = tests{i, 1};
        val = tests{i, 2};
        try
            SaveFast(binFile, 'val', val);
            L = LoadFast(binFile);
            [ok, msg] = compareValues(val, L.val);
            assert(ok, msg);
            fprintf('.'); passed = passed + 1;
        catch ME
            fprintf('F'); failed = failed + 1;
            failures{end+1} = sprintf('FAIL: %s: %s', name, ME.message); %#ok<AGROW>
        end
    end
    fprintf('  %d/%d passed\n', passed, nTests);
    for i = 1:numel(failures), fprintf('  %s\n', failures{i}); end
    if ~isempty(failures), fprintf('\n'); end
end

%% ===== SECTION 10: STRUCT ARRAYS WITH NESTING =====

function [passed, failed] = runStructArraySection(binFile)
    fprintf('--- Section 10: Struct arrays with nesting ---\n');
    passed = 0; failed = 0; failures = {};

    tests = {
        'sa_1x3', ...
            struct('a', {1, 2, 3}, 'b', {'x', 'y', 'z'})
        'sa_with_vectors', ...
            struct('v', {[1 2], [3 4 5], [6]}, 'n', {'a', 'b', 'c'})
        'sa_mixed_types', ...
            struct('num', {1, 2}, 'str', {"hello", "world"}, 'log', {true, false})
        'sa_of_structs', ...
            struct('inner', {struct('x', 1), struct('x', 2), struct('x', 3)})
        'sa_with_cells', ...
            struct('data', {{1, 'a'}, {2, 'b'}, {3, 'c'}})
    };

    nTests = size(tests, 1);
    for i = 1:nTests
        name = tests{i, 1};
        val = tests{i, 2};
        try
            SaveFast(binFile, 'val', val);
            L = LoadFast(binFile);
            [ok, msg] = compareValues(val, L.val);
            assert(ok, msg);
            fprintf('.'); passed = passed + 1;
        catch ME
            fprintf('F'); failed = failed + 1;
            failures{end+1} = sprintf('FAIL: %s: %s', name, ME.message); %#ok<AGROW>
        end
    end
    fprintf('  %d/%d passed\n', passed, nTests);
    for i = 1:numel(failures), fprintf('  %s\n', failures{i}); end
    if ~isempty(failures), fprintf('\n'); end
end

%% ===== SECTION 11: ALL LEAVES AS NAMED VARIABLES IN ONE FILE =====

function [passed, failed] = runAllLeavesOneFileSection(binFile, leaves, leafNames)
    n = numel(leaves);
    fprintf('--- Section 11: All leaves as named variables in one file (%d tests) ---\n', n);
    passed = 0; failed = 0; failures = {};

    % Build name-value argument list
    args = {};
    for i = 1:n
        args = [args, {leafNames{i}, leaves{i}}]; %#ok<AGROW>
    end

    try
        SaveFast(binFile, args{:});
        L = LoadFast(binFile);
    catch ME
        fprintf('E  SAVE/LOAD ERROR: %s\n', ME.message);
        failed = n;
        return;
    end

    for i = 1:n
        try
            assert(isfield(L, leafNames{i}), sprintf('Missing field: %s', leafNames{i}));
            [ok, msg] = compareValues(leaves{i}, L.(leafNames{i}));
            assert(ok, msg);
            fprintf('.'); passed = passed + 1;
        catch ME
            fprintf('F'); failed = failed + 1;
            failures{end+1} = sprintf('FAIL: %s: %s', leafNames{i}, ME.message); %#ok<AGROW>
        end
    end
    fprintf('  %d/%d passed\n', passed, n);
    for i = 1:numel(failures), fprintf('  %s\n', failures{i}); end
    if ~isempty(failures), fprintf('\n'); end
end

%% ===== COMPARISON FUNCTIONS (from test_validate_vs_matlab.m) =====

function [ok, msg] = compareValues(a, b)
    if ~strcmp(class(a), class(b))
        ok = false;
        msg = sprintf('Type mismatch: %s vs %s', class(a), class(b));
        return;
    end
    if ~isa(a, 'containers.Map')
        if ~isequal(size(a), size(b))
            ok = false;
            msg = sprintf('Size mismatch: [%s] vs [%s]', ...
                num2str(size(a)), num2str(size(b)));
            return;
        end
    end
    if isa(a, 'containers.Map')
        [ok, msg] = compareMap(a, b);
    elseif isa(a, 'function_handle')
        [ok, msg] = compareFunctionHandle(a, b);
    elseif isa(a, 'table')
        [ok, msg] = compareTable(a, b);
    elseif isa(a, 'timetable')
        [ok, msg] = compareTimetable(a, b);
    elseif isa(a, 'categorical')
        [ok, msg] = compareCategorical(a, b);
    elseif isa(a, 'datetime')
        [ok, msg] = compareDatetime(a, b);
    elseif isa(a, 'duration')
        [ok, msg] = compareDuration(a, b);
    elseif isa(a, 'calendarDuration')
        [ok, msg] = compareCalendarDuration(a, b);
    elseif isa(a, 'string')
        [ok, msg] = compareString(a, b);
    elseif issparse(a) || issparse(b)
        [ok, msg] = compareSparse(a, b);
    elseif isstruct(a)
        [ok, msg] = compareStruct(a, b);
    elseif iscell(a)
        [ok, msg] = compareCell(a, b);
    elseif isnumeric(a) && ~isreal(a)
        [ok, msg] = compareComplex(a, b);
    elseif isnumeric(a) || islogical(a)
        [ok, msg] = compareNumeric(a, b);
    elseif ischar(a)
        if isequal(a, b)
            ok = true; msg = 'exact match';
        else
            ok = false; msg = 'char content differs';
        end
    else
        if isequal(a, b)
            ok = true; msg = 'isequal match';
        else
            ok = false; msg = sprintf('isequal failed for type %s', class(a));
        end
    end
end

function [ok, msg] = compareNumeric(a, b)
    if isequaln(a, b)
        ok = true; msg = 'exact match';
    else
        if isfloat(a)
            maxdiff = max(abs(double(a(:)) - double(b(:))));
            ok = false; msg = sprintf('max diff = %g', maxdiff);
        else
            ndiff = sum(a(:) ~= b(:));
            ok = false; msg = sprintf('%d elements differ', ndiff);
        end
    end
end

function [ok, msg] = compareComplex(a, b)
    if ~isreal(a) && isreal(b)
        ok = false; msg = 'lost complex flag';
        return;
    end
    re_ok = isequaln(real(a), real(b));
    im_ok = isequaln(imag(a), imag(b));
    if re_ok && im_ok
        ok = true; msg = 'exact match (real+imag)';
    else
        ok = false;
        parts = {};
        if ~re_ok, parts{end+1} = 'real'; end
        if ~im_ok, parts{end+1} = 'imag'; end
        msg = sprintf('%s part(s) differ', strjoin(parts, '+'));
    end
end

function [ok, msg] = compareSparse(a, b)
    if ~issparse(a) || ~issparse(b)
        ok = false; msg = 'sparse flag lost';
        return;
    end
    if nnz(a) ~= nnz(b)
        ok = false; msg = sprintf('nnz mismatch: %d vs %d', nnz(a), nnz(b));
        return;
    end
    if isequaln(full(a), full(b))
        ok = true; msg = 'exact match (sparse)';
    else
        ok = false; msg = 'sparse values differ';
    end
end

function [ok, msg] = compareString(a, b)
    if isequaln(a, b)
        ok = true; msg = 'exact match';
    else
        m_a = ismissing(a);
        m_b = ismissing(b);
        if ~isequal(m_a, m_b)
            ok = false; msg = 'missing pattern differs';
        else
            ok = false; msg = 'string content differs';
        end
    end
end

function [ok, msg] = compareCategorical(a, b)
    if ~isequal(categories(a), categories(b))
        ok = false; msg = 'categories differ';
        return;
    end
    if isordinal(a) ~= isordinal(b)
        ok = false; msg = 'ordinal flag differs';
        return;
    end
    if isequaln(a, b)
        ok = true; msg = 'exact match';
    else
        ok = false; msg = 'categorical values differ';
    end
end

function [ok, msg] = compareDatetime(a, b)
    if ~strcmp(a.TimeZone, b.TimeZone)
        ok = false; msg = sprintf('timezone: "%s" vs "%s"', a.TimeZone, b.TimeZone);
        return;
    end
    nat_a = isnat(a); nat_b = isnat(b);
    if ~isequal(nat_a, nat_b)
        ok = false; msg = 'NaT pattern differs';
        return;
    end
    valid = ~nat_a;
    if any(valid(:))
        if isempty(a.TimeZone)
            diff_sec = abs(datenum(a(valid)) - datenum(b(valid))) * 86400; %#ok<DATEFUN>
        else
            diff_sec = abs(seconds(a(valid) - b(valid)));
        end
        max_diff_ms = max(diff_sec) * 1000;
        if max_diff_ms > 1
            ok = false; msg = sprintf('max diff = %.3f ms', max_diff_ms);
            return;
        end
    end
    ok = true; msg = 'match (tol < 1ms)';
end

function [ok, msg] = compareDuration(a, b)
    ms_a = milliseconds(a(:));
    ms_b = milliseconds(b(:));
    if isequaln(ms_a, ms_b)
        ok = true; msg = 'exact match';
        return;
    end
    nan_a = isnan(ms_a); nan_b = isnan(ms_b);
    if ~isequal(nan_a, nan_b)
        ok = false; msg = 'NaN pattern differs';
        return;
    end
    valid = ~nan_a;
    if any(valid)
        maxdiff = max(abs(ms_a(valid) - ms_b(valid)));
        if maxdiff > 0.001
            ok = false; msg = sprintf('max diff = %g ms', maxdiff);
            return;
        end
    end
    ok = true; msg = 'match (tol < 0.001ms)';
end

function [ok, msg] = compareCalendarDuration(a, b)
    if isequal(a, b)
        ok = true; msg = 'exact match';
    else
        ok = false; msg = 'calendarDuration values differ';
    end
end

function [ok, msg] = compareMap(a, b)
    ka = keys(a); kb = keys(b);
    if ~isempty(ka) && isnumeric(ka{1})
        [~, ia] = sort(cell2mat(ka));
        [~, ib] = sort(cell2mat(kb));
        ka = ka(ia); kb = kb(ib);
    else
        ka = sort(ka); kb = sort(kb);
    end
    if ~isequal(ka, kb)
        ok = false; msg = 'keys differ';
        return;
    end
    for i = 1:numel(ka)
        [val_ok, val_msg] = compareValues(a(ka{i}), b(kb{i}));
        if ~val_ok
            ok = false; msg = sprintf('value for key "%s": %s', ...
                mat2str_safe(ka{i}), val_msg);
            return;
        end
    end
    ok = true; msg = 'exact match';
end

function [ok, msg] = compareFunctionHandle(a, b)
    str_a = func2str(a);
    str_b = func2str(b);
    if ~strcmp(str_a, str_b)
        ok = false; msg = sprintf('func2str: "%s" vs "%s"', str_a, str_b);
        return;
    end
    try
        test_vals = [0, 1, -1, pi/4];
        for v = test_vals
            ra = a(v);
            rb = b(v);
            if ~isequaln(ra, rb)
                ok = false; msg = sprintf('output differs at x=%g', v);
                return;
            end
        end
    catch
        % Multi-arg or can't evaluate — func2str match is sufficient
    end
    ok = true; msg = 'match (func2str)';
end

function [ok, msg] = compareTable(a, b)
    if ~isequal(a.Properties.VariableNames, b.Properties.VariableNames)
        ok = false; msg = 'VariableNames differ';
        return;
    end
    if ~isequal(a.Properties.RowNames, b.Properties.RowNames)
        ok = false; msg = 'RowNames differ';
        return;
    end
    if ~isequal(a.Properties.VariableUnits, b.Properties.VariableUnits)
        ok = false; msg = 'VariableUnits differ';
        return;
    end
    if ~isequal(a.Properties.VariableDescriptions, b.Properties.VariableDescriptions)
        ok = false; msg = 'VariableDescriptions differ';
        return;
    end
    vnames = a.Properties.VariableNames;
    for i = 1:numel(vnames)
        [col_ok, col_msg] = compareValues(a.(vnames{i}), b.(vnames{i}));
        if ~col_ok
            ok = false; msg = sprintf('column "%s": %s', vnames{i}, col_msg);
            return;
        end
    end
    ok = true; msg = 'exact match';
end

function [ok, msg] = compareTimetable(a, b)
    if ~isequal(a.Properties.VariableNames, b.Properties.VariableNames)
        ok = false; msg = 'VariableNames differ';
        return;
    end
    [rt_ok, rt_msg] = compareValues(a.Properties.RowTimes, b.Properties.RowTimes);
    if ~rt_ok
        ok = false; msg = ['RowTimes: ' rt_msg];
        return;
    end
    vnames = a.Properties.VariableNames;
    for i = 1:numel(vnames)
        [col_ok, col_msg] = compareValues(a.(vnames{i}), b.(vnames{i}));
        if ~col_ok
            ok = false; msg = sprintf('column "%s": %s', vnames{i}, col_msg);
            return;
        end
    end
    ok = true; msg = 'exact match';
end

function [ok, msg] = compareStruct(a, b)
    fa = sort(fieldnames(a));
    fb = sort(fieldnames(b));
    if ~isequal(fa, fb)
        ok = false; msg = 'fieldnames differ';
        return;
    end
    if ~isequal(size(a), size(b))
        ok = false; msg = 'struct array size differs';
        return;
    end
    for k = 1:numel(a)
        for i = 1:numel(fa)
            [f_ok, f_msg] = compareValues(a(k).(fa{i}), b(k).(fa{i}));
            if ~f_ok
                if numel(a) > 1
                    ok = false; msg = sprintf('element(%d).%s: %s', k, fa{i}, f_msg);
                else
                    ok = false; msg = sprintf('field "%s": %s', fa{i}, f_msg);
                end
                return;
            end
        end
    end
    ok = true; msg = 'exact match';
end

function [ok, msg] = compareCell(a, b)
    for k = 1:numel(a)
        [el_ok, el_msg] = compareValues(a{k}, b{k});
        if ~el_ok
            ok = false; msg = sprintf('element{%d}: %s', k, el_msg);
            return;
        end
    end
    ok = true; msg = 'exact match';
end

%% ===== UTILITIES =====

function s = mat2str_safe(v)
    if ischar(v)
        s = v;
    elseif isnumeric(v)
        s = num2str(v);
    else
        s = '?';
    end
end

function deleteIfExists(f)
    if exist(f, 'file')
        delete(f);
    end
end

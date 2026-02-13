function results = test_validate_vs_matlab()
% TEST_VALIDATE_VS_MATLAB - Validate SaveFast/LoadFast against MATLAB save/load
%
% Systematically saves every supported type with both MATLAB's native
% save/load and SaveFast/LoadFast, then compares the results to catch
% any data corruption, precision loss, or metadata loss.
%
% Usage:
%   results = test_validate_vs_matlab()  % Run all tests, return results
%   test_validate_vs_matlab()            % Run all tests, display summary
%
% Output:
%   results - Struct with per-field match status and diagnostics

    fprintf('\n');
    fprintf('=================================================================\n');
    fprintf('  MATLAB vs SaveFast/LoadFast VALIDATION\n');
    fprintf('=================================================================\n');
    fprintf('Start time: %s\n', datestr(now));
    fprintf('=================================================================\n\n');

    % Build master test data
    data = buildTestData();
    fields = fieldnames(data);

    matFile = [tempname '.mat'];
    binFile = [tempname '.bin'];
    cleanup = onCleanup(@() deleteFiles(matFile, binFile)); %#ok<NASGU>

    % Results tracking
    n = numel(fields);
    mat_ok = zeros(n, 1);
    fast_ok = zeros(n, 1);
    mat_msg = cell(n, 1);
    fast_msg = cell(n, 1);

    fprintf('%-30s %-20s %-20s\n', 'Field', 'MATLAB save/load', 'SaveFast/LoadFast');
    fprintf('%s\n', repmat('-', 1, 70));

    for i = 1:n
        fn = fields{i};
        original = data.(fn);

        % --- MATLAB save/load ---
        try
            val = original; %#ok<NASGU>
            save(matFile, 'val', '-v7.3');
            tmp = load(matFile, 'val');
            loaded_mat = tmp.val;
            [mat_ok(i), mat_msg{i}] = compareValues(original, loaded_mat);
        catch ME
            mat_ok(i) = 0;
            mat_msg{i} = ['ERROR: ' ME.message];
        end

        % --- SaveFast/LoadFast ---
        try
            SaveFast(binFile, 'val', original);
            tmp2 = LoadFast(binFile);
            loaded_fast = tmp2.val;
            [fast_ok(i), fast_msg{i}] = compareValues(original, loaded_fast);
        catch ME
            fast_ok(i) = 0;
            fast_msg{i} = ['ERROR: ' ME.message];
        end

        % Print row
        mat_status = statusStr(mat_ok(i), mat_msg{i});
        fast_status = statusStr(fast_ok(i), fast_msg{i});
        fprintf('%-30s %-20s %-20s\n', fn, mat_status, fast_status);
    end

    % Summary
    fprintf('%s\n', repmat('-', 1, 70));
    fprintf('Summary: %d/%d MATLAB matches, %d/%d SaveFast matches\n', ...
        sum(mat_ok), n, sum(fast_ok), n);
    fprintf('=================================================================\n\n');

    % Print details for mismatches
    any_mismatch = false;
    for i = 1:n
        if ~mat_ok(i) || ~fast_ok(i)
            if ~any_mismatch
                fprintf('MISMATCH DETAILS:\n');
                any_mismatch = true;
            end
            fn = fields{i};
            if ~mat_ok(i)
                fprintf('  [%s] MATLAB: %s\n', fn, mat_msg{i});
            end
            if ~fast_ok(i)
                fprintf('  [%s] SaveFast: %s\n', fn, fast_msg{i});
            end
        end
    end
    if any_mismatch
        fprintf('\n');
    end

    % Build results struct
    results.fields = fields;
    results.mat_ok = mat_ok;
    results.fast_ok = fast_ok;
    results.mat_msg = mat_msg;
    results.fast_msg = fast_msg;
    results.mat_pass = sum(mat_ok);
    results.fast_pass = sum(fast_ok);
    results.total = n;

    if nargout == 0
        clear results;
    end
end

%% ===== Test Data Builder =====

function data = buildTestData()
    % Build a struct with one field per type/scenario

    % Use a fixed rng seed for reproducibility
    rng(42);

    % --- A. Numeric scalars & arrays ---
    data.dbl_scalar = pi;
    data.dbl_matrix = rand(50, 50);
    data.dbl_3d = rand(4, 5, 3);
    data.sgl_matrix = single(rand(10));
    data.i8 = int8(randi([-128 127], 1, 5));
    data.u8 = uint8(randi([0 255], 1, 5));
    data.i16 = int16(randi([-1000 1000], 1, 5));
    data.u16 = uint16(randi([0 65535], 1, 5));
    data.i32 = int32(randi([-1e5 1e5], 1, 5));
    data.u64 = uint64(randi([0 1e9], 1, 5));

    % --- B. Special numeric values ---
    data.special_vals = [NaN, Inf, -Inf, 0, -0, eps, realmin, realmax];
    data.special_single = single([NaN, Inf, -Inf, 0]);

    % --- C. Complex ---
    data.cplx_double = randn(5) + 1i*randn(5);
    data.cplx_single = single(randn(3) + 1i*randn(3));
    data.cplx_nan = [1+2i, NaN+3i, Inf-1i];

    % --- D. Sparse ---
    data.sp_double = sprandn(100, 100, 0.05);
    data.sp_logical = sparse(logical(eye(20)));
    data.sp_complex = sparse([1 3], [2 4], [1+1i 2-2i], 5, 5);
    data.sp_empty = sparse(10, 10);

    % --- E. Char & Logical ---
    data.char_vec = ['Hello, World! Special: tab' char(9) 'here'];
    data.char_matrix = ['abc'; 'def'; 'ghi'];
    data.logical_mat = logical(randi([0 1], 5, 5));

    % --- F. String ---
    data.str_array = ["alpha", "beta"; "gamma", "delta"];
    data.str_missing = ["hello", string(missing), ""];
    data.str_scalar = "single string";

    % --- G. Categorical ---
    data.cat_plain = categorical({'a', 'b', 'c', 'a', 'b'});
    data.cat_ordinal = categorical({'low', 'med', 'high'}, ...
        {'low', 'med', 'high'}, 'Ordinal', true);
    data.cat_undefined = categorical({'x', '', 'y'}, {'x', 'y'});

    % --- H. Datetime ---
    data.dt_no_tz = datetime(2024, 1, 1:5);
    data.dt_with_tz = datetime(2024, 6, 15, 12, 0, 0, 'TimeZone', 'America/New_York');
    data.dt_nat = [datetime(2024, 1, 1), NaT, datetime(2024, 12, 31)];

    % --- I. Duration & CalendarDuration ---
    data.dur_vec = hours(1:5) + minutes(30);
    data.dur_nan = [hours(1), duration(NaN, 0, 0), minutes(-30)];
    data.caldur_vec = calmonths(1:3) + caldays(15) + hours(6);

    % --- J. containers.Map ---
    data.map_char_keys = containers.Map({'a', 'b', 'c'}, {1, 'two', [3 4 5]});
    data.map_num_keys = containers.Map({1, 2, 3}, {'x', 'y', 'z'});

    % --- K. function_handle ---
    data.fh_builtin = @sin;
    data.fh_anon = @(x, y) x.^2 + y.^2;

    % --- L. Struct & Cell ---
    data.struct_simple = struct('a', 1, 'b', 'text', 'c', rand(3));
    data.struct_array = struct('x', {1, 2, 3}, 'y', {'a', 'b', 'c'});
    data.cell_mixed = {1, 'two', [3; 4; 5], {6, 7}, struct('a', 8)};
    data.cell_2d = reshape({1, 2, 3, 4, 5, 6}, 2, 3);

    % --- M. Table & Timetable ---
    data.tbl_mixed = buildMixedTable();
    data.tbl_empty = table(double.empty(0, 1), cell(0, 1), ...
        'VariableNames', {'Num', 'Txt'});
    tt_times = datetime(2024, 1, 1) + hours(0:4)';
    data.tt_basic = timetable(tt_times, (1:5)', randn(5, 1), ...
        'VariableNames', {'A', 'B'});

    % --- N. Deep nesting ---
    data.nested_struct = struct('level1', ...
        struct('level2', ...
            struct('data', rand(3), ...
                   'cell', {{1+2i, 'text', int32(42)}})));
    data.nested_cell = {randn(2)+1i*randn(2), sparse(eye(3)), ...
        "hello", categorical({'a'}), datetime(2024,1,1), @sin, ...
        table([1;2], {'a';'b'}, 'VariableNames', {'N','S'})};
    data.nested_table_in_struct = struct( ...
        'info', 'metadata', ...
        'tbl', table(["x";"y";"z"], categorical({'a';'b';'c'}), ...
            'VariableNames', {'Str', 'Cat'}));
end

function tbl = buildMixedTable()
    Num = (1:5)';
    Str = ["a"; "b"; "c"; "d"; "e"];
    Cat = categorical({'x'; 'y'; 'x'; 'y'; 'x'});
    Log = logical([1; 0; 1; 0; 1]);
    tbl = table(Num, Str, Cat, Log, ...
        'RowNames', {'r1', 'r2', 'r3', 'r4', 'r5'});
    tbl.Properties.VariableUnits = {'m', '', '', ''};
    tbl.Properties.VariableDescriptions = {'count', 'label', 'group', 'flag'};
end

%% ===== Comparison Logic =====

function [ok, msg] = compareValues(a, b)
    % Compare two values with type-appropriate logic
    % Returns [match_bool, description_string]

    % Type check
    if ~strcmp(class(a), class(b))
        ok = false;
        msg = sprintf('Type mismatch: %s vs %s', class(a), class(b));
        return;
    end

    % Size check (except containers.Map which has no meaningful size)
    if ~isa(a, 'containers.Map')
        if ~isequal(size(a), size(b))
            ok = false;
            msg = sprintf('Size mismatch: [%s] vs [%s]', ...
                num2str(size(a)), num2str(size(b)));
            return;
        end
    end

    % Type-specific comparison
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
        % Fallback
        if isequal(a, b)
            ok = true; msg = 'isequal match';
        else
            ok = false; msg = sprintf('isequal failed for type %s', class(a));
        end
    end
end

% --- Numeric ---
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

% --- Complex ---
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

% --- Sparse ---
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

% --- String ---
function [ok, msg] = compareString(a, b)
    % Use isequaln: missing (like NaN) should compare as equal
    if isequaln(a, b)
        ok = true; msg = 'exact match';
    else
        % Check missing pattern
        m_a = ismissing(a);
        m_b = ismissing(b);
        if ~isequal(m_a, m_b)
            ok = false; msg = 'missing pattern differs';
        else
            ok = false; msg = 'string content differs';
        end
    end
end

% --- Categorical ---
function [ok, msg] = compareCategorical(a, b)
    if ~isequal(categories(a), categories(b))
        ok = false; msg = 'categories differ';
        return;
    end
    if isordinal(a) ~= isordinal(b)
        ok = false; msg = 'ordinal flag differs';
        return;
    end
    % Use isequaln: <undefined> (like NaN) should compare as equal
    if isequaln(a, b)
        ok = true; msg = 'exact match';
    else
        ok = false; msg = 'categorical values differ';
    end
end

% --- Datetime ---
function [ok, msg] = compareDatetime(a, b)
    % Check timezone
    if ~strcmp(a.TimeZone, b.TimeZone)
        ok = false; msg = sprintf('timezone: "%s" vs "%s"', a.TimeZone, b.TimeZone);
        return;
    end
    % Check NaT pattern
    nat_a = isnat(a);
    nat_b = isnat(b);
    if ~isequal(nat_a, nat_b)
        ok = false; msg = 'NaT pattern differs';
        return;
    end
    % Compare valid elements
    valid = ~nat_a;
    if any(valid(:))
        % Need timezone for datenum comparison
        if isempty(a.TimeZone)
            diff_sec = abs(datenum(a(valid)) - datenum(b(valid))) * 86400;
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

% --- Duration ---
function [ok, msg] = compareDuration(a, b)
    ms_a = milliseconds(a(:));
    ms_b = milliseconds(b(:));
    if isequaln(ms_a, ms_b)
        ok = true; msg = 'exact match';
        return;
    end
    % NaN-aware comparison
    nan_a = isnan(ms_a);
    nan_b = isnan(ms_b);
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

% --- CalendarDuration ---
function [ok, msg] = compareCalendarDuration(a, b)
    if isequal(a, b)
        ok = true; msg = 'exact match';
    else
        ok = false; msg = 'calendarDuration values differ';
    end
end

% --- containers.Map ---
function [ok, msg] = compareMap(a, b)
    ka = keys(a);
    kb = keys(b);
    % Sort keys — cell of char uses sort(), cell of numeric needs sortrows
    if ~isempty(ka) && isnumeric(ka{1})
        [~, ia] = sort(cell2mat(ka));
        [~, ib] = sort(cell2mat(kb));
        ka = ka(ia);
        kb = kb(ib);
    else
        ka = sort(ka);
        kb = sort(kb);
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

% --- function_handle ---
function [ok, msg] = compareFunctionHandle(a, b)
    str_a = func2str(a);
    str_b = func2str(b);
    if ~strcmp(str_a, str_b)
        ok = false; msg = sprintf('func2str: "%s" vs "%s"', str_a, str_b);
        return;
    end
    % Test evaluation at sample points
    try
        % Try single-arg evaluation
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

% --- Table ---
function [ok, msg] = compareTable(a, b)
    % Variable names
    if ~isequal(a.Properties.VariableNames, b.Properties.VariableNames)
        ok = false; msg = 'VariableNames differ';
        return;
    end
    % Row names
    if ~isequal(a.Properties.RowNames, b.Properties.RowNames)
        ok = false; msg = 'RowNames differ';
        return;
    end
    % Units
    if ~isequal(a.Properties.VariableUnits, b.Properties.VariableUnits)
        ok = false; msg = 'VariableUnits differ';
        return;
    end
    % Descriptions
    if ~isequal(a.Properties.VariableDescriptions, b.Properties.VariableDescriptions)
        ok = false; msg = 'VariableDescriptions differ';
        return;
    end
    % Compare each column
    vnames = a.Properties.VariableNames;
    for i = 1:numel(vnames)
        col_a = a.(vnames{i});
        col_b = b.(vnames{i});
        [col_ok, col_msg] = compareValues(col_a, col_b);
        if ~col_ok
            ok = false; msg = sprintf('column "%s": %s', vnames{i}, col_msg);
            return;
        end
    end
    ok = true; msg = 'exact match';
end

% --- Timetable ---
function [ok, msg] = compareTimetable(a, b)
    % Variable names
    if ~isequal(a.Properties.VariableNames, b.Properties.VariableNames)
        ok = false; msg = 'VariableNames differ';
        return;
    end
    % Row times
    [rt_ok, rt_msg] = compareValues(a.Properties.RowTimes, b.Properties.RowTimes);
    if ~rt_ok
        ok = false; msg = ['RowTimes: ' rt_msg];
        return;
    end
    % Compare each column
    vnames = a.Properties.VariableNames;
    for i = 1:numel(vnames)
        col_a = a.(vnames{i});
        col_b = b.(vnames{i});
        [col_ok, col_msg] = compareValues(col_a, col_b);
        if ~col_ok
            ok = false; msg = sprintf('column "%s": %s', vnames{i}, col_msg);
            return;
        end
    end
    ok = true; msg = 'exact match';
end

% --- Struct ---
function [ok, msg] = compareStruct(a, b)
    fa = sort(fieldnames(a));
    fb = sort(fieldnames(b));
    if ~isequal(fa, fb)
        ok = false; msg = 'fieldnames differ';
        return;
    end
    % Check struct array size
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

% --- Cell ---
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

%% ===== Utilities =====

function s = statusStr(ok, msg)
    if ok
        s = ['MATCH (' msg ')'];
        % Truncate if too long
        if length(s) > 25
            s = 'MATCH';
        end
    else
        s = 'MISMATCH';
    end
end

function s = mat2str_safe(v)
    if ischar(v)
        s = v;
    elseif isnumeric(v)
        s = num2str(v);
    else
        s = '?';
    end
end

function deleteFiles(varargin)
    for i = 1:nargin
        if exist(varargin{i}, 'file')
            delete(varargin{i});
        end
    end
end

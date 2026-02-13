function SaveFast(filepath, varargin)
% SAVEFAST - Fast binary save with multi-variable support
%
% Syntax:
%   SaveFast(filepath, var1)                    % Single variable (backward compatible)
%   SaveFast(filepath, 'name1', var1, ...)      % Named variables
%   SaveFast(filepath, struct('a', 1, 'b', 2))  % Save struct fields as variables
%
% Examples:
%   % Single variable (backward compatible)
%   data = rand(1000);
%   SaveFast('test.bin', data);
%
%   % Multiple named variables
%   A = magic(5);
%   B = 'hello';
%   SaveFast('test.bin', 'A', A, 'B', B);
%
%   % Save struct fields
%   results = struct('data', rand(100), 'params', struct('alpha', 0.1));
%   SaveFast('test.bin', results);
%
% Supported types:
%   - All numeric types (double, single, int8-64, uint8-64)
%   - Complex numeric arrays (bit flag 0x40)
%   - Sparse arrays (bit flag 0x80, combinable with 0x40)
%   - char, logical, struct, cell
%   - string, table, timetable, categorical
%   - datetime, duration, calendarDuration
%   - containers.Map, function_handle
%   - gpuArray and tall arrays (transparently gathered)
%
% File Format Version 3:
%   Header:
%     - Magic bytes: [77 70] ('MF' for MatFast)
%     - Version: 3
%     - Number of variables
%     - Timestamp
%   Type ID encoding:
%     - Bits 0-5: base type ID (0-63)
%     - Bit 6 (64):  complex flag (numeric types 1-10 only)
%     - Bit 7 (128): sparse flag  (numeric types 1-10 + logical 12)
%
% See also: LoadFast, SaveMatFast, LoadMatFast

% Version: 3.0
% Date: 2026-02-12

    %% Input Validation
    if nargin < 2
        error('SaveFast:NotEnoughInputs', ...
            'Usage: SaveFast(filepath, data) or SaveFast(filepath, ''name1'', var1, ...)');
    end

    if ~ischar(filepath) && ~isstring(filepath)
        error('SaveFast:InvalidFilepath', 'Filepath must be a string');
    end

    filepath = char(filepath);

    %% Parse inputs into variable names and values
    [varNames, varValues] = parseInputs(varargin{:});
    numVars = length(varNames);

    if numVars == 0
        error('SaveFast:NoVariables', 'No variables to save');
    end

    %% Check file writability
    [fdir, ~, ~] = fileparts(filepath);
    if ~isempty(fdir) && ~exist(fdir, 'dir')
        error('SaveFast:DirectoryNotFound', ...
            'Directory does not exist: %s', fdir);
    end

    % Check if file exists and is writable
    if exist(filepath, 'file')
        ftest = fopen(filepath, 'a');
        if ftest == -1
            error('SaveFast:PermissionDenied', ...
                'Cannot write to file (check permissions): %s', filepath);
        end
        fclose(ftest);
    end

    %% Open file for writing
    [fp, errormsg] = fopen(filepath, 'w');
    if ~isempty(errormsg) || fp == -1
        error('SaveFast:FileOpenError', ...
            'Failed to open file for writing: %s\nError: %s', filepath, errormsg);
    end

    % Ensure file is closed on error
    cleanupObj = onCleanup(@() fclose(fp));

    try
        %% Write File Header
        % Magic bytes: 'MF' (MatFast)
        fwrite(fp, uint8([77 70]), 'uint8');

        % Format version
        fwrite(fp, uint16(3), 'uint16');

        % Number of variables
        fwrite(fp, uint32(numVars), 'uint32');

        % Timestamp
        timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
        fwrite(fp, uint8(length(timestamp)), 'uint8');
        fwrite(fp, timestamp, 'char*1');

        %% Write each variable
        for i = 1:numVars
            varName = varNames{i};
            varData = varValues{i};

            % Write variable name
            if length(varName) > 255
                error('SaveFast:NameTooLong', ...
                    'Variable name exceeds 255 chars: %s', varName);
            end

            fwrite(fp, uint8(length(varName)), 'uint8');
            fwrite(fp, varName, 'char*1');

            % Write variable data (recursive)
            writeElement(fp, varData);
        end

    catch ME
        % File will be closed by cleanupObj
        rethrow(ME);
    end

    % File closed automatically by cleanupObj
end

%% Helper Functions

function [names, values] = parseInputs(varargin)
    % Parse variable inputs into names and values

    names = {};
    values = {};

    if nargin == 1
        % Single variable case
        data = varargin{1};

        if isstruct(data) && numel(data) == 1
            % Struct: save each field as a variable
            fnames = fieldnames(data);
            names = fnames';
            values = cell(1, length(fnames));
            for i = 1:length(fnames)
                values{i} = data.(fnames{i});
            end
        else
            % Single unnamed variable (backward compatible)
            names = {'data'};
            values = {data};
        end

    elseif mod(nargin, 2) == 0
        % Paired name-value arguments
        for i = 1:2:nargin
            if ~ischar(varargin{i}) && ~isstring(varargin{i})
                error('SaveFast:InvalidVariableName', ...
                    'Variable names must be strings. Argument %d is %s', ...
                    i, class(varargin{i}));
            end

            varName = char(varargin{i});

            % Validate variable name
            if ~isvarname(varName)
                error('SaveFast:InvalidVariableName', ...
                    'Invalid MATLAB variable name: %s', varName);
            end

            names{end+1} = varName; %#ok<AGROW>
            values{end+1} = varargin{i+1}; %#ok<AGROW>
        end

    else
        error('SaveFast:InvalidArguments', ...
            'Arguments must be pairs of (name, value) or a single struct/variable');
    end
end

function type_id = getBaseTypeId(tp)
    % Return base type ID for a class name, or -1 if not a base type
    switch tp
        case 'double',  type_id = 1;
        case 'single',  type_id = 2;
        case 'int8',    type_id = 3;
        case 'uint8',   type_id = 4;
        case 'int16',   type_id = 5;
        case 'uint16',  type_id = 6;
        case 'int32',   type_id = 7;
        case 'uint32',  type_id = 8;
        case 'int64',   type_id = 9;
        case 'uint64',  type_id = 10;
        case 'char',    type_id = 11;
        case 'logical', type_id = 12;
        case 'struct',  type_id = 20;
        case 'cell',    type_id = 21;
        otherwise,       type_id = -1;
    end
end

function writeElement(fp, data)
    % Recursive element writer with type dispatch

    % Transparent conversions: gpuArray and tall
    if isa(data, 'gpuArray')
        data = gather(data);
    end
    if isa(data, 'tall')
        data = gather(data);
    end

    % Dispatch by type — check special types first
    if isa(data, 'string')
        writeString(fp, data);
        return;
    elseif isa(data, 'table')
        writeTable(fp, data);
        return;
    elseif isa(data, 'timetable')
        writeTimetable(fp, data);
        return;
    elseif isa(data, 'categorical')
        writeCategorical(fp, data);
        return;
    elseif isa(data, 'datetime')
        writeDatetime(fp, data);
        return;
    elseif isa(data, 'duration')
        writeDuration(fp, data);
        return;
    elseif isa(data, 'calendarDuration')
        writeCalendarDuration(fp, data);
        return;
    elseif isa(data, 'containers.Map')
        writeContainersMap(fp, data);
        return;
    elseif isa(data, 'function_handle')
        writeFunctionHandle(fp, data);
        return;
    end

    % Empty array
    if isempty(data)
        writeEmptyElement(fp, class(data), size(data));
        return;
    end

    tp = class(data);
    base_id = getBaseTypeId(tp);

    if base_id == -1
        error('SaveFast:UnsupportedType', ...
            'Unsupported type: %s. Cannot serialize arbitrary class objects.', tp);
    end

    % Handle sparse arrays
    if issparse(data)
        writeSparse(fp, data, base_id);
        return;
    end

    % Handle complex numeric arrays
    if isnumeric(data) && ~isreal(data)
        writeComplex(fp, data, base_id);
        return;
    end

    % Standard types
    dims = size(data);
    ndims_val = numel(dims);

    % Write type ID
    fwrite(fp, uint8(base_id), 'uint8');

    % Write dimensions
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');

    % Write data based on type
    if base_id < 20
        % Simple types
        if base_id == 12 % logical
            fwrite(fp, uint8(data(:)), 'uint8');
        else
            fwrite(fp, data(:), tp);
        end

    elseif base_id == 20
        % Struct
        writeStruct(fp, data);

    elseif base_id == 21
        % Cell array
        writeCell(fp, data);
    end
end

function writeEmptyElement(fp, className, dims)
    % Write empty element marker with dimensions preserved
    fwrite(fp, uint8(0), 'uint8'); % Type ID 0 = empty
    fwrite(fp, uint8(length(className)), 'uint8');
    fwrite(fp, className, 'char*1');
    ndims_val = numel(dims);
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');
end

function writeStruct(fp, data)
    % Write struct array

    field_names = fieldnames(data);
    num_fields = numel(field_names);

    % Write number of fields
    fwrite(fp, uint32(num_fields), 'uint32');

    % Write field names
    for i = 1:num_fields
        fname = field_names{i};
        fwrite(fp, uint8(length(fname)), 'uint8');
        fwrite(fp, fname, 'char*1');
    end

    % Write data (iterate through struct array elements)
    num_elements = numel(data);
    for k = 1:num_elements
        for f = 1:num_fields
            writeElement(fp, data(k).(field_names{f}));
        end
    end
end

function writeCell(fp, data)
    % Write cell array

    num_elements = numel(data);
    for k = 1:num_elements
        writeElement(fp, data{k});
    end
end

function writeComplex(fp, data, base_id)
    % Write complex numeric array: [type_id|0x40][ndims][dims][real][imag]
    type_id = bitor(base_id, 64); % Set bit 6
    dims = size(data);
    ndims_val = numel(dims);
    tp = class(data);

    fwrite(fp, uint8(type_id), 'uint8');
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');

    % Write real and imaginary parts separately
    fwrite(fp, real(data(:)), tp);
    fwrite(fp, imag(data(:)), tp);
end

function writeSparse(fp, data, base_id)
    % Write sparse array: [type_id|0x80][ndims][dims][nnz][rows][cols][values]
    % If also complex, type_id|0x80|0x40
    is_complex = ~isreal(data);
    type_id = bitor(base_id, 128); % Set bit 7
    if is_complex
        type_id = bitor(type_id, 64); % Also set bit 6
    end

    dims = size(data);
    ndims_val = numel(dims);

    fwrite(fp, uint8(type_id), 'uint8');
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');

    % Extract triplets
    [rows, cols, vals] = find(data);
    nz = numel(vals);
    tp = class(data);

    fwrite(fp, uint64(nz), 'uint64');
    fwrite(fp, uint64(rows(:)), 'uint64');
    fwrite(fp, uint64(cols(:)), 'uint64');

    if is_complex
        fwrite(fp, real(vals(:)), tp);
        fwrite(fp, imag(vals(:)), tp);
    else
        if base_id == 12 % sparse logical
            fwrite(fp, uint8(vals(:)), 'uint8');
        else
            fwrite(fp, vals(:), tp);
        end
    end
end

function writeString(fp, data)
    % Write string array: [30][ndims][dims] then per element [is_missing][len][chars]
    fwrite(fp, uint8(30), 'uint8');
    dims = size(data);
    ndims_val = numel(dims);
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');

    for k = 1:numel(data)
        if ismissing(data(k))
            fwrite(fp, uint8(1), 'uint8'); % is_missing = true
        else
            fwrite(fp, uint8(0), 'uint8'); % is_missing = false
            chars = char(data(k));
            fwrite(fp, uint32(length(chars)), 'uint32');
            if ~isempty(chars)
                fwrite(fp, chars, 'char*1');
            end
        end
    end
end

function writeTable(fp, data)
    % Write table: [31][nrows][nvars][varnames][has_rownames][rownames]
    %              [has_units][units][has_descr][descr] then each column
    fwrite(fp, uint8(31), 'uint8');
    fwrite(fp, uint64(height(data)), 'uint64');

    varNames = data.Properties.VariableNames;
    nvars = numel(varNames);
    fwrite(fp, uint32(nvars), 'uint32');

    % Write variable names as cell via writeElement
    writeElement(fp, varNames);

    % Row names
    rowNames = data.Properties.RowNames;
    if ~isempty(rowNames)
        fwrite(fp, uint8(1), 'uint8');
        writeElement(fp, rowNames);
    else
        fwrite(fp, uint8(0), 'uint8');
    end

    % Variable units
    units = data.Properties.VariableUnits;
    if ~isempty(units)
        fwrite(fp, uint8(1), 'uint8');
        writeElement(fp, units);
    else
        fwrite(fp, uint8(0), 'uint8');
    end

    % Variable descriptions
    descr = data.Properties.VariableDescriptions;
    if ~isempty(descr)
        fwrite(fp, uint8(1), 'uint8');
        writeElement(fp, descr);
    else
        fwrite(fp, uint8(0), 'uint8');
    end

    % Write each column
    for i = 1:nvars
        writeElement(fp, data{:, i});
    end
end

function writeTimetable(fp, data)
    % Write timetable: [32][nrows][nvars][varnames][rowtimes] then like table
    fwrite(fp, uint8(32), 'uint8');
    fwrite(fp, uint64(height(data)), 'uint64');

    varNames = data.Properties.VariableNames;
    nvars = numel(varNames);
    fwrite(fp, uint32(nvars), 'uint32');

    % Write variable names
    writeElement(fp, varNames);

    % Write row times
    writeElement(fp, data.Properties.RowTimes);

    % Variable units
    units = data.Properties.VariableUnits;
    if ~isempty(units)
        fwrite(fp, uint8(1), 'uint8');
        writeElement(fp, units);
    else
        fwrite(fp, uint8(0), 'uint8');
    end

    % Variable descriptions
    descr = data.Properties.VariableDescriptions;
    if ~isempty(descr)
        fwrite(fp, uint8(1), 'uint8');
        writeElement(fp, descr);
    else
        fwrite(fp, uint8(0), 'uint8');
    end

    % Write each column
    for i = 1:nvars
        writeElement(fp, data{:, i});
    end
end

function writeCategorical(fp, data)
    % Write categorical: [33][ndims][dims][is_ordinal][is_protected][num_cats][cat names][codes]
    fwrite(fp, uint8(33), 'uint8');
    dims = size(data);
    ndims_val = numel(dims);
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');

    fwrite(fp, uint8(isordinal(data)), 'uint8');
    fwrite(fp, uint8(isprotected(data)), 'uint8');

    cats = categories(data);
    num_cats = numel(cats);
    fwrite(fp, uint32(num_cats), 'uint32');

    % Write category names
    for i = 1:num_cats
        catName = cats{i};
        fwrite(fp, uint32(length(catName)), 'uint32');
        fwrite(fp, catName, 'char*1');
    end

    % Write codes (category indices, 0 = undefined)
    % Use double(data) trick: undefined -> NaN, category i -> i
    % But we need integer codes. Use ismember approach.
    codes = zeros(numel(data), 1, 'uint32');
    for i = 1:num_cats
        codes(data == cats{i}) = i;
    end
    % undefined elements remain 0
    fwrite(fp, codes, 'uint32');
end

function writeDatetime(fp, data)
    % Write datetime: [34][ndims][dims][has_tz][tz_str][format_len][format_str][posixtime doubles]
    fwrite(fp, uint8(34), 'uint8');
    dims = size(data);
    ndims_val = numel(dims);
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');

    % Timezone
    tz = data.TimeZone;
    if ~isempty(tz)
        fwrite(fp, uint8(1), 'uint8');
        fwrite(fp, uint8(length(tz)), 'uint8');
        fwrite(fp, tz, 'char*1');
    else
        fwrite(fp, uint8(0), 'uint8');
    end

    % Format
    fmt = data.Format;
    if isempty(fmt)
        fmt = '';
    end
    fwrite(fp, uint8(length(fmt)), 'uint8');
    if ~isempty(fmt)
        fwrite(fp, fmt, 'char*1');
    end

    % Data as posixtime (seconds since epoch as doubles)
    % For datetime without timezone, we need to set a temporary timezone
    if isempty(tz)
        % Convert via datenum to preserve values without timezone ambiguity
        fwrite(fp, datenum(data), 'double');
    else
        fwrite(fp, posixtime(data), 'double');
    end
end

function writeDuration(fp, data)
    % Write duration: [35][ndims][dims][format_len][format_str][milliseconds doubles]
    fwrite(fp, uint8(35), 'uint8');
    dims = size(data);
    ndims_val = numel(dims);
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');

    % Format
    fmt = data.Format;
    if isempty(fmt)
        fmt = '';
    end
    fwrite(fp, uint8(length(fmt)), 'uint8');
    if ~isempty(fmt)
        fwrite(fp, fmt, 'char*1');
    end

    % Data as milliseconds (double)
    fwrite(fp, milliseconds(data(:)), 'double');
end

function writeCalendarDuration(fp, data)
    % Write calendarDuration: [36][ndims][dims][months][days][millis]
    fwrite(fp, uint8(36), 'uint8');
    dims = size(data);
    ndims_val = numel(dims);
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');

    % Split into components
    [y, m, d, t] = split(data, {'years', 'months', 'days', 'time'});
    months_total = int32(y(:) * 12 + m(:));
    days_total = int32(d(:));
    millis = milliseconds(t(:));

    fwrite(fp, months_total, 'int32');
    fwrite(fp, days_total, 'int32');
    fwrite(fp, millis, 'double');
end

function writeContainersMap(fp, data)
    % Write containers.Map: [37][num_entries][key_type_str][value_type_str][keys cell][values cell]
    fwrite(fp, uint8(37), 'uint8');

    k = keys(data);
    v = values(data);
    num_entries = uint32(length(k));
    fwrite(fp, num_entries, 'uint32');

    % Key type
    kt = data.KeyType;
    fwrite(fp, uint8(length(kt)), 'uint8');
    fwrite(fp, kt, 'char*1');

    % Value type
    vt = data.ValueType;
    fwrite(fp, uint8(length(vt)), 'uint8');
    fwrite(fp, vt, 'char*1');

    % Keys and values as cell arrays
    writeElement(fp, k);
    writeElement(fp, v);
end

function writeFunctionHandle(fp, data)
    % Write function_handle: [38][str_len][func_str]
    fwrite(fp, uint8(38), 'uint8');

    funcStr = func2str(data);
    fwrite(fp, uint16(length(funcStr)), 'uint16');
    fwrite(fp, funcStr, 'char*1');
end

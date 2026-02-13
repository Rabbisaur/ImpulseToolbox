function varargout = LoadFast(filepath, varargin)
% LOADFAST - Fast binary load with multi-variable support
%
% Syntax:
%   data = LoadFast(filepath)                    % Load all variables as struct
%   [var1, var2] = LoadFast(filepath)           % Load all variables as separate outputs
%   data = LoadFast(filepath, 'var1', 'var2')   % Load specific variables
%   LoadFast(filepath)                          % Load into workspace (like load)
%
% Examples:
%   % Load all variables as struct
%   S = LoadFast('test.bin');
%   % Access: S.A, S.B, etc.
%
%   % Load specific variables
%   [A, B] = LoadFast('test.bin');
%
%   % Load into workspace
%   LoadFast('test.bin');  % Variables appear in caller workspace
%
%   % Load specific variables by name
%   data = LoadFast('test.bin', 'A', 'B');
%
% Supported types:
%   - All numeric types (double, single, int8-64, uint8-64)
%   - Complex numeric arrays, sparse arrays
%   - char, logical, struct, cell
%   - string, table, timetable, categorical
%   - datetime, duration, calendarDuration
%   - containers.Map, function_handle
%
% See also: SaveFast, SaveMatFast, LoadMatFast

% Version: 3.0
% Date: 2026-02-12

    %% Input Validation
    if ~ischar(filepath) && ~isstring(filepath)
        error('LoadFast:InvalidFilepath', 'Filepath must be a string');
    end

    filepath = char(filepath);

    if ~exist(filepath, 'file')
        error('LoadFast:FileNotFound', 'File does not exist: %s', filepath);
    end

    %% Open file
    [fp, errormsg] = fopen(filepath, 'r');
    if ~isempty(errormsg) || fp == -1
        error('LoadFast:FileOpenError', ...
            'Failed to open file: %s\nError: %s', filepath, errormsg);
    end

    % Ensure file is closed on error
    cleanupObj = onCleanup(@() fclose(fp));

    try
        %% Detect file format version
        magic = fread(fp, 2, 'uint8');

        if isequal(magic, [77; 70])  % 'MF'
            % Version 2 or 3 format
            [varNames, varValues] = readFormatV2V3(fp);

        elseif isequal(magic, [1; 1])
            % Version 1 format (original SaveMatFast)
            fseek(fp, 0, 'bof');  % Rewind
            [varNames, varValues] = readFormatV1(fp);

        else
            error('LoadFast:InvalidFormat', ...
                'Unrecognized file format. Magic bytes: [%d %d]', magic(1), magic(2));
        end

    catch ME
        rethrow(ME);
    end

    % File closed automatically by cleanupObj

    %% Filter requested variables
    if nargin > 1
        requestedVars = varargin;
        [varNames, varValues] = filterVariables(varNames, varValues, requestedVars);
    end

    %% Return data based on output arguments
    if nargout == 0
        % Load into caller workspace (like MATLAB load)
        for i = 1:length(varNames)
            assignin('caller', varNames{i}, varValues{i});
        end

    elseif nargout == 1
        % Return as struct
        S = struct();
        for i = 1:length(varNames)
            S.(varNames{i}) = varValues{i};
        end
        varargout{1} = S;

    else
        % Return as separate outputs
        if nargout > length(varNames)
            error('LoadFast:TooManyOutputs', ...
                'File contains %d variables but %d outputs requested', ...
                length(varNames), nargout);
        end

        varargout = varValues(1:nargout);
    end
end

%% Format Readers

function [varNames, varValues] = readFormatV2V3(fp)
    % Read Version 2/3 format (multi-variable)

    % Read version
    formatVersion = fread(fp, 1, 'uint16');
    if formatVersion ~= 2 && formatVersion ~= 3
        warning('LoadFast:VersionMismatch', ...
            'Expected version 2 or 3, got version %d', formatVersion);
    end

    % Read number of variables
    numVars = fread(fp, 1, 'uint32');

    % Read timestamp (informational only)
    timestamp_len = fread(fp, 1, 'uint8');
    timestamp = fread(fp, timestamp_len, '*char')'; %#ok<NASGU>

    % Read variables
    varNames = cell(1, numVars);
    varValues = cell(1, numVars);

    for i = 1:numVars
        % Read variable name
        name_len = fread(fp, 1, 'uint8');
        varName = fread(fp, name_len, '*char')';
        varNames{i} = varName;

        % Read variable data
        varValues{i} = readElement(fp);
    end
end

function [varNames, varValues] = readFormatV1(fp)
    % Read Version 1 format (single variable, original SaveMatFast)

    % Read header
    header.identificationBytes = fread(fp, 2, 'uint8');
    header.headersize = fread(fp, 1, 'uint8');
    header.dateBytes = fread(fp, 14, '*char');
    header.dataSize = fread(fp, 1, 'double');
    header.elementSize = fread(fp, 1, 'uint8');
    header.elementType = strtrim(fread(fp, 6, '*char')');
    header.dataNumberDimensions = fread(fp, 1, 'double');
    header.dataDimSize = fread(fp, header.dataNumberDimensions, 'double')';
    header.commentsize = fread(fp, 1, 'uint8');
    header.comment = fread(fp, header.commentsize, '*char');

    % Seek to end of header
    fseek(fp, header.headersize, 'bof');

    % Read data
    if strcmp(header.elementType, 'struct') || strcmp(header.elementType, 'cell')
        % Complex type
        data = readElement(fp);
        if ~isempty(data)
            data = reshape(data, header.dataDimSize);
        end

    else
        % Simple type with checksum
        checksum = fread(fp, 2, 'uint8');
        if checksum(1) ~= 1 || checksum(2) ~= 1
            error('LoadFast:ChecksumMismatch', ...
                'Version 1 format checksum failed');
        end

        % Calculate elements
        if header.dataSize > 0 && header.elementSize > 0
            num_elements = header.dataSize / header.elementSize;
        else
            num_elements = prod(header.dataDimSize);
        end

        % Map element type to proper fread precision, preserving class.
        % 'logical' (7 chars) is truncated to 'logica' by the 6-byte header.
        elemType = header.elementType;
        if strcmp(elemType, 'logica') || strcmp(elemType, 'logical')
            data = logical(fread(fp, num_elements, 'uint8=>uint8'));
        elseif strcmp(elemType, 'char')
            data = fread(fp, num_elements, '*char');
        else
            data = fread(fp, num_elements, ['*' elemType]);
        end

        data = reshape(data, header.dataDimSize);
    end

    % Return single variable with default name
    varNames = {'data'};
    varValues = {data};
end

function data = readElement(fp)
    % Recursive element reader with bit-flag decoding

    % Read type ID byte
    raw_type_id = fread(fp, 1, 'uint8');

    if isempty(raw_type_id) || raw_type_id == 0
        % Empty element
        if raw_type_id == 0
            % Read class name
            class_len = fread(fp, 1, 'uint8');
            className = fread(fp, class_len, '*char')';
            % Read dimensions
            ndims_val = fread(fp, 1, 'uint8');
            dims = fread(fp, ndims_val, 'double')';
            % Return empty array of that class with correct shape
            if strcmp(className, 'cell')
                data = cell(dims);
            elseif strcmp(className, 'struct')
                data = reshape(struct(), dims);
            else
                data = cast(zeros(dims), className);
            end
        else
            data = [];
        end
        return;
    end

    % Decode bit flags
    is_complex = bitand(raw_type_id, 64) > 0;   % bit 6
    is_sparse  = bitand(raw_type_id, 128) > 0;  % bit 7
    base_id    = bitand(raw_type_id, 63);        % bits 0-5

    % Dispatch on base type
    if is_sparse
        data = readSparse(fp, base_id, is_complex);
    elseif is_complex && base_id >= 1 && base_id <= 10
        data = readComplex(fp, base_id);
    elseif base_id >= 1 && base_id <= 12
        data = readSimple(fp, base_id);
    elseif base_id == 20
        data = readStructType(fp);
    elseif base_id == 21
        data = readCellType(fp);
    elseif base_id == 30
        data = readString(fp);
    elseif base_id == 31
        data = readTable(fp);
    elseif base_id == 32
        data = readTimetable(fp);
    elseif base_id == 33
        data = readCategorical(fp);
    elseif base_id == 34
        data = readDatetime(fp);
    elseif base_id == 35
        data = readDuration(fp);
    elseif base_id == 36
        data = readCalendarDuration(fp);
    elseif base_id == 37
        data = readContainersMap(fp);
    elseif base_id == 38
        data = readFunctionHandle(fp);
    else
        error('LoadFast:UnknownTypeID', 'Unknown type ID: %d (base: %d)', raw_type_id, base_id);
    end
end

%% Type-specific readers

function tp = baseIdToType(base_id)
    types = {'double','single','int8','uint8','int16','uint16', ...
             'int32','uint32','int64','uint64','char','logical'};
    if base_id >= 1 && base_id <= 12
        tp = types{base_id};
    else
        error('LoadFast:InvalidBaseId', 'Invalid base type ID: %d', base_id);
    end
end

function data = readSimple(fp, base_id)
    % Read simple numeric/char/logical array
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);
    tp = baseIdToType(base_id);

    if base_id == 11  % char
        data = fread(fp, num_elements, '*char');
    elseif base_id == 12  % logical
        data = logical(fread(fp, num_elements, 'uint8=>uint8'));
    else
        data = fread(fp, num_elements, ['*' tp]);
    end
    data = reshape(data, dims);
end

function data = readComplex(fp, base_id)
    % Read complex numeric array
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);
    tp = baseIdToType(base_id);

    re = fread(fp, num_elements, ['*' tp]);
    im = fread(fp, num_elements, ['*' tp]);
    data = complex(re, im);
    data = reshape(data, dims);
end

function data = readSparse(fp, base_id, is_complex)
    % Read sparse array
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    tp = baseIdToType(base_id);

    nz = fread(fp, 1, 'uint64');

    if nz == 0
        if base_id == 12
            data = sparse(dims(1), dims(2));
            data = data ~= 0; % sparse logical
        else
            data = sparse(dims(1), dims(2));
        end
        return;
    end

    rows = fread(fp, nz, 'uint64');
    cols = fread(fp, nz, 'uint64');

    if is_complex
        re = fread(fp, nz, ['*' tp]);
        im = fread(fp, nz, ['*' tp]);
        vals = complex(double(re), double(im));
    else
        if base_id == 12 % sparse logical
            vals = logical(fread(fp, nz, 'uint8=>uint8'));
        else
            vals = fread(fp, nz, ['*' tp]);
        end
    end

    data = sparse(double(rows), double(cols), vals, dims(1), dims(2));
end

function data = readStructType(fp)
    % Read struct array
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);

    % Read field names
    num_fields = fread(fp, 1, 'uint32');
    field_names = cell(1, num_fields);
    for i = 1:num_fields
        name_len = fread(fp, 1, 'uint8');
        field_names{i} = fread(fp, name_len, '*char')';
    end

    % Read values
    all_values = cell(num_elements, num_fields);
    for k = 1:num_elements
        for f = 1:num_fields
            all_values{k, f} = readElement(fp);
        end
    end

    % Convert to struct
    data = cell2struct(all_values, field_names, 2);
    data = reshape(data, dims);
end

function data = readCellType(fp)
    % Read cell array
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);

    data = cell(dims);
    for k = 1:num_elements
        data{k} = readElement(fp);
    end
end

function data = readString(fp)
    % Read string array
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);

    data = strings(dims);
    for k = 1:num_elements
        is_missing = fread(fp, 1, 'uint8');
        if is_missing
            data(k) = string(missing);
        else
            len = fread(fp, 1, 'uint32');
            if len > 0
                chars = fread(fp, len, '*char')';
                data(k) = string(chars);
            else
                data(k) = "";
            end
        end
    end
end

function data = readTable(fp)
    % Read table
    nrows = fread(fp, 1, 'uint64');
    nvars = fread(fp, 1, 'uint32');

    % Read variable names
    varNames = readElement(fp);

    % Row names
    has_rownames = fread(fp, 1, 'uint8');
    if has_rownames
        rowNames = readElement(fp);
    else
        rowNames = {};
    end

    % Variable units
    has_units = fread(fp, 1, 'uint8');
    if has_units
        units = readElement(fp);
    else
        units = {};
    end

    % Variable descriptions
    has_descr = fread(fp, 1, 'uint8');
    if has_descr
        descr = readElement(fp);
    else
        descr = {};
    end

    % Read columns
    cols = cell(1, nvars);
    for i = 1:nvars
        cols{i} = readElement(fp);
    end

    % Construct table
    data = table(cols{:}, 'VariableNames', varNames);

    if ~isempty(rowNames)
        data.Properties.RowNames = rowNames;
    end
    if ~isempty(units)
        data.Properties.VariableUnits = units;
    end
    if ~isempty(descr)
        data.Properties.VariableDescriptions = descr;
    end
end

function data = readTimetable(fp)
    % Read timetable
    nrows = fread(fp, 1, 'uint64'); %#ok<NASGU>
    nvars = fread(fp, 1, 'uint32');

    % Read variable names
    varNames = readElement(fp);

    % Read row times
    rowTimes = readElement(fp);

    % Variable units
    has_units = fread(fp, 1, 'uint8');
    if has_units
        units = readElement(fp);
    else
        units = {};
    end

    % Variable descriptions
    has_descr = fread(fp, 1, 'uint8');
    if has_descr
        descr = readElement(fp);
    else
        descr = {};
    end

    % Read columns
    cols = cell(1, nvars);
    for i = 1:nvars
        cols{i} = readElement(fp);
    end

    % Construct timetable
    data = timetable(cols{:}, 'RowTimes', rowTimes, 'VariableNames', varNames);

    if ~isempty(units)
        data.Properties.VariableUnits = units;
    end
    if ~isempty(descr)
        data.Properties.VariableDescriptions = descr;
    end
end

function data = readCategorical(fp)
    % Read categorical
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);

    is_ordinal = fread(fp, 1, 'uint8');
    is_protected = fread(fp, 1, 'uint8');

    num_cats = fread(fp, 1, 'uint32');
    cats = cell(1, num_cats);
    for i = 1:num_cats
        name_len = fread(fp, 1, 'uint32');
        cats{i} = fread(fp, name_len, '*char')';
    end

    codes = fread(fp, num_elements, 'uint32');

    % Reconstruct categorical
    % Map codes to category names; code 0 = undefined
    cellData = cell(dims);
    for k = 1:num_elements
        if codes(k) == 0
            cellData{k} = '';  % will become <undefined>
        else
            cellData{k} = cats{codes(k)};
        end
    end

    if is_ordinal
        data = categorical(cellData, cats, 'Ordinal', true);
    else
        data = categorical(cellData, cats);
    end

    if is_protected && ~is_ordinal
        data = setcats(data, cats);
        data = categorical(data, cats, 'Ordinal', false);
        % Protected but not ordinal: use addcats/setcats via protected categorical
        % Actually, protected categorical is created via ordinal or specific constructor
        % For now, just create with the categories specified
    end

    data = reshape(data, dims);
end

function data = readDatetime(fp)
    % Read datetime
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);

    % Timezone
    has_tz = fread(fp, 1, 'uint8');
    if has_tz
        tz_len = fread(fp, 1, 'uint8');
        tz = fread(fp, tz_len, '*char')';
    else
        tz = '';
    end

    % Format
    fmt_len = fread(fp, 1, 'uint8');
    if fmt_len > 0
        fmt = fread(fp, fmt_len, '*char')';
    else
        fmt = '';
    end

    % Read numeric data
    numData = fread(fp, num_elements, 'double');

    % Reconstruct
    if isempty(tz)
        % Saved as datenum
        data = datetime(numData, 'ConvertFrom', 'datenum');
    else
        % Saved as posixtime
        data = datetime(numData, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
    end

    if ~isempty(fmt)
        data.Format = fmt;
    end

    data = reshape(data, dims);
end

function data = readDuration(fp)
    % Read duration
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);

    % Format
    fmt_len = fread(fp, 1, 'uint8');
    if fmt_len > 0
        fmt = fread(fp, fmt_len, '*char')';
    else
        fmt = '';
    end

    % Read milliseconds
    ms = fread(fp, num_elements, 'double');
    data = milliseconds(ms);

    if ~isempty(fmt)
        data.Format = fmt;
    end

    data = reshape(data, dims);
end

function data = readCalendarDuration(fp)
    % Read calendarDuration
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);

    months_total = fread(fp, num_elements, 'int32');
    days_total = fread(fp, num_elements, 'int32');
    millis = fread(fp, num_elements, 'double');

    data = calmonths(double(months_total)) + caldays(double(days_total)) + milliseconds(millis);
    data = reshape(data, dims);
end

function data = readContainersMap(fp)
    % Read containers.Map
    num_entries = fread(fp, 1, 'uint32');

    % Key type
    kt_len = fread(fp, 1, 'uint8');
    keyType = fread(fp, kt_len, '*char')'; %#ok<NASGU>

    % Value type
    vt_len = fread(fp, 1, 'uint8');
    valueType = fread(fp, vt_len, '*char')'; %#ok<NASGU>

    % Keys and values
    k = readElement(fp);
    v = readElement(fp);

    if num_entries == 0
        data = containers.Map();
    else
        data = containers.Map(k, v);
    end
end

function data = readFunctionHandle(fp)
    % Read function_handle
    str_len = fread(fp, 1, 'uint16');
    funcStr = fread(fp, str_len, '*char')';
    data = str2func(funcStr);
end

%% Helper Functions

function [names, values] = filterVariables(allNames, allValues, requestedNames)
    % Filter to requested variable names

    names = {};
    values = {};

    for i = 1:length(requestedNames)
        reqName = requestedNames{i};
        idx = find(strcmp(allNames, reqName), 1);

        if isempty(idx)
            warning('LoadFast:VariableNotFound', ...
                'Variable "%s" not found in file', reqName);
        else
            names{end+1} = allNames{idx}; %#ok<AGROW>
            values{end+1} = allValues{idx}; %#ok<AGROW>
        end
    end

    if isempty(names)
        error('LoadFast:NoVariablesLoaded', ...
            'None of the requested variables were found');
    end
end

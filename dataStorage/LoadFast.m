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
% Features:
%   - Multi-variable support (like MATLAB load)
%   - Backward compatible with original LoadMatFast (version 1)
%   - Can load into struct or workspace
%   - Selective variable loading
%   - Format version detection
%
% See also: SaveFast, SaveMatFast, LoadMatFast

% Version: 2.0
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
            % Version 2 format
            [varNames, varValues] = readFormatV2(fp);

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

function [varNames, varValues] = readFormatV2(fp)
    % Read Version 2 format (multi-variable)

    % Read version
    formatVersion = fread(fp, 1, 'uint16');
    if formatVersion ~= 2
        warning('LoadFast:VersionMismatch', ...
            'Expected version 2, got version %d', formatVersion);
    end

    % Read number of variables
    numVars = fread(fp, 1, 'uint32');

    % Read timestamp (informational only)
    timestamp_len = fread(fp, 1, 'uint8');
    timestamp = fread(fp, timestamp_len, '*char')';
    % (timestamp available but not used)

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
    % Recursive element reader

    % Read type ID
    type_id = fread(fp, 1, 'uint8');

    if isempty(type_id) || type_id == 0
        % Empty element
        if type_id == 0
            % Read class name
            class_len = fread(fp, 1, 'uint8');
            className = fread(fp, class_len, '*char')';
            % Return empty array of that class
            data = feval(className, []);
        else
            data = [];
        end
        return;
    end

    % Read dimensions
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);

    % Type mapping
    typeMap = {1,'double'; 2,'single'; 3,'int8'; 4,'uint8'; 5,'int16'; ...
               6,'uint16'; 7,'int32'; 8,'uint32'; 9,'int64'; 10,'uint64'; ...
               11,'char'; 12,'logical'; 20,'struct'; 21,'cell'};

    tp = '';
    for i = 1:size(typeMap, 1)
        if typeMap{i,1} == type_id
            tp = typeMap{i,2};
            break;
        end
    end

    if isempty(tp)
        error('LoadFast:UnknownTypeID', 'Unknown type ID: %d', type_id);
    end

    % Read data based on type
    if type_id < 20
        % Simple type
        if type_id == 11  % char
            data = fread(fp, num_elements, '*char');
        elseif type_id == 12  % logical
            data = fread(fp, num_elements, 'uint8=>uint8');
            data = logical(data);
        else
            data = fread(fp, num_elements, ['*' tp]);
        end
        data = reshape(data, dims);

    elseif type_id == 20
        % Struct
        data = readStruct(fp, dims);

    elseif type_id == 21
        % Cell
        data = readCell(fp, dims);
    end
end

function data = readStruct(fp, dims)
    % Read struct array

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

function data = readCell(fp, dims)
    % Read cell array

    num_elements = prod(dims);
    data = cell(dims);

    for k = 1:num_elements
        data{k} = readElement(fp);
    end
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

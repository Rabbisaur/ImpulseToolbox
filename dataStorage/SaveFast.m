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
% Features:
%   - Multi-variable support (like MATLAB save)
%   - Backward compatible with original SaveMatFast
%   - Format versioning for future compatibility
%   - Better error handling and validation
%   - Supports all MATLAB data types (numeric, char, logical, struct, cell)
%
% File Format Version 2:
%   Header:
%     - Magic bytes: [77 70] ('MF' for MatFast)
%     - Version: 2
%     - Number of variables
%     - Timestamp
%   For each variable:
%     - Name length + name
%     - Data (using recursive serialization)
%
% See also: LoadFast, SaveMatFast, LoadMatFast

% Version: 2.0
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
        fwrite(fp, uint16(2), 'uint16');

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
                warning('SaveFast:NameTooLong', ...
                    'Variable name truncated to 255 chars: %s', varName);
                varName = varName(1:255);
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

function writeElement(fp, data)
    % Recursive element writer with type dispatch

    if isempty(data)
        % Empty array
        writeEmptyElement(fp, class(data));
        return;
    end

    tp = class(data);
    dims = size(data);
    ndims_val = numel(dims);

    % Type ID mapping
    typeMap = containers.Map(...
        {'double','single','int8','uint8','int16','uint16','int32','uint32','int64','uint64','char','logical','struct','cell'}, ...
        {1,      2,       3,     4,      5,      6,       7,      8,       9,      10,      11,    12,       20,      21});

    if ~isKey(typeMap, tp)
        warning('SaveFast:UnsupportedType', ...
            'Unsupported type: %s. Saving as empty.', tp);
        writeEmptyElement(fp, tp);
        return;
    end

    type_id = typeMap(tp);

    % Write type ID
    fwrite(fp, uint8(type_id), 'uint8');

    % Write dimensions
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');

    % Write data based on type
    if type_id < 20
        % Simple types
        if type_id == 12 % logical
            fwrite(fp, uint8(data(:)), 'uint8');
        else
            fwrite(fp, data(:), tp);
        end

    elseif type_id == 20
        % Struct
        writeStruct(fp, data);

    elseif type_id == 21
        % Cell array
        writeCell(fp, data);
    end
end

function writeEmptyElement(fp, className)
    % Write empty element marker
    fwrite(fp, uint8(0), 'uint8'); % Type ID 0 = empty
    fwrite(fp, uint8(length(className)), 'uint8');
    fwrite(fp, className, 'char*1');
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

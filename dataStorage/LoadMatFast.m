function [datamat, header] = LoadMatFast(filepath)

% this function defines a binary file format to store large matrix much
% faster than save function in the matlab
% use with SaveMatFast
%
% Updates:
% - Supports struct and cell arrays via recursive binary deserialization.

if ~exist(filepath,'file')
    error([filepath,' does not exist.'])
end

[fp, errormsg] = fopen(filepath,'r');
if ~isempty(errormsg)
    error(['Failed to open file ', filepath])
end

% read out header
header.identificationBytes = fread(fp,2,'uint8');
header.headersize = fread(fp,1,'uint8');
header.dateBytes = fread(fp,14,'*char');
header.dataSize = fread(fp,1,'double');
header.elementSize = fread(fp,1,'uint8');
header.elementType = fread(fp,6,'*char')';
% Clean up padding
header.elementType = strtrim(header.elementType);

header.dataNumberDimensions = fread(fp,1,'double');
header.dataDimSize = fread(fp,header.dataNumberDimensions,'double')'; % Ensure row vector
header.commentsize = fread(fp,1,'uint8');
header.comment = fread(fp,header.commentsize,'*char');

% Check if we are at the end of the header
currPos = ftell(fp);
if currPos ~= header.headersize
    % Warning or adjustment if header calculation differs (legacy compatibility)
    fseek(fp, header.headersize, -1);
end

% read data
if strcmp(header.elementType, 'struct') || strcmp(header.elementType, 'cell')
    % Recursive Load
    datamat = read_element(fp);
    
    % Reshape root object (recursive loader handles internal reshaping, but root needs it)
    if ~isempty(datamat)
        datamat = reshape(datamat, header.dataDimSize);
    end
else
    % Legacy/Simple Load
    % check checksum
    checksum = double(fread(fp,2,'uint8'));
    if checksum(1)~=1 || checksum(2) ~= 1
        fclose(fp);
        error('Check sum mismatch!')
    end

    % Calculate total elements
    if header.dataSize > 0 && header.elementSize > 0
        num_elements = header.dataSize / header.elementSize;
    else
        num_elements = prod(header.dataDimSize);
    end
    
    datamat = fread(fp, num_elements, header.elementType);
    
    datamat = reshape(datamat, header.dataDimSize);
    
    if strcmp(header.elementType, 'logical')
         datamat = logical(datamat);
    end
end

fclose(fp);
end

function data = read_element(fp)
    % Read Type ID
    type_id = fread(fp, 1, 'uint8');
    
    if isempty(type_id)
        data = [];
        return;
    end
    
    if type_id == 0 
        data = []; 
        return; 
    end
    
    % Read Dims
    ndims_val = fread(fp, 1, 'uint8');
    dims = fread(fp, ndims_val, 'double')';
    num_elements = prod(dims);
    
    switch type_id
        case 1,  tp = 'double'; sz=8;
        case 2,  tp = 'single'; sz=4;
        case 3,  tp = 'int8';   sz=1;
        case 4,  tp = 'uint8';  sz=1;
        case 5,  tp = 'int16';  sz=2;
        case 6,  tp = 'uint16'; sz=2;
        case 7,  tp = 'int32';  sz=4;
        case 8,  tp = 'uint32'; sz=4;
        case 9,  tp = 'int64';  sz=8;
        case 10, tp = 'uint64'; sz=8;
        case 11, tp = 'char';   sz=1;
        case 12, tp = 'logical';sz=1;
        case 20, tp = 'struct';
        case 21, tp = 'cell';
        otherwise, tp = 'unknown';
    end
    
    if type_id < 20
        % Simple Type
        if type_id == 11
            data = fread(fp, num_elements, '*char');
        elseif type_id == 12
            data = fread(fp, num_elements, 'uint8=>uint8');
            data = logical(data);
        else
            data = fread(fp, num_elements, ['*' tp]);
        end
        data = reshape(data, dims);
        
    elseif type_id == 20 % Struct
        % Read Fields
        num_fields = fread(fp, 1, 'uint32');
        field_names = cell(1, num_fields);
        for i = 1:num_fields
            nm_len = fread(fp, 1, 'uint8');
            field_names{i} = fread(fp, nm_len, '*char')';
        end
        
        % Pre-allocate struct
        % Initialize with last field to preallocate
        data = struct(); 
        % We can't easily preallocate a struct array with unknown values without initialization.
        % To be fast, we might need to handle the array loop carefully.
        % If we just loop and assign data(k).field = val, MATLAB handles it but it might be slow for huge arrays.
        % However, "struct" serialization in this format is: Array Order -> Field Order
        
        % Efficient Pre-allocation strategy:
        % Use cell arrays to collect data, then converting to struct is usually faster than incremental growth.
        
        % Collecting data:
        % We need to read (num_elements * num_fields) items.
        % They are stored: Element 1 (Field 1, Field 2...), Element 2...
        
        total_items = num_elements * num_fields;
        all_values = cell(num_elements, num_fields);
        
        for k = 1:num_elements
            for f = 1:num_fields
                all_values{k, f} = read_element(fp);
            end
        end
        
        % Construct Struct
        data = cell2struct(all_values, field_names, 2);
        data = reshape(data, dims);
        
    elseif type_id == 21 % Cell
        data = cell(dims);
        for k = 1:num_elements
            data{k} = read_element(fp);
        end
    end
end

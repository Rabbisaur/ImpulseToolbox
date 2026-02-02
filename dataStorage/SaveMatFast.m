function SaveMatFast(filepath,datamat,comment)

% this function defines a binary file format to store large matrix much
% faster than save function in the matlab
% use with LoadMatFast
%
% Updates:
% - Supports struct and cell arrays via recursive binary serialization.
% - Optimized for large embedded numerical matrices.
%
% Usage: SaveMatFast(filepath, variable, comment)

if nargin == 2
    comment = 'nil';
end

if ~ischar(comment)
    comment = char(comment);
end

if numel(comment) > 42
    comment = comment(1:42);
end
if numel(comment) < 42
    tmp = repmat(' ',1,42);
    tmp(1:numel(comment)) = comment;
    comment = tmp;
end

% header 64 bytes
% identification 2 bytes
header.identificationBytes = [1 1];

% date bytes yymmddhhmmss 14 bytes char
header.dateBytes = datestr(now,'yyyymmddHHMMSS');

datatype = class(datamat);
% Handle mapping of class names to header element types (max 6 chars)
switch datatype
    case {'double','single','int8','uint8','int16','uint16','int32','uint32','int64','uint64','char','logical'}
        header.elementType = pad_type(datatype, 6);
        header.elementSize = get_element_size(datatype);
        is_complex_type = false;
        header.dataSize = numel(datamat) * header.elementSize;
        
    case 'struct'
        header.elementType = pad_type('struct', 6);
        header.elementSize = 0; % Unknown/Variable
        is_complex_type = true;
        header.dataSize = 0; % Unknown/Variable
        
    case 'cell'
        header.elementType = pad_type('cell', 6);
        header.elementSize = 0; % Unknown/Variable
        is_complex_type = true;
        header.dataSize = 0; % Unknown/Variable
        
    otherwise
        error(['Unsupported data type: ' datatype]);
end

header.dataNumberDimensions = numel(size(datamat)); 
header.dataDimSize = size(datamat); 

% comment 1 bytes
header.commentsize = numel(comment);
% comment max 42 char
header.comment = comment;

% headersize calculation
% Fixed parts:
% id(2) + headersize(1) + date(14) + dataSize(8) + elementSize(1) + elementType(6) + ndims(8)
% + dimArray(ndims*8) + commentSize(1) + comment(42)
header.headersize = 2 + 1 + 14 + 8 + 1 + 6 + 8 + header.dataNumberDimensions*8 + 1 + 42;

[fp,errormsg]= fopen(filepath,'w');
if ~isempty(errormsg)
    error(['Error openning file: ', errormsg])
end

% write header
fwrite(fp,uint8(header.identificationBytes),'uint8');
fwrite(fp,uint8(header.headersize),'uint8');
fwrite(fp,char(header.dateBytes),'char*1');
fwrite(fp,double(header.dataSize),'double');
fwrite(fp,uint8(header.elementSize),'uint8');
fwrite(fp,char(header.elementType),'char*1');
fwrite(fp,double(header.dataNumberDimensions),'double');
fwrite(fp,double(header.dataDimSize),'double');

fwrite(fp,uint8(header.commentsize),'uint8');
fwrite(fp,char(header.comment),'char*1');

% Check position
currPos = ftell(fp);
if currPos~=header.headersize
    fclose(fp);
    error(['currPos=',num2str(currPos), ' but header.headersize=', num2str(header.headersize)])
end

% write data
if is_complex_type
    % Recursively write structure/cell
    write_element(fp, datamat);
else
    % write checksum for simple types (Legacy compatibility)
    fwrite(fp,uint8([1 1]),'uint8');
    % Fast write for simple matrix
    if strcmp(datatype, 'logical')
         fwrite(fp, uint8(datamat), 'uint8');
    else
         fwrite(fp, datamat, datatype);
    end
end

fclose(fp);
end

%% Helper Functions

function type_str = pad_type(str, len)
    type_str = repmat(' ', 1, len);
    n = min(length(str), len);
    type_str(1:n) = str(1:n);
end

function sz = get_element_size(type)
    switch type
        case {'char','int8','uint8','logical'}, sz = 1;
        case {'int16','uint16'}, sz = 2;
        case {'int32','uint32','single'}, sz = 4;
        case {'int64','uint64','double'}, sz = 8;
        otherwise, sz = 1; 
    end
end

function write_element(fp, data)
    % Element Header:
    % Type ID (1 byte)
    % ndims (1 byte) -> optimized from 8 bytes, assuming < 255 dims
    % dims (ndims * 8 bytes double)
    
    tp = class(data);
    dims = size(data);
    ndims_val = numel(dims);
    
    % Type Mapping
    % 1=double, 2=single, 3=int8, 4=uint8, 5=int16, 6=uint16, 
    % 7=int32, 8=uint32, 9=int64, 10=uint64, 11=char, 12=logical
    % 20=struct, 21=cell
    
    switch tp
        case 'double', type_id = 1;
        case 'single', type_id = 2;
        case 'int8',   type_id = 3;
        case 'uint8',  type_id = 4;
        case 'int16',  type_id = 5;
        case 'uint16', type_id = 6;
        case 'int32',  type_id = 7;
        case 'uint32', type_id = 8;
        case 'int64',  type_id = 9;
        case 'uint64', type_id = 10;
        case 'char',   type_id = 11;
        case 'logical',type_id = 12;
        case 'struct', type_id = 20;
        case 'cell',   type_id = 21;
        otherwise
            warning(['Skipping unsupported type inside struct/cell: ' tp]);
            type_id = 0; % Unknown
    end
    
    fwrite(fp, uint8(type_id), 'uint8');
    
    if type_id == 0
        return; 
    end
    
    fwrite(fp, uint8(ndims_val), 'uint8');
    fwrite(fp, double(dims), 'double');
    
    % Write Body
    if type_id < 20
        % Simple Type: Fast Write
        if type_id == 12 % logical stored as uint8
            fwrite(fp, uint8(data), 'uint8');
        else
            fwrite(fp, data, tp);
        end
    elseif type_id == 20 % Struct
        % Write Field Names
        field_names = fieldnames(data);
        num_fields = numel(field_names);
        fwrite(fp, uint32(num_fields), 'uint32');
        
        for i = 1:num_fields
            fname = field_names{i};
            fwrite(fp, uint8(length(fname)), 'uint8');
            fwrite(fp, fname, 'char*1');
        end
        
        % Write Values (Field-major or Element-major? Element-major matches struct indexing)
        % For optimal compatibility with massive arrays, we should iterate fields then elements?
        % No, data(i).field is standard.
        % Actually, MATLAB structs are column-major arrays.
        % We will iterate linearly through the struct array.
        % For huge struct arrays with small fields, this is overhead. 
        % But for a single struct with huge fields (common case), this is fine.
        
        % Optimization: Structure of Arrays vs Array of Structures.
        % MATLAB 'struct' is Structure of Arrays internally if fields are scalar? No.
        % Let's stick to simple serialization: Linear index of struct array, then each field.
        
        num_elements = numel(data);
        for k = 1:num_elements
            for f = 1:num_fields
                val = data(k).(field_names{f});
                write_element(fp, val);
            end
        end
        
    elseif type_id == 21 % Cell
        num_elements = numel(data);
        for k = 1:num_elements
            write_element(fp, data{k});
        end
    end
end
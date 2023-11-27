function SaveMatFast(filepath,datamat,comment)

% this function defines a binary file format to store large matrix much
% faster than save function in the matlab
% use with LoadMatFast

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
    tmp = repmat(' ',42,1);
    tmp(1:numel(comment)) = comment;
    comment = tmp;
end

% header 64 bytes
% identification 2 bytes
header.identificationBytes = [1 1];

% date bytes yymmddhhmmss 14 bytes char
header.dateBytes = datestr(datetime,'yyyymmddHHMMSS');
% data size in bytes, double, 8 bytes * numel (8 bytes)
datatype = class(datamat);
switch datatype
    case 'char'
        header.elementSize = 1;
    case 'uint8'
        header.elementSize = 1;
    case 'int8'
        header.elementSize = 1;
    case 'int16'
        header.elementSize = 2;
    case 'uint16'
        header.elementSize = 2;
    case 'int32'
        header.elementSize = 4;
    case 'uint32'
        header.elementSize = 4; 
    case 'int64'
        header.elementSize = 8;
    case 'uint64'
        header.elementSize = 8;
    case 'single'
        header.elementSize = 4;
    case 'double'
        header.elementSize = 8;
end
header.dataSize = numel(datamat)*header.elementSize;
header.dataNumberDimensions = numel(size(datamat)); % 8 bytes
header.dataDimSize = size(datamat); % n by 8 bytes
% comment 1 bytes
header.commentsize = numel(comment);
% comment max 42 char
header.comment = comment;

% headersize 1 byte
header.headersize = 2 +1+ 14 + 8 + 8 + 1 + 6 + header.dataNumberDimensions*8 + 1 + 42;

[fp,errormsg]= fopen(filepath,'w');
if ~isempty(errormsg)
    error(['Error openning file: ', errormsg])
end

header.elementType = repmat(' ',6,1);
elementType = class(datamat);
header.elementType(1:numel(elementType)) = elementType;

% write header
fwrite(fp,uint8(header.identificationBytes),'uint8'); % 2
fwrite(fp,uint8(header.headersize),'uint8'); % 1
fwrite(fp,char(header.dateBytes),'char'); % 14
fwrite(fp,double(header.dataSize),'double'); % 8
fwrite(fp,uint8(header.elementSize),'uint8'); % 1
fwrite(fp,char(header.elementType),'char'); % 6
fwrite(fp,double(header.dataNumberDimensions),'double'); % header.dataNumberDimensions *8
fwrite(fp,double(header.dataDimSize),'double');

fwrite(fp,uint8(header.commentsize),'uint8');
fwrite(fp,char(header.comment),'char');

% fseek(fp,header.headersize,-1);
currPos = ftell(fp);

if currPos~=header.headersize
    fclose(fp);
    error(['currPos=',num2str(currPos), ' but header.headersize=', num2str(header.headersize)])
end

% write checksum
fwrite(fp,uint8([1 1]),'uint8');

% write data
fwrite(fp,datamat,elementType);

fclose(fp);
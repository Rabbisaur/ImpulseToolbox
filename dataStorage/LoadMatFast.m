function [datamat, header] = LoadMatFast(filepath)

% this function defines a binary file format to store large matrix much
% faster than save function in the matlab
% use with SaveMatFast

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
idx = header.elementType == ' ';
header.elementType(idx) = '';
header.dataNumberDimensions = fread(fp,1,'double');
header.dataDimSize = fread(fp,header.dataNumberDimensions,'double');
header.commentsize = fread(fp,1,'uint8');
header.comment = fread(fp,header.commentsize,'*char');

% read data
fseek(fp,header.headersize,-1);
% check checksum
checksum = double(fread(fp,2,'uint8'));

if checksum(1)~=1 || checksum(2) ~= 1
    fclose(fp);
    error('Check sum mismatch!')
end

datamat = fread(fp,header.dataSize,header.elementType);


datamat = reshape(datamat,header.dataDimSize');

fclose(fp);

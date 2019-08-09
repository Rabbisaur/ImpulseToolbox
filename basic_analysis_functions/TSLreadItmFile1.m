function items = TSLreadItmFile(path)

% path = '/home/wangfeng/storage/projects/ncscmovie/stimulus/NCSCMOV.ITM';

if exist(path,'file') ~= 2
    error(['ITM file ',path,' does not exist, check spelling.'])
end

% open file read only
fid = fopen(path,'r');

tline = fgetl(fid);
while tline ~= -1
%     disp(tline)
    numinline = textscan(tline,'%d');
    numinline = numinline{1};
    if ~isempty(numinline)
        itmID = numinline(1);
    end
%     strinline = textscan(tline,'%s');
%     itmName = strinline{1}(end);
    tline = fgetl(fid);
end

NumConditions = itmID;
itmName = cell(NumConditions,1);

fseek(fid,0,-1);

tline = fgetl(fid);
while tline ~= -1
%     disp(tline)
    numinline = textscan(tline,'%d');
    numinline = numinline{1};
    if ~isempty(numinline)
        itmID = numinline(1);
    end
    if itmID > 0
        strinline = textscan(tline,'%s');
        itmName{itmID} = strinline{1}{end};
    end
    tline = fgetl(fid);
end
fclose(fid);

items.numitems = NumConditions;
items.itmName = itmName;
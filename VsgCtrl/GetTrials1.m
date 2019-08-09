function [b_Trialidx] = GetTrials(sMbmpath,cMbmvect)
sMatDir = sMbmpath(1:end-4);
% load Expmark
FileName = [sMatDir '/' 'Expmark.mat'];
if exist(FileName,'file')
    load(FileName);
else
    fprintf('%s donesn''t exist!\nCannot Proceed!\n',FileName);
    error('Expmark.mat not found!');
end
% load MbmInfo
FileName = [sMatDir '/' 'mbminfo.mat'];
if exist(FileName,'file') == 2
    load(FileName);
else
    fprintf('%s donesn''t exist!\nCannot Proceed!\n',FileName);
    error('mbminfo.mat not found!');
end
% get valid trials, no abortion & have vsgtrigger
b_ValidTrials = ~(Expmark(3,:) > 0  | Expmark(5,:) == 0);
% get trial idx taken mbm into consideration
b_Trialidx = b_ValidTrials;
for j = 1:numel(cMbmvect)
    if ~isempty(cMbmvect{j})
%        idx = MbmInfo(1,b_ValidTrials)<0;
        idx = MbmInfo(1,:)<0;
        for i = 1:numel(cMbmvect{j})
%            idx = idx | MbmInfo(i,b_ValidTrials) == cMbmvect{j}(i);
        idx = idx | MbmInfo(j,:) == cMbmvect{j}(i);
        end
        b_Trialidx = idx & b_Trialidx;
    end
end
return
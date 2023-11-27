function explog = VnCReadLogFile(filepath)
% filepath = 'F:\RF_mapping\091118_B1\091118_data\091118_B1_runstim.log';
fp = fopen(filepath,'r');
explog.name = fgetl(fp);
thisline = fgetl(fp);

% Read out parameters
while thisline(1) ~= '#'
thisline = fgetl(fp);
end
numpara = 0;
while 1
thisline = fgetl(fp);
if thisline(1) == '*'
    break
end
numpara = numpara + 1;
paraname = textscan(thisline,'%s');
paraname = paraname{1};
idx = find(thisline == '=');
paravalue = textscan(thisline((idx+1):end),'%f','delimiter',',');
paravalue = paravalue{1};

cmd = ['explog.para.', paraname{1}, '= paravalue;'];
eval(cmd);
end
thisline = fgetl(fp);
tmp = textscan(thisline,'%s','delimiter',';');
tmp = tmp{1};
tmp(end) = [];
for thisCell = 1:numel(tmp)
    idx = tmp{thisCell} == '*' | tmp{thisCell} == '#';
    tmp{thisCell}(idx) = [];
end
fieldnames = tmp;
thisTrial = 0;
while ~feof(fp)
    thisTrial = thisTrial + 1;
    thisline = fgetl(fp);
    idx = strfind(thisline,'Error = ');
    idx2 = find(thisline == ';',1,'last');
    explog.error{thisTrial} = thisline(idx:idx2);
    tmp = textscan(thisline,'%f','delimiter',';');
    tmp = tmp{1};

    if numel(tmp) == numel(fieldnames)
        for thisCell = 1:numel(tmp)
            cmd = ['explog.',fieldnames{thisCell},'(', num2str(thisTrial), ')=',num2str(tmp(thisCell)),';'];
            eval(cmd);
        end
    else
        error('Number of data does not match number of fields')
    end
end

fclose(fp);
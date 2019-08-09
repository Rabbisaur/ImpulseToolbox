function data = readRFdata(filepath)
% filepath = './fovea10by10by10.txt';

fp = fopen(filepath,'r');

linenum = 0;
beginFlag = false;
while(~feof(fp))
    thisline = fgetl(fp);
    if ~beginFlag
        thisline = textscan(thisline,'%s');
        if thisline{1}{1} == '#'
            beginFlag = 1;
        end
    else
        linenum = linenum + 1;
        thisline = textscan(thisline,'%f');
        thisline = thisline{1};
        ElecID{linenum} = thisline(1);
        ArrayCoord(:,linenum) = [thisline(2);thisline(3);thisline(4)];
        WorldCoord(:,linenum) = [thisline(5);thisline(6);thisline(7)];
        VoxelCoord(:,linenum) = [thisline(8);thisline(9);thisline(10)];
        pRF(:,linenum) = [thisline(11);thisline(12);thisline(13)];
    end
end
fclose(fp);
data.ElecID = ElecID;
data.ArrayCoord = ArrayCoord;
data.WorldCoord = WorldCoord;
data.VoxelCoord = VoxelCoord;
data.pRF = pRF;
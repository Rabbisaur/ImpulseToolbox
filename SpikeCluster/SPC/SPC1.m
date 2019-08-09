function [Class Model] = SPC(Input, Temprature, SWCycles, KNearNeighb)
% Maple to Matlab for Superparamagnetic Clustering
%   Usage: [Class Model] = SPC(Input, Temprature, SWCycles, KNearNeighb)

% Write data file in disk
Input =double(Input);
Input = zscore(Input, 0, 1);
DataFileName = 'SPC.dat';
save(DataFileName, 'Input', '-ascii');
% Write parameter file in disk
[NumRow NumCol] = size(Input);
MinTemp = min(Temprature);
MaxTemp = max(Temprature);
TempStepVal = (MaxTemp-MinTemp)/(numel(Temprature)-1);
ParaFileName = 'SPC.para';
FileID = fopen(ParaFileName, 'wt');
fprintf(FileID, 'NumberOfPoints: %s\n', num2str(NumRow));
fprintf(FileID, 'DataFile: %s\n', DataFileName);
fprintf(FileID, 'OutFile: %s\n', 'SPC');
fprintf(FileID, 'Dimensions: %s\n', num2str(NumCol));
fprintf(FileID, 'MinTemp: %s\n', num2str(MinTemp));
fprintf(FileID, 'MaxTemp: %s\n', num2str(MaxTemp));
fprintf(FileID, 'TempStep: %s\n', num2str(TempStepVal));
fprintf(FileID, 'SWCycles: %s\n', num2str(SWCycles));
fprintf(FileID, 'KNearestNeighbours: %s\n', num2str(KNearNeighb));
fprintf(FileID, 'MSTree|\n');
fprintf(FileID, 'DirectedGrowth|\n');
fprintf(FileID, 'SaveSuscept|\n');
fprintf(FileID, 'WriteLables|\n');
fprintf(FileID, 'WriteCorFile~\n');
fprintf(FileID, 'ForceRandomSeed: %s\n', num2str(randi(10000, 1)));   
fclose(FileID);
% Clustering
System = computer;
switch System
    case {'PCWIN';'PCWIN64'}
        ExePath = [pwd '\SPCWin.exe'];
        if exist(ExePath, 'file')==2
            delete(ExePath);
        end
        ExePath = which('SPCWin.exe');
        copyfile(ExePath, pwd, 'f');
        dos(['SPCWin.exe ' ParaFileName]);
    case {'MAC'}
        ExePath = [pwd '/SPCMac.exe'];
        if exist(ExePath, 'file')==2
            delete(ExePath);
        end
        ExePath = which('SPCMac.exe');
        copyfile(ExePath, pwd, 'f');
	    unix(['./SPCMac.exe ' ParaFileName]);
   case {'MACI';'MACI64'}
       ExePath = [pwd '/SPCMaci.exe'];
       if exist(ExePath, 'file')==2
           delete(ExePath);
       end
       ExePath = which('SPCMaci.exe');
       copyfile(ExePath, pwd, 'f');
	   unix(['./SPCMaci.exe ' ParaFileName]);
    case {'GLNX86';'GLNXA64';'GLNXI64'}
       ExePath = [pwd '/SPCLinux.exe'];
       if exist(ExePath, 'file')==2
           delete(ExePath);
       end
       ExePath = which('SPCLinux.exe');
       copyfile(ExePath, pwd, 'f');
	   unix(['./SPCLinux.exe ' ParaFileName]);
end
Class = load('SPC.dg_01.lab');
Model = load('SPC.dg_01');
% Delete all the temporary files
delete SPC*.*
clc;
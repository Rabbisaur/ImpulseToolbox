function [IsoI, IsoDist, LRatio] = IsoI(Wav, Unit, Fs)
% Maple to Matlab for the Isolation Information quality measures, the Isolation
% Distance and LRatio quality (Neymotin et al., J Neurosci 2011)
%   Usage: [IsoI IsoDist LRatio] = IsoI(Wav, Unit, Fs)

% Write data file in disk
[NumRow, NumCol] = size(Wav);
WavexName = 'IsoIwavex.txt';
FileID = fopen(WavexName, 'wt');
fprintf(FileID, '%s\n', num2str(1));
fprintf(FileID, '%s\n', num2str(NumCol));
fprintf(FileID, '%s\n', num2str(Fs));
fclose(FileID);
Wav = double([Unit Wav]);
save(WavexName, 'Wav', '-ascii', '-append');
% Computing isolation quality
System = computer;
switch System
    case {'PCWIN';'PCWIN64'}
        Wave2fPath = [pwd '\IsoIwave2f.exe'];
        IsoRatPath = [pwd '\IsoIisorat.exe'];
        IsoIPath = [pwd '\IsoIisoi.exe'];
        if exist(Wave2fPath, 'file')==2 && exist(IsoRatPath, 'file')==2 && exist(IsoIPath, 'file')==2
            delete(Wave2fPath);
            delete(IsoRatPath);
            delete(IsoIPath);
        end
        copyfile(which('IsoIwave2f.exe'), pwd, 'f');
        copyfile(which('IsoIisorat.exe'), pwd, 'f');
        copyfile(which('IsoIisoi.exe'), pwd, 'f');
        dos('IsoIwave2f IsoIwavex.txt IsoIwavefeat.txt');
        dos('IsoIisorat IsoIwavefeat.txt IsoIwave_IsoDLRat.txt');
        dos('IsoIisoi IsoIwavefeat.txt IsoIwave_IsoI.txt');
    case {'MAC'}
    case {'MACI';'MACI64'}
    case {'GLNX86';'GLNXA64';'GLNXI64'}
end
Class = load('SPC.dg_01.lab');
Model = load('SPC.dg_01');



% Delete all the temporary files
delete IsoI*.* _iso*
clc;
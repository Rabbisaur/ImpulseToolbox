clc, clear
%------Configure
Path = './test.gsq';

SeqLength = 1000;

CenterXLow = 3;
CenterXHi = 3;

CenterYLow = 3;
CenterYHi = 3;

SizeXLow = 5;
SizeXHi = 5;

SizeYLow = 5;
SizeYHi = 5;

SpFreqLow = 0.1;
SpFreqHi = 8.0;

TempFreqLow = 0;
TempFreqHi = 4;

OrientationLow = 0;
OrientationHi = 360;

PhaseLow = 0;
PhaseHi = 360;

ContrastLow = 0;
ContrastHi = 0.8;

MeanLumLow = 0.1;
MeanLumHi = 0.9;

BGLumLow = 0.5;
BGLumHi = 0.5;

FramesLow = 5;
FramesHi = 5;


%------Generation of Sequence
Seq.CenterX = (CenterXHi-CenterXLow)*rand(SeqLength,1)+ CenterXLow;
Seq.CenterY = (CenterYHi-CenterYLow)*rand(SeqLength,1)+ CenterYLow;
Seq.SizeX = (SizeXHi-SizeXLow)*rand(SeqLength,1)+ SizeXLow;
Seq.SizeY = (SizeYHi-SizeYLow)*rand(SeqLength,1)+ SizeYLow;
Seq.SpFreq = (SpFreqHi - SpFreqLow)*rand(SeqLength,1)+ SpFreqLow;
Seq.TempFreq = (TempFreqHi-TempFreqLow)*rand(SeqLength,1)+ TempFreqLow;
Seq.Orientation = (OrientationHi - OrientationLow)*rand(SeqLength,1)+ OrientationLow;
Seq.Phase = (PhaseHi - PhaseLow)*rand(SeqLength,1)+ PhaseLow;
Seq.Contrast = (ContrastHi - ContrastLow)*rand(SeqLength,1)+ ContrastLow;
Seq.MeanLum = (MeanLumHi- MeanLumLow)*rand(SeqLength,1)+ MeanLumLow;
Seq.BGLum = (BGLumHi- BGLumLow)* rand(SeqLength,1)+ BGLumLow;
Seq.Frames = round((FramesHi- FramesLow)* rand(SeqLength,1)+ FramesLow);

mSeq = struct2array(Seq);

%-------Write Sequence File
fileID = fopen(Path,'w');
fprintf(fileID,'RTSGratings Sequence File\n');
fprintf(fileID,'CenterX, CenterY, SizeX, SizeY, SpFreq, TempFreq, Orientation, Phase, Contrast, MeanLum, BGLum, Frames;\n');
fprintf(fileID,'SequenceLegnth = %d\n',SeqLength);
fprintf(fileID,'[DATA Begin]\n');
for i=1:SeqLength
    fprintf(fileID,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f;\n',mSeq(i,1),mSeq(i,2),mSeq(i,3),mSeq(i,4),mSeq(i,5),mSeq(i,6),mSeq(i,7),mSeq(i,8),mSeq(i,9),mSeq(i,10),mSeq(i,11),mSeq(i,12));
end
fprintf(fileID,'[DATA End]\n');
fclose(fileID);
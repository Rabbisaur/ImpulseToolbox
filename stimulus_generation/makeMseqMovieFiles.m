
n = 16;
R0 = 1;
T0 = 57;
c = 32;
r = 32;
p = 2^n/c/r;
framespertrial = p;
machineframespermovieframe = 1;
machinepixelspermoviepixel = 2;
filename = './moviefiles/mseq';

mseq = genmsequence(n,R0,T0);
% mseq = -mseq;
mseq(mseq == -1) = 0; % rescale to [0,1];

M = genmsequenceMovie(n,c,r,mseq);
makeCortexMovie(M,framespertrial,machineframespermovieframe,machinepixelspermoviepixel,filename)
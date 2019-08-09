function makeCortexMovie(M,framespertrial,machineframespermovieframe,machinepixelspermoviepixel,filename)

numdigits = 5;

frames = size(M,3);
if size(M,1)*machinepixelspermoviepixel ~= round(size(M,1)*machinepixelspermoviepixel) ||...
        size(M,2)*machinepixelspermoviepixel ~= round(size(M,2)*machinepixelspermoviepixel)
    error('machinepixelspermoviepixel must be integar')
end

Numtrials = ceil(frames/framespertrial);
framesinlasttrial = mod(frames,framespertrial);
if framesinlasttrial == 0
    framesinlasttrial = framespertrial;
end

for i = 1:(Numtrials-1)
    Mthistrial = M(:,:,((i-1)*framespertrial+1):i*framespertrial);
    MachineM = zeros(size(M,1)*machinepixelspermoviepixel,size(M,2)*machinepixelspermoviepixel,framespertrial);
    for j = 1:framespertrial
        MachineM(:,:,j) = imresize(Mthistrial(:,:,j),machinepixelspermoviepixel,'nearest');
    end
    MachineM = round(MachineM * 127) + 128;
    MachineM = uint8(MachineM);
    numzeros = numdigits - (floor(log10(i))+1);
    trialfilename = [filename,repmat('0',1,numzeros),num2str(i),'.ctx'];
    savecx_movie(trialfilename,MachineM,machineframespermovieframe);
end
% the last trial
Mthistrial = M(:,:,(end-framesinlasttrial+1):end);
MachineM = zeros(size(M,1)*machinepixelspermoviepixel,size(M,2)*machinepixelspermoviepixel,framesinlasttrial);
for j = 1:framesinlasttrial
    MachineM(:,:,j) = imresize(Mthistrial(:,:,j),machinepixelspermoviepixel,'nearest');
end
MachineM = round(MachineM * 127) + 128;
MachineM = uint8(MachineM);
numzeros = numdigits - (floor(log10(Numtrials))+1);
trialfilename = [filename,repmat('0',1,numzeros),num2str(Numtrials),'.ctx'];
savecx_movie(trialfilename,MachineM,machineframespermovieframe);
function savecx_movie(filename,M,machineframepermovieframe,notes)
switch nargin
    case 3
        notes = '';
    case 4
        % do nothing
    otherwise
end

movieframes = size(M,3);
machineframes = movieframes * machineframepermovieframe;

MachineM = zeros(size(M,1),size(M,2),machineframes);
MachineM = uint8(MachineM);
for thismovieframe = 1:movieframes
    for thismachineframe = 1:machineframepermovieframe
        MachineM(:,:,((thismovieframe-1)*machineframepermovieframe+thismachineframe)) = ...
            M(:,:,thismovieframe);
    end
end

%       dmns=[depth, x, y, nframes], in which
%               depth,          bitmap depth (1,2,4, or 8)
%               x,              x dimension of the image
%               y,              y dimension of the image
%               nframes,        number of frames in the movie
dmns = [8, size(MachineM,1),size(MachineM,2),machineframes];
savecx_movie_firstframe(filename, notes, dmns, MachineM(:,:,1))
for i = 2:machineframes
    savecx_movie_succframe(filename, notes, dmns, MachineM(:,:,i))
end

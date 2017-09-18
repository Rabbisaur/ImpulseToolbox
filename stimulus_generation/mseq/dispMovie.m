function dispMovie(M,framespertrial,t1,t2)
close all
frames = size(M,3);
trials = floor(frames/framespertrial);
figure
blankscreen = zeros(size(M,1),size(M,2));
for j = 1:trials
    imagesc(blankscreen);
    pause(t2)
    for i = 1:framespertrial
        imagesc(M(:,:,(j-1)*framespertrial+i))
        pause(t1)
    end
    imagesc(blankscreen);
    pause(t2)
end
function M = genmsequenceMovie(n,c,r,mseq)
frames = length(mseq);
p = 2^n/r/c;
M = zeros(r,c,frames);
noiselength = length(mseq);

for i = 0:c-1
    for j = 0:r-1
        for k = 1:frames
            idx = k+p*i + p*c*j;
            idx = mod(idx,noiselength);
            if idx == 0;
                idx = noiselength;
            end
            M(i+1,j+1,k) = mseq(idx);
        end
    end
end
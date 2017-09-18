clear,clc
n = 26;
seq.n = n;
seqV = zeros(2^n-1,1);
seqR0 = zeros(2^n-1,1);
seqT = zeros(2^n-1,1);

parfor T0 = 1:(2^n-1)
    R0 = 1;
    b = genmsequence(n,R0,T0);
    seqV(T0)= sum(b) == -1;
    seqT(T0)=T0;
    seqR0(T0) = R0;
    disp(['T = ', num2str(T0),'; V = ', num2str(seqV(T0))])
end

seq.L = seqV;
seq.T = seqT;
seq.R0 = seqR0;
save(['./seq',num2str(n), '.mat'],'seq')
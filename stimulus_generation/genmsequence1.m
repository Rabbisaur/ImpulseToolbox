function b = genmsequence(n,R0,T0)

R = zeros(1,n);
R1=R0;
for j = 1:n
    R(j) = floor(R1 / 2^(n-j));
    R1 = R1 - R(j)* 2^(n-j);
end
T = zeros(1,n);
T1=T0;
for j = 1:n
    T(j) = floor(T1 / 2^(n-j));
    T1 = T1 - T(j)* 2^(n-j);
end
R = logical(R);
T = logical(T);
% b = R(end);
% t = 0;
% i = 1;
% while(1)
%     i = i+1;
%     P = mod(sum(R & T),2);
%     R = circshift(R,[0,1]);
%     R(1) = P;
%     b(i) = R(end);
%     %     disp(R)
%     if mod(i-n,2) == 0 && i-n > 2
%         t = issame(b(n+(1:(i-n)/2)), b(n+(i-n)/2 +1:(i)));
%     end
%     if t > 0
%         L = (i-n)/2;
%         b = b(n+1:n+1+L-1);
%         
%         break;
%     end
% end
for i = 1:n
    P = mod(sum(R & T),2);
    R = circshift(R,[0,1]);
    R(1) = P;
end
b = zeros(2^n-1,1);
for i = 1:2^n-1
    P = mod(sum(R & T),2);
    R = circshift(R,[0,1]);
    R(1) = P;
    b(i) = R(end);
end
b = double(b);
b = 1 - 2 * b;
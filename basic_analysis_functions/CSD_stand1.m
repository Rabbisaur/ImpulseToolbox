function [CSD] = CSD_stand(mLFP, h, bVaknin, bFilter)
% mLFP is the LFP averaged across trials, (continuous data in electrode x times form)  
% bVaknin can be 0 or 1, to determine whether to abort the the first and the
% last electrode
% bFilter can be 0 or 1, if bFilter is 1, filter the result with 3 point
% h is inter contact distance in meters

if nargin<3
    bVaknin=1; bFilter=0;
elseif nargin<4
    bFilter=0;
end
% if exist([filepath 'mLFP.mat'], 'file')~=2
%     error([filepath 'mLFP does NOT exist, run getCSD first']);
% end
% load([filepath 'mLFP']);
% load([filepath 'tLFP']);

% h=100e-6;
cond = 0.300;

%%%%if the file is mLFP,  average LFP first, then compute the CSD
[m1, m2]= size(mLFP);
if bVaknin
  LFP=zeros(m1+2,m2);
  LFP(1,:)=mLFP(1,:);
  LFP(2:m1+1,:)=mLFP;
  LFP(m1+2,:)=mLFP(m1,:);
else
  LFP=mLFP;
end
CSD = -cond*D1(m1+2,h)*LFP*1e-6;
%%%%%% if file is tLFP, compute the CSD first, then average together
% [n1,n2,n3]=size(tLFP);
% tCSD=zeros(n1,n2,n3);
% for i=1:n1
%     if bVaknin
%       LFP_t=zeros(n2+2,n3);
%       LFP_t(1,:)=squeeze(tLFP(i,1,:));
%       LFP_t(2:n2+1,:)=squeeze(tLFP(i,:,:));
%       LFP_t(n2+2,:)=squeeze(tLFP(i,n2,:));
%     else
%       LFP_t = squeeze(tLFP(i,:,:));
%     end
%     tCSD(i,:,:) = -cond*D1(n2+2,h)*LFP_t*1e-6;
% end
% CSD_t=squeeze(mean(tCSD));
%%% filter the original CSD matrix if needed 
b0=0.54;
b1=0.23;
[n1,n2]=size(CSD); 
if bFilter %filter iCSD (does not change size of CSD matrix)         
  CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
  CSD_add(n1+2,:)=zeros(1,n2);
  CSD_add(2:n1+1,:)=CSD;        %CSD_add has n1+2 rows
  CSD = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows
end

% if bFilter %filter iCSD (does not change size of CSD matrix)         
%   CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
%   CSD_add(n1+2,:)=zeros(1,n2);
%   CSD_add(2:n1+1,:)=CSD_t;        %CSD_add has n1+2 rows
%   CSD_t = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows
% end

%%%%% substract the baseline of before stimulus onset
CSD_base= mean(CSD(:,1:100),2);
% CSD_base_t=mean(CSD_t(:,1:100),2);
for j=1:n2
   CSD(:,j)=CSD(:,j)-CSD_base; 
%    CSD_t(:,j)=CSD_t(:,j)-CSD_base_t;
end

 function out = D1(N,h)
%The matrix form of the standard double derivative formula, called D1 in
%Freeman and Nicholson (1975).
%N: number of electrodes
%h: inter-contact distance.
%out is a (N-2)x(N) matrix.
if nargin < 2, h = 100e-6 ;end;
if nargin < 1, N = 16 ;end;

for i=1:N-2
    for j=1:N
        if (i == j-1)
            out(i,j) = -2/h^2;
        elseif (abs(i-j+1) == 1)
            out(i,j) = 1/h^2;
        else
            out(i,j) = 0;
        end;
    end;
end;

function out = S_general(N,b0,b1)
%S = S_general(N,b0,b1)
%This is the three point filter matrix.
%Returns matrix of size (N-2)x(N),
%which represents a three point "spatial noise" filter with mid
%coeffescient b0 (will be normalized by the function) and neighbouring
%coeffescients b1 (will also be normalized).
%Default filter has b0 = 2 and b1 = 1 and number_of_electrodes = 18.
if nargin < 1, N = 18; end
if nargin < 3, b1 = 1; b0 = 2; end;
c = b0 + 2*b1;
out = zeros(N-2,N);
for i=1:N-2
    for j=1:N
        if (i == j-1)
            out(i,j) = b0/c;
        elseif (abs(i-j+1) == 1)
            out(i,j) = b1/c;
        end;
    end;
end;


function [Output, AlignedID] = Aligning(Input, Method)
%Align waveforms
% Usage:
%   [Output AlignedID] = Aligning(Input, Method)
% Inputs:
%   Input      a waveform in each row
%   Method     alignment methods, 1/2/3-none/global negative/postive
% Outputs:
%   Output	   interpolated waveforms
%   AlignedID  index of waveforms adjusted
% Copyright 2012-2012 Minggui Chen, Beijing Normal University
% Revision: 1.0 Date: 2012/2/9 23:59:59

%% Check INPUT
if ~isnumeric(Input) || size(Input, 2)<2
    error('INPUT must be a numeric matrix with 2/more coloumns.');
end
switch nargin
    case 1
        Method = 1;
    case 2
    otherwise
        error('ALIGNING accepts 1/2 input arguments.');
end
%% Align the input
warning('off');
[Row, Col] = size(Input);
Output = nan([Row Col], 'single');
switch Method
    case 1  %no further alignment
        Output = Input;
        AlignedID = false([size(Output, 1) 1]);
        return;
    case 2  %global Minimum
        [~, RefID] = min(Input, [], 2);
    case 3  %global Maximum
        [~, RefID] = max(Input, [], 2);
    case 4  %Global Absolute Maximum
        [~, RefID] = max(abs(Input), [], 2);
    case 5  %global Minimum slope
        SmoothInput = Smoothing(Input, 'uniform', -round(Col/40):round(Col/40));
        Slope = diff(SmoothInput, 1, 2);
        [~, Index] = sort(Slope, 2, 'ascend');
        RefID = round(median(Index(:, 1:max(5, round(Col*0.05))), 2));
    case 6  %Global Maximum Slope
        SmoothInput = Smoothing(Input, 'uniform', -round(Col/40):round(Col/40));
        Slope = diff(SmoothInput, 1, 2);
        [~, Index] = sort(Slope, 2, 'descend');
        RefID = round(median(Index(:, 1:max(5, round(Col*0.05))), 2));
    case 7  %Global Absolute Slope
        SmoothInput = Smoothing(Input, 'uniform', -round(Col/40):round(Col/40));
        Slope = diff(SmoothInput, 1, 2);
        [~, Index] = sort(abs(Slope), 2, 'descend');
        RefID = round(median(Index(:, 1:max(5, round(Col*0.05))), 2));
    case 8  %Multi-peak Energy Comparison(Chan, Neurocomputing2010)
        [~, MaxABSIndex] = max(abs(Input), [], 2);
        RefID = nan([Row 1]);
        Energy = Input.^2;
        Energy(:,1) = Energy(:,1)-Input(:,1+1).^2;
        Energy(:,2:end-1) = Energy(:,2:end-1)-Input(:,1:end-2).*Input(:,3:end);
        Energy(:,end) = Energy(:,end)-Input(:,end-1).^2;
        Energy(Energy<=0) = 0;
        Energy = sqrt(Energy);
        Slope = diff(Input, 1, 2);
        for i = 1:Row
            if Input(i,MaxABSIndex(i))<=mean(Input(i,:))
                PreIndex = find(Slope(i,1:MaxABSIndex(i)-1)>0, 1, 'last')+1;
                PostIndex = find(Slope(i,MaxABSIndex(i)+1:end)<0, 1, 'first');
            else
                PreIndex = find(Slope(i,1:MaxABSIndex(i)-1)<0, 1, 'last')+1;
                PostIndex = find(Slope(i,MaxABSIndex(i)+1:end)>0, 1, 'first');
            end
            PostIndex = PostIndex+MaxABSIndex(i);
            if isempty(PreIndex)
                PreIndex = 1;
            end
            if isempty(PostIndex)
                PostIndex = Col;
            end
            if Energy(i,PreIndex)>=Energy(i,PostIndex) && ...
                    Energy(i,MaxABSIndex(i))-Energy(i,PreIndex)<=Energy(i,PreIndex)-Energy(i,PostIndex)              
                RefID(i) = PreIndex;
            else
                RefID(i) = MaxABSIndex(i);
            end
        end
    otherwise
        error('No current METHOD was found.');
end
IDTable = tabulate(RefID);
[~, MaxID] = max(IDTable(:,3));
AlignID = IDTable(MaxID,1);
for i = 1:Row
    CurrShift = AlignID-RefID(i);
    Output(i,max(1,1+CurrShift):min(Col,Col+CurrShift)) = Input(i,max(1,1-CurrShift):min(Col,Col-CurrShift));
end
if Method~=1
    AlignedID = logical(sum(isnan(Output), 2));
    Output = InterpNAN(Output);
end



function Output = InterpNAN(Output)
[Row, Col] = size(Output);
for i = [1 Col]
    Output(isnan(Output(:,i)),i) = median(Output(~isnan(Output(:,i))));
end
for i = 1:Row
    NANID = isnan(Output(i,:));
    if sum(NANID)>0
        Output(i,:) = interp1(find(~NANID), Output(i,~NANID), 1:Col, 'cubic');
    end
end
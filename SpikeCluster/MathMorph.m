function [Output] = MathMorph(Input, Element, Operation)
%Apply mathmatical morphology operation
% Usage:
%   [Output] = MathMorph(Input, Element, Operation)
% Inputs:
%   Input       m by n, an input signal in each row
%   Element     structuring element less than n columns
%   Operation   morphlogical operations, erosion/dilation/opening/closing
% Outputs:
%   Output	    morphologically analyzed output
% Copyright 2012-2012 Minggui Chen, Beijing Normal University
% Revision: 1.0 Date: 2012/2/16 23:59:59

%% Check inputs
assert(nargin==3, 'MMORPHOLOGY accepts 3 input arguments.');
if ~isnumeric(Input) || numel(size(Input))>2
    error('INPUT must be a numeric matrix less than 2 dimensions.');
end
if ~isnumeric(Element) || numel(size(Element))>2
    error('ELEMENT must be a numeric matrix less than 2 dimensions.');
end
if size(Input, 2)==1
    Input = Input';
end
if size(Element, 2)==1
    Element = Element';
end
[InputRow InputCol] = size(Input);
[EleRow EleCol] = size(Element);
if EleRow==1
    Element = ones([InputRow 1])*Element;
else
    assert(EleRow==InputRow, 'ELEMENT and INPUT must have equal number of rows.');
end
assert(EleCol<=InputCol, 'The number of columns in ELEMENT must be small than INPUT.');
assert(ischar(Operation), 'OPERATION must be a string.');
%% Morphological processing
switch lower(Operation)
    case 'erosion'
        Output = MErosion(Input, Element);
    case 'dilation'
        Output = MDilation(Input, Element);
    case 'opening'
        Output = MDilation(MErosion(Input, Element), Element);
    case 'closing'
        Output = MErosion(MDilation(Input, Element), Element);
    otherwise
        error('No current operation was found.');
end

%% Erosion
function Output = MErosion(Input, Element)
[InputRow InputCol] = size(Input);
[~, EleCol] = size(Element);
Output = nan([InputRow InputCol-EleCol+1], 'single');
for i = 1:InputCol-EleCol+1
    LoIndex = max(1, 1-i);
    HiIndex = min(EleCol, InputCol-i);
    Index = LoIndex:HiIndex;
    Output(:,i) = min(Input(:,i+Index)-Element(:,Index), [], 2);
end
%% Dilation
function Output = MDilation(Input, Element)
[InputRow InputCol] = size(Input);
[~, EleCol] = size(Element);
Output = nan([InputRow InputCol+EleCol-1], 'single');
for i = 0:InputCol+EleCol-2
    LoIndex = max(0, i+1-InputCol);
    HiIndex = min(EleCol-1, i);
    Index = LoIndex:HiIndex;
    Output(:,i+1) = max(Input(:,i-Index+1)+Element(:,Index+1), [], 2);
end
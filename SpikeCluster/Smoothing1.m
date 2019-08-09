function output = Smoothing(input,method,varargin)
% Smoothing each row with the specified method
%   Usage: output = Smoothing(input,method,varargin)
%     
% Inputs
%   input - a vector/2-D matrix, smoothing each row
%   method - Gaussian/Uniform
%   varargin
%     'uniform':
%        win - total kernel window vector
%     'gaussian':
%        sigma - the std of gauusian kernel
%        win - total kernel window vector
%
% External Includes:
%
% Outputs:
%   output - the weighted average of all points within WIN
%
% History:
%   create: 5/31/2011, MGChen@BNU
%   flexible window: 7/9/2011, MGChen@BNU

%% check input vector
input = single(input);
if isempty(input)
    output = input;
    return;
end
if ~isreal(input) && ~isempty(input) || ~ismatrix(input)
    error('INPUT must be a 2-D matrix.');
end
bTransposed = 0;
if isvector(input) && size(input,1)>1
    bTransposed = 1;
    input = input';
end

%% window
inputlen = size(input, 2);
switch lower(method)
    case {'uniform';'uni'}
        if nargin ~= 2+1
            error('Smoothing(UNIFORM) accepts 1 variable argument only.');
        end
        win = varargin{1};  
        rawkernel = win.^0;
    case {'gaussian';'gauss'}
        if nargin ~= 2+2
            error('Smoothing(GAUSSIAN) accepts 2 variable arguments only.');
        end
        sigma = varargin{1};
        win = varargin{2};
        if isequal(win,0)
            rawkernel = 1;
        else            
            rawkernel = exp(-win.^2/sigma^2/2);
        end
        %disp(['Gauss smooth wil result in an epiphenomena of osccilation.']);
    otherwise
        error('No such method was found.');
end
winlen = max(abs(win));
kernel = zeros(1,2*winlen+1);
kernel(win+winlen+1) = rawkernel;
kernel = kernel/sum(kernel);  %set the mean to be equal
if sum(~isnan(kernel)) == 0
    error('KERNEL has been set to be all zero-weighted.');
end

padded(:, (winlen+1):(winlen+inputlen)) = input;
padded(:, 1:winlen) = input(:,1)*ones(1,winlen);
padded(:, (length(padded)+1):(length(padded)+winlen)) = input(:,inputlen)*ones(1,winlen);

output = conv2(padded, fliplr(kernel));  %numel(padded)+numel(kernel)-1
output = output(:, (2*winlen+1):(2*winlen+inputlen));

if bTransposed
    output = output';
end
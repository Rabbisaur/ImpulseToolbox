function hfile = hdr_InhomogeneityCorrect(hfile, opts)
% HDR::InhomogeneityCorrect  - attempt automatic inhomogeneity correction
%
% FORMAT:       [hdr = ] hdr.InhomogeneityCorrect([opts])
%
% Input fields:
%
%       opts        optional struct with settings
%        .mask      either 3D uint8/logical data or HDR object with preseg
%                   if omitted, try automatic mask detection
%        .model     either of 'log', {'mult'}
%        .numpasses number of passes, default 3 (valid: 1 through 5)
%        .order     polynomial order, default 3 (valid: 2 through 7)
%        .xmask     use mask in conjunction with autodetected mask
%
% Using: pmbfilter.

% Version:  v0.9d
% Build:    14082616
% Date:     Aug-26 2014, 4:43 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2011, 2012, 2014, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% neuroelf library
global ne_methods;

% argument check
if nargin < 1 || ...
    numel(hfile) ~= 1 || ...
   ~xffisobject(hfile, true, 'hdr')
    error( ...
        'xff:BadArgument', ...
        'Invalid call to ''%s''.', ...
        mfilename ...
    );
end
bc = xffgetcont(hfile.L);
if any(bc.ImgDim.DataType == [128, 1536, 2048, 2304])
    error( ...
        'xff:BadArgument', ...
        'Not defined for complex datatypes.' ...
    );
end
if ndims(bc.VoxelData) ~= 3
    error( ...
        'xff:BadArgument', ...
        'Invalid call to ''%s''.', ...
        mfilename ...
    );
end
if nargin < 2 || ...
    numel(opts) ~= 1 || ...
   ~isstruct(opts)
    opts = struct;
end
if ~isfield(opts, 'mask') || ...
    isempty(opts.mask)
    opts.mask = [];
end
if numel(opts.mask) == 1 && ...
    xffisobject(opts.mask, true, 'hdr')
    mbc = xffgetcont(opts.mask.L);
    opts.mask = mbc.VoxelData(:, :, :);
    if ~isequal(size(opts.mask), size(bc.VoxelData))
        opts.mask = [];
    end
    if ~isa(opts.mask, 'uint8')
        opts.mask = uint8([]);
    end
    if ~isempty(opts.mask)
        opts.mask(opts.mask < 226) = 0;
        opts.mask(opts.mask > 225) = opts.mask(opts.mask > 225) - 225;
        um = unique(opts.mask(:) + 1);
        ur = uint8(1:max(um));
        ur(uo) = 1:numel(um);
        opts.mask = ur(opts.mask);
    end
end
if ~isfield(opts, 'model') || ...
   ~ischar(opts.model) || ...
   ~any(strcmpi(opts.model(:)', {'l', 'log', 'm', 'mult'}))
    opts.model = 'mult';
else
    opts.model = lower(opts.model(1));
    if opts.model == 'l'
        opts.model = 'log';
    else
        opts.model = 'mult';
    end
end
if ~isfield(opts, 'numpasses') || ...
    numel(opts.numpasses) ~= 1 || ...
   ~isa(opts.numpasses, 'double') || ...
    isnan(opts.numpasses) || ...
   ~any((1:5) == opts.numpasses)
    opts.numpasses = 3;
end
if ~isfield(opts, 'order') || ...
    numel(opts.order) ~= 1 || ...
   ~isa(opts.order, 'double') || ...
    isnan(opts.order) || ...
   ~any((2:7) == opts.order)
    opts.order = 3;
end
if ~isfield(opts, 'xmask') || ...
    numel(opts.xmask) ~= 1 || ...
   ~islogical(opts.xmask)
    opts.xmask = false;
end

% get data
vd = aft_GetVolume(hfile, 1);

% apply correction (pre-filter)
for pc = 1:(opts.numpasses-1)
    vd = ne_methods.pmbfilter(vd, opts.order, opts.mask, struct('xmask', opts.xmask));
end

% apply final pass
vd = ne_methods.pmbfilter(vd, opts.order, opts.mask, struct( ...
    'bcutoff', 0.1, ...
    'cmask',   true, ...
    'robust',  true, ...
    'xmask',   opts.xmask));
vd(vd < 0) = 0;

% set back
bc.ImgDim.DataType = 16;
bc.ImgDim.BitsPerPixel = 32;
bc.ImgDim.ScalingSlope = 1;
bc.ImgDim.ScalingIntercept = 0;
bc.VoxelData = single(vd);

% set in output
xffsetcont(hfile.L, bc);

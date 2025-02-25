function hfile = bless(hfile, lev)
% xff::bless  - make object non-unwindable
%
% FORMAT:       [obj = ] bless(obj [, lev]);
%
% Input fields:
%
%       obj         valid xff object
%       lev         optional levels, if not given up to base WS
%
% Output fields:
%
%       obj         copy of input object (for convenience)

% Version:  v0.9d
% Build:    14061918
% Date:     Jun-19 2014, 6:28 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010, 2014, Jochen Weber
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

% global storage
global xffclup;
global xffcont;

% argument check
if nargin < 1 || ...
    any(~xffisobject(hfile, true))
    return;
end
if nargin < 2 || ...
    isa(lev, 'double') || ...
    numel(lev) ~= 1 || ...
    isnan(lev) || ...
    lev < 1
    lev = Inf;
else
    lev = fix(lev);
end

% iterate over object
sfile = struct(hfile);
for sc = 1:numel(sfile)

    % remove stack entirely
    ifile = find(xffclup == sfile(sc).L);
    if isinf(lev)
        xffcont(ifile).U = {};

    % remove part of stack
    else
        U = xffcont(ifile).U;
        if ~isempty(U)
            U(1:min(numel(U), lev)) = [];
            xffcont(ifile).U = U;
        end
    end
end

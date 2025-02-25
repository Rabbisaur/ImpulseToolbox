function s = subsref(p, S)
% xprogress::subsref  - support for .Progress calls
%
% FORMAT:       pbar.Progress(...);
%
% Input fields:
%
%       pbar        xprogress object

% Version:  v0.9c
% Build:    11050319
% Date:     May-02 2011, 6:17 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, Jochen Weber
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

% argument check
if nargin < 2 || ...
    numel(p) ~= 1 || ...
   ~isstruct(S) || ...
   ~isfield(S, 'type') || ...
   ~isfield(S, 'subs') || ...
    isempty(S) || ...
   ~strcmp(S(1).type, '.')
    error( ...
        'xprogress:BadCall', ...
        'Invalid call to xprogress::subsref.' ...
    );
end

% preset output
s = [];

% allow cell with one char call type
ssubs = S(1).subs(:)';
if iscell(ssubs) && ...
    numel(ssubs) == 1
    ssubs = ssubs{1};
end

% disallow non-char addressing
if ~ischar(ssubs) || ...
    isempty(ssubs)
    error( ...
        'xprogress:BadCall', ...
        'Invalid call to xprogress:subsref.' ...
    );
end

% what call
switch (lower(ssubs(:)'))
    case {'close'}
        xprogress(p, 'close', []);
    case {'progress'}
        if numel(S) ~= 2 || ...
           ~strcmp(S(2).type, '()')
            error( ...
                'xprogress:BadCall', ...
                'Invalid .Progress call to object.' ...
            );
        end
        xprogress(p, S(2).subs{:});
    otherwise
        try
            sp = struct(p);
            s = subsref(sp, S(2:end));
        catch ne_eo;
            rethrow(ne_eo);
        end
end

% FUNCTION ne_cm_setstype: set statistics type
function ne_cm_setstype(varargin)

% Version:  v1.0
% Build:    14091811
% Date:     Sep-18 2014, 11:45 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2014, Jochen Weber
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

% global variable
global ne_gcfg;
hFig = ne_gcfg.h.CM.CMFig;

% nothing to do without useful input
if nargin < 3 || ...
   ~ischar(varargin{3}) || ...
    isempty(varargin{3}) || ...
   ~any(strcmpi(varargin{3}(:)', {'ols', 'olsrob', 'rob'}))
    return;
end

% update UI settings
switch lower(varargin{3}(:)')
    case {'ols'}
        hFig.RadioGroupSetOne('Stats', 1);
        hFig.SetGroupEnabled('Robust', 'off');
    case {'olsrob'}
        hFig.RadioGroupSetOne('Stats', 3);
        hFig.SetGroupEnabled('Robust', 'on');
    case {'rob'}
        hFig.RadioGroupSetOne('Stats', 2);
        hFig.SetGroupEnabled('Robust', 'on');
end

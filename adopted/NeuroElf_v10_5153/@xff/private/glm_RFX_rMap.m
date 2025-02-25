function map = glm_RFX_rMap(hfile, c, r, mapopts)
% GLM::RFX_rMap  - calculate a second-level r contrast map
%
% FORMAT:       map = glm.RFX_rMap(c, r [, mapopts])
%
% Input fields:
%
%       c           NxC contrast vector
%       r           SxR regression data (or VOI object)
%       mapopts     structure with optional fields
%        .allrs     boolean flag, treat multiple Rs as one model (false)
%        .bbox      bounding box for VMP if necessary (default: full MNI)
%        .bvcomp    BV-compatible map names (length restriction, true)
%        .cnames    1xC cell array with contrast names
%        .const     also create (t-) map of constant term
%        .estfwhm   estimate smoothness and store in Map.RunTimeVars (true)
%        .groups    Gx2 cell array with group names and subject IDs
%        .imeth     interpolation method if necessary (default: 'cubic')
%        .meanr     boolean flag, remove mean from map (added as cov)
%        .meanrmsk  mask to get mean from (object or XxYxZ logical)
%        .minnum    minumum number of subjects to compute (2 * sqrt(N))
%        .names     1xN cell array with map names
%        .rnames    1xR cell array with regressor names
%        .rank      flag, rank-transform data before regression
%        .robust    flag, use robust regression in addition to OLS
%        .robwmaps  create robust-weight summary maps (default: false)
%        .smk       smoothing kernel for data (prior to regression, 0)
%        .subsel    subject selection (otherwise all subjects)
%        .swmaps    weight map per subject and regression (default: false)
%        .thresh    1x2 threshold (lower, upper), as p-values!
%        .voiidx    index into VOI list (only used if r is a VOI object)
%
% Output fields:
%
%       map         MAP/VMP/SMP object with maps
%
% Using: correlinvtstat, findfirst, fitrobustbisquare_img, lsqueeze,
%        newnatresvmp, ranktrans, resestsmooth, robustt, sdist, smoothdata3,
%        ztrans.

% Version:  v1.0
% Build:    14100210
% Date:     Oct-02 2014, 10:45 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% Copyright (c) 2010, 2011, 2014, Jochen Weber
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
correlinvtstat = ne_methods.correlinvtstat;
findfirst      = ne_methods.findfirst;
limitrangec    = ne_methods.limitrangec;
ranktrans      = ne_methods.ranktrans;
sdist          = ne_methods.sdist;

% argument check
if nargin < 3 || ...
    numel(hfile) ~= 1 || ...
   ~xffisobject(hfile, true, 'glm') || ...
   ~isa(c, 'double') || ...
    isempty(c) || ...
    any(isnan(c(:)) | isinf(c(:))) || ...
   ((~isa(r, 'double') || ...
     isempty(r) || ...
     any(isinf(r(:)))) && ...
    (numel(r) ~= 1 || ...
    ~xffisobject(r, true, 'voi')))
    error( ...
        'xff:BadArgument', ...
        'Invalid call to %s.', ...
        mfilename ...
    );
end
sbc = xffgetscont(hfile.L);
bc = sbc.C;
bcrtv = bc.RunTimeVars;
glmfile = sbc.F;
glmid = sbc.C.RunTimeVars.xffID;
if isempty(glmfile)
    glmfile = glmid;
end
if bc.ProjectTypeRFX == 0 && ...
    bc.SeparatePredictors ~= 2
    error( ...
        'xff:BadArgument', ...
        'Invalid call to %s.', ...
        mfilename ...
    );
end
isrfx = (bc.ProjectTypeRFX > 0);
ffxspred = glm_SubjectPredictors(hfile);
ffxsubs = glm_Subjects(hfile);
if isrfx
    numsubs = numel(bc.GLMData.Subject);
    numspred = size(bc.GLMData.Subject(1).BetaMaps, ...
        ndims(bc.GLMData.Subject(1).BetaMaps));
else
    ffxpred = bc.Predictor;
    ffxpred = {ffxpred(:).Name2};
    ffxpred = ffxpred(:);
    numsubs = numel(ffxsubs);
    numspred = numel(ffxspred) + 1;
end
if ~any(bc.ProjectType == [1, 2])
    error( ...
        'xff:Unsupported', ...
        'RFX correlation maps of FMRs are not yet supported.' ...
    );
end
if numsubs < 3
    error( ...
        'xff:BadArgument', ...
        'Invalid RFX GLM object.' ...
    );
end
if bc.ProjectType == 1
    if isrfx
        msz = size(bc.GLMData.RFXGlobalMap);
    else
        msz = size(bc.GLMData.MCorrSS);
    end
else
    if isrfx
        msz = numel(bc.GLMData.RFXGlobalMap);
    else
        msz = numel(bc.GLMData.MCorrSS);
    end
end
if nargin < 4 || ...
   ~isstruct(mapopts) || ...
    numel(mapopts) ~= 1
    mapopts = struct;
end
if ~isfield(mapopts, 'allrs') || ...
   ~islogical(mapopts.allrs) || ...
    numel(mapopts.allrs) ~= 1
    mapopts.allrs = false;
end
if ~isfield(mapopts, 'bbox') || ...
   ~isa(mapopts.bbox, 'double') || ...
   ~isequal(size(mapopts.bbox), [2, 3]) || ...
    any(isinf(mapopts.bbox(:)) | isnan(mapopts.bbox(:)) | ...
        mapopts.bbox(:) < 0 | mapopts.bbox(:) > 256 | mapopts.bbox(:) ~= fix(mapopts.bbox(:))) || ...
    any(diff(mapopts.bbox) < 0)
    mapopts.bbox = [44, 38, 44; 242, 194, 212];
end
if ~isfield(mapopts, 'bvcomp') || ...
   ~islogical(mapopts.bvcomp) || ...
    numel(mapopts.bvcomp) ~= 1
    mapopts.bvcomp = true;
end
if ~isfield(mapopts, 'cnames') || ...
   ~iscell(mapopts.cnames)
    mapopts.cnames = {};
end
if ~isfield(mapopts, 'const') || ...
   ~islogical(mapopts.const) || ...
    numel(mapopts.const) ~= 1
    mapopts.const = false;
end
if ~isfield(mapopts, 'estfwhm') || ...
   ~islogical(mapopts.estfwhm) || ...
    numel(mapopts.estfwhm) ~= 1
    mapopts.estfwhm = true;
end
if ~isfield(mapopts, 'groups') || ...
   ~iscell(mapopts.groups) || ...
    size(mapopts.groups, 2) ~= 2
    ngrp = 0;
    if ~isfield(mapopts, 'subsel') || ...
       ~isa(mapopts.subsel, 'double') || ...
        isempty(mapopts.subsel) || ...
        any(isinf(mapopts.subsel(:)) | isnan(mapopts.subsel(:)))
        ga = ones(numsubs, 1);
    else
        mapopts.subsel = mapopts.subsel(:)';
        mapopts.subsel(mapopts.subsel < 1 | mapopts.subsel > numsubs) = [];
        ga = zeros(numsubs, 1);
        ga(unique(round(mapopts.subsel))) = 1;
    end
    mapopts.groups = [];
    gamx = 1;
else
    ga = zeros(numsubs, 1);
    for gc = 1:size(mapopts.groups, 1)
        if ~ischar(mapopts.groups{gc, 1}) || ...
            isempty(mapopts.groups{gc, 1}) || ...
            isempty(mapopts.groups{gc, 2}) || ...
            any(isinf(mapopts.groups{gc, 2}(:)) | isnan(mapopts.groups{gc, 2}(:))) || ...
            any(mapopts.groups{gc, 2}(:) < 1 | mapopts.groups{gc, 2}(:) > numsubs) || ...
            any(ga(round(mapopts.groups{gc, 2}(:))) > 0)
            error( ...
                'xff:BadArgument', ...
                'Invalid group assignment.' ...
            );
        end
        mapopts.groups{gc, 2} = unique(round(mapopts.groups{gc, 2}(:)));
        ga(mapopts.groups{gc, 2}) = gc;
    end
    ngrp = size(mapopts.groups, 1);
    nval = sum(ga > 0) - ngrp;
    if nval < 3
        error( ...
            'xff:BadArgument', ...
            'Too few subjects for grouping.' ...
        );
    end
    gas = ga(ga > 0);
    gag = tril(ones(ngrp));
    gag(1:(ngrp + 1):end) = 0;
    gamx = sum(gag(:));
    gag(gag > 0) = 1:gamx;
end
gax = find(ga > 0);
if ~isfield(mapopts, 'imeth') || ...
   ~ischar(mapopts.imeth) || ...
    isempty(mapopts.imeth) || ...
    isempty(regexpi(mapopts.imeth(:)', '^(cubic|lanczos\d|linear|nearest)$'))
    mapopts.imeth = 'cubic';
else
    mapopts.imeth = lower(mapopts.imeth(:)');
end
if ~isfield(mapopts, 'meanr') || ...
   ~islogical(mapopts.meanr) || ...
    numel(mapopts.meanr) ~= 1
    mapopts.meanr = false;
end
if isfield(mapopts, 'meanrmsk') && ...
    numel(mapopts.meanrmsk) == 1 && ...
    xffisobject(mapopts.meanrmsk, true, 'msk')
    mbc = xffgetcont(mapopts.meanrmsk.L);
    if numel(mbc.Mask) == prod(msz)
        mapopts.meanrmsk = ne_methods.lsqueeze(mbc.Mask > 0);
    else
        mapopts.meanrmsk = [];
    end
elseif isfield(mapopts, 'meanrmsk') && ...
    islogical(mapopts.meanrmsk) && ...
    numel(mapopts.meanrmsk) == prod(msz)
    mapopts.meanrmsk = ne_methods.lsqueeze(mapopts.meanrmsk);
else
    mapopts.meanrmsk = [];
end
if isempty(mapopts.meanrmsk) && ...
    mapopts.meanr
    mapopts.meanrmsk = all(bc.GLMData.Subject(1).BetaMaps ~= 0, ...
        ndims(bc.GLMData.Subject(1).BetaMaps));
    for sc = 1:numsubs
        mapopts.meanrmsk = (mapopts.meanrmsk & ...
            all(bc.GLMData.Subject(1).BetaMaps ~= 0, ...
            ndims(bc.GLMData.Subject(1).BetaMaps)));
    end
    mapopts.meanrmsk = ne_methods.lsqueeze(mapopts.meanrmsk);
end
if ~isfield(mapopts, 'minnum') || ...
   ~isa(mapopts.minnum, 'double') || ...
    numel(mapopts.minnum) ~= 1 || ...
    isinf(mapopts.minnum) || ...
    isnan(mapopts.minnum) || ...
    mapopts.minnum < 0
    mapopts.minnum = 0;
end
if ~isfield(mapopts, 'names') || ...
   ~iscell(mapopts.names) || ...
    isempty(mapopts.names)
    mapopts.names = {};
end
if ~isfield(mapopts, 'rank') || ...
   ~islogical(mapopts.rank) || ...
    numel(mapopts.rank) ~= 1 || ...
   ~mapopts.rank
    mapopts.rank = false;
    ranktxt = '';
else
    ranktxt = 'Rank-';
end
if ~isfield(mapopts, 'rnames') || ...
   ~iscell(mapopts.rnames)
    mapopts.rnames = {};
end
if ~isfield(mapopts, 'robust') || ...
   ~islogical(mapopts.robust) || ...
    numel(mapopts.robust) ~= 1
    mapopts.robust = false;
end
if ~isfield(mapopts, 'robwmaps') || ...
   ~islogical(mapopts.robwmaps) || ...
    numel(mapopts.robwmaps) ~= 1
    mapopts.robwmaps = false;
end
if ~isfield(mapopts, 'smk') || ...
    numel(mapopts.smk) ~= 1 || ...
   ~isa(mapopts.smk, 'double') || ...
    isinf(mapopts.smk) || ...
    isnan(mapopts.smk) || ...
    mapopts.smk <= (0.5 * bc.Resolution)
    mapopts.smk = 0;
else
    mapopts.smk = min(mapopts.smk, 8 * bc.Resolution);
end
if ~isfield(mapopts, 'subsel') || ...
   ~isa(mapopts.subsel, 'double') || ...
    isempty(mapopts.subsel) || ...
    any(isinf(mapopts.subsel(:)) | isnan(mapopts.subsel(:))) || ...
    numel(unique(round(mapopts.subsel(:)))) ~= numel(mapopts.subsel) || ...
    any(mapopts.subsel(:) < 1 | mapopts.subsel(:) > numsubs)
    mapopts.subsel = 1:numsubs;
else
    mapopts.subsel = round(mapopts.subsel(:)');
end
if ~isfield(mapopts, 'swmaps') || ...
   ~islogical(mapopts.swmaps) || ...
    numel(mapopts.swmaps) ~= 1
    mapopts.swmaps = false;
end
if ~isfield(mapopts, 'voiidx') || ...
   ~isa(mapopts.voiidx, 'double') || ...
    numel(mapopts.voiidx) ~= 1 || ...
    isinf(mapopts.voiidx) || ...
    isnan(mapopts.voiidx) || ...
    mapopts.voiidx < 1
    mapopts.voiidx = 1;
else
    mapopts.voiidx = floor(mapopts.voiidx);
end
if numel(r) == 1
    try
        rbc = xffgetcont(r.L);
        rbc = numel(rbc.VOI);
        if size(c, 2) ~= numspred && ...
            size(c, 2) ~= (numspred - 1)
            ct = c';
        else
            ct = c;
        end
        r = glm_VOIBetas(hfile, r, ...
            struct('c', ct, 'vl', min(mapopts.voiidx, rbc)));
    catch ne_eo;
        rethrow(ne_eo);
    end
end
if any(size(r) == numsubs) && ...
    numel(mapopts.subsel) ~= numsubs
    rsr = repmat({':'}, 1, ndims(r));
    rsr{findfirst(size(r) == numsubs)} = mapopts.subsel;
    r = r(rsr{:});
end
if ~any(size(r) == numel(mapopts.subsel))
    error( ...
        'xff:BadArgument', ...
        'Correlation regressors must match in size with number of subjects.' ...
    );
end
rsdim = findfirst(size(r) == numel(mapopts.subsel));
nanr = isnan(r);
for dc = 1:ndims(r)
    if dc == rsdim
        continue;
    end
    nanr = any(nanr, dc);
end
if any(nanr)
    rsr = repmat({':'}, 1, ndims(r));
    rsr{rsdim} = find(~nanr);
    r = r(rsr{:});
    mapopts.subsel = mapopts.subsel(~nanr);
end
if ~isfield(mapopts, 'thresh') || ...
   ~isa(mapopts.thresh, 'double') || ...
    numel(mapopts.thresh) ~= 2 || ...
    any(isinf(mapopts.thresh) | isnan(mapopts.thresh) | mapopts.thresh <= 0 | mapopts.thresh >= 0.5)
    mapopts.thresh = [0.005, 0.0001];
else
    mapopts.thresh = -sort(-mapopts.thresh(:)');
end
subsel = mapopts.subsel;
numsubs = numel(subsel);
if mapopts.minnum == 0
    mapopts.minnum = ceil(2 * sqrt(numsubs));
elseif mapopts.minnum < 1
    mapopts.minnum = ceil(mapopts.minnum * numsubs);
end
mapopts.minnum = min(numsubs, mapopts.minnum);
if mapopts.meanr
    mvm = zeros(numsubs, 1);
end
thresh = mapopts.thresh;
if size(r, 2) == numel(mapopts.subsel) && ...
    size(r, 1) ~= numel(mapopts.subsel)
    r = r';
end
if size(c, 1) == 1 && ...
   (size(c, 2) == (numspred - 1) || ...
    size(c, 2) == numspred)
    c = c';
end
if size(c, 1) == (numspred - 1)
    c(end+1,:) = 0;
end
if size(c, 1) ~= numspred
    error( ...
        'xff:BadArgument', ...
        'Contrast vector must span all conditions.' ...
    );
end
nummaps = size(c, 2);
numrs = size(r, 2);
if mapopts.allrs
    nval = numsubs - (1 + numrs);
else
    nval = numsubs - 2;
end
if numel(mapopts.cnames) ~= nummaps
    sprednames = glm_SubjectPredictors(hfile);
    mapopts.cnames = cell(1, nummaps);
    for cc = 1:nummaps
        con = find(c(:, cc) > 0);
        coff = find(c(:, cc) < 0);
        conn = cell(1, numel(con));
        coffn = cell(1, numel(coff));
        for occ = 1:numel(con)
            if c(con(occ), cc) ~= 1
                conn{occ} = sprintf('%g * %s', sprednames{con(occ)});
            else
                conn{occ} = sprednames{con(occ)};
            end
        end
        for occ = 1:numel(coff)
            if c(coff(occ), cc) ~= -1
                coffn{occ} = sprintf('%g * %s', sprednames{coff(occ)});
            else
                coffn{occ} = sprednames{coff(occ)};
            end
        end
        if ~isempty(con)
            if numel(con) > 1
                connc = sprintf('%s + ', conn{:});
                connc = sprintf('(%s)', connc(1:end-3));
            else
                connc = conn{1};
            end
        else
            connc = 'Baseline';
        end
        if ~isempty(coff)
            if numel(coff) > 1
                coffnc = sprintf('%s + ', coffn{:});
                coffnc = sprintf('(%s)', coffnc(1:end-3));
            else
                coffnc = coffn{1};
            end
        else
            coffnc = '';
        end
        if ~isempty(coffnc)
            mapopts.cnames{cc} = sprintf('%s > %s', connc, coffnc);
        else
            mapopts.cnames{cc} = connc;
        end
    end
end
if numel(mapopts.rnames) ~= numrs
    mapopts.rnames = cell(1, numrs);
    for cc = 1:numrs
        mapopts.rnames = sprintf('reg%02d', cc);
    end
end
nummapst = nummaps * (numrs + mapopts.const) * ...
    (1 + mapopts.meanr) * (1 + mapopts.robust) * ...
    (1 + double(mapopts.robust) .* ...
    (2 * double(mapopts.robwmaps) + numsubs * double(mapopts.swmaps)));
if bc.ProjectType == 1
    subsa = {':', ':', ':', []};
    subsr = {':', ':', ':', []};
else
    subsa = {':', []};
    subsr = {':', []};
end
rpma = [msz, 1];
if numel(mapopts.names) ~= (nummaps * numrs)
    mapopts.names = cell(1, nummaps * numrs);
end
for cc = 1:numel(mapopts.names)
    if ~ischar(mapopts.names{cc})
        ccr = 1 + mod(cc - 1, numrs);
        ccc = 1 + round((cc - ccr) / numrs);
        mapopts.names{cc} = sprintf('%sCorr: (%s, %s)', ranktxt, ...
            mapopts.cnames{ccc}, mapopts.rnames{ccr});
    end
    if mapopts.bvcomp && ...
        numel(mapopts.names{cc}) > 96
        mapopts.names{cc} = ...
            [mapopts.names{cc}(1:46) ' ... ' mapopts.names{cc}(end-45:end)];
    else
        mapopts.names{cc} = mapopts.names{cc}(:)';
    end
end
if isrfx
    szmap = size(bc.GLMData.RFXGlobalMap);
else
    szmap = size(bc.GLMData.MCorrSS);
end

% create map container
copymaps = true;
switch (bc.ProjectType)

    % VTC/VMP
    case {1}
        map = bless(ne_methods.newnatresvmp(), 1);
        mapc = xffgetcont(map.L);
        mapc.Resolution = bc.Resolution;
        if ~isfield(bcrtv, 'SubjectSPMsn') || ...
           ~isstruct(bcrtv.SubjectSPMsn) || ...
            isempty(fieldnames(bcrtv.SubjectSPMsn))
            if isrfx
                szmap = size(bc.GLMData.RFXGlobalMap);
            else
                szmap = size(bc.GLMData.MCorrSS);
            end
            mapc.XStart = bc.XStart;
            mapc.XEnd = bc.XEnd;
            mapc.YStart = bc.YStart;
            mapc.YEnd = bc.YEnd;
            mapc.ZStart = bc.ZStart;
            mapc.ZEnd = bc.ZEnd;
        else
            mapopts.bbox(2, :) = mapopts.bbox(1, :) + mapc.Resolution .* ...
                ceil((mapopts.bbox(2, :) - mapopts.bbox(1, :)) ./ mapc.Resolution - 0.01);
            mapc.XStart = mapopts.bbox(1, 1);
            mapc.XEnd = mapopts.bbox(2, 1);
            mapc.YStart = mapopts.bbox(1, 2);
            mapc.YEnd = mapopts.bbox(2, 2);
            mapc.ZStart = mapopts.bbox(1, 3);
            mapc.ZEnd = mapopts.bbox(2, 3);
            szmap = round((mapopts.bbox(2, :) - mapopts.bbox(1, :)) ./ mapc.Resolution);
            copymaps = false;
            sbbox = struct('BBox', mapopts.bbox, 'ResXYZ', mapc.Resolution);
            rpma = round(diff(mapopts.bbox) ./ mapc.Resolution);
            msz = rpma;
        end
        mapc.RunTimeVars.AutoSave = true;
        mapc.RunTimeVars.TrfPlus = bcrtv.TrfPlus;
        mapf = 'VMPData';

    % MTC/SMP
    case {2}
        map = bless(xff('new:smp'), 1);
        mapc = xffgetcont(map.L);
        mapc.NrOfVertices = bc.NrOfVertices;
        mapf = 'SMPData';
end

% set some common fields
mapc.NrOfMaps = nummapst;
mapc.Map.Type = 2;
mapc.Map.LowerThreshold = ...
    correlinvtstat(-sdist('tinv', thresh(1), nval), numsubs);
mapc.Map.UpperThreshold = ...
    correlinvtstat(-sdist('tinv', thresh(2), nval), numsubs);
mapc.Map.DF1 = nval;
mapc.Map.DF2 = 0;
mapc.Map.NrOfFDRThresholds = 0;
mapc.Map.FDRThresholds = zeros(0, 3);
mapc.Map.(mapf) = single(zeros(rpma));

% replicate
mapc.Map = mapc.Map(1, ones(1, nummapst));

% rank-transform
if mapopts.rank
    r = ranktrans(r, 1);
end

% what models?
if mapopts.allrs
    micc = 1:numrs:(nummaps * numrs);
    numrsi = numrs;
else
    micc = 1:(nummaps * numrs);
    numrsi = 1;
end

% computation
conmaps = zeros([msz, numsubs]);
tmc = 1;
lcc = 0;
for icc = 1:numel(micc)

    % which contrast and regressor
    cr = 1 + mod(micc(icc) - 1, numrs);
    cc = 1 + round((micc(icc) - cr) / numrs);

    % zero out and fill conmaps if necessary
    if lcc ~= cc
        conmaps(:) = 0;
        keepsubs = true(1, numsubs);

        % fill contrast maps
        for pc = 1:numspred
            if c(pc, cc) ~= 0
                subsr{end} = pc;
                subsrs = struct('type', '()', 'subs', {subsr});
                for sc = 1:numsubs
                    if ~keepsubs(sc)
                        continue;
                    end
                    subsa{end} = sc;
                    subsas = struct('type', '()', 'subs', {subsa});
                    if isrfx
                        if copymaps
                            conmaps = subsasgn(conmaps, subsas, ...
                                subsref(conmaps, subsas) + c(pc, cc) .* ...
                                subsref(bc.GLMData.Subject(subsel(sc)).BetaMaps, subsrs));
                        else
                            conmaps = subsasgn(conmaps, subsas, ...
                                subsref(conmaps, subsas) + c(pc, cc) .* ...
                                aft_SampleBVBox(hfile, sbbox, ...
                                    (sc - 1) * numspred + pc, mapopts.imeth));
                        end
                    else
                        keepsubi = findfirst(~cellfun('isempty', regexpi(ffxpred, ...
                            sprintf('^subject\\s+%s:\\s*%s', ...
                            ffxsubs{gax(sc)}, ffxspred{pc}))));
                        if ~isempty(keepsubi)
                            conmaps = subsasgn(conmaps, subsas, ...
                                subsref(conmaps, subsas) + c(pc, cc) .* ...
                                subsref(bc.GLMData.Subject(subsel(sc)).BetaMaps, subsrs));
                        else
                            conmaps = subsasgn(conmaps, subsas, 0);
                            keepsubs(sc) = false;
                        end
                    end
                end
            end
        end

        % check maps
        for sc = 1:numsubs
            if ~keepsubs(sc)
                continue;
            end
            subsa{end} = sc;
            subsas = struct('type', '()', 'subs', {subsa});
            cmtest = ne_methods.lsqueeze(subsref(conmaps, subsas));
            if all(isnan(cmtest) | isinf(cmtest) | cmtest == 0)
                conmaps = subsasgn(conmaps, subsas, 0);
                keepsubs(sc) = false;
            end
        end
        keepsubsn = sum(keepsubs);
        if mapopts.allrs
            nval = keepsubsn - (1 + numrs);
        else
            nval = keepsubsn - 2;
        end

        % minumum criterion not met? then set to 0 for all subjects
        conmapsa = (~isinf(conmaps));
        conmapsa = conmapsa & (~isnan(conmaps));
        conmapsa = conmapsa & (conmaps ~= 0);
        conmapsa = sum(conmapsa, numel(subsr));
        conmaps(repmat(conmapsa < mapopts.minnum, [ones(1, numel(msz)), numsubs])) = 0;

        % rank transform?
        if mapopts.rank
            subsa{end} = keepsubs;
            subsas = struct('type', '()', 'subs', {subsa});
            conmaps = subsasgn(conmaps, subsas, ...
                ranktrans(subsref(conmaps, subsas), numel(subsr), ...
                struct('meancenter', true, 'nozero', true)));
        else
            conmaps(isinf(conmaps) | isnan(conmaps)) = 0;
        end

        % keep track!
        lcc = cc;
    end

    % generate design matrix/ces
    if mapopts.allrs
        X = [ne_methods.ztrans(r(keepsubs, :)), ones(keepsubsn, 1)];
    else
        X = [ne_methods.ztrans(r(keepsubs, cr)), ones(keepsubsn, 1)];
    end
    iXX = pinv(X' * X);
    if mapopts.meanr
        for sc = 1:keepsubsn
            subsa{end} = sc;
            subsas = struct('type', '()', 'subs', {subsa});
            mv = subsref(conmaps, subsas);
            mvm(sc) = mean(mv(mapopts.meanrmsk));
        end
        if mapopts.rank
            Xm = [X, ne_methods.ztrans(ranktrans(mvm(keepsubs), 1))];
        else
            Xm = [X, ne_methods.ztrans(mvm(keepsubs))];
        end
        iXXm = pinv(Xm' * Xm);
    end

    % smooth data
    if mapopts.smk > 0 && ...
        ndims(conmaps) == 4
        resim = (isinf(conmaps) | isnan(conmaps) | conmaps == 0);
        conmaps = ne_methods.smoothdata3(conmaps, mapopts.smk / mapc.Resolution);
        conmaps(resim) = 0;
    end

    % set additional data
    artv = struct( ...
        'SourceGLM',   glmfile, ...
        'SourceGLMID', glmid, ...
        'Contrast',    c(:, lcc), ...
        'FWHMResEst',  [], ...
        'FWHMResImg',  [], ...
        'GlobSigMap',  [], ...
        'MeanRem',     mapopts.meanr, ...
        'Regressors',  r, ...
        'RFXGLM',      true, ...
        'Robust',      mapopts.robust, ...
        'SubPreds',    {ffxspred}, ...
        'SubSel',      subsel(:));

    % OLS computations first
    subsa{end} = keepsubs;
    subsas = struct('type', '()', 'subs', {subsa});
    conmaps = subsref(conmaps, subsas);
    betas = iXX * X' * reshape(conmaps, prod(msz), keepsubsn)';
    resim = conmaps - reshape((X * betas)', [msz, keepsubsn]);
    stder = sqrt(sum(resim .^ 2, ndims(conmaps))) .* sqrt(1 / nval);
    if mapopts.estfwhm
        [artv.FWHMResEst, artv.FWHMResImg] = ...
            ne_methods.resestsmooth(resim, bc.Resolution);
    end

    % first maps
    for irc = 1:numrsi
        tmap = betas(irc, :)' ./ (stder(:) .* sqrt(iXX(irc, irc)));
        tmap(isinf(tmap) | isnan(tmap)) = 0;
        mapc.Map(tmc).(mapf) = reshape(single(correlinvtstat( ...
            tmap, keepsubsn)), [msz, 1]);
        mapc.Map(tmc).Name = mapopts.names{icc+irc-1};
        if mapopts.allrs
            artv.Regressors = r(:, irc);
        else
            artv.Regressors = r(:, cr);
        end
        mapc.Map(tmc).DF1 = nval;
        mapc.Map(tmc).RunTimeVars = artv;
        tmc = tmc + 1;
    end
    if mapopts.const
        tmap = betas(numrsi+1, :)' ./ (stder(:) .* sqrt(iXX(numrsi+1, numrsi+1)));
        tmap(isinf(tmap) | isnan(tmap)) = 0;
        mapc.Map(tmc).(mapf) = reshape(single(tmap), [msz, 1]);
        mapc.Map(tmc).Type = 1;
        mapc.Map(tmc).Name = sprintf('%s (intercept-t)', mapopts.names{icc});
        mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval);
        mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval);
        artv.Regressors = [];
        mapc.Map(tmc).DF1 = nval;
        mapc.Map(tmc).RunTimeVars = artv;
        tmc = tmc + 1;
    end

    % with mean removed ?
    if mapopts.meanr
        betas = iXXm * Xm' * reshape(conmaps, prod(msz), keepsubsn)';
        resim = conmaps - reshape((Xm * betas)', [msz, keepsubsn]);
        stder = sqrt(sum(resim .^ 2, ndims(conmaps))) .* sqrt(1 / (nval - 1));
        if mapopts.estfwhm
            [artv.FWHMResEst, artv.FWHMResImg] = ...
                ne_methods.resestsmooth(resim, bc.Resolution);
        end
        artv.GlobSigMap = single(reshape(mapopts.meanrmsk, rpma));

        % first maps
        for irc = 1:numrsi
            tmap = betas(irc, :)' ./ (stder(:) .* sqrt(iXXm(irc, irc)));
            tmap(isinf(tmap) | isnan(tmap)) = 0;
            mapc.Map(tmc).(mapf) = reshape(single(correlinvtstat( ...
                tmap, keepsubsn)), [msz, 1]);
            mapc.Map(tmc).Name = sprintf('%s (mean-rem)', mapopts.names{icc+irc-1});
            mapc.Map(tmc).DF1 = nval - 1;
            mapc.Map(tmc).LowerThreshold = ...
                correlinvtstat(-sdist('tinv', thresh(1), nval - 1), nval + 1);
            mapc.Map(tmc).UpperThreshold = ...
                correlinvtstat(-sdist('tinv', thresh(2), nval - 1), nval + 1);
            if mapopts.allrs
                artv.Regressors = r(:, irc);
            else
                artv.Regressors = r(:, cr);
            end
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if mapopts.const
            tmap = betas(numrsi+1, :)' ./ (stder(:) .* sqrt(iXXm(numrsi+1, numrsi+1)));
            tmap(isinf(tmap) | isnan(tmap)) = 0;
            mapc.Map(tmc).(mapf) = reshape(single(tmap), [msz, 1]);
            mapc.Map(tmc).Type = 1;
            mapc.Map(tmc).DF1 = nval - 1;
            mapc.Map(tmc).Name = sprintf('%s (intercept-t, mean-rem)', mapopts.names{icc});
            mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval - 1);
            mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval - 1);
            artv.Regressors = [];
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        artv.GlobSigMap = [];
    end

    % compute robust stats
    if mapopts.robust

        % perform fit
        [b, w] = ne_methods.fitrobustbisquare_img(X, conmaps);
        if mapopts.estfwhm
            ptc = zeros(size(conmaps));
            for bmc = 1:size(X, 2)
                ptc = ptc + repmat(b(:, :, :, bmc), [1, 1, 1, size(X, 1)]) .* ...
                    repmat(reshape(X(:, bmc), [1, 1, 1, size(X, 1)]), szmap);
            end
            ptc = w .* ptc + (1 - w) .* conmaps;
            [artv.FWHMResEst, artv.FWHMResImg] = ...
                ne_methods.resestsmooth(conmaps - ptc, bc.Resolution);
        end
        rsicv = zeros(1, size(X, 2));
        for irc = 1:numrsi
            rsicv(:) = 0;
            rsicv(irc) = 1;
            rt = ne_methods.robustt(X, conmaps, b, w, rsicv);
            rt(isinf(rt) | isnan(rt)) = 0;
            corm = correlinvtstat(rt, keepsubsn);
            corm(isinf(corm) | isnan(corm)) = 0;
            mapc.Map(tmc).(mapf) = single(corm);
            mapc.Map(tmc).Name = sprintf('%s (robust)', mapopts.names{icc+irc-1});
            if mapopts.allrs
                artv.Regressors = r(:, irc);
            else
                artv.Regressors = r(:, cr);
            end
            mapc.Map(tmc).DF1 = nval;
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if mapopts.const
            rsicv(:) = 0;
            rsicv(numrsi+1) = 1;
            rt = ne_methods.robustt(X, conmaps, b, w, rsicv);
            rt(isinf(rt) | isnan(rt)) = 0;
            mapc.Map(tmc).Type = 1;
            mapc.Map(tmc).(mapf) = single(rt);
            mapc.Map(tmc).Name = sprintf('%s (robust, intercept-t)', mapopts.names{icc});
            mapc.Map(tmc).LowerThreshold = -sdist('tinv', thresh(1), nval);
            mapc.Map(tmc).UpperThreshold = -sdist('tinv', thresh(2), nval);
            artv.Regressors = [];
            mapc.Map(tmc).DF1 = nval;
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if isfield(artv, 'FWHMResEst')
            artv.FWHMResEst = [];
            artv.FWHMResImg = [];
        end
        if mapopts.robwmaps
            mapc.Map(tmc).Type = 145;
            mapc.Map(tmc).Name = sprintf('%s (mean robust weight)', mapopts.names{icc});
            mapc.Map(tmc).LowerThreshold = 0.75;
            mapc.Map(tmc).UpperThreshold = 1;
            mapc.Map(tmc).(mapf) = single((1 / size(w, ndims(w))) .* sum(w, ndims(w)));
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
            mapc.Map(tmc).Name = sprintf('%s (weighted number of outliers)', mapopts.names{icc});
            mapc.Map(tmc).Type = 146;
            mapc.Map(tmc).LowerThreshold = 1;
            mapc.Map(tmc).UpperThreshold = 0.5 * size(w, ndims(w));
            mapc.Map(tmc).(mapf) = single(-3 .* sum( ...
                limitrangec(w - 2/3, -1/3, 0, -1/6), ndims(w)));
            mapc.Map(tmc).RunTimeVars = artv;
            tmc = tmc + 1;
        end
        if mapopts.swmaps
            for irc = find(keepsubs(:)')
                mapc.Map(tmc).Name = sprintf('%s (%s outlier)', ...
                    mapopts.names{icc}, ffxsubs{subsel(irc)});
                mapc.Map(tmc).Type = 144;
                mapc.Map(tmc).LowerThreshold = 0.25;
                mapc.Map(tmc).UpperThreshold = 1;
                mapc.Map(tmc).(mapf) = single(1 - w(:, :, :, irc));
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
        end
        if mapopts.meanr
            [b, w] = ne_methods.fitrobustbisquare_img(Xm, conmaps);
            artv.GlobSigMap = single(reshape(mapopts.meanrmsk, rpma));
            if mapopts.estfwhm
                ptc = zeros(size(conmaps));
                for bmc = 1:size(Xm, 2)
                    ptc = ptc + repmat(b(:, :, :, bmc), [1, 1, 1, size(Xm, 1)]) .* ...
                        repmat(reshape(Xm(:, bmc), [1, 1, 1, size(Xm, 1)]), szmap);
                end
                ptc = w .* ptc + (1 - w) .* conmaps;
                [artv.FWHMResEst, artv.FWHMResImg] = ...
                    ne_methods.resestsmooth(conmaps - ptc, bc.Resolution);
            end
            rsicv = zeros(1, size(Xm, 2));
            for irc = 1:numrsi
                rsicv(:) = 0;
                rsicv(irc) = 1;
                rt = ne_methods.robustt(Xm, conmaps, b, w, rsicv);
                rt(isinf(rt) | isnan(rt)) = 0;
                corm = correlinvtstat(rt, keepsubsn);
                corm(isinf(corm) | isnan(corm)) = 0;
                mapc.Map(tmc).(mapf) = single(corm);
                mapc.Map(tmc).Name = sprintf('%s (robust, mean-rem)', mapopts.names{icc+irc-1});
                mapc.Map(tmc).DF1 = nval - 1;
                mapc.Map(tmc).LowerThreshold = ...
                    correlinvtstat(-sdist('tinv', thresh(1), nval - 1), nval + 1);
                mapc.Map(tmc).UpperThreshold = ...
                    correlinvtstat(-sdist('tinv', thresh(2), nval - 1), nval + 1);
                if mapopts.allrs
                    artv.Regressors = r(:, irc);
                else
                    artv.Regressors = r(:, cr);
                end
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
            if mapopts.const
                rsicv(:) = 0;
                rsicv(numrsi+1) = 1;
                rt = ne_methods.robustt(Xm, conmaps, b, w, rsicv);
                rt(isinf(rt) | isnan(rt)) = 0;
                mapc.Map(tmc).Type = 1;
                mapc.Map(tmc).(mapf) = single(rt);
                mapc.Map(tmc).Name = sprintf('%s (robust, intercept-t, mean-rem)', mapopts.names{icc});
                mapc.Map(tmc).DF1 = nval - 1;
                mapc.Map(tmc).LowerThreshold = ...
                    -sdist('tinv', thresh(1), nval - 1);
                mapc.Map(tmc).UpperThreshold = ...
                    -sdist('tinv', thresh(2), nval - 1);
                artv.Regressors = [];
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
            if mapopts.robwmaps
                mapc.Map(tmc).Type = 145;
                mapc.Map(tmc).Name = sprintf('%s (mean robust weight, w/o mean)', mapopts.names{icc});
                mapc.Map(tmc).LowerThreshold = 0.75;
                mapc.Map(tmc).UpperThreshold = 1;
                mapc.Map(tmc).(mapf) = single((1 / size(w, ndims(w))) .* sum(w, ndims(w)));
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
                mapc.Map(tmc).Name = sprintf('%s (weighted number of outliers, w/o mean)', mapopts.names{icc});
                mapc.Map(tmc).Type = 146;
                mapc.Map(tmc).LowerThreshold = 1;
                mapc.Map(tmc).UpperThreshold = 0.5 * size(w, ndims(w));
                mapc.Map(tmc).(mapf) = single(-3 .* sum( ...
                    limitrangec(w - 2/3, -1/3, 0, -1/6), ndims(w)));
                mapc.Map(tmc).RunTimeVars = artv;
                tmc = tmc + 1;
            end
            if mapopts.swmaps
                for irc = 1:find(keepsubs(:)')
                    mapc.Map(tmc).Name = sprintf('%s (%s outlier, w/o mean)', ...
                        mapopts.names{icc}, ffxsubs{subsel(irc)});
                    mapc.Map(tmc).Type = 144;
                    mapc.Map(tmc).LowerThreshold = 0.25;
                    mapc.Map(tmc).UpperThreshold = 1;
                    mapc.Map(tmc).(mapf) = single(1 - w(:, :, :, irc));
                    mapc.Map(tmc).RunTimeVars = artv;
                    tmc = tmc + 1;
                end
            end
            artv.GlobSigMap = [];
        end
    end
end
if tmc <= numel(mapc.Map)
    mapc.Map(tmc:end) = [];
end

% put back
xffsetcont(map.L, mapc);

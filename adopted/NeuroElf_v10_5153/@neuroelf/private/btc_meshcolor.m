% FUNCTION btc_meshcolor: colorize mesh accordingly
function btc_meshcolor(srf, recolor, source, target, varargin)

% Version:  v1.0
% Build:    15040309
% Date:     Apr-03 2015, 9:39 AM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/

% Copyright (c) 2010 - 2014, 2015, Jochen Weber
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

% if invalid SRF, return
if ~isxff(srf, 'srf')
    return;
end
rerender = (nargin < 5);

% get surface content
srfsc = getscont(srf);
srfc = srfsc.C;
srfh = srfsc.H;

% recolor default: off
if nargin < 2 || ...
   ~islogical(recolor) || ...
    numel(recolor) ~= 1
    recolor = false;
end
if nargin < 4 || ...
   ~ischar(source) || ...
   ~isfield(ne_gcfg.cc, source(:)') || ...
   ~isstruct(ne_gcfg.cc.(source(:)')) || ...
   ~isfield(ne_gcfg.cc.(source(:)'), 'Config') || ...
   ~isfield(ne_gcfg.cc.(source(:)').Config, 'sattype') || ...
   ~strcmpi(ne_gcfg.cc.(source(:)').Config.sattype, 'surf') || ...
    numel(target) ~= 1 || ...
  (~isa(target, 'double') && ...
   ~isa(target, 'matlab.graphics.primitive.Patch')) || ...
   ~ishandle(target) || ...
   ~strcmpi(get(target, 'Type'), 'patch') || ...
   (size(get(target, 'FaceVertexCData'), 1) ~= srfc.NrOfVertices && ...
    ~isempty(srfc.TriangleVertex))

    % don't do anything if more than 2 arguments!
    if nargin > 2
        return;
    end

    % get main window config
    source = '';
    fcfg = ne_gcfg.fcfg;

    % check whether any stats selected
    hsrf = srfh.Surface;
    if isfield(srfh, 'Stats')
        sths = srfh.Stats;
        stvar = sths{1};
        stvix = sths{2};
    else
        stvar = [];
        stvix = [];
    end
else
    fcfg = ne_gcfg.cc.(source(:)').Config;
    stvar = fcfg.SurfStatsVar;
    if ~isxff(stvar, {'cmp', 'glm', 'mtc', 'smp'})
        stvar = [];
        stvix = [];
    else
        stvix = fcfg.SurfStatsVarIdx;
    end
    hsrf = target;
end
midx = min(size(srfh.VertexMorphMeshes, 1), max(0, fcfg.srfcfg.time));

% get configuration
scfg = ne_gcfg.c.ini.Surface;
sbar = ne_gcfg.c.ini.Statistics.ShowThreshBars;

% get handles
ptio = srfh.SurfTIO;

% get color information (what color type and RGB info)
vcol = srfc.VertexColor(:, 1);
vrgb = srfc.VertexColor(:, 2:4);

% deleted vertices?
delv = false;
delc = (vcol(:, 1) >= 10512);
if any(delc)
    delv = true;
end

% recolor needed?
if numel(ptio.Layer) < 1
    recolor = true;
end

% recolor base
if recolor

    % get uint8 and alpha per-vertex vectors
    ctio = ptio.Rendered;
    atio = ones(size(ctio, 1), 1);

    % basic colors
    ccnv = [srfc.ConvexRGBA; srfc.ConcaveRGBA];

    % find vertices that need RGB coloring
    rgbc = find(isnan(vcol(:, 1)));

    % and already fill those with their RGB code
    ctio(rgbc, 1, :) = vrgb(rgbc, :);

    % which vertices get the "base color"
    basec = (~isnan(vcol) & (vcol < 2));

    % fill those also
    ctio(basec, 1, :) = round(255 .* ccnv(vcol(basec) + 1, 1:3));

    % for now the alpha setting of base colors is unheeded (as is in BV)
    % atio(basec) = ccnv(vcol(basec) + 1, 4);

    % set RGB vertices as if colored by base
    if ~isempty(rgbc)
        vcol(rgbc) = 0;
    end

    % find indexed (painted) colors and deleted
    idxc = find((vcol >= 10000) & (vcol <= 10511));

    % any indexed colors
    if ~isempty(idxc)

        % get RGB color table
        idxcol = [cat(3, scfg.IndexConvexColors{:}); cat(3, scfg.IndexConcaveColors{:})];

        % get base index
        idx = 0.5 .* (vcol(idxc) - 9998);

        % and which of the pair
        idxv = idx - floor(idx);
        idx = floor(idx);

        % fix to a maximum (of 11, as in BV standard configuration)
        idx = min(numel(scfg.IndexConcaveColors), idx);

        % then compute final index
        idxiv = 1 + double(idxv ~= 0) + (6 .* (idx - 1));

        % get RGB values for those indices (apart two indices each)
        ctio(idxc, 1, 1) = idxcol(idxiv);
        ctio(idxc, 1, 2) = idxcol(idxiv + 2);
        ctio(idxc, 1, 3) = idxcol(idxiv + 4);
    end

    % set to first layer
    atio(any(isnan(srfc.VertexCoordinate), 2)) = 0;
    setlayer(ptio, 1, ctio, atio);
end

% remove stats and other layers
if numel(ptio.Layer) > 1 && ...
    rerender
    dellayer(ptio, 2:numel(ptio.Layer));
end
srfh.SurfStatsBars = {zeros(256, 0, 3)};
sbari = 1;
sbarsz = fcfg.SurfBarSize;
sbarl = sbarsz(1) - 1;

% alpha blending overrider
mapca = zeros(srfc.NrOfVertices, 1);

% match ?
if ~isempty(stvix) && ...
    isxff(stvar, {'cmp', 'glm', 'mtc', 'smp'}) && ...
    srfc.NrOfVertices == stvar.NrOfVertices && ...
    rerender

    % for SMPs
    sttyp = stvar.FileType;
    if any(strcmp(sttyp, {'mtc', 'smp'}))

        % get maps
        maps = stvar.Map(stvix);
        strtv = stvar.RunTimeVars;

        % and iterate over maps
        for mc = 1:numel(maps)

            % for SMP
            if strcmp(sttyp, 'smp')
                % get data
                mapd = maps(mc).SMPData(:);

                % cluster ?
                if maps(mc).EnableClusterCheck && ...
                    isempty(maps(mc).SMPDataCT)
                    stvar.ClusterTable(srf, stvix(mc));
                    maps(mc) = stvar.Map(stvix(mc));
                end
                if maps(mc).EnableClusterCheck && ...
                   ~isempty(maps(mc).SMPDataCT)
                    mapd = mapd .* maps(mc).SMPDataCT;
                end

            % for MTC
            else

                % compute indices
                ncnds = strtv.NrOfConditions;
                ntcpc = strtv.NrOfTCsPerCondition;
                nvptc = strtv.NrOfVolumesPerTC;
                smvi = strtv.SubMapVol + (stvix(mc) - 1) * ntcpc * nvptc;
                cthr = strtv.ConditionThresholds(stvix(mc), 2, :);

                % interpolate map data
                mapic = [Inf, Inf; smvi, 1; ncnds * ntcpc * nvptc + 1, 1; size(stvar.MTCData)];
                if strtv.SubMapVol < 2 || ...
                    strtv.SubMapVol > (nvptc - 1)
                    simeth = 'linear';
                else
                    simeth = 'cubic';
                end
                if strtv.SubMapVol <= nvptc
                    mapd = flexinterpn_method(stvar.MTCData, mapic, simeth);
                else
                    mapd = zeros(1, stvar.NrOfVertices);
                end
                mapd = mapd(:);

                % also get error metric as alpha value
                mapic(2) = mapic(2) + nvptc;
                mapc = flexinterpn_method(stvar.MTCData, mapic, simeth);
                mapc = mapc(:);
                if strcmpi(strtv.TCNames{2}, 'sd')
                    mapc = (1 / sqrt(strtv.NrOfConditionOnsets(stvix(mc)))) .* mapc;
                end
                if lower(strtv.TCNames{2}(1)) == 's'
                    mapc = abs(mapd) ./ abs(mapc);
                    mapc(isnan(mapc)) = 0;
                end
                mapc = limitrangec((1 / max(cthr(2) - cthr(1), sqrt(eps))) .* ...
                    (mapc(:) - cthr(1)), 0, 1, 0);
            end

            % get thresholds
            ltr = maps(mc).LowerThreshold;
            utr = maps(mc).UpperThreshold;
            if numel(ltr) ~= 1 || ...
               ~isa(ltr, 'double') || ...
                isinf(ltr) || ...
                isnan(ltr) || ...
                ltr <= 0
                ltr = sqrt(eps);
                maps(mc).LowerThreshold = ltr;
            end
            if numel(utr) ~= 1 || ...
               ~isa(utr, 'double') || ...
                isinf(utr) || ...
                isnan(utr) || ...
                utr <= ltr
                utr = ltr + sqrt(eps);
                maps(mc).UpperThreshold = utr;
            end
            tls = maps(mc).ShowPositiveNegativeFlag;
            if numel(tls) ~= 1 || ...
               ~isa(tls, 'double') || ...
               ~any((0:3) == tls)
                tls = 3;
                maps(mc).ShowPositiveNegativeFlag = tls;
            end

            % threshold data
            mapd = threshmapc(double(mapd), ltr, utr, tls);

            % get colors
            lut = maps(mc).OverlayColors;
            if isempty(lut)
                lut = ne_gcfg.lut.Colors;
            end
            if sbar && ...
                mod(size(lut, 1), 2) == 0 && ...
                any(tls == [1, 2, 3])
                xlut = reshape(lut, [round(0.5 * size(lut, 1)), 2, 3]);
                xlut = flexinterpn(limitrangec((1 / 255) .* xlut, 0, 1, 0), ...
                    [Inf, Inf, Inf; 1, 1, 1; (size(xlut, 1) - 1) / sbarl, 1, 1; size(xlut)]);
                if tls == 1
                    srfh.SurfStatsBars{sbari} = xlut(:, 1, :);
                    sbari = sbari + 1;
                elseif tls == 2
                    srfh.SurfStatsBars{sbari} = xlut(:, 2, :);
                    sbari = sbari + 1;
                elseif tls == 3
                    srfh.SurfStatsBars{sbari} = xlut;
                    sbari = sbari + 1;
                end
            end

            % transparency setting
            if strcmp(sttyp, 'smp')
                alpha = maps(mc).TransColorFactor;
                if ~isa(alpha, 'double') || ...
                    numel(alpha) ~= 1 || ...
                    isinf(alpha) || ...
                    isnan(alpha) || ...
                    alpha > 1 || ...
                    alpha < -100
                    alpha = 1;
                    maps(mc).TransColorFactor = alpha;
                end

                % no transparency
                if alpha == 1

                    % simply show all values ~= 0
                    mapc = double(mapd ~= 0);

                % for positive alpha, set values ~= 0 to that value
                elseif alpha >= 0
                    mapc = alpha .* double(mapd ~= 0);

                % and for negative values
                else

                    % multiply with (absolute) stats value
                    mapc = min(1, abs(alpha .* mapd));
                end

            % restrict alpha to super-threshold data
            else
                mapc(mapd == 0) = 0;
            end

            % override alpha
            mapca = max(mapca, limitrangec(abs(mapd), 0, 1, 0));

            % set RGB as layer
            setlayer(ptio, 1 + mc, threshlutc(mapd, lut), mapc);
        end
    end
end

% join layers
if fcfg.join && ...
    numel(ptio.Layer) > 2

    % apply "stats-map join" (BV-like display) to the three images
    joinlayers(ptio, 2:numel(ptio.Layer));
end

% render colors
if rerender
    render(ptio);
end

% set coloring and limits
srfealpha = 0;
if midx == 0 || ...
   ~recolor
    set(hsrf, 'CData', (1 / 255) .* double(ptio.Rendered));
    if ~isempty(srfc.TriangleVertex)
        if ~isfield(srfc.RunTimeVars, 'VertexAlpha') || ...
           ~isa(srfc.RunTimeVars.VertexAlpha, 'double') || ...
            numel(srfc.RunTimeVars.VertexAlpha) ~= numel(mapca)
            srfgalpha = srfh.SurfProps{4};
        else
            srfgalpha = srfc.RunTimeVars.VertexAlpha;
        end
    else
        srfealpha = srfh.SurfProps{4};
        srfgalpha = 0;
    end
elseif midx == round(midx)
    set(hsrf, 'CData', (1 / 255) .* srfh.VertexMorphMeshes{midx, 4});
    if ~isempty(srfc.TriangleVertex)
        srfgalpha = srfh.VertexMorphMeshes{midx, 5}{4};
    else
        srfealpha = srfh.VertexMorphMeshes{midx, 5}{4};
        srfgalpha = 0;
    end
else
    bidx = floor(midx);
    tmrp = midx - bidx;
    bmrp = 1 - tmrp;
    if bidx == 0
        c1 = double(ptio.Rendered);
        if ~isfield(srfc.RunTimeVars, 'VertexAlpha') || ...
           ~isa(srfc.RunTimeVars.VertexAlpha, 'double') || ...
            numel(srfc.RunTimeVars.VertexAlpha) ~= numel(mapca)
            a1 = srfh.SurfProps{4};
        else
            a1 = srfc.RunTimeVars.VertexAlpha;
        end
    else
        c1 = srfh.VertexMorphMeshes{bidx, 4};
        a1 = srfh.VertexMorphMeshes{bidx, 5}{4};
    end
    c2 = srfh.VertexMorphMeshes{bidx+1, 4};
    a2 = srfh.VertexMorphMeshes{bidx+1, 5}{4};
    set(hsrf, 'CData', (1 / 255) .* (bmrp .* c1 + tmrp .* c2));
    if ~isempty(srfc.TriangleVertex)
        srfgalpha = bmrp .* a1 + tmrp .* a2;
    else
        srfealpha = srfh.SurfProps{4};
        srfgalpha = 0;
    end
end

% check rendered if needed
if delv && ...
   ~strcmpi(fcfg.renderer, 'opengl')
    disp('Warning: invisible faces require OpenGL rendering!');
    delv = false;
end

% visibility
if ~isa(srfgalpha, 'double')
    srfgalpha = 1;
end

% deleted vertices
if delv

    % we need to convert the triangles to a vertex -> triangle list!
    [vnei, bnei, vtri] = ...
        mesh_trianglestoneighbors(numel(delc), srfc.TriangleVertex);

    % and then create a triangle-based alpha data!
    if numel(srfgalpha) == 1
        ta = srfgalpha .* ones(size(srfc.TriangleVertex, 1), 1);
    else
        ta = mean(srfgalpha) .* ones(size(srfc.TriangleVertex, 1), 1);
    end

    % override
    ta = limitrangec(ta + srfgalpha .* mean(reshape( ...
        mapca(srfc.TriangleVertex(:), 1), size(srfc.TriangleVertex)), 2), 0, 1, 0);

    % and set deleted (with all three vertices deleted!) to 0
    ta(histcount(cat(2, vtri{delc}), 1, numel(ta), 1) > 0) = 0;

    % then set as FaceVertexAlphaData and FaceAlpha to flat
    set(hsrf, 'FaceVertexAlphaData', ta);
    set(hsrf, 'FaceAlpha', 'flat');
    if isempty(source)
        srfh.SurfProps{4} = 'flat';
        srfh.SurfProps{6} = ta;
    end

% otherwise
else

    % set to flat alpha
    if numel(srfgalpha) == 1 && ...
        all(mapca == 0)
        set(hsrf, 'FaceAlpha', srfgalpha);
        if isempty(source)
            srfh.SurfProps{6} = [];
        end
    else

        % mix correctly
        if numel(srfgalpha) == 1
            ta = srfgalpha .* ones(numel(mapca), 1);
        else
            ta = srfgalpha;
        end
        ta = limitrangec(ta + (1 - srfgalpha) .* mapca, 0, 1, 0);

        % then set
        set(hsrf, 'FaceVertexAlphaData', ta);
        set(hsrf, 'FaceAlpha', 'interp', 'AlphaDataMapping', 'none');
        if isempty(source)
            srfh.SurfProps{6} = ta;
        end
    end
end

% set markers
if lower(srfh.SurfProps{7}(1)) == 'n'
    set(hsrf, 'Marker', '.', 'MarkerFaceColor', 'none', 'MarkerSize', 1, ...
        'MarkerEdgeColor', 'none');
else
    set(hsrf, 'Marker', '.', 'MarkerFaceColor', 'none', 'MarkerSize', 1, ...
        'MarkerEdgeColor', 'flat');
end
if srfealpha > 0
    set(hsrf, 'EdgeAlpha', srfealpha);
end

% write back
srf.SetHandle('SurfProps', srfh.SurfProps);
srf.SetHandle('SurfStatsBars', srfh.SurfStatsBars);

function [main_rise_inds, is_sim, ind_bounds] = detect_main_voltage_rise(mean_waveform, par)
% detect_main_voltage_rise    

% Author: Ariel Tankus.
% Created: 21.08.2008.
% Copyright (C) Ariel Tankus, 2010.  All rights reserved.


% Minimal voltage rise between two successive sample points.  When this
% threshold is crossed, the main rise begins:
%min_voltage_rise = 1.5;
min_voltage_rise = 1.2;
eps_volt = 0.1;      % epsilon volts; this is used instead of 0, to avoid
                     % small fluctuations.
min_d2 = 0.3;        % Minimal value of second derivative, above which
                     % the rise starts.
min_d2_perc_from_max = 10;   % If min_d2 is never exceeded, use this
                             % threshold as percentage of maximal second
                             % derivative value.

max_mean_waveform = max(mean_waveform);
is_sim_waveform = (max_mean_waveform < 3);
if (is_sim_waveform)
    % Simulator files, reduce req. thresholds:
    % (simulator "voltage" is between [-1, 1]).
%    min_voltage_rise = min_voltage_rise ./ 10;
%    eps_volt = eps_volt ./ 10;
    mean_waveform = mean_waveform .* 100;
end
if (nargout >= 2)
    is_sim = is_sim_waveform;
end
is_sim = 0;  %minggui chen

d = diff(mean_waveform);      % derivative.
% The main voltage rise end at par.w_pre, because all spikes are aligned
% according to this time.
% Detect start of rise:  (+1: to get index of center of valley.)
%ind = find((d(1:(par.w_pre-1)) <= 0) & (d(2:par.w_pre) > 0), 1, 'last') + 1;

ind_v_rise = find(d > min_voltage_rise, 1, 'first') + 1;
if ((isempty(ind_v_rise)) || (ind_v_rise >= par.w_pre))
    % no main rise detected.  waveform suspected as not-a-spike.
    main_rise_inds = [NaN, par.w_pre];
    return;
end

% refine detection of start:
% Find the first rise higher than eps_volt:
ind = find((d(1:(ind_v_rise-1)) <= eps_volt) & ...
           (d(2:ind_v_rise) > eps_volt), 1, 'last') + 1;
if (isempty(ind))
    ind = 1;
end
% Find the first time when the second detivative starts rising:
dd = diff(d);
ind2 = find(dd(ind:(ind_v_rise-1)) >= min_d2, 1, 'first') + 1;   %+1: to take
                                                                 %central point.
if (isempty(ind2))
    % necessarily: for all i: dd(i) < min_d2.  Use a relative threshold,
    % w.r.t the max dd value:
    th = min_d2_perc_from_max ./ 100 * max(dd(ind:(ind_v_rise-1)));
    ind2 = find(dd(ind:(ind_v_rise-1)) >= th, 1, 'first') + 1;
end
ind2 = ind2 + ind - 1;  % make ind2 relative to beginning of spike.

% The extrinsic curvature in the 2D plane:
curv = dd ./ ((1+d(1:(end-1)).^2).^(3/2));
[max_c, max_ind] = max(curv(max(ind-1, 1):(ind_v_rise-1)));
max_ind = max_ind + ind - 1;    % make max_ind relative to beginning of spike.

%if (isempty(ind))
%    % necessarily, all values of d(1:ind_v_rise) > eps_volt, because
%    % d(ind_v_rise) > min_voltage_rise, and min_voltage_rise > eps_volt.
%    % In this case we begin from the minimum of d:
%    [min_val, ind] = min(d(1:ind_v_rise));
%    ind = ind(1);          % ensure only one minimum exists.
%end
%main_rise_inds = [ind, par.w_pre];

main_rise_inds = [max_ind, par.w_pre];
%main_rise_inds = [ind_v_rise, par.w_pre];

if (nargout >= 3)
    ind_bounds = [ind, ind_v_rise];
end
